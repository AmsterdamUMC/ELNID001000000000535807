# workflow/scripts/atac_deseq2.R
# Snakemake rule: atac_deseq2
# Runs DESeq2 differential accessibility analysis on merged ATAC-seq counts.
#
# REVISION NOTE
#   The fatigue models (continuous MFI, categorical fatigue group, HC
#   comparison) use the same code path as before. All new behavior is inside
#   `if (!is_fatigue_model)` blocks or inside the non-HC status branch.
#   New optional params (all NULL by default):
#     frip_adjust         TRUE/FALSE. Default TRUE. Non-fatigue models only.
#     adjust_for_group    TRUE/FALSE. Default FALSE. Non-HC status contrast only.
#     pool_types          TRUE/FALSE. Default FALSE. Within-type status contrast.

log_file <- snakemake@log[[1]]
con      <- file(log_file, open = "wt")
sink(con, type = "message")
sink(con, type = "output")

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(readr)
  library(DESeq2)
  library(apeglm)
  library(readxl)
  library(csaw)
})

source(snakemake@input[["utils"]])

# ── Parameters ────────────────────────────────────────────────────────────────
min_row_sums   <- snakemake@params$min_row_sums
deseq_coef     <- snakemake@params$deseq_coef
subset_type    <- snakemake@params$subset_type    # "CD", "UC", or NULL
subset_status  <- snakemake@params$subset_status  # "Active", "Remission", or NULL
padj_threshold <- snakemake@params$padj_threshold
cols_to_remove <- snakemake@params$cols_to_remove
continuous     <- isTRUE(snakemake@params$continuous)
status_contrast  <- isTRUE(snakemake@params$status_contrast)

hc_comparison    <- isTRUE(snakemake@params$hc_comparison)
hc_fatigue_level <- snakemake@params$hc_fatigue_level
if (is.null(hc_fatigue_level) || !nzchar(hc_fatigue_level)) hc_fatigue_level <- "Low"

## Fatigue models keep the old code path. Every new rule below applies only
## when this is FALSE.
is_fatigue_model <- continuous || hc_comparison ||
  (!status_contrast && grepl("Fatigue", deseq_coef))

## Read a TRUE/FALSE param. An absent param gives the default.
param_flag <- function(name, default) {
  v <- snakemake@params[[name]]
  if (is.null(v)) return(default)
  v <- as.logical(v)
  if (length(v) != 1 || is.na(v)) stop("param '", name, "' must be TRUE or FALSE")
  v
}

frip_adjust        <- param_flag("frip_adjust", TRUE)
adjust_for_group   <- param_flag("adjust_for_group", FALSE)
pool_types         <- param_flag("pool_types", FALSE)

has_subset_type <- !is.null(subset_type) && subset_type != "None"

## Pooled within-type status contrast: keep CD and UC together and estimate
## Active vs Remission inside each type with a Type:Status term. The dispersion
## is then estimated from all samples, not from the 12 samples of one type.
pool_types_active <- status_contrast && !hc_comparison && pool_types && has_subset_type

if (is_fatigue_model && (!frip_adjust || adjust_for_group || pool_types)) {
  warning("frip_adjust, adjust_for_group and pool_types ",
          "are ignored in fatigue models.")
}

# helper
## hc_fatigue_level: "Low" or "High". Selects which fatigue group pairs with
## No_Fatigue (HC) for this comparison. Ignored when hc_comparison is FALSE.
select_samples <- function(metadata, full_counts, subset_type, hc_comparison,
                           hc_fatigue_level = "Low") {
  if (hc_comparison) {
    keep_types <- if (!is.null(subset_type) && subset_type != "None") {
      c("HC", subset_type)
    } else {
      c("HC", "CD", "UC")
    }
    message("HC comparison: keeping Type in ", paste(keep_types, collapse = ", "))
    metadata <- metadata[metadata$Type %in% keep_types, ]
    
    fatigue_group <- paste0(hc_fatigue_level, "_Fatigue")
    message("HC comparison: keeping Status in {Healthy, Remission}, ",
            "group in {No_Fatigue, ", fatigue_group, "}")
    metadata <- metadata[metadata$Status %in% c("Healthy", "Remission") &
                           metadata$group  %in% c("No_Fatigue", fatigue_group), ]
    
  } else if (!is.null(subset_type) && subset_type != "None") {
    message("Subsetting to Type == ", subset_type)
    metadata <- metadata[metadata$Type == subset_type, ]
    
  } else {
    message("Removing Healthy Controls")
    metadata <- metadata[metadata$Type != "HC", ]
  }
  
  metadata$Status <- droplevels(as.factor(metadata$Status))
  metadata$group  <- droplevels(as.factor(metadata$group))
  full_counts <- full_counts[, rownames(metadata)]
  list(metadata = metadata, full_counts = full_counts)
}
# ── Load inputs ───────────────────────────────────────────────────────────────
message("Loading merged count table...")
full_counts <- readRDS(snakemake@input$counts_rds)

metadata <- read_excel(snakemake@input$metadata_file) %>%
  dplyr::filter(batch != "no_ATACseq", !Patient_ID %in% cols_to_remove) %>%
  dplyr::mutate(
    MFI_scaled  = as.numeric(scale(MFI_totalscore)),
    FRiP_scaled = as.numeric(scale(FRiP)),   # library-quality covariate
    sex         = as.factor(sex),
    group       = factor(group),
    Status      = factor(Status)
  ) %>%
  tibble::column_to_rownames("Patient_ID")

if (any(is.na(metadata$FRiP_scaled)))
  stop("NA in FRiP for: ",
       paste(rownames(metadata)[is.na(metadata$FRiP_scaled)], collapse = ", "))

print(rownames(metadata))
## In pooled mode, keep all IBD samples. The type is selected by the coefficient.
sel_type    <- if (pool_types_active) NULL else subset_type
sel         <- select_samples(metadata, full_counts, sel_type, hc_comparison, hc_fatigue_level)
metadata    <- sel$metadata
full_counts <- sel$full_counts

# ── Subset by disease status (Active / Remission) ────────────────────────────
# Not applied for status_contrast or HC comparisons — Status is the main
# factor there, or already restricted by select_samples().
stratified_by_status <- !is.null(subset_status) && subset_status != "None"
if (stratified_by_status && !status_contrast && !hc_comparison) {
  message("Subsetting to Status == ", subset_status)
  metadata    <- metadata[metadata$Status == subset_status, ]
  full_counts <- full_counts[, rownames(metadata)]
}

## ── Design formula ────────────────────────────────────────────────────────────
## Built as a term list, not as nested branches. Each covariate is added only
## if it can be estimated in the current subset. The function reports what it
## kept, what it dropped, and how many residual degrees of freedom remain.
##
## Rules (fatigue models):
##   * MFI_scaled or group or Status is the term of interest.
##   * FRiP_scaled is always a covariate: it explains count variance. It is not
##     a confounder of MFI (no correlation with MFI in this cohort), so it is
##     included for precision, not for confounding control.
##   * Type enters every pooled model. Disease type is a baseline variable and
##     the fatigue effect may differ by disease.
##   * Status is a common cause, not a mediator. Disease activity can cause
##     fatigue AND change chromatin:
##         activity -> fatigue        activity -> chromatin
##     For Status to be a mediator, fatigue would have to cause activity, which
##     is not the assumed direction. A common cause is a confounder, so the
##     default is to adjust for it. The MFI term then measures the part of
##     fatigue that is NOT explained by disease activity, which is the question
##     of interest: fatigue that persists in remission.
##     The cost: if fatigue is partly caused by activity, adjusting removes the
##     part of the signal that runs through activity, so the effect shrinks.
##     Set adjust_for_status = FALSE to fit the unadjusted version and compare.
##     In subsets that hold one Status only, the term is dropped automatically.
##   * batch enters only when the subset spans more than one batch. Note that
##     batch and Status overlap in the CD cohort, so they compete for the same
##     variance and the Status coefficient is not interpretable there.
##   * sex enters only when both levels are present with at least 2 samples
##     each. A 5/1 split spends a parameter and estimates nothing.
##
## Rules (status contrast, non-HC):
##   * Status is the exposure. In PBMCs, activity changes cell composition, and
##     composition changes FRiP:  activity -> composition -> FRiP.
##     FRiP is then DOWNSTREAM of Status, not a common cause. Adjusting for it
##     removes part of the real Status effect. Not adjusting keeps any
##     technical quality difference. Fit both (frip_adjust TRUE and FALSE).
##   * Fatigue group is also downstream of activity (activity -> fatigue).
##     It is not a confounder of Status, so it is not in the model by default.

## ----------------------------------------------------------------- config ----

## Snakemake param, with a default. Fit both ways and compare the MFI results
## before you decide which to report.
adjust_for_status <- tryCatch(
  as.logical(snakemake@params$adjust_for_status),
  error = function(e) NA
)
if (length(adjust_for_status) != 1 || is.na(adjust_for_status)) {
  adjust_for_status <- TRUE
}

## ---------------------------------------------------------------- helpers ----

## Columns that may appear in a design. Only these are touched: converting all
## 131 metadata columns to factors would change sample IDs, dates and free text.
design_cols <- intersect(c("Status", "group", "Type", "batch", "sex"),
                         names(metadata))

## Empty factor levels build all-zero model matrix columns and make the model
## rank-deficient. Drop them once, here.
metadata[design_cols] <- lapply(metadata[design_cols], function(x) {
  if (is.character(x)) x <- factor(x)
  if (is.factor(x))    x <- droplevels(x)
  x
})

## Re-level without losing values. factor(x, levels = ...) silently turns any
## value not in `levels` into NA, which drops the sample from the model matrix
## but not from the count matrix. Fail instead.
set_ref <- function(x, ref) {
  x <- droplevels(as.factor(x))
  if (!ref %in% levels(x)) {
    stop("reference level '", ref, "' not present. levels: ",
         paste(levels(x), collapse = ", "))
  }
  stats::relevel(x, ref = ref)
}

## Samples with NA in any design column are dropped by model.matrix but stay in
## the counts, which misaligns the two. Catch it before the model is built.
assert_complete <- function(md, cols, label = "") {
  cols <- intersect(cols, names(md))
  bad  <- cols[vapply(md[cols], anyNA, logical(1))]
  if (length(bad)) {
    stop(label, ": NA values in design column(s): ",
         paste(bad, collapse = ", "),
         ". Fix the metadata or drop those samples upstream.")
  }
  invisible(TRUE)
}

## A factor is usable if it has >= 2 levels and every level has >= min_n rows.
usable_factor <- function(x, min_n = 2) {
  if (is.null(x))    return(FALSE)
  if (!is.factor(x)) return(FALSE)
  tb <- table(x)
  length(tb) >= 2 && all(tb >= min_n)
}

## Build a formula from a term vector, keeping only terms that keep the model
## full rank. Terms are added in order, so put the term of interest first.
build_design <- function(terms, md, label = "") {
  terms <- terms[!is.na(terms) & nzchar(terms)]
  assert_complete(md, terms, label)
  kept    <- character(0)
  dropped <<- character(0)          # recorded in design_record
  for (tm in terms) {
    trial <- c(kept, tm)
    f  <- stats::as.formula(paste("~", paste(trial, collapse = " + ")))
    mm <- try(stats::model.matrix(f, md), silent = TRUE)
    if (inherits(mm, "try-error")) {
      message("  dropped (error): ", tm)
      dropped <<- c(dropped, tm); next
    }
    if (qr(mm)$rank < ncol(mm)) {
      message("  dropped (aliased): ", tm)
      dropped <<- c(dropped, tm); next
    }
    if (ncol(mm) >= nrow(md)) {
      message("  dropped (no residual df left): ", tm)
      dropped <<- c(dropped, tm); next
    }
    kept <- trial
  }
  f  <- stats::as.formula(paste("~", paste(kept, collapse = " + ")))
  mm <- stats::model.matrix(f, md)
  rdf <- nrow(md) - ncol(mm)
  message(sprintf("  %s: n = %d, params = %d, residual df = %d",
                  label, nrow(md), ncol(mm), rdf))
  message("  design: ", paste(deparse(f), collapse = ""))
  if (rdf < 4) {
    warning(label, ": only ", rdf, " residual df. Coefficients from this ",
            "model are unstable and depend on covariate choice. Do not ",
            "report a peak count without a sensitivity check.")
  }
  f
}

## ------------------------------------------------ non-fatigue diagnostics ----

## Variance inflation of each model column. A high value on the term of
## interest means a covariate carries the same information, so the test has
## little power. VIF > 5 is a warning, VIF > 10 is a serious problem.
design_vif <- function(f, md) {
  mm <- stats::model.matrix(f, md)
  x  <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]
  if (ncol(x) < 2) return(stats::setNames(rep(1, ncol(x)), colnames(x)))
  vif <- vapply(seq_len(ncol(x)), function(j) {
    r2 <- summary(stats::lm(x[, j] ~ x[, -j, drop = FALSE]))$r.squared
    1 / max(1 - r2, 1e-12)
  }, numeric(1))
  stats::setNames(vif, colnames(x))
}

## Association of FRiP with the term of interest.
## Spearman for numeric terms, Wilcoxon and R^2 for two-level factors.
qc_association <- function(md, term, cols) {
  cols <- intersect(cols, names(md))
  g    <- md[[term]]
  rows <- lapply(cols, function(q) {
    v <- suppressWarnings(as.numeric(md[[q]]))
    if (all(is.na(v))) return(NULL)
    if (is.factor(g) && nlevels(g) == 2) {
      p  <- suppressWarnings(stats::wilcox.test(v ~ g)$p.value)
      r2 <- summary(stats::lm(v ~ g))$r.squared
      data.frame(column = q, test = "wilcoxon", p = p, r2 = r2)
    } else if (is.numeric(g)) {
      ct <- suppressWarnings(stats::cor.test(v, g, method = "spearman"))
      data.frame(column = q, test = "spearman", p = ct$p.value,
                 r2 = unname(ct$estimate)^2)
    } else NULL
  })
  do.call(rbind, rows)
}

## Cross-tables of the term of interest with each factor covariate. A covariate
## level that holds one term level only means full overlap there.
factor_overlap <- function(md, term, covs) {
  covs <- intersect(covs, names(md))
  covs <- covs[vapply(md[covs], is.factor, logical(1)) & covs != term]
  lapply(stats::setNames(covs, covs), function(cv) {
    tb <- table(md[[cv]], md[[term]])
    message("  ", term, " x ", cv, ":")
    print(tb)
    one_level <- rowSums(tb > 0) == 1
    if (any(one_level)) {
      warning(term, " and ", cv, " overlap: level(s) ",
              paste(rownames(tb)[one_level], collapse = ", "),
              " contain one ", term, " level only. The ", term,
              " effect in these samples cannot be separated from ", cv, ".")
    }
    as.data.frame(tb)
  })
}

## --------------------------------------------------------------- context ----

is_pooled <- !has_subset_type || pool_types_active

## Covariates offered to every model, in order of priority. The builder keeps
## the ones that survive. Sex is last, so it is the first term sacrificed when
## degrees of freedom run out.
covars <- c(
  if (is_pooled)                     "Type",
  if (usable_factor(metadata$batch)) "batch",
  if (usable_factor(metadata$sex))   "sex"
)

message("\n[design] subset = ",
        if (pool_types_active) paste0("pooled (", subset_type, " coefficient)")
        else if (is_pooled) "pooled" else subset_type,
        " | adjust_for_status = ", adjust_for_status)

## Status is added only where it varies. usable_factor() returns FALSE in
## status-stratified subsets, so the term disappears there without a branch.
status_term <- if (adjust_for_status && usable_factor(metadata$Status)) {
  "Status"
} else {
  NULL
}

## Non-fatigue terms only.
frip_term <- if (frip_adjust) "FRiP_scaled" else NULL

## ----------------------------------------------------------- the formula ----

term_of_interest <- NULL   # used by the non-fatigue diagnostics

if (status_contrast) {
  metadata$Status <- set_ref(metadata$Status, "Remission")
  term_of_interest <- "Status"
  
  if (hc_comparison) {
    ## Unchanged HC path. group holds only {No_Fatigue, <level>_Fatigue}, so
    ## the reference must match hc_fatigue_level.
    group_ref <- paste0(hc_fatigue_level, "_Fatigue")
    metadata$group <- set_ref(metadata$group, group_ref)
    design_formula <- build_design(
      c("Status", "group", "FRiP_scaled", covars), metadata, "status_contrast")
    
  } else {
    ## Fatigue group is downstream of activity. Add it only on request.
    group_term <- NULL
    if (adjust_for_group) {
      metadata$group <- set_ref(metadata$group, "Low_Fatigue")
      group_term <- "group"
    }
    
    if (pool_types_active) {
      metadata$Type <- set_ref(metadata$Type, "CD")
      assert_complete(metadata, c("Type", "Status"), "status_contrast_pooled")
      cell_tb <- table(metadata$Type, metadata$Status)
      message("  Type x Status:"); print(cell_tb)
      if (any(cell_tb < 2)) {
        stop("pool_types needs >= 2 samples in every Type x Status cell.")
      }
      ## ~ Type + Type:Status gives one Active-vs-Remission coefficient per type.
      expected_coef <- paste0("Type", subset_type, ".StatusActive")
      if (!identical(deseq_coef, expected_coef)) {
        stop("pool_types: set deseq_coef to '", expected_coef, "', not '",
             deseq_coef, "'.")
      }
      design_formula <- build_design(
        c("Type", "Type:Status", group_term, frip_term,
          setdiff(covars, "Type")),
        metadata, "status_contrast_pooled")
      if (!grepl("Type:Status", paste(deparse(design_formula), collapse = ""))) {
        stop("pool_types: the Type:Status term was dropped. Check the log.")
      }
    } else {
      design_formula <- build_design(
        c("Status", group_term, frip_term, covars),
        metadata, "status_contrast")
    }
  }
  
} else if (continuous) {
  ## stratified_by_status no longer changes the formula: it selects the
  ## samples upstream, and must not also change the model.
  design_formula <- build_design(
    c("MFI_scaled", "FRiP_scaled", status_term, covars),
    metadata, "continuous")
  
} else if (grepl("Fatigue", deseq_coef)) {
  metadata$group <- set_ref(metadata$group, "Low_Fatigue")
  ## Same reasoning as the continuous branch: fatigue group is the exposure,
  ## and Status is a common cause of exposure and outcome.
  design_formula <- build_design(
    c("group", "FRiP_scaled", status_term, covars),
    metadata, "categorical_fatigue")
  
} else {
  metadata$Type <- set_ref(metadata$Type, "CD")
  term_of_interest <- "Type"
  ## Type is the term of interest. Status is a covariate here regardless of
  ## adjust_for_status, because that switch concerns the fatigue models only.
  ## Status is kept: the question is CD vs UC at equal disease activity.
  design_formula <- build_design(
    c("Type", setdiff(covars, "Type"), frip_term,
      if (usable_factor(metadata$Status)) "Status" else NULL),
    metadata, "type_contrast")
}

## ------------------------------------------------------------ final check ----

mm <- stats::model.matrix(design_formula, metadata)
stopifnot(qr(mm)$rank == ncol(mm))
stopifnot(nrow(metadata) > ncol(mm))

## model.matrix silently drops rows with NA. If it did, the metadata no longer
## lines up with the count matrix columns and every result would be wrong.
if (nrow(mm) != nrow(metadata)) {
  stop("model.matrix returned ", nrow(mm), " rows for ", nrow(metadata),
       " samples. Some design column contains NA.")
}

## Non-fatigue diagnostics: collinearity, QC association, factor overlap.
diag_record <- NULL
if (!is_fatigue_model && !is.null(term_of_interest)) {
  message("\n[diagnostics] term of interest = ", term_of_interest)
  vif <- design_vif(design_formula, metadata)
  message("  VIF: ", paste(sprintf("%s = %.1f", names(vif), vif), collapse = " | "))
  if (any(vif > 5)) {
    warning("VIF > 5 for: ", paste(names(vif)[vif > 5], collapse = ", "),
            ". These columns share information. Compare fits with and ",
            "without the correlated covariate before you interpret peak counts.")
  }
  qc_assoc <- qc_association(metadata, term_of_interest, "FRiP")
  if (!is.null(qc_assoc)) {
    message("  FRiP vs ", term_of_interest, ":")
    print(qc_assoc)
    if (any(qc_assoc$p < 0.05)) {
      warning("FRiP differs by ", term_of_interest, ". Adjusting for it can ",
              "remove real signal; not adjusting can keep technical signal. ",
              "Run both (frip_adjust TRUE and FALSE).")
    }
  }
  overlap <- factor_overlap(metadata, term_of_interest,
                            c("batch", "sex", "Type", "group"))
  diag_record <- list(vif = vif, qc_association = qc_assoc,
                      factor_overlap = overlap)
}

## Record the design so every result file can be traced back to its model.
design_record <- list(
  subset            = if (pool_types_active) paste0("pooled_", subset_type)
  else if (is_pooled) "pooled" else subset_type,
  formula           = paste(deparse(design_formula), collapse = ""),
  adjust_for_status = adjust_for_status,
  n                 = nrow(metadata),
  n_params          = ncol(mm),
  residual_df       = nrow(metadata) - ncol(mm),
  dropped_terms     = dropped,
  batch_table       = if ("batch" %in% names(metadata))
    as.data.frame(table(metadata$batch)) else NULL,
  sex_table         = if ("sex" %in% names(metadata))
    as.data.frame(table(metadata$sex)) else NULL,
  status_table      = if ("Status" %in% names(metadata))
    as.data.frame(table(metadata$Status)) else NULL
)
if (!is_fatigue_model) {
  design_record$frip_adjust        <- frip_adjust
  design_record$adjust_for_group   <- adjust_for_group
  design_record$pool_types         <- pool_types_active
  design_record$diagnostics        <- diag_record
}
saveRDS(design_record, snakemake@output$design_record)
# ── Build DESeq2 dataset ──────────────────────────────────────────────────────
#make sure order of rownames(metadata) is same as colnames(full_counts)
message("Building DESeq2 dataset...")
dds <- DESeqDataSetFromMatrix(
  countData = full_counts[, rownames(metadata)],
  colData   = metadata,
  design    = design_formula
)

if (is_fatigue_model) {
  # Filter low-count peaks and sex chromosomes (unchanged)
  keep <- rowSums(counts(dds)) >= min_row_sums
  dds  <- dds[keep, ]
  dds  <- dds[!grepl("chr[XY]", rownames(dds)), ]
} else {
  ## Remove unplaced contigs, chrM and sex chromosomes BEFORE the fit, so they
  ## do not enter size factors, dispersion trend or the FDR correction.
  dds <- dds[!grepl("chr[XYM]|chrUn|random|_alt|EBV", rownames(dds),
                    ignore.case = TRUE), ]
  ## Group-aware filter: a peak must have enough counts in at least as many
  ## samples as the smallest group. A total-count filter keeps peaks that are
  ## open in one or two samples only.
  keep <- edgeR::filterByExpr(counts(dds),
                              design = stats::model.matrix(design_formula,
                                                           as.data.frame(colData(dds))))
  message(sprintf("filterByExpr: keeping %d of %d peaks", sum(keep), length(keep)))
  dds <- dds[keep, ]
}

colSums(counts(dds))
# ── Run DESeq2 ────────────────────────────────────────────────────────────────
message("Running DESeq2...")
set.seed(123)

use_loess_norm <- isTRUE(snakemake@params$loess_norm)   # absent param -> FALSE, current behavior

if (use_loess_norm) {
  message("Applying loess normalization offsets (csaw) for intensity-dependent bias...")
  suppressPackageStartupMessages(library(csaw))
  cnt <- counts(dds)
  se  <- SummarizedExperiment(assays = list(counts = cnt))
  se$totals <- colSums(cnt)                          # <- what csaw was missing
  se <- normOffsets(se, se.out = TRUE)
  nf <- exp(assay(se, "offset"))                     # log-offsets -> multiplicative
  nf <- nf / exp(rowMeans(log(nf)))                  # row geomean = 1 (DESeq2 requirement)
  normalizationFactors(dds) <- nf
} else {
  dds <- estimateSizeFactors(dds, type = "poscounts")
}

dds <- DESeq(dds, fitType = "local", minReplicatesForReplace = Inf)
vsd <- vst(dds, blind = TRUE)

if (!deseq_coef %in% resultsNames(dds)) {
  stop("deseq_coef '", deseq_coef, "' not in resultsNames(dds): ",
       paste(resultsNames(dds), collapse = ", "))
}

# apeglm shrinkage requires a factor coefficient; for continuous MFI use "normal"
if (continuous) {
  res <- lfcShrink(dds, coef = deseq_coef, type = "normal")
} else {
  res <- lfcShrink(dds, coef = deseq_coef, type = "apeglm")
}

filtered_res <- as.data.frame(res) %>%
  tibble::rownames_to_column("rowname") %>%
  dplyr::filter(!grepl("chrUn|random", rowname, ignore.case = TRUE))

if (!is_fatigue_model) {
  ## lfcShrink keeps the p-values of the UNSHRUNK Wald test but replaces the
  ## LFC. A volcano plot of shrunk LFC vs Wald p then shows significant peaks
  ## near LFC = 0. Save the unshrunk values next to the shrunk ones.
  ##   log2FoldChange      shrunk (use for ranking by effect, MA plots)
  ##   log2FoldChange_MLE  unshrunk (use for volcano x-axis)
  ##   stat                Wald statistic (use for GSEA / fgsea ranking)
  res_mle <- as.data.frame(results(dds, name = deseq_coef))
  idx <- match(filtered_res$rowname, rownames(res_mle))
  filtered_res$log2FoldChange_MLE <- res_mle$log2FoldChange[idx]
  filtered_res$lfcSE_MLE          <- res_mle$lfcSE[idx]
  filtered_res$stat               <- res_mle$stat[idx]
  
  sig <- which(filtered_res$padj < padj_threshold)
  message(sprintf(
    "Shrinkage QC: %d significant peaks | %d with |shrunk LFC| < 0.1 | %d with NA shrunk LFC",
    length(sig),
    sum(abs(filtered_res$log2FoldChange[sig]) < 0.1, na.rm = TRUE),
    sum(is.na(filtered_res$log2FoldChange[sig]))))
}

frac_neg   <- mean(filtered_res$log2FoldChange < 0, na.rm = TRUE)
med_lfc    <- median(filtered_res$log2FoldChange, na.rm = TRUE)
hi_bm      <- filtered_res$baseMean > quantile(filtered_res$baseMean, 0.75, na.rm = TRUE)
med_lfc_hi <- median(filtered_res$log2FoldChange[hi_bm], na.rm = TRUE)
message(sprintf(
  "Skew QC: frac negative LFC = %.3f | median LFC = %.4f | median LFC (top baseMean quartile) = %.4f",
  frac_neg, med_lfc, med_lfc_hi))

# ── Outputs ───────────────────────────────────────────────────────────────────
message("Saving outputs...")
saveRDS(dds, snakemake@output$dds_rds)
saveRDS(vsd, snakemake@output$vsd_rds)
write_csv(filtered_res, snakemake@output$differential_csv)

vsd_df <- as.data.frame(assay(vsd)) %>% rownames_to_column("region")
write_csv(vsd_df, snakemake@output$vsd_counts_csv)

plot_pca(vsd, "group",  snakemake@output$group_pca_png)
plot_pca(vsd, "sex",    snakemake@output$sex_pca_png)
plot_pca(vsd, "Type",   snakemake@output$type_pca_png)
plot_pca(vsd, "Status", snakemake@output$status_pca_png)

message(
  "DESeq2 complete. ",
  sum(filtered_res$padj < padj_threshold, na.rm = TRUE),
  " peaks significant at padj < ", padj_threshold
)