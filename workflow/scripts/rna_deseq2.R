log_file <- snakemake@log[[1]]
con      <- file(log_file, open = "wt")
sink(con, type = "message")
sink(con, type = "output")

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(readr)
  library(DESeq2)
})

# ── Parameters ────────────────────────────────────────────────────────────────
min_counts      <- snakemake@params$min_counts
min_samples     <- snakemake@params$min_samples
deseq_coef      <- snakemake@params$deseq_coef
subset_type     <- snakemake@params$subset_type    # "CD", "UC", or NULL
subset_status   <- snakemake@params$subset_status  # "Active", "Remission", or NULL
continuous      <- isTRUE(snakemake@params$continuous)
status_contrast <- isTRUE(snakemake@params$status_contrast)
hc_comparison   <- isTRUE(snakemake@params$hc_comparison)

## Whether to adjust the fatigue term for disease activity. Status is a common
## cause of fatigue and of expression, so it is a confounder, not a mediator,
## and the default is to adjust. Fit both ways and compare before reporting.
adjust_for_status <- tryCatch(as.logical(snakemake@params$adjust_for_status),
                              error = function(e) NA)
if (length(adjust_for_status) != 1 || is.na(adjust_for_status)) {
  adjust_for_status <- TRUE
}

## Quality covariate. There is no RNA-seq equivalent of FRiP: DESeq2 size
## factors already absorb sequencing depth, so do not substitute library size.
## If a real quality metric exists in the metadata, name it here. For this
## cohort the strongest candidate is globin fraction, because the whole-blood
## libraries were made without globin depletion, so globin content changes how
## much of each library is informative. Other options: exonic rate, RIN,
## duplication rate. Leave NULL to fit without one.
qc_covar <- snakemake@params$qc_covar               # e.g. "globin_frac_scaled"
if (is.null(qc_covar) || identical(qc_covar, "None")) qc_covar <- NULL

## ── Params shared with the ATAC script (non-fatigue models only) ────────────
## Fatigue models and HC comparisons keep the old code path.
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

## adjust_for_group: fatigue group as a covariate in the non-HC status
##   contrast. Default FALSE, because fatigue is downstream of activity.
## pool_types: fit CD and UC together with ~ Type + Type:Status and read the
##   Status coefficient of one type. Default FALSE.
adjust_for_group <- param_flag("adjust_for_group", FALSE)
pool_types       <- param_flag("pool_types", FALSE)

has_subset_type   <- !is.null(subset_type) && subset_type != "None"
pool_types_active <- status_contrast && !hc_comparison && pool_types && has_subset_type

## frip_adjust is an ATAC param. It has no effect here, so the _noFRiP entries
## in COMPARISONS give the same RNA-seq result as the entries without the suffix.
if (!is.null(snakemake@params[["frip_adjust"]])) {
  message("Note: frip_adjust is an ATAC-seq param and has no effect in this script.")
}

# ── Helper ────────────────────────────────────────────────────────────────────
strip_ensembl_version <- function(ids) sub("\\..*", "", ids)

# ── Load inputs ───────────────────────────────────────────────────────────────
message("Loading inputs...")
ct_clean <- readRDS(snakemake@input$counts_rds)
metadata <- readRDS(snakemake@input$metadata_rds)
G_list   <- readRDS(snakemake@input$glist_rds)

stopifnot(all(rownames(metadata) %in% colnames(ct_clean)))

# ── Subset samples ────────────────────────────────────────────────────────────
# Three cases, because the HC comparisons need the controls kept, and the
# ordinary comparisons need them removed.
stratified_by_status <- !is.null(subset_status) && subset_status != "None"

if (hc_comparison) {
  ## Healthy controls versus the quietest patients: remission AND low fatigue.
  ## `subset_type` names which disease to compare against; NULL means both.
  ## Do NOT filter on Type first — that would remove every control.
  message("HC comparison: keeping HC plus ",
          if (is.null(subset_type) || subset_type == "None") "CD and UC"
          else subset_type,
          " (Remission, Low_Fatigue only)")
  
  is_hc      <- metadata$Type == "HC"
  is_patient <- if (!is.null(subset_type) && subset_type != "None") {
    metadata$Type == subset_type
  } else {
    metadata$Type != "HC"
  }
  keep_rows <- which(
    is_hc |
      (is_patient &
         metadata$Status == "Remission" &
         metadata$group  == "Low_Fatigue")
  )
  metadata <- metadata[keep_rows, ]
  
  n_hc <- sum(metadata$Type == "HC")
  message("  HC: ", n_hc, " | patients: ", nrow(metadata) - n_hc)
  if (n_hc < 2 || nrow(metadata) - n_hc < 2) {
    stop("HC comparison needs at least 2 samples per side. Check that ",
         "controls are coded Type == 'HC', Status == 'Healthy' and ",
         "group == 'No_Fatigue'.")
  }
  
} else {
  if (has_subset_type && !pool_types_active) {
    message("Subsetting to Type == ", subset_type)
    metadata <- metadata[which(metadata$Type == subset_type), ]
  } else {
    ## Pooled mode also comes here: keep CD and UC, the coefficient selects the type.
    message("Removing Healthy Controls")
    metadata <- metadata[which(metadata$Type != "HC"), ]
  }
  
  # Not applied for status_contrast comparisons — Status is the main factor
  if (stratified_by_status && !status_contrast) {
    message("Subsetting to Status == ", subset_status)
    metadata <- metadata[which(metadata$Status == subset_status), ]
  }
}

ct_clean <- ct_clean[, rownames(metadata), drop = FALSE]
stopifnot(identical(rownames(metadata), colnames(ct_clean)))
message("Samples after subsetting: ", nrow(metadata))

# ── Filter sex chromosome genes ───────────────────────────────────────────────
autosomal <- G_list$ensembl_gene_id[
  !G_list$chromosome_name %in% c("X", "Y", "chrX", "chrY")
]

# ── Design formula ────────────────────────────────────────────────────────────
# Built as a term list, not as nested branches. Each covariate is added only if
# it can be estimated in the current subset. The builder reports what it kept,
# what it dropped, and how many residual degrees of freedom remain.
#
# Rules:
#   * MFI_scaled, group or Status is the term of interest and goes first.
#   * qc_covar, if given, is a precision covariate. It explains count variance;
#     it is not there to remove confounding.
#   * Type enters every pooled model. Disease type is a baseline variable and
#     the fatigue effect may differ by disease.
#   * Status is a common cause, not a mediator:
#         activity -> fatigue        activity -> expression
#     For Status to be a mediator, fatigue would have to cause activity, which
#     is not the assumed direction. So the default is to adjust for it. The MFI
#     term then measures the part of fatigue NOT explained by disease activity,
#     which is the question of interest. In status-stratified subsets the term
#     is constant and is dropped automatically.
#   * In the non-HC status contrast, Status is the exposure and fatigue group
#     is downstream of it (activity -> fatigue). Group is not a confounder
#     there, so it enters only when adjust_for_group = TRUE.
#   * batch enters only when the subset spans more than one batch.
#   * sex enters only when both levels are present with at least 2 samples
#     each. A 5/1 split spends a parameter and estimates nothing. Note that
#     sex-chromosome genes are already removed above, so sex here adjusts for
#     autosomal expression differences only.

## Columns that may appear in a design. Only these are touched: converting the
## whole metadata to factors would change sample IDs, dates and free text.
design_cols <- intersect(c("Status", "group", "Type", "batch", "sex"),
                         names(metadata))

## Empty factor levels build all-zero model matrix columns and make the model
## rank-deficient. Subsetting above is exactly what creates them. Drop once.
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
## NA is deliberately NOT checked here: a design column with NA is a metadata
## problem, and assert_complete() must raise it as an error. Silently dropping
## such a term would hide the problem instead.
usable_factor <- function(x, min_n = 2) {
  if (is.null(x))    return(FALSE)
  if (!is.factor(x)) return(FALSE)
  tb <- table(x)
  length(tb) >= 2 && all(tb >= min_n)
}

dropped <- character(0)

## Build a formula from a term vector, keeping only terms that keep the model
## full rank. Terms are added in order, so put the term of interest first.
build_design <- function(terms, md, label = "") {
  terms <- terms[!is.na(terms) & nzchar(terms)]
  assert_complete(md, terms, label)
  kept    <- character(0)
  dropped <<- character(0)
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
  f   <- stats::as.formula(paste("~", paste(kept, collapse = " + ")))
  mm  <- stats::model.matrix(f, md)
  rdf <- nrow(md) - ncol(mm)
  message(sprintf("  %s: n = %d, params = %d, residual df = %d",
                  label, nrow(md), ncol(mm), rdf))
  message("  design: ", paste(deparse(f), collapse = ""))
  if (rdf < 4) {
    warning(label, ": only ", rdf, " residual df. Coefficients from this ",
            "model are unstable and depend on covariate choice. Do not ",
            "report a gene count without a sensitivity check.")
  }
  f
}

is_pooled <- !has_subset_type || pool_types_active

## Covariates offered to every model, in priority order. Sex is last, so it is
## the first term sacrificed when degrees of freedom run out.
covars <- c(
  if (is_pooled)                     "Type",
  if (usable_factor(metadata$batch)) "batch",
  if (usable_factor(metadata$sex))   "sex"
)

message("\n[design] subset = ",
        if (pool_types_active) paste0("pooled (", subset_type, " coefficient)")
        else if (is_pooled) "pooled" else subset_type,
        " | status = ", if (stratified_by_status) subset_status else "both",
        " | adjust_for_status = ", adjust_for_status,
        " | qc_covar = ", if (is.null(qc_covar)) "none" else qc_covar)

## Status is added only where it varies. usable_factor() returns FALSE in
## status-stratified subsets, so the term disappears there without a branch.
status_term <- if (adjust_for_status && usable_factor(metadata$Status)) {
  "Status"
} else {
  NULL
}

if (status_contrast) {
  metadata$Status <- set_ref(metadata$Status, "Remission")
  
  if (hc_comparison) {
    ## Unchanged HC path. In an HC comparison, group is No_Fatigue for every
    ## control and Low_Fatigue for every patient, so it marks exactly the same
    ## split as Status, and build_design() drops group as aliased. The guard
    ## below covers subsets where group has only one level.
    group_term <- if (usable_factor(metadata$group)) {
      metadata$group <- set_ref(metadata$group, "Low_Fatigue")
      "group"
    } else {
      message("  group has one usable level in this subset, left out")
      NULL
    }
    design_formula <- build_design(
      c("Status", group_term, qc_covar, covars), metadata, "status_contrast")
    
  } else {
    ## Fatigue group is downstream of activity. Add it only on request.
    group_term <- NULL
    if (adjust_for_group) {
      if (usable_factor(metadata$group)) {
        metadata$group <- set_ref(metadata$group, "Low_Fatigue")
        group_term <- "group"
      } else {
        message("  group has one usable level in this subset, left out")
      }
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
        c("Type", "Type:Status", group_term, qc_covar, setdiff(covars, "Type")),
        metadata, "status_contrast_pooled")
      if (!grepl("Type:Status", paste(deparse(design_formula), collapse = ""))) {
        stop("pool_types: the Type:Status term was dropped. Check the log.")
      }
    } else {
      design_formula <- build_design(
        c("Status", group_term, qc_covar, covars), metadata, "status_contrast")
    }
  }
  
} else if (continuous) {
  ## stratified_by_status no longer changes the formula: it selects the samples
  ## upstream, and must not also change the model.
  design_formula <- build_design(
    c("MFI_scaled", qc_covar, status_term, covars), metadata, "continuous")
  
} else if (grepl("Fatigue", deseq_coef)) {
  metadata$group <- set_ref(metadata$group, "Low_Fatigue")
  design_formula <- build_design(
    c("group", qc_covar, status_term, covars), metadata, "categorical_fatigue")
  
} else {
  metadata$Type <- set_ref(metadata$Type, "CD")
  ## Type is the term of interest. Status is a covariate here regardless of
  ## adjust_for_status, which concerns the fatigue models only.
  design_formula <- build_design(
    c("Type", setdiff(covars, "Type"), qc_covar,
      if (usable_factor(metadata$Status)) "Status" else NULL),
    metadata, "type_contrast")
}

## model.matrix silently drops rows with NA. If it did, the metadata no longer
## lines up with the count matrix columns and every result would be wrong.
mm <- stats::model.matrix(design_formula, metadata)
stopifnot(qr(mm)$rank == ncol(mm))
stopifnot(nrow(metadata) > ncol(mm))
if (nrow(mm) != nrow(metadata)) {
  stop("model.matrix returned ", nrow(mm), " rows for ", nrow(metadata),
       " samples. Some design column contains NA.")
}

# ── Build DESeq2 dataset ──────────────────────────────────────────────────────
message("Building DESeq2 dataset...")
dds <- DESeqDataSetFromMatrix(
  countData = ct_clean[, rownames(metadata), drop = FALSE],
  colData   = metadata,
  design    = design_formula
)

keep <- rowSums(counts(dds) >= min_counts) >= min_samples
dds  <- dds[keep, ]
dds  <- dds[rownames(dds) %in% autosomal, ]
message("Genes kept: ", nrow(dds))

# ── Run DESeq2 ────────────────────────────────────────────────────────────────
message("Running DESeq2...")
set.seed(123)
dds <- DESeq(dds)
vsd <- vst(dds, blind = TRUE)

## The coefficient name must exist, or results() returns the wrong contrast or
## errors deep in the call. Check it against the fitted model.
if (!deseq_coef %in% resultsNames(dds)) {
  stop("deseq_coef '", deseq_coef, "' not in resultsNames(dds): ",
       paste(resultsNames(dds), collapse = ", "))
}
res <- results(dds, name = deseq_coef)

res_df <- as.data.frame(res) %>%
  arrange(padj) %>%
  mutate(ensembl_gene_id = strip_ensembl_version(rownames(.))) %>%
  left_join(G_list, by = "ensembl_gene_id") %>%
  dplyr::filter(!is.na(hgnc_symbol))

# ── Record the design ─────────────────────────────────────────────────────────
design_record <- list(
  subset            = if (pool_types_active) paste0("pooled_", subset_type)
  else if (is_pooled) "pooled" else subset_type,
  subset_status     = if (stratified_by_status) subset_status else "both",
  hc_comparison     = hc_comparison,
  formula           = paste(deparse(design_formula), collapse = ""),
  deseq_coef        = deseq_coef,
  adjust_for_status = adjust_for_status,
  qc_covar          = if (is.null(qc_covar)) NA_character_ else qc_covar,
  dropped_terms     = dropped,
  n                 = nrow(metadata),
  n_params          = ncol(mm),
  residual_df       = nrow(metadata) - ncol(mm),
  n_genes           = nrow(dds),
  median_lfcSE      = median(res$lfcSE, na.rm = TRUE),
  batch_table       = if ("batch" %in% names(metadata))
    as.data.frame(table(metadata$batch)) else NULL,
  sex_table         = if ("sex" %in% names(metadata))
    as.data.frame(table(metadata$sex)) else NULL,
  status_table      = if ("Status" %in% names(metadata))
    as.data.frame(table(metadata$Status)) else NULL
)
if (!is_fatigue_model) {
  design_record$adjust_for_group <- adjust_for_group
  design_record$pool_types       <- pool_types_active
}

# ── Save outputs ──────────────────────────────────────────────────────────────
message("Saving outputs...")
saveRDS(dds, snakemake@output$dds_rds)
saveRDS(vsd, snakemake@output$vsd_rds)
saveRDS(metadata, snakemake@output$metadata)
write_csv(res_df, snakemake@output$results_csv)
saveRDS(design_record, snakemake@output$design_record)

n_sig <- sum(res_df$padj < snakemake@params$padj_threshold, na.rm = TRUE)
message("DESeq2 complete. ", n_sig, " genes significant at padj < ",
        snakemake@params$padj_threshold,
        " | median lfcSE = ", signif(design_record$median_lfcSE, 3))