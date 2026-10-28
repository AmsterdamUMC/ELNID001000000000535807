# Run this once after creating the epicatac conda environment:
# conda activate epicatac
# Rscript workflow/scripts/install_epicatac.R

# if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
# remotes::install_github("GfellerLab/EPIC-ATAC", upgrade = "never")
# message("EPIC-ATAC installed successfully.")

# workflow/scripts/atac_celltype.R
# Estimates PBMC cell type proportions per sample using EPIC-ATAC,
# stratified by the comparison's relevant grouping variable.
# Produces boxplots (categorical) or scatter plots (continuous/status).
log_file <- snakemake@log[[1]]
con      <- file(log_file, open = "wt")
sink(con, type = "message")
sink(con, type = "output")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(scales)
  library(readxl)
  library(EPICATAC)
})

# ── Parameters ────────────────────────────────────────────────────────────────
comparison    <- snakemake@params$comparison
subset_type   <- snakemake@params$subset_type
subset_status <- snakemake@params$subset_status
continuous    <- isTRUE(snakemake@params$continuous)
status_contrast <- isTRUE(snakemake@params$status_contrast)
cols_to_remove  <- snakemake@params$cols_to_remove
hc_comparison    <- isTRUE(snakemake@params$hc_comparison)
hc_fatigue_level <- snakemake@params$hc_fatigue_level
if (is.null(hc_fatigue_level) || !nzchar(hc_fatigue_level)) hc_fatigue_level <- "Low"

# ── Determine x-axis grouping variable from comparison type ──────────────────
if (status_contrast) {
  group_var   <- "Status"
  x_label     <- "Disease Status"
} else if (continuous) {
  group_var   <- "MFI_scaled"
  x_label     <- "MFI Score (scaled)"
} else {
  group_var   <- "group"
  x_label     <- "Fatigue Group"
}

# ── Load metadata ─────────────────────────────────────────────────────────────
message("Loading metadata...")
metadata <- read_excel(snakemake@input$metadata_file) %>%
  dplyr::filter(batch != "no_ATACseq",
                !Patient_ID %in% cols_to_remove) %>%
  dplyr::mutate(MFI_scaled = as.numeric(scale(MFI_totalscore)))


# Apply type subset
if (hc_comparison) {
  keep_types <- if (!is.null(subset_type) && subset_type != "None") {
    c("HC", subset_type)
  } else {
    c("HC", "CD", "UC")
  }
  message("HC comparison: keeping Type in ", paste(keep_types, collapse = ", "))
  metadata <- metadata %>% dplyr::filter(Type %in% keep_types)
} else if (!is.null(subset_type) && subset_type != "None") {
  metadata <- metadata %>% dplyr::filter(Type == subset_type)
} else {
  metadata <- metadata %>% dplyr::filter(Type != "HC")
}

# Apply status subset (fatigue comparisons stratified by status)
if (!is.null(subset_status) && subset_status != "None" && !status_contrast) {
  metadata <- metadata %>% dplyr::filter(Status == subset_status)
}

# Restrict HC comparison to Healthy/Remission, No_Fatigue/Low_Fatigue only
if (hc_comparison) {
  fatigue_group <- paste0(hc_fatigue_level, "_Fatigue")
  message("HC comparison: keeping Status in {Healthy, Remission}, ",
          "group in {No_Fatigue, ", fatigue_group, "}")
  metadata <- metadata %>%
    dplyr::filter(Status %in% c("Healthy", "Remission"),
                  group  %in% c("No_Fatigue", fatigue_group))
}

message("Samples in this comparison: ", nrow(metadata))

# ── Load and prepare count matrix ─────────────────────────────────────────────
message("Loading count matrix...")
full_counts <- readRDS(snakemake@input$counts_rds)
full_counts <- full_counts[, metadata$Patient_ID]

# Convert rownames to chr:start-end format expected by EPIC-ATAC
rownames(full_counts) <- gsub("_", ":", rownames(full_counts))
rownames(full_counts) <- sub("^([^:]+:[^:]+):", "\\1-", rownames(full_counts))

# Keep standard chromosomes only
standard_chrs <- grepl("^chr([0-9]+|X|Y):[0-9]+-[0-9]+$", rownames(full_counts))
full_counts   <- full_counts[standard_chrs, ]
message("Peaks after chr filter: ", nrow(full_counts))

# ── Normalise ─────────────────────────────────────────────────────────────────
normalized <- EPICATAC:::get_TPMlike_counts(full_counts)

# ── Run EPIC-ATAC (with and without neutrophils) ──────────────────────────────
message("Running EPIC-ATAC (full reference)...")
out_full <- EPIC_ATAC(bulk = normalized, reference = atacRef_PBMC, ATAC = TRUE)

message("Running EPIC-ATAC (neutrophil-excluded reference)...")
custom_ref <- atacRef_PBMC
custom_ref$refProfiles     <- custom_ref$refProfiles[,
                                                     colnames(custom_ref$refProfiles) != "Neutrophils"]
custom_ref$refProfiles.var <- custom_ref$refProfiles.var[,
                                                         colnames(custom_ref$refProfiles.var) != "Neutrophils"]
custom_ref$sigGenes <- custom_ref$sigGenes[
  custom_ref$sigGenes %in% rownames(custom_ref$refProfiles)]
out_no_neut <- EPIC_ATAC(bulk = normalized, reference = custom_ref, ATAC = TRUE)

# ── Helper: build long data frame ─────────────────────────────────────────────
make_long <- function(epic_out, meta) {
  fracs <- as.data.frame(epic_out[["cellFractions"]]) %>%
    tibble::rownames_to_column("Patient_ID") %>%
    dplyr::left_join(
      meta %>% dplyr::select(Patient_ID, group, Type, Status, MFI_scaled),
      by = "Patient_ID"
    ) %>%
    tidyr::pivot_longer(
      cols      = -c(Patient_ID, group, Type, Status, MFI_scaled),
      names_to  = "cell_type",
      values_to = "fraction"
    )
  fracs
}

# ── Helper: make plot ─────────────────────────────────────────────────────────
make_plot <- function(df_long, group_var, x_label, title_suffix) {
  if (continuous || status_contrast && group_var == "Status") {
    # scatter for continuous MFI; also scatter for status (2 points per cell type)
  }
  
  if (continuous) {
    stat_df <- df_long %>%
      group_by(cell_type) %>%
      group_modify(~ {
        ok <- which(!is.na(.x$fraction) & !is.na(.x$MFI_scaled))   # guard NA columns
        ct <- suppressWarnings(cor.test(.x$fraction[ok], .x$MFI_scaled[ok], method = "spearman"))
        r2 <- summary(lm(fraction ~ MFI_scaled, data = .x[ok, ]))$r.squared
        tibble(rho = unname(ct$estimate), p = ct$p.value, r2 = r2,
               label = sprintf("rho=%.2f  p=%.3f\nR\u00b2=%.2f",
                               unname(ct$estimate), ct$p.value, r2))
      }) %>%
      ungroup()
    p <- ggplot(df_long, aes(x = .data[[group_var]], y = fraction,
                             label = Patient_ID)) +
      geom_point(size = 2.5, color = "#08306b", alpha = 0.8) +
      geom_smooth(method = "lm", se = TRUE, color = "#525252", linewidth = 0.8) +
      geom_text_repel(size = 2.5, max.overlaps = Inf) +
      facet_wrap(~ cell_type, scales = "free_y", nrow = 2) +
      scale_y_continuous(labels = percent_format()) +
      labs(x = x_label, y = "Cell fraction",
           title = paste0("Cell Type Proportions: ", title_suffix)) +
      theme_bw() +
      theme(strip.text = element_text(size = 10),
            plot.title = element_text(hjust = 0.5, face = "bold")) +
      geom_text(data = stat_df,
                aes(x = -Inf, y = Inf, label = label),
                hjust = -0.05, vjust = 1.15, size = 3, lineheight = 0.9,
                inherit.aes = FALSE)
  } else {
    # categorical: boxplot with jitter (group or Status on x)
    fill_vals <- if (group_var == "group") {
      c("High_Fatigue" = "#533AB7", "Low_Fatigue" = "#1D9E75", "No_Fatigue" = "#7FA6D6")
    } else {
      c("Active" = "#BB0C00", "Remission" = "#00AFBB", "Healthy" = "#4DAF4A")
    }
    p <- ggplot(df_long, aes(x = .data[[group_var]], y = fraction,
                             fill = .data[[group_var]], label = Patient_ID)) +
      geom_boxplot(alpha = 0.7, outlier.shape = NA) +
      geom_jitter(width = 0.1, size = 2, show.legend = FALSE) +
      geom_text_repel(size = 2.5, max.overlaps = Inf, show.legend = FALSE) +
      facet_wrap(~ cell_type, scales = "free_y", nrow = 2) +
      scale_fill_manual(values = fill_vals) +
      scale_y_continuous(labels = percent_format()) +
      labs(x = NULL, y = "Cell fraction", fill = NULL,
           title = paste0("Cell Type Proportions: ", title_suffix)) +
      theme_bw() +
      theme(legend.position = "bottom",
            strip.text      = element_text(size = 10),
            plot.title      = element_text(hjust = 0.5, face = "bold"))
  }
  p
}

# ── Comparison label for plot titles ─────────────────────────────────────────
comparison_label <- switch(comparison,
                           "CD_fatigue"                      = "High vs Low Fatigue (CD)",
                           "UC_fatigue"                      = "High vs Low Fatigue (UC)",
                           "UC_vs_CD"                        = "UC vs CD",
                           "CD_fatigue_Active"               = "High vs Low Fatigue (CD, Active)",
                           "CD_fatigue_Remission"            = "High vs Low Fatigue (CD, Remission)",
                           "UC_fatigue_Active"               = "High vs Low Fatigue (UC, Active)",
                           "UC_fatigue_Remission"            = "High vs Low Fatigue (UC, Remission)",
                           "UC_vs_CD_Active"                 = "UC vs CD (Active)",
                           "UC_vs_CD_Remission"              = "UC vs CD (Remission)",
                           "CD_fatigue_continuous"           = "MFI Score (CD)",
                           "UC_fatigue_continuous"           = "MFI Score (UC)",
                           "CD_fatigue_continuous_Active"    = "MFI Score (CD, Active)",
                           "CD_fatigue_continuous_Remission" = "MFI Score (CD, Remission)",
                           "UC_fatigue_continuous_Active"    = "MFI Score (UC, Active)",
                           "UC_fatigue_continuous_Remission" = "MFI Score (UC, Remission)",
                           "CD_Active_vs_Remission"          = "Active vs Remission (CD)",
                           "UC_Active_vs_Remission"          = "Active vs Remission (UC)",
                           "Active_vs_Remission"             = "Active vs Remission (pooled)",
                           comparison
)

# ── Generate plots and CSVs ───────────────────────────────────────────────────
for (version in c("full", "no_neutrophils")) {
  epic_out <- if (version == "full") out_full else out_no_neut
  df_long  <- make_long(epic_out, metadata)
  
  # Save CSV
  out_csv <- if (version == "full") {
    snakemake@output$fractions_csv
  } else {
    snakemake@output$fractions_no_neut_csv
  }
  df_long %>%
    dplyr::select(Patient_ID, cell_type, fraction, group, Type, Status, MFI_scaled) %>%
    write.csv(out_csv, row.names = FALSE)
  message("Fractions saved to: ", out_csv)
  
  # Save plot
  out_png <- if (version == "full") {
    snakemake@output$plot_png
  } else {
    snakemake@output$plot_no_neut_png
  }
  title_suffix <- paste0(comparison_label,
                         if (version == "no_neutrophils") " (no neutrophils)" else "")
  p <- make_plot(df_long, group_var, x_label, title_suffix)
  ggsave(out_png, plot = p, width = 14, height = 8, dpi = 300, bg = "white")
  message("Plot saved to: ", out_png)
}

message("Cell type deconvolution complete.")