# workflow/scripts/atac_annotate.R
# Snakemake rule: atac_annotate
# Annotates DESeq2 peaks with ChIPseeker and generates volcano plot.

log_file <- snakemake@log[[1]]
con      <- file(log_file, open = "wt")
sink(con, type = "message")
sink(con, type = "output")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(GenomicRanges)
  library(ChIPseeker)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  library(DESeq2)
})

source(file.path(snakemake@scriptdir, "utils.R"))

# ── Load inputs ───────────────────────────────────────────────────────────────
message("Loading inputs...")
dds           <- readRDS(snakemake@input$dds_rds)
gr_counts     <- readRDS(snakemake@input$gr_rds)
filtered_res  <- read_csv(snakemake@input$differential_csv, show_col_types = FALSE)
deseq_coef   <- snakemake@params$deseq_coef
comparison        <- snakemake@params$comparison
continuous        <- isTRUE(snakemake@params$continuous)

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
                           comparison   # fallback
)

# ── Annotate peaks ────────────────────────────────────────────────────────────
message("Annotating peaks with ChIPseeker...")
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

gr_names       <- paste0(seqnames(gr_counts), "_", start(gr_counts), "_", end(gr_counts))
names(gr_counts) <- gr_names
gr_tested      <- gr_counts[rownames(dds)]

peak_annot <- annotatePeak(gr_tested, TxDb = txdb, annoDb = "org.Hs.eg.db")
anno_df    <- as.data.frame(peak_annot) %>%
  mutate(rowname = paste0(seqnames, "_", start, "_", end)) %>%
  right_join(filtered_res, by = "rowname")

# ── Outputs ───────────────────────────────────────────────────────────────────
write_csv(anno_df, snakemake@output$annotated_csv)

volcano_plot(
  df         = anno_df,
  gene_col   = "SYMBOL",
  title      = "Accessibility",
  plot_path  = snakemake@output$volcano_png,
  comparison_title = comparison_label,
  lfc_thr    = snakemake@params$lfc_threshold,
  pval_thr   = snakemake@params$pval_threshold,
  padj_thr   = snakemake@params$padj_threshold
)

message("Annotation complete. Annotated peaks: ", nrow(anno_df))
