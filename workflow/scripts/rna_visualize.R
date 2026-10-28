# workflow/scripts/rna_visualize.R
# Snakemake rule: rna_visualize
# Generates PCA plot, heatmap of top variable genes, and volcano plot.

log_file <- snakemake@log[[1]]
con      <- file(log_file, open = "wt")
sink(con, type = "message")
sink(con, type = "output")

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(readr)
  library(matrixStats)
  library(RColorBrewer)
  library(DESeq2)
})

source(file.path(snakemake@scriptdir, "utils.R"))

# ── Parameters ────────────────────────────────────────────────────────────────
top_var_peaks   <- snakemake@params$top_var_peaks
grouping_factor <- snakemake@params$grouping_factor
deseq_coef      <- snakemake@params$deseq_coef

# ── Load inputs ───────────────────────────────────────────────────────────────
message("Loading inputs...")
dds      <- readRDS(snakemake@input$dds_rds)
vsd      <- readRDS(snakemake@input$vsd_rds)
res_df   <- read_csv(snakemake@input$results_csv, show_col_types = FALSE)
G_list   <- readRDS(snakemake@input$glist_rds)
metadata <- readRDS(snakemake@input$metadata_rds)

# ── PCA ───────────────────────────────────────────────────────────────────────
message("Generating PCA plot...")
plot_pca(vsd, grouping_factor, snakemake@output$pca_png)

# ── Heatmap of top variable genes ─────────────────────────────────────────────
message("Generating heatmap...")
count_mat <- data.frame(assay(dds))
peak_var  <- rowVars(as.matrix(count_mat))
top_peaks <- order(peak_var, decreasing = TRUE)[seq_len(top_var_peaks)]

mat_top <- count_mat[top_peaks, ] %>%
  rownames_to_column("ensembl_gene_id") %>%
  left_join(dplyr::select(G_list, ensembl_gene_id, hgnc_symbol),
            by = "ensembl_gene_id") %>%
  dplyr::filter(!is.na(hgnc_symbol), hgnc_symbol != "") %>%
  mutate(variance = rowVars(dplyr::select(., where(is.numeric)) %>% as.matrix())) %>%
  group_by(hgnc_symbol) %>%
  slice_max(variance, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  dplyr::select(-ensembl_gene_id, -variance) %>%
  column_to_rownames("hgnc_symbol")

grouping_vector <- factor(metadata[[grouping_factor]])

groups <- levels(grouping_vector)
n_groups <- length(groups)
pal <- brewer.pal(3, "Set2")
group_cols <- setNames(pal[seq_len(n_groups)], groups)

col_side <- group_cols[grouping_vector]
coul     <- colorRampPalette(brewer.pal(8, "PiYG"))(25)

png(snakemake@output$heatmap_png, width = 6, height = 8, units = "in", res = 300)
heatmap(as.matrix(mat_top), col = coul, ColSideColors = col_side, labRow = NA, Colv = NA)
dev.off()
message("Heatmap saved to: ", snakemake@output$heatmap_png)

# ── Volcano ───────────────────────────────────────────────────────────────────
message("Generating volcano plot...")
volcano_plot(
  df        = res_df,
  gene_col  = "hgnc_symbol",
  title     = "Expression",
  plot_path = snakemake@output$volcano_png,
  comparison_title = deseq_coef,
  lfc_thr   = snakemake@params$lfc_threshold,
  pval_thr  = snakemake@params$pval_threshold,
  padj_thr  = snakemake@params$padj_threshold
)

message("Visualization complete.")
