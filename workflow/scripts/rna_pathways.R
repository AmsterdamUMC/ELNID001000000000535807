# workflow/scripts/rna_pathways.R
# Snakemake rule: rna_pathways
# Runs fgsea GO:BP enrichment and generates pathway dot plots.

log_file <- snakemake@log[[1]]
con      <- file(log_file, open = "wt")
sink(con, type = "message")
sink(con, type = "output")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(forcats)
  library(msigdbr)
  library(fgsea)
})

source(file.path(snakemake@scriptdir, "utils.R"))

# ── Parameters ────────────────────────────────────────────────────────────────
pval_thr         <- snakemake@params$pval_threshold
padj_thr         <- snakemake@params$padj_threshold
fgsea_min_size   <- snakemake@params$fgsea_min_size
fgsea_max_size   <- snakemake@params$fgsea_max_size
nes_thr          <- snakemake@params$fgsea_nes_threshold
go_pattern       <- snakemake@params$go_pattern

# ── Load inputs ───────────────────────────────────────────────────────────────
message("Loading DE results...")
res_df <- read_csv(snakemake@input$results_csv, show_col_types = FALSE)

# ── Build gene rankings ────────────────────────────────────────────────────────
message("Building GSEA rankings...")
rankings <- res_df %>%
  dplyr::filter(!is.na(hgnc_symbol), hgnc_symbol != "",
                !is.na(pvalue), !is.na(log2FoldChange)) %>%
  group_by(hgnc_symbol) %>%
  slice_min(pvalue, n = 1, with_ties = FALSE) %>%
  mutate(score = -log10(pvalue) * sign(log2FoldChange)) %>%
  ungroup() %>%
  dplyr::filter(is.finite(score)) %>%
  arrange(desc(score)) %>%
  { setNames(.$score, .$hgnc_symbol) }

message("Ranked ", length(rankings), " genes.")

# ── Load GO:BP gene sets ───────────────────────────────────────────────────────
message("Loading MSigDB GO:BP pathways...")
msigdbr_gobp <- msigdbr(species = "Homo sapiens") %>%
  dplyr::filter(gs_subcollection == "GO:BP")
pathways <- split(msigdbr_gobp$gene_symbol, msigdbr_gobp$gs_name)

# ── Run fgsea ─────────────────────────────────────────────────────────────────
message("Running fgsea...")
set.seed(123)
fgsea_res <- fgsea(
  pathways = pathways,
  stats    = rankings,
  minSize  = fgsea_min_size,
  maxSize  = fgsea_max_size
) %>%
  mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ", "))

write.csv(fgsea_res, snakemake@output$fgsea_csv, row.names = FALSE, quote = FALSE)
message("Full fgsea results saved: ", nrow(fgsea_res), " pathways tested.")

# ── Filter significant results ────────────────────────────────────────────────
sig_paths <- fgsea_res %>%
  dplyr::filter(pval < pval_thr, abs(NES) > nes_thr) %>%
  arrange(NES) %>%
  mutate(pathway = gsub("GOBP_", "", pathway),
         pathway = fct_reorder(pathway, NES))

message(nrow(sig_paths), " significant pathways (pval < ", pval_thr,
        ", |NES| > ", nes_thr, ").")

# ── Dot plot: all significant pathways ────────────────────────────────────────
pathway_dotplot(
  pathways  = sig_paths %>% dplyr::filter(padj < padj_thr),
  title     = "Pathways Enriched in High Fatigue: RNA-seq",
  plot_path = snakemake@output$dotplot_top_png,
  x_col     = "NES",
  y_col     = "pathway",
  color_col = "pval",
  size_col  = "size"
)

# ── Dot plot: topic-specific subset ───────────────────────────────────────────
specific_paths <- sig_paths %>%
  dplyr::filter(grepl(go_pattern, pathway, ignore.case = TRUE)) %>%
  mutate(pathway = fct_reorder(as.character(pathway), NES))

pathway_dotplot(
  pathways  = specific_paths,
  title     = "Selected Pathways: High Fatigue RNA-seq",
  plot_path = snakemake@output$dotplot_sel_png,
  x_col     = "NES",
  y_col     = "pathway",
  color_col = "pval",
  size_col  = "size"
)

message("RNA pathway enrichment complete.")
