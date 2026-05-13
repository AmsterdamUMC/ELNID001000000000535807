# workflow/scripts/atac_pathways.R
# Snakemake rule: atac_pathways
# Runs GREAT pathway enrichment for High and Medium fatigue differential peaks.

log_file <- snakemake@log[[1]]
con      <- file(log_file, open = "wt")
sink(con, type = "message")
sink(con, type = "output")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(GenomicRanges)
  library(rGREAT)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(forcats)
})

source(file.path(snakemake@scriptdir, "utils.R"))

# ── Parameters ────────────────────────────────────────────────────────────────
great_padj_binom <- snakemake@params$great_padj_binom
great_padj_hyper <- snakemake@params$great_padj_hyper
lfc_thr          <- snakemake@params$lfc_threshold
pval_thr         <- snakemake@params$pval_threshold

# ── Load inputs ───────────────────────────────────────────────────────────────
message("Loading inputs...")
gr_counts    <- readRDS(snakemake@input$gr_rds)
diff_res     <- read_csv(snakemake@input$differential_csv, show_col_types = FALSE)

# ── Helper: run GREAT for one direction ───────────────────────────────────────
run_great_pathway <- function(gr_counts, diff_res, nm, out_csv) {
  filtered <- diff_res %>%
    dplyr::filter(abs(log2FoldChange) > lfc_thr, pvalue < pval_thr) %>%
    tidyr::separate(rowname, into = c("chr", "start", "end"),
                    sep = "_", remove = FALSE) %>%
    dplyr::mutate(start = as.integer(start), end = as.integer(end))

  if (nm == "hf") {
    filtered <- dplyr::filter(filtered, log2FoldChange > 0)
  } else {
    filtered <- dplyr::filter(filtered, log2FoldChange < 0)
  }

  if (nrow(filtered) == 0) {
    message("No significant peaks for condition: ", nm)
    write.csv(data.frame(), out_csv, row.names = FALSE)
    return(data.frame())
  }

  diff_gr <- makeGRangesFromDataFrame(
    filtered, seqnames.field = "chr",
    start.field = "start", end.field = "end", ignore.strand = TRUE
  )
  seqlevelsStyle(diff_gr) <- seqlevelsStyle(gr_counts)

  res <- great(diff_gr, background = gr_counts,
               "c5.go.bp", "TxDb.Hsapiens.UCSC.hg38.knownGene")
  tb  <- getEnrichmentTable(res)

  sig_tb <- tb %>%
    dplyr::filter(p_adjust <= great_padj_binom, p_adjust_hyper <= great_padj_hyper) %>%
    dplyr::arrange(desc(fold_enrichment)) %>%
    mutate(id = gsub("GOBP_", "", id),
           id = factor(id, levels = id))

  write.csv(sig_tb, out_csv, quote = FALSE, row.names = FALSE)
  message("GREAT (", nm, "): ", nrow(sig_tb), " significant pathways saved to ", out_csv)
  sig_tb
}

# ── Run for both directions ───────────────────────────────────────────────────
sig_hf <- run_great_pathway(gr_counts, diff_res, "hf", snakemake@output$great_hf_csv)
sig_mf <- run_great_pathway(gr_counts, diff_res, "mf", snakemake@output$great_mf_csv)

# ── Dot plots ─────────────────────────────────────────────────────────────────
make_great_dotplot <- function(sig_tb, label, plot_path) {
  if (nrow(sig_tb) == 0) {
    message("No pathways to plot for: ", label)
    return(invisible(NULL))
  }
  filtered <- sig_tb %>%
    dplyr::filter(observed_gene_hits >= 5, gene_set_size >= 10) %>%
    dplyr::slice_max(fold_enrichment, n = 30)

  pathway_dotplot(
    pathways  = filtered,
    title     = paste("GREAT Pathways:", label),
    plot_path = plot_path,
    x_col     = "fold_enrichment",
    y_col     = "id",
    color_col = "p_value",
    size_col  = "observed_gene_hits"
  )
}

make_great_dotplot(sig_hf, "High Fatigue (more accessible)",  snakemake@output$dotplot_hf)
make_great_dotplot(sig_mf, "Medium Fatigue (more accessible)", snakemake@output$dotplot_mf)

message("ATAC pathway enrichment complete.")
