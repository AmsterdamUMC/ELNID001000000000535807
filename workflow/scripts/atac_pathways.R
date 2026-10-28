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

source(file.path(snakemake@params$script_dir, "utils.R"))

# ── Parameters ────────────────────────────────────────────────────────────────
great_padj_binom <- snakemake@params$great_padj_binom
great_padj_hyper <- snakemake@params$great_padj_hyper
lfc_thr          <- snakemake@params$lfc_threshold
pval_thr         <- snakemake@params$pval_threshold
padj_thr         <- snakemake@params$padj_threshold
comparison       <- snakemake@params$comparison
database         <- snakemake@params$database

# Direction labels depend on comparison type
if (comparison == "UC_vs_CD") {
  up_label   <- "UC"
  down_label <- "CD"
} else {
  up_label   <- "High Fatigue"
  down_label <- "Low Fatigue"
}

# ── Load inputs ───────────────────────────────────────────────────────────────
message("Loading inputs for comparison: ", comparison)
gr_counts <- readRDS(snakemake@input$gr_rds)
diff_res  <- read_csv(snakemake@input$differential_csv, show_col_types = FALSE)

# ── Helper: run GREAT for one direction ──────────────────────────────────────
run_great_pathway <- function(gr_counts, diff_res, direction, label, out_csv, mode = c("lfc_pval", "padj")) {
  mode <- match.arg(mode)
  
  keep <- if (mode == "padj") {
    which(!is.na(diff_res$padj) & diff_res$padj < padj_thr)
  } else {
    which(!is.na(diff_res$pvalue) &
            abs(diff_res$log2FoldChange) > lfc_thr &
            diff_res$pvalue < pval_thr)
  }
  
  filtered <- diff_res[keep, ] %>%
    tidyr::separate(rowname, into = c("chr", "start", "end"),
                    sep = "_", remove = FALSE) %>%
    dplyr::mutate(start = as.integer(start), end = as.integer(end))
  
  filtered <- if (direction == "up") {
    dplyr::filter(filtered, log2FoldChange > 0)
  } else {
    dplyr::filter(filtered, log2FoldChange < 0)
  }
  
  if (nrow(filtered) == 0) {
    message("No significant peaks for: ", label, " [", mode, "]")
    write.csv(data.frame(), out_csv, row.names = FALSE)
    return(invisible(data.frame()))
  }
  
  diff_gr <- makeGRangesFromDataFrame(
    filtered, seqnames.field = "chr",
    start.field = "start", end.field = "end", ignore.strand = TRUE
  )
  seqlevelsStyle(diff_gr) <- seqlevelsStyle(gr_counts)
  
  res   <- great(diff_gr, background = gr_counts,
                 database, "TxDb.Hsapiens.UCSC.hg38.knownGene")
  tb    <- getEnrichmentTable(res)
  tb <- tb %>%dplyr::arrange(desc(fold_enrichment)) %>%
    mutate(
      id = gsub("GOBP_", "", id), id = gsub("KEGG_", "", id),
      id = factor(id, levels = id)
    )
  
  sig_tb <- tb %>%
    dplyr::filter(p_adjust <= great_padj_binom, p_adjust_hyper <= great_padj_hyper) %>%
    dplyr::arrange(desc(fold_enrichment)) %>%
    mutate(
      id = gsub("GOBP_", "", id), id = gsub("KEGG_", "", id),
      id = factor(id, levels = id)
    )
  
  write.csv(tb, out_csv, quote = FALSE, row.names = FALSE)
  message("GREAT (", label, "): ", nrow(sig_tb), " significant pathways.")
  sig_tb
}

# ── Run for both directions ───────────────────────────────────────────────────
sig_up   <- run_great_pathway(gr_counts, diff_res, "up",   up_label,   snakemake@output$great_up_csv)
sig_down <- run_great_pathway(gr_counts, diff_res, "down", down_label, snakemake@output$great_down_csv)

# ── Dot plots ────────────────────────────────────────────────────────────────
make_great_dotplot <- function(sig_tb, label, plot_path, mode) {
  sig_lab <- if (mode == "padj") paste0("FDR<", padj_thr)
  else paste0("p<", pval_thr, ", |LFC|>", lfc_thr)
  if (nrow(sig_tb) == 0) {
    message("No pathways to plot for: ", plot_path, " — writing empty placeholder")
    p <- ggplot() +
      labs(title = paste("No significant pathways:", label)) +
      theme_void() +
      theme(plot.title = element_text(size = 12, hjust = 0.5))
    ggsave(plot_path, plot = p, width = 8, height = 4, dpi = 150)
    return(invisible(NULL))
  }
  filtered <- sig_tb %>%
    dplyr::filter(observed_gene_hits >= 5, gene_set_size >= 10) %>%
    dplyr::slice_max(fold_enrichment, n = 30)
  
  pathway_dotplot(
    pathways  = filtered,
    title     = paste("GREAT ", database, "Pathways:", label, "-", gsub("_", " ",
                                                          sub("^[^_]+_", "", comparison))),
    plot_path = plot_path,
    x_col     = "fold_enrichment",
    y_col     = "id",
    color_col = "p_value",
    size_col  = "observed_gene_hits"
  )
}

for (m in c("lfc_pval", "padj")) {
  sfx <- if (m == "padj") "_padj" else ""
  up   <- run_great_pathway(gr_counts, diff_res, "up",   up_label,
                            snakemake@output[[paste0("great_up",   sfx, "_csv")]], m)
  down <- run_great_pathway(gr_counts, diff_res, "down", down_label,
                            snakemake@output[[paste0("great_down", sfx, "_csv")]], m)
  make_great_dotplot(up,   up_label,   snakemake@output[[paste0("dotplot_up",  sfx, "_png")]], m)
  make_great_dotplot(down, down_label, snakemake@output[[paste0("dotplot_dn",  sfx, "_png")]], m)
}

message("ATAC pathway enrichment complete.")