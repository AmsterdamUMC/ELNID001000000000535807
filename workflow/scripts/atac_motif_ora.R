# workflow/scripts/atac_motif_ora.R
# ORA on TFs whose motifs are enriched in DA peaks (up vs down),
# using HOMER knownResults.txt from atac_homer_da rule.
log_file <- snakemake@log[[1]]
con      <- file(log_file, open = "wt")
sink(con, type = "message")
sink(con, type = "output")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

source(file.path(snakemake@params$script_dir, "utils.R"))

# ── Parameters ────────────────────────────────────────────────────────────────
pval_thr    <- snakemake@params$pval_thr
min_gs_size <- snakemake@params$min_gs_size
max_gs_size <- snakemake@params$max_gs_size
comparison  <- snakemake@params$comparison

comparison_label <- switch(comparison,
                           "CD_fatigue"           = "High vs Low Fatigue (CD)",
                           "UC_fatigue"           = "High vs Low Fatigue (UC)",
                           "UC_vs_CD"             = "UC vs CD",
                           comparison
)

if (grepl("UC_vs_CD", comparison)) {
  up_label   <- "UC"
  down_label <- "CD"
} else {
  up_label   <- "High Fatigue"
  down_label <- "Low Fatigue"
}

# ── Helper: parse knownResults.txt and return significant TF names ────────────
parse_known_motifs <- function(path, pval_thr) {
  message("Parsing: ", path)
  tryCatch({
    df <- read_tsv(path, show_col_types = FALSE) %>%
      dplyr::rename(
        motif_name = `Motif Name`,
        pvalue     = `P-value`,
        padj       = `q-value (Benjamini)`
      ) %>%
      dplyr::filter(pvalue < pval_thr) %>%
      dplyr::mutate(TF_name = gsub("[(/].*", "", motif_name))
    
    message(nrow(df), " significant motifs at pval_thr < ", pval_thr)
    unique(df$TF_name)
    print(unique(df$TF_name))
  }, error = function(e) {
    message("Failed to parse ", path, ": ", e$message)
    character(0)
  })
}

# ── Helper: run ORA and save outputs ─────────────────────────────────────────
run_ora <- function(tfs, label, out_csv, out_png) {
  if (length(tfs) == 0) {
    message("Skipping ORA for ", label, " — no significant motifs")
    write.csv(data.frame(), out_csv, row.names = FALSE)
    p <- ggplot() +
      labs(title = paste("No significant motifs:", label)) +
      theme_void() +
      theme(plot.title = element_text(size = 12, hjust = 0.5))
    ggsave(out_png, plot = p, width = 8, height = 4, dpi = 150)
    return(invisible(NULL))
  }
  
  message("Running ORA on ", length(tfs), " TFs for: ", label)
  
  enrich_go <- tryCatch(
    clusterProfiler::enrichGO(
      gene          = tfs,
      OrgDb         = org.Hs.eg.db,
      keyType       = "SYMBOL",
      ont           = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff  = pval_thr,
      minGSSize     = min_gs_size,
      maxGSSize     = max_gs_size,
      readable      = TRUE
    ),
    error = function(e) {
      message("enrichGO failed for ", label, ": ", e$message)
      NULL
    }
  )
  
  if (is.null(enrich_go) || nrow(as.data.frame(enrich_go)) == 0) {
    message("No enriched GO terms for: ", label)
    write.csv(data.frame(), out_csv, row.names = FALSE)
    p <- ggplot() +
      labs(title = paste("No enriched GO terms:", label)) +
      theme_void() +
      theme(plot.title = element_text(size = 12, hjust = 0.5))
    ggsave(out_png, plot = p, width = 8, height = 4, dpi = 150)
    return(invisible(NULL))
  }
  
  write.csv(as.data.frame(enrich_go), out_csv, row.names = FALSE)
  message(nrow(as.data.frame(enrich_go)), " enriched GO terms saved to ", out_csv)
  
  p <- dotplot(enrich_go, showCategory = 20) +
    labs(title = paste("Motif TF GO Enrichment:", label, "-", comparison_label)) +
    theme(plot.title = element_text(size = 12))
  ggsave(out_png, plot = p, width = 10, height = 8, dpi = 300)
  message("Dot plot saved to ", out_png)
}

# ── Parse HOMER results and run ORA ──────────────────────────────────────────
tfs_up   <- parse_known_motifs(snakemake@input$homer_up_txt,   pval_thr)
tfs_down <- parse_known_motifs(snakemake@input$homer_down_txt, pval_thr)

run_ora(tfs_up,   up_label,   snakemake@output$ora_up_csv,   snakemake@output$dotplot_up_png)
run_ora(tfs_down, down_label, snakemake@output$ora_down_csv, snakemake@output$dotplot_down_png)

message("Motif ORA complete.")