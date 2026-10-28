log_file <- snakemake@log[[1]]
con      <- file(log_file, open = "wt")
sink(con, type = "message")
sink(con, type = "output")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
})

# ── Parameters ────────────────────────────────────────────────────────────────
lfc_thr    <- snakemake@params$lfc_thr
pval_thr   <- snakemake@params$pval_thr
genome     <- snakemake@params$genome       # e.g. "hg38"
up_dir     <- snakemake@output$up_dir
down_dir   <- snakemake@output$down_dir

dir.create(up_dir,   recursive = TRUE, showWarnings = FALSE)
dir.create(down_dir, recursive = TRUE, showWarnings = FALSE)

# ── Load differential peaks ───────────────────────────────────────────────────
message("Loading differential peaks...")
diff_res <- read_csv(snakemake@input$differential_csv, show_col_types = FALSE)

# ── Helper: filter peaks and write BED file ───────────────────────────────────
write_bed <- function(diff_res, direction, lfc_thr, pval_thr, bed_path) {
  filtered <- diff_res %>%
    dplyr::filter(pvalue < pval_thr, abs(log2FoldChange) > lfc_thr) %>%
    tidyr::separate(rowname, into = c("chr", "start", "end"),
                    sep = "_", remove = FALSE) %>%
    dplyr::mutate(start = as.integer(start), end = as.integer(end))
  
  filtered <- if (direction == "up") {
    dplyr::filter(filtered, log2FoldChange > 0)
  } else {
    dplyr::filter(filtered, log2FoldChange < 0)
  }
  
  if (nrow(filtered) == 0) {
    message("No significant peaks for direction: ", direction)
    return(invisible(FALSE))
  }
  
  bed <- filtered %>%
    dplyr::select(chr, start, end) %>%
    dplyr::mutate(
      name  = paste0(chr, "_", start, "_", end),
      score = 0,
      strand = "."
    )
  
  write.table(bed, bed_path,
              sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE)
  message(nrow(bed), " peaks written to ", bed_path)
  return(invisible(TRUE))
}

# ── Helper: run HOMER findMotifsGenome.pl ─────────────────────────────────────
run_homer <- function(bed_path, out_dir, genome) {
  if (!file.exists(bed_path) || file.size(bed_path) == 0) {
    message("Skipping HOMER for empty BED: ", bed_path)
    return(invisible(NULL))
  }
  cmd <- paste(
    "findMotifsGenome.pl",
    bed_path,
    genome,
    out_dir,
    "-size given",
    "-mask",
    "-p", snakemake@threads
  )
  message("Running HOMER: ", cmd)
  ret <- system(cmd)
  if (ret != 0) stop("HOMER failed with exit code: ", ret)
}

# ── Write BED files ───────────────────────────────────────────────────────────
up_bed   <- file.path(up_dir,   "peaks_up.bed")
down_bed <- file.path(down_dir, "peaks_down.bed")

write_bed(diff_res, "up",   lfc_thr, pval_thr, up_bed)
write_bed(diff_res, "down", lfc_thr, pval_thr, down_bed)

# ── Run HOMER ─────────────────────────────────────────────────────────────────
run_homer(up_bed,   up_dir,   genome)
run_homer(down_bed, down_dir, genome)

message("HOMER motif analysis complete.")