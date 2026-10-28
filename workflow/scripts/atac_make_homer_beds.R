# workflow/scripts/atac_make_homer_beds.R
# Called by atac_homer_da shell rule.
# Reads DA peaks + merged GR, writes 3 BED files (up, down, background),
# then prints their paths to stdout for the shell to capture.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(GenomicRanges)
})

args     <- commandArgs(trailingOnly = TRUE)
da_csv   <- args[1]
gr_rds   <- args[2]
lfc_thr  <- as.numeric(args[3])
pval_thr <- as.numeric(args[4])
out_dir  <- args[5]   # directory to write BED files into (persistent, not /tmp)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

da <- read_csv(da_csv, show_col_types = FALSE)
gr <- readRDS(gr_rds)

# ── Background: all consensus peaks (0-based BED) ────────────────────────────
bg_bed <- file.path(out_dir, "background.bed")
as.data.frame(gr) %>%
  dplyr::select(seqnames, start, end) %>%
  dplyr::mutate(start = start - 1L) %>%
  write.table(bg_bed, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE)

# ── Helper: peak ID "chr:start-end" → 0-based BED ────────────────────────────
write_bed <- function(peak_ids, path) {
  if (length(peak_ids) == 0) {
    file.create(path)
    return(invisible(NULL))
  }
  coords <- do.call(rbind, strsplit(peak_ids, "_"))
  coords <- matrix(coords, ncol = 3)   # ← forces matrix shape even for 1 row
  data.frame(
    chr   = coords[, 1],
    start = as.integer(coords[, 2]) - 1L,
    end   = as.integer(coords[, 3])
  ) %>%
    write.table(path, sep = "\t", quote = FALSE,
                row.names = FALSE, col.names = FALSE)
}

up_peaks   <- da %>%
  dplyr::filter(log2FoldChange >  lfc_thr, pvalue < pval_thr) %>%
  pull(rowname)
down_peaks <- da %>%
  dplyr::filter(log2FoldChange < -lfc_thr, pvalue < pval_thr) %>%
  pull(rowname)

message("Up peaks:   ", length(up_peaks))
message("Down peaks: ", length(down_peaks))


up_bed   <- file.path(out_dir, "peaks_up.bed")
down_bed <- file.path(out_dir, "peaks_down.bed")
write_bed(up_peaks,   up_bed)
write_bed(down_peaks, down_bed)

# Print paths to stdout for shell to capture
cat(up_bed, "\n")
cat(down_bed, "\n")
cat(bg_bed, "\n")
cat(length(up_peaks), "\n")
cat(length(down_peaks), "\n")