# workflow/scripts/atac_merge_counts.R
# Snakemake rule: atac_merge_counts
# Merges two batch-specific ATAC-seq peak count tables into a single
# GRanges-based count matrix by reducing overlapping peaks.

log_file <- snakemake@log[[1]]
con      <- file(log_file, open = "wt")
sink(con, type = "message")
sink(con, type = "output")

suppressPackageStartupMessages({
  library(dplyr)
  library(GenomicRanges)
  library(tibble)
  library(readr)
})

# ── Parameters from config ────────────────────────────────────────────────────
count_columns     <- snakemake@params$count_columns
new_count_columns <- snakemake@params$new_count_columns
cols_to_remove    <- snakemake@params$cols_to_remove

# ── Helper functions ──────────────────────────────────────────────────────────
aggregate_counts <- function(revmap, mcols_df, cols) {
  sapply(cols, function(col) {
    col_counts <- mcols_df[[col]]
    vapply(revmap, function(ix) sum(col_counts[ix], na.rm = TRUE), numeric(1))
  })
}

rename_and_filter_mcols <- function(mcols_df) {
  colnames(mcols_df) <- new_count_columns
  mcols_df <- as.data.frame(mcols_df) %>%
    dplyr::select(-dplyr::any_of(cols_to_remove))
  mcols_df <- mcols_df[, order(colnames(mcols_df))]
  mcols_df
}

# ── Main ──────────────────────────────────────────────────────────────────────
message("Reading count tables...")
ct1 <- read.table(snakemake@input$file1, sep = "\t", header = TRUE)
ct2 <- read.table(snakemake@input$file2, sep = "\t", header = TRUE)

gr1 <- makeGRangesFromDataFrame(ct1, keep.extra.columns = TRUE)
gr2 <- makeGRangesFromDataFrame(ct2, keep.extra.columns = TRUE)
gr_all <- c(gr1, gr2)

message("Reducing overlapping peaks...")
gr_union <- GenomicRanges::reduce(gr_all, with.revmap = TRUE)

all_mcols <- as.data.frame(mcols(gr_all))
agg       <- aggregate_counts(gr_union$revmap, all_mcols, count_columns)
mcols(gr_union)[count_columns] <- agg
mcols(gr_union)$revmap         <- NULL
mcols(gr_union)                <- rename_and_filter_mcols(mcols(gr_union))

message("Building full count table...")
full_ct <- as.data.frame(gr_union) %>%
  mutate(region_id = paste0(seqnames, "_", start, "_", end)) %>%
  column_to_rownames("region_id") %>%
  dplyr::select(-seqnames, -start, -end, -width, -strand)

message("Saving outputs...")
saveRDS(full_ct,   snakemake@output$counts_rds)
saveRDS(gr_union,  snakemake@output$gr_rds)
message("Done. Merged ", nrow(full_ct), " regions.")
