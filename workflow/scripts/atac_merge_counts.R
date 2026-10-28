# workflow/scripts/atac_merge_counts.R
log_file <- snakemake@log[[1]]
con      <- file(log_file, open = "wt")
sink(con, type = "message")
sink(con, type = "output")

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(GenomicRanges)
  library(tibble)
})

# ── Parameters ────────────────────────────────────────────────────────────────
count_beds     <- snakemake@input$count_beds
sample_map     <- snakemake@params$sample_map     # long_name → short_name
cols_to_remove <- snakemake@params$cols_to_remove

# ── Load per-sample count BED files ───────────────────────────────────────────
message("Reading ", length(count_beds), " count files...")

count_list <- map(count_beds, function(path) {
  #long_name <- sub("\\.bam$", "", sub("_counts\\.bed$", "", basename(path)))
  long_name <- sub("_counts\\.bed$", "", basename(path))
  df <- read.table(path, header = FALSE, sep = "\t",
                   col.names = c("chr", "start", "end", long_name))
  df
})

# ── Assemble matrix ───────────────────────────────────────────────────────────
message("Assembling count matrix...")

coords       <- count_list[[1]][, c("chr", "start", "end")]
count_matrix <- bind_cols(coords, map(count_list, ~ .x[, 4, drop = FALSE]))

# ── Rename: long BAM stem → short sample ID ───────────────────────────────────
message("Renaming columns...")

# R mangles column names when reading: digits->X prefix, hyphens->dots
# Apply the same mangling to sample_map keys so they match
names(sample_map) <- sub("\\.bam$", "", names(sample_map))  # strip .bam
names(sample_map) <- make.names(names(sample_map))           # apply same R mangling

# Verify keys now match column names
unmatched <- setdiff(names(count_matrix)[-(1:3)], names(sample_map))
if (length(unmatched) > 0) {
  warning("These columns had no match in sample_map: ", 
          paste(unmatched, collapse=", "))
}

colnames(count_matrix) <- dplyr::recode(colnames(count_matrix), !!!sample_map)

message("Columns after renaming: ", paste(colnames(count_matrix), collapse=", "))

# ── Remove unwanted samples (uses short names) ────────────────────────────────
if (length(cols_to_remove) > 0) {
  message("Removing: ", paste(cols_to_remove, collapse = ", "))
  count_matrix <- count_matrix %>% dplyr::select(-any_of(cols_to_remove))
}

# ── Sort sample columns alphabetically ────────────────────────────────────────
coord_cols  <- c("chr", "start", "end")
sample_cols <- sort(setdiff(colnames(count_matrix), coord_cols))
count_matrix <- count_matrix[, c(coord_cols, sample_cols)]

# ── Build GRanges ─────────────────────────────────────────────────────────────
message("Building GRanges object...")
gr_union <- makeGRangesFromDataFrame(count_matrix, keep.extra.columns = TRUE)

# ── Build count table ─────────────────────────────────────────────────────────
message("Building count table...")
full_ct <- count_matrix %>%
  mutate(region_id = paste0(chr, "_", start, "_", end)) %>%
  column_to_rownames("region_id") %>%
  dplyr::select(-chr, -start, -end)

# ── Save ──────────────────────────────────────────────────────────────────────
message("Saving outputs...")
saveRDS(full_ct,  snakemake@output$counts_rds)
saveRDS(gr_union, snakemake@output$gr_rds)
message("Done. ", nrow(full_ct), " peaks x ", ncol(full_ct), " samples.")