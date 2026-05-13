# workflow/scripts/rna_prepare.R
# Snakemake rule: rna_prepare
# Loads featureCounts table, filters metadata, fetches BioMart gene annotations.

log_file <- snakemake@log[[1]]
con      <- file(log_file, open = "wt")
sink(con, type = "message")
sink(con, type = "output")

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(readxl)
  library(biomaRt)
})

# ── Parameters ────────────────────────────────────────────────────────────────
min_mfi_score   <- snakemake@params$min_mfi_score    # 40
high_mfi_cutoff <- snakemake@params$high_mfi_cutoff  # 80
low_mfi_cutoff  <- snakemake@params$low_mfi_cutoff   # 60

# ── Helper: strip Ensembl version suffix ──────────────────────────────────────
strip_ensembl_version <- function(ids) sub("\\..*", "", ids)

# ── Helper: read and reformat featureCounts table ─────────────────────────────
read_count_table <- function(path, sample_ids) {
  ct <- read.table(path, sep = "\t", header = TRUE)
  samples_regex <- paste0("\\b(", paste(sample_ids, collapse = "|"), ")\\b")
  ct_filtered   <- ct[, !grepl("output", colnames(ct)) |
                         grepl(samples_regex, colnames(ct))]
  colnames(ct_filtered)[7:ncol(ct_filtered)] <- sample_ids
  rownames(ct_filtered) <- ct_filtered$Geneid
  ct_filtered[1:6]      <- NULL
  ct_filtered           <- ct_filtered[, sample_ids]
  ct_filtered
}

# ── Load and filter metadata ──────────────────────────────────────────────────
message("Loading metadata...")
metadata <- read_excel(snakemake@input$metadata_file) %>%
  dplyr::filter(Type == "CD", MFI_totalscore >= min_mfi_score) %>%
  dplyr::filter(MFI_totalscore >= high_mfi_cutoff | MFI_totalscore < low_mfi_cutoff) %>%
  dplyr::mutate(
    MFI_scaled = as.numeric(scale(MFI_totalscore)),
    Seks       = as.factor(Seks),
    group      = factor(
      ifelse(MFI_totalscore < high_mfi_cutoff, "Low_fatigue", "High_fatigue"),
      levels = c("Low_fatigue", "High_fatigue")
    )
  ) %>%
  tibble::column_to_rownames("Patient_ID")

message("Retained ", nrow(metadata), " samples after metadata filtering.")

# ── Load count table ──────────────────────────────────────────────────────────
message("Reading count table...")
ct       <- read_count_table(snakemake@input$counts_raw, rownames(metadata))
ct_clean <- ct
rownames(ct_clean) <- strip_ensembl_version(rownames(ct_clean))

# ── Fetch BioMart annotations ─────────────────────────────────────────────────
message("Connecting to Ensembl BioMart...")
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

G_list <- getBM(
  filters    = "ensembl_gene_id",
  attributes = c("ensembl_gene_id", "hgnc_symbol", "chromosome_name", "entrezgene_id"),
  values     = rownames(ct_clean),
  mart       = mart
)
message("Fetched annotations for ", nrow(G_list), " genes.")

# ── Save outputs ──────────────────────────────────────────────────────────────
saveRDS(ct_clean,  snakemake@output$counts_rds)
saveRDS(metadata,  snakemake@output$metadata_rds)
saveRDS(G_list,    snakemake@output$glist_rds)
message("Preparation complete.")
