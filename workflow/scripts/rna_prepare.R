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
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  
})

cols_to_remove    <- snakemake@params$cols_to_remove
# ── Helper: strip Ensembl version suffix ──────────────────────────────────────
strip_ensembl_version <- function(ids) sub("\\..*", "", ids)

# ── Helper: read and reformat featureCounts table ─────────────────────────────
read_count_table <- function(path, sample_ids) {
  ct <- read.table(path, sep = "\t", header = TRUE, check.names = FALSE)
  
  # drop HC from the requested sample list
  #sample_ids <- sample_ids[!grepl("^HC", sample_ids)]
  
  # extract sample id from columns like output/rnaseq/bam/CD1/CD1_filtered.bam
  is_bam <- grepl("_filtered\\.bam$", colnames(ct))
  ids <- sub(".*[/.]([^/.]+)_filtered\\.bam$", "\\1", colnames(ct))
  colnames(ct)[is_bam] <- ids[is_bam]
  
  # sanity check before subsetting
  missing <- setdiff(sample_ids, colnames(ct))
  if (length(missing)) stop("No count columns for: ", paste(missing, collapse = ", "))
  
  meta_cols <- c("Geneid", "Chr", "Start", "End", "Strand", "Length")
  ct_filtered <- ct[, c(meta_cols, sample_ids)]
  rownames(ct_filtered) <- ct_filtered$Geneid
  ct_filtered <- ct_filtered[, sample_ids]
  ct_filtered
}

# ── Load and filter metadata ──────────────────────────────────────────────────
message("Loading metadata...")
metadata <- read_excel(snakemake@input$metadata_file) %>%
  dplyr::filter(Type == "HC" | batch != "no_ATACseq", !Patient_ID %in% cols_to_remove) %>%
  dplyr::mutate(
    MFI_scaled = as.numeric(scale(MFI_totalscore)),
    sex        = as.factor(sex),
    group      = factor(group),
    Status     = factor(Status)
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
gene_ids <- rownames(ct_clean)

G_list <- data.frame(
  ensembl_gene_id = gene_ids,
  hgnc_symbol = unname(mapIds(
    org.Hs.eg.db,
    keys = gene_ids,
    keytype = "ENSEMBL",
    column = "SYMBOL",
    multiVals = "first"
  )),
  chromosome_name = unname(mapIds(
    org.Hs.eg.db,
    keys = gene_ids,
    keytype = "ENSEMBL",
    column = "CHR",
    multiVals = "first"
  )),
  entrezgene_id = unname(mapIds(
    org.Hs.eg.db,
    keys = gene_ids,
    keytype = "ENSEMBL",
    column = "ENTREZID",
    multiVals = "first"
  )),
  stringsAsFactors = FALSE
)

message(
  "Fetched annotations for ",
  sum(!is.na(G_list$hgnc_symbol)),
  " / ",
  length(gene_ids),
  " genes."
)
message("Fetched annotations for ", nrow(G_list), " genes.")

# ── Save outputs ──────────────────────────────────────────────────────────────
saveRDS(ct_clean,  snakemake@output$counts_rds)
saveRDS(metadata,  snakemake@output$metadata_rds)
saveRDS(G_list,    snakemake@output$glist_rds)
message("Preparation complete.")
