# workflow/scripts/rna_deseq2.R
# Snakemake rule: rna_deseq2
# Runs DESeq2 differential expression analysis for High vs Low fatigue.

log_file <- snakemake@log[[1]]
con      <- file(log_file, open = "wt")
sink(con, type = "message")
sink(con, type = "output")

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(readr)
  library(DESeq2)
})

# ── Parameters ────────────────────────────────────────────────────────────────
min_counts  <- snakemake@params$min_counts
min_samples <- snakemake@params$min_samples
deseq_coef  <- snakemake@params$deseq_coef

# ── Helper: strip Ensembl version suffix ──────────────────────────────────────
strip_ensembl_version <- function(ids) sub("\\..*", "", ids)

# ── Load inputs ───────────────────────────────────────────────────────────────
message("Loading inputs...")
ct_clean <- readRDS(snakemake@input$counts_rds)
metadata <- readRDS(snakemake@input$metadata_rds)
G_list   <- readRDS(snakemake@input$glist_rds)

# ── Filter sex chromosome genes ───────────────────────────────────────────────
autosomal <- G_list$ensembl_gene_id[
  !G_list$chromosome_name %in% c("X", "Y", "chrX", "chrY")
]

# ── Build DESeq2 dataset ──────────────────────────────────────────────────────
message("Building DESeq2 dataset...")
dds <- DESeqDataSetFromMatrix(
  countData = ct_clean[, rownames(metadata)],
  colData   = metadata,
  design    = ~ group + Seks
)

# Filter lowly expressed genes
keep <- rowSums(counts(dds) >= min_counts) >= min_samples
dds  <- dds[keep, ]

# Remove sex chromosome genes
dds <- dds[rownames(dds) %in% autosomal, ]

# ── Run DESeq2 ────────────────────────────────────────────────────────────────
message("Running DESeq2...")
set.seed(123)
dds <- DESeq(dds)
vsd <- vst(dds, blind = TRUE)

res <- results(dds, name = deseq_coef)
res_df <- as.data.frame(res) %>%
  arrange(padj) %>%
  mutate(ensembl_gene_id = strip_ensembl_version(rownames(.))) %>%
  left_join(G_list, by = "ensembl_gene_id") %>%
  dplyr::filter(!is.na(hgnc_symbol))

# ── Save outputs ──────────────────────────────────────────────────────────────
saveRDS(dds,    snakemake@output$dds_rds)
saveRDS(vsd,    snakemake@output$vsd_rds)
write_csv(res_df, snakemake@output$results_csv)

n_sig <- sum(res_df$padj < snakemake@params$padj_threshold, na.rm = TRUE)
message("DESeq2 complete. ", n_sig, " genes significant at padj < ",
        snakemake@params$padj_threshold)
