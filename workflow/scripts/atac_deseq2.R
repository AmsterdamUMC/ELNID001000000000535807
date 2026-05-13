# workflow/scripts/atac_deseq2.R
# Snakemake rule: atac_deseq2
# Runs DESeq2 differential accessibility analysis on merged ATAC-seq counts.

log_file <- snakemake@log[[1]]
con      <- file(log_file, open = "wt")
sink(con, type = "message")
sink(con, type = "output")

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(readr)
  library(DESeq2)
  library(apeglm)
})

source(file.path(snakemake@scriptdir, "utils.R"))

# ── Parameters ────────────────────────────────────────────────────────────────
min_row_sums  <- snakemake@params$min_row_sums
deseq_coef    <- snakemake@params$deseq_coef

# ── Load inputs ───────────────────────────────────────────────────────────────
message("Loading merged count table...")
full_counts   <- readRDS(snakemake@input$counts_rds)
metadata      <- read.csv(snakemake@input$metadata_file, row.names = 1)

# Ensure factor levels are set correctly
metadata$group <- factor(metadata$group, levels = c("Medium_fatigue", "High_fatigue"))
metadata$batch <- factor(metadata$batch)
metadata$sex   <- factor(metadata$sex)
metadata$MFI_scaled <- scale(metadata$MFI)

# ── Build DESeq2 dataset ──────────────────────────────────────────────────────
message("Building DESeq2 dataset...")
dds <- DESeqDataSetFromMatrix(
  countData = full_counts[, rownames(metadata)],
  colData   = metadata,
  design    = ~ group + batch
)

# Filter low-count peaks and sex chromosomes
keep <- rowSums(counts(dds)) >= min_row_sums
dds  <- dds[keep, ]
dds  <- dds[!grepl("chr[XY]", rownames(dds)), ]

# ── Run DESeq2 ────────────────────────────────────────────────────────────────
message("Running DESeq2...")
set.seed(123)
dds <- DESeq(dds, fitType = "local", minReplicatesForReplace = Inf)
vsd <- vst(dds, blind = TRUE)

res <- lfcShrink(dds, coef = deseq_coef, type = "apeglm")

# Filter random/unmapped chromosomes
filtered_res <- as.data.frame(res) %>%
  tibble::rownames_to_column("rowname") %>%
  dplyr::filter(!grepl("chrUn|random", rowname, ignore.case = TRUE))

# ── Outputs ───────────────────────────────────────────────────────────────────
message("Saving outputs...")
saveRDS(dds, snakemake@output$dds_rds)
saveRDS(vsd, snakemake@output$vsd_rds)
write_csv(filtered_res, snakemake@output$differential_csv)

vsd_df <- as.data.frame(assay(vsd)) %>% rownames_to_column("region")
write_csv(vsd_df, snakemake@output$vsd_counts_csv)

plot_pca(vsd, "group", snakemake@output$pca_png)

message("DESeq2 complete. ", sum(filtered_res$padj < snakemake@params$padj_threshold, na.rm = TRUE),
        " peaks significant at padj < ", snakemake@params$padj_threshold)
