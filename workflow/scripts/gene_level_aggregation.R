#!/usr/bin/env Rscript
# Aggregate peak-level DESeq2 results to gene-level calls using csaw::combineTests()
# (Simes combination of peak p-values) + inverse-variance weighted mean logFC for
# an accompanying effect-size summary.
#
# Restricted to promoter peaks: distal intergenic peaks have unreliable nearest-gene
# assignment and should not be pooled into a gene-level call.
#
# Usage (from Snakemake or CLI):
#   Rscript gene_level_aggregation.R <annotated_peaks.csv> <comparison_name> <output.csv>
suppressPackageStartupMessages({
  library(csaw)
  library(dplyr)
  library(S4Vectors)
})

peaks_path      <- snakemake@input$annotated_peaks
comparison_name <- snakemake@params$comparison
out_path        <- snakemake@output$aggragation_csv

stopifnot(file.exists(peaks_path))

peaks <- read.csv(peaks_path, stringsAsFactors = FALSE)

required_cols <- c("annotation", "SYMBOL", "log2FoldChange", "lfcSE", "pvalue")
stopifnot(all(required_cols %in% colnames(peaks)))

is_promoter <- grepl("^Promoter", peaks$annotation)
has_symbol  <- !is.na(peaks$SYMBOL) & nzchar(peaks$SYMBOL)
peaks <- peaks[which(is_promoter & has_symbol), ]

assert_complete <- function(df, cols) {
  for (col in cols) {
    if (any(is.na(df[[col]]))) {
      stop(sprintf("NA values found in column '%s' after filtering", col))
    }
  }
}
print(peaks %>% dplyr::filter(is.na(pvalue)))
peaks <- peaks %>% dplyr::filter(!is.na(pvalue))
assert_complete(peaks, c("log2FoldChange", "lfcSE", "pvalue", "SYMBOL"))

ids <- factor(peaks$SYMBOL)

tab <- DataFrame(
  PValue = peaks$pvalue,
  logFC  = peaks$log2FoldChange
)

# Simes combination of peak-level p-values into one p-value per gene.
combined <- combineTests(ids, tab, pval.col = "PValue", fc.col = "logFC")
combined <- as.data.frame(combined)
combined$SYMBOL <- rownames(combined)

# Inverse-variance weighted mean logFC per gene: combineTests() reports
# direction/significance but not a single effect-size number, so this
# pairs with it for the dot-matrix color scale.
wmean <- peaks %>%
  group_by(SYMBOL) %>%
  summarise(
    n_peaks              = n(),
    wmean_log2FC         = weighted.mean(log2FoldChange, w = 1 / (lfcSE^2)),
    frac_same_direction  = mean(sign(log2FoldChange) == sign(weighted.mean(log2FoldChange, w = 1 / (lfcSE^2)))),
    .groups = "drop"
  )

gene_level <- combined %>%
  left_join(wmean, by = "SYMBOL") %>%
  mutate(comparison = comparison_name) %>%
  select(SYMBOL, comparison, n_peaks, num.tests, direction, PValue, FDR,
         wmean_log2FC, frac_same_direction)

write.csv(gene_level, out_path, row.names = FALSE)

cat(sprintf("[%s] %d genes tested, %d at FDR < 0.05, %d at nominal p < 0.01\n",
            comparison_name, nrow(gene_level),
            sum(gene_level$FDR < 0.05, na.rm = TRUE),
            sum(gene_level$PValue < 0.01, na.rm = TRUE)))