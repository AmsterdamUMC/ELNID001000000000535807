# workflow/scripts/utils.R
# Shared plotting helpers used by all pipeline scripts.
# Sourced explicitly by each script; does NOT use global CONFIG.
# All thresholds must be passed as arguments.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(RColorBrewer)
  library(vegan)
})

#' Volcano plot: differential accessibility or expression
#' @param df          Results data frame (log2FoldChange, pvalue, padj, gene column)
#' @param gene_col    Name of the gene/symbol column in df
#' @param title       Axis label: "Accessibility" or "Expression"
#' @param plot_path   Output PNG path
#' @param lfc_thr     |log2FC| threshold (default 1)
#' @param pval_thr    p-value threshold (default 0.05)
#' @param padj_thr    adjusted p-value threshold for point shape (default 0.05)
volcano_plot <- function(df, gene_col, title, plot_path, comparison_title,
                         lfc_thr  = 1.0,
                         pval_thr = 0.05,
                         padj_thr = 0.05) {
  df <- df %>%
    mutate(
      StatSig = case_when(
        log2FoldChange >  lfc_thr & pvalue < pval_thr ~ "Upregulated",
        log2FoldChange < -lfc_thr & pvalue < pval_thr ~ "Downregulated",
        TRUE ~ "NS"
      ),
      StatSig = factor(StatSig, levels = c("Upregulated", "Downregulated", "NS")),
      shape   = factor(
        dplyr::case_when(
          is.na(padj)        ~ "NS",
          padj < padj_thr    ~ "Significant",
          TRUE               ~ "NS"
        ),
        levels = c("Significant", "NS")
      ),
      delabel = if_else(
        padj < padj_thr & StatSig != "NS" & !is.na(.data[[gene_col]]),
        .data[[gene_col]], NA_character_
      )
    )
  p <- ggplot(df, aes(x = log2FoldChange, y = -log10(pvalue),
                      color = StatSig, label = delabel, shape = shape)) +
    geom_hline(yintercept = -log10(pval_thr), color = "gray50", linetype = "dashed") +
    geom_vline(xintercept = c(-lfc_thr, lfc_thr), color = "gray50", linetype = "dashed") +
    geom_point(size = 2, alpha = 0.7) +
    geom_text_repel(max.overlaps = Inf, size = 3) +
    scale_color_manual(
      values = c("Upregulated" = "#00AFBB", "Downregulated" = "#BB0C00", "NS" = "gray70"),
      name   = title
    ) +
    scale_shape_manual(
      values = c("Significant" = 16, "NS" = 1),
      name   = paste0("padj < ", padj_thr)
    ) +
    labs(
      x     = "log\u2082(fold change)",
      y     = expression(-Log[10]("p-value")),
      title = paste0("Differential ", title, ": ", comparison_title)
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position    = "top",
          panel.grid.minor   = element_blank(),
           plot.title         = element_text(hjust = 0.5, face = "bold")) #+
    # coord_cartesian(xlim = c(-5, 5))
  
  ggsave(plot_path, p, width = 10, height = 7, dpi = 300, bg = "white")
  message("Volcano plot saved to: ", plot_path)
}

#' PCA plot from a DESeq2 VST object
#' @param vsd_counts  DESeq2 VST SummarizedExperiment
#' @param color_by    colData column to colour by
#' @param plot_path   Output PNG path
plot_pca <- function(vsd_counts, color_by, plot_path) {
  
  meta <- as.data.frame(colData(vsd_counts))
  
  # check grouping
  if (is.null(meta[[color_by]]) || length(unique(meta[[color_by]])) < 2) {
    
    message("Skipping PCA: <2 groups in ", color_by)
    
    p <- ggplot() +
      theme_void() +
      annotate("text", x = 0, y = 0,
               label = paste0("No variation in ", color_by))
    
    ggsave(plot_path, p, width = 6, height = 5, units = "in")
    return(invisible(NULL))
  }
  
  mat     <- assay(vsd_counts)
  pca     <- prcomp(t(mat), scale. = TRUE)
  pcaData <- data.frame(pca$x, meta)
  
  pct <- round(100 * summary(pca)$importance[2, 1:2], 1)
  
  form <- as.formula(paste("dist(t(mat)) ~", color_by))
  
  adonis_res <- adonis2(form, data = meta)
  
  pval <- signif(adonis_res$`Pr(>F)`[1], 3)
  r2   <- round(adonis_res$R2[1], 3)
  
  p <- ggplot(pcaData, aes(PC1, PC2, color = !!sym(color_by))) +
    geom_point(size = 3) +
    stat_ellipse(aes(group = !!sym(color_by)),
                 type = "norm", linewidth = 1, alpha = 0.2) +
    geom_text_repel(aes(label = rownames(pcaData)), size = 3) +
    labs(
      x = paste0("PC1 (", pct[1], "%)"),
      y = paste0("PC2 (", pct[2], "%)"),
      subtitle = paste0("PERMANOVA: R² = ", r2, ", p = ", pval)
    ) +
    scale_color_brewer(palette = "Set2") +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  ggsave(plot_path, p, width = 6, height = 5, units = "in")
  message("PCA plot saved to: ", plot_path)
}
#' Dot plot for pathway enrichment results (fgsea or GREAT)
#' @param pathways  Data frame of enrichment results
#' @param title     Plot title
#' @param plot_path Output PNG path
#' @param x_col     Column for x-axis (e.g. fold_enrichment or NES)
#' @param y_col     Column for y-axis labels
#' @param color_col Column for colour aesthetic
#' @param size_col  Column for size aesthetic
pathway_dotplot <- function(pathways, title, plot_path,
                            x_col     = "fold_enrichment",
                            y_col     = "id",
                            color_col = "p_value",
                            size_col  = "observed_gene_hits") {
  if (nrow(pathways) == 0) {
    message("No pathways to plot for: ", plot_path, " — writing empty placeholder")
    p <- ggplot() +
      labs(title = paste("No significant pathways:", title)) +
      theme_void() +
      theme(plot.title = element_text(size = 12, hjust = 0.5))
    ggsave(plot_path, plot = p, width = 8, height = 4, dpi = 150)
    return(invisible(NULL))
  }
  pathways[[y_col]] <- factor(pathways[[y_col]], levels = pathways[[y_col]][order(pathways[[x_col]], decreasing = FALSE)])
  p <- ggplot(pathways, aes(
    x     = .data[[x_col]],
    y     = .data[[y_col]],
    color = .data[[color_col]],
    size  = .data[[size_col]]
  )) +
    geom_point() +
    labs(title = title, x = x_col, y = NULL) +
    theme_bw() +
    theme(axis.text.x  = element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y  = element_text(size = 9),
          axis.title   = element_text(size = 12),
          legend.title = element_text(size = 10),
          legend.text  = element_text(size = 8),
          plot.title   = element_text(size = 14, hjust = 0.5))
  
  plot_height <- max(6, nrow(pathways) * 0.15)
  ggsave(plot_path, plot = p, dpi = 400, width = 12, height = plot_height, limitsize = FALSE)
  message("Pathway plot saved to: ", plot_path)
}