# workflow/rules/rnaseq.smk
# ──────────────────────────────────────────────────────────────────────────────
# RNA-seq rules:
#   rna_prepare  →  rna_deseq2  →  rna_visualize  →  rna_pathways
# ──────────────────────────────────────────────────────────────────────────────

RNA_OUT = config["output_dir"] + "/rna"

# ── Step 1: Load counts + metadata, fetch BioMart annotations ─────────────────
rule rna_prepare:
    input:
        counts_raw    = config["rna"]["counts_raw"],
        metadata_file = config["rna"]["metadata_file"],
    output:
        counts_rds   = RNA_OUT + "/counts_clean.rds",
        metadata_rds = RNA_OUT + "/metadata.rds",
        glist_rds    = RNA_OUT + "/gene_annotations.rds",
    params:
        min_mfi_score   = 40,   # MFI_totalscore >= 40 (filter in metadata)
        high_mfi_cutoff = 80,   # >= 80 → High_fatigue
        low_mfi_cutoff  = 60,   # < 60  → Low_fatigue
    log:
        RNA_OUT + "/logs/rna_prepare.log",
    conda:
        "/net/beegfs/groups/tytgat/ELNID001000000000535807/analysis/workflow/envs/rna.yaml"
    script:
        "../scripts/rna_prepare.R"


# ── Step 2: DESeq2 differential expression ────────────────────────────────────
rule rna_deseq2:
    input:
        counts_rds   = RNA_OUT + "/counts_clean.rds",
        metadata_rds = RNA_OUT + "/metadata.rds",
        glist_rds    = RNA_OUT + "/gene_annotations.rds",
    output:
        dds_rds    = RNA_OUT + "/dds.rds",
        vsd_rds    = RNA_OUT + "/vsd.rds",
        results_csv= RNA_OUT + "/differentialAnalysis_RNAseq_HighVSLow.csv",
    params:
        min_counts    = config["rna_min_counts"],
        min_samples   = config["rna_min_samples"],
        lfc_threshold = config["lfc_threshold"],
        pval_threshold= config["pval_threshold"],
        padj_threshold= config["rna_padj_threshold"],
        deseq_coef    = "group_High_fatigue_vs_Low_fatigue",
    log:
        RNA_OUT + "/logs/rna_deseq2.log",
    conda:
        "/net/beegfs/groups/tytgat/ELNID001000000000535807/analysis/workflow/envs/rna.yaml"
    script:
        "../scripts/rna_deseq2.R"


# ── Step 3: QC visualizations (PCA, heatmap, volcano) ────────────────────────
rule rna_visualize:
    input:
        dds_rds      = RNA_OUT + "/dds.rds",
        vsd_rds      = RNA_OUT + "/vsd.rds",
        results_csv  = RNA_OUT + "/differentialAnalysis_RNAseq_HighVSLow.csv",
        glist_rds    = RNA_OUT + "/gene_annotations.rds",
        metadata_rds = RNA_OUT + "/metadata.rds",
    output:
        pca_png     = RNA_OUT + "/pca.png",
        heatmap_png = RNA_OUT + "/heatmap.png",
        volcano_png = RNA_OUT + "/volcano_plot.png",
    params:
        top_var_peaks  = config["rna_top_var_peaks"],
        lfc_threshold  = config["lfc_threshold"],
        pval_threshold = config["pval_threshold"],
        padj_threshold = config["rna_padj_threshold"],
    log:
        RNA_OUT + "/logs/rna_visualize.log",
    conda:
        "/net/beegfs/groups/tytgat/ELNID001000000000535807/analysis/workflow/envs/rna.yaml"
    script:
        "../scripts/rna_visualize.R"


# ── Step 4: GSEA pathway enrichment ───────────────────────────────────────────
rule rna_pathways:
    input:
        results_csv = RNA_OUT + "/differentialAnalysis_RNAseq_HighVSLow.csv",
    output:
        fgsea_csv         = RNA_OUT + "/enrichedPathways_HighFatigue.csv",
        dotplot_top_png   = RNA_OUT + "/pathwayEnrichment_topPathways.png",
        dotplot_sel_png   = RNA_OUT + "/pathwayEnrichment_selectedPathways.png",
    params:
        pval_threshold      = config["pval_threshold"],
        padj_threshold      = config["rna_padj_threshold"],
        fgsea_min_size      = config["fgsea_min_size"],
        fgsea_max_size      = config["fgsea_max_size"],
        fgsea_nes_threshold = config["fgsea_nes_threshold"],
        go_pattern          = config["go_pattern"],
    log:
        RNA_OUT + "/logs/rna_pathways.log",
    conda:
        "/net/beegfs/groups/tytgat/ELNID001000000000535807/analysis/workflow/envs/rna.yaml"
    script:
        "../scripts/rna_pathways.R"
