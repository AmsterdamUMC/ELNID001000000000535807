# workflow/rules/atac.smk
# ──────────────────────────────────────────────────────────────────────────────
# ATAC-seq rules:
#   atac_merge_counts  →  atac_deseq2  →  atac_annotate  →  atac_pathways
# ──────────────────────────────────────────────────────────────────────────────

ATAC_OUT = config["output_dir"] + "/atac"

# ── Step 1: Merge overlapping peak count tables from two batches ───────────────
rule atac_merge_counts:
    input:
        file1 = config["atac"]["count_file1"],
        file2 = config["atac"]["count_file2"],
    output:
        counts_rds = ATAC_OUT + "/merged_counts.rds",
        gr_rds     = ATAC_OUT + "/merged_gr.rds",
    params:
        count_columns     = config["atac_count_columns"],
        new_count_columns = config["atac_new_count_columns"],
        cols_to_remove    = config["atac_cols_to_remove"],
    log:
        ATAC_OUT + "/logs/atac_merge_counts.log",
    conda:
        "/net/beegfs/groups/tytgat/ELNID001000000000535807/analysis/workflow/envs/atac.yaml"
    script:
        "../scripts/atac_merge_counts.R"


# ── Step 2: DESeq2 differential accessibility ─────────────────────────────────
rule atac_deseq2:
    input:
        counts_rds    = ATAC_OUT + "/merged_counts.rds",
        metadata_file = config["atac"]["metadata_file"],
    output:
        dds_rds           = ATAC_OUT + "/dds.rds",
        vsd_rds           = ATAC_OUT + "/vsd.rds",
        differential_csv  = ATAC_OUT + "/differential_peaks.csv",
        vsd_counts_csv    = ATAC_OUT + "/vsd_counts.csv",
        pca_png           = ATAC_OUT + "/pca_group.png",
    params:
        min_row_sums  = config["atac_min_row_sums"],
        lfc_threshold = config["lfc_threshold"],
        pval_threshold= config["pval_threshold"],
        padj_threshold= config["padj_threshold"],
        deseq_coef    = "group_High_fatigue_vs_Medium_fatigue",
    log:
        ATAC_OUT + "/logs/atac_deseq2.log",
    conda:
        "/net/beegfs/groups/tytgat/ELNID001000000000535807/analysis/workflow/envs/atac.yaml"
    script:
        "../scripts/atac_deseq2.R"


# ── Step 3: Peak annotation and volcano plot ──────────────────────────────────
rule atac_annotate:
    input:
        dds_rds          = ATAC_OUT + "/dds.rds",
        gr_rds           = ATAC_OUT + "/merged_gr.rds",
        differential_csv = ATAC_OUT + "/differential_peaks.csv",
    output:
        annotated_csv = ATAC_OUT + "/annotated_peaks.csv",
        volcano_png   = ATAC_OUT + "/volcano_plot.png",
    params:
        lfc_threshold  = config["lfc_threshold"],
        pval_threshold = config["pval_threshold"],
        padj_threshold = config["padj_threshold"],
    log:
        ATAC_OUT + "/logs/atac_annotate.log",
    conda:
        "/net/beegfs/groups/tytgat/ELNID001000000000535807/analysis/workflow/envs/atac.yaml"
    script:
        "../scripts/atac_annotate.R"


# ── Step 4: GREAT pathway enrichment (High fatigue & Medium fatigue) ──────────
rule atac_pathways:
    input:
        gr_rds           = ATAC_OUT + "/merged_gr.rds",
        differential_csv = ATAC_OUT + "/differential_peaks.csv",
    output:
        great_hf_csv  = ATAC_OUT + "/GREATpathwayEnrichment_hf.csv",
        great_mf_csv  = ATAC_OUT + "/GREATpathwayEnrichment_mf.csv",
        dotplot_hf    = ATAC_OUT + "/GREATpathways_hf.png",
        dotplot_mf    = ATAC_OUT + "/GREATpathways_mf.png",
    params:
        great_padj_binom = config["great_padj_binom"],
        great_padj_hyper = config["great_padj_hyper"],
        lfc_threshold    = config["lfc_threshold"],
        pval_threshold   = config["pval_threshold"],
    log:
        ATAC_OUT + "/logs/atac_pathways.log",
    conda:
        "/net/beegfs/groups/tytgat/ELNID001000000000535807/analysis/workflow/envs/atac.yaml"
    script:
        "../scripts/atac_pathways.R"
