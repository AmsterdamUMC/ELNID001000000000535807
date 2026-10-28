# workflow/rules/rnaseq.smk
# ──────────────────────────────────────────────────────────────────────────────
# RNA-seq rules:
#   rna_prepare  →  rna_deseq2  →  rna_visualize  →  rna_pathways
# ──────────────────────────────────────────────────────────────────────────────

RNA_OUT = config["output_dir"] + "/rna"

rule rna_prepare:
    input:
        counts_raw    = config["rna"]["counts_raw"],
        metadata_file = config["metadata_file"],
    output:
        counts_rds   = RNA_OUT + "/counts_clean.rds",
        metadata_rds = RNA_OUT + "/metadata.rds",
        glist_rds    = RNA_OUT + "/gene_annotations.rds",
    params:
      cols_to_remove    = config["rna_samples_to_remove"],
          
    log:
        "logs/rna_prepare.log",
    conda:
        "/net/beegfs/groups/tytgat/ELNID001000000000535807/analysis/workflow/envs/rna.yaml"
    script:
        "../scripts/rna_prepare.R"


rule rna_deseq2:
    input:
        counts_rds   = RNA_OUT + "/counts_clean.rds",
        metadata_rds = RNA_OUT + "/metadata.rds",
        glist_rds    = RNA_OUT + "/gene_annotations.rds",
    output:
        dds_rds     = RNA_OUT + "/{comparison}/dds.rds",
        vsd_rds     = RNA_OUT + "/{comparison}/vsd.rds",
        metadata    = RNA_OUT + "/{comparison}/metadata.rds",
        results_csv = RNA_OUT + "/{comparison}/differential_expression.csv",
        design_records   = RNA_OUT + "/{comparison}/design_record.rds",
    params:
        script_dir     = "scripts",
        min_counts     = config["rna_min_counts"],
        min_samples    = config["rna_min_samples"],
        lfc_threshold  = lambda wc: config["continuous_lfc_threshold"] if COMPARISONS[wc.comparison]["continuous"] else config["lfc_threshold"],
        padj_threshold = config["rna_padj_threshold"],
        deseq_coef     = lambda wc: COMPARISONS[wc.comparison]["coef"],
        subset_type    = lambda wc: COMPARISONS[wc.comparison]["subset"],
        subset_status  = lambda wc: COMPARISONS[wc.comparison]["subset_status"],  # "Active", "Remission", or None
        continuous     = lambda wc: COMPARISONS[wc.comparison]["continuous"],
        status_contrast  = lambda wc: COMPARISONS[wc.comparison]["status_contrast"],
        hc_comparison = lambda w: COMPARISONS[w.comparison].get("hc_comparison", False),
        adjust_for_group = lambda wc: COMPARISONS[wc.comparison].get("adjust_for_group"),
        pool_types       = lambda wc: COMPARISONS[wc.comparison].get("pool_types"),
    log:
        "logs/rna_deseq2_{comparison}.log",
    conda:
        "/net/beegfs/groups/tytgat/ELNID001000000000535807/analysis/workflow/envs/rna.yaml"
    script:
        "../scripts/rna_deseq2.R"


rule rna_visualize:
    input:
        dds_rds     = RNA_OUT + "/{comparison}/dds.rds",
        vsd_rds     = RNA_OUT + "/{comparison}/vsd.rds",
        results_csv = RNA_OUT + "/{comparison}/differential_expression.csv",
        glist_rds   = RNA_OUT + "/gene_annotations.rds",
        metadata_rds= RNA_OUT + "/{comparison}/metadata.rds",
    output:
        pca_png     = RNA_OUT + "/{comparison}/pca.png",
        heatmap_png = RNA_OUT + "/{comparison}/heatmap.png",
        volcano_png = RNA_OUT + "/{comparison}/volcano_plot.png",
    params:
        script_dir      = "scripts",
        top_var_peaks   = config["rna_top_var_peaks"],
        lfc_threshold  = lambda wc: config["continuous_lfc_threshold"] if COMPARISONS[wc.comparison]["continuous"] else config["lfc_threshold"],
        pval_threshold  = config["pval_threshold"],
        padj_threshold  = config["padj_threshold"],
        deseq_coef      = lambda wc: COMPARISONS[wc.comparison]["coef"],
        comparison      = lambda wc: wc.comparison,
        grouping_factor = lambda wc: COMPARISONS[wc.comparison]["grouping_factor"],
    log:
        "logs/rna_visualize_{comparison}.log",
    conda:
        "/net/beegfs/groups/tytgat/ELNID001000000000535807/analysis/workflow/envs/rna.yaml"
    script:
        "../scripts/rna_visualize.R"


rule rna_pathways:
    input:
        results_csv = RNA_OUT + "/{comparison}/differential_expression.csv",
        utils       = "../scripts/utils.R"
    output:
        fgsea_csv       = RNA_OUT + "/{comparison}/enrichedPathways_{database}.csv",
        dotplot_top_png = RNA_OUT + "/{comparison}/pathwayEnrichment_{database}_topPathways.png",
        #dotplot_sel_png = RNA_OUT + "/{comparison}/pathwayEnrichment_selectedPathways.png",
    params:
        script_dir          = "scripts",
        padj_threshold      = config["rna_padj_threshold"],
        fgsea_min_size      = config["fgsea_min_size"],
        fgsea_max_size      = config["fgsea_max_size"],
        fgsea_nes_threshold = config["fgsea_nes_threshold"],
        go_pattern          = config["go_pattern"],
        comparison          = lambda wc: wc.comparison,
        database            = lambda wc: RNA_DATABASES[wc.database],
    log:
        "logs/rna_pathways_{comparison}_{database}.log",
    conda:
        "/net/beegfs/groups/tytgat/ELNID001000000000535807/analysis/workflow/envs/rna.yaml"
    script:
        "../scripts/rna_pathways.R"
