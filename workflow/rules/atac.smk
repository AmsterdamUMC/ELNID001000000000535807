# workflow/rules/atac.smk
from pathlib import Path
import re

ATAC_OUT        = config["output_dir"] + "/atac"
ATAC_BAM_DIRS   = config["atac"]["bam_dir"]
ATAC_SAMPLE_MAP = config["atac_sample_map"]          # long_name → short_name
ATAC_SAMPLES    = [s.replace(".bam", "") for s in ATAC_SAMPLE_MAP.keys()]
ATAC_PEAK_DIRS  = config["atac"]["peak_dirs"]        # 3 batch peak dirs
def find_bam(sample, bam_dirs):
    """Search all bam_dirs for {sample}.bam, return first match."""
    for d in bam_dirs:
        p = Path(d) / f"{sample}.bam"
        if p.exists():
            return str(p)
    raise ValueError(f"BAM not found for '{sample}' in: {bam_dirs}")

def find_peaks(peak_dirs):
    """Collect all narrowPeak files across all peak dirs."""
    peaks = []
    for d in peak_dirs:
        peaks += list(Path(d).glob("*.narrowPeak"))
    if not peaks:
        raise ValueError(f"No .narrowPeak files found in: {peak_dirs}")
    return [str(p) for p in peaks]


# ── Step 0a: Build union peak set from all batches ────────────────────────────
rule atac_union_peaks:
    input:
        narrowpeaks = find_peaks(ATAC_PEAK_DIRS),
        blacklist   = config["blacklist"],
    output:
        union_peaks = ATAC_OUT + "/union_peaks.bed",
    log:
        "logs/atac_union_peaks.log",
    conda:
        "/net/beegfs/groups/tytgat/ELNID001000000000535807/analysis/workflow/envs/atac.yaml"
    shell:
        """
        (
        cat {input.narrowpeaks} \
          | grep -v "^#" \
          | cut -f1-3 \
          | sort -k1,1 -k2,2n \
          | bedtools merge -i stdin \
          | bedtools intersect -a stdin -b {input.blacklist} -v \
          > {output.union_peaks}

        echo "Union peaks after blacklist filter: $(wc -l < {output.union_peaks})"
        ) > {log} 2>&1
        """


# ── Step 0b: Count reads per peak per BAM ─────────────────────────────────────
rule atac_count_reads:
    input:
        bam   = lambda wc: find_bam(wc.sample, ATAC_BAM_DIRS),
        bai   = lambda wc: find_bam(wc.sample, ATAC_BAM_DIRS) + ".bai",
        peaks = ATAC_OUT + "/union_peaks.bed",
    output:
        counts = ATAC_OUT + "/sample_counts/{sample}_counts.bed",
    resources:
        mem_mb = 190000
    log:
        "logs/atac_count_{sample}.log",
    conda:
        "/net/beegfs/groups/tytgat/ELNID001000000000535807/analysis/workflow/envs/atac.yaml"
    shell:
        """
        (
        bedtools coverage \
          -a {input.peaks} \
          -b {input.bam} \
          -counts \
          > {output.counts}

        echo "Counted $(wc -l < {output.counts}) peaks for {wildcards.sample}"
        ) > {log} 2>&1
        """

# ── Cell type deconvolution (EPIC-ATAC) ──────────────────────────────────────
rule atac_celltype:
    input:
        counts_rds    = ATAC_OUT + "/merged_counts.rds",
        metadata_file = config["metadata_file"],
    output:
        fractions_csv        = ATAC_OUT + "/{comparison}/celltype_fractions.csv",
        fractions_no_neut_csv= ATAC_OUT + "/{comparison}/celltype_fractions_no_neutrophils.csv",
        plot_png             = ATAC_OUT + "/{comparison}/celltype_proportions.png",
        plot_no_neut_png     = ATAC_OUT + "/{comparison}/celltype_proportions_no_neutrophils.png",
    params:
        comparison      = lambda wc: wc.comparison,
        subset_type     = lambda wc: COMPARISONS[wc.comparison]["subset"],
        subset_status   = lambda wc: COMPARISONS[wc.comparison]["subset_status"],
        continuous      = lambda wc: COMPARISONS[wc.comparison]["continuous"],
        status_contrast = lambda wc: COMPARISONS[wc.comparison]["status_contrast"],
        cols_to_remove  = config["atac_samples_to_remove"],
        hc_comparison  = lambda wc: COMPARISONS[wc.comparison].get("hc_comparison", False),
        hc_fatigue_level = lambda wc: COMPARISONS[wc.comparison].get("hc_fatigue_level", "Low"),
    log:
        "logs/atac_celltype_{comparison}.log",
    conda:
        "/net/beegfs/groups/tytgat/ELNID001000000000535807/analysis/workflow/envs/epicatac.yaml"
    script:
        "../scripts/atac_celltype.R"
 

# ── Step 1: Assemble count matrix ─────────────────────────────────────────────
rule atac_merge_counts:
    input:
        count_beds = expand(
            ATAC_OUT + "/sample_counts/{sample}_counts.bed",
            sample=ATAC_SAMPLES
        ),
    output:
        counts_rds = ATAC_OUT + "/merged_counts.rds",
        gr_rds     = ATAC_OUT + "/merged_gr.rds",
    params:
        sample_map     = ATAC_SAMPLE_MAP,
        cols_to_remove = config["atac_samples_to_remove"],
    log:
        "logs/atac_merge_counts.log",
    conda:
        "/net/beegfs/groups/tytgat/ELNID001000000000535807/analysis/workflow/envs/atac.yaml"
    script:
        "../scripts/atac_merge_counts.R"

# ── Step 2: DESeq2 differential accessibility ─────────────────────────────────
rule atac_deseq2:
    input:
        counts_rds    = ATAC_OUT + "/merged_counts.rds",
        metadata_file = config["metadata_file"],
        utils         = "scripts/utils.R"
    output:
        dds_rds          = ATAC_OUT + "/{comparison}/dds.rds",
        vsd_rds          = ATAC_OUT + "/{comparison}/vsd.rds",
        design_records   = ATAC_OUT + "/{comparison}/design_record.rds",
        differential_csv = ATAC_OUT + "/{comparison}/differential_peaks.csv",
        vsd_counts_csv   = ATAC_OUT + "/{comparison}/vsd_counts.csv",
        group_pca_png    = ATAC_OUT + "/{comparison}/pca_group.png",
        sex_pca_png      = ATAC_OUT + "/{comparison}/pca_sex.png",
        type_pca_png     = ATAC_OUT + "/{comparison}/pca_type.png",
        status_pca_png   = ATAC_OUT + "/{comparison}/pca_status.png",
    params:
        cols_to_remove = config["atac_samples_to_remove"],
        min_row_sums   = config["atac_min_row_sums"],
        padj_threshold = config["padj_threshold"],
        lfc_threshold  = lambda wc: config["continuous_lfc_threshold"] if COMPARISONS[wc.comparison]["continuous"] else config["lfc_threshold"],
        deseq_coef     = lambda wc: COMPARISONS[wc.comparison]["coef"],
        subset_type    = lambda wc: COMPARISONS[wc.comparison]["subset"],
        subset_status  = lambda wc: COMPARISONS[wc.comparison]["subset_status"],
        continuous     = lambda wc: COMPARISONS[wc.comparison]["continuous"],
        status_contrast  = lambda wc: COMPARISONS[wc.comparison]["status_contrast"],
        loess_norm     = config["atac_loess_norm"],
        hc_comparison  = lambda wc: COMPARISONS[wc.comparison].get("hc_comparison", False),
        hc_fatigue_level = lambda wc: COMPARISONS[wc.comparison].get("hc_fatigue_level", "Low"),
        frip_adjust      = lambda wc: COMPARISONS[wc.comparison].get("frip_adjust"),
        adjust_for_group = lambda wc: COMPARISONS[wc.comparison].get("adjust_for_group"),
        pool_types       = lambda wc: COMPARISONS[wc.comparison].get("pool_types"),

    log:
        "logs/atac_deseq2_{comparison}.log",
    conda:
        "/net/beegfs/groups/tytgat/ELNID001000000000535807/analysis/workflow/envs/atac.yaml"
    script:
        "../scripts/atac_deseq2.R"




# ── Step 2b: HOMER known motif enrichment on DA peaks ────────────────────────
rule atac_homer_da:
    input:
        differential_csv = ATAC_OUT + "/{comparison}/differential_peaks.csv",
        gr_rds           = ATAC_OUT + "/merged_gr.rds",
        beds_script      = "scripts/atac_make_homer_beds.R",
    output:
        homer_up_txt   = ATAC_OUT + "/{comparison}/homer_up/knownResults.txt",
        homer_down_txt = ATAC_OUT + "/{comparison}/homer_down/knownResults.txt",
    params:
        lfc_threshold  = lambda wc: config["continuous_lfc_threshold"] if COMPARISONS[wc.comparison]["continuous"] else config["lfc_threshold"],
        pval_threshold = config["pval_threshold"],
        genome         = "hg38",
        beds_dir       = ATAC_OUT + "/{comparison}/homer_beds",
        out_up         = ATAC_OUT + "/{comparison}/homer_up",
        out_down       = ATAC_OUT + "/{comparison}/homer_down",
    threads: 4
    log:
        "logs/atac_homer_da_{comparison}.log",
    conda:
        "/net/beegfs/groups/tytgat/ELNID001000000000535807/analysis/workflow/envs/homer.yaml"
    shell:
        """
        set -euo pipefail
        # Install hg38 into the conda env's HOMER if not already present
        if [ ! -d "$CONDA_PREFIX/share/homer/data/genomes/hg38" ]; then
          echo "Installing hg38 into HOMER..." >> {log}
          cd $CONDA_PREFIX/share/homer && \
          perl configureHomer.pl -install hg38 >> {log} 2>&1
          cd -
        fi
        # Step 1: extract BED files via dedicated R script
        mkdir -p {params.beds_dir}
        Rscript {input.beds_script} \\
            {input.differential_csv} \\
            {input.gr_rds} \\
            {params.lfc_threshold} \\
            {params.pval_threshold} \\
            {params.beds_dir} \\
            > {params.beds_dir}/bed_paths.txt 2>> {log}

        up_bed=$(sed -n '1p' {params.beds_dir}/bed_paths.txt | tr -d '\\n')
        down_bed=$(sed -n '2p' {params.beds_dir}/bed_paths.txt | tr -d '\\n')
        bg_bed=$(sed -n '3p' {params.beds_dir}/bed_paths.txt | tr -d '\\n')
        n_up=$(sed -n '4p' {params.beds_dir}/bed_paths.txt | tr -d '\\n')
        n_down=$(sed -n '5p' {params.beds_dir}/bed_paths.txt | tr -d '\\n')

        echo "Up peaks: $n_up   Down peaks: $n_down" >> {log}

        mkdir -p {params.out_up} {params.out_down}

        # Step 2: run HOMER on up peaks (skip if empty)
        if [ "$n_up" -gt 0 ]; then
            findMotifsGenome.pl \\
                $up_bed {params.genome} {params.out_up} \\
                -bg $bg_bed -p {threads} -nomotif >> {log} 2>&1
        else
            echo "No up peaks — creating empty knownResults.txt" >> {log}
            touch {output.homer_up_txt}
        fi

        # Step 3: run HOMER on down peaks (skip if empty)
        if [ "$n_down" -gt 0 ]; then
            findMotifsGenome.pl \\
                $down_bed {params.genome} {params.out_down} \\
                -bg $bg_bed -p {threads} -nomotif >> {log} 2>&1
        else
            echo "No down peaks — creating empty knownResults.txt" >> {log}
            touch {output.homer_down_txt}
        fi
        """

# ── Step 3: Peak annotation and volcano plot ──────────────────────────────────
rule atac_annotate:
    input:
        dds_rds          = ATAC_OUT + "/{comparison}/dds.rds",
        gr_rds           = ATAC_OUT + "/merged_gr.rds",
        differential_csv = ATAC_OUT + "/{comparison}/differential_peaks.csv",
    output:
        annotated_csv = ATAC_OUT + "/{comparison}/annotated_peaks.csv",
        volcano_png   = ATAC_OUT + "/{comparison}/volcano_plot.png",
    params:
        lfc_threshold  = lambda wc: config["continuous_lfc_threshold"] if COMPARISONS[wc.comparison]["continuous"] else config["lfc_threshold"],
        pval_threshold = config["pval_threshold"],
        padj_threshold = config["padj_threshold"],
        comparison     = lambda wc: wc.comparison,
        deseq_coef     = lambda wc: COMPARISONS[wc.comparison]["coef"],
    log:
        "logs/{comparison}_atac_annotate.log",
    conda:
        "/net/beegfs/groups/tytgat/ELNID001000000000535807/analysis/workflow/envs/atac.yaml"
    script:
        "../scripts/atac_annotate.R"


# ── Step 4: GREAT pathway enrichment ─────────────────────────────────────────
rule atac_pathways:
    input:
        gr_rds           = ATAC_OUT + "/merged_gr.rds",
        differential_csv = ATAC_OUT + "/{comparison}/differential_peaks.csv",
    output:
        great_up_csv   = ATAC_OUT + "/{comparison}/GREAT_{database}_up.csv",
        great_down_csv = ATAC_OUT + "/{comparison}/GREAT_{database}_down.csv",
        dotplot_up_png = ATAC_OUT + "/{comparison}/GREAT_{database}_up_dotplot.png",
        dotplot_dn_png = ATAC_OUT + "/{comparison}/GREAT_{database}_down_dotplot.png",
        great_up_padj_csv   = ATAC_OUT + "/{comparison}/GREAT_{database}_up_padj.csv",
        great_down_padj_csv = ATAC_OUT + "/{comparison}/GREAT_{database}_down_padj.csv",
        dotplot_up_padj_png = ATAC_OUT + "/{comparison}/GREAT_{database}_up_dotplot_padj.png",
        dotplot_dn_padj_png = ATAC_OUT + "/{comparison}/GREAT_{database}_down_dotplot_padj.png",
    params:
        script_dir       = "scripts",
        great_padj_binom = config["great_padj_binom"],
        great_padj_hyper = config["great_padj_hyper"],
        lfc_threshold    = lambda wc: config["continuous_lfc_threshold"] if COMPARISONS[wc.comparison]["continuous"] else config["lfc_threshold"],
        pval_threshold   = config["pval_threshold"],
        padj_threshold   = config["padj_threshold"],
        continuous       = lambda wc: COMPARISONS[wc.comparison]["continuous"],
        comparison       = lambda wc: wc.comparison,
        deseq_coef       = lambda wc: COMPARISONS[wc.comparison]["coef"],
        database         = lambda wc: ATAC_DATABASES[wc.database],
    log:
        "logs/{comparison}_atac_pathways{database}.log",
    conda:
        "/net/beegfs/groups/tytgat/ELNID001000000000535807/analysis/workflow/envs/atac.yaml"
    script:
        "../scripts/atac_pathways.R"


# ── Step 5: Motif ORA using pre-computed HOMER knownResults.txt ───────────────
rule atac_motif_ora:
    input:
        homer_up_txt   = ATAC_OUT + "/{comparison}/homer_up/knownResults.txt",
        homer_down_txt = ATAC_OUT + "/{comparison}/homer_down/knownResults.txt",
    output:
        ora_up_csv       = ATAC_OUT + "/{comparison}/motif_ORA_up.csv",
        ora_down_csv     = ATAC_OUT + "/{comparison}/motif_ORA_down.csv",
        dotplot_up_png   = ATAC_OUT + "/{comparison}/motif_ORA_up_dotplot.png",
        dotplot_down_png = ATAC_OUT + "/{comparison}/motif_ORA_down_dotplot.png",
    params:
        script_dir  = "scripts",
        comparison  = lambda wc: wc.comparison,
        homer_padj  = config["homer_padj_threshold"],
        pval_thr    = config["pval_threshold"],
        min_gs_size = config.get("ora_min_gs_size", 10),
        max_gs_size = config.get("ora_max_gs_size", 500),
    log:
        "logs/atac_motif_ora_{comparison}.log",
    conda:
        "/net/beegfs/groups/tytgat/ELNID001000000000535807/analysis/workflow/envs/motif.yaml"
    script:
        "../scripts/atac_motif_ora.R"
        
rule gene_level_aggregation:
    input:  
      annotated_peaks = ATAC_OUT + "/{comparison}/annotated_peaks.csv"
    output: 
      aggragation_csv = ATAC_OUT + "/{comparison}/gene_level_aggregation.csv"
    params:
        comparison  = "{comparison}",
    log:
        "logs/atac_gene_level_aggregation_{comparison}.log",
    conda:
        "/net/beegfs/groups/tytgat/ELNID001000000000535807/analysis/workflow/envs/atac.yaml"
    script:  
      "../scripts/gene_level_aggregation.R"
