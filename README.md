# Fatigue Study Multi-Omics Snakemake Pipeline - ELNID001000000000535807

Modular Snakemake pipeline for ATAC-seq differential accessibility and
RNA-seq differential expression analysis comparing fatigue groups.

## Directory structure

```
ELNID001000000000535807/analysis/workflow/
├── config/
│   ├── config.yaml          # All paths and thresholds — edit this first
│   └── atac_metadata.csv    # ATAC-seq sample metadata (6 samples used in DESeq2)
├── Snakefile            # Master file; declares rule all
├── rules/
│   ├── atac.smk         # ATAC-seq rules (4 steps)
│   └── rnaseq.smk       # RNA-seq rules (4 steps)
├── envs/
│   ├── base.yaml         # Basic env
│   ├── atac.yaml         # ATAC-seq env 
│   └── rna.yaml       # RNA-seq env 
├── scripts/
    ├── utils.R                # Shared plotting helpers (volcano, PCA, dotplot)
    ├── atac_merge_counts.R    # Step 1: merge two-batch peak count tables
    ├── atac_deseq2.R          # Step 2: DESeq2 differential accessibility
    ├── atac_annotate.R        # Step 3: ChIPseeker annotation + volcano
    ├── atac_pathways.R        # Step 4: GREAT GO:BP enrichment
    ├── rna_prepare.R          # Step 1: load counts + BioMart annotations
    ├── rna_deseq2.R           # Step 2: DESeq2 differential expression
    ├── rna_visualize.R        # Step 3: PCA, heatmap, volcano
    └── rna_pathways.R         # Step 4: fgsea GO:BP enrichment
└── results/                 # All outputs written here (created automatically)
```

## Quickstart

### 1. Edit config
Open `config/config.yaml` and update the `atac` and `rna` input paths to
point to your data files. Thresholds can be adjusted there too.

### 2. Dry run (check the DAG without executing)
```bash
snakemake -n --snakefile workflow/Snakefile
```

### 3. Run locally
```bash
# if needed, run this first:
# snakemake --use-conda --conda-create-envs-only --cores 2 --conda-frontend conda
snakemake --snakefile workflow/Snakefile --cores 4 --conda-frontend conda --use-conda 
```

### 4. Run on a cluster (SLURM example)
```bash
snakemake --snakefile workflow/Snakefile \
  --executor slurm \
  --default-resources slurm_partition=normal mem_mb=32000 \
  --jobs 20
```

## Pipeline DAG

```
ATAC-seq
  atac_merge_counts  →  atac_deseq2  →  atac_annotate
                                      →  atac_pathways

RNA-seq
  rna_prepare  →  rna_deseq2  →  rna_visualize
                              →  rna_pathways
```

## Key outputs

| File | Description |
|------|-------------|
| `results/atac/differential_peaks.csv` | DESeq2 results (apeglm LFC shrinkage) |
| `results/atac/annotated_peaks.csv` | Peaks annotated with ChIPseeker |
| `results/atac/GREATpathwayEnrichment_hf.csv` | GREAT results, High fatigue peaks |
| `results/atac/GREATpathwayEnrichment_mf.csv` | GREAT results, Medium fatigue peaks |
| `results/rna/differentialAnalysis_RNAseq_HighVSLow.csv` | DESeq2 DE results |
| `results/rna/enrichedPathways_HighFatigue.csv` | Full fgsea results |

## Notes

- The ATAC pipeline drops Low fatigue samples from the DESeq2 contrast
  (High vs Medium). The Low fatigue columns are listed in `atac_cols_to_remove`
  in `config.yaml` and can be reinstated if needed.
- BioMart queries require an internet connection during `rna_prepare`.
- All thresholds (LFC, p-value, padj, GREAT, fgsea) are centralised in
  `config/config.yaml` — no need to edit individual R scripts.
