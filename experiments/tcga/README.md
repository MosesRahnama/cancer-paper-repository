# TCGA Validation: Hypothesis 1

Empirical validation of Hypothesis 1 from the boundary-logic framework using TCGA RNA-seq data.

## Data Source

- **ISB Cancer Genomics Cloud** BigQuery: `isb-cgc-bq.TCGA.RNAseq_hg38_gdc_current`
- **Cohorts**: TCGA-SKCM (melanoma, 473 samples) and TCGA-LUAD (lung adenocarcinoma, 589 samples)
- **Expression**: TPM (transcripts per million) and raw STAR counts

## Target Genes

| Category | Genes |
|---|---|
| Immune checkpoint | CD274 (PD-L1), PDCD1LG2 (PD-L2), PDCD1 (PD-1) |
| MHC-I | HLA-A, HLA-B, HLA-C, B2M |
| Gap junctions | GJA1 (Cx43), GJB2 (Cx26), GJA5 (Cx40), GJB6 (Cx30) |
| Circadian clock | ARNTL (BMAL1), CLOCK, PER1, PER2, CRY1, CRY2 |
| Differentiation | MYC, TP53, CDH1 (E-cadherin), VIM (Vimentin) |

## Scripts

- `tcga_extract.py` — Queries BigQuery and exports expression matrices to GCS
- `tcga_analysis.py` — Runs Spearman correlations and generates figures

## Usage

```bash
# Requires: google-cloud-bigquery, google-cloud-storage, pandas, scipy, matplotlib
# Set GOOGLE_APPLICATION_CREDENTIALS before running

python experiments/tcga/tcga_extract.py    # Extract data
python experiments/tcga/tcga_analysis.py   # Run analysis
```

## Key Results

Circadian coherence (CV across clock genes) negatively correlates with PD-L1 in both cohorts (SKCM: rho=-0.381, p=9.6e-18; LUAD: rho=-0.303, p=1.1e-12), supporting the boundary-framework prediction that immune checkpoint engagement and temporal decoherence are linked.
