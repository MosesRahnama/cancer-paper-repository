# TCGA Cancer Boundary-Logic Analysis

Complete reproducible pipeline for extracting and analyzing TCGA cancer genomics data to validate the boundary-logic failure hypothesis.

## Overview

This analysis validates that cancer progression involves systematic failures in cellular boundary mechanisms—specifically temporal coherence (circadian regulation), immune distinction (MHC-I presentation), and communication (gap junctions). The key finding is that boundary failure is not monolithic: tumors exhibit distinct failure modes with different clinical outcomes.

**Key Results:**
- **n=3,611 tumor samples** across 6 cancer types (plus 309 normal samples; 3,920 total rows)
- **100% replication**: Circadian CV is negatively associated with PD-L1 across all 6 cancer types (all FDR q < 0.05)
- **Two distinct boundary-failure modes identified:**
  - **Active Masking**: High PD-L1 + locked circadian clock (high BMAL1, low PER1) → better survival
  - **Decoherence**: Low PD-L1 + low MHC-I → worse survival, genuine loss of coherence
- **Survival association (patient-level deduplicated; BH across all 12 survival tests)**: Active Masking vs Decoherence is significant in SKCM (p=0.0011, q=0.0132) and directionally consistent but nominal in LUAD (p=0.049, q=0.294)

## Prerequisites

### 1. Google Cloud Platform Setup

You need access to the ISB-CGC BigQuery datasets. This requires:

1. **GCP Account** with BigQuery API enabled
2. **Service Account Credentials**:
   - Create a service account in your GCP project
   - Grant `BigQuery Data Viewer` and `BigQuery Job User` roles
   - Download JSON credentials file
   - Place at `infra/secrets/gcp-credentials.json` (relative to repo root)
   - OR set environment variable: `export GOOGLE_APPLICATION_CREDENTIALS=/path/to/credentials.json`

3. **GCS Bucket** (optional, for cloud backup):
   - Update `PROJECT_ID` and `GCS_BUCKET` in [tcga_config.py](tcga_config.py) if using your own bucket

### 2. Python Environment

**Required packages:**
```bash
pip install google-cloud-bigquery google-cloud-storage pandas numpy scipy matplotlib seaborn
```

**Tested versions:**
- Python 3.10+
- google-cloud-bigquery >= 3.11.0
- pandas >= 2.0.0
- scipy >= 1.11.0 (required for `scipy.stats.ecdf` and `CensoredData`)
- matplotlib >= 3.7.0
- seaborn >= 0.12.0

## Data Sources

All data is extracted from **ISB-CGC BigQuery public datasets** hosted at `isb-cgc-bq.TCGA.*`:

1. **RNAseq_hg38_gdc_current**: RNA-seq expression data (TPM values) for 21 target genes
2. **clinical_gdc_current**: Patient demographics and vital status
3. **clinical_diagnosis_gdc_current**: Tumor staging and diagnosis details

### Target Genes (n=21)

**Immune checkpoints** (3 genes):
- CD274 (PD-L1), PDCD1LG2 (PD-L2), PDCD1 (PD-1)

**MHC Class I** (4 genes):
- HLA-A, HLA-B, HLA-C, B2M (β2-microglobulin)

**Gap junctions** (4 genes):
- GJA1 (Cx43), GJB2 (Cx26), GJA5 (Cx40), GJB6 (Cx30)

**Circadian clock** (6 genes):
- ARNTL (BMAL1), CLOCK, PER1, PER2, CRY1, CRY2

**Differentiation/control** (4 genes):
- MYC, TP53, CDH1 (E-cadherin), VIM (Vimentin)

### Cancer Types (n=6)

| Project ID  | Cancer Type                      | Tumor Samples | Normal Samples |
|-------------|----------------------------------|---------------|----------------|
| TCGA-SKCM   | Skin Cutaneous Melanoma          | 472           | 1              |
| TCGA-LUAD   | Lung Adenocarcinoma              | 530           | 59             |
| TCGA-BRCA   | Breast Invasive Carcinoma        | 1113          | 113            |
| TCGA-COAD   | Colon Adenocarcinoma             | 473           | 41             |
| TCGA-HNSC   | Head and Neck Squamous Cell      | 522           | 44             |
| TCGA-LUSC   | Lung Squamous Cell Carcinoma     | 503           | 51             |
| **Total**   |                                  | **3611**      | **309**        |

## Pipeline Execution

### Quick Start (Full Pipeline)

```bash
cd experiments/tcga

# 1. Extract data from BigQuery (~5 minutes)
python tcga_extract_expanded.py

# 2. Run all core analyses (each ~1-2 minutes)
python tcga_multicancer.py         # Multi-cancer correlation validation
python tcga_tumor_normal.py        # Tumor vs normal comparison
python tcga_immune_subtype.py      # Active Masking vs Decoherence classification
python tcga_survival.py            # Kaplan-Meier survival analysis
python tcga_stage_analysis.py      # Stage-stratified analysis + FDR correction

# 3. OPTIONAL: Run robustness checks with covariate control (for peer review)
# Requires optional covariate file (see script #8 documentation)
python robustness_check.py         # Multivariable Cox models with stage + microenvironment controls

# 4. Sync manuscript artifacts + emit provenance manifest
python sync_manuscript_artifacts.py --rewrite-tex

# All results saved to current directory + figures/ subdirectory
```

### Detailed Script Documentation

---

#### 1. **tcga_config.py**
**Type**: Configuration module (imported by all other scripts)
**Purpose**: Shared configuration, gene lists, and utility functions
**Key functions**:
- `get_bq_client()`: Initialize BigQuery client with credentials
- `compute_circadian_cv(df)`: Calculate circadian coherence (CV across 6 clock genes)
- `classify_boundary_failure(df)`: Classify tumors into Active_Masking/Decoherence/Mixed
- `normalize_stage(stage_str)`: Normalize AJCC stage labels to I/II/III/IV

**No execution required** (library module only)

---

#### 2. **tcga_extract_expanded.py**
**Type**: Data extraction
**Runtime**: ~5 minutes
**BigQuery cost**: ~$0.10 (processes ~50 GB)
**Outputs**:
- `tcga_expanded_tpm.csv` (3,920 samples × 25 columns: metadata + 21 genes)
- `tcga_clinical.csv` (3,646 patients × 12 columns: demographics, vital status, staging)

**What it does:**
1. Queries `TCGA.RNAseq_hg38_gdc_current` for 6 cancer types
2. Pivots from long format (gene × sample rows) to wide format (sample rows, gene columns)
3. Queries `TCGA.clinical_gdc_current` + `TCGA.clinical_diagnosis_gdc_current` for survival and staging
4. Deduplicates clinical records (keeps first diagnosis per patient)
5. Uploads to GCS bucket (optional, if configured)

**Key BigQuery query pattern:**
```sql
SELECT
  project_short_name,
  case_barcode,
  sample_barcode,
  sample_type_name,
  MAX(IF(gene_name = 'CD274', tpm_unstranded, NULL)) AS CD274,
  MAX(IF(gene_name = 'ARNTL', tpm_unstranded, NULL)) AS ARNTL,
  ...
FROM `isb-cgc-bq.TCGA.RNAseq_hg38_gdc_current`
WHERE project_short_name IN ('TCGA-SKCM', 'TCGA-LUAD', ...)
  AND gene_name IN ('CD274', 'ARNTL', 'PER1', ...)
GROUP BY project_short_name, case_barcode, sample_barcode, sample_type_name
```

**Run this first** before any analysis scripts.

---

#### 3. **tcga_multicancer.py**
**Type**: Hypothesis 1 validation (multi-cancer replication)
**Runtime**: ~1 minute
**Inputs**: `tcga_expanded_tpm.csv`
**Outputs**:
- `tcga_multicancer_correlations.csv` (correlation results for all 6 cancer types)
- **Figures**:
  - `multicancer_correlation_heatmap.png` (rows=comparisons, cols=cancer types)
  - `multicancer_circadian_cv_forest.png` (forest plot: circadian CV vs PD-L1)

**What it does:**
1. For each cancer type, computes Spearman correlations:
   - Circadian CV vs Checkpoint genes (CD274, PDCD1LG2)
   - Circadian CV vs MHC-I composite (mean of HLA-A, HLA-B, B2M)
   - Circadian CV vs Gap junction composite
2. Applies Benjamini-Hochberg FDR correction across all tests
3. Generates heatmap showing effect sizes (rho) and significance

**Key finding**: Circadian CV vs PD-L1 shows a consistent **negative** correlation in all 6 cancer types (100% replication; all FDR q < 0.05).

---

#### 4. **tcga_tumor_normal.py**
**Type**: Paired tumor vs normal comparison
**Runtime**: ~1 minute
**Inputs**: `tcga_expanded_tpm.csv`
**Outputs**:
- `tumor_normal_comparison.csv` (statistical test results)
- **Figures**:
  - `tumor_normal_circadian_cv_paired.png` (paired comparison per cancer type)
  - `tumor_normal_circadian_cv_summary.png` (summary boxplot across all types)

**What it does:**
1. Matches tumor and normal samples by case_barcode
2. Compares circadian CV, gap junction expression, PD-L1 expression
3. Uses **Wilcoxon signed-rank test** for paired samples
4. Uses **Mann-Whitney U test** on project-level unpaired tumor vs normal distributions (independent-sample check)

**Key finding**: Tumors have **lower** circadian CV than normals in 4/5 cancer types → clock is "locked" not "broken".

---

#### 5. **tcga_immune_subtype.py**
**Type**: Boundary-failure mode classification
**Runtime**: ~1 minute
**Inputs**: `tcga_expanded_tpm.csv`
**Outputs**:
- `immune_subtype_comparison.csv` (Kruskal-Wallis test results)
- **Figures**:
  - `immune_subtype_circadian_cv.png` (violin plots: circadian CV by subtype)
  - `immune_subtype_boundary_scatter.png` (PD-L1 vs circadian CV scatter, colored by subtype)

**What it does:**
1. Classifies each tumor using within-cancer-type medians:
   - **Active_Masking**: CD274 > median AND ARNTL > median AND PER1 < median
   - **Decoherence**: CD274 < median AND B2M < median
   - **Mixed**: everything else
2. Compares circadian CV, gap junction expression, clock gene ratios across groups
3. Uses **Kruskal-Wallis test** + pairwise Mann-Whitney U with FDR correction

**Key finding**: Active Masking tumors have **significantly different** circadian patterns than Decoherence (all p < 10^-10).

---

#### 6. **tcga_survival.py**
**Type**: Kaplan-Meier survival analysis
**Runtime**: ~2 minutes
**Inputs**: `tcga_expanded_tpm.csv`, `tcga_clinical.csv`
**Outputs**:
- `survival_logrank_results.csv` (log-rank test p-values)
- **Figures**:
  - `survival_circadian_quartile.png` (KM curves stratified by circadian CV quartiles)
  - `survival_boundary_failure.png` (KM curves: Active Masking vs Decoherence vs Mixed)

**What it does:**
1. Filters to tumor samples and retains one tumor sample per case (Primary > Recurrent > Metastatic when duplicates exist)
2. Joins expression data with clinical outcomes (vital_status, days_to_death, days_to_last_follow_up)
3. Computes survival time: days_to_death if dead, else days_to_last_follow_up (censored)
4. Stratifies patients by:
   - Circadian CV quartiles (Q1=most coherent, Q4=least coherent)
   - Boundary-failure mode (Active Masking vs Decoherence)
5. Estimates KM curves using `scipy.stats.ecdf` with `CensoredData.right_censored()`
6. Tests group differences using `scipy.stats.logrank()`

**Key finding**:
- After enforcing one tumor sample per case, Active Masking vs Decoherence is:
- SKCM: log-rank p = 0.0011, BH q = 0.0132 across all 12 survival tests (significant)
- LUAD: log-rank p = 0.049, BH q = 0.294 across all 12 survival tests (nominal only)

**Statistical note**: Uses `CensoredData.right_censored(all_times, is_censored_mask)` where `is_censored_mask[i] = True` if patient i is censored (alive or lost to follow-up).

---

#### 7. **tcga_stage_analysis.py**
**Type**: Stage-stratified analysis + comprehensive FDR correction
**Runtime**: ~1 minute
**Inputs**: `tcga_expanded_tpm.csv`, `tcga_clinical.csv`
**Outputs**:
- `stage_analysis_results.csv` (Spearman correlation: circadian CV vs stage)
- `master_fdr_results.csv` (comprehensive FDR correction across ALL analyses)
- **Figures**:
  - `stage_circadian_cv_boxplot.png` (circadian CV by stage, per cancer type)

**What it does:**
1. Normalizes AJCC stage labels to I/II/III/IV using `normalize_stage()`
2. Tests if circadian CV increases with stage (ordinal Spearman correlation)
3. **Collects ALL p-values** from multicancer, tumor-normal, immune subtype, and survival analyses
4. Applies comprehensive **Benjamini-Hochberg FDR correction** to the full set of tests

**Key finding**: Circadian CV has no FDR-significant association with stage in any cohort; boundary mode appears early and is maintained through progression. (Some non-circadian markers show modest stage associations.)

---

#### 8. **robustness_check.py** (OPTIONAL)
**Type**: Multivariable survival analysis with covariate control
**Runtime**: ~3-5 minutes
**Inputs**:
- `tcga_expanded_tpm.csv` (required)
- `tcga_clinical.csv` (required)
- `tcga_purity_immune_covariates.csv` (OPTIONAL - see below)

**Outputs**:
- `robustness_cox_coefficients.csv` (covariate-level Cox coefficients, HRs, CIs, p-values)
- `robustness_primary_tests.csv` (pre-specified survival robustness tests with BH q-values)
- `robustness_ph_diagnostics.csv` (proportional hazards diagnostics)
- `robustness_partial_correlations.csv` (optional partial-correlation controls; empty header-only file if covariates unavailable)

**What it does:**
1. Runs multivariable Cox proportional hazards models for Active Masking vs Decoherence survival
2. Controls for tumor stage (categorical: I/II/III/IV/missing)
3. Optionally controls for purity and microenvironment covariates if provided
4. Tests proportional hazards assumption via covariate × log(time) interaction screens
5. Tests continuous PD-L1 × circadian clock interaction models
6. Runs a sensitivity model adding continuous PD-L1 + B2M to test boundary mode beyond class components
7. Applies Benjamini-Hochberg FDR correction across primary robustness tests

**Interpretation note:** PH interaction screens in current outputs indicate broad non-proportionality across several covariates. Treat Cox hazard ratios as time-averaged associations unless followed by explicit time-varying or stratified survival models.

**Optional Covariate File (Purity/Microenvironment):**

This script can optionally incorporate tumor purity and microenvironment estimates as covariates. To use this feature, create a CSV file with the following structure:

**File name** (one of these, checked in order):
- `tcga_purity_immune_covariates.csv` (recommended)
- `tcga_purity_covariates.csv`
- `tcga_purity.csv`

**Required columns:**
- `case_barcode` or `bcr_patient_barcode` or `submitter_id` (patient identifier)

**Optional columns** (any recognized variant):
- **Tumor purity**: `purity`, `tumor_purity`, `cpe`, `estimate_purity`, `absolute_purity`, `consensus_purity`
- **Lymphoid infiltration**: `lymphocyte_infiltration`, `lymphocyte_score`, `immune_lymphoid`
- **Myeloid infiltration**: `myeloid_infiltration`, `myeloid_score`, `immune_myeloid`, `macrophage_regulation`
- **Leukocyte fraction**: `leukocyte_fraction`, `leukocyte_score`, `immune_leukocyte`
- **Stromal fraction**: `stromal_fraction`, `stromal_score`
- **IFN-γ response**: `ifn_gamma_response`, `ifng_response`, `ifn-gamma_response`
- **TGF-β response**: `tgfb_response`, `tgf_beta_response`, `tgf-beta_response`
- **Project ID**: `project_short_name`, `project_id`, `cohort`

**Where to get these estimates:**
- **ESTIMATE algorithm**: Yoshihara et al. 2013 - estimates tumor purity and immune/stromal scores from expression data
- **ABSOLUTE**: Carter et al. 2012 - absolute tumor purity from copy number and allelic fractions
- **xCell/CIBERSORT**: Cell-type deconvolution tools for immune infiltration estimates
- **Pre-computed TCGA estimates**: Available from Thorsson et al. 2018 Immunity paper (PanCancer immune landscape)

**Note**: If no covariate file is found, the script runs without optional controls (stage-only adjustment). This is **acceptable** because:
1. Our primary results (circadian CV vs PD-L1 correlation, Active Masking subtype definition) do not depend on survival modeling
2. Stage-adjusted Cox models provide an initial clinical-covariate check
3. Optional purity/microenvironment controls are a sensitivity check, not a primary requirement

**When to run this script:**
- For manuscript peer review requiring multivariable survival models
- To test robustness of survival findings to stage and microenvironment confounders
- To explore continuous PD-L1 × clock interactions beyond categorical Active Masking/Decoherence subtypes

**When NOT to run this script:**
- For initial exploratory analysis (use `tcga_survival.py` instead - faster, simpler)
- If you don't have covariate data and reviewers haven't requested it
- For paper figures (main figures come from `tcga_survival.py`)

---

---

#### 9. **sensitivity_threshold_sweep.py**
**Type**: Subtype classification robustness analysis
**Runtime**: ~30 seconds
**Inputs**: `tcga_expanded_tpm.csv`, `tcga_clinical.csv`
**Outputs**:
- `threshold_sensitivity_results.csv` (log-rank p-values across 11 percentile thresholds)
- `results/threshold_sensitivity_sweep.png` (bar chart: -log10(p) vs threshold percentile)

**What it does:**
1. Varies the percentile cutoff used for AM/DC classification from 25th to 75th in steps of 5
2. Re-classifies all tumors at each threshold
3. Re-runs log-rank survival test (AM vs DC) for SKCM and LUAD at each threshold
4. Reports whether the melanoma survival signal is robust to threshold choice

**Key finding**: SKCM survival signal is significant from the 25th to 65th percentile (9/11 thresholds, p < 0.05). The median (50th) is actually conservative — the strongest signal occurs at the 35th-40th percentile (p = 0.0005). This demonstrates that the survival finding is not an artifact of the specific threshold chosen.

---

#### 10. **sensitivity_rmst.py**
**Type**: PH-assumption-free survival analysis
**Runtime**: ~2 minutes (includes 1,000 bootstrap iterations)
**Inputs**: `tcga_expanded_tpm.csv`, `tcga_clinical.csv`
**Outputs**:
- `rmst_results.csv` (RMST estimates and bootstrap 95% CIs at 4 time horizons)
- `results/rmst_am_vs_dc.png` (bar chart: RMST difference with CI per cohort)

**What it does:**
1. Computes Restricted Mean Survival Time (RMST) — area under the KM curve up to a time horizon tau — for AM and DC groups separately
2. Tests RMST difference (AM minus DC) at tau = 36, 60, 84, 120 months
3. Constructs bootstrap 95% confidence intervals (n=1,000 resamples, seed=42)
4. Reports results for SKCM and LUAD

**Why it matters**: Cox models in the robustness analysis show PH violation fractions of 0.64-1.0. RMST requires NO proportional hazards assumption. It answers: "On average, how many more months do AM patients survive than DC patients over this time window?"

**Key finding**: Melanoma AM patients survive 24.2 months longer than DC patients at the 10-year horizon (bootstrap 95% CI: 10.8-36.8 months). Significant at every time horizon tested (36, 60, 84, 120 months). LUAD shows a 7.8-month advantage at 5 years (CI: 1.3-14.2).

---

#### 11. **sensitivity_immune_residualization.py**
**Type**: Stromal contamination sensitivity analysis
**Runtime**: ~30 seconds
**Inputs**: `tcga_expanded_tpm.csv`, `tcga_purity_immune_covariates.csv`
**Outputs**:
- `immune_residualization_results.csv` (R2 of immune effects on CV, residual correlations, purity-tertile results)
- `results/immune_residualization_scatter.png` (residual CV vs PD-L1 scatter after immune removal)
- `results/purity_stratified_correlations.png` (CV vs PD-L1 within purity tertiles)

**What it does:**
1. **Continuous residualization**: Regresses circadian CV on immune covariates (leukocyte fraction, lymphocyte infiltration, IFN-gamma response, purity, stromal fraction) via OLS, then correlates the residual CV with PD-L1
2. **Purity stratification**: Splits tumors into purity tertiles and re-computes the CV-PD-L1 correlation within each stratum

**Why it matters**: The primary critique of the "locked clock" finding is that it could reflect immune cell infiltration (T-cells expressing clock genes) rather than tumor-intrinsic restructuring. This analysis quantifies how much of the circadian CV variance is explained by immune infiltration, and whether the CV-PD-L1 coupling persists after removing those effects.

**Key finding**: Immune covariates explain only 0.7-13.8% of circadian CV variance (R2). After removing immune effects, the CV-PD-L1 correlation remains significant in all four cohorts tested (SKCM, LUAD, BRCA, HNSC). The correlation also persists within all three purity tertiles for SKCM, LUAD, and BRCA.

---

#### 12. **composite_observability_index.py**
**Type**: Independent metric validation (breaks tautology concern)
**Runtime**: ~30 seconds
**Inputs**: `tcga_expanded_tpm.csv`, `tcga_clinical.csv`
**Outputs**:
- `observability_index_results.csv` (AM/DC separation, survival, PD-L1 correlation)
- `results/observability_index_survival.png` (KM curves: high vs low observability)

**What it does:**
1. Constructs a composite "observability index" from features **independent of the AM/DC classification rule**:
   - MHC-I: mean log2(HLA-A, HLA-B, HLA-C) — NOT B2M (used in DC definition)
   - Gap junctions: mean log2(GJA1, GJB2, GJA5, GJB6) — not used in AM/DC at all
   - Orthogonal clock: CV of CLOCK, PER2, CRY1, CRY2 — NOT ARNTL or PER1 (used in AM definition)
2. Z-scores each component within cohort and averages with equal weights
3. Tests whether this index: (a) separates AM from DC, (b) predicts survival, (c) correlates with PD-L1

**Why it matters**: If cancer = loss of distinction and therapy = restoration of distinction, the framework risks circularity unless distinction is independently operationalized. This index uses orthogonal features to show the boundary-failure concept captures real biology, not just definitional artifacts.

**Key finding**: The index separates AM from DC with p = 1.7x10^-10 (SKCM) and p = 1.5x10^-20 (LUAD), and correlates with PD-L1 at rho = +0.48 (SKCM) and +0.54 (LUAD). Directional survival trend (p = 0.063 SKCM, p = 0.083 LUAD) — not FDR-significant but consistent, expected given weaker/orthogonal features.

---

#### 13. **run_external_validation.py**
**Type**: External cohort replication (non-TCGA)
**Runtime**: ~10 seconds (after data download)
**Inputs**: Pre-downloaded GEO data in `geo_cache/` (see below)
**Outputs**:
- `external_validation_results.csv` (replication statistics)
- `results/external_validation_geo.png` (3-panel: CV vs PD-L1, response by subtype, CV by response)

**Data acquisition (prerequisite):**
The script `external_validation_geo.py` downloads the raw data from GEO:
```bash
python external_validation_geo.py  # Downloads GSE91061 FPKM matrix + metadata
python run_external_validation.py   # Runs the validation analysis
```
This creates `geo_cache/` containing:
- `GSE91061_fpkm.csv.gz` (16 MB, FPKM expression matrix, 22,187 genes x 109 samples)
- `GSE91061_family.xml` (sample metadata with response annotations)
- `gse91061_target_expression.csv` (extracted 15 target genes x 109 samples)
- `gse91061_sample_meta.csv` (patient/timepoint parsed from column names)

**Note**: `geo_cache/` is `.gitignore`d. Raw GEO downloads must be regenerated by running `external_validation_geo.py`.

**What it does:**
1. Loads pre-treatment melanoma samples from GSE91061 (Riaz et al. 2017, nivolumab)
2. Maps Entrez gene IDs to our target gene symbols (15/15 found)
3. Computes circadian CV and tests correlation with PD-L1
4. Classifies AM/DC using within-cohort medians
5. Tests immunotherapy response rate by boundary-failure subtype (Fisher exact test)

**Key finding**: The circadian CV-PD-L1 negative correlation replicates in this independent cohort (rho = -0.283, p = 0.049). This is the 7th cohort showing the same direction (7/7 negative). The AM-vs-DC response rate comparison (Fisher p = 1.0) is underpowered (n=10 AM) and non-significant, as expected for this sample size.

---

### Utility Scripts (Discovery and Schema Inspection)

These scripts were used during development to discover the BigQuery schema. **Not required for reproduction**, but included for transparency:

- **tcga_discover.py**: Lists all TCGA tables and inspects RNAseq schema
- **tcga_clinical_discover.py**: Discovers clinical_gdc_current schema
- **tcga_clinical_discover2.py**: Discovers diagnosis and follow-up table names
- **tcga_clinical_discover3.py**: Inspects diagnosis and follow-up schemas with sample rows
- **tcga_schema.py**: General schema exploration utility

You can ignore these unless you want to understand how the BigQuery schema was discovered.

---

## File Descriptions

### Data Files (CSV)

| File | Rows | Columns | Description |
|------|------|---------|-------------|
| `tcga_expanded_tpm.csv` | 3,920 | 25 | Expression data (21 genes) for 6 cancer types |
| `tcga_clinical.csv` | 3,646 | 12 | Clinical data with survival outcomes |
| `tcga_expression_tpm.csv` | 1,062 | 25 | *Legacy*: Original 2-cancer extraction (SKCM+LUAD only) |
| `tcga_expression_counts.csv` | 1,062 | 25 | *Legacy*: Original 2-cancer extraction (count data) |
| `hypothesis1_correlations.csv` | 34 | 6 | *Legacy 2-cohort output* from `tcga_analysis.py` (SKCM+LUAD) |
| `tcga_multicancer_correlations.csv` | 102 | 8 | Detailed multi-cancer correlation results |
| `immune_subtype_comparison.csv` | 6 | 18 | Subtype CV contrasts (Kruskal + pairwise Mann-Whitney with FDR) |
| `tumor_normal_comparison.csv` | 35 | varies | Tumor-normal comparisons (paired Wilcoxon + unpaired Mann-Whitney with FDR) |
| `stage_analysis_results.csv` | 48 | 6 | Stage-trend and stage-group tests across markers/cancers |
| `survival_logrank_results.csv` | 12 | 10 | Log-rank tests with within-family and global FDR columns |
| `master_fdr_results.csv` | 244 | 9 | Comprehensive FDR-corrected p-values across all analyses |
| `tcga_purity_immune_covariates.csv` | varies | 3-10 | **OPTIONAL**: Tumor purity and microenvironment estimates for robustness checks (user-provided; see script #8 documentation for accepted columns) |
| `robustness_cox_coefficients.csv` | varies | 13 | Multivariable Cox coefficients (HR, CI, p-value per covariate) |
| `robustness_primary_tests.csv` | varies | varies | Primary robustness tests with BH correction and PH-screen caution columns |
| `robustness_ph_diagnostics.csv` | varies | varies | Proportional hazards interaction diagnostics (generated by `robustness_check.py`) |
| `robustness_partial_correlations.csv` | varies | 7 | Optional partial-correlation sensitivity tests with available controls |
| `threshold_sensitivity_results.csv` | 22 | 6 | Log-rank p-values for AM-vs-DC survival across 11 percentile thresholds (SKCM + LUAD) |
| `rmst_results.csv` | 8 | 10 | RMST estimates and bootstrap 95% CIs at 4 time horizons (SKCM + LUAD) |
| `immune_residualization_results.csv` | ~16 | 8 | Immune R2, residual correlations, and purity-tertile results for 4 cohorts |
| `observability_index_results.csv` | 6 | 6 | Composite observability index: AM/DC separation, survival, PD-L1 correlation |
| `external_validation_results.csv` | 1 | 11 | GSE91061 external validation: CV-PD-L1 replication and response analysis |
| `control_budget_test_results.csv` | varies | varies | MYC vs coherence-proxy data for budget-escape empirical test (generated by `test_control_budget.py`) |
| `evidence_ledger_prespecified.csv` | varies | varies | Source CSV for the pre-specified evidence table (`tables/evidence_table_prespecified.tex`) |

### Figures (PNG, 300 DPI)

| File | Type | Description |
|------|------|-------------|
| `hypothesis1_TCGA_SKCM.png` | Scatter | *Legacy*: SKCM-only circadian CV vs PD-L1 |
| `hypothesis1_TCGA_LUAD.png` | Scatter | *Legacy*: LUAD-only circadian CV vs PD-L1 |
| `multicancer_correlation_heatmap.png` | Heatmap | Effect sizes (rho) for 3 hypotheses × 6 cancer types |
| `multicancer_circadian_cv_forest.png` | Forest plot | Circadian CV vs PD-L1 correlation across cancer types |
| `tumor_normal_circadian_cv_paired.png` | Paired plot | Paired tumor-normal circadian CV comparison |
| `tumor_normal_circadian_cv_summary.png` | Boxplot | Summary of tumor vs normal across all cancer types |
| `immune_subtype_circadian_cv.png` | Violin | Circadian CV distribution by boundary-failure mode |
| `immune_subtype_boundary_scatter.png` | Scatter | PD-L1 vs circadian CV, colored by Active Masking/Decoherence |
| `survival_circadian_quartile.png` | KM curves | Survival stratified by circadian CV quartiles |
| `survival_boundary_failure.png` | KM curves | Survival stratified by Active Masking vs Decoherence |
| `stage_circadian_cv_boxplot.png` | Boxplot | Circadian CV by tumor stage (I/II/III/IV) |
| `control_budget_pareto_frontier.png` | Scatter | Pareto frontier overlay for budget-escape test (generated alongside `control_budget_combined.png` by `test_control_budget.py`) |

**Sensitivity analysis figures** (generated by scripts 9-13, written to `results/`):

| File | Type | Description |
|------|------|-------------|
| `threshold_sensitivity_sweep.png` | Bar chart | -log10(p) across percentile thresholds for AM/DC survival |
| `rmst_am_vs_dc.png` | Bar chart | RMST difference (AM-DC months) with bootstrap 95% CI |
| `immune_residualization_scatter.png` | Scatter | Residual CV vs PD-L1 after immune-fraction removal |
| `purity_stratified_correlations.png` | Scatter | CV vs PD-L1 within tumor-purity tertiles |
| `observability_index_survival.png` | KM curves | Survival by independent observability index |
| `external_validation_geo.png` | 3-panel | External validation: GSE91061 nivolumab melanoma |

## Key Statistical Methods

1. **Spearman rank correlation** (`scipy.stats.spearmanr`): Non-parametric correlation for non-normal distributions
2. **Benjamini-Hochberg FDR correction** (`scipy.stats.false_discovery_control`): Controls false discovery rate across 244 multiple tests
3. **Wilcoxon signed-rank test** (`scipy.stats.wilcoxon`): Paired comparison for matched tumor-normal samples
4. **Mann-Whitney U test** (`scipy.stats.mannwhitneyu`): Unpaired comparison for independent groups
5. **Kruskal-Wallis H test** (`scipy.stats.kruskal`): Multi-group comparison (Active Masking vs Decoherence vs Mixed)
6. **Kaplan-Meier estimation**: `scipy.stats.ecdf` with `CensoredData.right_censored()` for survival curves
7. **Log-rank test** (`scipy.stats.logrank`): Compare survival curves between groups

## Reproducibility Notes

### Why BigQuery Instead of Raw Files?

**Problem**: TCGA raw RNA-seq files are massive:
- ~60,000 genes per sample
- ~10,000 samples per cancer type
- Total download: **50-100 GB** per cancer type

**Solution**: ISB-CGC BigQuery public datasets allow **server-side filtering**:
- Query only 21 target genes (not all 60,000)
- Filter to 6 cancer types
- Pivot to wide format on the server
- Download only **~2 MB** of final data

**BigQuery cost**: ~$0.10 for the entire pipeline (processes ~50 GB of data, but you only download a few MB).

### Reproducing Exact Results

All random seeds are fixed where applicable. However, note:
- BigQuery table `*_current` datasets are **periodically updated** by ISB-CGC
- To ensure exact reproducibility, you can use versioned tables (e.g., `*_r35` for release 35) instead of `*_current`
- Current results based on data accessed: **February 2026**

### Active Masking Interpretation

The key conceptual shift in this analysis:

**Previous interpretation (incorrect)**:
- "Tumors have broken/destroyed circadian clocks"
- Lower CV = less coherence

**Correct interpretation (this analysis)**:
- **Active Masking tumors have LOWER circadian CV** (more coherent!)
- The clock is **locked in growth phase** (high BMAL1, low PER1), not broken
- This "locked" state correlates with high PD-L1 → active immune evasion
- These patients **survive longer** because the clock still exists and can potentially be reset

**Decoherence tumors**:
- Higher circadian CV (genuinely incoherent)
- Low PD-L1 + low MHC-I → invisible to immune system
- Worse survival outcomes

This interpretation is consistent with:
1. Tumors having **lower** CV than normals (locked ≠ broken)
2. Active Masking surviving **longer** than Decoherence (locked clock can be reset)
3. Combination therapy prediction (anti-PD-1 + differentiation agent works best for Active Masking)

## Citation

If you use this pipeline or data, please cite:

> Rahnama, M. (2026). *Cancer as Boundary Logic Failure: A Multi-Cancer Validation of Temporal Coherence Loss and Immune Evasion*. Manuscript in preparation.

## License

Repository code is covered by the repository's custom non-commercial research license (`LICENSE` at repo root). TCGA data access/usage remains subject to NCI/ISB-CGC terms.

## Contact

For questions about the pipeline or data:
- GitHub: https://github.com/MosesRahnama/cancer-paper-repository
- Issues: https://github.com/MosesRahnama/cancer-paper-repository/issues
