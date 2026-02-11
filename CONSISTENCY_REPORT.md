# Internal Consistency Check: Budget Escape Pivot

**Date:** 2026-02-11 (updated)
**Paper:** `paper/Cancer_As_Boundary_Logic_Failure.tex`
**Previous review:** 2026-02-10 (identified 5 "reduced investment" passages)
**Current review:** 2026-02-11 (full re-audit of manuscript, code, figures, data)

---

## Executive Summary

All five "reduced investment" passages flagged in the 2026-02-10 report have been **corrected**. The Budget Escape narrative is now internally consistent from Abstract through Conclusion. No `B_ctrl` remnants. No "reduced investment" language anywhere in the manuscript.

**Status:** ✅ **LANGUAGE CONSISTENCY RESOLVED**

Three discrepancy classes were tracked during the 2026-02-11 audit:

1. Stage figure/text FDR scope mismatch (Section 3) — resolved
2. Tumor–normal circadian-CV $p$ value drift between manuscript/figure/CSV (Section 4) — resolved
3. Survival BH-FDR scope inconsistency across text/table/plots (Section 5) — resolved

**Status:** ✅ **PRIMARY CONSISTENCY ISSUES RESOLVED**

---

## ✅ PREVIOUSLY FLAGGED ISSUES — ALL RESOLVED

### 1. Abstract (formerly Line 55)
- **Old:** "reduced differentiation investment"
- **Current:** "metabolic expansion (Warburg effect), proliferative burden, and tumor-immune inflammatory signaling produce systemic metabolic stress, including cachexia"
- **Status:** ✅ FIXED

### 2. Introduction / Distinction Section (formerly Line 193)
- **Old:** "investment in stable identity-bearing differentiation is reduced"
- **Current:** "the cost of stable identity-bearing differentiation is often externalized to the host"
- **Status:** ✅ FIXED

### 3. Introduction / Distinction Section (formerly Line 196)
- **Old:** "reduced local differentiation investment"
- **Current:** "metabolic expansion (Warburg), direct host resource burden, and tumor-immune inflammatory catabolism"
- **Status:** ✅ FIXED

### 4. Warburg Section (formerly Line 346)
- **Old:** "reduced investment in differentiation-specific regulatory order"
- **Current:** "a mechanism to support high-entropy proliferative states without sacrificing regulatory order...Malignancy expands the metabolic budget to support both proliferation and active evasion (Active Masking), increasing organism-level burden"
- **Status:** ✅ FIXED

### 5. Cachexia Section (formerly Line 354)
- **Old:** "reduced differentiation investment"
- **Current:** "local proliferative advantage with expanded metabolic capacity is coupled to global host burden and inflammatory wasting"
- **Status:** ✅ FIXED

### 6. Limitations Section (formerly Line 823, flagged as MINOR)
- **Old:** "reduced investment in identity-forming information"
- **Current:** "distinguishing energetic replication cost from decoupled identity-forming processes"
- **Status:** ✅ FIXED

---

## ✅ NARRATIVE FLOW VERIFICATION

The Budget Escape story now flows correctly through the manuscript:

1. **Section 6 (Control Budget):** Introduces homeostatic budget constraint → data contradicts simple trade-off → introduces Budget Escape ($C_G + C_S > B_{\text{homeostasis}}$) → supported by new Figure (control_budget_combined.png)
2. **Section 8 (Warburg):** Warburg as metabolic autonomy and expansion, not reduced investment
3. **Section 8 (Cachexia):** Three-mechanism model (substrate burden, budget escape, inflammatory catabolism) — no "reduced investment" language
4. **Section 10 (TCGA):** Budget Escape listed as Implication #3 with 5/6 positive-correlation evidence
5. **Conclusion:** "forced distinction" continuum, no trade-off language

---

## ✅ MANUSCRIPT STRUCTURAL INTEGRITY CHECKS (LaTeX)

Read-only integrity checks on `paper/Cancer_As_Boundary_Logic_Failure.tex` (including `\\input{tables/evidence_table_prespecified.tex}`):

- Figures: 23 `\\includegraphics{...}` calls; all referenced image files exist in-repo.
- Tables: 5 `table` environments; the pre-specified evidence table label `tab:evidence-table` is defined in the `tables/evidence_table_prespecified.tex` input file.
- Cross-references: 0 missing `\\ref{}` / `\\eqref{}` targets; 0 duplicate `\\label{}` keys.
- Citations: 44 unique `\\cite{}` keys; all have corresponding `\\bibitem{}` entries; 0 uncited bibitems.

---

## ✅ RESOLVED: STAGE FIGURE–TEXT FDR SCOPE ALIGNMENT

### Resolution Status

The manuscript (Section 10, Stage-Stratified Analysis, ~line 745) states:

> "11/12 circadian-stage tests were non-significant after global FDR correction, with one borderline signal (TCGA-HNSC Spearman trend, **q ≈ 0.050**)."

This q ≈ 0.050 value comes from the **global BH-FDR** across all ~200+ tests in `master_fdr_results.csv` (computed by `tcga_stage_analysis.py`). In that global pool, the HNSC circadian Spearman p=0.0287 gets q=0.04998 — essentially borderline FDR-significant.

Within-family BH correction **across only the 6 circadian Spearman tests** yields **q=0.172** for HNSC, while the **global** BH-FDR across the full pooled analysis set yields **q≈0.050** (as reported in `master_fdr_results.csv` and the pre-specified evidence table).

### Why This Matters

These are both valid corrections but applied to different test pools. The manuscript language explicitly says “global FDR correction”, so the stage Spearman summary plot should report the global q-values (and, if desired, optionally note the within-family q-values for transparency).

### Applied Fix

`generate_manuscript_figures.py` now includes a `_load_global_stage_qvalues()` helper that pulls the global q-values from `master_fdr_results.csv`. The current `results/stage_circadian_spearman_summary.png` shows HNSC at approximately `q=0.05`, matching manuscript wording and correction scope.

---

## ✅ RESOLVED: TUMOR–NORMAL P-VALUE / FIGURE SYNC

### Resolution Status

The manuscript Tumor-vs-Matched-Normal subsection (Section 10, ~lines 649–676), manuscript-facing figure (`results/tumor_normal_circadian_cv_paired.png`), and `experiments/tcga/tumor_normal_comparison.csv` are now aligned for circadian-CV paired Wilcoxon values.

For `Circadian_CV` paired Wilcoxon tests, the CSV currently contains:

| Cohort | Manuscript/Figure $p$ | Current CSV `wilcoxon_p` |
|---|---:|---:|
| TCGA-LUAD | $2.28\\times10^{-6}$ | $2.28\\times10^{-6}$ |
| TCGA-BRCA | $5.13\\times10^{-9}$ | $5.13\\times10^{-9}$ |
| TCGA-HNSC | $7.0\\times10^{-7}$ | $7.02\\times10^{-7}$ (matches) |
| TCGA-LUSC | $8.6\\times10^{-7}$ | $8.61\\times10^{-7}$ (matches) |
| TCGA-COAD | “not significant” | $0.564$ (still non-significant) |

### Applied Fix

- Copied latest tumor-normal figure artifacts from `experiments/tcga/figures/` into `results/`.
- Updated manuscript paired-Wilcoxon circadian-CV values for LUAD and BRCA to match the current CSV output.

---

## ✅ RESOLVED: SURVIVAL q-VALUE SCOPE CONSISTENCY

### Resolution Status

`experiments/tcga/survival_logrank_results.csv` contains two different BH-FDR corrections:
- `q_value_test_family`: BH within each test family (6 AM-vs-DC tests; 6 Q1-vs-Q4 tests)
- `q_value_all_survival_tests`: BH across all 12 survival tests (6 AM-vs-DC + 6 Q1-vs-Q4)

The manuscript narrative text, pre-specified evidence table, and survival summary plot now consistently use `q_value_all_survival_tests` (BH across all 12 survival tests).

### Concrete Example (AM vs DC)

| Cohort | $p$ | q (within AM-vs-DC family, m=6) | q (all survival tests, m=12) |
|---|---:|---:|---:|
| TCGA-SKCM | 0.0011 | 0.0066 | 0.0132 |
| TCGA-LUAD | 0.0489 | 0.147 | 0.294 |

### Applied Fix

- Updated Abstract, intro preview, survival results, and implications text to:
  - SKCM AM-vs-DC: $q=0.0132$
  - LUAD AM-vs-DC: $q=0.294$
- Updated survival captions to explicitly state BH correction across all 12 survival tests.
- Kept the evidence table aligned to the same q-scope.

---

## ✅ RESOLVED: ARTIFACT SYNC + CANONICAL MANUSCRIPT FIGURE DIRECTORY

The repository currently contains duplicate PNG basenames under both:
- `results/` (used by the manuscript via `\\includegraphics{results/...}`), and
- `experiments/tcga/figures/` (direct outputs of several TCGA scripts).

### Applied Fix

- Added deterministic sync script: `experiments/tcga/sync_manuscript_artifacts.py`
- Script behavior:
  - Parses manuscript `\\includegraphics{...}` entries
  - Syncs shared TCGA basenames between `results/` and `experiments/tcga/figures/`
  - Resolves conflicts by deterministic policy (newer mtime wins; tie-break to `results/`)
  - Optionally rewrites manuscript paths to canonical `results/...` form (`--rewrite-tex`)
- Executed sync + rewrite:
  - `python experiments/tcga/sync_manuscript_artifacts.py --rewrite-tex`
  - `control_budget_combined.png` is now referenced from `results/` (canonicalized)
  - Bytewise-different duplicate PNG basenames between `results/` and `experiments/tcga/figures/`: **0**

### Provenance Contract Implemented

The sync script now emits a provenance manifest:
- `paper/figure_provenance.csv`

For each manuscript figure, it records:
- `figure_path`
- `source_inputs`
- `generator_script`
- `figure_timestamp_utc`
- `figure_sha256`

---

## ✅ DATA–TEXT CROSS-CHECKS (all verified)

| Claim in Manuscript | Data Source | Match? |
|---|---|---|
| Circadian CV–PD-L1: ρ = −0.381 (SKCM) | `multicancer_correlations.csv` | ✅ |
| Circadian CV–PD-L1: ρ = −0.303 (LUAD) | `multicancer_correlations.csv` | ✅ |
| 100% replication across 6 cohorts (all FDR q<0.05) | `multicancer_correlations.csv` | ✅ |
| B2M–PD-L1: ρ = +0.760 (SKCM) | `hypothesis1_correlations.csv` | ✅ |
| AM-vs-DC survival SKCM: p=0.0011, q=0.0132 (all 12 tests) | `survival_logrank_results.csv` | ✅ |
| AM-vs-DC survival LUAD: p=0.049, q=0.294 (all 12 tests) | `survival_logrank_results.csv` | ✅ |
| Circadian CV quartiles: no FDR-significant survival | `survival_logrank_results.csv` | ✅ |
| Tumor < Normal circadian CV in 4/5 types | `tumor_normal_comparison.csv` | ✅ |
| Cox adjusted AM-vs-DC SKCM: HR=0.507, q=0.069 | `robustness_primary_tests.csv` | ✅ |
| Cox adjusted AM-vs-DC LUAD: HR=0.288, q=0.045 | `robustness_primary_tests.csv` | ✅ |
| PD-L1×clock interaction SKCM: HR=0.217, q=0.045 | `robustness_primary_tests.csv` | ✅ |
| Therapy table: ATRA → cancer=0, control → cancer=376 | `scenario_comparison.json` | ✅ |

---

## ✅ FIGURE–DATA CROSS-CHECKS (all verified)

| Figure | Data Source | Match? |
|---|---|---|
| `survival_test_qvalue_summary.png` | `survival_logrank_results.csv` | ✅ |
| `robustness_primary_hr_forest.png` | `robustness_primary_tests.csv` | ✅ |
| `immune_subtype_prevalence_stacked.png` | `immune_subtype_comparison.csv` | ✅ |
| `multicancer_circadian_cv_forest.png` | `multicancer_correlations.csv` | ✅ |
| `multicancer_correlation_heatmap.png` | `multicancer_correlations.csv` | ✅ |
| `survival_boundary_failure.png` | `survival_logrank_results.csv` | ✅ |
| `tumor_normal_circadian_cv_paired.png` | `tumor_normal_comparison.csv` | ✅ synced to current CSV (Section 4 resolved) |
| `stage_circadian_spearman_summary.png` | `master_fdr_results.csv` | ✅ global q-values (HNSC q≈0.050) |

---

## ✅ CODE QUALITY NOTES

- `tcga_extract.py` has been refactored to use `tcga_config.py` helpers (no more hardcoded credentials in the script itself)
- `tcga_tumor_normal.py` updated with proper sample deduplication (priority: Primary > Recurrent > Metastatic) and separated paired (Wilcoxon) from unpaired (Mann-Whitney) tests
- `test_therapy.py` has a new test (`test_differentiation_therapy_enforces_growth_arrest_even_with_high_global_division`) validating per-cell growth arrest
- `generate_manuscript_figures.py` is a clean figure-generation script reading from CSV outputs
- `sync_manuscript_artifacts.py` adds deterministic manuscript-figure sync and provenance generation (`paper/figure_provenance.csv`)
- Old root-level `.tex` copies removed; canonical manuscript now in `paper/` directory

---

## ⚠️ INTERPRETATION / REPRODUCIBILITY RISKS (NON-CLINICAL)

These are not “errors” but important constraints on interpretation and reproducibility:

- **Bulk RNA-seq CV is a proxy:** “Circadian CV” here is within-sample cross-gene dispersion in bulk RNA-seq, not a direct measure of rhythmic phase coherence; the manuscript correctly notes this limitation.
- **Survival is natural-history association:** TCGA lacks immunotherapy response endpoints and treatment is heterogeneous; subtype–survival associations are not evidence of checkpoint blockade efficacy.
- **Definitional coupling in subtype analyses:** The Active Masking definition includes ARNTL/PER1, so downstream CV differences are internal-consistency evidence, not independent validation.
- **Multiple-testing scope must be explicit:** Several outputs provide both within-family and global q-values; mixing them without explicit scope creates avoidable inconsistencies (see Sections 3 and 5).
- **Cox PH sensitivity:** The robustness pipeline includes PH diagnostics; if broad non-proportionality is present, Cox coefficients should be treated as time-averaged associations.

---

## SUMMARY

| Category | Status |
|---|---|
| "Reduced investment" language | ✅ All 5+1 passages fixed |
| Budget Escape narrative flow | ✅ Consistent Abstract → Conclusion |
| Manuscript structural integrity (LaTeX) | ✅ Refs/cites/figure paths resolve |
| Data–text numerical accuracy | ✅ Synced to current CSV outputs and consistent q-scope reporting |
| Figure–data consistency | ✅ Stage and tumor-normal figures aligned with stated data sources |
| Code refactoring | ✅ Improved since last review |

**Action required:** none (blocking consistency issues resolved).

---

*Previous review (2026-02-10) identified 5 "reduced investment" passages; all were corrected before the 2026-02-11 audit. The stale Feb 10 report has been removed from this document to avoid confusion.*
