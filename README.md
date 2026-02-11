# Cancer as Boundary Logic Failure

**Author:** Moses Rahnama (Mina Analytics)
**Paper:** *Cancer as Boundary Logic Failure: A Computational and Information-Theoretic Framework for the Internal Self-Referential Malignancy*
**Manuscript:** [`paper/Cancer_As_Boundary_Logic_Failure.tex`](paper/Cancer_As_Boundary_Logic_Failure.tex)

---

## What This Repository Contains

A complete reproducible research package for the boundary logic failure theory of cancer, comprising:

1. **The manuscript** (`paper/`) with all LaTeX source, tables, and figure references
2. **TCGA empirical analysis pipeline** (`experiments/tcga/`) — extraction, correlation, survival, robustness, and sensitivity analyses across 6 cancer types (n=3,611 tumors)
3. **Sensitivity and validation analyses** — threshold robustness, RMST, immune residualization, composite observability index, and external validation in an independent immunotherapy cohort
4. **Toy simulation framework** (`src/`, `examples/`, `run_simulation.py`) — tissue graph model with four therapeutic operators
5. **All generated outputs** (`results/`) — figures, CSVs, and JSON used in the manuscript

**This is a theoretical framework with empirical support, not a clinical tool.**

---

## Key Empirical Findings

| Finding | Evidence | Location |
|---------|----------|----------|
| Circadian CV–PD-L1 negative coupling | 6/6 TCGA cohorts, all FDR q < 0.05 (ρ: −0.125 to −0.381) | Section 10.1 |
| External replication | 7th cohort (GSE91061 nivolumab melanoma): ρ = −0.283, p = 0.049 | Section 10.9 |
| "Locked, not broken" clock | Tumors show lower circadian CV than matched normals in 4/5 cancer types | Section 10.4 |
| Active Masking vs Decoherence survival | SKCM: log-rank p = 0.0011, BH q = 0.0132; RMST +24 months at 10 years | Section 10.6 |
| Threshold robustness | SKCM signal significant at 9/11 percentile thresholds (25th–65th) | Section 10.6 |
| Stromal confound survived | Residual CV–PD-L1 significant in 4/4 cohorts after immune residualization | Section 10.8 |
| Independent observability index | Separates AM from DC at p < 10⁻¹⁰ using orthogonal features | Section 10.9 |
| Budget Escape (not trade-off) | Positive proliferation–coherence correlation in 5/6 cancer types | Section 6 |

---

## Repository Structure

```
cancer-paper-repository/
  paper/
    Cancer_As_Boundary_Logic_Failure.tex   # Main manuscript (LaTeX)
    tables/
      evidence_table_prespecified.tex      # Pre-specified endpoint audit table

  experiments/tcga/
    README.md                              # Detailed pipeline documentation (15 scripts)
    tcga_config.py                         # Shared configuration and utility functions

    # --- Data Extraction (requires GCP BigQuery) ---
    tcga_extract.py                        # Original 2-cohort extraction
    tcga_extract_expanded.py               # 6-cohort extraction → tcga_expanded_tpm.csv
    tcga_clinical_discover*.py             # Schema discovery utilities

    # --- Core Analyses ---
    tcga_multicancer.py                    # Cross-cancer correlation validation
    tcga_analysis.py                       # Legacy 2-cohort analysis
    tcga_tumor_normal.py                   # Tumor vs matched normal comparison
    tcga_immune_subtype.py                 # AM/DC boundary-failure classification
    tcga_survival.py                       # KM survival + log-rank tests
    tcga_stage_analysis.py                 # Stage-stratified analysis + global FDR
    robustness_check.py                    # Multivariable Cox models with covariates

    # --- Sensitivity Analyses (new) ---
    sensitivity_threshold_sweep.py         # AM/DC threshold robustness (25th–75th pct)
    sensitivity_rmst.py                    # PH-free RMST survival analysis
    sensitivity_immune_residualization.py   # Stromal confound: residualization + purity strata
    composite_observability_index.py       # Independent observability metric
    run_external_validation.py             # GSE91061 nivolumab melanoma validation
    external_validation_geo.py             # GEO data download utility

    # --- Manuscript Figure Generation ---
    generate_manuscript_figures.py         # Summary figures from CSV outputs
    sync_manuscript_artifacts.py           # Deterministic figure sync + provenance
    convert_thorsson_to_covariates.py      # Immune covariate preprocessing
    test_control_budget.py                 # Control-budget empirical stress test

    # --- Data Files (CSV) ---
    tcga_expanded_tpm.csv                  # 3,920 samples × 25 columns
    tcga_clinical.csv                      # 3,646 patients × 12 columns
    tcga_purity_immune_covariates.csv      # Thorsson immune landscape covariates
    tcga_multicancer_correlations.csv      # Cross-cancer Spearman results
    tumor_normal_comparison.csv            # Paired Wilcoxon + unpaired MW results
    immune_subtype_comparison.csv          # Subtype CV contrasts
    survival_logrank_results.csv           # Log-rank tests with FDR
    stage_analysis_results.csv             # Stage-trend tests
    master_fdr_results.csv                 # Comprehensive FDR across all analyses
    robustness_primary_tests.csv           # Cox robustness tests with FDR
    threshold_sensitivity_results.csv      # Threshold sweep results
    rmst_results.csv                       # RMST estimates + bootstrap CIs
    immune_residualization_results.csv     # Stromal confound sensitivity
    observability_index_results.csv        # Composite index results
    external_validation_results.csv        # GSE91061 replication statistics

  src/
    core/
      cell.py                              # Cell class: state, alignment, energy
      tissue.py                            # Tissue graph: nodes, edges, dynamics
      metrics.py                           # Coherence (κ), entropy, signaling density
    simulations/
      engine.py                            # Simulation wrapper with history
      therapy.py                           # 4 therapeutic operators
      parameter_profiles.py                # Calibrated parameter presets
    visualization/
      plots.py                             # Plotting utilities

  examples/
    run_full_scenario.py                   # Multi-therapy comparison scenario

  tests/
    test_core_dynamics.py                  # Core transition checks
    test_therapy.py                        # Therapy operator checks + growth arrest
    test_parameter_profiles.py             # Parameter conversion validation

  results/                                 # All manuscript-facing outputs
    README.md                              # Artifact index
    *.png                                  # 29 figures referenced by the manuscript
    *.json                                 # Simulation scenario outputs

  CONSISTENCY_REPORT.md                    # Internal consistency audit trail
  CITATION.cff                             # Citation metadata
  CONTRIBUTING.md                          # Contribution protocol
  LICENSE                                  # Non-commercial research license
```

---

## Reproduction

### Full Pipeline (from scratch)

```bash
# 1. Setup
python -m venv .venv && .venv\Scripts\Activate.ps1  # Windows
pip install -e ".[dev]"
pip install google-cloud-bigquery google-cloud-storage GEOparse

# 2. TCGA data extraction (requires GCP BigQuery credentials)
cd experiments/tcga
python tcga_extract_expanded.py

# 3. Core analyses
python tcga_multicancer.py
python tcga_tumor_normal.py
python tcga_immune_subtype.py
python tcga_survival.py
python tcga_stage_analysis.py
python robustness_check.py

# 4. Sensitivity analyses (uses existing CSV outputs, no GCP needed)
python sensitivity_threshold_sweep.py
python sensitivity_rmst.py
python sensitivity_immune_residualization.py
python composite_observability_index.py

# 5. External validation (downloads from GEO, no GCP needed)
python external_validation_geo.py
python run_external_validation.py

# 6. Manuscript figures
python generate_manuscript_figures.py

# 7. Simulation (from repo root)
cd ../..
python examples/run_full_scenario.py
python run_simulation.py
python run_phase2_extractions.py

# 8. Tests
pytest
```

### Without GCP (using committed CSV data)

Steps 4–8 above work without GCP credentials because the extracted CSV files are committed to the repository. Steps 2–3 require BigQuery access to ISB-CGC public datasets.

---

## Data Sources

| Source | Access | Used For |
|--------|--------|----------|
| ISB-CGC BigQuery (`isb-cgc-bq.TCGA.*`) | Public (GCP account required) | RNA-seq expression (21 genes × 6 cancer types) |
| TCGA clinical tables | Public (GCP account required) | Survival, staging, demographics |
| Thorsson et al. 2018 (*Immunity*) | Pre-computed, committed as CSV | Immune landscape covariates (purity, infiltration) |
| GEO GSE91061 (Riaz et al. 2017) | Public (internet required) | External validation (nivolumab melanoma) |

---

## Statistical Methods

1. **Spearman rank correlation** — non-parametric correlation across 6 cancer types
2. **Benjamini-Hochberg FDR** — multiple-testing correction across 244 tests (global) and within families
3. **Wilcoxon signed-rank** — paired tumor vs matched normal comparisons
4. **Mann-Whitney U** — unpaired group comparisons and purity-stratified tests
5. **Kruskal-Wallis H** — multi-group subtype comparisons
6. **Kaplan-Meier + log-rank** — survival analysis with `scipy.stats.ecdf` and `CensoredData`
7. **Multivariable Cox PH** — clinical covariate adjustment with PH diagnostics
8. **Restricted Mean Survival Time (RMST)** — PH-free survival comparison with bootstrap CIs
9. **OLS residualization** — immune-fraction confound sensitivity
10. **Fisher exact test** — response-rate comparison in external validation

---

## Simulation Framework

The toy simulation models tissue as a labeled graph. Cells are nodes; edges represent regulatory, temporal, and adhesion coupling. Cancer emerges when termination predicates fail and self-referential duplication (the `rec_succ` analogy) dominates.

| Cell State | Division Rate | Energy Cost | Organism Alignment |
|------------|--------------|-------------|-------------------|
| HEALTHY | 0.01 | 1.0 | 0.95 |
| PRE_CANCER | 0.01 | 1.0 | 0.50 |
| CANCER | 0.20 | 2.7 | 0.05 |

Four therapeutic operators implement the "forced distinction" design principle:
- **Differentiation therapy** (ATRA-like): full state restoration
- **Forced distinction**: partial observability increase
- **Checkpoint inhibitor** (anti-PD-1-like): immune visibility restoration
- **Bioelectric reprogramming** (Levin-like): adhesion and alignment restoration

**Disclaimer:** The simulation demonstrates logical failure modes, not molecular causation. Parameters are illustrative, not clinically calibrated.

---

## Consistency and Audit Trail

- [`CONSISTENCY_REPORT.md`](CONSISTENCY_REPORT.md): tracks all data–text–figure cross-checks, resolved discrepancies, and FDR scope decisions
- [`paper/figure_provenance.csv`](paper/figure_provenance.csv): SHA-256 hashes and generator scripts for every manuscript figure
- [`experiments/tcga/README.md`](experiments/tcga/README.md): detailed documentation for all 15 analysis scripts with inputs, outputs, and key findings

---

## Citation

```bibtex
@article{rahnama2026boundary,
  title={Cancer as Boundary Logic Failure: A Computational and
         Information-Theoretic Framework for the Internal
         Self-Referential Malignancy},
  author={Rahnama, Moses},
  year={2026},
  note={Manuscript in preparation.
        Repository: github.com/MosesRahnama/cancer-paper-repository}
}
```

See `CITATION.cff` for machine-readable citation metadata.

---

## License

Custom non-commercial research license. Non-commercial use (research, teaching, internal evaluation) is permitted with attribution. Commercial/clinical use requires a separate license. See [`LICENSE`](LICENSE) for full terms.
