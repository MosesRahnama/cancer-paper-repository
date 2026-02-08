# Cancer as Boundary Logic Failure: Simulation Framework

**Author:** Moses Rahnama (Mina Analytics)
**Paper:** *Cancer as Boundary Logic Failure: A Computational and Information-Theoretic Framework for the Internal Self-Referential Malignancy*

---

## Repository Overview

A simulation framework implementing the boundary logic failure model of cancer. The model represents tissue as a labeled graph where cells are nodes connected by regulatory, temporal, and adhesion edges. Cancer emerges when a cell loses informational coupling with the organism and enters a self-referential duplication loop analogous to `rec_succ` in computability theory.

The framework allows you to:
- Simulate cancer emergence from healthy tissue
- Track tissue coherence (kappa) as cancer progresses
- Apply four distinct therapeutic interventions
- Visualize coherence collapse, cancer growth, and energy depletion
- Compare therapy outcomes

**This is a theoretical model, not a clinical tool.** It demonstrates logical failure modes, not molecular causation.
Legacy extraction artifacts are intentionally excluded from this published repository.

---

## Quick Start

```bash
# from repo root
python -m venv .venv
# Windows PowerShell
.venv\Scripts\Activate.ps1

# editable install with test tools
pip install -e ".[dev]"

# Run manuscript-linked scenarios (writes to results/)
python examples/run_full_scenario.py
python run_simulation.py
```

---

## Repository Structure

```
cancer-paper-repository/
  README.md                          # This file
  pyproject.toml                     # Package metadata/build config
  requirements.txt                   # Python dependencies
  requirements-dev.txt               # Dev/test dependency set
  CONTRIBUTING.md                    # Contribution protocol
  RELEASE_CHECKLIST.md               # Pre-release verification steps
  CITATION.cff                       # Citation metadata
  LICENSE                            # Rights and usage terms
  .gitignore                         # Git ignore rules
  run_simulation.py                  # Simple demo script
  run_phase2_extractions.py          # Phase-2 extracted operators + multi-seed sweep

  src/
    core/
      cell.py                        # Cell class with state, identity alignment, energy cost
      tissue.py                      # Tissue graph: nodes, edges, mutation, progression, division
      metrics.py                     # Coherence (kappa), entropy, signaling density, alignment

    simulations/
      engine.py                      # Simulation wrapper with history recording
      therapy.py                     # 4 therapeutic operators + apply_therapy_to_all_cancer
      parameter_profiles.py          # Evidence-aware run parameter presets

    visualization/
      plots.py                       # Coherence-vs-cancer plot, tissue network, state distribution

  examples/
    run_full_scenario.py             # Compares 5 therapy strategies, saves plots + JSON

  tests/
    test_parameter_profiles.py       # Parameter conversion and validation
    test_core_dynamics.py            # Core transition + pre-cancer accumulation checks
    test_therapy.py                  # Therapy operator behavior checks

  results/
    README.md                        # Artifact index and interpretation scope
    scenario_comparison.json         # Numerical results from therapy comparison
    scenario_comparison_metadata.json # Parameter profile + source metadata
    therapy_comparison.png           # Multi-therapy coherence + cancer + pre-cancer plot
    kappa_trajectory.png             # Single-therapy coherence trajectory
```

---

## Core Concepts

### Cell States
| State | Meaning | Division Rate | Energy Cost |
|-------|---------|--------------|-------------|
| HEALTHY | Organism-aligned, boundary-coupled | 0.01 | 1.0 |
| PRE_CANCER | Boundary weakening, partial decoupling | 0.01 | 1.0 |
| CANCER | Self-referential, fully decoupled (rec_succ) | 0.20 | 2.7 (illustrative Warburg proxy) |

### Edge Types
| Type | Label | Meaning |
|------|-------|---------|
| REGULATORY | SIG | Growth inhibition, contact inhibition, differentiation cues |
| LINEAGE | DIV | Parent-child division history |
| TEMPORAL | TMP | Circadian/hormonal timing sync |
| ADHESION | HAR | Gap junction communication, tissue adhesion |

### Coherence Metric (kappa)

```
kappa = w_sig * rho_sig + w_tmp * rho_tmp + w_identity * (1 - H)
```

Where:
- `rho_sig` = regulatory signaling density
- `rho_tmp` = temporal coordination density
- `H` = mean identity entropy across cells
- Weights: (0.4, 0.2, 0.4)

| kappa | Interpretation |
|-------|---------------|
| 0.9-1.0 | Healthy tissue |
| 0.7-0.9 | Minor perturbation |
| 0.5-0.7 | Pre-cancerous |
| 0.3-0.5 | Early malignancy |
| 0.1-0.3 | Advanced malignancy |
| 0.0-0.1 | Terminal / Metastatic |

### Cancer Transform (what happens when a cell becomes cancerous)

1. Identity alignment drops: p_align -> 0.05
2. Division rate increases: 0.01 -> 0.20
3. Energy cost increases: 1.0 -> 2.7 (illustrative Warburg burden proxy)
4. Regulatory and temporal edges are severed (boundary collapse)
5. Cancer daughters receive NO regulatory or temporal edges (rec_succ without halt)

---

## Therapeutic Operators

| Therapy | What It Does | Biological Analog |
|---------|-------------|-------------------|
| `forced_distinction` | Partially restores organism alignment + temporal coupling | Neoantigen presentation, making cancer "visible" |
| `differentiation_therapy` | Fully restores healthy state, halts division, reconnects edges | ATRA in acute promyelocytic leukemia (>90% CR) |
| `checkpoint_inhibitor` | Partially restores organism alignment (immune visibility) | Anti-PD-1/PD-L1 antibodies |
| `bioelectric_reprogramming` | Restores alignment + adhesion, reduces division rate | Ion channel manipulation (Levin's work) |

---

## Latest Results

From `examples/run_full_scenario.py` (100 cells, seed=42, calibrated monthly profile):

| Therapy | Pre-therapy kappa | Post-therapy kappa | Cancer cells remaining |
|---------|------------------|--------------------|----------------------|
| No Therapy (Control) | 0.418 | 0.418 | 376 |
| Forced Distinction | 0.418 | 0.231 | 376 |
| **Differentiation Therapy (ATRA)** | **0.418** | **0.441** | **0** |
| Checkpoint Inhibitor | 0.418 | 0.232 | 376 |
| Bioelectric Reprogramming | 0.418 | 0.199 | 376 |

This run uses `calibrated-monthly-oed` parameters (`p_mutation=0.002`, `p_progression=0.0092`, `p_healthy_division=0.01`) and logs profile metadata in `results/scenario_comparison_metadata.json`.

In this illustrative parameterization (single seed, fixed operator settings), differentiation therapy is the only intervention that eliminates all cancer cells. This should not be read as a universal efficacy ranking; it shows that stronger identity restoration in the toy model can reach a halt state, while partial distinction-restoration operators improve observability without full clearance.
The cancer-state `energy_cost=2.7` setting is an illustrative modeling choice for qualitative stress testing, not an empirically fitted physiological constant.

---

## Phase 2 Extraction Runner

Run the organized high-value extraction script:

```bash
python run_phase2_extractions.py
```

Optional stressor mode with unmasking:

```bash
python run_phase2_extractions.py --seeds 40 --apply-unmasking --unmask-severity 0.5
```

Outputs are written to `results/` (JSON + summary plots).

---

## Reproducibility and Quality Gates

```bash
# unit tests
pytest

# script compile smoke checks
python -m py_compile run_simulation.py
python -m py_compile run_phase2_extractions.py
python -m py_compile examples/run_full_scenario.py
```

CI workflow is included at `.github/workflows/ci.yml`.

---

## Foundational References

1. Rahnama, M. (2026). *Cancer as Boundary Logic Failure*. [This paper]
2. Rahnama, M. (2025). *Strong Normalization for the Safe Fragment of a Self-Referential Operator Kernel*. arXiv:2512.00081
3. Rahnama, M. (2025). *The Observation Principle, Laws of Entanglement, Laws of Information*. Zenodo. https://zenodo.org/records/17861545
4. Rahnama, M. (2026). *Thermodynamic Constraints on Measurement Events: A Boundary Framework for Classical Information*. Zenodo. https://doi.org/10.5281/zenodo.18445561

---

## Limitations

- This is NOT a predictive cancer model
- Does NOT capture genetics, metabolism, angiogenesis, immune editing, clonal heterogeneity, or spatial constraints
- The coherence metric kappa is not calibrated against clinical data
- Simulation parameters are illustrative, not empirically derived
- Demonstrates logical failure modes, not molecular causation
- Several parameters (including cancer `energy_cost=2.7`) are illustrative placeholders pending calibration/sensitivity studies

---

## License

All rights reserved by the author. See `LICENSE`.
