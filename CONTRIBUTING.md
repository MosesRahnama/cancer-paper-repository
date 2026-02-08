# Contributing

## Scope
This repository is a companion simulation framework for the manuscript. Keep contributions aligned with:
- Reproducibility of published figures/results
- Transparent assumptions and parameter provenance
- Clear separation of exploratory code vs manuscript-linked code

## Development Setup
```bash
python -m venv .venv
. .venv/Scripts/activate  # Windows PowerShell: .venv\Scripts\Activate.ps1
pip install -e ".[dev]"
```

## Required Checks
```bash
pytest
python run_simulation.py
python examples/run_full_scenario.py
```

## Change Guidelines
- Preserve default behavior of manuscript-facing scripts unless a change is explicitly justified.
- Document new parameters in `README.md` and include rationale/citation where applicable.
- For any new output artifact, update `results/README.md`.
- Avoid committing `__pycache__` and transient local artifacts.

## Pull Request Notes
Include:
1. What changed.
2. Why it changed.
3. Validation steps and exact commands run.

Maintainers: release procedure is documented in `.github/maintainers/RELEASE_CHECKLIST.md`.
