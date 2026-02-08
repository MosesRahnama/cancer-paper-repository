# Release Checklist

1. Validate environment:
   - `pip install -e ".[dev]"`
2. Run quality gates:
   - `pytest`
   - `python -m py_compile run_simulation.py run_phase2_extractions.py examples/run_full_scenario.py`
3. Regenerate default manuscript artifacts:
   - `python examples/run_full_scenario.py`
   - `python run_simulation.py`
   - `python run_phase2_extractions.py`
4. Verify artifact set in `results/` matches `results/README.md`.
5. Confirm citation metadata in `CITATION.cff` (especially repository URL).
6. Tag version and publish.
