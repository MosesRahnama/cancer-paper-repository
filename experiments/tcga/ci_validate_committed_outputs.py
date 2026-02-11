"""
CI validator for committed TCGA artifacts and extraction version locking.

This script is intentionally deterministic:
- Uses only committed CSV/cache files (no BigQuery access).
- Refuses network downloads by default where applicable.
- Re-runs analysis scripts and checks expected outputs + key invariants.
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
from typing import Iterable

import pandas as pd


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(os.path.dirname(SCRIPT_DIR))
RESULTS_DIR = os.path.join(REPO_ROOT, "results")
PYTHON = sys.executable


def run_cmd(args: list[str], cwd: str | None = None, expect_success: bool = True) -> subprocess.CompletedProcess:
    proc = subprocess.run(
        args,
        cwd=cwd or SCRIPT_DIR,
        text=True,
        capture_output=True,
        check=False,
    )
    if proc.stdout:
        print(proc.stdout, end="" if proc.stdout.endswith("\n") else "\n")
    if proc.stderr:
        print(proc.stderr, end="" if proc.stderr.endswith("\n") else "\n")

    if expect_success and proc.returncode != 0:
        raise RuntimeError(f"Command failed ({proc.returncode}): {' '.join(args)}")
    if not expect_success and proc.returncode == 0:
        raise RuntimeError(f"Command unexpectedly succeeded: {' '.join(args)}")
    return proc


def assert_files_exist(paths: Iterable[str]) -> None:
    missing = [p for p in paths if not os.path.exists(p)]
    if missing:
        raise FileNotFoundError("Missing required files:\n" + "\n".join(missing))


def validate_external_results(csv_path: str) -> None:
    df = pd.read_csv(csv_path)
    if len(df) != 1:
        raise ValueError(f"Expected one row in {csv_path}, found {len(df)}")
    row = df.iloc[0]

    expected = {
        "n_pre_treatment": 49,
        "cv_pdl1_rho": -0.2831,
        "cv_pdl1_p": 0.048748,
        "fisher_p": 1.0,
    }
    if int(row["n_pre_treatment"]) != expected["n_pre_treatment"]:
        raise ValueError(f"Unexpected n_pre_treatment: {row['n_pre_treatment']}")
    if abs(float(row["cv_pdl1_rho"]) - expected["cv_pdl1_rho"]) > 5e-4:
        raise ValueError(f"Unexpected cv_pdl1_rho: {row['cv_pdl1_rho']}")
    if abs(float(row["cv_pdl1_p"]) - expected["cv_pdl1_p"]) > 5e-4:
        raise ValueError(f"Unexpected cv_pdl1_p: {row['cv_pdl1_p']}")
    if abs(float(row["fisher_p"]) - expected["fisher_p"]) > 1e-9:
        raise ValueError(f"Unexpected fisher_p: {row['fisher_p']}")


def validate_release_locking() -> None:
    with tempfile.TemporaryDirectory(prefix="tcga-ci-") as tmpdir:
        manifest = os.path.join(tmpdir, "tcga_extract_expanded_manifest.json")

        # Default dry-run must succeed and emit a pinned release manifest.
        run_cmd(
            [
                PYTHON,
                "tcga_extract_expanded.py",
                "--dry-run",
                "--output-dir",
                tmpdir,
            ]
        )
        if not os.path.exists(manifest):
            raise FileNotFoundError(f"Expected manifest not found: {manifest}")

        with open(manifest, "r", encoding="utf-8") as f:
            payload = json.load(f)
        release_tag = str(payload.get("tcga_release_tag", ""))
        if release_tag == "current":
            raise ValueError("Default extraction release tag must be pinned, not 'current'.")
        if not release_tag.startswith("r"):
            raise ValueError(f"Unexpected default release tag: {release_tag}")

        # current without explicit opt-in must fail.
        run_cmd(
            [
                PYTHON,
                "tcga_extract_expanded.py",
                "--dry-run",
                "--release-tag",
                "current",
                "--output-dir",
                tmpdir,
            ],
            expect_success=False,
        )

        # current with explicit opt-in must pass.
        run_cmd(
            [
                PYTHON,
                "tcga_extract_expanded.py",
                "--dry-run",
                "--release-tag",
                "current",
                "--allow-current",
                "--output-dir",
                tmpdir,
            ]
        )


def main() -> None:
    print("=" * 72)
    print("CI validation: committed TCGA outputs + extraction version lock")
    print("=" * 72)

    required_inputs = [
        os.path.join(SCRIPT_DIR, "tcga_expanded_tpm.csv"),
        os.path.join(SCRIPT_DIR, "tcga_clinical.csv"),
        os.path.join(SCRIPT_DIR, "tcga_purity_immune_covariates.csv"),
        os.path.join(SCRIPT_DIR, "geo_cache", "gse91061_target_expression.csv"),
        os.path.join(SCRIPT_DIR, "geo_cache", "gse91061_sample_meta.csv"),
    ]
    assert_files_exist(required_inputs)
    print("Input artifact presence check: PASS")

    run_cmd([PYTHON, "external_validation_geo.py", "--no-download"])
    run_cmd([PYTHON, "run_external_validation.py"])
    run_cmd([PYTHON, "sensitivity_threshold_sweep.py"])
    run_cmd([PYTHON, "sensitivity_rmst.py"])
    run_cmd([PYTHON, "sensitivity_immune_residualization.py"])
    run_cmd([PYTHON, "composite_observability_index.py"])

    required_outputs = [
        os.path.join(SCRIPT_DIR, "external_validation_results.csv"),
        os.path.join(SCRIPT_DIR, "threshold_sensitivity_results.csv"),
        os.path.join(SCRIPT_DIR, "rmst_results.csv"),
        os.path.join(SCRIPT_DIR, "immune_residualization_results.csv"),
        os.path.join(SCRIPT_DIR, "observability_index_results.csv"),
        os.path.join(RESULTS_DIR, "external_validation_geo.png"),
        os.path.join(RESULTS_DIR, "threshold_sensitivity_sweep.png"),
        os.path.join(RESULTS_DIR, "rmst_am_vs_dc.png"),
        os.path.join(RESULTS_DIR, "immune_residualization_scatter.png"),
        os.path.join(RESULTS_DIR, "purity_stratified_correlations.png"),
        os.path.join(RESULTS_DIR, "observability_index_survival.png"),
    ]
    assert_files_exist(required_outputs)
    print("Output artifact presence check: PASS")

    validate_external_results(os.path.join(SCRIPT_DIR, "external_validation_results.csv"))
    print("External validation numeric invariants: PASS")

    validate_release_locking()
    print("Extraction version-lock behavior: PASS")

    print("=" * 72)
    print("TCGA CI validation passed.")
    print("=" * 72)


if __name__ == "__main__":
    main()
