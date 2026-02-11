"""
Synchronize manuscript figure artifacts and emit a provenance manifest.

Goals:
1) Keep manuscript-facing figures in `results/`.
2) Eliminate stale duplicate drift for shared basenames between:
   - `results/`
   - `experiments/tcga/figures/`
3) Write a provenance CSV with figure path, source inputs, generator script,
   timestamp, and SHA-256 hash for auditability.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import re
import shutil
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent.parent
PAPER_TEX = REPO_ROOT / "paper" / "Cancer_As_Boundary_Logic_Failure.tex"
RESULTS_DIR = REPO_ROOT / "results"
TCGA_FIGURES_DIR = SCRIPT_DIR / "figures"
PROVENANCE_CSV = REPO_ROOT / "paper" / "figure_provenance.csv"


INCLUDEGRAPHICS_PATTERN = re.compile(
    r"\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}"
)


@dataclass
class SyncEvent:
    basename: str
    source: str
    copied_to_results: bool
    copied_to_figures: bool


FIGURE_METADATA: dict[str, dict[str, str]] = {
    "control_budget_combined.png": {
        "source_inputs": "experiments/tcga/control_budget_test_results.csv",
        "generator_script": "experiments/tcga/test_control_budget.py",
    },
    "multicancer_correlation_heatmap.png": {
        "source_inputs": "experiments/tcga/tcga_multicancer_correlations.csv",
        "generator_script": "experiments/tcga/tcga_multicancer.py",
    },
    "multicancer_circadian_cv_forest.png": {
        "source_inputs": "experiments/tcga/tcga_multicancer_correlations.csv",
        "generator_script": "experiments/tcga/tcga_multicancer.py",
    },
    "hypothesis1_TCGA_SKCM.png": {
        "source_inputs": "experiments/tcga/hypothesis1_correlations.csv",
        "generator_script": "experiments/tcga/tcga_analysis.py",
    },
    "hypothesis1_TCGA_LUAD.png": {
        "source_inputs": "experiments/tcga/hypothesis1_correlations.csv",
        "generator_script": "experiments/tcga/tcga_analysis.py",
    },
    "tumor_normal_circadian_cv_paired.png": {
        "source_inputs": "experiments/tcga/tumor_normal_comparison.csv",
        "generator_script": "experiments/tcga/tcga_tumor_normal.py",
    },
    "tumor_normal_circadian_cv_summary.png": {
        "source_inputs": "experiments/tcga/tumor_normal_comparison.csv",
        "generator_script": "experiments/tcga/tcga_tumor_normal.py",
    },
    "immune_subtype_circadian_cv.png": {
        "source_inputs": "experiments/tcga/immune_subtype_comparison.csv",
        "generator_script": "experiments/tcga/tcga_immune_subtype.py",
    },
    "immune_subtype_boundary_scatter.png": {
        "source_inputs": "experiments/tcga/immune_subtype_comparison.csv",
        "generator_script": "experiments/tcga/tcga_immune_subtype.py",
    },
    "immune_subtype_prevalence_stacked.png": {
        "source_inputs": "experiments/tcga/immune_subtype_comparison.csv",
        "generator_script": "experiments/tcga/generate_manuscript_figures.py",
    },
    "survival_boundary_failure.png": {
        "source_inputs": "experiments/tcga/survival_logrank_results.csv",
        "generator_script": "experiments/tcga/tcga_survival.py",
    },
    "survival_circadian_quartile.png": {
        "source_inputs": "experiments/tcga/survival_logrank_results.csv",
        "generator_script": "experiments/tcga/tcga_survival.py",
    },
    "survival_test_qvalue_summary.png": {
        "source_inputs": "experiments/tcga/survival_logrank_results.csv",
        "generator_script": "experiments/tcga/generate_manuscript_figures.py",
    },
    "stage_circadian_cv_boxplot.png": {
        "source_inputs": (
            "experiments/tcga/stage_analysis_results.csv;"
            "experiments/tcga/master_fdr_results.csv"
        ),
        "generator_script": "experiments/tcga/tcga_stage_analysis.py",
    },
    "stage_circadian_spearman_summary.png": {
        "source_inputs": (
            "experiments/tcga/stage_analysis_results.csv;"
            "experiments/tcga/master_fdr_results.csv"
        ),
        "generator_script": "experiments/tcga/generate_manuscript_figures.py",
    },
    "robustness_primary_hr_forest.png": {
        "source_inputs": "experiments/tcga/robustness_primary_tests.csv",
        "generator_script": "experiments/tcga/generate_manuscript_figures.py",
    },
    "therapy_comparison.png": {
        "source_inputs": "results/scenario_comparison.json",
        "generator_script": "examples/run_full_scenario.py",
    },
    "kappa_trajectory.png": {
        "source_inputs": "results/scenario_comparison.json",
        "generator_script": "examples/run_full_scenario.py",
    },
    "coherence_vs_cancer.png": {
        "source_inputs": "results/scenario_comparison.json",
        "generator_script": "run_simulation.py",
    },
    "state_distribution.png": {
        "source_inputs": "results/scenario_comparison.json",
        "generator_script": "run_simulation.py",
    },
    "tissue_network.png": {
        "source_inputs": "results/scenario_comparison.json",
        "generator_script": "run_simulation.py",
    },
    "phase2_baseline_kappa_summary.png": {
        "source_inputs": "results/phase2_baseline_seed_sweep.json",
        "generator_script": "run_phase2_extractions.py",
    },
    "phase2_unmasking_cancer_summary.png": {
        "source_inputs": "results/phase2_unmasking_seed_sweep.json",
        "generator_script": "run_phase2_extractions.py",
    },
}


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _same_file(a: Path, b: Path) -> bool:
    return a.exists() and b.exists() and _sha256(a) == _sha256(b)


def _select_sync_source(results_path: Path, figures_path: Path) -> Path | None:
    if results_path.exists() and figures_path.exists():
        if _same_file(results_path, figures_path):
            return results_path
        r_mtime = results_path.stat().st_mtime
        f_mtime = figures_path.stat().st_mtime
        if r_mtime >= f_mtime:
            return results_path
        return figures_path
    if results_path.exists():
        return results_path
    if figures_path.exists():
        return figures_path
    return None


def _copy_if_needed(src: Path, dst: Path) -> bool:
    if dst.exists() and _same_file(src, dst):
        return False
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)
    return True


def _parse_includegraphics_paths(tex_path: Path) -> list[str]:
    content = tex_path.read_text(encoding="utf-8")
    return [m.group(1).replace("\\", "/") for m in INCLUDEGRAPHICS_PATTERN.finditer(content)]


def _rewrite_tex_to_results(tex_path: Path) -> int:
    content = tex_path.read_text(encoding="utf-8")
    new_content = re.sub(
        r"(\\includegraphics(?:\[[^\]]*\])?\{)experiments/tcga/figures/([^}]+\})",
        r"\1results/\2",
        content,
    )
    if new_content == content:
        return 0
    tex_path.write_text(new_content, encoding="utf-8", newline="\n")
    return 1


def _sync_shared_manuscript_figures(include_paths: list[str]) -> list[SyncEvent]:
    events: list[SyncEvent] = []
    for include_path in include_paths:
        if not (
            include_path.startswith("results/")
            or include_path.startswith("experiments/tcga/figures/")
        ):
            continue

        basename = Path(include_path).name
        results_path = RESULTS_DIR / basename
        figures_path = TCGA_FIGURES_DIR / basename
        sync_src = _select_sync_source(results_path, figures_path)
        if sync_src is None:
            continue

        copied_to_results = _copy_if_needed(sync_src, results_path)
        copied_to_figures = False
        # Keep mirrors in lockstep only for TCGA figure namespace.
        if figures_path.exists() or include_path.startswith("experiments/tcga/figures/"):
            copied_to_figures = _copy_if_needed(sync_src, figures_path)

        src_name = (
            "results"
            if sync_src.resolve() == results_path.resolve()
            else "experiments/tcga/figures"
        )
        events.append(
            SyncEvent(
                basename=basename,
                source=src_name,
                copied_to_results=copied_to_results,
                copied_to_figures=copied_to_figures,
            )
        )
    return events


def _resolve_figure_file(include_path: str) -> Path:
    include_path = include_path.replace("\\", "/")
    if include_path.startswith("results/"):
        return REPO_ROOT / include_path
    if include_path.startswith("experiments/tcga/figures/"):
        return REPO_ROOT / include_path
    # Fallback for other relative paths in TeX.
    return PAPER_TEX.parent / include_path


def _write_provenance(include_paths: list[str], output_csv: Path) -> int:
    rows = []
    for include_path in include_paths:
        figure_file = _resolve_figure_file(include_path)
        if not figure_file.exists():
            continue
        stat = figure_file.stat()
        timestamp = datetime.fromtimestamp(stat.st_mtime, tz=timezone.utc).isoformat()
        basename = figure_file.name
        metadata = FIGURE_METADATA.get(
            basename,
            {"source_inputs": "UNKNOWN", "generator_script": "UNKNOWN"},
        )
        rows.append(
            {
                "figure_path": include_path,
                "source_inputs": metadata["source_inputs"],
                "generator_script": metadata["generator_script"],
                "figure_timestamp_utc": timestamp,
                "figure_sha256": _sha256(figure_file),
            }
        )

    output_csv.parent.mkdir(parents=True, exist_ok=True)
    with output_csv.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "figure_path",
                "source_inputs",
                "generator_script",
                "figure_timestamp_utc",
                "figure_sha256",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)
    return len(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--rewrite-tex",
        action="store_true",
        help="Rewrite manuscript figure paths from experiments/tcga/figures/* to results/*.",
    )
    parser.add_argument(
        "--provenance-out",
        type=Path,
        default=PROVENANCE_CSV,
        help=f"Provenance CSV output path (default: {PROVENANCE_CSV}).",
    )
    args = parser.parse_args()

    rewritten = _rewrite_tex_to_results(PAPER_TEX) if args.rewrite_tex else 0
    include_paths = _parse_includegraphics_paths(PAPER_TEX)
    sync_events = _sync_shared_manuscript_figures(include_paths)
    n_rows = _write_provenance(include_paths, args.provenance_out)

    print("Manuscript artifact sync completed.")
    if rewritten:
        print("  - Rewrote manuscript figure paths to results/: yes")
    else:
        print("  - Rewrote manuscript figure paths to results/: no")
    print(f"  - Synced figure entries inspected: {len(sync_events)}")
    for evt in sync_events:
        print(
            f"    * {evt.basename}: source={evt.source}, "
            f"results_updated={evt.copied_to_results}, figures_updated={evt.copied_to_figures}"
        )
    print(f"  - Wrote provenance rows: {n_rows}")
    print(f"  - Provenance file: {args.provenance_out}")


if __name__ == "__main__":
    main()
