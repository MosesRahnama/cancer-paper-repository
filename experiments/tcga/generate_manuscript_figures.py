"""
Generate additional manuscript figures from existing TCGA analysis outputs.

This script reads CSV outputs already produced by the TCGA workflows and writes
publication-ready summary plots into the repository-level `results/` folder so
they can be referenced directly by the paper TeX source.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.transforms import blended_transform_factory
import numpy as np
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent.parent
RESULTS_DIR = REPO_ROOT / "results"

PROJECT_ORDER = [
    "TCGA-SKCM",
    "TCGA-LUAD",
    "TCGA-BRCA",
    "TCGA-COAD",
    "TCGA-HNSC",
    "TCGA-LUSC",
]

PROJECT_LABELS = {
    "TCGA-SKCM": "Melanoma",
    "TCGA-LUAD": "Lung Adeno",
    "TCGA-BRCA": "Breast",
    "TCGA-COAD": "Colon",
    "TCGA-HNSC": "Head & Neck",
    "TCGA-LUSC": "Lung Squam.",
}


def _bh_fdr(p_values: pd.Series) -> pd.Series:
    """Benjamini-Hochberg correction with NaN pass-through."""
    p = pd.to_numeric(p_values, errors="coerce")
    out = pd.Series(np.nan, index=p.index, dtype=float)
    valid = p.dropna()
    if valid.empty:
        return out

    order = np.argsort(valid.values)
    ranked = valid.values[order]
    n = len(ranked)
    q = ranked * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0.0, 1.0)

    q_series = pd.Series(q, index=valid.index[order])
    out.loc[q_series.index] = q_series
    return out


def _read_csv(name: str) -> pd.DataFrame:
    path = SCRIPT_DIR / name
    if not path.exists():
        raise FileNotFoundError(f"Missing required file: {path}")
    return pd.read_csv(path)


def _load_global_stage_qvalues() -> dict[str, float]:
    """
    Return global BH-FDR q-values for circadian stage Spearman tests.

    These come from master_fdr_results.csv where FDR is computed across the
    full analysis pool, matching manuscript language ("global FDR correction").
    """
    path = SCRIPT_DIR / "master_fdr_results.csv"
    if not path.exists():
        return {}

    df = pd.read_csv(path)
    required = {"source", "cancer_type", "test_description", "q_value"}
    if not required.issubset(df.columns):
        return {}

    stage = df[df["source"] == "stage_analysis"].copy()
    stage = stage[stage["test_description"] == "Circadian CV_spearman_stage_trend"].copy()
    if stage.empty:
        return {}

    stage["q_value"] = pd.to_numeric(stage["q_value"], errors="coerce")
    stage = stage.dropna(subset=["q_value"])
    return dict(zip(stage["cancer_type"], stage["q_value"]))


def plot_survival_test_qvalue_summary() -> Path:
    df = _read_csv("survival_logrank_results.csv")
    df = df[df["test"].isin(["logrank_AM_vs_DC", "logrank_Q1_vs_Q4"])].copy()
    if df.empty:
        raise RuntimeError("No survival tests available for plotting.")

    test_labels = {
        "logrank_AM_vs_DC": "Active Masking vs Decoherence",
        "logrank_Q1_vs_Q4": "Circadian CV Q1 vs Q4",
    }

    df["project"] = df["cancer_type"]
    df["project"] = pd.Categorical(df["project"], categories=PROJECT_ORDER, ordered=True)
    df = df.sort_values(["project", "test"])

    if "q_value_all_survival_tests" in df.columns:
        df["q_value"] = pd.to_numeric(df["q_value_all_survival_tests"], errors="coerce")
    else:
        df["q_value"] = _bh_fdr(df["p_value"])

    q_floor = 1e-300
    df["neglog10_q"] = -np.log10(np.clip(df["q_value"], q_floor, 1.0))

    pivot = df.pivot(index="project", columns="test", values="neglog10_q")
    pivot = pivot.reindex(PROJECT_ORDER)

    x = np.arange(len(PROJECT_ORDER))
    width = 0.38

    fig, ax = plt.subplots(figsize=(12.5, 5.5))
    bars1 = ax.bar(
        x - width / 2,
        pivot.get("logrank_AM_vs_DC", pd.Series(index=pivot.index, dtype=float)).fillna(0.0),
        width=width,
        label=test_labels["logrank_AM_vs_DC"],
        color="#2c7fb8",
        alpha=0.9,
    )
    bars2 = ax.bar(
        x + width / 2,
        pivot.get("logrank_Q1_vs_Q4", pd.Series(index=pivot.index, dtype=float)).fillna(0.0),
        width=width,
        label=test_labels["logrank_Q1_vs_Q4"],
        color="#fdae61",
        alpha=0.9,
    )

    threshold = -np.log10(0.05)
    ax.axhline(threshold, color="black", linestyle="--", linewidth=1.0, alpha=0.7)
    ax.text(
        len(PROJECT_ORDER) - 0.3,
        threshold + 0.03,
        "FDR q = 0.05",
        ha="right",
        va="bottom",
        fontsize=9,
    )

    for bars, test_name in [(bars1, "logrank_AM_vs_DC"), (bars2, "logrank_Q1_vs_Q4")]:
        sub = df[df["test"] == test_name].set_index("project")
        for b, proj in zip(bars, PROJECT_ORDER):
            q = sub["q_value"].get(proj, np.nan)
            if pd.notna(q) and q < 0.05:
                ax.text(
                    b.get_x() + b.get_width() / 2,
                    b.get_height() + 0.04,
                    "*",
                    ha="center",
                    va="bottom",
                    fontsize=13,
                    fontweight="bold",
                )

    ax.set_xticks(x)
    ax.set_xticklabels([PROJECT_LABELS[p] for p in PROJECT_ORDER], rotation=0)
    ax.set_ylabel(r"$-\log_{10}(\mathrm{FDR}\ q)$")
    ax.set_title("Survival Evidence Across Cohorts: Boundary Mode vs Circadian Quartiles")
    ax.legend(frameon=False, loc="upper right")
    ax.set_ylim(0, max(1.8, float(np.nanmax(df["neglog10_q"]) + 0.8)))
    ax.grid(axis="y", alpha=0.25)
    ax.text(
        0.01,
        0.97,
        "* = FDR-significant (q < 0.05)",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=8.8,
        bbox={"facecolor": "white", "alpha": 0.8, "edgecolor": "none"},
    )
    fig.tight_layout()

    out = RESULTS_DIR / "survival_test_qvalue_summary.png"
    fig.savefig(out, dpi=300, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    return out


def plot_robustness_hr_forest() -> Path:
    df = _read_csv("robustness_primary_tests.csv").copy()
    keep_tests = {
        "am_vs_dc_adjusted": "AM vs DC (adjusted)",
        "am_vs_dc_plus_pdl1_b2m": "AM vs DC + PD-L1/B2M",
        "pdl1_x_clock_interaction": "PD-L1 x Clock interaction",
        "high_pdl1_clock_effect": "Clock effect (high PD-L1 stratum)",
    }
    df = df[df["test_name"].isin(keep_tests)].copy()
    if df.empty:
        raise RuntimeError("No robustness test rows available for plotting.")

    df["project"] = pd.Categorical(df["project"], categories=PROJECT_ORDER, ordered=True)
    test_order = list(keep_tests.keys())
    df["test_name"] = pd.Categorical(df["test_name"], categories=test_order, ordered=True)
    df = df.sort_values(["project", "test_name"]).reset_index(drop=True)

    df["row_label"] = (
        df["project"].map(PROJECT_LABELS).astype(str)
        + " | "
        + df["test_name"].map(keep_tests).astype(str)
    )
    df["q_value"] = pd.to_numeric(df.get("q_value", np.nan), errors="coerce")
    df["p_value"] = pd.to_numeric(df["p_value"], errors="coerce")
    df["hazard_ratio"] = pd.to_numeric(df["hazard_ratio"], errors="coerce")
    df["ci95_low"] = pd.to_numeric(df["ci95_low"], errors="coerce")
    df["ci95_high"] = pd.to_numeric(df["ci95_high"], errors="coerce")

    y_pos = np.arange(len(df))[::-1]
    colors = np.where(df["q_value"] < 0.05, "#1b9e77", "#7f7f7f")
    x = df["hazard_ratio"].to_numpy()
    xerr_low = np.maximum(x - df["ci95_low"].to_numpy(), 1e-9)
    xerr_high = np.maximum(df["ci95_high"].to_numpy() - x, 1e-9)

    fig, ax = plt.subplots(figsize=(14.2, 7.2))
    for i in range(len(df)):
        ax.errorbar(
            x[i],
            y_pos[i],
            xerr=np.array([[xerr_low[i]], [xerr_high[i]]]),
            fmt="o",
            color=colors[i],
            ecolor=colors[i],
            elinewidth=1.5,
            capsize=3,
            markersize=6,
            alpha=0.95,
        )

    ax.axvline(1.0, color="black", linestyle="--", linewidth=1.0, alpha=0.8)
    ax.set_xscale("log")
    ax.set_yticks(y_pos)
    ax.set_yticklabels(df["row_label"])
    ax.tick_params(axis="y", labelsize=8, pad=5)
    ax.set_xlabel("Hazard Ratio (log scale, 95% CI)")
    ax.set_title("Clinical Robustness Models: Hazard Ratios by Test and Cohort")

    x_min = max(1e-3, float(np.nanmin(df["ci95_low"]) * 0.75))
    x_max = float(np.nanmax(df["ci95_high"]) * 1.35)
    ax.set_xlim(x_min, x_max)

    right_text = blended_transform_factory(ax.transAxes, ax.transData)
    for i, (_, row) in enumerate(df.iterrows()):
        q_txt = f"{row['q_value']:.3f}" if pd.notna(row["q_value"]) else "NA"
        p_txt = f"{row['p_value']:.3g}" if pd.notna(row["p_value"]) else "NA"
        y_offset = 0.16 if i % 2 == 0 else -0.16
        ax.text(
            1.01,
            y_pos[i] + y_offset,
            f"p={p_txt}, q={q_txt}",
            transform=right_text,
            ha="left",
            va="center",
            fontsize=8.5,
            color="#2f2f2f",
            clip_on=False,
        )

    ax.grid(axis="x", which="both", alpha=0.2)
    ax.text(
        0.01,
        0.98,
        "Green points are FDR-significant (q < 0.05).",
        transform=ax.transAxes,
        fontsize=9,
        ha="left",
        va="top",
        bbox={"facecolor": "white", "alpha": 0.85, "edgecolor": "none"},
    )
    fig.subplots_adjust(left=0.30, right=0.82, top=0.90, bottom=0.12)

    out = RESULTS_DIR / "robustness_primary_hr_forest.png"
    fig.savefig(out, dpi=300, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    return out


def plot_immune_subtype_prevalence() -> Path:
    df = _read_csv("immune_subtype_comparison.csv").copy()
    required = [
        "project",
        "n_total",
        "n_active_masking",
        "n_decoherence",
        "n_mixed",
    ]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise RuntimeError(f"Missing subtype columns: {missing}")

    df = df[df["project"].isin(PROJECT_ORDER)].copy()
    df["project"] = pd.Categorical(df["project"], categories=PROJECT_ORDER, ordered=True)
    df = df.sort_values("project")

    pct_am = 100.0 * df["n_active_masking"] / df["n_total"]
    pct_dc = 100.0 * df["n_decoherence"] / df["n_total"]
    pct_mx = 100.0 * df["n_mixed"] / df["n_total"]

    x = np.arange(len(df))
    fig, ax = plt.subplots(figsize=(11.5, 5.6))
    ax.bar(x, pct_am, color="#1b9e77", label="Active Masking")
    ax.bar(x, pct_dc, bottom=pct_am, color="#d95f02", label="Decoherence")
    ax.bar(x, pct_mx, bottom=pct_am + pct_dc, color="#7570b3", label="Mixed")

    for i, n_total in enumerate(df["n_total"].astype(int)):
        ax.text(i, 102.0, f"n={n_total}", ha="center", va="bottom", fontsize=9)

    ax.set_xticks(x)
    ax.set_xticklabels([PROJECT_LABELS[p] for p in df["project"]])
    ax.set_ylabel("Subtype prevalence (%)")
    ax.set_ylim(0, 108)
    ax.set_title("Boundary-Failure Subtype Composition by Cancer Type")
    ax.legend(frameon=False, loc="upper right")
    ax.grid(axis="y", alpha=0.2)
    fig.tight_layout()

    out = RESULTS_DIR / "immune_subtype_prevalence_stacked.png"
    fig.savefig(out, dpi=300, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    return out


def plot_stage_circadian_spearman() -> Path:
    df = _read_csv("stage_analysis_results.csv").copy()
    df = df[
        (df["metric"] == "Circadian CV")
        & (df["test"] == "spearman_stage_trend")
        & (df["cancer_type"].isin(PROJECT_ORDER))
    ].copy()
    if df.empty:
        raise RuntimeError("No circadian stage Spearman rows available.")

    df["project"] = pd.Categorical(df["cancer_type"], categories=PROJECT_ORDER, ordered=True)
    df = df.sort_values("project")
    df["p_value"] = pd.to_numeric(df["p_value"], errors="coerce")
    df["statistic"] = pd.to_numeric(df["statistic"], errors="coerce")
    global_q = _load_global_stage_qvalues()
    if global_q:
        df["q_value"] = df["cancer_type"].map(global_q)
    else:
        df["q_value"] = _bh_fdr(df["p_value"])

    y = np.arange(len(df))[::-1]
    colors = np.where(df["q_value"] < 0.05, "#2c7fb8", "#9e9e9e")

    fig, ax = plt.subplots(figsize=(10.8, 4.8))
    ax.barh(y, df["statistic"], color=colors, alpha=0.9)
    ax.axvline(0.0, color="black", linewidth=1.0, alpha=0.8)

    ax.set_yticks(y)
    ax.set_yticklabels([PROJECT_LABELS[p] for p in df["project"]])
    ax.set_xlabel("Spearman rho (stage vs circadian CV)")
    ax.set_title("Stage Association of Circadian CV Across Cancer Types")

    x_min = float(min(-0.2, np.nanmin(df["statistic"]) - 0.08))
    x_max = float(max(0.2, np.nanmax(df["statistic"]) + 0.12))
    ax.set_xlim(x_min, x_max)

    for yi, (_, row) in zip(y, df.iterrows()):
        p_txt = f"{row['p_value']:.3g}" if pd.notna(row["p_value"]) else "NA"
        q_txt = f"{row['q_value']:.3g}" if pd.notna(row["q_value"]) else "NA"
        ax.text(
            x_max - 0.01,
            yi,
            f"p={p_txt}, q={q_txt}",
            ha="right",
            va="center",
            fontsize=8.5,
        )

    ax.text(
        0.01,
        0.03,
        "Blue bars: q < 0.05 after global BH correction (master_fdr_results.csv).",
        transform=ax.transAxes,
        fontsize=8.8,
        ha="left",
        va="bottom",
    )
    ax.grid(axis="x", alpha=0.2)
    fig.tight_layout()

    out = RESULTS_DIR / "stage_circadian_spearman_summary.png"
    fig.savefig(out, dpi=300, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    return out


def main() -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    outputs = [
        plot_survival_test_qvalue_summary(),
        plot_robustness_hr_forest(),
        plot_immune_subtype_prevalence(),
        plot_stage_circadian_spearman(),
    ]
    print("Generated manuscript figures:")
    for out in outputs:
        print(f"  - {out}")


if __name__ == "__main__":
    main()
