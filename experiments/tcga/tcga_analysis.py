"""
Hypothesis 1 Correlation Analysis
===================================
Tests the boundary-logic prediction:
  PD-L1 expression should inversely correlate with tissue-coherence proxies
  (gap junction density and circadian gene coherence).

Data: TCGA-SKCM (melanoma) + TCGA-LUAD (lung adenocarcinoma)
Source: ISB-CGC BigQuery, extracted to gs://kernel-o6/tcga_cancer_boundary/data/
"""
import argparse
import os
import sys
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(__file__))
from tcga_config import DATA_DIR, FIGURE_DIR, upload_to_gcs  # noqa: E402

OUTPUT_DIR = FIGURE_DIR

# ── Define gene groups ─────────────────────────────────────────────────────
CHECKPOINT = ["CD274", "PDCD1LG2", "PDCD1"]
MHC_I = ["HLA_A", "HLA_B", "HLA_C", "B2M"]
GAP_JUNCTION = ["GJA1", "GJB2", "GJA5", "GJB6"]
CIRCADIAN = ["ARNTL", "CLOCK", "PER1", "PER2", "CRY1", "CRY2"]
DIFFERENTIATION = ["CDH1", "VIM", "MYC"]


def log_transform(series):
    """Log2(TPM + 1) transform."""
    return np.log2(series + 1)


def compute_circadian_coherence(df_subset):
    """
    Compute per-sample circadian coherence as coefficient of variation
    across circadian genes. Lower CV = more coherent expression.
    """
    circ_cols = [c for c in CIRCADIAN if c in df_subset.columns]
    circ_vals = df_subset[circ_cols].apply(log_transform)
    cv = circ_vals.std(axis=1) / circ_vals.mean(axis=1)
    return cv


def run_correlations(df_subset, project_name):
    """Run Spearman correlations for Hypothesis 1."""
    results = []

    pdl1 = log_transform(df_subset["CD274"])

    # PD-L1 vs each gap junction gene
    for gj in GAP_JUNCTION:
        if gj in df_subset.columns:
            gj_expr = log_transform(df_subset[gj])
            mask = pdl1.notna() & gj_expr.notna()
            r, p = stats.spearmanr(pdl1[mask], gj_expr[mask])
            results.append({
                "project": project_name,
                "comparison": f"CD274 vs {gj}",
                "category": "PD-L1 vs Gap Junction",
                "rho": r, "p_value": p, "n": mask.sum()
            })

    # PD-L1 vs each circadian gene
    for circ in CIRCADIAN:
        if circ in df_subset.columns:
            circ_expr = log_transform(df_subset[circ])
            mask = pdl1.notna() & circ_expr.notna()
            r, p = stats.spearmanr(pdl1[mask], circ_expr[mask])
            results.append({
                "project": project_name,
                "comparison": f"CD274 vs {circ}",
                "category": "PD-L1 vs Circadian",
                "rho": r, "p_value": p, "n": mask.sum()
            })

    # PD-L1 vs circadian coherence (CV)
    circ_cv = compute_circadian_coherence(df_subset)
    mask = pdl1.notna() & circ_cv.notna()
    r, p = stats.spearmanr(pdl1[mask], circ_cv[mask])
    results.append({
        "project": project_name,
        "comparison": "CD274 vs Circadian CV",
        "category": "PD-L1 vs Circadian Coherence",
        "rho": r, "p_value": p, "n": mask.sum()
    })

    # PD-L1 vs MHC-I components
    for mhc in MHC_I:
        if mhc in df_subset.columns:
            mhc_expr = log_transform(df_subset[mhc])
            mask = pdl1.notna() & mhc_expr.notna()
            r, p = stats.spearmanr(pdl1[mask], mhc_expr[mask])
            results.append({
                "project": project_name,
                "comparison": f"CD274 vs {mhc}",
                "category": "PD-L1 vs MHC-I",
                "rho": r, "p_value": p, "n": mask.sum()
            })

    # Gap junction composite vs MHC-I composite
    gj_cols = [c for c in GAP_JUNCTION if c in df_subset.columns]
    mhc_cols = [c for c in MHC_I if c in df_subset.columns]
    gj_mean = df_subset[gj_cols].apply(log_transform).mean(axis=1)
    mhc_mean = df_subset[mhc_cols].apply(log_transform).mean(axis=1)
    mask = gj_mean.notna() & mhc_mean.notna()
    r, p = stats.spearmanr(gj_mean[mask], mhc_mean[mask])
    results.append({
        "project": project_name,
        "comparison": "GapJunction_mean vs MHC_I_mean",
        "category": "Gap Junction vs MHC-I",
        "rho": r, "p_value": p, "n": mask.sum()
    })

    # E-cadherin vs VIM (EMT proxy)
    if "CDH1" in df_subset.columns and "VIM" in df_subset.columns:
        cdh1 = log_transform(df_subset["CDH1"])
        vim = log_transform(df_subset["VIM"])
        mask = cdh1.notna() & vim.notna()
        r, p = stats.spearmanr(cdh1[mask], vim[mask])
        results.append({
            "project": project_name,
            "comparison": "CDH1 vs VIM",
            "category": "EMT axis",
            "rho": r, "p_value": p, "n": mask.sum()
        })

    return pd.DataFrame(results)


def plot_hypothesis1(df_subset, project_name):
    """Generate Figure for Hypothesis 1."""
    fig, axes = plt.subplots(2, 3, figsize=(16, 10))
    fig.suptitle(
        f"Hypothesis 1: PD-L1 vs Tissue Coherence Proxies — {project_name}",
        fontsize=14, fontweight="bold"
    )

    pdl1 = log_transform(df_subset["CD274"])

    panels = [
        ("GJA1", "Connexin-43 (GJA1)", axes[0, 0]),
        ("GJB2", "Connexin-26 (GJB2)", axes[0, 1]),
        ("ARNTL", "BMAL1 (ARNTL)", axes[0, 2]),
        ("CLOCK", "CLOCK", axes[1, 0]),
        ("PER1", "PER1", axes[1, 1]),
    ]

    for gene, label, ax in panels:
        if gene not in df_subset.columns:
            continue
        y = log_transform(df_subset[gene])
        mask = pdl1.notna() & y.notna()
        r, p = stats.spearmanr(pdl1[mask], y[mask])

        ax.scatter(pdl1[mask], y[mask], alpha=0.35, s=12, c="steelblue", edgecolors="none")
        ax.set_xlabel("log2(CD274 + 1)")
        ax.set_ylabel(f"log2({label} + 1)")

        sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
        ax.set_title(f"rho ={r:.3f}, p = {p:.2e} {sig}", fontsize=10)

        # Add trend line
        z = np.polyfit(pdl1[mask], y[mask], 1)
        xline = np.linspace(pdl1[mask].min(), pdl1[mask].max(), 100)
        ax.plot(xline, np.polyval(z, xline), "r-", alpha=0.6, linewidth=1.5)

    # Circadian coherence panel
    ax = axes[1, 2]
    circ_cv = compute_circadian_coherence(df_subset)
    mask = pdl1.notna() & circ_cv.notna()
    r, p = stats.spearmanr(pdl1[mask], circ_cv[mask])
    ax.scatter(pdl1[mask], circ_cv[mask], alpha=0.35, s=12, c="darkorange", edgecolors="none")
    ax.set_xlabel("log2(CD274 + 1)")
    ax.set_ylabel("Circadian CV (lower = more coherent)")
    sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
    ax.set_title(f"rho ={r:.3f}, p = {p:.2e} {sig}", fontsize=10)
    z = np.polyfit(pdl1[mask], circ_cv[mask], 1)
    xline = np.linspace(pdl1[mask].min(), pdl1[mask].max(), 100)
    ax.plot(xline, np.polyval(z, xline), "r-", alpha=0.6, linewidth=1.5)

    plt.tight_layout()
    out_path = os.path.join(OUTPUT_DIR, f"hypothesis1_{project_name.replace('-', '_')}.png")
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out_path}")
    return out_path


def main():
    parser = argparse.ArgumentParser(description="Run 2-cohort Hypothesis 1 correlation analysis.")
    parser.add_argument(
        "--data-dir",
        default=DATA_DIR,
        help="Directory containing tcga_expression_tpm.csv and output CSV.",
    )
    parser.add_argument(
        "--figure-dir",
        default=FIGURE_DIR,
        help="Directory to write figure PNG files.",
    )
    parser.add_argument(
        "--no-gcs",
        action="store_true",
        help="Skip GCS upload and keep outputs local only.",
    )
    args = parser.parse_args()

    os.makedirs(args.figure_dir, exist_ok=True)

    df = pd.read_csv(os.path.join(args.data_dir, "tcga_expression_tpm.csv"))
    print(f"Loaded {len(df)} samples")
    print(f"Projects: {df['project_short_name'].value_counts().to_dict()}")

    global OUTPUT_DIR
    OUTPUT_DIR = args.figure_dir

    all_results = []

    for proj in ["TCGA-SKCM", "TCGA-LUAD"]:
        subset = df[df["project_short_name"] == proj].copy()
        # Focus on tumor samples for primary analysis
        tumor = subset[subset["sample_type_name"].str.contains("Tumor|Metastatic", case=False, na=False)]
        print(f"\n{'='*60}")
        print(f"  {proj}: {len(tumor)} tumor samples")
        print(f"{'='*60}")

        # Run correlations
        res = run_correlations(tumor, proj)
        all_results.append(res)

        # Print results
        for _, row in res.iterrows():
            sig = "***" if row["p_value"] < 0.001 else "**" if row["p_value"] < 0.01 else "*" if row["p_value"] < 0.05 else "ns"
            print(f"  {row['comparison']:35s}  rho={row['rho']:+.3f}  p={row['p_value']:.2e}  n={row['n']:4d}  {sig}")

        # Generate plot
        plot_hypothesis1(tumor, proj)

    # Combine and save
    results_df = pd.concat(all_results, ignore_index=True)
    results_path = os.path.join(args.data_dir, "hypothesis1_correlations.csv")
    results_df.to_csv(results_path, index=False)
    print(f"\nAll correlations saved to: {results_path}")

    if not args.no_gcs:
        try:
            upload_to_gcs(results_path, "analysis/hypothesis1_correlations.csv")
            for fig_file in os.listdir(args.figure_dir):
                if fig_file.endswith(".png"):
                    local_path = os.path.join(args.figure_dir, fig_file)
                    upload_to_gcs(local_path, f"figures/{fig_file}")
        except Exception as e:
            print(f"  GCS upload failed (non-fatal): {e}")

    print("\nDone.")


if __name__ == "__main__":
    main()
