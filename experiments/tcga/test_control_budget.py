"""
Test the Control Budget hypothesis: C_G(G) + C_S(S) <= B_ctrl

Approach: Pareto Frontier Analysis
- Define G (Proliferation Score) and S (Stability/Differentiation Score)
- Test for inverse boundary constraint (high G + high S should be impossible)
- Use quantile regression to detect the upper boundary
"""
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from scipy.stats import spearmanr
import warnings
warnings.filterwarnings('ignore')

# Add current directory to path for config import
sys.path.insert(0, os.path.dirname(__file__))
from tcga_config import DATA_DIR, FIGURE_DIR, ALL_PROJECTS, log_transform, setup_plotting

# ============================================================================
# Define Metrics
# ============================================================================

def compute_proliferation_score(df):
    """
    Compute proliferation score from available genes.

    Our current gene set includes:
    - MYC: master proliferation regulator
    - TP53: when lost/mutated, removes growth checkpoint

    We'll use MYC as primary proliferation proxy (log-transformed TPM).
    Note: Ideally would include MKI67, PCNA, TOP2A, FOXM1 for complete score.
    """
    # For now, use MYC as proliferation proxy
    if 'MYC' in df.columns:
        prolif = log_transform(df['MYC'])
    else:
        raise ValueError("MYC not found in dataset")

    return prolif


def compute_differentiation_score(df, cancer_type=None):
    """
    Compute differentiation/stability score.

    Differentiation = organized, stable, identity-bearing state

    Proxies from current gene set:
    - CDH1 (E-cadherin): epithelial differentiation marker, cell adhesion
    - Gap junctions (GJA1, GJB2, GJA5, GJB6): intercellular coupling, tissue organization
    - Circadian coherence: inverse of circadian CV (organized temporal state)

    Higher differentiation = higher CDH1 + higher gap junction coupling + lower circadian CV
    """
    components = []

    # 1. E-cadherin (differentiation marker)
    if 'CDH1' in df.columns:
        components.append(log_transform(df['CDH1']))

    # 2. Gap junction coupling (average of available connexins)
    gj_genes = ['GJA1', 'GJB2', 'GJA5', 'GJB6']
    gj_present = [g for g in gj_genes if g in df.columns]
    if gj_present:
        gj_mean = df[gj_present].apply(log_transform).mean(axis=1)
        components.append(gj_mean)

    # 3. Circadian coherence (inverse of CV - more coherent = more differentiated/stable)
    circ_genes = ['ARNTL', 'CLOCK', 'PER1', 'PER2', 'CRY1', 'CRY2']
    circ_present = [g for g in circ_genes if g in df.columns]
    if circ_present and len(circ_present) >= 4:
        circ_vals = df[circ_present].apply(log_transform)
        circ_cv = circ_vals.std(axis=1) / circ_vals.mean(axis=1)
        # Invert CV: high coherence (low CV) = high differentiation
        # Normalize to similar scale as other components
        circ_coherence = -circ_cv  # Negative CV, so high coherence = high value
        components.append(circ_coherence)

    if not components:
        raise ValueError("No differentiation markers found in dataset")

    # Combine components (mean of normalized scores)
    diff_score = pd.DataFrame({f'comp_{i}': c for i, c in enumerate(components)})

    # Z-score normalize each component
    diff_score_norm = (diff_score - diff_score.mean()) / diff_score.std()

    # Average
    diff_final = diff_score_norm.mean(axis=1)

    return diff_final


def compute_total_load(prolif, diff):
    """
    Total resource demand = normalized proliferation + normalized differentiation.

    If budget constraint exists, high total_load should correlate with metabolic stress.
    """
    # Z-score normalize both
    prolif_norm = (prolif - prolif.mean()) / prolif.std()
    diff_norm = (diff - diff.mean()) / diff.std()

    # Total load
    total = prolif_norm + diff_norm

    return total


def quantile_regression_boundary(x, y, quantile=0.95):
    """
    Fit quantile regression to detect upper boundary.

    Use numpy percentile-based approach since we don't have statsmodels.

    Returns: slope, intercept of the boundary line
    """
    # Sort by x
    sorted_idx = np.argsort(x)
    x_sorted = x[sorted_idx]
    y_sorted = y[sorted_idx]

    # Sliding window to find 95th percentile at each x position
    window_size = max(50, int(len(x) * 0.1))  # At least 50 points or 10% of data

    x_boundary = []
    y_boundary = []

    # Create bins along x-axis
    n_bins = 20
    x_bins = np.linspace(x.min(), x.max(), n_bins + 1)

    for i in range(n_bins):
        mask = (x >= x_bins[i]) & (x < x_bins[i+1])
        if mask.sum() > 10:  # Need at least 10 points
            y_quantile = np.percentile(y[mask], quantile * 100)
            x_center = (x_bins[i] + x_bins[i+1]) / 2
            x_boundary.append(x_center)
            y_boundary.append(y_quantile)

    if len(x_boundary) < 3:
        return None, None, None, None

    x_boundary = np.array(x_boundary)
    y_boundary = np.array(y_boundary)

    # Fit linear regression to boundary points
    slope, intercept, r_value, p_value, std_err = stats.linregress(x_boundary, y_boundary)

    return slope, intercept, p_value, (x_boundary, y_boundary)


# ============================================================================
# Main Analysis
# ============================================================================

def main():
    print("="*80)
    print("CONTROL BUDGET HYPOTHESIS TEST")
    print("="*80)
    print("\nTesting: C_G(G) + C_S(S) <= B_ctrl")
    print("Prediction: High proliferation + High differentiation is IMPOSSIBLE")
    print("Expected: Negative boundary in G vs S scatter plot\n")

    # Load data
    expr_file = os.path.join(DATA_DIR, "tcga_expanded_tpm.csv")
    if not os.path.exists(expr_file):
        print(f"ERROR: {expr_file} not found")
        print("Run tcga_extract_expanded.py first")
        return

    df = pd.read_csv(expr_file)
    print(f"Loaded {len(df)} samples from {expr_file}")

    # Filter to tumor samples only
    tumor_mask = df['sample_type_name'].str.contains('Tumor|Metastatic', case=False, na=False)
    df_tumor = df[tumor_mask].copy()
    print(f"Filtered to {len(df_tumor)} tumor samples\n")

    # Setup plotting
    setup_plotting()

    # Results storage
    results = []

    # ========================================================================
    # Analysis 1: Pareto Frontier per Cancer Type
    # ========================================================================

    print("="*80)
    print("ANALYSIS 1: PARETO FRONTIER TEST (Per Cancer Type)")
    print("="*80)

    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    axes = axes.flatten()

    for idx, cancer in enumerate(ALL_PROJECTS):
        df_cancer = df_tumor[df_tumor['project_short_name'] == cancer].copy()

        if len(df_cancer) < 50:
            print(f"\nSkipping {cancer}: too few samples ({len(df_cancer)})")
            continue

        print(f"\n{cancer} (n={len(df_cancer)})")
        print("-" * 40)

        # Compute scores
        try:
            prolif = compute_proliferation_score(df_cancer)
            diff = compute_differentiation_score(df_cancer, cancer)
        except Exception as e:
            print(f"  ERROR: {e}")
            continue

        # Test for negative correlation (should exist if budget constrains both)
        rho, p_val = spearmanr(prolif, diff)
        print(f"  Spearman correlation (G vs S): rho = {rho:.3f}, p = {p_val:.2e}")

        # Quantile regression on 95th percentile
        slope, intercept, p_boundary, boundary_points = quantile_regression_boundary(
            prolif.values, diff.values, quantile=0.95
        )

        if slope is not None:
            print(f"  95th percentile boundary slope: {slope:.3f} (p = {p_boundary:.2e})")
            if slope < 0 and p_boundary < 0.05:
                print(f"  [+] NEGATIVE BOUNDARY DETECTED: Budget constraint supported")
            else:
                print(f"  [-] No clear negative boundary")

        # Store results
        results.append({
            'cancer_type': cancer,
            'n_samples': len(df_cancer),
            'correlation_rho': rho,
            'correlation_p': p_val,
            'boundary_slope': slope,
            'boundary_p': p_boundary,
            'budget_supported': (slope is not None and slope < 0 and p_boundary < 0.05)
        })

        # Plot
        ax = axes[idx]
        ax.scatter(prolif, diff, alpha=0.5, s=20, c='steelblue', edgecolors='none')

        # Add boundary line if detected
        if slope is not None and boundary_points is not None:
            x_bound, y_bound = boundary_points
            ax.plot(x_bound, y_bound, 'r--', linewidth=2, alpha=0.7,
                   label=f'95% boundary\nslope={slope:.2f}')

            # Extend line for visualization
            x_range = np.array([prolif.min(), prolif.max()])
            y_pred = slope * x_range + intercept
            ax.plot(x_range, y_pred, 'r-', linewidth=1, alpha=0.3)

        ax.set_xlabel('Proliferation Score (G)', fontsize=10)
        ax.set_ylabel('Differentiation Score (S)', fontsize=10)
        ax.set_title(f'{cancer}\nρ={rho:.2f}, slope={slope:.2f}' if slope else f'{cancer}\nρ={rho:.2f}',
                    fontsize=11)
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)

    # Remove empty subplots
    for idx in range(len(ALL_PROJECTS), len(axes)):
        fig.delaxes(axes[idx])

    plt.tight_layout()
    output_file = os.path.join(FIGURE_DIR, "control_budget_pareto_frontier.png")
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\nSaved: {output_file}")
    plt.close()

    # ========================================================================
    # Analysis 2: Combined Analysis (All Cancer Types)
    # ========================================================================

    print("\n" + "="*80)
    print("ANALYSIS 2: COMBINED PARETO FRONTIER (All Cancers)")
    print("="*80)

    # Compute scores for all tumors
    prolif_all = compute_proliferation_score(df_tumor)
    diff_all = compute_differentiation_score(df_tumor)

    # Overall correlation
    rho_all, p_all = spearmanr(prolif_all, diff_all)
    print(f"\nOverall Spearman correlation: rho = {rho_all:.3f}, p = {p_all:.2e}")

    # Quantile regression
    slope_all, intercept_all, p_boundary_all, boundary_all = quantile_regression_boundary(
        prolif_all.values, diff_all.values, quantile=0.95
    )

    if slope_all is not None:
        print(f"95th percentile boundary slope: {slope_all:.3f} (p = {p_boundary_all:.2e})")
        if slope_all < 0 and p_boundary_all < 0.05:
            print("[+] NEGATIVE BOUNDARY DETECTED: Strong evidence for budget constraint")
        else:
            print("[-] No clear negative boundary")

    # Plot combined
    fig, ax = plt.subplots(1, 1, figsize=(10, 8))

    # Color by cancer type
    for cancer in ALL_PROJECTS:
        mask = df_tumor['project_short_name'] == cancer
        ax.scatter(prolif_all[mask], diff_all[mask],
                  label=cancer.replace('TCGA-', ''),
                  alpha=0.6, s=30, edgecolors='none')

    # Add boundary
    if slope_all is not None and boundary_all is not None:
        x_bound, y_bound = boundary_all
        ax.plot(x_bound, y_bound, 'k--', linewidth=3, alpha=0.8,
               label=f'95% boundary (slope={slope_all:.2f})')

    ax.set_xlabel('Proliferation Score (G)', fontsize=12)
    ax.set_ylabel('Differentiation/Stability Score (S)', fontsize=12)
    ax.set_title(f'Control Budget: Pareto Frontier Test\n' +
                f'n={len(df_tumor)}, ρ={rho_all:.3f}, boundary slope={slope_all:.2f}',
                fontsize=13, weight='bold')
    ax.legend(fontsize=9, loc='best')
    ax.grid(True, alpha=0.3)

    output_file = os.path.join(FIGURE_DIR, "control_budget_combined.png")
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\nSaved: {output_file}")
    plt.close()

    # ========================================================================
    # Analysis 3: Total Load vs Metabolic Stress Test
    # ========================================================================

    print("\n" + "="*80)
    print("ANALYSIS 3: TOTAL LOAD vs METABOLIC STRESS")
    print("="*80)
    print("\nTesting: High (G+S) -> Metabolic stress/burden")
    print("Proxies for stress: VIM (mesenchymal/stress), low CDH1 (loss of organization)")

    # Compute total load
    total_load = compute_total_load(prolif_all, diff_all)

    # Metabolic stress proxies from available genes
    if 'VIM' in df_tumor.columns:
        vim = log_transform(df_tumor['VIM'])
        rho_vim, p_vim = spearmanr(total_load, vim)
        print(f"\nTotal Load vs VIM (stress/mesenchymal): rho = {rho_vim:.3f}, p = {p_vim:.2e}")

    # High load should correlate with PD-L1 (immune evasion = stress response)
    if 'CD274' in df_tumor.columns:
        pdl1 = log_transform(df_tumor['CD274'])
        rho_pdl1, p_pdl1 = spearmanr(total_load, pdl1)
        print(f"Total Load vs PD-L1 (immune stress): rho = {rho_pdl1:.3f}, p = {p_pdl1:.2e}")

    # ========================================================================
    # Summary Table
    # ========================================================================

    print("\n" + "="*80)
    print("SUMMARY: CONTROL BUDGET EVIDENCE")
    print("="*80)

    df_results = pd.DataFrame(results)

    print("\nPer-Cancer Results:")
    print(df_results.to_string(index=False))

    n_supported = df_results['budget_supported'].sum()
    n_total = len(df_results)

    print(f"\n{'='*80}")
    print(f"BUDGET CONSTRAINT SUPPORTED: {n_supported}/{n_total} cancer types")
    print(f"Overall evidence: {'STRONG' if n_supported >= 4 else 'MODERATE' if n_supported >= 2 else 'WEAK'}")
    print(f"{'='*80}")

    # Save results
    results_file = os.path.join(DATA_DIR, "control_budget_test_results.csv")
    df_results.to_csv(results_file, index=False)
    print(f"\nResults saved to: {results_file}")


if __name__ == "__main__":
    main()
