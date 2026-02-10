"""
Full Scenario: Healthy -> Cancer -> Therapy -> Recovery
========================================================

Generates plots and saves simulation results to results/ folder.

Usage: python examples/run_full_scenario.py
"""

import sys
import os
import json
import argparse

# Ensure src is on path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from simulations.engine import Simulation
from simulations.parameter_profiles import CALIBRATED_MONTHLY_PROFILE
from simulations.therapy import (
    apply_therapy_to_all_cancer,
    differentiation_therapy,
    forced_distinction,
    checkpoint_inhibitor,
    bioelectric_reprogramming,
    combination_checkpoint_differentiation,
)

RESULTS_DIR = os.path.join(os.path.dirname(__file__), "..", "results")
os.makedirs(RESULTS_DIR, exist_ok=True)


def run_scenario(
    therapy_name,
    therapy_fn=None,
    seed=42,
    run_kwargs=None,
    n_cells=100,
    energy_budget=10000.0,
    pre_steps=200,
    post_steps=100,
):
    """Run a complete scenario with a specific therapy."""
    if run_kwargs is None:
        run_kwargs = CALIBRATED_MONTHLY_PROFILE.run_kwargs()
    sim = Simulation(n_cells=n_cells, energy_budget=energy_budget, seed=seed)

    # Phase 1: Cancer progression
    sim.run(steps=pre_steps, **run_kwargs)
    pre_therapy = dict(sim.history[-1])

    # Phase 2: Apply therapy
    if therapy_fn is None:
        treated = 0
    else:
        treated = apply_therapy_to_all_cancer(sim.tissue, therapy_fn)

    # Phase 3: Post-therapy evolution
    sim.run(steps=post_steps, **run_kwargs)
    post_therapy = dict(sim.history[-1])

    return {
        "therapy": therapy_name,
        "treated_cells": treated,
        "pre_therapy": pre_therapy,
        "post_therapy": post_therapy,
        "kappa_series": sim.kappa_series(),
        "cancer_series": sim.cancer_series(),
        "pre_cancer_series": [h["pre_cancer"] for h in sim.history],
        "energy_series": sim.energy_series(),
    }


def main():
    parser = argparse.ArgumentParser(
        description="Run full therapy comparison scenario and generate plots."
    )
    parser.add_argument("--seed", type=int, default=42, help="Random seed.")
    parser.add_argument(
        "--n-cells", type=int, default=100, help="Initial number of healthy cells."
    )
    parser.add_argument(
        "--energy-budget", type=float, default=10000.0, help="Initial global energy budget."
    )
    parser.add_argument(
        "--pre-steps",
        type=int,
        default=200,
        help="Pre-therapy simulation steps.",
    )
    parser.add_argument(
        "--post-steps",
        type=int,
        default=100,
        help="Post-therapy simulation steps.",
    )
    args = parser.parse_args()

    print("=" * 60)
    print("  Cancer Boundary Logic: Full Scenario Comparison")
    print("=" * 60)

    run_kwargs = CALIBRATED_MONTHLY_PROFILE.run_kwargs()
    print(
        "Parameter profile: "
        f"{CALIBRATED_MONTHLY_PROFILE.name} "
        f"({CALIBRATED_MONTHLY_PROFILE.step_unit})"
    )
    print(
        "Run parameters: "
        f"p_mutation={run_kwargs['p_mutation']}, "
        f"p_progression={run_kwargs['p_progression']:.4f}, "
        f"p_healthy_division={run_kwargs['p_healthy_division']}"
    )
    print(
        "Run config: "
        f"seed={args.seed}, n_cells={args.n_cells}, "
        f"energy_budget={args.energy_budget}, "
        f"pre_steps={args.pre_steps}, post_steps={args.post_steps}"
    )

    therapies = [
        ("No Therapy (Control)", None),
        ("Forced Distinction", forced_distinction),
        ("Differentiation Therapy (ATRA)", differentiation_therapy),
        ("Checkpoint Inhibitor (anti-PD-1)", checkpoint_inhibitor),
        ("Bioelectric Reprogramming (Levin)", bioelectric_reprogramming),
        ("Combination (Anti-PD-1 + ATRA)", combination_checkpoint_differentiation),
    ]

    all_results = []

    for name, fn in therapies:
        print(f"\n--- {name} ---")
        result = run_scenario(
            name,
            fn,
            seed=args.seed,
            run_kwargs=run_kwargs,
            n_cells=args.n_cells,
            energy_budget=args.energy_budget,
            pre_steps=args.pre_steps,
            post_steps=args.post_steps,
        )
        pre = result["pre_therapy"]
        post = result["post_therapy"]
        print(f"  Pre-therapy:  kappa={pre['kappa']:.3f}, cancer={pre['cancer']}")
        print(f"  Treated:      {result['treated_cells']} cells")
        print(f"  Post-therapy: kappa={post['kappa']:.3f}, cancer={post['cancer']}")
        all_results.append(result)

    # Save numerical results
    results_file = os.path.join(RESULTS_DIR, "scenario_comparison.json")
    serializable = []
    for r in all_results:
        serializable.append({
            "therapy": r["therapy"],
            "treated_cells": r["treated_cells"],
            "pre_therapy_kappa": r["pre_therapy"]["kappa"],
            "pre_therapy_cancer": r["pre_therapy"]["cancer"],
            "post_therapy_kappa": r["post_therapy"]["kappa"],
            "post_therapy_cancer": r["post_therapy"]["cancer"],
        })
    with open(results_file, "w") as f:
        json.dump(serializable, f, indent=2)
    print(f"\nResults saved to {results_file}")

    metadata_file = os.path.join(RESULTS_DIR, "scenario_comparison_metadata.json")
    with open(metadata_file, "w") as f:
        json.dump(
            {
                "parameter_profile": {
                    "name": CALIBRATED_MONTHLY_PROFILE.name,
                    "step_unit": CALIBRATED_MONTHLY_PROFILE.step_unit,
                    "seed": args.seed,
                    "n_cells": args.n_cells,
                    "energy_budget": args.energy_budget,
                    "pre_steps": args.pre_steps,
                    "post_steps": args.post_steps,
                    **run_kwargs,
                    "source_refs": CALIBRATED_MONTHLY_PROFILE.source_refs,
                    "note": CALIBRATED_MONTHLY_PROFILE.note,
                }
            },
            f,
            indent=2,
        )
    print(f"Metadata saved to {metadata_file}")

    # Generate comparison plot
    try:
        import matplotlib
        matplotlib.use("Agg")  # headless
        import matplotlib.pyplot as plt

        fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(12, 10), sharex=True)

        for r in all_results:
            ax1.plot(r["kappa_series"], label=r["therapy"], linewidth=1.5)
            ax2.plot(r["cancer_series"], label=r["therapy"], linewidth=1.5)
            ax3.plot(r["pre_cancer_series"], label=r["therapy"], linewidth=1.5)

        ax1.axvline(x=args.pre_steps, color="black", linestyle="--", alpha=0.5, label="Therapy applied")
        ax2.axvline(x=args.pre_steps, color="black", linestyle="--", alpha=0.5)
        ax3.axvline(x=args.pre_steps, color="black", linestyle="--", alpha=0.5)

        ax1.set_ylabel("Coherence (kappa)")
        ax1.set_ylim(-0.05, 1.05)
        ax1.legend(loc="lower left", fontsize=8)
        ax1.set_title("Therapy Comparison: Coherence Recovery")
        ax1.grid(True, alpha=0.3)

        ax2.set_ylabel("Cancer Cell Count")
        ax2.legend(loc="upper left", fontsize=8)
        ax2.grid(True, alpha=0.3)

        ax3.set_ylabel("Pre-Cancer Cell Count")
        ax3.set_xlabel("Time Steps")
        ax3.legend(loc="upper left", fontsize=8)
        ax3.grid(True, alpha=0.3)

        plt.tight_layout()
        plot_file = os.path.join(RESULTS_DIR, "therapy_comparison.png")
        plt.savefig(plot_file, dpi=150, bbox_inches="tight")
        print(f"Plot saved to {plot_file}")
        plt.close()

        # Individual kappa trajectory
        fig, ax = plt.subplots(figsize=(10, 5))
        best = next(
            (r for r in all_results if r["therapy"] == "Differentiation Therapy (ATRA)"),
            all_results[0],
        )
        ax.plot(best["kappa_series"], "b-", linewidth=2, label="Coherence (kappa)")
        ax.axhline(y=0.5, color="orange", linestyle="--", alpha=0.5, label="Pre-cancer threshold")
        ax.axhline(y=0.3, color="red", linestyle="--", alpha=0.5, label="Malignancy threshold")
        ax.axvline(x=args.pre_steps, color="black", linestyle="--", alpha=0.5, label="Therapy applied")
        ax.set_xlabel("Time Steps")
        ax.set_ylabel("Tissue Coherence (kappa)")
        ax.set_title(f"Coherence Trajectory: {best['therapy']}")
        ax.set_ylim(-0.05, 1.05)
        ax.legend()
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        kappa_file = os.path.join(RESULTS_DIR, "kappa_trajectory.png")
        plt.savefig(kappa_file, dpi=150, bbox_inches="tight")
        print(f"Plot saved to {kappa_file}")
        plt.close()

    except ImportError:
        print("Install matplotlib for plots: pip install matplotlib")

    # Print summary table
    print("\n" + "=" * 70)
    print(f"{'Therapy':<35} {'Pre-kappa':>10} {'Post-kappa':>11} {'Cancer':>8}")
    print("-" * 70)
    for r in serializable:
        print(
            f"{r['therapy']:<35} "
            f"{r['pre_therapy_kappa']:>10.3f} "
            f"{r['post_therapy_kappa']:>11.3f} "
            f"{r['post_therapy_cancer']:>8}"
        )
    print("=" * 70)


if __name__ == "__main__":
    main()
