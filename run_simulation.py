"""
Cancer as Boundary Logic Failure -- Full Simulation Demo
=========================================================

Moses Rahnama, 2026
Companion code for: Cancer_As_Boundary_Logic_Failure.tex

This script demonstrates:
  1. Healthy tissue initialization with coherence ~1.0
  2. Stochastic cancer emergence via mutation + progression
  3. Coherence collapse as cancer cells decouple (rec_succ without halt)
  4. Energy budget depletion (cachexia / Warburg effect proxy)
  5. Therapeutic intervention via forced distinction
  6. Post-therapy coherence recovery

Run:  python run_simulation.py
Deps: pip install matplotlib networkx numpy
"""

import sys
import os
import argparse

# Add src to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "src"))
RESULTS_DIR = os.path.join(os.path.dirname(__file__), "results")
os.makedirs(RESULTS_DIR, exist_ok=True)

from simulations.engine import Simulation
from simulations.parameter_profiles import CALIBRATED_MONTHLY_PROFILE
from simulations.therapy import (
    apply_therapy_to_all_cancer,
    forced_distinction,
)

def main():
    parser = argparse.ArgumentParser(
        description="Run single-scenario cancer boundary logic simulation."
    )
    parser.add_argument("--seed", type=int, default=42, help="Random seed.")
    parser.add_argument(
        "--n-cells", type=int, default=100, help="Initial number of healthy cells."
    )
    parser.add_argument(
        "--energy-budget", type=float, default=10000.0, help="Initial global energy budget."
    )
    parser.add_argument(
        "--pre-steps", type=int, default=200, help="Pre-therapy simulation steps."
    )
    parser.add_argument(
        "--post-steps", type=int, default=100, help="Post-therapy simulation steps."
    )
    parser.add_argument(
        "--skip-plots", action="store_true", help="Skip plotting output files."
    )
    args = parser.parse_args()

    print("=" * 60)
    print("  Cancer as Boundary Logic Failure")
    print("  Simulation Framework")
    print("  Moses Rahnama, 2026")
    print("=" * 60)

    run_kwargs = CALIBRATED_MONTHLY_PROFILE.run_kwargs()
    print("\nParameter profile:")
    print(f"  name={CALIBRATED_MONTHLY_PROFILE.name}")
    print(f"  step_unit={CALIBRATED_MONTHLY_PROFILE.step_unit}")
    print(
        "  "
        f"p_mutation={run_kwargs['p_mutation']}, "
        f"p_progression={run_kwargs['p_progression']:.4f}, "
        f"p_healthy_division={run_kwargs['p_healthy_division']}"
    )
    print(
        "  "
        f"seed={args.seed}, n_cells={args.n_cells}, energy_budget={args.energy_budget}, "
        f"pre_steps={args.pre_steps}, post_steps={args.post_steps}"
    )

    # Phase 1: Cancer progression
    print(f"\n--- Phase 1: Cancer Progression ({args.pre_steps} steps) ---")
    sim = Simulation(n_cells=args.n_cells, energy_budget=args.energy_budget, seed=args.seed)
    sim.run(steps=args.pre_steps, **run_kwargs)
    print(sim.summary())

    # Phase 2: Apply therapy
    cancer_count = sim.history[-1]["cancer"]
    if cancer_count > 0:
        print(f"\n--- Phase 2: Applying Forced Distinction Therapy ---")
        treated = apply_therapy_to_all_cancer(
            sim.tissue, forced_distinction
        )
        print(f"  Treated {treated} cancer cells")

        # Phase 3: Post-therapy evolution
        print(f"\n--- Phase 3: Post-Therapy Evolution ({args.post_steps} steps) ---")
        sim.run(steps=args.post_steps, **run_kwargs)
        print(sim.summary())
    else:
        print("\nNo cancer emerged in this run. Try a different seed.")

    # Visualize
    if args.skip_plots:
        print("\n--- Plot generation skipped (--skip-plots) ---")
    else:
        try:
            from visualization.plots import (
                plot_coherence_vs_cellcount,
                plot_state_distribution,
                plot_tissue_network,
            )
            print("\n--- Generating Plots ---")
            plot_coherence_vs_cellcount(
                sim,
                save_path=os.path.join(RESULTS_DIR, "coherence_vs_cancer.png"),
            )
            plot_state_distribution(
                sim,
                save_path=os.path.join(RESULTS_DIR, "state_distribution.png"),
            )
            plot_tissue_network(
                sim.tissue,
                save_path=os.path.join(RESULTS_DIR, "tissue_network.png"),
            )
        except ImportError:
            print("\nInstall matplotlib and networkx for visualizations:")
            print("  pip install matplotlib networkx")
        except Exception as e:
            print(f"\nVisualization skipped (headless mode): {e}")

    print("\nDone.")

if __name__ == "__main__":
    main()
