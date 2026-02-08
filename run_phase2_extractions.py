"""
Phase 2 extracted operators and seed-sweep runner.

This script ports high-value extracted ideas into the organized companion
codebase without replacing the modular engine:
1) Harmonic-aware coherence variant (kappa_har)
2) Regulatory unmasking operator (pre-cancer stressor)
3) Force re-entanglement therapy operator
4) Multi-seed therapy sweep with aggregate JSON + plots

Usage examples:
  python run_phase2_extractions.py
  python run_phase2_extractions.py --seeds 40 --apply-unmasking --unmask-severity 0.5
"""

import argparse
import json
import os
import random
import sys
from collections import defaultdict
from statistics import fmean, pstdev

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


HERE = os.path.dirname(__file__)
SRC_DIR = os.path.join(HERE, "src")
RESULTS_DIR = os.path.join(HERE, "results")
os.makedirs(RESULTS_DIR, exist_ok=True)
if SRC_DIR not in sys.path:
    sys.path.insert(0, SRC_DIR)

from core.cell import CellState
from core.metrics import compute_identity_entropy
from core.tissue import RelationType
from simulations.engine import Simulation
from simulations.parameter_profiles import CALIBRATED_MONTHLY_PROFILE
from simulations.therapy import (
    apply_therapy_to_all_cancer,
    bioelectric_reprogramming,
    checkpoint_inhibitor,
    differentiation_therapy,
    forced_distinction,
)


def compute_harmonic_kappa(tissue, w_sig=0.4, w_har=0.2, w_identity=0.4):
    """
    Harmonic-aware coherence variant extracted from support material.

    kappa_har = clip( w_sig * rho_sig + w_har * rho_har + w_identity * (1 - H) )
    where:
      rho_sig = fraction of regulatory edges among all edges
      rho_har = adhesion-edge density per node (clipped to [0,1])
      H       = mean identity entropy
    """
    total_edges = len(tissue.edges)
    if total_edges == 0 or not tissue.cells:
        return 0.0

    sig_count = sum(1 for _, _, r in tissue.edges if r == RelationType.REGULATORY)
    har_count = sum(1 for _, _, r in tissue.edges if r == RelationType.ADHESION)

    rho_sig = sig_count / total_edges
    rho_har = min(1.0, har_count / max(1, len(tissue.cells)))
    h = compute_identity_entropy(tissue)

    kappa = w_sig * rho_sig + w_har * rho_har + w_identity * (1.0 - h)
    return max(0.0, min(1.0, kappa))


def apply_unmasking(tissue, severity, rng, max_fraction=0.12):
    """
    Experimental stress operator inspired by extracted "constraint removal".

    Applies a pre-cancerous tilt to a small fraction of healthy cells and prunes
    one regulatory edge per targeted cell.
    """
    healthy_ids = [cid for cid, c in tissue.cells.items() if c.state == CellState.HEALTHY]
    if not healthy_ids:
        return 0

    n_targets = max(1, int(len(healthy_ids) * max_fraction * max(0.0, min(1.0, severity))))
    targets = rng.sample(healthy_ids, min(n_targets, len(healthy_ids)))

    for cid in targets:
        cell = tissue.cells[cid]
        cell.state = CellState.PRE_CANCER
        cell.p_align = max(0.2, cell.p_align - 0.35 * severity)
        cell.perspective = {
            "self": min(0.95, 0.4 + 0.3 * severity),
            "organism": max(0.05, 0.6 - 0.3 * severity),
        }
        cell.division_rate = min(0.2, cell.division_rate * (1.0 + 0.75 * severity))

        for edge in list(tissue.edges):
            s, t, rel = edge
            if cid in (s, t) and rel == RelationType.REGULATORY:
                tissue.edges.discard(edge)
                break

    return len(targets)


def force_reentanglement(tissue, target_id):
    """
    Experimental therapy operator inspired by extracted re-entanglement notes.

    Unlike differentiation therapy, this does not force healthy state; it
    restores observability and coupling pressure.
    """
    if target_id not in tissue.cells:
        return

    cell = tissue.cells[target_id]
    if cell.state != CellState.CANCER:
        return

    cell.p_align = max(cell.p_align, 0.45)
    cell.perspective = {"self": 0.55, "organism": 0.45}
    cell.division_rate = max(0.03, cell.division_rate * 0.35)

    tissue.edges.add(("clock", target_id, RelationType.TEMPORAL))

    for other_id, other in tissue.cells.items():
        if other_id == target_id or other.state != CellState.HEALTHY:
            continue
        tissue.edges.add((target_id, other_id, RelationType.REGULATORY))
        tissue.edges.add((other_id, target_id, RelationType.REGULATORY))
        tissue.edges.add((target_id, other_id, RelationType.ADHESION))
        break


def run_trial(seed, therapy_name, therapy_fn, n_cells, energy_budget, apply_unmasking_flag, unmask_severity):
    sim = Simulation(n_cells=n_cells, energy_budget=energy_budget, seed=seed)
    rng = random.Random(seed + 10007)
    run_kwargs = CALIBRATED_MONTHLY_PROFILE.run_kwargs()

    unmasked_cells = 0
    if apply_unmasking_flag:
        unmasked_cells = apply_unmasking(sim.tissue, severity=unmask_severity, rng=rng)

    sim.run(steps=200, **run_kwargs)
    pre = dict(sim.history[-1])
    pre_harmonic = compute_harmonic_kappa(sim.tissue)

    if therapy_fn is None:
        treated = 0
    else:
        treated = apply_therapy_to_all_cancer(sim.tissue, therapy_fn)

    sim.run(steps=100, **run_kwargs)
    post = dict(sim.history[-1])
    post_harmonic = compute_harmonic_kappa(sim.tissue)

    return {
        "seed": seed,
        "therapy": therapy_name,
        "treated_cells": treated,
        "unmasked_cells": unmasked_cells,
        "pre_kappa": pre["kappa"],
        "post_kappa": post["kappa"],
        "pre_kappa_harmonic": pre_harmonic,
        "post_kappa_harmonic": post_harmonic,
        "pre_cancer": pre["cancer"],
        "post_cancer": post["cancer"],
    }


def summarize(records):
    grouped = defaultdict(list)
    for r in records:
        grouped[r["therapy"]].append(r)

    summary = {}
    for therapy, rows in grouped.items():
        summary[therapy] = {
            "n": len(rows),
            "post_kappa_mean": fmean(r["post_kappa"] for r in rows),
            "post_kappa_sd": pstdev(r["post_kappa"] for r in rows) if len(rows) > 1 else 0.0,
            "post_kappa_harmonic_mean": fmean(r["post_kappa_harmonic"] for r in rows),
            "post_kappa_harmonic_sd": (
                pstdev(r["post_kappa_harmonic"] for r in rows) if len(rows) > 1 else 0.0
            ),
            "post_cancer_mean": fmean(r["post_cancer"] for r in rows),
            "post_cancer_sd": pstdev(r["post_cancer"] for r in rows) if len(rows) > 1 else 0.0,
            "treated_cells_mean": fmean(r["treated_cells"] for r in rows),
        }
    return summary


def plot_summary(summary, results_dir, tag):
    therapy_order = [
        "No Therapy (Control)",
        "Forced Distinction",
        "Force Re-Entanglement (Extracted)",
        "Differentiation Therapy (ATRA)",
        "Checkpoint Inhibitor (anti-PD-1)",
        "Bioelectric Reprogramming (Levin)",
    ]
    labels = [t for t in therapy_order if t in summary]
    if not labels:
        return

    x = list(range(len(labels)))

    base_mean = [summary[t]["post_kappa_mean"] for t in labels]
    base_sd = [summary[t]["post_kappa_sd"] for t in labels]
    har_mean = [summary[t]["post_kappa_harmonic_mean"] for t in labels]
    har_sd = [summary[t]["post_kappa_harmonic_sd"] for t in labels]

    fig, ax = plt.subplots(figsize=(13, 6))
    width = 0.38
    ax.bar([i - width / 2 for i in x], base_mean, width=width, yerr=base_sd, capsize=4, label="kappa")
    ax.bar(
        [i + width / 2 for i in x],
        har_mean,
        width=width,
        yerr=har_sd,
        capsize=4,
        label="kappa_harmonic",
    )
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=18, ha="right")
    ax.set_ylabel("Post-therapy coherence")
    ax.set_title(f"Phase 2 Seed Sweep ({tag}): base vs harmonic coherence")
    ax.grid(True, axis="y", alpha=0.3)
    ax.legend()
    plt.tight_layout()
    out1 = os.path.join(results_dir, f"phase2_{tag}_kappa_summary.png")
    plt.savefig(out1, dpi=150, bbox_inches="tight")
    plt.close(fig)

    cancer_mean = [summary[t]["post_cancer_mean"] for t in labels]
    cancer_sd = [summary[t]["post_cancer_sd"] for t in labels]
    fig, ax = plt.subplots(figsize=(13, 5))
    ax.bar(x, cancer_mean, yerr=cancer_sd, capsize=4, color="tomato")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=18, ha="right")
    ax.set_ylabel("Post-therapy cancer cell count")
    ax.set_title(f"Phase 2 Seed Sweep ({tag}): residual cancer burden")
    ax.grid(True, axis="y", alpha=0.3)
    plt.tight_layout()
    out2 = os.path.join(results_dir, f"phase2_{tag}_cancer_summary.png")
    plt.savefig(out2, dpi=150, bbox_inches="tight")
    plt.close(fig)

    return out1, out2


def run_sweep(seeds, n_cells, energy_budget, apply_unmasking_flag, unmask_severity):
    therapies = [
        ("No Therapy (Control)", None),
        ("Forced Distinction", forced_distinction),
        ("Force Re-Entanglement (Extracted)", force_reentanglement),
        ("Differentiation Therapy (ATRA)", differentiation_therapy),
        ("Checkpoint Inhibitor (anti-PD-1)", checkpoint_inhibitor),
        ("Bioelectric Reprogramming (Levin)", bioelectric_reprogramming),
    ]

    records = []
    for seed in range(seeds):
        for therapy_name, therapy_fn in therapies:
            records.append(
                run_trial(
                    seed=seed,
                    therapy_name=therapy_name,
                    therapy_fn=therapy_fn,
                    n_cells=n_cells,
                    energy_budget=energy_budget,
                    apply_unmasking_flag=apply_unmasking_flag,
                    unmask_severity=unmask_severity,
                )
            )
    return records


def main():
    parser = argparse.ArgumentParser(description="Run phase-2 extracted operators and seed sweep.")
    parser.add_argument("--seeds", type=int, default=20, help="Number of random seeds to run.")
    parser.add_argument("--n-cells", type=int, default=100, help="Initial healthy cells.")
    parser.add_argument("--energy-budget", type=float, default=10000.0, help="Initial energy budget.")
    parser.add_argument(
        "--apply-unmasking",
        action="store_true",
        help="Apply extracted unmasking stressor before progression.",
    )
    parser.add_argument(
        "--unmask-severity",
        type=float,
        default=0.4,
        help="Unmasking severity in [0,1], used only with --apply-unmasking.",
    )
    args = parser.parse_args()
    run_kwargs = CALIBRATED_MONTHLY_PROFILE.run_kwargs()

    records = run_sweep(
        seeds=args.seeds,
        n_cells=args.n_cells,
        energy_budget=args.energy_budget,
        apply_unmasking_flag=args.apply_unmasking,
        unmask_severity=args.unmask_severity,
    )
    summary = summarize(records)
    tag = "unmasking" if args.apply_unmasking else "baseline"
    plot_paths = plot_summary(summary, RESULTS_DIR, tag)

    payload = {
        "config": {
            "seeds": args.seeds,
            "n_cells": args.n_cells,
            "energy_budget": args.energy_budget,
            "apply_unmasking": args.apply_unmasking,
            "unmask_severity": args.unmask_severity,
            "parameter_profile_name": CALIBRATED_MONTHLY_PROFILE.name,
            "step_unit": CALIBRATED_MONTHLY_PROFILE.step_unit,
            **run_kwargs,
            "source_refs": CALIBRATED_MONTHLY_PROFILE.source_refs,
        },
        "summary": summary,
        "records": records,
    }
    json_path = os.path.join(RESULTS_DIR, f"phase2_{tag}_seed_sweep.json")
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)

    print("=" * 72)
    print(f"Phase 2 sweep complete ({tag})")
    print(f"Records: {len(records)}")
    print(f"JSON:    {json_path}")
    if plot_paths:
        print(f"Plots:   {plot_paths[0]}")
        print(f"         {plot_paths[1]}")
    print("=" * 72)


if __name__ == "__main__":
    main()
