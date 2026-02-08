"""
Simulation engine for the Boundary Logic Failure cancer model.

Wraps Tissue with history recording, parameter presets, and
therapy injection points.

Reference: Cancer_As_Boundary_Logic_Failure.tex, Section 11
"""

from typing import Dict, List, Optional
from core.cell import CellState
from core.tissue import Tissue
from core.metrics import compute_coherence, compute_organism_alignment


class Simulation:
    """Run a tissue simulation and record metrics at each step."""

    def __init__(
        self,
        n_cells: int = 100,
        energy_budget: float = 10000.0,
        seed: Optional[int] = None,
    ):
        import random
        if seed is not None:
            random.seed(seed)

        self.tissue = Tissue(n_cells=n_cells, energy_budget=energy_budget)
        self.step_count = 0

        # History
        self.history: List[Dict] = []

    def step(self, **kwargs):
        """Advance one timestep and record metrics."""
        self.tissue.step(**kwargs)
        self.step_count += 1
        self._record()

    def run(self, steps: int = 200, **kwargs):
        """Run for multiple steps."""
        for _ in range(steps):
            self.step(**kwargs)

    def _record(self):
        counts = self.tissue.count_by_state()
        self.history.append({
            "time": self.step_count,
            "kappa": compute_coherence(self.tissue),
            "alignment": compute_organism_alignment(self.tissue),
            "total_cells": self.tissue.total_cells(),
            "healthy": counts.get(CellState.HEALTHY, 0),
            "pre_cancer": counts.get(CellState.PRE_CANCER, 0),
            "cancer": counts.get(CellState.CANCER, 0),
            "energy": self.tissue.energy_budget,
        })

    # --- Convenience accessors ---

    def kappa_series(self) -> List[float]:
        return [h["kappa"] for h in self.history]

    def cancer_series(self) -> List[int]:
        return [h["cancer"] for h in self.history]

    def energy_series(self) -> List[float]:
        return [h["energy"] for h in self.history]

    def summary(self) -> str:
        if not self.history:
            return "No simulation data."
        h = self.history[-1]
        kappa = h["kappa"]

        if kappa >= 0.9:
            status = "Healthy tissue"
        elif kappa >= 0.7:
            status = "Minor perturbation"
        elif kappa >= 0.5:
            status = "Pre-cancerous"
        elif kappa >= 0.3:
            status = "Early malignancy"
        elif kappa >= 0.1:
            status = "Advanced malignancy"
        else:
            status = "Terminal / Metastatic"

        lines = [
            f"=== Simulation t={h['time']} ===",
            f"  Cells:   {h['total_cells']} (healthy={h['healthy']}, "
            f"pre-cancer={h['pre_cancer']}, cancer={h['cancer']})",
            f"  Kappa:   {kappa:.3f}  ({status})",
            f"  Energy:  {h['energy']:.1f}",
        ]
        return "\n".join(lines)
