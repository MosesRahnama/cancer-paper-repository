import random

from core.cell import CellState
from core.tissue import RelationType, Tissue
from simulations.engine import Simulation


def test_precancer_progression_severs_regulatory_and_temporal_edges():
    tissue = Tissue(n_cells=6, energy_budget=1000.0)
    target = "cell_0"
    tissue.cells[target].transform_to_precancer()

    random.seed(123)
    tissue.step(p_mutation=0.0, p_progression=1.0, p_healthy_division=0.0)

    assert tissue.cells[target].state == CellState.CANCER
    assert all(
        not (target in (s, t) and r in (RelationType.REGULATORY, RelationType.TEMPORAL))
        for s, t, r in tissue.edges
    )


def test_model_can_accumulate_precancer_population_with_slow_progression():
    sim = Simulation(n_cells=80, energy_budget=8000.0, seed=11)
    sim.run(
        steps=60,
        p_mutation=0.02,
        p_progression=0.002,
        p_healthy_division=0.01,
    )

    peak_precancer = max(h["pre_cancer"] for h in sim.history)
    assert peak_precancer > 0
