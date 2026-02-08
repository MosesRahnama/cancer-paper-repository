from core.cell import CellState
from core.tissue import Tissue
from simulations.therapy import (
    apply_therapy_to_all_cancer,
    differentiation_therapy,
    forced_distinction,
)


def _make_cancer_tissue(n_cells=8, n_cancer=3):
    tissue = Tissue(n_cells=n_cells, energy_budget=1000.0)
    for idx in range(n_cancer):
        tissue.cells[f"cell_{idx}"].transform_to_cancer()
    return tissue


def test_apply_therapy_counts_only_cancer_cells():
    tissue = _make_cancer_tissue(n_cells=8, n_cancer=3)
    treated = apply_therapy_to_all_cancer(tissue, forced_distinction)
    assert treated == 3


def test_differentiation_therapy_restores_healthy_state():
    tissue = _make_cancer_tissue(n_cells=5, n_cancer=1)
    treated = apply_therapy_to_all_cancer(tissue, differentiation_therapy)

    assert treated == 1
    assert tissue.cells["cell_0"].state == CellState.HEALTHY
    assert tissue.cells["cell_0"].division_rate <= 0.01
