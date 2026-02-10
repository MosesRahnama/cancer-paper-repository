"""
Therapeutic intervention operators for the Boundary Logic Failure model.

Each function acts on a Tissue, modifying cancer cells to restore
boundary observability. These implement the "forced distinction"
principle from Section 12 of the paper.

Reference: Cancer_As_Boundary_Logic_Failure.tex, Section 12
"""

from core.cell import Cell, CellState
from core.tissue import Tissue, RelationType


def forced_distinction(tissue: Tissue, target_id: str):
    """
    Core therapeutic principle: restore informational distinguishability.

    Increases organism alignment so the immune system (the body's
    Born Rule) can distinguish the cancer cell from self.
    """
    if target_id not in tissue.cells:
        return
    cell = tissue.cells[target_id]
    if cell.state != CellState.CANCER:
        return

    cell.p_align = 0.4  # Partially restored -- now distinguishable
    cell.perspective = {"self": 0.5, "organism": 0.5}
    # Re-add temporal coupling (clock sync restored)
    tissue.edges.add(("clock", target_id, RelationType.TEMPORAL))


def differentiation_therapy(tissue: Tissue, target_id: str):
    """
    ATRA-like differentiation therapy.

    Forces cancer cells to complete the arrested differentiation
    program. Does NOT kill -- redirects into organism's boundary.
    >90% complete remission in APL.
    """
    if target_id not in tissue.cells:
        return
    cell = tissue.cells[target_id]
    if cell.state != CellState.CANCER:
        return

    cell.state = CellState.HEALTHY
    cell.p_align = 0.85
    cell.division_rate = 0.0  # Growth arrest
    cell.energy_cost = 1.0    # Normal metabolism
    cell.perspective = {"self": 0.15, "organism": 0.85}

    # Restore signaling edges to nearby healthy cells
    for other_id, other in tissue.cells.items():
        if other_id != target_id and other.state == CellState.HEALTHY:
            tissue.edges.add((target_id, other_id, RelationType.REGULATORY))
            tissue.edges.add((other_id, target_id, RelationType.REGULATORY))
            break  # Connect to at least one neighbor

    tissue.edges.add(("clock", target_id, RelationType.TEMPORAL))


def checkpoint_inhibitor(tissue: Tissue, target_id: str):
    """
    Anti-PD-1/PD-L1 checkpoint inhibitor.

    Removes the false "self" signal so the immune system can
    detect the cancer cell. In boundary logic: restoring
    distinguishability without changing cell state.
    """
    if target_id not in tissue.cells:
        return
    cell = tissue.cells[target_id]
    if cell.state != CellState.CANCER:
        return

    # Make cancer cell partially visible (immune can now "see" it)
    cell.perspective["organism"] = 0.3
    cell.p_align = 0.3


def bioelectric_reprogramming(tissue: Tissue, target_id: str):
    """
    Ion channel / gap junction manipulation.

    Restores collective voltage patterns that suppress tumorigenesis.
    Based on Levin's work: hyperpolarized cells resist tumor
    formation despite active oncogene expression.
    """
    if target_id not in tissue.cells:
        return
    cell = tissue.cells[target_id]
    if cell.state != CellState.CANCER:
        return

    cell.p_align = 0.6
    cell.perspective = {"self": 0.35, "organism": 0.65}
    cell.division_rate = 0.05  # Reduced but not halted

    # Restore adhesion edges (gap junction reconnection)
    for other_id, other in tissue.cells.items():
        if other_id != target_id and other.state == CellState.HEALTHY:
            tissue.edges.add((target_id, other_id, RelationType.ADHESION))
            break

    tissue.edges.add(("clock", target_id, RelationType.TEMPORAL))


def combination_checkpoint_differentiation(tissue: Tissue, target_id: str):
    """
    Combination therapy: checkpoint inhibitor (anti-PD-1) + differentiation (ATRA).

    Two-step mechanism:
      1. Restore observability (checkpoint inhibitor) -- so the cell can
         "hear" organism-level stop signals.
      2. Restore identity (differentiation) -- redirect cell back into
         the organism's boundary program.

    Rationale: checkpoint inhibition alone creates a stalemate (the immune
    system can see the cell, but the cell keeps proliferating).
    Differentiation alone may fail if cells can't receive the signal.
    The combination is predicted to be synergistic -- restoring the
    communication channel first, then sending the corrective signal.
    """
    if target_id not in tissue.cells:
        return
    cell = tissue.cells[target_id]
    if cell.state != CellState.CANCER:
        return

    # Step 1: Restore observability (forced distinction / checkpoint inhibitor)
    checkpoint_inhibitor(tissue, target_id)

    # Step 2: Restore identity (differentiation therapy)
    differentiation_therapy(tissue, target_id)


def apply_therapy_to_all_cancer(tissue: Tissue, therapy_fn):
    """Apply a therapy function to every cancer cell in the tissue."""
    cancer_ids = [
        cid for cid, cell in tissue.cells.items()
        if cell.state == CellState.CANCER
    ]
    for cid in cancer_ids:
        therapy_fn(tissue, cid)
    return len(cancer_ids)
