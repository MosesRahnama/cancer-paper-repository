"""
Tissue network simulation for the Boundary Logic Failure cancer model.

A tissue is modeled as a labeled graph T = (V, E), where each node is a
Cell and each edge carries a label indicating the type of coupling.

Edge types:
  - REGULATORY (SIG): growth inhibition, differentiation cues, contact inhibition
  - LINEAGE (DIV): division history (parent-child)
  - TEMPORAL (TMP): shared timing fields (circadian, hormonal)
  - ADHESION (HAR): tissue adhesion, gap junction communication

Reference: Cancer_As_Boundary_Logic_Failure.tex, Section 9
"""

import random
from enum import Enum
from typing import Dict, List, Set, Tuple, Optional

from .cell import Cell, CellState


class RelationType(Enum):
    REGULATORY = "SIG"   # Signaling / contact inhibition
    LINEAGE = "DIV"      # Division history (parent-child)
    TEMPORAL = "TMP"     # Circadian / hormonal sync
    ADHESION = "HAR"     # Tissue adhesion / gap junctions


class Tissue:
    """A tissue region modeled as a labeled directed graph.

    Attributes:
        cells: Dict mapping cell ID to Cell object.
        edges: Set of (source_id, target_id, RelationType) tuples.
        energy_budget: Global energy pool. Depletes as cells consume.
        time: Current simulation timestep.
    """

    def __init__(self, n_cells: int = 100, energy_budget: float = 10000.0):
        self.cells: Dict[str, Cell] = {}
        self.edges: Set[Tuple[str, str, RelationType]] = set()
        self.energy_budget = energy_budget
        self.time = 0
        self._next_id = 0

        # Initialize healthy tissue
        for i in range(n_cells):
            cid = f"cell_{i}"
            self.cells[cid] = Cell(id=cid)
            self._next_id = i + 1

            # Chain regulatory edges (nearest-neighbor signaling)
            if i > 0:
                prev = f"cell_{i - 1}"
                self.edges.add((prev, cid, RelationType.REGULATORY))
                self.edges.add((cid, prev, RelationType.REGULATORY))

            # Temporal coupling (all cells synced to organism clock)
            self.edges.add(("clock", cid, RelationType.TEMPORAL))

        # Add adhesion edges (every 3rd neighbor for tissue structure)
        for i in range(n_cells):
            for j in [i + 2, i + 3]:
                if j < n_cells:
                    a, b = f"cell_{i}", f"cell_{j}"
                    self.edges.add((a, b, RelationType.ADHESION))

    def get_neighbors(self, cell_id: str, edge_type: Optional[RelationType] = None) -> List[str]:
        """Get neighboring cell IDs, optionally filtered by edge type."""
        neighbors = []
        for s, t, r in self.edges:
            if edge_type and r != edge_type:
                continue
            if s == cell_id and t in self.cells:
                neighbors.append(t)
            elif t == cell_id and s in self.cells:
                neighbors.append(s)
        return neighbors

    def regulatory_degree(self, cell_id: str) -> int:
        """Count regulatory (signaling) edges for a cell."""
        return len(self.get_neighbors(cell_id, RelationType.REGULATORY))

    def step(
        self,
        p_mutation: float = 0.001,
        p_progression: float = 0.05,
        p_healthy_division: float = 0.01,
        max_healthy_neighbors: int = 4,
    ):
        """Advance the simulation by one timestep.

        1. Each cell consumes energy from global budget.
        2. Random mutation: HEALTHY -> PRE_CANCER (small probability).
        3. Progression: PRE_CANCER -> CANCER (boundary collapse).
        4. Division: healthy cells respect contact inhibition;
           cancer cells do not.
        """
        self.time += 1
        new_cells: Dict[str, Cell] = {}

        for cid, cell in list(self.cells.items()):
            # Energy consumption
            if self.energy_budget < cell.energy_cost:
                continue
            self.energy_budget -= cell.energy_cost

            initial_state = cell.state

            # --- Mutation (HEALTHY -> PRE_CANCER) ---
            if initial_state == CellState.HEALTHY and random.random() < p_mutation:
                cell.transform_to_precancer()

            # --- Progression (PRE_CANCER -> CANCER) ---
            elif initial_state == CellState.PRE_CANCER and random.random() < p_progression:
                cell.transform_to_cancer()
                # Sever regulatory and temporal edges (boundary collapse)
                self.edges = {
                    (s, t, r)
                    for s, t, r in self.edges
                    if not (
                        cid in (s, t)
                        and r in (RelationType.REGULATORY, RelationType.TEMPORAL)
                    )
                }

            # --- Division ---
            # Healthy-cell division can be swept independently via p_healthy_division.
            division_prob = (
                p_healthy_division
                if cell.state == CellState.HEALTHY
                else cell.division_rate
            )
            if random.random() < division_prob:
                if cell.state == CellState.HEALTHY:
                    # Contact inhibition: only divide if below neighbor threshold
                    if self.regulatory_degree(cid) < max_healthy_neighbors:
                        daughter = self._divide(cid, cell)
                        if daughter:
                            new_cells[daughter.id] = daughter
                else:
                    # Cancer cells: no contact inhibition (rec_succ without halt)
                    daughter = self._divide(cid, cell)
                    if daughter:
                        new_cells[daughter.id] = daughter

        self.cells.update(new_cells)

    def _divide(self, parent_id: str, parent: Cell) -> Optional[Cell]:
        """Execute cell division. Returns the daughter cell."""
        daughter_id = f"cell_{self._next_id}"
        self._next_id += 1
        daughter = parent.create_daughter(daughter_id)
        self.cells[daughter_id] = daughter

        # Add lineage edge
        self.edges.add((parent_id, daughter_id, RelationType.LINEAGE))

        # Healthy daughters get regulatory edges to parent neighbors
        if daughter.state == CellState.HEALTHY:
            for neighbor in self.get_neighbors(parent_id, RelationType.REGULATORY)[:2]:
                self.edges.add((daughter_id, neighbor, RelationType.REGULATORY))
            self.edges.add(("clock", daughter_id, RelationType.TEMPORAL))

        # Cancer daughters get NO regulatory or temporal edges
        # They only connect via lineage (the rec_succ step)
        return daughter

    def run(self, steps: int = 200, **kwargs):
        """Run the simulation for a given number of steps."""
        for _ in range(steps):
            self.step(**kwargs)

    # --- State queries ---

    def count_by_state(self) -> Dict[CellState, int]:
        """Count cells in each state."""
        counts = {s: 0 for s in CellState}
        for cell in self.cells.values():
            counts[cell.state] += 1
        return counts

    def total_cells(self) -> int:
        return len(self.cells)

    def edge_counts(self) -> Dict[RelationType, int]:
        """Count edges by type."""
        counts = {r: 0 for r in RelationType}
        for _, _, r in self.edges:
            counts[r] += 1
        return counts
