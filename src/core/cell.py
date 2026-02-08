"""
Cell representation for the Boundary Logic Failure cancer model.

Each cell carries:
  - A discrete state (HEALTHY, PRE_CANCER, CANCER)
  - An identity-alignment scalar p_align in [0,1] (organism vs self)
  - An energetic cost proxy epsilon > 0
  - A division-rate proxy r >= 0

Reference: Cancer_As_Boundary_Logic_Failure.tex, Section 9 (Toy Formalism)
"""

from enum import Enum
from dataclasses import dataclass, field
from typing import Dict


class CellState(Enum):
    HEALTHY = "HEALTHY"
    PRE_CANCER = "PRE_CANCER"
    CANCER = "CANCER"


@dataclass
class Cell:
    """A single cell in the tissue network.

    Attributes:
        id: Unique identifier for this cell.
        state: Current state (HEALTHY, PRE_CANCER, CANCER).
        p_align: Identity alignment scalar [0, 1].
            1.0 = fully organism-aligned.
            0.0 = fully self-aligned (cancerous).
        energy_cost: Metabolic cost per timestep. Cancer cells draw
            more from the global budget (Warburg effect proxy).
        division_rate: Probability of division per timestep.
            Healthy cells are constrained by contact inhibition;
            cancer cells are not.
        perspective: Decomposition of identity between self and organism.
            Used for computing identity entropy H.
    """

    id: str
    state: CellState = CellState.HEALTHY
    p_align: float = 1.0
    energy_cost: float = 1.0
    division_rate: float = 0.01
    perspective: Dict[str, float] = field(
        default_factory=lambda: {"self": 0.1, "organism": 0.9}
    )

    def transform_to_precancer(self):
        """Initiate boundary weakening. The cell begins losing
        organism-level coupling but retains partial alignment."""
        self.state = CellState.PRE_CANCER
        self.p_align = 0.5
        self.perspective = {"self": 0.4, "organism": 0.6}

    def transform_to_cancer(self):
        """Complete boundary collapse. The cell enters the
        self-referential loop (rec_succ without termination).

        Effects:
          - Identity collapses to pure self-reference
          - Division rate increases (no contact inhibition)
          - Energy cost increases (Warburg effect)
        """
        self.state = CellState.CANCER
        self.p_align = 0.05
        self.division_rate = 0.2
        self.energy_cost = 2.7  # Illustrative Warburg burden proxy (not fitted)
        self.perspective = {"self": 0.95, "organism": 0.05}

    def create_daughter(self, daughter_id: str) -> "Cell":
        """Division: create a daughter cell inheriting the parent state.

        For healthy cells, the daughter acquires a new identity
        (differentiation). For cancer cells, the daughter remains
        insufficiently distinguishable from the parent for robust
        organism-level control (duplication with weakened distinction).
        """
        daughter = Cell(
            id=daughter_id,
            state=self.state,
            p_align=self.p_align,
            energy_cost=self.energy_cost,
            division_rate=self.division_rate,
            perspective=dict(self.perspective),
        )
        return daughter
