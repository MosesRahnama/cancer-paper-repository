"""
Evidence-aware parameter profiles for simulation runs.

These profiles map literature-level annual progression rates to per-step
probabilities. They do not claim lesion-specific molecular fidelity; they
provide transparent, reproducible parameterization for the toy model.
"""

from dataclasses import dataclass
from typing import Dict, Tuple


def annual_to_step_probability(annual_rate: float, steps_per_year: int = 12) -> float:
    """Convert annual risk to per-step Bernoulli probability."""
    if annual_rate < 0.0 or annual_rate > 1.0:
        raise ValueError("annual_rate must be in [0, 1]")
    if steps_per_year <= 0:
        raise ValueError("steps_per_year must be positive")
    return 1.0 - (1.0 - annual_rate) ** (1.0 / steps_per_year)


@dataclass(frozen=True)
class ParameterProfile:
    name: str
    step_unit: str
    p_mutation: float
    p_progression: float
    p_healthy_division: float
    source_refs: Tuple[str, ...]
    note: str

    def run_kwargs(self) -> Dict[str, float]:
        return {
            "p_mutation": self.p_mutation,
            "p_progression": self.p_progression,
            "p_healthy_division": self.p_healthy_division,
        }


AGGRESSIVE_DEMO_PROFILE = ParameterProfile(
    name="aggressive-demo",
    step_unit="abstract step",
    p_mutation=0.001,
    p_progression=0.05,
    p_healthy_division=0.01,
    source_refs=(),
    note="Legacy demo profile retained for backwards compatibility.",
)


CALIBRATED_MONTHLY_PROFILE = ParameterProfile(
    name="calibrated-monthly-oed",
    step_unit="1 step = 1 month",
    p_mutation=0.002,
    p_progression=annual_to_step_probability(0.105, steps_per_year=12),
    p_healthy_division=0.01,
    source_refs=(
        "10.1155/2015/854636",
        "10.1101/cshperspect.a026542",
        "10.1126/science.aab4082",
    ),
    note=(
        "Progression anchored to pooled oral epithelial dysplasia malignant "
        "transformation rate (~10.5% annual) converted to monthly hazard."
    ),
)
