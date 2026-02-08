"""
Coherence and health metrics for the Boundary Logic Failure model.

The coherence metric kappa summarizes three tissue-level properties:
  (i)   density of regulatory signaling
  (ii)  temporal alignment
  (iii) identity dispersion

Reference: Cancer_As_Boundary_Logic_Failure.tex, Section 9 (Eq. 6-7)
"""

import math
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .tissue import Tissue


def compute_identity_entropy(tissue: "Tissue") -> float:
    """Compute mean Shannon entropy of identity distribution across cells.

    Each cell has p_align in [0,1]. We treat (p_align, 1-p_align) as a
    binary distribution and compute its entropy. High entropy = cells are
    uncertain about identity. Low entropy = cells are committed (either
    to organism or to self).

    Returns:
        H: Mean entropy in [0, 1] (normalized by log2).
    """
    if not tissue.cells:
        return 0.0

    total_entropy = 0.0
    for cell in tissue.cells.values():
        p = cell.p_align
        # Clamp to avoid log(0)
        p = max(1e-10, min(1.0 - 1e-10, p))
        h = -(p * math.log2(p) + (1 - p) * math.log2(1 - p))
        total_entropy += h

    return total_entropy / len(tissue.cells)


def compute_signaling_density(tissue: "Tissue") -> float:
    """Fraction of edges that are regulatory (signaling).

    High density = strong cross-cell communication.
    Low density = decoupled, isolated cells.
    """
    from .tissue import RelationType

    total = len(tissue.edges)
    if total == 0:
        return 0.0

    sig_count = sum(1 for _, _, r in tissue.edges if r == RelationType.REGULATORY)
    return sig_count / total


def compute_temporal_density(tissue: "Tissue") -> float:
    """Fraction of edges that are temporal (circadian/hormonal sync).

    Measures how much of the tissue is coupled to organism-level
    timing signals.
    """
    from .tissue import RelationType

    total = len(tissue.edges)
    if total == 0:
        return 0.0

    tmp_count = sum(1 for _, _, r in tissue.edges if r == RelationType.TEMPORAL)
    return tmp_count / total


def compute_coherence(
    tissue: "Tissue",
    w_sig: float = 0.4,
    w_tmp: float = 0.2,
    w_identity: float = 0.4,
) -> float:
    """Compute tissue coherence score kappa in [0, 1].

    kappa = w_sig * rho_sig + w_tmp * rho_tmp + w_identity * (1 - H)

    where:
      rho_sig  = regulatory signaling density
      rho_tmp  = temporal coordination density
      H        = mean identity entropy

    High kappa (~1.0) = healthy, coordinated tissue.
    Low kappa (~0.0)  = decoupled, malignant tissue.

    Args:
        tissue: The Tissue object.
        w_sig: Weight for signaling density.
        w_tmp: Weight for temporal density.
        w_identity: Weight for identity coherence (1 - H).

    Returns:
        kappa: Coherence score clipped to [0, 1].
    """
    rho_sig = compute_signaling_density(tissue)
    rho_tmp = compute_temporal_density(tissue)
    H = compute_identity_entropy(tissue)

    kappa = w_sig * rho_sig + w_tmp * rho_tmp + w_identity * (1.0 - H)
    return max(0.0, min(1.0, kappa))


def compute_organism_alignment(tissue: "Tissue") -> float:
    """Mean organism alignment across all cells.

    Simple average of p_align values. Quick health indicator.
    """
    if not tissue.cells:
        return 0.0
    return sum(c.p_align for c in tissue.cells.values()) / len(tissue.cells)
