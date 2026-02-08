"""
Visualization functions for the Boundary Logic Failure model.

Generates:
- Coherence vs cancer count dual-axis plots
- Tissue network graphs (color-coded by state)
- State distribution bar charts
"""

from typing import List, Optional, TYPE_CHECKING

if TYPE_CHECKING:
    from simulations.engine import Simulation
    from core.tissue import Tissue


def plot_coherence_vs_cellcount(
    sim: "Simulation",
    title: str = "Cancer Boundary Logic: Coherence Collapse vs Tumor Growth",
    save_path: Optional[str] = None,
    show: bool = False,
):
    """
    Dual-axis plot: kappa (blue) and cancer count (red).

    The crossing point where kappa decline and cancer count acceleration
    meet represents the "Event Horizon of the Disease" -- where tumor
    entropy production exceeds host dissipation capacity.
    """
    import matplotlib.pyplot as plt

    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(12, 10), sharex=True)

    steps = range(len(sim.history))
    kappa = sim.kappa_series()
    cancer = sim.cancer_series()
    energy = sim.energy_series()

    # Coherence
    ax1.plot(steps, kappa, "b-", linewidth=2, label="Coherence (kappa)")
    ax1.axhline(y=0.5, color="orange", linestyle="--", alpha=0.5, label="Pre-cancer threshold")
    ax1.axhline(y=0.3, color="red", linestyle="--", alpha=0.5, label="Malignancy threshold")
    ax1.set_ylabel("Tissue Coherence (kappa)")
    ax1.set_ylim(-0.05, 1.05)
    ax1.legend(loc="upper right")
    ax1.set_title(title)
    ax1.grid(True, alpha=0.3)

    # Cancer count
    ax2.plot(steps, cancer, "r-", linewidth=2, label="Cancer cells")
    ax2.set_ylabel("Cancer Cell Count")
    ax2.legend(loc="upper left")
    ax2.grid(True, alpha=0.3)

    # Energy budget (cachexia)
    ax3.plot(steps, energy, "g-", linewidth=2, label="Energy budget")
    ax3.set_ylabel("Energy Budget")
    ax3.set_xlabel("Time Steps")
    ax3.legend(loc="upper right")
    ax3.grid(True, alpha=0.3)

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"Saved to {save_path}")
    if show:
        plt.show()
    else:
        plt.close(fig)


def plot_tissue_network(
    tissue: "Tissue",
    title: str = "Tissue Network",
    save_path: Optional[str] = None,
    show: bool = False,
):
    """
    Network graph visualization. Nodes colored by state:
    green=healthy, yellow=pre-cancer, red=cancer, blue=clock.
    """
    import matplotlib.pyplot as plt
    import networkx as nx
    from core.cell import CellState

    G = nx.Graph()

    color_map = {
        CellState.HEALTHY: "green",
        CellState.PRE_CANCER: "yellow",
        CellState.CANCER: "red",
    }

    for cid, cell in tissue.cells.items():
        G.add_node(cid, color=color_map.get(cell.state, "gray"), size=200)

    # Add clock node
    if any(s == "clock" or t == "clock" for s, t, _ in tissue.edges):
        G.add_node("clock", color="blue", size=600)

    for s, t, r in tissue.edges:
        if s in G.nodes and t in G.nodes:
            G.add_edge(s, t)

    pos = nx.spring_layout(G, seed=42)
    colors = [G.nodes[n].get("color", "gray") for n in G.nodes()]
    sizes = [G.nodes[n].get("size", 200) for n in G.nodes()]

    fig = plt.figure(figsize=(12, 8))
    nx.draw(
        G, pos,
        node_color=colors,
        node_size=sizes,
        with_labels=False,
        alpha=0.8,
        edge_color="lightgray",
    )
    plt.title(f"{title}\nGreen=Healthy  Yellow=PreCancer  Red=Cancer  Blue=Clock")

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"Saved to {save_path}")
    if show:
        plt.show()
    else:
        plt.close(fig)


def plot_state_distribution(
    sim: "Simulation",
    title: str = "Cell State Distribution Over Time",
    save_path: Optional[str] = None,
    show: bool = False,
):
    """Stacked area chart showing healthy/pre-cancer/cancer proportions."""
    import matplotlib.pyplot as plt

    steps = range(len(sim.history))
    healthy = [h["healthy"] for h in sim.history]
    pre_cancer = [h["pre_cancer"] for h in sim.history]
    cancer = [h["cancer"] for h in sim.history]

    fig = plt.figure(figsize=(10, 6))
    plt.stackplot(
        steps, healthy, pre_cancer, cancer,
        labels=["Healthy", "Pre-Cancer", "Cancer"],
        colors=["green", "gold", "red"],
        alpha=0.7,
    )
    plt.xlabel("Time Steps")
    plt.ylabel("Cell Count")
    plt.title(title)
    plt.legend(loc="upper left")
    plt.grid(True, alpha=0.3)

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"Saved to {save_path}")
    if show:
        plt.show()
    else:
        plt.close(fig)
