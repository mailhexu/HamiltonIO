"""Minimal example: plot |H_ij| vs distance for SrMnO3 Wannier90 TB.

Usage:
    uv run python examples/wannier/hoppings_vs_distance_SrMnO3.py

This uses data from tests/data/SrMnO3_wannier/ and the WannierHam helpers to:
    1. Load the tight-binding model.
    2. Compute unique hoppings annotated with Wannier-center distances.
    3. Plot:
       - Scatter of |H_ij| vs distance for all unique hoppings.
       - Binned average |H_ij| vs distance.
"""

import os

import matplotlib.pyplot as plt
from ase.io import read

from HamiltonIO.wannier import WannierHam


def main() -> None:
    # Project root (this file is examples/wannier/...)
    base = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
    data_dir = os.path.join(base, "tests", "data", "SrMnO3_wannier")
    prefix = "SrMnO3"

    # Read atoms from .pwi file using ASE
    pwi_file = os.path.join(data_dir, "SrMnO3.pwi")
    atoms = read(pwi_file, format="espresso-in")

    # Read Wannier TB model, passing atoms
    ham = WannierHam.read_from_wannier_dir(data_dir, prefix, atoms=atoms)

    # Use cell from underlying atoms object for Cartesian conversion
    cell = atoms.get_cell()

    # Distance-resolved unique hoppings (Wannier centers, unique R/i/j convention)
    entries = ham.distance_resolved_hoppings(cell=cell)

    # 1) Scatter: all unique hoppings
    distances = [e["distance"] for e in entries]
    mags = [abs(e["hopping"]) for e in entries]

    plt.figure(figsize=(3.4, 2.5))
    plt.scatter(distances, mags, s=3, alpha=0.6)
    plt.xlabel("Distance (Å)", fontsize=9)
    plt.ylabel(r"|$H_{ij}$| (eV)", fontsize=9)
    plt.yscale("log")
    plt.tick_params(labelsize=8)
    plt.tight_layout()
    plt.savefig("hoppings_scatter.pdf", dpi=300, bbox_inches="tight")

    # 2) Binned average |H_ij|
    centers, avg = WannierHam.bin_hoppings_by_distance(entries, dr=0.2)

    plt.figure(figsize=(3.4, 2.5))
    plt.plot(centers, avg, marker="o", linestyle="-", linewidth=1, markersize=3)
    plt.xlabel("Distance (Å)", fontsize=9)
    plt.ylabel(r"Average |$H_{ij}$| (eV)", fontsize=9)
    plt.tick_params(labelsize=8)
    plt.tight_layout()
    plt.savefig("hoppings_binned.pdf", dpi=300, bbox_inches="tight")

    plt.show()


if __name__ == "__main__":
    main()
