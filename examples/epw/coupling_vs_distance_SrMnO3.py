"""Minimal example: plot |g| vs distance for SrMnO3 EPW coupling.

Usage:
    uv run python examples/epw/coupling_vs_distance_SrMnO3.py

This uses SrMnO3 EPW data and the Epmat helpers to:
    1. Load the electron-phonon coupling matrix elements.
    2. Compute couplings annotated with distances (both Rk and Rg metrics).
    3. Plot:
       - Scatter of |g| vs distance for Rk-based (WF-WF) distances.
       - Scatter of |g| vs distance for Rg-based (WF-atom) distances.
       - Binned average |g| vs distance for both metrics.
"""

import os

import matplotlib.pyplot as plt

from HamiltonIO.epw.epwparser import Epmat, read_crystal_fmt


def main() -> None:
    # Project root (this file is examples/epw/...)
    base = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
    data_dir = os.path.join(base, "HamiltonIO", "epw", "test", "up")

    # Read crystal structure
    crystal = read_crystal_fmt(os.path.join(data_dir, "crystal.fmt"))

    # Read EPW data
    ep = Epmat()
    ep.crystal = crystal
    ep.read(path=data_dir, prefix="SrMnO3", epmat_ncfile="epmat.nc")

    # Use cell from crystal structure (convert to Cartesian, alat is in Bohr)
    from ase.units import Bohr

    cell = crystal.at.reshape(3, 3) * crystal.alat * Bohr

    # Analyze mode 0 (first phonon mode)
    imode = 0
    print(f"Analyzing mode {imode}")
    print(f"Crystal has {crystal.natom} atoms, {crystal.nmode} modes")
    print(f"Epmat has {ep.nwann} Wannier functions")

    # 1) Rk-based distances (electron WF to WF)
    print("\nExtracting Rk-based couplings (WF-WF distances)...")
    entries_Rk = ep.distance_resolved_couplings_Rk(imode=imode, cell=cell)
    print(f"Found {len(entries_Rk)} non-zero coupling elements")

    distances_Rk = [e["distance"] for e in entries_Rk]
    mags_Rk = [abs(e["coupling"]) for e in entries_Rk]

    # 2) Rg-based distances (electron WF to atomic displacement)
    print("\nExtracting Rg-based couplings (WF-atom distances)...")
    entries_Rg = ep.distance_resolved_couplings_Rg(imode=imode, cell=cell)
    print(f"Found {len(entries_Rg)} non-zero coupling elements")

    distances_Rg = [e["distance"] for e in entries_Rg]
    mags_Rg = [abs(e["coupling"]) for e in entries_Rg]

    # 3) Plot Rk-based scatter
    plt.figure(figsize=(7, 5))

    plt.subplot(2, 2, 1)
    plt.scatter(distances_Rk, mags_Rk, s=3, alpha=0.6, label="Rk-based")
    plt.xlabel("Distance (Å)", fontsize=9)
    plt.ylabel(r"|$g$| (eV/Å)", fontsize=9)
    plt.yscale("log")
    plt.title(f"Rk-based (WF-WF) - Mode {imode}", fontsize=9)
    plt.tick_params(labelsize=8)
    plt.legend(fontsize=8)

    # 4) Plot Rg-based scatter
    plt.subplot(2, 2, 2)
    plt.scatter(distances_Rg, mags_Rg, s=3, alpha=0.6, color="orange", label="Rg-based")
    plt.xlabel("Distance (Å)", fontsize=9)
    plt.ylabel(r"|$g$| (eV/Å)", fontsize=9)
    plt.yscale("log")
    plt.title(f"Rg-based (WF-atom) - Mode {imode}", fontsize=9)
    plt.tick_params(labelsize=8)
    plt.legend(fontsize=8)

    # 5) Binned Rk-based
    centers_Rk, avg_Rk = Epmat.bin_couplings_by_distance(entries_Rk, dr=0.2)

    plt.subplot(2, 2, 3)
    plt.plot(centers_Rk, avg_Rk, marker="o", linestyle="-", linewidth=1, markersize=3)
    plt.xlabel("Distance (Å)", fontsize=9)
    plt.ylabel(r"Average |$g$| (eV/Å)", fontsize=9)
    plt.title(f"Binned Rk-based - Mode {imode}", fontsize=9)
    plt.tick_params(labelsize=8)

    # 6) Binned Rg-based
    centers_Rg, avg_Rg = Epmat.bin_couplings_by_distance(entries_Rg, dr=0.2)

    plt.subplot(2, 2, 4)
    plt.plot(
        centers_Rg,
        avg_Rg,
        marker="o",
        linestyle="-",
        linewidth=1,
        markersize=3,
        color="orange",
    )
    plt.xlabel("Distance (Å)", fontsize=9)
    plt.ylabel(r"Average |$g$| (eV/Å)", fontsize=9)
    plt.title(f"Binned Rg-based - Mode {imode}", fontsize=9)
    plt.tick_params(labelsize=8)

    plt.tight_layout()
    plt.savefig("coupling_vs_distance.pdf", dpi=300, bbox_inches="tight")
    print("\nPlot saved to coupling_vs_distance.pdf")


if __name__ == "__main__":
    main()
