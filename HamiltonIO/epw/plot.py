"""Plotting utilities for EPW data."""

import os
from .epwparser import Epmat, read_crystal_fmt
from ase.units import Bohr


def plot_epw_distance(
    ax,
    path,
    imode,
    distance_type="Rk",
    epmat_ncfile="epmat.nc",
    crystal_fmt_file="crystal.fmt",
    ylim=(1e-5, 10),
    **scatter_kwargs,
):
    """
    Plot electron-phonon coupling matrix elements vs. distance.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The matplotlib axes object to plot on.
    path : str
        Directory containing the EPW files.
    imode : int
        The phonon mode index to analyze.
    distance_type : str, optional
        The type of distance to plot, either 'Rk' (WF-WF) or 'Rg' (WF-atom).
        Default is 'Rk'.
    epmat_ncfile : str, optional
        Name of the NetCDF file. Default is 'epmat.nc'.
    crystal_fmt_file : str, optional
        Name of the crystal format file. Default is 'crystal.fmt'.
    ylim : tuple, optional
        The y-axis limits for the plot. Default is (1e-5, 10).
    **scatter_kwargs : dict
        Additional keyword arguments to pass to `ax.scatter`.
    """
    # Load data
    crystal = read_crystal_fmt(os.path.join(path, crystal_fmt_file))
    ep = Epmat()
    ep.crystal = crystal
    ep.read(
        path=path, prefix="prefix", epmat_ncfile=epmat_ncfile
    )  # prefix is not used when ncfile is present

    # Get cell
    cell = crystal.at.reshape(3, 3) * crystal.alat * Bohr

    # Calculate distances
    if distance_type == "Rk":
        entries = ep.distance_resolved_couplings_Rk(imode=imode, cell=cell)
    elif distance_type == "Rg":
        entries = ep.distance_resolved_couplings_Rg(imode=imode, cell=cell)
    else:
        raise ValueError("distance_type must be 'Rk' or 'Rg'")

    distances = [e["distance"] for e in entries]
    magnitudes = [abs(e["coupling"]) for e in entries]

    # Default scatter settings
    defaults = {"s": 5, "alpha": 0.6}
    defaults.update(scatter_kwargs)

    # Plot
    ax.scatter(distances, magnitudes, **defaults)
    ax.set_xlabel("Distance (Å)")
    ax.set_ylabel(r"|$g$| (eV/Å)")
    ax.set_yscale("log")
    ax.set_ylim(ylim)

    return ax
