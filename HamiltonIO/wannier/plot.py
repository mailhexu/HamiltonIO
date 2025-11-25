"""Plotting utilities for Wannier90 data."""

import os

from ase.io import read

from .wannier_hamiltonian import WannierHam


def plot_wannier_distance(
    ax,
    path,
    prefix,
    structure_file=None,
    structure_format=None,
    ylim=(1e-5, 10),
    **scatter_kwargs,
):
    """
    Plot Wannier hopping matrix elements vs. distance.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The matplotlib axes object to plot on.
    path : str
        Directory containing the Wannier90 files.
    prefix : str
        Prefix for Wannier90 files (e.g., "wannier90" for wannier90_hr.dat).
    structure_file : str, optional
        Structure file to read. If None, looks for common files.
    structure_format : str, optional
        Format of structure file (e.g., 'vasp', 'cif', 'espresso-in').
        If None, will be inferred from file extension.
    ylim : tuple, optional
        The y-axis limits for the plot. Default is (1e-5, 10).
    **scatter_kwargs : dict
        Additional keyword arguments to pass to `ax.scatter`.

    Returns
    -------
    ax : matplotlib.axes.Axes
        The configured matplotlib axes.
    """
    # Load structure
    if structure_file is None:
        # Try common structure file names
        for fname in ["POSCAR", "structure.cif", "input.pwi", "structure.vasp"]:
            fpath = os.path.join(path, fname)
            if os.path.exists(fpath):
                structure_file = fname
                break
        if structure_file is None:
            raise FileNotFoundError(
                "No structure file found. Please specify structure_file."
            )

    atoms = read(os.path.join(path, structure_file), format=structure_format)

    # Load Wannier Hamiltonian
    ham = WannierHam.read_from_wannier_dir(path, prefix, atoms=atoms)

    # Get distance-resolved hoppings
    from ase import Atoms

    if isinstance(atoms, Atoms):
        cell = atoms.get_cell()
    else:
        cell = atoms[0].get_cell()
    entries = ham.distance_resolved_hoppings(cell)

    # Extract data
    distances = [e["distance"] for e in entries]
    magnitudes = [abs(e["hopping"]) for e in entries]

    # Default scatter settings
    defaults = {"s": 5, "alpha": 0.6}
    defaults.update(scatter_kwargs)

    # Plot
    ax.scatter(distances, magnitudes, **defaults)
    ax.set_xlabel("Distance (Å)")
    ax.set_ylabel("|t| (eV)")
    ax.set_yscale("log")
    ax.set_ylim(ylim)

    return ax
