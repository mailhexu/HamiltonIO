"""Band-structure path generation, calculation, and plotting utilities."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
from ase import Atoms
from ase.cell import Cell


def _display_label(label: str) -> str:
    if label == "G":
        return r"$\Gamma$"
    return label


def _cell_and_pbc_from_structure(structure: Any) -> tuple[Cell, np.ndarray]:
    if isinstance(structure, Atoms):
        return structure.cell, np.asarray(structure.pbc, dtype=bool)
    if isinstance(structure, Cell):
        return structure, np.array([True, True, True], dtype=bool)
    return Cell(structure), np.array([True, True, True], dtype=bool)


def _resolve_structure(model: Any, structure: Any | None) -> Any:
    if structure is not None:
        return structure
    for attr in ("atoms", "cell"):
        if hasattr(model, attr):
            value = getattr(model, attr)
            if value is not None:
                return value
    raise ValueError(
        "No structure or cell was provided. Pass an ASE Atoms, ASE Cell, or 3x3 "
        "cell-like object, or set model.atoms/model.cell."
    )


def _group_band_path(bandpath, eps: float = 1e-8, shift: float = 0.15):
    xs, special_xs, labels = bandpath.get_linear_kpoint_axis()
    kpts = np.asarray(bandpath.kpts, dtype=float)

    discontinuities = xs[1:] - xs[:-1] < eps
    segments = [0] + list(np.where(discontinuities)[0] + 1) + [len(xs)]

    xcoords = []
    kpt_segments = []
    for iseg, (start, end) in enumerate(zip(segments[:-1], segments[1:])):
        xcoords.append(np.asarray(xs[start:end], dtype=float) + iseg * shift)
        kpt_segments.append(kpts[start:end])

    special_xs = np.asarray(special_xs, dtype=float).copy()
    duplicate_special = special_xs[1:] - special_xs[:-1] < eps
    for idx in np.where(duplicate_special)[0] + 1:
        special_xs[idx:] += shift

    if len(xcoords) == 1:
        xcoords = xcoords[0]

    return xcoords, kpt_segments, special_xs, labels


@dataclass
class BandPathData:
    """K-point path data for band-structure calculations."""

    kpoints: np.ndarray
    xcoords: np.ndarray | list[np.ndarray]
    kpath_labels: list[tuple[float, str]]
    special_points: dict[str, np.ndarray]
    path: str | None = None

    def __post_init__(self):
        self.kpoints = np.asarray(self.kpoints, dtype=float)
        self.special_points = {
            str(key): np.asarray(value, dtype=float)
            for key, value in self.special_points.items()
        }
        if isinstance(self.xcoords, list):
            self.xcoords = [np.asarray(x, dtype=float) for x in self.xcoords]
        else:
            self.xcoords = np.asarray(self.xcoords, dtype=float)
        self.kpath_labels = [(float(x), str(label)) for x, label in self.kpath_labels]


@dataclass
class ElectronicBandStructure:
    """Electronic band-structure data and plotting helper."""

    energies: np.ndarray
    kpoints: np.ndarray
    xcoords: np.ndarray | list[np.ndarray]
    kpath_labels: list[tuple[float, str]]
    special_points: dict[str, np.ndarray]
    fermi_energy: float | None = None
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self):
        self.energies = np.asarray(self.energies, dtype=float)
        self.kpoints = np.asarray(self.kpoints, dtype=float)
        self.special_points = {
            str(key): np.asarray(value, dtype=float)
            for key, value in self.special_points.items()
        }
        if isinstance(self.xcoords, list):
            self.xcoords = [np.asarray(x, dtype=float) for x in self.xcoords]
        else:
            self.xcoords = np.asarray(self.xcoords, dtype=float)
        self.kpath_labels = [(float(x), str(label)) for x, label in self.kpath_labels]

    def plot(
        self,
        ax=None,
        filename: str | Path | None = None,
        show: bool = False,
        shift: float = 0.0,
        fermi_energy: float | None = None,
        linewidth: float = 1.5,
        color: str = "blue",
        linestyle: str = "-",
        ylabel: str | None = None,
        ylim: tuple[float, float] | None = None,
        **kwargs,
    ):
        """Plot the band structure and return the Matplotlib axes."""
        import matplotlib.pyplot as plt

        if ax is None:
            _, ax = plt.subplots(constrained_layout=True)

        effective_fermi = self.fermi_energy if fermi_energy is None else fermi_energy
        energy_shift = float(shift)
        if effective_fermi is not None:
            energy_shift += float(effective_fermi)
        plotted = self.energies - energy_shift

        if isinstance(self.xcoords, list):
            start = 0
            for x in self.xcoords:
                stop = start + len(x)
                for band in plotted[start:stop].T:
                    ax.plot(
                        x,
                        band,
                        linewidth=linewidth,
                        color=color,
                        linestyle=linestyle,
                        **kwargs,
                    )
                start = stop
            xmin = self.xcoords[0][0]
            xmax = self.xcoords[-1][-1]
        else:
            for band in plotted.T:
                ax.plot(
                    self.xcoords,
                    band,
                    linewidth=linewidth,
                    color=color,
                    linestyle=linestyle,
                    **kwargs,
                )
            xmin = self.xcoords[0]
            xmax = self.xcoords[-1]

        ax.set_xlim(float(xmin), float(xmax))
        if ylim is None:
            ymin = float(np.min(plotted))
            ymax = float(np.max(plotted))
            pad = 0.05 * (ymax - ymin if ymax > ymin else 1.0)
            ylim = (ymin - pad, ymax + pad)
        ax.set_ylim(*ylim)

        tick_positions = [x for x, _ in self.kpath_labels]
        tick_labels = [label for _, label in self.kpath_labels]
        ax.set_xticks(tick_positions)
        ax.set_xticklabels(tick_labels)
        ax.vlines(
            tick_positions,
            ymin=ylim[0],
            ymax=ylim[1],
            color="black",
            linewidth=linewidth / 5,
        )
        ax.set_xlabel("Wave vector")
        if ylabel is None:
            ylabel = (
                "Energy - $E_F$ (eV)" if effective_fermi is not None else "Energy (eV)"
            )
        ax.set_ylabel(ylabel)

        if filename is not None:
            ax.figure.savefig(filename, dpi=300, bbox_inches="tight")
        if show:
            plt.show()
        return ax

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-serializable dictionary."""
        xcoords = (
            [x.tolist() for x in self.xcoords]
            if isinstance(self.xcoords, list)
            else self.xcoords.tolist()
        )
        return {
            "energies": self.energies.tolist(),
            "kpoints": self.kpoints.tolist(),
            "xcoords": xcoords,
            "kpath_labels": [(float(x), label) for x, label in self.kpath_labels],
            "special_points": {k: v.tolist() for k, v in self.special_points.items()},
            "fermi_energy": self.fermi_energy,
            "metadata": self.metadata,
        }

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "ElectronicBandStructure":
        """Build an electronic band structure from serialized data."""
        return cls(**data)


def make_band_path(
    structure: Any,
    path: str | None = None,
    npoints: int = 300,
    special_points: dict[str, Any] | None = None,
    eps: float = 2e-4,
) -> BandPathData:
    """Generate a high-symmetry band path using ASE."""
    cell, pbc = _cell_and_pbc_from_structure(structure)
    bandpath = cell.bandpath(
        path=path,
        npoints=npoints,
        special_points=special_points,
        eps=eps,
        pbc=pbc,
    )
    xcoords, kpt_segments, special_xs, labels = _group_band_path(bandpath)
    kpoints = np.concatenate(kpt_segments, axis=0)
    kpath_labels = [
        (float(x), _display_label(label)) for x, label in zip(special_xs, labels)
    ]
    return BandPathData(
        kpoints=kpoints,
        xcoords=xcoords,
        kpath_labels=kpath_labels,
        special_points=bandpath.special_points,
        path=getattr(bandpath, "path", path),
    )


def _extract_eigenvalues(result: Any, source: str) -> np.ndarray:
    if source == "solve_all":
        if isinstance(result, tuple):
            evals = result[0]
        else:
            evals = result
    else:
        if not isinstance(result, tuple) or len(result) < 3:
            raise TypeError("HS_and_eigen() must return a tuple containing eigenvalues")
        evals = result[2]
    evals = np.asarray(evals, dtype=float)
    if evals.ndim != 2:
        raise ValueError(
            f"Expected eigenvalues with shape (nkpts, nbands), got {evals.shape}"
        )
    return evals


def _solve_eigenvalues(model: Any, kpoints: np.ndarray, solver_kwargs: dict[str, Any]):
    if hasattr(model, "solve_all"):
        return _extract_eigenvalues(
            model.solve_all(kpoints, **solver_kwargs), "solve_all"
        )
    if hasattr(model, "HS_and_eigen"):
        return _extract_eigenvalues(
            model.HS_and_eigen(kpoints, **solver_kwargs), "HS_and_eigen"
        )
    raise TypeError(
        "Band structure calculation requires a model with solve_all(kpts, ...) "
        "or HS_and_eigen(kpts, ...)."
    )


def calculate_band_structure(
    model: Any,
    structure: Any | None = None,
    path: str | None = None,
    npoints: int = 300,
    special_points: dict[str, Any] | None = None,
    fermi_energy: float | None = None,
    solver_kwargs: dict[str, Any] | None = None,
) -> ElectronicBandStructure:
    """Calculate electronic bands along an ASE-generated high-symmetry path."""
    solver_kwargs = {} if solver_kwargs is None else dict(solver_kwargs)
    resolved_structure = _resolve_structure(model, structure)
    band_path = make_band_path(
        resolved_structure, path=path, npoints=npoints, special_points=special_points
    )
    energies = _solve_eigenvalues(model, band_path.kpoints, solver_kwargs)
    return ElectronicBandStructure(
        energies=energies,
        kpoints=band_path.kpoints,
        xcoords=band_path.xcoords,
        kpath_labels=band_path.kpath_labels,
        special_points=band_path.special_points,
        fermi_energy=fermi_energy,
        metadata={
            "model_class": model.__class__.__name__,
            "path": band_path.path,
            "npoints": npoints,
            "units": {"energies": "eV"},
            "solver_kwargs": solver_kwargs,
        },
    )


def plot_band_structure(
    model: Any,
    structure: Any | None = None,
    filename: str | Path | None = None,
    show: bool = False,
    ax=None,
    **kwargs,
):
    """Calculate and plot a band structure.

    Returns
    -------
    tuple
        ``(band_structure, ax)`` where ``band_structure`` is an
        :class:`ElectronicBandStructure` and ``ax`` is the Matplotlib axes.
    """
    calculate_keys = {
        "path",
        "npoints",
        "special_points",
        "fermi_energy",
        "solver_kwargs",
    }
    calculate_kwargs = {
        key: kwargs.pop(key) for key in list(kwargs) if key in calculate_keys
    }
    plot_kwargs = kwargs
    if "fermi_energy" in calculate_kwargs:
        plot_kwargs["fermi_energy"] = calculate_kwargs["fermi_energy"]
    bands = calculate_band_structure(model, structure=structure, **calculate_kwargs)
    ax = bands.plot(ax=ax, filename=filename, show=show, **plot_kwargs)
    return bands, ax
