from pathlib import Path

import matplotlib
import numpy as np
import pytest
from ase import Atoms

matplotlib.use("Agg")

from HamiltonIO.bandstructure import (  # noqa: E402
    ElectronicBandStructure,
    calculate_band_structure,
    make_band_path,
    plot_band_structure,
)


class SolveAllModel:
    def __init__(self):
        self.calls = []

    def solve_all(self, kpts, scale=1.0):
        self.calls.append((np.array(kpts), scale))
        kpts = np.asarray(kpts)
        energies = scale * np.column_stack([kpts[:, 0], kpts[:, 0] + kpts[:, 1] + 1.0])
        evecs = np.zeros((len(kpts), 2, 2), dtype=complex)
        return energies, evecs


class HSEModel:
    def HS_and_eigen(self, kpts, offset=0.0):
        kpts = np.asarray(kpts)
        hams = np.zeros((len(kpts), 2, 2), dtype=complex)
        evals = np.column_stack([kpts[:, 2] + offset, kpts[:, 2] + offset + 2.0])
        evecs = np.zeros((len(kpts), 2, 2), dtype=complex)
        return hams, None, evals, evecs


class UnsupportedModel:
    pass


@pytest.fixture
def cubic_atoms():
    return Atoms("H", cell=[1.0, 1.0, 1.0], pbc=True)


def test_make_band_path_automatic(cubic_atoms):
    path = make_band_path(cubic_atoms, npoints=12)

    assert path.kpoints.shape == (12, 3)
    assert len(path.kpath_labels) >= 2
    assert path.kpath_labels[0][1] == r"$\Gamma$"
    assert "G" in path.special_points
    assert isinstance(path.path, str)


def test_make_band_path_respects_atoms_pbc():
    atoms = Atoms("H", cell=[1.0, 1.0, 10.0], pbc=[True, True, False])

    path = make_band_path(atoms, npoints=12)

    assert path.path == "MGXM"
    assert "R" not in path.special_points


def test_make_band_path_explicit_path(cubic_atoms):
    path = make_band_path(cubic_atoms, path="GX", npoints=8)

    assert path.kpoints.shape == (8, 3)
    assert path.kpath_labels[0][1] == r"$\Gamma$"
    assert path.kpath_labels[-1][1] == "X"
    np.testing.assert_allclose(path.kpoints[0], [0.0, 0.0, 0.0])
    np.testing.assert_allclose(path.kpoints[-1], path.special_points["X"])


def test_make_band_path_segmented_path(cubic_atoms):
    path = make_band_path(cubic_atoms, path="GX,MR", npoints=10)

    assert isinstance(path.xcoords, list)
    assert len(path.xcoords) == 2
    assert sum(len(x) for x in path.xcoords) == len(path.kpoints)
    assert path.kpath_labels[0][1] == r"$\Gamma$"
    assert path.kpath_labels[-1][1] == "R"


def test_electronic_band_structure_preserves_raw_energies():
    energies = np.array([[1.0, 2.0], [3.0, 4.0]])
    bands = ElectronicBandStructure(
        energies=energies,
        kpoints=np.zeros((2, 3)),
        xcoords=np.array([0.0, 1.0]),
        kpath_labels=[(0.0, r"$\Gamma$"), (1.0, "X")],
        special_points={"G": np.zeros(3), "X": np.array([0.0, 0.5, 0.0])},
        fermi_energy=1.5,
    )

    np.testing.assert_allclose(bands.energies, energies)


def test_electronic_band_structure_serialization_roundtrip():
    bands = ElectronicBandStructure(
        energies=np.array([[1.0, 2.0], [3.0, 4.0]]),
        kpoints=np.zeros((2, 3)),
        xcoords=np.array([0.0, 1.0]),
        kpath_labels=[(0.0, r"$\Gamma$"), (1.0, "X")],
        special_points={"G": np.zeros(3), "X": np.array([0.0, 0.5, 0.0])},
        fermi_energy=1.5,
        metadata={"path": "GX"},
    )

    restored = ElectronicBandStructure.from_dict(bands.to_dict())

    np.testing.assert_allclose(restored.energies, bands.energies)
    np.testing.assert_allclose(restored.kpoints, bands.kpoints)
    assert restored.kpath_labels == bands.kpath_labels
    assert restored.fermi_energy == 1.5
    assert restored.metadata["path"] == "GX"


def test_calculate_band_structure_uses_solve_all(cubic_atoms):
    model = SolveAllModel()
    bands = calculate_band_structure(
        model, cubic_atoms, path="GX", npoints=6, solver_kwargs={"scale": 2.0}
    )

    assert bands.energies.shape == (6, 2)
    assert len(model.calls) == 1
    assert model.calls[0][1] == 2.0
    np.testing.assert_allclose(bands.energies[:, 0], 2.0 * bands.kpoints[:, 0])


def test_calculate_band_structure_falls_back_to_hs_and_eigen(cubic_atoms):
    bands = calculate_band_structure(
        HSEModel(), cubic_atoms, path="GX", npoints=5, solver_kwargs={"offset": 3.0}
    )

    assert bands.energies.shape == (5, 2)
    np.testing.assert_allclose(bands.energies[:, 0], bands.kpoints[:, 2] + 3.0)


def test_calculate_band_structure_rejects_unsupported_model(cubic_atoms):
    with pytest.raises(TypeError, match="solve_all.*HS_and_eigen"):
        calculate_band_structure(UnsupportedModel(), cubic_atoms, path="GX", npoints=5)


def test_calculate_band_structure_resolves_model_atoms(cubic_atoms):
    model = SolveAllModel()
    model.atoms = cubic_atoms

    bands = calculate_band_structure(model, path="GX", npoints=5)

    assert bands.kpoints.shape == (5, 3)
    assert bands.metadata["model_class"] == "SolveAllModel"


def test_public_package_exports():
    import HamiltonIO

    assert HamiltonIO.ElectronicBandStructure is ElectronicBandStructure
    assert HamiltonIO.make_band_path is make_band_path
    assert HamiltonIO.calculate_band_structure is calculate_band_structure
    assert HamiltonIO.plot_band_structure is plot_band_structure


def test_plot_saves_file_and_shifts_display_only(tmp_path):
    filename = tmp_path / "bands.png"
    energies = np.array([[1.0, 2.0], [2.0, 3.0], [3.0, 4.0]])
    bands = ElectronicBandStructure(
        energies=energies.copy(),
        kpoints=np.zeros((3, 3)),
        xcoords=np.array([0.0, 0.5, 1.0]),
        kpath_labels=[(0.0, r"$\Gamma$"), (1.0, "X")],
        special_points={"G": np.zeros(3), "X": np.array([0.0, 0.5, 0.0])},
    )

    ax = bands.plot(filename=filename, fermi_energy=2.0)

    assert filename.exists()
    assert ax.get_ylabel() == "Energy - $E_F$ (eV)"
    assert [tick.get_text() for tick in ax.get_xticklabels()] == [r"$\Gamma$", "X"]
    np.testing.assert_allclose(bands.energies, energies)


def test_plot_uses_stored_fermi_energy_by_default(tmp_path):
    filename = tmp_path / "stored_fermi.png"
    bands = ElectronicBandStructure(
        energies=np.array([[1.0, 2.0], [2.0, 3.0]]),
        kpoints=np.zeros((2, 3)),
        xcoords=np.array([0.0, 1.0]),
        kpath_labels=[(0.0, r"$\Gamma$"), (1.0, "X")],
        special_points={"G": np.zeros(3), "X": np.array([0.0, 0.5, 0.0])},
        fermi_energy=1.0,
    )

    ax = bands.plot(filename=filename)

    assert filename.exists()
    assert ax.get_ylabel() == "Energy - $E_F$ (eV)"
    np.testing.assert_allclose(ax.lines[0].get_ydata(), [0.0, 1.0])


def test_plot_segmented_path(tmp_path):
    filename = tmp_path / "segmented.png"
    bands = ElectronicBandStructure(
        energies=np.array([[0.0, 1.0], [1.0, 2.0], [3.0, 4.0], [4.0, 5.0]]),
        kpoints=np.zeros((4, 3)),
        xcoords=[np.array([0.0, 1.0]), np.array([1.15, 2.15])],
        kpath_labels=[(0.0, r"$\Gamma$"), (1.0, "X"), (1.15, "M"), (2.15, "R")],
        special_points={"G": np.zeros(3), "X": np.array([0.0, 0.5, 0.0])},
    )

    ax = bands.plot(filename=filename)

    assert filename.exists()
    assert len(ax.lines) == 4
    np.testing.assert_allclose(ax.lines[0].get_xdata(), [0.0, 1.0])
    np.testing.assert_allclose(ax.lines[2].get_xdata(), [1.15, 2.15])


def test_plot_band_structure_one_call(cubic_atoms, tmp_path):
    filename = Path(tmp_path) / "one_call.png"

    bands, ax = plot_band_structure(
        SolveAllModel(), cubic_atoms, path="GX", npoints=5, filename=filename
    )

    assert filename.exists()
    assert bands.energies.shape == (5, 2)
    assert ax.get_xlabel() == "Wave vector"


def test_plot_band_structure_forwards_matplotlib_kwargs(cubic_atoms, tmp_path):
    filename = Path(tmp_path) / "kwargs.png"

    _, ax = plot_band_structure(
        SolveAllModel(),
        cubic_atoms,
        path="GX",
        npoints=5,
        filename=filename,
        alpha=0.25,
        marker="o",
    )

    assert filename.exists()
    assert ax.lines[0].get_alpha() == 0.25
    assert ax.lines[0].get_marker() == "o"
