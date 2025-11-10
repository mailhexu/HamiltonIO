"""Tests for distance-resolved coupling analysis in EPW electron-phonon data.

Tests cover:
    - Distance computation from Wannier centers (Rk-based)
    - Distance computation from atomic positions (Rg-based)
    - Binning of coupling magnitudes by distance
    - Correctness of mode-to-atom mapping
"""

from pathlib import Path

import numpy as np
import pytest
from ase.units import Bohr

from HamiltonIO.epw.epwparser import Epmat, read_crystal_fmt


@pytest.fixture
def srmo3_epmat():
    """Load SrMnO3 EPW data from test data."""
    # Find test data directory
    project_root = Path(__file__).parent.parent
    data_dir = project_root / "HamiltonIO" / "epw" / "test" / "up"

    # Read crystal structure
    crystal = read_crystal_fmt(str(data_dir / "crystal.fmt"))

    # Read EPW data
    ep = Epmat()
    ep.crystal = crystal
    ep.read(path=str(data_dir), prefix="SrMnO3", epmat_ncfile="epmat.nc")

    # Get cell in Angstrom
    cell = crystal.at.reshape(3, 3) * crystal.alat * Bohr

    return ep, crystal, cell


def test_epmat_loading(srmo3_epmat):
    """Test that EPW data loads correctly."""
    ep, crystal, cell = srmo3_epmat

    # Check basic properties
    assert ep.nwann > 0, "Should have Wannier functions"
    assert crystal.natom > 0, "Should have atoms"
    assert crystal.nmode == crystal.natom * 3, "Should have 3 modes per atom"

    # Check that data is loaded
    assert hasattr(ep, "epmatfile") or hasattr(ep, "epmat_wann"), (
        "Should have matrix data"
    )
    assert ep.nRk > 0, "Should have k-space R vectors"
    assert ep.nRg > 0, "Should have g-space R vectors"

    # Check cell shape
    assert cell.shape == (3, 3), "Cell should be 3x3 matrix"


def test_distance_resolved_couplings_Rk(srmo3_epmat):
    """Test Rk-based distance computation (WF to WF)."""
    ep, crystal, cell = srmo3_epmat
    imode = 0

    # Get distance-resolved couplings
    entries = ep.distance_resolved_couplings_Rk(imode=imode, cell=cell)

    # Should have couplings
    assert len(entries) > 0, "Should have non-zero coupling elements"

    # Check structure of each entry
    for entry in entries:
        assert "imode" in entry
        assert "Rg" in entry
        assert "Rk" in entry
        assert "i" in entry
        assert "j" in entry
        assert "distance" in entry
        assert "coupling" in entry

        # Check types
        assert entry["imode"] == imode
        assert len(entry["Rg"]) == 3
        assert len(entry["Rk"]) == 3
        assert isinstance(entry["i"], int)
        assert isinstance(entry["j"], int)
        assert isinstance(entry["distance"], float)
        assert entry["distance"] >= 0, "Distance should be non-negative"
        assert isinstance(entry["coupling"], complex)

    # Distances should be reasonable
    distances = np.array([e["distance"] for e in entries])
    assert distances.min() >= 0
    assert distances.max() < 100, "Distances should be reasonable (< 100 Å)"
    assert len(np.unique(distances)) > 1, "Should have multiple distinct distances"


def test_distance_resolved_couplings_Rg(srmo3_epmat):
    """Test Rg-based distance computation (WF to atom)."""
    ep, crystal, cell = srmo3_epmat
    imode = 0

    # Get distance-resolved couplings
    entries = ep.distance_resolved_couplings_Rg(imode=imode, cell=cell)

    # Should have couplings
    assert len(entries) > 0, "Should have non-zero coupling elements"

    # Check structure of each entry
    for entry in entries:
        assert "imode" in entry
        assert "atom_index" in entry
        assert "Rg" in entry
        assert "Rk" in entry
        assert "i" in entry
        assert "j" in entry
        assert "distance" in entry
        assert "coupling" in entry

        # Check types
        assert entry["imode"] == imode
        assert isinstance(entry["atom_index"], int)
        assert 0 <= entry["atom_index"] < crystal.natom, "Atom index should be valid"
        assert len(entry["Rg"]) == 3
        assert len(entry["Rk"]) == 3
        assert isinstance(entry["i"], int)
        assert isinstance(entry["j"], int)
        assert isinstance(entry["distance"], float)
        assert entry["distance"] >= 0, "Distance should be non-negative"
        assert isinstance(entry["coupling"], complex)

    # Distances should be reasonable
    distances = np.array([e["distance"] for e in entries])
    assert distances.min() >= 0
    assert distances.max() < 100, "Distances should be reasonable (< 100 Å)"


def test_mode_to_atom_mapping(srmo3_epmat):
    """Test that mode-to-atom mapping is correct."""
    ep, crystal, cell = srmo3_epmat

    # Test several modes
    for imode in range(
        min(9, crystal.nmode)
    ):  # Test first 9 modes (3 atoms × 3 directions)
        entries = ep.distance_resolved_couplings_Rg(imode=imode, cell=cell)

        # All entries for this mode should map to the same atom
        atom_indices = set(e["atom_index"] for e in entries)
        assert len(atom_indices) == 1, f"Mode {imode} should map to single atom"

        # Verify the mapping: atom_index = imode // 3
        expected_atom = imode // 3
        actual_atom = atom_indices.pop()
        assert actual_atom == expected_atom, (
            f"Mode {imode} should map to atom {expected_atom}, got {actual_atom}"
        )


def test_bin_couplings_by_distance(srmo3_epmat):
    """Test binning of coupling magnitudes by distance."""
    ep, crystal, cell = srmo3_epmat
    imode = 0

    entries = ep.distance_resolved_couplings_Rk(imode=imode, cell=cell)

    # Bin with default dr=0.1 Å
    centers, avg_mag = Epmat.bin_couplings_by_distance(entries, dr=0.1)

    # Should produce bins
    assert len(centers) > 0, "Should have bin centers"
    assert len(avg_mag) > 0, "Should have average magnitudes"
    assert len(centers) == len(avg_mag), "Should have same number of centers and values"

    # Bin centers should be monotonically increasing
    assert np.all(np.diff(centers) > 0), "Bin centers should be increasing"

    # Average magnitudes should be non-negative
    assert np.all(avg_mag >= 0), "Average magnitudes should be non-negative"

    # Try different bin width
    centers2, avg_mag2 = Epmat.bin_couplings_by_distance(entries, dr=0.5)
    assert len(centers2) < len(centers), "Larger dr should give fewer bins"


def test_bin_couplings_empty():
    """Test binning with empty input."""
    centers, avg_mag = Epmat.bin_couplings_by_distance([], dr=0.1)

    assert len(centers) == 0, "Empty input should give no bins"
    assert len(avg_mag) == 0, "Empty input should give no values"


def test_Rk_vs_Rg_consistency(srmo3_epmat):
    """Test that Rk and Rg methods give same coupling values but different distances."""
    ep, crystal, cell = srmo3_epmat
    imode = 0

    entries_Rk = ep.distance_resolved_couplings_Rk(imode=imode, cell=cell)
    entries_Rg = ep.distance_resolved_couplings_Rg(imode=imode, cell=cell)

    # Should have same number of entries
    assert len(entries_Rk) == len(entries_Rg), (
        "Rk and Rg methods should return same number of entries"
    )

    # Create lookup dictionaries
    Rk_lookup = {(e["Rg"], e["Rk"], e["i"], e["j"]): e for e in entries_Rk}
    Rg_lookup = {(e["Rg"], e["Rk"], e["i"], e["j"]): e for e in entries_Rg}

    # Check that coupling values match
    matches = 0
    dist_diffs = 0
    for key in Rk_lookup:
        if key in Rg_lookup:
            matches += 1
            # Coupling values should be identical
            assert np.isclose(
                Rk_lookup[key]["coupling"], Rg_lookup[key]["coupling"], rtol=1e-10
            ), f"Coupling mismatch for {key}"

            # Distances should generally be different (unless special case)
            if Rk_lookup[key]["distance"] != Rg_lookup[key]["distance"]:
                dist_diffs += 1

    assert matches > 0, "Should have matching entries"
    # Most distances should differ between Rk and Rg
    assert dist_diffs > matches * 0.5, "Most Rk and Rg distances should differ"


def test_onsite_Rk_distances(srmo3_epmat):
    """Test that Rk=0, i=j couplings have near-zero distance."""
    ep, crystal, cell = srmo3_epmat
    imode = 0

    entries = ep.distance_resolved_couplings_Rk(imode=imode, cell=cell)

    # Find onsite couplings (Rk=0, i=j)
    onsite = [e for e in entries if e["Rk"] == (0, 0, 0) and e["i"] == e["j"]]

    # Should have some onsite terms
    assert len(onsite) > 0, "Should have onsite diagonal elements"

    # Onsite distances should be very small (essentially zero)
    for entry in onsite:
        assert entry["distance"] < 1e-6, (
            f"Onsite Rk distance should be ~0, got {entry['distance']}"
        )


def test_wannier_centers_available(srmo3_epmat):
    """Test that Wannier centers are available in crystal structure."""
    ep, crystal, cell = srmo3_epmat

    # Should have Wannier centers
    assert hasattr(crystal, "w_centers"), "Crystal should have w_centers"
    assert crystal.w_centers is not None, "w_centers should not be None"

    # Check shape
    w_centers = np.asarray(crystal.w_centers).reshape(ep.nwann, 3)
    assert w_centers.shape == (ep.nwann, 3), (
        f"w_centers should be ({ep.nwann}, 3), got {w_centers.shape}"
    )


def test_atomic_positions_available(srmo3_epmat):
    """Test that atomic positions are available in crystal structure."""
    ep, crystal, cell = srmo3_epmat

    # Should have atomic positions
    assert hasattr(crystal, "tau"), "Crystal should have tau (atomic positions)"
    assert crystal.tau is not None, "tau should not be None"

    # Check shape
    tau = np.asarray(crystal.tau).reshape(crystal.natom, 3)
    assert tau.shape == (crystal.natom, 3), (
        f"tau should be ({crystal.natom}, 3), got {tau.shape}"
    )
