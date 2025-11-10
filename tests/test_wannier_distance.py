"""Tests for distance-resolved hopping analysis in Wannier90 Hamiltonian.

Tests cover:
    - Correctness of iter_unique_hoppings (symmetry reduction)
    - Distance computation from Wannier centers
    - Binning of hopping magnitudes by distance
"""

import os
from pathlib import Path

import numpy as np
import pytest
from ase.io import read

from HamiltonIO.wannier import WannierHam


@pytest.fixture
def srmo3_ham():
    """Load SrMnO3 Wannier TB model from test data."""
    test_dir = Path(__file__).parent / "data" / "SrMnO3_wannier"
    pwi_file = test_dir / "SrMnO3.pwi"
    atoms = read(pwi_file, format="espresso-in")
    ham = WannierHam.read_from_wannier_dir(test_dir, "SrMnO3", atoms=atoms)
    return ham, atoms


def test_iter_unique_hoppings(srmo3_ham):
    """Test that iter_unique_hoppings yields unique hopping terms."""
    ham, _ = srmo3_ham

    # Collect all unique hoppings
    unique_hops = list(ham.iter_unique_hoppings())

    # Should have at least some hoppings
    assert len(unique_hops) > 0

    # Check structure: (R, i, j, value)
    for R, i, j, val in unique_hops:
        assert len(R) == 3, "R should be 3D"
        assert isinstance(i, (int, np.integer)), "i should be integer"
        assert isinstance(j, (int, np.integer)), "j should be integer"
        assert isinstance(val, (complex, np.complexfloating)), "value should be complex"

    # Check symmetry convention: only positive R (first non-zero component >= 0)
    for R, i, j, val in unique_hops:
        nz = np.nonzero(R)[0]
        if len(nz) > 0:
            assert R[nz[0]] >= 0, (
                f"First non-zero R component should be >= 0, got R={R}"
            )

    # Check R=0 convention: only upper triangle (i <= j)
    r0_hops = [(i, j) for R, i, j, _ in unique_hops if np.all(np.array(R) == 0)]
    for i, j in r0_hops:
        assert i <= j, f"For R=0, should have i <= j, got i={i}, j={j}"


def test_distance_resolved_hoppings(srmo3_ham):
    """Test distance computation for hopping matrix elements."""
    ham, atoms = srmo3_ham
    cell = atoms.get_cell()

    # Get distance-resolved hoppings
    entries = ham.distance_resolved_hoppings(cell)

    # Should have hoppings
    assert len(entries) > 0

    # Check structure of each entry
    for entry in entries:
        assert "R" in entry
        assert "i" in entry
        assert "j" in entry
        assert "distance" in entry
        assert "hopping" in entry

        # Check types
        assert len(entry["R"]) == 3
        assert isinstance(entry["i"], int)
        assert isinstance(entry["j"], int)
        assert isinstance(entry["distance"], float)
        assert entry["distance"] >= 0, "Distance should be non-negative"
        assert isinstance(entry["hopping"], complex)

    # Distances should be reasonable (not all zero, not all huge)
    distances = np.array([e["distance"] for e in entries])
    assert distances.min() >= 0
    assert distances.max() < 100, "Distances should be reasonable (< 100 Å)"
    assert len(np.unique(distances)) > 1, "Should have multiple distinct distances"


def test_distance_resolved_hoppings_onsite(srmo3_ham):
    """Test that R=0, i=j hoppings have near-zero distance."""
    ham, atoms = srmo3_ham
    cell = atoms.get_cell()

    entries = ham.distance_resolved_hoppings(cell)

    # Find onsite hoppings (R=0, i=j)
    onsite = [e for e in entries if e["R"] == (0, 0, 0) and e["i"] == e["j"]]

    # Should have some onsite terms
    assert len(onsite) > 0, "Should have onsite diagonal elements"

    # Onsite distances should be very small (essentially zero)
    for entry in onsite:
        assert entry["distance"] < 1e-6, (
            f"Onsite distance should be ~0, got {entry['distance']}"
        )


def test_bin_hoppings_by_distance(srmo3_ham):
    """Test binning of hopping magnitudes by distance."""
    ham, atoms = srmo3_ham
    cell = atoms.get_cell()

    entries = ham.distance_resolved_hoppings(cell)

    # Bin with default dr=0.1 Å
    centers, avg_mag = WannierHam.bin_hoppings_by_distance(entries, dr=0.1)

    # Should produce bins
    assert len(centers) > 0
    assert len(avg_mag) > 0
    assert len(centers) == len(avg_mag)

    # Bin centers should be monotonically increasing
    assert np.all(np.diff(centers) > 0), "Bin centers should be increasing"

    # Average magnitudes should be non-negative
    assert np.all(avg_mag >= 0), "Average magnitudes should be non-negative"

    # Try different bin width
    centers2, avg_mag2 = WannierHam.bin_hoppings_by_distance(entries, dr=0.5)
    assert len(centers2) < len(centers), "Larger dr should give fewer bins"


def test_bin_hoppings_empty():
    """Test binning with empty input."""
    centers, avg_mag = WannierHam.bin_hoppings_by_distance([], dr=0.1)

    assert len(centers) == 0
    assert len(avg_mag) == 0


def test_distance_mapping_consistency(srmo3_ham):
    """Test that Hamiltonian indices map consistently to distances."""
    ham, atoms = srmo3_ham
    cell = atoms.get_cell()

    entries = ham.distance_resolved_hoppings(cell)

    # Group by (R, i, j) to check uniqueness
    seen = set()
    for entry in entries:
        key = (entry["R"], entry["i"], entry["j"])
        assert key not in seen, (
            f"Duplicate hopping found: R={entry['R']}, i={entry['i']}, j={entry['j']}"
        )
        seen.add(key)

    # Verify distances are computed from Wannier positions
    positions_frac = np.asarray(ham.positions)

    # Check a few entries manually
    for entry in entries[:10]:  # Check first 10 for speed
        R = np.array(entry["R"])
        i = entry["i"]
        j = entry["j"]

        # Compute distance manually
        d_frac = positions_frac[j] + R - positions_frac[i]
        d_cart = d_frac @ cell
        expected_dist = float(np.linalg.norm(d_cart))

        assert np.isclose(entry["distance"], expected_dist, rtol=1e-6), (
            f"Distance mismatch for R={R}, i={i}, j={j}"
        )


def test_wannier_centers_exist(srmo3_ham):
    """Test that Wannier centers are available."""
    ham, _ = srmo3_ham

    # Should have positions attribute
    assert hasattr(ham, "positions")
    assert ham.positions is not None
    assert len(ham.positions) == ham.nbasis

    # Positions should be fractional coordinates in [0, 1)
    positions = np.asarray(ham.positions)
    assert positions.shape == (ham.nbasis, 3)
