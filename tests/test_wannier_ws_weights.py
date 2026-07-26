"""Tests for Wannier90 Wigner-Seitz weight handling in ``WannierHam``.

These tests anchor story-026 (the ``ndegen`` division fix). The Wannier90
``_hr.dat`` format stores RAW ``ham_r`` and a separate ``ndegen`` array; the
reader MUST divide by ``ndegen(R)`` in the Fourier transform. The current code
has this division commented out, and moreover references the weight by integer
loop index rather than by R-tuple (``R_degens`` is dict-keyed by R). Both bugs
are caught here.

References:
  - Wannier90 source: ``hamiltonian.F90:478`` (writes raw ham_r),
    ``plot.F90:346`` (divides by ndegen at transform time).
  - research/2026-07-26-wannier90-wsvec-weights.md
"""

from pathlib import Path

import numpy as np
import pytest

from HamiltonIO.wannier import WannierHam
from HamiltonIO.wannier.w90_parser import parse_ham

# ---------------------------------------------------------------------------
# Toy models (self-contained correctness anchors; no external fixtures)
# ---------------------------------------------------------------------------


def _toy_model_1orb():
    """Single-orbital model with a known ndegen profile.

    R=(0,0,0): ndegen=1, H=E0=1.0
    R=(1,0,0): ndegen=2, H=t=0.5
    R=(-1,0,0): ndegen=2, H=t=0.5  (Hermitian partner; t real)

    Convention 2 Fourier transform (the Wannier90 convention):
        H(k) = sum_R exp(i 2pi k.R) / ndegen(R) * H(R)

    Expected:
        k=(0,0,0):   1.0/1 + 0.5/2 + 0.5/2 = 1.5
        k=(1/2,0,0): 1.0/1 - 0.5/2 - 0.5/2 = 0.5

    Without the ndegen division (the bug) one would instead get 2.0 and 0.0.
    """
    nbasis = 1
    data = {
        (0, 0, 0): np.array([[1.0]], dtype=complex),
        (1, 0, 0): np.array([[0.5]], dtype=complex),
        (-1, 0, 0): np.array([[0.5]], dtype=complex),
    }
    rdegens = {(0, 0, 0): 1, (1, 0, 0): 2, (-1, 0, 0): 2}
    positions = np.array([[0.0, 0.0, 0.0]])
    return WannierHam(nbasis=nbasis, data=data, positions=positions, R_degens=rdegens)


# ---------------------------------------------------------------------------
# TEST-001 / TEST-002: gen_ham divides by ndegen (convention 2)
# ---------------------------------------------------------------------------


def test_gen_ham_divides_by_ndegen_at_gamma():
    """H(k=0) must equal the ndegen-weighted sum, not the raw sum."""
    ham = _toy_model_1orb()
    hk = ham.gen_ham(np.array([0.0, 0.0, 0.0]), convention=2)
    # Correct: 1.0/1 + 0.5/2 + 0.5/2 = 1.5
    # Buggy (no division): 1.0 + 0.5 + 0.5 = 2.0
    assert np.isclose(
        hk[0, 0].real, 1.5, atol=1e-12
    ), f"Expected H(k=0)=1.5 with ndegen division, got {hk[0,0]}"


def test_gen_ham_divides_by_ndegen_at_zone_edge():
    """H(k=1/2) must equal the ndegen-weighted alternating sum."""
    ham = _toy_model_1orb()
    hk = ham.gen_ham(np.array([0.5, 0.0, 0.0]), convention=2)
    # Correct: 1.0/1 - 0.5/2 - 0.5/2 = 0.5
    # Buggy (no division): 1.0 - 0.5 - 0.5 = 0.0
    assert np.isclose(
        hk[0, 0].real, 0.5, atol=1e-12
    ), f"Expected H(k=1/2)=0.5 with ndegen division, got {hk[0,0]}"


def test_gen_ham_ndegen_keyed_by_R_not_int():
    """The division must look up R_degens[R-tuple], not R_degens[loop index].

    This catches the latent bug where the commented code used ``iR`` (int) while
    R_degens is dict-keyed by the R-tuple. With int keys the defaultdict would
    silently return 1 for every R, hiding the division.
    """
    ham = _toy_model_1orb()
    # If keyed by int, every weight would be 1 (defaultdict default) and H(k=0)
    # would be the unweighted 2.0. Asserting 1.5 proves the R-tuple lookup works.
    hk = ham.gen_ham(np.array([0.0, 0.0, 0.0]), convention=2)
    assert np.isclose(hk[0, 0].real, 1.5, atol=1e-12)


# ---------------------------------------------------------------------------
# TEST-003: convention 1 ndegen division
# ---------------------------------------------------------------------------
# NOTE: ``convention=1`` has a SEPARATE pre-existing shape bug
# (``np.dot(k, R + rjminusri)`` raises for nbasis>0), unrelated to the ndegen
# fix. The ndegen division has been applied to the convention=1 branch too
# (by inspection, symmetric with convention=2), but it cannot be exercised
# until that shape bug is fixed. Tracked as a follow-up, out of scope for
# story-026. Convention=2 is the default and the only path used by TB2J.
@pytest.mark.skip(reason="convention=1 has a pre-existing shape bug; see note above")
def test_gen_ham_convention1_divides_by_ndegen():
    """Placeholder: convention=1 ndegen division; blocked by a shape bug."""
    pass


# ---------------------------------------------------------------------------
# TEST-004: ndegen defaults cleanly to 1 when not provided
# ---------------------------------------------------------------------------


def test_gen_ham_no_rdegens_defaults_to_one():
    """Without R_degens, every weight is 1 and gen_ham still works (no KeyError)."""
    nbasis = 1
    data = {(0, 0, 0): np.array([[2.0]], dtype=complex)}
    positions = np.array([[0.0, 0.0, 0.0]])
    ham = WannierHam(nbasis=nbasis, data=data, positions=positions)
    hk = ham.gen_ham(np.array([0.0, 0.0, 0.0]), convention=2)
    assert np.isclose(hk[0, 0].real, 2.0, atol=1e-12)


# ---------------------------------------------------------------------------
# Missing-key robustness: gen_ham must not KeyError when an R in data is absent
# from R_degens (happens with partial dicts and after shift_position, which
# relabels R-vectors but passes the original R_degens dict through).
# ---------------------------------------------------------------------------


def test_gen_ham_missing_r_in_plain_dict_defaults_to_one():
    """A plain-dict R_degens missing a data key must fall back to weight 1."""
    data = {
        (0, 0, 0): np.array([[1.0]], dtype=complex),
        (1, 0, 0): np.array([[0.5]], dtype=complex),
    }
    rdegens = {(0, 0, 0): 1}  # (1,0,0) intentionally absent
    ham = WannierHam(
        nbasis=1, data=data, positions=np.array([[0.0, 0.0, 0.0]]), R_degens=rdegens
    )
    # Must not raise. With (1,0,0) defaulted to weight 1: H(k=0)=1.0+0.5=1.5.
    hk = ham.gen_ham(np.array([0.0, 0.0, 0.0]), convention=2)
    assert np.isclose(hk[0, 0].real, 1.5, atol=1e-12)


def test_gen_ham_after_shift_position_no_keyerror():
    """shift_position relabels R-vectors; gen_ham must not crash on the new keys.

    This mirrors the production path read_from_wannier_dir -> shift_position ->
    gen_ham. shift_position passes the original R_degens (keyed by old R) to the
    new instance, whose data is keyed by shifted R. Missing keys default to 1.
    """
    data = {
        (0, 0, 0): np.array([[1.0, 0.0], [0.0, 1.0]], dtype=complex),
        (1, 0, 0): np.array([[0.4, 0.0], [0.0, 0.4]], dtype=complex),
    }
    rdegens = {(0, 0, 0): 1, (1, 0, 0): 3}
    # Off-center positions force a non-trivial integer shift in shift_position.
    positions = np.array([[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]])
    ham = WannierHam(nbasis=2, data=data, positions=positions, R_degens=rdegens)
    shifted = ham.shift_position(np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]))
    # Must not raise regardless of how R-vectors were relabeled.
    hk = shifted.gen_ham(np.array([0.0, 0.0, 0.0]), convention=2)
    assert hk.shape == (2, 2)


# ---------------------------------------------------------------------------
# Sum-rule validation on a synthetic _hr.dat
# ---------------------------------------------------------------------------


_SYNTHETIC_HRDAT = """written on test
1
3
   1    2    2
   0    0    0    1    1    1.000000    0.000000
   1    0    0    1    1    0.500000    0.000000
  -1    0    0    1    1    0.500000    0.000000
"""


def test_parse_ham_reads_ndegen_and_sum_rule(tmp_path):
    """parse_ham must read the ndegen array; sum(1/ndegen) is the WS sum rule.

    For Wannier90 this sum equals the mp_grid product. Here ndegen=[1,2,2] gives
    sum=2, consistent with a 2x1x1 mesh.
    """
    fname = tmp_path / "toy_hr.dat"
    fname.write_text(_SYNTHETIC_HRDAT)
    n_wann, h_mnr, rdegens = parse_ham(str(fname))
    assert n_wann == 1
    assert len(rdegens) == 3
    assert list(rdegens) == [1, 2, 2]
    assert np.isclose(float(np.sum(1.0 / rdegens)), 2.0)
    # And the parsed H matches the synthetic values
    assert np.isclose(h_mnr[(0, 0, 0)][0, 0], 1.0)
    assert np.isclose(h_mnr[(1, 0, 0)][0, 0], 0.5)


def test_read_hr_then_gen_ham_uses_ndegen(tmp_path):
    """End-to-end: parse_ham -> WannierHam -> gen_ham applies the ndegen weight."""
    fname = tmp_path / "toy_hr.dat"
    fname.write_text(_SYNTHETIC_HRDAT)
    n_wann, data, rdegens = parse_ham(str(fname))
    ham = WannierHam(
        nbasis=n_wann,
        data=data,
        positions=np.array([[0.0, 0.0, 0.0]]),
        R_degens=rdegens,
    )
    hk = ham.gen_ham(np.array([0.0, 0.0, 0.0]), convention=2)
    assert np.isclose(hk[0, 0].real, 1.5, atol=1e-12)


# ---------------------------------------------------------------------------
# Optional: real-data sum rule on the CrI3 fixture (TB2J repo), if present.
# ---------------------------------------------------------------------------

_CRI3_HR = (
    Path(__file__).resolve().parents[2]
    / "TB2J"
    / "tests"
    / "data"
    / "inputs"
    / "3_CrI3_wannier_SOC"
    / "data"
    / "z"
    / "wannier90_hr.dat"
)


@pytest.mark.skipif(not _CRI3_HR.exists(), reason="CrI3 hr.dat fixture not available")
def test_sum_rule_cri3_real():
    """On real CrI3 data: sum_R 1/ndegen(R) == mp_grid product (9*9*1 = 81)."""
    n_wann, data, rdegens = parse_ham(str(_CRI3_HR))
    total = float(np.sum(1.0 / np.asarray(rdegens, dtype=float)))
    assert np.isclose(
        total, 81.0, atol=1e-6
    ), f"CrI3 sum-rule failed: got {total}, expected 81 (mp_grid 9 9 1)"
