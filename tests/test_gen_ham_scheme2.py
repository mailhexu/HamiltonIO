"""Tests for scheme-2 (per-orbital-pair ``_wsvec.dat``) Fourier transform.

Story-028. Scheme 2 applies the Wannier90 ``use_ws_distance`` correction: for
each orbital pair ``(i, j)`` and lattice vector ``R``, the phase is averaged
over the WS image translations ``T`` and divided by
``ndegen(R) * N_T(i,j,R)``:

    H_ij(k) = sum_R sum_{T in Tlist(i,j,R)}
              exp(i 2pi k.(R+T)) / (ndegen(R) * N_T) * H_ij(R)

Source: Wannier90 ``plot.F90:335-340``, ``postw90_common.F90:723-727``.
Reference: research/2026-07-26-wannier90-wsvec-weights.md, ADR-001.
"""

from pathlib import Path

import numpy as np
import pytest

from HamiltonIO.wannier import WannierHam

# ---------------------------------------------------------------------------
# Independent brute-force reference (derived from the math, not the impl)
# ---------------------------------------------------------------------------


def _scheme2_reference(ham, k):
    """Compute H(k) with the per-pair WS formula, independent of gen_ham."""
    Hk = np.zeros((ham.nbasis, ham.nbasis), dtype=complex)
    shifts = ham.ws_shifts or {}
    for R, mat in ham.data.items():
        ndegen = ham.R_degens.get(R, 1)
        rblock = shifts.get(R, {})
        for i in range(ham.nbasis):
            for j in range(ham.nbasis):
                tlist = rblock.get((i, j), [(0, 0, 0)])
                n_t = len(tlist)
                acc = 0.0 + 0.0j
                for t in tlist:
                    rt = np.array(R, dtype=float) + np.array(t, dtype=float)
                    acc += np.exp(2.0j * np.pi * np.dot(k, rt))
                acc /= ndegen * n_t
                Hk[i, j] += mat[i, j] * acc
    return Hk


def _toy_ws_model():
    """2-orbital model with a mixed wsvec profile.

    R=(0,0,0): ndegen=1, ws default (N_T=1, T=0).
    R=(1,0,0): ndegen=2, ws per pair:
        (0,0): [(0,0,0), (0,1,0)]   N_T=2
        (0,1): [(0,0,0)]            N_T=1
        (1,0): [(0,0,0)]            N_T=1
        (1,1): [(0,0,0), (0,-1,0)]  N_T=2
    """
    nbasis = 2
    data = {
        (0, 0, 0): np.array([[1.0, 0.2], [0.2, 0.8]], dtype=complex),
        (1, 0, 0): np.array([[0.3, 0.1], [0.1, 0.25]], dtype=complex),
    }
    rdegens = {(0, 0, 0): 1, (1, 0, 0): 2}
    ws_shifts = {
        (1, 0, 0): {
            (0, 0): np.array([[0, 0, 0], [0, 1, 0]]),
            (0, 1): np.array([[0, 0, 0]]),
            (1, 0): np.array([[0, 0, 0]]),
            (1, 1): np.array([[0, 0, 0], [0, -1, 0]]),
        },
    }
    positions = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
    ham = WannierHam(
        nbasis=nbasis,
        data=data,
        positions=positions,
        R_degens=rdegens,
        use_ws=True,
        ws_shifts=ws_shifts,
    )
    return ham


# ---------------------------------------------------------------------------
# TEST-001: scheme-2 formula correctness
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "k",
    [
        np.array([0.0, 0.0, 0.0]),
        np.array([0.5, 0.0, 0.0]),
        np.array([0.0, 0.5, 0.0]),
        np.array([0.25, 0.75, 0.0]),
        np.array([0.3, -0.2, 0.5]),
    ],
)
def test_gen_ham_scheme2_matches_reference(k):
    """gen_ham(use_ws=True) must match the independent brute-force reference."""
    ham = _toy_ws_model()
    hk = ham.gen_ham(k, convention=2)
    expected = _scheme2_reference(ham, k)
    assert np.allclose(
        hk, expected, atol=1e-12
    ), f"scheme-2 mismatch at k={k}:\n got      {hk}\n expected {expected}"


# ---------------------------------------------------------------------------
# TEST-002: use_ws=False stays on scheme 1 (regression)
# ---------------------------------------------------------------------------


def test_gen_ham_scheme1_when_use_ws_false():
    """With use_ws=False, gen_ham ignores ws_shifts and uses scheme 1."""
    ham = _toy_ws_model()
    ham.use_ws = False
    k = np.array([0.25, 0.5, 0.0])
    hk = ham.gen_ham(k, convention=2)
    # Scheme-1 reference: phase = exp(ikR)/ndegen, no T sum.
    expected = np.zeros((ham.nbasis, ham.nbasis), dtype=complex)
    for R, mat in ham.data.items():
        phase = np.exp(2.0j * np.pi * np.dot(k, R)) / ham.R_degens[R]
        expected += mat * phase
    assert np.allclose(hk, expected, atol=1e-12)


# ---------------------------------------------------------------------------
# TEST-003: missing ws entry for an (R,i,j) falls back to scheme-1 element
# ---------------------------------------------------------------------------


def test_gen_ham_scheme2_missing_pair_defaults_to_scheme1():
    """If ws_shifts lacks an (R,i,j), that element uses T=[0], N_T=1."""
    nbasis = 1
    data = {(1, 0, 0): np.array([[0.4]], dtype=complex)}
    rdegens = {(1, 0, 0): 2}
    # ws_shifts present but missing the (0,0) pair for this R.
    ws_shifts = {(1, 0, 0): {}}
    ham = WannierHam(
        nbasis=nbasis,
        data=data,
        positions=np.array([[0.0, 0.0, 0.0]]),
        R_degens=rdegens,
        use_ws=True,
        ws_shifts=ws_shifts,
    )
    k = np.array([0.0, 0.0, 0.0])
    hk = ham.gen_ham(k, convention=2)
    # Fallback: exp(0)/(2*1) * 0.4 = 0.2.
    assert np.isclose(hk[0, 0], 0.2, atol=1e-12)


# ---------------------------------------------------------------------------
# TEST-004: default attributes on a plain WannierHam (no ws args)
# ---------------------------------------------------------------------------


def test_wannierham_defaults_no_ws():
    """A WannierHam built without ws args has use_ws=False, ws_shifts=None."""
    ham = WannierHam(
        nbasis=1,
        data={(0, 0, 0): np.array([[1.0]], dtype=complex)},
        positions=np.array([[0.0, 0.0, 0.0]]),
    )
    assert ham.use_ws is False
    assert ham.ws_shifts is None


# ---------------------------------------------------------------------------
# Auto-detect helper (tested directly; read_from_wannier_dir calls it)
# ---------------------------------------------------------------------------


from HamiltonIO.wannier.wannier_hamiltonian import _autodetect_wsvec  # noqa: E402

_TRUE_WSVEC = (
    "## written on 1Jan2026 at 00:00:00 with use_ws_distance=.true.\n"
    "   0    0    0    1    1\n   1\n   0    0    0\n"
)
_FALSE_WSVEC = (
    "## written on 1Jan2026 at 00:00:00 with use_ws_distance=.false.\n"
    "   0    0    0    1    1\n   1\n   0    0    0\n"
)


def test_autodetect_wsvec_present_true(tmp_path):
    (tmp_path / "wannier90_wsvec.dat").write_text(_TRUE_WSVEC)
    use_ws, shifts = _autodetect_wsvec(str(tmp_path), "wannier90")
    assert use_ws is True
    assert (0, 0, 0) in shifts


def test_autodetect_wsvec_present_false(tmp_path, recwarn):
    (tmp_path / "wannier90_wsvec.dat").write_text(_FALSE_WSVEC)
    use_ws, shifts = _autodetect_wsvec(str(tmp_path), "wannier90")
    assert use_ws is False
    assert shifts is None


def test_autodetect_wsvec_absent(tmp_path):
    use_ws, shifts = _autodetect_wsvec(str(tmp_path), "wannier90")
    assert use_ws is False
    assert shifts is None


# ---------------------------------------------------------------------------
# Integration: read_from_wannier_dir auto-detects wsvec on the CrI3 fixture
# ---------------------------------------------------------------------------

_CRI3_DIR = (
    Path(__file__).resolve().parents[2]
    / "TB2J"
    / "tests"
    / "data"
    / "inputs"
    / "3_CrI3_wannier_SOC"
    / "data"
    / "z"
)


@pytest.mark.skipif(not _CRI3_DIR.exists(), reason="CrI3 fixture not available")
def test_read_from_wannier_dir_autodetects_wsvec_cri3():
    from HamiltonIO.wannier.w90_parser import parse_atoms

    atoms = parse_atoms(str(_CRI3_DIR / "wannier90.win"))
    ham = WannierHam.read_from_wannier_dir(
        path=str(_CRI3_DIR), prefix="wannier90", atoms=atoms
    )
    assert ham.use_ws is True
    assert ham.ws_shifts is not None
    assert len(ham.ws_shifts) > 0
    # The model must be solvable (gen_ham runs end-to-end).
    hk = ham.gen_ham(np.array([0.0, 0.0, 0.0]), convention=2)
    assert hk.shape == (ham.nbasis, ham.nbasis)


# ---------------------------------------------------------------------------
# shift_position propagates ws state (re-keys to shifted R)
# ---------------------------------------------------------------------------


def test_shift_position_propagates_ws_state():
    """shift_position must carry use_ws/ws_shifts onto the new model."""
    ham = _toy_ws_model()
    # Force a non-trivial integer shift by giving different reference positions.
    ref = np.array([[0.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    # Set positions so round(pos - ref) is nonzero for orbital 1.
    ham._positions = np.array([[0.0, 0.0, 0.0], [0.0, 0.2, 0.0]])
    shifted = ham.shift_position(ref)
    assert shifted.use_ws is ham.use_ws
    assert shifted.ws_shifts is not None
    # gen_ham must run without KeyError on the shifted model.
    hk = shifted.gen_ham(np.array([0.1, 0.2, 0.3]), convention=2)
    assert hk.shape == (ham.nbasis, ham.nbasis)
