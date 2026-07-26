"""Tests for the Wannier90 ``_wsvec.dat`` parser (story-027).

Format (confirmed from Wannier90 ``ws_distance.F90:283-313`` ``ws_write_vec`` and
from local CrI3 / SrMnO3 samples):

    ## written on {date} at {time} with use_ws_distance={.true.|.false.}
    {Rx} {Ry} {Rz} {iw} {jw}      # block header (1-based iw, jw)
    {N_T}                          # number of image translations
    {Tx} {Ty} {Tz}                 # N_T lines of integer shifts
    ...                            # next block, order: irpt, iw, jw

The shift ``T`` is ``irdist_ws - irvec`` (the supercell translation only); the
full lattice vector for the Fourier phase is ``R + T``.
"""

from pathlib import Path

import numpy as np
import pytest

from HamiltonIO.wannier.wsvec_parser import parse_wsvec

# ---------------------------------------------------------------------------
# Synthetic files
# ---------------------------------------------------------------------------

_TRUE_HEADER = "## written on 1Jan2026 at 00:00:00 with use_ws_distance=.true.\n"
_FALSE_HEADER = "## written on 1Jan2026 at 00:00:00 with use_ws_distance=.false.\n"


def _block(R, iw, jw, shifts):
    """Render one wsvec block. shifts is a list of (Tx,Ty,Tz)."""
    lines = [f"  {R[0]:>3d} {R[1]:>3d} {R[2]:>3d} {iw:>3d} {jw:>3d}\n"]
    lines.append(f"  {len(shifts):>3d}\n")
    for s in shifts:
        lines.append(f"  {s[0]:>3d} {s[1]:>3d} {s[2]:>3d}\n")
    return "".join(lines)


def _synthetic_true_file():
    """2 wannier, 1 R-point, mixed N_T."""
    return (
        _TRUE_HEADER
        + _block((0, 0, 0), 1, 1, [(0, 0, 0)])
        + _block((0, 0, 0), 1, 2, [(0, 0, 0), (1, 0, 0)])
        + _block((0, 0, 0), 2, 1, [(0, 0, 0), (-1, 0, 0)])
        + _block((0, 0, 0), 2, 2, [(0, 0, 0)])
    )


def _synthetic_false_file():
    """2 wannier, 1 R-point, trivial (N_T=1, T=0)."""
    return (
        _FALSE_HEADER
        + _block((0, 0, 0), 1, 1, [(0, 0, 0)])
        + _block((0, 0, 0), 1, 2, [(0, 0, 0)])
        + _block((0, 0, 0), 2, 1, [(0, 0, 0)])
        + _block((0, 0, 0), 2, 2, [(0, 0, 0)])
    )


# ---------------------------------------------------------------------------
# TEST-001: header .true. + structure
# ---------------------------------------------------------------------------


def test_parse_wsvec_header_true(tmp_path):
    fname = tmp_path / "wsvec.dat"
    fname.write_text(_synthetic_true_file())
    result = parse_wsvec(str(fname))
    assert result["use_ws_distance"] is True
    shifts = result["shifts"]
    # One R key
    assert list(shifts.keys()) == [(0, 0, 0)]
    rblock = shifts[(0, 0, 0)]
    # 0-based (i,j) keys
    assert set(rblock.keys()) == {(0, 0), (0, 1), (1, 0), (1, 1)}
    # (0,0): single shift (0,0,0)
    assert np.array_equal(rblock[(0, 0)], np.array([[0, 0, 0]]))
    # (0,1): two shifts, order preserved
    assert rblock[(0, 1)].shape == (2, 3)
    assert np.array_equal(rblock[(0, 1)], np.array([[0, 0, 0], [1, 0, 0]]))


# ---------------------------------------------------------------------------
# TEST-002: header .false. -> trivial, warn
# ---------------------------------------------------------------------------


def test_parse_wsvec_header_false_trivial(tmp_path, recwarn):
    fname = tmp_path / "wsvec.dat"
    fname.write_text(_synthetic_false_file())
    result = parse_wsvec(str(fname))
    assert result["use_ws_distance"] is False
    # Every block must be N_T=1, T=(0,0,0)
    for R, block in result["shifts"].items():
        for (i, j), arr in block.items():
            assert arr.shape == (1, 3), f"{R},{i},{j} not trivial"
            assert np.array_equal(arr, np.array([[0, 0, 0]]))


# ---------------------------------------------------------------------------
# TEST-003: CrI3 real fixture anchor
# ---------------------------------------------------------------------------

_CRI3_WSVEC = (
    Path(__file__).resolve().parents[2]
    / "TB2J"
    / "tests"
    / "data"
    / "inputs"
    / "3_CrI3_wannier_SOC"
    / "data"
    / "z"
    / "wannier90_wsvec.dat"
)


@pytest.mark.skipif(not _CRI3_WSVEC.exists(), reason="CrI3 wsvec fixture not available")
def test_parse_wsvec_cri3_anchor():
    result = parse_wsvec(str(_CRI3_WSVEC))
    assert result["use_ws_distance"] is True
    # R=(1,4,0): orbital pair (0,0) has N_T=2 with shifts (0,-9,0) and (0,0,0).
    block = result["shifts"][(1, 4, 0)]
    arr = block[(0, 0)]
    assert arr.shape == (2, 3), f"expected N_T=2, got {arr.shape}"
    expected = np.array([[0, -9, 0], [0, 0, 0]])
    assert np.array_equal(arr, expected), f"shifts mismatch: {arr}"


# ---------------------------------------------------------------------------
# TEST-004: SrMnO3 real fixture anchor
# ---------------------------------------------------------------------------

_SRMO3_WSVEC = (
    Path(__file__).resolve().parents[2]
    / "superexchange_MFT"
    / "data"
    / "SrMnO3"
    / "Wannier_SrMnO3_NM"
    / "wannier90_wsvec.dat"
)


@pytest.mark.skipif(
    not _SRMO3_WSVEC.exists(), reason="SrMnO3 wsvec fixture not available"
)
def test_parse_wsvec_srmno3_anchor():
    result = parse_wsvec(str(_SRMO3_WSVEC))
    assert result["use_ws_distance"] is True
    # Corner R=(-2,-2,-2): pair (0,0) has N_T=8 (the 8 cube-corner images).
    block = result["shifts"][(-2, -2, -2)]
    arr = block[(0, 0)]
    assert arr.shape == (8, 3), f"expected N_T=8, got {arr.shape}"
    # The 8 corner shifts of a cube with mp_grid multiples of 4.
    assert np.array_equal(arr[0], np.array([0, 0, 0]))


# ---------------------------------------------------------------------------
# TEST-005: corrupt file -> raise
# ---------------------------------------------------------------------------


def test_parse_wsvec_corrupt_truncated(tmp_path):
    """A block header with no following N_T line must raise."""
    fname = tmp_path / "wsvec.dat"
    fname.write_text(_TRUE_HEADER + _block((0, 0, 0), 1, 1, [(0, 0, 0)]))  # ok
    # Append a header with no N_T
    content = fname.read_text() + "  0  0  0    1    2\n"
    fname.write_text(content)
    with pytest.raises((ValueError, IndexError, EOFError)):
        parse_wsvec(str(fname))


def test_parse_wsvec_bad_header_line(tmp_path):
    """A block header that doesn't parse as 5 ints must raise."""
    fname = tmp_path / "wsvec.dat"
    fname.write_text(_TRUE_HEADER + "garbage not a header\n")
    with pytest.raises((ValueError, IndexError)):
        parse_wsvec(str(fname))


# ---------------------------------------------------------------------------
# validate_ws_weights: sum-rule + wsvec N_T statistics (story-030)
# ---------------------------------------------------------------------------

from HamiltonIO.wannier.wsvec_parser import validate_ws_weights  # noqa: E402

_HRDAT_NDEGEN_1_2_2 = """written on test
1
3
   1    2    2
   0    0    0    1    1    1.000000    0.000000
   1    0    0    1    1    0.500000    0.000000
  -1    0    0    1    1    0.500000    0.000000
"""


def test_validate_ws_weights_sum_rule_pass(tmp_path):
    """ndegen=[1,2,2] -> sum(1/ndegen)=2 == mp_grid product 2*1*1."""
    hr = tmp_path / "toy_hr.dat"
    hr.write_text(_HRDAT_NDEGEN_1_2_2)
    res = validate_ws_weights(str(hr), mp_grid=(2, 1, 1))
    assert res["sum_rule_ok"] is True
    assert np.isclose(res["sum_rule_value"], 2.0)
    assert res["sum_rule_expected"] == 2
    assert res["ws_stats"] is None


def test_validate_ws_weights_sum_rule_fail(tmp_path):
    """Tampered ndegen -> sum rule mismatch detected."""
    bad = """written on test
1
3
   1    1    1
   0    0    0    1    1    1.000000    0.000000
   1    0    0    1    1    0.500000    0.000000
  -1    0    0    1    1    0.500000    0.000000
"""
    hr = tmp_path / "bad_hr.dat"
    hr.write_text(bad)
    res = validate_ws_weights(str(hr), mp_grid=(2, 1, 1))
    # sum(1/1+1/1+1/1)=3 != 2
    assert res["sum_rule_ok"] is False
    assert np.isclose(res["sum_rule_value"], 3.0)


def test_validate_ws_weights_with_wsvec_stats(tmp_path):
    """With a wsvec file, report per-pair N_T statistics."""
    hr = tmp_path / "toy_hr.dat"
    hr.write_text(_HRDAT_NDEGEN_1_2_2)
    ws = tmp_path / "toy_wsvec.dat"
    # 1 R, 1 orbital -> 1 block with N_T=3.
    ws.write_text(
        _TRUE_HEADER + _block((1, 0, 0), 1, 1, [(0, 0, 0), (0, 1, 0), (0, 0, 1)])
    )
    res = validate_ws_weights(str(hr), mp_grid=(2, 1, 1), wsvec_path=str(ws))
    assert res["ws_stats"] is not None
    assert res["ws_stats"]["use_ws_distance"] is True
    assert res["ws_stats"]["min"] == 3
    assert res["ws_stats"]["max"] == 3
    assert res["ws_stats"]["mean"] == 3.0


_CRI3_HR_VAL = (
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
_CRI3_WS_VAL = Path(_CRI3_HR_VAL).with_name("wannier90_wsvec.dat")


@pytest.mark.skipif(not _CRI3_HR_VAL.exists(), reason="CrI3 fixture not available")
def test_validate_ws_weights_cri3_real():
    """Real CrI3: sum rule == 81 (mp_grid 9 9 1); wsvec N_T stats are finite."""
    res = validate_ws_weights(str(_CRI3_HR_VAL), mp_grid=(9, 9, 1))
    assert res["sum_rule_ok"] is True
    assert np.isclose(res["sum_rule_value"], 81.0)
    if _CRI3_WS_VAL.exists():
        res2 = validate_ws_weights(
            str(_CRI3_HR_VAL), mp_grid=(9, 9, 1), wsvec_path=str(_CRI3_WS_VAL)
        )
        assert res2["ws_stats"]["use_ws_distance"] is True
        assert res2["ws_stats"]["min"] >= 1
        assert res2["ws_stats"]["max"] >= 2  # CrI3 has some N_T=2 blocks
