"""Tests for EPW Wigner-Seitz metadata preservation in Epmat (story 1).

EPW pre-bakes WS weights into .epmatwp, so the per-pair ndegen in wigner.fmt is
informational/validation metadata (not a physics input). These tests check that
Epmat.read_Rvectors preserves the full per-pair ndegen arrays and surfaces the
use_ws status, instead of flattening to 1D and discarding the info.

Reference: specs/research/2026-07-26-epw-wigner-seitz.md
"""

import numpy as np

from HamiltonIO.epw.epwparser import Epmat
from HamiltonIO.epw.wigner import WignerData


def _make_wigner_data(nrr=5, dims=1, dims2=1):
    """Build a small WignerData with known ndegen values."""
    irvec = np.zeros((nrr, 3), dtype=int)
    for i in range(nrr):
        irvec[i] = [i - 2, 0, 0]
    wslen = np.ones(nrr, dtype=float)
    ndegen_k = np.ones((nrr, dims, dims), dtype=int)
    ndegen_q = np.ones((nrr, dims2, dims2), dtype=int)
    ndegen_g = np.ones((nrr, dims, dims2), dtype=int)
    # Set a nontrivial pattern so we can distinguish per-pair values
    for ir in range(nrr):
        if dims > 1:
            ndegen_k[ir] = (ir + 1) % 4 + 1
        if dims2 > 1:
            ndegen_q[ir] = (ir + 2) % 3 + 1
        ndegen_g[ir] = (ir + 1) % 5 + 1
    return WignerData(
        dims=dims,
        dims2=dims2,
        nrr_k=nrr,
        irvec_k=irvec.copy(),
        ndegen_k=ndegen_k,
        wslen_k=wslen.copy(),
        nrr_q=nrr,
        irvec_q=irvec.copy(),
        ndegen_q=ndegen_q,
        wslen_q=wslen.copy(),
        nrr_g=nrr,
        irvec_g=irvec.copy(),
        ndegen_g=ndegen_g,
        wslen_g=wslen.copy(),
    )


def _write_wigner_fmt(tmp_path, wd, name="wigner.fmt"):
    fname = tmp_path / name
    wd.to_file(str(fname))
    return str(tmp_path), name


# ---------------------------------------------------------------------------
# use_ws=.false. (dims=1, dims2=1): the common case (k444q444 SrMnO3)
# ---------------------------------------------------------------------------


def test_epmat_use_ws_false_preserves_full_and_1d(tmp_path):
    wd = _make_wigner_data(nrr=5, dims=1, dims2=1)
    path, fname = _write_wigner_fmt(tmp_path, wd)
    ep = Epmat()
    ep.read_Rvectors(path, fname=fname)
    assert ep.use_ws is False
    assert ep.dims == 1 and ep.dims2 == 1
    # Full per-pair arrays have the (nR, dims, dims) shape.
    assert ep.ndegen_k_full.shape == (5, 1, 1)
    assert ep.ndegen_q_full.shape == (5, 1, 1)
    assert ep.ndegen_g_full.shape == (5, 1, 1)
    # 1D backwards-compat arrays are still present.
    assert ep.ndegen_k.shape == (5,)
    assert np.array_equal(ep.ndegen_k, ep.ndegen_k_full[:, 0, 0])


# ---------------------------------------------------------------------------
# use_ws=.true. (dims=nbndsub, dims2=nat): per-pair ndegen preserved
# ---------------------------------------------------------------------------


def test_epmat_use_ws_true_preserves_per_pair_ndegen(tmp_path):
    wd = _make_wigner_data(nrr=5, dims=2, dims2=3)
    path, fname = _write_wigner_fmt(tmp_path, wd)
    ep = Epmat()
    ep.read_Rvectors(path, fname=fname)
    assert ep.use_ws is True
    assert ep.dims == 2 and ep.dims2 == 3
    # Per-pair shapes preserved.
    assert ep.ndegen_k_full.shape == (5, 2, 2)
    assert ep.ndegen_q_full.shape == (5, 3, 3)
    assert ep.ndegen_g_full.shape == (5, 2, 3)
    # The full arrays match what we wrote (per-pair, not flattened).
    np.testing.assert_array_equal(ep.ndegen_k_full, wd.ndegen_k)
    np.testing.assert_array_equal(ep.ndegen_g_full, wd.ndegen_g)


def test_epmat_use_ws_true_1d_is_representative(tmp_path):
    """The 1D ndegen_* (backwards compat) is the (0,0) pair slice of _full."""
    wd = _make_wigner_data(nrr=5, dims=2, dims2=3)
    path, fname = _write_wigner_fmt(tmp_path, wd)
    ep = Epmat()
    ep.read_Rvectors(path, fname=fname)
    assert ep.ndegen_k.shape == (5,)
    np.testing.assert_array_equal(ep.ndegen_k, ep.ndegen_k_full[:, 0, 0])


def test_epmat_roundtrip_wigner_fmt(tmp_path):
    """Writing then reading wigner.fmt preserves the per-pair ndegen."""
    wd = _make_wigner_data(nrr=4, dims=2, dims2=2)
    path, fname = _write_wigner_fmt(tmp_path, wd)
    ep = Epmat()
    ep.read_Rvectors(path, fname=fname)
    # irvec and ndegen survive the round-trip
    np.testing.assert_array_equal(ep.Rk, wd.irvec_k)
    np.testing.assert_array_equal(ep.ndegen_g_full, wd.ndegen_g)


# ---------------------------------------------------------------------------
# validate_epw_ws_weights (story 2)
# ---------------------------------------------------------------------------

from HamiltonIO.epw.epwparser import validate_epw_ws_weights  # noqa: E402


def test_validate_use_ws_false_global_sum_rule(tmp_path):
    """use_ws=False: global sum_R 1/ndegen == mp_grid product on the 1D array."""
    # ndegen all = 1 over 8 R-points -> sum = 8.
    wd = _make_wigner_data(nrr=8, dims=1, dims2=1)
    path, fname = _write_wigner_fmt(tmp_path, wd)
    ep = Epmat()
    ep.read_Rvectors(path, fname=fname)
    res = validate_epw_ws_weights(ep, mp_grid=(2, 2, 2))
    assert res["use_ws"] is False
    assert res["k"]["sum_rule_ok"] is True
    assert np.isclose(res["k"]["sum_rule_value"], 8.0)
    assert res["k"]["sum_rule_expected"] == 8


def test_validate_reports_ndegen_stats(tmp_path):
    """The helper reports per-channel ndegen stats (min/max/mean, zero fraction)."""
    wd = _make_wigner_data(nrr=5, dims=2, dims2=2)
    # Zero out one pair at R=0 to create a zeroed entry
    wd.ndegen_g[0, 0, 0] = 0
    path, fname = _write_wigner_fmt(tmp_path, wd)
    ep = Epmat()
    ep.read_Rvectors(path, fname=fname)
    res = validate_epw_ws_weights(ep)
    assert res["use_ws"] is True
    for ch in ("k", "q", "g"):
        assert "min" in res[ch] and "max" in res[ch] and "mean" in res[ch]
        assert "zero_fraction" in res[ch]
    # g has at least one zero entry
    assert res["g"]["zero_fraction"] > 0


def test_validate_detects_bad_sum_rule(tmp_path):
    """A tampered ndegen (all 1s but claiming mp_grid 3x3x3) fails the sum rule."""
    wd = _make_wigner_data(nrr=5, dims=1, dims2=1)
    path, fname = _write_wigner_fmt(tmp_path, wd)
    ep = Epmat()
    ep.read_Rvectors(path, fname=fname)
    # sum(1/1 * 5) = 5 != 27
    res = validate_epw_ws_weights(ep, mp_grid=(3, 3, 3))
    assert res["k"]["sum_rule_ok"] is False
