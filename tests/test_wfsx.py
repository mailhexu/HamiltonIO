#!/usr/bin/env python3
"""
Unit tests for the SIESTA WFSX wavefunction parser.

Fixtures:
    si_prim_sz_6k.selected.WFSX
        Scalar (nspin=1) primitive 2-atom fcc Si cell, a=5.43 Ang, 8 bands,
        8 AOs, 6 k-points along Gamma-X: (0,0,0), (0.1,0,0), ... (0.5,0,0)
        in primitive fractional coordinates.
    si_sc_pso_gamma.selected.WFSX
        Spinor (nspin=4-family) 8-atom conventional Si cell at Gamma with
        Spin.Orbit on (scalar Si pseudo, so the SO splitting is ~0).
"""

from pathlib import Path

import numpy as np
import pytest
import sisl

from HamiltonIO.siesta.wfsx import SiestaWFSXParser, WFSXData

DATA_DIR = Path(__file__).parent / "data"
SCALAR_WFSX = DATA_DIR / "si_prim_sz_6k.selected.WFSX"
SPINOR_WFSX = DATA_DIR / "si_sc_pso_gamma.selected.WFSX"

# a=5.43 Ang fcc primitive cell, rows
FCC_CELL = 5.43 * np.array([[0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]])


@pytest.fixture(scope="module")
def scalar_data() -> WFSXData:
    return SiestaWFSXParser(str(SCALAR_WFSX), cell=FCC_CELL).read()


@pytest.fixture(scope="module")
def spinor_data() -> WFSXData:
    return SiestaWFSXParser(str(SPINOR_WFSX), cell=np.eye(3)).read()


def test_scalar_sizes_and_kpoints(scalar_data):
    d = scalar_data
    assert d.kpoints.shape == (6, 3)
    assert d.kpoints_cart.shape == (6, 3)
    assert d.eigenvalues.shape == (6, 8)
    assert d.coefficients.shape == (6, 8, 8)
    assert np.iscomplexobj(d.coefficients)
    assert d.is_spinor is False
    assert d.norb == 8
    np.testing.assert_allclose(d.kpoints[0], 0.0, atol=1e-12)
    np.testing.assert_allclose(d.kpoints[1], (0.1, 0.0, 0.0), atol=1e-8)


def test_scalar_eigenvalues(scalar_data):
    # SIESTA Fermi-shifted values as stored (eV). Reference values are quoted
    # to 3 decimals (third element differs by ~2e-4 from -4.367), so the
    # rounded trio is checked at 1e-3; the tight check below pins 1e-6.
    assert np.allclose(
        scalar_data.eigenvalues[0][:3],
        [-16.634, -4.367, -4.367],
        atol=1e-3,
    )
    # tight values of the stored doubles
    np.testing.assert_allclose(
        scalar_data.eigenvalues[0][:3],
        [-16.633740, -4.366985, -4.366803],
        atol=1e-6,
    )


def test_spinor(spinor_data):
    d = spinor_data
    # 8 atoms x 4 PAO = no_u = 32 orbitals; each state stores the two spinor
    # components of every orbital => norb = 2 * no_u = 64 columns, and the
    # non-collinear Hamiltonian is 64x64 => 64 bands at Gamma.
    assert d.is_spinor is True
    assert d.norb == 64
    assert d.kpoints.shape == (1, 3)
    assert d.eigenvalues.shape == (1, 64)
    assert d.coefficients.shape == (1, 64, 64)
    assert np.iscomplexobj(d.coefficients)
    np.testing.assert_allclose(d.kpoints_cart[0], 0.0, atol=1e-12)


def test_spinor_interleaved_components(spinor_data):
    # With SOC ~ 0 the Gamma eigenstates come in Kramers pairs, each almost
    # pure in one spin channel. SIESTA stores the two spinor components of
    # orbital i in consecutive slots (2i, 2i+1), so a pure-spin band has
    # weight only on a single parity of the column index.
    d = spinor_data
    c = d.coefficients[0]
    weight_even = (np.abs(c) ** 2)[:, 0::2].sum(axis=1)
    weight_odd = (np.abs(c) ** 2)[:, 1::2].sum(axis=1)
    assert np.all((weight_even < 1e-8) | (weight_odd < 1e-8))
    assert np.all(weight_even + weight_odd > 0.1)


def test_kpoints_cart_roundtrip():
    # kpoints_cart must reproduce sisl's stored per-k Cartesian k (1/Ang, 2pi)
    for path in (SCALAR_WFSX, SPINOR_WFSX):
        d = SiestaWFSXParser(str(path)).read()
        sile = sisl.get_sile(str(path))
        stored = np.array([st.info["k"] for st in sile.yield_eigenstate()])
        np.testing.assert_allclose(d.kpoints_cart, stored, atol=1e-12)


def test_no_cell_gives_nan_kpoints():
    d = SiestaWFSXParser(str(SCALAR_WFSX)).read()
    assert d.kpoints.shape == (6, 3)
    assert np.all(np.isnan(d.kpoints))
    # everything else is unaffected
    assert d.kpoints_cart.shape == (6, 3)
    assert d.coefficients.shape == (6, 8, 8)
