"""OpenMX ``.scfout`` parser tests on committed Si diamond fixtures.

The fixtures (``tests/data/openmx_si_{prim,sc,sc_p}.*``) were produced
with OpenMX 3.9 (GGA-PBE, Si7.0-s2p2d1 / P7.0-s2p2d1, SpinPolarization
off): a 2-atom fcc primitive cell on a 4x4x4 k-grid and 8-atom
conventional cells (pristine Si8 and substitutional Si7P) on 2x2x2
grids. Physics oracles: OpenMX's own eigenvalue blocks in the ``.out``
files and direct real-space diagonalization.
"""

import os
import re

import numpy as np
import pytest
from ase.units import Ha
from scipy.linalg import eigh

DATA = os.path.join(os.path.dirname(__file__), "data")

PRIM_SCFOUT = os.path.join(DATA, "openmx_si_prim.scfout")
PRIM_OUT = os.path.join(DATA, "openmx_si_prim.out")
SC_SCFOUT = os.path.join(DATA, "openmx_si_sc.scfout")
SC_P_SCFOUT = os.path.join(DATA, "openmx_si_sc_p.scfout")

B_DIAMOND = np.array([[-1, 1, 1], [1, -1, 1], [1, 1, -1]])


def _parse(path):
    from HamiltonIO.openmx import OpenmxParser

    return OpenmxParser(path)


def _out_mesh_eigenvalues(path, max_k=2):
    """``[(k, evals_eV), ...]`` from an OpenMX ``.out`` kloop dump."""
    lines = open(path).read().splitlines()
    res, i = [], 0
    while i < len(lines) - 3 and len(res) < max_k:
        if lines[i].strip().startswith("k1="):
            k = (
                float(lines[i].split("k1=")[1].split()[0]),
                float(lines[i].split("k2=")[1].split()[0]),
                float(lines[i].split("k3=")[1].split()[0]),
            )
            vals, j = [], i + 2
            while j < len(lines):
                m = re.match(r"\s*(\d+)\s+([-\d.Ee+]+)\s+([-\d.Ee+]+)\s*$", lines[j])
                if not m:
                    break
                vals.append(float(m.group(2)))
                j += 1
            res.append((k, np.array(vals) * Ha))
            i = j
        i += 1
    return res


def test_prim_shapes_and_structure():
    parser = _parse(PRIM_SCFOUT)
    data = parser.data
    assert data.natom == 2
    assert data.spin_p_switch == 0
    assert data.norb == 26  # Si7.0-s2p2d1: 2s + 6p + 5d per atom
    assert data.symbols() == ["Si", "Si"]
    assert len(data.Rlist) == 343  # symmetric {-3..3}^3 image cube
    assert (0, 0, 0) in {tuple(int(v) for v in R) for R in data.Rlist}
    atoms = parser.atoms
    assert atoms.get_chemical_formula() == "Si2"
    assert np.allclose(atoms.cell.lengths(), [2.715 * np.sqrt(2)] * 3, atol=1e-3)


def test_sc_and_dopant_shapes():
    sc = _parse(SC_SCFOUT)
    assert sc.data.natom == 8
    assert sc.data.norb == 104
    assert sc.atoms.get_chemical_formula() == "Si8"
    assert np.allclose(sc.atoms.cell.lengths(), [5.43] * 3, atol=1e-3)
    dop = _parse(SC_P_SCFOUT)
    assert dop.atoms.get_chemical_symbols() == ["Si"] * 4 + ["P"] + ["Si"] * 3
    assert dop.data.norb == 104  # P carries the same 13-orbital basis
    assert dop.data.valence_electrons == pytest.approx(33.0)  # 7*4 + 5
    # far image shells may be empty; the R=0 self-overlap must be exact
    assert dop.data.S.shape == (len(dop.data.Rlist), 104, 104)
    assert np.abs(dop.data.S[0].diagonal() - 1.0).max() < 1e-4


def test_units_and_fermi_level():
    parser = _parse(PRIM_SCFOUT)
    assert parser.efermi == pytest.approx(parser.data.chem_p * Ha)
    assert parser.efermi == pytest.approx(-0.1029 * Ha, abs=1e-3)
    # cell vectors are converted from Bohr to Angstrom
    sc = _parse(SC_SCFOUT)
    a_Bohr = 5.43 / 0.529177210903
    assert np.allclose(sc.data.tv, a_Bohr * np.eye(3), atol=1e-3)


def test_eigenvalues_match_openmx_output():
    """H(k) diagonalization reproduces OpenMX's own k-mesh eigenvalues."""
    parser = _parse(PRIM_SCFOUT)
    model = parser.get_model()
    for k, evals_out in _out_mesh_eigenvalues(PRIM_OUT, max_k=2):
        H, S, E, V = model.HS_and_eigen(np.atleast_2d(np.asarray(k)))
        assert np.abs(np.sort(E[0]) - evals_out).max() < 2e-3  # eV


def test_gamma_matches_direct_real_space_diagonalization():
    model = _parse(SC_SCFOUT).get_model()
    H_sum = sum(np.asarray(block) for block in model.HR)
    S_sum = sum(np.asarray(block) for block in model.SR)
    e_direct = eigh(H_sum, S_sum, eigvals_only=True)
    H, S, E, V = model.HS_and_eigen(np.atleast_2d(np.zeros(3)))
    assert np.abs(E[0] - e_direct).max() < 1e-8


def test_translation_phase_convention():
    """Blocks obey H(-R) = H(R)^dag and H(k) stays Hermitian off-grid."""
    model = _parse(SC_SCFOUT).get_model()
    Rdict = {tuple(int(v) for v in R): i for i, R in enumerate(model.Rlist)}
    dev = 0.0
    for T, iR in list(Rdict.items())[:40]:
        mT = tuple(-v for v in T)
        if mT in Rdict:
            dev = max(
                dev,
                np.abs(
                    np.asarray(model.HR)[iR] - np.asarray(model.HR)[Rdict[mT]].conj().T
                ).max(),
            )
    assert dev < 1e-7
    H, S, E, V = model.HS_and_eigen(np.atleast_2d([0.13, 0.27, 0.41]))
    assert np.abs(H[0] - H[0].conj().T).max() < 1e-7
    assert np.abs(S[0] - S[0].conj().T).max() < 1e-7


def test_spin_api_on_unpolarized_data():
    from HamiltonIO.openmx import OpenmxWrapper

    parser = _parse(SC_SCFOUT)
    model = parser.get_model()
    assert isinstance(model, OpenmxWrapper)
    assert model.nbasis == 104
    assert model.nspin == 1
    # unpolarized data has no spin channels to select
    with pytest.raises(ValueError):
        parser.get_model(spin=0)


def test_efermi_matches_output_chemp():
    """efermi is ChemP converted Hartree -> eV (sc fixture: -0.09902 Ha)."""
    model = _parse(SC_SCFOUT).get_model()
    assert model.efermi == pytest.approx(-0.09902050121809 * Ha, abs=1e-6)
    prim = _parse(PRIM_SCFOUT)
    assert prim.efermi == pytest.approx(-0.10287990806711 * Ha, abs=1e-6)
    dop = _parse(SC_P_SCFOUT)
    assert dop.efermi == pytest.approx(-0.08645988543752 * Ha, abs=1e-6)


def test_get_models_unpolarized_returns_single():
    from HamiltonIO.openmx import OpenmxWrapper

    models = _parse(SC_P_SCFOUT).get_models()
    assert isinstance(models, OpenmxWrapper)
    assert models.nbasis == 104
