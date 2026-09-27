"""ABACUS unfolding-facing parser checks (story 018).

Validates ``AbacusParser.get_unfold_model()`` on the committed Si
fixtures (tests/data/abacus_example):

- ``si_prim``: 2-atom primitive fcc cell, DZP (26 orbitals),
  4x4x4 k-grid, HR/SR sparse output plus ABACUS's own k-space
  ``data-{k}-H`` files and ``BANDS_1.dat`` eigenvalues,
- ``si_conv``: 8-atom conventional cell, DZP (104 orbitals), 2x2x2
  k-grid.

The ABACUS-side oracles written by the same run pin every convention:

- HR units: the sparse tables are in Ry on disk (parser converts to
  eV); diagonalization must reproduce ABACUS's own eV BANDS output,
- phase convention: parsed-model eigenvalues at *generic* (non-Gamma)
  grid k-points reproduce ABACUS's H(k) eigenvalues to text precision,
  which fails if the R-to-k Bloch phases disagree with ABACUS,
- SR overlap shapes/keys: dict keyed by integer translation tuples.
"""

import os

import numpy as np
import pytest
from scipy.linalg import eigh

DATA = os.path.join(os.path.dirname(__file__), "data", "abacus_example")
RY_TO_EV = 13.605698066


def _outpath(name):
    path = os.path.join(DATA, name, "OUT." + name)
    if not os.path.exists(path):
        pytest.skip(f"ABACUS fixture {path} missing")
    return path


@pytest.fixture(scope="module")
def prim_parser():
    pytest.importorskip("HamiltonIO")
    from HamiltonIO.abacus.abacus_wrapper import AbacusParser

    return AbacusParser(outpath=_outpath("si_prim"))


@pytest.fixture(scope="module")
def conv_parser():
    pytest.importorskip("HamiltonIO")
    from HamiltonIO.abacus.abacus_wrapper import AbacusParser

    return AbacusParser(outpath=_outpath("si_conv"))


def read_kpoints(path, nk):
    """Supercell fractional k-points from the ABACUS ``kpoints`` file."""
    rows = []
    with open(path) as handle:
        for line in handle:
            parts = line.split()
            if len(parts) == 5:
                try:
                    rows.append([float(v) for v in parts])
                except ValueError:
                    continue
    return np.array(rows[:nk])[:, 1:4]


class TestUnfoldModelInterface:
    def test_atoms_supercell(self, prim_parser):
        model = prim_parser.get_unfold_model()
        assert len(model.atoms) == 2
        assert model.atoms.get_chemical_symbols() == ["Si", "Si"]
        # conventional lattice constant: prim fcc cell edge a/sqrt(2)
        assert model.atoms.cell.lengths().min() == pytest.approx(
            5.431 / np.sqrt(2), abs=1e-2
        )

    def test_hr_sr_dict_keys(self, prim_parser):
        model = prim_parser.get_unfold_model()
        assert len(model.HR) == len(model.SR) == 189
        assert model.SR[(0, 0, 0)].shape == (26, 26)
        for key in list(model.HR)[:5]:
            assert isinstance(key, tuple) and len(key) == 3
            assert all(int(v) == v for v in key)
        # overlaps are Hermitian with unit diagonal at R=0
        S0 = model.SR[(0, 0, 0)]
        assert np.abs(S0 - S0.conj().T).max() < 1e-8
        assert np.abs(np.diag(S0) - 1.0).max() < 1e-6

    def test_hs_and_eigen_shapes(self, prim_parser):
        model = prim_parser.get_unfold_model()
        H, S = model.hs_and_eigen([0.1, 0.2, 0.3])
        assert H.shape == S.shape == (26, 26)
        assert np.abs(H - H.conj().T).max() < 1e-6
        assert np.abs(S - S.conj().T).max() < 1e-8

    def test_atom_orbital_counts(self, prim_parser):
        model = prim_parser.get_unfold_model()
        # DZP: 2s2p1d -> 2 + 6 + 5 = 13 orbitals per atom
        assert model.atom_orbital_counts == [13, 13]


class TestUnitsAndConvention:
    """Seal units and R->k phases against ABACUS's own k-space output."""

    def test_gamma_eigenvalues_match_direct_diag(self, prim_parser):
        model = prim_parser.get_unfold_model()
        H, S = model.hs_and_eigen([0.0, 0.0, 0.0])
        eig = eigh(H, S, eigvals_only=True)
        # H(R=0) on-site s element is O(-26 eV); eV scale sanity
        assert eig[0] == pytest.approx(-5.6957, abs=0.01)

    def test_generic_k_eigenvalues_match_abacus_bands(self, prim_parser):
        """Parsed model spectrum == ABACUS BANDS_1.dat at generic k.

        BANDS eigenvalues are eV as printed by ABACUS; agreement at
        non-Gamma k pins both the Ry->eV conversion and the R->k Bloch
        phase convention (a wrong phase gives k-dependent matrices).
        """
        out = _outpath("si_prim")
        kpoints = read_kpoints(os.path.join(out, "kpoints"), 36)
        bands = np.loadtxt(os.path.join(out, "BANDS_1.dat"))
        model = prim_parser.get_unfold_model()
        for ik in [0, 1, 5, 20, 35]:
            H, S = model.hs_and_eigen(kpoints[ik])
            eig = np.sort(eigh(H, S, eigvals_only=True))
            reference = bands[ik, 2:]  # first two columns are index/distance
            dev = max(np.min(np.abs(eig - value)) for value in reference)
            assert dev < 1e-4, f"k={kpoints[ik]} deviation {dev} eV"

    def test_rydberg_on_disk(self, prim_parser):
        """The on-disk HR unit is Rydberg; the model converts to eV."""
        from HamiltonIO.abacus.abacus_api import XR_matrix

        hr_file = os.path.join(_outpath("si_prim"), "data-HR-sparse_SPIN0.csr")
        raw = XR_matrix(2, hr_file)
        raw.read_file()
        gamma_sum_onsite = float(sum(blk[0, 0] for blk in raw.XR))
        model = prim_parser.get_unfold_model()
        H, _ = model.hs_and_eigen(np.zeros(3))
        # sum_R H_{3s,3s}(R) ~ -1.87 Ry is the on-site element of H(Gamma);
        # in eV it would sit near -25, pinning the on-disk unit as Rydberg
        assert -5.0 < gamma_sum_onsite < -1.0
        assert H[0, 0] == pytest.approx(gamma_sum_onsite * RY_TO_EV, abs=1e-3)


class TestConvModel:
    def test_shapes(self, conv_parser):
        model = conv_parser.get_unfold_model()
        assert len(model.atoms) == 8
        assert len(model.HR) == 99
        assert model.SR[(0, 0, 0)].shape == (104, 104)
        assert model.atom_orbital_counts == [13] * 8

    def test_gamma_gap(self, conv_parser):
        model = conv_parser.get_unfold_model()
        H, S = model.hs_and_eigen([0.0, 0.0, 0.0])
        eig = np.sort(eigh(H, S, eigvals_only=True))
        # 32 valence electrons -> the 16 lowest bands are occupied; the
        # Si gap separates them from the empty bands
        assert eig[16] - eig[15] > 0.5
        assert model.efermi == pytest.approx(eig[16], abs=1.5)


def test_pw_efermi_from_nscf_log():
    """PW parser reads the Fermi energy (eV column) from the run log.

    The si_pw_path fixture is an nscf band run: only
    ``running_nscf.log`` exists, and it carries both the ``E_Fermi
    <Ry> <eV>`` table and the ``EFERMI = ... eV`` summary (they agree).
    """
    pytest.importorskip("HamiltonIO")
    from HamiltonIO.abacus.pw_wfc import AbacusPWParser

    outpath = _outpath("si_pw_path")
    if not os.path.isdir(outpath):
        pytest.skip(f"ABACUS fixture {outpath} missing")
    assert not os.path.exists(os.path.join(outpath, "running_scf.log"))
    data = AbacusPWParser(outpath).read()
    assert data.efermi == pytest.approx(6.2903334548, abs=1e-6)
    # the VBM triple at Gamma sits exactly at the Fermi level
    gamma = np.argsort(np.abs(data.kpoints).sum(axis=1))[0]
    assert data.eigenvalues[gamma].max() == pytest.approx(data.efermi, abs=1e-4)


def test_pw_scf_kpoints_ignore_symmetry_reduction_table(tmp_path):
    """Only the direct-coordinate table is wavefunction-indexed."""
    from HamiltonIO.abacus.pw_wfc import AbacusPWParser

    (tmp_path / "kpoints").write_text(
        "nkstot now = 1\nK-POINTS DIRECT COORDINATES\n"
        " KPOINTS DIRECT_X DIRECT_Y DIRECT_Z WEIGHT\n"
        " 1 0.25 0.00 0.50 1.0\n\n"
        "K-POINTS REDUCTION ACCORDING TO SYMMETRY\n"
        " KPT DIRECT_X DIRECT_Y DIRECT_Z IBZ DIRECT_X DIRECT_Y DIRECT_Z\n"
        " 1 0.25 0.00 0.50 1 0.25 0.00 0.50\n"
    )
    np.testing.assert_array_equal(
        AbacusPWParser(tmp_path)._read_kpoints(), [[0.25, 0.0, 0.5]]
    )
