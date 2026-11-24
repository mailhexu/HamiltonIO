"""SislParser(sparse=True): H(k) and S(k) from the sparse matrices, no dense H(R).

The sparse mode must give the same H(k) and S(k) as the default dense mode for
unpolarized, collinear and spin-orbit Hamiltonians, with and without Loewdin
orthogonalization, while not building the real-space arrays.

The tests write small random Hamiltonians with sisl (TSHS + fdf geometry), so no
SIESTA run is needed.
"""

import numpy as np
import pytest

sisl = pytest.importorskip("sisl")

from HamiltonIO.siesta.sisl_wrapper import SislParser  # noqa: E402

_KPTS = [(0.0, 0.0, 0.0), (0.5, 0.0, 0.0), (0.21, -0.37, 0.0)]


def _hermitian(H):
    try:
        return (H + H.transpose(conjugate=True)) / 2
    except TypeError:  # sisl < 0.15.3
        return (H + H.transpose(hermitian=True)) / 2


def _write_siesta(path, spin, seed=0):
    """Random Hermitian non-orthogonal Hamiltonian with a 3x3x1 supercell."""
    rng = np.random.default_rng(seed)
    orbs = [sisl.AtomicOrbital("2s", R=2.5), sisl.AtomicOrbital("2pz", R=2.5)]
    geom = sisl.Geometry(
        [[0.0, 0.0, 0.0], [1.4, 0.6, 0.2]],
        atoms=[sisl.Atom(6, orbs), sisl.Atom(8, orbs)],
        lattice=sisl.Lattice([3.0, 3.2, 10.0], nsc=[3, 3, 1]),
    )
    H = sisl.Hamiltonian(geom, spin=sisl.Spin(spin), orthogonal=False)
    for i in range(geom.no):
        for j in range(geom.no_s):
            H[i, j] = list(rng.standard_normal(H.dim - 1)) + [
                0.05 * rng.standard_normal()
            ]
    H = _hermitian(H)
    for i in range(geom.no):
        H[i, i, H.S_idx] = 1.0
    H.write(str(path / "siesta.TSHS"))
    geom.write(str(path / "siesta.fdf"))
    with open(path / "siesta.fdf", "a") as f:
        f.write("\nSystemLabel siesta\n")
    return str(path / "siesta.fdf")


def _models(fdf, **kwargs):
    m = SislParser(fdf_fname=fdf, **kwargs).get_model()
    return m if isinstance(m, tuple) else (m,)


@pytest.mark.parametrize("orth", [False, True])
@pytest.mark.parametrize("spin", ["unpolarized", "polarized", "spin-orbit"])
def test_sparse_hk_matches_dense(tmp_path, spin, orth):
    fdf = _write_siesta(tmp_path, spin)
    dense = _models(fdf, orth=orth)
    sparse = _models(fdf, orth=orth, sparse=True)
    assert len(dense) == len(sparse)
    for md, ms in zip(dense, sparse):
        assert ms.HR is None and ms.SR is None
        assert ms.nbasis == md.nbasis
        for k in _KPTS:
            Hd, Sd = md.gen_ham(k)
            Hs, Ss = ms.gen_ham(k)
            np.testing.assert_allclose(Hs, Hd, atol=1e-12)
            if orth:
                assert Sd is None and Ss is None
            else:
                np.testing.assert_allclose(Ss, Sd, atol=1e-12)
            np.testing.assert_allclose(ms.HSE_k(k)[2], md.HSE_k(k)[2], atol=1e-10)


def test_sparse_rejects_split_soc(tmp_path):
    fdf = _write_siesta(tmp_path, "spin-orbit")
    with pytest.raises(ValueError, match="dense H"):
        SislParser(fdf_fname=fdf, read_H_soc=True, sparse=True).get_model()


@pytest.mark.parametrize("spin", ["unpolarized", "polarized", "spin-orbit"])
def test_sparse_real_space_accessors(tmp_path, spin):
    fdf = _write_siesta(tmp_path, spin)
    for dense, sparse in zip(_models(fdf), _models(fdf, sparse=True)):
        for R in [(0, 0, 0), (1, 0, 0), (-1, 1, 0)]:
            np.testing.assert_allclose(sparse.get_hamR(R), dense.get_hamR(R))
        np.testing.assert_allclose(sparse.get_H0(), dense.get_H0())
        dense_blocks = dense.get_intra_atomic_blocks()
        sparse_blocks = sparse.get_intra_atomic_blocks()
        assert dense_blocks.keys() == sparse_blocks.keys()
        for atom in dense_blocks:
            np.testing.assert_allclose(
                sparse_blocks[atom]["H_full"], dense_blocks[atom]["H_full"]
            )
        assert sparse.HR is None and sparse.SR is None


def test_sparse_noncollinear_real_space_fourier_transform(tmp_path):
    # The dense parser currently cannot read sisl's four-component NC layout.
    model, = _models(_write_siesta(tmp_path, "non-collinear"), sparse=True)
    for k in _KPTS:
        reconstructed = sum(
            model.get_hamR(R) * np.exp(2j * np.pi * np.dot(k, R))
            for R in model.Rlist
        )
        np.testing.assert_allclose(reconstructed, model.gen_ham(k)[0], atol=1e-12)
    assert model.HR is None and model.SR is None
