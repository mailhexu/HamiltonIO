"""GPAW readers for band unfolding: LCAO model tables and planewave eigendata.

Two routes into the unfolding package's interfaces (story 017):

``GpawLcaoModel``
    Real-space LCAO Hamiltonian/overlap tables built from a converged
    LCAO calculation (live calculator, or a ``mode='all'`` ``.gpw``
    restart). Exposes the backend interface of
    ``unfolding.lcao_unfolder``: ``.atoms`` (the supercell), ``.HR`` /
    ``.SR`` dicts keyed by integer translation tuples, and
    ``hs_and_eigen(k) -> (H, S)`` at any supercell fractional k-point.

``GpawPWParser``
    Planewave eigendata from a ``mode='all'`` ``.gpw`` restart,
    mirroring the SIESTA :class:`~HamiltonIO.siesta.wfsx.WFSXData`
    surface: fractional k-points, eigenvalues (eV), and per-k expansion
    coefficients on the stored plane-wave grids together with their
    integer G vectors (consumable by ``unfolding.pw_unfolder``).

Conventions (pinned against gpaw 26.7.0 with the default new
``gpaw.new`` calculators; see tests/test_gpaw_unfold.py):

- k-points are supercell fractional coordinates; GPAW plane waves carry
  the ``e^{-2 pi i k.f}``-family Bloch phases, matching HamiltonIO
  convention 2 on the same fractional coordinates.
- LCAO: GPAW's own transform (``gpaw.lcao.tightbinding``) is
  ``A(k) = sum_R A(R) e^{-2 pi i k.R}``; convention 2 wants
  ``A(k) = sum_R A(R) e^{+2 pi i k.R}``. The stored convention-2 block
  at translation ``T`` is therefore the plain inverse lattice Fourier
  transform over the full grid. Even-grid Nyquist translations are
  distributed equally over their positive and negative real-space images;
  this leaves sampled matrices unchanged and preserves Hermiticity at
  generic k-points. The resulting tables follow HamiltonIO convention 2.
- The k-grid must be a full (symmetry-off) **Gamma-centered**
  Monkhorst-Pack grid (``kpts={'size': (n, n, n), 'gamma': True}``):
  ``get_lcao_hamiltonian`` is evaluated at every stored k-point, the
  real-space transform divides by ``prod(N_c)``, and the Hermiticity of
  the derived tables relies on the grid containing Gamma. GPAW's
  default even-sized grids are half-shifted (``(2, 2, 2)`` stores
  ``{+-0.25}**3``) and yield non-Hermitian tables because the LCAO
  matrices carry orbital position phases; such files are rejected with
  a remedy.
- H is returned in eV, S dimensionless, both as full Hermitian matrices
  (``get_lcao_hamiltonian`` returns eV and already symmetrizes the
  legacy triangle storage; the new calculator returns full matrices).
- The matrices carry all PAW contributions: ``S_kMM`` is the PAW
  overlap (pseudo-AO overlap plus augmentation ``dS``) and ``H_kMM``
  includes the ``dH`` projector terms, so unfolding weights built from
  these tables need no additional PAW correction at the matrix level.
  (The projector-only weights of a full all-electron PW reconstruction
  are a separate research item; not needed here.)
- Planewave coefficients are pseudo-wavefunction coefficients `c_nG`,
  renormalized so `sum_G |c_nG|^2 = 1` per band. GPAW normalizes in the
  PAW S metric instead; pseudo-L2 norm can differ from one. Per-band
  scaling preserves relative pseudo coset fractions, but those are not
  physical PAW all-electron weights. A PAW-consistent projection needs
  augmentation overlaps between the primitive bra and the supercell state.
  Integer G vectors are recovered exactly from
  ``desc.G_plus_k_Gv = (g + kpt) @ (2 pi * icell)``.
"""

from dataclasses import dataclass
from typing import Tuple

import numpy as np
from ase.units import Ha

from HamiltonIO.gpaw.gpaw_api import get_hr_sr_from_calc

__all__ = [
    "GpawLcaoModel",
    "GpawPWParser",
    "GpawPWData",
]


def _hermitize(mat: np.ndarray) -> np.ndarray:
    """Return the Hermitian part of ``mat`` (keeps Hermitian input exact)."""
    return 0.5 * (mat + np.conj(mat.swapaxes(-1, -2)))


def _grid_dims(ibzk_qc: np.ndarray) -> np.ndarray:
    """Infer the Monkhorst-Pack dimensions of a full (symmetry-off) grid."""
    kpts = np.atleast_2d(np.asarray(ibzk_qc, dtype=float))
    dims = [len(np.unique(np.round(kpts[:, i], 10))) for i in range(3)]
    dims = np.asarray(dims, dtype=int)
    if np.prod(dims) != len(kpts):
        raise ValueError(
            f"k-point set with {len(kpts)} points is not a full "
            f"Monkhorst-Pack grid (unique counts per axis {dims.tolist()}); "
            "the LCAO real-space transform requires symmetry='off'"
        )
    return dims


def _translation_set(dims: np.ndarray) -> np.ndarray:
    """Real-space cell set of gpaw TightBinding, negated (see module doc)."""
    cells = np.indices(tuple(int(n) for n in dims)).reshape(3, -1).T
    half = dims[np.newaxis, :] // 2
    cells = ((cells + half) % dims) - half
    return -cells


class GpawLcaoModel:
    """Real-space LCAO tables of one GPAW spin channel (unfolding model).

    Attributes
    ----------
    atoms : ase.Atoms
        The supercell the calculation was run on (calculator detached).
    HR, SR : dict
        ``{T_tuple: (nao, nao) complex array}`` over the supercell
        translations resolved by the k-grid, HamiltonIO convention 2
        (``H(k) = sum_T H[T] e^{+2 pi i k.T}``), H in eV.
    hs_and_eigen : callable
        ``hs_and_eigen(k) -> (H, S)`` at a supercell fractional k-point.
    """

    def __init__(
        self, atoms, HR, SR, kgrid, kpoints, efermi=None, nspin=1, spin=0, nbands=None
    ):
        self._atoms = atoms
        self.HR = HR
        self.SR = SR
        self._T = np.array(sorted(HR), dtype=float)
        self._H_arr = np.array([HR[tuple(t)] for t in self._T])
        self._S_arr = np.array([SR[tuple(t)] for t in self._T])
        self.kgrid = tuple(int(n) for n in kgrid)
        self.kpoints = np.asarray(kpoints, dtype=float)
        self.efermi = efermi
        self.nspin = int(nspin)
        self.spin = int(spin)
        self.nbands = nbands

    # -- construction -------------------------------------------------------

    @classmethod
    def from_calc(cls, calc, spin: int = 0) -> "GpawLcaoModel":
        """Extract the tables from a converged LCAO calculator.

        Works with legacy and new (``gpaw.new``) calculators; the data
        must live on one MPI rank (run serially or collect first).
        """
        H_skMM, S_kMM = get_hr_sr_from_calc(calc)
        if H_skMM is None:
            raise RuntimeError(
                "get_lcao_hamiltonian returned None (non-master rank); "
                "run the extraction on a single process"
            )
        H_skMM = np.asarray(H_skMM)
        S_kMM = np.asarray(S_kMM)
        if H_skMM.ndim != 4:
            raise ValueError(
                f"expected H_skMM with shape (nspin, nk, nao, nao), got "
                f"{H_skMM.shape}"
            )
        nspin, nk, nao, nao2 = H_skMM.shape
        if nao != nao2 or S_kMM.shape != (nk, nao, nao):
            raise ValueError(
                f"inconsistent matrix shapes: H {H_skMM.shape}, S {S_kMM.shape}"
            )
        if not 0 <= spin < nspin:
            raise ValueError(f"spin={spin} out of range for nspin={nspin}")

        ibzk_qc = np.asarray(calc.get_ibz_k_points(), dtype=float)
        if len(ibzk_qc) != nk:
            raise ValueError(f"{nk} H matrices but {len(ibzk_qc)} IBZ k-points")
        dims = _grid_dims(ibzk_qc)
        if not np.any(np.all(np.abs(ibzk_qc) < 1e-8, axis=1)):
            raise ValueError(
                "the k-grid must be Gamma-centered (GPAW's default even "
                "grids are half-shifted, e.g. kpts=(2, 2, 2) gives "
                "{+-0.25}**3); pass kpts={'size': (n, n, n), 'gamma': True}. "
                "The real-space transform of half-shifted grids is not "
                "Hermitian because GPAW's LCAO matrices carry orbital "
                "position phases."
            )

        herm = np.abs(H_skMM[spin] - np.conj(np.swapaxes(H_skMM[spin], 1, 2))).max()
        if herm > 1e-8:
            raise ValueError(
                f"k-space Hamiltonian not Hermitian (max deviation {herm:.3e}); "
                "the .gpw restart likely lacks the wavefunctions - write it "
                "with mode='all'"
            )

        # Convention-2 blocks: (1/nk) sum_j H(k_j) e^{-2 pi i k_j.T} over the
        # negated TightBinding cell set (module docstring).
        translations = _translation_set(dims)
        phase = np.exp(-2j * np.pi * (ibzk_qc @ translations.T))  # (nk, nT)
        H_T = np.tensordot(phase, H_skMM[spin].astype(np.complex128), axes=(0, 0))
        S_T = np.tensordot(phase, S_kMM.astype(np.complex128), axes=(0, 0))
        H_T /= nk
        S_T /= nk

        # Even grids store a single Nyquist translation (+N/2), which is
        # indistinguishable from -N/2 on the sampled grid. Splitting its
        # Fourier coefficient equally between the two images preserves all
        # sampled matrices and restores H[-T] = H[T]^dagger at generic k.
        from itertools import product

        HR, SR = {}, {}
        for i, T in enumerate(translations):
            images = [
                (-int(t), int(t)) if n % 2 == 0 and abs(t) == n // 2 else (int(t),)
                for t, n in zip(T, dims)
            ]
            scale = 1.0 / np.prod([len(axis) for axis in images])
            for image in product(*images):
                HR[image] = H_T[i] * scale
                SR[image] = S_T[i] * scale

        atoms = calc.atoms.copy()
        atoms.calc = None
        efermi = None
        try:
            efermi = float(calc.get_fermi_level())
        except Exception:  # pragma: no cover - pickle-only paths
            pass
        return cls(atoms, HR, SR, dims, ibzk_qc, efermi=efermi, nspin=nspin, spin=spin)

    @classmethod
    def from_file(cls, gpw_path, spin: int = 0) -> "GpawLcaoModel":
        """Read the tables from a ``mode='all'`` LCAO ``.gpw`` restart."""
        from gpaw import restart

        atoms, calc = restart(str(gpw_path), txt=None)
        model = cls.from_calc(calc, spin=spin)
        model._atoms = atoms.copy()
        model._atoms.calc = None
        return model

    # -- model interface ----------------------------------------------------

    @property
    def atoms(self):
        return self._atoms

    def hs_and_eigen(self, k):
        """Return ``(H, S)`` at supercell fractional ``k`` (convention 2).

        Sampled matrices are reconstructed exactly. The paired Nyquist
        shells keep the unsampled interpolant Hermitian as well.
        """
        k = np.asarray(k, dtype=float).reshape(3)
        phase = np.exp(2j * np.pi * (k @ self._T.T))
        H = np.tensordot(phase, self._H_arr, axes=(0, 0))
        S = np.tensordot(phase, self._S_arr, axes=(0, 0))
        return _hermitize(H), _hermitize(S)

    def eigenvalues(self, kpoints):
        """Band energies (eV) at each supercell fractional k-point."""
        from scipy.linalg import eigh

        kpoints = np.atleast_2d(np.asarray(kpoints, dtype=float))
        return np.array(
            [eigh(*self.hs_and_eigen(k), eigvals_only=True) for k in kpoints]
        )


@dataclass
class GpawPWData:
    """Aggregated planewave contents of a GPAW ``.gpw`` restart.

    Attributes
    ----------
    kpoints : (nk, 3) float ndarray
        Supercell fractional k-points of the stored (symmetry-off) grid.
    kpoints_cart : (nk, 3) float ndarray
        Cartesian k-points in 1/Angstrom including the 2*pi factor,
        ``k_cart = k_frac @ 2*pi*inv(cell).T`` (WFSXData convention).
    eigenvalues : (nk, nband) float ndarray
        Eigenvalues in eV, GPAW-native (absolute, not Fermi-shifted).
    coefficients : tuple of (nband, npw_k) complex ndarrays
        Pseudo-wavefunction coefficients on the stored plane-wave grid
        of each k-point, renormalized to ``sum_G |c|**2 = 1`` per band.
    gvecs : tuple of (npw_k, 3) int ndarrays
        Integer supercell reciprocal-lattice vectors of each k-point's
        plane waves (same row order as the coefficients).
    efermi : float
        Fermi level in eV of the restarted calculation.
    cell : (3, 3) float ndarray
        Supercell lattice vectors in Angstrom (rows).
    nspin, spin : int
        Number of collinear spin channels of the calculation and the
        channel that was read.
    nbands : int
        Number of stored bands per k-point.
    """

    kpoints: np.ndarray
    kpoints_cart: np.ndarray
    eigenvalues: np.ndarray
    coefficients: Tuple[np.ndarray, ...]
    gvecs: Tuple[np.ndarray, ...]
    efermi: float
    cell: np.ndarray
    nspin: int
    spin: int
    nbands: int

    @property
    def npw(self):
        """Plane-wave count of each stored k-point."""
        return [len(g) for g in self.gvecs]


class GpawPWParser:
    """Parser for planewave GPAW restarts (unfolding eigendata source).

    Requires a ``.gpw`` written with ``mode='all'`` (wavefunctions
    included) from a symmetry-off run, and reads it back through
    ``gpaw.restart``. Serial reads only (the stored wavefunctions live
    on one rank).
    """

    def __init__(self, gpw_path: str, spin: int = 0):
        self.gpw_path = str(gpw_path)
        self.spin = int(spin)

    def read(self) -> GpawPWData:
        from gpaw import restart

        atoms, calc = restart(self.gpw_path, txt=None)
        ibzwfs = _get_ibzwfs(calc)
        ibz = ibzwfs.ibz
        kpt_kc = np.asarray(ibz.kpt_kc, dtype=float)
        cell = np.asarray(atoms.cell, dtype=float)
        bcart = 2.0 * np.pi * np.linalg.inv(cell).T

        kpoints, eig_list, coeff_list, gvec_list = [], [], [], []
        for wfs in ibzwfs:
            if wfs.spin != self.spin:
                continue
            ik = int(wfs.k)
            kpoints.append(kpt_kc[ik])
            eig_list.append(np.asarray(wfs.eig_n, dtype=float) * Ha)

            psit = wfs.psit_nX
            if getattr(psit, "data", None) is None:
                raise ValueError(
                    f"{self.gpw_path} carries no wavefunction data; rewrite it "
                    "with calc.write(name, mode='all')"
                )
            data = np.asarray(psit.data)
            if data.ndim != 2:
                raise NotImplementedError(
                    "only scalar (spinor-free) planewave states are supported; "
                    f"got coefficient array of shape {data.shape}"
                )
            norms = np.einsum("ng,ng->n", data, data.conj()).real
            if np.any(norms <= 0):
                raise ValueError(f"non-positive coefficient norm at k-point {ik}")
            coeff_list.append(data / np.sqrt(norms)[:, None])
            gvec_list.append(_integer_g_vectors(psit.desc))

        if not kpoints:
            raise ValueError(
                f"no wavefunctions for spin={self.spin} in {self.gpw_path}"
            )
        kpoints = np.asarray(kpoints, dtype=float)
        eigenvalues = np.asarray(eig_list, dtype=float)
        nbands_set = {e.shape[0] for e in eigenvalues}
        if len(nbands_set) != 1:
            raise ValueError(f"varying band counts across k-points: {nbands_set}")
        nspin = int(getattr(ibzwfs, "nspins", 1))
        efermi = float(calc.get_fermi_level())
        return GpawPWData(
            kpoints=kpoints,
            kpoints_cart=kpoints @ bcart,
            eigenvalues=eigenvalues,
            coefficients=tuple(coeff_list),
            gvecs=tuple(gvec_list),
            efermi=efermi,
            cell=cell,
            nspin=nspin,
            spin=self.spin,
            nbands=int(eigenvalues.shape[1]),
        )


def _get_ibzwfs(calc):
    """Return the new-calculator wavefunction collection."""
    dft = getattr(calc, "dft", None)
    if dft is None:
        raise NotImplementedError(
            "GpawPWParser needs a gpaw >= 25 new-generation calculator "
            "(gpaw.new); legacy calculators are not supported"
        )
    return dft.ibzwfs


def _integer_g_vectors(desc) -> np.ndarray:
    """Integer SC reciprocal vectors of a new-style PW descriptor's states.

    ``desc.G_plus_k_Gv`` stores ``(g + kpt) @ (2 pi * icell)`` in
    bohr^-1; inverting that relation recovers the integer ``g`` exactly.
    """
    gk = np.asarray(desc.G_plus_k_Gv, dtype=float)
    kpt = np.asarray(desc.kpt_c, dtype=float)
    bmat = np.asarray(desc.icell, dtype=float) * 2.0 * np.pi
    g = np.linalg.solve(bmat.T, (gk - kpt @ bmat).T).T
    gi = np.rint(g)
    if np.abs(g - gi).max() > 1e-6:
        raise ValueError(
            "could not recover integer G vectors from the PW descriptor "
            f"(max deviation {np.abs(g - gi).max():.3e})"
        )
    return np.asarray(gi, dtype=int)
