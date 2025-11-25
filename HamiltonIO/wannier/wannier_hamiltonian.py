import copy
import os
from collections import defaultdict

import numpy as np
from scipy.io import netcdf_file
from scipy.linalg import eigh
from scipy.sparse import csr_matrix

from HamiltonIO.hamiltonian import Hamiltonian

from .utils import auto_assign_basis_name
from .w90_parser import parse_atoms, parse_ham, parse_tb, parse_xyz


class WannierHam(Hamiltonian):
    def __init__(
        self,
        nbasis,
        data=None,
        R_degens=None,
        positions=None,
        sparse=False,
        ndim=3,
        nspin=1,
        double_site_energy=2.0,
    ):
        """
        :param nbasis: number of basis.
        :param data: a dictionary of {R: matrix}. R is a tuple, matrix is a nbasis*nbasis matrix.
        :param R_degens: degeneracy of R.
        :param positions: reduced positions.
        :param sparse: Bool, whether to use a sparse matrix.
        :param ndim: number of dimensions.
        :param nspin: number of spins.
        """
        if data is not None:
            self.data = data
        else:
            self.data = defaultdict(lambda: np.zeros((nbasis, nbasis), dtype=complex))
        self._nbasis = nbasis
        self._nspin = nspin
        self._norb = nbasis // nspin
        self._ndim = ndim
        if R_degens is not None:
            self.R_degens = R_degens
        else:
            self.R_degens = np.ones(len(self.data.keys()), dtype=int)
        if positions is None:
            self._positions = np.zeros((nbasis, self.ndim))
        else:
            self._positions = positions
        self.prepare_phase_rjri()
        self.sparse = sparse
        self.double_site_energy = double_site_energy
        if sparse:
            self._matrix = csr_matrix
        self.atoms = None
        self.R2kfactor = 2.0j * np.pi
        self.k2Rfactor = -2.0j * np.pi
        self.is_siesta = False
        self.is_orthogonal = True
        self._name = "Wannier"

    def set_atoms(self, atoms):
        self.atoms = atoms

    @property
    def nspin(self):
        return self._nspin

    @property
    def norb(self):
        """
        norb: number of orbitals, if spin/spinor, norb=nbasis/2
        """
        return self._norb

    @property
    def nbasis(self):
        return self._nbasis

    @property
    def ndim(self):
        return self._ndim

    @property
    def xcart(self):
        raise NotImplementedError()

    @property
    def xred(self):
        return self._positions

    @property
    def positions(self):
        return self._positions

    @property
    def onsite_energies(self):
        return self.data[(0, 0, 0)].diagonal() * 2

    @property
    def hoppings(self):
        """
        The hopping parameters, not including any onsite energy.
        """
        data = copy.deepcopy(self.data)
        np.fill_diagonal(data[(0, 0, 0)], 0.0)
        return data

    @staticmethod
    def read_from_wannier_dir(path, prefix, atoms=None, nls=True, groupby=None):
        """
        read tight binding model from a wannier function directory.
        :param path: path
        :param prefix: prefix to the wannier files, often wannier90, or wannier90_up, or wannier90_dn for vasp.
        """
        if atoms is None:
            atoms = parse_atoms(os.path.join(path, prefix + ".win"))
        cell = atoms.get_cell()
        tb_fname = os.path.join(path, prefix + "_tb.dat")
        hr_fname = os.path.join(path, prefix + "_hr.dat")
        if os.path.exists(tb_fname):
            xcart, nbasis, data, R_degens = parse_tb(fname=tb_fname)
            xred = cell.scaled_positions(xcart)
        else:
            nbasis, data, R_degens = parse_ham(fname=hr_fname)
            xyz_fname = os.path.join(path, prefix + "_centres.xyz")
            if os.path.exists(xyz_fname):
                has_xyz = True
                xcart, _, _ = parse_xyz(fname=xyz_fname)
                xred = cell.scaled_positions(xcart)
            else:
                raise FileNotFoundError(f"The file {xyz_fname} does not exist.")

        if groupby == "spin":
            # error message if nbasis is not even.
            norb = nbasis // 2
            xtmp = copy.deepcopy(xred)
            if has_xyz:
                xred[::2] = xtmp[:norb]
                xred[1::2] = xtmp[norb:]
            for key, val in data.items():
                dtmp = copy.deepcopy(val)
                data[key][::2, ::2] = dtmp[:norb, :norb]
                data[key][::2, 1::2] = dtmp[:norb, norb:]
                data[key][1::2, ::2] = dtmp[norb:, :norb]
                data[key][1::2, 1::2] = dtmp[norb:, norb:]
        if has_xyz:
            ind, positions = auto_assign_basis_name(xred, atoms)
        m = WannierHam(nbasis=nbasis, data=data, positions=xred, R_degens=R_degens)
        if has_xyz:
            nm = m.shift_position(positions)
        else:
            nm = m
        nm.set_atoms(atoms)
        return nm

    @staticmethod
    def load_banddownfolder(path, prefix, atoms=None, nls=True, groupby="spin"):
        from banddownfolder.scdm.lwf import LWF

        lwf = LWF.load_nc(fname=os.path.join(path, f"{prefix}.nc"))
        nbasis = lwf.nwann
        nspin = 1
        positions = lwf.wann_centers
        ndim = lwf.ndim
        H_mnR = defaultdict(lambda: np.zeros((nbasis, nbasis), dtype=complex))

        for iR, R in enumerate(lwf.Rlist):
            R = tuple(R)
            val = lwf.HwannR[iR]
            if np.linalg.norm(R) < 0.001:
                H_mnR[R] = val / 2.0
                # H_mnR[R] -= np.diag(np.diag(val) / 2.0)
            else:
                H_mnR[R] = val / 2.0
        m = WannierHam(nbasis, data=H_mnR, nspin=nspin, ndim=ndim, positions=positions)
        m.atoms = atoms
        return m

    def gen_ham(self, k, convention=2):
        """
        generate hamiltonian matrix at k point.
        H_k( i, j)=\sum_R H_R(i, j)^phase.
        There are two conventions,
        first:
        phase =e^{ik(R+rj-ri)}. often better used for berry phase.
        second:
        phase= e^{ikR}. We use the first convention here.

        :param k: kpoint
        :param convention: 1 or 2.
        """
        Hk = np.zeros((self.nbasis, self.nbasis), dtype="complex")
        if convention == 2:
            for iR, (R, mat) in enumerate(self.data.items()):
                phase = np.exp(self.R2kfactor * np.dot(k, R))  # / self.R_degens[iR]
                H = mat * phase
                Hk += H + H.conjugate().T
        elif convention == 1:
            for iR, (R, mat) in enumerate(self.data.items()):
                phase = (
                    np.exp(self.R2kfactor * np.dot(k, R + self.rjminusri))
                    # / self.R_degens[iR]
                )
                H = mat * phase
                Hk += H + H.conjugate().T
        else:
            raise ValueError("convention should be either 1 or 2.")
        return Hk

    def solve(self, k, convention=2):
        Hk = self.gen_ham(k, convention=convention)
        return eigh(Hk)

    def HSE_k(self, kpt, convention=2):
        H = self.gen_ham(tuple(kpt), convention=convention)
        S = None
        evals, evecs = eigh(H)
        return H, S, evals, evecs

    def solve_all(self, kpts, convention=2):
        _, _, evals, evecs = self.HS_and_eigen(kpts, convention=convention)
        return evals, evecs

    def HS_and_eigen(self, kpts, convention=2):
        """
        calculate eigens for all kpoints.
        :param kpts: list of k points.
        """
        nk = len(kpts)
        hams = np.zeros((nk, self.nbasis, self.nbasis), dtype=complex)
        evals = np.zeros((nk, self.nbasis), dtype=float)
        evecs = np.zeros((nk, self.nbasis, self.nbasis), dtype=complex)
        for ik, k in enumerate(kpts):
            hams[ik], S, evals[ik], evecs[ik] = self.HSE_k(
                tuple(k), convention=convention
            )
        return hams, None, evals, evecs

    def prepare_phase_rjri(self):
        """
        The matrix P: P(i, j) = r(j)-r(i)
        """
        self.rjminusri = self.xred[None, :, :] - self.xred[:, None, :]

    def to_sparse(self):
        for key, val in self.data:
            self.data[key] = self._matrix(val)

    @property
    def Rlist(self):
        return list(self.data.keys())

    @property
    def nR(self):
        """
        number of R
        """
        return len(self.Rlist)

    @property
    def site_energies(self):
        """
        on site energies.
        """
        return self.data[(0, 0, 0)].diagonal() * 2

    @property
    def ham_R0(self):
        """
        return hamiltonian at R=0. Note that the data is halfed for R=0.
        """
        return self.data[(0, 0, 0)] + self.data[(0, 0, 0)].T.conj()

    def get_hamR(self, R):
        """
        return the hamiltonian at H(i, j) at R.
        """
        nzR = np.nonzero(R)[0]
        if len(nzR) != 0 and R[nzR[0]] < 0:
            newR = tuple(-np.array(R))
            return self.data[newR].T.conj()
        elif len(nzR) == 0:
            newR = R
            mat = self.data[newR]
            return mat + self.data[(0, 0, 0)].T.conj()
        else:
            newR = R
            return self.data[newR]

    @staticmethod
    def _positive_R_mat(R, mat):
        nzR = np.nonzero(R)[0]
        if len(nzR) != 0 and R[nzR[0]] < 0:
            newR = tuple(np.array(-R))
            newmat = mat.T.conj()
        elif len(nzR) == 0:
            newR = R
            newmat = (mat + mat.T.conj()) / 2.0
        else:
            newR = R
            newmat = mat
        return newR, newmat

    def _to_positive_R(self):
        """
        make all the R positive.
        t(i, j, R) = t(j, i, -R).conj() if R is negative.
        """
        new_WannierHam = WannierHam(self.nbasis, sparse=self.sparse)
        for R, mat in self.data:
            newR, newmat = self._positive_R_mat(R, mat)
            new_WannierHam[newR] += newmat
        return new_WannierHam

    def shift_position(self, rpos):
        """
        shift the positions of basis set to near reference positions.
        E.g. reduced position 0.8, with refernce 0.0 will goto -0.2.
        This can move the wannier functions to near the ions.
        """
        pos = self.positions
        shift = np.zeros((self.nbasis, self.ndim), dtype="int")
        shift[:, :] = np.round(pos - rpos)
        newpos = copy.deepcopy(pos)
        for i in range(self.nbasis):
            newpos[i] = pos[i] - shift[i]
        d = WannierHam(
            self.nbasis, ndim=self.ndim, nspin=self.nspin, R_degens=self.R_degens
        )
        d._positions = newpos

        for R, v in self.data.items():
            for i in range(self.nbasis):
                for j in range(self.nbasis):
                    sR = tuple(np.array(R) - shift[i] + shift[j])
                    nzR = np.nonzero(sR)[0]
                    if len(nzR) != 0 and sR[nzR[0]] < 0:
                        newR = tuple(-np.array(sR))
                        d.data[newR][j, i] += v[i, j].conj()
                    elif len(nzR) == 0:
                        newR = sR
                        d.data[newR][i, j] += v[i, j] * 0.5
                        d.data[newR][j, i] += v[i, j].conj() * 0.5
                    else:
                        d.data[sR][i, j] += v[i, j]
        return d

    def save(self, fname):
        """
        Save model into a netcdf file.
        :param fname: filename.
        """
        # from netCDF4 import Dataset
        # root = Dataset(fname, 'w', format="NETCDF4")
        root = netcdf_file(fname, mode="w")
        root.createDimension("nR", self.nR)
        root.createDimension("ndim", self.ndim)
        root.createDimension("nbasis", self.nbasis)
        root.createDimension("nspin", self.nspin)
        root.createDimension("natom", len(self.atoms))
        R = root.createVariable("R", "i4", ("nR", "ndim"))
        data_real = root.createVariable("data_real", "f8", ("nR", "nbasis", "nbasis"))
        data_imag = root.createVariable("data_imag", "f8", ("nR", "nbasis", "nbasis"))
        positions = root.createVariable("positions", "f8", ("nbasis", "ndim"))

        if self.atoms is not None:
            atom_numbers = root.createVariable("atom_numbers", "i4", ("natom",))
            atom_xred = root.createVariable("atom_xred", "f8", ("natom", "ndim"))
            atom_cell = root.createVariable("atom_cell", "f8", ("ndim", "ndim"))

        atom_cell.unit = "Angstrom"
        positions.unit = "1"
        data_real.unit = "eV"
        data_imag.unit = "eV"

        R[:] = np.array(self.Rlist)
        d = np.array(tuple(self.data.values()))
        data_real[:] = np.real(d)
        data_imag[:] = np.imag(d)
        positions[:] = np.array(self.positions)

        if self.atoms is not None:
            atom_numbers[:] = np.array(self.atoms.get_atomic_numbers())
            atom_xred[:] = np.array(self.atoms.get_scaled_positions(wrap=False))
            atom_cell[:] = np.array(self.atoms.get_cell())
        root.close()

    @staticmethod
    def from_tbmodel(model):
        """
        translate from a tbmodel type tight binding model
        """
        ret = WannierHam(nbasis=model.size)
        for R, v in model.hop.items():
            ret.data[R] = v
        ret._positions = np.reshape(model.pos, (model.size, model.dim))
        return ret

    @staticmethod
    def from_tbmodel_hdf5(fname):
        """
        load model from a hdf5 file. It uses the tbmodel parser.
        """

        from tbmodels import Model

        m = Model.from_hdf5_file(fname)
        ret = WannierHam(nbasis=m.size)
        for R, v in m.hop.items():
            ret.data[R] = v
        ret.positions = np.reshape(m.pos, (m.size, m.ndim))
        return ret

    def to_spin_polarized(self, order=1):
        """
        repeat to get spin polarized.
        order =1 : orb1_up, orb1_dn, orb2_up, orb2_dn...
        order =2 : orb1_up, orb2_up, ... orb1_dn, orb2_dn...
        """
        ret = WannierHam(self.nbasis * 2)
        if self.positions is None:
            ret.positions = None
        else:
            ret.positions = np.repeat(self.positions, 2, axis=0)
        for R, mat in self.data.items():
            if order == 1:
                ret.data[R][::2, ::2] = mat
                ret.data[R][1::2, 1::2] = mat
            elif order == 2:
                ret.data[R][: self.norb, : self.norb] = mat
                ret.data[R][self.norb :, self.norb :] = mat
        return ret

    def validate(self):
        # make sure all R are 3d.
        for R in self.data.keys():
            if len(R) != self.ndim:
                raise ValueError("Dimension of R should be ndim %s" % (self.ndim))

    def iter_unique_hoppings(self):
        """Yield unique hopping matrix elements.

        Iterates over symmetry-reduced hopping terms using the convention
        t(i, j, R) = t(j, i, -R)^*.

        - Only R with first non-zero component >= 0 are used.
        - For R == (0, 0, 0), only i <= j are returned.

        Yields:
            tuple: (R, i, j, value)
        """
        for R, mat in self.data.items():
            nz = np.nonzero(R)[0]
            if len(nz) and R[nz[0]] < 0:
                # skip negative R; covered by corresponding positive R term
                continue

            if len(nz) == 0:
                # onsite / R=0: keep upper triangle including diagonal
                for i in range(self.nbasis):
                    for j in range(i, self.nbasis):
                        val = mat[i, j]
                        if val != 0:
                            yield R, i, j, val
            else:
                # offsite positive-R matrix elements: take all i, j
                for i in range(self.nbasis):
                    for j in range(self.nbasis):
                        val = mat[i, j]
                        if val != 0:
                            yield R, i, j, val

    def distance_resolved_hoppings(self, cell, use_absolute=True):
        """Return unique hoppings annotated with Wannier-center distances.

        Args:
            cell (array-like): 3x3 lattice vectors in Angstrom.
            use_absolute (bool): If True (default), use absolute distance
                between Wannier centers r_j + R and r_i.

        Returns:
            list[dict]: Each item has keys
                "R", "i", "j", "distance", "hopping".
        """
        positions_frac = np.asarray(self.positions)
        cell = np.asarray(cell)

        results = []
        for R, i, j, hij in self.iter_unique_hoppings():
            # real-space vector from i to j includes lattice translation R
            d_frac = positions_frac[j] + np.asarray(R) - positions_frac[i]
            d_cart = d_frac @ cell
            if use_absolute:
                dist = float(np.linalg.norm(d_cart))
            else:
                dist = d_cart
            results.append(
                {
                    "R": tuple(R),
                    "i": int(i),
                    "j": int(j),
                    "distance": dist,
                    "hopping": complex(hij),
                }
            )
        return results

    @staticmethod
    def bin_hoppings_by_distance(entries, dr=0.1):
        """Bin hopping magnitudes as a function of distance.

        Args:
            entries (Iterable[dict]): Output of ``distance_resolved_hoppings``.
            dr (float): Bin width in Angstrom.

        Returns:
            tuple[np.ndarray, np.ndarray]: (bin_centers, avg_abs_h)
        """
        if not entries:
            return np.array([]), np.array([])

        distances = np.array([e["distance"] for e in entries], dtype=float)
        mags = np.abs([e["hopping"] for e in entries])

        dmax = float(distances.max())
        nbins = int(np.ceil(dmax / dr)) or 1
        edges = np.linspace(0.0, nbins * dr, nbins + 1)

        sums = np.zeros(nbins, dtype=float)
        counts = np.zeros(nbins, dtype=int)

        inds = np.clip((distances / dr).astype(int), 0, nbins - 1)
        for idx, mag in zip(inds, mags):
            sums[idx] += mag
            counts[idx] += 1

        with np.errstate(divide="ignore", invalid="ignore"):
            avg = np.where(counts > 0, sums / counts, 0.0)

        centers = (edges[:-1] + edges[1:]) * 0.5
        return centers, avg


def merge_tbmodels_spin(tbmodel_up, tbmodel_dn):
    """
    Merge a spin up and spin down model to one spinor model.
    """
    tbmodel = WannierHam(
        nbasis=tbmodel_up.nbasis * 2,
        data=None,
        positions=np.vstack([tbmodel_up.positions, tbmodel_dn.positions]),
        sparse=False,
        ndim=tbmodel_up.ndim,
        nspin=2,
        double_site_energy=2.0,
    )
    norb = tbmodel.norb
    for R in tbmodel_up.data:
        tbmodel.data[R][:norb, :norb] = tbmodel_up.data[R][:, :]
        tbmodel.data[R][norb:, norb:] = tbmodel_dn.data[R][:, :]
    return tbmodel
