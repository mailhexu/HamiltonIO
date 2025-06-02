import os
from collections import defaultdict

import numpy as np
from ase.atoms import Atoms
from scipy.linalg import eigh
from TB2J.utils import symbol_number

from HamiltonIO.lcao_hamiltonian import LCAOHamiltonian
from HamiltonIO.mathutils.lowdin import Lowdin
from HamiltonIO.mathutils.rotate_spin import (
    rotate_spinor_matrix,
)
from HamiltonIO.model.kR_convert import R_to_onek
from HamiltonIO.siesta.mysiesta_nc import MySiestaNC

try:
    import sisl
except Exception as e:
    print(e)
    print("Sisl is not installed, please install it by 'pip install sisl'")


class SiestaHamiltonian(LCAOHamiltonian):
    def __init__(
        self,
        HR,
        SR,
        Rlist,
        nbasis,
        atoms=None,
        orbs=None,
        nspin=1,
        HR_nosoc=None,
        HR_soc=None,
        nel=None,
        orth=False,
    ):
        """
        Hamiltonian object for Siesta
        params:
        HR: numpy array. index (iR, iorb1, iorb2)
             Real space Hamiltonian
        SR: numpy array. index (iR, iorb1, iorb2)
            Real space overlap matrix
        Rlist: list of numpy array. index (iR, 3)
            Real space lattice vectors
        nbasis: int
            number of basis
        atoms: Atoms object
            Atoms object
        orbs: list of str
            list of orbital names
        nspin: int
            number of spin
        HR_nosoc: numpy array. index (iR, iorb1, iorb2)
            Real space Hamiltonian without SOC
        HR_soc: numpy array. index (iR, iorb1, iorb2)
            Real space Hamiltonian with SOC
        orth: bool
            whether the Hamiltonian is orthogonalized with Lowdin's method
        """

        super().__init__(
            HR=HR,
            SR=SR,
            Rlist=Rlist,
            nbasis=nbasis,
            atoms=atoms,
            orbs=orbs,
            nspin=nspin,
            HR_nosoc=HR_nosoc,
            HR_soc=HR_soc,
            nel=nel,
            orth=orth,
        )
        self._name = "SIESTA"
        self.is_orthogonal = False
        self.R2kfactor = 2j * np.pi


class SislParser:
    """
    Siesta Hamiltonian parser using sisl as backend
    """

    def __init__(self, fdf_fname=None, ispin=None, read_H_soc=False, orth=False):
        self.fdf = sisl.get_sile(fdf_fname)
        self.ham = self.fdf.read_hamiltonian()
        self.spin = self.ham.spin
        self.ispin = ispin
        self.read_H_soc = read_H_soc
        self.orth = orth

    def get_model(self):
        self.read_spin(self.ispin)
        self.geom, self.atoms = self.read_geometry(fdf=self.fdf)
        orbs, nbasis, norb = self.read_basis()
        Rlist = self.read_Rlist(self.geom)
        self.nR = len(Rlist)
        HR, SR = self.read_HS(self.nR, norb)
        nel = self.read_nel(fdf=self.fdf)
        HR_soc = None
        HR_nosoc = None

        so_strength = self.read_so_strength(self.fdf)
        if self.read_H_soc:
            from TB2J.pauli import chargepart, spinpart

            HR_soc = self.read_HR_soc(self.fdf)
            HR_nosoc = HR - HR_soc * so_strength
            # move the charge part of SOC to HR_nosoc
            for iR, _ in enumerate(Rlist):
                HR_nosoc[iR] += chargepart(HR_soc[iR])
                HR_soc[iR] = spinpart(HR_soc[iR])
            # HR_nosoc += chargepart(HR_soc)
            # HR_soc = spinpart(HR_soc)
            # print("max of real part of HR_soc", np.max(np.abs(np.real(HR_soc))))
            # print("max of imag part of HR_soc", np.max(np.abs(np.imag(HR_soc))))
            # print("max of real part of HR_nosoc", np.max(np.abs(np.real(HR_nosoc))))
            # print("max of imag part of HR_nosoc", np.max(np.abs(np.imag(HR_nosoc))))
            # print("max of real part of HR", np.max(np.abs(np.real(HR))))
            # print("max of imag part of HR", np.max(np.abs(np.imag(HR))))
        else:
            HR_soc = None
            HR_nosoc = None
        if self.spintype == "spinor" or self.spintype == "nonpolarized":
            model = SiestaHamiltonian(
                HR=HR,
                SR=SR,
                Rlist=Rlist,
                nbasis=nbasis,
                nspin=1,
                orbs=orbs,
                HR_nosoc=HR_nosoc,
                HR_soc=HR_soc,
                nel=nel,
                atoms=self.atoms,
                orth=self.orth,
            )
            return model
        else:
            model_up = SiestaHamiltonian(
                HR=HR[0],
                SR=SR,
                Rlist=Rlist,
                nbasis=nbasis,
                nspin=1,
                orbs=orbs,
                HR_nosoc=HR_nosoc,
                HR_soc=HR_soc,
                nel=nel,
                atoms=self.atoms,
                orth=self.orth,
            )
            model_down = SiestaHamiltonian(
                HR=HR[1],
                SR=SR,
                Rlist=Rlist,
                nbasis=nbasis,
                nspin=1,
                orbs=orbs,
                HR_nosoc=HR_nosoc,
                HR_soc=HR_soc,
                nel=nel,
                atoms=self.atoms,
                orth=self.orth,
            )
            model = (model_up, model_down)
            if self.ispin is not None:
                return model[self.ispin]
            else:
                return model
        return model

    def read_Rlist(self, fdf):
        return self.geom.sc_off

    def read_so_strength(self, fdf):
        return fdf.get("Spin.OrbitStrength", default=1.0)

    def read_HS(self, nR, norb, dense=True):
        smat = np.asarray(self.ham.tocsr(self.ham.S_idx).todense())
        SR = np.reshape(smat, (self.norb, self.nR, self.norb)).transpose((1, 0, 2))
        if self.spintype == "spinor":
            # SR = np.kron(SR, np.eye(2))
            # SR = np.kron(np.eye(2), SR)
            SR2 = np.zeros((nR, self.norb * 2, self.norb * 2), dtype=complex)
            SR2[:, ::2, ::2] = SR
            SR2[:, 1::2, 1::2] = SR
            SR = SR2
        HR = self._get_HR_all(dense=dense)

        return HR, SR

    def read_spin(self, ispin):
        self.spintype = None
        if self.ham.spin.is_colinear:
            if ispin == "merge":
                self.spintype = "spinor"
            else:
                self.spintype = "colinear"
        elif self.ham.spin.is_spinorbit or self.spin.is_noncolinear:
            self.spintype = "spinor"
        else:
            self.spintype = "nonpolarized"

        if ispin == "up":
            ispin = 0
        elif ispin == "down":
            ispin = 1
        if ispin not in [None, 0, 1, "merge"]:
            raise ValueError("ispin should be None/0/1, but is %s" % ispin)
        self.ispin = ispin
        if not self.ham.spin.is_colinear:
            self.ispin = None
        return self.ispin

    def read_basis(self):
        self.orbs = []
        self.orb_dict = defaultdict(lambda: [])

        sdict = list(symbol_number(self.atoms).keys())
        _atoms = self.ham._geometry._atoms
        if self.ham.spin.is_colinear and (self.ispin in [0, 1]):
            for ia, a in enumerate(_atoms):
                symnum = sdict[ia]
                try:
                    orb_names = [f"{symnum}|{x.name()}|up" for x in a.orbital]
                except Exception:
                    orb_names = [f"{symnum}|{x.name()}|up" for x in a.orbitals]
                    raise Warning(
                        "It seems you are using an acient version of sisl. Please consider updating sisl to the latest version"
                    )
                self.orbs += orb_names
                self.orb_dict[ia] += orb_names
            self.norb = len(self.orbs)
            self.nbasis = self.norb
        elif (
            self.ham.spin.is_spinorbit
            or self.ham.spin.is_noncolinear
            or self.ispin == "merge"
        ):
            # set the orbitals
            for ia, a in enumerate(_atoms):
                symnum = sdict[ia]
                orb_names = []
                for x in a.orbitals:
                    for spin in {"up", "down"}:
                        orb_names.append(f"{symnum}|{x.name()}|{spin}")
                self.orbs += orb_names
                self.orb_dict[ia] += orb_names
            self.norb = len(self.orbs) // 2
            self.nbasis = len(self.orbs)
        else:  # non-polarized, spin=None
            for ia, a in enumerate(_atoms):
                symnum = sdict[ia]
                orb_names = []
                for x in a.orbitals:
                    orb_names.append(f"{symnum}|{x.name()}")
                self.orbs += orb_names
                self.orb_dict[ia] += orb_names
            self.norb = len(self.orbs)
            self.nbasis = len(self.orbs)

        return self.orbs, self.nbasis, self.norb

    def read_geometry(self, fdf):
        geom = fdf.read_geometry()
        _atoms = geom._atoms
        atomic_numbers = []
        self.positions = geom.xyz
        self.cell = np.array(geom.lattice.cell)
        for ia, a in enumerate(_atoms):
            atomic_numbers.append(a.Z % 200)
            # %200 is for synthetic atoms and ghost atoms
        self.atoms = Atoms(
            numbers=atomic_numbers, cell=self.cell, positions=self.positions
        )
        return geom, self.atoms

    def get_nc_path(self, fdf):
        fname = fdf.dir_file(fdf.get("System_Label", default="siesta") + ".nc")
        return fname

    def get_so_strength(self, fdf):
        return fdf.get("Spin_strength", default=1.0)

    def read_nel(self, fdf):
        fname = self.get_nc_path(fdf)
        nc = MySiestaNC(fname)
        self.nel = nc.read_qtot()

        try:
            fname = self.get_nc_path(fdf)
            nc = MySiestaNC(fname)
            self.nel = nc.read_qtot()
            # self.nel = fdf.read_qtot()[0]
        except Exception as e:
            print(e)
            self.nel = None
        return self.nel

    # def read_from_sisl(self, fdf_fname, read_H_soc=False):
    #    if fdf_fname is not None:
    #        fname = self.get_nc_path(fdf)
    #        self.ham = fdf.read_hamiltonian()

    def read_HR_soc(self, fdf):
        # put the soc part of the hamiltonian into HR_soc
        fname = self.get_nc_path(fdf)
        nc = MySiestaNC(fname)
        ham_soc = nc.read_soc_hamiltonian()
        # ham_soc = fdf.read_soc_hamiltonian()
        mat = ham_soc._csr.todense()
        HR_soc = self.convert_sisl_to_spinorham(mat)
        return HR_soc

    def _get_HR_all_nonpolarized(self, dense=True):
        """
        Get the Hamiltonian matrix in real space
        """
        mat = self.ham._csr.todense()
        self._HR = mat.reshape((self.norb, self.nR, self.norb, 2))[
            :, :, :, 0
        ].transpose((1, 0, 2))
        return self._HR

    def _get_HR_all_colinear(self, dense=True, ispin=None):
        """
        Get the Hamiltonian matrix in real space
        """
        # mat = self.ham._csr.todense()
        norb, norb_sc, ndspin = self.ham._csr.shape
        ndspin = ndspin
        # HRs = np.zeros((2, self.nR, self.norb, self.norb), dtype=float)

        # if False:
        #    for ispin in range(2):
        #        for iR in range(self.nR):
        #            for iorb in range(self.norb):
        #                for jorb in range(self.norb):
        #                    HRs[ispin, iR, iorb, jorb] = self.ham[iorb, jorb, ispin, iR]

        dmat = self.ham._csr.todense()
        mat = dmat.reshape((self.norb, self.nR, self.norb, 3))
        self._HR = mat.transpose((3, 1, 0, 2))[:2, :, :, :]
        if ispin is not None:
            return self._HR[ispin]
        else:
            return self._HR

    def _get_HR_all_SOC(self, dense=True):
        """
        Get the Hamiltonian matrix in real space
        """
        mat = self.ham._csr.todense()
        HR = self.convert_sisl_to_spinorham(mat)
        return HR

    def convert_sisl_to_spinorham(self, mat):
        norb, norb_sc, ndspin = mat.shape
        # nbasis = norb * 2
        mat = mat.reshape((norb, self.nR, norb, ndspin)).transpose((1, 0, 2, 3))
        HRs = np.zeros((self.nR, norb * 2, norb * 2), dtype=complex)
        # up-up:
        HRs[:, ::2, ::2] = mat[:, :, :, 0] + 1j * mat[:, :, :, 4]
        # down-down:
        HRs[:, 1::2, 1::2] = mat[:, :, :, 1] + 1j * mat[:, :, :, 5]
        # up-down:
        HRs[:, ::2, 1::2] = mat[:, :, :, 2] + 1j * mat[:, :, :, 3]
        # down-up:
        HRs[:, 1::2, ::2] = mat[:, :, :, 6] + 1j * mat[:, :, :, 7]
        return HRs

    def _get_HR_all(self, dense=True):
        if self.ham.spin.is_colinear:
            return self._get_HR_all_colinear(dense=dense)
        elif self.ham.spin.is_spinorbit or self.ham.spin.is_noncolinear:
            return self._get_HR_all_SOC(dense=dense)
        else:
            return self._get_HR_all_nonpolarized(dense=dense)

    @property
    def HR(self):
        if self.split_soc:
            _HR = np.zeros_like(self.HR_soc)
            for iR, _ in enumerate(self.Rlist):
                theta, phi = self.soc_rotation_angle
                # _HR[iR] = (
                #    rotate_Matrix_from_z_to_spherical(self.HR_nosoc[iR], theta, phi)
                #    + self.HR_soc[iR]
                # )
                # _HR[iR] = rotate_Matrix_from_z_to_spherical(self.HR_soc[iR], theta, phi)
                _HR[iR] = (
                    rotate_spinor_matrix(self.HR_nosoc[iR], theta, phi)
                    + self.HR_soc[iR]
                )
            return _HR
        else:
            return self.get_HR_all()

    def view_info(self):
        print(self.orb_dict)
        print(self.atoms)

    def solve_with_sisl(self, k, convention=2):
        if convention == 1:
            gauge = "r"
        elif convention == 2:
            gauge = "R"
        if self.ispin in [0, 1]:
            evals, evecs = self.ham.eigh(
                k=k, spin=self.ispin, eigvals_only=False, gauge=gauge
            )
        elif self.ispin is None:
            evals, evecs = self.ham.eigh(k=k, eigvals_only=False, gauge=gauge)
        elif self.ispin == "merge":
            evals0, evecs0 = self.ham.eigh(k=k, spin=0, eigvals_only=False, gauge=gauge)
            evals1, evecs1 = self.ham.eigh(k=k, spin=1, eigvals_only=False, gauge=gauge)
            evals = np.zeros(self.nbasis, dtype=float)
            evecs = np.zeros((self.nbasis, self.nbasis), dtype=complex)
            evals[::2] = evals0
            evals[1::2] = evals1
            evecs[::2, ::2] = evecs0
            evecs[1::2, 1::2] = evecs1
        return evals, evecs

    def solve(self, k, convention=2):
        # if convention == 1:
        #    gauge = "r"
        # elif convention == 2:
        #    gauge = "R"
        if self.ispin in [0, 1]:
            evals, evecs = eigh(
                self.Hk(k, convention=convention), self.Sk(k, convention=convention)
            )
        elif self.ispin is None:
            evals, evecs = eigh(
                self.Hk(k, convention=convention), self.Sk(k, convention=convention)
            )
        elif self.ispin == "merge":
            evals0, evecs0 = eigh(
                self.Hk(k, spin=0, convention=convention),
                self.Sk(k, convention=convention),
            )
            evals1, evecs1 = eigh(
                self.Hk(k, spin=1, convention=convention),
                self.Sk(k, convention=convention),
            )
            evals, evals = self.merge_spin(evals0, evals1, evecs0, evecs1)
        return evals, evecs

    def merge_spin(self, evals0, evals1, evecs0, evecs1):
        evals = np.zeros(self.nbasis, dtype=float)
        evecs = np.zeros((self.nbasis, self.nbasis), dtype=complex)
        evals[::2] = evals0
        evals[1::2] = evals1
        evecs[::2, ::2] = evecs0
        evecs[1::2, 1::2] = evecs1
        return evals, evecs

    def Hk(self, k, convention=2, ispin=None):
        if ispin is None:
            ispin = self.ispin
        if ispin is None:
            HR = self.HR
            H = R_to_onek(k, self.Rlist, HR)
        elif ispin in [0, 1]:
            HR = self._get_HR_all_colinear(dense=True, ispin=ispin)
            H = R_to_onek(k, self.Rlist, HR)
        elif ispin == "merge":
            H = np.zeros((self.nbasis, self.nbasis), dtype="complex")
            H[::2, ::2] = self.Hk(k, spin=0, convention=convention)
            H[1::2, 1::2] = self.Hk(k, spin=1, convention=convention)
        return H

    def eigen(self, k, convention=2):
        return self.solve(k)

    def gen_ham(self, k, convention=2):
        return self.Hk(k, convention=convention)

    def Sk(self, k, convention=2):
        if convention == 1:
            pass
        elif convention == 2:
            pass
        S0 = self.ham.Sk(k, gauge="R", format="dense")
        if self.ispin is None:
            S = S0
        elif self.ispin in [0, 1]:
            S = S0
        elif self.ispin == "merge":
            S = np.zeros((self.nbasis, self.nbasis), dtype="complex")
            S[::2, ::2] = S0
            S[1::2, 1::2] = S0
        return S

    def solve_all(self, kpts, orth=False):
        evals = []
        evecs = []
        for ik, k in enumerate(kpts):
            if orth:
                S = self.Sk(k)
                Smh = Lowdin(S)
                H = self.gen_ham(k)
                Horth = Smh.T.conj() @ H @ Smh
                evalue, evec = eigh(Horth)
            else:
                evalue, evec = self.solve(k)
            evals.append(evalue)
            evecs.append(evec)
        return np.array(evals, dtype=float), np.array(evecs, dtype=complex, order="C")

    def HSE_k(self, k, convention=2):
        Hk = self.Hk(k, convention=convention)
        Sk = self.Sk(k, convention=convention)
        evalue, evec = eigh(Hk, Sk)
        return Hk, Sk, evalue, evec

    def HS_and_eigen(self, kpts, convention=2):
        nkpts = len(kpts)
        self.evals = np.zeros((nkpts, self.nbasis), dtype=float)
        self.nkpts = nkpts
        if not self._use_cache:
            self.evecs = np.zeros((nkpts, self.nbasis, self.nbasis), dtype=complex)
            self.H = np.zeros((nkpts, self.nbasis, self.nbasis), dtype=complex)
            self.S = np.zeros((nkpts, self.nbasis, self.nbasis), dtype=complex)
        else:
            self._prepare_cache()

        for ik, k in enumerate(kpts):
            if self._use_cache:
                self.evecs = np.memmap(
                    os.path.join(self.cache_path, "evecs.dat"),
                    mode="w+",
                    shape=(nkpts, self.nbasis, self.nbasis),
                    dtype=complex,
                )
                self.H = np.memmap(
                    os.path.join(self.cache_path, "H.dat"),
                    mode="w+",
                    shape=(nkpts, self.nbasis, self.nbasis),
                    dtype=complex,
                )
                self.S = np.memmap(
                    os.path.join(self.cache_path, "H.dat"),
                    mode="w+",
                    shape=(nkpts, self.nbasis, self.nbasis),
                    dtype=complex,
                )
            Hk = self.Hk(k, convention=convention)
            Sk = self.Sk(k, convention=convention)
            evalue, evec = eigh(Hk, Sk)
            self.H[ik] = Hk
            self.S[ik] = Sk
            self.evals[ik] = evalue
            self.evecs[ik] = evec
            if self._use_cache:
                del self.evecs
                del self.H
                del self.S
        return self.H, self.S, self.evals, self.evecs

    def _prepare_cache(self, path="./TB2J_results/cache"):
        if not os.path.exists(path):
            os.makedirs(path)
        else:
            os.remove(path)
        self.cache_path = path

    def get_HSE(self, kpt, convention=2):
        Hk = self.Hk(kpt, convention=convention)
        Sk = self.Sk(kpt, convention=convention)
        evalue, evec = self.solve(kpt, convention=convention)
        return Hk, Sk, evalue, evec

    def get_fermi_level(self):
        return 0.0


class SislWFSXWrapper(SislParser):
    """Wrapper for retrieving eigenvalues and eigenvectors from siesta WFSX file

    Parameters
    ----------
    geom : sisl.Geometry
        the geometry containing information about atomic positions and orbitals
    wfsx_sile: sisl.io.siesta.wfsxSileSiesta
        file reader for WFSX file
    spin : sisl.physics.Spin
        spin object carrying information about spin configurations and spin component.
    ispin : None or int
        index of spin channel to be considered. Only takes effect for collinear spin calculations (UP: 0, DOWN: 1).
        (default: None)
    shift_fermi: None or float
        energy shift to be applied to all energies. If `None` no shift is applied. (default: None)
    """

    def __init__(
        self, sisl_hamiltonian, geom, wfsx_sile, spin, ispin=None, shift_fermi=None
    ):
        # super().__init__(geom, spin=spin, ispin=ispin, shift_fermi=shift_fermi)
        super().__init__(sisl_hamiltonian, geom=geom, spin=spin)
        self.geom = geom
        self.wfsx_sile = wfsx_sile
        self.read_all()

    def read_all(self):
        """Read all wavefunctions, eigenvalues, and k-points from WFSX file."""
        evals = []
        evecs = []
        wfsx_kpts = []

        def change_gauge(k, evec):
            """Change the eigenvector gauge"""
            phase = np.dot(
                self.geom.xyz[self.geom.o2a(np.arange(self.geom.no)), :],
                np.dot(k, self.geom.rcell),
            )
            if self.spin.has_noncolinear:
                # for NC/SOC we have a 2x2 spin-box per orbital
                phase = np.repeat(phase, 2)
            # r -> R gauge tranformation.
            return evec * np.exp(1j * phase).reshape(1, -1)

        # Read wavefunctions and eigenvalues
        for wfc in self.wfsx_sile.yield_eigenstate(parent=self.geom):
            wfsx_kpts.append(wfc.info["k"])
            evals.append(wfc.c)
            # To get the same eigvecs as eigh returns we need to transpose the
            # array and change the gauge from 'r' to 'R'
            evecs.append(change_gauge(wfc.info["k"], wfc.state).T)

        # If any k-point occurs more than once in the WaveFuncKPoints block,
        # we discard the duplicates
        is_duplicate = self._is_duplicate(wfsx_kpts)
        self.wfsx_kpts = wfsx_kpts[~is_duplicate]
        self.evals = np.array(evals, dtype=float)[~is_duplicate]
        if self.shift_fermi is not None:
            self.evals += self.shift_fermi
        self.evecs = np.array(evecs, dtype=np.complex64, order="C")[~is_duplicate]

    def _is_duplicate(self, array):
        # TODO: Move into utils
        # Find all matches (i,j): array[i] == array[j]
        matches = np.all(np.isclose(array[None, :], array[:, None]), axis=-1)
        # Remove double counting of matches: (i,j) and (j,i)
        # Also, remove matches of elements with themselves: (i,i)
        matches = np.triu(matches, 1)

        # Finally determine elements which are duplicates
        return np.any(matches, axis=0)

    def find_all(self, kpts):
        """Find the correct eigenvectors and eigenvalues and sort them
        to match the order in kpts.

        Parameters
        ----------
        kpts : list of float (3,) or (nk, 3)
            list of k points

        Returns
        -------
        evals : list of float (nk, n)
            list of eigenvalues for every requested k point
        evecs :
            list of eiegnvector for every requested k point
        """
        kpts = np.asarray(kpts)
        sort_idx = np.where(
            np.all(np.isclose(self.wfsx_kpts[None, :], kpts[:, None]), axis=-1)
        )[1]
        if len(sort_idx) < len(kpts):
            # k-point not found
            raise ValueError(
                f"{self.__class__.__name__}._read_all unable to find at least one "
                "required k point in '{self.wfsx_sile.file}'. Please, ensure that "
                "all k points listed below are included:\n"
                + "\n".join([str(k) for k in kpts])
            )
        if not np.all(np.isclose(self.wfsx_kpts[sort_idx], kpts)):
            # raise ValueError(       wfsx_kpts = np.asarray(wfsx_kpts)
            pass
