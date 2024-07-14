import os
import numpy as np
import copy
from collections import defaultdict
from dataclasses import dataclass


@dataclass
class BaseHamiltonian:
    is_orthogonal: bool = True
    _name: str = "Base Hamiltonian"
    R2kfactor: float = 2 * np.pi
    _nspin: int = 1
    norb: int = 0
    nbasis: int = 0


class Hamiltonian(BaseHamiltonian):
    """
    Abstract class for tight-binding-like Hamiltonian.
    """

    def __init__(
        self,
        R2kfactor=2 * np.pi,
        nspin=1,
        norb=0,
        is_orthogonal=True,
        _name="Generic Hamiltonian",
    ):
        #: :math:`\alpha` used in :math:`H(k)=\sum_R  H(R) \exp( \alpha k \cdot R)`,
        #: Should be :math:`2\pi i` or :math:`-2\pi i`
        self.is_orthogonal = is_orthogonal

        self.R2kfactor = R2kfactor

        #: number of spin. 1 for collinear, 2 for spinor.
        self._nspin = nspin

        #:number of orbitals. Each orbital can have two spins.
        self.norb = norb

        #: nbasis=nspin*norb
        self.nbasis = nspin * norb

        #: The array of cartesian coordinate of all basis. shape:nbasis,3
        self.xcart = None

        #: The array of cartesian coordinate of all basis. shape:nbasis,3
        self.xred = None

        #: The order of the spinor basis.
        #: 1: orb1_up, orb2_up,  ... orb1_down, orb2_down,...
        #: 2: orb1_up, orb1_down, orb2_up, orb2_down,...
        self._name = "Generic Hamiltonian"

    @property
    def Rlist(self):
        """
        list of R vectors.
        """
        return np.asarray(self._Rlist)

    @Rlist.setter
    def Rlist(self, value):
        self._Rlist = value

    @property
    def nR(self):
        """
        number of R vectors.
        """
        return len(self._Rlist)

    @property
    def orb_names(self):
        """
        list of orbital names.
        """
        return self._orb_names

    @property
    def orb_idx(self):
        """
        dictionary of orbital index, the keys are the orbital names, the keys are the indices.
        """
        return self._orb_idx

    @property
    def name(self):
        """
        name of the Hamiltonian. e.g. "Wannier90", "SIESTA"
        """
        return self._name

    def get_HR(self, R=None, iR=None, ispin=0, dense=True):
        """
        get the Hamiltonian H(R), array of shape (nbasis, nbasis)
        parameters:
        =================
        R: array-like
            the R vector.
        iR: int
            the index of the R vector.
        ispin: int
            the index of the spin. For collinear system: 0 for up, 1 for down, None for both.
              For non-spin-polarized system, or non-collinear/SOC system, ispin is ignored.
        dense: bool
            if True, convert to dense matrix if the default format is sparse.

        Returns:
        =================
        H: array of shape (nbasis, nbasis)
        """
        raise NotImplementedError()

    def get_all_HR(self, dense=True, ispin=None):
        """
        get the Hamiltonian H(R) for all R vectors.
        parameters:
        =================
        dense: bool
            if True, convert to dense matrix if the default format is sparse.
        ispin: int or None
            For collinear system: 0 for up, 1 for down, None for both.
            For others, ispin is ignored.
        Returns:
        =================
        H:  array of shape ( nR, nbasis, nbasis)
          if ispin is None, and kspin=2, H.shape=(nR, nbasis, nbasis, 2)

        """
        raise NotImplementedError()

    def get_orbs(self):
        """
        returns the orbitals.
        """
        raise NotImplementedError()

    def HSE_k(self, kpt):
        raise NotImplementedError()

    def HS_and_eigen(self, kpts):
        """
        get Hamiltonian, overlap matrices, eigenvalues, eigen vectors for all kpoints.

        :param:

        * kpts: list of k points.

        :returns:

        * H, S, eigenvalues, eigenvectors for all kpoints
        * H: complex array of shape (nkpts, nbasis, nbasis)
        * S: complex array of shape (nkpts, nbasis, nbasis). S=None if the basis set is orthonormal.
        * evals: complex array of shape (nkpts, nbands)
        * evecs: complex array of shape (nkpts, nbasis, nbands)
        """
        raise NotImplementedError()

    def get_Hk(self, kpt, ispin=None, dense=True):
        """
        get the Hamiltonian H(k), array of shape (nbasis, nbasis)
        parameters:
        =================
        kpt: array-like
            the k vector.
        ispin: int
            the index of the spin. For collinear system: 0 for up, 1 for down, None for both.
              For non-spin-polarized system, or non-collinear/SOC system, ispin is ignored.
        dense: bool
            if True, convert to dense matrix if the default format is sparse.

        Returns:
        =================
        H: array of shape (nbasis, nbasis)
        """
        raise NotImplementedError()

    def get_Sk(self, kpt, ispin=None, dense=True):
        """
        get the overlap matrix S(k), array of shape (nbasis, nbasis)
        parameters:
        =================
        kpt: array-like
            the k vector.
        ispin: int
            the index of the spin. For collinear system: 0 for up, 1 for down, None for both.
              For non-spin-polarized system, or non-collinear/SOC system, ispin is ignored.
        dense: bool
            if True, convert to dense matrix if the default format is sparse.

        Returns:
        =================
        S: array of shape (nbasis, nbasis)
        """
        raise NotImplementedError()

    def get_HS_and_eigen_k(self, kpt, ispin=None, dense=True):
        """
        get Hamiltonian, overlap matrices, eigenvalues, eigen vectors for a kpoint.
        parameters:
        =================
        kpt: array-like
            the k vector.
        ispin: int
            the index of the spin. For collinear system: 0 for up, 1 for down, None for both.
              For non-spin-polarized system, or non-collinear/SOC system, ispin is ignored.
        dense: bool
            if True, convert to dense matrix if the default format is sparse.

        returns:
        =================
        H, S, eigenvalues, eigenvectors
        H: array of shape (nbasis, nbasis)
        S: array of shape (nbasis, nbasis). S=None if the basis set is orthonormal.
        evals: array of shape (nbands,)
        evecs: array of shape (nbasis, nbands)
        """
        raise NotImplementedError()
