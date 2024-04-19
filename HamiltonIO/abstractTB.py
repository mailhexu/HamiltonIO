import os
import numpy as np
import copy
from scipy.linalg import eigh
from scipy.sparse import csr_matrix
from scipy.io import netcdf_file
from collections import defaultdict


class AbstractTB:
    def __init__(self, R2kfactor, nspin, norb):
        #: :math:`\alpha` used in :math:`H(k)=\sum_R  H(R) \exp( \alpha k \cdot R)`,
        #: Should be :math:`2\pi i` or :math:`-2\pi i`
        self.is_siesta = False
        self.is_orthogonal = True
        self.R2kfactor = R2kfactor

        #: number of spin. 1 for collinear, 2 for spinor.
        self.nspin = nspin

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

        self._name = None

    @property
    def name(self):
        return self._name

    def get_hamR(self, R):
        """
        get the Hamiltonian H(R), array of shape (nbasis, nbasis)
        """
        raise NotImplementedError()

    def get_orbs(self):
        """
        returns the orbitals.
        """
        raise NotImplementedError()

    def HSE(self, kpt):
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



