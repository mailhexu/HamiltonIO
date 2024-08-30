#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The abacus wrapper
"""

import numpy as np
from scipy.linalg import eigh

from HamiltonIO.hamiltonian import Hamiltonian
from HamiltonIO.mathutils.lowdin import Lowdin_symmetric_orthonormalization
from HamiltonIO.mathutils.rotate_spin import (
    rotate_spinor_matrix_einsum_R,
)
from HamiltonIO.model.kR_convert import R_to_k, R_to_onek
from HamiltonIO.model.occupations import Occupations


class LCAOHamiltonian(Hamiltonian):
    def __init__(
        self,
        HR,
        SR,
        Rlist,
        nbasis=None,
        orbs=None,
        atoms=None,
        nspin=1,
        HR_soc=None,
        HR_nosoc=None,
        nel=None,
        so_strength=1.0,
        orth=False,
    ):
        self.R2kfactor = 2j * np.pi
        self.is_orthogonal = False
        self.split_soc = False
        self._name = "LCAOHamiltonian"
        self._HR = HR
        self.SR = SR
        self.Rlist = Rlist
        self.nbasis = nbasis
        self.orbs = orbs
        self.atoms = atoms
        self.nspin = nspin
        self.norb = nbasis * nspin
        self.nel = nel
        if HR_soc is not None:
            self.set_HR_soc(HR_soc=HR_soc, HR_nosoc=HR_nosoc, HR_full=HR)
        self.soc_rotation_angle = [0.0, 0.0]
        self.so_strength = so_strength
        self.orth = orth
        if orth:
            self.is_orthogonal = True
        self.H0 = None

    @property
    def Rlist(self):
        return self._Rlist

    @Rlist.setter
    def Rlist(self, Rlist):
        self._Rlist = Rlist
        self._build_Rdict()

    def get_Ridx(self, R):
        """
        Get the index of R in the Rlist
        """
        return self.Rdict[tuple(R)]

    # def get_H0(self):
    #    return self._get_H0(*self.soc_rotation_angle)

    def get_H0(self):
        R0 = self.get_Ridx((0, 0, 0))
        if self.split_soc:
            if self.H0 is None:
                theta, phi = self.soc_rotation_angle
                self.H0 = rotate_spinor_matrix_einsum_R(self.HR_nosoc, theta, phi)[R0]
            return self.H0
        else:
            return self.HR[R0]

    def get_HR_soc(self, R):
        return self.HR_soc[self.get_Ridx(R)]

    def get_Hk_soc(self, kpts):
        return R_to_k(kpts, self.Rlist, self.HR_soc)

    def set_HR_soc(self, HR_soc=None, HR_nosoc=None, HR_full=None):
        self.split_soc = True
        self.HR_soc = HR_soc
        self.HR_nosoc = HR_nosoc
        self.HR_full = HR_full
        if HR_soc is None:
            self.HR_soc = HR_full - HR_nosoc
        if HR_nosoc is None:
            self.HR_nosoc = HR_full - HR_soc
        if HR_full is None:
            self.HR_full = HR_soc + HR_nosoc

    def set_so_strength(self, so_strength):
        self.so_strength = so_strength

    def set_Hsoc_rotation_angle(self, angle):
        """
        Set the rotation angle for SOC part of Hamiltonian
        """
        self.soc_rotation_angle = angle

    @property
    def HR(self):
        if self.split_soc:
            theta, phi = self.soc_rotation_angle
            _HR = self.get_HR_from_soc(theta, phi, self.so_strength)
            return _HR
        else:
            return self._HR

    def get_HR_from_soc(self, theta, phi, so_strength):
        HR = (
            rotate_spinor_matrix_einsum_R(self.HR_nosoc, theta, phi)
            + self.HR_soc * so_strength
        )
        return HR

    @HR.setter
    def set_HR(self, HR):
        self._HR = HR

    def _build_Rdict(self):
        if hasattr(self, "Rdict"):
            pass
        else:
            self.Rdict = {}
            for iR, R in enumerate(self.Rlist):
                self.Rdict[tuple(R)] = iR

    def get_hamR(self, R):
        return self.HR[self.Rdict[tuple(R)]]

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
        # Hk = np.zeros((self.nbasis, self.nbasis), dtype="complex")
        # Sk = np.zeros((self.nbasis, self.nbasis), dtype="complex")
        if convention == 2:
            # for iR, R in enumerate(self.Rlist):
            #    phase = np.exp(self.R2kfactor * np.dot(k, R))
            #    H = self.HR[iR] * phase
            #    # Hk += H + H.conjugate().T
            #    Hk += H
            #    S = self.SR[iR] * phase
            #    # Sk += S + S.conjugate().T
            #    Sk += S
            #    # Hk = (Hk + Hk.conj().T)/2
            #    # Sk = (Sk + Sk.conj().T)/2
            Hk = R_to_onek(k, self.Rlist, self.HR)
            Sk = R_to_onek(k, self.Rlist, self.SR)
            if self.orth:
                Hk = Lowdin_symmetric_orthonormalization(Hk, Sk)
                Sk = None
        elif convention == 1:
            # TODO: implement the first convention (the r convention)
            raise NotImplementedError("convention 1 is not implemented yet.")
            pass
        else:
            raise ValueError("convention should be either 1 or 2.")
        return Hk, Sk

    def solve(self, k, convention=2):
        Hk, Sk = self.gen_ham(k, convention=convention)
        return eigh(Hk, Sk)

    def solve_all(self, kpts, convention=2):
        nk = len(kpts)
        evals = np.zeros((nk, self.nbasis), dtype=float)
        evecs = np.zeros((nk, self.nbasis, self.nbasis), dtype=complex)
        for ik, k in enumerate(kpts):
            evals[ik], evecs[ik] = self.solve(k, convention=convention)
        return evals, evecs

    def HSE_k(self, kpt, convention=2):
        H, S = self.gen_ham(tuple(kpt), convention=convention)
        evals, evecs = eigh(H, S)
        return H, S, evals, evecs

    def HS_and_eigen(self, kpts, convention=2):
        """
        calculate eigens for all kpoints.
        :param kpts: list of k points.
        """
        nk = len(kpts)
        hams = np.zeros((nk, self.nbasis, self.nbasis), dtype=complex)
        if self.orth:
            Ss = None
        else:
            Ss = np.zeros((nk, self.nbasis, self.nbasis), dtype=complex)
        evals = np.zeros((nk, self.nbasis), dtype=float)
        evecs = np.zeros((nk, self.nbasis, self.nbasis), dtype=complex)
        if self.orth:
            for ik, k in enumerate(kpts):
                hams[ik], _, evals[ik], evecs[ik] = self.HSE_k(
                    tuple(k), convention=convention
                )
        else:
            for ik, k in enumerate(kpts):
                hams[ik], Ss[ik], evals[ik], evecs[ik] = self.HSE_k(
                    tuple(k), convention=convention
                )
        return hams, Ss, evals, evecs

    def get_fermi_energy(self, evals, width=0.01, kweights=None, nspin=2):
        occ = Occupations(nel=self.nel, wk=kweights, nspin=nspin)
        efermi = occ.efermi(evals)
        return efermi
