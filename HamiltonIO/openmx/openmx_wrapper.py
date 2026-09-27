"""OpenMX ``.scfout`` model: exposes HR/SR tables and ``H(k)`` solvers.

Mirrors the SIESTA/ABACUS wrappers (subclass of
:class:`~HamiltonIO.lcao_hamiltonian.LCAOHamiltonian`): the real-space
tables are ``(nR, norb, norb)`` arrays over the integer image
translations ``Rlist`` and the k-space assembly is convention 2,
``M(k) = sum_R exp(+2 pi i k.R) M(R)`` with ``k`` in fractional
reciprocal coordinates. This is exactly OpenMX's own phase
(``EigenValue_Problem.c``: ``exp(i 2 pi k.(l1,l2,l3))`` with the same
``atv_ijk`` translations), so eigenvalues match OpenMX band output.
"""

from __future__ import annotations

import numpy as np
from ase.atoms import Atoms
from ase.units import Bohr, Ha

from HamiltonIO.lcao_hamiltonian import LCAOHamiltonian

from .scfout import Scfout, read_scfout


class OpenmxWrapper(LCAOHamiltonian):
    """HamiltonIO model for one spin channel of an OpenMX calculation.

    For unpolarized runs the single wrapper carries the (norb, norb)
    matrices; collinear runs yield one wrapper per channel with
    ``nspin=1``; non-collinear runs yield one wrapper with the spinor
    dimension ``2*norb``.
    """

    def __init__(
        self,
        HR,
        SR,
        Rlist,
        nbasis=None,
        atoms=None,
        nspin=1,
        orbs=None,
        nel=None,
        efermi=None,
        data=None,
    ):
        super().__init__(
            HR=HR,
            SR=SR,
            Rlist=Rlist,
            nbasis=nbasis,
            atoms=atoms,
            orbs=orbs,
            nspin=nspin,
            nel=nel,
        )
        self._name = "OpenMX"
        self.efermi = efermi
        self.data = data


class OpenmxParser:
    """Parser for OpenMX ``.scfout`` files.

    Parameters
    ----------
    scfout : path-like
        The ``.scfout`` binary (Hartree/Bohr inside, converted here).
    atoms : ase.Atoms, optional
        Override the structure reconstructed from the file.
    """

    def __init__(self, scfout=None, atoms=None):
        if scfout is None:
            raise ValueError("scfout path is required")
        self.scfout_path = scfout
        self.data: Scfout = read_scfout(scfout)
        self._atoms = atoms

    # -- structure ----------------------------------------------------------

    @property
    def atoms(self) -> Atoms:
        """Supercell atoms (the OpenMX cell) with cell in Angstrom."""
        if self._atoms is None:
            data = self.data
            self._atoms = Atoms(
                symbols=data.symbols(),
                positions=data.positions * Bohr,
                cell=data.tv * Bohr,
                pbc=True,
            )
        return self._atoms

    @property
    def efermi(self) -> float:
        """Fermi level in eV (OpenMX ``ChemP`` in Hartree)."""
        return self.data.chem_p * Ha

    @property
    def n_electrons(self) -> float:
        return self.data.valence_electrons

    # -- models -------------------------------------------------------------

    def get_model(self, spin=None) -> OpenmxWrapper:
        """Build the HamiltonIO model.

        ``spin=None`` returns the only model (spin 0 and spin 3). For
        collinear runs pass ``spin=0``/``spin=1`` to select the channel;
        the returned model is single-channel (``nspin=1``), mirroring
        the ABACUS/SIESTA wrappers.
        """
        data = self.data
        Rlist = data.Rlist
        SR = data.S
        eV_per_Ha = float(Ha)
        if data.spin_p_switch == 1:
            if spin not in (0, 1):
                raise ValueError(
                    "collinear OpenMX data: pass spin=0 (up) or spin=1 (down)"
                )
            model = OpenmxWrapper(
                HR=np.asarray(data.H[spin]) * eV_per_Ha,
                SR=SR,
                Rlist=Rlist,
                nbasis=data.norb,
                atoms=self.atoms,
                nspin=1,
                nel=data.valence_electrons / 2,
                efermi=self.efermi,
                data=data,
            )
            return model
        if data.spin_p_switch != 1 and spin is not None:
            raise ValueError(
                "only collinear (SpinPolarization on) OpenMX data has spin "
                "channels; pass spin=None"
            )
        return OpenmxWrapper(
            HR=np.asarray(data.H) * eV_per_Ha,
            SR=SR,
            Rlist=Rlist,
            nbasis=data.norb * (2 if data.spin_p_switch == 3 else 1),
            atoms=self.atoms,
            nspin=1,
            nel=data.valence_electrons,
            efermi=self.efermi,
            data=data,
        )

    def get_models(self):
        """All models: one for spin 0/3, ``(up, down)`` for collinear."""
        if self.data.spin_p_switch == 1:
            return self.get_model(0), self.get_model(1)
        return self.get_model()
