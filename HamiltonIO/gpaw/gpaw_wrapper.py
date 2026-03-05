#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GPAW Wrapper for HamiltonIO
"""

from typing import Optional, Tuple, Union

import numpy as np

from HamiltonIO.gpaw.gpaw_api import (
    convert_k_to_r_space_manual,
    get_hr_sr_from_calc,
    read_gpaw_calculator,
    read_gpaw_pickle,
)
from HamiltonIO.gpaw.orbital_api import GPAWOrbital, parse_gpaw_orbital
from HamiltonIO.lcao_hamiltonian import LCAOHamiltonian


class GPAWWrapper(LCAOHamiltonian):
    """
    GPAW wrapper class extending LCAOHamiltonian

    This class provides an interface between GPAW's LCAO output and TB2J.

    Parameters
    ----------
    HR : np.ndarray
        Real-space Hamiltonian matrices
        Shape: (nR, nao, nao) for nspin=1
               (nR, nao, nao) for nspin=2 (model_up or model_dn)
    SR : np.ndarray
        Real-space Overlap matrices
        Shape: (nR, nao, nao)
    Rlist : np.ndarray
        List of R vectors (cell indices)
        Shape: (nR, 3)
    nbasis : int
        Number of atomic orbitals
    atoms : ASE Atoms object, optional
        Atomic structure
    nspin : int, optional
        Spin polarization (1=non-polarized, 2=collinear)
    orbs : List[GPAWOrbital], optional
        List of orbital descriptors
    nel : float, optional
        Number of electrons
    """

    def __init__(
        self,
        HR: np.ndarray,
        SR: np.ndarray,
        Rlist: np.ndarray,
        nbasis: int,
        atoms=None,
        nspin: int = 1,
        orbs: Optional[list] = None,
        nel: Optional[float] = None,
    ):
        super().__init__(
            HR=HR,
            SR=SR,
            Rlist=Rlist,
            nbasis=nbasis,
            atoms=atoms,
            nspin=nspin,
            orbs=orbs,
        )
        self._name = "GPAW"
        self.R2kfactor = 2j * np.pi
        self._build_Rdict()
        self._positions = None
        # Store k-space data for direct access (if available)
        self.H_k = None  # k-space Hamiltonian
        self.S_k = None  # k-space Overlap
        self.kpts = None  # k-points list

    @property
    def positions(self):
        """
        Get orbital positions (in Angstroms)

        For GPAW LCAO, orbitals are centered on atoms.
        Each atom has multiple orbitals at the same position.

        Returns
        -------
        np.ndarray
            Orbital positions, shape: (nbasis, 3)
        """
        if self._positions is None:
            if self.atoms is not None:
                # For LCAO, all orbitals on the same atom share the atom's position
                # We need to map each orbital to its atom's position
                natoms = len(self.atoms)
                atom_positions = self.atoms.get_positions()

                if self.orbs is not None:
                    # Use orbital information to map orbitals to atom positions
                    self._positions = np.array(
                        [atom_positions[orb.iatom] for orb in self.orbs]
                    )
                else:
                    # Fallback: assume equal number of orbitals per atom
                    # This is approximate but should work for basic cases
                    norbs_per_atom = self.nbasis // natoms
                    if self.nbasis % natoms != 0:
                        norbs_per_atom = 1  # Minimum assumption
                    self._positions = np.repeat(atom_positions, norbs_per_atom, axis=0)[
                        : self.nbasis
                    ]
            else:
                # No atoms data, return zeros
                self._positions = np.zeros((self.nbasis, 3))

        return self._positions

    def _build_Rdict(self):
        """
        Build dictionary mapping R vectors to indices

        Creates self.Rdict such that self.Rdict[tuple(R)] = iR
        """
        self.Rdict = {}
        for iR, R in enumerate(self.Rlist):
            self.Rdict[tuple(R)] = iR

    def get_hamR(self, R: Union[tuple, list, np.ndarray]) -> np.ndarray:
        """
        Get the Hamiltonian matrix for a given R vector

        Parameters
        ----------
        R : tuple, list, or np.ndarray
            R vector (3 integers)

        Returns
        -------
        np.ndarray
            Hamiltonian matrix H[R]
            Shape: (nao, nao)
        """
        return self.HR[self.Rdict[tuple(R)]]

    def set_k_data(self, H_k, S_k, kpts):
        """
        Set k-space Hamiltonian and Overlap data directly

        This is used when the data is already in k-space format
        (e.g., from GPAW's get_lcao_hamiltonian).

        Parameters
        ----------
        H_k : np.ndarray
            k-space Hamiltonian
            Shape: (nspins, nkpts, nao, nao) or (nkpts, nao, nao)
        S_k : np.ndarray
            k-space Overlap
            Shape: (nkpts, nao, nao)
        kpts : np.ndarray
            k-points in fractional coordinates
            Shape: (nkpts, 3)
        """
        self.H_k = H_k
        self.S_k = S_k
        self.kpts = kpts

    def gen_ham(self, k, convention=2):
        """
        Generate Hamiltonian matrix at k point.

        For GPAW, if k-space data is available (H_k), return the Hamiltonian
        for the closest k-point. Otherwise, use the standard R→k transform.

        Parameters
        ----------
        k : array-like
            k-point (3 fractional coordinates)
        convention : int, optional
            Convention for phase factor (default: 2)

        Returns
        -------
        Hk : np.ndarray
            k-space Hamiltonian at k-point
            Shape: (nao, nao)
        Sk : np.ndarray or None
            k-space Overlap at k-point (None if orthogonal)
            Shape: (nao, nao)
        """
        k = np.asarray(k)

        # Check if we have k-space data
        if self.H_k is not None and self.kpts is not None:
            # Find closest k-point in our list
            distances = np.linalg.norm(self.kpts - k, axis=1)
            ik = np.argmin(distances)

            # For LCAO k-space data, we return the closest k-point
            if self.H_k.ndim == 3 and self.H_k.shape[0] == len(self.kpts):
                # H_k has shape (nkpts, nao, nao)
                Hk = self.H_k[ik]
            else:
                # Fallback: try to handle other shapes
                Hk = self.H_k[ik]

            Sk = self.S_k[ik] if self.S_k is not None else None

            if self.orth:
                from HamiltonIO.mathutils.lowdin import (
                    Lowdin_symmetric_orthonormalization,
                )

                Hk = Lowdin_symmetric_orthonormalization(Hk, Sk)
                Sk = None

            return Hk, Sk

        # Fall back to standard R→k transform
        return super().gen_ham(k, convention)


class GPAWParser:
    """
    Parser for GPAW LCAO calculations

    This class coordinates reading GPAW calculation output and creating
    GPAWWrapper instances.

    Parameters
    ----------
    gpw_file : str, optional
        Path to the .gpw file containing the converged GPAW calculation
    calc : GPAW calculator object, optional
        Pre-loaded GPAW calculator object
    pickle_file : str, optional
        Path to the pickle file created by dump_hamiltonian()
    """

    def __init__(
        self,
        gpw_file: Optional[str] = None,
        calc=None,
        pickle_file: Optional[str] = None,
    ):
        self.gpw_file = gpw_file
        self.calc = calc
        self.pickle_file = pickle_file
        self.atoms = None
        self.basis = None
        self.efermi = None
        self.nel = None

        # Read data
        if calc is not None:
            self.calc = calc
        elif gpw_file is not None:
            self.read_calc(gpw_file)
        elif pickle_file is not None:
            self.read_pickle(pickle_file)
        else:
            raise ValueError("Must provide either calc, gpw_file, or pickle_file")

    def read_calc(self, gpw_file: str):
        """
        Load GPAW calculator from .gpw file

        Parameters
        ----------
        gpw_file : str
            Path to the .gpw file
        """
        self.calc = read_gpaw_calculator(gpw_file)
        self.gpw_file = gpw_file
        self.atoms = self.calc.atoms
        return self.calc

    def read_pickle(self, pickle_file: str):
        """
        Load data from pickle file

        Parameters
        ----------
        pickle_file : str
            Path to the pickle file
        """
        self.pickle_data = read_gpaw_pickle(pickle_file)
        self.pickle_file = pickle_file
        # Extract available data from pickle
        if "H_skMM" in self.pickle_data:
            self.H_skMM = self.pickle_data["H_skMM"]
        if "S_kMM" in self.pickle_data:
            self.S_kMM = self.pickle_data["S_kMM"]
        if "atoms_data" in self.pickle_data:
            # Reconstruct atoms object from atoms_data
            from ase import Atoms

            self.atoms = Atoms(**self.pickle_data["atoms_data"])

    def read_basis(self) -> list:
        """
        Parse orbital information from GPAW calculator or reconstruct from atoms

        Returns
        -------
        List[GPAWOrbital]
            List of orbital descriptors
        """
        if self.calc is not None:
            self.basis = parse_gpaw_orbital(self.calc)
        elif self.atoms is not None:
            # Reconstruct basis from atoms (simplified version)
            # We need to determine the number of orbitals per atom
            # This depends on the basis set used in GPAW

            # For LCAO, GPAW uses atomic orbitals based on the setup
            # We'll create a simple basis based on atomic numbers
            # This is a placeholder - real implementation needs GPAW setup info

            from ase.data import atomic_numbers

            symbols = self.atoms.get_chemical_symbols()

            # Simple mapping: count orbitals per element
            # This should be replaced with actual GPAW basis info
            element_orbitals = {
                "H": [("s", 0)],
                "He": [("s", 0)],
                "Li": [("s", 0)],
                "Be": [("s", 0)],
                "B": [("s", 0), ("p", 0)],
                "C": [("s", 0), ("p", 0)],
                "N": [("s", 0), ("p", 0)],
                "O": [("s", 0), ("p", 0)],
                "Fe": [("s", 0), ("p", 0), ("d", 0)],
                # Add more elements as needed
            }

            basis = []
            for iatom, symbol in enumerate(symbols):
                if symbol in element_orbitals:
                    for orb_name, m in element_orbitals[symbol]:
                        orb = GPAWOrbital(
                            iatom=iatom,
                            element=symbol,
                            n=0,
                            l=0,
                            m=m,
                            z=atomic_numbers[symbol],
                            M_index=len(basis),
                        )
                        basis.append(orb)
                else:
                    # Default to s orbital for unknown elements
                    orb = GPAWOrbital(
                        iatom=iatom,
                        element=symbol,
                        n=0,
                        l=0,
                        m=0,
                        z=atomic_numbers[symbol],
                        M_index=len(basis),
                    )
                    basis.append(orb)

            self.basis = basis
        else:
            raise ValueError(
                "Calculator object or atoms data not available for reading basis"
            )
        return self.basis

    def read_efermi(self) -> float:
        """
        Read Fermi energy from GPAW calculator or pickle

        Returns
        -------
        float
            Fermi energy in eV
        """
        if self.calc is not None:
            self.efermi = self.calc.get_fermi_level()
        elif hasattr(self, "pickle_data") and "calc_data" in self.pickle_data:
            # For pickle files, we may not have efermi stored
            # Calculate from H eigenvalues or use default
            # For now, return 0.0 as placeholder
            self.efermi = 0.0
        else:
            raise ValueError(
                "Calculator object or pickle data not available for reading efermi"
            )
        return self.efermi

    def read_nel(self) -> float:
        """
        Read number of electrons from GPAW calculator or pickle

        Returns
        -------
        float
            Number of electrons
        """
        if self.calc is not None:
            self.nel = self.calc.get_number_of_electrons()
        elif hasattr(self, "pickle_data") and self.atoms is not None:
            # Calculate from atomic numbers
            from ase.data import atomic_numbers

            self.nel = sum(
                atomic_numbers[symbol] for symbol in self.atoms.get_chemical_symbols()
            )
        else:
            raise ValueError(
                "Calculator object or atoms data not available for reading nel"
            )
        return self.nel

    def read_hsr(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray, int]:
        """
        Read Hamiltonian and Overlap matrices from GPAW calculator or pickle

        Returns
        -------
        H_skMM : np.ndarray
            k-space Hamiltonian in eV
        S_kMM : np.ndarray
            k-space Overlap matrix
        Rlist : np.ndarray
            List of R vectors (placeholder for k-space)
        nbasis : int
            Number of atomic orbitals
        """
        # Check if we have pickle data first
        if hasattr(self, "pickle_data") and "H_skMM" in self.pickle_data:
            H_skMM = self.pickle_data["H_skMM"]
            S_kMM = self.pickle_data["S_kMM"]
        elif self.calc is None:
            raise ValueError(
                "Calculator object or pickle data not available for reading HSR"
            )
        else:
            H_skMM, S_kMM = get_hr_sr_from_calc(self.calc)

        # Determine nbasis
        if H_skMM.ndim == 4:
            _, _, nbasis, _ = H_skMM.shape
        else:
            _, nbasis, _ = H_skMM.shape

        # For k-space, Rlist is just a placeholder (Gamma point only)
        Rlist = np.array([[0, 0, 0]])

        return H_skMM, S_kMM, Rlist, nbasis

    def get_models(self) -> Union[GPAWWrapper, Tuple[GPAWWrapper, GPAWWrapper]]:
        """
        Create GPAWWrapper instance(s) from the parsed data

        Returns
        -------
        model : GPAWWrapper or (model_up, model_dn)
            For nspin=1: single GPAWWrapper
            For nspin=2: tuple of (model_up, model_dn)

        Notes
        -----
        TB2J uses real-space hopping matrices (HR, SR, Rlist) directly for
        exchange calculation. Therefore, we perform Fourier transform using
        GPAW's TightBinding class when gpw_file is provided.

        For accurate results:
        - Use dense k-mesh (≥7×7×7)
        - Disable point group symmetry: symmetry={'point_group': False}
        - This ensures IBZ ≈ full BZ for proper Fourier transform
        """
        # Read all necessary data
        self.read_basis()
        self.read_efermi()
        self.read_nel()

        # Get k-space Hamiltonian and Overlap
        H_skMM, S_kMM, Rlist, nbasis = self.read_hsr()

        # Determine spin
        if H_skMM.ndim == 4:
            nspin = 2  # Spin-polarized
            nspins, nkpts, nao, _ = H_skMM.shape
        else:
            nspin = 1  # Non-polarized
            nkpts, nao, _ = H_skMM.shape

        # Get k-points
        if hasattr(self, "pickle_data") and "calc_data" in self.pickle_data:
            ibzk_qc = self.pickle_data["calc_data"].get("ibzk_kc")
        elif self.calc is not None:
            ibzk_qc = self.calc.get_ibz_k_points()
        else:
            raise ValueError(
                "K-point data not available. Please provide either gpw_file or "
                "pickle_file with calc_data containing k-point data."
            )

        # Use manual Fourier transform if calc is available
        # TB2J needs real-space matrices for exchange calculation
        # Note: GPAW's TightBinding.bloch_to_real_space uses time-reversal symmetry
        # which TB2J doesn't use, so we use manual transform instead
        if self.calc is not None:
            print("Using manual k→R Fourier transform (no time-reversal symmetry)...")

            # Get k-point data for manual transform
            ibzk_qc = self.calc.get_ibz_k_points()
            # Get k-point weights (not used in transform, but needed for completeness)
            weight_k = self.calc.wfs.kd.weight_k
            cell = self.atoms.cell

            # Determine Rmax based on number of k-points
            # For proper Fourier transform, we need nR ≈ nkpts
            nkpts = len(ibzk_qc)
            # Find Rmax such that (2*Rmax+1)^3 ≈ nkpts
            import math

            Rmax = int(math.ceil((nkpts ** (1 / 3) - 1) / 2))
            print(f"  Auto-detected Rmax={Rmax} for {nkpts} k-points")

            HR, SR, Rlist, nbasis = convert_k_to_r_space_manual(
                H_skMM, S_kMM, ibzk_qc, weight_k, cell, Rmax=Rmax
            )
            print(f"  Generated {len(Rlist)} R vectors")
            print(f"  HR shape: {HR.shape}")
            print(f"  SR shape: {SR.shape}")

            # Verify matrices are valid
            # S at R=0 should be positive definite
            iR0 = np.argmin(np.linalg.norm(Rlist, axis=1))
            S_evals = np.linalg.eigvalsh(SR[iR0].real)
            if not np.all(S_evals > 0):
                raise ValueError(
                    f"Overlap matrix at R=0 is not positive definite. "
                    f"Min eigenvalue: {S_evals.min():.2e}"
                )
            print(
                f"  S at R=0 is positive definite (min eigenvalue: {S_evals.min():.2e})"
            )

            # Clean up small imaginary parts from numerical noise
            HR = np.real_if_close(HR, tol=1e-8)
            SR = np.real_if_close(SR, tol=1e-8)

            # Symmetrize HR matrices to restore crystal symmetry
            # This is necessary because with symmetry='off', the Fourier transform
            # does not preserve symmetry equivalence between HR matrices
            try:
                from HamiltonIO.gpaw.symmetrize_hr import symmetrize_hr

                Rlist_tuples = [tuple(R) for R in Rlist]

                if HR.ndim == 4:
                    # Spin-polarized case: symmetrize each spin channel separately
                    for ispin in range(HR.shape[0]):
                        HR[ispin] = symmetrize_hr(
                            HR[ispin], Rlist_tuples, self.atoms, verbose=False
                        )
                else:
                    # Non-polarized case
                    HR = symmetrize_hr(HR, Rlist_tuples, self.atoms, verbose=False)

                print("  HR matrices symmetrized using crystal symmetry")
            except ImportError:
                print(
                    "  Warning: symmetrization module not available. HR matrices not symmetrized."
                )
                print("  Exchange parameters may have broken symmetry equivalence.")

            # Create wrappers with real-space matrices
            if nspin == 2:
                # HR has shape (nspins, nR, nao, nao)
                model_up = GPAWWrapper(
                    HR=HR[0],
                    SR=SR,
                    Rlist=Rlist,
                    nbasis=nbasis,
                    atoms=self.atoms,
                    nspin=1,
                    orbs=self.basis,
                    nel=self.nel,
                )
                model_dn = GPAWWrapper(
                    HR=HR[1],
                    SR=SR,
                    Rlist=Rlist,
                    nbasis=nbasis,
                    atoms=self.atoms,
                    nspin=1,
                    orbs=self.basis,
                    nel=self.nel,
                )

                # Also set k-space data for gen_ham interpolation
                # This ensures gen_ham uses k-space data instead of R→k transform
                model_up.set_k_data(H_skMM[0], S_kMM, ibzk_qc)
                model_dn.set_k_data(H_skMM[1], S_kMM, ibzk_qc)

                model_up.efermi = self.efermi
                model_dn.efermi = self.efermi
                return model_up, model_dn
            else:
                model = GPAWWrapper(
                    HR=HR,
                    SR=SR,
                    Rlist=Rlist,
                    nbasis=nbasis,
                    atoms=self.atoms,
                    nspin=nspin,
                    orbs=self.basis,
                    nel=self.nel,
                )

                # Also set k-space data for gen_ham interpolation
                model.set_k_data(H_skMM, S_kMM, ibzk_qc)

                model.efermi = self.efermi
                return model

        # Fallback: Use k-space interpolation if no calc available
        print("Warning: Using k-space interpolation (calc not available).")
        print("Exchange parameters may be less accurate.")
        HR = np.zeros((1, nao, nao), dtype=complex)
        SR = np.zeros((1, nao, nao), dtype=complex)

        if nspin == 2:
            HR[0] = np.mean(H_skMM, axis=1).mean(axis=0)
            SR[0] = np.mean(S_kMM, axis=0)

            model_up = GPAWWrapper(
                HR=HR,
                SR=SR,
                Rlist=np.array([[0, 0, 0]]),
                nbasis=nbasis,
                atoms=self.atoms,
                nspin=1,
                orbs=self.basis,
                nel=self.nel,
            )
            model_dn = GPAWWrapper(
                HR=HR,
                SR=SR,
                Rlist=np.array([[0, 0, 0]]),
                nbasis=nbasis,
                atoms=self.atoms,
                nspin=1,
                orbs=self.basis,
                nel=self.nel,
            )

            model_up.set_k_data(H_skMM[0], S_kMM, ibzk_qc)
            model_dn.set_k_data(H_skMM[1], S_kMM, ibzk_qc)

            model_up.efermi = self.efermi
            model_dn.efermi = self.efermi
            return model_up, model_dn
        else:
            HR[0] = np.mean(H_skMM, axis=0)
            SR[0] = np.mean(S_kMM, axis=0)

            model = GPAWWrapper(
                HR=HR,
                SR=SR,
                Rlist=np.array([[0, 0, 0]]),
                nbasis=nbasis,
                atoms=self.atoms,
                nspin=nspin,
                orbs=self.basis,
                nel=self.nel,
            )

            model.set_k_data(H_skMM, S_kMM, ibzk_qc)

            model.efermi = self.efermi
            return model

    def get_basis(self):
        """
        Get basis information in symbol-number format

        Returns
        -------
        list or tuple
            For nspin=1: list of (atom_index, orbital_symbol, spin) tuples
            For nspin=2: tuple of (basis_up, basis_dn)
        """
        from HamiltonIO.utils import symbol_number_list

        if self.basis is None:
            self.read_basis()

        slist = symbol_number_list(self.atoms)

        # Determine spin from calculator
        if self.calc is not None:
            nspin = self.calc.wfs.nspins
        elif hasattr(self, "H_skMM") and self.H_skMM.ndim == 4:
            nspin = 2
        else:
            nspin = 1

        if nspin == 2:
            basis_up = []
            basis_dn = []
            for b in self.basis:
                basis_up.append((slist[b.iatom], b.sym, "up"))
                basis_dn.append((slist[b.iatom], b.sym, "down"))
            return basis_up, basis_dn
        else:
            basis = []
            for b in self.basis:
                basis.append((slist[b.iatom], b.sym, "up"))
            return basis
