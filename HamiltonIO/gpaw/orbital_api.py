#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Parser for GPAW orbital information
"""

from dataclasses import dataclass
from typing import List


@dataclass
class GPAWOrbital:
    """
    Orbital descriptor for GPAW basis functions

    Attributes
    ----------
    iatom : int
        Atom index (0-based)
    element : str
        Element symbol (e.g., "Fe", "O")
    n : int
        Principal quantum number
    l : int
        Angular momentum quantum number (0=s, 1=p, 2=d, 3=f)
    m : int
        Magnetic quantum number (-l, ..., +l)
    z : int
        Zeta function index (for multiple-zeta basis sets)
    M_index : int
        Global orbital index in the Hamiltonian matrix
    """

    iatom: int = 0
    element: str = ""
    n: int = 0
    l: int = 0
    m: int = 0
    z: int = 0
    M_index: int = 0

    @property
    def l_symbol(self) -> str:
        """Return the angular momentum symbol (s, p, d, f)"""
        symbols = ["s", "p", "d", "f"]
        if self.l < len(symbols):
            return symbols[self.l]
        return f"l{self.l}"

    @property
    def sym(self) -> str:
        """Return orbital symbol (e.g., '3d', '2p')"""
        return f"{self.n}{self.l_symbol}"


def parse_gpaw_orbital(calc) -> List[GPAWOrbital]:
    """
    Parse orbital information from GPAW calculator object

    Parameters
    ----------
    calc : GPAW calculator object
        Converged GPAW calculator with LCAO basis set

    Returns
    -------
    List[GPAWOrbital]
        List of GPAWOrbital objects describing each atomic orbital

    Notes
    -----
    GPAW orbital mapping pattern (current API):
    - calc.wfs.basis_functions.M_a[a] gives the starting orbital index for atom a
    - calc.wfs.setups[a].nao gives the number of orbitals for atom a
    - calc.wfs.setups[a].basis_functions_J is a list of Spline objects

    Each Spline has:
    - l: angular momentum quantum number (0=s, 1=p, 2=d, 3=f)

    Note: Current GPAW versions (24.x+) use Spline objects which only have
    angular momentum l, not principal quantum number n. We infer n from
    typical values for each element.
    """
    # Default principal quantum numbers for common elements by l value
    # Format: {element: {l: n}}
    DEFAULT_N_VALUES = {
        "H": {0: 1},  # 1s
        "He": {0: 1},  # 1s
        "Li": {0: 2, 1: 2},  # 2s, 2p
        "Be": {0: 2, 1: 2},  # 2s, 2p
        "B": {0: 2, 1: 2},  # 2s, 2p
        "C": {0: 2, 1: 2},  # 2s, 2p
        "N": {0: 2, 1: 2},  # 2s, 2p
        "O": {0: 2, 1: 2},  # 2s, 2p
        "F": {0: 2, 1: 2},  # 2s, 2p
        "Ne": {0: 2, 1: 2},  # 2s, 2p
        "Na": {0: 3, 1: 3, 2: 3},  # 3s, 3p, 3d
        "Mg": {0: 3, 1: 3, 2: 3},  # 3s, 3p, 3d
        "Al": {0: 3, 1: 3, 2: 3},  # 3s, 3p, 3d
        "Si": {0: 3, 1: 3, 2: 3},  # 3s, 3p, 3d
        "P": {0: 3, 1: 3, 2: 3},  # 3s, 3p, 3d
        "S": {0: 3, 1: 3, 2: 3},  # 3s, 3p, 3d
        "Cl": {0: 3, 1: 3, 2: 3},  # 3s, 3p, 3d
        "Ar": {0: 3, 1: 3, 2: 3},  # 3s, 3p, 3d
        "K": {0: 4, 1: 4, 2: 3},  # 4s, 4p, 3d
        "Ca": {0: 4, 1: 4, 2: 3},  # 4s, 4p, 3d
        "Sc": {0: 4, 1: 4, 2: 3},  # 4s, 4p, 3d
        "Ti": {0: 4, 1: 4, 2: 3},  # 4s, 4p, 3d
        "V": {0: 4, 1: 4, 2: 3},  # 4s, 4p, 3d
        "Cr": {0: 4, 1: 4, 2: 3},  # 4s, 4p, 3d
        "Mn": {0: 4, 1: 4, 2: 3},  # 4s, 4p, 3d
        "Fe": {0: 4, 1: 4, 2: 3},  # 4s, 4p, 3d
        "Co": {0: 4, 1: 4, 2: 3},  # 4s, 4p, 3d
        "Ni": {0: 4, 1: 4, 2: 3},  # 4s, 4p, 3d
        "Cu": {0: 4, 1: 4, 2: 3},  # 4s, 4p, 3d
        "Zn": {0: 4, 1: 4, 2: 3},  # 4s, 4p, 3d
    }

    # Fallback to a reasonable default based on l value
    FALLBACK_N = {0: 3, 1: 3, 2: 3, 3: 4}  # Default to n=3 for s,p,d and n=4 for f

    def get_default_n(element: str, l: int) -> int:
        """Get default principal quantum number for element and l value"""
        if element in DEFAULT_N_VALUES and l in DEFAULT_N_VALUES[element]:
            return DEFAULT_N_VALUES[element][l]
        return FALLBACK_N.get(l, 3)

    orbs = []

    # Get number of atoms from the calculator
    natoms = len(calc.atoms)

    # Iterate over each atom
    for a in range(natoms):
        # Get element symbol
        element = calc.atoms[a].symbol

        # Get starting orbital index for this atom
        # M_a may not exist in some GPAW versions, compute it manually
        if hasattr(calc.wfs.basis_functions, "M_a"):
            M_start = calc.wfs.basis_functions.M_a[a]
        else:
            # Compute M_start manually: sum of nao for previous atoms
            M_start = 0
            for i in range(a):
                M_start += calc.wfs.setups[i].nao

        # Get the Setup object for this atom
        setup = calc.wfs.setups[a]

        # Iterate over basis functions for this atom
        M_offset = 0  # Track orbital offset within this atom
        for z, bf in enumerate(setup.basis_functions_J):
            # In current GPAW, bf is a Spline object with only 'l' attribute
            # The 'n' (principal quantum number) is not stored
            # We infer it from typical values for the element
            n = get_default_n(element, bf.l)

            # For each angular momentum l, there are (2*l+1) m values
            for m in range(-bf.l, bf.l + 1):
                orb = GPAWOrbital(
                    iatom=a,
                    element=element,
                    n=n,
                    l=bf.l,
                    m=m,
                    z=z,
                    M_index=M_start + M_offset,
                )
                orbs.append(orb)
                M_offset += 1

    return orbs


def bset_to_symnum_type(bset, atoms):
    """
    Convert the basis set to symbol number type for compatibility

    Parameters
    ----------
    bset : List[GPAWOrbital]
        List of GPAWOrbital objects
    atoms : ASE Atoms object
        Atomic structure

    Returns
    -------
    List[Tuple]
        List of (atom_index, orbital_symbol, spin) tuples
    """
    from HamiltonIO.utils import symbol_number_list

    slist = symbol_number_list(atoms)
    result = []
    for b in bset:
        result.append((slist[b.iatom], b.sym, "up"))  # Default to spin up
    return result
