#!/usr/bin/env python3
"""
Symmetrize HR matrices using crystal symmetry operations.

This module provides functions to symmetrize real-space Hamiltonian matrices
by averaging over symmetry-equivalent (i, j, R) combinations.
"""

from typing import List, Tuple

import numpy as np

try:
    from spglib import get_symmetry_dataset

    HAS_SPGLIB = True
except ImportError:
    HAS_SPGLIB = False
    print("Warning: spglib not available. Symmetrization will not work.")


def get_symmetry_operations(atoms) -> Tuple[List[np.ndarray], List[np.ndarray]]:
    """
    Get symmetry operations for a crystal structure.

    Parameters
    ----------
    atoms : ase.Atoms
        Crystal structure

    Returns
    -------
    rotations : list of np.ndarray
        Rotation matrices (3x3) for each symmetry operation
    translations : list of np.ndarray
        Translation vectors (3,) for each symmetry operation
    """
    if not HAS_SPGLIB:
        return [], []

    # spglib expects: (lattice, positions, numbers)
    # NOTE: positions must be in fractional coordinates!
    lattice = atoms.get_cell().T
    # Convert to fractional coordinates if needed
    cell = atoms.get_cell()
    positions_cart = atoms.get_positions()
    positions_frac = np.dot(positions_cart, np.linalg.inv(cell))
    numbers = atoms.get_atomic_numbers()

    dataset = get_symmetry_dataset((lattice, positions_frac, numbers))

    rotations = dataset.rotations
    translations = dataset.translations

    return rotations, translations


def apply_symmetry_to_bond(
    iatom: int,
    jatom: int,
    R: Tuple[int, int, int],
    rotation: np.ndarray,
    translation: np.ndarray,
    atoms,
    pos_frac: np.ndarray = None,
) -> List[Tuple[int, int, Tuple[int, int, int]]]:
    """
    Apply a symmetry operation to a bond (iatom, jatom, R).

    IMPORTANT: Only returns bonds where the atom mapping preserves the bond type.
    For a bond (i, j, R), we require that atom i maps to atom i' and atom j maps to atom j'.
    Self-bonds (i==j) will only map to other self-bonds, and inter-atomic bonds (i!=j)
    will only map to other inter-atomic bonds.

    Parameters
    ----------
    iatom : int
        Index of first atom
    jatom : int
        Index of second atom
    R : tuple of 3 ints
        Unit cell translation vector
    rotation : np.ndarray
        3x3 rotation matrix
    translation : np.ndarray
        3-element translation vector (in fractional coordinates)
    atoms : ase.Atoms
        Crystal structure
    pos_frac : np.ndarray, optional
        Pre-computed positions in fractional coordinates

    Returns
    -------
    equivalent_bonds : list of tuples
        List of (iatom', jatom', R') that are equivalent to the input bond
    """
    cell = atoms.get_cell()
    natoms = len(atoms)

    # Get positions in fractional coordinates
    if pos_frac is None:
        positions = atoms.get_positions()
        pos_frac = np.dot(positions, np.linalg.inv(cell))

    # Position of atom i
    pos_i = pos_frac[iatom]

    # Position of atom j in cell R
    pos_j_R = pos_frac[jatom] + np.array(R)

    # Apply rotation to the relative position
    rel_pos = pos_j_R - pos_i
    rel_pos_rotated = np.dot(rotation, rel_pos)

    # Position of atom i after rotation and translation
    pos_i_rotated = np.dot(rotation, pos_i) + translation

    # Target position for atom j
    # We need: rel_pos_rotated = pos_j'_R' - pos_i_rotated
    # => pos_j'_R' = pos_i_rotated + rel_pos_rotated
    pos_j_target = pos_i_rotated + rel_pos_rotated

    equivalent_bonds = []

    # For inter-atomic bonds (iatom != jatom), we need to find (iatom', jatom', R')
    # where atom i' is the image of atom i, and atom j' is the image of atom j
    # For self-bonds (iatom == jatom), we need to find (iatom', iatom', R')
    # where atom i' is the image of atom i

    # Find which atom position maps to pos_i_rotated
    for iatom_new in range(natoms):
        pos_i_new = pos_frac[iatom_new]
        # Check if atom i' is the image of atom i under this symmetry operation
        if np.allclose(pos_i_new, pos_i_rotated, atol=1e-6):
            # Found iatom' - now find jatom' and R'
            for jatom_new in range(natoms):
                pos_j_new = pos_frac[jatom_new]

                # R' is the integer part of (pos_j_target - pos_j_new)
                R_new = pos_j_target - pos_j_new
                R_new_rounded = np.round(R_new).astype(int)

                # Check if this is a valid mapping
                if np.allclose(R_new, R_new_rounded, atol=1e-6):
                    # Verify the mapping is consistent
                    pos_check = pos_j_new + R_new_rounded
                    if np.allclose(pos_j_target, pos_check, atol=1e-6):
                        # Verify the bond type is preserved:
                        # Self-bonds (i==j) should only map to self-bonds (i'==j')
                        # Inter-atomic bonds (i!=j) should only map to inter-atomic bonds (i'!=j')
                        if (iatom == jatom) == (iatom_new == jatom_new):
                            equivalent_bonds.append(
                                (iatom_new, jatom_new, tuple(R_new_rounded))
                            )

    return equivalent_bonds


def symmetrize_hr_full(
    HR: np.ndarray, Rlist: List[Tuple[int, int, int]], atoms, verbose: bool = False
) -> np.ndarray:
    """
    Symmetrize HR matrices by averaging over symmetry-equivalent (i,j,R) combinations.

    This is the FULL symmetrization that considers atom indices, not just R vectors.

    Parameters
    ----------
    HR : np.ndarray
        Real-space Hamiltonian matrices
        Shape: (nR, nao, nao) or (nspins, nR, nao, nao)
    Rlist : list of tuples
        List of R vectors corresponding to HR array
    atoms : ase.Atoms
        Crystal structure
    verbose : bool
        Print progress information

    Returns
    -------
    HR_sym : np.ndarray
        Symmetrized HR matrices (same shape as input)
    """
    if not HAS_SPGLIB:
        print("Warning: spglib not available. Returning HR as-is.")
        return HR

    # Get symmetry operations
    rotations, translations = get_symmetry_operations(atoms)

    if len(rotations) == 0:
        print("Warning: No symmetry operations found. Returning HR as-is.")
        return HR

    if verbose:
        print(f"Symmetrizing HR using {len(rotations)} symmetry operations...")
        print("Considering full (i,j,R) symmetrization with atom permutations...")

    # Check if HR has spin dimension
    if HR.ndim == 4:
        nspins, nR, nao, _ = HR.shape
        has_spin = True
    else:
        nR, nao, _ = HR.shape
        nspins = 1
        has_spin = False
        HR = HR[np.newaxis, ...]

    # Get orbital information
    # We need to know which orbitals belong to which atoms
    # For now, assume equal number of orbitals per atom
    norbs_per_atom = nao // len(atoms)

    # Get positions in fractional coordinates (for consistency with symmetry operations)
    cell = atoms.get_cell()
    positions = atoms.get_positions()
    pos_frac = np.dot(positions, np.linalg.inv(cell))

    if verbose:
        print(f"  Number of atoms: {len(atoms)}")
        print(f"  Number of orbitals per atom: {norbs_per_atom}")

    # Create lookup dictionary for R vectors
    R_to_idx = {R: i for i, R in enumerate(Rlist)}

    # Initialize symmetrized HR
    HR_sym = np.zeros_like(HR)

    # Track which (i,j,R) have been processed
    processed = set()

    # Build mapping of (iatom, jatom, R) to HR blocks
    # For each bond (iatom, jatom, R), we need to:
    # 1. Find all equivalent bonds
    # 2. Average the corresponding HR blocks
    # 3. Assign the averaged HR block to all equivalent bonds

    # Process each R vector
    for iR, R in enumerate(Rlist):
        # Process each pair of atoms
        for iatom in range(len(atoms)):
            for jatom in range(len(atoms)):
                bond_key = (iatom, jatom, R)

                if bond_key in processed:
                    continue

                # Find all equivalent bonds
                equivalent_bonds = set()
                equivalent_bonds.add(bond_key)

                for rotation, translation in zip(rotations, translations):
                    equiv = apply_symmetry_to_bond(
                        iatom, jatom, R, rotation, translation, atoms, pos_frac
                    )
                    for e in equiv:
                        equivalent_bonds.add(e)

                # Mark all as processed
                processed.update(equivalent_bonds)

                # Get indices of equivalent R vectors that exist in our Rlist
                valid_equiv = []
                for iatom_e, jatom_e, R_e in equivalent_bonds:
                    if R_e in R_to_idx:
                        iR_e = R_to_idx[R_e]
                        valid_equiv.append((iatom_e, jatom_e, iR_e))

                if len(valid_equiv) == 0:
                    if verbose and R == (0, 0, 0):
                        print(
                            f"  Warning: No equivalent bond found for ({iatom},{jatom},{R})"
                        )
                    # Use original HR block
                elif len(valid_equiv) == 1:
                    # No averaging needed
                    iatom_e, jatom_e, iR_e = valid_equiv[0]
                    orb_i = slice(iatom * norbs_per_atom, (iatom + 1) * norbs_per_atom)
                    orb_j = slice(jatom * norbs_per_atom, (jatom + 1) * norbs_per_atom)
                    orb_i_e = slice(
                        iatom_e * norbs_per_atom, (iatom_e + 1) * norbs_per_atom
                    )
                    orb_j_e = slice(
                        jatom_e * norbs_per_atom, (jatom_e + 1) * norbs_per_atom
                    )

                    for ispin in range(nspins):
                        HR_sym[ispin, iR, orb_i, orb_j] = HR[
                            ispin, iR_e, orb_i_e, orb_j_e
                        ]
                else:
                    # Average over equivalent bonds
                    if (
                        verbose
                        and len(equivalent_bonds) > 1
                        and R == (0, 0, 0)
                        and iatom == 0
                        and jatom == 1
                    ):
                        print(
                            f"  Bond ({iatom},{jatom},{R}): averaging over {len(valid_equiv)} equivalent bonds"
                        )

                    # Collect all HR blocks for equivalent bonds
                    HR_blocks = []
                    for iatom_e, jatom_e, iR_e in valid_equiv:
                        orb_i_e = slice(
                            iatom_e * norbs_per_atom, (iatom_e + 1) * norbs_per_atom
                        )
                        orb_j_e = slice(
                            jatom_e * norbs_per_atom, (jatom_e + 1) * norbs_per_atom
                        )
                        HR_blocks.append(HR[:, iR_e, orb_i_e, orb_j_e])

                    # Average
                    HR_avg = np.mean(HR_blocks, axis=0)

                    # Assign to all equivalent positions
                    for iatom_e, jatom_e, iR_e in valid_equiv:
                        orb_i = slice(
                            iatom * norbs_per_atom, (iatom + 1) * norbs_per_atom
                        )
                        orb_j = slice(
                            jatom * norbs_per_atom, (jatom + 1) * norbs_per_atom
                        )
                        orb_i_e = slice(
                            iatom_e * norbs_per_atom, (iatom_e + 1) * norbs_per_atom
                        )
                        orb_j_e = slice(
                            jatom_e * norbs_per_atom, (jatom_e + 1) * norbs_per_atom
                        )

                        for ispin in range(nspins):
                            HR_sym[ispin, iR_e, orb_i_e, orb_j_e] = HR_avg[ispin]

    if not has_spin:
        HR_sym = HR_sym[0]

    return HR_sym


def symmetrize_hr(
    HR: np.ndarray, Rlist: List[Tuple[int, int, int]], atoms, verbose: bool = False
) -> np.ndarray:
    """
    Symmetrize HR matrices (wrapper that calls full symmetrization).

    Parameters
    ----------
    HR : np.ndarray
        Real-space Hamiltonian matrices
    Rlist : list of tuples
        List of R vectors
    atoms : ase.Atoms
        Crystal structure
    verbose : bool
        Print progress

    Returns
    -------
    HR_sym : np.ndarray
        Symmetrized HR matrices
    """
    # Use full (i,j,R) symmetrization
    return symmetrize_hr_full(HR, Rlist, atoms, verbose=verbose)


if __name__ == "__main__":
    """Test the full symmetrization function"""
    import sys

    sys.path.insert(0, "HamiltonIO")

    from HamiltonIO.gpaw.gpaw_wrapper import GPAWParser

    print("=" * 80)
    print("Testing Full (i,j,R) HR Symmetrization")
    print("=" * 80)

    # Load data
    parser = GPAWParser(
        gpw_file="fe_bcc_no_symmetry.gpw",
        pickle_file="fe_bcc_no_symmetry_hamiltonian.pkl",
    )
    model_up, model_dn = parser.get_models()
    atoms = parser.atoms

    # Get HR and Rlist
    HR = model_up.HR
    Rlist = [tuple(R) for R in model_up.Rlist]

    print("\nOriginal HR:")
    print(f"  Shape: {HR.shape}")
    print(f"  Number of R vectors: {len(Rlist)}")

    # Check the bond (atom0 -> atom1, R=0,0,0) before symmetrization
    # Assuming 9 orbitals per atom (18 total / 2 atoms)
    norbs_per_atom = HR.shape[1] // len(atoms)
    print(f"  Orbitals per atom: {norbs_per_atom}")

    orb_0 = slice(0, norbs_per_atom)
    orb_1 = slice(norbs_per_atom, 2 * norbs_per_atom)

    print("\nBefore symmetrization (atom0->atom1 blocks):")
    test_R = [(0, 0, 0), (-1, 0, 0), (0, -1, 0), (-1, -1, 0)]
    for R in test_R:
        if R in Rlist:
            iR = Rlist.index(R)
            hr_block = HR[iR, orb_0, orb_1]
            trace = np.trace(hr_block.real)
            norm = np.linalg.norm(hr_block)
            print(f"  R={R}: trace={trace:10.3f}, norm={norm:11.3f}")

    # Symmetrize
    print("\nSymmetrizing HR with full (i,j,R) method...")
    HR_sym = symmetrize_hr_full(HR, Rlist, atoms, verbose=True)

    print("\nAfter symmetrization (atom0->atom1 blocks):")
    for R in test_R:
        if R in Rlist:
            iR = Rlist.index(R)
            hr_block = HR_sym[iR, orb_0, orb_1]
            trace = np.trace(hr_block.real)
            norm = np.linalg.norm(hr_block)
            print(f"  R={R}: trace={trace:10.3f}, norm={norm:11.3f}")

    print("\n" + "=" * 80)
    print("CONCLUSION:")
    print("=" * 80)
    print("Full (i,j,R) symmetrization complete!")
    print("HR matrix blocks for equivalent bonds should now have same trace/norm.")
