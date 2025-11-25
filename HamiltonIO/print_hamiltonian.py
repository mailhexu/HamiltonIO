#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Hamiltonian printing utilities for intra-atomic analysis.

Provides functions to print and analyze intra-atomic (on-site) Hamiltonian blocks
with support for Pauli decomposition and SOC splitting.
"""

import numpy as np
from HamiltonIO.mathutils.pauli import pauli_block_all


def print_intra_atomic_hamiltonian(
    hamiltonian,
    atom_indices=None,
    output_file=None,
    pauli_decomp=True,
    show_matrix=False,
):
    """
    Print intra-atomic (on-site, R=(0,0,0)) Hamiltonian decomposition per atom.

    Parameters:
        hamiltonian: LCAOHamiltonian object
            Hamiltonian object with intra-atomic block extraction capability
        atom_indices: list of int or None
            Atom indices to print. If None, print all atoms.
        output_file: str or None
            If specified, write output to this file instead of stdout.
        pauli_decomp: bool
            If True and nspin==2, apply Pauli decomposition (I, σx, σy, σz).
        show_matrix: bool
            If True, print full matrix elements. Otherwise show summary statistics.

    For spinor systems (nspin=2), decomposes into:
        - I component (charge/scalar part)
        - σx component (spin-x coupling)
        - σy component (spin-y coupling)
        - σz component (spin-z coupling)

    For SOC-split systems (split_soc=True), shows decomposition into:
        - Non-SOC (scalar) part
        - SOC (spin-dependent) part
    """
    import sys

    # Open output file if specified
    if output_file:
        f = open(output_file, "w")
    else:
        f = sys.stdout

    try:
        # Get intra-atomic blocks
        blocks = hamiltonian.get_intra_atomic_blocks(atom_indices=atom_indices)

        # Get atom symbols if available
        atom_symbols = None
        if hasattr(hamiltonian, "atoms") and hamiltonian.atoms is not None:
            atom_symbols = hamiltonian.atoms.get_chemical_symbols()

        # Print header
        f.write("=" * 80 + "\n")
        f.write("Intra-Atomic Hamiltonian Analysis (R=(0,0,0))\n")
        f.write("=" * 80 + "\n")
        f.write(f"System: {hamiltonian._name}\n")
        f.write(f"Total basis functions: {hamiltonian.nbasis}\n")
        f.write(
            f"Spin configuration: {'Spinor (nspin=2)' if hamiltonian.nspin == 2 else 'Collinear/Non-polarized (nspin=1)'}\n"
        )
        f.write(f"SOC splitting: {'Yes' if hamiltonian.split_soc else 'No'}\n")
        f.write("=" * 80 + "\n\n")

        # Print each atom
        for iatom in sorted(blocks.keys()):
            atom_data = blocks[iatom]
            symbol = atom_symbols[iatom] if atom_symbols else "Unknown"
            norb_atom = len(atom_data["orbital_indices"])

            f.write(f"\n{'=' * 80}\n")
            f.write(f"Atom {iatom} ({symbol})\n")
            f.write(f"{'=' * 80}\n")
            f.write(f"Orbital indices: {atom_data['orbital_indices']}\n")
            f.write(f"Number of basis functions: {norb_atom}\n")

            # Print orbital labels if available
            if hasattr(hamiltonian, "orbs") and hamiltonian.orbs:
                orbital_labels = [
                    hamiltonian.orbs[i]
                    for i in atom_data["orbital_indices"]
                    if i < len(hamiltonian.orbs)
                ]
                if orbital_labels:
                    f.write("Orbital labels:\n")
                    for idx, label in zip(atom_data["orbital_indices"], orbital_labels):
                        f.write(f"  [{idx}] {label}\n")

            f.write(f"{'-' * 80}\n")

            # Process full Hamiltonian
            print_hamiltonian_block(
                f,
                atom_data["H_full"],
                "Full Hamiltonian",
                pauli_decomp=pauli_decomp,
                show_matrix=show_matrix,
            )

            # Process SOC decomposition if available
            if hamiltonian.split_soc:
                f.write(f"\n{'-' * 80}\n")
                f.write("SOC Decomposition:\n")
                f.write(f"{'-' * 80}\n")

                print_hamiltonian_block(
                    f,
                    atom_data["H_nosoc"],
                    "Non-SOC Part",
                    pauli_decomp=pauli_decomp,
                    show_matrix=show_matrix,
                )

                f.write("\n")
                print_hamiltonian_block(
                    f,
                    atom_data["H_soc"],
                    "SOC Part",
                    pauli_decomp=pauli_decomp,
                    show_matrix=show_matrix,
                )

                # Verify sum
                if atom_data["H_nosoc"] is not None and atom_data["H_soc"] is not None:
                    H_reconstructed = atom_data["H_nosoc"] + atom_data["H_soc"]
                    diff_norm = np.linalg.norm(H_reconstructed - atom_data["H_full"])
                    f.write(
                        f"\nVerification: ||H_full - (H_nosoc + H_soc)|| = {diff_norm:.2e}\n"
                    )

        f.write(f"\n{'=' * 80}\n")
        f.write("End of Intra-Atomic Hamiltonian Analysis\n")
        f.write("=" * 80 + "\n")

    finally:
        if output_file:
            f.close()


def print_hamiltonian_block(f, H, title, pauli_decomp=True, show_matrix=False):
    """
    Print a Hamiltonian block with optional Pauli decomposition.

    Parameters:
        f: file object to write to
        H: numpy array, Hamiltonian matrix block
        title: str, title for this block
        pauli_decomp: bool, whether to apply Pauli decomposition
        show_matrix: bool, whether to show full matrix
    """
    f.write(f"\n{title}:\n")
    f.write(f"  Shape: {H.shape}\n")

    # Basic statistics
    f.write(f"  Trace: {np.trace(H):.6f}\n")
    f.write(f"  Frobenius norm: {np.linalg.norm(H):.6f}\n")
    f.write(f"  Max abs element: {np.max(np.abs(H)):.6f}\n")

    # Show matrix if requested - print FULL matrices
    if show_matrix:
        # (a) Real part of original matrix
        f.write("\n  (a) Real part of H:\n")
        for row in H.real:
            f.write("    " + " ".join(f"{val:10.4f}" for val in row) + "\n")

        # (b) Imaginary part of original matrix
        if np.any(np.abs(np.imag(H)) > 1e-10):
            f.write("\n  (b) Imaginary part of H:\n")
            for row in H.imag:
                f.write("    " + " ".join(f"{val:10.4f}" for val in row) + "\n")

    # Pauli decomposition for spinor systems (check matrix shape, not nspin)
    # Spinor matrices have even dimension (2x the number of orbitals)
    if pauli_decomp and H.shape[0] % 2 == 0 and H.shape[0] > 1:
        f.write("\n  Pauli Decomposition (H = I⊗σ₀ + Mx⊗σₓ + My⊗σᵧ + Mz⊗σᵤ):\n")
        f.write("  Note: The Pauli components are matrices of shape (norb, norb)\n")
        f.write("        where norb = nbasis/2 (half the original size)\n")
        MI, Mx, My, Mz = pauli_block_all(H)

        f.write("\n    Statistics:\n")
        f.write(
            f"      I (charge):   Trace={np.trace(MI):10.6f}, Norm={np.linalg.norm(MI):10.6f}, Max|elem|={np.max(np.abs(MI)):10.6f}\n"
        )
        f.write(
            f"      σₓ (spin-x):  Trace={np.trace(Mx):10.6f}, Norm={np.linalg.norm(Mx):10.6f}, Max|elem|={np.max(np.abs(Mx)):10.6f}\n"
        )
        f.write(
            f"      σᵧ (spin-y):  Trace={np.trace(My):10.6f}, Norm={np.linalg.norm(My):10.6f}, Max|elem|={np.max(np.abs(My)):10.6f}\n"
        )
        f.write(
            f"      σᵤ (spin-z):  Trace={np.trace(Mz):10.6f}, Norm={np.linalg.norm(Mz):10.6f}, Max|elem|={np.max(np.abs(Mz)):10.6f}\n"
        )

        # Spin magnitude
        spin_norm = np.sqrt(
            np.linalg.norm(Mx) ** 2 + np.linalg.norm(My) ** 2 + np.linalg.norm(Mz) ** 2
        )
        f.write(f"      Total spin norm: {spin_norm:.6f}\n")

        # Print full Pauli component matrices if show_matrix is True
        if show_matrix:
            # (c) I component
            f.write("\n  (c) I (Identity/Charge) component - Real part:\n")
            for row in MI.real:
                f.write("    " + " ".join(f"{val:10.4f}" for val in row) + "\n")
            if np.any(np.abs(np.imag(MI)) > 1e-10):
                f.write("\n  (c) I (Identity/Charge) component - Imaginary part:\n")
                for row in MI.imag:
                    f.write("    " + " ".join(f"{val:10.4f}" for val in row) + "\n")

            # (d) σₓ component
            f.write("\n  (d) σₓ (Spin-X) component - Real part:\n")
            for row in Mx.real:
                f.write("    " + " ".join(f"{val:10.4f}" for val in row) + "\n")
            if np.any(np.abs(np.imag(Mx)) > 1e-10):
                f.write("\n  (d) σₓ (Spin-X) component - Imaginary part:\n")
                for row in Mx.imag:
                    f.write("    " + " ".join(f"{val:10.4f}" for val in row) + "\n")

            # (e) σᵧ component
            f.write("\n  (e) σᵧ (Spin-Y) component - Real part:\n")
            for row in My.real:
                f.write("    " + " ".join(f"{val:10.4f}" for val in row) + "\n")
            if np.any(np.abs(np.imag(My)) > 1e-10):
                f.write("\n  (e) σᵧ (Spin-Y) component - Imaginary part:\n")
                for row in My.imag:
                    f.write("    " + " ".join(f"{val:10.4f}" for val in row) + "\n")

            # (f) σᵤ component
            f.write("\n  (f) σᵤ (Spin-Z) component - Real part:\n")
            for row in Mz.real:
                f.write("    " + " ".join(f"{val:10.4f}" for val in row) + "\n")
            if np.any(np.abs(np.imag(Mz)) > 1e-10):
                f.write("\n  (f) σᵤ (Spin-Z) component - Imaginary part:\n")
                for row in Mz.imag:
                    f.write("    " + " ".join(f"{val:10.4f}" for val in row) + "\n")
