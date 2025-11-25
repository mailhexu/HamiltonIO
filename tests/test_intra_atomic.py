#!/usr/bin/env python3
"""
Unit tests for intra-atomic Hamiltonian decomposition functionality.

Purpose:
    Test the extraction and decomposition of intra-atomic (R=(0,0,0)) Hamiltonian blocks.

How to run:
    uv run python -m pytest tests/test_intra_atomic.py -v
"""

import numpy as np
import pytest


def create_mock_hamiltonian(natoms=2, norb_per_atom=3, nspin=1, with_soc=False):
    """
    Create a mock LCAOHamiltonian for testing.

    Parameters:
        natoms: number of atoms
        norb_per_atom: number of orbitals per atom
        nspin: 1 for collinear, 2 for spinor
        with_soc: if True, create SOC-split Hamiltonian
    """
    from HamiltonIO.lcao_hamiltonian import LCAOHamiltonian
    from ase.atoms import Atoms

    # Create mock atoms
    symbols = ["Fe"] * natoms
    positions = [[0, 0, i * 3.0] for i in range(natoms)]
    atoms = Atoms(symbols=symbols, positions=positions, cell=[10, 10, 10])

    nbasis = natoms * norb_per_atom * nspin
    nR = 3  # Just a few R vectors including (0,0,0)

    # Create Rlist with (0,0,0) at index 0
    Rlist = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]])

    # Create random HR and SR
    HR = np.random.rand(nR, nbasis, nbasis) + 1j * np.random.rand(nR, nbasis, nbasis)
    SR = np.random.rand(nR, nbasis, nbasis) + 1j * np.random.rand(nR, nbasis, nbasis)

    # Make Hermitian
    for iR in range(nR):
        HR[iR] = (HR[iR] + HR[iR].conj().T) / 2
        SR[iR] = (SR[iR] + SR[iR].conj().T) / 2

    # Create mock orbitals (ABACUS-style for simplicity)
    from dataclasses import dataclass

    @dataclass
    class MockOrb:
        iatom: int
        sym: str

    orbs = []
    for iatom in range(natoms):
        for iorb in range(norb_per_atom * nspin):
            orbs.append(MockOrb(iatom=iatom, sym=f"orb_{iorb}"))

    if with_soc:
        # Split into non-SOC and SOC parts
        HR_nosoc = HR * 0.8  # 80% non-SOC
        HR_soc = HR * 0.2  # 20% SOC
        ham = LCAOHamiltonian(
            HR=None,
            SR=SR,
            Rlist=Rlist,
            nbasis=nbasis,
            atoms=atoms,
            orbs=orbs,
            nspin=nspin,
            HR_nosoc=HR_nosoc,
            HR_soc=HR_soc,
        )
    else:
        ham = LCAOHamiltonian(
            HR=HR,
            SR=SR,
            Rlist=Rlist,
            nbasis=nbasis,
            atoms=atoms,
            orbs=orbs,
            nspin=nspin,
        )

    return ham


def test_orbital_atom_mapping_collinear():
    """Test orbital-to-atom mapping for collinear system."""
    ham = create_mock_hamiltonian(natoms=2, norb_per_atom=3, nspin=1)
    orb_atom_map = ham._get_orbital_atom_mapping()

    assert len(orb_atom_map) == 2  # 2 atoms
    assert len(orb_atom_map[0]) == 3  # 3 orbitals per atom
    assert len(orb_atom_map[1]) == 3
    assert orb_atom_map[0] == [0, 1, 2]
    assert orb_atom_map[1] == [3, 4, 5]


def test_orbital_atom_mapping_spinor():
    """Test orbital-to-atom mapping for spinor system."""
    ham = create_mock_hamiltonian(natoms=2, norb_per_atom=3, nspin=2)
    orb_atom_map = ham._get_orbital_atom_mapping()

    assert len(orb_atom_map) == 2  # 2 atoms
    assert len(orb_atom_map[0]) == 6  # 3 orbitals * 2 spins per atom
    assert len(orb_atom_map[1]) == 6


def test_get_intra_atomic_blocks_collinear():
    """Test extraction of intra-atomic blocks for collinear system."""
    ham = create_mock_hamiltonian(natoms=2, norb_per_atom=3, nspin=1)
    blocks = ham.get_intra_atomic_blocks()

    assert len(blocks) == 2  # 2 atoms
    assert 0 in blocks
    assert 1 in blocks

    # Check shape of block for atom 0
    assert blocks[0]["H_full"].shape == (3, 3)
    assert blocks[0]["H_nosoc"] is None  # No SOC split
    assert blocks[0]["H_soc"] is None
    assert blocks[0]["orbital_indices"] == [0, 1, 2]


def test_get_intra_atomic_blocks_with_soc():
    """Test extraction of intra-atomic blocks with SOC decomposition."""
    ham = create_mock_hamiltonian(natoms=2, norb_per_atom=3, nspin=2, with_soc=True)
    blocks = ham.get_intra_atomic_blocks()

    assert len(blocks) == 2

    # Check SOC decomposition
    assert blocks[0]["H_nosoc"] is not None
    assert blocks[0]["H_soc"] is not None

    # Verify sum: H_full = H_nosoc + H_soc
    H_reconstructed = blocks[0]["H_nosoc"] + blocks[0]["H_soc"]
    assert np.allclose(H_reconstructed, blocks[0]["H_full"], rtol=1e-10)


def test_get_intra_atomic_blocks_filtered():
    """Test extraction with atom filtering."""
    ham = create_mock_hamiltonian(natoms=3, norb_per_atom=2, nspin=1)
    blocks = ham.get_intra_atomic_blocks(atom_indices=[0, 2])

    assert len(blocks) == 2
    assert 0 in blocks
    assert 2 in blocks
    assert 1 not in blocks


def test_print_intra_atomic_hamiltonian_collinear():
    """Test printing functionality for collinear system."""
    import tempfile
    import os
    from HamiltonIO.print_hamiltonian import print_intra_atomic_hamiltonian

    ham = create_mock_hamiltonian(natoms=2, norb_per_atom=2, nspin=1)

    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".txt") as f:
        output_file = f.name

    try:
        # Print to file
        print_intra_atomic_hamiltonian(
            ham, output_file=output_file, pauli_decomp=False, show_matrix=False
        )

        # Check file was created and has content
        assert os.path.exists(output_file)
        with open(output_file, "r") as f:
            content = f.read()
            assert "Intra-Atomic Hamiltonian Analysis" in content
            assert "Atom 0" in content
            assert "Atom 1" in content
            assert "Fe" in content  # Element symbol
    finally:
        if os.path.exists(output_file):
            os.remove(output_file)


def test_print_intra_atomic_hamiltonian_spinor_with_pauli():
    """Test printing with Pauli decomposition for spinor system."""
    import tempfile
    import os
    from HamiltonIO.print_hamiltonian import print_intra_atomic_hamiltonian

    ham = create_mock_hamiltonian(natoms=1, norb_per_atom=2, nspin=2)

    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".txt") as f:
        output_file = f.name

    try:
        # Print with Pauli decomposition
        print_intra_atomic_hamiltonian(
            ham, output_file=output_file, pauli_decomp=True, show_matrix=False
        )

        # Check Pauli components are mentioned
        with open(output_file, "r") as f:
            content = f.read()
            assert "Pauli Decomposition" in content
            assert "I (charge) component" in content or "charge" in content.lower()
            assert "spin-x" in content.lower() or "σₓ" in content or "σx" in content
    finally:
        if os.path.exists(output_file):
            os.remove(output_file)


def test_print_intra_atomic_with_soc_split():
    """Test printing with SOC decomposition."""
    import tempfile
    import os
    from HamiltonIO.print_hamiltonian import print_intra_atomic_hamiltonian

    ham = create_mock_hamiltonian(natoms=1, norb_per_atom=2, nspin=2, with_soc=True)

    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".txt") as f:
        output_file = f.name

    try:
        print_intra_atomic_hamiltonian(
            ham, output_file=output_file, pauli_decomp=True, show_matrix=False
        )

        with open(output_file, "r") as f:
            content = f.read()
            assert "SOC Decomposition" in content
            assert "Non-SOC Part" in content
            assert "SOC Part" in content
            assert "Verification" in content  # Check for sum verification
    finally:
        if os.path.exists(output_file):
            os.remove(output_file)


if __name__ == "__main__":
    # Run tests
    pytest.main([__file__, "-v"])
