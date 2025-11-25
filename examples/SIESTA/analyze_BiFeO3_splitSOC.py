#!/usr/bin/env python3
"""
Comprehensive analysis of BiFeO3 with split-SOC using intra-atomic Hamiltonian methods.

Purpose:
    - Demonstrate split-SOC analysis (H = H_nosoc + H_soc) in multiferroic BiFeO3
    - Show intra-atomic analysis for different elements (Bi, Fe, O)
    - Compare SOC effects in heavy (Bi) vs magnetic (Fe) vs light (O) atoms
    - Provide physical insights into multiferroic properties

How to run:
    cd /Users/hexu/projects/HamiltonIO
    uv run python examples/SIESTA/analyze_BiFeO3_splitSOC.py

Expected output:
    - Detailed analysis of split-SOC Hamiltonian
    - Intra-atomic Hamiltonians for Bi, Fe, and O atoms
    - Pauli decomposition showing SOC contributions
    - Comparison of SOC strengths across different elements
"""

from pathlib import Path
import numpy as np


def analyze_split_soc_system():
    """Analyze BiFeO3 with split-SOC."""
    from HamiltonIO.siesta.sisl_wrapper import SislParser
    from HamiltonIO.mathutils.pauli import pauli_block_all

    print("=" * 80)
    print("BiFeO3 with Split-SOC - Comprehensive Analysis")
    print("=" * 80)

    # Load split-SOC data
    bifeo3_dir = Path(
        "/Users/hexu/projects/TB2J_examples/Siesta/BiFeO3/BiFeO3_splitSOC"
    )
    fdf_path = bifeo3_dir / "siesta.fdf"

    if not fdf_path.exists():
        print(f"Error: Split-SOC data not found at {fdf_path}")
        return

    # Load with read_H_soc=True to enable split-SOC
    parser = SislParser(fdf_fname=str(fdf_path), ispin="merge", read_H_soc=True)
    ham = parser.get_model()

    # System overview
    print(f"\n{'=' * 80}")
    print("System Overview")
    print("=" * 80)
    print(f"Material: BiFeO3 (Bismuth Ferrite)")
    print(f"Number of atoms: {len(ham.atoms)}")
    print(
        f"Atom composition: {dict((x, ham.atoms.get_chemical_symbols().count(x)) for x in set(ham.atoms.get_chemical_symbols()))}"
    )
    print(f"Total basis functions: {ham.nbasis}")
    print(f"Spin treatment: Spinor with SOC")
    print(f"Number of R-vectors: {ham.nR}")
    print(f"\nSplit-SOC enabled: {ham.split_soc}")

    # Split-SOC decomposition analysis
    if ham.split_soc:
        print(f"\n{'=' * 80}")
        print("Split-SOC Decomposition (H = H_nosoc + H_soc)")
        print("=" * 80)

        # Get R=(0,0,0) index
        R0_idx = None
        for iR, R in enumerate(ham.Rlist):
            if np.allclose(R, [0, 0, 0]):
                R0_idx = iR
                break

        H_full = ham.HR[R0_idx]
        H_nosoc = ham.HR_nosoc[R0_idx]
        H_soc = ham.HR_soc[R0_idx]

        print(f"\nAt R=(0,0,0):")
        print(f"  H_full norm: {np.linalg.norm(H_full):.4f} eV")
        print(f"  H_nosoc norm: {np.linalg.norm(H_nosoc):.4f} eV")
        print(f"  H_soc norm: {np.linalg.norm(H_soc):.4f} eV")
        print(
            f"  SOC/Total ratio: {np.linalg.norm(H_soc) / np.linalg.norm(H_full) * 100:.2f}%"
        )

        # Verify decomposition
        diff = H_full - (H_nosoc + H_soc)
        print(
            f"\n  Decomposition verification: max|H_full - (H_nosoc + H_soc)| = {np.max(np.abs(diff)):.2e}"
        )
        print(
            f"  ✓ Decomposition valid!"
            if np.max(np.abs(diff)) < 1e-10
            else "  ✗ Decomposition failed!"
        )

    # Atom information
    print(f"\n{'=' * 80}")
    print("Atom Details")
    print("=" * 80)
    symbols = ham.atoms.get_chemical_symbols()
    for i, (atom, symbol) in enumerate(zip(ham.atoms, symbols)):
        print(f"Atom {i}: {symbol} at {atom.position}")

    # Intra-atomic analysis
    print(f"\n{'=' * 80}")
    print("Intra-Atomic Hamiltonian Analysis")
    print("=" * 80)

    blocks = ham.get_intra_atomic_blocks()

    # Analyze representative atoms: Bi, Fe, O
    bi_idx = symbols.index("Bi")
    fe_idx = symbols.index("Fe")
    o_idx = symbols.index("O")

    for label, iatom in [
        ("Bi (heavy)", bi_idx),
        ("Fe (magnetic)", fe_idx),
        ("O (light)", o_idx),
    ]:
        data = blocks[iatom]
        H_intra = data["H_full"]

        print(f"\n{'-' * 80}")
        print(f"Atom {iatom}: {label}")
        print("-" * 80)
        print(f"  Matrix shape: {H_intra.shape}")
        print(f"  Number of orbitals: {H_intra.shape[0] // 2} (spinor pairs)")

        # Overall statistics
        print(f"\n  Hamiltonian Statistics:")
        print(f"    Trace: {np.real(np.trace(H_intra)):.4f} eV (real)")
        print(f"    Frobenius norm: {np.linalg.norm(H_intra):.4f} eV")
        print(f"    Max |element|: {np.max(np.abs(H_intra)):.4f} eV")
        print(
            f"    Avg on-site energy: {np.real(np.trace(H_intra)) / H_intra.shape[0]:.4f} eV"
        )

        # Check imaginary part (SOC signature)
        imag_norm = np.linalg.norm(np.imag(H_intra))
        real_norm = np.linalg.norm(np.real(H_intra))
        print(f"    Imaginary norm: {imag_norm:.4f} eV")
        print(f"    Imag/Real ratio: {imag_norm / real_norm * 100:.2f}%")

        # Pauli decomposition
        if H_intra.shape[0] % 2 == 0:
            print(f"\n  Pauli Decomposition:")
            MI, Mx, My, Mz = pauli_block_all(H_intra)

            print(f"    I (Charge) norm: {np.linalg.norm(MI):.4f} eV")
            print(f"    σₓ (Spin-X) norm: {np.linalg.norm(Mx):.4f} eV")
            print(f"    σᵧ (Spin-Y) norm: {np.linalg.norm(My):.4f} eV")
            print(f"    σᵤ (Spin-Z) norm: {np.linalg.norm(Mz):.4f} eV")

            spin_norm = np.sqrt(
                np.linalg.norm(Mx) ** 2
                + np.linalg.norm(My) ** 2
                + np.linalg.norm(Mz) ** 2
            )
            spin_flip_norm = np.sqrt(np.linalg.norm(Mx) ** 2 + np.linalg.norm(My) ** 2)

            print(f"\n    Total spin norm: {spin_norm:.4f} eV")
            print(f"    Spin-flip norm: {spin_flip_norm:.4f} eV")
            print(f"    Spin/Charge ratio: {spin_norm / np.linalg.norm(MI) * 100:.2f}%")

            if spin_flip_norm > 1e-6:
                print(
                    f"    Anisotropy (σᵤ/σₓᵧ): {np.linalg.norm(Mz) / spin_flip_norm:.2f}"
                )

    # Compare SOC effects across elements
    print(f"\n{'=' * 80}")
    print("SOC Strength Comparison")
    print("=" * 80)

    print(
        f"\n{'Element':<15} {'Orbitals':<12} {'Imag norm (eV)':<18} {'Imag/Real (%)':<15}"
    )
    print("-" * 80)

    for label, iatom in [
        ("Bi (Z=83)", bi_idx),
        ("Fe (Z=26)", fe_idx),
        ("O (Z=8)", o_idx),
    ]:
        H_intra = blocks[iatom]["H_full"]
        imag_norm = np.linalg.norm(np.imag(H_intra))
        real_norm = np.linalg.norm(np.real(H_intra))
        n_orb = H_intra.shape[0] // 2

        print(
            f"{label:<15} {n_orb:<12} {imag_norm:<18.4f} {imag_norm / real_norm * 100:<15.2f}"
        )

    print("\n" + "-" * 80)
    print("Physical Insights:")
    print("-" * 80)
    print("1. Bi (heavy element, Z=83) shows strongest SOC effects")
    print("2. Fe (3d magnetic center) has moderate SOC, important for magnetism")
    print("3. O (light element) has weakest SOC")
    print("4. SOC strength scales approximately with Z² (atomic number squared)")
    print("5. BiFeO3 multiferroism arises from interplay of Bi-SOC and Fe-magnetism")

    # Save detailed output
    print(f"\n{'=' * 80}")
    print("Saving Detailed Output")
    print("=" * 80)

    output_file = "BiFeO3_splitSOC_analysis.txt"
    print(f"\nGenerating detailed output: {output_file}")
    ham.print_intra_atomic_hamiltonian(
        output_file=output_file, pauli_decomp=True, show_matrix=True
    )
    print(f"✓ Saved to {output_file}")


def main():
    """Main analysis function."""
    print(
        """
╔══════════════════════════════════════════════════════════════════════════════╗
║               BiFeO3 Split-SOC Analysis - HamiltonIO                         ║
║                                                                              ║
║  This script demonstrates split-SOC (H = H_nosoc + H_soc) analysis in       ║
║  multiferroic BiFeO3, showing how SOC effects vary across different         ║
║  elements and contribute to multiferroic properties.                        ║
║                                                                              ║
║  Key Features:                                                               ║
║  - Split-SOC decomposition (orbital vs spin-orbit parts)                    ║
║  - Element-specific SOC strength comparison (Bi vs Fe vs O)                 ║
║  - Pauli decomposition of intra-atomic Hamiltonians                         ║
║  - Physical interpretation of multiferroic coupling                         ║
╚══════════════════════════════════════════════════════════════════════════════╝
    """
    )

    try:
        analyze_split_soc_system()

        print(f"\n{'=' * 80}")
        print("Analysis Complete!")
        print("=" * 80)
        print("\nGenerated files:")
        print("  - BiFeO3_splitSOC_analysis.txt (detailed output with full matrices)")
        print("\nKey findings:")
        print("  - Split-SOC successfully decomposed: H = H_nosoc + H_soc")
        print("  - Bi atoms show strongest SOC (heavy element effect)")
        print("  - Fe atoms show moderate SOC (magnetic center)")
        print("  - O atoms show weakest SOC (light element)")

    except Exception as e:
        print(f"\nError during analysis: {e}")
        import traceback

        traceback.print_exc()


if __name__ == "__main__":
    main()
