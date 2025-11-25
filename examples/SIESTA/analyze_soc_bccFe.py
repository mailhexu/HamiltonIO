#!/usr/bin/env python3
"""
Comprehensive analysis of BCC Fe with SOC using intra-atomic Hamiltonian methods.

Purpose:
    - Demonstrate detailed analysis of spin-orbit coupling effects in bccFe
    - Show how to extract and interpret Pauli decomposition
    - Compare SOC vs non-SOC calculations
    - Provide physical insights into magnetic properties

How to run:
    cd examples/SIESTA
    uv run python analyze_soc_bccFe.py

Expected output:
    - Detailed analysis of intra-atomic Hamiltonian
    - Pauli decomposition showing SOC contributions
    - Physical interpretation of spin components
    - Comparison with non-SOC calculation
"""

from pathlib import Path
import numpy as np


def analyze_soc_system():
    """Analyze BCC Fe with SOC."""
    from HamiltonIO.siesta.sisl_wrapper import SislParser
    from HamiltonIO.mathutils.pauli import pauli_block_all

    print("=" * 80)
    print("BCC Fe with SOC - Comprehensive Analysis")
    print("=" * 80)

    # Load SOC data
    script_dir = Path(__file__).parent
    fdf_path = script_dir / "bccFe/SOC/siesta.fdf"

    if not fdf_path.exists():
        print(f"Error: SOC data not found at {fdf_path}")
        return

    parser = SislParser(fdf_fname=str(fdf_path), ispin="merge")
    ham = parser.get_model()

    # System overview
    print(f"\n{'=' * 80}")
    print("System Overview")
    print("=" * 80)
    print(f"Material: BCC Iron (Fe)")
    print(f"Number of atoms: {len(ham.atoms)}")
    print(f"Total basis functions: {ham.nbasis}")
    print(f"Spin treatment: {'Spinor (SOC)' if ham.nbasis > 10 else 'Collinear'}")
    print(f"Number of R-vectors: {ham.nR}")

    # Orbital structure
    print(f"\n{'=' * 80}")
    print("Orbital Structure (Spinor Pairs)")
    print("=" * 80)

    if hasattr(ham, "orbs") and ham.orbs:
        # Group orbitals by type
        s_orbs = [orb for orb in ham.orbs if "|s" in orb or "|S" in orb]
        p_orbs = [orb for orb in ham.orbs if "|p" in orb or "|P" in orb]
        d_orbs = [orb for orb in ham.orbs if "|d" in orb or "|D" in orb]

        print(f"\nTotal orbitals: {len(ham.orbs)}")
        print(f"  s orbitals: {len(s_orbs)}")
        print(f"  p orbitals: {len(p_orbs)}")
        print(f"  d orbitals: {len(d_orbs)}")

        print("\ns orbitals (spin pairs):")
        for orb in s_orbs:
            idx = ham.orbs.index(orb)
            print(f"  [{idx:2d}] {orb}")

        print("\np orbitals (spin pairs):")
        for orb in p_orbs:
            idx = ham.orbs.index(orb)
            print(f"  [{idx:2d}] {orb}")

        print("\nd orbitals (spin pairs):")
        for orb in d_orbs:
            idx = ham.orbs.index(orb)
            print(f"  [{idx:2d}] {orb}")

    # Intra-atomic analysis
    print(f"\n{'=' * 80}")
    print("Intra-Atomic Hamiltonian (R=(0,0,0))")
    print("=" * 80)

    blocks = ham.get_intra_atomic_blocks()

    for iatom, data in blocks.items():
        H_full = data["H_full"]

        print(f"\nAtom {iatom} ({ham.atoms[iatom].symbol}):")
        print(f"  Matrix shape: {H_full.shape}")
        print(f"  Matrix type: {'Complex' if np.iscomplexobj(H_full) else 'Real'}")

        # Overall statistics
        print(f"\n  Overall Hamiltonian Statistics:")
        print(f"    Trace: {np.trace(H_full):.6f}")
        print(f"    Frobenius norm: {np.linalg.norm(H_full):.6f}")
        print(f"    Max |element|: {np.max(np.abs(H_full)):.6f}")
        print(
            f"    Average on-site energy: {np.trace(H_full) / H_full.shape[0]:.6f} eV"
        )

        # Check for imaginary parts (signature of SOC)
        imag_norm = np.linalg.norm(np.imag(H_full))
        print(f"    Imaginary part norm: {imag_norm:.6f}")
        print(f"    → {'SOC present' if imag_norm > 1e-6 else 'No significant SOC'}")

        # Pauli decomposition
        if H_full.shape[0] % 2 == 0:
            print(f"\n  Pauli Decomposition (H = I⊗σ₀ + Mx⊗σₓ + My⊗σᵧ + Mz⊗σᵤ):")
            MI, Mx, My, Mz = pauli_block_all(H_full)

            # I component (charge)
            print(f"\n    I (Identity/Charge) Component:")
            print(f"      Shape: {MI.shape}")
            print(f"      Trace: {np.trace(MI):.6f} eV")
            print(f"      Norm: {np.linalg.norm(MI):.6f} eV")
            print(f"      Max |element|: {np.max(np.abs(MI)):.6f} eV")
            print(f"      → Spin-independent orbital energies")

            # σₓ component
            print(f"\n    σₓ (Spin-X) Component:")
            print(f"      Trace: {np.trace(Mx):.6f} eV (should be ≈0)")
            print(f"      Norm: {np.linalg.norm(Mx):.6f} eV")
            print(f"      Max |element|: {np.max(np.abs(Mx)):.6f} eV")
            print(f"      → Spin-flip coupling in x-direction")

            # σᵧ component
            print(f"\n    σᵧ (Spin-Y) Component:")
            print(f"      Trace: {np.trace(My):.6f} eV (should be ≈0)")
            print(f"      Norm: {np.linalg.norm(My):.6f} eV")
            print(f"      Max |element|: {np.max(np.abs(My)):.6f} eV")
            print(f"      → Spin-flip coupling in y-direction (complex)")

            # σᵤ component
            print(f"\n    σᵤ (Spin-Z) Component:")
            print(f"      Trace: {np.trace(Mz):.6f} eV")
            print(f"      Norm: {np.linalg.norm(Mz):.6f} eV")
            print(f"      Max |element|: {np.max(np.abs(Mz)):.6f} eV")
            print(f"      → Spin splitting along z-axis")

            # Spin analysis
            spin_norm = np.sqrt(
                np.linalg.norm(Mx) ** 2
                + np.linalg.norm(My) ** 2
                + np.linalg.norm(Mz) ** 2
            )
            spin_flip_norm = np.sqrt(np.linalg.norm(Mx) ** 2 + np.linalg.norm(My) ** 2)

            print(f"\n    Spin Analysis:")
            print(f"      Total spin norm: {spin_norm:.6f} eV")
            print(f"      Spin-flip norm (√(σₓ² + σᵧ²)): {spin_flip_norm:.6f} eV")
            print(
                f"      Spin/Charge ratio: {spin_norm / np.linalg.norm(MI):.4f} ({spin_norm / np.linalg.norm(MI) * 100:.2f}%)"
            )
            print(
                f"      Anisotropy (σᵤ/σₓᵧ): {np.linalg.norm(Mz) / spin_flip_norm:.2f}"
            )
            print(
                f"      → {'Uniaxial' if np.linalg.norm(Mz) > 2 * spin_flip_norm else 'Isotropic'} spin character"
            )

    # Print formatted output to file
    print(f"\n{'=' * 80}")
    print("Saving Detailed Output")
    print("=" * 80)

    output_file = "bccFe_SOC_detailed_analysis.txt"
    print(f"\nGenerating detailed output file: {output_file}")
    ham.print_intra_atomic_hamiltonian(
        output_file=output_file, pauli_decomp=True, show_matrix=True
    )
    print(f"✓ Saved to {output_file}")


def compare_soc_vs_nonsoc():
    """Compare SOC vs non-SOC calculations."""
    from HamiltonIO.siesta.sisl_wrapper import SislParser

    print(f"\n\n{'=' * 80}")
    print("Comparison: SOC vs Non-SOC")
    print("=" * 80)

    script_dir = Path(__file__).parent

    # Load both systems
    fdf_nosoc = script_dir / "bccFe/nonpolarized/siesta.fdf"
    fdf_soc = script_dir / "bccFe/SOC/siesta.fdf"

    if not fdf_nosoc.exists() or not fdf_soc.exists():
        print("Error: Data files not found")
        return

    parser_nosoc = SislParser(fdf_fname=str(fdf_nosoc), ispin=None)
    ham_nosoc = parser_nosoc.get_model()

    parser_soc = SislParser(fdf_fname=str(fdf_soc), ispin="merge")
    ham_soc = parser_soc.get_model()

    # Get intra-atomic blocks
    blocks_nosoc = ham_nosoc.get_intra_atomic_blocks()
    blocks_soc = ham_soc.get_intra_atomic_blocks()

    print("\n" + "-" * 80)
    print("System Comparison")
    print("-" * 80)
    print(f"{'Property':<30} {'Non-SOC':>20} {'SOC':>20}")
    print("-" * 80)
    print(f"{'Basis functions':<30} {ham_nosoc.nbasis:>20} {ham_soc.nbasis:>20}")
    print(f"{'Matrix type':<30} {'Real':>20} {'Complex':>20}")
    print(f"{'Number of R-vectors':<30} {ham_nosoc.nR:>20} {ham_soc.nR:>20}")

    # Compare atom 0
    H_nosoc = blocks_nosoc[0]["H_full"]
    H_soc = blocks_soc[0]["H_full"]

    print(f"\n{'Atom 0 Hamiltonian':<30} {'Non-SOC':>20} {'SOC':>20}")
    print("-" * 80)
    print(f"{'Shape':<30} {str(H_nosoc.shape):>20} {str(H_soc.shape):>20}")
    print(
        f"{'Trace (eV)':<30} {np.trace(H_nosoc):>20.2f} {np.real(np.trace(H_soc)):>20.2f}"
    )
    print(
        f"{'Norm (eV)':<30} {np.linalg.norm(H_nosoc):>20.2f} {np.linalg.norm(H_soc):>20.2f}"
    )
    print(
        f"{'Max |element| (eV)':<30} {np.max(np.abs(H_nosoc)):>20.2f} {np.max(np.abs(H_soc)):>20.2f}"
    )

    # Imaginary parts
    imag_norm = np.linalg.norm(np.imag(H_soc))
    print(f"{'Imaginary norm (eV)':<30} {0.0:>20.2f} {imag_norm:>20.2f}")

    print("\n" + "=" * 80)
    print("Key Observations:")
    print("=" * 80)
    print("1. SOC doubles the basis size (spin-up and spin-down coupled)")
    print("2. SOC introduces complex matrix elements (imaginary parts)")
    print("3. Total energy (trace) approximately doubles with SOC")
    print("4. Imaginary parts indicate spin-orbit coupling strength")
    print(f"   → SOC strength: ~{imag_norm:.2f} eV")


def main():
    """Main analysis function."""
    print("""
╔══════════════════════════════════════════════════════════════════════════════╗
║                   BCC Fe SOC Analysis - HamiltonIO                           ║
║                                                                              ║
║  This script demonstrates comprehensive analysis of spin-orbit coupling     ║
║  effects in BCC iron using intra-atomic Hamiltonian methods.                ║
║                                                                              ║
║  See SOC_ANALYSIS.md for detailed physical interpretation.                  ║
╚══════════════════════════════════════════════════════════════════════════════╝
    """)

    try:
        # Analyze SOC system
        analyze_soc_system()

        # Compare with non-SOC
        compare_soc_vs_nonsoc()

        print(f"\n{'=' * 80}")
        print("Analysis Complete!")
        print("=" * 80)
        print("\nGenerated files:")
        print("  - bccFe_SOC_detailed_analysis.txt (detailed output)")
        print("\nFor more information, see:")
        print("  - SOC_ANALYSIS.md (physical interpretation)")
        print("  - examples/SIESTA/print_intra_atomic_Fe.py (basic usage)")

    except Exception as e:
        print(f"\nError during analysis: {e}")
        import traceback

        traceback.print_exc()


if __name__ == "__main__":
    main()
