#!/usr/bin/env python3
"""
Example script demonstrating intra-atomic Hamiltonian analysis for SIESTA bccFe.

Purpose:
    - Load SIESTA Hamiltonian from bccFe calculations
    - Extract and print intra-atomic (on-site, R=(0,0,0)) Hamiltonian blocks per atom
    - Show SOC/non-SOC decomposition if available
    - Apply Pauli decomposition for spinor systems to understand spin components

How to run:
    cd examples/SIESTA
    uv run python print_intra_atomic_Fe.py

Expected behavior:
    - Prints intra-atomic Hamiltonian decomposition for each Fe atom
    - For SOC calculations: shows separate non-SOC and SOC contributions
    - For spinor systems: shows Pauli decomposition (I, σx, σy, σz components)

Physical interpretation:
    - I (identity) component: charge/scalar part of the Hamiltonian
    - σx, σy, σz components: spin-dependent couplings in x, y, z directions
    - In SOC systems: non-SOC part is mostly diagonal (crystal field),
      SOC part couples spin states
"""

from pathlib import Path


def main():
    # Import here to avoid issues if sisl is not installed
    from HamiltonIO.siesta.sisl_wrapper import SislParser

    # Get script directory to build relative paths
    script_dir = Path(__file__).parent

    # Example 1: Non-polarized bccFe (no spin)
    print("=" * 80)
    print("Example 1: bccFe non-polarized calculation")
    print("=" * 80)

    fdf_path = script_dir / "bccFe/nonpolarized/siesta.fdf"
    if fdf_path.exists():
        parser = SislParser(fdf_fname=str(fdf_path), ispin=None)
        ham = parser.get_model()

        print(f"\nLoaded Hamiltonian: {ham._name}")
        print(f"Number of atoms: {len(ham.atoms)}")
        print(f"Number of basis functions: {ham.nbasis}")
        print(f"nspin: {ham.nspin}")

        # Print intra-atomic Hamiltonian analysis
        print("\n" + "=" * 80)
        print("Intra-Atomic Hamiltonian Analysis")
        print("=" * 80)
        ham.print_intra_atomic_hamiltonian(pauli_decomp=False, show_matrix=False)

        # Optionally save to file
        output_file = "bccFe_nonpolarized_intra_atomic.txt"
        print(f"\nSaving detailed output to: {output_file}")
        ham.print_intra_atomic_hamiltonian(
            output_file=output_file, pauli_decomp=False, show_matrix=True
        )
    else:
        print(f"Data file not found: {fdf_path}")
        print("Please ensure SIESTA example data exists in examples/SIESTA/bccFe/")

    # Example 2: bccFe with SOC (spinor calculation)
    print("\n\n" + "=" * 80)
    print("Example 2: bccFe with SOC (spinor calculation)")
    print("=" * 80)

    fdf_path_soc = script_dir / "bccFe/SOC/siesta.fdf"
    if fdf_path_soc.exists():
        try:
            parser_soc = SislParser(fdf_fname=str(fdf_path_soc), ispin="merge")
            ham_soc = parser_soc.get_model()

            print(f"\nLoaded Hamiltonian: {ham_soc._name}")
            print(f"Number of atoms: {len(ham_soc.atoms)}")
            print(f"Number of basis functions: {ham_soc.nbasis}")
            print(f"nspin: {ham_soc.nspin}")

            # Print intra-atomic Hamiltonian with Pauli decomposition
            print("\n" + "=" * 80)
            print("Intra-Atomic Hamiltonian Analysis with Pauli Decomposition")
            print("=" * 80)
            print("\nPhysical interpretation:")
            print("  I (identity): Charge/scalar part (diagonal energy levels)")
            print("  σx: Spin-x coupling (off-diagonal in spin space)")
            print("  σy: Spin-y coupling (off-diagonal in spin space, imaginary)")
            print("  σz: Spin-z coupling (diagonal in spin space, spin splitting)")
            print("=" * 80)

            ham_soc.print_intra_atomic_hamiltonian(pauli_decomp=True, show_matrix=False)

            # Save to file with full matrix
            output_file_soc = "bccFe_SOC_intra_atomic.txt"
            print(f"\nSaving detailed output with matrices to: {output_file_soc}")
            ham_soc.print_intra_atomic_hamiltonian(
                output_file=output_file_soc, pauli_decomp=True, show_matrix=True
            )

        except Exception as e:
            print(f"Error processing SOC data: {e}")
            import traceback

            traceback.print_exc()
            print("This may require specific SIESTA/sisl configuration for SOC")
    else:
        print(f"SOC data file not found: {fdf_path_soc}")

    # Example 3: Filter specific atoms
    print("\n\n" + "=" * 80)
    print("Example 3: Analyzing specific atoms only")
    print("=" * 80)

    if fdf_path.exists():
        parser = SislParser(fdf_fname=str(fdf_path), ispin=None)
        ham = parser.get_model()

        print("\nAnalyzing only atom 0:")
        ham.print_intra_atomic_hamiltonian(
            atom_indices=[0], pauli_decomp=False, show_matrix=True
        )

    print("\n" + "=" * 80)
    print("Analysis complete!")
    print("=" * 80)
    print("\nKey insights from intra-atomic analysis:")
    print("  1. On-site energies show crystal field splitting")
    print(
        "  2. For SOC systems: compare non-SOC (orbital) vs SOC (spin-orbit) contributions"
    )
    print("  3. Pauli decomposition reveals spin texture and magnetic anisotropy")
    print(
        "  4. Use this for debugging DFT calculations and understanding local physics"
    )


if __name__ == "__main__":
    main()
