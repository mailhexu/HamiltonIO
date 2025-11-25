#!/usr/bin/env python3
"""
HamiltonIO SIESTA Command Line Interface

Provides command-line tools for analyzing SIESTA Hamiltonians.
"""

import argparse
import sys
from pathlib import Path


def get_version():
    """Get the HamiltonIO version."""
    try:
        from .. import __version__

        return __version__
    except ImportError:
        return "0.2.5"  # Fallback version from pyproject.toml


def analyze_intra_atomic(
    fdf_file,
    output_file=None,
    atoms=None,
    pauli_decomp=True,
    show_matrix=False,
    read_soc=False,
    ispin=None,
):
    """
    Analyze intra-atomic (on-site) Hamiltonian for SIESTA calculations.

    Parameters:
        fdf_file: str or Path
            Path to SIESTA .fdf file
        output_file: str or None
            Output file path. If None, print to stdout
        atoms: list of int or None
            List of atom indices to analyze. If None, analyze all atoms
        pauli_decomp: bool
            Whether to perform Pauli decomposition for spinor systems
        show_matrix: bool
            Whether to show full matrix elements
        read_soc: bool
            Whether to read split-SOC data (H = H_nosoc + H_soc)
        ispin: str or None
            Spin channel: None, 'up', 'down', or 'merge'

    Returns:
        SiestaHamiltonian: The loaded Hamiltonian object
    """
    from .sisl_wrapper import SislParser
    from HamiltonIO.print_hamiltonian import print_intra_atomic_hamiltonian

    # Validate input
    fdf_path = Path(fdf_file)
    if not fdf_path.exists():
        raise FileNotFoundError(f"FDF file not found: {fdf_path}")

    # Load Hamiltonian
    parser = SislParser(fdf_fname=str(fdf_path), ispin=ispin, read_H_soc=read_soc)
    ham = parser.get_model()

    # Perform analysis
    print_intra_atomic_hamiltonian(
        ham,
        atom_indices=atoms,
        output_file=output_file,
        pauli_decomp=pauli_decomp,
        show_matrix=show_matrix,
    )

    return ham


def intra_atomic_command(args):
    """Handle the intra-atomic analysis command."""
    print("=== HamiltonIO SIESTA Intra-Atomic Analysis ===")
    print(f"Version: {get_version()}")
    print()

    # Validate input file
    fdf_path = Path(args.fdf)
    if not fdf_path.exists():
        print(f"✗ Error: FDF file not found: {fdf_path}")
        sys.exit(1)

    print(f"✓ Input FDF file: {fdf_path.absolute()}")

    # Check for NetCDF file
    nc_file = fdf_path.parent / (fdf_path.stem + ".nc")
    if not nc_file.exists():
        print(f"✗ Error: NetCDF file not found: {nc_file}")
        print(
            "  Please ensure the SIESTA calculation has completed and generated the .nc file."
        )
        sys.exit(1)

    print(f"✓ NetCDF file: {nc_file}")
    print()

    # Parse atom indices
    atoms = None
    if args.atoms:
        try:
            atoms = [int(x) for x in args.atoms.split(",")]
            print(f"  Analyzing atoms: {atoms}")
        except ValueError:
            print(f"✗ Error: Invalid atom indices: {args.atoms}")
            print("  Please provide comma-separated integers (e.g., '0,1,2')")
            sys.exit(1)
    else:
        print("  Analyzing: all atoms")

    # Other options
    print(f"  Split-SOC: {'enabled' if args.split_soc else 'disabled'}")
    print(f"  Pauli decomposition: {'enabled' if args.pauli else 'disabled'}")
    print(f"  Show matrices: {'yes' if args.show_matrix else 'no'}")
    if args.ispin:
        print(f"  Spin channel: {args.ispin}")
    if args.output:
        print(f"  Output file: {args.output}")
    print()

    # Run analysis
    try:
        print("Loading SIESTA Hamiltonian...")
        ham = analyze_intra_atomic(
            fdf_file=args.fdf,
            output_file=args.output,
            atoms=atoms,
            pauli_decomp=args.pauli,
            show_matrix=args.show_matrix,
            read_soc=args.split_soc,
            ispin=args.ispin,
        )

        print()
        print("=" * 80)
        print("✓ Analysis complete!")
        print("=" * 80)
        print()
        print(f"System: {ham._name}")
        print(f"Number of atoms: {len(ham.atoms) if ham.atoms else 'N/A'}")
        print(f"Basis functions: {ham.nbasis}")
        print(f"Split-SOC: {ham.split_soc}")
        print()

        if args.output:
            output_path = Path(args.output)
            if output_path.exists():
                size = output_path.stat().st_size
                if size < 1024:
                    size_str = f"{size} B"
                elif size < 1024**2:
                    size_str = f"{size / 1024:.1f} KB"
                else:
                    size_str = f"{size / (1024**2):.1f} MB"
                print(f"Output written to: {output_path.absolute()}")
                print(f"File size: {size_str}")
            else:
                print(f"⚠ Warning: Output file not found: {output_path}")

    except FileNotFoundError as e:
        print(f"✗ Error: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"✗ Error during analysis: {e}")
        import traceback

        traceback.print_exc()
        sys.exit(1)


def main():
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        prog="hamiltonio-siesta",
        description="HamiltonIO SIESTA Hamiltonian analysis tools",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Analyze intra-atomic Hamiltonian for all atoms
  %(prog)s intra-atomic siesta.fdf

  # Analyze specific atoms with split-SOC
  %(prog)s intra-atomic siesta.fdf --atoms 0,1,2 --split-soc

  # Save to file with full matrices
  %(prog)s intra-atomic siesta.fdf -o analysis.txt --show-matrix

  # Analyze with spin-orbit coupling (SOC/non-collinear)
  %(prog)s intra-atomic siesta.fdf --ispin merge --pauli

For more information, visit: https://github.com/mailhexu/HamiltonIO
        """,
    )

    parser.add_argument(
        "--version",
        action="version",
        version=f"HamiltonIO SIESTA CLI {get_version()}",
    )

    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    # Intra-atomic analysis command
    intra_parser = subparsers.add_parser(
        "intra-atomic",
        help="Analyze intra-atomic (on-site) Hamiltonian",
        description="Extract and analyze the on-site (R=0) Hamiltonian blocks for each atom.",
    )
    intra_parser.add_argument("fdf", help="Path to SIESTA .fdf file")
    intra_parser.add_argument(
        "-o", "--output", help="Output file path (default: print to stdout)"
    )
    intra_parser.add_argument(
        "-a",
        "--atoms",
        help="Comma-separated atom indices to analyze (e.g., '0,1,2'). Default: all atoms",
    )
    intra_parser.add_argument(
        "--no-pauli",
        dest="pauli",
        action="store_false",
        default=True,
        help="Disable Pauli decomposition for spinor systems",
    )
    intra_parser.add_argument(
        "--show-matrix",
        action="store_true",
        help="Show full matrix elements (can be large for big systems)",
    )
    intra_parser.add_argument(
        "--split-soc",
        action="store_true",
        help="Read split-SOC data (H = H_nosoc + H_soc). "
        "Requires SIESTA calculation with SOC.split.SR.SO = True",
    )
    intra_parser.add_argument(
        "--ispin",
        choices=["up", "down", "merge"],
        help="Spin channel for collinear calculations. "
        "Use 'merge' for SOC/non-collinear as spinor system",
    )

    # Parse arguments
    args = parser.parse_args()

    # Handle commands
    if args.command == "intra-atomic":
        intra_atomic_command(args)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
