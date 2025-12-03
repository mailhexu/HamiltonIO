#!/usr/bin/env python3
"""
HamiltonIO ABACUS Command Line Interface

Provides command-line tools for analyzing ABACUS Hamiltonians.
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
        return "0.3.0"  # Fallback version from pyproject.toml


def analyze_intra_atomic(
    outpath=None,
    outpath_nosoc=None,
    outpath_soc=None,
    output_file=None,
    atoms=None,
    pauli_decomp=True,
    show_matrix=False,
    spin=None,
    binary=False,
):
    """
    Analyze intra-atomic (on-site) Hamiltonian for ABACUS calculations.

    Parameters:
        outpath: str or Path
            Path to ABACUS OUT.* directory (for single calculation)
        outpath_nosoc: str or Path
            Path to ABACUS OUT.* directory with SOC=0 (for split-SOC)
        outpath_soc: str or Path
            Path to ABACUS OUT.* directory with SOC=1 (for split-SOC)
        output_file: str or None
            Output file path. If None, print to stdout
        atoms: list of int or None
            List of atom indices to analyze. If None, analyze all atoms
        pauli_decomp: bool
            Whether to perform Pauli decomposition for spinor systems
        show_matrix: bool
            Whether to show full matrix elements
        spin: str or None
            Spin configuration: 'non-polarized', 'collinear', 'noncollinear', or None (auto-detect)
        binary: bool
            Whether to read binary CSR files

    Returns:
        AbacusWrapper: The loaded Hamiltonian object
    """
    from .abacus_wrapper import AbacusParser, AbacusSplitSOCParser
    from HamiltonIO.print_hamiltonian import print_intra_atomic_hamiltonian

    # Determine if split-SOC or single calculation
    if outpath_nosoc is not None and outpath_soc is not None:
        # Split-SOC mode
        outpath_nosoc = Path(outpath_nosoc)
        outpath_soc = Path(outpath_soc)

        if not outpath_nosoc.exists():
            raise FileNotFoundError(f"OUT directory not found: {outpath_nosoc}")
        if not outpath_soc.exists():
            raise FileNotFoundError(f"OUT directory not found: {outpath_soc}")

        parser = AbacusSplitSOCParser(
            outpath_nosoc=str(outpath_nosoc),
            outpath_soc=str(outpath_soc),
            binary=binary,
        )
        ham = parser.parse()

    elif outpath is not None:
        # Single calculation mode
        outpath = Path(outpath)
        if not outpath.exists():
            raise FileNotFoundError(f"OUT directory not found: {outpath}")

        parser = AbacusParser(outpath=str(outpath), spin=spin, binary=binary)
        result = parser.get_models()

        # Handle collinear case (returns tuple)
        if isinstance(result, tuple):
            ham = result
            print(
                "\nNote: Collinear calculation detected. Analyzing both spin components."
            )
        else:
            ham = result

    else:
        raise ValueError(
            "Must provide either 'outpath' or both 'outpath_nosoc' and 'outpath_soc'"
        )

    # Perform analysis
    if ham is not None:
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
    print("=== HamiltonIO ABACUS Intra-Atomic Analysis ===")
    print(f"Version: {get_version()}")
    print()

    # Validate input
    if args.outpath_nosoc and args.outpath_soc:
        print(f"✓ Split-SOC mode")
        print(f"  OUT.* (nosoc): {args.outpath_nosoc}")
        print(f"  OUT.* (soc):   {args.outpath_soc}")
        mode = "split-soc"
    elif args.outpath:
        print(f"✓ Single calculation mode")
        print(f"  OUT.* directory: {args.outpath}")
        mode = "single"
    else:
        print(
            "✗ Error: Must provide either --outpath or both --outpath-nosoc and --outpath-soc"
        )
        sys.exit(1)

    # Parse atom indices
    atoms = None
    if args.atoms:
        try:
            atoms = [int(x.strip()) for x in args.atoms.split(",")]
            print(f"\n  Analyzing atoms: {atoms}")
        except ValueError:
            print(f"✗ Error: Invalid atom indices: {args.atoms}")
            sys.exit(1)
    else:
        print(f"\n  Analyzing atoms: all")

    print(f"  Pauli decomposition: {'enabled' if not args.no_pauli else 'disabled'}")
    print(f"  Show matrices: {'yes' if args.show_matrix else 'no'}")
    print(f"  Binary format: {'yes' if args.binary else 'no'}")

    if args.output:
        print(f"  Output file: {args.output}")
    print()

    try:
        print("Loading ABACUS Hamiltonian...")

        if mode == "split-soc":
            ham = analyze_intra_atomic(
                outpath_nosoc=args.outpath_nosoc,
                outpath_soc=args.outpath_soc,
                output_file=args.output,
                atoms=atoms,
                pauli_decomp=not args.no_pauli,
                show_matrix=args.show_matrix,
                binary=args.binary,
            )
        else:
            ham = analyze_intra_atomic(
                outpath=args.outpath,
                output_file=args.output,
                atoms=atoms,
                pauli_decomp=not args.no_pauli,
                show_matrix=args.show_matrix,
                spin=args.spin,
                binary=args.binary,
            )

        # Summary
        print()
        print("=" * 80)
        print("✓ Analysis complete!")
        print("=" * 80)
        print()
        print(f"System: ABACUS")
        if ham is not None:
            # If tuple (collinear), use first component for summary
            h_summary = ham[0] if isinstance(ham, tuple) else ham

            if hasattr(h_summary, "atoms") and h_summary.atoms is not None:
                print(f"Number of atoms: {len(h_summary.atoms)}")
            if hasattr(h_summary, "nbasis"):
                print(f"Basis functions: {h_summary.nbasis}")
            if hasattr(h_summary, "split_soc"):
                print(f"Split-SOC: {h_summary.split_soc}")
        print()

        if args.output:
            import os

            file_size = os.path.getsize(args.output)
            print(f"Output written to: {Path(args.output).absolute()}")
            print(f"File size: {file_size / 1024:.1f} KB")

    except FileNotFoundError as e:
        print(f"✗ Error: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"✗ Error during analysis: {e}")
        import traceback

        traceback.print_exc()
        sys.exit(1)


def main():
    """Main entry point for hamiltonio-abacus CLI."""
    parser = argparse.ArgumentParser(
        prog="hamiltonio-abacus",
        description="HamiltonIO ABACUS analysis tools",
    )

    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {get_version()}"
    )

    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    # intra-atomic subcommand
    intra_parser = subparsers.add_parser(
        "intra-atomic",
        help="Extract and analyze the on-site (R=0) Hamiltonian blocks for each atom.",
    )

    # Input options
    input_group = intra_parser.add_argument_group("Input options")
    input_group.add_argument(
        "--outpath",
        type=str,
        help="Path to ABACUS OUT.* directory (for single calculation)",
    )
    input_group.add_argument(
        "--outpath-nosoc",
        type=str,
        help="Path to OUT.* directory with SOC=0 (for split-SOC analysis)",
    )
    input_group.add_argument(
        "--outpath-soc",
        type=str,
        help="Path to OUT.* directory with SOC=1 (for split-SOC analysis)",
    )

    # Output options
    output_group = intra_parser.add_argument_group("Output options")
    output_group.add_argument(
        "-o",
        "--output",
        type=str,
        help="Output file path (default: print to stdout)",
    )

    # Analysis options
    analysis_group = intra_parser.add_argument_group("Analysis options")
    analysis_group.add_argument(
        "-a",
        "--atoms",
        type=str,
        help="Comma-separated atom indices to analyze (e.g., '0,1,2'). Default: all atoms",
    )
    analysis_group.add_argument(
        "--no-pauli",
        action="store_true",
        help="Disable Pauli decomposition for spinor systems",
    )
    analysis_group.add_argument(
        "--show-matrix",
        action="store_true",
        help="Show full matrix elements (can be large for big systems)",
    )
    analysis_group.add_argument(
        "--spin",
        type=str,
        choices=["non-polarized", "collinear", "noncollinear"],
        help="Spin configuration (auto-detected if not specified)",
    )
    analysis_group.add_argument(
        "--binary",
        action="store_true",
        help="Read binary CSR format files",
    )

    intra_parser.set_defaults(func=intra_atomic_command)

    # Parse arguments
    args = parser.parse_args()

    if args.command is None:
        parser.print_help()
        sys.exit(0)

    # Execute command
    args.func(args)


if __name__ == "__main__":
    main()
