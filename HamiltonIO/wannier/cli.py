#!/usr/bin/env python3
"""
HamiltonIO Wannier90 Command Line Interface

Provides command-line tools for analyzing Wannier90 tight-binding models.
"""

import argparse
import os
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
    path=".",
    prefix="wannier90",
    structure_file=None,
    output_file=None,
    atoms=None,
    show_matrix=False,
):
    """
    Analyze intra-atomic (on-site) Hamiltonian for Wannier90 calculations.

    Parameters:
        path: str or Path
            Path to Wannier90 directory
        prefix: str
            Wannier90 file prefix
        structure_file: str or None
            Structure file name (e.g., POSCAR, structure.cif)
        output_file: str or None
            Output file path. If None, print to stdout
        atoms: list of int or None
            List of atom indices to analyze. If None, analyze all atoms
        show_matrix: bool
            Whether to show full matrix elements

    Returns:
        WannierHam: The loaded Hamiltonian object
    """
    from HamiltonIO.print_hamiltonian import print_intra_atomic_hamiltonian

    from .wannier_hamiltonian import WannierHam

    # Load Hamiltonian
    input_path = Path(path)

    # Resolve structure file path if provided
    if structure_file is not None:
        if os.path.isabs(structure_file):
            structure_path = structure_file
        else:
            structure_path = str(input_path / structure_file)
    else:
        structure_path = None

    ham = WannierHam.read_from_wannier_dir(
        path=str(input_path),
        prefix=prefix,
        structure_file=structure_path,
    )

    # Perform analysis
    print_intra_atomic_hamiltonian(
        ham,
        atom_indices=atoms,
        output_file=output_file,
        pauli_decomp=False,
        show_matrix=show_matrix,
    )

    return ham


def intra_atomic_command(args):
    """Handle the intra-atomic analysis command."""
    print("=== HamiltonIO Wannier90 Intra-Atomic Analysis ===")
    print(f"Version: {get_version()}")
    print()

    # Validate input directory
    input_path = Path(args.path)
    if not input_path.exists():
        print(f"Error: Directory not found: {input_path}")
        sys.exit(1)

    print(f"Input directory: {input_path.absolute()}")
    print(f"Wannier90 prefix: {args.prefix}")

    if args.structure_file is not None:
        print(f"Structure file: {args.structure_file}")
    print()

    # Check Wannier90 file
    hr_file = input_path / f"{args.prefix}_hr.dat"
    if not hr_file.exists():
        print(f"Error: Wannier90 HR file not found: {hr_file}")
        print(f"Please ensure {args.prefix}_hr.dat exists in the directory.")
        sys.exit(1)

    # Parse atom indices
    atoms = None
    if args.atoms:
        try:
            atoms = [int(x) for x in args.atoms.split(",")]
            print(f"Analyzing atoms: {atoms}")
        except ValueError:
            print(f"Error: Invalid atom indices: {args.atoms}")
            print("Please provide comma-separated integers (e.g., '0,1,2')")
            sys.exit(1)
    else:
        print("Analyzing: all atoms")

    print(f"Show matrices: {'yes' if args.show_matrix else 'no'}")
    if args.output:
        print(f"Output file: {args.output}")
    print()

    try:
        print("Loading Wannier90 Hamiltonian...")
        ham = analyze_intra_atomic(
            path=args.path,
            prefix=args.prefix,
            structure_file=args.structure_file,
            output_file=args.output,
            atoms=atoms,
            show_matrix=args.show_matrix,
        )

        print()
        print("=" * 80)
        print("Analysis complete!")
        print("=" * 80)
        print()
        print(f"System: {ham._name}")
        print(f"Basis functions: {ham.nbasis}")
        if ham.atoms is not None:
            print(f"Number of atoms: {len(ham.atoms)}")

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
                print(f"Warning: Output file not found: {output_path}")

    except FileNotFoundError as e:
        print(f"Error: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"Error during analysis: {e}")
        import traceback

        traceback.print_exc()
        sys.exit(1)


def distance_analysis(args):
    """Analyze Wannier hopping matrix elements vs. distance."""
    import matplotlib

    matplotlib.use("Agg")  # Non-interactive backend
    import matplotlib.pyplot as plt

    print("=== HamiltonIO Wannier90 Distance Analysis ===")
    print(f"Version: {get_version()}")
    print()

    # Validate input directory
    input_path = Path(args.path)
    if not input_path.exists():
        print(f" Error: Directory not found: {input_path}")
        sys.exit(1)

    print(f" Input directory: {input_path.absolute()}")
    print(f" Wannier90 prefix: {args.prefix}")
    print(f" Structure file: {args.structure_file or 'auto-detect'}")
    print(f" Output file: {args.output}")
    print()

    # Check Wannier90 file
    hr_file = input_path / f"{args.prefix}_hr.dat"
    if not hr_file.exists():
        print(f" Error: Wannier90 HR file not found: {hr_file}")
        print(f" Please ensure {args.prefix}_hr.dat exists in the directory.")
        sys.exit(1)

    # Load and plot
    print(" Loading Wannier90 Hamiltonian...")
    try:
        from .plot import plot_wannier_distance

        fig, ax = plt.subplots(figsize=(6, 5))

        plot_wannier_distance(
            ax=ax,
            path=str(input_path),
            prefix=args.prefix,
            structure_file=args.structure_file,
            structure_format=args.structure_format,
            ylim=args.ylim,
        )

        ax.set_title(f"Wannier Hopping vs Distance: {args.prefix}")

        plt.tight_layout()
        plt.savefig(args.output, dpi=150)

        output_size = Path(args.output).stat().st_size
        print(" Plot created successfully!")
        print(f"   Output file: {Path(args.output).absolute()}")
        if output_size < 1024:
            print(f"   File size: {output_size} B")
        elif output_size < 1024**2:
            print(f"   File size: {output_size / 1024:.1f} KB")
        else:
            print(f"   File size: {output_size / (1024**2):.1f} MB")

    except FileNotFoundError as e:
        print(f" Error: {e}")
        sys.exit(1)
    except Exception as e:
        print(f" Error creating plot: {e}")
        import traceback

        traceback.print_exc()
        sys.exit(1)


def main():
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        prog="hamiltonio-wannier",
        description="HamiltonIO Wannier90 tight-binding model tools",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Analyze intra-atomic Hamiltonian
  %(prog)s intra-atomic --path ./wannier_calc --prefix wannier90

  # Analyze with external structure file (when .win lacks cell/positions)
  %(prog)s intra-atomic -p ./ -n wannier90 --structure-file POSCAR

  # Analyze specific atoms with full matrices
  %(prog)s intra-atomic -p ./ -n wannier90 --atoms 0,1,2 --show-matrix -o onsite.txt

  # Analyze hopping vs distance (auto-detect structure file)
  %(prog)s distance --path ./wannier_calc --prefix wannier90

  # Specify structure file and format
  %(prog)s distance -p ./ -n wannier90 -s POSCAR -f vasp

  # Custom output and y-axis limits
  %(prog)s distance -p ./ -n material -o hopping.pdf --ylim 1e-4 10

For more information, visit: https://github.com/mailhexu/HamiltonIO
        """,
    )

    parser.add_argument(
        "--version",
        action="version",
        version=f"HamiltonIO Wannier90 CLI {get_version()}",
    )

    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    # Intra-atomic analysis command
    intra_parser = subparsers.add_parser(
        "intra-atomic",
        help="Analyze intra-atomic (on-site) Hamiltonian",
        description="Extract and analyze the on-site (R=0) Hamiltonian blocks for each atom.",
    )
    intra_parser.add_argument(
        "-p",
        "--path",
        default=".",
        help="Path to Wannier90 directory (default: current directory)",
    )
    intra_parser.add_argument(
        "-n",
        "--prefix",
        default="wannier90",
        help="Wannier90 file prefix (default: wannier90)",
    )
    intra_parser.add_argument(
        "-s",
        "--structure-file",
        help="Structure file name (e.g., POSCAR, structure.cif). "
        "Use if .win file lacks cell and atomic positions.",
    )
    intra_parser.add_argument(
        "-a",
        "--atoms",
        help="Comma-separated atom indices to analyze (e.g., '0,1,2'). Default: all atoms",
    )
    intra_parser.add_argument(
        "-o",
        "--output",
        help="Output file path (default: print to stdout)",
    )
    intra_parser.add_argument(
        "--show-matrix",
        action="store_true",
        help="Show full matrix elements (can be large for big systems)",
    )

    # Distance analysis command
    distance_parser = subparsers.add_parser(
        "distance",
        help="Analyze Wannier hopping matrix elements vs. distance",
    )
    distance_parser.add_argument(
        "-p",
        "--path",
        default=".",
        help="Path to Wannier90 directory (default: current directory)",
    )
    distance_parser.add_argument(
        "-n",
        "--prefix",
        default="wannier90",
        help="Wannier90 file prefix (default: wannier90)",
    )
    distance_parser.add_argument(
        "-s",
        "--structure-file",
        help="Structure file name (e.g., POSCAR, structure.cif). "
        "Auto-detected if not specified.",
    )
    distance_parser.add_argument(
        "-f",
        "--structure-format",
        help="Structure file format (e.g., vasp, cif, espresso-in). "
        "Auto-detected if not specified.",
    )
    distance_parser.add_argument(
        "-o",
        "--output",
        default="wannier_distance.pdf",
        help="Output plot filename (default: wannier_distance.pdf)",
    )
    distance_parser.add_argument(
        "--ylim",
        nargs=2,
        type=float,
        default=[1e-5, 10],
        metavar=("MIN", "MAX"),
        help="Y-axis limits (default: 1e-5 10)",
    )

    # Parse arguments
    args = parser.parse_args()

    # Handle commands
    if args.command == "intra-atomic":
        intra_atomic_command(args)
    elif args.command == "distance":
        distance_analysis(args)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
