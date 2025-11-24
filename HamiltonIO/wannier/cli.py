#!/usr/bin/env python3
"""
HamiltonIO Wannier90 Command Line Interface

Provides command-line tools for analyzing Wannier90 tight-binding models.
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
    if args.command == "distance":
        distance_analysis(args)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
