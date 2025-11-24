#!/usr/bin/env python3
"""
HamiltonIO EPW Command Line Interface

Provides command-line tools for converting and analyzing EPW
(Electron-Phonon Wannier) data.
"""

import argparse
import sys
import time
from pathlib import Path

from .epwparser import save_epmat_to_nc, Epmat, read_crystal_fmt
from .plot import plot_epw_distance


def get_version():
    """Get the HamiltonIO version."""
    try:
        from .. import __version__

        return __version__
    except ImportError:
        return "0.2.5"  # Fallback version from pyproject.toml


def validate_epw_files(path, prefix):
    """Validate that all required EPW files exist and return their info."""
    required_files = ["epwdata.fmt", "wigner.fmt", f"{prefix}.epmatwp"]
    file_info = []
    missing_files = []

    for filename in required_files:
        filepath = Path(path) / filename
        if filepath.exists():
            size = filepath.stat().st_size
            file_info.append((filename, size))
        else:
            missing_files.append(filename)

    return file_info, missing_files


def format_file_size(size_bytes):
    """Format file size in human-readable format."""
    if size_bytes < 1024:
        return f"{size_bytes} B"
    elif size_bytes < 1024**2:
        return f"{size_bytes / 1024:.1f} KB"
    elif size_bytes < 1024**3:
        return f"{size_bytes / (1024**2):.1f} MB"
    else:
        return f"{size_bytes / (1024**3):.1f} GB"


def convert_to_netcdf(args):
    """Convert EPW binary files to NetCDF format."""
    print("=== HamiltonIO EPW to NetCDF Converter ===")
    print(f"Version: {get_version()}")
    print()

    # Validate input directory
    input_path = Path(args.path)
    if not input_path.exists():
        print(f" Error: Directory not found: {input_path}")
        sys.exit(1)

    if not input_path.is_dir():
        print(f" Error: Path is not a directory: {input_path}")
        sys.exit(1)

    print(f" Input directory: {input_path.absolute()}")
    print(f" EPW file prefix: {args.prefix}")
    print(f" Output file: {args.output}")
    if args.dry_run:
        print(" Mode: Dry run (no conversion will be performed)")
    print()

    # Validate required files
    print(" Checking required files...")
    file_info, missing_files = validate_epw_files(args.path, args.prefix)

    total_size = 0
    for filename, size in file_info:
        print(f" {filename}: {format_file_size(size)}")
        total_size += size

    print(f"\n Total input size: {format_file_size(total_size)}")

    if missing_files:
        print("\n Missing required files:")
        for filename in missing_files:
            print(f"   - {filename}")
        print("\nPlease ensure all required EPW output files are present.")
        sys.exit(1)

    # Check output file
    output_path = Path(args.path) / args.output
    if output_path.exists():
        print(f"\n  Warning: Output file exists: {output_path}")
        if not args.force:
            response = input("Overwrite? (y/N): ")
            if response.lower() != "y":
                print("Conversion cancelled.")
                sys.exit(0)
        else:
            print("Overwriting existing file (force mode).")

    # Dry run mode - just check files and exit
    if args.dry_run:
        print("\n Dry run completed - all required files are present and valid.")
        print("   Run without --dry-run to perform the actual conversion.")
        return

    # Perform conversion
    print(f"\n Converting {args.prefix}.epmatwp to {args.output}...")
    print("   This may take a while for large files...")

    start_time = time.time()
    try:
        save_epmat_to_nc(path=args.path, prefix=args.prefix, ncfile=args.output)
        conversion_time = time.time() - start_time

        if output_path.exists():
            output_size = output_path.stat().st_size
            print(" Conversion completed successfully!")
            print(f"   Output file: {output_path.absolute()}")
            print(f"   File size: {format_file_size(output_size)}")
            print(f"   Conversion time: {conversion_time:.2f} seconds")
        else:
            print(" Error: Conversion failed - output file not created")
            sys.exit(1)

    except KeyboardInterrupt:
        print("\n Conversion cancelled by user")
        sys.exit(1)
    except Exception as e:
        print("\n Error: Conversion failed")
        print(f"   Details: {str(e)}")
        sys.exit(1)


def distance_analysis(args):
    """Analyze electron-phonon coupling vs. distance."""
    import matplotlib

    matplotlib.use("Agg")  # Non-interactive backend
    import matplotlib.pyplot as plt
    from ase.units import Bohr

    print("=== HamiltonIO EPW Distance Analysis ===")
    print(f"Version: {get_version()}")
    print()

    # Validate input directory
    input_path = Path(args.path)
    if not input_path.exists():
        print(f" Error: Directory not found: {input_path}")
        sys.exit(1)

    print(f" Input directory: {input_path.absolute()}")
    print(f" NetCDF file: {args.epmat_ncfile}")
    print(f" Crystal file: {args.crystal_fmt_file}")
    print(f" Mode index: {args.imode}")
    print(f" Distance type: {args.distance_type}")
    print(f" Output file: {args.output}")
    print()

    # Check required files
    ncfile_path = input_path / args.epmat_ncfile
    crystal_path = input_path / args.crystal_fmt_file

    if not ncfile_path.exists():
        print(f" Error: NetCDF file not found: {ncfile_path}")
        print(" Run 'hamiltonio-epw epw_to_nc' first to convert EPW data.")
        sys.exit(1)

    if not crystal_path.exists():
        print(f" Error: Crystal file not found: {crystal_path}")
        sys.exit(1)

    # Load data
    print(" Loading EPW data...")
    try:
        crystal = read_crystal_fmt(str(crystal_path))
        ep = Epmat()
        ep.crystal = crystal
        ep.read(path=str(input_path), prefix="prefix", epmat_ncfile=args.epmat_ncfile)
    except Exception as e:
        print(f" Error loading data: {e}")
        sys.exit(1)

    # Get cell
    cell = crystal.at.reshape(3, 3) * crystal.alat * Bohr

    # Calculate distances
    print(f" Calculating {args.distance_type}-based distances for mode {args.imode}...")
    try:
        if args.distance_type == "Rk":
            entries = ep.distance_resolved_couplings_Rk(imode=args.imode, cell=cell)
        elif args.distance_type == "Rg":
            entries = ep.distance_resolved_couplings_Rg(imode=args.imode, cell=cell)
        else:
            print(
                f" Error: Invalid distance_type '{args.distance_type}'."
                " Use 'Rk' or 'Rg'."
            )
            sys.exit(1)
    except Exception as e:
        print(f" Error calculating distances: {e}")
        sys.exit(1)

    print(f" Found {len(entries)} non-zero coupling elements")

    # Create plot
    print(" Creating plot...")
    fig, ax = plt.subplots(figsize=(6, 5))

    try:
        plot_epw_distance(
            ax=ax,
            path=str(input_path),
            imode=args.imode,
            distance_type=args.distance_type,
            epmat_ncfile=args.epmat_ncfile,
            crystal_fmt_file=args.crystal_fmt_file,
            ylim=args.ylim,
        )

        # Add title
        distance_label = "WF-WF" if args.distance_type == "Rk" else "WF-atom"
        ax.set_title(f"Mode {args.imode}: {distance_label} distance")

        plt.tight_layout()
        plt.savefig(args.output, dpi=150)

        output_size = Path(args.output).stat().st_size
        print(" Plot created successfully!")
        print(f"   Output file: {Path(args.output).absolute()}")
        print(f"   File size: {format_file_size(output_size)}")

    except Exception as e:
        print(f" Error creating plot: {e}")
        sys.exit(1)


def main():
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        prog="hamiltonio-epw",
        description="HamiltonIO EPW (Electron-Phonon Wannier) tools",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Convert EPW data to NetCDF
  %(prog)s epw_to_nc --path ./epw_data --prefix material --output epmat.nc
  %(prog)s epw_to_nc -p /path/to/files -n my_material -o output.nc
  %(prog)s epw_to_nc --path ./ --prefix test --force

  # Analyze distance-resolved coupling
  %(prog)s distance --path ./epw_data --imode 0 --distance-type Rk
  %(prog)s distance -p ./ -m 5 -t Rg -o mode5_Rg.pdf
  %(prog)s distance -m 0 --ylim 1e-4 10

For more information, visit: https://github.com/mailhexu/HamiltonIO
        """,
    )

    parser.add_argument(
        "--version", action="version", version=f"HamiltonIO EPW CLI {get_version()}"
    )

    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    # Convert command
    convert_parser = subparsers.add_parser(
        "epw_to_nc", help="Convert EPW binary files to NetCDF format"
    )
    convert_parser.add_argument(
        "-p",
        "--path",
        default=".",
        help="Path to directory containing EPW files (default: current directory)",
    )
    convert_parser.add_argument(
        "-n", "--prefix", default="epw", help="EPW file prefix (default: epw)"
    )
    convert_parser.add_argument(
        "-o",
        "--output",
        default="epmat.nc",
        help="Output NetCDF filename (default: epmat.nc)",
    )
    convert_parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing output file without prompting",
    )
    convert_parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Check files and show conversion info without performing conversion",
    )

    # Distance analysis command
    distance_parser = subparsers.add_parser(
        "distance",
        help="Analyze electron-phonon coupling vs. distance",
    )
    distance_parser.add_argument(
        "-p",
        "--path",
        default=".",
        help="Path to directory containing EPW files (default: current directory)",
    )
    distance_parser.add_argument(
        "-m",
        "--imode",
        type=int,
        required=True,
        help="Phonon mode index to analyze",
    )
    distance_parser.add_argument(
        "-t",
        "--distance-type",
        choices=["Rk", "Rg"],
        default="Rk",
        help="Distance type: 'Rk' for WF-WF, 'Rg' for WF-atom (default: Rk)",
    )
    distance_parser.add_argument(
        "-o",
        "--output",
        default="epw_distance.pdf",
        help="Output plot filename (default: epw_distance.pdf)",
    )
    distance_parser.add_argument(
        "--epmat-ncfile",
        default="epmat.nc",
        help="NetCDF file containing EPW data (default: epmat.nc)",
    )
    distance_parser.add_argument(
        "--crystal-fmt-file",
        default="crystal.fmt",
        help="Crystal format file (default: crystal.fmt)",
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
    if args.command == "epw_to_nc":
        convert_to_netcdf(args)
    elif args.command == "distance":
        distance_analysis(args)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
