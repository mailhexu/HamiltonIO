#!/usr/bin/env python3
"""
HamiltonIO EPW Command Line Interface

Provides command-line tools for converting and analyzing EPW (Electron-Phonon Wannier) data.
"""

import argparse
import sys
import time
from pathlib import Path

from .epwparser import save_epmat_to_nc


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


def main():
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        prog="hamiltonio-epw",
        description="HamiltonIO EPW (Electron-Phonon Wannier) tools",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s epw_to_nc --path ./epw_data --prefix material --output epmat.nc
  %(prog)s epw_to_nc -p /path/to/files -n my_material -o output.nc
  %(prog)s epw_to_nc --path ./ --prefix test --force

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

    # Parse arguments
    args = parser.parse_args()

    # Handle commands
    if args.command == "epw_to_nc":
        convert_to_netcdf(args)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
