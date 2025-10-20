def cli_convert_to_netcdf():
    """Simple CLI for converting EPW data to NetCDF format."""
    import argparse
    import os
    import sys

    parser = argparse.ArgumentParser(description="Convert EPW data to NetCDF format")
    parser.add_argument(
        "--convert-netcdf", action="store_true", help="Convert EPW data to NetCDF"
    )
    parser.add_argument(
        "-p",
        "--path",
        default=".",
        help="Path to EPW files (default: current directory)",
    )
    parser.add_argument(
        "-n", "--prefix", default="epw", help="EPW file prefix (default: epw)"
    )
    parser.add_argument(
        "-o",
        "--output",
        default="epmat.nc",
        help="Output NetCDF file (default: epmat.nc)",
    )

    args = parser.parse_args()

    if not args.convert_netcdf:
        parser.print_help()
        return

    print("=== EPW to NetCDF Converter ===")

    # Check if directory exists
    if not os.path.exists(args.path):
        print(f"❌ Directory not found: {args.path}")
        sys.exit(1)

    # Check required files
    required_files = ["epwdata.fmt", "wigner.fmt", f"{args.prefix}.epmatwp"]
    missing = []

    for file in required_files:
        filepath = os.path.join(args.path, file)
        if os.path.exists(filepath):
            size = os.path.getsize(filepath)
            print(f"✅ {file}: {size:,} bytes")
        else:
            print(f"❌ {file}: Not found")
            missing.append(file)

    if missing:
        print(f"\\nMissing files: {missing}")
        sys.exit(1)

    # Check if output exists
    output_path = os.path.join(args.path, args.output)
    if os.path.exists(output_path):
        print(f"⚠️  Output file exists: {args.output}")
        response = input("Overwrite? (y/N): ")
        if response.lower() != "y":
            sys.exit(0)

    # Perform conversion
    print(f"\\n🔄 Converting {args.prefix}.epmatwp to {args.output}...")

    try:
        from epwparser import save_epmat_to_nc

        save_epmat_to_nc(path=args.path, prefix=args.prefix, ncfile=args.output)

        if os.path.exists(output_path):
            size = os.path.getsize(output_path)
            print("✅ Conversion completed!")
            print(f"   Output: {args.output} ({size / 1024**2:.1f} MB)")
        else:
            print("❌ Conversion failed - output file not created")
            sys.exit(1)

    except Exception as e:
        print(f"❌ Conversion failed: {e}")
        sys.exit(1)
