#!/usr/bin/env python3
"""
Simple example usage of the WignerData class for EPW calculations.

This script demonstrates basic file I/O operations with Wigner-Seitz data.
"""

from wigner import WignerData


def main():
    """Main example function."""
    print("=== WignerData Simple Example ===\n")

    # 1. Read data from file
    print("1. Reading Wigner-Seitz data from file...")
    data = WignerData.from_file("test/up/wigner.fmt")
    print("   ✅ Data loaded successfully!")
    print(data.summary())
    print()

    # 2. Write data to a new file
    print("2. Writing data to new file...")
    output_file = "example_output.fmt"
    data.to_file(output_file)
    print(f"   ✅ Data written to '{output_file}'")
    print()

    # 3. Verify round-trip by reading the written file
    print("3. Verifying round-trip...")
    data_copy = WignerData.from_file(output_file)
    print("   ✅ Round-trip successful!")
    print(data_copy.summary())
    print()

    print("=== Example completed successfully! ===")


if __name__ == "__main__":
    main()