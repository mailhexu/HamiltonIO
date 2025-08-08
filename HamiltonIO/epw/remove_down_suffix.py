#!/usr/bin/env python3
"""
Script to remove '.down' suffix from filenames in a directory.

This script will rename files by removing '.down' from their names.
For example:
- 'crystal.down.fmt' -> 'crystal.fmt'
- 'epwdata.down.fmt' -> 'epwdata.fmt'
- 'wigner.down.fmt' -> 'wigner.fmt'

Usage:
    python remove_down_suffix.py [directory] [options]

Options:
    -d, --directory DIR    Directory to process (default: current directory)
    -n, --dry-run         Show what would be renamed without actually doing it
    -v, --verbose         Show detailed output
    -r, --recursive       Process subdirectories recursively
    -h, --help           Show this help message
"""

import os
import sys
import argparse
import glob
from pathlib import Path


def remove_down_from_filename(filepath, dry_run=False, verbose=False):
    """
    Remove '.down' from a filename.
    
    Parameters
    ----------
    filepath : str or Path
        Path to the file to rename
    dry_run : bool, optional
        If True, only show what would be done without actually renaming
    verbose : bool, optional
        If True, show detailed output
        
    Returns
    -------
    bool
        True if file was renamed (or would be renamed in dry-run), False otherwise
    """
    filepath = Path(filepath)
    
    # Check if filename contains '.down'
    if '.down' not in filepath.name:
        if verbose:
            print(f"   Skip: {filepath.name} (no '.down' found)")
        return False
    
    # Create new filename by removing '.down'
    new_name = filepath.name.replace('.down', '')
    new_filepath = filepath.parent / new_name
    
    # Check if target file already exists
    if new_filepath.exists():
        print(f"   ⚠️  Warning: Target file already exists: {new_filepath}")
        print(f"      Skipping: {filepath}")
        return False
    
    if dry_run:
        print(f"   Would rename: {filepath.name} -> {new_name}")
        return True
    else:
        try:
            filepath.rename(new_filepath)
            if verbose:
                print(f"   ✅ Renamed: {filepath.name} -> {new_name}")
            else:
                print(f"   {filepath.name} -> {new_name}")
            return True
        except Exception as e:
            print(f"   ❌ Error renaming {filepath}: {e}")
            return False


def process_directory(directory, dry_run=False, verbose=False, recursive=False):
    """
    Process all files in a directory to remove '.down' from filenames.
    
    Parameters
    ----------
    directory : str or Path
        Directory to process
    dry_run : bool, optional
        If True, only show what would be done
    verbose : bool, optional
        If True, show detailed output
    recursive : bool, optional
        If True, process subdirectories recursively
        
    Returns
    -------
    tuple
        (total_files_processed, files_renamed)
    """
    directory = Path(directory)
    
    if not directory.exists():
        print(f"❌ Directory does not exist: {directory}")
        return 0, 0
    
    if not directory.is_dir():
        print(f"❌ Not a directory: {directory}")
        return 0, 0
    
    print(f"Processing directory: {directory.absolute()}")
    if dry_run:
        print("🔍 DRY RUN MODE - No files will be actually renamed")
    print()
    
    total_files = 0
    renamed_files = 0
    
    # Get file pattern based on recursive option
    if recursive:
        pattern = "**/*"
        files = directory.glob(pattern)
        files = [f for f in files if f.is_file()]  # Filter only files
    else:
        files = [f for f in directory.iterdir() if f.is_file()]
    
    # Filter files that contain '.down'
    down_files = [f for f in files if '.down' in f.name]
    
    if not down_files:
        print("No files with '.down' found in the directory.")
        return 0, 0
    
    print(f"Found {len(down_files)} file(s) with '.down' in the name:")
    print()
    
    for filepath in sorted(down_files):
        total_files += 1
        if remove_down_from_filename(filepath, dry_run, verbose):
            renamed_files += 1
    
    return total_files, renamed_files


def main():
    """Main function to handle command line arguments and execute the script."""
    parser = argparse.ArgumentParser(
        description="Remove '.down' suffix from filenames in a directory",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Process current directory
    python remove_down_suffix.py
    
    # Process specific directory
    python remove_down_suffix.py -d /path/to/directory
    
    # Dry run to see what would be renamed
    python remove_down_suffix.py -n
    
    # Process recursively with verbose output
    python remove_down_suffix.py -r -v
    
    # Process test/down directory
    python remove_down_suffix.py -d test/down
        """
    )
    
    parser.add_argument(
        '-d', '--directory',
        default='.',
        help='Directory to process (default: current directory)'
    )
    
    parser.add_argument(
        '-n', '--dry-run',
        action='store_true',
        help='Show what would be renamed without actually doing it'
    )
    
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Show detailed output'
    )
    
    parser.add_argument(
        '-r', '--recursive',
        action='store_true',
        help='Process subdirectories recursively'
    )
    
    args = parser.parse_args()
    
    print("=== Remove '.down' Suffix from Filenames ===")
    print()
    
    try:
        total, renamed = process_directory(
            args.directory,
            dry_run=args.dry_run,
            verbose=args.verbose,
            recursive=args.recursive
        )
        
        print()
        print("=== Summary ===")
        print(f"Total files processed: {total}")
        print(f"Files renamed: {renamed}")
        
        if args.dry_run and renamed > 0:
            print()
            print("To actually rename the files, run the command without --dry-run")
        elif renamed > 0:
            print("✅ Operation completed successfully!")
        
    except KeyboardInterrupt:
        print("\n❌ Operation cancelled by user")
        sys.exit(1)
    except Exception as e:
        print(f"❌ Unexpected error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()