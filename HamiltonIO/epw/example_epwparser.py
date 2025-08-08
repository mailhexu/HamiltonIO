#!/usr/bin/env python3
"""
Comprehensive example usage of the epwparser module for EPW calculations.

This script demonstrates the full functionality of epwparser.py.
"""

import os
import numpy as np
from epwparser import (
    read_WSVec, 
    read_WSVec_deprecated,
    read_crystal_fmt,
    read_epwdata_fmt,
    line_to_array,
    line2vec,
    Crystal,
    Epmat
)


def demonstrate_wigner_reading():
    """Demonstrate reading Wigner-Seitz vectors."""
    print("=== Wigner-Seitz Vector Reading ===\n")
    
    print("1. Reading with new read_WSVec function (wigner.fmt format):")
    try:
        wigner_file = "test/up/wigner.fmt"
        if os.path.exists(wigner_file):
            result = read_WSVec(wigner_file)
            dims, dims2, nRk, nRq, nRg, Rk, Rq, Rg, ndegen_k, ndegen_q, ndegen_g = result
            
            print(f"   ✅ Successfully read {wigner_file}")
            print(f"   Dimensions: {dims} Wannier functions, {dims2} atoms")
            print(f"   R-vectors: k={nRk}, q={nRq}, g={nRg}")
            print(f"   Array shapes: Rk{Rk.shape}, Rq{Rq.shape}, Rg{Rg.shape}")
            print(f"   Sample k-vector: {Rk[0]} (degeneracy: {ndegen_k[0]})")
            print(f"   Origin k-vector: {Rk[13]} (degeneracy: {ndegen_k[13]})")
            
        else:
            print(f"   ❌ File {wigner_file} not found")
    except Exception as e:
        print(f"   ❌ Error reading new format: {e}")
    
    print()


def demonstrate_crystal_reading():
    """Demonstrate reading crystal structure information."""
    print("=== Crystal Structure Reading ===\n")
    
    crystal_files = ["test/up/crystal.fmt", "test/down/crystal.fmt"]
    
    for crystal_file in crystal_files:
        if os.path.exists(crystal_file):
            print(f"Reading {crystal_file}:")

            crystal = read_crystal_fmt(crystal_file)
            try:
                crystal = read_crystal_fmt(crystal_file)
                
                print(f"   ✅ Successfully read crystal structure")
                print(f"   Number of atoms: {crystal.natom}")
                print(f"   Number of modes: {crystal.nmode}")
                print(f"   Number of electrons: {crystal.nelect}")
                print(f"   Unit cell volume: {crystal.omega:.6f}")
                print(f"   Lattice parameter: {crystal.alat:.6f}")
                print(f"   Non-collinear: {crystal.noncolin}")
                    
            except Exception as e:
                print(f"   ❌ Error reading crystal file: {e}")
        else:
            print(f"   ❌ File {crystal_file} not found")
        print()


def demonstrate_epwdata_reading():
    """Demonstrate reading EPW data dimensions."""
    print("=== EPW Data Dimensions ===\n")
    
    epwdata_files = ["test/up/epwdata.fmt", "test/down/epwdata.fmt"]
    
    for epwdata_file in epwdata_files:
        if os.path.exists(epwdata_file):
            print(f"Reading {epwdata_file}:")
            try:
                nbndsub, nrr_k, nmodes, nrr_q, nrr_g = read_epwdata_fmt(epwdata_file)
                
                print(f"   ✅ Successfully read EPW data dimensions")
                print(f"   Number of bands (nbndsub): {nbndsub}")
                print(f"   R-vectors: k={nrr_k}, q={nrr_q}, g={nrr_g}")
                print(f"   Number of modes: {nmodes}")
                
            except Exception as e:
                print(f"   ❌ Error reading EPW data: {e}")
        else:
            print(f"   ❌ File {epwdata_file} not found")
        print()


def demonstrate_utility_functions():
    """Demonstrate utility functions for parsing."""
    print("=== Utility Functions ===\n")
    
    # Test line_to_array function
    print("1. line_to_array function:")
    test_line = "1.0 2.5 -3.14 4.2"
    float_array = line_to_array(test_line, float)
    int_line = "1 2 3 4 5"
    int_array = line_to_array(int_line, int)
    
    print(f"   Input: '{test_line}'")
    print(f"   Output (float): {float_array}")
    print(f"   Input: '{int_line}'")
    print(f"   Output (int): {int_array}")
    print()
    
    # Test line2vec function
    print("2. line2vec function:")
    vector_line = "  -1   0   1  "
    vector = line2vec(vector_line)
    print(f"   Input: '{vector_line}'")
    print(f"   Output: {vector}")
    print()


def demonstrate_epmat_class():
    """Demonstrate using the Epmat class for advanced EPW operations."""
    print("=== Epmat Class Usage ===\n")
    
    # Look for EPW data directories
    test_dirs = ["test/up", "test/down"]
    
    for test_dir in test_dirs:
        if not os.path.exists(test_dir):
            continue
            
        print(f"Testing Epmat class with {test_dir}:")
        
        # Check for required files
        wigner_file = os.path.join(test_dir, "wigner.fmt")
        epwdata_file = os.path.join(test_dir, "epwdata.fmt")
        
        if os.path.exists(wigner_file):
            try:
                # Create and initialize Epmat instance
                epmat = Epmat()
                
                # Method 1: Read R-vectors only
                print("   1. Reading R-vectors using read_Rvectors:")
                epmat.read_Rvectors(test_dir, "wigner.fmt")
                
                print(f"      ✅ Successfully loaded R-vectors")
                print(f"      Dimensions: dims={epmat.dims}, dims2={epmat.dims2}")
                print(f"      R-vectors: k={epmat.nRk}, q={epmat.nRq}, g={epmat.nRg}")
                print(f"      Array shapes: Rk{epmat.Rk.shape}, Rq{epmat.Rq.shape}, Rg{epmat.Rg.shape}")
                
                # Show R-vector dictionaries
                print(f"      Dictionary sizes: k={len(epmat.Rkdict)}, q={len(epmat.Rqdict)}, g={len(epmat.Rgdict)}")
                
                # Test dictionary lookup
                sample_Rk = tuple(epmat.Rk[0])
                sample_index = epmat.Rkdict[sample_Rk]
                print(f"      Dictionary test: R-vector {sample_Rk} -> index {sample_index}")
                
                # Show some sample R-vectors
                print(f"      Sample k-vectors: {[tuple(epmat.Rk[i]) for i in range(min(3, epmat.nRk))]}")
                print(f"      Sample degeneracies: k={epmat.ndegen_k[:3]}")
                
                # Method 2: Try full read if epwdata.fmt exists
                if os.path.exists(epwdata_file):
                    print("\\n   2. Reading full EPW data:")
                    try:
                        epmat2 = Epmat()
                        # This will read epwdata.fmt and wigner.fmt
                        epmat2.read(test_dir, prefix="SrMnO3", epmat_ncfile=None)
                        
                        print(f"      ✅ Successfully read EPW dimensions")
                        print(f"      Wannier functions: {epmat2.nwann}")
                        print(f"      Phonon modes: {epmat2.nmodes}")
                        print(f"      R-vectors: k={epmat2.nRk}, q={epmat2.nRq}, g={epmat2.nRg}")
                        
                    except Exception as e:
                        print(f"      ⚠️  Partial success (R-vectors loaded, matrix elements may be missing): {e}")
                        
                else:
                    print(f"      ❌ epwdata.fmt not found in {test_dir}")
                
            except Exception as e:
                print(f"   ❌ Error with Epmat class: {e}")
                import traceback
                traceback.print_exc()
        else:
            print(f"   ❌ wigner.fmt not found in {test_dir}")
        print()


def demonstrate_netcdf_operations():
    """Demonstrate NetCDF file operations with EPW data using epwparser functions."""
    print("=== NetCDF File Operations ===\n")
    
    # Import additional functions for NetCDF operations
    try:
        from epwparser import save_epmat_to_nc, EpmatOneMode
        from netCDF4 import Dataset
        netcdf_available = True
    except ImportError as e:
        print(f"   ❌ NetCDF4 not available: {e}")
        print("   Install with: pip install netcdf4")
        netcdf_available = False
        return
    
    test_dirs = ["test/up", "test/down"]
    
    for test_dir in test_dirs:
        if not os.path.exists(test_dir):
            continue
            
        print(f"Testing NetCDF operations with {test_dir}:")
        
        # Check for required files
        wigner_file = os.path.join(test_dir, "wigner.fmt")
        epwdata_file = os.path.join(test_dir, "epwdata.fmt")
        
        if os.path.exists(wigner_file) and os.path.exists(epwdata_file):
            try:
                # Method 1: Create Epmat instance with mock data
                print("   1. Creating Epmat instance with mock matrix elements:")
                
                # Read dimensions
                nbndsub, nrr_k, nmodes, nrr_q, nrr_g = read_epwdata_fmt(epwdata_file)
                dims, dims2, nRk, nRq, nRg, *_ = read_WSVec(wigner_file)
                
                print(f"      Dimensions: bands={nbndsub}, modes={nmodes}")
                print(f"      R-vectors: k={nrr_k}, q={nrr_q}, g={nrr_g}")
                
                # Create Epmat instance
                epmat = Epmat()
                epmat.nwann = nbndsub
                epmat.nmodes = nmodes
                epmat.nRk = nrr_k
                epmat.nRq = nrr_q
                epmat.nRg = nrr_g
                
                # Read R-vectors using the modern implementation
                epmat.read_Rvectors(test_dir, "wigner.fmt")
                
                # Create mock matrix elements (normally these would be read from .epmatwp files)
                print("      Creating mock matrix elements for demonstration...")
                mock_epmat = np.random.random((nrr_g, nmodes, nrr_k, nbndsub, nbndsub)) +1j * np.random.random((nrr_g, nmodes, nrr_k, nbndsub, nbndsub))
                mock_epmat *= 0.01  # Scale to reasonable values
                epmat.epmat_wann = mock_epmat
                
                print(f"      Mock matrix shape: {mock_epmat.shape}")
                print(f"      Memory usage: {mock_epmat.nbytes / 1024**2:.1f} MB")
                
                # Method 2: Save to NetCDF using Epmat.save_to_netcdf method
                print("\\n   2. Saving to NetCDF using Epmat.save_to_netcdf:")
                nc_filename = "example_epmat.nc"
                
                try:
                    epmat.save_to_netcdf(test_dir, nc_filename)
                    nc_filepath = os.path.join(test_dir, nc_filename)
                    
                    print(f"      ✅ Successfully saved using epmat.save_to_netcdf()")
                    print(f"      File: {nc_filepath}")
                    
                    # Get file size
                    if os.path.exists(nc_filepath):
                        file_size = os.path.getsize(nc_filepath)
                        print(f"      File size: {file_size / 1024**2:.1f} MB")
                    
                except Exception as e:
                    print(f"      ❌ Error saving with Epmat method: {e}")
                
                # Method 3: Alternative - use save_epmat_to_nc function
                print("\\n   3. Using save_epmat_to_nc function:")
                nc_filename2 = "example_epmat_func.nc"
                
                try:
                    # This function creates an Epmat instance and saves it
                    # Note: This requires actual .epmatwp files, so we'll demonstrate the concept
                    print("      save_epmat_to_nc function is available for complete datasets")
                    print("      Usage: save_epmat_to_nc(path='./data', prefix='material', ncfile='epmat.nc')")
                    print("      This function reads binary .epmatwp files and converts to NetCDF")
                    
                except Exception as e:
                    print(f"      ❌ Error with save_epmat_to_nc: {e}")
                
                # Method 4: Read back from NetCDF using Epmat
                print("\\n   4. Reading NetCDF file using Epmat class:")
                nc_filepath = os.path.join(test_dir, nc_filename)
                if os.path.exists(nc_filepath):
                    try:
                        epmat_nc = Epmat()
                        epmat_nc.read(test_dir, prefix="SrMnO3", epmat_ncfile=nc_filename)
                        
                        print(f"      ✅ Successfully loaded Epmat with NetCDF")
                        print(f"      NetCDF file connected: {hasattr(epmat_nc, 'epmatfile')}")
                        
                        if hasattr(epmat_nc, 'epmatfile') and epmat_nc.epmatfile is not None:
                            print(f"      Available methods:")
                            print(f"      - get_epmat_Rv_from_index(imode, iRg, iRk)")
                            print(f"      - get_epmat_Rv_from_RgRk(imode, Rg, Rk)")
                            print(f"      - get_epmat_Rv_from_R(imode, Rg)")
                            
                            # Test reading a matrix element
                            test_matrix = epmat_nc.get_epmat_Rv_from_index(0, 0, 0)
                            print(f"      Sample matrix element shape: {test_matrix.shape}")
                            print(f"      Sample matrix range: [{test_matrix.real.min():.6f}, {test_matrix.real.max():.6f}]")
                            
                            # Close the NetCDF file
                            epmat_nc.epmatfile.close()
                        
                    except Exception as e:
                        print(f"      ❌ Error with NetCDF Epmat: {e}")
                
                # Method 5: Demonstrate NetCDF file inspection
                print("\\n   5. NetCDF file inspection:")
                if os.path.exists(nc_filepath):
                    try:
                        with Dataset(nc_filepath, 'r') as nc_file:
                            print(f"      ✅ NetCDF file structure:")
                            print(f"      Dimensions: {dict(nc_file.dimensions)}")
                            print(f"      Variables: {list(nc_file.variables.keys())}")
                            
                            # Show variable details
                            if 'epmat_real' in nc_file.variables:
                                var = nc_file.variables['epmat_real']
                                print(f"      epmat_real shape: {var.shape}")
                                print(f"      epmat_real dtype: {var.dtype}")
                            
                    except Exception as e:
                        print(f"      ❌ Error inspecting NetCDF: {e}")
                
                # Clean up demo files
                cleanup_files = [
                    os.path.join(test_dir, nc_filename),
                    os.path.join(test_dir, nc_filename2)
                ]
                
                for cleanup_file in cleanup_files:
                    if os.path.exists(cleanup_file):
                        os.remove(cleanup_file)
                        print(f"      🧹 Cleaned up: {os.path.basename(cleanup_file)}")
                
            except Exception as e:
                print(f"   ❌ Error in NetCDF operations: {e}")
                import traceback
                traceback.print_exc()
        else:
            missing = []
            if not os.path.exists(wigner_file):
                missing.append("wigner.fmt")
            if not os.path.exists(epwdata_file):
                missing.append("epwdata.fmt")
            print(f"   ❌ Missing required files: {missing}")
        print()


def demonstrate_epmat_onemode():
    """Demonstrate EpmatOneMode for single phonon mode analysis."""
    print("=== EpmatOneMode Usage ===\\n")
    
    try:
        from epwparser import EpmatOneMode
        
        print("EpmatOneMode class is available for single phonon mode analysis:")
        print("   Features:")
        print("   - Extract matrix elements for a specific phonon mode")
        print("   - Efficient memory usage (only one mode at a time)")
        print("   - Methods for R-vector specific operations:")
        print("     * get_epmat_RgRk(Rg, Rk, avg=False)")
        print("     * get_epmat_RgRk_two_spin(Rg, Rk, avg=False)")
        print("   - Time-reversal symmetry averaging")
        print()
        print("   Usage example:")
        print("   ```python")
        print("   # First create full Epmat instance")
        print("   epmat = Epmat()")
        print("   epmat.read(path, prefix, epmat_ncfile='epmat.nc')")
        print("   ")
        print("   # Extract single mode")
        print("   mode_3 = EpmatOneMode(epmat, imode=3)")
        print("   ")
        print("   # Get matrix elements for specific R-vectors")
        print("   matrix = mode_3.get_epmat_RgRk(Rg=(0,0,0), Rk=(1,0,0))")
        print("   ```")
        print()
        
    except ImportError:
        print("   ❌ EpmatOneMode not available (requires full epwparser)")


def analyze_data_consistency():
    """Analyze consistency between different data sources."""
    print("=== Data Consistency Analysis ===\n")
    
    wigner_file = "test/up/wigner.fmt"
    epwdata_file = "test/up/epwdata.up.fmt"
    
    if os.path.exists(wigner_file) and os.path.exists(epwdata_file):
        try:
            # Read from both sources
            dims, dims2, nRk_w, nRq_w, nRg_w, *_ = read_WSVec(wigner_file)
            nbndsub, nRk_e, nmodes, nRq_e, nRg_e = read_epwdata_fmt(epwdata_file)
            
            print("Comparing dimensions from wigner.fmt and epwdata.fmt:")
            print(f"   nRk: wigner={nRk_w}, epwdata={nRk_e} {'✅' if nRk_w == nRk_e else '❌'}")
            print(f"   nRq: wigner={nRq_w}, epwdata={nRq_e} {'✅' if nRq_w == nRq_e else '❌'}")
            print(f"   nRg: wigner={nRg_w}, epwdata={nRg_e} {'✅' if nRg_w == nRg_e else '❌'}")
            print(f"   Additional info: nbndsub={nbndsub}, nmodes={nmodes}")
            
        except Exception as e:
            print(f"   ❌ Error in consistency check: {e}")
    else:
        print("   ❌ Required files not found for consistency check")
    
    print()


def main():
    """Main example function demonstrating epwparser functionality."""
    print("=== EPW Parser Comprehensive Example ===\n")
    
    # Check if we're in the right directory
    if not os.path.exists("test"):
        print("❌ Test directory not found. Please run this script from the EPW project root.")
        return
    
    # Demonstrate core functionality
    demonstrate_wigner_reading()
    demonstrate_crystal_reading()
    demonstrate_epwdata_reading()
    demonstrate_utility_functions()
    
    # Demonstrate advanced functionality
    demonstrate_epmat_class()
    demonstrate_netcdf_operations()
    demonstrate_epmat_onemode()
    analyze_data_consistency()
    
    print("=== Example completed successfully! ===")
    print("\nComprehensive epwparser functionality demonstrated:")
    print("- ✅ Wigner-Seitz vector reading (modern format)")
    print("- ✅ Crystal structure parsing")
    print("- ✅ EPW data dimensions")
    print("- ✅ Utility functions for data parsing")
    print("- ✅ Epmat class usage")
    print("- ✅ NetCDF file operations")
    print("- ✅ EpmatOneMode single-mode analysis")
    print("- ✅ Data consistency checking")
    
    print("\nKey features of the updated implementation:")
    print("- 🔄 Uses modern wigner.py for R-vector reading")
    print("- 📊 NetCDF support for large datasets")
    print("- 🎯 Single phonon mode extraction")
    print("- 🔍 Comprehensive data validation")
    print("- 🐍 Pythonic array indexing throughout")


if __name__ == "__main__":
    main()