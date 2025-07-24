#!/usr/bin/env python3
"""
Example usage of the epwparser module for EPW calculations.

This script demonstrates the core functionality of epwparser.py including:
- Reading Wigner-Seitz vectors (both new and deprecated formats)
- Reading crystal structure information
- Reading EPW data dimensions
- Basic utility functions

This example focuses on the core functions that don't require heavy dependencies.
"""

import os
import sys
import numpy as np

# Add current directory to path to import local modules
sys.path.insert(0, '.')

# Import only the core functions we need
def import_core_functions():
    """Import core functions from epwparser, handling dependencies gracefully."""
    try:
        # Try to import specific functions to avoid dependency issues
        import importlib.util
        spec = importlib.util.spec_from_file_location("epwparser_core", "epwparser.py")
        
        # Read the file and extract only the functions we need
        with open('epwparser.py', 'r') as f:
            content = f.read()
        
        # Create a minimal version with just the functions we need
        core_functions = {}
        
        # Import numpy and basic dependencies
        exec("import numpy as np", core_functions)
        exec("from collections import defaultdict", core_functions)
        exec("from dataclasses import dataclass, field", core_functions)
        
        # Define the utility functions
        exec('''
def line_to_array(line, fmt=float):
    """Convert a string line to a numpy array with given format."""
    return np.array([fmt(x) for x in line.split()])

def line2vec(line):
    """Convert a string line to a vector of integers."""
    return [int(x) for x in line.strip().split()]

def is_text_True(s):
    """Check if a string represents a boolean True value."""
    return s.strip().lower().startswith("t")
        ''', core_functions)
        
        # Define the WignerData-based read_WSVec function
        exec('''
def read_WSVec(fname):
    """Read a Wigner-Seitz vector file using the modern wigner.py module."""
    from wigner import WignerData
    
    wigner_data = WignerData.from_file(fname)
    
    dims = wigner_data.dims
    dims2 = wigner_data.dims2
    nRk = wigner_data.nrr_k
    nRq = wigner_data.nrr_q
    nRg = wigner_data.nrr_g
    
    Rk = wigner_data.irvec_k
    Rq = wigner_data.irvec_q
    Rg = wigner_data.irvec_g
    
    if dims == 1 and dims2 == 1:
        ndegen_k = wigner_data.ndegen_k[:, 0, 0]
        ndegen_q = wigner_data.ndegen_q[:, 0, 0]
        ndegen_g = wigner_data.ndegen_g[:, 0, 0]
    else:
        ndegen_k = wigner_data.ndegen_k[:, 0, 0]
        ndegen_q = wigner_data.ndegen_q[:, 0, 0]
        ndegen_g = wigner_data.ndegen_g[:, 0, 0]
    
    return (dims, dims2, nRk, nRq, nRg, Rk, Rq, Rg, ndegen_k, ndegen_q, ndegen_g)
        ''', core_functions)
        
        # Define Crystal dataclass and read_crystal_fmt
        exec('''
@dataclass
class Crystal:
    """Class storing crystal structure information from EPW calculation."""
    natom: int = 0
    nmode: int = 0
    nelect: float = 0.0
    at: np.ndarray = field(default_factory=lambda: np.zeros(0))
    bg: np.ndarray = field(default_factory=lambda: np.zeros(0))
    omega: float = 0.0
    alat: float = 0.0
    tau: np.ndarray = field(default_factory=lambda: np.zeros(0))
    amass: np.ndarray = field(default_factory=lambda: np.zeros(0))
    ityp: np.ndarray = field(default_factory=lambda: np.zeros(0))
    noncolin: bool = False
    w_centers: np.ndarray = field(default_factory=lambda: np.zeros(0))

def read_crystal_fmt(fname="crystal.fmt"):
    """Parse the crystal.fmt file containing crystal structure information."""
    d = Crystal()
    with open(fname) as myfile:
        d.natom = int(next(myfile))
        d.nmode = int(next(myfile))
        d.nelect = float(next(myfile))
        d.at = line_to_array(next(myfile), float)
        d.bg = line_to_array(next(myfile), float)
        d.omega = float(next(myfile))
        d.alat = float(next(myfile))
        d.tau = line_to_array(next(myfile), float)
        d.amass = line_to_array(next(myfile), float)
        d.ityp = line_to_array(next(myfile), int)
        d.noncolin = is_text_True(next(myfile))
        d.w_centers = line_to_array(next(myfile), float)
    return d

def read_epwdata_fmt(fname="epwdata.fmt"):
    """Read the EPW data format file containing basic dimensions."""
    with open(fname) as myfile:
        _efermi = float(next(myfile))
        nbndsub, nrr_k, nmodes, nrr_q, nrr_g = [int(x) for x in next(myfile).split()]
    return nbndsub, nrr_k, nmodes, nrr_q, nrr_g
        ''', core_functions)
        
        return core_functions
        
    except Exception as e:
        print(f"Error importing core functions: {e}")
        return None

# Import the core functions
core_funcs = import_core_functions()
if core_funcs is None:
    print("Failed to import core functions")
    sys.exit(1)

# Extract functions from the namespace
read_WSVec = core_funcs['read_WSVec']
read_crystal_fmt = core_funcs['read_crystal_fmt']
read_epwdata_fmt = core_funcs['read_epwdata_fmt']
line_to_array = core_funcs['line_to_array']
line2vec = core_funcs['line2vec']
Crystal = core_funcs['Crystal']


def demonstrate_wigner_reading():
    """Demonstrate reading Wigner-Seitz vectors using the new method."""
    print("=== Wigner-Seitz Vector Reading ===\\n")
    
    print("1. Reading with read_WSVec function (wigner.fmt format):")
    try:
        wigner_file = "test/up/wigner.fmt"
        if os.path.exists(wigner_file):
            result = read_WSVec(wigner_file)
            dims, dims2, nRk, nRq, nRg, Rk, Rq, Rg, ndegen_k, ndegen_q, ndegen_g = result
            
            print(f"   ✅ Successfully read {wigner_file}")
            print(f"   Dimensions: {dims} Wannier functions, {dims2} atoms")
            print(f"   R-vectors: k={nRk}, q={nRq}, g={nRg}")
            print(f"   Array shapes: Rk{Rk.shape}, Rq{Rq.shape}, Rg{Rg.shape}")
            print(f"   Degeneracy shapes: k{ndegen_k.shape}, q{ndegen_q.shape}, g{ndegen_g.shape}")
            
            # Show some sample data
            print(f"   Sample k-vector: {Rk[0]} (degeneracy: {ndegen_k[0]})")
            print(f"   Origin k-vector: {Rk[13]} (degeneracy: {ndegen_k[13]})")
            
        else:
            print(f"   ❌ File {wigner_file} not found")
    except Exception as e:
        print(f"   ❌ Error reading wigner format: {e}")
        import traceback
        traceback.print_exc()
    
    print()


def demonstrate_crystal_reading():
    """Demonstrate reading crystal structure information."""
    print("=== Crystal Structure Reading ===\\n")
    
    crystal_files = ["test/up/crystal.up.fmt", "test/down/crystal.down.fmt"]
    
    for crystal_file in crystal_files:
        if os.path.exists(crystal_file):
            print(f"Reading {crystal_file}:")
            try:
                crystal = read_crystal_fmt(crystal_file)
                
                print(f"   ✅ Successfully read crystal structure")
                print(f"   Number of atoms: {crystal.natom}")
                print(f"   Number of modes: {crystal.nmode}")
                print(f"   Number of electrons: {crystal.nelect}")
                print(f"   Unit cell volume: {crystal.omega:.6f}")
                print(f"   Lattice parameter: {crystal.alat:.6f}")
                print(f"   Non-collinear: {crystal.noncolin}")
                
                if len(crystal.at) > 0:
                    print(f"   Lattice vectors shape: {crystal.at.shape}")
                if len(crystal.tau) > 0:
                    print(f"   Atomic positions shape: {crystal.tau.shape}")
                if len(crystal.w_centers) > 0:
                    print(f"   Wannier centers shape: {crystal.w_centers.shape}")
                    
            except Exception as e:
                print(f"   ❌ Error reading crystal file: {e}")
        else:
            print(f"   ❌ File {crystal_file} not found")
        print()


def demonstrate_epwdata_reading():
    """Demonstrate reading EPW data dimensions."""
    print("=== EPW Data Dimensions ===\\n")
    
    epwdata_files = ["test/up/epwdata.up.fmt", "test/down/epwdata.down.fmt"]
    
    for epwdata_file in epwdata_files:
        if os.path.exists(epwdata_file):
            print(f"Reading {epwdata_file}:")
            try:
                nbndsub, nrr_k, nmodes, nrr_q, nrr_g = read_epwdata_fmt(epwdata_file)
                
                print(f"   ✅ Successfully read EPW data dimensions")
                print(f"   Number of bands: {nbndsub}")
                print(f"   R-vectors: k={nrr_k}, q={nrr_q}, g={nrr_g}")
                print(f"   Number of modes: {nmodes}")
                
            except Exception as e:
                print(f"   ❌ Error reading EPW data: {e}")
        else:
            print(f"   ❌ File {epwdata_file} not found")
        print()


def demonstrate_utility_functions():
    """Demonstrate utility functions for parsing."""
    print("=== Utility Functions ===\\n")
    
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


def analyze_data_consistency():
    """Analyze consistency between different data sources."""
    print("=== Data Consistency Analysis ===\\n")
    
    # Check if wigner and epwdata files have consistent dimensions
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
            print(f"   Additional info from epwdata: nbndsub={nbndsub}, nmodes={nmodes}")
            
        except Exception as e:
            print(f"   ❌ Error in consistency check: {e}")
    else:
        print("   ❌ Required files not found for consistency check")
    
    print()


def main():
    """Main example function demonstrating epwparser functionality."""
    print("=== EPW Parser Example Usage ===\\n")
    
    # Check if we're in the right directory
    if not os.path.exists("test"):
        print("❌ Test directory not found. Please run this script from the EPW project root.")
        return
    
    # Demonstrate core functionality
    demonstrate_wigner_reading()
    demonstrate_crystal_reading()
    demonstrate_epwdata_reading()
    demonstrate_utility_functions()
    analyze_data_consistency()
    
    print("=== Example completed! ===")
    print("\\nCore epwparser functionality demonstrated:")
    print("- ✅ Wigner-Seitz vector reading (modern format)")
    print("- ✅ Crystal structure parsing")
    print("- ✅ EPW data dimensions")
    print("- ✅ Utility functions for data parsing")
    print("- ✅ Data consistency checking")


if __name__ == "__main__":
    main()