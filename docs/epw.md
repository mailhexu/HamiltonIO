# EPW Interface Documentation

The EPW (Electron-Phonon Wannier) interface in HamiltonIO provides tools for reading, processing, and analyzing electron-phonon coupling data from Quantum ESPRESSO's EPW code.

## Overview

The EPW module handles:
- **Wigner-Seitz vectors** for electrons, phonons, and electron-phonon interactions
- **EPW matrix elements** in both binary and NetCDF formats
- **Crystal structure** information from EPW calculations
- **Electron-phonon coupling** matrices and transformations

## Core Components

### WignerData Class

The `WignerData` class handles Wigner-Seitz grid information for EPW calculations.

```python
from HamiltonIO.epw.wigner import WignerData

# Read Wigner-Seitz data
data = WignerData.from_file("wigner.fmt")
print(data.summary())

# Access data
print(f"Electrons: {data.nrr_k} WS vectors")
print(f"Phonons: {data.nrr_q} WS vectors") 
print(f"Electron-phonon: {data.nrr_g} WS vectors")

# Write to file
data.to_file("output_wigner.fmt")
```

**Key Attributes:**
- `dims`: Number of Wannier functions
- `dims2`: Number of atoms
- `irvec_k/q/g`: WS vectors for electrons/phonons/electron-phonon
- `ndegen_k/q/g`: WS degeneracies
- `wslen_k/q/g`: WS vector lengths

### EPW Parser

The `epwparser` module provides comprehensive parsing of EPW output files.

```python
from HamiltonIO.epw.epwparser import Epmat, read_crystal_fmt

# Read crystal structure
crystal = read_crystal_fmt("crystal.fmt")
print(f"Atoms: {crystal.natom}, Modes: {crystal.nmode}")

# Read EPW matrix elements
epmat = Epmat()
epmat.read(path="./", prefix="material", epmat_ncfile="epmat.nc")
epmat.print_info()

# Extract single mode data
mode_data = EpmatOneMode(epmat, imode=0)
matrix = mode_data.get_epmat_RgRk(Rg=(0,0,0), Rk=(0,0,0))
```

### File Format Support

#### Binary Format (.epmatwp)
```python
from HamiltonIO.epw.epwparser import read_epmatwp

# Read binary EPW matrix elements
epmat_data = read_epmatwp("material.epmatwp", path="./")
print(f"Shape: {epmat_data.shape}")  # (nRg, nmodes, nRk, nwann, nwann)
```

#### NetCDF Format
```python
from HamiltonIO.epw.epwparser import save_epmat_to_nc

# Convert binary to NetCDF
save_epmat_to_nc(path="./", prefix="material", ncfile="epmat.nc")
```

## Command Line Tools

### EPW to NetCDF Converter

Convert EPW binary files to NetCDF format for easier analysis:

```bash
python -m HamiltonIO.epw.epwparser --convert-netcdf -p /path/to/epw -n material -o epmat.nc
```

**Options:**
- `-p, --path`: Directory containing EPW files
- `-n, --prefix`: EPW file prefix
- `-o, --output`: Output NetCDF filename

**Required Files:**
- `epwdata.fmt`: Basic dimensions
- `wigner.fmt`: Wigner-Seitz vectors
- `{prefix}.epmatwp`: EPW matrix elements

## Data Structures

### Crystal Structure
```python
@dataclass
class Crystal:
    natom: int          # Number of atoms
    nmode: int          # Number of phonon modes
    nelect: float       # Number of electrons
    at: np.ndarray      # Lattice vectors
    bg: np.ndarray      # Reciprocal lattice vectors
    omega: float        # Unit cell volume
    tau: np.ndarray     # Atomic positions
    amass: np.ndarray   # Atomic masses
    w_centers: np.ndarray  # Wannier centers
```

### EPW Matrix Elements
```python
@dataclass
class Epmat:
    crystal: Crystal
    nwann: int          # Number of Wannier functions
    nmodes: int         # Number of phonon modes
    nRk/nRq/nRg: int   # Number of R vectors
    Rk/Rq/Rg: np.ndarray  # R vector lists
    ndegen_k/q/g: np.ndarray  # Degeneracies
```

## Usage Examples

### Basic Workflow

```python
from HamiltonIO.epw.epwparser import Epmat
from HamiltonIO.epw.wigner import WignerData

# 1. Read Wigner-Seitz data
wigner = WignerData.from_file("wigner.fmt")
print(wigner.summary())

# 2. Load EPW matrix elements
epmat = Epmat()
epmat.read(path="./", prefix="SrMnO3", epmat_ncfile="epmat.nc")

# 3. Extract specific mode
mode_0 = EpmatOneMode(epmat, imode=0)

# 4. Get matrix elements for specific R vectors
g_matrix = mode_0.get_epmat_RgRk(Rg=(0,0,0), Rk=(1,0,0))
print(f"Matrix shape: {g_matrix.shape}")
```

### File Conversion

```python
from HamiltonIO.epw.epwparser import save_epmat_to_nc

# Convert binary EPW data to NetCDF
save_epmat_to_nc(
    path="/path/to/epw/data",
    prefix="material_name", 
    ncfile="epmat.nc"
)
```

### Data Analysis

```python
# Find origin vectors
k_orig, q_orig, g_orig = wigner.find_origin_indices()

# Access specific R vector data
origin_matrix = epmat.get_epmat_Rv_from_RgRk(
    imode=0, 
    Rg=(0,0,0), 
    Rk=(0,0,0)
)

# Time-reversal symmetry analysis
matrix_avg = mode_0.get_epmat_RgRk(
    Rg=(0,0,0), 
    Rk=(1,0,0), 
    avg=True  # Average with time-reversed
)
```

## File Formats

### wigner.fmt
Contains Wigner-Seitz vectors and degeneracies:
```
nrr_k nrr_q nrr_g dims dims2
# For each k-vector:
ix iy iz wslen
ndegen_matrix...
```

### epwdata.fmt
Basic EPW dimensions:
```
efermi
nwann nrr_k nmodes nrr_q nrr_g
```

### crystal.fmt
Crystal structure information:
```
natom
nmode
nelect [nbndskip]
lattice_vectors
reciprocal_vectors
volume
alat
atomic_positions
masses
types
flags
wannier_centers
cutoff_length
```

## Performance Notes

- **NetCDF format** is recommended for large datasets (faster I/O, compression)
- **Binary format** (.epmatwp) is the native EPW output but slower to read
- Use `EpmatOneMode` for memory-efficient single-mode analysis
- Wigner-Seitz data is typically small and fast to load

## Integration with HamiltonIO

The EPW interface integrates with the broader HamiltonIO ecosystem:

```python
# Convert to standard Hamiltonian format
from HamiltonIO.lcao_hamiltonian import LCAOHamiltonian

# EPW data can be used to construct electron-phonon Hamiltonians
# for transport and spectroscopy calculations
```

## Error Handling

Common issues and solutions:

- **Missing files**: Ensure all required EPW output files are present
- **Format mismatch**: Use `read_WSVec()` for new format, `read_WSVec_deprecated()` for old
- **Memory issues**: Use NetCDF format and process modes individually
- **Unit conversion**: Matrix elements are automatically converted to eV/Å units