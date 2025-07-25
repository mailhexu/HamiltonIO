# SIESTA Interface Documentation

The SIESTA interface in HamiltonIO provides comprehensive tools for reading and analyzing electronic structure data from SIESTA DFT calculations using the sisl library as backend.

## Overview

The SIESTA interface handles:
- **SIESTA output files** (`.fdf`, `.nc`, `.HSX`, `.WFSX`)
- **Spin configurations** (non-polarized, collinear, non-collinear, spin-orbit)
- **Real-space Hamiltonians** and overlap matrices
- **Electronic structure analysis** with k-space transformations

## Core Components

### SiestaHamiltonian Class

The main class for handling SIESTA Hamiltonian data, extending `LCAOHamiltonian`.

```python
from HamiltonIO.siesta import SiestaHam

# Read from SIESTA calculation
model = SiestaHam("siesta.fdf", spin=0)  # spin-up channel

# Basic properties
print(f"Number of basis functions: {model.nbasis}")
print(f"Number of R vectors: {len(model.Rlist)}")
print(f"Spin configuration: {model.nspin}")

# Get Hamiltonian and overlap matrices
H_k = model.get_Hk([0, 0, 0])  # k-space Hamiltonian
S_k = model.get_Sk([0, 0, 0])  # k-space overlap matrix

# Solve eigenvalue problem
H, S, evals, evecs = model.HSE_k([0, 0, 0])
```

### SislParser Class

Low-level parser using sisl backend for flexible data extraction.

```python
from HamiltonIO.siesta import SislParser

# Initialize parser
parser = SislParser(
    fdf_fname="siesta.fdf",
    ispin=0,  # spin channel (0=up, 1=down, None=both)
    read_H_soc=True,  # read spin-orbit coupling
    orth=False  # don't orthogonalize
)

# Get model
model = parser.get_model()
```

## Spin Configurations

### Non-Polarized Systems
```python
# Single spin channel
model = SiestaHam("siesta.fdf")
print(f"Spin type: non-polarized")
print(f"Basis functions: {model.nbasis}")
```

### Collinear Spin-Polarized
```python
# Separate spin channels
model_up = SiestaHam("siesta.fdf", spin=0)    # spin-up
model_down = SiestaHam("siesta.fdf", spin=1)  # spin-down

# Or get both channels
parser = SislParser("siesta.fdf", ispin=None)
model_up, model_down = parser.get_model()
```

### Non-Collinear/Spin-Orbit
```python
# Spinor representation
model = SiestaHam("siesta.fdf", read_H_soc=True)
print(f"Basis functions (spinor): {model.nbasis}")  # 2 × orbitals

# Access SOC components
if hasattr(model, 'HR_soc'):
    print("SOC Hamiltonian available")
    print(f"Max SOC strength: {model.get_max_Hsoc_abs()}")
```

## File Format Support

### Input Files
- **`.fdf`**: Main SIESTA input file
- **`.nc`**: NetCDF output with Hamiltonian data
- **`.HSX`**: Hamiltonian and overlap matrices
- **`.WFSX`**: Wavefunctions (for SislWFSXWrapper)

### Reading Examples
```python
# From FDF file (recommended)
model = SiestaHam("calculation.fdf", spin=0)

# Direct NetCDF access
from HamiltonIO.siesta.mysiesta_nc import MySiestaNC
nc = MySiestaNC("siesta.nc")
qtot = nc.read_qtot()  # total charge
H_soc = nc.read_soc_hamiltonian()  # SOC Hamiltonian
```

## Advanced Features

### Spin-Orbit Coupling Analysis

```python
# Read with SOC splitting
parser = SislParser("siesta.fdf", read_H_soc=True)
model = parser.get_model()

# Access SOC and non-SOC parts
if model.split_soc:
    print("SOC Hamiltonian split available")
    
    # Rotate SOC quantization axis
    model.set_Hsoc_rotation_angle([np.pi/4, np.pi/6])  # theta, phi
    
    # Adjust SOC strength
    model.set_so_strength(0.5)  # 50% SOC strength
```

### Orbital Analysis

```python
# Get orbital information
orbs = model.orbs
print("Orbital labels:")
for i, orb in enumerate(orbs):
    print(f"  {i}: {orb}")

# Orbital dictionary by atom
parser = SislParser("siesta.fdf")
parser.get_model()
print("Orbitals by atom:")
for atom_idx, orb_list in parser.orb_dict.items():
    print(f"  Atom {atom_idx}: {orb_list}")
```

### k-Space Calculations

```python
# Single k-point
k_point = [0.25, 0.25, 0.0]
H, S, evals, evecs = model.HSE_k(k_point)

# Multiple k-points
from ase.dft.kpoints import monkhorst_pack
kpts = monkhorst_pack([4, 4, 4])

H_all, S_all, evals_all, evecs_all = model.HS_and_eigen(kpts)
print(f"Eigenvalues shape: {evals_all.shape}")  # (nk, nbasis)
```

### Real-Space Analysis

```python
# Get all R-space matrices
HR_all = model.get_all_HR()  # All R vectors
SR_all = model.get_SR_all()  # Overlap matrices

# Specific R vector
R_vector = (1, 0, 0)
HR_specific = model.get_HR(R=R_vector)
```

## Data Structures

### SiestaHamiltonian Properties
```python
# Dimensions
model.nbasis      # Number of basis functions
model.norb        # Number of orbitals
model.nspin       # Number of spin channels
model.nel         # Number of electrons

# Spatial data
model.atoms       # ASE Atoms object
model.Rlist       # R-vector list
model.positions   # Atomic positions
model.cell        # Unit cell

# Electronic structure
model.orbs        # Orbital labels
model.efermi      # Fermi energy (if available)
```

### Spin-Orbit Properties
```python
# SOC-specific attributes (when read_H_soc=True)
model.HR_soc      # SOC part of Hamiltonian
model.HR_nosoc    # Non-SOC part
model.split_soc   # Boolean flag
model.soc_rotation_angle  # [theta, phi] for quantization axis
model.so_strength # SOC strength factor
```

## Usage Examples

### Basic Electronic Structure

```python
from HamiltonIO.siesta import SiestaHam
import numpy as np

# 1. Load SIESTA calculation
model = SiestaHam("siesta.fdf", spin=0)

# 2. Calculate band structure along high-symmetry path
k_path = np.array([
    [0.0, 0.0, 0.0],  # Gamma
    [0.5, 0.0, 0.0],  # X
    [0.5, 0.5, 0.0],  # M
    [0.0, 0.0, 0.0]   # Gamma
])

bands = []
for k in k_path:
    evals, evecs = model.solve(k)
    bands.append(evals)

bands = np.array(bands)
print(f"Band structure shape: {bands.shape}")
```

### Spin-Polarized Analysis

```python
# Compare spin channels
model_up = SiestaHam("siesta.fdf", spin=0)
model_down = SiestaHam("siesta.fdf", spin=1)

# Density of states at Gamma point
evals_up, _ = model_up.solve([0, 0, 0])
evals_down, _ = model_down.solve([0, 0, 0])

print("Spin splitting at Gamma:")
for i in range(min(len(evals_up), len(evals_down))):
    splitting = evals_up[i] - evals_down[i]
    print(f"  Band {i}: {splitting:.3f} eV")
```

### Spin-Orbit Coupling Analysis

```python
# Load with SOC
model = SiestaHam("siesta.fdf", read_H_soc=True)

# Analyze SOC strength
if model.split_soc:
    soc_strength = model.get_max_Hsoc_abs()
    print(f"Maximum SOC matrix element: {soc_strength:.3f} eV")
    
    # Compare with/without SOC
    model.set_so_strength(0.0)  # Turn off SOC
    evals_nosoc, _ = model.solve([0, 0, 0])
    
    model.set_so_strength(1.0)  # Full SOC
    evals_soc, _ = model.solve([0, 0, 0])
    
    print("SOC effect on eigenvalues:")
    for i in range(min(len(evals_nosoc), len(evals_soc))):
        shift = evals_soc[i] - evals_nosoc[i]
        print(f"  Level {i}: {shift:.3f} eV")
```

### Wavefunction Analysis (WFSX)

```python
from HamiltonIO.siesta import SislWFSXWrapper
import sisl

# Read geometry and wavefunction data
geom = sisl.get_sile("siesta.fdf").read_geometry()
wfsx = sisl.get_sile("siesta.WFSX")
ham = sisl.get_sile("siesta.fdf").read_hamiltonian()

# Create wrapper
wrapper = SislWFSXWrapper(
    sisl_hamiltonian=ham,
    geom=geom,
    wfsx_sile=wfsx,
    spin=ham.spin
)

# Access wavefunctions
kpts = [[0, 0, 0], [0.25, 0, 0]]
evals, evecs = wrapper.find_all(kpts)
```

## Integration with HamiltonIO

### Convert to Standard Format

```python
# Convert to LCAOHamiltonian
from HamiltonIO.lcao_hamiltonian import LCAOHamiltonian

# Extract SIESTA data
HR = model.get_all_HR()
SR = model.get_SR_all() 
Rlist = model.Rlist

# Create LCAO model
lcao_model = LCAOHamiltonian(
    HR=HR,
    SR=SR,
    Rlist=Rlist,
    nbasis=model.nbasis,
    atoms=model.atoms,
    nspin=model.nspin
)
```

### Fermi Level Analysis

```python
# Calculate Fermi energy
from ase.dft.kpoints import monkhorst_pack

kpts = monkhorst_pack([8, 8, 8])
evals, _ = model.solve_all(kpts)

# Get Fermi level (if nel is available)
if hasattr(model, 'nel') and model.nel is not None:
    efermi = model.get_fermi_energy(evals, kweights=None)
    print(f"Fermi energy: {efermi:.3f} eV")
```

## Performance Notes

- **sisl backend**: Leverages optimized sisl library for file I/O
- **Memory management**: Large systems benefit from sparse matrix storage
- **k-point parallelization**: Use `solve_all()` for efficient multiple k-points
- **SOC calculations**: Spinor representation doubles memory requirements

## Common Issues

- **File paths**: Ensure `.nc` file is in same directory as `.fdf`
- **Spin indexing**: `spin=0` (up), `spin=1` (down), `spin=None` (both)
- **SOC compatibility**: Requires SIESTA compiled with SOC support
- **Basis ordering**: Orbital labels follow SIESTA conventions
- **Units**: Energies in eV, positions in Bohr (converted automatically)
- **Memory**: Large systems may require chunked k-point calculations