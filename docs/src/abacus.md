# ABACUS Interface Documentation

The ABACUS interface in HamiltonIO provides tools for reading and analyzing electronic structure data from ABACUS DFT calculations, supporting both collinear and non-collinear spin configurations.

## Overview

The ABACUS interface handles:
- **ABACUS output files** (sparse CSR format, STRU files, log files)
- **Spin configurations** (non-polarized, collinear, non-collinear/SOC)
- **Real-space Hamiltonians** and overlap matrices in sparse format
- **Orbital basis information** and atomic structure

## Core Components

### AbacusWrapper Class

Main class for handling ABACUS Hamiltonian data, extending `LCAOHamiltonian`.

```python
from HamiltonIO.abacus import AbacusWrapper

# Create wrapper with parsed data
model = AbacusWrapper(
    HR=HR_matrix,
    SR=SR_matrix, 
    Rlist=R_vectors,
    nbasis=num_basis,
    atoms=atoms_object,
    nspin=1
)

# Basic properties
print(f"Number of basis functions: {model.nbasis}")
print(f"Number of R vectors: {len(model.Rlist)}")
print(f"Fermi energy: {model.efermi} eV")
```

### AbacusParser Class

Parser for ABACUS output files with automatic spin detection.

```python
from HamiltonIO.abacus import AbacusParser

# Initialize parser
parser = AbacusParser(
    outpath="./OUT.material/",  # ABACUS output directory
    spin=None,  # Auto-detect from log
    binary=False  # Text format CSR files
)

# Get models based on spin configuration
if parser.spin == "collinear":
    model_up, model_down = parser.get_models()
elif parser.spin == "noncollinear":
    model = parser.get_models()
```

## File Format Support

### Required ABACUS Output Files

#### Hamiltonian and Overlap (CSR format)
- **`data-HR-sparse_SPIN0.csr`**: Hamiltonian (spin-up or non-collinear)
- **`data-HR-sparse_SPIN1.csr`**: Hamiltonian (spin-down, collinear only)
- **`data-SR-sparse_SPIN0.csr`**: Overlap matrix

#### Structure and Metadata
- **`STRU`** or **`Stru`**: Atomic structure file
- **`running_scf.log`**: SCF log with Fermi energy and electron count
- **`Orbital`**: Orbital basis information

### Reading Examples

```python
# Automatic parsing
parser = AbacusParser(outpath="./OUT.Fe/")

# Check spin configuration
print(f"Detected spin: {parser.spin}")
print(f"Number of electrons: {parser.nel}")
print(f"Fermi energy: {parser.efermi} eV")

# Get atomic structure
atoms = parser.atoms
print(f"Chemical formula: {atoms.get_chemical_formula()}")
```

## Spin Configurations

### Non-Polarized Systems
```python
parser = AbacusParser(outpath="./OUT.material/", spin="non-polarized")
# Single model returned
model = parser.get_models()
```

### Collinear Spin-Polarized
```python
parser = AbacusParser(outpath="./OUT.Fe/", spin="collinear")
model_up, model_down = parser.get_models()

print(f"Spin-up basis: {model_up.nbasis}")
print(f"Spin-down basis: {model_down.nbasis}")

# Compare eigenvalues
evals_up, _ = model_up.solve([0, 0, 0])
evals_down, _ = model_down.solve([0, 0, 0])
```

### Non-Collinear/Spin-Orbit
```python
parser = AbacusParser(outpath="./OUT.material/", spin="noncollinear")
model = parser.get_models()

print(f"Spinor basis functions: {model.nbasis}")  # 2 × orbitals
print(f"Spin channels: {model.nspin}")  # Should be 2
```

## Advanced Features

### Split SOC Analysis

For calculations with separate SOC and non-SOC parts:

```python
from HamiltonIO.abacus import AbacusSplitSOCParser

# Parse both SOC and non-SOC calculations
parser = AbacusSplitSOCParser(
    outpath_nosoc="./OUT.nosoc/",
    outpath_soc="./OUT.soc/",
    binary=False
)

model = parser.parse()

# Access SOC components
print(f"SOC splitting available: {model.split_soc}")
if hasattr(model, 'HR_soc'):
    print("SOC Hamiltonian matrix available")
    print("Non-SOC Hamiltonian matrix available")
```

### Orbital Basis Information

```python
# Get orbital basis details
basis_info = parser.basis
print("Orbital basis:")
for i, orb in enumerate(basis_info):
    print(f"  {i}: Atom {orb.iatom}, Symbol {orb.sym}")

# Get formatted basis names
if parser.spin == "collinear":
    basis_up, basis_down = parser.get_basis()
    print("Spin-up orbitals:", basis_up[:5])  # First 5
    print("Spin-down orbitals:", basis_down[:5])
```

### Binary File Support

```python
# For binary CSR files (faster I/O)
parser = AbacusParser(
    outpath="./OUT.material/",
    binary=True  # Read binary format
)

models = parser.get_models()
```

## Data Structures

### AbacusWrapper Properties
```python
# Basic dimensions
model.nbasis      # Number of basis functions
model.norb        # Number of orbitals
model.nspin       # Number of spin channels
model.nel         # Number of electrons

# Electronic structure
model.efermi      # Fermi energy from SCF
model.atoms       # ASE Atoms object
model.basis       # Orbital basis information

# Hamiltonian data
model.HR          # Real-space Hamiltonian
model.SR          # Real-space overlap matrix
model.Rlist       # R-vector list
```

### File Format Details

#### CSR Sparse Matrix Format
```
Matrix Dimension: <nbasis>
Matrix number: <nR>
<R_x> <R_y> <R_z> <nnz>
<data_values>
<column_indices>
<row_pointers>
```

#### STRU File Format
```
ATOMIC_SPECIES
<element> <mass> <pseudopotential>

NUMERICAL_ORBITAL
<orbital_file>

LATTICE_CONSTANT
<lattice_constant>

LATTICE_VECTORS
<a1_x> <a1_y> <a1_z>
<a2_x> <a2_y> <a2_z>
<a3_x> <a3_y> <a3_z>

ATOMIC_POSITIONS
<coordinate_type>
<element>
<magnetism>
<number_of_atoms>
<x> <y> <z> <fix_x> <fix_y> <fix_z>
```

## Usage Examples

### Basic Electronic Structure

```python
from HamiltonIO.abacus import AbacusParser
import numpy as np

# 1. Parse ABACUS calculation
parser = AbacusParser(outpath="./OUT.material/")
model = parser.get_models()

# 2. Calculate band structure
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
print(f"Band energies shape: {bands.shape}")
```

### Collinear Spin Analysis

```python
# Parse collinear calculation
parser = AbacusParser(outpath="./OUT.Fe/")
model_up, model_down = parser.get_models()

# Compare spin channels at Gamma point
evals_up, _ = model_up.solve([0, 0, 0])
evals_down, _ = model_down.solve([0, 0, 0])

print("Spin splitting analysis:")
for i in range(min(10, len(evals_up), len(evals_down))):
    splitting = evals_up[i] - evals_down[i]
    print(f"  Band {i}: {splitting:.3f} eV")

# Magnetic moment analysis
total_moment = np.sum(evals_up < parser.efermi) - np.sum(evals_down < parser.efermi)
print(f"Total magnetic moment: {total_moment:.1f} μB")
```

### Non-Collinear/SOC Analysis

```python
# Parse non-collinear calculation
parser = AbacusParser(outpath="./OUT.soc/")
model = parser.get_models()

print(f"Spinor system with {model.nbasis} basis functions")

# Analyze eigenvalue spectrum
evals, evecs = model.solve([0, 0, 0])
print(f"Eigenvalue range: {evals.min():.3f} to {evals.max():.3f} eV")

# Check for Kramers degeneracy
kramers_pairs = []
for i in range(0, len(evals)-1, 2):
    degeneracy = abs(evals[i+1] - evals[i])
    kramers_pairs.append(degeneracy)
    
avg_degeneracy = np.mean(kramers_pairs)
print(f"Average Kramers degeneracy: {avg_degeneracy:.6f} eV")
```

### Real-Space Analysis

```python
# Examine hopping parameters
parser = AbacusParser(outpath="./OUT.material/")
model = parser.get_models()

# Get R-space matrices
HR_all = model.HR  # All R-vectors
print(f"Number of R vectors: {len(model.Rlist)}")

# Analyze nearest-neighbor hoppings
R_neighbors = [(1,0,0), (0,1,0), (0,0,1)]
for R in R_neighbors:
    R_idx = model.get_Ridx(R)
    H_R = HR_all[R_idx]
    max_hopping = np.max(np.abs(H_R))
    print(f"R={R}: max hopping = {max_hopping:.3f} eV")
```

### Structure Analysis

```python
# Analyze atomic structure
atoms = parser.atoms
print(f"Chemical formula: {atoms.get_chemical_formula()}")
print(f"Space group: {atoms.info.get('spacegroup', 'Unknown')}")

# Orbital composition
basis_info = parser.basis
orbital_count = {}
for orb in basis_info:
    symbol = atoms[orb.iatom].symbol
    if symbol not in orbital_count:
        orbital_count[symbol] = 0
    orbital_count[symbol] += 1

print("Orbitals per element:")
for element, count in orbital_count.items():
    print(f"  {element}: {count} orbitals")
```

## Integration with HamiltonIO

### Convert to Standard Format

```python
# ABACUS models are already LCAOHamiltonian subclasses
from HamiltonIO.lcao_hamiltonian import LCAOHamiltonian

# Direct usage
model = parser.get_models()
print(f"Model type: {type(model)}")  # AbacusWrapper(LCAOHamiltonian)

# Access parent class methods
H_k = model.get_Hk([0, 0, 0])
S_k = model.get_Sk([0, 0, 0])
```

### k-Space Calculations

```python
# Efficient multiple k-point calculations
from ase.dft.kpoints import monkhorst_pack

kpts = monkhorst_pack([6, 6, 6])
H_all, S_all, evals_all, evecs_all = model.HS_and_eigen(kpts)

print(f"k-points: {len(kpts)}")
print(f"Eigenvalues shape: {evals_all.shape}")
print(f"Energy range: {evals_all.min():.3f} to {evals_all.max():.3f} eV")
```

## Performance Notes

- **Sparse matrices**: ABACUS uses CSR format for memory efficiency
- **Binary files**: Use `binary=True` for faster I/O with large systems
- **Memory management**: Non-collinear systems require 4× memory for spinors
- **k-point calculations**: Vectorized operations for multiple k-points

## Common Issues

- **File paths**: Ensure output directory contains all required files
- **Spin detection**: Check `running_scf.log` for correct nspin value
- **Units**: Energies converted from Ry to eV automatically
- **Basis ordering**: Follows ABACUS orbital ordering conventions
- **SOC calculations**: Requires non-collinear ABACUS calculation
- **Memory**: Large systems benefit from binary format and chunked calculations