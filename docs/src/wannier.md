# Wannier90 Interface Documentation

The Wannier90 interface in HamiltonIO provides tools for reading and manipulating tight-binding Hamiltonians from Wannier90 calculations.

## Overview

The Wannier interface handles:
- **Wannier90 output files** (`_hr.dat`, `_tb.dat`, `.win`, `.xyz`)
- **Tight-binding Hamiltonians** in real space with R-vector indexing
- **Wannier function centers** and orbital assignments
- **k-space transformations** and eigenvalue calculations

## Core Components

### WannierHam Class

The `WannierHam` class is the main interface for Wannier90 tight-binding models.

```python
from HamiltonIO.wannier import WannierHam

# Read from Wannier90 directory
model = WannierHam.read_from_wannier_dir(
    path="./wannier_calc/", 
    prefix="wannier90",
    atoms=atoms
)

# Access basic properties
print(f"Number of basis functions: {model.nbasis}")
print(f"Number of R vectors: {model.nR}")
print(f"R vectors: {model.Rlist}")

# Get Hamiltonian at specific R
H_R0 = model.get_hamR((0, 0, 0))  # On-site terms
H_R1 = model.get_hamR((1, 0, 0))  # Nearest neighbor

# Generate k-space Hamiltonian
k_point = [0.0, 0.0, 0.0]  # Gamma point
H_k = model.gen_ham(k_point)

# Solve eigenvalue problem
eigenvals, eigenvecs = model.solve(k_point)
```

### File Parsers

#### Hamiltonian Parser (`_hr.dat`)
```python
from HamiltonIO.wannier.w90_parser import parse_ham

# Parse Hamiltonian file
n_wann, H_mnR, R_degens = parse_ham("wannier90_hr.dat")

print(f"Number of Wannier functions: {n_wann}")
print(f"R-space Hamiltonian keys: {list(H_mnR.keys())}")
print(f"Degeneracies: {R_degens}")
```

#### Tight-Binding Parser (`_tb.dat`)
```python
from HamiltonIO.wannier.w90_parser import parse_tb

# Parse tight-binding file (includes centers)
centers, n_wann, H_mnR, R_degens = parse_tb("wannier90_tb.dat")

print(f"Wannier centers shape: {centers.shape}")
```

#### Structure Parser (`.win`)
```python
from HamiltonIO.wannier.w90_parser import parse_atoms, parse_cell

# Parse atomic structure
atoms = parse_atoms("wannier90.win")
print(f"Chemical symbols: {atoms.get_chemical_symbols()}")

# Parse unit cell
cell = parse_cell("wannier90.win")
print(f"Lattice vectors:\n{cell}")
```

#### Wannier Centers (`.xyz`)
```python
from HamiltonIO.wannier.w90_parser import parse_xyz

# Parse Wannier centers and atomic positions
wann_pos, atom_symbols, atom_pos = parse_xyz("wannier90.xyz")

print(f"Wannier positions shape: {wann_pos.shape}")
print(f"Atomic symbols: {atom_symbols}")
```

## Advanced Features

### Spin-Orbit Coupling Support

```python
# Read spin-orbit coupled system
model = WannierHam.read_from_wannier_dir(
    path="./soc_calc/",
    prefix="wannier90",
    groupby="spin"  # Reorder for spinor basis
)

# Convert to spin-polarized representation
model_spin = model.to_spin_polarized(order=1)  # Interleaved spin
```

### Position Shifting

```python
# Shift Wannier functions to atomic positions
reference_positions = atoms.get_scaled_positions()
shifted_model = model.shift_position(reference_positions)
```

### Basis Assignment

```python
from HamiltonIO.wannier.utils import auto_assign_basis_name

# Automatically assign orbital names
basis_dict, shifted_pos = auto_assign_basis_name(
    positions=model.positions,
    atoms=atoms,
    max_distance=0.1
)

print("Basis assignments:")
for name, index in basis_dict.items():
    print(f"  {name}: {index}")
```

## Data Structures

### WannierHam Properties
```python
# Basic dimensions
model.nbasis      # Number of basis functions
model.norb        # Number of orbitals (nbasis/nspin)
model.nspin       # Number of spin channels
model.nR          # Number of R vectors

# Spatial information
model.positions   # Wannier function positions
model.Rlist       # List of R vectors
model.atoms       # ASE Atoms object

# Hamiltonian data
model.data        # Dictionary {R: H(R)} 
model.onsite_energies  # Diagonal elements at R=0
model.hoppings    # Off-diagonal hopping terms
```

### File Format Support

#### `_hr.dat` Format
```
Header line
n_wann
n_R
degeneracy_weights (15 per line)
R_x R_y R_z m n H_real H_imag
...
```

#### `_tb.dat` Format
```
n_wann
n_R
Wannier_centers (3*n_wann values)
R_vectors and Hamiltonian_elements
...
```

## Usage Examples

### Basic Workflow

```python
from HamiltonIO.wannier import WannierHam
from ase.io import read

# 1. Read atomic structure
atoms = read("POSCAR")  # or from .win file

# 2. Load Wannier model
model = WannierHam.read_from_wannier_dir(
    path="./wannier90_calc/",
    prefix="wannier90",
    atoms=atoms
)

# 3. Calculate band structure
k_path = [[0,0,0], [0.5,0,0], [0.5,0.5,0], [0,0,0]]
bands = []

for k in k_path:
    evals, evecs = model.solve(k)
    bands.append(evals)

bands = np.array(bands)
```

### k-Space Analysis

```python
# Generate k-point mesh
from ase.dft.kpoints import monkhorst_pack

kpts = monkhorst_pack([8, 8, 8])
all_evals, all_evecs = model.solve_all(kpts)

print(f"Eigenvalues shape: {all_evals.shape}")  # (nk, nbasis)
print(f"Eigenvectors shape: {all_evecs.shape}")  # (nk, nbasis, nbasis)
```

### Real-Space Analysis

```python
# Examine hopping parameters
R_neighbors = [(1,0,0), (0,1,0), (0,0,1), (-1,0,0), (0,-1,0), (0,0,-1)]

for R in R_neighbors:
    if R in model.data:
        H_R = model.get_hamR(R)
        max_hopping = np.max(np.abs(H_R))
        print(f"R={R}: max hopping = {max_hopping:.3f} eV")
```

### File I/O Operations

```python
# Save model to NetCDF
model.save("wannier_model.nc")

# Load model from NetCDF
loaded_model = WannierHam.load_MyTB("wannier_model.nc")

# Convert from other formats
from tbmodels import Model
tbmodel = Model.from_hdf5_file("model.hdf5")
wannier_model = WannierHam.from_tbmodel(tbmodel)
```

## Integration with HamiltonIO

The Wannier interface integrates seamlessly with the broader HamiltonIO ecosystem:

```python
# Convert to LCAO format
from HamiltonIO.lcao_hamiltonian import LCAOHamiltonian

# Extract data for LCAO format
HR_data = np.array([model.data[R] for R in model.Rlist])
SR_data = None  # Wannier functions are orthogonal
Rlist = np.array(model.Rlist)

lcao_model = LCAOHamiltonian(
    HR=HR_data,
    SR=SR_data, 
    Rlist=Rlist,
    nbasis=model.nbasis,
    atoms=model.atoms,
    orth=True  # Orthogonal basis
)
```

## Performance Notes

- **Memory efficiency**: Use sparse matrices for large systems
- **k-point calculations**: Vectorized operations for multiple k-points
- **R-space storage**: Dictionary format allows efficient sparse storage
- **File formats**: `_tb.dat` includes centers, `_hr.dat` is more compact

## Common Issues

- **Gauge consistency**: Ensure consistent gauge between Wannier90 and analysis
- **Spin ordering**: Check spin ordering for SOC calculations (`groupby` parameter)
- **Position wrapping**: Use `shift_position()` to move functions near atoms
- **Unit conventions**: Energies in eV, positions in reduced coordinates
- **R-vector symmetry**: Hamiltonian uses time-reversal: H(R) = H†(-R)