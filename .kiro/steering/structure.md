# Project Structure

## Root Directory Layout
```
HamiltonIO/                 # Main package directory
├── __init__.py            # Package exports: Hamiltonian, SiestaHam, WannierHam
├── hamiltonian.py         # Abstract base classes (BaseHamiltonian, Hamiltonian)
├── lcao_hamiltonian.py    # LCAO implementation (LCAOHamiltonian)
├── orbital.py             # Orbital-related utilities
├── utils.py               # General utilities
└── [subpackages]/         # DFT code-specific implementations
```

## Core Subpackages

### DFT Code Interfaces
- **`siesta/`** - SIESTA DFT code integration
- **`wannier/`** - Wannier90 tight-binding interface  
- **`abacus/`** - ABACUS DFT code support
- **`gpaw/`** - GPAW integration
- **`lawaf/`** - LAWAF (electron wannier) support
- **`epw/`** - EPW (Electron-Phonon Wannier) tools

### Utilities and Models
- **`mathutils/`** - Mathematical utilities (k-R conversion, Pauli matrices, spin rotation)
- **`model/`** - Physical models (density matrix, Fermi level, occupations)
- **`format/`** - File format handlers (sparse matrices)
- **`output/`** - Output generation utilities
- **`wantibexos/`** - Specialized analysis tools

## Key Architecture Components

### Base Classes (hamiltonian.py)
- `BaseHamiltonian` - Dataclass with basic properties
- `Hamiltonian` - Abstract interface defining standard methods:
  - `get_HR()`, `get_Hk()` - Real/k-space Hamiltonian access
  - `get_SR()`, `get_Sk()` - Overlap matrix access  
  - `HS_and_eigen()` - Eigenvalue calculations

### Implementation Pattern
Each DFT code follows consistent structure:
```
code_name/
├── __init__.py           # Exports main class
├── [code]_wrapper.py     # Main implementation class
├── [code]_parser.py      # File parsing utilities
└── utils.py              # Code-specific utilities
```

## Configuration and Build
- **`pyproject.toml`** - Modern Python packaging configuration
- **`.ruff.toml`** - Linter/formatter settings
- **`.pre-commit-config.yaml`** - Git hooks configuration
- **`docs/`** - Sphinx documentation source
- **`examples/`** - Usage examples organized by DFT code
- **`tests/`** - Test suite (minimal structure currently)

## Naming Conventions
- **Classes**: PascalCase (e.g., `SiestaHam`, `WannierHam`)
- **Methods**: snake_case with descriptive prefixes:
  - `get_*` for data retrieval
  - `set_*` for configuration
  - `_build_*` for internal construction
- **Properties**: snake_case (e.g., `norb`, `nbasis`, `Rlist`)
- **Private attributes**: Leading underscore (e.g., `_name`, `_nspin`)

## Data Flow Pattern
1. **Input**: DFT code output files → Code-specific parser
2. **Processing**: Raw data → Hamiltonian/Overlap matrices in R-space
3. **Transformation**: R-space → k-space via Fourier transform
4. **Analysis**: Eigenvalue solving, band structure, DOS calculations