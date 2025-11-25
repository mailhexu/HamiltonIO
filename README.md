# HamiltonIO

A library for reading and analyzing Hamiltonian files from various DFT codes (ABACUS, SIESTA, EPW/Quantum ESPRESSO, Wannier90).

## Features

- **Multi-code support**: ABACUS, SIESTA, EPW, Wannier90
- **Spin configurations**: Non-polarized, collinear, non-collinear, spin-orbit coupling
- **Real and k-space analysis**: Band structures, DOS, Hamiltonian decomposition
- **Command-line tools**: Easy-to-use CLI for common analysis tasks
- **Python API**: Full programmatic access for custom workflows

## Installation

```bash
pip install hamiltonIO
```

Or install from source:

```bash
git clone https://github.com/mailhexu/HamiltonIO.git
cd HamiltonIO
pip install .
```

## Quick Start

### Command Line Tools

HamiltonIO provides CLI tools for quick analysis:

```bash
# Convert EPW files to NetCDF
hamiltonio-epw epw_to_nc --path ./data --prefix material --output epmat.nc

# Analyze ABACUS intra-atomic Hamiltonians
hamiltonio-abacus intra-atomic --outpath OUT.ABACUS

# Analyze SIESTA intra-atomic Hamiltonians  
hamiltonio-siesta intra-atomic siesta.fdf
```

For detailed CLI usage, see [CLI Documentation](docs/src/cli.md).

### Python Library

```python
# ABACUS
from HamiltonIO.abacus import AbacusParser
parser = AbacusParser(outpath="./OUT.material/")
model = parser.get_models()
evals, evecs = model.solve([0, 0, 0])  # Solve at Gamma point

# SIESTA
from HamiltonIO.siesta import SiestaHam
model = SiestaHam("siesta.fdf")
evals, evecs = model.solve([0, 0, 0])

# EPW
from HamiltonIO.epw import Epmat
epmat = Epmat()
epmat.read(path="./", prefix="material", epmat_ncfile="epmat.nc")

# Wannier90
from HamiltonIO.wannier import Wannier90Hamiltonian
model = Wannier90Hamiltonian("wannier90_hr.dat")
evals, evecs = model.solve([0, 0, 0])
```

## Documentation

Detailed documentation for each interface:

- [**Installation**](docs/src/install.md) - Installation methods and requirements
- [**CLI Tools**](docs/src/cli.md) - Command-line interface usage
- [**ABACUS**](docs/src/abacus.md) - ABACUS interface and examples
- [**SIESTA**](docs/src/siesta.md) - SIESTA interface and examples
- [**EPW**](docs/src/epw.md) - EPW/Quantum ESPRESSO interface
- [**Wannier90**](docs/src/wannier.md) - Wannier90 interface
- [**Distance Analysis**](docs/src/distance_analysis.md) - Hopping distance analysis

## Key Features by Code

### ABACUS
- Automatic spin detection (non-polarized, collinear, noncollinear)
- Binary and text CSR format support
- Split-SOC analysis via finite difference
- Intra-atomic Hamiltonian decomposition

### SIESTA
- NetCDF and HSX format support
- Split-SOC from native SIESTA output
- Collinear and non-collinear calculations
- Integration with sisl library

### EPW (Quantum ESPRESSO)
- Binary to NetCDF conversion
- Wigner-Seitz grid handling
- Electron-phonon coupling matrices
- Crystal structure from EPW data

### Wannier90
- HR format parsing
- Distance-resolved hopping analysis
- Band interpolation
- Orbital projections

## Common Workflows

### Intra-Atomic Analysis

```bash
# ABACUS: Analyze atoms 0, 1, 2 with full matrices
hamiltonio-abacus intra-atomic --outpath OUT.Fe \
    --atoms 0,1,2 --show-matrix -o analysis.txt

# SIESTA: Split-SOC analysis
hamiltonio-siesta intra-atomic siesta.fdf \
    --split-soc --show-matrix -o soc_analysis.txt
```

### File Conversion

```bash
# EPW: Convert binary to NetCDF (faster I/O)
hamiltonio-epw epw_to_nc --path ./epw_data \
    --prefix material --output epmat.nc
```

### Band Structure

```python
import numpy as np
from HamiltonIO.abacus import AbacusParser

# Load model
parser = AbacusParser(outpath="./OUT.material/")
model = parser.get_models()

# High-symmetry path
k_path = np.array([
    [0.0, 0.0, 0.0],  # Gamma
    [0.5, 0.0, 0.0],  # X
    [0.5, 0.5, 0.0],  # M
])

# Calculate bands
bands = []
for k in k_path:
    evals, _ = model.solve(k)
    bands.append(evals)
```

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is licensed under the terms specified in the [LICENSE](LICENSE) file.

## Citation

If you use HamiltonIO in your research, please cite:

```bibtex
@software{hamiltonIO,
  author = {Xu, He},
  title = {HamiltonIO: A Library for DFT Hamiltonian I/O},
  url = {https://github.com/mailhexu/HamiltonIO},
  year = {2024}
}
```

## Contact

For questions, bug reports, or feature requests:
- GitHub Issues: https://github.com/mailhexu/HamiltonIO/issues
- Email: mailhexu@gmail.com

## Acknowledgments

HamiltonIO builds upon and integrates with:
- [sisl](https://github.com/zerothi/sisl) - SIESTA interface
- [ASE](https://wiki.fysik.dtu.dk/ase/) - Atomic simulation environment
- [NumPy](https://numpy.org/) - Numerical computing
