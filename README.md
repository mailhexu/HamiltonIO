# HamiltonIO
A library of IO of Hamiltonian files of DFT codes

## Installation
```bash
pip install hamiltonIO
```

## Usage

### Command Line Tools

HamiltonIO includes command-line tools for converting data files:

**EPW File Converter:**
```bash
# Convert EPW binary files to NetCDF format
hamiltonio-epw epw_to_nc --path /path/to/epw/files --prefix material --output epmat.nc

# Check files without converting
hamiltonio-epw epw_to_nc --dry-run --path ./data --prefix test
```

### Python Library

See the examples in the documentation for detailed Python API usage.

