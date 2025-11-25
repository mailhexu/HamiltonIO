# HamiltonIO
A library of IO of Hamiltonian files of DFT codes

## Installation
```bash
pip install hamiltonIO
```

## Usage

### Command Line Tools

HamiltonIO includes command-line tools for data file conversion and analysis:

**EPW File Converter:**
```bash
# Convert EPW binary files to NetCDF format
hamiltonio-epw epw_to_nc --path /path/to/epw/files --prefix material --output epmat.nc

# Check files without converting
hamiltonio-epw epw_to_nc --dry-run --path ./data --prefix test
```

**ABACUS Hamiltonian Analysis:**
```bash
# Analyze intra-atomic (on-site) Hamiltonians
hamiltonio-abacus intra-atomic --outpath OUT.ABACUS

# Analyze specific atoms with full matrix output
hamiltonio-abacus intra-atomic --outpath OUT.ABACUS --atoms 0,1,2 --show-matrix

# Analyze split-SOC Hamiltonian (H = H_nosoc + H_soc) using finite difference
# Requires two calculations: one with SOC strength = 0 and one with SOC strength = 1
hamiltonio-abacus intra-atomic --outpath-nosoc soc0/OUT.ABACUS --outpath-soc soc1/OUT.ABACUS --show-matrix

# Save analysis to file
hamiltonio-abacus intra-atomic --outpath OUT.ABACUS -o analysis.txt
```

**SIESTA Hamiltonian Analysis:**
```bash
# Analyze intra-atomic (on-site) Hamiltonians
hamiltonio-siesta intra-atomic siesta.fdf

# Analyze specific atoms with full matrix output
hamiltonio-siesta intra-atomic siesta.fdf --atoms 0,1,2 --show-matrix

# Analyze split-SOC Hamiltonian (H = H_nosoc + H_soc)
hamiltonio-siesta intra-atomic siesta.fdf --split-soc --show-matrix

# Save analysis to file
hamiltonio-siesta intra-atomic siesta.fdf -o analysis.txt

# Analyze specific spin channel (for collinear calculations)
hamiltonio-siesta intra-atomic siesta.fdf --ispin up
```

### Python Library

See the examples in the documentation for detailed Python API usage.

