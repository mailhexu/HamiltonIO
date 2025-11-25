# Command Line Interface (CLI)

HamiltonIO provides command-line tools for analyzing and converting DFT calculation data. All CLI tools are automatically installed with the package.

## Installation

```bash
pip install hamiltonIO
```

After installation, the following commands will be available:
- `hamiltonio-epw` - EPW file converter and analysis
- `hamiltonio-abacus` - ABACUS Hamiltonian analysis
- `hamiltonio-siesta` - SIESTA Hamiltonian analysis

## General Usage

All CLI tools support:
- `--help` - Show detailed help and options
- `--version` - Display version information

## EPW CLI

### EPW to NetCDF Converter

Convert EPW binary files to NetCDF format for easier analysis and faster I/O.

**Basic Usage:**
```bash
hamiltonio-epw epw_to_nc --path /path/to/epw --prefix material --output epmat.nc
```

**Options:**
- `-p, --path DIR` - Directory containing EPW files (default: current directory)
- `-n, --prefix PREFIX` - EPW file prefix (default: epw)
- `-o, --output FILE` - Output NetCDF filename (default: epmat.nc)
- `--force` - Overwrite existing output file without prompting
- `--dry-run` - Check files without converting

**Required Files:**
The converter validates that these files exist:
- `epwdata.fmt` - Basic dimensions
- `wigner.fmt` - Wigner-Seitz vectors
- `{prefix}.epmatwp` - EPW matrix elements

**Examples:**

```bash
# Convert with custom prefix
hamiltonio-epw epw_to_nc --prefix graphene --output graphene_epmat.nc

# Check files before converting
hamiltonio-epw epw_to_nc --path ./data --prefix test --dry-run

# Batch conversion
for prefix in mat1 mat2 mat3; do
    hamiltonio-epw epw_to_nc --path ./data --prefix $prefix \
        --output "${prefix}_epmat.nc" --force
done
```

## ABACUS CLI

### Intra-Atomic Hamiltonian Analysis

Analyze intra-atomic (on-site, R=(0,0,0)) Hamiltonian blocks with Pauli decomposition and SOC splitting.

**Basic Usage:**
```bash
hamiltonio-abacus intra-atomic --outpath OUT.ABACUS
```

**Options:**
- `--outpath DIR` - ABACUS OUT.* directory (single calculation)
- `--outpath-nosoc DIR` - OUT.* directory with SOC=0 (for split-SOC)
- `--outpath-soc DIR` - OUT.* directory with SOC=1 (for split-SOC)
- `-o, --output FILE` - Save analysis to file (default: stdout)
- `--atoms LIST` - Comma-separated atom indices (e.g., "0,1,2")
- `--show-matrix` - Display full matrix elements
- `--no-pauli` - Disable Pauli decomposition
- `--spin TYPE` - Spin configuration: non-polarized, collinear, noncollinear
- `--binary` - Read binary CSR files

**Examples:**

```bash
# Basic analysis
hamiltonio-abacus intra-atomic --outpath OUT.Fe

# Analyze specific atoms with full output
hamiltonio-abacus intra-atomic --outpath OUT.Fe --atoms 0,1,2 --show-matrix

# Split-SOC analysis (requires two calculations)
hamiltonio-abacus intra-atomic \
    --outpath-nosoc soc0/OUT.ABACUS \
    --outpath-soc soc1/OUT.ABACUS \
    --show-matrix

# Save to file
hamiltonio-abacus intra-atomic --outpath OUT.Fe -o analysis.txt
```

**Split-SOC Analysis:**

To analyze H = H_nosoc + H_soc using finite difference:

1. Run ABACUS with SOC strength = 0
2. Run ABACUS with SOC strength = 1
3. Use `--outpath-nosoc` and `--outpath-soc` options

```bash
hamiltonio-abacus intra-atomic \
    --outpath-nosoc ./soc_strength_0/OUT.material \
    --outpath-soc ./soc_strength_1/OUT.material
```

**Output Format:**

The analysis provides:
- System information (basis size, spin configuration)
- Per-atom blocks with orbital indices
- Hamiltonian statistics (trace, norm, max element)
- Pauli decomposition (I, σₓ, σᵧ, σᵤ components)
- SOC decomposition (H_nosoc, H_soc if available)
- Verification of H_full = H_nosoc + H_soc

## SIESTA CLI

### Intra-Atomic Hamiltonian Analysis

Analyze intra-atomic (on-site, R=(0,0,0)) Hamiltonian blocks from SIESTA calculations.

**Basic Usage:**
```bash
hamiltonio-siesta intra-atomic siesta.fdf
```

**Options:**
- `FDF_FILE` - Path to SIESTA .fdf file (required)
- `-o, --output FILE` - Save analysis to file (default: stdout)
- `--atoms LIST` - Comma-separated atom indices (e.g., "0,1,2")
- `--show-matrix` - Display full matrix elements
- `--no-pauli` - Disable Pauli decomposition
- `--split-soc` - Read split-SOC data (H = H_nosoc + H_soc)
- `--ispin SPIN` - Spin channel: up, down, or merge (for collinear)

**Examples:**

```bash
# Basic analysis
hamiltonio-siesta intra-atomic siesta.fdf

# Analyze specific atoms with matrices
hamiltonio-siesta intra-atomic siesta.fdf --atoms 0,1,2 --show-matrix

# Split-SOC analysis
hamiltonio-siesta intra-atomic siesta.fdf --split-soc --show-matrix

# Save to file
hamiltonio-siesta intra-atomic siesta.fdf -o analysis.txt

# Analyze spin-up channel only (collinear)
hamiltonio-siesta intra-atomic siesta.fdf --ispin up
```

**Required Files:**
- `.fdf` file - SIESTA input file
- `.nc` file - NetCDF output (same basename as .fdf)

**Split-SOC Analysis:**

SIESTA can write separate SOC and non-SOC Hamiltonian components. Use `--split-soc` to analyze:
- H_nosoc: Scalar relativistic part
- H_soc: Spin-orbit coupling part
- H_full = H_nosoc + H_soc

**Output Format:**

Similar to ABACUS CLI:
- System information with orbital labels
- Per-atom Hamiltonian blocks
- Statistics and decompositions
- SOC verification if split-SOC is used

## Common Workflows

### 1. Quick Analysis

```bash
# ABACUS
hamiltonio-abacus intra-atomic --outpath OUT.material

# SIESTA  
hamiltonio-siesta intra-atomic calculation.fdf
```

### 2. Detailed Analysis with Output File

```bash
# ABACUS with full matrices
hamiltonio-abacus intra-atomic --outpath OUT.Fe \
    --show-matrix -o fe_analysis.txt

# SIESTA with specific atoms
hamiltonio-siesta intra-atomic siesta.fdf \
    --atoms 0,1 --show-matrix -o atom_analysis.txt
```

### 3. SOC Analysis

```bash
# ABACUS finite difference method
hamiltonio-abacus intra-atomic \
    --outpath-nosoc soc0/OUT.material \
    --outpath-soc soc1/OUT.material \
    -o soc_analysis.txt

# SIESTA split-SOC
hamiltonio-siesta intra-atomic siesta.fdf \
    --split-soc -o soc_analysis.txt
```

### 4. Batch Processing

```bash
# Process multiple ABACUS calculations
for dir in OUT.*; do
    name=$(basename $dir)
    hamiltonio-abacus intra-atomic --outpath $dir -o "${name}_analysis.txt"
done

# Process multiple SIESTA calculations
for fdf in *.fdf; do
    name=$(basename $fdf .fdf)
    hamiltonio-siesta intra-atomic $fdf -o "${name}_analysis.txt"
done
```

## Troubleshooting

### EPW CLI

**Error: Missing files**
- Ensure all required files (epwdata.fmt, wigner.fmt, *.epmatwp) are present
- Use `--dry-run` to check file availability

**Error: Cannot overwrite**
- Use `--force` to overwrite existing output
- Or remove the existing output file manually

### ABACUS CLI

**Error: OUT directory not found**
- Verify the path to OUT.* directory
- Check that ABACUS calculation has completed

**Error: Missing files**
- Required: `data-HR-sparse_SPIN*.csr`, `data-SR-sparse_SPIN0.csr`, `STRU`, `running_scf.log`
- Ensure ABACUS was run with proper output settings

**Error: Cannot detect spin**
- Manually specify with `--spin` option
- Check `running_scf.log` for nspin value

### SIESTA CLI

**Error: NetCDF file not found**
- Ensure .nc file exists with same basename as .fdf
- Check that SIESTA calculation has completed
- Verify SIESTA was compiled with NetCDF support

**Error: Cannot read split-SOC**
- Ensure SIESTA was run with SOC splitting enabled
- Check that both H_nosoc and H_soc are in .nc file
- Not all SIESTA versions support split-SOC output

## Tips

1. **Use dry-run for EPW**: Always check files before conversion
2. **Save output to files**: Use `-o` for large analyses
3. **Start without --show-matrix**: Full matrices can be very large
4. **Filter atoms**: Use `--atoms` to focus on specific atoms
5. **Batch processing**: Use shell loops for multiple calculations

## Getting Help

For detailed help on any command:

```bash
hamiltonio-epw --help
hamiltonio-epw epw_to_nc --help
hamiltonio-abacus intra-atomic --help
hamiltonio-siesta intra-atomic --help
```

For issues and bug reports, visit:
https://github.com/mailhexu/HamiltonIO/issues
