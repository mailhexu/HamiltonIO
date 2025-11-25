# Distance-Resolved Matrix Element Analysis

This tutorial demonstrates how to analyze Wannier hopping matrix elements and electron-phonon coupling matrix elements as a function of distance in real space.

## Overview

Distance-resolved analysis is essential for understanding:
- **Localization** of electronic states (hopping decay)
- **Electron-phonon coupling locality** (coupling range)
- **Convergence** of Wannier representations
- **Physical interpretation** of tight-binding and electron-phonon models

## Wannier90: Hopping vs Distance

### Basic Usage

```python
from HamiltonIO.wannier import WannierHam
from ase.io import read
import matplotlib.pyplot as plt

# Load Wannier Hamiltonian
atoms = read("structure.cif")
ham = WannierHam.read_from_wannier_dir("./wannier90_calc", "wannier90", atoms=atoms)

# Get distance-resolved hoppings
cell = atoms.get_cell()
entries = ham.distance_resolved_hoppings(cell)

# Extract distances and magnitudes
distances = [e["distance"] for e in entries]
magnitudes = [abs(e["hopping"]) for e in entries]

# Plot scatter
plt.scatter(distances, magnitudes, s=3, alpha=0.6)
plt.xlabel("Distance (Å)")
plt.ylabel("|t| (eV)")
plt.yscale("log")
plt.title("Hopping Matrix Elements vs Distance")
plt.savefig("hoppings_vs_distance.pdf")
```

### Understanding the Output

Each entry in `distance_resolved_hoppings()` contains:
- `R`: Lattice vector (tuple of 3 integers)
- `i`, `j`: Wannier function indices
- `distance`: Real-space distance in Angstroms (scalar)
- `hopping`: Complex hopping matrix element t(R, i, j)

The distance is computed as:
```
d = |r_j + R - r_i|
```
where `r_i` and `r_j` are Wannier function centers in fractional coordinates, and `R` is the lattice vector.

### Binned Average

For cleaner visualization, bin the hoppings by distance:

```python
# Bin hoppings with 0.1 Å resolution
bin_centers, avg_magnitudes = WannierHam.bin_hoppings_by_distance(entries, dr=0.1)

plt.plot(bin_centers, avg_magnitudes, marker="o", linestyle="-")
plt.xlabel("Distance (Å)")
plt.ylabel("Average |t| (eV)")
plt.yscale("log")
plt.title("Binned Hopping vs Distance")
plt.savefig("hoppings_binned.pdf")
```

### Complete Example

See `examples/wannier/hoppings_vs_distance_SrMnO3.py` for a full working example:

```python
"""Plot hopping matrix elements vs distance for SrMnO3 Wannier TB model."""

import matplotlib.pyplot as plt
from ase.io import read
from HamiltonIO.wannier import WannierHam

# Load structure and Wannier model
atoms = read("SrMnO3.pwi", format="espresso-in")
ham = WannierHam.read_from_wannier_dir("./", "SrMnO3", atoms=atoms)

# Get distance-resolved hoppings
cell = atoms.get_cell()
entries = ham.distance_resolved_hoppings(cell)

# Extract data
distances = [e["distance"] for e in entries]
mags = [abs(e["hopping"]) for e in entries]

# Create 2-panel plot
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))

# Panel 1: Scatter plot
ax1.scatter(distances, mags, s=3, alpha=0.6)
ax1.set_xlabel("Distance (Å)")
ax1.set_ylabel("|t| (eV)")
ax1.set_yscale("log")
ax1.set_title("All Hopping Elements")

# Panel 2: Binned average
centers, avg_mag = WannierHam.bin_hoppings_by_distance(entries, dr=0.1)
ax2.plot(centers, avg_mag, marker="o", linestyle="-", linewidth=1, markersize=3)
ax2.set_xlabel("Distance (Å)")
ax2.set_ylabel("Average |t| (eV)")
ax2.set_yscale("log")
ax2.set_title("Binned Average")

plt.tight_layout()
plt.savefig("hoppings_vs_distance.pdf")
print(f"Found {len(entries)} hopping elements")
print(f"Distance range: {min(distances):.2f} - {max(distances):.2f} Å")
```

### Using the CLI Tool

HamiltonIO provides a command-line tool for quick Wannier distance analysis:

```bash
# Analyze hopping vs distance (auto-detect structure file)
hamiltonio-wannier distance --path ./wannier_calc --prefix wannier90

# Specify structure file explicitly
hamiltonio-wannier distance -p ./ -n wannier90 -s POSCAR -f vasp

# Custom output and y-axis limits
hamiltonio-wannier distance -p ./ -n material -o hopping.pdf --ylim 1e-4 10
```

**CLI Options:**
- `-p, --path`: Path to Wannier90 directory (default: current directory)
- `-n, --prefix`: Wannier90 file prefix (default: wannier90)
- `-s, --structure-file`: Structure filename (auto-detected if omitted)
- `-f, --structure-format`: Structure format (vasp, cif, espresso-in, etc.)
- `-o, --output`: Output plot filename (default: wannier_distance.pdf)
- `--ylim MIN MAX`: Y-axis limits (default: 1e-5 10)

**Example workflow:**
```bash
cd /path/to/wannier_calculation
hamiltonio-wannier distance -n mymat -o hoppings.pdf --ylim 1e-3 5
```

### Using the Plotting API

For more control and integration into existing scripts:

```python
import matplotlib.pyplot as plt
from HamiltonIO.wannier import plot_wannier_distance

# Create figure
fig, ax = plt.subplots(figsize=(6, 5))

# Plot using the API
plot_wannier_distance(
    ax=ax,
    path="./wannier_calc",
    prefix="wannier90",
    structure_file="POSCAR",
    structure_format="vasp",
    ylim=(1e-4, 10),
    color="steelblue",  # Pass matplotlib scatter kwargs
    s=3,
    alpha=0.6,
)

ax.set_title("Wannier Hopping vs Distance")
plt.tight_layout()
plt.savefig("plot.pdf")
```

See `examples/wannier/plot_api_example.py` for a complete working example.

## EPW: Electron-Phonon Coupling vs Distance

### Two Distance Metrics

EPW electron-phonon coupling matrices `g(Rg, Rk, mode, i, j)` involve **two distance metrics**:

1. **Rk-based distance** (electron WF to WF): Distance between electron Wannier functions
   ```
   d_Rk = |w_j + Rk - w_i|
   ```
   - Analogous to hopping distance in tight-binding models
   - Measures electron delocalization

2. **Rg-based distance** (electron WF to atom): Distance between electron WF and atomic displacement
   ```
   d_Rg = |τ_atom + Rg - w_i|
   ```
   - Measures locality of electron-phonon interaction
   - `τ_atom` is position of atom associated with the phonon mode

### Basic Usage: Rk-based

```python
from HamiltonIO.epw.epwparser import Epmat, read_crystal_fmt
from ase.units import Bohr
import matplotlib.pyplot as plt

# Load EPW data
crystal = read_crystal_fmt("crystal.fmt")
ep = Epmat()
ep.crystal = crystal
ep.read(path="./", prefix="material", epmat_ncfile="epmat.nc")

# Get cell in Cartesian coordinates (Angstrom)
cell = crystal.at.reshape(3, 3) * crystal.alat * Bohr

# Get distance-resolved couplings (Rk-based: WF-to-WF)
imode = 0  # First phonon mode
entries_Rk = ep.distance_resolved_couplings_Rk(imode=imode, cell=cell)

# Extract data
distances = [e["distance"] for e in entries_Rk]
magnitudes = [abs(e["coupling"]) for e in entries_Rk]

# Plot
plt.scatter(distances, magnitudes, s=3, alpha=0.6)
plt.xlabel("Distance (Å)")
plt.ylabel("|g| (eV/Å)")
plt.yscale("log")
plt.title(f"Electron-Phonon Coupling vs WF Distance - Mode {imode}")
plt.savefig("coupling_Rk_vs_distance.pdf")
```

### Basic Usage: Rg-based

```python
# Get distance-resolved couplings (Rg-based: WF-to-atom)
entries_Rg = ep.distance_resolved_couplings_Rg(imode=imode, cell=cell)

# Extract data
distances = [e["distance"] for e in entries_Rg]
magnitudes = [abs(e["coupling"]) for e in entries_Rg]

# Plot
plt.scatter(distances, magnitudes, s=3, alpha=0.6, color="orange")
plt.xlabel("Distance (Å)")
plt.ylabel("|g| (eV/Å)")
plt.yscale("log")
plt.title(f"Electron-Phonon Coupling vs Atom Distance - Mode {imode}")
plt.savefig("coupling_Rg_vs_distance.pdf")
```

### Understanding EPW Output

**Rk-based entries** (`distance_resolved_couplings_Rk`):
- `imode`: Phonon mode index
- `Rg`, `Rk`: Lattice vectors for phonon and electron
- `i`, `j`: Wannier function indices
- `distance`: Distance between WF i and WF j+Rk
- `coupling`: Complex coupling matrix element g

**Rg-based entries** (`distance_resolved_couplings_Rg`):
- `imode`: Phonon mode index
- `atom_index`: Index of atom associated with this mode (imode // 3)
- `Rg`, `Rk`: Lattice vectors
- `i`, `j`: Wannier function indices  
- `distance`: Distance between WF i and atom+Rg
- `coupling`: Complex coupling matrix element g

### Binned Analysis

```python
# Bin Rk-based couplings
centers_Rk, avg_Rk = Epmat.bin_couplings_by_distance(entries_Rk, dr=0.2)

# Bin Rg-based couplings
centers_Rg, avg_Rg = Epmat.bin_couplings_by_distance(entries_Rg, dr=0.2)

# Compare both metrics
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))

ax1.plot(centers_Rk, avg_Rk, marker="o", label="Rk (WF-WF)")
ax1.set_xlabel("Distance (Å)")
ax1.set_ylabel("Average |g| (eV/Å)")
ax1.set_title("Rk-based Distance")
ax1.legend()

ax2.plot(centers_Rg, avg_Rg, marker="o", color="orange", label="Rg (WF-atom)")
ax2.set_xlabel("Distance (Å)")
ax2.set_ylabel("Average |g| (eV/Å)")
ax2.set_title("Rg-based Distance")
ax2.legend()

plt.tight_layout()
plt.savefig("coupling_comparison.pdf")
```

### Complete Example

See `examples/epw/coupling_vs_distance_SrMnO3.py` for a full working example:

```python
"""Plot electron-phonon coupling vs distance for SrMnO3."""

import os
import matplotlib.pyplot as plt
import numpy as np
from HamiltonIO.epw.epwparser import Epmat, read_crystal_fmt
from ase.units import Bohr

# Load crystal structure
data_dir = "./epw_output"
crystal = read_crystal_fmt(os.path.join(data_dir, "crystal.fmt"))

# Read EPW data
ep = Epmat()
ep.crystal = crystal
ep.read(path=data_dir, prefix="SrMnO3", epmat_ncfile="epmat.nc")

# Get cell (convert to Cartesian Angstrom)
cell = crystal.at.reshape(3, 3) * crystal.alat * Bohr

# Analyze mode 0
imode = 0
print(f"Analyzing mode {imode}")
print(f"Crystal has {crystal.natom} atoms, {crystal.nmode} modes")
print(f"Epmat has {ep.nwann} Wannier functions")

# Get both distance metrics
entries_Rk = ep.distance_resolved_couplings_Rk(imode=imode, cell=cell)
entries_Rg = ep.distance_resolved_couplings_Rg(imode=imode, cell=cell)

# Extract data
distances_Rk = [e["distance"] for e in entries_Rk]
mags_Rk = [abs(e["coupling"]) for e in entries_Rk]

distances_Rg = [e["distance"] for e in entries_Rg]
mags_Rg = [abs(e["coupling"]) for e in entries_Rg]

# Create 4-panel plot
plt.figure(figsize=(12, 10))

# Panel 1: Rk scatter
plt.subplot(2, 2, 1)
plt.scatter(distances_Rk, mags_Rk, s=3, alpha=0.6)
plt.xlabel("Distance (Å)")
plt.ylabel("|g| (eV/Å)")
plt.yscale("log")
plt.title(f"Rk-based (WF-WF) - Mode {imode}")

# Panel 2: Rg scatter
plt.subplot(2, 2, 2)
plt.scatter(distances_Rg, mags_Rg, s=3, alpha=0.6, color="orange")
plt.xlabel("Distance (Å)")
plt.ylabel("|g| (eV/Å)")
plt.yscale("log")
plt.title(f"Rg-based (WF-atom) - Mode {imode}")

# Panel 3: Rk binned
centers_Rk, avg_Rk = Epmat.bin_couplings_by_distance(entries_Rk, dr=0.2)
plt.subplot(2, 2, 3)
plt.plot(centers_Rk, avg_Rk, marker="o", linestyle="-", linewidth=1, markersize=3)
plt.xlabel("Distance (Å)")
plt.ylabel("Average |g| (eV/Å)")
plt.title(f"Binned Rk-based - Mode {imode}")

# Panel 4: Rg binned
centers_Rg, avg_Rg = Epmat.bin_couplings_by_distance(entries_Rg, dr=0.2)
plt.subplot(2, 2, 4)
plt.plot(centers_Rg, avg_Rg, marker="o", linestyle="-", linewidth=1, 
         markersize=3, color="orange")
plt.xlabel("Distance (Å)")
plt.ylabel("Average |g| (eV/Å)")
plt.title(f"Binned Rg-based - Mode {imode}")

plt.tight_layout()
plt.savefig("coupling_vs_distance.pdf")
print(f"\nPlot saved to coupling_vs_distance.pdf")
```

### Using the CLI Tool

HamiltonIO provides a command-line tool for quick distance analysis without writing Python code:

```bash
# Convert EPW data to NetCDF first (if not already done)
hamiltonio-epw epw_to_nc --path ./epw_data --prefix material

# Analyze Rk-based (WF-WF) distance for mode 0
hamiltonio-epw distance --path ./epw_data --imode 0 --distance-type Rk -o mode0_Rk.pdf

# Analyze Rg-based (WF-atom) distance for mode 5
hamiltonio-epw distance -p ./epw_data -m 5 -t Rg -o mode5_Rg.pdf

# Custom y-axis limits
hamiltonio-epw distance -p ./ -m 0 --ylim 1e-4 10 -o custom_ylim.pdf
```

**CLI Options:**
- `-p, --path`: Directory containing EPW files (default: current directory)
- `-m, --imode`: Phonon mode index (required)
- `-t, --distance-type`: Either 'Rk' (WF-WF) or 'Rg' (WF-atom), default: Rk
- `-o, --output`: Output plot filename (default: epw_distance.pdf)
- `--epmat-ncfile`: NetCDF filename (default: epmat.nc)
- `--crystal-fmt-file`: Crystal format file (default: crystal.fmt)
- `--ylim MIN MAX`: Y-axis limits (default: 1e-5 10)

**Example workflow:**
```bash
# 1. Convert binary EPW data to NetCDF
cd /path/to/epw_calculation
hamiltonio-epw epw_to_nc --prefix mymat --output epmat.nc

# 2. Analyze first 3 modes with both distance types
hamiltonio-epw distance -m 0 -t Rk -o mode0_Rk.pdf
hamiltonio-epw distance -m 0 -t Rg -o mode0_Rg.pdf
hamiltonio-epw distance -m 1 -t Rk -o mode1_Rk.pdf
hamiltonio-epw distance -m 1 -t Rg -o mode1_Rg.pdf
hamiltonio-epw distance -m 2 -t Rk -o mode2_Rk.pdf
hamiltonio-epw distance -m 2 -t Rg -o mode2_Rg.pdf
```

### Using the Plotting API

For more control and integration into existing scripts, use the `plot_epw_distance` function:

```python
import matplotlib.pyplot as plt
from HamiltonIO.epw import plot_epw_distance

# Create figure
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))

# Plot Rk-based distances
plot_epw_distance(
    ax=ax1,
    path="./epw_data",
    imode=0,
    distance_type="Rk",
    ylim=(1e-4, 10),
    color="blue",  # Pass matplotlib scatter kwargs directly
    s=5,
    alpha=0.6,
)
ax1.set_title("Rk-based (WF-WF)")

# Plot Rg-based distances
plot_epw_distance(
    ax=ax2,
    path="./epw_data",
    imode=0,
    distance_type="Rg",
    ylim=(1e-4, 10),
    color="green",
)
ax2.set_title("Rg-based (WF-atom)")

plt.tight_layout()
plt.savefig("comparison.pdf")
```

See `examples/epw/plot_api_example.py` for a complete working example.

## Interpreting Results

### Wannier Hopping Decay

Expected behavior:
- **Short-range (~2-4 Å)**: Strong hoppings between nearest-neighbor orbitals
- **Medium-range (4-8 Å)**: Exponential decay
- **Long-range (>8 Å)**: Negligible hoppings (noise level)

**Troubleshooting:**
- If hoppings don't decay → Wannier functions are not well-localized
- If decay is too fast → May need more Wannier functions or better initial projections
- Large scatter at all distances → Check Wannier gauge consistency

### Electron-Phonon Coupling Decay

Expected behavior:

**Rk-based (WF-WF distance):**
- Tracks electron delocalization
- Similar decay to hopping matrix elements
- Useful for understanding electronic contribution to coupling

**Rg-based (WF-atom distance):**
- Tracks phonon locality
- Typically decays faster than Rk-based
- Directly related to physical electron-phonon interaction range

**Key insights:**
- If Rg-based coupling is short-ranged → Local electron-phonon interaction
- If Rk-based coupling is long-ranged → Delocalized electronic states contribute
- Compare different phonon modes to identify which modes couple most strongly

## Advanced Analysis

### Mode-by-Mode Analysis

```python
# Analyze multiple modes
modes_to_analyze = [0, 1, 2, 3, 4, 5]  # First 6 modes

fig, axes = plt.subplots(2, 3, figsize=(15, 10))

for idx, imode in enumerate(modes_to_analyze):
    ax = axes.flat[idx]
    
    # Get Rg-based distances (WF-to-atom)
    entries = ep.distance_resolved_couplings_Rg(imode=imode, cell=cell)
    centers, avg_mag = Epmat.bin_couplings_by_distance(entries, dr=0.2)
    
    ax.plot(centers, avg_mag, marker="o", linewidth=1, markersize=3)
    ax.set_xlabel("Distance (Å)")
    ax.set_ylabel("Average |g| (eV/Å)")
    ax.set_title(f"Mode {imode}")
    ax.set_yscale("log")

plt.tight_layout()
plt.savefig("coupling_by_mode.pdf")
```

### Atom-Resolved Analysis

```python
# Group modes by atom (each atom has 3 modes: x, y, z)
natoms = crystal.natom

for iatom in range(natoms):
    modes = [iatom * 3, iatom * 3 + 1, iatom * 3 + 2]  # x, y, z for this atom
    
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    
    for idx, imode in enumerate(modes):
        entries = ep.distance_resolved_couplings_Rg(imode=imode, cell=cell)
        centers, avg_mag = Epmat.bin_couplings_by_distance(entries, dr=0.2)
        
        axes[idx].plot(centers, avg_mag, marker="o")
        axes[idx].set_xlabel("Distance (Å)")
        axes[idx].set_ylabel("|g| (eV/Å)")
        axes[idx].set_title(f"Atom {iatom} - Mode {imode}")
        axes[idx].set_yscale("log")
    
    plt.tight_layout()
    plt.savefig(f"coupling_atom{iatom}.pdf")
```

### Convergence Tests

Check if your Wannier representation includes enough R vectors:

```python
# Get maximum distances
entries = ham.distance_resolved_hoppings(cell)
distances = np.array([e["distance"] for e in entries])
magnitudes = np.array([abs(e["hopping"]) for e in entries])

# Find distance where hopping magnitude drops below threshold
threshold = 1e-4  # eV
max_significant_distance = distances[magnitudes > threshold].max()

print(f"Maximum R vector distance: {distances.max():.2f} Å")
print(f"Distance where |t| > {threshold}: {max_significant_distance:.2f} Å")

if max_significant_distance > distances.max() * 0.8:
    print("WARNING: May need larger R-space grid (increase EPW/Wannier90 cutoff)")
else:
    print("R-space grid appears adequate")
```

## API Reference

### Wannier90 Methods

**`WannierHam.distance_resolved_hoppings(cell, use_absolute=True)`**
- `cell`: 3×3 lattice vectors in Angstrom
- `use_absolute`: If True, return scalar distances; if False, return 3D vectors
- Returns: List of dicts with keys `["R", "i", "j", "distance", "hopping"]`

**`WannierHam.bin_hoppings_by_distance(entries, dr=0.1)` (static method)**
- `entries`: Output from `distance_resolved_hoppings()`
- `dr`: Bin width in Angstroms
- Returns: `(bin_centers, average_magnitudes)` as numpy arrays

### EPW Methods

**`Epmat.distance_resolved_couplings_Rk(imode, cell, use_absolute=True)`**
- `imode`: Phonon mode index (0-based)
- `cell`: 3×3 lattice vectors in Angstrom
- `use_absolute`: If True, return scalar distances; if False, return 3D vectors
- Returns: List of dicts with WF-to-WF distances

**`Epmat.distance_resolved_couplings_Rg(imode, cell, use_absolute=True)`**
- `imode`: Phonon mode index (0-based)
- `cell`: 3×3 lattice vectors in Angstrom  
- `use_absolute`: If True, return scalar distances; if False, return 3D vectors
- Returns: List of dicts with WF-to-atom distances

**`Epmat.bin_couplings_by_distance(entries, dr=0.1)` (static method)**
- `entries`: Output from `distance_resolved_couplings_*`
- `dr`: Bin width in Angstroms
- Returns: `(bin_centers, average_magnitudes)` as numpy arrays

## Common Issues

**Q: Why are there multiple entries at the same distance?**
A: Multiple orbital pairs (i, j) or R vectors can give the same distance due to crystal symmetry.

**Q: My couplings/hoppings don't decay at large distances**
A: Check that:
- Wannier functions are well-localized (check spread in Wannier90 output)
- R-space grid is large enough
- Data is from a converged calculation

**Q: What's the difference between Rk and Rg distances for EPW?**
A: 
- **Rk distance**: Separation between electron Wannier functions (analogous to hopping)
- **Rg distance**: Separation between electron WF and atomic displacement (physical coupling range)

**Q: Should I use binned or scatter plots?**
A: 
- **Scatter**: Shows all data, useful for identifying outliers
- **Binned**: Cleaner, shows average trends, better for comparing systems

**Q: How do I choose the bin width `dr`?**
A: Start with 0.1-0.2 Å. Adjust based on:
- Data density (more data → smaller bins possible)
- Distance range of interest
- Crystal lattice constant (should resolve nearest-neighbor shells)

## See Also

- [Wannier90 Documentation](wannier.md) - Full Wannier90 interface reference
- [EPW Documentation](epw.md) - Complete EPW interface guide
- Wannier90 User Guide: http://wannier.org
- EPW Documentation: https://docs.epw-code.org/
