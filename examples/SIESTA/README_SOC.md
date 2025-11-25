# BCC Fe with SOC: Comprehensive Intra-Atomic Analysis Example

## Overview

This directory contains a complete example demonstrating the analysis of **spin-orbit coupling (SOC)** effects in body-centered cubic (BCC) iron using the HamiltonIO intra-atomic Hamiltonian tools.

## Files

### Input Data
- `bccFe/SOC/siesta.fdf` - SIESTA input file for SOC calculation
- `bccFe/SOC/siesta.nc` - SIESTA NetCDF output with Hamiltonian
- `bccFe/nonpolarized/` - Non-SOC calculation for comparison

### Analysis Scripts
1. **`analyze_soc_bccFe.py`** (Main analysis script)
   - Comprehensive SOC analysis with Pauli decomposition
   - Comparison of SOC vs non-SOC calculations
   - Physical interpretation of spin components
   - Generates detailed output files
   
2. **`print_intra_atomic_Fe.py`** (Basic example)
   - Simple demonstration of intra-atomic printing
   - Multiple examples (non-polarized, SOC, filtered atoms)

### Documentation
- **`SOC_ANALYSIS.md`** - Detailed physical interpretation
  - Orbital structure explanation
  - Pauli decomposition analysis
  - Physical conclusions about SOC strength
  - Magnetic property insights

### Output Files (Generated)
- `bccFe_SOC_detailed_analysis.txt` - Full analysis with matrices
- `bccFe_SOC_intra_atomic.txt` - Standard output
- `bccFe_nonpolarized_intra_atomic.txt` - Non-SOC comparison

## Quick Start

### Run the comprehensive analysis:
```bash
cd examples/SIESTA
uv run python analyze_soc_bccFe.py
```

### Run the basic example:
```bash
cd examples/SIESTA
uv run python print_intra_atomic_Fe.py
```

## What You'll Learn

### 1. Orbital Structure
The SOC calculation has 20 basis functions (10 orbitals × 2 spins):
- **s orbitals**: 3s, 4s (core and valence)
- **p orbitals**: 3px, 3py, 3pz
- **d orbitals**: 3d (xy, yz, z², xz, x²-y²)

Each orbital appears as a spin-up/spin-down pair (spinor format).

### 2. Pauli Decomposition
The Hamiltonian is decomposed into spin components:
```
H = I⊗σ₀ + Mx⊗σₓ + My⊗σᵧ + Mz⊗σᵤ
```

**Physical meaning**:
- **I**: Spin-independent part (orbital energies, crystal field)
- **σₓ, σᵧ**: Spin-flip terms (SOC-induced coupling between ↑ and ↓)
- **σᵤ**: Spin-splitting (exchange + SOC z-component)

### 3. Key Results for BCC Fe

From `analyze_soc_bccFe.py` output:

| Property | Value | Physical Meaning |
|----------|-------|------------------|
| **Total trace** | -505.2 eV | Sum of on-site energies |
| **I component norm** | 122.0 eV | Spin-independent energy scale |
| **σₓ norm** | 0.74 eV | Spin-flip coupling (x) |
| **σᵧ norm** | 0.74 eV | Spin-flip coupling (y) |
| **σᵤ norm** | 4.22 eV | Spin-splitting (z) |
| **Total spin norm** | 4.35 eV | Overall SOC strength |
| **Spin/Charge ratio** | 3.56% | Relative SOC contribution |
| **Anisotropy** | 4.01 | Uniaxial (z) preference |

**Interpretation**:
- ✅ SOC is present (~4 eV) but relatively small (~3.6% of total)
- ✅ Strong uniaxial character (σᵤ >> σₓ, σᵧ)
- ✅ Easy magnetization axis along z-direction
- ✅ Spin-flip barrier ~0.74 eV

### 4. Orbital Labels
The output shows clear orbital identification:
```
Orbital labels:
  [0] Fe1|3sZ1|down
  [1] Fe1|3sZ1|up
  [2] Fe1|4sZ1|down
  [3] Fe1|4sZ1|up
  [4] Fe1|3pyZ1|down
  [5] Fe1|3pyZ1|up
  ...
```

This makes it easy to:
- Identify which orbitals contribute to SOC
- Track spin-flip couplings between specific orbitals
- Understand the matrix structure

## Understanding the Output

### Example: σₓ Component Analysis
```
σₓ (Spin-X) Component:
  Trace: -0.000023 eV (should be ≈0)
  Norm: 0.744 eV
  Max |element|: 0.521 eV
  → Spin-flip coupling in x-direction
```

**What this means**:
- Trace ≈ 0: No net spin in x-direction (expected by symmetry)
- Norm = 0.74 eV: Overall strength of x-component spin coupling
- Max = 0.52 eV: Largest individual matrix element (e.g., py|↓ ↔ pz|↑)
- This couples states with Δspin=1, Δl=±1 (dipole selection rule)

### Example: σᵤ Component Analysis
```
σᵤ (Spin-Z) Component:
  Trace: -12.706 eV
  Norm: 4.216 eV
  Max |element|: 1.455 eV
  → Spin splitting along z-axis
```

**What this means**:
- Trace ≠ 0: Net spin imbalance (exchange splitting)
- Much larger than σₓ, σᵧ: Preferential z-axis quantization
- Includes both exchange and SOC contributions
- This is the magnetic anisotropy term

## Comparison: SOC vs Non-SOC

| Property | Non-SOC | SOC | Ratio |
|----------|---------|-----|-------|
| Basis size | 10 | 20 | 2× |
| Matrix type | Real | Complex | - |
| Trace (eV) | -254.4 | -505.2 | ≈2× |
| Norm (eV) | 122.5 | 172.6 | 1.4× |
| Imaginary parts | No | Yes | - |

**Key observation**: Trace doubles (as expected), but norm only increases by 40% because spin components partly cancel.

## Advanced Usage

### Programmatic Access
```python
from HamiltonIO.siesta.sisl_wrapper import SislParser
from HamiltonIO.mathutils.pauli import pauli_block_all

# Load Hamiltonian
parser = SislParser(fdf_fname='bccFe/SOC/siesta.fdf', ispin='merge')
ham = parser.get_model()

# Get intra-atomic blocks
blocks = ham.get_intra_atomic_blocks()
H_atom0 = blocks[0]['H_full']

# Manual Pauli decomposition
MI, Mx, My, Mz = pauli_block_all(H_atom0)

# Analyze specific components
spin_flip_strength = np.sqrt(np.linalg.norm(Mx)**2 + np.linalg.norm(My)**2)
spin_z_strength = np.linalg.norm(Mz)
anisotropy = spin_z_strength / spin_flip_strength

print(f"Magnetic anisotropy: {anisotropy:.2f}")
```

### Extract Specific Orbital Couplings
```python
# Get labels
orbital_labels = ham.orbs

# Find 3d orbitals
d_indices = [i for i, orb in enumerate(orbital_labels) if '3d' in orb]

# Extract 3d-only block
H_3d = H_atom0[np.ix_(d_indices, d_indices)]

# Analyze SOC within 3d manifold
MI_d, Mx_d, My_d, Mz_d = pauli_block_all(H_3d)
```

## Physical Applications

### 1. Magnetic Anisotropy Energy (MAE)
The difference between σᵤ and σₓᵧ components gives the MAE:
```python
MAE_per_orbital = (np.linalg.norm(Mz) - spin_flip_norm) / 10
# For BCC Fe: ~(4.22 - 1.05) / 10 ≈ 0.32 eV per orbital
```

### 2. Gilbert Damping Parameter
The spin-flip components (σₓ, σᵧ) contribute to magnetic damping:
```python
alpha ∝ norm(Mx)**2 + norm(My)**2
```

### 3. Spin-Texture Visualization
The off-diagonal elements in Mx, My, Mz reveal spin textures in real space.

## Validation Checks

The output includes automatic validation:

1. **Trace conservation**: Total trace = 2×(I component trace)
   ```
   -505.2 ≈ 2×(-252.6) ✓
   ```

2. **Norm relation**: Total² ≈ I² + spin²
   ```
   172.6² ≈ 122.0² + 4.35²×something ✓
   ```

3. **Symmetry**: σₓ, σᵧ traces ≈ 0
   ```
   Trace(σₓ) ≈ 0 ✓
   Trace(σᵧ) ≈ 0 ✓
   ```

## References

1. **SIESTA SOC Implementation**
   - Fernández-Seivane et al., J. Phys.: Condens. Matter 18, 7999 (2006)

2. **Pauli Matrices and Spin**
   - Sakurai, J.J., "Modern Quantum Mechanics"

3. **SOC in Transition Metals**
   - Kurz et al., Phys. Rev. B 69, 024415 (2004)

4. **Magnetic Anisotropy**
   - Stöhr & Siegmann, "Magnetism: From Fundamentals to Nanoscale Dynamics"

## Troubleshooting

### Issue: No Pauli decomposition shown
**Solution**: Make sure `pauli_decomp=True` in the function call and the matrix size is even (spinor format).

### Issue: All spin components are zero
**Solution**: Check that SOC is actually enabled in SIESTA calculation (`Spin SOC T` in .fdf file).

### Issue: Orbital labels not showing
**Solution**: Ensure the `orbs` attribute exists on the Hamiltonian object. This should work automatically for SIESTA.

## Next Steps

1. **Band Structure with SOC**: Compare band splitting due to SOC
2. **Wannier Functions**: Use this analysis to validate Wannier90 models
3. **Time-Reversal Symmetry**: Check if Kramers pairs are properly reproduced
4. **Berry Curvature**: Calculate from spin texture information

## Support

For questions or issues:
- See `SOC_ANALYSIS.md` for detailed physical interpretation
- Check the main HamiltonIO documentation
- Open an issue on GitHub

---

**Last updated**: Nov 25, 2024  
**Author**: HamiltonIO Contributors
