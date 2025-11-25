# BCC Fe with SOC: Intra-Atomic Hamiltonian Analysis

## System Overview

**Material**: Body-Centered Cubic (BCC) Iron  
**Calculation Type**: Spin-Orbit Coupling (SOC) included  
**Code**: SIESTA DFT calculation

### System Parameters
- **Atoms**: 1 Fe atom at origin
- **Total basis functions**: 20 (10 orbitals × 2 spin components)
- **Basis set**: DZP (Double-Zeta Polarized)
- **nspin**: 1 (spinor format, not collinear)
- **Number of R-vectors**: 343

## Orbital Structure

The 20 basis functions consist of 10 orbitals, each with spin-up and spin-down components:

### s orbitals (4 total)
- `[0-1]` Fe1|3sZ1 (down/up) - Core 3s orbital
- `[2-3]` Fe1|4sZ1 (down/up) - Valence 4s orbital

### p orbitals (6 total)
- `[4-5]` Fe1|3pyZ1 (down/up)
- `[6-7]` Fe1|3pzZ1 (down/up)
- `[8-9]` Fe1|3pxZ1 (down/up)

### d orbitals (10 total)
- `[10-11]` Fe1|3dxyZ1 (down/up)
- `[12-13]` Fe1|3dyzZ1 (down/up)
- `[14-15]` Fe1|3dz2Z1 (down/up)
- `[16-17]` Fe1|3dxzZ1 (down/up)
- `[18-19]` Fe1|3dx2-y2Z1 (down/up)

**Note**: The ordering shows orbitals are stored in spinor format (down/up pairs).

## Intra-Atomic Hamiltonian (R=(0,0,0))

### Overall Statistics
- **Matrix Shape**: 20×20 (complex)
- **Trace**: -505.19 eV (total on-site energy)
- **Frobenius Norm**: 172.58 eV
- **Max |element|**: 85.19 eV (diagonal element, likely core 3s orbital)
- **Complex**: Yes (imaginary parts present due to SOC)

### Physical Interpretation

#### 1. Energy Scale
The trace of -505 eV represents the sum of all on-site energies. Divided by 10 orbitals gives approximately **-50.5 eV per orbital** (average). This is reasonable for Fe with:
- Deep core 3s: ~-85 eV
- Valence 4s: ~-12 eV  
- 3p orbitals: ~-50 eV
- 3d orbitals: near Fermi level

#### 2. Spin-Orbit Coupling Effects

The **imaginary parts** in the Hamiltonian are a signature of SOC. These arise from:
```
H_SOC ∝ L·S
```
where L is orbital angular momentum and S is spin. This couples states with different spin and orbital angular momentum.

## Pauli Decomposition Analysis

The Hamiltonian is decomposed as:
```
H = I⊗σ₀ + Mx⊗σₓ + My⊗σᵧ + Mz⊗σᵤ
```

### Component Analysis

#### I (Identity/Charge) Component
- **Trace**: -252.60 eV (half of total, as expected)
- **Norm**: 121.96 eV
- **Max |element|**: 83.86 eV
- **Physical meaning**: Spin-independent part of the Hamiltonian
  - Contains orbital energies without spin effects
  - Crystal field splitting
  - Coulomb interactions

#### σₓ (Spin-X) Component  
- **Trace**: -0.000023 eV (≈ 0, as expected by symmetry)
- **Norm**: 0.744 eV
- **Max |element|**: 0.521 eV
- **Physical meaning**: Spin-flip coupling in x-direction
  - Connects |↑⟩ and |↓⟩ states
  - Induced by SOC: L×∇V · S has x-component
  - Small but non-zero due to SOC

#### σᵧ (Spin-Y) Component
- **Trace**: -0.000033 eV (≈ 0, as expected)
- **Norm**: 0.744 eV  
- **Max |element|**: 0.521 eV
- **Physical meaning**: Spin-flip coupling in y-direction
  - Complex/imaginary contributions
  - Equal magnitude to σₓ (isotropic SOC)
  - Small coupling (~0.5 eV max)

#### σᵤ (Spin-Z) Component
- **Trace**: -12.71 eV (non-zero, breaking spin symmetry)
- **Norm**: 4.22 eV
- **Max |element|**: 1.46 eV
- **Physical meaning**: Spin-dependent diagonal terms
  - Zeeman-like splitting between |↑⟩ and |↓⟩
  - Exchange splitting + SOC z-component
  - **Much larger than σₓ, σᵧ** → preferential z-axis quantization

### Spin Norm Analysis

**Total spin norm**: 4.35 eV  
**Spin/Charge ratio**: 0.0356 (3.56%)

This means:
- SOC effects are **~3.6% of the total Hamiltonian magnitude**
- Dominated by σᵤ component (4.22 eV out of 4.35 eV)
- Spin-flip terms (σₓ, σᵧ) are small: ~0.74 eV each
- **Conclusion**: System has strong uniaxial spin character (z-direction)

## Matrix Structure Insights

### Diagonal Blocks (Real Part)
Looking at the first 10×10 block (real part):

1. **3s orbitals** (rows/cols 0-3):
   - Large negative diagonal: -85.19 eV (3s|down), -82.53 eV (3s|up)
   - **Spin splitting**: ~2.7 eV between up/down
   - Small off-diagonal: 1.69 eV (3s-4s coupling)

2. **4s orbitals** (rows/cols 2-3):
   - Diagonal: -12.75 eV (down), -12.05 eV (up)
   - **Spin splitting**: ~0.7 eV

3. **3p orbitals** (rows/cols 4-9):
   - Diagonal: ~-52 eV (down), ~-49 eV (up)
   - **Spin splitting**: ~3 eV
   - Very small inter-orbital coupling

4. **3d orbitals** (rows/cols 10-19):
   - Near Fermi level (not shown in first 10×10)
   - Expected to show larger SOC effects

### Off-Diagonal Structure (Imaginary Part)

The imaginary parts show **spin-orbit coupling** between orbitals:

**Example**: Row 4-9 (3p orbitals) have imaginary elements ~±0.52 eV
- These couple py, pz, px orbitals with different spins
- Pattern: `iL×∇V` operator connects p orbitals with ΔL = ±1
- Example: `3py|down` ↔ `3pz|up` via Lₓ operator

## Physical Conclusions

### 1. Spin-Orbit Coupling Strength
- **SOC energy scale**: ~4 eV (from spin norm)
- **Relative strength**: 3.6% of total Hamiltonian
- **Type**: Predominantly uniaxial (z-axis preferred)

### 2. Magnetic Properties
- **Easy axis**: z-direction (σᵤ >> σₓ, σᵧ)
- **Spin-flip barrier**: ~0.74 eV (from σₓ, σᵧ)
- **Exchange splitting**: ~2-3 eV (from σᵤ trace and diagonal differences)

### 3. Orbital Character
- **s orbitals**: Minimal SOC (spherical symmetry)
- **p orbitals**: Moderate SOC (~0.5 eV imaginary parts)
- **d orbitals**: Expected to have strongest SOC (not fully shown in 10×10 block)

### 4. Validation
✅ **Trace conservation**: Total trace = 2×(I trace) = 2×(-252.6) ≈ -505 eV  
✅ **Norm relation**: Total norm² ≈ (I norm)² + (spin norm)²  
✅ **Symmetry**: σₓ, σᵧ have zero trace (no net spin in x,y)  
✅ **Complex structure**: Imaginary parts present only in off-diagonal (SOC)

## Comparison: SOC vs Non-SOC

| Property | Non-SOC (10×10) | SOC (20×20) |
|----------|----------------|-------------|
| Basis size | 10 | 20 (10×2 spins) |
| Matrix type | Real | Complex |
| Spin treatment | Separate up/down | Coupled spinor |
| SOC effects | None | σₓ, σᵧ, σᵤ ≠ 0 |
| Trace | -254.4 eV | -505.2 eV (≈2×) |
| Imaginary parts | No | Yes (SOC coupling) |

## Recommendations for Further Analysis

1. **Orbital-resolved DOS**: Project onto specific orbital characters (s, p, d)
2. **Spin texture**: Analyze σₓ, σᵧ, σᵤ spatial distribution
3. **Band structure**: Compare SOC vs non-SOC bands
4. **Magnetic anisotropy energy**: Calculate from σᵤ - (σₓ, σᵧ) difference
5. **Wannierization**: Use this as validation for Wannier90 tight-binding

## Usage Example

```python
from HamiltonIO.siesta.sisl_wrapper import SislParser
from pathlib import Path

# Load SOC data
parser = SislParser(fdf_fname='examples/SIESTA/bccFe/SOC/siesta.fdf', ispin='merge')
ham = parser.get_model()

# Analyze intra-atomic Hamiltonian
ham.print_intra_atomic_hamiltonian(
    atom_indices=None,       # All atoms
    output_file='soc_analysis.txt',
    pauli_decomp=True,       # Show Pauli decomposition
    show_matrix=True         # Show matrix elements
)

# Programmatic access
blocks = ham.get_intra_atomic_blocks()
for iatom, data in blocks.items():
    H_full = data['H_full']
    # Your analysis here
```

## References

- SIESTA Manual: https://departments.icmab.es/leem/siesta/
- Pauli matrices and SOC: Relativistic Quantum Mechanics textbooks
- SOC in DFT: Kurz et al., PRB 69, 024415 (2004)
