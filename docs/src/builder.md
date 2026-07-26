# Tight-Binding Builder

`HamiltonIO.builder.TightBindingBuilder` builds small orthogonal tight-binding models from structures, orbital assignments, onsite terms, and hopping terms. It returns a `TightBindingModelData` object containing the generated `LCAOHamiltonian` plus orbital metadata helpers.

```python
import numpy as np
from HamiltonIO.builder import TightBindingBuilder

builder = TightBindingBuilder.from_data(
    cell=np.eye(3) * 5.0,
    positions=[[0.0, 0.0, 0.0], [1.8, 0.0, 0.0]],
    symbols=["Fe", "O"],
)
builder.add_orbitals(atom=0, labels=["s"])
builder.add_orbitals(atom=1, labels=["px", "py", "pz"])
builder.set_onsite(species="Fe", orbital="s", value=0.0)
builder.set_onsite(species="O", shell=1, value=-3.0)
builder.add_hopping((0, "s"), (1, "px"), R=(0, 0, 0), value=-1.5)

model_data = builder.build(nspin=2, nel=4)
ham = model_data.hamiltonian
```

The spinless orbital order is atom-major: atom order from the structure, then the user-provided orbital label order for each atom. For `nspin=2`, HamiltonIO uses spin-interleaved order: `[orb0_up, orb0_dn, orb1_up, orb1_dn, ...]`.

Supported selectors are:

- `atom=<int>` for all orbitals on one atom.
- `species=<str>` for all orbitals on atoms with a chemical symbol.
- `orbital=<str>` for a label such as `s`, `px`, or `dxy`.
- `shell=<int>` for angular momentum `l`.
- `(atom_index, orbital_label)` for one atom-local label.
- `index=<int>` for one spinless orbital index.

Hopping terms use integer lattice translations `R`. By default, `add_hopping` also inserts the Hermitian counterpart, so a term `H_R(i,j)=t` creates `H_-R(j,i)=conj(t)`. Disable this only when explicitly assembling both directions yourself.

The first builder version is intentionally orthogonal. It sets the overlap to identity at `R=(0, 0, 0)` and zero elsewhere. Slater-Koster parameterization, sp-bond superexchange defaults, and Hubbard workflows belong in downstream examples such as TBUpy, not in the generic HamiltonIO builder.
