"""Minimal TightBindingBuilder example."""

import numpy as np

from HamiltonIO.builder import TightBindingBuilder

builder = TightBindingBuilder.from_data(
    cell=np.eye(3) * 5.0,
    positions=[[0.0, 0.0, 0.0], [1.8, 0.0, 0.0]],
    symbols=["Fe", "O"],
)
builder.add_orbitals(atom=0, labels=["s"])
builder.add_orbitals(atom=1, labels=["px"])
builder.set_onsite((0, "s"), value=0.0)
builder.set_onsite((1, "px"), value=-3.0)
builder.add_hopping((0, "s"), (1, "px"), R=(0, 0, 0), value=-1.5)

model_data = builder.build(nspin=1, nel=2)
print(model_data.hamiltonian.get_H0())
