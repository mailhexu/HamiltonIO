from dataclasses import dataclass
import numpy as np
from typing import Union
from ase.atoms import Atoms


@dataclass
class Orbital:
    name: str = "Unamed"
    n: int = 0
    l : int = 0
    m : int = 0
    atom : Atoms = None
    xred : np.ndarray = None
    spin: Union[int, None] = None

    def __str__(self):
        return f"{self.name}"

    @property
    def cell(self):
        return self.atom.get_cell()

    @property
    def xcart(self):
        return self.cell.get_cartesian_coords(self.xred)


def test_orbital():
    from ase.build import bulk
    from ase.io import read
    from ase.spacegroup import crystal

    atoms = bulk("Si")
    orbital = Orbital(name="1s", n=1, l=0, m=0, atom=atoms, xred=[0.0, 0.0, 0.0])
    print(orbital)
    print(orbital.cell)
    print(orbital.xcart)

    atoms = crystal("Si", [(0, 0, 0)], spacegroup=227, cellpar=[5.43, 5.43, 5.43, 90, 90, 90])
    orbital = Orbital(name="1s", n=1, l=0, m=0, atom=atoms, xred=[0.0, 0.0, 0.0])
    print(orbital)
    print(orbital.cell)
    print(orbital.xcart)

test_orbital()
