from types import SimpleNamespace

import numpy as np
from ase.units import Ry, eV

from HamiltonIO.lcao_hamiltonian import LCAOHamiltonian
from HamiltonIO.siesta.mysiesta_nc import MySiestaNC


class _Variable:
    unit = "Ry"

    def __init__(self, values):
        self.values = np.asarray(values, dtype=float)
        self.shape = self.values.shape

    def __getitem__(self, key):
        return self.values[key]


class _FakeHamiltonian:
    def __init__(self, nnzs):
        self._csr = SimpleNamespace(_D=np.zeros((nnzs, 8), dtype=float))

    def transpose(self, *, spin, sort):
        assert spin is False
        assert sort is True
        return self


class _FakeSiestaNC:
    def __init__(self, real, imaginary):
        self.hamiltonian = _FakeHamiltonian(real.shape[1])
        self.groups = {
            "SPARSE": SimpleNamespace(
                variables={
                    "ReH_so": _Variable(real),
                    "ImH_so": _Variable(imaginary),
                }
            )
        }

    def _r_class_spin(self, _hamiltonian_class, **_kwargs):
        return self.hamiltonian


def test_lcao_solve_uses_generalized_hermitian_eigensolver():
    model = LCAOHamiltonian(
        HR=np.diag([2.0, 6.0])[None, ...],
        SR=np.diag([1.0, 2.0])[None, ...],
        Rlist=np.zeros((1, 3), dtype=int),
        nbasis=2,
    )

    eigenvalues, _ = model.solve([0.0, 0.0, 0.0])

    np.testing.assert_allclose(eigenvalues, [2.0, 3.0])


def test_mr487_conjugated_h21_component_uses_negative_imaginary_part():
    real = np.arange(8, dtype=float).reshape(4, 2)
    imaginary = np.arange(10, 18, dtype=float).reshape(4, 2)
    reader = _FakeSiestaNC(real, imaginary)

    result = MySiestaNC.read_soc_hamiltonian(reader)

    expected = np.column_stack(
        (
            real[0],
            real[1],
            real[2],
            -imaginary[2],
            imaginary[0],
            imaginary[1],
            real[3],
            -imaginary[3],
        )
    )
    np.testing.assert_allclose(result._csr._D, expected * Ry / eV)
