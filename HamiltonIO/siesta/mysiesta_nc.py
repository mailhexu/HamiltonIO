from sisl.io.siesta._help import _mat_spin_convert
from sisl.io.sile import SileError
from sisl.io.siesta.siesta_nc import ncSileSiesta
from sisl.physics.hamiltonian import Hamiltonian
from ase.units import Ry, eV


class MySiestaNC(ncSileSiesta):
    def read_soc_hamiltonian(self, **kwargs) -> Hamiltonian:
        """Returns a spin-orbit coupling Hamiltonian from the underlying NetCDF file"""
        try:
            H = self._read_class_spin(Hamiltonian, **kwargs)
        except AttributeError as E:
            H = self._r_class_spin(Hamiltonian, **kwargs)

        sp = self.groups["SPARSE"]
        if sp.variables["H_so"].unit != "Ry":
            raise SileError(
                f"{self}.read_soc_hamiltonian requires the stored matrix to be in Ry!"
            )

        for i in range(len(H.spin)):
            H._csr._D[:, i] = sp.variables["H_so"][i, :] * Ry / eV

        # fix siesta specific notation
        # _mat_spin_convert(H)
        H._csr._D[:, 3] *= -1
        # H._csr._D[:, 7] *= -1
        return H.transpose(spin=False, sort=kwargs.get("sort", True))

    def read_qtot(self):
        """Returns the total charge of the system"""
        return self._value("Qtot")[:][0]


def test_mysieta_nc():
    # Create a new instance of the MySiestaNC class
    sile = MySiestaNC(
        "/home/hexu/projects/TB2J_examples/Siesta/BiFeO3/BiFeO3_splitSOC/siesta.nc"
    )
    # Read the total charge of the system
    Qtot = sile.read_qtot()
    # Read the spin-orbit coupling Hamiltonian
    H = sile.read_soc_hamiltonian()
    print(Qtot)


if __name__ == "__main__":
    test_mysieta_nc()
