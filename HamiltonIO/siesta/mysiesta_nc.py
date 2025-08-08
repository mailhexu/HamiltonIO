try:
    from sisl.io.siesta._help import _mat_spin_convert as _mat_siesta2sisl  # noqa: F401
except ImportError:
    from sisl.io.siesta._help import _mat_siesta2sisl  # noqa: F401


from ase.units import Ry, eV
from sisl.io.siesta.siesta_nc import ncSileSiesta
from sisl.io.sile import SileError
from sisl.physics.hamiltonian import Hamiltonian


class MySiestaNC(ncSileSiesta):
    def read_soc_hamiltonian(self, **kwargs) -> Hamiltonian:
        """Returns a spin-orbit coupling Hamiltonian from the underlying NetCDF file"""
        #H = self._read_class_spin(Hamiltonian, **kwargs)
        try:
            #H = #self._read_class_spin(Hamiltonian, **kwargs)
            H = self._r_class_spin(Hamiltonian, **kwargs)
        except AttributeError:
            H = self.read_hamiltonian()

        sp = self.groups["SPARSE"]
        if sp.variables["H_so"].unit != "Ry":
            raise SileError(
                f"{self}.read_soc_hamiltonian requires the stored matrix to be in Ry!"
            )

        for i in range(len(H.spin)):
            H._csr._D[:, i] = sp.variables["H_so"][i, :] * Ry / eV

        # fix siesta specific notation
        _mat_siesta2sisl(H)
        #H._csr._D[:, 3] *= -1
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
    print(f"Total charge: {Qtot}")
    print(H)


if __name__ == "__main__":
    test_mysieta_nc()
