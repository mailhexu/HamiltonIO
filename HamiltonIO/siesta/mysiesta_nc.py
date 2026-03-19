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
        # H = self._read_class_spin(Hamiltonian, **kwargs)
        try:
            # H = #self._read_class_spin(Hamiltonian, **kwargs)
            H = self._r_class_spin(Hamiltonian, **kwargs)
        except AttributeError:
            H = self.read_hamiltonian()

        sp = self.groups["SPARSE"]
        if "H_so" in sp.variables:
            if sp.variables["H_so"].unit != "Ry":
                raise SileError(
                    f"{self}.read_soc_hamiltonian requires the stored matrix to be in Ry!"
                )
            # H_so has shape (n_spin_components, nnzs)
            # For SOC, there are 8 spin components
            n_spin = sp.variables["H_so"].shape[0]
            for i in range(n_spin):
                H._csr._D[:, i] = sp.variables["H_so"][i, :] * Ry / eV
        elif "ReH_so" in sp.variables and "ImH_so" in sp.variables:
            re = sp.variables["ReH_so"]
            im = sp.variables["ImH_so"]
            if re.unit != "Ry" or im.unit != "Ry":
                raise SileError(
                    f"{self}.read_soc_hamiltonian requires the stored matrix to be in Ry!"
                )
            # ReH_so and ImH_so are stored with Fortran shape (nnzs, spin_cmplx=4),
            # which in Python/NumPy (row-major) appears as (spin_cmplx=4, nnzs).
            # The 4 complex spin components follow CVEC4 ordering:
            #   cvec4(0) = H11,  cvec4(1) = H22,  cvec4(2) = H12,  cvec4(3) = H21
            #
            # Reconstruct VEC8 (8 real spin components) per m_spin_conventions.F90:
            #   D[:,0] = real(H11)   = re[0,:]
            #   D[:,1] = real(H22)   = re[1,:]
            #   D[:,2] = real(H12)   = re[2,:]
            #   D[:,3] = -imag(H12)  = -im[2,:]   ← minus sign per VEC8 convention
            #   D[:,4] = imag(H11)   = im[0,:]
            #   D[:,5] = imag(H22)   = im[1,:]
            #   D[:,6] = real(H21)   = re[3,:]
            #   D[:,7] = imag(H21)   = im[3,:]
            scale = Ry / eV
            H._csr._D[:, 0] = re[0, :] * scale  # real(H11)
            H._csr._D[:, 1] = re[1, :] * scale  # real(H22)
            H._csr._D[:, 2] = re[2, :] * scale  # real(H12)
            H._csr._D[:, 3] = -im[2, :] * scale  # -imag(H12)
            H._csr._D[:, 4] = im[0, :] * scale  # imag(H11)
            H._csr._D[:, 5] = im[1, :] * scale  # imag(H22)
            H._csr._D[:, 6] = re[3, :] * scale  # real(H21)
            H._csr._D[:, 7] = im[3, :] * scale  # imag(H21)
        else:
            raise SileError(
                f"{self}.read_soc_hamiltonian could not find H_so or ReH_so/ImH_so "
                f"variables in the SPARSE group!"
            )

        # fix siesta specific notation
        _mat_siesta2sisl(H)
        # H._csr._D[:, 3] *= -1
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
