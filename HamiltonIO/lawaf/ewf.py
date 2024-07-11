import numpy as np
from ase import Atoms
from scipy.linalg import eigh
from hamiltonIO.hamiltonian import Hamiltonian
from lawaf.mathutils.kR_convert import R_to_onek


class EWF(Hamiltonian):
    """
    LWF
    elements:
        Rlist: list of R-vectors. (nR, 3)
        wannR: Wannier functions in real space, (nR, nbasis, nwann)
        HwannR: total Hamiltonian in real space (nR, nwann, nwann)
        kpts: k-points
        kweights: weights of k-points
        wann_centers: centers of Wannier functions.
        wann_names: names of Wannier functions.
    """

    Rlist: np.ndarray = None
    Rdeg: np.ndarray = None
    wannR: np.ndarray = None
    HwannR: np.ndarray = None
    kpts: np.ndarray = None
    kweights: np.ndarray = None
    wann_centers: np.ndarray = None
    wann_names: list = None
    atoms: Atoms = None

    def __init__(
        self,
        Rlist,
        Rdeg,
        wannR,
        HwannR,
        kpts,
        kweights,
        wann_centers,
        wann_names,
        atoms=None,
    ):
        super().__init__(
            _name="LaWaF Electron Wannier",
        )

    def __post_init__(self):
        self.nR, self.nbasis, self.nwann = self.wannR.shape
        self.natoms = self.nbasis // 3
        self.nR = self.Rlist.shape[0]
        # self.check_normalization()

    def save_pickle(self, filename):
        """
        save the LWF to pickle file.
        """
        import pickle

        with open(filename, "wb") as f:
            pickle.dump(self, f)

    @classmethod
    def load_pickle(cls, filename):
        """
        load the LWF from pickle file.
        """
        import pickle

        with open(filename, "rb") as f:
            return pickle.load(f)

    def write_to_netcdf(self, filename):
        """
        write the LWF to netcdf file.
        """
        import xarray as xr

        print(f"wann_masses: {self.wann_masses}")
        ds = xr.Dataset(
            {
                "factor": self.factor,
                "Rlist": (["nR", "dim"], self.Rlist),
                "Rdeg": (["nR"], self.Rdeg),
                "wannR": (
                    ["ncplx", "nR", "nbasis", "nwann"],
                    np.stack([np.real(self.wannR), np.imag(self.wannR)], axis=0),
                ),
                "Hwann_R": (
                    ["ncplx", "nR", "nwann", "nwann"],
                    np.stack([np.real(self.HwannR), np.imag(self.HwannR)], axis=0),
                ),
                "kpts": (["nkpt", "dim"], self.kpts),
                "kweights": (["nkpt"], self.kweights),
                "wann_centers": (["nwann", "dim"], self.wann_centers),
                "wann_names": (["nwann"], self.wann_names),
            }
        )
        ds.to_netcdf(filename, group="wannier", mode="w")

        atoms = self.atoms
        if atoms is not None:
            ds2 = xr.Dataset(
                {
                    "positions": (["natom", "dim"], atoms.get_positions()),
                    "masses": (["natom"], atoms.get_masses()),
                    "cell": (["dim", "dim"], atoms.get_cell()),
                    "atomic_numbers": (["natom"], atoms.get_atomic_numbers()),
                }
            )
            ds2.to_netcdf(filename, group="atoms", mode="a")

    @classmethod
    def load_from_netcdf(cls, filename):
        """
        load the LWF from netcdf file.
        """
        import xarray as xr

        ds = xr.open_dataset(filename, group="wannier")
        wannR = ds["wannR"].values[0] + 1j * ds["wannR"].values[1]
        HwannR = ds["Hwann_R"].values[0] + 1j * ds["Hwann_R"].values[1]

        ds_atoms = xr.open_dataset(filename, group="atoms")
        atoms = Atoms(
            positions=ds_atoms["positions"].values,
            masses=ds_atoms["masses"].values,
            cell=ds_atoms["cell"].values,
            atomic_numbers=ds_atoms["atomic_numbers"].values,
        )

        return cls(
            Rlist=ds["Rlist"].values,
            Rdeg=ds["Rdeg"].values,
            wannR=wannR,
            HwannR=HwannR,
            kpts=ds["kpts"].values,
            kweights=ds["kweights"].values,
            wann_centers=ds["wann_centers"].values,
            wann_names=ds["wann_names"].values,
        )

    def remove_phase(self, Hk, k):
        """
        remove the phase of the R-vector
        """
        self.dr = self.wann_centers[None, :, :] - self.wann_centers[:, None, :]
        phase = np.exp(-2.0j * np.pi * np.einsum("ijk, k->ij", self.dr, k))
        return Hk * phase

    def check_normalization(self):
        """
        check the normalization of the LWF.
        """
        self.wann_norm = np.sum(self.wannR * self.wannR.conj(), axis=(0, 1)).real
        print(f"Norm of Wannier functions: {self.wann_norm}")

    def get_Hk(self, kpt):
        """
        get the Hamiltonian at k-point.
        """
        Hk = R_to_onek(kpt, self.Rlist, self.HwannR)
        return Hk

    def solve_k(self, kpt):
        """
        solve the Hamiltonian at k-point with NAC.
        """
        # if np.linalg.norm(kpt) < 1e-6:
        #    Hk = self.get_Hk_noNAC(kpt)
        # else:
        Hk = self.get_Hk(kpt)
        evals, evecs = eigh(Hk)
        return evals, evecs

    def solve_all(self, kpts):
        """
        solve the Hamiltonian at all k-points with NAC.
        """
        evals = []
        evecs = []
        for k in kpts:
            e, v = self.solve_k(k)
            evals.append(e)
            evecs.append(v)
        return np.array(evals), np.array(evecs)

    def HS_and_eigen(self, kpt):
        Hk = self.get_Hk(kpt)
        evals, evecs = eigh(Hk)
        return Hk, None, evals, evecs
