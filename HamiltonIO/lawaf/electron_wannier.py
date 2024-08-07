import numpy as np
from ase import Atoms
from scipy.linalg import eigh
from HamiltonIO.hamiltonian import Hamiltonian
from HamiltonIO.mathutils.kR_convert import R_to_onek


class LawafHamiltonian(Hamiltonian):
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
    def __init__(
        self,
        Rlist=None,
        Rdeg=None,
        wannR=None,
        HwannR=None,
        SwannR=None,
        wann_centers=None,
        wann_names=None,
        atoms=None,
        kpts=None,
        kweights=None,
        is_orthogonal =True,
    ):
        if wannR is not None:
            norb = wannR.shape[2]
        else:
            norb = HwannR.shape[1]
        super().__init__(
            _name="LaWaF Electron Wannier",
            is_orthogonal=is_orthogonal, 
            R2kfactor=2 * np.pi,
            nspin=1,
            #norb=wannR.shape[2],
        )
        self.Rlist = Rlist
        self.Rdeg = Rdeg
        self.wannR = wannR
        self.HwannR = HwannR
        self.SwannR = SwannR
        self.wann_centers = wann_centers
        self.wann_names = wann_names
        self.atoms = atoms
        self.kpts = kpts
        self.kweights = kweights
        self._Rdict = {tuple(R): i for i, R in enumerate(self.Rlist)}
        # self.check_normalization()
        if wannR is not None:
            self.nwann = self.wannR.shape[2]
        else:
            self.nwann = self.HwannR.shape[1]


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

    def write_to_txt(self, filename):
        """
        write the LWF to txt file.
        """
        wannR_fname = filename + "_wannR.txt"
        with open(wannR_fname, "w") as f:
            for iR, R in enumerate(self.Rlist):
                f.write(f"# R = {R}\n")
                for i in range(self.nwann):
                    f.write(f"# {self.wann_names[i]}:\n")
                    for j in range(self.nbasis):
                        f.write(f"{j:} {self.wannR[iR, j, i]}\n")
                    f.write("\n")
                f.write("\n")


    def write_to_netcdf(self, filename):
        """
        write the LWF to netcdf file.
        """
        import xarray as xr

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

    def get_HR(self, R=None, iR=None, ispin=0):
        """
        get the Hamiltonian at R-vector.
        """
        if iR is None:
            iR = self._Rdict[tuple(R)]
        H = self.HwannR[iR]
        return H

    def get_all_HR(self):
        return self.HwannR

    def get_Hk(self, kpt):
        """
        get the Hamiltonian at k-point.
        """
        Hk = R_to_onek(kpt, self.Rlist, self.HwannR)
        return Hk

    def get_Sk(self, kpt):
        """
        get the overlap matrix at k-point.
        """
        if self.is_orthogonal:
            Sk = None
        else:
            Sk = R_to_onek(kpt, self.Rlist, self.SwannR)
        return Sk

    def solve_k(self, kpt):
        """
        solve the Hamiltonian at k-point with NAC.
        """
        # if np.linalg.norm(kpt) < 1e-6:
        #    Hk = self.get_Hk_noNAC(kpt)
        # else:
        Hk = self.get_Hk(kpt)
        Sk = self.get_Sk(kpt)
        evals, evecs = eigh(Hk, Sk)
        return evals, evecs

    def solve_all(self, kpts):
        """
        solve the Hamiltonian at all k-points.
        """
        evals = []
        evecs = []
        for k in kpts:
            e, v = self.solve_k(k)
            evals.append(e)
            evecs.append(v)
        return np.array(evals), np.array(evecs)

    def HS_and_eigen(self, kpts):
        Hks=[]
        if self.is_orthogonal:
            Sks = None
        else:
            Sks=[]
        evals=[]
        evecs=[]
        for kpt in kpts:
            Hk = self.get_Hk(kpt)
            Hks.append(Hk)
            if not self.is_orthogonal:
                Sk = self.get_Sk(kpt)
                Sks.append(Sk)
            evals, evecs = eigh(Hk)
            evals.append(evals)
            evecs.append(evecs)
        Hks = np.array(Hks)
        evals = np.array(evals)
        evecs = np.array(evecs)
        return Hk, Sks, evals, evecs

    def HSE_k(self, kpt):
        Hk = self.get_Hk(kpt)
        Sk = self.get_Sk(kpt)
        evals, evecs = eigh(Hk, Sk)
        return Hk, Sk, evals, evecs

    def get_HS_and_eigen_k(self, kpt):
        Hk = self.get_Hk(kpt)
        Sk = self.get_Sk(kpt)
        evals, evecs = eigh(Hk, Sk)
        return Hk, Sk, evals, evecs
