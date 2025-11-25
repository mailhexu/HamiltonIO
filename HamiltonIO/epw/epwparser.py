"""
Parser for EPW (Electron-Phonon Wannier) files.

This module provides functionality to parse and handle various EPW-related files
generated from Quantum Espresso's EPW code. It includes parsers for:
- WSVec files containing Wigner-Seitz grid information
- Crystal format files containing structural information
- EPW data format files
- EPW matrix elements

The module provides both low-level parsing functions and high-level classes
for handling the EPW data structure.
"""

import os
from collections import defaultdict
from dataclasses import dataclass, field

import numpy as np
from ase.units import Bohr, Ry
from netCDF4 import Dataset

from HamiltonIO.wannier.w90_parser import parse_ham


def line_to_array(line, fmt=float):
    """Convert a string line to a numpy array with given format.

    Parameters
    ----------
    line : str
        The input string containing space-separated values
    fmt : type, optional
        The format to convert elements to, by default float

    Returns
    -------
    numpy.ndarray
        Array containing the converted values
    """
    return np.array([fmt(x) for x in line.split()])


def line2vec(line):
    """Convert a string line to a vector of integers.

    Parameters
    ----------
    line : str
        The input string containing space-separated integers

    Returns
    -------
    list
        List of integers parsed from the input line
    """
    return [int(x) for x in line.strip().split()]


def read_WSVec_deprecated(fname):
    """Read a Wigner-Seitz vector file (deprecated format).

    This function reads the old WSVec format with "=" separators.
    Use read_WSVec() for the new wigner.fmt format.

    Parameters
    ----------
    fname : str
        Path to the WSVec file

    Returns
    -------
    tuple
        Contains:
        - dims : int
            First dimension parameter
        - dims2 : int
            Second dimension parameter
        - nRk : int
            Number of R vectors for k-space
        - nRq : int
            Number of R vectors for q-space
        - nRg : int
            Number of R vectors for g-space
        - Rk : numpy.ndarray
            R vectors in k-space
        - Rq : numpy.ndarray
            R vectors in q-space
        - Rg : numpy.ndarray
            R vectors in g-space
        - ndegen_k : numpy.ndarray
            Degeneracy of k-space vectors
        - ndegen_q : numpy.ndarray
            Degeneracy of q-space vectors
        - ndegen_g : numpy.ndarray
            Degeneracy of g-space vectors
    """
    with open(fname) as myfile:
        lines = myfile.readlines()

    (dims, dims2, nRk, nRq, nRg) = tuple(int(line.split("=")[1]) for line in lines[:5])

    Rk = np.zeros((nRk, 3), dtype=int)
    Rq = np.zeros((nRq, 3), dtype=int)
    Rg = np.zeros((nRg, 3), dtype=int)
    ndegen_k = np.zeros(nRk, dtype=int)
    ndegen_q = np.zeros(nRq, dtype=int)
    ndegen_g = np.zeros(nRg, dtype=int)

    start = 6
    end = start + nRk
    for i, line in enumerate(lines[start:end]):
        Rk[i] = line2vec(line)

    start = end + 1
    end = start + nRq
    for i, line in enumerate(lines[start:end]):
        Rq[i] = line2vec(line)

    start = end + 1
    end = start + nRg
    for i, line in enumerate(lines[start:end]):
        Rg[i] = line2vec(line)

    start = end + 1
    end = start + nRk
    for i, line in enumerate(lines[start:end]):
        ndegen_k[i] = int(line.strip())

    start = end + 1
    end = start + nRq
    for i, line in enumerate(lines[start:end]):
        ndegen_q[i] = int(line.strip())

    start = end + 1
    end = start + nRg
    for i, line in enumerate(lines[start:end]):
        ndegen_g[i] = int(line.strip())

    return (dims, dims2, nRk, nRq, nRg, Rk, Rq, Rg, ndegen_k, ndegen_q, ndegen_g)


def read_WSVec(fname):
    """Read a Wigner-Seitz vector file using the modern wigner.py module.

    This function uses the WignerData class to read the new wigner.fmt format
    and returns data in the same format as the deprecated function for compatibility.

    Parameters
    ----------
    fname : str
        Path to the wigner.fmt file

    Returns
    -------
    tuple
        Contains:
        - dims : int
            Number of Wannier functions
        - dims2 : int
            Number of atoms
        - nRk : int
            Number of R vectors for k-space
        - nRq : int
            Number of R vectors for q-space
        - nRg : int
            Number of R vectors for g-space
        - Rk : numpy.ndarray
            R vectors in k-space, shape (nRk, 3)
        - Rq : numpy.ndarray
            R vectors in q-space, shape (nRq, 3)
        - Rg : numpy.ndarray
            R vectors in g-space, shape (nRg, 3)
        - ndegen_k : numpy.ndarray
            Degeneracy of k-space vectors (flattened for compatibility)
        - ndegen_q : numpy.ndarray
            Degeneracy of q-space vectors (flattened for compatibility)
        - ndegen_g : numpy.ndarray
            Degeneracy of g-space vectors (flattened for compatibility)
    """
    from HamiltonIO.epw.wigner import WignerData

    # Read data using the modern WignerData class
    wigner_data = WignerData.from_file(fname)

    # Extract data in the format expected by the old interface
    dims = wigner_data.dims
    dims2 = wigner_data.dims2
    nRk = wigner_data.nrr_k
    nRq = wigner_data.nrr_q
    nRg = wigner_data.nrr_g

    # R vectors are already in the correct Pythonic format (nR, 3)
    Rk = wigner_data.irvec_k
    Rq = wigner_data.irvec_q
    Rg = wigner_data.irvec_g

    # For compatibility with the old interface, we need to flatten the degeneracy arrays
    # The old format expected 1D arrays, but the new format has multi-dimensional arrays
    # We'll take the diagonal elements for the case where dims=dims2=1
    if dims == 1 and dims2 == 1:
        ndegen_k = wigner_data.ndegen_k[:, 0, 0]
        ndegen_q = wigner_data.ndegen_q[:, 0, 0]
        ndegen_g = wigner_data.ndegen_g[:, 0, 0]
    else:
        # For more complex cases, we might need to handle this differently
        # For now, we'll take the first element of each degeneracy matrix
        ndegen_k = wigner_data.ndegen_k[:, 0, 0]
        ndegen_q = wigner_data.ndegen_q[:, 0, 0]
        ndegen_g = wigner_data.ndegen_g[:, 0, 0]

    return (dims, dims2, nRk, nRq, nRg, Rk, Rq, Rg, ndegen_k, ndegen_q, ndegen_g)


@dataclass
class Crystal:
    """Class storing crystal structure information from EPW calculation.

    Attributes
    ----------
    natom : int
        Number of atoms in the unit cell
    nmode : int
        Number of phonon modes
    nelect : float
        Number of electrons
    at : numpy.ndarray
        Real space lattice vectors
    bg : numpy.ndarray
        Reciprocal space lattice vectors
    omega : float
        Unit cell volume
    alat : float
        Lattice parameter
    tau : numpy.ndarray
        Atomic positions
    amass : numpy.ndarray
        Atomic masses
    ityp : numpy.ndarray
        Atomic types
    noncolin : bool
        Whether the calculation is non-collinear
    w_centers : numpy.ndarray
        Wannier function centers
    """

    natom: int = 0
    nmode: int = 0
    nelect: float = 0.0
    nbandskip: int = 0
    at: np.ndarray = field(default_factory=lambda: np.zeros(0))
    bg: np.ndarray = field(default_factory=lambda: np.zeros(0))
    omega: float = 0.0
    alat: float = 0.0
    tau: np.ndarray = field(default_factory=lambda: np.zeros(0))
    amass: np.ndarray = field(default_factory=lambda: np.zeros(0))
    ityp: np.ndarray = field(default_factory=lambda: np.zeros(0))
    noncolin: bool = False
    do_cutoff_2D_epw: bool = False
    w_centers: np.ndarray = field(default_factory=lambda: np.zeros(0))
    L: float = 0.0


def is_text_True(s):
    """Check if a string represents a boolean True value.

    Parameters
    ----------
    s : str
        Input string to check

    Returns
    -------
    bool
        True if string starts with 't' or 'T', False otherwise
    """
    return s.strip().lower().startswith("t")


def read_crystal_fmt(fname="crystal.fmt"):
    """Parse the crystal.fmt file containing crystal structure information.

    Parameters
    ----------
    fname : str, optional
        Path to the crystal.fmt file, by default "crystal.fmt"

    Returns
    -------
    Crystal
        Crystal object containing the parsed structural information
    """
    """
    parser to the crystal.fmt file
    see line 114 (qe version 6.8) epw_write in ephwann_shuffle.f90.
    """
    d = Crystal()
    with open(fname) as myfile:
        d.natom = int(next(myfile))
        d.nmode = int(next(myfile))
        s = next(myfile).strip().split()
        d.nelect = float(s[0])
        if len(s) > 1:
            d.nbndskip = int(s[1])
        d.at = line_to_array(next(myfile), float)
        d.bg = line_to_array(next(myfile), float)
        d.omega = float(next(myfile))
        d.alat = float(next(myfile))
        d.tau = line_to_array(next(myfile), float)
        d.amass = line_to_array(next(myfile), float)
        d.ityp = line_to_array(next(myfile), int)
        d.noncolin = is_text_True(next(myfile))
        d.do_cutoff_2D_epw = is_text_True(next(myfile))
        d.w_centers = line_to_array(next(myfile), float)
        d.L = float(next(myfile))
    return d


def read_epwdata_fmt(fname="epwdata.fmt"):
    """Read the EPW data format file containing basic dimensions.

    Parameters
    ----------
    fname : str, optional
        Path to the epwdata.fmt file, by default "epwdata.fmt"

    Returns
    -------
    tuple
        Contains:
        - nbndsub : int
            Number of bands
        - nrr_k : int
            Number of R vectors for k-space
        - nmodes : int
            Number of phonon modes
        - nrr_q : int
            Number of R vectors for q-space
        - nrr_g : int
            Number of R vectors for g-space
    """
    with open(fname) as myfile:
        _efermi = float(next(myfile))
        nbndsub, nrr_k, nmodes, nrr_q, nrr_g = [int(x) for x in next(myfile).split()]
    return nbndsub, nrr_k, nmodes, nrr_q, nrr_g


def read_epmatwp(fname="./sic.epmatwp", path="./"):
    """Read the EPW matrix elements file.

    Parameters
    ----------
    fname : str, optional
        Name of the epmatwp file, by default "./sic.epmatwp"
    path : str, optional
        Path to the directory containing the file, by default './'

    Returns
    -------
    numpy.ndarray
        5D array containing the EPW matrix elements with shape
        (nrr_g, nmodes, nrr_k, nbndsub, nbndsub)
    """
    nbndsub, nrr_k, nmodes, nrr_q, nrr_g = read_epwdata_fmt(
        os.path.join(path, "epwdata.fmt")
    )
    # Is this H or H.T?
    mat = np.fromfile(os.path.join(path, fname), dtype=complex)
    return np.reshape(mat, (nrr_g, nmodes, nrr_k, nbndsub, nbndsub), order="C")


@dataclass
class EpmatOneMode:
    """Class for handling EPW matrix elements for a single phonon mode.

    This class provides functionality to store and manipulate electron-phonon
    matrix elements for a specific phonon mode.

    Attributes
    ----------
    nwann : int
        Number of Wannier functions
    nRk : int
        Number of R vectors in k-space
    nRq : int
        Number of R vectors in q-space
    nRg : int
        Number of R vectors in g-space
    Rk : numpy.ndarray
        R vectors in k-space
    Rq : numpy.ndarray
        R vectors in q-space
    Rg : numpy.ndarray
        R vectors in g-space
    ndegen_k : numpy.ndarray
        Degeneracy of k-space vectors
    ndegen_q : numpy.ndarray
        Degeneracy of q-space vectors
    ndegen_g : numpy.ndarray
        Degeneracy of g-space vectors
    """

    nwann: int = 0
    nRk: int = None
    nRq: int = None
    nRg: int = None
    Rk: np.ndarray = field(default=None)
    Rq: np.ndarray = field(default=None)
    Rg: np.ndarray = field(default=None)
    ndegen_k: np.ndarray = field(default=None)
    ndegen_q: np.ndarray = field(default=None)
    ndegen_g: np.ndarray = field(default=None)

    def __init__(self, epmat, imode, close_nc=False):
        """Initialize EPW matrix elements for a single mode.

        Parameters
        ----------
        epmat : Epmat
            EPW matrix elements object containing data for all modes
        imode : int
            Index of the phonon mode to extract
        close_nc : bool, optional
            Whether to close the netCDF file after reading, by default False
        """
        self.nwann = epmat.nwann
        self.nRk = epmat.nRk
        self.nRq = epmat.nRq
        self.nRg = epmat.nRg
        self.Rk = epmat.Rk
        self.Rq = epmat.Rq
        self.Rg = epmat.Rg

        self.Rkdict = epmat.Rkdict
        self.Rqdict = epmat.Rqdict
        self.Rgdict = epmat.Rgdict

        self.ndegen_k = epmat.ndegen_k
        self.ndegen_q = epmat.ndegen_q
        self.ndegen_g = epmat.ndegen_g

        self.data = np.zeros(
            (self.nRg, self.nRk, self.nwann, self.nwann), dtype=complex
        )
        for Rg, iRg in self.Rgdict.items():
            dreal = epmat.epmatfile.variables["epmat_real"][iRg, imode, :, :, :] * (
                Ry / Bohr
            )
            dimag = epmat.epmatfile.variables["epmat_imag"][iRg, imode, :, :, :] * (
                Ry / Bohr
            )
            self.data[iRg] = dreal + 1.0j * dimag
            # self.data = np.swapaxes(self.data, 2, 3)

        if close_nc:
            epmat.epmatfile.close()

    def get_epmat_RgRk(self, Rg, Rk, avg=False):
        """Get EPW matrix elements for given R vectors in k and g space.

        Parameters
        ----------
        Rg : tuple
            R vector in g-space
        Rk : tuple
            R vector in k-space
        avg : bool, optional
            Whether to average with the time-reversed counterpart, by default False

        Returns
        -------
        numpy.ndarray
            Complex array of EPW matrix elements
        """
        iRg = self.Rgdict[tuple(Rg)]
        iRk = self.Rkdict[tuple(Rk)]
        ret = np.copy(self.data[iRg, iRk, :, :])

        if avg:
            Rg2 = tuple(np.array(Rg) - np.array(Rk))
            Rk2 = tuple(-np.array(Rk))
            if Rg2 in self.Rgdict and Rk in self.Rkdict:
                iRg2 = self.Rgdict[tuple(Rg2)]
                iRk2 = self.Rkdict[tuple(Rk2)]
                ret += self.data[iRg2, iRk2, :, :].T
                ret /= 2.0
        return ret

    def get_epmat_RgRk_two_spin(self, Rg, Rk, avg=False):
        """Get EPW matrix elements for two-spin case.

        Parameters
        ----------
        Rg : tuple
            R vector in g-space
        Rk : tuple
            R vector in k-space
        avg : bool, optional
            Whether to average with the time-reversed counterpart, by default False

        Returns
        -------
        numpy.ndarray
            Complex array of EPW matrix elements for two-spin case
        """
        H = self.get_epmat_RgRk(Rg, Rk, avg=avg)
        ret = np.zeros((self.nwann * 2, self.nwann * 2), dtype=complex)
        ret[::2, ::2] = H
        ret[1::2, 1::2] = H
        return ret


@dataclass
class Epmat:
    """Class for handling the complete set of EPW matrix elements.

    This class provides functionality to read, store, and manipulate
    the full set of electron-phonon matrix elements for all modes.

    Attributes
    ----------
    crystal : Crystal
        Crystal structure information
    nwann : int
        Number of Wannier functions
    nmodes : int
        Number of phonon modes
    nRk : int
        Number of R vectors in k-space
    nRq : int
        Number of R vectors in q-space
    nRg : int
        Number of R vectors in g-space
    Rk : numpy.ndarray
        R vectors in k-space
    Rq : numpy.ndarray
        R vectors in q-space
    Rg : numpy.ndarray
        R vectors in g-space
    ndegen_k : numpy.ndarray
        Degeneracy of k-space vectors
    ndegen_q : numpy.ndarray
        Degeneracy of q-space vectors
    ndegen_g : numpy.ndarray
        Degeneracy of g-space vectors
    Hwann : numpy.ndarray
        Wannier Hamiltonian
    Rlist : numpy.ndarray
        List of R vectors
    epmat_wann : numpy.ndarray
        EPW matrix elements in Wannier gauge
    epmat_ncfile : Dataset
        NetCDF file containing EPW matrix elements
    """

    crystal: Crystal = field(default_factory=Crystal)
    nwann: int = 0
    nmodes: int = 0
    nRk: int = None
    nRq: int = None
    nRg: int = None
    Rk: np.ndarray = field(default=None)
    Rq: np.ndarray = field(default=None)
    Rg: np.ndarray = field(default=None)
    ndegen_k: np.ndarray = field(default=None)
    ndegen_q: np.ndarray = field(default=None)
    ndegen_g: np.ndarray = field(default=None)
    Hwann: np.ndarray = field(default=None)
    Rlist: np.ndarray = field(default=None)
    epmat_wann: np.ndarray = field(default=None)
    epmat_ncfile: Dataset = field(default=None)

    Rk_dict: dict = field(default_factory=dict)
    Rq_dict: dict = field(default_factory=dict)
    Rg_dict: dict = field(default_factory=dict)

    def read_Rvectors(self, path, fname="wigner.fmt"):
        """Read R vectors from a Wigner-Seitz vector file.

        This method uses the modern wigner.py implementation to read
        Wigner-Seitz vectors in the new format.

        Parameters
        ----------
        path : str
            Directory containing the file
        fname : str, optional
            Name of the Wigner-Seitz vector file, by default "wigner.fmt"
            For legacy format, use "WSVecDeg.dat" and read_WSVec_deprecated()
        """
        fullfname = os.path.join(path, fname)

        # Use the modern read_WSVec function which leverages wigner.py
        (
            dims,
            dims2,
            self.nRk,
            self.nRq,
            self.nRg,
            self.Rk,
            self.Rq,
            self.Rg,
            self.ndegen_k,
            self.ndegen_q,
            self.ndegen_g,
        ) = read_WSVec(fullfname)

        # Create dictionaries for fast R-vector lookup
        self.Rkdict = {tuple(self.Rk[i]): i for i in range(self.nRk)}
        self.Rqdict = {tuple(self.Rq[i]): i for i in range(self.nRq)}
        self.Rgdict = {tuple(self.Rg[i]): i for i in range(self.nRg)}

        # Store dimensions for reference
        self.dims = dims
        self.dims2 = dims2

    def read_Wannier_Hamiltonian(self, path, fname):
        """Read Wannier Hamiltonian from file.

        Parameters
        ----------
        path : str
            Directory containing the file
        fname : str
            Name of the Hamiltonian file
        """
        nwann, HR = parse_ham(fname=os.path.join(path, fname))
        nR = len(HR)
        self.Rlist = np.array(list(HR.keys()), dtype=int)
        self.Hwann = np.zeros((nR, nwann, nwann), dtype=complex)
        for i, H in enumerate(HR.values()):
            self.Hwann[i] = H

    def read_epmat(self, path, prefix):
        """Read EPW matrix elements from binary file.

        Parameters
        ----------
        path : str
            Directory containing the file
        prefix : str
            Prefix for the EPW files
        """
        mat = np.fromfile(os.path.join(path, f"{prefix}.epmatwp"), dtype=complex)
        # mat = mat.reshape((self.nRg, self.nmodes, self.nRk,
        #                   self.nwann, self.nwann), order='F')
        mat.shape = (self.nRg, self.nmodes, self.nRk, self.nwann, self.nwann)
        self.epmat_wann = mat

    def read(self, path, prefix, epmat_ncfile=None):
        """Read all necessary EPW data.

        Parameters
        ----------
        path : str
            Directory containing the files
        prefix : str
            Prefix for the EPW files
        epmat_ncfile : str, optional
            Name of netCDF file containing EPW matrix elements, by default None
        """
        # self.crystal = read_crystal_fmt(
        #    fname=os.path.join(path, 'crystal.fmt'))
        # self.read_Wannier_Hamiltonian(path,  f"{prefix}_hr.dat")
        (self.nwann, self.nRk, self.nmodes, self.nRq, self.nRg) = read_epwdata_fmt(
            os.path.join(path, "epwdata.fmt")
        )
        self.read_Rvectors(path, "wigner.fmt")
        if epmat_ncfile is not None:
            self.epmatfile = Dataset(os.path.join(path, epmat_ncfile), "r")
        else:
            self.read_epmat(path, prefix)

    def get_epmat_Rv_from_index(self, imode, iRg, iRk):
        """Get EPW matrix elements for given mode and R vector indices.

        Parameters
        ----------
        imode : int
            Mode index
        iRg : int
            Index of R vector in g-space
        iRk : int
            Index of R vector in k-space

        Returns
        -------
        numpy.ndarray
            Complex array of EPW matrix elements
        """
        dreal = self.epmatfile.variables["epmat_real"][iRg, imode, iRk, :, :]
        dimag = self.epmatfile.variables["epmat_imag"][iRg, imode, iRk, :, :]
        return (dreal + 1.0j * dimag) * (Ry / Bohr)

    def get_epmat_Rv_from_RgRk(self, imode, Rg, Rk):
        """Get EPW matrix elements for given mode and R vectors.

        Parameters
        ----------
        imode : int
            Mode index
        Rg : tuple
            R vector in g-space
        Rk : tuple
            R vector in k-space

        Returns
        -------
        numpy.ndarray
            Complex array of EPW matrix elements
        """
        iRg = self.Rgdict[tuple(Rg)]
        iRk = self.Rkdict[tuple(Rk)]
        dreal = self.epmatfile.variables["epmat_real"][iRg, imode, iRk, :, :]
        dimag = self.epmatfile.variables["epmat_imag"][iRg, imode, iRk, :, :]
        return (dreal + 1.0j * dimag) * (Ry / Bohr)

    def get_epmat_Rv_from_R(self, imode, Rg):
        """Get EPW matrix elements for given mode and g-space R vector.

        Parameters
        ----------
        imode : int
            Mode index
        Rg : tuple
            R vector in g-space

        Returns
        -------
        dict
            Dictionary mapping k-space R vectors to EPW matrix elements
        """
        iRg = self.Rgdict[tuple(Rg)]
        g = defaultdict(lambda: np.zeros((self.nwann, self.nwann)))
        for Rk in self.Rkdict:
            iRk = self.Rkdict[tuple(Rk)]
            g[tuple(Rk)] = self.get_epmat_Rv_from_index(imode, iRg, iRk)
        return g

    def distance_resolved_couplings_Rk(self, imode, cell, use_absolute=True):
        """Return electron-phonon couplings annotated with Rk-based distances.

        This method extracts g(imode, Rg, Rk, i, j) matrix elements and computes
        the distance between electron Wannier functions (WF to WF distance).

        Parameters
        ----------
        imode : int
            Phonon mode index
        cell : array-like
            3x3 lattice vectors in Angstrom (real-space lattice)
        use_absolute : bool, optional
            If True (default), return scalar distance norm.
            If False, return 3D distance vector.

        Returns
        -------
        list[dict]
            Each item contains:
                - "imode": phonon mode index
                - "Rg": tuple of Rg vector
                - "Rk": tuple of Rk vector
                - "i": Wannier function index i
                - "j": Wannier function index j
                - "distance": distance between WF i and WF j+Rk
                - "coupling": complex coupling matrix element g
        """
        # Reshape w_centers to (nwann, 3) if needed
        w_centers = np.asarray(self.crystal.w_centers).reshape(self.nwann, 3)
        cell = np.asarray(cell)

        results = []
        # Iterate over all Rg vectors
        for Rg in self.Rgdict.keys():
            iRg = self.Rgdict[Rg]
            # Iterate over all Rk vectors
            for Rk in self.Rkdict.keys():
                iRk = self.Rkdict[Rk]
                # Get coupling matrix for this (imode, Rg, Rk)
                g_mat = self.get_epmat_Rv_from_index(imode, iRg, iRk)

                # Iterate over Wannier function pairs
                for i in range(self.nwann):
                    for j in range(self.nwann):
                        g_ij = g_mat[i, j]
                        if g_ij != 0:  # Only store non-zero couplings
                            # Distance from WF i to WF j+Rk (electron hopping-like)
                            d_frac = w_centers[j] + np.asarray(Rk) - w_centers[i]
                            d_cart = d_frac @ cell
                            if use_absolute:
                                dist = float(np.linalg.norm(d_cart))
                            else:
                                dist = d_cart

                            results.append(
                                {
                                    "imode": int(imode),
                                    "Rg": tuple(Rg),
                                    "Rk": tuple(Rk),
                                    "i": int(i),
                                    "j": int(j),
                                    "distance": dist,
                                    "coupling": complex(g_ij),
                                }
                            )
        return results

    def distance_resolved_couplings_Rg(self, imode, cell, use_absolute=True):
        """Return electron-phonon couplings annotated with Rg-based distances.

        This method extracts g(imode, Rg, Rk, i, j) matrix elements and computes
        the distance between electron Wannier function and atomic displacement.

        Parameters
        ----------
        imode : int
            Phonon mode index (0-based)
        cell : array-like
            3x3 lattice vectors in Angstrom (real-space lattice)
        use_absolute : bool, optional
            If True (default), return scalar distance norm.
            If False, return 3D distance vector.

        Returns
        -------
        list[dict]
            Each item contains:
                - "imode": phonon mode index
                - "atom_index": index of atom associated with this mode
                - "Rg": tuple of Rg vector
                - "Rk": tuple of Rk vector
                - "i": Wannier function index i
                - "j": Wannier function index j
                - "distance": distance between WF i and atom+Rg
                - "coupling": complex coupling matrix element g
        """
        # Reshape positions
        w_centers = np.asarray(self.crystal.w_centers).reshape(self.nwann, 3)
        tau = np.asarray(self.crystal.tau).reshape(self.crystal.natom, 3)
        cell = np.asarray(cell)

        # Map mode to atom: each atom has 3 modes (x, y, z displacements)
        atom_index = imode // 3

        results = []
        # Iterate over all Rg vectors
        for Rg in self.Rgdict.keys():
            iRg = self.Rgdict[Rg]
            # Iterate over all Rk vectors
            for Rk in self.Rkdict.keys():
                iRk = self.Rkdict[Rk]
                # Get coupling matrix for this (imode, Rg, Rk)
                g_mat = self.get_epmat_Rv_from_index(imode, iRg, iRk)

                # Iterate over Wannier function pairs
                for i in range(self.nwann):
                    for j in range(self.nwann):
                        g_ij = g_mat[i, j]
                        if g_ij != 0:  # Only store non-zero couplings
                            # Distance from WF i to atom+Rg
                            # (electron-phonon coupling locality)
                            d_frac = tau[atom_index] + np.asarray(Rg) - w_centers[i]
                            d_cart = d_frac @ cell
                            if use_absolute:
                                dist = float(np.linalg.norm(d_cart))
                            else:
                                dist = d_cart

                            results.append(
                                {
                                    "imode": int(imode),
                                    "atom_index": int(atom_index),
                                    "Rg": tuple(Rg),
                                    "Rk": tuple(Rk),
                                    "i": int(i),
                                    "j": int(j),
                                    "distance": dist,
                                    "coupling": complex(g_ij),
                                }
                            )
        return results

    @staticmethod
    def bin_couplings_by_distance(entries, dr=0.1):
        """Bin coupling magnitudes as a function of distance.

        Parameters
        ----------
        entries : list[dict]
            Output of ``distance_resolved_couplings_Rk`` or
            ``distance_resolved_couplings_Rg``.
        dr : float, optional
            Bin width in Angstrom, by default 0.1

        Returns
        -------
        tuple[np.ndarray, np.ndarray]
            (bin_centers, avg_abs_coupling) where avg_abs_coupling is the
            average |g| in each distance bin.
        """
        if not entries:
            return np.array([]), np.array([])

        distances = np.array([e["distance"] for e in entries], dtype=float)
        mags = np.abs([e["coupling"] for e in entries])

        dmax = float(distances.max())
        nbins = int(np.ceil(dmax / dr)) or 1
        bin_edges = np.linspace(0, dmax, nbins + 1)
        bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

        bin_indices = np.digitize(distances, bin_edges) - 1
        bin_indices = np.clip(bin_indices, 0, nbins - 1)

        avg_vals = np.zeros(nbins)
        for i in range(nbins):
            mask = bin_indices == i
            if mask.any():
                avg_vals[i] = mags[mask].mean()

        return bin_centers, avg_vals

    def print_info(self):
        """Print basic information about the EPW matrix elements."""
        print(
            f"nwann: {self.nwann}\n nRk:{self.nRk}, "
            f"nmodes:{self.nmodes}, nRg:{self.nRg}"
        )

    def save_to_netcdf(self, path, fname):
        """Save EPW matrix elements to a netCDF file.

        Parameters
        ----------
        path : str
            Directory to save the file in
        fname : str
            Name of the output netCDF file
        """
        root = Dataset(os.path.join(path, fname), "w")
        root.createDimension("nwann", self.nwann)
        root.createDimension("nRk", self.nRk)
        root.createDimension("nRq", self.nRq)
        root.createDimension("nRg", self.nRg)
        root.createDimension("nmodes", self.nmodes)
        root.createDimension("three", 3)
        root.createVariable(
            "epmat_real",
            "double",
            dimensions=("nRg", "nmodes", "nRk", "nwann", "nwann"),
            zlib=False,
        )
        root.createVariable(
            "epmat_imag",
            "double",
            dimensions=("nRg", "nmodes", "nRk", "nwann", "nwann"),
            zlib=False,
        )

        root.variables["epmat_real"][:] = self.epmat_wann.real
        root.variables["epmat_imag"][:] = self.epmat_wann.imag
        print("ncfile written")
        root.close()


def save_epmat_to_nc(path, prefix, ncfile="epmat.nc"):
    """Save EPW matrix elements to a NetCDF file.

    Parameters
    ----------
    path : str
        Path to the directory containing EPW files
    prefix : str
        Prefix for the EPW files
    ncfile : str, optional
        Name of the output NetCDF file, by default 'epmat.nc'
    """
    ep = Epmat()
    ep.read(path=path, prefix=prefix)
    ep.save_to_netcdf(path=path, fname=ncfile)


def test():
    """Test function for saving EPW matrix elements to netCDF format.

    Uses a sample SiC dataset to demonstrate the file conversion process.
    """
    path = "/home/hexu/projects/epw_test/sic_small_kq/NM/epw"
    save_epmat_to_nc(path=path, prefix="sic", ncfile="epmat.nc")


def test_read_data():
    """Test function for reading and manipulating EPW data.

    Demonstrates loading EPW data from a netCDF file and performing
    various operations on the matrix elements, including time-reversal
    symmetry checks.
    """
    path = "/home/hexu/projects/SrMnO3/spinphon_data/epw555_10"
    prefix = "SrMnO3"
    # path = "/home/hexu/projects/epw_test/sic_small_kq/NM/epw"
    # prefix = 'sic'
    ep = Epmat()
    ep.read(path=path, prefix=prefix, epmat_ncfile="epmat.nc")
    # ep.get_epmat_Rv_from_index(0, 0)
    # d = ep.get_epmat_Rv_from_R(0, (0, 0, 0))
    # print(d[(0, 0, 2)].imag)

    ep1mode = EpmatOneMode(ep, imode=3)
    Rg = (0, 0, 0)
    Rk = (0, 1, 0)
    Rg2 = tuple(np.array(Rg) - np.array(Rk))
    Rk2 = tuple(-np.array(Rk))
    dv1 = ep1mode.get_epmat_RgRk(Rg=Rg, Rk=Rk, avg=False).real
    dv2 = ep1mode.get_epmat_RgRk(Rg=Rg2, Rk=Rk2, avg=False).real.T
    print(dv1)
    print("-" * 10)
    print(dv2)
    print("-" * 10)
    print(dv1 - dv2)
