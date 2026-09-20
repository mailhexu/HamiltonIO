"""SIESTA WFSX wavefunction-file parser.

Reads the wavefunction coefficients stored by SIESTA in ``.WFSX`` files
(both the old Fortran-unformatted and the netCDF variants are handled by
sisl underneath) and aggregates them into a single :class:`WFSXData` object.

Conventions of the data produced here (matching what SIESTA stores on disk):

- ``kpoints_cart`` is the k-point exactly as stored in the file: Cartesian
  coordinates in inverse Angstrom *including* the 2*pi factor,
  ``k_cart = k_frac @ 2*pi*inv(cell).T``.
- ``kpoints`` are fractional reciprocal coordinates of the cell passed to
  :class:`SiestaWFSXParser` via ``k_frac = k_cart @ cell.T / (2*pi)``. The
  WFSX file itself does not store the lattice, so without ``cell`` the
  fractional k-points are filled with NaN.
- ``eigenvalues`` are the energies as stored (eV, Fermi-shifted by SIESTA).
- ``coefficients`` are in the gauge stored by SIESTA (sisl labels this the
  ``"orbital"`` gauge): ``coefficients[ik, iband, i]`` is the expansion
  coefficient of stored orbital slot ``i`` for band ``iband`` at k-point
  ``ik``.

For scalar (nspin=1) files each state has ``no_u`` coefficients (one per
orbital). For spinor (non-collinear, nspin=4/8) files each state stores the
two spinor components of every orbital in *consecutive* slots, so
``norb == 2*no_u`` and ``coefficients[ik, iband, 2*i]`` /
``coefficients[ik, iband, 2*i + 1]`` are the two components of orbital
``i``.
"""

import warnings
from dataclasses import dataclass

import numpy as np

try:
    import sisl
except Exception as e:  # pragma: no cover
    print(e)
    print("Sisl is not installed, please install it by 'pip install sisl'")


@dataclass
class WFSXData:
    """Aggregated contents of a SIESTA WFSX file.

    Attributes
    ----------
    kpoints : (nk, 3) float ndarray
        Fractional reciprocal coordinates of the cell passed to the parser;
        all NaN when no cell was given.
    kpoints_cart : (nk, 3) float ndarray
        Cartesian k-points in 1/Angstrom with the 2*pi factor included,
        exactly as stored in the file.
    eigenvalues : (nk, nbands) float ndarray
        Eigenvalues in eV as stored (SIESTA writes them Fermi-shifted).
    coefficients : (nk, nbands, norb) complex ndarray
        Wavefunction coefficients in the SIESTA-stored gauge.
    is_spinor : bool
        True when the file comes from a non-collinear (nspin=4/8) run.
    norb : int
        Number of AO coefficients stored per state. For scalar files this is
        the number of orbitals ``no_u``; for spinor files it is ``2*no_u``
        (the two spinor components of each orbital are stored in consecutive
        slots ``2i`` and ``2i+1``).
    """

    kpoints: np.ndarray
    kpoints_cart: np.ndarray
    eigenvalues: np.ndarray
    coefficients: np.ndarray
    is_spinor: bool
    norb: int


class SiestaWFSXParser:
    """Parser for SIESTA WFSX wavefunction files backed by sisl."""

    def __init__(self, wfsx_path: str, cell=None):
        """
        Parameters
        ----------
        wfsx_path : str
            Path to the WFSX file.
        cell : (3, 3) array_like, optional
            Lattice vectors in Angstrom as rows. Needed only to convert the
            stored Cartesian k-points to fractional coordinates; without it
            ``WFSXData.kpoints`` is filled with NaN.
        """
        self.wfsx_path = str(wfsx_path)
        self.cell = None if cell is None else np.asarray(cell, dtype=float)

    def read(self) -> WFSXData:
        """Parse the file once and return the aggregated data."""
        sile = sisl.get_sile(self.wfsx_path)
        nspin, _, nk, _ = sile.read_sizes()
        is_spinor = nspin in (4, 8)
        if nspin == 2:
            raise NotImplementedError(
                "Collinear spin-polarized WFSX files (nspin=2) store two "
                "eigenstate sets per k-point; not supported by this parser."
            )

        kpoints_cart = np.empty((nk, 3), dtype=float)
        kpoints_cart[:] = np.nan
        eigenvalues = []
        coefficients = []
        # We deliberately use the raw stored Cartesian k (no parent lattice):
        # the fractional reduction is done here via `cell`. Silence sisl's
        # per-k warning about not being able to make that conversion itself.
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore", message=".*cannot convert stored k-points.*"
            )
            for ik, st in enumerate(sile.yield_eigenstate()):
                kpoints_cart[ik] = st.info["k"]
                eigenvalues.append(np.asarray(st.eig, dtype=float))
                # EigenstateElectron.state is (nwf, norb) = (n_bands, norb)
                coefficients.append(np.asarray(st.state, dtype=np.complex128))

        eigenvalues = np.asarray(eigenvalues, dtype=float)
        coefficients = np.asarray(coefficients, dtype=np.complex128)

        if self.cell is not None:
            kpoints = kpoints_cart @ self.cell.T / (2.0 * np.pi)
        else:
            kpoints = np.full(kpoints_cart.shape, np.nan)

        return WFSXData(
            kpoints=kpoints,
            kpoints_cart=kpoints_cart,
            eigenvalues=eigenvalues,
            coefficients=coefficients,
            is_spinor=is_spinor,
            norb=coefficients.shape[-1],
        )
