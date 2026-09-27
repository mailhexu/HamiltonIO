"""ABACUS planewave wavefunction/energy file parser.

Reads the per-k planewave wavefunction files written by ABACUS
(``out_wfc_pw 1`` -> ``OUT.${suffix}/WAVEFUNC${K}.txt``) together with
the band eigenvalues (``out_band 1`` -> ``BANDS_1.dat``) and the k-point
list (``OUT.${suffix}/kpoints``), and aggregates them into a single
:class:`AbacusPWData` object mirroring
:class:`HamiltonIO.siesta.wfsx.WFSXData`.

Conventions of the data produced here (matching what ABACUS stores):

- ``kpoints`` are supercell fractional reciprocal coordinates of the
  cell that produced the run (taken from the ``kpoints`` file, so they
  are exact; the Cartesian value stored in the WAVEFUNC header is kept
  as ``kpoints_cart_tpiba`` in ABACUS's internal ``2*pi/alat`` units).
- ``gvecs[ik]`` holds integer Miller indices w.r.t. the supercell
  reciprocal lattice for every stored plane wave.  ABACUS writes FFT-box
  indices (``0 <= idx < n``); they are unfolded to Miller indices with
  the standard wrap ``idx > n/2 -> idx - n``, with the wavefunction FFT
  dimensions parsed from ``running_scf.log``.
- ``coefficients[ik]`` is the ``(nband, npw_k)`` complex coefficient
  array of the stored gauge (``psi_nk = sum_G c(G) e^{i(k+G).r}``;
  coefficients are dimensionless, bands orthonormal).
- ``eigenvalues[ik]`` are the energies in eV as written by
  ``out_band``.

Spin-polarized (nspin=2) runs store two WAVEFUNC files per k-point
(``WAVEFUNC${K}_S${spin}.txt`` in some versions); this parser reads the
unpolarized single-channel layout only, which is what scalar unfolding
needs.
"""

from dataclasses import dataclass
from pathlib import Path

import numpy as np

_WFC_HEADER_TOKENS = 11


@dataclass
class AbacusPWData:
    """Aggregated contents of an ABACUS planewave output directory.

    Attributes
    ----------
    kpoints : (nk, 3) float ndarray
        Supercell fractional k-points, exactly as listed in ``kpoints``.
    kpoints_cart_tpiba : (nk, 3) float ndarray
        Cartesian k-points in ABACUS internal units of ``2*pi/alat``,
        as stored in the WAVEFUNC headers.
    eigenvalues : (nk, nband) float ndarray
        Band energies in eV (from ``BANDS_1.dat``).
    gvecs : list of (npw_k, 3) int ndarray
        Per-k integer Miller indices of the stored plane waves.
    coefficients : list of (nband, npw_k) complex ndarray
        Per-k planewave coefficients (ABACUS-stored gauge).
    fft_grid : (3,) int ndarray
        Wavefunction FFT dimensions used to unfold the stored indices.
    """

    kpoints: np.ndarray
    kpoints_cart_tpiba: np.ndarray
    eigenvalues: np.ndarray
    gvecs: list
    coefficients: list
    fft_grid: np.ndarray


class AbacusPWParser:
    """Parser for ABACUS planewave output directories."""

    def __init__(self, outpath: str):
        self.outpath = Path(outpath)

    # -- public API --------------------------------------------------------

    def read(self) -> AbacusPWData:
        """Parse the output directory once and return the aggregated data."""
        wfc_files = self._find_wfc_files()
        fft_grid = self._read_fft_grid()
        kpoints = self._read_kpoints()
        kpoints_cart = np.zeros((len(wfc_files), 3), dtype=float)
        gvecs, coefficients = [], []
        nband = None
        for ik, fname in enumerate(wfc_files):
            kcart, ng, nband_k, gvec, coeff = _read_wfc_file(fname, fft_grid)
            kpoints_cart[ik] = kcart
            gvecs.append(gvec)
            coefficients.append(coeff)
            if nband is None:
                nband = nband_k
            elif nband_k != nband:
                raise ValueError(f"{fname}: nband={nband_k} differs from {nband}")
        eigenvalues = self._read_bands(nband)
        if eigenvalues.shape[0] != len(wfc_files):
            raise ValueError(
                f"{self.outpath}: {len(wfc_files)} WAVEFUNC files but "
                f"{eigenvalues.shape[0]} k-point rows in BANDS_1.dat"
            )
        if kpoints.shape[0] != len(wfc_files):
            raise ValueError(
                f"{self.outpath}: {len(wfc_files)} WAVEFUNC files but "
                f"{kpoints.shape[0]} k-points in the kpoints file"
            )
        return AbacusPWData(
            kpoints=kpoints,
            kpoints_cart_tpiba=kpoints_cart,
            eigenvalues=eigenvalues,
            gvecs=gvecs,
            coefficients=coefficients,
            fft_grid=np.asarray(fft_grid, dtype=int),
        )

    # -- pieces ------------------------------------------------------------

    def _find_wfc_files(self):
        files = sorted(
            self.outpath.glob("WAVEFUNC*.txt"),
            key=lambda p: int(
                "".join(c for c in p.stem.split("WAVEFUNC")[-1] if c.isdigit()) or 0
            ),
        )
        if not files:
            raise FileNotFoundError(
                f"no WAVEFUNC*.txt files in {self.outpath} "
                "(run ABACUS with out_wfc_pw 1)"
            )
        return files

    def _read_fft_grid(self):
        import re

        wfc_pattern = (
            r"fft grid for wave functions\s*=\s*\[\s*(\d+),\s*(\d+),\s*(\d+)\s*\]"
        )
        rho_pattern = (
            r"fft grid for charge/potential\s*=\s*\[\s*(\d+),\s*(\d+),\s*(\d+)\s*\]"
        )
        pattern = None
        for patterns in ([wfc_pattern], [wfc_pattern, rho_pattern]):
            for name in ("running_scf.log", "running_nscf.log"):
                log = self.outpath / name
                if not log.exists():
                    continue
                text = log.read_text()
                for pat in patterns:
                    match = re.search(pat, text)
                    if match is not None:
                        pattern = match
                        break
                if pattern is not None:
                    break
            if pattern is not None:
                break
        if pattern is None:
            raise ValueError(
                f"could not find the wavefunction FFT grid dimensions in {self.outpath}"
            )
        return tuple(int(g) for g in pattern.groups())

    def _read_kpoints(self):
        """Supercell fractional k-points from the ``kpoints`` file."""
        fname = self.outpath / "kpoints"
        rows = []
        started = False
        with open(fname) as handle:
            for line in handle:
                parts = line.split()
                if not started:
                    started = "K-POINTS DIRECT COORDINATES" in line
                    continue
                if not parts:
                    if rows:
                        break
                    continue
                if len(parts) < 4:
                    continue
                try:
                    rows.append([float(v) for v in parts[1:4]])
                except ValueError:
                    continue
        if not rows:
            raise ValueError(f"no k-points parsed from {fname}")
        return np.asarray(rows, dtype=float)

    def _read_bands(self, nband):
        """Eigenvalues (eV) from ``BANDS_1.dat``; shape (nk, nband)."""
        fname = self.outpath / "BANDS_1.dat"
        rows = []
        with open(fname) as handle:
            for line in handle:
                parts = line.split()
                if len(parts) < 3:
                    continue
                try:
                    values = [float(v) for v in parts[2 : 2 + nband]]
                except ValueError:
                    continue
                if len(values) == nband:
                    rows.append(values)
        if not rows:
            raise FileNotFoundError(
                f"{fname} missing or empty (run ABACUS with out_band 1)"
            )
        return np.asarray(rows, dtype=float)


def _read_wfc_file(fname, fft_grid):
    """Parse one WAVEFUNC${K}.txt file.

    Returns ``(k_cart_tpiba, ng, nband, gvec_miller, coefficients)``.
    """
    with open(fname) as handle:
        header = handle.readline()
        if "Kpoint" not in header:
            raise ValueError(f"{fname}: unexpected header {header!r}")
        values = handle.readline().split()
        if len(values) < _WFC_HEADER_TOKENS:
            raise ValueError(f"{fname}: short header line {values!r}")
        k_cart = np.array(
            [float(values[2]), float(values[3]), float(values[4])], dtype=float
        )
        ng = int(values[6])
        nband = int(values[7])
        line = handle.readline()
        while "<Reciprocal Lattice Vector>" not in line:
            if not line:
                raise ValueError(f"{fname}: missing <Reciprocal Lattice Vector>")
            line = handle.readline()
        recip = np.array(
            [[float(v) for v in handle.readline().split()[:3]] for _ in range(3)]
        )
        if not np.isfinite(recip).all():
            raise ValueError(f"{fname}: unreadable reciprocal lattice block")
        line = handle.readline()
        while line.strip() != "<G vectors>":
            if not line:
                raise ValueError(f"{fname}: missing <G vectors> block")
            line = handle.readline()
        gvec = np.zeros((ng, 3), dtype=int)
        filled = 0
        while filled < ng:
            parts = line.split()
            if len(parts) >= 3:
                try:
                    gvec[filled] = [int(v) for v in parts[:3]]
                    filled += 1
                except ValueError:
                    pass
            line = handle.readline()
            if not line:
                raise ValueError(f"{fname}: truncated <G vectors> block")
        coeff = np.zeros((nband, ng), dtype=complex)
        for ib in range(nband):
            # Skip everything (blank lines and the closing duplicate
            # "< Band ib >" of the previous band) up to the band marker.
            marker = ""
            while f"< Band {ib + 1} >" not in marker:
                marker = handle.readline()
                if not marker:
                    raise ValueError(f"{fname}: expected < Band {ib + 1} > marker")
            got = 0
            while got < ng:
                line = handle.readline()
                if not line:
                    raise ValueError(f"{fname}: truncated band {ib + 1}")
                parts = line.split()
                if len(parts) % 2:
                    parts = parts[:-1]
                if not parts:
                    continue
                vals = np.asarray([float(v) for v in parts])
                pairs = len(vals) // 2
                if not pairs:
                    continue
                take = min(pairs, ng - got)
                coeff[ib, got : got + take] = (
                    vals[0 : 2 * take : 2] + 1j * vals[1 : 2 * take : 2]
                )
                got += take
    # FFT-box index -> Miller index (standard wrap: idx > n/2 <-> idx - n)
    miller = np.asarray(gvec, dtype=float)
    for axis in range(3):
        n = fft_grid[axis]
        half = n // 2
        mask = gvec[:, axis] > half
        miller[mask, axis] = gvec[mask, axis] - n
    return k_cart, ng, nband, miller.astype(int), coeff
