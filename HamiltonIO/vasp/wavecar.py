"""Read standard full-sphere VASP WAVECAR pseudo plane waves and energies.

The licensed WAVECAR is supplied by the caller. pymatgen handles VASP's
Fortran records; this module owns the backend-to-unfolder data conversion.
Gamma-half and noncollinear layouts are rejected rather than silently
misinterpreted as full-sphere scalar coefficients.
"""

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class VaspPWData:
    lattice: np.ndarray  # row vectors, Angstrom
    kpoints: np.ndarray  # fractional reciprocal coordinates
    gvecs: tuple[np.ndarray, ...]
    coefficients: tuple[np.ndarray, ...]  # per k (nspin, nband, 1, npw)
    eigenvalues: tuple[np.ndarray, ...]  # per k (nspin, nband), eV
    fermi_energy: float
    atom_symbols: tuple[str, ...]
    atom_positions: np.ndarray  # supercell fractional positions


def read_wavecar(path, poscar, *, kpoint_indices=None, band_indices=None):
    """Read selected standard-wavefunction k-points and bands.

    ``kpoint_indices`` and ``band_indices`` are zero-based. Selecting a
    subset reduces the returned data but pymatgen itself opens the complete
    WAVECAR internally. The lattice/positions come from the supplied POSCAR.
    """
    from pymatgen.io.vasp.inputs import Poscar
    from pymatgen.io.vasp.outputs import Wavecar

    wave = Wavecar(str(path))
    if wave.vasp_type != "std":
        raise ValueError(
            "only full-sphere scalar WAVECAR (vasp_type='std') is supported"
        )
    structure = Poscar.from_file(str(poscar)).structure
    lattice = np.asarray(wave.a, dtype=float)
    if not np.allclose(lattice, structure.lattice.matrix, atol=1e-5):
        raise ValueError("WAVECAR and POSCAR lattices differ")
    kindex = (
        range(int(wave.nk))
        if kpoint_indices is None
        else [int(i) for i in kpoint_indices]
    )
    bindex = (
        range(int(wave.nb)) if band_indices is None else [int(i) for i in band_indices]
    )
    kindex, bindex = tuple(kindex), tuple(bindex)
    if (
        not kindex
        or not bindex
        or min(kindex) < 0
        or max(kindex) >= wave.nk
        or min(bindex) < 0
        or max(bindex) >= wave.nb
    ):
        raise ValueError("kpoint_indices or band_indices outside WAVECAR range")
    gvecs, coefficients, eigenvalues = [], [], []
    for ik in kindex:
        g = np.asarray(wave.Gpoints[ik], dtype=float)
        if not np.allclose(g, np.rint(g), atol=1e-8):
            raise ValueError(
                "WAVECAR reciprocal vectors are not integer Miller indices"
            )
        gvecs.append(np.asarray(np.rint(g), dtype=int))
        coefficients.append(
            np.asarray(
                [
                    [
                        np.asarray(wave.coeffs[spin][ik][ib], dtype=complex)
                        for ib in bindex
                    ]
                    for spin in range(wave.spin)
                ]
            )[:, :, None, :]
        )
        eigenvalues.append(
            np.asarray(
                [
                    wave.band_energy[spin][ik][list(bindex), 0]
                    for spin in range(wave.spin)
                ]
            )
        )
    return VaspPWData(
        lattice,
        np.asarray([wave.kpoints[ik] for ik in kindex]),
        tuple(gvecs),
        tuple(coefficients),
        tuple(eigenvalues),
        float(wave.efermi),
        tuple(str(site.specie.symbol) for site in structure),
        np.asarray(structure.frac_coords),
    )
