"""Read ABINIT PAW WFKs and their matching JTH PAW XML datasets.

Wavefunction coefficients/eigenvalues are parsed by :mod:`.wfk`; the XML
reader is invoked here, never in the unfolding engine. ABINIT writes no
cprj into an ETSF WFK, so projector overlaps are computed from the
recorded coefficients and the matching XML projector functions.
"""

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from ase.data import chemical_symbols

from .wfk import WFKData, _require_netcdf4, read_wfk


@dataclass(frozen=True)
class AbinitPawData:
    wavefunctions: WFKData
    positions: np.ndarray  # reduced positions in the WFK cell
    species: tuple[str, ...]
    datasets: dict  # species -> pypao PawPseudo
    kweights: np.ndarray


def read_paw_wfk(path, datasets) -> AbinitPawData:
    """Read a full-G PAW WFK and matching ``{symbol: PAW-XML path}``.

    Reject norm-conserving WFKs, absent species datasets, or unsupported
    spinors before any projection is attempted.
    """
    from pypao.libpsp.formats.pawxml import read_paw_xml

    wfk = read_wfk(path)
    if wfk.nspinor != 1:
        raise ValueError(
            "PAW projector projection currently requires scalar WFK coefficients"
        )
    Dataset = _require_netcdf4()
    with Dataset(path) as nc:
        if int(nc.variables["usepaw"][:]) != 1:
            raise ValueError("WFK is not a PAW calculation (usepaw != 1)")
        positions = np.asarray(nc.variables["reduced_atom_positions"][:], dtype=float)
        atom_type = np.asarray(nc.variables["atom_species"][:], dtype=int)
        atomic_numbers = np.asarray(nc.variables["atomic_numbers"][:], dtype=float)
        kweights = np.asarray(nc.variables["kpoint_weights"][:], dtype=float)
    if positions.shape != (len(atom_type), 3) or kweights.shape != (len(wfk.kpoints),):
        raise ValueError("invalid PAW WFK atomic positions or kpoint weights")
    species = tuple(
        chemical_symbols[int(round(atomic_numbers[t - 1]))] for t in atom_type
    )
    missing = set(species) - set(datasets)
    if missing:
        raise ValueError(f"missing matching PAW XML datasets for {sorted(missing)}")
    pseudos = {symbol: read_paw_xml(Path(datasets[symbol])) for symbol in set(species)}
    return AbinitPawData(wfk, positions, species, pseudos, kweights)
