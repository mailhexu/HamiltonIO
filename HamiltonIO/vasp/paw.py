"""Read VASP PAW POTCAR projector and partial-wave radial tables.

Only standard PAW datasets with tabulated reciprocal-space nonlocal
projectors and ``PAW radial sets`` are accepted. The POTCAR is licensed
VASP input; callers must supply their own file. Nothing here redistributes it.
"""

import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.integrate import simpson


@dataclass(frozen=True)
class VaspPawDataset:
    symbol: str
    angular_momenta: np.ndarray  # one per radial projector
    q_grid: np.ndarray  # 1/Angstrom, radial reciprocal grid
    projectors_q: np.ndarray  # (nradial, 100) reciprocal tabulation
    radial_grid: np.ndarray  # Bohr
    partial_ae: np.ndarray  # u=rR
    partial_pseudo: np.ndarray  # u=rR
    overlap_correction: np.ndarray  # full (l,m) onsite dS


def _numbers(text):
    try:
        return np.array([float(token.replace("D", "E")) for token in text.split()])
    except ValueError as exc:
        raise ValueError("malformed numeric block in POTCAR") from exc


def _radial_sections(text):
    """Return labelled numeric blocks after ``PAW radial sets``."""
    sections = {}
    current = None
    for line in text.splitlines():
        stripped = line.strip()
        if stripped in ("grid", "pseudo wavefunction", "ae wavefunction"):
            current = stripped
            sections.setdefault(current, []).append([])
        elif (
            current
            and stripped
            and all(
                re.fullmatch(r"[+-]?(?:\d*\.\d+|\d+\.?\d*)(?:[eEdD][+-]?\d+)?", word)
                for word in stripped.split()
            )
        ):
            sections[current][-1].extend(_numbers(stripped))
        elif current:
            current = None
    return sections


def read_potcar_paw(path, symbol):
    """Parse the requested species' PAW projectors and overlap correction.

    VASP's 100-point reciprocal tables use ``q_i=i*gmax/100`` in 1/Angstrom.
    PAW radial partial waves are ``u=rR`` on the supplied logarithmic Bohr
    grid. The monopole correction is the radial AE-minus-pseudo overlap.
    """
    text = Path(path).read_text()
    chunks = text.split("End of Dataset")
    matching = [
        chunk
        for chunk in chunks
        if chunk.strip().splitlines()
        and chunk.strip().splitlines()[0].split()[1].split("_")[0] == symbol
    ]
    if len(matching) != 1:
        raise ValueError(f"expected exactly one PAW POTCAR dataset for {symbol!r}")
    content = matching[0]
    if "PAW radial sets" not in content or "Non local Part" not in content:
        raise ValueError("POTCAR lacks PAW radial sets and nonlocal projector tables")
    nonlocal_text, radial_text = content.split("PAW radial sets", 1)
    header, *blocks = nonlocal_text.split("Non local Part")
    match = re.search(r"([0-9.]+)\s+[TF]\s*$", header)
    if not match:
        raise ValueError("POTCAR nonlocal projector cutoff not found")
    gmax = float(match.group(1))
    angular, qtables = [], []
    for block in blocks:
        parts = block.split("Reciprocal Space Part")
        descriptor = parts[0].split()
        if len(descriptor) < 3:
            raise ValueError("invalid POTCAR nonlocal projector header")
        l, n = map(int, descriptor[:2])
        if len(parts) != n + 1:
            raise ValueError(
                "POTCAR radial projector count disagrees with nonlocal header"
            )
        angular.extend([l] * n)
        for pair in parts[1:]:
            if "Real Space Part" not in pair:
                raise ValueError("POTCAR projector lacks real-space companion")
            reciprocal, _ = pair.split("Real Space Part", 1)
            q = _numbers(reciprocal)
            if len(q) != 100:
                raise ValueError("expected 100 reciprocal projector samples")
            qtables.append(q)
    first, *rest = radial_text.strip().splitlines()
    radial = _radial_sections("\n".join(rest))
    ngrid = int(first.split()[0])
    grid = np.asarray(radial["grid"][0])
    ae = np.asarray(radial["ae wavefunction"])
    pseudo = np.asarray(radial["pseudo wavefunction"])
    if (
        grid.shape != (ngrid,)
        or ae.shape != pseudo.shape
        or ae.shape != (len(angular), ngrid)
    ):
        raise ValueError("POTCAR radial partial-wave dimensions mismatch")
    # Angular harmonics are orthonormal. dS couples radial channels with
    # equal (l,m), and vanishes between different angular channels.
    labels = [(n, l, m) for n, l in enumerate(angular) for m in range(-l, l + 1)]
    dS = np.zeros((len(labels), len(labels)))
    for i, (ni, li, mi) in enumerate(labels):
        for j, (nj, lj, mj) in enumerate(labels[: i + 1]):
            if li == lj and mi == mj:
                dS[i, j] = dS[j, i] = simpson(
                    ae[ni] * ae[nj] - pseudo[ni] * pseudo[nj],
                    x=grid,
                )
    return VaspPawDataset(
        symbol,
        np.asarray(angular),
        np.arange(100) * (gmax / 100),
        np.asarray(qtables),
        grid,
        ae,
        pseudo,
        dS,
    )
