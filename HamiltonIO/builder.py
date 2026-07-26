"""User-facing builders for small tight-binding Hamiltonians."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping

import numpy as np
from ase import Atoms

from HamiltonIO.lcao_hamiltonian import LCAOHamiltonian

ORBITAL_LM: dict[str, tuple[int, int]] = {
    "s": (0, 0),
    "px": (1, -1),
    "py": (1, 0),
    "pz": (1, 1),
    "dz2": (2, 0),
    "d3z2-r2": (2, 0),
    "d3z^2-r^2": (2, 0),
    "dx2-y2": (2, 1),
    "dx^2-y^2": (2, 1),
    "dxy": (2, 2),
    "dyz": (2, 3),
    "dxz": (2, 4),
}


@dataclass(frozen=True)
class OrbitalSpec:
    """Minimal orbital metadata used by the builder and downstream consumers."""

    iatom: int
    label: str
    l: int
    m: int

    @property
    def atom(self) -> int:
        return self.iatom

    @property
    def sym(self) -> str:
        return self.label


@dataclass(frozen=True)
class _Selector:
    atom: int | None = None
    species: str | None = None
    orbital: str | None = None
    shell: int | None = None
    index: int | None = None


@dataclass(frozen=True)
class _OnsiteTerm:
    selector: _Selector
    value: complex
    mode: str = "set"


@dataclass(frozen=True)
class _HoppingTerm:
    source: _Selector
    target: _Selector
    R: tuple[int, int, int]
    value: complex
    hermitian: bool = True


@dataclass
class TightBindingModelData:
    """Result returned by :class:`TightBindingBuilder`."""

    hamiltonian: LCAOHamiltonian
    orbitals: tuple[OrbitalSpec, ...]
    symbols: tuple[str, ...]
    orbital_index: dict[tuple[int, str], int]

    def select_orbitals(self, selector: Any = None, **kwargs: Any) -> np.ndarray:
        norm = _normalize_selector(selector, **kwargs)
        return np.asarray(
            _resolve_selector(norm, self.orbitals, self.symbols), dtype=int
        )

    def shell_groups(self) -> dict[tuple[int, int], np.ndarray]:
        groups: dict[tuple[int, int], list[int]] = {}
        for idx, orb in enumerate(self.orbitals):
            groups.setdefault((orb.iatom, orb.l), []).append(idx)
        return {key: np.asarray(value, dtype=int) for key, value in groups.items()}

    def hubbard_dict(
        self, spec: Mapping[Any, Mapping[str, float]]
    ) -> dict[str, dict[str, float]]:
        """Convert compact ``(species, shell)`` keys to TBUpy-style dictionaries."""

        result: dict[str, dict[str, float]] = {}
        for key, params in spec.items():
            if isinstance(key, tuple) and len(key) == 2:
                species, shell = key
                entry = dict(params)
                entry["L"] = int(shell)
                result[str(species)] = entry
            elif isinstance(key, str):
                result[key] = dict(params)
            else:
                raise ValueError(
                    "hubbard_dict keys must be species strings or (species, shell) tuples"
                )
        return result

    def initial_density(
        self, nel: float | None = None, moment: float = 0.0
    ) -> np.ndarray:
        """Build a simple diagonal collinear density for downstream SCF starts."""

        nspin = int(getattr(self.hamiltonian, "nspin", 1))
        if nspin not in (1, 2):
            raise ValueError("initial_density supports nspin=1 or nspin=2")
        norb = len(self.orbitals)
        nel_value = getattr(self.hamiltonian, "nel", None) if nel is None else nel
        if nel_value is None:
            raise ValueError(
                "initial_density requires nel or a Hamiltonian built with nel"
            )
        nel_eff = float(nel_value)
        if nspin == 1:
            rho = np.zeros((1, norb, norb), dtype=complex)
            rho[0, np.arange(norb), np.arange(norb)] = nel_eff / norb
            return rho
        rho = np.zeros((2, norb, norb), dtype=complex)
        base = nel_eff / (2.0 * norb)
        rho[0, np.arange(norb), np.arange(norb)] = base + 0.5 * moment
        rho[1, np.arange(norb), np.arange(norb)] = base - 0.5 * moment
        return rho


class TightBindingBuilder:
    """Fluent builder for orthogonal tight-binding ``LCAOHamiltonian`` objects."""

    def __init__(self, atoms: Atoms):
        self.atoms = atoms.copy()
        self._orbital_labels: dict[int, list[str]] = {}
        self._onsite_terms: list[_OnsiteTerm] = []
        self._hopping_terms: list[_HoppingTerm] = []

    @classmethod
    def from_data(
        cls,
        *,
        cell: Any,
        positions: Any,
        symbols: list[str] | tuple[str, ...],
        pbc: bool | tuple[bool, bool, bool] = True,
        scaled_positions: bool = False,
    ) -> "TightBindingBuilder":
        if scaled_positions:
            atoms = Atoms(
                symbols=list(symbols), scaled_positions=positions, cell=cell, pbc=pbc
            )
        else:
            atoms = Atoms(
                symbols=list(symbols), positions=positions, cell=cell, pbc=pbc
            )
        return cls(atoms)

    @classmethod
    def from_atoms(cls, atoms: Atoms) -> "TightBindingBuilder":
        return cls(atoms)

    def add_orbitals(
        self, *, atom: int, labels: list[str] | tuple[str, ...]
    ) -> "TightBindingBuilder":
        atom = int(atom)
        if atom < 0 or atom >= len(self.atoms):
            raise ValueError(f"atom index {atom} is out of range")
        labels = [str(label) for label in labels]
        if len(set(labels)) != len(labels):
            raise ValueError(f"duplicate orbital labels for atom {atom}: {labels}")
        for label in labels:
            if label not in ORBITAL_LM:
                raise ValueError(f"unsupported orbital label {label!r}")
        if atom in self._orbital_labels:
            raise ValueError(f"orbitals for atom {atom} have already been assigned")
        self._orbital_labels[atom] = labels
        return self

    def set_onsite(
        self,
        selector: Any = None,
        *,
        value: complex,
        mode: str = "set",
        **kwargs: Any,
    ) -> "TightBindingBuilder":
        if mode not in ("set", "add"):
            raise ValueError("onsite mode must be 'set' or 'add'")
        self._onsite_terms.append(
            _OnsiteTerm(_normalize_selector(selector, **kwargs), complex(value), mode)
        )
        return self

    def add_hopping(
        self,
        source: Any,
        target: Any,
        *,
        R: Any,
        value: complex,
        hermitian: bool = True,
    ) -> "TightBindingBuilder":
        self._hopping_terms.append(
            _HoppingTerm(
                _normalize_selector(source),
                _normalize_selector(target),
                _normalize_R(R),
                complex(value),
                bool(hermitian),
            )
        )
        return self

    def build(
        self, *, nspin: int = 1, nel: float | None = None
    ) -> TightBindingModelData:
        if nspin not in (1, 2):
            raise ValueError(f"nspin must be 1 or 2, got {nspin}")
        orbitals = self._build_orbitals()
        if not orbitals:
            raise ValueError("at least one orbital must be assigned")
        symbols = tuple(self.atoms.get_chemical_symbols())
        orbital_index = {
            (orb.iatom, orb.label): idx for idx, orb in enumerate(orbitals)
        }

        Rlist = self._build_Rlist()
        R_index = {R: idx for idx, R in enumerate(Rlist)}
        norb = len(orbitals)
        nbasis = norb * nspin
        HR = np.zeros((len(Rlist), nbasis, nbasis), dtype=complex)
        SR = np.zeros_like(HR)
        SR[R_index[(0, 0, 0)]] = np.eye(nbasis, dtype=complex)

        onsite_values: dict[int, complex] = {}
        for term in self._onsite_terms:
            for iorb in _resolve_selector(term.selector, orbitals, symbols):
                if term.mode == "set":
                    if iorb in onsite_values and not np.isclose(
                        onsite_values[iorb], term.value
                    ):
                        raise ValueError(
                            f"conflicting onsite values for orbital index {iorb}"
                        )
                    onsite_values[iorb] = term.value
                else:
                    onsite_values[iorb] = onsite_values.get(iorb, 0.0) + term.value
        for iorb, value in onsite_values.items():
            for spin in range(nspin):
                idx = iorb * nspin + spin
                HR[R_index[(0, 0, 0)], idx, idx] = value

        for term in self._hopping_terms:
            source = _resolve_one(term.source, orbitals, symbols, "hopping source")
            target = _resolve_one(term.target, orbitals, symbols, "hopping target")
            self._insert_hopping(HR, R_index, nspin, source, target, term.R, term.value)
            if term.hermitian:
                self._insert_hopping(
                    HR,
                    R_index,
                    nspin,
                    target,
                    source,
                    _neg_R(term.R),
                    np.conjugate(term.value),
                )

        _validate_hermitian(HR, Rlist)
        hio_orbs = _spin_duplicate_orbitals(orbitals, nspin)
        ham = LCAOHamiltonian(
            HR=HR,
            SR=SR,
            Rlist=np.asarray(Rlist, dtype=int),
            nbasis=nbasis,
            orbs=hio_orbs,
            atoms=self.atoms.copy(),
            nspin=nspin,
            nel=nel,
        )
        return TightBindingModelData(
            hamiltonian=ham,
            orbitals=tuple(orbitals),
            symbols=symbols,
            orbital_index=orbital_index,
        )

    def _build_orbitals(self) -> list[OrbitalSpec]:
        orbitals: list[OrbitalSpec] = []
        for iatom in range(len(self.atoms)):
            for label in self._orbital_labels.get(iatom, []):
                l_val, m_val = ORBITAL_LM[label]
                orbitals.append(OrbitalSpec(iatom=iatom, label=label, l=l_val, m=m_val))
        return orbitals

    def _build_Rlist(self) -> list[tuple[int, int, int]]:
        Rlist = [(0, 0, 0)]
        for term in self._hopping_terms:
            for R in (term.R, _neg_R(term.R) if term.hermitian else None):
                if R is not None and R not in Rlist:
                    Rlist.append(R)
        return Rlist

    @staticmethod
    def _insert_hopping(
        HR: np.ndarray,
        R_index: dict[tuple[int, int, int], int],
        nspin: int,
        source: int,
        target: int,
        R: tuple[int, int, int],
        value: complex,
    ) -> None:
        iR = R_index[R]
        for spin in range(nspin):
            HR[iR, source * nspin + spin, target * nspin + spin] += value


def _normalize_selector(selector: Any = None, **kwargs: Any) -> _Selector:
    if selector is not None and kwargs:
        raise ValueError(
            "selector positional form cannot be combined with selector keywords"
        )
    if isinstance(selector, _Selector):
        return selector
    if isinstance(selector, tuple) and len(selector) == 2:
        return _Selector(atom=int(selector[0]), orbital=str(selector[1]))
    if isinstance(selector, int):
        return _Selector(index=int(selector))
    if selector is not None:
        raise ValueError(f"unsupported selector form {selector!r}")
    allowed = {"atom", "species", "orbital", "shell", "index"}
    unknown = set(kwargs) - allowed
    if unknown:
        raise ValueError(f"unknown selector keys: {sorted(unknown)}")
    return _Selector(
        atom=None if kwargs.get("atom") is None else int(kwargs["atom"]),
        species=None if kwargs.get("species") is None else str(kwargs["species"]),
        orbital=None if kwargs.get("orbital") is None else str(kwargs["orbital"]),
        shell=None if kwargs.get("shell") is None else int(kwargs["shell"]),
        index=None if kwargs.get("index") is None else int(kwargs["index"]),
    )


def _resolve_selector(
    selector: _Selector,
    orbitals: list[OrbitalSpec] | tuple[OrbitalSpec, ...],
    symbols: tuple[str, ...],
) -> list[int]:
    if selector.index is not None:
        if selector.index < 0 or selector.index >= len(orbitals):
            raise ValueError(f"orbital index {selector.index} is out of range")
        matches = [selector.index]
    else:
        matches = []
        for idx, orb in enumerate(orbitals):
            if selector.atom is not None and orb.iatom != selector.atom:
                continue
            if selector.species is not None and symbols[orb.iatom] != selector.species:
                continue
            if selector.orbital is not None and orb.label != selector.orbital:
                continue
            if selector.shell is not None and orb.l != selector.shell:
                continue
            matches.append(idx)
    if not matches:
        raise ValueError(f"selector {selector} matched no orbitals")
    return matches


def _resolve_one(
    selector: _Selector,
    orbitals: list[OrbitalSpec],
    symbols: tuple[str, ...],
    label: str,
) -> int:
    matches = _resolve_selector(selector, orbitals, symbols)
    if len(matches) != 1:
        raise ValueError(
            f"{label} selector {selector} matched {len(matches)} orbitals; expected exactly one"
        )
    return matches[0]


def _normalize_R(R: Any) -> tuple[int, int, int]:
    arr = np.asarray(R)
    if arr.shape != (3,):
        raise ValueError("lattice translation R must have shape (3,)")
    if not np.allclose(arr, np.rint(arr)):
        raise ValueError("lattice translation R must contain integers")
    return tuple(int(x) for x in arr)


def _neg_R(R: tuple[int, int, int]) -> tuple[int, int, int]:
    return (-R[0], -R[1], -R[2])


def _spin_duplicate_orbitals(
    orbitals: list[OrbitalSpec], nspin: int
) -> list[OrbitalSpec]:
    if nspin == 1:
        return list(orbitals)
    duplicated: list[OrbitalSpec] = []
    for orb in orbitals:
        duplicated.extend([orb, orb])
    return duplicated


def _validate_hermitian(HR: np.ndarray, Rlist: list[tuple[int, int, int]]) -> None:
    R_index = {R: idx for idx, R in enumerate(Rlist)}
    for R, iR in R_index.items():
        neg = _neg_R(R)
        if neg not in R_index:
            raise ValueError(f"missing Hermitian counterpart R={neg}")
        diff = HR[iR] - HR[R_index[neg]].conjugate().T
        if not np.allclose(diff, 0.0):
            raise ValueError(f"Hamiltonian is not Hermitian for R={R}")


__all__ = ["ORBITAL_LM", "OrbitalSpec", "TightBindingBuilder", "TightBindingModelData"]
