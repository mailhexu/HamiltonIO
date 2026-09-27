#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The abacus wrapper
"""

import os
import re
from pathlib import Path

import numpy as np

from HamiltonIO.abacus.abacus_api import XR_matrix, read_csr_new_format
from HamiltonIO.abacus.orbital_api import parse_abacus_orbital
from HamiltonIO.abacus.stru_api import read_abacus
from HamiltonIO.lcao_hamiltonian import LCAOHamiltonian
from HamiltonIO.mathutils.pauli import chargepart, spinpart
from HamiltonIO.model.kR_convert import R_to_onek
from HamiltonIO.utils import symbol_number_list

SPLIT_SR_ATOL = 1e-12

_ORBLETTERS = {"S": 0, "P": 1, "D": 2, "F": 3, "G": 4, "H": 5, "I": 6}


def _parse_orbital_zeta_counts(orb_file):
    """Return ``[(l, nz), ...]`` parsed from an ABACUS ``.orb`` header."""
    import re

    counts = {}
    with open(orb_file) as handle:
        for line in handle:
            if line.startswith("SUMMARY"):
                break
            match = re.match(r"\s*Number of ([SPDFGI])orbital\s*-->\s*(\d+)", line)
            if match:
                counts[_ORBLETTERS[match.group(1)]] = int(match.group(2))
    if not counts:
        raise ValueError(f"no orbital zeta counts found in {orb_file}")
    return sorted(counts.items())


def _align_csr_blocks(
    reference_basis,
    reference_rlist,
    basis,
    rlist,
    matrix,
    label,
    allow_missing_zero=False,
):
    if basis != reference_basis:
        raise ValueError(f"{label} basis {basis} does not match {reference_basis}")
    reference = [tuple(r) for r in reference_rlist]
    candidate = [tuple(r) for r in rlist]
    if len(reference) != len(set(reference)) or len(candidate) != len(set(candidate)):
        raise ValueError(f"{label} contains duplicate R vectors")
    reference_set = set(reference)
    candidate_set = set(candidate)
    if allow_missing_zero:
        if not candidate_set <= reference_set:
            raise ValueError(
                f"{label} R sets contain translations outside the reference"
            )
        aligned = np.zeros((len(reference),) + matrix.shape[1:], dtype=matrix.dtype)
        positions = {r: i for i, r in enumerate(candidate)}
        for index, r in enumerate(reference):
            if r in positions:
                aligned[index] = matrix[positions[r]]
        return aligned
    if reference_set != candidate_set:
        raise ValueError(f"{label} R sets do not match the reference")
    positions = {r: i for i, r in enumerate(candidate)}
    return matrix[[positions[r] for r in reference]]


class AbacusWrapper(LCAOHamiltonian):
    def __init__(
        self,
        HR,
        SR,
        Rlist,
        nbasis=None,
        atoms=None,
        nspin=1,
        basis=None,
        HR_soc=None,
        HR_nosoc=None,
        nel=None,
    ):
        super().__init__(
            HR,
            SR,
            Rlist,
            nbasis=nbasis,
            atoms=atoms,
            nspin=nspin,
            orbs=None,
            HR_soc=HR_soc,
            HR_nosoc=HR_nosoc,
            nel=nel,
        )
        self._name = "ABACUS"


class AbacusUnfoldModel:
    """Unfolding-facing view of an ABACUS LCAO Hamiltonian.

    Exposes the minimal model interface consumed by the band-unfolding
    package (see ``unfolding.lcao_unfolder.HamiltonIOModel``):

    - ``atoms``: the ASE Atoms object of the supercell (Angstrom),
    - ``HR``/``SR``: dicts ``{R_tuple: (norb, norb) ndarray}`` keyed by
      integer supercell-lattice translations (HR in eV),
    - ``hs_and_eigen(k)``: ``(H, S)`` at one supercell fractional
      k-point, built with HamiltonIO convention 2
      (``H(k) = sum_R H(R) exp(2 pi i k.R)``),
    - ``atom_orbital_counts``: per-atom orbital counts for building
      :class:`~unfolding.mapping.RelabelMap` orbital counts.

    A collinear (nspin=2) run must select one spin channel before
    construction; noncollinear (nspin=4) tables are passed through as
    stored (spinor-block-doubled orbital ordering).
    """

    def __init__(self, atoms, HR, SR, Rlist, basis=None, efermi=None):
        self.atoms = atoms
        keys = [tuple(int(v) for v in R) for R in Rlist]
        if len(keys) != len(set(keys)):
            raise ValueError("duplicate R vectors in Rlist")
        self._keys = keys
        self._HR = np.asarray(HR)
        self._SR = np.asarray(SR)
        self.basis = list(basis) if basis is not None else []
        self.efermi = efermi

    @property
    def HR(self):
        return {R: self._HR[i] for i, R in enumerate(self._keys)}

    @property
    def SR(self):
        return {R: self._SR[i] for i, R in enumerate(self._keys)}

    @property
    def Rlist(self):
        return list(self._keys)

    @property
    def atom_orbital_counts(self):
        """Per-supercell-atom orbital counts (atom-list order)."""
        counts = np.zeros(len(self.atoms), dtype=int)
        for orb in self.basis:
            counts[orb.iatom] += 1
        return [int(c) for c in counts]

    def hs_and_eigen(self, k):
        """Return ``(H, S)`` at a single supercell fractional k-point."""
        k = np.asarray(k, dtype=float).reshape(3)
        H = R_to_onek(k, self._keys, self._HR)
        S = R_to_onek(k, self._keys, self._SR)
        return H, S


class AbacusParser:
    def __init__(self, spin=None, outpath=None, binary=False):
        if outpath is None:
            raise ValueError("outpath is required")
        self.outpath = Path(outpath)
        if spin is None:
            self.spin = self.read_spin()
        else:
            self.spin = spin
        self.binary = binary
        # read the information
        self.read_atoms()
        if (Path(self.outpath) / "running_scf.log").exists():
            self.efermi = self.read_efermi()
            self.nel = self.read_nel()
        else:
            self.efermi = None
            self.nel = None
        if self.spin in ["non-polarized", "collinear"]:
            self.basis = self.read_basis()
        elif self.spin == "noncollinear":
            self.basis = self.read_basis(nspin=2)

    def read_spin(self):
        log = Path(self.outpath) / "running_scf.log"
        if not log.exists():
            raise FileNotFoundError(
                f"{log} not found; pass spin=('non-polarized'|'collinear'|"
                "'noncollinear') explicitly to work from the CSR tables alone"
            )
        with open(str(log)) as myfile:
            for line in myfile:
                if line.strip().startswith("nspin"):
                    nspin = int(line.strip().split()[-1])
                    if nspin == 1:
                        return "non-polarized"
                    elif nspin == 2:
                        return "collinear"
                    elif nspin == 4:
                        return "noncollinear"
                    else:
                        raise ValueError("nspin should be either 1 or 4.")

    def read_atoms(self):
        path1 = str(Path(self.outpath) / "../STRU")
        path2 = str(Path(self.outpath) / "../Stru")
        if os.path.exists(path1):
            self.atoms = read_abacus(path1, verbose=True)
        elif os.path.exists(path2):
            self.atoms = read_abacus(path2, verbose=True)
        else:
            raise Exception("The STRU or Stru file cannot be found.")
        return self.atoms

    def read_basis(self, nspin=1):
        fname = Path(self.outpath) / "Orbital"
        if fname.exists():
            self.basis = parse_abacus_orbital(fname, nspin=nspin)
        else:
            self.basis = self.read_basis_from_stru(nspin=nspin)
        return self.basis

    def _resolve_orb_file(self, orb_file):
        """Locate a STRU-referenced ``.orb`` file (paths are relative)."""
        p = Path(orb_file)
        if p.is_absolute() and p.exists():
            return p
        run_root = Path(self.outpath).parent
        candidates = [run_root / p]
        echo = Path(self.outpath) / "INPUT"
        if echo.exists():
            for line in echo.read_text().splitlines():
                parts = line.split()
                if len(parts) >= 2 and parts[0] == "orbital_dir":
                    directory = (
                        line.split(None, 1)[1].split("#")[0].strip().rstrip(" =")
                    )
                    candidates.append(run_root / directory / p)
                    candidates.append(Path(directory) / p)
                    break
        for candidate in candidates:
            if candidate.exists():
                return candidate
        raise FileNotFoundError(f"orbital file {orb_file!r} not found near {run_root}")

    def read_basis_from_stru(self, nspin=1):
        """Synthesize the basis list from the STRU orbital files.

        ABACUS <=3.10 does not write ``OUT.${suffix}/Orbital``.  The
        per-species zeta counts are parsed from the ``.orb`` files
        referenced by STRU and expanded in ABACUS orbital order
        (atom-major, then l, zeta, m), so orbital indexing matches the
        Hamiltonian tables.
        """
        orb_files = self.atoms.info.get("basis")
        if not orb_files:
            raise FileNotFoundError(
                f"neither {Path(self.outpath) / 'Orbital'} nor STRU orbital "
                "information is available"
            )
        from HamiltonIO.abacus.orbital_api import AbacusOrbital

        orb_files = {sym: self._resolve_orb_file(f) for sym, f in orb_files.items()}

        nz_by_l = {}
        for symbol, orb_file in orb_files.items():
            nz_by_l[symbol] = _parse_orbital_zeta_counts(orb_file)
        orbs = []
        for symbol in orb_files:  # ATOMIC_SPECIES block order
            for iatom, sym in enumerate(self.atoms.get_chemical_symbols()):
                if sym != symbol:
                    continue
                for l, nz in nz_by_l[symbol]:
                    for z in range(1, nz + 1):
                        for m in range(-l, l + 1):
                            for spin in range(nspin):
                                orbs.append(
                                    AbacusOrbital(
                                        iatom=iatom,
                                        sym=f"l{l}m{m}Z{z}",
                                        spin=spin,
                                        element=symbol,
                                        l=l,
                                        m=m,
                                        z=z,
                                    )
                                )
        return orbs

    def read_HSR_collinear(self, binary=None):
        p = Path(self.outpath)
        binary = self.binary if binary is None else binary
        filenames, is_new = self._select_csr_set(
            p,
            ("hrs1_nao.csr", "hrs2_nao.csr", "srs1_nao.csr"),
            (
                "data-HR-sparse_SPIN0.csr",
                "data-HR-sparse_SPIN1.csr",
                "data-SR-sparse_SPIN0.csr",
            ),
            binary,
        )
        HR_filename = filenames[:2]
        SR_filename = filenames[2]
        if is_new:
            if binary:
                raise ValueError(
                    "binary=True is unsupported for new-format text CSR files"
                )
            nbasis, Rlist, HR_up = read_csr_new_format(HR_filename[0], is_complex=False)
            nbasis_dn, Rlist_dn, HR_dn = read_csr_new_format(
                HR_filename[1], is_complex=False
            )
            nbasis_sr, Rlist_sr, SR = read_csr_new_format(SR_filename, is_complex=False)
            HR_dn = _align_csr_blocks(
                nbasis, Rlist, nbasis_dn, Rlist_dn, HR_dn, "HR down"
            )
            SR = _align_csr_blocks(
                nbasis, Rlist, nbasis_sr, Rlist_sr, SR, "SR", allow_missing_zero=True
            )
            return nbasis, Rlist, HR_up * 13.605698066, HR_dn * 13.605698066, SR
        nbasis, Rlist, HR_up = self._read_legacy_csr(HR_filename[0], 2, binary)
        nbasis_dn, Rlist_dn, HR_dn = self._read_legacy_csr(HR_filename[1], 2, binary)
        nbasis_sr, Rlist_sr, SR = self._read_legacy_csr(SR_filename, 2, binary)
        HR_dn = _align_csr_blocks(nbasis, Rlist, nbasis_dn, Rlist_dn, HR_dn, "HR down")
        SR = _align_csr_blocks(
            nbasis, Rlist, nbasis_sr, Rlist_sr, SR, "SR", allow_missing_zero=True
        )
        HR_up *= 13.605698066
        HR_dn *= 13.605698066
        return nbasis, Rlist, HR_up, HR_dn, SR

    def Read_HSR_noncollinear(self, binary=None):
        p = Path(self.outpath)
        binary = self.binary if binary is None else binary
        filenames, is_new = self._select_csr_set(
            p,
            ("hrs1_nao.csr", "srs1_nao.csr"),
            ("data-HR-sparse_SPIN0.csr", "data-SR-sparse_SPIN0.csr"),
            binary,
        )
        HR_filename, SR_filename = filenames
        if is_new:
            if binary:
                raise ValueError(
                    "binary=True is unsupported for new-format text CSR files"
                )
            _Ry_to_eV = 13.605698066
            nbasis, Rlist, HR = read_csr_new_format(HR_filename, is_complex=True)
            nbasis_sr, Rlist_s, SR = read_csr_new_format(SR_filename, is_complex=True)
            HR *= _Ry_to_eV
            SR = _align_csr_blocks(
                nbasis, Rlist, nbasis_sr, Rlist_s, SR, "SR", allow_missing_zero=True
            )
        else:
            nbasis, Rlist, HR = self._read_legacy_csr(HR_filename, 4, binary)
            nbasis_sr, Rlist_sr, SR = self._read_legacy_csr(SR_filename, 4, binary)
            HR *= 13.605698066
            SR = _align_csr_blocks(
                nbasis, Rlist, nbasis_sr, Rlist_sr, SR, "SR", allow_missing_zero=True
            )
        return nbasis, Rlist, HR, SR

    def read_HSR_nonpolarized(self, binary=None):
        """Read R-space HR/SR for a non-polarized (nspin=1) run.

        Returns ``(nbasis, Rlist, HR, SR)`` with HR in eV.
        """
        p = Path(self.outpath)
        binary = self.binary if binary is None else binary
        filenames, is_new = self._select_csr_set(
            p,
            ("hrs1_nao.csr", "srs1_nao.csr"),
            ("data-HR-sparse_SPIN0.csr", "data-SR-sparse_SPIN0.csr"),
            binary,
        )
        HR_filename, SR_filename = filenames
        if is_new:
            if binary:
                raise ValueError(
                    "binary=True is unsupported for new-format text CSR files"
                )
            nbasis, Rlist, HR = read_csr_new_format(HR_filename, is_complex=False)
            nbasis_sr, Rlist_sr, SR = read_csr_new_format(SR_filename, is_complex=False)
            SR = _align_csr_blocks(
                nbasis, Rlist, nbasis_sr, Rlist_sr, SR, "SR", allow_missing_zero=True
            )
        else:
            nbasis, Rlist, HR = self._read_legacy_csr(HR_filename, 2, binary)
            nbasis_sr, Rlist_sr, SR = self._read_legacy_csr(SR_filename, 2, binary)
            SR = _align_csr_blocks(
                nbasis, Rlist, nbasis_sr, Rlist_sr, SR, "SR", allow_missing_zero=True
            )
        HR = np.asarray(HR) * 13.605698066
        return nbasis, Rlist, HR, SR

    @staticmethod
    def _read_legacy_csr(filename, nspin, binary):
        matrix = XR_matrix(nspin, filename)
        if binary:
            matrix.read_file_binary()
        else:
            matrix.read_file()
        return matrix.basis_num, matrix.R_direct_coor, matrix.XR

    @staticmethod
    def _find_csr(path, new_name, old_name):
        # ABACUS v3.11+ renamed CSR files (hrs1_nao.csr etc.) vs v3.10-LTS (data-HR-sparse_SPIN0.csr)
        new_p = path / new_name
        old_p = path / old_name
        if new_p.exists():
            return str(new_p)
        if old_p.exists():
            return str(old_p)
        raise FileNotFoundError(f"Neither {new_name} nor {old_name} found in {path}")

    @staticmethod
    def _select_csr_set(path, new_names, legacy_names, binary):
        new_paths = tuple(str(path / name) for name in new_names)
        legacy_paths = tuple(str(path / name) for name in legacy_names)
        new_present = [Path(filename).exists() for filename in new_paths]
        legacy_present = [Path(filename).exists() for filename in legacy_paths]
        if binary:
            if all(legacy_present):
                return legacy_paths, False
            if any(new_present):
                raise ValueError(
                    "binary=True is unsupported for new-format text CSR files"
                )
            raise FileNotFoundError(f"legacy binary CSR files are incomplete in {path}")
        if all(new_present):
            return new_paths, True
        if all(legacy_present) and not any(new_present):
            return legacy_paths, False
        if any(new_present) or any(legacy_present):
            raise ValueError(f"mixed new/legacy CSR file set in {path}")
        raise FileNotFoundError(f"No CSR file set found in {path}")

    def get_models(self):
        if self.spin == "non-polarized":
            nbasis, Rlist, HR, SR = self.read_HSR_nonpolarized()
            model = AbacusWrapper(
                HR=HR, SR=SR, Rlist=Rlist, nbasis=nbasis, nspin=1, atoms=self.atoms
            )
            model.efermi = self.efermi
            model.basis = self.get_basis()
            model.orbs = self.basis
            model.atoms = self.atoms
            return model
        elif self.spin == "collinear":
            nbasis, Rlist, HR_up, HR_dn, SR = self.read_HSR_collinear()
            model_up = AbacusWrapper(
                HR=HR_up, SR=SR, Rlist=Rlist, nbasis=nbasis, nspin=1, atoms=self.atoms
            )
            model_dn = AbacusWrapper(
                HR=HR_dn, SR=SR, Rlist=Rlist, nbasis=nbasis, nspin=1, atoms=self.atoms
            )
            model_up.efermi = self.efermi
            model_dn.efermi = self.efermi
            model_up.basis, model_dn.basis = self.get_basis()
            model_up.orbs = self.basis
            model_dn.orbs = self.basis
            model_up.atoms = self.atoms
            model_dn.atoms = self.atoms
            return model_up, model_dn
        elif self.spin == "noncollinear":
            nbasis, Rlist, HR, SR = self.Read_HSR_noncollinear()
            model = AbacusWrapper(
                HR=HR, SR=SR, Rlist=Rlist, nbasis=nbasis, nspin=2, atoms=self.atoms
            )
            model.efermi = self.efermi
            model.basis = self.get_basis()
            model.orbs = self.basis
            model.atoms = self.atoms
            return model

    def read_efermi(self):
        fname = str(Path(self.outpath) / "running_scf.log")
        efermi = None
        number = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
        with open(fname, "r") as myfile:
            for line in myfile:
                assignment = re.match(r"\s*(?:EFERMI|E_Fermi)\s*=", line)
                if assignment:
                    match = re.fullmatch(
                        rf"\s*(?:EFERMI|E_Fermi)\s*=\s*({number})(?:\s+eV)?\s*",
                        line,
                    )
                    if not match:
                        raise ValueError(f"EFERMI in {fname} must have a numeric value")
                    efermi = float(match.group(1))
                    continue
                table_start = re.match(r"\s*E_Fermi\b(.*)", line)
                if table_start and re.match(r"\s*[-+.\d]", table_start.group(1)):
                    match = re.fullmatch(
                        rf"\s*E_Fermi\s+({number})\s+({number})\s*", line
                    )
                    if not match:
                        raise ValueError(
                            f"EFERMI in {fname} must have numeric Ry and eV columns"
                        )
                    efermi = float(match.group(2))
        if efermi is None:
            raise ValueError(f"EFERMI not found in the {str(fname)}  file.")
        return efermi

    def read_nel(self):
        """
        Reading the number of electrons from the scf log file.
        """
        fname = str(Path(self.outpath) / "running_scf.log")
        nel = None
        with open(fname, "r") as myfile:
            for line in myfile:
                if "number of electrons" in line:
                    nel = float(line.split()[-1])
        if nel is None:
            raise ValueError(f"number of electron not found in the {str(fname)}  file.")
        return nel

    def get_unfold_model(self, spin=None):
        """Return the :class:`AbacusUnfoldModel` for band unfolding.

        Parameters
        ----------
        spin : {"up", "down", None}
            Channel selection for collinear (nspin=2) runs. ``None``
            defaults to the up channel. Ignored for non-polarized and
            noncollinear runs.
        """
        if self.spin == "non-polarized":
            nbasis, Rlist, HR, SR = self.read_HSR_nonpolarized()
            return AbacusUnfoldModel(
                self.atoms, HR, SR, Rlist, basis=self.basis, efermi=self.efermi
            )
        elif self.spin == "collinear":
            nbasis, Rlist, HR_up, HR_dn, SR = self.read_HSR_collinear()
            if spin in (None, "up"):
                HR = HR_up
            elif spin == "down":
                HR = HR_dn
            else:
                raise ValueError(f"unknown spin channel: {spin!r}")
            return AbacusUnfoldModel(
                self.atoms, HR, SR, Rlist, basis=self.basis, efermi=self.efermi
            )
        elif self.spin == "noncollinear":
            nbasis, Rlist, HR, SR = self.Read_HSR_noncollinear()
            return AbacusUnfoldModel(
                self.atoms, HR, SR, Rlist, basis=self.basis, efermi=self.efermi
            )
        raise ValueError(f"unknown spin: {self.spin!r}")

    def get_basis(self):
        slist = symbol_number_list(self.atoms)
        if self.spin == "non-polarized":
            return [(slist[b.iatom], b.sym) for b in self.basis]
        elif self.spin == "collinear":
            basis_up = []
            basis_dn = []
            for b in self.basis:
                basis_up.append((slist[b.iatom], b.sym, "up"))
                basis_dn.append((slist[b.iatom], b.sym, "down"))
            return basis_up, basis_dn
        elif self.spin == "noncollinear":
            basis = []
            for b in self.basis:
                basis.append((slist[b.iatom], b.sym, b.spin_symbol))
            return basis


class AbacusSingleStepSOCParser:
    """Read ABACUS output with single-step SOC separation (out_mat_hs2_soc=1)."""

    def __init__(self, outpath=None, binary=False):
        if outpath is None:
            raise ValueError("outpath is required")
        self.outpath = Path(outpath)
        self.binary = binary
        self.parser = AbacusParser(outpath=self.outpath, binary=binary)

    def parse(self):
        nbasis, Rlist, HR, SR = self.parser.Read_HSR_noncollinear()
        p = Path(self.outpath)
        soc_file = p / "hrs1-soc_nao.csr"
        Ry_to_eV = 13.605698066
        if not soc_file.exists():
            raise FileNotFoundError(f"SOC Hamiltonian file not found: {soc_file}")
        nbasis_soc, Rlist_soc, HR_soc = read_csr_new_format(
            str(soc_file), is_complex=True
        )
        HR_soc = _align_csr_blocks(
            nbasis, Rlist, nbasis_soc, Rlist_soc, HR_soc, "SOC Hamiltonian"
        )
        HR_soc *= Ry_to_eV
        HR_soc = np.array([spinpart(HR) for HR in HR_soc])

        model = AbacusWrapper(
            HR=HR,
            SR=SR,
            Rlist=Rlist,
            nbasis=nbasis,
            nspin=2,
            atoms=self.parser.atoms,
            HR_soc=HR_soc,
            HR_nosoc=HR - HR_soc,
            nel=self.parser.nel,
        )
        model.efermi = self.parser.efermi
        model.basis = self.parser.get_basis()
        model.orbs = self.parser.basis
        model.atoms = self.parser.atoms
        return model


class AbacusSplitSOCParser:
    """
    Abacus parser with Hamiltonian split to SOC and non-SOC parts
    """

    def __init__(self, outpath_nosoc=None, outpath_soc=None, binary=False):
        if outpath_nosoc is None or outpath_soc is None:
            raise ValueError("both SOC output paths are required")
        self.outpath_nosoc = Path(outpath_nosoc)
        self.outpath_soc = Path(outpath_soc)
        self.binary = binary
        self.parser_nosoc = AbacusParser(outpath=self.outpath_nosoc, binary=binary)
        self.parser_soc = AbacusParser(outpath=self.outpath_soc, binary=binary)
        spin1 = self.parser_nosoc.read_spin()
        spin2 = self.parser_soc.read_spin()
        if spin1 != "noncollinear" or spin2 != "noncollinear":
            raise ValueError("Spin should be noncollinear")

    def parse(self):
        nbasis, Rlist, HR_nosoc, SR = self.parser_nosoc.Read_HSR_noncollinear()
        nbasis2, Rlist2, HR2, SR2 = self.parser_soc.Read_HSR_noncollinear()
        HR2 = _align_csr_blocks(nbasis, Rlist, nbasis2, Rlist2, HR2, "SOC Hamiltonian")
        SR2 = _align_csr_blocks(
            nbasis, Rlist, nbasis2, Rlist2, SR2, "SOC overlap", allow_missing_zero=True
        )
        if SR.size == 0 or SR2.size == 0:
            raise ValueError("SOC and non-SOC overlap domains must not be empty")
        if not np.isfinite(SR).all() or not np.isfinite(SR2).all():
            raise ValueError("SOC and non-SOC overlap values must be finite")
        max_difference = float(np.max(np.abs(SR2 - SR)))
        if not np.isfinite(max_difference):
            raise ValueError("SOC and non-SOC overlap difference must be finite")
        if max_difference > SPLIT_SR_ATOL:
            raise ValueError(
                f"SOC and non-SOC overlap differs: max difference {max_difference:.3e} exceeds {SPLIT_SR_ATOL:.1e}"
            )
        if self.parser_nosoc.get_basis() != self.parser_soc.get_basis():
            raise ValueError("SOC and non-SOC basis descriptors do not match")
        HR_soc = HR2 - HR_nosoc
        for iR, _ in enumerate(Rlist):
            spart, cpart = spinpart(HR_soc[iR]), chargepart(HR_soc[iR])
            HR_nosoc[iR] += cpart
            HR_soc[iR] = spart

        model = AbacusWrapper(
            HR=None,
            SR=SR,
            Rlist=Rlist,
            nbasis=nbasis,
            nspin=2,
            HR_soc=HR_soc,
            HR_nosoc=HR_nosoc,
            nel=self.parser_nosoc.nel,
        )
        model.efermi = self.parser_soc.efermi
        model.basis = self.parser_nosoc.basis
        model.orbs = self.parser_nosoc.basis
        model.atoms = self.parser_nosoc.atoms
        return model


def test_abacus_wrapper_collinear():
    outpath = "/Users/hexu/projects/TB2J_abacus/abacus-tb2j-master/abacus_example/case_Fe/1_no_soc/OUT.Fe"
    parser = AbacusParser(outpath=outpath, spin=None, binary=False)
    atoms = parser.read_atoms()
    # atoms=parser.read_atoms_out()
    # parser.read_HSR_collinear()
    model_up, model_dn = parser.get_models()
    H, S, E, V = model_up.HSE_k([0, 0, 0])
    # print(H.diagonal().real)
    # print(model_up.get_HR0().diagonal().real)
    print(parser.efermi)
    return model_up, model_dn, parser, atoms, H, S, E, V


# def test_abacus_wrapper_ncl():
#    outpath = "/Users/hexu/projects/TB2J_abacus/abacus-tb2j-master/abacus_example/case_Fe/2_soc/OUT.Fe"
#
#    parser = AbacusParser(outpath=outpath, spin=None, binary=False)
#    #atoms = parser.read_atoms()
#    #model = parser.get_models()
#    #H, S, E, V = model.HSE_k([0, 0, 0])
#    #print(parser.efermi)
#    retrun parser


if __name__ == "__main__":
    # test_abacus_wrapper()
    # test_abacus_wrapper_ncl()
    pass
