"""Pure-python reader for OpenMX ``*.scfout`` binary files.

The layout follows ``openmx3.9/source/read_scfout.c`` (SCFOUT version 3):
files are native-endian (little-endian on all architectures OpenMX
supports); the reader validates the header and refuses byte-swapped
files rather than guessing.

Storage is block-per-neighbour-pair: for atom ``ct`` (1-based in the
original C arrays) and neighbour index ``h``, the block
``Hks[spin][ct][h]`` has shape
``(Total_NumOrbs[ct], Total_NumOrbs[natn[ct][h]])`` where the column
atom ``natn[ct][h]`` sits in the periodic image indexed by
``ncn[ct][h]`` into the ``atv``/``atv_ijk`` image tables. Energies are
in Hartree, lengths in Bohr. Image translations are keyed by the
integer triples ``atv_ijk[Rn][1:4]``; the image set is the symmetric
cube ``{-CpyCell..CpyCell}**3`` (``truncation.c: Generation_ATV``).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

SCFOUT_VERSION = 3


@dataclass
class Scfout:
    """Parsed contents of one ``.scfout`` file.

    Attributes
    ----------
    natom, spin_p_switch, version, solver, tcpy_cell :
        Header scalars. ``spin_p_switch`` is 0 (unpolarized), 1
        (collinear) or 3 (non-collinear).
    total_numorbs : (natom,) int
        Number of local orbitals per atom.
    fnan : (natom,) int
        Number of first-neighbour entries per atom.
    natn, ncn : list of int arrays
        ``natn[ct][h]`` (0-based) global index of neighbour ``h`` of
        atom ``ct``; ``ncn[ct][h]`` its image index into ``atv_ijk``.
    atv_ijk : (nR, 3) int
        Integer image translations ``R`` (row 0 is ``(0, 0, 0)``).
    tv : (3, 3) float
        Cell vectors in Bohr (rows).
    positions : (natom, 3) float
        Cartesian atomic coordinates in Bohr.
    H : (nR, norb, norb) for spin 0, (nspin, nR, norb, norb) real for
        spin 1, and (nR, 2*norb, 2*norb) complex for spin 3 (assembled
        non-collinear blocks, orbital-major spinor order).
    S : (nR, norb, norb) float
        Stored overlap blocks (single-spin dimension).
    chem_p, valence_electrons, e_temp, total_spin_s : float
        Fermi level (Hartree), valence electrons, smearing (K), spin S.
    input_lines : list of str
        The echoed OpenMX input file stored at the tail of the file.
    """

    natom: int = 0
    spin_p_switch: int = 0
    version: int = SCFOUT_VERSION
    tcpy_cell: int = 0
    solver: int = 0
    total_numorbs: np.ndarray = None
    fnan: np.ndarray = None
    natn: list = field(default_factory=list)
    ncn: list = field(default_factory=list)
    atv_ijk: np.ndarray = None
    tv: np.ndarray = None
    positions: np.ndarray = None
    H: np.ndarray = None
    S: np.ndarray = None
    chem_p: float = 0.0
    e_temp: float = 0.0
    valence_electrons: float = 0.0
    total_spin_s: float = 0.0
    input_lines: list = field(default_factory=list)

    @property
    def norb(self) -> int:
        """Single-spin orbital count (spinor half-dimension for spin 3)."""
        return int(np.sum(self.total_numorbs))

    @property
    def Rlist(self) -> np.ndarray:
        """Integer image translations, one row per stored shell."""
        return np.asarray(self.atv_ijk)

    def symbols(self) -> list:
        """Element symbols parsed from the echoed input file."""
        return parse_symbols(self.input_lines, self.natom)


def parse_symbols(input_lines, natom) -> list:
    """Extract per-atom element symbols from echoed OpenMX input lines.

    Species labels are resolved through ``Definition.of.Atomic.Species``
    (``label  basis-spec  vps``): the element symbol is the leading
    alphabetic prefix of the basis spec (e.g. ``Si7.0-s2p2d1 -> Si``).
    """
    label_to_symbol = {}
    in_def = False
    in_atoms = False
    atoms = []
    for line in input_lines:
        low = line.lower()
        if "<definition.of.atomic.species" in low:
            in_def = True
            continue
        if in_def:
            if "definition.of.atomic.species>" in low:
                in_def = False
                continue
            tok = line.split()
            if len(tok) >= 2:
                label = tok[0]
                spec = tok[1]
                sym = spec[:2] if spec[:2].isalpha() else spec[:1]
                label_to_symbol[label] = sym
            continue
        if "<atoms.speciesandcoordinates" in low:
            in_atoms = True
            continue
        if in_atoms:
            if "atoms.speciesandcoordinates>" in low:
                in_atoms = False
                continue
            tok = line.split()
            if len(tok) >= 5 and tok[0].isdigit():
                label = tok[1]
                atoms.append(label_to_symbol.get(label, label[:2]))
    if len(atoms) != natom:
        raise ValueError(
            f"parsed {len(atoms)} species lines but scfout has {natom} atoms"
        )
    return atoms


def assemble_ncol(H_re, H_im):
    """Assemble the non-collinear ``(2n, 2n)`` block from OpenMX parts.

    ``H_re`` is ``(4, n, n)``: ``a = Hks[0]`` (Re H_aa), ``b = Hks[1]``
    (Re H_bb), ``c = Hks[2]`` (Re H_ab), ``d = Hks[3]`` (Im H_ab).
    ``H_im`` is ``(3, n, n)``: the ``iHks`` imaginary SOC+U parts for
    aa, bb, ab. Spinor basis is orbital-major (orb1_up, orb1_dn, ...),
    matching ``EigenValue_Problem.c``.
    """
    n = H_re.shape[-1]
    out = np.zeros((2 * n, 2 * n), dtype=complex)
    out[0::2, 0::2] = H_re[0] + 1j * H_im[0]
    out[1::2, 1::2] = H_re[1] + 1j * H_im[1]
    ab = H_re[2] + 1j * (H_re[3] + H_im[2])
    out[0::2, 1::2] = ab
    out[1::2, 0::2] = ab.conj()
    return out


def read_scfout(path) -> Scfout:
    """Read one OpenMX ``.scfout`` file."""
    with open(path, "rb") as fh:
        raw = fh.read()
    return parse_scfout_bytes(raw)


def parse_scfout_bytes(raw: bytes) -> Scfout:
    off = 0

    def take(dtype, count, shape=None):
        nonlocal off
        item = np.dtype(dtype).itemsize
        arr = np.frombuffer(raw, dtype=dtype, count=count, offset=off)
        off += item * count
        if shape is not None:
            return arr.reshape(shape)
        return arr

    # --- header -----------------------------------------------------------
    head = take("<i4", 6)
    atomnum = int(head[0])
    code = int(head[1])
    if not 0 <= code <= 15:  # read_scfout.c endianness guard
        raise ValueError(
            "unsupported .scfout header (endianness or version mismatch): "
            f"atomnum={atomnum}, code={code}"
        )
    spin_p_switch = code % 4
    version = code // 4
    if version != SCFOUT_VERSION:
        raise ValueError(
            f".scfout version {version} not supported (need {SCFOUT_VERSION})"
        )
    tcpy_cell = int(head[5])

    order_max = int(take("<i4", 1)[0])

    nR = tcpy_cell + 1
    take("<f8", 4 * nR, (nR, 4))  # atv (Bohr) -- keys come from atv_ijk
    atv_ijk = take("<i4", 4 * nR, (nR, 4))[:, 1:4].astype(int)

    Total_NumOrbs = take("<i4", atomnum)
    FNAN = take("<i4", atomnum)

    natn, ncn = [], []
    for a in range(atomnum):
        natn.append(take("<i4", int(FNAN[a]) + 1))
    for a in range(atomnum):
        ncn.append(take("<i4", int(FNAN[a]) + 1))

    tv = take("<f8", 12).reshape(3, 4)[:, 1 : 3 + 1]  # cell vectors (Bohr)
    take("<f8", 12).reshape(3, 4)  # rtv (reciprocal lattice, Bohr^-1)
    positions = take("<f8", 4 * atomnum, (atomnum, 4))[:, 1:4]

    nspin_blocks = spin_p_switch + 1

    def read_stack():
        """One full neighbour-block sweep, in file order."""
        blocks = []
        for ct in range(atomnum):
            for h in range(int(FNAN[ct]) + 1):
                gh = int(natn[ct][h]) - 1  # to 0-based
                n = int(Total_NumOrbs[ct])
                m = int(Total_NumOrbs[gh])
                blocks.append(take("<f8", n * m).reshape(n, m))
        return blocks

    H_blocks = [read_stack() for _ in range(nspin_blocks)]
    iH_blocks = [read_stack() for _ in range(3)] if spin_p_switch == 3 else None

    OLP_blocks = read_stack()
    for _ in range(3 * order_max):  # OLPpo (position operator)
        read_stack()
    for _ in range(3):  # OLPmo (momentum operator)
        read_stack()
    for _ in range(nspin_blocks):  # DM
        read_stack()
    for _ in range(2):  # iDM
        read_stack()

    solver = int(take("<i4", 1)[0])
    d10 = take("<f8", 10)
    chem_p, e_temp = float(d10[0]), float(d10[1])
    valence_electrons, total_spin_s = float(d10[8]), float(d10[9])

    num_lines = int(take("<i4", 1)[0])
    input_lines = []
    for _ in range(num_lines):
        line = take("V1", 256).tobytes()
        input_lines.append(line.split(b"\x00")[0].decode("latin-1"))

    # --- assemble dense real-space tables ---------------------------------
    norb = int(Total_NumOrbs.sum())
    S_dense = np.zeros((nR, norb, norb), dtype=float)
    H_dense = [np.zeros((nR, norb, norb), dtype=float) for _ in range(nspin_blocks)]
    iH_dense = (
        [np.zeros((nR, norb, norb), dtype=float) for _ in range(3)]
        if spin_p_switch == 3
        else None
    )

    Rindex = {}
    for iR in range(nR):
        Rindex[tuple(int(v) for v in atv_ijk[iR])] = iR

    orb_off = np.concatenate([[0], np.cumsum(Total_NumOrbs)])

    b = 0
    for ct in range(atomnum):
        for h in range(int(FNAN[ct]) + 1):
            gh = int(natn[ct][h]) - 1
            iR = Rindex[tuple(int(v) for v in atv_ijk[int(ncn[ct][h])])]
            sl = np.s_[orb_off[ct] : orb_off[ct + 1], orb_off[gh] : orb_off[gh + 1]]
            S_dense[iR][sl] = OLP_blocks[b]
            for spin in range(nspin_blocks):
                H_dense[spin][iR][sl] = H_blocks[spin][b]
            if iH_dense is not None:
                for spin in range(3):
                    iH_dense[spin][iR][sl] = iH_blocks[spin][b]
            b += 1

    if spin_p_switch == 3:
        H_out = np.stack(
            [
                assemble_ncol(
                    np.stack([H_dense[s][iR] for s in range(4)]),
                    np.stack([iH_dense[s][iR] for s in range(3)]),
                )
                for iR in range(nR)
            ]
        )
    elif spin_p_switch == 1:
        H_out = np.stack(H_dense)  # (2, nR, norb, norb)
    else:
        H_out = H_dense[0]  # (nR, norb, norb)

    return Scfout(
        natom=atomnum,
        spin_p_switch=spin_p_switch,
        version=version,
        tcpy_cell=tcpy_cell,
        solver=solver,
        total_numorbs=np.asarray(Total_NumOrbs, dtype=int),
        fnan=np.asarray(FNAN, dtype=int),
        natn=natn,
        ncn=ncn,
        atv_ijk=atv_ijk,
        tv=tv,
        positions=positions,
        H=H_out,
        S=S_dense,
        chem_p=chem_p,
        e_temp=e_temp,
        valence_electrons=valence_electrons,
        total_spin_s=total_spin_s,
        input_lines=input_lines,
    )
