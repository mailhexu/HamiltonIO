"""Read ABINIT 9/10 ETSF netCDF WFK wavefunctions and eigenvalues."""

from __future__ import annotations

import os
import re
from dataclasses import dataclass

import numpy as np

HARTREE_TO_EV = 27.211386245988
_SUPPORTED_ABINIT_MAJORS = {9, 10}


def _decode_char(value) -> str:
    """Decode an ETSF fixed-width character variable without netCDF helpers."""
    array = np.asarray(value)
    if array.dtype.kind == "S":
        raw = b"".join(array.reshape(-1).tolist())
        return raw.decode("ascii", errors="replace").rstrip("\x00 ")
    if array.dtype.kind == "U":
        return "".join(array.reshape(-1).tolist()).rstrip("\x00 ")
    return str(array.item()).rstrip("\x00 ")


def _require_netcdf4():
    try:
        from netCDF4 import Dataset
    except ModuleNotFoundError as exc:
        if exc.name != "netCDF4":
            raise
        raise ImportError(
            "netCDF4 is required for ABINIT WFK unfolding. "
            "Install netCDF4 in the HamiltonIO environment"
        ) from exc
    return Dataset


@dataclass(frozen=True)
class WFKData:
    """Validated ABINIT WFK arrays in Hartree, independent of any unfolder."""

    kpoints: np.ndarray
    gvecs: tuple[np.ndarray, ...]
    coefficients: tuple[np.ndarray, ...]
    eigenvalues: tuple[np.ndarray, ...]
    rprimd: np.ndarray
    codvsn: str
    istwfk: np.ndarray
    fermi_energy: float | None = None
    usepaw: int = 0

    def __post_init__(self):
        def frozen(array, dtype):
            value = np.array(array, dtype=dtype, copy=True)
            value.setflags(write=False)
            return value

        kpoints = frozen(self.kpoints, float)
        if kpoints.ndim != 2 or kpoints.shape[1] != 3 or len(kpoints) == 0:
            raise ValueError("kpoints must have nonzero shape (nk, 3)")
        if not (
            len(self.gvecs)
            == len(self.coefficients)
            == len(self.eigenvalues)
            == len(kpoints)
        ):
            raise ValueError("WFK arrays need one entry per k-point")
        gvecs, coefficients, eigenvalues = [], [], []
        shape = None
        for ik, (g, c, e) in enumerate(
            zip(self.gvecs, self.coefficients, self.eigenvalues)
        ):
            g, c, e = frozen(g, int), frozen(c, complex), frozen(e, float)
            if g.ndim != 2 or g.shape[1] != 3 or not len(g):
                raise ValueError(f"invalid G vectors at k-point {ik}")
            if (
                c.ndim != 4
                or e.ndim != 2
                or c.shape[:2] != e.shape
                or c.shape[-1] != len(g)
            ):
                raise ValueError(
                    f"invalid coefficient/eigenvalue shape at k-point {ik}"
                )
            if shape is not None and c.shape[:3] != shape:
                raise ValueError(
                    f"inconsistent spin/band/spinor dimensions at k-point {ik}"
                )
            shape = c.shape[:3]
            gvecs.append(g)
            coefficients.append(c)
            eigenvalues.append(e)
        rprimd, istwfk = frozen(self.rprimd, float), frozen(self.istwfk, int)
        if rprimd.shape != (3, 3) or istwfk.shape != (len(kpoints),):
            raise ValueError("invalid primitive_vectors or istwfk shape")
        if not isinstance(self.codvsn, str) or not self.codvsn.strip():
            raise ValueError("codvsn must be non-empty")
        for key, val in (
            ("kpoints", kpoints),
            ("gvecs", tuple(gvecs)),
            ("coefficients", tuple(coefficients)),
            ("eigenvalues", tuple(eigenvalues)),
            ("rprimd", rprimd),
            ("istwfk", istwfk),
        ):
            object.__setattr__(self, key, val)
        if self.fermi_energy is not None:
            object.__setattr__(self, "fermi_energy", float(self.fermi_energy))

    @property
    def nspinor(self):
        return self.coefficients[0].shape[2]

    @property
    def nspin(self):
        return self.coefficients[0].shape[0]


def _wfk_variable(nc, name):
    try:
        return nc.variables[name]
    except KeyError as exc:
        raise ValueError(f"not an ETSF ABINIT WFK: missing variable {name!r}") from exc


def read_wfk(path) -> WFKData:
    """Read an ABINIT 9/10 ETSF netCDF WFK into backend-free eigendata.

    Header and metadata arrays are small. Basis and band arrays are read one
    k-point at a time, with padded WFK axes sliced to
    `number_of_coefficients` and `number_of_states` before materializing
    them. Fortran WFKs, half-sphere storage, and real-only coefficient storage
    are rejected deliberately because their reconstruction semantics are
    outside the first adapter contract.
    """
    if not os.path.isfile(path):
        raise FileNotFoundError(f"ABINIT WFK not found: {path!s}")
    Dataset = _require_netcdf4()
    try:
        nc = Dataset(path, "r")
    except (OSError, RuntimeError) as exc:
        raise ValueError(
            f"cannot read {path!s} as a netCDF ABINIT WFK; "
            "rerun ABINIT with iomode 3 to produce an ETSF netCDF WFK"
        ) from exc

    try:
        kpoints = np.asarray(
            _wfk_variable(nc, "reduced_coordinates_of_kpoints")[:], dtype=float
        )
        npw = np.asarray(_wfk_variable(nc, "number_of_coefficients")[:], dtype=int)
        nstates = np.asarray(_wfk_variable(nc, "number_of_states")[:], dtype=int)
        istwfk = np.asarray(_wfk_variable(nc, "istwfk")[:], dtype=int)
        usepaw = int(nc.variables["usepaw"][:]) if "usepaw" in nc.variables else 0
        gvec_var = _wfk_variable(nc, "reduced_coordinates_of_plane_waves")
        coeff_var = _wfk_variable(nc, "coefficients_of_wavefunctions")
        eig_var = _wfk_variable(nc, "eigenvalues")
        rprimd = np.asarray(_wfk_variable(nc, "primitive_vectors")[:], dtype=float)
        codvsn = _decode_char(_wfk_variable(nc, "codvsn")[:])
        fermi = (
            float(nc.variables["fermi_energy"][:])
            if "fermi_energy" in nc.variables
            else None
        )

        match = re.match(r"\s*(\d+)", codvsn)
        if match is None or int(match.group(1)) not in _SUPPORTED_ABINIT_MAJORS:
            supported = ", ".join(str(v) for v in sorted(_SUPPORTED_ABINIT_MAJORS))
            raise ValueError(
                f"unsupported ABINIT WFK major version {codvsn!r}; supported ABINIT WFK major versions: {supported}"
            )
        if kpoints.ndim != 2 or kpoints.shape[1] != 3:
            raise ValueError(
                "not an ETSF ABINIT WFK: invalid reduced_coordinates_of_kpoints shape"
            )
        if npw.shape != (len(kpoints),) or istwfk.shape != (len(kpoints),):
            raise ValueError(
                "not an ETSF ABINIT WFK: inconsistent per-k npw or istwfk dimensions"
            )
        if np.any(istwfk != 1):
            bad = np.where(istwfk != 1)[0].tolist()
            raise ValueError(
                f"WFK has istwfk={istwfk.tolist()} (bad k-point indices {bad}); "
                "rerun ABINIT with istwfk 1 for full-G planewave unfolding"
            )

        coeff_shape = tuple(coeff_var.shape)
        if len(coeff_shape) != 6:
            raise ValueError(
                "not an ETSF ABINIT WFK: invalid coefficient array dimensions"
            )
        nspin, nk, mband, nspinor, max_npw, ncomplex = coeff_shape
        if ncomplex != 2:
            raise ValueError(
                "WFK uses real coefficient storage; this adapter requires "
                "complex real/imag coefficient pairs"
            )
        eig_shape = tuple(eig_var.shape)
        if nk != len(kpoints) or len(eig_shape) != 3 or eig_shape[:2] != (nspin, nk):
            raise ValueError(
                "not an ETSF ABINIT WFK: inconsistent coefficient/eigenvalue dimensions"
            )
        if nstates.shape != (nspin, nk):
            raise ValueError("not an ETSF ABINIT WFK: invalid number_of_states shape")
        gvec_shape = tuple(gvec_var.shape)
        if gvec_shape != (nk, max_npw, 3):
            raise ValueError(
                "not an ETSF ABINIT WFK: invalid reduced plane-wave dimensions"
            )

        gvecs: list[np.ndarray] = []
        coefficients: list[np.ndarray] = []
        eigenvalues: list[np.ndarray] = []
        for ik in range(nk):
            if not np.all(nstates[:, ik] == nstates[0, ik]):
                raise ValueError(
                    "WFK has spin-dependent nband; this adapter requires common nband"
                )
            nband = int(nstates[0, ik])
            nplane = int(npw[ik])
            if not (
                0 < nband <= mband and nband <= eig_shape[2] and 0 < nplane <= max_npw
            ):
                raise ValueError(f"WFK has invalid nband/npw at k-point {ik}")
            raw = np.asarray(coeff_var[:, ik, :nband, :, :nplane, :], dtype=float)
            coefficients.append(raw[..., 0] + 1j * raw[..., 1])
            eigenvalues.append(np.asarray(eig_var[:, ik, :nband], dtype=float))
            gvecs.append(np.asarray(gvec_var[ik, :nplane, :], dtype=int))

        return WFKData(
            kpoints=kpoints,
            gvecs=tuple(gvecs),
            coefficients=tuple(coefficients),
            eigenvalues=tuple(eigenvalues),
            rprimd=rprimd,
            codvsn=codvsn,
            istwfk=istwfk,
            fermi_energy=fermi,
            usepaw=usepaw,
        )
    finally:
        nc.close()
