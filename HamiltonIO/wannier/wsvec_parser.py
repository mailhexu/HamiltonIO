"""Parser for Wannier90 ``_wsvec.dat`` (per-orbital-pair Wigner-Seitz images).

Written by Wannier90 when ``write_wsvec=.true.`` (see ``ws_distance.F90``,
``ws_write_vec``). For each orbital pair ``(i, j)`` and lattice vector ``R`` it
lists the supercell translations ``T`` that bring Wannier function ``j`` at
``R`` into the Wigner-Seitz cell of ``i`` at the origin.

File format::

    ## written on {date} at {time} with use_ws_distance={.true.|.false.}
    {Rx} {Ry} {Rz} {iw} {jw}     # block header (1-based iw, jw)
    {N_T}                         # number of image translations
    {Tx} {Ty} {Tz}                # N_T lines of integer shifts
    ...                           # next block, loop order: irpt, iw, jw

The shift ``T`` is ``irdist_ws - irvec`` (the supercell translation only); the
full lattice vector for the Fourier phase is ``R + T``. When
``use_ws_distance=.false.`` every block is the trivial ``N_T=1, T=(0,0,0)``.

Reference: research/2026-07-26-wannier90-wsvec-weights.md, ADR-001.
"""

import warnings

import numpy as np


def parse_wsvec(filename):
    """Parse a Wannier90 ``_wsvec.dat`` file.

    Returns a dict with:
      - ``use_ws_distance``: bool from the file header.
      - ``shifts``: ``{R_tuple: {(i, j): ndarray(shape=(N_T, 3), dtype=int)}}``
        with 0-based orbital indices ``i, j``.

    Raises ValueError on a truncated/malformed file. Emits a ``UserWarning`` when
    the header says ``use_ws_distance=.false.`` (the file is the trivial
    identity and carries no per-pair correction).
    """
    with open(filename) as f:
        lines = f.readlines()

    if not lines:
        raise ValueError(f"{filename}: empty wsvec file")

    header = lines[0].strip()
    if "use_ws_distance" not in header:
        raise ValueError(
            f"{filename}: missing 'use_ws_distance' marker in header: {header!r}"
        )
    use_ws_distance = ".true." in header

    if not use_ws_distance:
        warnings.warn(
            f"{filename}: wsvec header says use_ws_distance=.false.; "
            "file is trivial (all N_T=1, T=(0,0,0)). Using global ndegen only.",
            UserWarning,
            stacklevel=2,
        )

    shifts = {}
    idx = 1
    n = len(lines)
    while idx < n:
        raw = lines[idx].strip()
        idx += 1
        if not raw:
            continue
        parts = raw.split()
        # Block header must be 5 integers: Rx Ry Rz iw jw
        if len(parts) < 5:
            raise ValueError(
                f"{filename}: line {idx} is not a valid block header: {raw!r}"
            )
        try:
            rx, ry, rz, iw, jw = (int(p) for p in parts[:5])
        except ValueError as exc:
            raise ValueError(
                f"{filename}: line {idx} block header not integer: {raw!r}"
            ) from exc
        R = (rx, ry, rz)
        i = iw - 1
        j = jw - 1

        # N_T count line
        nt = _read_next_int(lines, idx, filename)
        idx += 1

        if nt < 1:
            raise ValueError(
                f"{filename}: block (R={R}, i={i}, j={j}) has N_T={nt} < 1"
            )

        tvecs = np.zeros((nt, 3), dtype=int)
        for k in range(nt):
            tline, idx = _read_next_nonblank(lines, idx, filename)
            tparts = tline.split()
            if len(tparts) < 3:
                raise ValueError(
                    f"{filename}: shift line {idx} needs 3 ints: {tline!r}"
                )
            try:
                tvecs[k] = [int(tparts[0]), int(tparts[1]), int(tparts[2])]
            except ValueError as exc:
                raise ValueError(
                    f"{filename}: shift line {idx} not integer: {tline!r}"
                ) from exc

        shifts.setdefault(R, {})[(i, j)] = tvecs

    return {"use_ws_distance": use_ws_distance, "shifts": shifts}


def _read_next_int(lines, idx, filename):
    line, new_idx = _read_next_nonblank(lines, idx, filename)
    try:
        return int(line.split()[0])
    except (ValueError, IndexError) as exc:
        raise ValueError(
            f"{filename}: expected integer count at line {new_idx}: {line!r}"
        ) from exc


def _read_next_nonblank(lines, idx, filename):
    while idx < len(lines):
        line = lines[idx].strip()
        idx += 1
        if line:
            return line, idx
    raise EOFError(f"{filename}: unexpected EOF while reading block data")


def validate_ws_weights(hr_path, mp_grid, wsvec_path=None):
    """Validate Wannier90 Wigner-Seitz weights (story-030).

    Checks the WS sum rule on ``_hr.dat``:
        sum_R 1/ndegen(R) == Nk1 * Nk2 * Nk3
    (Wannier90 ``hamiltonian.F90:604-607``), and, when a ``_wsvec.dat`` is
    provided, reports per-pair image-count (``N_T``) statistics.

    Args:
        hr_path: path to ``{prefix}_hr.dat``.
        mp_grid: ``(Nk1, Nk2, Nk3)`` tuple (the Wannier90 ``mp_grid``).
        wsvec_path: optional path to ``{prefix}_wsvec.dat``.

    Returns:
        dict with ``sum_rule_ok`` (bool), ``sum_rule_value``,
        ``sum_rule_expected`` (mp_grid product), ``num_wann``, ``nrpts``, and
        ``ws_stats`` (``None`` if no wsvec; otherwise ``{use_ws_distance, n_blocks,
        min, max, mean}`` of the per-pair ``N_T`` counts).
    """
    # Local import to keep parse_wsvec importable without parse_ham at module load.
    from .w90_parser import parse_ham

    n_wann, _data, rdegens = parse_ham(hr_path)
    rdegens = np.asarray(rdegens, dtype=float)
    total = float(np.sum(1.0 / rdegens))
    expected = int(np.prod(mp_grid))
    result = {
        "sum_rule_value": total,
        "sum_rule_expected": expected,
        "sum_rule_ok": bool(np.isclose(total, expected, atol=1e-6)),
        "num_wann": n_wann,
        "nrpts": int(len(rdegens)),
        "ws_stats": None,
    }
    if wsvec_path is not None:
        ws = parse_wsvec(wsvec_path)
        if ws["use_ws_distance"]:
            nts = [
                len(arr) for rblock in ws["shifts"].values() for arr in rblock.values()
            ]
            result["ws_stats"] = {
                "use_ws_distance": True,
                "n_blocks": len(nts),
                "min": int(min(nts)) if nts else 0,
                "max": int(max(nts)) if nts else 0,
                "mean": float(np.mean(nts)) if nts else 0.0,
            }
        else:
            result["ws_stats"] = {"use_ws_distance": False}
    return result
