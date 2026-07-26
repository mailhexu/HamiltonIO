# -*- coding: utf-8 -*-
import os
import re
import struct
from typing import NoReturn

import numpy as np
from scipy.sparse import csr_matrix

MAX_DENSE_CSR_ELEMENTS = 20_000_000
MAX_NEW_CSR_FILE_BYTES = 512 * 1024 * 1024


def _validate_dense_dimensions(filename, basis_num, r_num):
    if basis_num <= 0 or r_num < 0:
        raise ValueError(
            f"{filename}: invalid CSR dimensions basis={basis_num}, R={r_num}"
        )
    elements = basis_num * basis_num * r_num
    if elements > MAX_DENSE_CSR_ELEMENTS:
        raise ValueError(
            f"{filename}: dense CSR element count {elements} exceeds MAX_DENSE_CSR_ELEMENTS={MAX_DENSE_CSR_ELEMENTS}"
        )


def _validate_block_data_size(filename, basis_num, data_size):
    if data_size < 0 or data_size > basis_num * basis_num:
        raise ValueError(
            f"{filename}: invalid CSR data_size={data_size} for basis={basis_num}"
        )


class XR_matrix:
    def __init__(self, nspin, XR_fileName):
        self.nspin = nspin
        self.XR_fileName = XR_fileName

    def read_file(self):
        with open(self.XR_fileName, "r") as fread:
            # Handle both old (v3.5.x–v3.10-LTS) and new (v3.11.0+) ABACUS CSR formats.
            # Old: "Matrix Dimension of H(R): N" / "Matrix number of H(R): N"
            # New: "N # number of localized basis" / "N # number of Bravais lattice vector R"
            # plus metadata header + comment lines (#) interspersed with data.
            basis_num = None
            R_num = None
            for line in fread:
                parts = line.split()
                if not parts:
                    continue
                # Old format: "Matrix Dimension of H(R): 54"
                if parts[0] == "Matrix" and basis_num is None:
                    basis_num = int(parts[-1])
                    continue
                # Old format: "Matrix number of H(R): 221"
                if parts[0] == "Matrix" and basis_num is not None:
                    R_num = int(parts[-1])
                    break
                # New format: "54 # number of localized basis"
                if (
                    len(parts) >= 5
                    and parts[1] == "#"
                    and "localized" in line
                    and "basis" in line
                ):
                    basis_num = int(parts[0])
                    continue
                # New format: "221 # number of Bravais lattice vector R"
                if (
                    len(parts) >= 2
                    and parts[1] == "#"
                    and "Bravais" in line
                    and "vector" in line
                ):
                    R_num = int(parts[0])
                    break

            if basis_num is None or R_num is None:
                raise ValueError(
                    f"Could not parse CSR header from {self.XR_fileName}. "
                    f"Expected old format ('Matrix Dimension of H(R): N') or "
                    f"new format ('N # number of localized basis')."
                )

            _validate_dense_dimensions(self.XR_fileName, basis_num, R_num)

            self.basis_num = basis_num
            self.R_num = R_num
            self.R_direct_coor = np.zeros([self.R_num, 3], dtype=int)
            if self.nspin != 4:
                self.XR = np.zeros(
                    [self.R_num, self.basis_num, self.basis_num], dtype=float
                )
            else:
                self.XR = np.zeros(
                    [self.R_num, self.basis_num, self.basis_num], dtype=complex
                )

            for iR in range(self.R_num):
                while True:
                    line = fread.readline()
                    if not line:
                        raise EOFError(
                            f"Unexpected EOF reading R-vector {iR} from {self.XR_fileName}"
                        )
                    parts = line.split()
                    if len(parts) < 4 or parts[0].startswith("#"):
                        continue
                    try:
                        int(parts[3])
                        break
                    except ValueError:
                        continue

                parts = line.split()
                self.R_direct_coor[iR, 0] = int(parts[0])
                self.R_direct_coor[iR, 1] = int(parts[1])
                self.R_direct_coor[iR, 2] = int(parts[2])
                data_size = int(parts[3])
                _validate_block_data_size(self.XR_fileName, self.basis_num, data_size)

                if self.nspin != 4:
                    data = np.zeros((data_size,), dtype=float)
                else:
                    data = np.zeros((data_size,), dtype=complex)

                indices = np.zeros((data_size,), dtype=int)
                indptr = np.zeros((self.basis_num + 1,), dtype=int)

                if data_size != 0:
                    if self.nspin != 4:
                        data_values = self._read_data_tokens(fread, data_size)
                        for index in range(data_size):
                            data[index] = float(data_values[index])
                    else:
                        complex_str = self._read_data_text(fread, data_size)
                        matches = re.findall(r"\(([^)]+)\)", complex_str)
                        if len(matches) < data_size:
                            raise ValueError(
                                f"Expected {data_size} complex values, got {len(matches)}"
                            )
                        for index in range(data_size):
                            parts_c = matches[index].split(",")
                            data[index] = complex(float(parts_c[0]), float(parts_c[1]))

                    indices_tokens = self._read_data_tokens(fread, data_size)
                    for index in range(data_size):
                        indices[index] = int(indices_tokens[index])

                    indptr_tokens = self._read_data_tokens(fread, self.basis_num + 1)
                    for index in range(self.basis_num + 1):
                        indptr[index] = int(indptr_tokens[index])

                    if not np.isfinite(data).all():
                        raise ValueError(
                            f"{self.XR_fileName}: CSR values must be finite"
                        )

                self.XR[iR] = csr_matrix(
                    (data, indices, indptr), shape=(self.basis_num, self.basis_num)
                ).toarray()

    @staticmethod
    def _read_data_tokens(fread, count):
        tokens = []
        while len(tokens) < count:
            line = fread.readline()
            if not line:
                raise EOFError("Unexpected EOF reading CSR data")
            parts = line.split()
            if not parts or parts[0].startswith("#"):
                continue
            tokens.extend(parts)
        return tokens[:count]

    @staticmethod
    def _read_data_text(fread, min_matches):
        text = ""
        while text.count("(") < min_matches:
            line = fread.readline()
            if not line:
                raise EOFError("Unexpected EOF reading CSR complex data")
            if line.strip().startswith("#"):
                continue
            text += " " + line
        return text

    def read_file_binary(self):
        with open(self.XR_fileName, "rb") as fread:
            self.basis_num = struct.unpack("i", fread.read(4))[0]
            self.R_num = struct.unpack("i", fread.read(4))[0]
            _validate_dense_dimensions(self.XR_fileName, self.basis_num, self.R_num)
            self.R_direct_coor = np.zeros([self.R_num, 3], dtype=int)
            if self.nspin != 4:
                self.XR = np.zeros(
                    [self.R_num, self.basis_num, self.basis_num], dtype=float
                )
            else:
                self.XR = np.zeros(
                    [self.R_num, self.basis_num, self.basis_num], dtype=complex
                )

            for iR in range(self.R_num):
                self.R_direct_coor[iR, 0] = struct.unpack("i", fread.read(4))[0]
                self.R_direct_coor[iR, 1] = struct.unpack("i", fread.read(4))[0]
                self.R_direct_coor[iR, 2] = struct.unpack("i", fread.read(4))[0]
                data_size = struct.unpack("i", fread.read(4))[0]
                _validate_block_data_size(self.XR_fileName, self.basis_num, data_size)

                if self.nspin != 4:
                    data = np.zeros((data_size,), dtype=float)
                else:
                    data = np.zeros((data_size,), dtype=complex)

                indices = np.zeros((data_size,), dtype=int)
                indptr = np.zeros((self.basis_num + 1,), dtype=int)

                if data_size != 0:
                    if self.nspin != 4:
                        for index in range(data_size):
                            data[index] = struct.unpack("d", fread.read(8))[0]
                    else:
                        for index in range(data_size):
                            real = struct.unpack("d", fread.read(8))[0]
                            imag = struct.unpack("d", fread.read(8))[0]
                            data[index] = complex(real, imag)

                    for index in range(data_size):
                        indices[index] = struct.unpack("i", fread.read(4))[0]

                    for index in range(self.basis_num + 1):
                        indptr[index] = struct.unpack("i", fread.read(4))[0]

                    if not np.isfinite(data).all():
                        raise ValueError(
                            f"{self.XR_fileName}: CSR values must be finite"
                        )

                self.XR[iR] = csr_matrix(
                    (data, indices, indptr), shape=(self.basis_num, self.basis_num)
                ).toarray()


def read_HR_SR(
    nspin=4,
    binary=False,
    HR_fileName="data-HR-sparse_SPIN0.csr",
    SR_fileName="data-SR-sparse_SPIN0.csr",
):
    """
    IN:
        nspin: int, different spins.
        binary: bool, whether the HR and SR matrices are binary files.
        HR_fileName: if nspin=1 or 4, str, HR file name;
                     if nspin=2, list or tuple, size=2, [HR_up_fileName, HR_dn_fileName].
        SR_fileName: str, SR file name.

    OUT:
        if nspin = 1 or 4:
            return basis_num, R_direct_coor, HR, SR
        elif nspin = 2:
            return basis_num, R_direct_coor, HR_up, HR_dn, SR

        basis_num: int, number of atomic orbital basis.
        R_direct_coor: numpy ndarray, shape=[R_num, 3], fractional coordinates (x, y, z) of the R-th primitive cell.
        HR or HR_up or HR_dn: numpy ndarray, shape=[R_num, basis_num, basis_num], HR matrix and unit is eV.
        SR: numpy ndarray, shape=[R_num, basis_num, basis_num], SR matrix.
    """
    Ry_to_eV = 13.605698066

    if nspin == 1 or nspin == 4:
        if not isinstance(HR_fileName, str):
            raise ValueError("The HR_fileName must be a str for nspin=1 or 4.")

        HR = XR_matrix(nspin, HR_fileName)
        SR = XR_matrix(nspin, SR_fileName)

        if binary:
            HR.read_file_binary()
            SR.read_file_binary()
        else:
            HR.read_file()
            SR.read_file()

        return HR.basis_num, HR.R_direct_coor, HR.XR * Ry_to_eV, SR.XR
    else:
        if not isinstance(HR_fileName, (list, tuple)):
            raise ValueError("The HR_fileName must be a list or a tuple for nspin=2.")

        if len(HR_fileName) != 2:
            raise ValueError("The size of the HR_fileName must be 2 for nspin=2.")

        HR_up = XR_matrix(nspin, HR_fileName[0])
        HR_dn = XR_matrix(nspin, HR_fileName[1])
        SR = XR_matrix(nspin, SR_fileName)

        if binary:
            HR_up.read_file_binary()
            HR_dn.read_file_binary()
            SR.read_file_binary()
        else:
            HR_up.read_file()
            HR_dn.read_file()
            SR.read_file()

        return (
            HR_up.basis_num,
            HR_up.R_direct_coor,
            HR_up.XR * Ry_to_eV,
            HR_dn.XR * Ry_to_eV,
            SR.XR,
        )


def read_csr_new_format(filename, is_complex=True):
    """Read the ABACUS 3.11+ text CSR format without accepting partial data."""
    filename = str(filename)
    if os.stat(filename).st_size > MAX_NEW_CSR_FILE_BYTES:
        raise ValueError(
            f"{filename}: CSR file size exceeds MAX_NEW_CSR_FILE_BYTES={MAX_NEW_CSR_FILE_BYTES}"
        )
    with open(filename, "r") as f:
        lines = f.readlines()

    def fail(message, line=None) -> NoReturn:
        location = f" at line {line + 1}" if line is not None else ""
        raise ValueError(f"{filename}{location}: {message}")

    basis_matches = []
    r_matches = []
    marker = None
    for line_number, line in enumerate(lines):
        text = line.strip()
        if text == "# CSR Format" or ("CSR Format" in text and text.startswith("#")):
            if marker is not None:
                fail("duplicate CSR Format marker", line_number)
            marker = line_number
        match = re.fullmatch(r"(\d+)\s+# number of localized basis", text)
        if match:
            basis_matches.append((line_number, int(match.group(1))))
        match = re.fullmatch(r"(\d+)\s+# number of Bravais lattice vector R", text)
        if match:
            r_matches.append((line_number, int(match.group(1))))
    if len(basis_matches) != 1:
        fail("expected exactly one localized basis header")
    if len(r_matches) != 1:
        fail("expected exactly one Bravais R-count header")
    if marker is None:
        fail("missing CSR Format marker")
    if basis_matches[0][0] > marker or r_matches[0][0] > marker:
        fail("CSR headers must precede the CSR Format marker")
    basis_num = basis_matches[0][1]
    R_num = r_matches[0][1]
    try:
        _validate_dense_dimensions(filename, basis_num, R_num)
    except ValueError as error:
        fail(str(error))

    def next_content(index):
        while index < len(lines):
            text = lines[index].strip()
            if text and not (text.startswith("#") and "CSR" not in text):
                break
            index += 1
        return index

    markers = {
        "values": ("# CSR values",),
        "column_indices": ("# CSR column_indices", "# CSR column indices"),
        "row_indptr": ("# CSR row_indptr", "# CSR row pointers"),
    }

    def section(index, name, next_name, block, allow_empty=False, limit=None):
        index = next_content(index)
        expected = markers[name]
        if index == len(lines) or lines[index].strip() not in expected:
            fail(
                f"R block {block}: expected {expected[0]}",
                index if index < len(lines) else None,
            )
        payload = []
        token_count = 0
        index += 1
        expected_next = markers[next_name]
        while True:
            index = next_content(index)
            if index == len(lines):
                fail(f"R block {block}: expected {expected_next[0]}")
            text = lines[index].strip()
            if text in expected_next:
                if not payload and not allow_empty:
                    fail(f"R block {block}: missing {name} data", index)
                return " ".join(payload), index
            if text.startswith("#"):
                fail(f"R block {block}: expected {expected_next[0]}", index)
            payload.append(text)
            token_count += len(text.split())
            if limit is not None and token_count > limit:
                fail(f"R block {block}: excess {name} tokens", index)
            index += 1

    def final_section(index, block):
        index = next_content(index)
        expected = markers["row_indptr"]
        if index == len(lines) or lines[index].strip() not in expected:
            fail(
                f"R block {block}: expected {expected}",
                index if index < len(lines) else None,
            )
        payload = []
        index += 1
        while len(payload) < basis_num + 1:
            index = next_content(index)
            if index == len(lines):
                fail(f"R block {block}: truncated row_indptr")
            text = lines[index].strip()
            if text.startswith("#"):
                fail(f"R block {block}: invalid row_indptr", index)
            payload.extend(text.split())
            if len(payload) > basis_num + 1:
                fail(f"R block {block}: invalid indptr", index)
            index += 1
        return payload, index

    dtype = complex if is_complex else float
    R_direct_coor = np.empty((R_num, 3), dtype=int)
    XR = np.zeros((R_num, basis_num, basis_num), dtype=dtype)
    index = marker + 1
    seen_r = set()
    float_pattern = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
    complex_pattern = re.compile(
        rf"\[?\s*\(\s*({float_pattern})\s*,\s*({float_pattern})\s*\)\s*\]?"
    )

    for block in range(R_num):
        index = next_content(index)
        if index == len(lines):
            fail(f"expected R block {block + 1} of {R_num}")
        parts = lines[index].split()
        if len(parts) != 4:
            fail(f"R block {block}: expected 'Rx Ry Rz nnz'", index)
        try:
            r = tuple(int(part) for part in parts[:3])
            nnz = int(parts[3])
        except ValueError:
            fail(f"R block {block}: R coordinates and nnz must be integers", index)
        if nnz < 0:
            fail(f"R block {block}: nnz must be non-negative", index)
        if r in seen_r:
            fail(f"R block {block}: duplicate R vector {r}", index)
        seen_r.add(r)
        R_direct_coor[block] = r
        index += 1
        if nnz == 0:
            section_start = next_content(index)
            if (
                section_start < len(lines)
                and lines[section_start].strip() in markers["values"]
            ):
                values_text, index = section(
                    index, "values", "column_indices", block, allow_empty=True
                )
                indices_text, index = section(
                    index, "column_indices", "row_indptr", block, allow_empty=True
                )
                indptr_tokens, index = final_section(index, block)
                if values_text or indices_text:
                    fail(f"R block {block}: zero nnz block has payload")
                try:
                    indptr = np.array(
                        [int(value) for value in indptr_tokens], dtype=int
                    )
                except ValueError:
                    fail(f"R block {block}: invalid row_indptr")
                if np.any(indptr):
                    fail(f"R block {block}: zero nnz block has nonzero indptr")
            continue

        values_text, index = section(
            index, "values", "column_indices", block, limit=nnz
        )
        if is_complex:
            matches = complex_pattern.findall(values_text)
            if complex_pattern.sub("", values_text).strip() or len(matches) != nnz:
                fail(f"R block {block}: expected exactly {nnz} complex values")
            data = np.array(
                [complex(float(real), float(imag)) for real, imag in matches]
            )
        else:
            try:
                data = np.array([float(value) for value in values_text.split()])
            except ValueError:
                fail(f"R block {block}: invalid values")
            if len(data) != nnz:
                fail(f"R block {block}: expected exactly {nnz} values")
        if not np.isfinite(data).all():
            fail(f"R block {block}: values must be finite")

        indices_text, index = section(
            index, "column_indices", "row_indptr", block, limit=nnz
        )
        indptr_tokens, index = final_section(index, block)
        try:
            indices = np.array(
                [int(value) for value in indices_text.split()], dtype=int
            )
        except ValueError:
            fail(f"R block {block}: invalid column indices")
        try:
            indptr = np.array([int(value) for value in indptr_tokens], dtype=int)
        except ValueError:
            fail(f"R block {block}: invalid row_indptr")
        if len(indices) != nnz or np.any(indices < 0) or np.any(indices >= basis_num):
            fail(f"R block {block}: invalid indices")
        if (
            len(indptr) != basis_num + 1
            or indptr[0] != 0
            or indptr[-1] != nnz
            or np.any(np.diff(indptr) < 0)
        ):
            fail(f"R block {block}: invalid indptr")
        XR[block] = csr_matrix(
            (data, indices, indptr), shape=(basis_num, basis_num)
        ).toarray()

    for line_number in range(index, len(lines)):
        text = lines[line_number].strip()
        if text and not text.startswith("#"):
            fail("extra CSR data after declared R blocks", line_number)
        if "CSR" in text:
            fail("extra CSR section after declared R blocks", line_number)
    return basis_num, R_direct_coor, XR
