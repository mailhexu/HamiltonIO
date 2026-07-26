import re
import struct
from pathlib import Path

import numpy as np
import pytest

from HamiltonIO.abacus import (
    AbacusSingleStepSOCParser,
    AbacusSplitSOCParser,
    abacus_api,
    abacus_wrapper,
)
from HamiltonIO.abacus.abacus_api import (
    MAX_DENSE_CSR_ELEMENTS,
    MAX_NEW_CSR_FILE_BYTES,
    XR_matrix,
    read_csr_new_format,
    read_HR_SR,
)
from HamiltonIO.abacus.abacus_wrapper import AbacusParser, _align_csr_blocks
from HamiltonIO.mathutils.pauli import chargepart

DATA_DIR = Path(__file__).parent / "data" / "abacus"
RY_TO_EV = 13.605698066


def _csr_text(values, r=(0, 0, 0)):
    value_text = " ".join(f"[({real},{imag})]" for real, imag in values)
    return f""" --- Ionic Step 1 ---
 # print H matrix in real space H(R)
 4 # number of spin directions
 1 # spin index
 4 # number of localized basis
 1 # number of Bravais lattice vector R

 user_defined_lattice
 1
 1 0 0
 0 1 0
 0 0 1
 Fe
 1
 Direct
 0 0 0

 #----------------------------------------------------------------------#
 #                               CSR Format                             #
 #----------------------------------------------------------------------#

 {r[0]} {r[1]} {r[2]} 4
 # CSR values
 {value_text}
 # CSR column_indices
 0 1 2 3
 # CSR row_indptr
 0 1 2 3 4
"""


def _strict_csr_text(basis=2, blocks=(((0, 0, 0), 2),)):
    body = [
        "# metadata",
        f"{basis} # number of localized basis",
        f"{len(blocks)} # number of Bravais lattice vector R",
        "",
        "# CSR Format",
    ]
    for r, nnz in blocks:
        body.append(f"{r[0]} {r[1]} {r[2]} {nnz}")
        if nnz:
            body.extend(
                [
                    "# CSR values",
                    " ".join("(1.0,0.0)" for _ in range(nnz)),
                    "# CSR column_indices",
                    " ".join(str(i % basis) for i in range(nnz)),
                    "# CSR row_indptr",
                    "0 " + " ".join("0" for _ in range(basis - 1)) + f" {nnz}",
                ]
            )
    return "\n".join(body) + "\n"


def _new_csr_text(basis, blocks, imaginary=0.0):
    body = [
        f"{basis} # number of localized basis",
        f"{len(blocks)} # number of Bravais lattice vector R",
        "# CSR Format",
    ]
    for r, value in blocks:
        body.extend(
            [
                f"{r[0]} {r[1]} {r[2]} {basis}",
                "# CSR values",
                " ".join(f"({value},{imaginary})" for _ in range(basis)),
                "# CSR column_indices",
                " ".join(str(i) for i in range(basis)),
                "# CSR row_indptr",
                " ".join(str(i) for i in range(basis + 1)),
            ]
        )
    return "\n".join(body) + "\n"


def _write_single_step_soc_case(tmp_path, soc_r=(0, 0, 0)):
    out = tmp_path / "OUT.Fe"
    out.mkdir(parents=True)
    (tmp_path / "STRU").write_text(
        """ATOMIC_SPECIES
Fe 55.845 Fe.upf

NUMERICAL_ORBITAL
Fe.orb

LATTICE_CONSTANT
1.0

LATTICE_VECTORS
1 0 0
0 1 0
0 0 1

ATOMIC_POSITIONS
Direct

Fe
0.0
1
0 0 0 1 1 1 mag 1
"""
    )
    (out / "running_scf.log").write_text(
        """number of electrons = 6
nspin = 4
EFERMI = 1.25 eV
"""
    )
    (out / "Orbital").write_text(
        """  #io    spec    l    m    z  sym
    0      Fe    0    0    1              s
    0      Fe    1    0    1             pz

"""
    )
    (out / "hrs1_nao.csr").write_text(
        _csr_text([(1.0, 0.0), (2.0, 0.0), (3.0, 0.0), (4.0, 0.0)])
    )
    (out / "hrs1-soc_nao.csr").write_text(
        _csr_text([(0.1, 0.0), (0.2, 0.0), (0.3, 0.0), (0.4, 0.0)], soc_r)
    )
    (out / "srs1_nao.csr").write_text(
        _csr_text([(1.0, 0.0), (1.0, 0.0), (1.0, 0.0), (1.0, 0.0)])
    )
    return out


CALCULATIONS = {
    "old_nspin2": {
        "out_dir": DATA_DIR / "old_nspin2" / "OUT.Fe",
        "nspin": 2,
        "basis": 2,
        "hr_files": ["data-HR-sparse_SPIN0.csr", "data-HR-sparse_SPIN1.csr"],
        "sr_file": "data-SR-sparse_SPIN0.csr",
    },
    "old_nspin4": {
        "out_dir": DATA_DIR / "old_nspin4" / "OUT.Fe",
        "nspin": 4,
        "basis": 4,
        "hr_files": ["data-HR-sparse_SPIN0.csr"],
        "sr_file": "data-SR-sparse_SPIN0.csr",
    },
    "new_nspin2": {
        "out_dir": DATA_DIR / "new_nspin2" / "OUT.Fe",
        "nspin": 2,
        "basis": 2,
        "hr_files": ["hrs1_nao.csr", "hrs2_nao.csr"],
        "sr_file": "srs1_nao.csr",
    },
    "new_nspin4": {
        "out_dir": DATA_DIR / "new_nspin4" / "OUT.Fe",
        "nspin": 4,
        "basis": 4,
        "hr_files": ["hrs1_nao.csr"],
        "sr_file": "srs1_nao.csr",
    },
}


class TestCSRParsing:
    """Test CSR file parsing for all ABACUS version × nspin combinations."""

    @pytest.mark.parametrize("calc_name", list(CALCULATIONS.keys()))
    def test_all_csr_files_parse(self, calc_name):
        """Every CSR file in each calculation directory must parse with correct shape and dtype."""
        info = CALCULATIONS[calc_name]
        expected_dtype = complex if info["nspin"] == 4 else float
        for csr_file in info["hr_files"] + [info["sr_file"]]:
            xr = XR_matrix(
                nspin=info["nspin"],
                XR_fileName=str(info["out_dir"] / csr_file),
            )
            xr.read_file()
            assert xr.basis_num == info["basis"]
            assert xr.R_num >= 1
            assert xr.XR.shape[1] == xr.XR.shape[2] == xr.basis_num
            assert xr.XR.dtype == expected_dtype
            assert np.count_nonzero(xr.XR) > 0

    @pytest.mark.parametrize("calc_name", ["old_nspin2", "new_nspin2"])
    def test_nspin2_hr_spin_channels_differ(self, calc_name):
        """Spin-up and spin-down HR must have different values."""
        info = CALCULATIONS[calc_name]
        assert len(info["hr_files"]) == 2
        xr_up = XR_matrix(
            nspin=2, XR_fileName=str(info["out_dir"] / info["hr_files"][0])
        )
        xr_up.read_file()
        xr_dn = XR_matrix(
            nspin=2, XR_fileName=str(info["out_dir"] / info["hr_files"][1])
        )
        xr_dn.read_file()
        assert not np.allclose(xr_up.XR[0], xr_dn.XR[0])

    @pytest.mark.parametrize("calc_name", ["old_nspin2", "new_nspin2"])
    def test_nspin2_hr_values(self, calc_name):
        """Verify specific HR values match the fixture data."""
        info = CALCULATIONS[calc_name]
        xr = XR_matrix(nspin=2, XR_fileName=str(info["out_dir"] / info["hr_files"][0]))
        xr.read_file()
        expected = np.array([[1.0, 0.3], [0.0, 2.0]])
        assert np.allclose(xr.XR[0], expected)

    @pytest.mark.parametrize("calc_name", ["old_nspin4", "new_nspin4"])
    def test_nspin4_complex_values(self, calc_name):
        """nspin=4 CSR must produce complex matrices."""
        info = CALCULATIONS[calc_name]
        xr = XR_matrix(nspin=4, XR_fileName=str(info["out_dir"] / info["hr_files"][0]))
        xr.read_file()
        assert xr.XR.dtype == complex
        nonzero = xr.XR[xr.XR != 0]
        assert len(nonzero) >= 1
        assert np.isclose(nonzero[0], complex(1.0, 0.0))

    @pytest.mark.parametrize("calc_name", list(CALCULATIONS.keys()))
    def test_r_direct_coor(self, calc_name):
        """R-vectors should be integer arrays with correct shape."""
        info = CALCULATIONS[calc_name]
        xr = XR_matrix(
            nspin=info["nspin"], XR_fileName=str(info["out_dir"] / info["hr_files"][0])
        )
        xr.read_file()
        assert xr.R_direct_coor.dtype == int
        assert xr.R_direct_coor.shape == (xr.R_num, 3)

    @pytest.mark.parametrize("calc_name", list(CALCULATIONS.keys()))
    def test_input_files_present(self, calc_name):
        """Each calculation directory must contain INPUT, STRU, KPT, and pseudopotential."""
        calc_dir = DATA_DIR / calc_name
        for f in [
            "INPUT",
            "STRU",
            "KPT",
            "Fe_ONCV_PBE_FR-1.0.upf",
            "Fe_gga_8au_100Ry_4s2p2d1f.orb",
        ]:
            assert (calc_dir / f).exists(), f"{f} missing from {calc_name}"


class TestLogRegex:
    """Test running_scf.log token parsing for old and new ABACUS versions."""

    @pytest.mark.parametrize("calc_name", list(CALCULATIONS.keys()))
    def test_nbands_regex(self, calc_name):
        log = (CALCULATIONS[calc_name]["out_dir"] / "running_scf.log").read_text()
        m = re.search(r"NBANDS\)?\s*=\s*(\d+)", log)
        assert m is not None
        assert int(m.group(1)) > 0

    @pytest.mark.parametrize("calc_name", list(CALCULATIONS.keys()))
    def test_electrons_regex(self, calc_name):
        log = (CALCULATIONS[calc_name]["out_dir"] / "running_scf.log").read_text()
        m = re.search(
            r"(?i)autoset(?:\s+the)?\s+number of electrons:?\s*=\s*(\d+)", log
        )
        assert m is not None
        assert int(m.group(1)) > 0

    @pytest.mark.parametrize("calc_name", list(CALCULATIONS.keys()))
    def test_efermi_regex(self, calc_name):
        log = (CALCULATIONS[calc_name]["out_dir"] / "running_scf.log").read_text()
        m = re.search(r"EFERMI\s*=\s*([0-9.]+)", log)
        assert m is not None


class TestFindCSR:
    """Test _find_csr filename auto-detection."""

    @pytest.mark.parametrize("calc_name", list(CALCULATIONS.keys()))
    def test_find_csr_in_real_fixture(self, calc_name):
        """_find_csr should find the correct CSR file in each fixture directory."""
        info = CALCULATIONS[calc_name]
        if "new" in calc_name:
            new_name, old_name = "hrs1_nao.csr", "data-HR-sparse_SPIN0.csr"
        else:
            new_name, old_name = "hrs1_nao.csr", "data-HR-sparse_SPIN0.csr"
        result = AbacusParser._find_csr(info["out_dir"], new_name, old_name)
        assert Path(result).exists()

    def test_prefers_new_filename(self, tmp_path):
        new_f = tmp_path / "hrs1_nao.csr"
        old_f = tmp_path / "data-HR-sparse_SPIN0.csr"
        new_f.write_text("new")
        old_f.write_text("old")
        result = AbacusParser._find_csr(
            tmp_path, "hrs1_nao.csr", "data-HR-sparse_SPIN0.csr"
        )
        assert result == str(new_f)

    def test_falls_back_to_old_filename(self, tmp_path):
        old_f = tmp_path / "data-HR-sparse_SPIN0.csr"
        old_f.write_text("old")
        result = AbacusParser._find_csr(
            tmp_path, "hrs1_nao.csr", "data-HR-sparse_SPIN0.csr"
        )
        assert result == str(old_f)

    def test_raises_when_neither_exists(self, tmp_path):
        with pytest.raises(FileNotFoundError, match="Neither"):
            AbacusParser._find_csr(tmp_path, "hrs1_nao.csr", "data-HR-sparse_SPIN0.csr")


class TestReadHRSR:
    """Test the high-level read_HR_SR function on fixture data."""

    @pytest.mark.parametrize("calc_name", ["old_nspin2", "new_nspin2"])
    def test_read_hr_sr_nspin2(self, calc_name):
        """read_HR_SR for nspin=2 returns (basis, Rlist, HR_up, HR_dn, SR)."""
        info = CALCULATIONS[calc_name]
        out = info["out_dir"]
        nbasis, Rlist, HR_up, HR_dn, SR = read_HR_SR(
            nspin=2,
            HR_fileName=[
                str(out / info["hr_files"][0]),
                str(out / info["hr_files"][1]),
            ],
            SR_fileName=str(out / info["sr_file"]),
        )
        assert nbasis == info["basis"]
        assert HR_up.shape == SR.shape
        assert HR_dn.shape == SR.shape
        assert not np.allclose(HR_up, HR_dn)

    @pytest.mark.parametrize("calc_name", ["old_nspin4", "new_nspin4"])
    def test_read_hr_sr_nspin4(self, calc_name):
        """read_HR_SR for nspin=4 returns (basis, Rlist, HR, SR) with complex HR."""
        info = CALCULATIONS[calc_name]
        out = info["out_dir"]
        nbasis, Rlist, HR, SR = read_HR_SR(
            nspin=4,
            HR_fileName=str(out / info["hr_files"][0]),
            SR_fileName=str(out / info["sr_file"]),
        )
        assert nbasis == info["basis"]
        assert HR.dtype == complex
        assert SR.dtype == complex


class TestStrictNewCSR:
    @pytest.mark.parametrize("data_size", [-1, 2])
    def test_legacy_text_rejects_invalid_block_data_size_before_payload(
        self, tmp_path, data_size
    ):
        filename = tmp_path / "legacy.csr"
        filename.write_text(
            f"Matrix Dimension of H(R): 1\nMatrix number of H(R): 1\n0 0 0 {data_size}\n"
        )
        with pytest.raises(ValueError, match="data_size"):
            XR_matrix(1, filename).read_file()

    @pytest.mark.parametrize("data_size", [-1, 2])
    def test_legacy_binary_rejects_invalid_block_data_size_before_payload(
        self, tmp_path, data_size
    ):
        filename = tmp_path / "legacy.bin"
        filename.write_bytes(struct.pack("8i", 1, 1, 0, 0, 0, data_size, 0, 0))
        with pytest.raises(ValueError, match="data_size"):
            XR_matrix(1, filename).read_file_binary()

    def test_rejects_oversized_file_before_opening(self, tmp_path, monkeypatch):
        filename = tmp_path / "missing.csr"

        class Stat:
            st_size = MAX_NEW_CSR_FILE_BYTES + 1

        monkeypatch.setattr(abacus_api.os, "stat", lambda path: Stat())
        with pytest.raises(ValueError, match="MAX_NEW_CSR_FILE_BYTES"):
            read_csr_new_format(filename)

    @pytest.mark.parametrize(
        "basis, r_count",
        [(MAX_DENSE_CSR_ELEMENTS + 1, 1), (1, MAX_DENSE_CSR_ELEMENTS + 1)],
    )
    def test_rejects_oversized_dense_declarations_before_allocation(
        self, tmp_path, monkeypatch, basis, r_count
    ):
        filename = tmp_path / "hrs1_nao.csr"
        filename.write_text(
            f"{basis} # number of localized basis\n{r_count} # number of Bravais lattice vector R\n# CSR Format\n"
        )
        monkeypatch.setattr(
            abacus_api.np,
            "zeros",
            lambda *args, **kwargs: pytest.fail("dense allocation attempted"),
        )
        with pytest.raises(ValueError, match="MAX_DENSE_CSR_ELEMENTS"):
            read_csr_new_format(filename)

    def test_legacy_rejects_nonfinite_values(self, tmp_path):
        filename = tmp_path / "legacy.csr"
        filename.write_text(
            "Matrix Dimension of H(R): 1\nMatrix number of H(R): 1\n0 0 0 1\nnan\n0\n0 1\n"
        )
        with pytest.raises(ValueError, match="finite"):
            XR_matrix(1, filename).read_file()

    @pytest.mark.parametrize(
        ("text", "match"),
        [
            ("2 # number of localized basis\n", "Bravais"),
            (_strict_csr_text().replace("# CSR Format", "# CSR values"), "CSR Format"),
            (_strict_csr_text().replace("2 # number", "two # number"), "basis"),
            (_strict_csr_text().replace("0 0 0 2", "0 0 0 -1"), "nnz"),
            (_strict_csr_text().replace("(1.0,0.0)", "(bad,0.0)", 1), "values"),
            (_strict_csr_text().replace("0 1", "0 2", 1), "indices"),
            (_strict_csr_text().replace("0 0 2", "1 0 2"), "indptr"),
            (
                _strict_csr_text().replace(
                    "# CSR row_indptr\n0 0 2", "# CSR row_indptr\n0 2"
                ),
                "indptr",
            ),
            (_strict_csr_text().replace("(1.0,0.0) (1.0,0.0)", "(1.0,0.0)"), "values"),
            (
                _strict_csr_text().replace(
                    "# CSR column_indices\n0 1", "# CSR column_indices\n0"
                ),
                "indices",
            ),
            (
                _strict_csr_text().replace("# CSR values", "# CSR column_indices", 1),
                "values",
            ),
            (_strict_csr_text() + "# CSR values\n", "extra"),
            (_strict_csr_text() + "1 0 0 0\n", "extra"),
        ],
    )
    def test_rejects_malformed_input(self, tmp_path, text, match):
        filename = tmp_path / "hrs1_nao.csr"
        filename.write_text(text)
        with pytest.raises(ValueError, match=match):
            read_csr_new_format(filename)

    @pytest.mark.parametrize(
        "text",
        [
            _strict_csr_text(blocks=(((0, 0, 0), 2),)).replace(
                "1 # number of Bravais", "2 # number of Bravais"
            ),
            _strict_csr_text(blocks=(((0, 0, 0), 2), ((0, 0, 0), 2))),
            _strict_csr_text(blocks=(((0, 0, 0), 2),))[:-5],
        ],
    )
    def test_rejects_wrong_count_duplicate_and_truncation(self, tmp_path, text):
        filename = tmp_path / "hrs1_nao.csr"
        filename.write_text(text)
        with pytest.raises(ValueError):
            read_csr_new_format(filename)

    def test_accepts_zero_nnz_and_legal_comments(self, tmp_path):
        filename = tmp_path / "hrs1_nao.csr"
        filename.write_text(_strict_csr_text(blocks=(((0, 0, 0), 0),)))
        basis, rlist, matrix = read_csr_new_format(filename)
        assert basis == 2
        assert np.array_equal(rlist, [[0, 0, 0]])
        assert not matrix.any()

    def test_accepts_zero_nnz_with_canonical_empty_sections(self, tmp_path):
        text = _strict_csr_text(blocks=(((0, 0, 0), 0),))
        text += "# CSR values\n# CSR column indices\n# CSR row pointers\n0 0 0\n"
        filename = tmp_path / "hrs1_nao.csr"
        filename.write_text(text)
        _, _, matrix = read_csr_new_format(filename)
        assert not matrix.any()

    @pytest.mark.parametrize("is_complex", [False, True])
    def test_accepts_multiline_section_payloads(self, tmp_path, is_complex):
        text = _strict_csr_text()
        values = (
            "(1.0,0.0)\n# legal payload comment\n(1.0,0.0)"
            if is_complex
            else "1.0\n# legal payload comment\n1.0"
        )
        text = text.replace("(1.0,0.0) (1.0,0.0)", values)
        text = text.replace("0 1\n# CSR row_indptr", "0\n1\n# CSR row_indptr")
        text = text.replace("# CSR row_indptr\n0 0 2", "# CSR row_indptr\n0\n0 2")
        filename = tmp_path / "hrs1_nao.csr"
        filename.write_text(text)
        _, _, matrix = read_csr_new_format(filename, is_complex=is_complex)
        assert matrix.any()

    @pytest.mark.parametrize("is_complex", [False, True])
    def test_accepts_canonical_abacus_markers_with_wrapped_payloads(
        self, tmp_path, is_complex
    ):
        values = "(1.0,0.0)\n(2.0,0.0)" if is_complex else "1.0\n2.0"
        text = _strict_csr_text()
        text = text.replace("(1.0,0.0) (1.0,0.0)", values)
        text = text.replace("# CSR column_indices", "# CSR column indices")
        text = text.replace("# CSR row_indptr", "# CSR row pointers")
        text = text.replace("0 1\n# CSR row pointers", "0\n1\n# CSR row pointers")
        text = text.replace("# CSR row pointers\n0 0 2", "# CSR row pointers\n0\n0 2")
        filename = tmp_path / "hrs1_nao.csr"
        filename.write_text(text)
        _, _, matrix = read_csr_new_format(filename, is_complex=is_complex)
        assert np.allclose(matrix[0, 1], [1, 2])

    @pytest.mark.parametrize(
        "replacement",
        [
            "(1.0,0.0)\n(1.0,0.0)\n(1.0,0.0)",
            "(1.0,0.0)\n# CSR column_indices",
        ],
    )
    def test_rejects_multiline_payload_overrun_and_truncation(
        self, tmp_path, replacement
    ):
        filename = tmp_path / "hrs1_nao.csr"
        filename.write_text(
            _strict_csr_text().replace("(1.0,0.0) (1.0,0.0)", replacement)
        )
        with pytest.raises(ValueError, match="values"):
            read_csr_new_format(filename)

    def test_rejects_excess_wrapped_indices_before_next_marker(self, tmp_path):
        filename = tmp_path / "hrs1_nao.csr"
        text = _strict_csr_text().replace(
            "# CSR column_indices\n0 1", "# CSR column_indices\n0\n1\n0"
        )
        filename.write_text(text)
        with pytest.raises(ValueError, match="indices"):
            read_csr_new_format(filename)

    def test_rejects_new_text_csr_in_binary_mode(self, tmp_path):
        out = _write_single_step_soc_case(tmp_path)
        with pytest.raises(ValueError, match="binary=True.*new-format"):
            AbacusParser(outpath=out, binary=True).Read_HSR_noncollinear()

    def test_public_split_soc_parsers_are_exported(self):
        assert AbacusSingleStepSOCParser is not None
        assert AbacusSplitSOCParser is not None


class TestCSRAlignment:
    def test_alignment_reorders_equal_unique_r_sets(self):
        rlist = np.array([[0, 0, 0], [1, 0, 0]])
        other_rlist = rlist[::-1]
        matrix = np.array([[[2.0]], [[1.0]]])
        aligned = _align_csr_blocks(1, rlist, 1, other_rlist, matrix, "SR")
        assert np.array_equal(aligned[:, 0, 0], [1.0, 2.0])

    def test_sr_subset_reorders_and_zero_fills_complex_blocks(self):
        rlist = np.array([[0, 0, 0], [1, 0, 0]])
        matrix = np.array([[[2 + 0.5j]]])
        aligned = _align_csr_blocks(
            1, rlist, 1, np.array([[1, 0, 0]]), matrix, "SR", allow_missing_zero=True
        )
        assert aligned.dtype == complex
        assert np.array_equal(aligned[:, 0, 0], [0, 2 + 0.5j])

    def test_default_hr_alignment_rejects_missing_r_blocks(self):
        with pytest.raises(ValueError, match="R sets"):
            _align_csr_blocks(
                1,
                np.array([[0, 0, 0], [1, 0, 0]]),
                1,
                np.array([[0, 0, 0]]),
                np.ones((1, 1, 1)),
                "HR",
            )

    @pytest.mark.parametrize(
        ("basis", "rlist", "match"),
        [
            (2, np.array([[0, 0, 0]]), "basis"),
            (1, np.array([[0, 0, 0], [0, 0, 0]]), "duplicate"),
            (1, np.array([[2, 0, 0]]), "R sets"),
        ],
    )
    def test_alignment_rejects_basis_and_r_mismatches(self, basis, rlist, match):
        with pytest.raises(ValueError, match=match):
            _align_csr_blocks(
                1,
                np.array([[0, 0, 0]]),
                basis,
                rlist,
                np.zeros((len(rlist), basis, basis)),
                "SR",
            )


class TestWrapperNewCSR:
    def _parser(self, outpath, binary=False):
        parser = AbacusParser.__new__(AbacusParser)
        parser.outpath = outpath
        parser.binary = binary
        return parser

    def _write_noncollinear_set(
        self, outpath, hr_blocks, sr_blocks, sr_basis=2, sr_imaginary=0.0
    ):
        outpath.mkdir()
        (outpath / "hrs1_nao.csr").write_text(_new_csr_text(2, hr_blocks))
        (outpath / "srs1_nao.csr").write_text(
            _new_csr_text(sr_basis, sr_blocks, sr_imaginary)
        )

    def test_new_noncollinear_reorders_and_preserves_complex_sr(self, tmp_path):
        out = tmp_path / "OUT"
        self._write_noncollinear_set(
            out,
            [((0, 0, 0), 1), ((1, 0, 0), 2)],
            [((1, 0, 0), 20), ((0, 0, 0), 10)],
            sr_imaginary=0.25,
        )
        _, rlist, _, sr = self._parser(out).Read_HSR_noncollinear()
        assert np.array_equal(rlist, [[0, 0, 0], [1, 0, 0]])
        assert np.allclose(sr[:, 0, 0], [10 + 0.25j, 20 + 0.25j])

    def test_new_noncollinear_sr_subset_reorders_and_zero_fills(self, tmp_path):
        out = tmp_path / "OUT"
        self._write_noncollinear_set(
            out, [((0, 0, 0), 1), ((1, 0, 0), 2)], [((1, 0, 0), 20)], sr_imaginary=0.25
        )
        _, rlist, _, sr = self._parser(out).Read_HSR_noncollinear()
        assert np.array_equal(rlist, [[0, 0, 0], [1, 0, 0]])
        assert np.allclose(sr[:, 0, 0], [0, 20 + 0.25j])

    @pytest.mark.parametrize(
        ("sr_blocks", "sr_basis", "match"),
        [([((0, 0, 0), 1)], 3, "basis"), ([((2, 0, 0), 1)], 2, "R sets")],
    )
    def test_new_noncollinear_rejects_sr_basis_and_r_mismatches(
        self, tmp_path, sr_blocks, sr_basis, match
    ):
        out = tmp_path / "OUT"
        self._write_noncollinear_set(out, [((0, 0, 0), 1)], sr_blocks, sr_basis)
        with pytest.raises(ValueError, match=match):
            self._parser(out).Read_HSR_noncollinear()

    def test_text_mode_rejects_mixed_new_and_legacy_set(self, tmp_path):
        out = tmp_path / "OUT"
        out.mkdir()
        (out / "hrs1_nao.csr").write_text(_new_csr_text(2, [((0, 0, 0), 1)]))
        (out / "data-SR-sparse_SPIN0.csr").write_text("legacy")
        with pytest.raises(ValueError, match="mixed"):
            self._parser(out).Read_HSR_noncollinear()

    def test_binary_mode_uses_legacy_files_when_new_files_coexist(
        self, tmp_path, monkeypatch
    ):
        out = tmp_path / "OUT"
        out.mkdir()
        for name in [
            "hrs1_nao.csr",
            "srs1_nao.csr",
            "data-HR-sparse_SPIN0.csr",
            "data-SR-sparse_SPIN0.csr",
        ]:
            (out / name).write_text("placeholder")
        calls = []

        class Reader:
            def __init__(self, nspin, filename):
                calls.append((filename, False))
                self.basis_num = 2
                self.R_direct_coor = np.array([[0, 0, 0]])
                self.XR = np.eye(2, dtype=complex)[None]

            def read_file_binary(self):
                calls[-1] = (calls[-1][0], True)

            def read_file(self):
                raise AssertionError("binary dispatch must not use text reader")

        monkeypatch.setattr(abacus_wrapper, "XR_matrix", Reader)
        self._parser(out, binary=True).Read_HSR_noncollinear()
        assert [Path(filename).name for filename, _ in calls] == [
            "data-HR-sparse_SPIN0.csr",
            "data-SR-sparse_SPIN0.csr",
        ]
        assert all(binary for _, binary in calls)


class TestFermiParsing:
    @pytest.mark.parametrize(
        "line, expected",
        [("EFERMI = -1.25e+1 eV", -12.5), ("E_Fermi = +.5E-2 eV", 0.005)],
    )
    def test_read_efermi_accepts_signed_scientific_values(
        self, tmp_path, line, expected
    ):
        (tmp_path / "running_scf.log").write_text(line)
        parser = AbacusParser.__new__(AbacusParser)
        parser.outpath = tmp_path
        assert parser.read_efermi() == expected

    def test_read_efermi_rejects_non_numeric_value(self, tmp_path):
        (tmp_path / "running_scf.log").write_text("EFERMI = unknown eV")
        parser = AbacusParser.__new__(AbacusParser)
        parser.outpath = tmp_path
        with pytest.raises(ValueError, match="numeric"):
            parser.read_efermi()

    @pytest.mark.parametrize("line", ["EFERMI = 1.2foo", "E_Fermi = 1.2foo"])
    def test_read_efermi_rejects_numeric_prefixes(self, tmp_path, line):
        (tmp_path / "running_scf.log").write_text(line)
        parser = AbacusParser.__new__(AbacusParser)
        parser.outpath = tmp_path
        with pytest.raises(ValueError, match="numeric"):
            parser.read_efermi()

    def test_read_efermi_ignores_unrelated_mentions(self, tmp_path):
        (tmp_path / "running_scf.log").write_text(
            "EFERMI convergence is enabled\nE_Fermi = 2.5 eV"
        )
        parser = AbacusParser.__new__(AbacusParser)
        parser.outpath = tmp_path
        assert parser.read_efermi() == 2.5

    @pytest.mark.parametrize(
        ("line", "expected"),
        [
            ("E_Fermi        -0.0262701204        -0.3574233252", -0.3574233252),
            ("E_Fermi        0.0203805901         0.2772921542", 0.2772921542),
        ],
    )
    def test_read_efermi_uses_ev_column_from_abacus_table(
        self, tmp_path, line, expected
    ):
        (tmp_path / "running_scf.log").write_text(line)
        parser = AbacusParser.__new__(AbacusParser)
        parser.outpath = tmp_path
        assert parser.read_efermi() == expected

    @pytest.mark.parametrize(
        "line",
        ["EFERMI = 1.2 eV trailing", "E_Fermi 0.1", "E_Fermi 0.1 1.2foo"],
    )
    def test_read_efermi_rejects_malformed_assignment_and_table_tokens(
        self, tmp_path, line
    ):
        (tmp_path / "running_scf.log").write_text(line)
        parser = AbacusParser.__new__(AbacusParser)
        parser.outpath = tmp_path
        with pytest.raises(ValueError, match="numeric"):
            parser.read_efermi()


class TestSingleStepSOCParser:
    def test_parse_split_soc_from_single_out_dir(self, tmp_path):
        out = _write_single_step_soc_case(tmp_path)

        model = AbacusSingleStepSOCParser(outpath=out).parse()

        expected_total = np.diag([1.0, 2.0, 3.0, 4.0]) * RY_TO_EV
        expected_soc = np.diag([-0.05, 0.05, -0.05, 0.05]) * RY_TO_EV
        assert model.split_soc is True
        assert model.HR_full is not None
        assert model.HR_soc is not None
        assert model.HR_nosoc is not None
        assert np.allclose(model.HR_full[0], expected_total)
        assert np.allclose(model.HR_soc[0], expected_soc)
        assert np.allclose(chargepart(model.HR_soc[0]), 0.0)
        assert np.allclose(model.HR_nosoc[0], expected_total - expected_soc)
        assert np.array_equal(model.Rlist, np.array([[0, 0, 0]]))
        assert getattr(model, "efermi") == 1.25
        assert model.nel == 6
        assert model.atoms is not None
        assert len(model.atoms) == 1
        assert getattr(model, "basis") == [
            ("Fe1", "sZ1", "up"),
            ("Fe1", "sZ1", "down"),
            ("Fe1", "pzZ1", "up"),
            ("Fe1", "pzZ1", "down"),
        ]

    def test_missing_soc_file_fails_clearly(self, tmp_path):
        out = _write_single_step_soc_case(tmp_path)
        (out / "hrs1-soc_nao.csr").unlink()

        with pytest.raises(FileNotFoundError, match="SOC Hamiltonian file not found"):
            AbacusSingleStepSOCParser(outpath=out).parse()

    def test_soc_rlist_must_match_total_rlist(self, tmp_path):
        out = _write_single_step_soc_case(tmp_path, soc_r=(1, 0, 0))

        with pytest.raises(ValueError, match="R sets"):
            AbacusSingleStepSOCParser(outpath=out).parse()


class TestSplitSOCParser:
    def _write_case(self, root, blocks, sr_blocks=None):
        out = _write_single_step_soc_case(root)
        (out / "hrs1_nao.csr").write_text(_new_csr_text(4, blocks))
        (out / "srs1_nao.csr").write_text(
            _new_csr_text(4, blocks if sr_blocks is None else sr_blocks)
        )
        return out

    def test_parse_reorders_soc_r_blocks(self, tmp_path):
        nosoc = self._write_case(
            tmp_path / "nosoc",
            [((0, 0, 0), 1), ((1, 0, 0), 2)],
            [((0, 0, 0), 10), ((1, 0, 0), 20)],
        )
        soc = self._write_case(
            tmp_path / "soc",
            [((1, 0, 0), 4), ((0, 0, 0), 3)],
            [((1, 0, 0), 20), ((0, 0, 0), 10)],
        )
        model = AbacusSplitSOCParser(nosoc, soc).parse()
        assert np.array_equal(model.Rlist, [[0, 0, 0], [1, 0, 0]])

    def test_parse_accepts_equal_subset_sr_blocks(self, tmp_path):
        nosoc = self._write_case(
            tmp_path / "nosoc", [((0, 0, 0), 1), ((1, 0, 0), 2)], [((1, 0, 0), 20)]
        )
        soc = self._write_case(
            tmp_path / "soc", [((1, 0, 0), 4), ((0, 0, 0), 3)], [((1, 0, 0), 20)]
        )
        model = AbacusSplitSOCParser(nosoc, soc).parse()
        assert np.array_equal(model.Rlist, [[0, 0, 0], [1, 0, 0]])

    def test_parse_rejects_mismatched_overlap_values(self, tmp_path):
        nosoc = self._write_case(
            tmp_path / "nosoc", [((0, 0, 0), 1)], [((0, 0, 0), 10)]
        )
        soc = self._write_case(tmp_path / "soc", [((0, 0, 0), 2)], [((0, 0, 0), 11)])
        with pytest.raises(ValueError, match="overlap differs.*max difference"):
            AbacusSplitSOCParser(nosoc, soc).parse()

    def test_parse_rejects_soc_r_mismatch(self, tmp_path):
        nosoc = self._write_case(tmp_path / "nosoc", [((0, 0, 0), 1)])
        soc = self._write_case(tmp_path / "soc", [((2, 0, 0), 1)])
        with pytest.raises(ValueError, match="R sets"):
            AbacusSplitSOCParser(nosoc, soc).parse()

    def test_parse_rejects_equal_size_different_basis_descriptors(self, tmp_path):
        nosoc = self._write_case(tmp_path / "nosoc", [((0, 0, 0), 1)])
        soc = self._write_case(tmp_path / "soc", [((0, 0, 0), 1)])
        (soc / "Orbital").write_text(
            """  #io    spec    l    m    z  sym
    0      Fe    1    0    1             pz
    0      Fe    0    0    1              s

"""
        )
        with pytest.raises(ValueError, match="basis descriptors"):
            AbacusSplitSOCParser(nosoc, soc).parse()
