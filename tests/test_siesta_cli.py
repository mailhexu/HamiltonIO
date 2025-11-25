#!/usr/bin/env python3
"""
Unit tests for SIESTA CLI functionality.

Purpose:
    Test the command-line interface and programmatic API for SIESTA intra-atomic analysis.

How to run:
    uv run python -m pytest tests/test_siesta_cli.py -v
"""

import os
import tempfile
import pytest
from pathlib import Path


@pytest.fixture
def bcc_fe_soc_path():
    """Path to BCC Fe SOC example."""
    base = Path(__file__).parent.parent
    fdf_path = base / "examples" / "SIESTA" / "bccFe" / "SOC" / "siesta.fdf"
    if not fdf_path.exists():
        pytest.skip("BCC Fe SOC example data not found")
    return str(fdf_path)


@pytest.fixture
def bcc_fe_nonpol_path():
    """Path to BCC Fe nonpolarized example."""
    base = Path(__file__).parent.parent
    fdf_path = base / "examples" / "SIESTA" / "bccFe" / "nonpolarized" / "siesta.fdf"
    if not fdf_path.exists():
        pytest.skip("BCC Fe nonpolarized example data not found")
    return str(fdf_path)


def test_analyze_intra_atomic_function_basic(bcc_fe_nonpol_path):
    """Test basic usage of analyze_intra_atomic function."""
    from HamiltonIO.siesta.cli import analyze_intra_atomic

    # Analyze without output file (should not raise error)
    ham = analyze_intra_atomic(
        fdf_file=bcc_fe_nonpol_path,
        atoms=None,
        read_soc=False,
        pauli_decomp=True,
        show_matrix=False,
        ispin="merge",
        output_file=None,
    )

    # Check returned Hamiltonian
    assert ham is not None
    assert ham.nbasis > 0
    assert len(ham.atoms) > 0


def test_analyze_intra_atomic_function_with_output(bcc_fe_nonpol_path):
    """Test analyze_intra_atomic with output file."""
    from HamiltonIO.siesta.cli import analyze_intra_atomic

    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".txt") as f:
        output_file = f.name

    try:
        ham = analyze_intra_atomic(
            fdf_file=bcc_fe_nonpol_path,
            atoms=[0],
            read_soc=False,
            pauli_decomp=True,
            show_matrix=True,
            ispin="merge",
            output_file=output_file,
        )

        # Check file was created
        assert os.path.exists(output_file)
        assert os.path.getsize(output_file) > 0

        # Check content
        with open(output_file, "r") as f:
            content = f.read()
            assert "Intra-Atomic Hamiltonian Analysis" in content
            assert "Atom 0" in content
            assert "Fe" in content

    finally:
        if os.path.exists(output_file):
            os.remove(output_file)


def test_analyze_intra_atomic_function_atom_filtering(bcc_fe_nonpol_path):
    """Test analyze_intra_atomic with atom filtering."""
    from HamiltonIO.siesta.cli import analyze_intra_atomic

    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".txt") as f:
        output_file = f.name

    try:
        ham = analyze_intra_atomic(
            fdf_file=bcc_fe_nonpol_path,
            atoms=[0],  # Only first atom
            read_soc=False,
            pauli_decomp=False,
            show_matrix=False,
            ispin="merge",
            output_file=output_file,
        )

        # Check output contains only atom 0
        with open(output_file, "r") as f:
            content = f.read()
            assert "Atom 0" in content
            # Should not have Atom 1 if there are multiple atoms
            if len(ham.atoms) > 1:
                assert "Atom 1" not in content

    finally:
        if os.path.exists(output_file):
            os.remove(output_file)


def test_analyze_intra_atomic_soc(bcc_fe_soc_path):
    """Test analyze_intra_atomic with SOC system."""
    from HamiltonIO.siesta.cli import analyze_intra_atomic

    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".txt") as f:
        output_file = f.name

    try:
        ham = analyze_intra_atomic(
            fdf_file=bcc_fe_soc_path,
            atoms=[0],
            read_soc=False,
            pauli_decomp=True,
            show_matrix=True,
            ispin="merge",
            output_file=output_file,
        )

        # Check content has Pauli decomposition
        with open(output_file, "r") as f:
            content = f.read()
            assert "Pauli Decomposition" in content
            # Should have spin components
            assert "spin" in content.lower() or "σ" in content

    finally:
        if os.path.exists(output_file):
            os.remove(output_file)


def test_cli_main_function_help():
    """Test CLI main function with --help."""
    from HamiltonIO.siesta.cli import main
    import sys

    # Save original argv
    original_argv = sys.argv

    try:
        # Test --help
        sys.argv = ["hamiltonio-siesta", "--help"]
        with pytest.raises(SystemExit) as excinfo:
            main()
        # --help should exit with code 0
        assert excinfo.value.code == 0

    finally:
        # Restore original argv
        sys.argv = original_argv


def test_cli_intra_atomic_command_help():
    """Test CLI intra-atomic subcommand with --help."""
    from HamiltonIO.siesta.cli import main
    import sys

    original_argv = sys.argv

    try:
        sys.argv = ["hamiltonio-siesta", "intra-atomic", "--help"]
        with pytest.raises(SystemExit) as excinfo:
            main()
        assert excinfo.value.code == 0

    finally:
        sys.argv = original_argv


def test_cli_version():
    """Test CLI --version flag."""
    from HamiltonIO.siesta.cli import main
    import sys

    original_argv = sys.argv

    try:
        sys.argv = ["hamiltonio-siesta", "--version"]
        with pytest.raises(SystemExit) as excinfo:
            main()
        # Version should exit with code 0
        assert excinfo.value.code == 0

    finally:
        sys.argv = original_argv


def test_cli_intra_atomic_command_basic(bcc_fe_nonpol_path):
    """Test CLI intra-atomic command basic usage."""
    from HamiltonIO.siesta.cli import main
    import sys

    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".txt") as f:
        output_file = f.name

    original_argv = sys.argv

    try:
        sys.argv = [
            "hamiltonio-siesta",
            "intra-atomic",
            bcc_fe_nonpol_path,
            "-o",
            output_file,
            "--atoms",
            "0",
        ]
        main()

        # Check output file
        assert os.path.exists(output_file)
        with open(output_file, "r") as f:
            content = f.read()
            assert "Atom 0" in content

    finally:
        sys.argv = original_argv
        if os.path.exists(output_file):
            os.remove(output_file)


def test_cli_invalid_fdf_file():
    """Test CLI with invalid FDF file."""
    from HamiltonIO.siesta.cli import main
    import sys

    original_argv = sys.argv

    try:
        sys.argv = [
            "hamiltonio-siesta",
            "intra-atomic",
            "/nonexistent/path/siesta.fdf",
        ]
        with pytest.raises(SystemExit):
            main()

    finally:
        sys.argv = original_argv


if __name__ == "__main__":
    # Run tests
    pytest.main([__file__, "-v"])
