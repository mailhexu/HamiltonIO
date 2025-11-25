#!/usr/bin/env python3
"""
Unit tests for ABACUS CLI functionality.

Purpose:
    Test the command-line interface and programmatic API for ABACUS intra-atomic analysis.

How to run:
    uv run python -m pytest tests/test_abacus_cli.py -v
"""

import os
import tempfile
import pytest
from pathlib import Path


@pytest.fixture
def fe_collinear_path():
    """Path to Fe collinear (no SOC) example."""
    out_path = Path(
        "/Users/hexu/projects/TB2J_abacus/abacus-tb2j-master/abacus_example/case_Fe/1_no_soc/OUT.Fe"
    )
    if not out_path.exists():
        pytest.skip("Fe collinear example data not found")
    return str(out_path)


@pytest.fixture
def fe_soc_path():
    """Path to Fe noncollinear (with SOC) example."""
    out_path = Path(
        "/Users/hexu/projects/TB2J_abacus/abacus-tb2j-master/abacus_example/case_Fe/2_soc/OUT.Fe"
    )
    if not out_path.exists():
        pytest.skip("Fe SOC example data not found")
    return str(out_path)


def test_analyze_intra_atomic_function_basic(fe_collinear_path):
    """Test basic usage of analyze_intra_atomic function."""
    from HamiltonIO.abacus.cli import analyze_intra_atomic

    # Analyze without output file (should not raise error)
    ham = analyze_intra_atomic(
        outpath=fe_collinear_path,
        atoms=None,
        pauli_decomp=True,
        show_matrix=False,
        spin=None,
        binary=False,
        output_file=None,
    )

    # Check returned Hamiltonian
    assert ham is not None
    assert ham.nbasis > 0
    assert len(ham.atoms) > 0


def test_analyze_intra_atomic_function_with_output(fe_collinear_path):
    """Test analyze_intra_atomic with output file."""
    from HamiltonIO.abacus.cli import analyze_intra_atomic

    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".txt") as f:
        output_file = f.name

    try:
        ham = analyze_intra_atomic(
            outpath=fe_collinear_path,
            atoms=[0],
            pauli_decomp=True,
            show_matrix=True,
            spin=None,
            binary=False,
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


def test_analyze_intra_atomic_function_atom_filtering(fe_collinear_path):
    """Test analyze_intra_atomic with atom filtering."""
    from HamiltonIO.abacus.cli import analyze_intra_atomic

    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".txt") as f:
        output_file = f.name

    try:
        ham = analyze_intra_atomic(
            outpath=fe_collinear_path,
            atoms=[0],  # Only first atom
            pauli_decomp=False,
            show_matrix=False,
            spin=None,
            binary=False,
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


def test_analyze_intra_atomic_noncollinear(fe_soc_path):
    """Test analyze_intra_atomic with noncollinear (SOC) system."""
    from HamiltonIO.abacus.cli import analyze_intra_atomic

    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".txt") as f:
        output_file = f.name

    try:
        ham = analyze_intra_atomic(
            outpath=fe_soc_path,
            atoms=[0],
            pauli_decomp=True,
            show_matrix=True,
            spin=None,
            binary=False,
            output_file=output_file,
        )

        # Check content has Pauli decomposition
        with open(output_file, "r") as f:
            content = f.read()
            assert "Pauli Decomposition" in content

    finally:
        if os.path.exists(output_file):
            os.remove(output_file)


def test_cli_main_function_help():
    """Test CLI main function with --help."""
    from HamiltonIO.abacus.cli import main
    import sys

    # Save original argv
    original_argv = sys.argv

    try:
        # Test --help
        sys.argv = ["hamiltonio-abacus", "--help"]
        with pytest.raises(SystemExit) as excinfo:
            main()
        # --help should exit with code 0
        assert excinfo.value.code == 0

    finally:
        # Restore original argv
        sys.argv = original_argv


def test_cli_intra_atomic_command_help():
    """Test CLI intra-atomic subcommand with --help."""
    from HamiltonIO.abacus.cli import main
    import sys

    original_argv = sys.argv

    try:
        sys.argv = ["hamiltonio-abacus", "intra-atomic", "--help"]
        with pytest.raises(SystemExit) as excinfo:
            main()
        assert excinfo.value.code == 0

    finally:
        sys.argv = original_argv


def test_cli_version():
    """Test CLI --version flag."""
    from HamiltonIO.abacus.cli import main
    import sys

    original_argv = sys.argv

    try:
        sys.argv = ["hamiltonio-abacus", "--version"]
        with pytest.raises(SystemExit) as excinfo:
            main()
        # Version should exit with code 0
        assert excinfo.value.code == 0

    finally:
        sys.argv = original_argv


def test_cli_intra_atomic_command_basic(fe_collinear_path):
    """Test CLI intra-atomic command basic usage."""
    from HamiltonIO.abacus.cli import main
    import sys

    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".txt") as f:
        output_file = f.name

    original_argv = sys.argv

    try:
        sys.argv = [
            "hamiltonio-abacus",
            "intra-atomic",
            "--outpath",
            fe_collinear_path,
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


def test_cli_invalid_outpath():
    """Test CLI with invalid OUT directory."""
    from HamiltonIO.abacus.cli import main
    import sys

    original_argv = sys.argv

    try:
        sys.argv = [
            "hamiltonio-abacus",
            "intra-atomic",
            "--outpath",
            "/nonexistent/path/OUT.Fe",
        ]
        with pytest.raises(SystemExit):
            main()

    finally:
        sys.argv = original_argv


def test_cli_no_input_error():
    """Test CLI with no input directories specified."""
    from HamiltonIO.abacus.cli import main
    import sys

    original_argv = sys.argv

    try:
        sys.argv = [
            "hamiltonio-abacus",
            "intra-atomic",
        ]
        with pytest.raises(SystemExit):
            main()

    finally:
        sys.argv = original_argv


if __name__ == "__main__":
    # Run tests
    pytest.main([__file__, "-v"])
