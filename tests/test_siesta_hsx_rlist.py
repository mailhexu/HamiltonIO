"""Regression: SislParser.read_Rlist for k-sampled .HSX archives.

`read_Rlist` must derive the real-space interaction shells from the
Hamiltonian's own lattice (`ham.lattice.sc_off`), not from the fdf
geometry (whose `sc_off` is always the single [[0, 0, 0]] row, so a
k-sampled .HSX with N stored R-blocks crashes `get_model` with a numpy
reshape ValueError).
"""

import os
from pathlib import Path

import numpy as np
import pytest

sisl = pytest.importorskip("sisl")

FIXTURE = Path(
    os.environ.get(
        "HAMILTONIO_SIESTA_HSX_FIXTURE",
        "/home/hexu/projects/unfolding_dev/unfolding/tests/data/si_example/si_sc.fdf",
    )
)


def test_read_rlist_k_sampled_hsx():
    if not FIXTURE.is_file():
        pytest.skip(f"k-sampled si_sc fixture not found: {FIXTURE}")
    from HamiltonIO.siesta.sisl_wrapper import SislParser

    parser = SislParser(str(FIXTURE))
    sc_off = parser.ham.lattice.sc_off
    assert len(sc_off) > 1, "fixture is not k-sampled (single R block)"

    model = parser.get_model()  # must not raise the reshape ValueError
    rlist = model.Rlist
    assert rlist.shape == (len(sc_off), 3)
    # the shells are exactly the Hamiltonian lattice offsets, R=0 included
    assert sorted(map(tuple, np.asarray(rlist))) == sorted(
        map(tuple, np.asarray(sc_off))
    )
