import numpy as np
import pytest

from HamiltonIO.builder import TightBindingBuilder


def _builder():
    builder = TightBindingBuilder.from_data(
        cell=np.eye(3) * 5.0,
        positions=np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]),
        symbols=["Fe", "O"],
    )
    builder.add_orbitals(atom=0, labels=["s"])
    builder.add_orbitals(atom=1, labels=["px", "py", "pz"])
    return builder


def test_orbital_order_and_selectors():
    data = _builder().build()

    assert [(orb.iatom, orb.label, orb.l) for orb in data.orbitals] == [
        (0, "s", 0),
        (1, "px", 1),
        (1, "py", 1),
        (1, "pz", 1),
    ]
    assert data.orbital_index[(1, "py")] == 2
    np.testing.assert_array_equal(data.select_orbitals(atom=1), [1, 2, 3])
    np.testing.assert_array_equal(data.select_orbitals(species="O", shell=1), [1, 2, 3])
    np.testing.assert_array_equal(data.select_orbitals((1, "px")), [1])
    np.testing.assert_array_equal(data.select_orbitals(index=0), [0])


def test_invalid_selectors_and_labels_raise():
    builder = TightBindingBuilder.from_data(
        cell=np.eye(3),
        positions=np.zeros((1, 3)),
        symbols=["H"],
    )
    with pytest.raises(ValueError, match="unsupported orbital label"):
        builder.add_orbitals(atom=0, labels=["bad"])
    builder.add_orbitals(atom=0, labels=["s"])
    data = builder.build()
    with pytest.raises(ValueError, match="matched no orbitals"):
        data.select_orbitals(species="X")
    with pytest.raises(ValueError, match="out of range"):
        data.select_orbitals(index=3)


def test_onsite_conflict_and_additive_shift():
    builder = _builder()
    builder.set_onsite(species="Fe", value=1.0)
    builder.set_onsite((0, "s"), value=0.5, mode="add")
    data = builder.build()

    H0 = data.hamiltonian.get_H0()
    assert H0[0, 0] == pytest.approx(1.5)

    conflict = _builder()
    conflict.set_onsite((0, "s"), value=1.0)
    conflict.set_onsite((0, "s"), value=2.0)
    with pytest.raises(ValueError, match="conflicting onsite"):
        conflict.build()


def test_hopping_hermitian_counterpart_and_complex_conjugation():
    builder = _builder()
    builder.add_hopping((0, "s"), (1, "px"), R=(1, 0, 0), value=1.0 + 2.0j)
    data = builder.build()
    ham = data.hamiltonian

    iR = ham.get_Ridx((1, 0, 0))
    iR_neg = ham.get_Ridx((-1, 0, 0))
    assert ham.HR[iR, 0, 1] == pytest.approx(1.0 + 2.0j)
    assert ham.HR[iR_neg, 1, 0] == pytest.approx(1.0 - 2.0j)
    assert ham.SR[ham.get_Ridx((0, 0, 0)), 0, 0] == pytest.approx(1.0)


def test_spin_duplication_uses_spin_interleaved_order():
    builder = _builder()
    builder.set_onsite(species="Fe", value=-1.0)
    builder.add_hopping((0, "s"), (1, "px"), R=(0, 0, 0), value=-0.5)
    data = builder.build(nspin=2, nel=4)
    ham = data.hamiltonian
    H0 = ham.get_H0()

    assert ham.nbasis == 8
    assert ham.nspin == 2
    assert len(ham.orbs) == 8
    assert ham.orbs[0].label == "s"
    assert ham.orbs[1].label == "s"
    assert H0[0, 0] == pytest.approx(-1.0)
    assert H0[1, 1] == pytest.approx(-1.0)
    assert H0[0, 2] == pytest.approx(-0.5)
    assert H0[1, 3] == pytest.approx(-0.5)
    assert H0[2, 0] == pytest.approx(-0.5)
    assert H0[3, 1] == pytest.approx(-0.5)


def test_hubbard_dict_and_initial_density_helpers():
    data = _builder().build(nspin=2, nel=4)

    assert data.hubbard_dict({("Fe", 0): {"U": 4.0, "J": 1.0}}) == {
        "Fe": {"U": 4.0, "J": 1.0, "L": 0}
    }
    rho = data.initial_density(moment=0.1)
    assert rho.shape == (2, 4, 4)
    assert np.trace(rho[0] + rho[1]) == pytest.approx(4.0)


def test_initial_density_requires_electron_count():
    data = _builder().build(nspin=2)

    with pytest.raises(ValueError, match="requires nel"):
        data.initial_density()
