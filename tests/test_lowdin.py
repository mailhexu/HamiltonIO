"""Tests for Lowdin's symmetric orthonormalization.

H' = S^(-1/2) H S^(-1/2) must have the generalized eigenvalues of (H, S), be Hermitian, and
S^(-1/2) must orthonormalize S, also for an ill-conditioned overlap matrix. For a
well-conditioned S the result equals the fractional matrix power S^(-1/2) H S^(-1/2).
"""

import numpy as np
import pytest
from scipy.linalg import eigh, fractional_matrix_power

from HamiltonIO.mathutils.lowdin import Lowdin, Lowdin_symmetric_orthonormalization


def _random_problem(n, smin, seed):
    """Random complex Hermitian H and Hermitian positive-definite S with smallest eigenvalue smin."""
    rng = np.random.default_rng(seed)
    a = rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n))
    H = a + a.conj().T
    q, _ = np.linalg.qr(rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n)))
    s = np.geomspace(smin, 2.0, n)
    S = (q * s) @ q.conj().T
    return H, S


@pytest.mark.parametrize("smin", [0.5, 1e-5])
def test_generalized_eigenvalues_and_orthonormality(smin):
    H, S = _random_problem(60, smin, seed=1)
    Hp = Lowdin_symmetric_orthonormalization(H, S)
    assert np.array_equal(Hp, Hp.conj().T)
    ref = eigh(H, S, eigvals_only=True)
    np.testing.assert_allclose(
        np.linalg.eigvalsh(Hp), ref, rtol=1e-9, atol=1e-9 * np.abs(ref).max()
    )
    X = Lowdin(S)
    np.testing.assert_allclose(X @ S @ X, np.eye(len(S)), atol=1e-9)


def test_matches_fractional_matrix_power():
    H, S = _random_problem(40, 0.3, seed=2)
    X = fractional_matrix_power(S, -0.5)
    np.testing.assert_allclose(
        Lowdin_symmetric_orthonormalization(H, S), X @ H @ X, rtol=1e-10, atol=1e-10
    )


def test_real_input():
    H, S = _random_problem(30, 0.2, seed=3)
    H, S = H.real + H.real.T, S.real
    Hp = Lowdin_symmetric_orthonormalization(H, S)
    np.testing.assert_allclose(
        np.linalg.eigvalsh(Hp), eigh(H, S, eigvals_only=True), atol=1e-10
    )
