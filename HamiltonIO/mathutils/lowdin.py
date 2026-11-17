import numpy as np
from scipy.linalg import inv, eigh

def Lowdin(S):
    """
    Calculate S^(-1/2).
    Which is used in lowind's symmetric orthonormalization.
    psi_prime = S^(-1/2) psi
    """
    eigval, eigvec = eigh(S)
    S_half=eigvec @ np.diag(np.sqrt(1.0 / eigval)) @ (eigvec.T.conj())
    return S_half


def Lowdin_symmetric_orthonormalization(H, S):
    """
    Lowdin's symmetric orthonormalization, H' = S^(-1/2) H S^(-1/2).

    S^(-1/2) is computed from the eigendecomposition of the Hermitian positive-definite
    overlap matrix (Lowdin). This is the matrix that the general fractional matrix power
    returns, but it is several times faster for large bases, and orthonormal and Hermitian
    to machine precision also for an ill-conditioned S. The anti-Hermitian parts of S and
    of the result, which are rounding noise of the input, are removed.
    """
    S_half = Lowdin(0.5 * (S + S.conj().T))
    H_prime = S_half @ H @ S_half
    return 0.5 * (H_prime + H_prime.conj().T)
