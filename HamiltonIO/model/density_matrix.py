import numpy as np
from HamiltonIO.occupations import Occupation, get_occupation


def get_density_matrix_k(
    evals=None, evecs=None, kweights=None, nel=None, width=0.1, S=None
):
    """
    Compute the density matrix in k-space.
    """
    occ = get_occupation(evals, kweights, nel, width=width)
    rho = np.einsum("kib, kb, kjb -> kij", evecs, occ, evecs.conj())
    for ik in range(len(rho)):
        rho[ik] = rho[ik] @ S[ik]
    return rho


def get_total_density_matrix(rho, kweights):
    """
    Compute the total density matrix."""
    return np.einsum("kij, k -> ij", rho, kweights)
