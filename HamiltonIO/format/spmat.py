"""
Sparse format of the Hamiltonian matrix and the overlap matrix.
"""
from dataclasses import dataclass
from scipy.sparse import csr_matrix, lil_matrix, spmatrix
from typing import Union

@dataclass
class M_sparse:
    norb: int
    norb_sc: int
    nR: int
    h: spmatrix
    R: Union[None, np.ndarray, list] = None
    
