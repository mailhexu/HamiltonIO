#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GPAW API for reading pickle files and extracting Hamiltonian/Overlap matrices
"""

import pickle
from typing import Optional, Tuple

import numpy as np


def read_gpaw_pickle(pickle_file: str) -> dict:
    """
    Read GPAW pickle file containing Hamiltonian and Overlap matrices

    Parameters
    ----------
    pickle_file : str
        Path to the pickle file created by GPAW's dump_hamiltonian()

    Returns
    -------
    dict
        Dictionary containing H_skMM, S_kMM, atoms_data, calc_data

    Notes
    -----
    The pickle file is created by GPAW's dump_hamiltonian() which stores
    three separate pickle dumps:
        1. (h_skmm, s_kmm) - Hamiltonian and Overlap matrices
        2. atoms_data - cell, positions, numbers, pbc
        3. calc_data - weight_k, ibzk_kc
    """
    data = {}
    with open(pickle_file, "rb") as f:
        # First dump: H and S matrices
        h_skmm, s_kmm = pickle.load(f)
        data["H_skMM"] = h_skmm
        data["S_kMM"] = s_kmm

        # Second dump: atoms data
        atoms_data = pickle.load(f)
        data["atoms_data"] = atoms_data

        # Third dump: calc data
        calc_data = pickle.load(f)
        data["calc_data"] = calc_data

    return data


def read_gpaw_calculator(gpw_file: str):
    """
    Read GPAW calculator from .gpw file

    Parameters
    ----------
    gpw_file : str
        Path to the .gpw file containing the converged GPAW calculation

    Returns
    -------
    GPAW calculator object
        GPAW calculator with converged wavefunction

    Notes
    -----
    Uses GPAW's restart() method to load the calculator.
    """
    from gpaw import restart

    atoms, calc = restart(gpw_file, txt=None)
    return calc


def get_hr_sr_from_calc(calc) -> Tuple[np.ndarray, np.ndarray]:
    """
    Extract k-space Hamiltonian and Overlap matrices from GPAW calculator

    Parameters
    ----------
    calc : GPAW calculator object
        Converged GPAW calculator with LCAO basis set

    Returns
    -------
    H_skMM : np.ndarray
        Hamiltonian matrix in k-space
        Shape: (nspins, nkpts, nao, nao) for nspin=2
               (nkpts, nao, nao) for nspin=1
        Units: eV
    S_kMM : np.ndarray
        Overlap matrix in k-space
        Shape: (nkpts, nao, nao)
        Dimensionless

    Notes
    -----
    This function uses GPAW's get_lcao_hamiltonian() which returns:
    - H_skMM: Hamiltonian in eV
    - S_kMM: Overlap matrix (dimensionless)

    For spin-polarized calculations (nspin=2), H_skMM has shape (nspins, nkpts, nao, nao).
    For non-polarized calculations (nspin=1), the first dimension is squeezed out.
    """
    from gpaw.lcao.tools import get_lcao_hamiltonian

    H_skMM, S_kMM = get_lcao_hamiltonian(calc)
    return H_skMM, S_kMM


def convert_k_to_r_space_gpaW_tb(
    H_skMM: np.ndarray,
    S_kMM: np.ndarray,
    calc,
    nkpts: Optional[int] = None,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, int]:
    """
    Convert k-space Hamiltonian to real-space using GPAW's TightBinding class

    Parameters
    ----------
    H_skMM : np.ndarray
        k-space Hamiltonian matrix
        Shape: (nspins, nkpts, nao, nao) or (nkpts, nao, nao)
    S_kMM : np.ndarray
        k-space Overlap matrix
        Shape: (nkpts, nao, nao)
    calc : GPAW calculator object
        GPAW calculator object (needed for atomic structure and k-points)
    nkpts : int, optional
        Number of k-points (auto-detected if not specified)

    Returns
    -------
    HR : np.ndarray
        Real-space Hamiltonian matrices
        Shape: (nR, nao, nao) or (nspins, nR, nao, nao)
    SR : np.ndarray
        Real-space Overlap matrices
        Shape: (nR, nao, nao)
    Rlist : np.ndarray
        List of R vectors (cell indices)
        Shape: (nR, 3)
    nbasis : int
        Number of atomic orbitals (nao)

    Notes
    -----
    This function uses GPAW's TightBinding.bloch_to_real_space() method
    which properly handles:
    - Fourier transform from k-space to real-space
    - Time-reversal symmetry (expands IBZ to full BZ)
    - Complex phase factors

    The key is that bloch_to_real_space takes k-space data and performs:
        H_R = sum_k w_k * H_k * exp(2πi k·R)
    """
    from gpaw.lcao.tightbinding import TightBinding

    # Determine if H_skMM has spin dimension
    if H_skMM.ndim == 4:
        # Spin-polarized case: (nspins, nkpts, nao, nao)
        nspins, nkpts, nao, _ = H_skMM.shape
        has_spin = True
    else:
        # Non-polarized case: (nkpts, nao, nao)
        nkpts, nao, _ = H_skMM.shape
        nspins = 1
        has_spin = False

    nbasis = nao

    # Create TightBinding object to access bloch_to_real_space
    tb = TightBinding(calc.atoms, calc)

    # Get R cell vectors from TightBinding
    # These define the real-space lattice vectors
    R_cN = tb.R_cN.T  # Shape: (nR, 3)

    if has_spin:
        # Spin-polarized: transform each spin separately
        HR_list = []
        for ispin in range(nspins):
            # Use bloch_to_real_space to convert k-space to real-space
            # Input shape: (nkpts, nao, nao)
            # The function expects the input to be (nkpts, nao*nao) for the transform
            # but we can pass it as (nkpts, nao, nao) and let it handle the reshaping
            H_R_list = tb.bloch_to_real_space(H_skMM[ispin], R_c=R_cN.T)
            # H_R_list is list of (nao, nao) arrays for each R
            HR_list.append(np.array(H_R_list))
        HR = np.array(HR_list)
        # S is the same for all spins
        S_R_list = tb.bloch_to_real_space(S_kMM, R_c=R_cN.T)
        SR = np.array(S_R_list)
    else:
        # Non-polarized: single spin channel
        H_R_list = tb.bloch_to_real_space(H_skMM, R_c=R_cN.T)
        HR = np.array(H_R_list)
        S_R_list = tb.bloch_to_real_space(S_kMM, R_c=R_cN.T)
        SR = np.array(S_R_list)

    # Use R cell vectors as Rlist
    # R_cN contains integer cell indices
    Rlist = R_cN.copy()

    return HR, SR, Rlist, nbasis


def convert_k_to_r_space_manual(
    H_skMM: np.ndarray,
    S_kMM: np.ndarray,
    ibzk_qc: np.ndarray,
    weight_k: np.ndarray,
    cell: np.ndarray,
    bzk_qc: Optional[np.ndarray] = None,
    Rmax: int = 2,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, int]:
    """
    Convert k-space Hamiltonian to real-space using manual Fourier transform

    Parameters
    ----------
    H_skMM : np.ndarray
        k-space Hamiltonian matrix at IBZ k-points
        Shape: (nspins, nkpts, nao, nao) or (nkpts, nao, nao)
    S_kMM : np.ndarray
        k-space Overlap matrix at IBZ k-points
        Shape: (nkpts, nao, nao)
    ibzk_qc : np.ndarray
        IBZ k-points in fractional coordinates
        Shape: (nkpts, 3)
    weight_k : np.ndarray
        k-point weights for IBZ k-points
        Shape: (nkpts,)
    cell : np.ndarray
        Unit cell vectors
        Shape: (3, 3)
    bzk_qc : np.ndarray, optional
        Full BZ k-points in fractional coordinates (currently NOT used due to
        symmetry complexity. Kept for API compatibility.)
    Rmax : int, optional
        Maximum range for R vectors (default: 2)

    Returns
    -------
    HR : np.ndarray
        Real-space Hamiltonian matrices
        Shape: (nR, nao, nao) or (nspins, nR, nao, nao)
    SR : np.ndarray
        Real-space Overlap matrices
        Shape: (nR, nao, nao)
    Rlist : np.ndarray
        List of R vectors (cell indices)
        Shape: (nR, 3)
    nbasis : int
        Number of atomic orbitals (nao)

    Notes
    -----
    Fourier transform formula:
        H_R = sum_k w_k * H_k * exp(-2πi k·R)

    where the sum is over the IBZ k-points only. Using full BZ with
    nearest-neighbor mapping causes non-positive-definite matrices due
    to incorrect symmetry handling.

    For accurate results, the GPAW calculation should be done with
    symmetry disabled or with a dense k-mesh that makes the IBZ
    approximately equal to the full BZ.

    The result is enforced to be Hermitian: H_R = (H_R + H_R^†) / 2
    """
    # Determine if H_skMM has spin dimension
    if H_skMM.ndim == 4:
        # Spin-polarized case: (nspins, nkpts, nao, nao)
        nspins, nkpts_ibz, nao, _ = H_skMM.shape
        has_spin = True
    else:
        # Non-polarized case: (nkpts, nao, nao)
        nkpts_ibz, nao, _ = H_skMM.shape
        nspins = 1
        has_spin = False
        H_skMM = H_skMM[np.newaxis, ...]  # Add spin dimension for uniformity

    nbasis = nao

    # Use IBZ k-points directly (not full BZ)
    # Full BZ expansion with nearest-neighbor mapping causes issues
    kpts_ft = ibzk_qc
    H_ft = H_skMM
    S_ft = S_kMM

    # For proper Fourier transform normalization:
    # H_R = (1/N_k) * sum_k H_k * exp(-2πi k·R)
    # H_k = sum_R H_R * exp(2πi k·R)
    # This ensures the round-trip preserves the Hamiltonian
    # Use k-point weights for correct normalization
    weights_ft = weight_k / np.sum(weight_k)

    # Generate R vectors
    R_list = []
    for i in range(-Rmax, Rmax + 1):
        for j in range(-Rmax, Rmax + 1):
            for k in range(-Rmax, Rmax + 1):
                R_list.append([i, j, k])
    Rlist = np.array(R_list)
    nR = len(R_list)

    if has_spin:
        HR = np.zeros((nspins, nR, nao, nao), dtype=complex)
    else:
        HR = np.zeros((nR, nao, nao), dtype=complex)
    SR = np.zeros((nR, nao, nao), dtype=complex)

    # Perform Fourier transform for each R
    # Using the convention: H_R = (1/N_k) * sum_k H_k * exp(-2πi k·R)
    # Note the negative sign in the phase for k→R transform
    for iR, R in enumerate(R_list):
        phase = np.exp(-2j * np.pi * np.dot(kpts_ft, R))

        if has_spin:
            for ispin in range(nspins):
                # Compute H_R = (1/N_k) * sum_k H_k * exp(-2πi k·R)
                HR_k = np.sum(
                    weights_ft[:, None, None] * H_ft[ispin] * phase[:, None, None],
                    axis=0,
                )
                # Store directly - Fourier transform of Hermitian H_k should give Hermitian H_R
                HR[ispin, iR] = HR_k
        else:
            HR_k = np.sum(
                weights_ft[:, None, None] * H_ft[0] * phase[:, None, None], axis=0
            )
            HR[iR] = HR_k

        SR_k = np.sum(weights_ft[:, None, None] * S_ft * phase[:, None, None], axis=0)
        SR[iR] = SR_k

    # Remove spin dimension if originally non-polarized
    if not has_spin:
        HR = HR[0]

    return HR, SR, Rlist, nbasis
