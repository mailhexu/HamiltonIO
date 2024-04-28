from pathlib import Path
import numpy as np
from HamiltonIO import SiestaHam

def test_SiestaHam_collinear():
    # Read the Hamiltonian from the Siesta output
    fdf_file = Path('bccFe/collinear/siesta.fdf')
    ham=SiestaHam(fdf_file, spin=0)

    # The list of R vectors
    Rlist   = ham.Rlist

    # The Hamiltonian in the format of H(spin, R, i, j)
    # where spin is index of spin 0/1 in the case of collinear spin.
    # R is the index of the R vector, and i, j are the orbital indices
    HRs = ham.get_HR_all()

    # THe overlap matrix in the format of S(R, i, j)
    SRs = ham.get_SR_all()

    # print the diagonal of the Hamilton
    #print("Diagonal of HR spin up, R=(0,0,0)", HRs[0, 0, :, :].diagonal())

    #print("Diagonal of overlap matrix at R=(0,0,0)", SRs[0, :, :].diagonal())

    # The Hk at k=(0,0,0)
    Hk0 = ham.Hk([0,0,0])

    # In the k-space, Hk at k=(0,0,0)
    Sk0= ham.Sk([0,0,0])

def test_SiestaHam_nonpolarized():
    # Read the Hamiltonian from the Siesta output
    fdf_file = Path('bccFe/nonpolarized/siesta.fdf')
    ham=SiestaHam(fdf_file, spin=0)

    # The list of R vectors
    Rlist   = ham.Rlist

    # The Hamiltonian in the format of H( R, i, j)
    HRs = ham.get_HR_all()

    # THe overlap matrix in the format of S(R, i, j)
    SRs = ham.get_SR_all()
    
    # Hk at k=(0,0,0)
    Hk0 = ham.Hk([0,0,0])

    # Sk at k=(0,0,0)
    Sk0=ham.Sk([0,0,0])


def test_SiestaHam_SOC():
    # Read the Hamiltonian from the Siesta output
    fdf_file = Path('bccFe/SOC/siesta.fdf')
    ham=SiestaHam(fdf_file)

    # The list of R vectors
    Rlist = ham.Rlist

    # The Hamiltonian in the format of H(R, i, j)
    HRs = ham.get_HR_all()
    
    # The overlap matrix in the format of S(R, i, j)
    SRs = ham.get_SR_all()

    #Hk at k=(0,0,0)
    Hk0 = ham.Hk([0,0,0])

    #Sk at k=(0,0,0)
    Sk0=ham.Sk([0,0,0])


if __name__ == '__main__':
    test_SiestaHam_collinear()
    test_SiestaHam_nonpolarized()
    test_SiestaHam_SOC()


