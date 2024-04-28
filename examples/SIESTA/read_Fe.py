from pathlib import Path
import numpy as np
from HamiltonIO import SiestaHam

def test_SiestaHam():
    fdf_file = Path('bccFe/collinear/siesta.fdf')
    ham=SiestaHam(fdf_file, spin=0)
    
    print(f"R   = {ham.Rlist}")
    print(ham.ham.spin)
    HRs = ham.get_HRs()
    SRs = ham.get_SRs()
    print(HRs[0, 0, :, :].diagonal())
    print("Sdiag:", SRs[0, :, :].diagonal())
    print("Sdiag:", SRs[0, :, :])
    print(np.sum(SRs, axis=0).diagonal())

    print(ham.Sk([0,0,0]).diagonal())

if __name__ == '__main__':
    test_SiestaHam()


