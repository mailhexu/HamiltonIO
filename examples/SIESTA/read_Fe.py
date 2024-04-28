from pathlib import Path
import numpy as np
from HamiltonIO import SiestaHam

def test_SiestaHam():
    fdf_file = Path('bccFe/siesta.fdf')
    ham=SiestaHam(fdf_file, spin=0)
    geom = ham.ham.geometry
    R = geom.lattice.sc_off
    mat=ham.ham._csr
    print("nsc:", ham.ham.nsc)
    print(f"shape: ", mat.shape)
    print(f"dtype: ", mat.dtype)
    print(f"R   = {R}")
    nrow, ncol, nspin = mat.shape
    #print(""nrow, ncol/nrow, 9*9*9)

    print(ham.Rlist)
    print(ham.ham.spin)

    HRs = ham.get_HRs()
    SRs = ham.get_SRs()
    print(HRs[0, 0, :, :].diagonal())
    print("Sdiag:", SRs[0, :, :].diagonal())
    print("Sdiag:", SRs[0, :, :])
    print(np.sum(SRs, axis=0).diagonal())

    print(ham.Sk([0,0,0]).diagonal())

    h = ham.ham
    for icell, cell in enumerate(ham.Rlist):
        for i in range(nrow):
            if icell == 0:
                print(h[i, i+icell*nrow, 0], end=' ')

if __name__ == '__main__':
    test_SiestaHam()


