from pathlib import Path
import numpy as np
from HamiltonIO import SiestaHam

def test_SiestaHam_collinear():
    fdf_file = Path('bccFe/collinear/siesta.fdf')
    ham=SiestaHam(fdf_file, spin=0)
    print(f"R   = {ham.Rlist}")
    print(ham.ham.spin)
    HRs = ham.get_HR_all()
    SRs = ham.get_SR_all()
    print(HRs[0, 0, :, :].diagonal())
    print("Sdiag:", SRs[0, :, :].diagonal())
    print(np.sum(SRs, axis=0).diagonal())

    print(ham.Sk([0,0,0]).diagonal())

def test_SiestaHam_nonpolarized():
    fdf_file = Path('bccFe/nonpolarized/siesta.fdf')
    ham=SiestaHam(fdf_file, spin=0)
    print(f"R   = {ham.Rlist}")
    print(ham.ham.spin)
    HRs = ham.get_HR_all()
    SRs = ham.get_SR_all()
    print(HRs[0, :, :].diagonal())
    print("Sdiag:", SRs[0, :, :].diagonal())
    print(np.sum(SRs, axis=0).diagonal())
    print(ham.Sk([0,0,0]).diagonal())


def test_SiestaHam_SOC():
    fdf_file = Path('bccFe/SOC/siesta.fdf')
    ham=SiestaHam(fdf_file)
    print(f"R   = {ham.Rlist}")
    print(ham.ham.spin)
    HRs = ham.get_HR_all()
    SRs = ham.get_SR_all()
    print(HRs[0, :, :].diagonal())
    print("Sdiag:", SRs[0, :, :].diagonal())
    print(np.sum(SRs, axis=0).diagonal())
    print(ham.Sk([0,0,0]).diagonal())

if __name__ == '__main__':
    #test_SiestaHam_collinear()
    #test_SiestaHam_nonpolarized()
    test_SiestaHam_SOC()


