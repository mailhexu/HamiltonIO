#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The abacus wrapper
"""

import os
from pathlib import Path

from TB2J.utils import symbol_number_list

from HamiltonIO.abacus.abacus_api import read_HR_SR
from HamiltonIO.abacus.orbital_api import parse_abacus_orbital
from HamiltonIO.abacus.stru_api import read_abacus
from HamiltonIO.lcao_hamiltonian import LCAOHamiltonian


class AbacusWrapper(LCAOHamiltonian):
    def __init__(
        self,
        HR,
        SR,
        Rlist,
        nbasis=None,
        atoms=None,
        nspin=1,
        basis=None,
        HR_soc=None,
        HR_nosoc=None,
        nel=None,
    ):
        super().__init__(
            HR,
            SR,
            Rlist,
            nbasis=nbasis,
            atoms=atoms,
            nspin=nspin,
            orbs=None,
            HR_soc=HR_soc,
            HR_nosoc=HR_nosoc,
            nel=nel,
        )
        self._name = "ABACUS"


class AbacusParser:
    def __init__(self, spin=None, outpath=None, binary=False):
        self.outpath = outpath
        if spin is None:
            self.spin = self.read_spin()
        else:
            self.spin = spin
        self.binary = binary
        # read the information
        self.read_atoms()
        self.efermi = self.read_efermi()
        self.nel = self.read_nel()
        if self.spin in ["non-polarized", "collinear"]:
            self.basis = self.read_basis()
        elif self.spin == "noncollinear":
            self.basis = self.read_basis(nspin=2)

    def read_spin(self):
        with open(str(Path(self.outpath) / "running_scf.log")) as myfile:
            for line in myfile:
                if line.strip().startswith("nspin"):
                    nspin = int(line.strip().split()[-1])
                    if nspin == 1:
                        return "non-polarized"
                    elif nspin == 2:
                        return "collinear"
                    elif nspin == 4:
                        return "noncollinear"
                    else:
                        raise ValueError("nspin should be either 1 or 4.")

    def read_atoms(self):
        path1 = str(Path(self.outpath) / "../STRU")
        path2 = str(Path(self.outpath) / "../Stru")
        if os.path.exists(path1):
            self.atoms = read_abacus(path1)
        elif os.path.exists(path2):
            self.atoms = read_abacus(path2)
        else:
            raise Exception("The STRU or Stru file cannot be found.")
        return self.atoms

    def read_basis(self, nspin=1):
        fname = str(Path(self.outpath) / "Orbital")
        self.basis = parse_abacus_orbital(fname, nspin=nspin)
        return self.basis

    def read_HSR_collinear(self, binary=None):
        p = Path(self.outpath)
        SR_filename = p / "data-SR-sparse_SPIN0.csr"
        HR_filename = [p / "data-HR-sparse_SPIN0.csr", p / "data-HR-sparse_SPIN1.csr"]
        nbasis, Rlist, HR_up, HR_dn, SR = read_HR_SR(
            nspin=2,
            binary=self.binary,
            HR_fileName=HR_filename,
            SR_fileName=SR_filename,
        )
        return nbasis, Rlist, HR_up, HR_dn, SR

    def Read_HSR_noncollinear(self, binary=None):
        p = Path(self.outpath)
        SR_filename = str(p / "data-SR-sparse_SPIN0.csr")
        HR_filename = str(p / "data-HR-sparse_SPIN0.csr")
        nbasis, Rlist, HR, SR = read_HR_SR(
            nspin=4,
            binary=self.binary,
            HR_fileName=HR_filename,
            SR_fileName=SR_filename,
        )
        return nbasis, Rlist, HR, SR

    def get_models(self):
        if self.spin == "collinear":
            nbasis, Rlist, HR_up, HR_dn, SR = self.read_HSR_collinear()
            model_up = AbacusWrapper(
                HR=HR_up, SR=SR, Rlist=Rlist, nbasis=nbasis, nspin=1, atoms=self.atoms
            )
            model_dn = AbacusWrapper(
                HR=HR_dn, SR=SR, Rlist=Rlist, nbasis=nbasis, nspin=1, atoms=self.atoms
            )
            model_up.efermi = self.efermi
            model_dn.efermi = self.efermi
            model_up.basis, model_dn.basis = self.get_basis()
            model_up.atoms = self.atoms
            model_dn.atoms = self.atoms
            return model_up, model_dn
        elif self.spin == "noncollinear":
            nbasis, Rlist, HR, SR = self.Read_HSR_noncollinear()
            model = AbacusWrapper(
                HR=HR, SR=SR, Rlist=Rlist, nbasis=nbasis, nspin=2, atoms=self.atoms
            )
            model.efermi = self.efermi
            model.basis = self.get_basis()
            model.atoms = self.atoms
            return model

    def read_efermi(self):
        """
        Reading the efermi from the scf log file.
        Search for the line EFERMI = xxxxx eV
        """
        fname = str(Path(self.outpath) / "running_scf.log")
        efermi = None
        with open(fname, "r") as myfile:
            for line in myfile:
                if "EFERMI" in line:
                    efermi = float(line.split()[2])
        if efermi is None:
            raise ValueError(f"EFERMI not found in the {str(fname)}  file.")
        return efermi

    def read_nel(self):
        """
        Reading the number of electrons from the scf log file.
        """
        fname = str(Path(self.outpath) / "running_scf.log")
        nel = None
        with open(fname, "r") as myfile:
            for line in myfile:
                if "number of electrons" in line:
                    nel = float(line.split()[-1])
        if nel is None:
            raise ValueError(f"number of electron not found in the {str(fname)}  file.")
        return nel

    def get_basis(self):
        slist = symbol_number_list(self.atoms)
        if self.spin == "collinear":
            basis_up = []
            basis_dn = []
            for b in self.basis:
                basis_up.append((slist[b.iatom], b.sym, "up"))
                basis_dn.append((slist[b.iatom], b.sym, "down"))
            return basis_up, basis_dn
        elif self.spin == "noncollinear":
            basis = []
            for b in self.basis:
                basis.append((slist[b.iatom], b.sym, "up"))
                basis.append((slist[b.iatom], b.sym, "down"))
            return basis


class AbacusSplitSOCParser:
    """
    Abacus parser with Hamiltonian split to SOC and non-SOC parts
    """

    def __init__(self, outpath_nosoc=None, outpath_soc=None, binary=False):
        self.outpath_nosoc = outpath_nosoc
        self.outpath_soc = outpath_soc
        self.binary = binary
        self.parser_nosoc = AbacusParser(outpath=outpath_nosoc, binary=binary)
        self.parser_soc = AbacusParser(outpath=outpath_soc, binary=binary)
        spin1 = self.parser_nosoc.read_spin()
        spin2 = self.parser_soc.read_spin()
        if spin1 != "noncollinear" or spin2 != "noncollinear":
            raise ValueError("Spin should be noncollinear")

    def parse(self):
        nbasis, Rlist, HR_nosoc, SR = self.parser_nosoc.Read_HSR_noncollinear()
        nbasis2, Rlist2, HR2, SR2 = self.parser_soc.Read_HSR_noncollinear()
        HR_soc = HR2 - HR_nosoc
        from TB2J.pauli import chargepart, spinpart

        for iR, _ in enumerate(Rlist):
            HR_soc[iR] = spinpart(HR_soc[iR])
            HR_nosoc[iR] += chargepart(HR_soc[iR])

        model = AbacusWrapper(
            HR=None,
            SR=SR,
            Rlist=Rlist,
            nbasis=nbasis,
            nspin=2,
            HR_soc=HR_soc,
            HR_nosoc=HR_nosoc,
            nel=self.parser_nosoc.nel,
        )
        model.efermi = self.parser_soc.efermi
        model.basis = self.parser_nosoc.basis
        model.atoms = self.parser_nosoc.atoms
        return model


def test_abacus_wrapper_collinear():
    outpath = "/Users/hexu/projects/TB2J_abacus/abacus-tb2j-master/abacus_example/case_Fe/1_no_soc/OUT.Fe"
    parser = AbacusParser(outpath=outpath, spin=None, binary=False)
    atoms = parser.read_atoms()
    # atoms=parser.read_atoms_out()
    # parser.read_HSR_collinear()
    model_up, model_dn = parser.get_models()
    H, S, E, V = model_up.HSE_k([0, 0, 0])
    # print(H.diagonal().real)
    # print(model_up.get_HR0().diagonal().real)
    print(parser.efermi)
    return model_up, model_dn, parser, atoms, H, S, E, V


# def test_abacus_wrapper_ncl():
#    outpath = "/Users/hexu/projects/TB2J_abacus/abacus-tb2j-master/abacus_example/case_Fe/2_soc/OUT.Fe"
#
#    parser = AbacusParser(outpath=outpath, spin=None, binary=False)
#    #atoms = parser.read_atoms()
#    #model = parser.get_models()
#    #H, S, E, V = model.HSE_k([0, 0, 0])
#    #print(parser.efermi)
#    retrun parser


if __name__ == "__main__":
    # test_abacus_wrapper()
    # test_abacus_wrapper_ncl()
    pass
