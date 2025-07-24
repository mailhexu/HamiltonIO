import numpy as np
from ase.build import molecule
from gpaw import GPAW, LCAO

# 1. Set up the system and calculator
# Create a water molecule using ASE
atoms = molecule('H2O')
# Set the cell size with sufficient vacuum
atoms.set_cell([6.0, 6.0, 6.0])
atoms.center()

# Initialize the GPAW calculator in LCAO mode
# Use the Double-zeta polarized (dzp) basis set
calc = GPAW(mode=LCAO(),
            xc='PBE',
            txt='h2o_lcao.txt')

atoms.calc = calc

# 2. Perform the self-consistent field calculation
energy = atoms.get_potential_energy()

# 3. Access the Kohn-Sham Hamiltonian and Overlap matrix
# The hamiltonian and overlap are stored per k-point and spin
# For a non-spin-polarized calculation without k-points, we have one set of matrices
h_skmm, s_kmm = calc.hamiltonian.get_hamiltonian_and_overlap()

# For a simple non-periodic system, we are interested in the Gamma point (k=0)
# and for a non-spin-polarized calculation, spin=0.
h_mm = h_skmm[0, 0]
s_mm = s_kmm[0]

# 4. Save the matrices to files
np.savetxt('hamiltonian.dat', h_mm)
np.savetxt('overlap.dat', s_mm)

print("Kohn-Sham Hamiltonian and overlap matrix have been saved to 'hamiltonian.dat' and 'overlap.dat'")
print(f"Dimensions of the Hamiltonian matrix: {h_mm.shape}")
print(f"Dimensions of the overlap matrix: {s_mm.shape}")

# To load the matrices back into numpy arrays:
# h_loaded = np.loadtxt('hamiltonian.dat')
# s_loaded = np.loadtxt('overlap.dat')
