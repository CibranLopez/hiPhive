"""
Calculate thermal properties within the harmonic approximation using a
hiPhive force constant potential in combination with phonopy.

Runs in approximately 30 seconds on an Intel Core i5-4670K CPU.
"""
import numpy as np
import matplotlib.pyplot as plt

from ase import Atoms
from ase.build import bulk
from hiphive import ForceConstantPotential

from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms

# parameters
dim = 5  # dimension in phonopy calculation
Nq = 51  # number of q-points along each segment of the path through the BZ
mesh = [32, 32, 32]  # q-point mesh for MSD calculation
temperatures = [300, 600, 900, 1200]  # temperatures for evaluating MSD

# set up phonopy
prim = bulk('Ni')
atoms_phonopy = PhonopyAtoms(symbols=prim.get_chemical_symbols(),
                             scaled_positions=prim.get_scaled_positions(),
                             cell=prim.cell)
phonopy = Phonopy(atoms_phonopy, supercell_matrix=dim*np.eye(3),
                  primitive_matrix=None)

# get force constants
fcp = ForceConstantPotential.read('fcc-nickel.fcp')
supercell = phonopy.get_supercell()
supercell = Atoms(cell=supercell.cell, numbers=supercell.numbers, pbc=True,
                  scaled_positions=supercell.get_scaled_positions())
fcs = fcp.get_force_constants(supercell)
print(fcs)
# access specific parts of the force constant matrix
fcs.print_force_constant((0, 1))
fcs.print_force_constant((10, 12))

# get mean square displacements
phonopy.set_force_constants(fcs.get_fc_array(order=2))
phonopy.set_mesh(mesh, is_eigenvectors=True, is_mesh_symmetry=False)
phonopy.set_thermal_displacements(temperatures=temperatures)
_, msds = phonopy.get_thermal_displacements()
msds = np.sum(msds, axis=1)  # sum up the MSD over x,y,z
for temperature, msd in zip(temperatures, msds):
    print('T = {:4d} K    MSD = {:.5f} A**2'.format(temperature, msd))


# FCC q-point paths
def get_band(q_start, q_stop, N):
    """ Return path between q_start and q_stop """
    return np.array([q_start + (q_stop-q_start)*i/(N-1) for i in range(N)])


G2X = get_band(np.array([0, 0, 0]), np.array([0.5, 0.5, 0]), Nq)
X2K2G = get_band(np.array([0.5, 0.5, 1.0]), np.array([0, 0, 0]), Nq)
G2L = get_band(np.array([0, 0, 0]), np.array([0.5, 0.5, 0.5]), Nq)
bands = [G2X, X2K2G, G2L]

# get phonon dispersion
phonopy.set_band_structure(bands)
qvecs, qnorms, freqs, _ = phonopy.get_band_structure()

# plot dispersion
fig = plt.figure()

kpts = [0.0, qnorms[0][-1], qnorms[1][-1], qnorms[2][-1]]
kpts_labels = ['$\\Gamma$', 'X', '$\\Gamma$', 'L']

plt.axvline(x=kpts[1], color='k', linewidth=0.9)
plt.axvline(x=kpts[2], color='k', linewidth=0.9)

for q, freq, in zip(qnorms, freqs):
    plt.plot(q, freq, color='b', linewidth=2.0)

plt.xlabel('Wave vector $\\vec{q}$', fontsize=14.0)
plt.ylabel('Frequency $\\omega$ (THz)', fontsize=14.0)
plt.xticks(kpts, kpts_labels, fontsize=14.0)
plt.xlim([0.0, qnorms[-1][-1]])
plt.ylim([0.0, 12.0])

plt.tight_layout()
plt.savefig('phonon_dispersion.pdf')
