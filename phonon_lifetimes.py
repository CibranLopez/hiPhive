"""
Calculate third order properties using a ForceConstantPotential and feeding
the resulting force constants to phono3py.

Runs in approximately 180 seconds on an Intel Core i5-4670K CPU.
"""
import h5py
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from ase import Atoms
from ase.io import write
from ase.build import bulk
from hiphive import ForceConstantPotential
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms

# parameters
dim = 5
mesh = 14
temperatures = [1200]

# get phono3py supercell
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

# write force constants and POSCAR
fcs.write_to_phonopy('fc2.hdf5')
fcs.write_to_phono3py('fc3.hdf5')
write('POSCAR', prim)

# call phono3py
phono3py_cmd = 'phono3py --dim="{0} {0} {0}" --fc2 --fc3 --br --mesh="'\
               '{1} {1} {1}" --ts="{2}"'.format(
                dim, mesh, ' '.join(str(T) for T in temperatures))
subprocess.call(phono3py_cmd, shell=True)

# collect phono3py data
with h5py.File('kappa-m{0}{0}{0}.hdf5'.format(mesh), 'r') as hf:
    temperatures = hf['temperature'][:]
    frequency = hf['frequency'][:]
    gamma = hf['gamma'][:]

# generate plot of lifetimes
ms = 4
fs = 14
plt.figure()
plt.plot(frequency.flatten(), gamma[0].flatten(), 'o', ms=ms)

plt.xlabel('Frequency (THz)', fontsize=fs)
plt.ylabel('Gamma (THz)', fontsize=fs)
plt.xlim(left=0)
plt.ylim(bottom=0)

plt.title('T={:d}K'.format(int(temperatures[0])))
plt.gca().tick_params(labelsize=fs)
plt.tight_layout()

plt.savefig('frequency_gamma.pdf')
