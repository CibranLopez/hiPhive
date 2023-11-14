"""
Basic example showing how to carry out an analysis of the ClusterSpace

Runs in approximately 4 minutes on an Intel Core i5-4670K CPU.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from ase.build import bulk
from hiphive import ClusterSpace

tol = 1e-4

# setup
max_cutoffs = [10, 10, 7]
atoms = bulk('Al', 'fcc', a=4.05)
cs = ClusterSpace(atoms, max_cutoffs, acoustic_sum_rules=False)

# get data
df = pd.DataFrame(cs.orbit_data)

# collect number of orbits, clusters and parameters
parameter_data = {}
for order in cs.cutoffs.orders:
    parameter_data[order] = []
    cutoffs = sorted(set(df.loc[df.order == order, 'maximum_distance']))
    for cutoff in cutoffs:
        selection = 'maximum_distance < {} and order == {}'.format(
            cutoff+tol, order)
        sub_orbits = df.query(selection)
        num_orbits = len(sub_orbits)
        num_clusters = sub_orbits.n_clusters.sum()
        num_params = sub_orbits.n_parameters.sum()
        parameter_data[order].append(
            [cutoff, num_orbits, num_clusters, num_params])
    parameter_data[order] = np.array(parameter_data[order])

# plot
fs = 14
lw = 2
ms = 9
fig = plt.figure(figsize=(14, 4))
ax1 = fig.add_subplot(1, 3, 1)
ax2 = fig.add_subplot(1, 3, 2)
ax3 = fig.add_subplot(1, 3, 3)
for i, ax in enumerate([ax1, ax2, ax3], start=1):
    for order, data in parameter_data.items():
        ax.semilogy(data[:, 0], data[:, i], '-o', ms=ms, lw=lw, label=order)
    ax.legend(loc=2,  title='Order', fontsize=fs)
    ax.set_xlim([0.0, np.max(max_cutoffs)])
    ax.set_ylim(bottom=0.9)
    ax.set_xlabel('Cutoff (A)', fontsize=fs)
    ax.tick_params(labelsize=fs)

ax1.set_ylabel('Number of orbits', fontsize=fs)
ax2.set_ylabel('Number of clusters', fontsize=fs)
ax3.set_ylabel('Number of parameters', fontsize=fs)


plt.tight_layout()
plt.savefig('clusters.svg')
