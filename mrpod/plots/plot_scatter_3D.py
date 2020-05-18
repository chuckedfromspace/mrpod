"""
Plot 3D scatters points prior and post POD.
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from matplotlib import rcParams
from mrpod import pod_modes
from mrpod.examples.scatter_points import scatter_3D

params = {
        'legend.fontsize': 18,
        'axes.labelsize': 18,
        'xtick.labelsize': 16,
        'ytick.labelsize': 16,
        'xtick.top': False,
        'ytick.right': False,
        'font.family': 'serif',
        'font.serif': "Times New Roman",  # or "Times"
        'mathtext.fontset': 'cm',  # 'stixsans', 'cm'
        'mathtext.rm': 'serif',  # 'sansserif', 'serif'
        'grid.color': 'k',
        'grid.linestyle': ':',
        'grid.linewidth': 0.5
        }
rcParams.update(params)

X = scatter_3D(250, R1=3.5, R2=1.5, angle=60)
corr_mat = X.T @ X
pod_results = pod_modes(X.T, num_of_modes=3)
proj_coeffs = pod_results['proj_coeffs']
eigenvec = pod_results['modes']
print(eigenvec)
fig = plt.figure(1, figsize=(8, 7))
ax = fig.gca(projection='3d')
coord_lim = 4

# ax.scatter(proj_coeffs[0,:], -proj_coeffs[1,:], c='b')
ax.scatter(X[0,:], X[1,:], X[2,:], c='b')
ax.set_xlim([-coord_lim, coord_lim])
ax.set_ylim([-coord_lim, coord_lim])
# ax.set_zlim([0, 1.2])
# ax.set_aspect(1)

# plot the eigenvectors

# ax.axis('off')

plt.show()
