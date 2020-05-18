"""
Plot 2D scatters points prior and post POD.
"""
import matplotlib.pyplot as plt
from matplotlib import rcParams
from mrpod import pod_modes
from mrpod.examples.scatter_points import scatter_circular

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

X = scatter_circular(250, R1=3.5, R2=1.5, angle=60)
corr_mat = X.T @ X
pod_results = pod_modes(X.T, num_of_modes=2)
proj_coeffs = pod_results['proj_coeffs']
eigenvec = pod_results['modes']

_, ax = plt.subplots(1, figsize=(5, 5))
coord_lim = 4

ax.scatter(X[0,:], X[1,:], c='b')
ax.set_xlim([-coord_lim, coord_lim])
ax.set_ylim([-coord_lim, coord_lim])
ax.set_aspect(1)

# plot the eigenvectors
ax.annotate('', xy=(eigenvec[0,0], eigenvec[0,1]), xytext=(0,0),
            arrowprops=dict(color='r', width=0.5, headlength=8, headwidth=5),
            color='r', ha='center')
ax.text(eigenvec[0, 0]+0.2, eigenvec[0, 1]+0.2, r'$v_1$', fontsize=16, color='r',
        ha='center', va='center')
ax.annotate('', xy=(-eigenvec[1,0], -eigenvec[1,1]), xytext=(0,0),
            arrowprops=dict(color='r', width=0.5, headlength=8, headwidth=5),
            color='r', ha='center')
ax.text(-eigenvec[1, 0]-0.2, -eigenvec[1, 1]+0.2, r'$v_2$', fontsize=16, color='r',
        ha='center', va='center')

# add fake axes
ax.hlines(0, -coord_lim, coord_lim, 'k')
ax.annotate('', xy=(coord_lim, 0), xytext=(coord_lim-0.2, 0),
            arrowprops=dict(color='k', width=0.5, headlength=8, headwidth=5),
            color='k', ha='center')
ax.vlines(0, -coord_lim, coord_lim, 'k')
ax.annotate('', xy=(0,coord_lim), xytext=(0,coord_lim-0.2),
            arrowprops=dict(color='k', width=0.5, headlength=8, headwidth=5),
            color='k', ha='center')
ax.text(coord_lim, -0.3, r'$x$', fontsize=16)
ax.text(0.2, coord_lim, r'$y$', fontsize=16)

# # projections
# ax.scatter(proj_coeffs[0,:], -proj_coeffs[1,:], c='b')
# ax.hlines(0, -coord_lim, coord_lim, 'k')
# ax.annotate('', xy=(coord_lim, 0), xytext=(coord_lim-0.2, 0),
#             arrowprops=dict(color='k', width=0.5, headlength=8, headwidth=5),
#             color='k', ha='center')
# ax.vlines(0, -coord_lim, coord_lim, 'k')
# ax.annotate('', xy=(0,coord_lim), xytext=(0,coord_lim-0.2),
#             arrowprops=dict(color='k', width=0.5, headlength=8, headwidth=5),
#             color='k', ha='center')
# ax.text(coord_lim, -0.3, r'$v_1$', fontsize=16)
# ax.text(0.2, coord_lim, r'$v_2$', fontsize=16)

ax.axis('off')

plt.show()
