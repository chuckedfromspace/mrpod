"""
Plot 3D scatters points prior and post POD.
"""
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from matplotlib import rcParams
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d
from mrpod import pod_modes
from mrpod.examples.scatter_points import scatter_3D

class Arrow3D(FancyArrowPatch):

    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0, 0), (0, 0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
        FancyArrowPatch.draw(self, renderer)

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
# eigvals = pod_results['eigvals']
# eigvals = eigvals/eigvals.sum()*100

# print(eigvals)

fig = plt.figure(1, figsize=plt.figaspect(0.6))
ax = fig.gca(projection='3d')
coord_lim = 4
x_axis = Arrow3D([-coord_lim, coord_lim], [0, 0], [0, 0],
                 mutation_scale=20, lw=1, arrowstyle="-|>", color="k")
y_axis = Arrow3D([0, 0], [-coord_lim, coord_lim], [0, 0],
                 mutation_scale=20, lw=1, arrowstyle="-|>", color="k")
z_axis = Arrow3D([0, 0], [0, 0], [0, 2],
                 mutation_scale=20, lw=1, arrowstyle="-|>", color="k")

ax.text(coord_lim, -0.3, 0, r'$x$', fontsize=16)
ax.text(0.2, coord_lim, 0, r'$y$', fontsize=16)
ax.text(0, 0, 2, r'$z$', fontsize=16)

ax.scatter(X[0,:], X[1,:], X[2,:], c='b')
ax.add_artist(x_axis)
ax.add_artist(y_axis)
ax.add_artist(z_axis)
ax.set_xlim([-coord_lim, coord_lim])
ax.set_ylim([-coord_lim, coord_lim])
ax.set_zlim([0, 2])
ax.set_xticklabels([])
ax.set_yticklabels([])
ax.set_zticklabels([])

# plot the eigenvectors
eigenvec1 = Arrow3D([0, eigenvec[0][0]], [0, eigenvec[0][1]], [0, eigenvec[0][2]],
                    mutation_scale=20, lw=1, arrowstyle="-|>", color="r")
eigenvec2 = Arrow3D([0, -eigenvec[1][0]], [0, -eigenvec[1][1]], [0, -eigenvec[1][2]],
                    mutation_scale=20, lw=1, arrowstyle="-|>", color="r")
eigenvec3 = Arrow3D([0, eigenvec[2][0]], [0, eigenvec[2][1]], [0, eigenvec[2][2]],
                    mutation_scale=20, lw=1, arrowstyle="-|>", color="r")
ax.add_artist(eigenvec1)
ax.add_artist(eigenvec2)
ax.add_artist(eigenvec3)

ax.text(*eigenvec[0], r'$\overrightarrow{v}_1$', color='r', fontsize=16)
ax.text(*-eigenvec[1], r'$\overrightarrow{v}_2$', color='r', fontsize=16)
ax.text(*eigenvec[2], r'$\overrightarrow{v}_3$', color='r', fontsize=16)
ax.view_init(30, -64)

# # projections
# _, ax = plt.subplots(1, figsize=(5, 5))
# ax.scatter(-proj_coeffs[1,:], proj_coeffs[2,:], c='b')
# ax.hlines(0, -coord_lim, coord_lim, 'k')
# ax.annotate('', xy=(coord_lim, 0), xytext=(coord_lim-0.2, 0),
#             arrowprops=dict(color='k', width=0.5, headlength=8, headwidth=5),
#             color='k', ha='center')
# ax.vlines(0, -coord_lim, coord_lim, 'k')
# ax.annotate('', xy=(0,coord_lim), xytext=(0,coord_lim-0.2),
#             arrowprops=dict(color='k', width=0.5, headlength=8, headwidth=5),
#             color='k', ha='center')
# ax.text(coord_lim, -0.3, r'$\overrightarrow{v}_1$', fontsize=16)
# ax.text(0.2, coord_lim, r'$\overrightarrow{v}_3$', fontsize=16)
# ax.axis('off')

# plt.savefig('fig_scatter_3D.png', bbox_inches='tight', dpi=150, transparent=False)
plt.show()
