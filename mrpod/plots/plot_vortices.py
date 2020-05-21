"""
Plot 2D vortices and create movies.
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import rcParams
from mrpod.examples.vortex_shedding import precessing_vortex_core, toroidal_vortex

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

# set up the canvas
x = np.arange(-2.9, 3.1, 0.05)
y = np.arange(0, 3, 0.05)
X, Y = np.meshgrid(x, y)
fig, ax = plt.subplots(1, 1, figsize=(6, 6))
ax.set_aspect('1')

def single_shot(func, vskip=3):
    """
    func is either precessing_vortex_core or toroidal_vortex
    """
    v_array = func
    v_x = v_array[0, :]
    v_y = v_array[1, :]
    ax.quiver(X[::vskip, ::vskip], Y[::vskip, ::vskip], v_x[::vskip, ::vskip],
              v_y[::vskip, ::vskip], units='width', scale=20)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xticks([-2, 0, 2])
    ax.set_yticks([0, 2])
    ax.set_ylabel(r'$y$')
    ax.set_xlabel(r'$x$')
    ax.set_xlim([-3.12, 3.12])
    ax.set_ylim([-0.2, 3])

    plt.savefig("fig_pvc_01.png", bbox_inches='tight', dpi=150, transparent=False)

def create_movie(func1, func2, vskip=3, f=500):
    t = np.arange(8)*1/(8*f)
    ims = []
    for _t in t:
        v_array = func1*np.cos(2*np.pi*f*_t) + func2*np.sin(2*np.pi*f*_t)
        v_x = v_array[0, :]
        v_y = v_array[1, :]
        Q = ax.quiver(X[::vskip, ::vskip], Y[::vskip, ::vskip], v_x[::vskip, ::vskip],
                      v_y[::vskip, ::vskip], units='width', scale=20)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_xticks([-2, 0, 2])
        ax.set_yticks([0, 2])
        ax.set_ylabel(r'$y$')
        ax.set_xlabel(r'$x$')
        ax.set_xlim([-3.12, 3.12])
        ax.set_ylim([-0.2, 3])
        ims.append([Q])

    ani = animation.ArtistAnimation(fig, ims, interval=150/f*500, blit=True)
    ani.save('mov_tv.gif', dpi=300, writer='imagemagick')

# plot a single frame
# single_shot(precessing_vortex_core(offset=-10))
# single_shot(toroidal_vortex(offset=0))

# create a movie
# create_movie(precessing_vortex_core(offset=-10), precessing_vortex_core(offset=0))
# create_movie(toroidal_vortex(offset=-10), toroidal_vortex(offset=0), f=1000)


plt.show()