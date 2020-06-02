"""
Plot 2D vortices and create movies.
"""
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.gridspec as gridspec
from matplotlib import rcParams
from mrpod.examples.vortex_shedding import (dataset_stat_mix, tv_1, tv_2,
                                            dataset_stationary, pvc_1, pvc_2)
from mrpod import pod_modes

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

def single_shot(func, vskip=3):
    """
    func is either precessing_vortex_core or toroidal_vortex
    """
    _, ax = plt.subplots(1, 1, figsize=(6, 6))
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
    ax.set_aspect('1')

    plt.savefig("fig_pvc_01.png", bbox_inches='tight', dpi=150, transparent=False)

def create_movie(func1, func2, vskip=3, f=500, add_noise=False):
    """
    func1 and func2 are spatially-shifted equivalences
    """
    np.random.seed(0)
    fig, ax = plt.subplots(1, 1, figsize=(6, 4))
    ax.set_aspect('1')
    t = np.arange(8)*1/(8*f)
    ims = []
    noise = np.zeros([2, 60, 120])
    for _t in t:
        if add_noise:
            noise = np.random.randn(2, 60, 120)*5
        v_array = func1*np.cos(2*np.pi*f*_t) + func2*np.sin(2*np.pi*f*_t) + noise
        v_x = v_array[0, :]
        v_y = v_array[1, :]
        Q = ax.quiver(X[::vskip, ::vskip], Y[::vskip, ::vskip], v_x[::vskip, ::vskip],
                      v_y[::vskip, ::vskip], units='width', scale=120)
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
    ani.save('mov_pvc_subnoise.gif', dpi=150, writer='imagemagick')

def dataset_pod(func1, func2, plot_modes=True, plot_phase=True, plot_ps=True,
                plot_movie=False, **kwargs):
    """
    Carry out POD on a synthesized dataset
    """
    v_array = dataset_stationary(func1, func2, **kwargs)
    # v_array = dataset_stat_mix()
    pod_results = pod_modes(v_array, num_of_modes=4, normalize_mode=True)
    proj_coeffs = pod_results['proj_coeffs']
    modes = pod_results['modes']
    eigvals = pod_results['eigvals']
    eigvals = eigvals/eigvals.sum()*100
    print(eigvals[:6])

    if plot_movie:
        fig, ax = plt.subplots(1, 1, figsize=(6, 4))
        ax.set_aspect('1')
        vskip = 3
        ims = []
        for i in range(21):
            # v_array = proj_coeffs[0, i]*modes[0, :] + proj_coeffs[1, i]*modes[1, :]
            v_x = v_array[i, 0, :]
            v_y = v_array[i, 1, :]
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

        ani = animation.ArtistAnimation(fig, ims, interval=80, blit=True)
        # ani.save('mov_mix.gif', dpi=150, writer='imagemagick')

    if plot_modes:
        plt.figure(1, figsize=(8, 4))
        gs1 = gridspec.GridSpec(1, 2)
        gs1.update(left=0.1, right=0.9, top=0.9, bottom=0.48, wspace=0.05)
        ax = np.empty([2], dtype='object')
        ax[0] = plt.subplot(gs1[0, 0])
        ax[1] = plt.subplot(gs1[0, 1])
        vskip = 3
        modes[1, :] = -1*modes[1, :]
        for i in [0, 1]:
            v_x = modes[i+0, 0, :]
            v_y = modes[i+0, 1, :]
            ax[i].quiver(X[::vskip, ::vskip], Y[::vskip, ::vskip], v_x[::vskip, ::vskip],
                         v_y[::vskip, ::vskip], units='width', scale=25)
            ax[i].spines['top'].set_visible(False)
            ax[i].spines['right'].set_visible(False)
            ax[i].set_xticks([-2, 0, 2])
            ax[i].set_yticks([0, 2])

            ax[i].set_xlabel(r'$x$')
            ax[i].set_xlim([-3.12, 3.12])
            ax[i].set_ylim([-0.2, 3])

        ax[0].set_ylabel(r'$y$')
        ax[1].set_yticks([])
        plt.savefig("fig_subnoise_modes.png", bbox_inches='tight', dpi=150, transparent=False)

    if plot_phase:
        _, ax = plt.subplots(1, figsize=(3, 3))
        ax.scatter(proj_coeffs[0,::3]/np.sqrt(eigvals[0]), 
                   -proj_coeffs[1,::3]/np.sqrt(eigvals[1]), c='b')
        # ax.set_xlim([-11, 11])
        # ax.set_ylim([-11, 11])
        # ax.set_xticks([])
        # ax.set_yticks([])
        ax.set_aspect('1')
        ax.set_xlabel(r'$a_1$')
        ax.set_ylabel(r'$a_2$')
        plt.savefig("fig_subnoise_phase.png", bbox_inches='tight', dpi=150, transparent=False)

    if plot_ps:
        _, ax = plt.subplots(1, figsize=(4, 2))
        line_colors = ['k', 'k--']
        for i in range(2):
            frq, Amp = signal.welch(proj_coeffs[i+0, :], 10000,
                                    window='boxcar', scaling='density', detrend=False)
            ax.plot(frq, np.log10(Amp), line_colors[i], label=r'$a_{%d}$' % (i+1, ))

        ax.set_xlim(0, 1500)
        ax.set_ylim(-3, 2.5)
        ax.set_xticks([0, 500, 1000, 1500])
        ax.set_yticks([-2, -1, 0, 1, 2,])
        ax.set_yticklabels(['0.01', '0.1', '1', '10', '100'], fontsize=12)
        ax.set_xlabel(r'$f$ [Hz]')
        ax.set_ylabel('PSD')
        ax.set_aspect('2.e2')
        ax.legend(loc='upper right', labelspacing=0.1,
                  handlelength=0.8, handletextpad=0.2, frameon=False)
        plt.savefig("fig_subnoise_ps.png", bbox_inches='tight', dpi=150, transparent=False)


def plot_time_shift(f=470):
    """
    """
    fig, ax = plt.subplots(1, 1, figsize=(12, 3))
    t = np.arange(150)*1e-4
    ax.plot(t, np.cos(2*np.pi*f*t), 'k', label=r'Re$(\overrightarrow{\Phi})$')
    ax.plot(t, np.sin(2*np.pi*f*t), 'k--', label=r'Im$(\overrightarrow{\Phi})$')
    ax.legend(loc=0, handletextpad=0.2, fontsize=14)
    ax.axis('off')
    ax.hlines(0, 0, 1.6e-2, 'k')
    ax.annotate('', xy=(1.6e-2, 0), xytext=(-2e-4, 0),
                arrowprops=dict(color='k', width=0.5, headlength=8, headwidth=5),
                color='k', ha='center')
    ax.vlines(0, -1.2, 1.2, 'k')
    ax.annotate('', xy=(0, 1.2), xytext=(0, 0),
                arrowprops=dict(color='k', width=0.5, headlength=8, headwidth=5),
                color='k', ha='center')
    ax.text(1.62e-2, -0.2, r'$t$', fontsize=16)
    ax.text(0.2e-3, 1.2, r'$Amp$', fontsize=16)
    ax.set_ylim([-1.2, 1.2])
    ax.set_xlim([-2e-4, 1.6e-2])
    plt.savefig("fig_time_shift.png", bbox_inches='tight', dpi=150, transparent=False)

# set up the canvas
x = np.arange(-2.9, 3.1, 0.05)
y = np.arange(0, 3, 0.05)
X, Y = np.meshgrid(x, y)



# plots
# plot_time_shift()

# single_shot(pvc_1())
# single_shot(toroidal_vortex(offset=0))

# create a movie
create_movie(pvc_1(), pvc_2(), f=470, add_noise=True)
# create_movie(tv_1(), tv_2(), f=470*2)

# dataset_pod(pvc_1(), pvc_2(), f=470, add_noise=True, noise_level=5)


plt.show()