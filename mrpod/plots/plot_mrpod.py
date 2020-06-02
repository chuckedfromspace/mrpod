import numpy as np
from scipy import signal
import pywt
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.gridspec as gridspec
from matplotlib import rcParams
from mrpod.examples.vortex_shedding import dataset_stat_mix, pvc_1, pvc_2, dataset_stationary
from mrpod import CompositeFilter, WaveletTransform, mrpod_detail_bundle
from mrpod.utils import pkl_load


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

dir_filterbank = "filters_sym12_j=8.pkl"


def dataset_mrpod(plot_modes=True, plot_phase=True, plot_ps=True):
    """
    Carry out MRPOD on a synthesized dataset
    """
    x = np.arange(-2.9, 3.1, 0.05)
    y = np.arange(0, 3, 0.05)
    X, Y = np.meshgrid(x, y)
    # v_array = dataset_stat_mix()
    v_array = dataset_stationary(pvc_1(), pvc_2(), f=470, add_noise=True, noise_level=0.4)
    n_scale = 0
    j_level = 3
    pod_results = mrpod_detail_bundle(v_array, js=[j_level], scales=[n_scale], seg=1, num_of_modes=10,
                                      full_path_filterbank=dir_filterbank)
    proj_coeffs = pod_results['proj_coeffs']
    modes = pod_results['modes']
    eigvals = pod_results['eigvals']
    eigvals = eigvals/eigvals.sum()*100
    print(eigvals[:6])
    
    if plot_modes:
        plt.figure(1, figsize=(8, 4))
        gs1 = gridspec.GridSpec(1, 2)
        gs1.update(left=0.1, right=0.9, top=0.9, bottom=0.48, wspace=0.05)
        ax = np.empty([2], dtype='object')
        ax[0] = plt.subplot(gs1[0, 0])
        ax[1] = plt.subplot(gs1[0, 1])
        vskip = 3
        # modes[1, :] = -1*modes[1, :]
        for i, _mode_n in enumerate([1, 0]):
            v_x = modes[_mode_n+0, 0, :]
            v_y = modes[_mode_n+0, 1, :]
            ax[i].quiver(X[::vskip, ::vskip], Y[::vskip, ::vskip], v_x[::vskip, ::vskip],
                         v_y[::vskip, ::vskip], units='width', scale=20)
            ax[i].spines['top'].set_visible(False)
            ax[i].spines['right'].set_visible(False)
            ax[i].set_xticks([-2, 0, 2])
            ax[i].set_yticks([0, 2])

            ax[i].set_xlabel(r'$x$')
            ax[i].set_xlim([-3.12, 3.12])
            ax[i].set_ylim([-0.2, 3])

        ax[0].set_ylabel(r'$y$')
        ax[1].set_yticks([])
        plt.savefig("fig_subnoise_mrpodmodes.png", bbox_inches='tight', dpi=150, transparent=False)

    if plot_phase:
        _, ax = plt.subplots(1, figsize=(3, 3))
        proj_coeffs = proj_coeffs[:, 10:-10]
        ax.scatter(proj_coeffs[0, ::2]/np.sqrt(eigvals[0]), 
                   proj_coeffs[1, ::2]/np.sqrt(eigvals[1]), c='b')
        # ax.set_xlim([-11, 11])
        # ax.set_ylim([-11, 11])
        ax.set_xticks([-6, 0, 6])
        ax.set_yticks([-6, 0, 6])
        ax.set_aspect('1')
        ax.set_xlabel(r'$a_1$')
        ax.set_ylabel(r'$a_2$')
        plt.savefig("fig_subnoise_mrpodphase.png", bbox_inches='tight', dpi=150, transparent=False)

    if plot_ps:
        w = pywt.Wavelet('sym8')
        w_filter = CompositeFilter(w.dec_lo)
        gain, f = w_filter.sqd_gain_fcn(j_level, n_scale, max_overlap=True)
        _, ax = plt.subplots(1, figsize=(4, 2))
        line_colors = ['k', 'k--']
        for i in range(2):
            frq, Amp = signal.welch(proj_coeffs[i+0, 5:-5], 10000,
                                    window='boxcar', scaling='density', detrend=False)
            ax.plot(frq, np.log10(Amp), line_colors[i], label=r'$a_{%d}$' % (i+1, ))

        # print()
        ax.plot(f*1e4, np.log10(gain**2), 'r')
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
        plt.savefig("fig_subnoise_mrpodps.png", bbox_inches='tight', dpi=150, transparent=False)

def filter_design():
    filterbank = pkl_load(dir_filterbank)
    _, ax = plt.subplots(1, figsize=(6, 3))
    
    # t = np.arange(400)*1e-4
    # Fourier_filter = np.sin(2*np.pi*500*t)
    # ax.plot(Fourier_filter, 'b')
    wavelet_filter = filterbank['10']
    to_pad = (400-len(wavelet_filter))//2
    wavelet_filter = np.pad(wavelet_filter, (to_pad, 400-to_pad))
    ax.plot(wavelet_filter, 'b')
    ax.hlines(0, 0, 400, 'k')
    ax.axis('off')
    ax.set_ylim([-0.2, 0.8])
    # ax.annotate('', xy=(coord_lim, 0), xytext=(coord_lim-0.2, 0),
    #             arrowprops=dict(color='k', width=0.5, headlength=8, headwidth=5),
    #             color='k', ha='center')



dataset_mrpod()
# filter_design()

plt.show()
