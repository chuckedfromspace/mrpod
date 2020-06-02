"""
Generate a damping wave.
"""
import numpy as np
from scipy import signal

def inverse_fft(f=470, bandwidth=(-15, 15), dN=1000, N=8000, fs=10000, add_noise=False):
    """
    Create a time-varying sinusoidal wave
    """
    df = (fs/2)/(N/2)
    f_c = int(f/df)
    shift = np.zeros(N) + 1e-2
    shift[1000+dN:4000+dN] = signal.windows.general_gaussian(3000, p=5.5, sig=3000/4)
    n = np.zeros((8000,), dtype=complex)
    n[f_c+bandwidth[0]:f_c+bandwidth[1]] = np.exp(1j*np.random.normal(0, 2*np.pi,
                                                  (bandwidth[1]-bandwidth[0],)))
    s = (np.fft.ifft(n)).real*110 * shift

    if add_noise:
        noise = np.random.uniform(-0.06, 0.06, N)
    else:
        noise = 0

    return s + noise

def signal_synth(N, f_p, fs=10000, A=2, num_of_cycles=120, delay=0, p=2.5, s=4, phase=0):
    '''
    N is the total length of the signal
    f_p is the central frequency of the signal
    '''
    cycles = num_of_cycles*np.int(1/f_p*fs)
    t_x = np.arange(cycles)/fs
    x = A*np.cos(2*np.pi*f_p*t_x+phase)*signal.windows.general_gaussian(cycles, p=p, sig=cycles/s)
    del_N = N-cycles
    x = np.pad(x, (del_N//2 + delay, del_N-del_N//2-delay), 'constant',
               constant_values=(0, 0))
    # t = np.append(t_x, np.arange(del_N-del_N//2-delay)/fs + t_x[-1])
    # t = np.append(np.arange(-del_N//2-delay, 0)/fs, t)

    return x
