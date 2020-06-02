"""
An example to generate a simple dataset of time-evolving vortex shedding
"""
import numpy as np
# from .damping_wave import signal_synth
# from ..utils import pkl_dump

def phantomGaussian(size, fwhm=3, center=None, scale=1, offset_x=0):
    """
    Create a Gaussian shape inside a square (soft edge shapes)

    Parameters
    ----------
    size: int
        Number of pixels of the side of a square.
    fwhm: int
        Number of pixels of the FWHM of the Gaussian shape.
    center: array_like
        Center of the shape.
    scale: 1, float, optional
        To rescale the size of the phantom shape.
    offset_x: 0, int, optional
        Additional offset to the x axis.

    Returns
    -------
    Output: 2D array 
        A Gaussian in a square.
    """
    x = np.arange(0, size, 1, float)
    R = size // 2
    y = x[:R,np.newaxis]

    if center is None:
        x0 = y0 = R
    else:
        x0 = center[0]
        y0 = center[1]
    output = np.exp(-4*np.log(2) * ((x-x0/scale-offset_x)**2 + (y-y0/scale)**2) / (fwhm/scale)**2)

    return output

def pvc_1(size=120, R=60, fac_x=0.951662, fac_y=0.9222):
    """Create a precessing vortex core with spatial shift
    """
    v_x_settings = [[-0.7, 35, 0, 0],
                    [-0.4, 20, 30, -2], [-0.4, 20, -30, -2],
                    [1.1, 32, 26, 28], [1.1, 32, -26, 28],
                    [0.4, 18, 48, 28], [0.4, 18, -48, 28],
                    [-0.6, 30, 50, 54], [-0.8, 30, 30, 48], [-0.6, 30, -50, 54], [-0.8, 30, -30, 48]]
    v_x_settings = np.array(v_x_settings)
    v_x_settings[:,2] = v_x_settings[:,2]*fac_x
    v_x_settings[:,3] = v_x_settings[:,3]*fac_x
    v_x = np.zeros([R, size])
    for i in v_x_settings:
        v_x += i[0]*phantomGaussian(size, fwhm=i[1], center=[i[2]+R, i[3]])

    v_y_settings = [[1., 30, 6, 0], [-1., 30, -6, 0], [0.5, 15, 3, 0], [-0.5, 15, -3, 0],
                    [-0.3, 10, 30, 2], [0.3, 10, -30, 2], [-0.5, 10, 30, 12], [0.5, 10, -30, 12],
                    [-1, 30, 25, 40], [1, 30, -25, 40],
                    [0.8, 35, 40, 70], [-0.8, 35, -40, 70],
                    [0.7, 10, 46, 35], [-0.7, 10, -46, 35],
                    [0.3, 10, 43, 20], [-0.3, 10, -43, 20],
                    [-0.1, 30, 60, 10], [0.1, 30, -60, 10]]
    v_y_settings = np.array(v_y_settings)
    v_y_settings[:,2] = v_y_settings[:,2]*fac_y
    v_y_settings[:,3] = v_y_settings[:,3]*fac_y
    v_y = np.zeros([R, size])
    for i in v_y_settings:
        v_y += i[0]*phantomGaussian(size, fwhm=i[1], center=[i[2]+R, i[3]])
    return np.array([v_x*1.078, v_y])

def tv_1(size=120, R=60, fac_x=0.9402, fac_y=0.755):
    """Create a toroidal vortex with spatial shift
    """
    v_x_settings = [[-0.7, 35, 0, 0],
                    [-0.4, 20, 30, -2], [-0.4, 20, -30, -2],
                    [1.1, 32, 26, 28], [1.1, 32, -26, 28],
                    [0.4, 18, 48, 28], [0.4, 18, -48, 28],
                    [-0.6, 30, 50, 54], [-0.8, 30, 30, 48], [-0.6, 30, -50, 54], [-0.8, 30, -30, 48]]
    v_x_settings = np.array(v_x_settings)
    v_x_settings[:,2] = v_x_settings[:,2]*fac_x
    v_x_settings[:,3] = v_x_settings[:,3]*fac_x
    v_x = np.zeros([R, size])
    for i in v_x_settings:
        v_x += i[0]*phantomGaussian(size, fwhm=i[1], center=[i[2], i[3]], scale=2, offset_x=R)

    v_y_settings = [[1., 30, 6, 0], [-1., 30, -6, 0], [0.5, 15, 3, 0], [-0.5, 15, -3, 0],
                    [-0.3, 10, 30, 2], [0.3, 10, -30, 2], [-0.5, 10, 30, 12], [0.5, 10, -30, 12],
                    [-1, 30, 25, 40], [1, 30, -25, 40],
                    [0.8, 35, 40, 70], [-0.8, 35, -40, 70],
                    [0.7, 10, 46, 35], [-0.7, 10, -46, 35],
                    [0.3, 10, 43, 20], [-0.3, 10, -43, 20],
                    [-0.1, 30, 60, 10], [0.1, 30, -60, 10]]
    v_y_settings = np.array(v_y_settings)
    v_y_settings[:,2] = v_y_settings[:,2]*fac_y
    v_y_settings[:,3] = v_y_settings[:,3]*fac_y
    v_y = np.zeros([R, size])
    for i in v_y_settings:
        v_y += i[0]*phantomGaussian(size, fwhm=i[1], center=[i[2], i[3]], scale=2, offset_x=R)
    return np.array([v_x*1.5766, v_y*1.489])*0.5

def pvc_2(size=120, R=60, fac_x=1, fac_y=1):
    """Create a precessing vortex core
    """
    v_x_settings = [[0.8, 15, 0, -7],
                    [0.4, 20, 38, -8], [0.4, 20, -38, -8],
                    [-0.7, 35, 17, 12], [-0.7, 35, -17, 12],
                    [-0.6, 20, 36, 12], [-0.6, 20, -36, 12],
                    [1.1, 32, 30, 38], [1.1, 32, -30, 38],
                    [-1, 30, 52, 70], [-1, 30, 34, 70], [-1, 30, -52, 70], [-1, 30, -34, 70]]
    v_x_settings = np.array(v_x_settings)
    v_x_settings[:,2] = v_x_settings[:,2]*fac_x
    v_x_settings[:,3] = v_x_settings[:,3]*fac_x
    v_x = np.zeros([R, size])
    for i in v_x_settings:
        v_x += i[0]*phantomGaussian(size, fwhm=i[1], center=[i[2]+R, i[3]])

    v_y_settings = [[1, 30, 18, 22], [-1, 30, -18, 22],
                    [-0.3, 10, 5, 0], [0.3, 10, -5, 0],
                    [-0.8, 15, 36, 18], [0.8, 15, -36, 18], [-0.6, 15, 32, 2], [0.6, 15, -32, 2],
                    [0.1, 25, 40, -2], [-0.1, 25, -40, -2],
                    [-1, 30, 30, 50], [1, 30, -30, 50],
                    [1, 35, 53, 64], [-1, 35, -53, 64],
                    [-0.3, 30, 75, 20], [0.3, 30, -75, 20]]
    v_y_settings = np.array(v_y_settings)
    v_y_settings[:,2] = v_y_settings[:,2]*fac_y
    v_y_settings[:,3] = v_y_settings[:,3]*fac_y
    v_y = np.zeros([R, size])
    for i in v_y_settings:
        v_y += i[0]*phantomGaussian(size, fwhm=i[1], center=[i[2]+R, i[3]])
    return np.array([v_x, v_y])

def tv_2(size=120, R=60, fac_x=1, fac_y=1):
    """Create a toroidal vortex
    """
    v_x_settings = [[0.8, 15, 0, -7],
                    [0.4, 20, 38, -8], [0.4, 20, -38, -8],
                    [-0.7, 35, 17, 12], [-0.7, 35, -17, 12],
                    [-0.6, 20, 36, 12], [-0.6, 20, -36, 12],
                    [1.1, 32, 30, 38], [1.1, 32, -30, 38],
                    [-1, 30, 52, 70], [-1, 30, 34, 70], [-1, 30, -52, 70], [-1, 30, -34, 70]]
    v_x_settings = np.array(v_x_settings)
    v_x_settings[:,2] = v_x_settings[:,2]*fac_x
    v_x_settings[:,3] = v_x_settings[:,3]*fac_x
    v_x = np.zeros([R, size])
    for i in v_x_settings:
        v_x += i[0]*phantomGaussian(size, fwhm=i[1], center=[i[2], i[3]], scale=2, offset_x=R)

    v_y_settings = [[1, 30, 18, 22], [-1, 30, -18, 22],
                    [-0.3, 10, 5, 0], [0.3, 10, -5, 0],
                    [-0.8, 15, 36, 18], [0.8, 15, -36, 18], [-0.6, 15, 32, 2], [0.6, 15, -32, 2],
                    [0.1, 25, 40, -2], [-0.1, 25, -40, -2],
                    [-1, 30, 30, 50], [1, 30, -30, 50],
                    [1, 35, 53, 64], [-1, 35, -53, 64]]
    v_y_settings = np.array(v_y_settings)
    v_y_settings[:,2] = v_y_settings[:,2]*fac_y
    v_y_settings[:,3] = v_y_settings[:,3]*fac_y
    v_y = np.zeros([R, size])
    for i in v_y_settings:
        v_y += i[0]*phantomGaussian(size, fwhm=i[1], center=[i[2], i[3]], scale=2, offset_x=R)
    return np.array([v_x, v_y])*0.5

def dataset_stationary(func1, func2, N=400, f=470, add_noise=False, noise_level=0.4):
    """
    Create a dataset using either PVC or TV, with a stationary wave-like motion
    """
    np.random.seed(0)
    t = np.arange(N)*1e-4
    v_array = []
    noise = np.zeros([2, 60, 120])
    for _t in t:
        if add_noise:
            noise = np.random.randn(2, 60, 120)*noise_level
        v_array.append(func1*np.cos(2*np.pi*f*_t) + func2*np.sin(2*np.pi*f*_t)+noise)

    return np.array(v_array)


def dataset_stat_mix():
    """
    Create a dataset using both PVC and TV, with a stationary wave-like motion
    """
    v_array = dataset_stationary(pvc_1(), pvc_2()) + 2*dataset_stationary(tv_1(), tv_2(), f=470*2)

    return v_array

# def dataset_custom(path_write, N=8000, add_noise=True):
#     """
#     Create a dataset containing both PVC and TV, both of which last only a small duration.
#     """
#     np.random.seed(0)
#     s_pvc1 = signal_synth(N, 470, A=1, delay=-1000)
#     s_pvc2 = signal_synth(N, 470, A=1, delay=-1000, phase=np.pi/2)
#     s_tv1 = signal_synth(N, 470*2, A=2, delay=-500, num_of_cycles=200)
#     s_tv2 = signal_synth(N, 470*2, A=2, delay=-500, num_of_cycles=200, phase=np.pi/2)
#     v_array = []
#     noise = np.zeros([2, 60, 120])
#     v_pvc1 = pvc_1()
#     v_pvc2 = pvc_2()
#     v_tv1 = tv_1()
#     v_tv2 = tv_2()
#     for _t in range(N):
#         if add_noise:
#             noise = np.random.randn(2, 60, 120)
#         _v = v_pvc1*s_pvc1[_t] + v_pvc2*s_pvc2[_t] + v_tv1*s_tv1[_t] + v_tv2*s_tv2[_t]
#         v_array.append(_v+noise)
#         print(_t)


#     pkl_dump(path_write, np.array(v_array))