import numpy as np
from .wavelet_transform import WaveletTransform
from .utils import pkl_dump

def ortho_check(v1, v2):
    """
    Check the orthogonality of two vectors.
    """
    v1 = np.reshape(v1, np.prod(np.shape(v1)), order='F')
    v2 = np.reshape(v2, np.prod(np.shape(v2)), order='F')

    return v1.T @ v2


def pod_eigendecomp(corr_mat, tol=1e-14):
    """
    Solving the eigenvalue problem of AX=Lambda.

    Parameters
    ----------
    corr_mat : ndarray
        Cross-correlation matrix with a dimension of NxN.

    Returns
    -------
    eigvals : 1d array
        Eigenvalues from the decomposition.
    proj_coeffs: ndarray
        Projection coefficients (temporal modes) of the shape of NxN.
    """
    Lambda, Psi = np.linalg.eigh(corr_mat)
    eigvals = Lambda[::-1]
    eigvals = eigvals[eigvals > tol]

    proj_coeffs = np.fliplr(Psi)[:, :len(eigvals)] @ np.diag(eigvals**(0.5))

    return eigvals, proj_coeffs.T

def pod_modes(data_array, pod_fcn=None, eigvals=None, proj_coeffs=None,
              num_of_modes=None, normalize_mode=True):
    """
    Calculate the POD modes based on the eigenvalue decomposition.

    Parameters
    ----------
    data_array : ndarray
        Data array arranged in a dimension of NxM0xM1x..., such as a time series (N) of multi-
        dimensional scalar/vector fields.
    pod_fcn : function object
        Function to solve the eigenvalue problem by taking in a pre-computed cross-correlation
        matrix. The output must conform to [eigvals, proj_coeffs]. If None, the built-in solver
        ``pod_eigendecomp`` is used.
    eigvals : 1d array
        Eigenvalues from the decomposition.
    proj_coeffs: ndarray
        Projection coefficients (temporal modes) of the shape of NxN.
    num_of_modes : None, int, optional
        Number of modes to be computed (from the most energetic Mode 1). By default, all the valid
        modes are computed.
    normalize_mode : True, bool, optional
        If True, the magnitudes of the modes will be normalized by their corresponding eingenvalues
        and the sample size.

    Returns
    -------
    dict
        A dictionary containing the eigenvalues, projection coefficients, POD modes, and the cross-
        correlation matrix.
    """
    # use the default pod function if none is supplied
    if pod_fcn is None:
        pod_fcn = pod_eigendecomp

    shape = np.shape(data_array)
    data_array = np.reshape(data_array, [shape[0], np.prod(shape[1:])], order='F')
    corr_mat = None

    if eigvals is None:
        # compute the cross-correlation matrix
        corr_mat = data_array @ data_array.T
        # solve the eigenvalue problem of corr_mat
        eigvals, proj_coeffs = pod_fcn(corr_mat)

    # make sure the desired number of modes does not exceed the valid modes
    if num_of_modes is None:
        num_of_modes = len(eigvals)
    else:
        num_of_modes = min(num_of_modes, len(eigvals))

    Omega = np.diag(eigvals[:num_of_modes]**0.5)
    modes = Omega @ proj_coeffs[:num_of_modes, :] @ data_array
    if normalize_mode:
        modes = np.diag(eigvals[:num_of_modes]**(-1)) @ modes
        modes = modes / (shape[0])**0.5*2
    modes = np.reshape(modes, [num_of_modes, *shape[1:]], order='F')

    dict_data = {'modes': modes,
                 'eigvals': eigvals,
                 'proj_coeffs': proj_coeffs,
                 'corr_mat': corr_mat}
    return dict_data



def mrpod_eigendecomp(corr_mat, js, scales, pod_fcn=None, reflect=False, **kwargs):
    """
    Apply eigenvalue decomposition to cross-correlation matrices reconstructed in specific
    bandpasses using ``WaveletTranform``.

    Parameters
    ----------
    corr_mat : ndarray
        Cross-correlation matrix with a dimension of NxN
    js : 1d arrray of int
        The correponding decomposition levels. Different j levels are admissible for "max-overlap".
    scales : 1d array of int
        All the scales included in the reconstruction. Discrete scales are admissible.
    pod_fcn : function object
        Function to solve the eigenvalue problem by taking in a pre-computed cross-correlation
        matrix. The output must conform to [eigvals, proj_coeffs]. If None, the built-in solver
        ``pod_eigendecomp`` is used.
    reflect : False, bool, optional
        If true, the corr_mat is padded symmetrically along both axes and is truncated after
        reconstruction and before eigenvalue decomposition. Such measure helps reducing the effect
        of uneven boundaries in the wavelet transform process.
    reflected : False, bool, optional
        If true, the supplied corr_mat has already been padded and is truncated before fed into the
        eigenvalue solver.

    Other Parameters
    ----------------
    **kwargs :
        Keyword arguments from ``WaveletTransform`` necessary to perform a wavelet transform.

    Returns
    -------
    eigvals : 1d array
        Eigenvalues from the decomposition
    proj_coeffs : ndarray
        Projection coefficients (temporal modes) of the shape of NxN
    K : ndarray
        Reconstructed cross-correlation matrix within the designed bandpasses.
    """
    # use the default pod function if none is supplied
    if pod_fcn is None:
        pod_fcn = pod_eigendecomp

    corr_mat = np.array(corr_mat)
    N = np.shape(corr_mat)[0]
    if reflect:
        corr_mat = np.pad(corr_mat, ((0, N)), 'symmetric')

    corr_mat_mra = WaveletTransform(X=corr_mat, **kwargs)
    K, T_oper = corr_mat_mra.detail_bundle_2d(js, scales)

    if reflect:
        K = K[:N, :N]

    eigvals, proj_coeffs = pod_fcn(K)

    return eigvals, proj_coeffs, K, T_oper

def mrpod_detail_bundle(data_array, *args, num_of_modes=50, seg=10, subtract_avg=False,
                        reflect=False, normalize_mode=True, full_path_write=None, **kwargs):
    """
    Computes MRPOD modes, eigvals and proj coeffs from a dataset arranged in an ndarray of the shape
    of NxM0xM1x..., with N being the dimension that will be wavelet transformed.

    Parameters
    ----------
    data_array : ndarray
        data array arranged in a dimension of NxM0xM1x..., such as a time series (N) of multi-
        dimensional scalar/vector fields.
    num_of_modes : 50, int, optional
        Number of modes to be computed (from the most energetic Mode 1). By default, all the valid
        modes are computed.
    seg : 10, int, optional
        In case of insufficient computer RAM, the ``data_array`` can be segmented along the M
        dimension of the reshaped NxM and pieced together after the filtering.
    subtract_avg : False, bool, optional
        If True, the ensemble average of the ``data_array`` along N is subtracted from the dataset.
    reflect : False, bool, optional
        If True, the cross-correlation matrix as well as the dataset is padded symmetrically along
        the sample axis.
    normalize_mode : True, bool, optional
        If True, the magnitudes of the modes will be normalized by their corresponding eingenvalues
        and the sample size.
    full_path_write : None, path obj, optional
        If provided, the output is saved locally in the specified directory.

    Other Parameters
    ----------------
    **kwargs :
        Keyword arguments from ``mrpod_eigendecomp()`` necessary to perform a wavelet transform.

    Returns
    -------
    dict
        A dictionary containing the eigenvalues, projection coefficients, POD modes, and the cross-
        correlation matrix.
    """
    shape = np.shape(data_array)
    if subtract_avg:
        avg = np.mean(data_array, axis=0)
        data_array = data_array - avg

    # reshape data array into an NxM matrix
    data_array = np.reshape(data_array, [shape[0], np.prod(shape[1:])], order='F')
    # compute the cross-correlation matrix
    corr_mat = data_array @ data_array.T
    # solve the eigenvalue problem of corr_mat
    eigvals, proj_coeffs, K, T_oper = mrpod_eigendecomp(corr_mat, *args, reflect=reflect,
                                                        **kwargs)
    # make sure the desired number of modes does not exceed the valid modes
    if num_of_modes is None:
        num_of_modes = len(eigvals)
    else:
        num_of_modes = min(num_of_modes, len(eigvals))

    # compute the modes
    T_oper = np.sum(T_oper, axis=(0, 1))
    Omega = np.diag(eigvals[:num_of_modes]**0.5)
    modes = np.zeros((num_of_modes, np.prod(shape[1:])))

    # divide the M dimension into managable segments
    seg_len = np.prod(shape[1:]) // seg
    for i in range(seg):
        data_seg = data_array[:, i*seg_len:(i+1)*seg_len]
        if reflect:
            data_seg = np.pad(data_seg, ((0, shape[0]), (0, 0)), 'symmetric')
        # compute
        B_seg = T_oper @ data_seg
        B_seg = B_seg[:shape[0], :]  # crop out the symmetric part
        modes_seg = Omega @ proj_coeffs[:num_of_modes, :] @ B_seg
        modes[:, i*seg_len:(i+1)*seg_len] = modes_seg

        print('Processing segmant %d' % i)
    data_array = []  # free up RAM

    if normalize_mode:
        modes = np.diag(eigvals[:num_of_modes]**(-1)) @ modes
        modes = modes / (shape[0])**0.5*2

    # convert NxM modes back into the original NxM0xM1... shape.
    modes = np.reshape(modes, [num_of_modes, *shape[1:]], order='F')

    dict_data = {'modes': modes,
                 'eigvals': eigvals,
                 'proj_coeffs': proj_coeffs[:num_of_modes, :],
                 'corr_mat': K}

    if full_path_write is not None:
        pkl_dump(full_path_write, dict_data)

    return dict_data
