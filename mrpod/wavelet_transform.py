import numpy as np
from .utils import pkl_dump, pkl_load, WAVELETS

def find_scale_index(level, x='0', y='1'):
    """
    Gray code order is used here to generate indices for the scales of each decomposition level in a
    wavelet packet transform. Either 'a' and 'd' (approximation and detail) or '0' and '1' are used.
    E.g., for a level=1 decomposition, the possible combination will be 'ad' and 'da'. For level=2,
    it will be 'aa', 'ad', 'da', 'dd'. The indices are also ordered by their corresponding frequency
    bandpasses.

    Parameters
    ----------
    level : int
        Dcomposition level for WPT.
    x : '0', str, optional
        First index in the gray code. 'a' is also generally used.
    y : '1', str, optional
        Second index in the gray code. 'd' is also generally used.

    Returns
    -------
    graycode_order : list
        The list of indices ordered by their corresponding frequencies.
    """

    graycode_order = [x, y]
    for _ in range(level - 1):
        graycode_order = [x + path for path in graycode_order] + \
                         [y + path for path in graycode_order[::-1]]
    return graycode_order

def transfer_fcn(coeff_filter, freq_domain):
    """
    Transfer function of a given wavelet filter.

    Parameters
    ----------
    coeff_filter : 1d array
        Coefficients of the filter
    freq_domain : 1d array
        Frequencies to be used for calculating the transfer function

    Returns
    -------
    T_fcn : 1d array
        Transfer function of the given filter
    """
    T_fcn = np.zeros(len(freq_domain))*(0+1j)
    for k in range(len(coeff_filter)):
        T_fcn += coeff_filter[k]*np.exp(-1j*2*np.pi*freq_domain*k)
    return T_fcn

def scale_to_frq(f_sample, j):
    """
    Convert scales at a given j to their corresponding center frequencies.

    Parameters
    ----------
    f_sample : float
        Sampling rate of the data.
    j : int
        The decomposition level.

    Returns
    -------
    frq : 1d array
        Center frequencies of all the scales at a given level j
    """
    frq = f_sample/2**(j+1)*(np.arange(2**j) + 0.5)

    return frq

def time_shift(w_j, L, j, scale):
    """
    Time shift the wavelet coefficient so that it matches the features in the original signal
    temporally. Only works if the half length of the wavelet is even and if the wavelet is of the
    LA type (symlet).

    Parameters
    ----------
    w_j : 1d array
        Wavelet coefficient.
    L : int
        Length of the wavelet used to calculate the wavelet coefficient.
    j : int
        Decomposition level.
    scale : int
        The specific scale (n) the wavelet coefficient corresponds to.

    Returns
    -------
    w_j : 1d array
        Wavelet coefficient corrected for its corresponding time shift
    """
    # compute the length of the composite wavelet at level j
    Lj = (2**j-1)*(L-1) + 1
    # find the c_j vector using the graycode
    c_j = find_scale_index(j, x='0', y='1')
    _s = int(0)
    # for the specified scale
    c_j_n = np.array(list(c_j[scale]), dtype='int')
    for _j in range(len(c_j[scale])):
        _s += int(c_j_n[_j]*2**_j)
    v_shift = round(-Lj/2 + (2**(j-1) - _s))
    w_j = np.roll(w_j, abs(v_shift))

    # # old code that seems to calculate all shifts for a given j. unnecessary.
    # S = np.zeros(2**j, dtype='int')
    # for _i in range(len(c_j)):  # len(c_j)==2**j
    #     # convert string to numpy array
    #     c_j[_i] = np.array(list(c_j[_i]), dtype='int')
    #     for _j in range(len(c_j[_i])):  # len(c_j[_i])==j
    #         # for each scale
    #         S[_i] += c_j[_i][_j]*2**_j
    # v_shift = -Lj/2 + (2**(j-1)-S)
    # v_shift = v_shift.astype('int')
    # w_j = np.roll(w_j, abs(v_shift[scale]))

    return w_j


class CompositeFilter():
    r"""Generate composite filter
    Composite filter generated from a given wavelet for a specific decomposition level and
    all the scales involved.

    """

    def __init__(self, g=None, wavelet='D(8)'):
        """
        Parameters
        ----------
        g : 1d array, optional
            Scaling (lowpass) filter for calculating the wavelet filter.
        wavelet : str, optional
            Name of  the scaling filter (g). Default is 'D(8)'. Built-in options are:
            'Haar', 'D(4)', 'D(8)', 'D(10)', 'D(12)', 'LA(8)', 'LA(12)', 'LA(16)'

            - 'D(8)': Daubechies standard wavelet with 8 elements.

            - 'LA(8)': Least-asymmetric variant of the Daubechies wavelet with 8 elements.

            Other wavelets can be found in the following sources:

            - *Wavelet Methods for Time Series Analysis* by Percival and Walden.

            - Import from the ``pywt`` package using ``pywt.Wavelet()``. Note that only half of the
              wavelet length is indicated in the name of the wavelets. E.g., the equivalent of
              'D(8)' and 'LA(8)' in ``pywt`` are 'db4' and 'sym4'.

        Notes
        -----
        The generated filters need to be normalized by :math:`2^{j/2}` for maximum overlap wavelet
        transforms. This is taken into account in ``WaveletTransform``.

        """
        if g is not None:
            self.g = g
        else:
            self.g = np.array(WAVELETS[wavelet])

        # length of the wavelet
        self.L = len(self.g)

        # calculate the wavelet (highpass) filter (h) paired to the given scaling filter (g)
        self.h = np.zeros(self.L)
        for _l in range(self.L):
            self.h[_l] = (-1)**_l*self.g[self.L-1-_l]

        # wavelet filters for the first decomposition level j=1
        self.filter_nodes = {'10': self.g, '11': self.h}

    def u_n(self, n):
        r"""
        Filter corresponds to the scale n of a level j, i.e. a selection of g or h.

        Parameters
        ----------
        n : int
            Index of the scale from a certain decomposition level j,
            :math:`n=0,1 \ldots, 2^{j-1}-1`.

        Returns
        -------
        u_n : 1d array
            Filter coefficients (either g or h).
        """
        if n % 4 in (0, 3):
            u_n = self.g
        elif n % 4 in (1, 2):
            u_n = self.h
        return u_n

    def filter_cascade(self, J):
        """
        Populate the collection of filters corresponding to all the j and n in a cascading scheme
        based on the specified target decomposition level J.

        Parameters
        ----------
        J : int
            Target decomposition level (or the maximum decomposition level). It should be larger
            than 1.
        """
        if J <= 1:
            raise ValueError('J should be larger than 1')

        for j in range(2, J+1):
            # determine the length of the composite filter
            Lj = (2**j-1)*(self.L-1) + 1

            for n in np.arange(2**j):  # all n-scales in level j
                # name of the node/subscale
                path = '%d%d' % (j, n)

                # compute the composite filter u_j
                u_j = np.zeros(Lj)
                _u_j = self.filter_nodes['%d%d' % (j-1, n//2)]
                _u = self.u_n(n)
                for l in range(Lj):  # all coefficients of u_j at n of j
                    for k in range(self.L):
                        index = l-2**(j-1)*k
                        if index < 0 or index >= len(_u_j):
                            u_j[l] += 0
                        else:
                            u_j[l] += _u[k]*_u_j[index]

                print('Processing the node: ' + path)
                self.filter_nodes.update({path: u_j})

    def sqd_gain_fcn(self, j, n, fs=None, max_overlap=False):
        r"""
        Squared gain function of a given wavelet filter at subscale n and decomposition level j.

        Parameters
        ----------
        j : int
            Decomposition level j
        n : int
            Index of the scale from a certain decomposition level j,
            :math:`n=0,1 \ldots, 2^{j-1}-1`.
        fs : None, 1d array, optional
            Frequencies to be used for caclulating the transfer function [0, 1/2).
        max_overlap : False, bool, optional
            If True, the output represent the squared gain from MODWT.

        Returns
        -------
        gain : 1d array
            Squared gain for the specified filter
        fs : 1d array
            Standard frequency range
        """
        if fs is None:
            fs = np.arange(0, 0.5, 1e-3)

        paths = find_scale_index(j, x='0', y='1')
        U_trans = 1
        for k in range(j):
            paths[n] = np.array(list(paths[n]), dtype='int')
            if paths[n][k] == 1:
                U_trans = U_trans*transfer_fcn(self.h, 2**k*fs)
            elif paths[n][k] == 0:
                U_trans = U_trans*transfer_fcn(self.g, 2**k*fs)
        gain = np.abs(U_trans)**2

        if max_overlap:
            gain = gain/2**j

        return gain, fs

    def max_J(self, N):
        """
        Maximum decomposition level for the data series with the given wavelet filter. Although
        there is technically no upper limit for MODWT, it is still recommended to apply the same
        criterion as in the case of DWT.

        Parameters
        ----------
        N : int
            Length of the data series

        Returns
        -------
        max_dec : int
            (Recommended) Maximum decomposition level for the data series with the given wavelet
            filter.
        """
        max_dec = np.log2(N/(self.L-1)+1)

        return int(np.floor(max_dec))

    def save_filterbank(self, full_path_write):
        """
        Save the computed filterbank to the specified local directory.

        Parameters
        ----------
        full_path_write : obj
            Full path to the file, recommended to use "./filterbank_[filtertype]_[j].pkl" to name
            the file.
        """
        pkl_dump(full_path_write, self.filter_nodes)


class WaveletTransform():
    """
    Basic class for discrete wavelet (packet) transform.
    """
    def __init__(self, X, j=None, mode='max-overlap', filterbank=None, full_path_filterbank=None):
        """
        Input either data or correlation matrix.

        Parameters
        ----------
        X : ndarray
            Input data. Two shapes are admissible: 
            - NxM with N being the sample size and M the total
            amount of physical coordinates (for 1d wavelet transform)
            - NxN with N being the sample size (for 2d wavelet transform)
        j : int
            Target decomposition level. A value is expected if 'standard' mode is used.
        mode : 'max-overlap', str. optional
            Two modes of discrete wavelet transform to choose from:
            - 'max-overlap', perform maximum overlap DWT
            - 'standard', perform DWT.
        filterbank : None, dict, optional
            Precalculated filterbank via ``CompositeFilter()``. It is also possible to input custom
            filters in the following format: {Node: Coefficients}. Node has to follow the naming
            convention '[j][0]'.
        full_path_filterbank : None, obj, optional
            Path to the filterbank pickle file created using ``save_filterbank()`` in
            ``CompositeFilter()``.
        """
        if mode == 'max-overlap':
            self.filter_fct = 1
        elif mode == 'standard':
            self.filter_fct = 0
            if j is None:
                raise ValueError('Provide the desired decomposition level j')
            else:
                self.j = j
        else:
            raise ValueError('Choose a valid DWT mode: max-overlap or standard.')

        if len(np.shape(X)) == 1:
            # convert a 0D array to 1D
            X = X.reshape(len(X), 1)
        elif len(np.shape(X)) > 2:
            raise ValueError('Input data needs to be converted into a 2-dimensional matrix first.')

        self.X = np.array(X)
        self.N = X.shape[0]
        self.M = X.shape[1]

        if full_path_filterbank is not None:
            self.filterbank = pkl_load(full_path_filterbank)
        elif filterbank is None:
            raise ValueError('Provide the precomputed filter bank or directory to the local file.')
        else:
            self.filterbank = filterbank

    @property
    def index_W(self):
        """
        Matrix conversion indices for data series and wavelet coefficient.
        """
        if self.filter_fct == 1:  # for MODWT
            index_t = np.broadcast_to(np.arange(self.N), (self.N, self.N)).T
            index_W = np.mod(np.add(index_t, -np.arange(self.N)), self.N)

        elif self.filter_fct == 0:  # for DWT
            if self.N % 2 != 0:
                raise Warning('The sample length of X is not integer multiples of 2**j. It will be'+
                              'truncated.')

            # the elements of X outside of integer multiples of 2**j will be truncated
            N_j = self.N // 2**self.j
            index_t = np.broadcast_to(np.arange(N_j), (self.N, N_j)).T
            index_W = 2**self.j*(index_t + 1)
            index_W = np.mod(np.add(index_W, -np.arange(self.N)), self.N)

        return index_W

    def filter_matrix(self, j, n):
        """
        Convert the composite filter u_j into its matrix form based on the sample size of N 

        Parameters
        ----------
        j : int
            The decomposition level.
        n : int
            The scale in level j.
        index_W : int
            Indices of the matrix form of the wavelet filter.

        Returns
        -------
        u_j_mat : ndarray
            Composite filter u_j in its matrix form
        """
        path = '%d%d' % (j, n)
        # for MODWT a factor of 2**(j/2) needs to be applied to the standard filterbank
        u_j = np.array(self.filterbank[path]) / (2**(j/2))**self.filter_fct
        # the elements of u_j is assumed zero outside of the length of Lj. In practice u_j is
        # shorter than N, its periodized version is derived by appending N-Lj zeros
        u_j_p = np.append(u_j, np.zeros(self.N-len(u_j)))
        # the index_W helps to select the correct elements of u_j_p to form u_j_mat directly
        u_j_mat = u_j_p[self.index_W]

        return u_j_mat

    def wavelet_coeff_1d(self, j, scales):
        """
        Decompose the input data with MODWPT at specified j and with desired scales, along the
        axis=0 dimension of the data (i.e., N in an NxM data). It is possible to combine multiple
        (discrete) scales to achieve desired filtering effect.

        Parameters
        ----------
        j : int
            The decomposition level.
        scales: list, 1d array of int
            All the scales included in the decomposition

        Returns
        -------
        W_j : ndarray
            Wavelet coefficients from decomposing input data. It has a dimension of len(scales)xNxM

        Notes
        -----
        1. Not recommended to set a large number of scales if M is very large
        2. Better to start with a single scale
        """
        W_j = np.zeros([len(scales), self.N, self.M])
        for i, n in enumerate(scales):
            u_j_mat = self.filter_matrix(j, n)
            W_j[i, :] = u_j_mat @ self.X

            print('Processing n=%d' % n)

        return W_j

    def power_spectrum_1d(self, j, scales):
        """
        Calculate the approximated power spectrum of the input data based on ``wavelet_coeff_1d``.

        Parameters
        ----------
        j : int
            The decomposition level.
        scales: list, 1d array of int
            All the scales included in the decomposition.

        Returns
        -------
        power_spectrum : 1d array
            Power spectrum at specific scales of j. It has a size of len(scales).
        """
        W_j = self.wavelet_coeff_1d(j, scales)
        power_spectrum = np.sum(W_j**2, axis=(1, 2))/self.N/self.M

        return power_spectrum

    def detail_bundle_1d(self, js, scales):
        """
        Reconstruct the input data with MODWT at specified j levels and with desired scales, along
        the axis=0 dimension of the data (i.e., N in an NxM data).

        Parameters
        ----------
        js : 1d array of int
            The correponding decomposition levels. Different j levels are admissible.
        scales: 1d array of int
            All the scales included in the reconstruction. Discrete scales are admissible.

        Returns
        -------
        B_details : ndarray
            Reconstructed detail bundle of input data with the dimension of NxM.

        Notes
        -----
        `js` and `scales` need to have the same length.
        """
        if self.filter_fct == 0:
            # since DWT doesn't allow mixed j levels
            js = int(np.ones(len(scales))*self.j)

        if len(js) != len(scales):
            raise ValueError('js and scales need to have the same length!')

        # construct the transfer operator T_oper
        T_oper = np.zeros([self.N, self.N])
        for j, n in zip(js, scales):
            u_j_mat = self.filter_matrix(j, n)
            T_oper += u_j_mat.T @ u_j_mat
            print('Processing j=%d and n=%d' % (j,n))

        B_details = T_oper @ self.X

        return B_details

    def wavelet_coeff_2d(self, j, scales):
        """
        Decompose the input data with MODWT at specified j and with desired scales, along both axes
        of the data dimension. Only data with a shape of NxN is admissible due to the precomputed
        Designed for handling cross-correlation matrices.

        Parameters
        ----------
        j : int
            The decomposition level.
        scales : list, 1d array of int
            All the scales included in the decomposition.

        Returns
        -------
        W_j : ndarray
            Wavelet coefficients from decomposing input data. It has a dimension of len(scales)xNxN.
        """
        if self.N != self.M:
            raise ValueError('Input data must have a square shape.')

        W_j = np.zeros([len(scales), self.N, self.N])
        for i, n in enumerate(scales):
            u_j_mat = self.filter_matrix(j, n)
            # perform wavelet transform along each axis of the input data sequentially.
            W_j[i, :] = u_j_mat @ self.X @ u_j_mat.T

            print('Processing n=%d' % n)

        return W_j

    def power_spectrum_2d(self, j, scales):
        """
        Calculate the power spectrum of the input data based on ``wavelet_coeff_2d``.

        Parameters
        ----------
        j : int
            The decomposition level.
        scales: list, 1d array of int
            All the scales included in the decomposition.

        Returns
        -------
        power_spectrum : 1d array
            Power spectrum at specific scales of j. It has a size of len(scales).
        """
        W_j = self.wavelet_coeff_2d(j, scales)
        power_spectrum = np.sum(W_j**2, axis=(1, 2))/self.N/self.N

        return power_spectrum

    def detail_bundle_2d(self, js, scales):
        """
        Reconstruct the input data (2d) with MODWPT at specified j and with desired scales.

        Parameters
        ----------
        js : 1d arrray of int
            The correponding decomposition levels. Different j levels are admissible.
        scales: 1d array of int
            All the scales included in the reconstruction. Discrete scales are admissible.

        Returns
        -------
        K_mat : ndarray
            Reconstructed input data with the dimension of NxN.
        T_mat: ndarray
            Transform matrix used to perform the transform.
        """
        if len(js) != len(scales):
            raise ValueError('js and scales need to have the same length!')
        if not self.N == self.M:
            raise ValueError('Input data must have a square shape')

        # construct the transfor operator T_oper
        T_oper = np.zeros([len(scales), 1, self.N, self.N])
        _count = 0
        for j, n in zip(js, scales):
            u_j_mat = self.filter_matrix(j, n)
            T_oper[_count, 0, :] = u_j_mat.T @ u_j_mat
            _count += 1

            print('Processing n=%d' % n)

        # only the first two axes of T_oper should be swapped in the transpose
        T_oper_transpose = np.transpose(T_oper, (1, 0, 2, 3))
        # reconstruct X along each its axis sequentially
        _K = T_oper @ self.X @ T_oper_transpose
        # sum over all scales used to obtain the final result
        K_mat = np.sum(_K, axis=(0, 1))

        return K_mat, T_oper
