Greetings
=========

Welcome to the documentation of the Python module ``mrpod`` (`GitHub depository
<https://github.com/chuckedfromspace/mrpod>`_) for performing
Multiresolution Proper Orthogonal Decomposition (MRPOD) of multi-dimensional
data series (vector and scalar) obtained in turbulent flows.

Why ``mrpod``
^^^^^^^^^^^^^

Data-driven, modal-decomposition techniques such as POD (also known as the
Principle Component Analysis, PCA) and Dynamic Mode Decomposition (DMD) have
been widely implemented to extract periodic, coherent structures in turbulent
flows. In reacting flows such as confined turbulent flames (in gas turbines),
however, we are often confronted with both periodic (e.g., hydrodynamic and
acoustic instabilities) and non-periodic dynamics (e.g., flame lift-off,
flashback, and bistability), which can coexist over a wide range of time scales
and may even interact with each other. 

In their most rudimentary forms, POD and DMD have proven insufficient at
separating these dynamics while resolving their temporal behaviors (such as
discontinuities) at the same time. On the other hand, by marrying the concept of
wavelet-based Multiresolution Analysis (MRA) with standard modal decompositions,
multiresolution DMD ([MRDMD]_) and multi-scale POD ([mPOD]_) have demonstrated
robust capabilities at identifying unsteady dynamics and discontinuities in time
series. 

Based on a similar concept, multiresolution POD ([MRPOD]_) has been developed by
combining Maximum-Overlap Discrete Wavelet Transform (MODWT) with conventional
snapshot POD. MRPOD has been successfully applied to time series of velocity
(vector) and scalar fields obtained by kHz-rate laser diagnostics
(e.g., Particle Image Velocimetry or PIV, Planar Laser-induced Fluorescence or
PLIF) in the so-called bistable turbulent swirl flame, a common phenomenon
encountered in gas turbines.

Bistable turbulent swirl flame
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Why MODWT
^^^^^^^^^

Why not ``pywt``
^^^^^^^^^^^^^^^^
Instead of using the existing Python library ``pywt`` to carry out wavelet
transform, a matrix-operation based routine was written from the ground up
specifically for more efficient 1-D and 2-D wavelet decomposition/reconstruction
of multi-dimensional data series stored in ndarrays. Although several commonly
used wavelet filters are built into ``mrpod``, the vast library of wavelet
filters in ``pywt`` should be taken advantage of when constructing custom 
composite filters using the filter-cascading method in ``mrpod``.

----------------

.. [MRDMD] Kutz, J., Fu, X., Brunton, S. Multiresolution dynamic mode
    decomposition. *SIAM Journal on Applied Dynamical Systems* 15 (2), 713-735,
    2016.

.. [mPOD] Mendez, M. A., Balabane, M., Buchlin, J. M. Multi-scale proper 
    orthogonal decomposition of complex fluid flows.
    *Journal of Fluid Mechanics* 870, 988-1036, 2019.

.. [MRPOD] Yin, Z., Stöhr, M. Time–Frequency Localisation of Intermittent
    Dynamics in a Bistable Turbulent Swirl Flame. *Journal of Fluid Mechanics*
    882, A30, 2020.