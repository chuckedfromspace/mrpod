Greetings
=========

Welcome to the documentation of the Python module ``mrpod`` for performing
Multiresolution Proper Orthogonal Decomposition (MRPOD) of multi-dimensional
time series.

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

Why MODWT
^^^^^^^^^

``mrpod`` was developed based primarily on the following two criteria:
    - Dynamics with various frequencies can be identified and adequately
    isolated.
    - Discontinuities in temporal behaviors can be properly resolved and align
    perfectly with the original data series for appropriate comparison.

Unlike the classical DWT, shift-invariant DWT such as MODWT is well-defined for
arbitrary sample sizes and is not sensitive to the "break-in" point in the time
series. Additionally, MODWT affords a more "square-looking" gain response and
hence a higher spectral isolation especially for cases with dynamics densely
packed in the frequency domain. Although MODWT sacrifices orthonormality, it can
still carry out an exact analysis of variance as well as a perfect
reconstruction of the time series.

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