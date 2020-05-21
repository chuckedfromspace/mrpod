"""
An example to generate a simple dataset of time-evolving vortex shedding
"""
import numpy as np

def phantomGaussian(size, fwhm=3, center=None, scale=1, offset_x=0):
    """
    Create a Gaussian shape inside a square (soft edge shapes)

    Parameters
    ----------
    size: int
        Number of pixels of the side of a square.
    fwhm: int
        Number of pixels of the FWHM of the Gaussian shape.
    center: array_like, [int, int]
        Center of the shape. Represent pixel locations.
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

def precessing_vortex_core(size=120, R=60, offset=0):
    """Create a precessing_vortex_core
    """
    v_y = (
        phantomGaussian(size, fwhm=30, center=[18+offset+R, 20+offset])  # righthand side
        - phantomGaussian(size, fwhm=10, center=[5+offset+R, 0+offset])*0.6
        - phantomGaussian(size, fwhm=25, center=[36+offset+R, 10+offset])*0.7
        - phantomGaussian(size, fwhm=30, center=[30+offset+R, 50+offset])
        + phantomGaussian(size, fwhm=35, center=[50+offset+R, 62+offset])
        - phantomGaussian(size, fwhm=30, center=[-18-offset+R, 20+offset])  # lefhand side
        + phantomGaussian(size, fwhm=10, center=[-5-offset+R, 0+offset])*0.6
        + phantomGaussian(size, fwhm=25, center=[-36-offset+R, 10+offset])*0.7
        + phantomGaussian(size, fwhm=30, center=[-30-offset+R, 50+offset])
        - phantomGaussian(size, fwhm=35, center=[-50-offset+R, 62+offset])
        )

    v_x = (
        - phantomGaussian(size, fwhm=30, center=[10+offset+R, 8+offset])  # righthand side
        - phantomGaussian(size, fwhm=15, center=[36+offset+R, 8+offset])
        + phantomGaussian(size, fwhm=35, center=[36+offset+R, 40+offset])
        - phantomGaussian(size, fwhm=45, center=[45+offset+R, 68+offset])
        - phantomGaussian(size, fwhm=30, center=[-10-offset+R, 8+offset])  # lefthand side
        - phantomGaussian(size, fwhm=15, center=[-36-offset+R, 8+offset])
        + phantomGaussian(size, fwhm=35, center=[-36-offset+R, 40+offset])
        - phantomGaussian(size, fwhm=45, center=[-45-offset+R, 68+offset])
        )

    return np.array([v_x, v_y])

def toroidal_vortex(size=120, R=60, offset=0, scale=2):
    """Create a toroidal vortex
    """
    v_y = (
        - phantomGaussian(size, fwhm=30, center=[18+offset, 20+offset], scale=scale, offset_x=R)  # righthand side
        + phantomGaussian(size, fwhm=10, center=[5+offset, 0+offset], scale=scale, offset_x=R)*0.6
        + phantomGaussian(size, fwhm=25, center=[36+offset, 10+offset], scale=scale, offset_x=R)*0.7
        + phantomGaussian(size, fwhm=30, center=[30+offset, 50+offset], scale=scale, offset_x=R)
        - phantomGaussian(size, fwhm=35, center=[50+offset, 62+offset], scale=scale, offset_x=R)
        + phantomGaussian(size, fwhm=40, center=[66+offset, 74+offset], scale=scale, offset_x=R)*0.5
        - phantomGaussian(size, fwhm=30, center=[50+offset, 94+offset], scale=scale, offset_x=R)*0.5
        - phantomGaussian(size, fwhm=30, center=[-18-offset, 20+offset], scale=scale, offset_x=R)  # lefhand side
        + phantomGaussian(size, fwhm=10, center=[-5-offset, 0+offset], scale=scale, offset_x=R)*0.6
        + phantomGaussian(size, fwhm=25, center=[-36-offset, 10+offset], scale=scale, offset_x=R)*0.7
        + phantomGaussian(size, fwhm=30, center=[-30-offset, 50+offset], scale=scale, offset_x=R)
        - phantomGaussian(size, fwhm=35, center=[-50-offset, 62+offset], scale=scale, offset_x=R)
        + phantomGaussian(size, fwhm=40, center=[-66-offset, 74+offset], scale=scale, offset_x=R)*0.5
        - phantomGaussian(size, fwhm=30, center=[-50-offset, 94+offset], scale=scale, offset_x=R)*0.5
        )

    v_x = (
        + phantomGaussian(size, fwhm=30, center=[10+offset, 8+offset], scale=scale, offset_x=R)  # righthand side
        + phantomGaussian(size, fwhm=15, center=[36+offset, 8+offset], scale=scale, offset_x=R)
        - phantomGaussian(size, fwhm=35, center=[36+offset, 40+offset], scale=scale, offset_x=R)
        + phantomGaussian(size, fwhm=45, center=[45+offset, 68+offset], scale=scale, offset_x=R)
        - phantomGaussian(size, fwhm=45, center=[55+offset, 96+offset], scale=scale, offset_x=R)*0.5
        - phantomGaussian(size, fwhm=30, center=[-10-offset, 8+offset], scale=scale, offset_x=R)  # lefthand side
        - phantomGaussian(size, fwhm=15, center=[-36-offset, 8+offset], scale=scale, offset_x=R)
        + phantomGaussian(size, fwhm=35, center=[-36-offset, 40+offset], scale=scale, offset_x=R)
        - phantomGaussian(size, fwhm=45, center=[-45-offset, 68+offset], scale=scale, offset_x=R)
        + phantomGaussian(size, fwhm=45, center=[-55-offset, 96+offset], scale=scale, offset_x=R)*0.5
        )

    return np.array([v_x, v_y])
