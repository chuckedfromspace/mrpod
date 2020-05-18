"""
Generate collections of scattered points in Cartesian coordinates.
"""
import numpy as np

def scatter_circular(size, R1=3, R2=None, angle=0, reset_random=True):
    """
    Generate a collection of random points distributed in a circle or elipse.

    Parameters
    ----------
    size : int
        Number of points
    R1 : float
        Radius of the circle
    R2 : float
        If not the same as R1, then the output is an elipse
    angle : 0, float, optional
        Angle to rotate the coordinates of all the points, should be [0, 180)

    Returns
    -------
    2d array
        A collection of scattered points with the shape of 2xsize.
    """
    if reset_random:
        np.random.seed(0)

    _a = np.random.random_sample(size)
    _b = np.random.random_sample(size)

    points = np.zeros((2, size))
    for i in range(size):
        a = min(_a[i], _b[i])
        b = max(_a[i], _b[i])
        points[0, i] = b*R1*np.cos(2*np.pi*a/b)
        points[1, i] = b*R2*np.sin(2*np.pi*a/b)

    radian = angle/180*np.pi
    rotate_mat = np.array([[np.cos(radian), -np.sin(radian)], [np.sin(radian), np.cos(radian)]])
    points = rotate_mat @ points

    return points

def scatter_3D(size, **kwargs):
    """
    Generate a collection of random points in a 3-D Cartesian coordinate.

    Parameters
    ----------
    size : int
        Number of points

    Other Parameters
    ----------------
    **kwargs
        This method takes keyward arguments from ``scatter_circle``

    Returns
    -------
    2d array
        A collection of scattered points with the shape of 3xsize.
    """
    points = scatter_circular(size, **kwargs)
    # _z = np.random.random_sample(size)
    _z = np.random.normal(0, 0.1, size)

    return np.vstack((points, _z))
