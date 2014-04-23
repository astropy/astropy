import numpy as np

from matplotlib.lines import Path

from astropy.coordinates.angle_utilities import angular_separation

# Tolerance for WCS round-tripping
ROUND_TRIP_TOL = 1e-1

# Tolerance for discontinuities relative to the median
DISCONT_FACTOR = 10.


def get_lon_lat_path(ax, transform, lon_lat):
    """
    Draw a curve, taking into account discontinuities.

    Parameters
    ----------
    ax : ~matplotlib.axes.Axes
        The axes in which to plot the grid
    transform : transformation class
        The transformation between the world and pixel coordinates
    lon_lat : `~numpy.ndarray`
        The longitude and latitude values along the curve, given as a (n,2)
        array.
    """

    # Get pixel limits
    # xlim = ax.get_xlim()
    # ylim = ax.get_ylim()

    # Transform line to pixel coordinates
    pixel = transform.inverted().transform(lon_lat)

    # In some spherical projections, some parts of the curve are 'behind' or
    # 'in front of' the plane of the image, so we find those by reversing the
    # transformation and finding points where the result is not consistent.

    lon_lat_check = transform.transform(pixel)

    sep = angular_separation(np.radians(lon_lat[:, 0]),
                             np.radians(lon_lat[:, 1]),
                             np.radians(lon_lat_check[:, 0]),
                             np.radians(lon_lat_check[:, 1]))

    sep[sep > np.pi] -= 2. * np.pi

    mask = np.abs(sep > ROUND_TRIP_TOL)

    # Mask values with invalid pixel positions
    mask = mask | np.isnan(pixel[:, 0]) | np.isnan(pixel[:, 1])

    # Mask values outside the viewport
    # This has now been disabled because it assumes specifically rectangular
    # axes, and also doesn't work if the coordinate direction is flipped.
    # outside = ((pixel[:, 0] < xlim[0]) | (pixel[:, 0] > xlim[-1]) |
    #            (pixel[:, 1] < ylim[0]) | (pixel[:, 1] > ylim[-1]))
    # mask[1:-1] = mask[1:-1] | (outside[2:] & outside[:-2])

    # We can now start to set up the codes for the Path.
    codes = np.zeros(lon_lat.shape[0], dtype=np.uint8)
    codes[:] = Path.LINETO
    codes[0] = Path.MOVETO
    codes[mask] = Path.MOVETO

    # Also need to move to point *after* a hidden value
    codes[1:][mask[:-1]] = Path.MOVETO

    # We now go through and search for discontinuities in the curve that would
    # be due to the curve going outside the field of view, invalid WCS values,
    # or due to discontinuities in the projection.

    # We start off by pre-computing the step in pixel coordinates from one
    # point to the next. The idea is to look for large jumps that might indicate
    # discontinuities.
    step = np.sqrt((pixel[1:, 0] - pixel[:-1, 0]) ** 2 +
                   (pixel[1:, 1] - pixel[:-1, 1]) ** 2)

    # We search for discontinuities by looking for places where the step
    # is larger by more than a given factor compared to the median
    # discontinuous = step > DISCONT_FACTOR * np.median(step)
    discontinuous = step[1:] > DISCONT_FACTOR * step[:-1]

    # Skip over discontinuities
    codes[2:][discontinuous] = Path.MOVETO

    # The above missed the first step, so check that too
    if step[0] > DISCONT_FACTOR * step[1]:
        codes[1] = Path.MOVETO

    # Create the path
    path = Path(pixel, codes=codes)

    # And add to the axes
    return path


def get_gridline_path(ax, transform, world):
    """
    Draw a grid line

    Parameters
    ----------
    ax : ~matplotlib.axes.Axes
        The axes in which to plot the grid
    transform : transformation class
        The transformation between the world and pixel coordinates
    world : `~numpy.ndarray`
        The world coordinates along the curve, given as a (n,2)
        array.
    """

    # Get pixel limits
    # xlim = ax.get_xlim()
    # ylim = ax.get_ylim()

    # Transform line to pixel coordinates
    pixel = transform.inverted().transform(world)

    # Mask values with invalid pixel positions
    mask = np.isnan(pixel[:, 0]) | np.isnan(pixel[:, 1])

    # Mask values outside the viewport
    # This has now been disabled because it assumes specifically rectangular
    # axes, and also doesn't work if the coordinate direction is flipped.
    # outside = ((pixel[:, 0] < xlim[0]) | (pixel[:, 0] > xlim[-1]) |
    #            (pixel[:, 1] < ylim[0]) | (pixel[:, 1] > ylim[-1]))
    # mask[1:-1] = mask[1:-1] | (outside[2:] & outside[:-2])

    # We can now start to set up the codes for the Path.
    codes = np.zeros(world.shape[0], dtype=np.uint8)
    codes[:] = Path.LINETO
    codes[0] = Path.MOVETO
    codes[mask] = Path.MOVETO

    # Also need to move to point *after* a hidden value
    codes[1:][mask[:-1]] = Path.MOVETO

    # We now go through and search for discontinuities in the curve that would
    # be due to the curve going outside the field of view, invalid WCS values,
    # or due to discontinuities in the projection.

    # Create the path
    path = Path(pixel, codes=codes)

    # And add to the axes
    return path
