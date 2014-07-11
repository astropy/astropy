import numpy as np


def reproject_image_2d(array, wcs_in, wcs_out, shape_out, mode='nearest'):
    """
    Reproject a 2D array from one WCS to another.

    Parameters
    ----------
    array : :class:`~numpy.ndarray`
        The array to reproject
    wcs_in : :class:`~astropy.wcs.WCS`
        The input WCS
    wcs_out : :class:`~astropy.wcs.WCS`
        The output WCS
    shape_out : tuple
        The shape of the output array
    """

    # Start off by defining pixel to pixel transformations

    def pixel_out_to_pixel_in(xp_out, yp_out):
        xw, yw = wcs_out.wcs_pix2world(xp_out, yp_out, 0)
        # TODO: currently assume coordinate frames match, fix this!
        xp_in, yp_in = wcs_in.wcs_world2pix(xw, yw, 0)
        return xp_in, yp_in

    if mode in ['nearest']:

        from scipy.ndimage import map_coordinates

        # Generate pixel coordinates of output image
        xp_out = np.arange(shape_out[1])
        yp_out = np.arange(shape_out[0])
        xp_out_grid, yp_out_grid = np.meshgrid(xp_out, yp_out)

        # Convert to pixel coordinates in input frmae
        xp_in_grid, yp_in_grid = pixel_out_to_pixel_in(xp_out_grid, yp_out_grid)

        # Interpolate values to new grid
        coordinates = [xp_in_grid.ravel(), yp_in_grid.ravel()]
        array_new = map_coordinates(array, coordinates, order=3).reshape(shape_out)

    else:

        raise ValueError("mode= should be one of 'nearest")

    return array_new
