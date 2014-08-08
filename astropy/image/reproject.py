import numpy as np

from ..wcs.utils import wcs_to_celestial_frame
from astropy.coordinates import UnitSphericalRepresentation
from .. import units as u


def reproject_image_2d(array, wcs_in, wcs_out, shape_out, mode='nearest', order=1):
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
    mode : {``'interpolation'``, ``'drizzle'``}
        The type of reprojection
    order : int
        The order of the interpolation (if ``mode`` is set to
        ``'interpolation'``). A value of ``0`` indicates nearest neighbor
        interpolation (the default).
    """

    # Find input/output frames
    frame_in = wcs_to_celestial_frame(wcs_in)
    frame_out = wcs_to_celestial_frame(wcs_out)

    # Defining pixel to pixel transformations

    def pixel_in_to_pixel_out(xp_in, yp_in):

        xw, yw = wcs_in.wcs_pix2world(xp_in, yp_in, 0)

        xw_unit_in = u.Unit(wcs_in.wcs.cunit[0])
        yw_unit_in = u.Unit(wcs_in.wcs.cunit[1])

        # TODO: for now assuming that coordinates are spherical, not
        # necessarily the case. Also assuming something about the order of the
        # arguments. Also assuming units.
        data = UnitSphericalRepresentation(xw * xw_unit_in,
                                           yw * yw_unit_in)
        coords_in = frame_in.realize_frame(data)
        coords_out = coords_in.transform_to(frame_out)

        xw_unit_out = u.Unit(wcs_out.wcs.cunit[0])
        yw_unit_out = u.Unit(wcs_out.wcs.cunit[1])

        xw = coords_out.spherical.lon.to(xw_unit_out).value
        yw = coords_out.spherical.lat.to(yw_unit_out).value

        xp_out, yp_out = wcs_out.wcs_world2pix(xw, yw, 0)
        return xp_out, yp_out

    def pixel_out_to_pixel_in(xp_out, yp_out):

        xw, yw = wcs_out.wcs_pix2world(xp_out, yp_out, 0)

        xw_unit_out = u.Unit(wcs_out.wcs.cunit[0])
        yw_unit_out = u.Unit(wcs_out.wcs.cunit[1])

        # TODO: for now assuming that coordinates are spherical, not
        # necessarily the case. Also assuming something about the order of the
        # arguments. Also assuming units.
        data = UnitSphericalRepresentation(xw * xw_unit_out,
                                           yw * yw_unit_out)
        coords_in = frame_out.realize_frame(data)
        coords_out = coords_in.transform_to(frame_in)

        xw_unit_in = u.Unit(wcs_in.wcs.cunit[0])
        yw_unit_in = u.Unit(wcs_in.wcs.cunit[1])

        xw = coords_out.spherical.lon.to(xw_unit_in).value
        yw = coords_out.spherical.lat.to(yw_unit_in).value

        xp_in, yp_in = wcs_in.wcs_world2pix(xw, yw, 0)
        return xp_in, yp_in

    if mode == 'interpolation':

        from scipy.ndimage import map_coordinates

        # Generate pixel coordinates of output image
        xp_out = np.arange(shape_out[1])
        yp_out = np.arange(shape_out[0])
        xp_out_grid, yp_out_grid = np.meshgrid(xp_out, yp_out)

        # Convert to pixel coordinates in input frmae
        xp_in_grid, yp_in_grid = pixel_out_to_pixel_in(xp_out_grid, yp_out_grid)

        # Interpolate values to new grid
        coordinates = [yp_in_grid.ravel(), xp_in_grid.ravel()]
        array_new = map_coordinates(array, coordinates,
                                    order=order, cval=np.nan,
                                    mode='constant').reshape(shape_out)

    elif mode == 'drizzle':

        # Generate pixel coordinates of input image pixel corners
        xp_in = np.linspace(0, array.shape[1], array.shape[1]+1) - 0.5
        yp_in = np.linspace(0, array.shape[0], array.shape[0]+1) - 0.5
        xp_in_grid, yp_in_grid = np.meshgrid(xp_out, yp_out)

        # Convert pixel coordainates to frame of reference of output image
        xp_out_grid, yp_out_grid = pixel_in_to_pixel_out(xp_in_grid, yp_in_grid)

        # Call drizzle algiorithm
        # TODO: magic here
        raise NotImplementedError("Drizzle mode not yet implemented")

    else:

        raise ValueError("mode= should be one of 'interpolation' or 'drizzle'")

    return array_new
