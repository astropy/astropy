# Functions/classes for WCSAxes related to APE14 WCSes

import numpy as np

from astropy.coordinates import SkyCoord, ICRS
from astropy import units as u
from astropy.wcs import WCS
from astropy.wcs.wcsapi import SlicedLowLevelWCS

from .frame import EllipticalFrame
from .transforms import CurvedTransform

__all__ = ['transform_coord_meta_from_wcs', 'WCSWorld2PixelTransform',
           'WCSPixel2WorldTransform']

IDENTITY = WCS(naxis=2)
IDENTITY.wcs.ctype = ["X", "Y"]
IDENTITY.wcs.crval = [0., 0.]
IDENTITY.wcs.crpix = [1., 1.]
IDENTITY.wcs.cdelt = [1., 1.]

# TODO: Make a ReversedPixelWCS thingy


def transform_coord_meta_from_wcs(wcs, frame_class, aslice=None):
    coord_meta = {}
    coord_meta['name'] = []
    coord_meta['type'] = []
    coord_meta['wrap'] = []
    coord_meta['unit'] = []
    coord_meta['format_unit'] = []

    invert_xy = False
    if aslice is not None:
        wcs_slice = list(aslice)
        wcs_slice[wcs_slice.index("x")] = slice(None)
        wcs_slice[wcs_slice.index("y")] = slice(None)
        wcs = SlicedLowLevelWCS(wcs, wcs_slice[::-1])

        invert_xy = aslice.index('x') < aslice.index('y')


    coord_meta['default_position'] = [''] * wcs.world_n_dim

    m = wcs.axis_correlation_matrix.copy()
    if invert_xy:
        m = m[:,::-1]

    for i, spine_name in enumerate(frame_class.spine_names[:2]):
        pos = np.nonzero(m[:, i % 2])[0]
        if len(pos) > 0:
            coord_meta['default_position'][pos[0]] = frame_class.spine_names[i]
            m[pos[0], :] = 0

    transform = WCSPixel2WorldTransform(wcs, invert_xy=invert_xy)

    for idx in range(wcs.world_n_dim):

        axis_type = wcs.world_axis_physical_types[idx]
        axis_unit = u.Unit(wcs.world_axis_units[idx])
        coord_wrap = None
        format_unit = axis_unit

        coord_type = 'scalar'

        if axis_type is not None:

            axis_type_split = axis_type.split('.')

            if "pos.helioprojective.lon" in axis_type:
                coord_wrap = 180.
                format_unit = u.arcsec
            elif "pos.helioprojective.lat" in axis_type:
                format_unit = u.arcsec
            elif "pos" in axis_type_split:
                if "lon" in axis_type_split:
                    coord_type = "longitude"
                elif "lat" in axis_type_split:
                    coord_type = "latitude"
                elif "ra" in axis_type_split:
                    coord_type = "longitude"
                    format_unit = u.hourangle
                elif "dec" in axis_type_split:
                    coord_type = "latitude"
                elif "alt" in axis_type_split:
                    coord_type = "longitude"
                elif "az" in axis_type_split:
                    coord_type = "latitude"
                elif "long" in axis_type_split:
                    coord_type = "longitude"

        coord_meta['type'].append(coord_type)
        coord_meta['wrap'].append(coord_wrap)
        coord_meta['format_unit'].append(format_unit)
        coord_meta['unit'].append(axis_unit)
        coord_meta['name'].append(axis_type or '')

    return transform, coord_meta


def wcsapi_to_celestial_frame(wcs):
    for cls, args, kwargs in wcs.world_axis_object_classes.values():
        if cls is SkyCoord:
            return kwargs.get('frame', ICRS())


class WCSWorld2PixelTransform(CurvedTransform):
    """
    WCS transformation from world to pixel coordinates
    """

    has_inverse = True
    frame_in = None

    def __init__(self, wcs, invert_xy=False):

        super().__init__()

        if wcs.pixel_n_dim != 2:
            raise ValueError('Only pixel_n_dim==2 is supported')

        self.wcs = wcs
        self.invert_xy = invert_xy

        self.frame_in = wcsapi_to_celestial_frame(wcs)

    def __eq__(self, other):
        # TODO: for performance, find a way to return True when they do match
        return False

    @property
    def input_dims(self):
        return self.wcs.world_n_dim

    def transform(self, world):

        # Convert to a list of arrays
        world = list(world.T)

        if len(world) != self.wcs.world_n_dim:
            raise ValueError(f"Expected {self.wcs.world_n_dim} world coordinates, got {len(world)} ")

        if len(world[0]) == 0:
            pixel = np.zeros((0, 2))
        else:
            pixel = self.wcs.world_to_pixel_values(*world)

        if self.invert_xy:
            pixel = pixel[::-1]

        pixel = np.array(pixel).T

        return pixel

    def inverted(self):
        """
        Return the inverse of the transform
        """
        return WCSPixel2WorldTransform(self.wcs, invert_xy=self.invert_xy)


class WCSPixel2WorldTransform(CurvedTransform):
    """
    WCS transformation from pixel to world coordinates
    """

    has_inverse = True

    def __init__(self, wcs, invert_xy=False):

        super().__init__()

        if wcs.pixel_n_dim != 2:
            raise ValueError('Only pixel_n_dim==2 is supported')

        self.wcs = wcs
        self.invert_xy = invert_xy

        self.frame_out = wcsapi_to_celestial_frame(wcs)

    def __eq__(self, other):
        # TODO: for performance, find a way to return True when they do match
        return False

    @property
    def output_dims(self):
        return self.wcs.world_n_dim

    def transform(self, pixel):

        # Convert to a list of arrays
        pixel = list(pixel.T)

        if len(pixel) != self.wcs.pixel_n_dim:
            raise ValueError(f"Expected {self.wcs.pixel_n_dim} world coordinates, got {len(pixel)} ")

        if self.invert_xy:
            pixel = pixel[::-1]

        if len(pixel[0]) == 0:
            world = np.zeros((0, self.wcs.world_n_dim))
        else:
            world = self.wcs.pixel_to_world_values(*pixel)

        # At the moment, one has to manually check that the transformation
        # round-trips, otherwise it should be considered invalid.
        # pixel_check = self.wcs.pixel_to_world_values(*world)
        # with np.errstate(invalid='ignore'):
        #     for ipix in range(len(pixel)):
        #         invalid = np.any(np.abs(pixel_check[ipix] - pixel[ipix]) > 1., axis=1)
        #     world[invalid] = np.nan

        world = np.array(world).T

        return world

    transform_non_affine = transform

    def inverted(self):
        """
        Return the inverse of the transform
        """
        return WCSWorld2PixelTransform(self.wcs, invert_xy=self.invert_xy)
