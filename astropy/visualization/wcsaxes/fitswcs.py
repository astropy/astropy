# Functions/classes for WCSAxes related to astropy.wcs.WCSLIB

import numpy as np

from astropy import units as u
from astropy.wcs import WCS
from astropy.wcs.utils import wcs_to_celestial_frame

from .transforms import Pixel2WorldTransform, World2PixelTransform

__all__ = ['transform_coord_meta_from_wcs', 'WCSWorld2PixelTransform',
           'WCSPixel2WorldTransform']

IDENTITY = WCS(naxis=2)
IDENTITY.wcs.ctype = ["X", "Y"]
IDENTITY.wcs.crval = [0., 0.]
IDENTITY.wcs.crpix = [1., 1.]
IDENTITY.wcs.cdelt = [1., 1.]


def transform_coord_meta_from_wcs(wcs, slice=None):

    transform = WCSPixel2WorldTransform(wcs, slice=slice)

    coord_meta = {}
    coord_meta['name'] = []
    coord_meta['type'] = []
    coord_meta['wrap'] = []
    coord_meta['unit'] = []
    coord_meta['format_unit'] = []

    for coord_index in range(wcs.wcs.naxis):

        ctype = wcs.wcs.ctype[coord_index]

        coord_wrap = None
        format_unit = None

        if ctype[:4] == 'RA--':
            coord_type = 'longitude'
            format_unit = u.hourangle
        elif ctype[:4] == 'HPLN':
            coord_type = 'longitude'
            format_unit = u.arcsec
            coord_wrap = 180.
        elif ctype[:4] == 'HPLT':
            coord_type = 'latitude'
            format_unit = u.arcsec
        elif ctype[:4] == 'HGLN':
            coord_type = 'latitude'
            coord_wrap = 180.
        elif ctype[1:4] == 'LON' or ctype[2:4] == 'LN':
            coord_type = 'longitude'
        elif ctype[:4] == 'DEC-' or ctype[1:4] == 'LAT' or ctype[2:4] == 'LT':
            coord_type = 'latitude'
        else:
            coord_type = 'scalar'

        coord_meta['type'].append(coord_type)
        coord_meta['wrap'].append(coord_wrap)
        coord_meta['format_unit'].append(format_unit)

        coord_meta['unit'].append(wcs.wcs.cunit[coord_index])
        coord_meta['name'].append(wcs.wcs.ctype[coord_index][:4].replace('-', ''))

    return transform, coord_meta


class WCSWorld2PixelTransform(World2PixelTransform):
    """
    WCS transformation from world to pixel coordinates
    """

    def __init__(self, wcs, slice=None):
        super().__init__()
        self.wcs = wcs
        if self.wcs.wcs.naxis > 2:
            if slice is None:
                raise ValueError("WCS has more than 2 dimensions, so ``slice`` should be set")
            elif len(slice) != self.wcs.wcs.naxis:
                raise ValueError("slice should have as many elements as WCS "
                                 "has dimensions (should be {})".format(self.wcs.wcs.naxis))
            else:
                self.slice = slice
                self.x_index = slice.index('x')
                self.y_index = slice.index('y')
        else:
            self.slice = None
        if wcs.has_celestial:
            self.frame_in = wcs_to_celestial_frame(wcs)

    def __eq__(self, other):
        return (isinstance(other, type(self)) and self.wcs == other.wcs and
                self.slice == other.slice)

    @property
    def input_dims(self):
        return self.wcs.wcs.naxis

    def transform(self, world):

        if world.shape[1] != self.wcs.wcs.naxis:
            raise ValueError("Second dimension of input values should "
                             "match number of WCS coordinates")

        if world.shape[0] == 0:
            pixel = np.zeros((0, 2))
        else:
            pixel = self.wcs.wcs_world2pix(world, 1) - 1

        if self.slice is None:
            return pixel
        else:
            return pixel[:, (self.x_index, self.y_index)]

    transform_non_affine = transform

    def inverted(self):
        return WCSPixel2WorldTransform(self.wcs, slice=self.slice)


class WCSPixel2WorldTransform(Pixel2WorldTransform):
    """
    WCS transformation from pixel to world coordinates
    """

    def __init__(self, wcs, slice=None):
        super().__init__()
        self.wcs = wcs
        self.slice = slice
        if self.slice is not None:
            self.x_index = slice.index('x')
            self.y_index = slice.index('y')
        if wcs.has_celestial:
            self.frame_out = wcs_to_celestial_frame(wcs)

    def __eq__(self, other):
        return (isinstance(other, type(self)) and self.wcs == other.wcs and
                self.slice == other.slice)

    @property
    def output_dims(self):
        return self.wcs.wcs.naxis

    def transform(self, pixel):

        if self.slice is None:
            pixel_full = pixel.copy()
        else:
            pixel_full = []
            for index in self.slice:
                if index == 'x':
                    pixel_full.append(pixel[:, 0])
                elif index == 'y':
                    pixel_full.append(pixel[:, 1])
                else:
                    pixel_full.append(index)
            pixel_full = np.array(np.broadcast_arrays(*pixel_full)).transpose()

        pixel_full += 1

        if pixel_full.shape[0] == 0:
            world = np.zeros((0, 2))
        else:
            world = self.wcs.wcs_pix2world(pixel_full, 1)

        # At the moment, one has to manually check that the transformation
        # round-trips, otherwise it should be considered invalid.
        pixel_check = self.wcs.wcs_world2pix(world, 1)
        with np.errstate(invalid='ignore'):
            invalid = np.any(np.abs(pixel_check - pixel_full) > 1., axis=1)
        world[invalid] = np.nan

        return world

    transform_non_affine = transform

    def inverted(self):
        return WCSWorld2PixelTransform(self.wcs, slice=self.slice)
