# Functions/classes for WCSAxes related to APE14 WCSes

import numpy as np

from astropy.coordinates import SkyCoord, ICRS, BaseCoordinateFrame
from astropy import units as u
from astropy.wcs import WCS
from astropy.wcs.wcsapi import SlicedLowLevelWCS

from .frame import RectangularFrame, EllipticalFrame, RectangularFrame1D
from .transforms import CurvedTransform

__all__ = ['transform_coord_meta_from_wcs', 'WCSWorld2PixelTransform',
           'WCSPixel2WorldTransform']

IDENTITY = WCS(naxis=2)
IDENTITY.wcs.ctype = ["X", "Y"]
IDENTITY.wcs.crval = [0., 0.]
IDENTITY.wcs.crpix = [1., 1.]
IDENTITY.wcs.cdelt = [1., 1.]


def transform_coord_meta_from_wcs(wcs, frame_class, slices=None):

    if slices is not None:
        slices = tuple(slices)

    if wcs.pixel_n_dim > 2:
        if slices is None:
            raise ValueError("WCS has more than 2 pixel dimensions, so "
                             "'slices' should be set")
        elif len(slices) != wcs.pixel_n_dim:
            raise ValueError("'slices' should have as many elements as WCS "
                             "has pixel dimensions (should be {})"
                             .format(wcs.pixel_n_dim))

    is_fits_wcs = isinstance(wcs, WCS)

    coord_meta = {}
    coord_meta['name'] = []
    coord_meta['type'] = []
    coord_meta['wrap'] = []
    coord_meta['unit'] = []
    coord_meta['visible'] = []
    coord_meta['format_unit'] = []

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
                coord_type = "longitude"
            elif "pos.helioprojective.lat" in axis_type:
                format_unit = u.arcsec
                coord_type = "latitude"
            elif "pos.heliographic.stonyhurst.lon" in axis_type:
                coord_wrap = 180.
                format_unit = u.arcsec
                coord_type = "longitude"
            elif "pos.heliographic.stonyhurst.lat" in axis_type:
                format_unit = u.arcsec
                coord_type = "latitude"
            elif "pos.heliographic.carrington.lon" in axis_type:
                coord_wrap = 360.
                format_unit = u.arcsec
                coord_type = "longitude"
            elif "pos.heliographic.carrington.lat" in axis_type:
                format_unit = u.arcsec
                coord_type = "latitude"
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

        # For FITS-WCS, for backward-compatibility, we need to make sure that we
        # provide aliases based on CTYPE for the name.
        if is_fits_wcs:
            if isinstance(wcs, WCS):
                alias = wcs.wcs.ctype[idx][:4].replace('-', '').lower()
            elif isinstance(wcs, SlicedLowLevelWCS):
                alias = wcs._wcs.wcs.ctype[wcs._world_keep[idx]][:4].replace('-', '').lower()
            name = (axis_type, alias) if axis_type else alias
        else:
            name = axis_type or ''

        coord_meta['name'].append(name)

    coord_meta['default_axislabel_position'] = [''] * wcs.world_n_dim
    coord_meta['default_ticklabel_position'] = [''] * wcs.world_n_dim
    coord_meta['default_ticks_position'] = [''] * wcs.world_n_dim

    invert_xy = False
    if slices is not None:
        wcs_slice = list(slices)
        wcs_slice[wcs_slice.index("x")] = slice(None)
        if 'y' in slices:
            wcs_slice[wcs_slice.index("y")] = slice(None)
            invert_xy = slices.index('x') > slices.index('y')
        wcs = SlicedLowLevelWCS(wcs, wcs_slice[::-1])
        world_keep = wcs._world_keep
    else:
        world_keep = list(range(wcs.world_n_dim))

    for i in range(len(coord_meta['type'])):
        coord_meta['visible'].append(i in world_keep)

    transform = WCSPixel2WorldTransform(wcs, invert_xy=invert_xy)

    m = wcs.axis_correlation_matrix.copy()
    if invert_xy:
        m = m[:, ::-1]

    if frame_class is RectangularFrame:

        for i, spine_name in enumerate('bltr'):
            pos = np.nonzero(m[:, i % 2])[0]
            if len(pos) > 0:
                index = world_keep[pos[0]]
                coord_meta['default_axislabel_position'][index] = spine_name
                coord_meta['default_ticklabel_position'][index] = spine_name
                coord_meta['default_ticks_position'][index] = spine_name
                m[pos[0], :] = 0

        # In the special and common case where the frame is rectangular and
        # we are dealing with 2-d WCS (after slicing), we show all ticks on
        # all axes for backward-compatibility.
        if len(world_keep) == 2:
            for index in world_keep:
                coord_meta['default_ticks_position'][index] = 'bltr'

    elif frame_class is RectangularFrame1D:

        for i, spine_name in enumerate('bt'):
            pos = np.nonzero(m[:, 0])[0]
            if len(pos) > 0:
                index = world_keep[pos[0]]
                coord_meta['default_axislabel_position'][index] = spine_name
                coord_meta['default_ticklabel_position'][index] = spine_name
                coord_meta['default_ticks_position'][index] = spine_name
                m[pos[0], :] = 0

        # In the special and common case where the frame is rectangular and
        # we are dealing with 2-d WCS (after slicing), we show all ticks on
        # all axes for backward-compatibility.
        if len(world_keep) == 1:
            for index in world_keep:
                coord_meta['default_ticks_position'][index] = 'bt'

    elif frame_class is EllipticalFrame:

        if 'longitude' in coord_meta['type']:
            lon_idx = coord_meta['type'].index('longitude')
            coord_meta['default_axislabel_position'][lon_idx] = 'h'
            coord_meta['default_ticklabel_position'][lon_idx] = 'h'
            coord_meta['default_ticks_position'][lon_idx] = 'h'

        if 'latitude' in coord_meta['type']:
            lat_idx = coord_meta['type'].index('latitude')
            coord_meta['default_axislabel_position'][lat_idx] = 'c'
            coord_meta['default_ticklabel_position'][lat_idx] = 'c'
            coord_meta['default_ticks_position'][lat_idx] = 'c'

    else:

        for i in range(len(coord_meta['type'])):
            if i in world_keep:
                index = world_keep[i]
                coord_meta['default_axislabel_position'][index] = frame_class.spine_names
                coord_meta['default_ticklabel_position'][index] = frame_class.spine_names
                coord_meta['default_ticks_position'][index] = frame_class.spine_names

    return transform, coord_meta


def wcsapi_to_celestial_frame(wcs):
    for cls, args, kwargs in wcs.world_axis_object_classes.values():
        if issubclass(cls, SkyCoord):
            return kwargs.get('frame', ICRS())
        elif issubclass(cls, BaseCoordinateFrame):
            return cls(**kwargs)


class WCSWorld2PixelTransform(CurvedTransform):
    """
    WCS transformation from world to pixel coordinates
    """

    has_inverse = True
    frame_in = None

    def __init__(self, wcs, invert_xy=False):

        super().__init__()

        if wcs.pixel_n_dim > 2:
            raise ValueError('Only pixel_n_dim =< 2 is supported')

        self.wcs = wcs
        self.invert_xy = invert_xy

        self.frame_in = wcsapi_to_celestial_frame(wcs)

    def __eq__(self, other):
        return (isinstance(other, type(self)) and self.wcs is other.wcs and
                self.invert_xy == other.invert_xy)

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

    transform_non_affine = transform

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

        if wcs.pixel_n_dim > 2:
            raise ValueError('Only pixel_n_dim =< 2 is supported')

        self.wcs = wcs
        self.invert_xy = invert_xy

        self.frame_out = wcsapi_to_celestial_frame(wcs)

    def __eq__(self, other):
        return (isinstance(other, type(self)) and self.wcs is other.wcs and
                self.invert_xy == other.invert_xy)

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

        if self.wcs.world_n_dim == 1:
            world = [world]

        # At the moment, one has to manually check that the transformation
        # round-trips, otherwise it should be considered invalid.
        pixel_check = self.wcs.world_to_pixel_values(*world)
        if self.wcs.pixel_n_dim == 1:
            pixel_check = [pixel_check]

        with np.errstate(invalid='ignore'):
            invalid = np.zeros(len(pixel[0]), dtype=bool)
            for ipix in range(len(pixel)):
                invalid |= np.abs(pixel_check[ipix] - pixel[ipix]) > 1.
            for iwrl in range(len(world)):
                world[iwrl][invalid] = np.nan

        world = np.array(world).T

        return world

    transform_non_affine = transform

    def inverted(self):
        """
        Return the inverse of the transform
        """
        return WCSWorld2PixelTransform(self.wcs, invert_xy=self.invert_xy)
