# This file includes the definition of a mix-in class that provides the low-
# and high-level WCS API to the astropy.wcs.WCS object. We keep this code
# isolated in this mix-in class to avoid making the main wcs.py file too
# long.

import numpy as np

from ... import units as u

from .low_level_api import BaseLowLevelWCS
from .high_level_api import HighLevelWCSMixin

__all__ = ['FITSWCSAPIMixin']

# Mapping from CTYPE axis name to UCD1

CTYPE_TO_UCD1 = {

    # Celestial coordinates
    'RA': 'pos.eq.ra',
    'DEC': 'pos.eq.dec',
    'GLON': 'pos.galactic.lon',
    'GLAT': 'pos.galactic.lat',
    'ELON': 'pos.ecliptic.lon',
    'ELAT': 'pos.ecliptic.lat',
    'TLON': 'pos.bodyrc.lon',
    'TLAT': 'pos.bodyrc.lat',
    'HPLT': 'custom:pos.helioprojective.lat',
    'HPLN': 'custom:pos.helioprojective.lon',

    # Spectral coordinates (WCS paper 3)
    'FREQ': 'em.freq',  # Frequency
    'ENER': 'em.energy',  # Energy
    'WAVN': 'em.wavenumber',  # Wavenumber
    'WAVE': 'em.wl',  # Vacuum wavelength
    'VRAD': 'spect.dopplerVeloc.radio',  # Radio velocity
    'VOPT': 'spect.dopplerVeloc.opt',  # Optical velocity
    'ZOPT': 'src.redshift',  # Redshift
    'AWAV': 'em.wl',  # Air wavelength
    'VELO': 'spect.dopplerVeloc',  # Apparent radial velocity
    'BETA': None,  # Beta factor (v/c)

    # Time coordinates (https://www.aanda.org/articles/aa/pdf/2015/02/aa24653-14.pdf)
    'TIME': 'time',
    'TAI': 'time',
    'TT': 'time',
    'TDT': 'time',
    'ET': 'time',
    'IAT': 'time',
    'UT1': 'time',
    'UTC': 'time',
    'GMT': 'time',
    'GPS': 'time',
    'TCG': 'time',
    'TCB': 'time',
    'TDB': 'time',
    'LOCAL': 'time'

    # UT() is handled separately in world_axis_physical_types

}


class FITSWCSAPIMixin(BaseLowLevelWCS, HighLevelWCSMixin):
    """
    A wrapper around the :class:`astropy.wcs.WCS` class that provides the
    low-level WCS API from APE 14.
    """

    @property
    def pixel_n_dim(self):
        return self.naxis

    @property
    def world_n_dim(self):
        return len(self.wcs.ctype)

    @property
    def world_axis_physical_types(self):
        types = []
        for axis_type in self.axis_type_names:
            if axis_type.startswith('UT('):
                types.append('time')
            else:
                types.append(CTYPE_TO_UCD1.get(axis_type, None))
        return types

    @property
    def world_axis_units(self):
        units = []
        for unit in self.wcs.cunit:
            if unit is None:
                unit = ''
            elif isinstance(unit, u.Unit):
                unit = unit.to_string(format='vounit')
            else:
                try:
                    unit = u.Unit(unit).to_string(format='vounit')
                except u.UnitsError:
                    unit = ''
            units.append(unit)
        return units

    @property
    def axis_correlation_matrix(self):

        # If there are any distortions present, we assume that there may be
        # correlations between all axes. Maybe if some distortions only apply
        # to the image plane we can improve this?
        for distortion_attribute in ('sip', 'det2im1', 'det2im2'):
            if getattr(self, distortion_attribute):
                return np.ones((self.n_world, self.n_pixel), dtype=bool)

        # Assuming linear world coordinates along each axis, the correlation
        # matrix would be given by whether or not the PC matrix is zero
        matrix = self.wcs.get_pc() != 0

        # We now need to check specifically for celestial coordinates since
        # these can assume correlations because of spherical distortions. For
        # each celestial coordinate we copy over the pixel dependencies from
        # the other celestial coordinates.
        celestial = (self.wcs.axis_types // 1000) % 10 == 2
        celestial_indices = np.nonzero(celestial)[0]
        for world1 in celestial_indices:
            for world2 in celestial_indices:
                if world1 != world2:
                    matrix[world1] |= matrix[world2]
                    matrix[world2] |= matrix[world1]

        return matrix

    def pixel_to_world_values(self, *pixel_arrays):
        return self.all_pix2world(*pixel_arrays, 0)

    def array_index_to_world_values(self, *indices):
        return self.all_pix2world(*indices[::-1], 0)

    def world_to_pixel_values(self, *world_arrays):
        return self.all_world2pix(*world_arrays, 0)

    def world_to_array_index_values(self, *world_arrays):
        return np.round(self.all_world2pix(*world_arrays, 0)[::-1]).astype(int)

    @property
    def world_axis_object_components(self):
        return _get_components_and_classes(self)[0]

    @property
    def world_axis_object_classes(self):
        return _get_components_and_classes(self)[1]

    @property
    def serialized_classes(self):
        return False


def _get_components_and_classes(wcs):

    # The aim of this function is to return whatever is needed for
    # world_axis_object_components and world_axis_object_classes. It's easier
    # to figure it out in one go and then return the values and let the
    # properties return part of it.

    # Avoid circular imports by importing here
    from ..utils import wcs_to_celestial_frame
    from ...coordinates import SkyCoord

    components = [None] * wcs.naxis
    classes = {}

    # Let's start off by checking whether the WCS has a pair of celestial
    # components

    if wcs.has_celestial:

        frame = wcs_to_celestial_frame(wcs)

        kwargs = {}
        kwargs['frame'] = frame
        kwargs['unit'] = u.deg

        classes['celestial'] = (SkyCoord, (), kwargs)

        components[wcs.wcs.lng] = ('celestial', 0, 'spherical.lon.degree')
        components[wcs.wcs.lat] = ('celestial', 1, 'spherical.lat.degree')

    # Fallback: for any remaining components that haven't been identified, just
    # return Quantity as the class to use

    for i in range(wcs.naxis):
        if components[i] is None:
            name = wcs.axis_type_names[i].lower()
            while name in classes:
                name += "_"
            classes[name] = (u.Quantity, (), {'unit': wcs.wcs.cunit[i]})
            components[i] = (name, 0, 'value')

    return components, classes
