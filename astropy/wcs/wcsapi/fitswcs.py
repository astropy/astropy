# This file includes the definition of a mix-in class that provides the low-
# and high-level WCS API to the astropy.wcs.WCS object. We keep this code
# isolated in this mix-in class to avoid making the main wcs.py file too
# long.

import warnings

import numpy as np

from astropy import units as u

from .low_level_api import BaseLowLevelWCS
from .high_level_api import HighLevelWCSMixin
from .sliced_low_level_wcs import SlicedLowLevelWCS

__all__ = ['custom_ctype_to_ucd_mapping', 'SlicedFITSWCS', 'FITSWCSAPIMixin']

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
    'HGLN': 'custom:pos.heliographic.stonyhurst.lon',
    'HGLT': 'custom:pos.heliographic.stonyhurst.lat',
    'CRLN': 'custom:pos.heliographic.carrington.lon',
    'CRLT': 'custom:pos.heliographic.carrington.lat',

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
    'BETA': 'custom:spect.doplerVeloc.beta',  # Beta factor (v/c)

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

    # UT() and TT() are handled separately in world_axis_physical_types

}

# Keep a list of additional custom mappings that have been registered. This
# is kept as a list in case nested context managers are used
CTYPE_TO_UCD1_CUSTOM = []


class custom_ctype_to_ucd_mapping:
    """
    A context manager that makes it possible to temporarily add new CTYPE to
    UCD1+ mapping used by :attr:`FITSWCSAPIMixin.world_axis_physical_types`.

    Parameters
    ----------
    mapping : dict
        A dictionary mapping a CTYPE value to a UCD1+ value

    Examples
    --------

    Consider a WCS with the following CTYPE::

        >>> from astropy.wcs import WCS
        >>> wcs = WCS(naxis=1)
        >>> wcs.wcs.ctype = ['SPAM']

    By default, :attr:`FITSWCSAPIMixin.world_axis_physical_types` returns `None`,
    but this can be overridden::

        >>> wcs.world_axis_physical_types
        [None]
        >>> with custom_ctype_to_ucd_mapping({'SPAM': 'food.spam'}):
        ...     wcs.world_axis_physical_types
        ['food.spam']
    """

    def __init__(self, mapping):
        CTYPE_TO_UCD1_CUSTOM.insert(0, mapping)
        self.mapping = mapping

    def __enter__(self):
        pass

    def __exit__(self, type, value, tb):
        CTYPE_TO_UCD1_CUSTOM.remove(self.mapping)


class SlicedFITSWCS(SlicedLowLevelWCS, HighLevelWCSMixin):
    pass


class FITSWCSAPIMixin(BaseLowLevelWCS, HighLevelWCSMixin):
    """
    A mix-in class that is intended to be inherited by the
    :class:`~astropy.wcs.WCS` class and provides the low- and high-level WCS API
    """

    @property
    def pixel_n_dim(self):
        return self.naxis

    @property
    def world_n_dim(self):
        return len(self.wcs.ctype)

    @property
    def array_shape(self):
        if self._naxis == [0, 0]:
            return None
        else:
            return tuple(self._naxis[::-1])

    @array_shape.setter
    def array_shape(self, value):
        if value is None:
            self._naxis = [0, 0]
        else:
            if len(value) != self.naxis:
                raise ValueError("The number of data axes, "
                                 "{}, does not equal the "
                                 "shape {}.".format(self.naxis, len(value)))
            self._naxis = list(value)[::-1]

    @property
    def pixel_shape(self):
        if self._naxis == [0, 0]:
            return None
        else:
            return tuple(self._naxis)

    @pixel_shape.setter
    def pixel_shape(self, value):
        if value is None:
            self._naxis = [0, 0]
        else:
            if len(value) != self.naxis:
                raise ValueError("The number of data axes, "
                                 "{}, does not equal the "
                                 "shape {}.".format(self.naxis, len(value)))
            self._naxis = list(value)

    @property
    def pixel_bounds(self):
        return self._pixel_bounds

    @pixel_bounds.setter
    def pixel_bounds(self, value):
        if value is None:
            self._pixel_bounds = value
        else:
            if len(value) != self.naxis:
                raise ValueError("The number of data axes, "
                                 "{}, does not equal the number of "
                                 "pixel bounds {}.".format(self.naxis, len(value)))
            self._pixel_bounds = list(value)

    @property
    def world_axis_physical_types(self):
        types = []
        # TODO: need to support e.g. TT(TAI)
        for ctype in self.wcs.ctype:
            if ctype.startswith(('UT(', 'TT(')):
                types.append('time')
            else:
                ctype_name = ctype.split('-')[0]
                for custom_mapping in CTYPE_TO_UCD1_CUSTOM:
                    if ctype_name in custom_mapping:
                        types.append(custom_mapping[ctype_name])
                        break
                else:
                    types.append(CTYPE_TO_UCD1.get(ctype_name, None))
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
    def world_axis_names(self):
        return list(self.wcs.cname)

    @property
    def axis_correlation_matrix(self):

        # If there are any distortions present, we assume that there may be
        # correlations between all axes. Maybe if some distortions only apply
        # to the image plane we can improve this?
        if self.has_distortion:
            return np.ones((self.world_n_dim, self.pixel_n_dim), dtype=bool)

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
        world = self.all_pix2world(*pixel_arrays, 0)
        return world[0] if self.world_n_dim == 1 else tuple(world)

    def array_index_to_world_values(self, *indices):
        world = self.all_pix2world(*indices[::-1], 0)
        return world[0] if self.world_n_dim == 1 else tuple(world)

    def world_to_pixel_values(self, *world_arrays):
        pixel = self.all_world2pix(*world_arrays, 0)
        return pixel[0] if self.pixel_n_dim == 1 else tuple(pixel)

    def world_to_array_index_values(self, *world_arrays):
        pixel_arrays = self.all_world2pix(*world_arrays, 0)[::-1]
        array_indices = tuple(np.asarray(np.floor(pixel + 0.5), dtype=np.int_) for pixel in pixel_arrays)
        return array_indices[0] if self.pixel_n_dim == 1 else array_indices

    @property
    def world_axis_object_components(self):
        return self._get_components_and_classes()[0]

    @property
    def world_axis_object_classes(self):
        return self._get_components_and_classes()[1]

    @property
    def serialized_classes(self):
        return False

    def _get_components_and_classes(self):

        # The aim of this function is to return whatever is needed for
        # world_axis_object_components and world_axis_object_classes. It's easier
        # to figure it out in one go and then return the values and let the
        # properties return part of it.

        # Since this method might get called quite a few times, we need to cache
        # it. We start off by defining a hash based on the attributes of the
        # WCS that matter here (we can't just use the WCS object as a hash since
        # it is mutable)
        wcs_hash = (self.naxis,
                    list(self.wcs.ctype),
                    list(self.wcs.cunit),
                    self.wcs.radesys,
                    self.wcs.equinox,
                    self.wcs.dateobs,
                    self.wcs.lng,
                    self.wcs.lat)

        # If the cache is present, we need to check that the 'hash' matches.
        if getattr(self, '_components_and_classes_cache', None) is not None:
            cache = self._components_and_classes_cache
            if cache[0] == wcs_hash:
                return cache[1]
            else:
                self._components_and_classes_cache = None

        # Avoid circular imports by importing here
        from astropy.wcs.utils import wcs_to_celestial_frame
        from astropy.coordinates import SkyCoord, EarthLocation
        from astropy.time.formats import FITS_DEPRECATED_SCALES
        from astropy.time import Time, TimeDelta

        components = [None] * self.naxis
        classes = {}

        # Let's start off by checking whether the WCS has a pair of celestial
        # components

        if self.has_celestial:

            try:
                frame = wcs_to_celestial_frame(self)
            except ValueError:
                # Some WCSes, e.g. solar, can be recognized by WCSLIB as being
                # celestial but we don't necessarily have frames for them.
                pass
            else:

                kwargs = {}
                kwargs['frame'] = frame
                kwargs['unit'] = u.deg

                classes['celestial'] = (SkyCoord, (), kwargs)

                components[self.wcs.lng] = ('celestial', 0, 'spherical.lon.degree')
                components[self.wcs.lat] = ('celestial', 1, 'spherical.lat.degree')

        # We can then make sure we correctly return Time objects where appropriate
        # (https://www.aanda.org/articles/aa/pdf/2015/02/aa24653-14.pdf)

        if 'time' in self.world_axis_physical_types:

            multiple_time = self.world_axis_physical_types.count('time') > 1

            for i in range(self.naxis):

                if self.world_axis_physical_types[i] == 'time':

                    if multiple_time:
                        name = f'time.{i}'
                    else:
                        name = 'time'

                    # Initialize delta
                    reference_time_delta = None

                    # Extract time scale
                    scale = self.wcs.ctype[i].lower()

                    if scale == 'time':
                        if self.wcs.timesys:
                            scale = self.wcs.timesys.lower()
                        else:
                            scale = 'utc'

                    # Drop sub-scales
                    if '(' in scale:
                        pos = scale.index('(')
                        scale, subscale = scale[:pos], scale[pos+1:-1]
                        warnings.warn(f'Dropping unsupported sub-scale '
                                      f'{subscale.upper()} from scale {scale.upper()}',
                                      UserWarning)

                    # TODO: consider having GPS as a scale in Time
                    # For now GPS is not a scale, we approximate this by TAI - 19s
                    if scale == 'gps':
                        reference_time_delta = TimeDelta(19, format='sec')
                        scale = 'tai'

                    elif scale.upper() in FITS_DEPRECATED_SCALES:
                        scale = FITS_DEPRECATED_SCALES[scale.upper()]

                    elif scale not in Time.SCALES:
                        raise ValueError(f'Unrecognized time CTYPE={self.wcs.ctype[i]}')

                    # Determine location
                    trefpos = self.wcs.trefpos.lower()

                    if trefpos.startswith('topocent'):
                        # Note that some headers use TOPOCENT instead of TOPOCENTER
                        if np.any(np.isnan(self.wcs.obsgeo[:3])):
                            warnings.warn('Missing or incomplete observer location '
                                          'information, setting location in Time to None',
                                          UserWarning)
                            location = None
                        else:
                            location = EarthLocation(*self.wcs.obsgeo[:3], unit=u.m)
                    elif trefpos == 'geocenter':
                        location = EarthLocation(0, 0, 0, unit=u.m)
                    elif trefpos == '':
                        location = None
                    else:
                        # TODO: implement support for more locations when Time supports it
                        warnings.warn(f"Observation location '{trefpos}' is not "
                                       "supported, setting location in Time to None", UserWarning)
                        location = None

                    reference_time = Time(np.nan_to_num(self.wcs.mjdref[0]),
                                          np.nan_to_num(self.wcs.mjdref[1]),
                                          format='mjd', scale=scale,
                                          location=location)

                    if reference_time_delta is not None:
                        reference_time = reference_time + reference_time_delta

                    def time_from_reference_and_offset(offset):
                        if isinstance(offset, Time):
                            return offset
                        return reference_time + TimeDelta(offset, format='sec')

                    def offset_from_time_and_reference(time):
                        return (time - reference_time).sec

                    classes[name] = (Time, (), {}, time_from_reference_and_offset)
                    components[i] = (name, 0, offset_from_time_and_reference)

        # Fallback: for any remaining components that haven't been identified, just
        # return Quantity as the class to use

        for i in range(self.naxis):
            if components[i] is None:
                name = self.wcs.ctype[i].split('-')[0].lower()
                if name == '':
                    name = 'world'
                while name in classes:
                    name += "_"
                classes[name] = (u.Quantity, (), {'unit': self.wcs.cunit[i]})
                components[i] = (name, 0, 'value')

        # Keep a cached version of result
        self._components_and_classes_cache = wcs_hash, (components, classes)

        return components, classes
