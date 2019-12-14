import pytest
import numpy as np

from astropy.coordinates import SkyCoord
from astropy.units import Quantity
from astropy.wcs import WCS
from astropy.wcs.wcsapi import BaseLowLevelWCS


@pytest.fixture
def spectral_1d_fitswcs():
    wcs = WCS(naxis=1)
    wcs.wcs.ctype = 'FREQ',
    wcs.wcs.cunit = 'Hz',
    wcs.wcs.cdelt = 3.e9,
    wcs.wcs.crval = 4.e9,
    wcs.wcs.crpix = 11.,
    wcs.wcs.cname == 'Frequency',
    return wcs


@pytest.fixture
def time_1d_fitswcs():
    wcs = WCS(naxis=1)
    wcs.wcs.ctype = 'TIME',
    wcs.wcs.mjdref = 30042,
    wcs.wcs.crval = 3.,
    wcs.wcs.crpix = 11.,
    wcs.wcs.cname == 'Frequency',
    return wcs


@pytest.fixture
def celestial_2d_fitswcs():
    wcs = WCS(naxis=2)
    wcs.wcs.ctype = 'RA---CAR', 'DEC--CAR'
    wcs.wcs.cunit = 'deg', 'deg'
    wcs.wcs.cdelt = -2., 2.
    wcs.wcs.crval = 4., 0.
    wcs.wcs.crpix = 6., 7.
    wcs.wcs.cname == 'Right Ascension', 'Declination'
    return wcs


class Spectral1DLowLevelWCS(BaseLowLevelWCS):

    @property
    def pixel_n_dim(self):
        return 1

    @property
    def world_n_dim(self):
        return 1

    @property
    def world_axis_physical_types(self):
        return ['em.freq']

    @property
    def world_axis_units(self):
        return ['Hz']

    def pixel_to_world_values(self, pixel_array):
        return np.asarray(pixel_array - 10) * 3e9 + 4e9

    def world_to_pixel_values(self, world_array):
        return np.asarray(world_array - 4e9) / 3e9 + 10

    @property
    def world_axis_object_components(self):
        return [('test', 0, 'value')]

    @property
    def world_axis_object_classes(self):
        return {'test': (Quantity, (), {'unit': 'Hz'})}


@pytest.fixture
def spectral_1d_ape14_wcs():
    return Spectral1DLowLevelWCS()


class Celestial2DLowLevelWCS(BaseLowLevelWCS):

    @property
    def pixel_n_dim(self):
        return 2

    @property
    def world_n_dim(self):
        return 2

    @property
    def world_axis_physical_types(self):
        return ['pos.eq.ra', 'pos.eq.dec']

    @property
    def world_axis_units(self):
        return ['deg', 'deg']

    def pixel_to_world_values(self, px, py):
        return (-(np.asarray(px) - 5.) * 2 + 4.,
                (np.asarray(py) - 6.) * 2)

    def world_to_pixel_values(self, wx, wy):
        return (-(np.asarray(wx) - 4.) / 2 + 5.,
                np.asarray(wy) / 2 + 6.)

    @property
    def world_axis_object_components(self):
        return [('test', 0, 'spherical.lon.degree'),
                ('test', 1, 'spherical.lat.degree')]

    @property
    def world_axis_object_classes(self):
        return {'test': (SkyCoord, (), {'unit': 'deg'})}


@pytest.fixture
def celestial_2d_ape14_wcs():
    return Celestial2DLowLevelWCS()
