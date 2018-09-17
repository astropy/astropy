from numpy.testing import assert_equal, assert_allclose

from .... import units as u
from ....units import Quantity
from ....coordinates import ICRS, Galactic, SkyCoord
from ...fitswcs_low_level_api import FITSLowLevelWCS
from ..base_high_level_api import HighLevelWCS

from ...tests.test_fitswcs_low_level_api import (WCS_SIMPLE_CELESTIAL,
                                                 WCS_SPECTRAL_CUBE)


def test_simple_celestial():

    llwcs = FITSLowLevelWCS(WCS_SIMPLE_CELESTIAL)
    hlwcs = HighLevelWCS(llwcs)

    coord = hlwcs.pixel_to_world(29, 39)
    assert isinstance(coord, SkyCoord)
    assert isinstance(coord.frame, ICRS)
    assert coord.ra.deg == 10
    assert coord.dec.deg == 20

    coord = hlwcs.numpy_index_to_world(39, 29)
    assert isinstance(coord, SkyCoord)
    assert isinstance(coord.frame, ICRS)
    assert coord.ra.deg == 10
    assert coord.dec.deg == 20

    coord = SkyCoord(10, 20, unit='deg', frame='icrs')

    x, y = hlwcs.world_to_pixel(coord)
    assert_allclose(x, 29.)
    assert_allclose(y, 39.)

    i, j = hlwcs.world_to_numpy_index(coord)
    assert_equal(i, 39)
    assert_equal(j, 29)

    # Check that if the coordinates are passed in a different frame things still
    # work properly

    coord_galactic = coord.galactic

    x, y = hlwcs.world_to_pixel(coord_galactic)
    assert_allclose(x, 29.)
    assert_allclose(y, 39.)

    i, j = hlwcs.world_to_numpy_index(coord_galactic)
    assert_equal(i, 39)
    assert_equal(j, 29)


def test_spectral_cube():

    llwcs = FITSLowLevelWCS(WCS_SPECTRAL_CUBE)
    hlwcs = HighLevelWCS(llwcs)

    coord, spec = hlwcs.pixel_to_world(29, 39, 44)
    assert isinstance(coord, SkyCoord)
    assert isinstance(coord.frame, Galactic)
    assert coord.l.deg == 25
    assert coord.b.deg == 10
    assert isinstance(spec, Quantity)
    assert spec.to_value(u.Hz) == 20

    coord, spec = hlwcs.numpy_index_to_world(44, 39, 29)
    assert isinstance(coord, SkyCoord)
    assert isinstance(coord.frame, Galactic)
    assert coord.l.deg == 25
    assert coord.b.deg == 10
    assert isinstance(spec, Quantity)
    assert spec.to_value(u.Hz) == 20

    coord = SkyCoord(25, 10, unit='deg', frame='galactic')
    spec = 20 * u.Hz

    x, y, z = hlwcs.world_to_pixel(coord, spec)
    assert_allclose(x, 29.)
    assert_allclose(y, 39.)
    assert_allclose(z, 44.)

    # Order of world coordinates shouldn't matter
    x, y, z = hlwcs.world_to_pixel(spec, coord)
    assert_allclose(x, 29.)
    assert_allclose(y, 39.)
    assert_allclose(z, 44.)

    i, j, k = hlwcs.world_to_numpy_index(coord, spec)
    assert_equal(i, 44)
    assert_equal(j, 39)
    assert_equal(k, 29)

    # Order of world coordinates shouldn't matter
    i, j, k = hlwcs.world_to_numpy_index(spec, coord)
    assert_equal(i, 44)
    assert_equal(j, 39)
    assert_equal(k, 29)
