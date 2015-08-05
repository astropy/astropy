from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from ...tests.helper import pytest, assert_quantity_allclose
from ... import units as u
from .. import Latitude, Longitude, EarthLocation, get_site, add_site

def test_get_site():
    # Compare to the IRAF observatory list available at:
    # http://tdc-www.harvard.edu/iraf/rvsao/bcvcorr/obsdb.html
    keck = get_site('keck')
    lon, lat, el = keck.to_geodetic()
    assert_quantity_allclose(lon, -1*Longitude('155:28.7', unit=u.deg),
                             atol=0.001*u.deg)
    assert_quantity_allclose(lat, Latitude('19:49.7', unit=u.deg),
                             atol=0.001*u.deg)
    assert_quantity_allclose(el, 4160*u.m, atol=1*u.m)

    keck = get_site('ctio')
    lon, lat, el = keck.to_geodetic()
    assert_quantity_allclose(lon, -1*Longitude('70.815', unit=u.deg),
                             atol=0.001*u.deg)
    assert_quantity_allclose(lat, Latitude('-30.16527778', unit=u.deg),
                             atol=0.001*u.deg)
    assert_quantity_allclose(el, 2215*u.m, atol=1*u.m)

def test_add_site():
    # Test observatory can be added and retrieved
    new_site_name = 'University of Washington'
    new_site_location = EarthLocation(-122.3080*u.deg, 47.6550*u.deg, 0*u.m)
    add_site(new_site_name, new_site_location)
    retrieved_location = get_site(new_site_name)
    assert retrieved_location == new_site_location

def test_bad_site():
    with pytest.raises(KeyError):
        get_site('nonexistent site')
