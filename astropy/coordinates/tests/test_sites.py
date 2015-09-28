from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from ...tests.helper import pytest, assert_quantity_allclose, remote_data
from ... import units as u
from .. import Longitude, Latitude, EarthLocation
from ..sites import get_site, get_site_names

def test_get_site():
    # Compare to the IRAF observatory list available at:
    # http://tdc-www.harvard.edu/iraf/rvsao/bcvcorr/obsdb.html
    keck = get_site('keck', online=False)
    lon, lat, el = keck.to_geodetic()
    assert_quantity_allclose(lon, -1*Longitude('155:28.7', unit=u.deg),
                             atol=0.001*u.deg)
    assert_quantity_allclose(lat, Latitude('19:49.7', unit=u.deg),
                             atol=0.001*u.deg)
    assert_quantity_allclose(el, 4160*u.m, atol=1*u.m)

    keck = get_site('ctio', online=False)
    lon, lat, el = keck.to_geodetic()
    assert_quantity_allclose(lon, -1*Longitude('70.815', unit=u.deg),
                             atol=0.001*u.deg)
    assert_quantity_allclose(lat, Latitude('-30.16527778', unit=u.deg),
                             atol=0.001*u.deg)
    assert_quantity_allclose(el, 2215*u.m, atol=1*u.m)

def test_get_site_names():
    names = get_site_names(show_aliases=True, online=False)
    assert 'keck' in names
    assert 'ctio' in names

@pytest.mark.xfail  # remove this when the data file gets uploaded
@remote_data
def test_get_site_online():
    keck = get_site('keck', online=True)
    lon, lat, el = keck.to_geodetic()
    assert_quantity_allclose(lon, -1*Longitude('155:28.7', unit=u.deg),
                             atol=0.001*u.deg)
    assert_quantity_allclose(lat, Latitude('19:49.7', unit=u.deg),
                             atol=0.001*u.deg)
    assert_quantity_allclose(el, 4160*u.m, atol=1*u.m)

    names = get_site_names(show_aliases=True, online=True)
    assert 'keck' in names
    assert 'ctio' in names


def test_bad_site():
    with pytest.raises(KeyError):
        get_site('nonexistent site', online=False)

@remote_data
# this will *try* the online so we have to make it remote_data, even though it
# falls back on the non-remote version
def test_with_EarthLocation():
    keckel = EarthLocation.of_site('keck')
    lon, lat, el = keckel.to_geodetic()
    assert_quantity_allclose(lon, -1*Longitude('155:28.7', unit=u.deg),
                             atol=0.001*u.deg)
    assert_quantity_allclose(lat, Latitude('19:49.7', unit=u.deg),
                             atol=0.001*u.deg)
    assert_quantity_allclose(el, 4160*u.m, atol=1*u.m)


    names = EarthLocation.get_site_names()
    assert 'keck' in names
    assert 'ctio' in names
