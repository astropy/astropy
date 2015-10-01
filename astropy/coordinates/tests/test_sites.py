from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from ...tests.helper import pytest, assert_quantity_allclose, remote_data
from ... import units as u
from .. import Longitude, Latitude, EarthLocation
from ..sites import get_builtin_sites, get_downloaded_sites, SiteRegistry

def test_builtin_sites():
    reg = get_builtin_sites()

    keck = reg.get_site('keck')
    lon, lat, el = keck.to_geodetic()
    assert_quantity_allclose(lon, -1*Longitude('155:28.7', unit=u.deg),
                             atol=0.001*u.deg)
    assert_quantity_allclose(lat, Latitude('19:49.7', unit=u.deg),
                             atol=0.001*u.deg)
    assert_quantity_allclose(el, 4160*u.m, atol=1*u.m)

    keck = reg.get_site('ctio')
    lon, lat, el = keck.to_geodetic()
    assert_quantity_allclose(lon, -1*Longitude('70.815', unit=u.deg),
                             atol=0.001*u.deg)
    assert_quantity_allclose(lat, Latitude('-30.16527778', unit=u.deg),
                             atol=0.001*u.deg)
    assert_quantity_allclose(el, 2215*u.m, atol=1*u.m)

    names = reg.get_site_names()
    assert 'keck' in names
    assert 'ctio' in names

    with pytest.raises(KeyError):
        reg.get_site('nonexistent site')

@pytest.mark.xfail  # remove this when the data file gets uploaded
@remote_data
def test_online_stes():
    reg = get_downloaded_sites()

    keck = reg.get_site('keck')
    lon, lat, el = keck.to_geodetic()
    assert_quantity_allclose(lon, -1*Longitude('155:28.7', unit=u.deg),
                             atol=0.001*u.deg)
    assert_quantity_allclose(lat, Latitude('19:49.7', unit=u.deg),
                             atol=0.001*u.deg)
    assert_quantity_allclose(el, 4160*u.m, atol=1*u.m)

    names = reg.get_site_names()
    assert 'keck' in names
    assert 'ctio' in names

    with pytest.raises(KeyError):
        reg.get_site('nonexistent site')


@remote_data
# this will *try* the online so we have to make it remote_data, even though it
# falls back on the non-remote version
def test_EarthLocation_basic():
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

    with pytest.raises(KeyError):
        EarthLocation.of_site('nonexistent site')


def test_EarthLocation_state_offline():
    EarthLocation._site_registry = None
    EarthLocation._get_site_registry(force_builtin=True)
    assert EarthLocation._site_registry is not None

    oldreg = EarthLocation._site_registry
    newreg = EarthLocation._get_site_registry()
    assert oldreg is newreg
    newreg = EarthLocation._get_site_registry(force_builtin=True)
    assert oldreg is not newreg


@pytest.mark.xfail  # remove this when the data file gets uploaded
@remote_data
def test_EarthLocation_state_online():
    EarthLocation._site_registry = None
    EarthLocation._get_site_registry(force_download=True)
    assert EarthLocation._site_registry is not None

    oldreg = EarthLocation._site_registry
    newreg = EarthLocation._get_site_registry()
    assert oldreg is newreg
    newreg = EarthLocation._get_site_registry(force_download=True)
    assert oldreg is not newreg


def test_registry():
    reg = SiteRegistry()

    assert len(reg.get_site_names()) == 0

    names = ['sitea', 'site A']
    loc = EarthLocation.from_geodetic(lat=1*u.deg, lon=2*u.deg,height=3*u.km)
    reg.add_site(names, loc)

    assert len(reg.get_site_names()) == 2

    loc1 = reg.get_site('SIteA')
    assert loc1 is loc

    loc2 = reg.get_site('sIte a')
    assert loc2 is loc
