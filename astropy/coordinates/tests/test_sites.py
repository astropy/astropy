import pytest

from astropy import units as u
from astropy.coordinates import EarthLocation, Latitude, Longitude
from astropy.coordinates.sites import (
    SiteRegistry,
    get_builtin_sites,
    get_downloaded_sites,
)
from astropy.tests.helper import assert_quantity_allclose
from astropy.units import allclose as quantity_allclose


def test_builtin_sites():
    reg = get_builtin_sites()

    greenwich = reg["greenwich"]
    lon, lat, el = greenwich.to_geodetic()
    assert_quantity_allclose(lon, Longitude("0:0:0", unit=u.deg), atol=10 * u.arcsec)
    assert_quantity_allclose(lat, Latitude("51:28:40", unit=u.deg), atol=1 * u.arcsec)
    assert_quantity_allclose(el, 46 * u.m, atol=1 * u.m)

    names = reg.names
    assert "greenwich" in names
    assert "example_site" in names

    with pytest.raises(
        KeyError,
        match="Site 'nonexistent' not in database. Use the 'names' attribute to see",
    ):
        reg["nonexistent"]


@pytest.mark.remote_data(source="astropy")
def test_online_sites():
    reg = get_downloaded_sites()

    keck = reg["keck"]
    lon, lat, el = keck.to_geodetic()
    assert_quantity_allclose(
        lon, -Longitude("155:28.7", unit=u.deg), atol=0.001 * u.deg
    )
    assert_quantity_allclose(lat, Latitude("19:49.7", unit=u.deg), atol=0.001 * u.deg)
    assert_quantity_allclose(el, 4160 * u.m, atol=1 * u.m)

    names = reg.names
    assert "keck" in names
    assert "ctio" in names

    # The JSON file contains `name` and `aliases` for each site, and astropy
    # should use names from both, but not empty strings [#12721].
    assert "" not in names
    assert "Royal Observatory Greenwich" in names

    with pytest.raises(
        KeyError,
        match="Site 'nonexistent' not in database. Use the 'names' attribute to see",
    ):
        reg["nonexistent"]

    with pytest.raises(
        KeyError,
        match="Site 'kec' not in database. Use the 'names' attribute to see available",
    ):
        reg["kec"]


@pytest.mark.remote_data(source="astropy")
# this will *try* the online so we have to make it remote_data, even though it
# could fall back on the non-remote version
def test_EarthLocation_basic():
    greenwichel = EarthLocation.of_site("greenwich")
    lon, lat, el = greenwichel.to_geodetic()
    assert_quantity_allclose(lon, Longitude("0:0:0", unit=u.deg), atol=10 * u.arcsec)
    assert_quantity_allclose(lat, Latitude("51:28:40", unit=u.deg), atol=1 * u.arcsec)
    assert_quantity_allclose(el, 46 * u.m, atol=1 * u.m)

    names = EarthLocation.get_site_names()
    assert "greenwich" in names
    assert "example_site" in names

    with pytest.raises(
        KeyError,
        match="Site 'nonexistent' not in database. Use EarthLocation.get_site_names",
    ):
        EarthLocation.of_site("nonexistent")


def test_EarthLocation_state_offline():
    EarthLocation._site_registry = None
    EarthLocation._get_site_registry(force_builtin=True)
    assert EarthLocation._site_registry is not None

    oldreg = EarthLocation._site_registry
    newreg = EarthLocation._get_site_registry()
    assert oldreg is newreg
    newreg = EarthLocation._get_site_registry(force_builtin=True)
    assert oldreg is not newreg


@pytest.mark.remote_data(source="astropy")
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

    assert len(reg.names) == 0

    names = ["sitea", "site A"]
    loc = EarthLocation.from_geodetic(lat=1 * u.deg, lon=2 * u.deg, height=3 * u.km)
    reg.add_site(names, loc)

    assert len(reg.names) == 2

    loc1 = reg["SIteA"]
    assert loc1 is loc

    loc2 = reg["sIte a"]
    assert loc2 is loc


def test_non_EarthLocation():
    """
    A regression test for a typo bug pointed out at the bottom of
    https://github.com/astropy/astropy/pull/4042
    """

    class EarthLocation2(EarthLocation):
        pass

    # This lets keeps us from needing to do remote_data
    # note that this does *not* mess up the registry for EarthLocation because
    # registry is cached on a per-class basis
    EarthLocation2._get_site_registry(force_builtin=True)

    el2 = EarthLocation2.of_site("greenwich")
    assert type(el2) is EarthLocation2
    assert el2.info.name == "Royal Observatory Greenwich"


def check_builtin_matches_remote(download_url=True):
    """
    This function checks that the builtin sites registry is consistent with the
    remote registry (or a registry at some other location).

    Note that current this is *not* run by the testing suite (because it
    doesn't start with "test", and is instead meant to be used as a check
    before merging changes in astropy-data)
    """
    builtin_registry = EarthLocation._get_site_registry(force_builtin=True)
    dl_registry = EarthLocation._get_site_registry(force_download=download_url)

    in_dl = {}
    matches = {}
    for name in builtin_registry.names:
        in_dl[name] = name in dl_registry
        if in_dl[name]:
            matches[name] = quantity_allclose(
                builtin_registry[name].geocentric, dl_registry[name].geocentric
            )
        else:
            matches[name] = False

    if not all(matches.values()):
        # this makes sure we actually see which don't match
        print("In builtin registry but not in download:")
        for name in in_dl:
            if not in_dl[name]:
                print("    ", name)
        print("In both but not the same value:")
        for name in matches:
            if not matches[name] and in_dl[name]:
                print(
                    "    ",
                    name,
                    "builtin:",
                    builtin_registry[name],
                    "download:",
                    dl_registry[name],
                )
        assert False, (
            "Builtin and download registry aren't consistent - failures printed to"
            " stdout"
        )


def test_meta_present():
    reg = get_builtin_sites()

    greenwich = reg["greenwich"]
    assert (
        greenwich.info.meta["source"]
        == "Ordnance Survey via http://gpsinformation.net/main/greenwich.htm and UNESCO"
    )
