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
from astropy.utils.exceptions import AstropyUserWarning


@pytest.fixture
def earthlocation_without_site_registry(monkeypatch):
    monkeypatch.setattr(EarthLocation, "_site_registry", None)


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


@pytest.mark.parametrize(
    "class_method,args",
    [(EarthLocation.get_site_names, []), (EarthLocation.of_site, ["greenwich"])],
)
def test_Earthlocation_refresh_cache_is_mandatory_kwarg(class_method, args):
    with pytest.raises(
        TypeError,
        match=(
            rf".*{class_method.__name__}\(\) takes [12] positional "
            "arguments? but [23] were given$"
        ),
    ):
        class_method(*args, False)


@pytest.mark.parametrize(
    "class_method,args",
    [(EarthLocation.get_site_names, []), (EarthLocation.of_site, ["greenwich"])],
)
@pytest.mark.parametrize("refresh_cache", [False, True])
def test_Earthlocation_refresh_cache(class_method, args, refresh_cache, monkeypatch):
    def get_site_registry_monkeypatched(force_download, force_builtin=False):
        assert force_download is refresh_cache
        return get_builtin_sites()

    monkeypatch.setattr(
        EarthLocation, "_get_site_registry", get_site_registry_monkeypatched
    )

    class_method(*args, refresh_cache=refresh_cache)


@pytest.mark.parametrize(
    "force_download,expectation",
    [
        (
            False,
            pytest.warns(
                AstropyUserWarning, match=r"use the option 'refresh_cache=True'\.$"
            ),
        ),
        (True, pytest.raises(OSError, match=r"^fail for test$")),
        ("url", pytest.raises(OSError, match=r"^fail for test$")),
    ],
)
@pytest.mark.parametrize(
    "class_method,args",
    [(EarthLocation.get_site_names, []), (EarthLocation.of_site, ["greenwich"])],
)
def test_EarthLocation_site_registry_connection_fail(
    force_download,
    expectation,
    class_method,
    args,
    earthlocation_without_site_registry,
    monkeypatch,
):
    def fail_download(*args, **kwargs):
        raise OSError("fail for test")

    monkeypatch.setattr(get_downloaded_sites, "__code__", fail_download.__code__)
    with expectation:
        class_method(*args, refresh_cache=force_download)


@pytest.mark.parametrize(
    "registry_kwarg",
    ["force_builtin", pytest.param("force_download", marks=pytest.mark.remote_data)],
)
def test_EarthLocation_state(earthlocation_without_site_registry, registry_kwarg):
    EarthLocation._get_site_registry(**{registry_kwarg: True})
    assert isinstance(EarthLocation._site_registry, SiteRegistry)

    oldreg = EarthLocation._site_registry
    assert oldreg is EarthLocation._get_site_registry()
    assert oldreg is not EarthLocation._get_site_registry(**{registry_kwarg: True})


def test_registry():
    reg = SiteRegistry()
    assert len(reg.names) == 0

    loc = EarthLocation.from_geodetic(lat=1 * u.deg, lon=2 * u.deg, height=3 * u.km)
    reg.add_site(["sitea", "site A"], loc)
    assert len(reg.names) == 2
    assert reg["SIteA"] is loc
    assert reg["sIte a"] is loc


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
    assert (
        get_builtin_sites()["greenwich"].info.meta["source"]
        == "Ordnance Survey via http://gpsinformation.net/main/greenwich.htm and UNESCO"
    )
