import shutil
from pathlib import Path
from threading import Lock

import pytest

from astropy import units as u
from astropy.config import set_temp_cache
from astropy.coordinates import EarthLocation, Latitude, Longitude, UnknownSiteException
from astropy.coordinates.sites import SiteRegistry
from astropy.tests.helper import assert_quantity_allclose

SITE_DATA_LOCK = Lock()
TEST_SITE_DATA_DIR = Path(__file__).with_name("data") / "sites"


@pytest.fixture(scope="function")
def local_site_data(monkeypatch):
    with SITE_DATA_LOCK:
        monkeypatch.setattr(EarthLocation, "_site_registry", None)
        with set_temp_cache(TEST_SITE_DATA_DIR):
            yield


@pytest.fixture(scope="function")
def remote_site_data(monkeypatch, tmp_path):
    with SITE_DATA_LOCK:
        monkeypatch.setattr(EarthLocation, "_site_registry", None)
        with set_temp_cache(tmp_path):
            yield


@pytest.fixture(scope="function")
def refreshable_site_data(monkeypatch, tmp_path):
    with SITE_DATA_LOCK:
        monkeypatch.setattr(EarthLocation, "_site_registry", None)
        shutil.copytree(TEST_SITE_DATA_DIR, tmp_path, dirs_exist_ok=True)
        with set_temp_cache(tmp_path):
            yield


@pytest.mark.usefixtures("local_site_data")
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

    with pytest.raises(TypeError, match="^site name None is not a 'str'$"):
        EarthLocation.of_site(None)


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


@pytest.mark.remote_data(source="astropy")
@pytest.mark.usefixtures("remote_site_data")
def test_Earthlocation_get_site_names_empty_cache():
    assert len(EarthLocation.get_site_names()) > 3


@pytest.mark.usefixtures("local_site_data")
def test_Earthlocation_get_site_names_no_refresh_cache():
    assert EarthLocation.get_site_names() == [
        "Royal Observatory Greenwich",
        "example_site",
        "greenwich",
    ]


@pytest.mark.remote_data(source="astropy")
@pytest.mark.usefixtures("refreshable_site_data")
def test_Earthlocation_get_site_names_refresh_cache():
    assert len(EarthLocation.get_site_names()) == 3
    # Regression test for #18572 - refresh_cache=True had no effect
    assert len(EarthLocation.get_site_names(refresh_cache=True)) > 3


@pytest.mark.remote_data(source="astropy")
@pytest.mark.usefixtures("remote_site_data")
def test_Earthlocation_of_site_empty_cache():
    keck = EarthLocation.of_site("keck")
    assert_quantity_allclose(keck.lon, -155.47833333 * u.deg)


@pytest.mark.usefixtures("local_site_data")
def test_Earthlocation_of_site_no_refresh_cache():
    with pytest.raises(UnknownSiteException, match="^Site 'keck' not in database"):
        EarthLocation.of_site("keck")


@pytest.mark.remote_data(source="astropy")
@pytest.mark.usefixtures("refreshable_site_data")
def test_Earthlocation_of_site_refresh_cache():
    with pytest.raises(UnknownSiteException, match="^Site 'keck' not in database"):
        EarthLocation.of_site("keck")
    # Regression test for #18572 - refresh_cache=True had no effect
    keck = EarthLocation.of_site("keck", refresh_cache=True)
    assert_quantity_allclose(keck.lon, -155.47833333 * u.deg)


def test_registry():
    reg = SiteRegistry()
    assert len(reg.names) == 0

    loc = EarthLocation.from_geodetic(lat=1 * u.deg, lon=2 * u.deg, height=3 * u.km)
    reg.add_site(["sitea", "site A"], loc)
    assert len(reg.names) == 2
    assert reg["SIteA"] is loc
    assert reg["sIte a"] is loc


@pytest.mark.usefixtures("local_site_data")
def test_non_EarthLocation():
    """
    A regression test for a typo bug pointed out at the bottom of
    https://github.com/astropy/astropy/pull/4042
    """

    class EarthLocation2(EarthLocation):
        pass

    el2 = EarthLocation2.of_site("greenwich")
    assert type(el2) is EarthLocation2
    assert el2.info.name == "Royal Observatory Greenwich"


@pytest.mark.usefixtures("local_site_data")
def test_meta_present():
    assert (
        EarthLocation.of_site("greenwich").info.meta["source"]
        == "Ordnance Survey via http://gpsinformation.net/main/greenwich.htm and UNESCO"
    )
