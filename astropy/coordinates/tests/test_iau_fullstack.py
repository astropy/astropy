# Licensed under a 3-clause BSD style license - see LICENSE.rst

import warnings

import erfa
import numpy as np
import pytest
from numpy import testing as npt

from astropy import units as u
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.coordinates.angle_utilities import golden_spiral_grid
from astropy.coordinates.builtin_frames import ICRS, AltAz
from astropy.coordinates.builtin_frames.utils import get_jd12
from astropy.time import Time
from astropy.utils import iers


# These fixtures are used in test_iau_fullstack
@pytest.fixture(scope="function")
def fullstack_icrs():
    rep = golden_spiral_grid(size=1000)
    return ICRS(rep)


@pytest.fixture(scope="function")
def fullstack_fiducial_altaz(fullstack_icrs):
    altazframe = AltAz(
        location=EarthLocation(lat=0 * u.deg, lon=0 * u.deg, height=0 * u.m),
        obstime=Time("J2000"),
    )
    with warnings.catch_warnings():  # Ignore remote_data warning
        warnings.simplefilter("ignore")
        result = fullstack_icrs.transform_to(altazframe)
    return result


@pytest.fixture(scope="function", params=["J2000.1", "J2010"])
def fullstack_times(request):
    return Time(request.param)


@pytest.fixture(
    scope="function",
    params=[(0, 0, 0), (23, 0, 0), (-70, 0, 0), (0, 100, 0), (23, 0, 3000)],
)
def fullstack_locations(request):
    value = request.param[0]
    return EarthLocation(lat=value * u.deg, lon=value * u.deg, height=value * u.m)


@pytest.fixture(
    scope="function",
    params=[
        (0 * u.bar, 0 * u.deg_C, 0, 1 * u.micron),
        (1 * u.bar, 0 * u.deg_C, 0 * u.one, 1 * u.micron),
        (1 * u.bar, 10 * u.deg_C, 0, 1 * u.micron),
        (1 * u.bar, 0 * u.deg_C, 50 * u.percent, 1 * u.micron),
        (1 * u.bar, 0 * u.deg_C, 0, 21 * u.cm),
    ],
)
def fullstack_obsconditions(request):
    return request.param


def _erfa_check(ira, idec, astrom):
    """
    This function does the same thing the astropy layer is supposed to do, but
    all in erfa
    """
    cra, cdec = erfa.atciq(ira, idec, 0, 0, 0, 0, astrom)
    az, zen, ha, odec, ora = erfa.atioq(cra, cdec, astrom)
    alt = np.pi / 2 - zen
    cra2, cdec2 = erfa.atoiq("A", az, zen, astrom)
    ira2, idec2 = erfa.aticq(cra2, cdec2, astrom)

    dct = locals()
    del dct["astrom"]
    return dct


def test_iau_fullstack(
    fullstack_icrs,
    fullstack_fiducial_altaz,
    fullstack_times,
    fullstack_locations,
    fullstack_obsconditions,
):
    """
    Test the full transform from ICRS <-> AltAz
    """

    # create the altaz frame
    altazframe = AltAz(
        obstime=fullstack_times,
        location=fullstack_locations,
        pressure=fullstack_obsconditions[0],
        temperature=fullstack_obsconditions[1],
        relative_humidity=fullstack_obsconditions[2],
        obswl=fullstack_obsconditions[3],
    )

    aacoo = fullstack_icrs.transform_to(altazframe)

    # compare aacoo to the fiducial AltAz - should always be different
    assert np.all(
        np.abs(aacoo.alt - fullstack_fiducial_altaz.alt) > 50 * u.milliarcsecond
    )
    assert np.all(
        np.abs(aacoo.az - fullstack_fiducial_altaz.az) > 50 * u.milliarcsecond
    )

    # if the refraction correction is included, we *only* do the comparisons
    # where altitude >5 degrees.  The SOFA guides imply that below 5 is where
    # where accuracy gets more problematic, and testing reveals that alt<~0
    # gives garbage round-tripping, and <10 can give ~1 arcsec uncertainty
    if fullstack_obsconditions[0].value == 0:
        # but if there is no refraction correction, check everything
        msk = slice(None)
        tol = 5 * u.microarcsecond
    else:
        msk = aacoo.alt > 5 * u.deg
        # most of them aren't this bad, but some of those at low alt are offset
        # this much.  For alt > 10, this is always better than 100 masec
        tol = 750 * u.milliarcsecond

    # now make sure the full stack round-tripping works
    icrs2 = aacoo.transform_to(ICRS())

    adras = np.abs(fullstack_icrs.ra - icrs2.ra)[msk]
    addecs = np.abs(fullstack_icrs.dec - icrs2.dec)[msk]
    assert np.all(
        adras < tol
    ), f"largest RA change is {np.max(adras.arcsec * 1000)} mas, > {tol}"
    assert np.all(
        addecs < tol
    ), f"largest Dec change is {np.max(addecs.arcsec * 1000)} mas, > {tol}"

    # check that we're consistent with the ERFA alt/az result
    iers_tab = iers.earth_orientation_table.get()
    xp, yp = u.Quantity(iers_tab.pm_xy(fullstack_times)).to_value(u.radian)
    lon = fullstack_locations.geodetic[0].to_value(u.radian)
    lat = fullstack_locations.geodetic[1].to_value(u.radian)
    height = fullstack_locations.geodetic[2].to_value(u.m)
    jd1, jd2 = get_jd12(fullstack_times, "utc")
    pressure = fullstack_obsconditions[0].to_value(u.hPa)
    temperature = fullstack_obsconditions[1].to_value(u.deg_C)
    # Relative humidity can be a quantity or a number.
    relative_humidity = u.Quantity(fullstack_obsconditions[2], u.one).value
    obswl = fullstack_obsconditions[3].to_value(u.micron)
    astrom, eo = erfa.apco13(
        jd1,
        jd2,
        fullstack_times.delta_ut1_utc,
        lon,
        lat,
        height,
        xp,
        yp,
        pressure,
        temperature,
        relative_humidity,
        obswl,
    )
    erfadct = _erfa_check(fullstack_icrs.ra.rad, fullstack_icrs.dec.rad, astrom)
    npt.assert_allclose(erfadct["alt"], aacoo.alt.radian, atol=1e-7)
    npt.assert_allclose(erfadct["az"], aacoo.az.radian, atol=1e-7)


def test_fiducial_roudtrip(fullstack_icrs, fullstack_fiducial_altaz):
    """
    Test the full transform from ICRS <-> AltAz
    """
    aacoo = fullstack_icrs.transform_to(fullstack_fiducial_altaz)

    # make sure the round-tripping works
    icrs2 = aacoo.transform_to(ICRS())
    npt.assert_allclose(fullstack_icrs.ra.deg, icrs2.ra.deg)
    npt.assert_allclose(fullstack_icrs.dec.deg, icrs2.dec.deg)


def test_future_altaz():
    """
    While this does test the full stack, it is mostly meant to check that a
    warning is raised when attempting to get to AltAz in the future (beyond
    IERS tables)
    """
    # this is an ugly hack to get the warning to show up even if it has already
    # appeared
    from astropy.coordinates.builtin_frames import utils
    from astropy.utils.exceptions import AstropyWarning

    if hasattr(utils, "__warningregistry__"):
        utils.__warningregistry__.clear()

    location = EarthLocation(lat=0 * u.deg, lon=0 * u.deg)
    t = Time("J2161")

    # check that these message(s) appear among any other warnings.  If tests are run with
    # --remote-data then the IERS table will be an instance of IERS_Auto which is
    # assured of being "fresh".  In this case getting times outside the range of the
    # table does not raise an exception.  Only if using IERS_B (which happens without
    # --remote-data, i.e. for all CI testing) do we expect another warning.
    with pytest.warns(
        AstropyWarning,
        match=r"Tried to get polar motions for " "times after IERS data is valid.*",
    ) as found_warnings:
        SkyCoord(1 * u.deg, 2 * u.deg).transform_to(AltAz(location=location, obstime=t))

    if isinstance(iers.earth_orientation_table.get(), iers.IERS_B):
        assert any(
            "(some) times are outside of range covered by IERS table." in str(w.message)
            for w in found_warnings
        )
