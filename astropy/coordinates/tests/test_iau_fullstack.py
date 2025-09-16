# Licensed under a 3-clause BSD style license - see LICENSE.rst

import erfa
import numpy as np
import pytest
from numpy import testing as npt

from astropy import units as u
from astropy.coordinates import EarthLocation, SkyCoord, golden_spiral_grid
from astropy.coordinates.builtin_frames import ICRS, AltAz
from astropy.coordinates.builtin_frames.utils import get_jd12
from astropy.tests.helper import assert_quantity_allclose
from astropy.time import Time
from astropy.utils import iers

ALTAZFRAME = AltAz(
    location=EarthLocation(lat=0 * u.deg, lon=0 * u.deg, height=0 * u.m),
    obstime=Time("J2000"),
)
DEFAULT_OBSCONDITIONS = {
    "pressure": 1 * u.bar,
    "temperature": 0 * u.deg_C,
    "relative_humidity": 0,
    "obswl": 1 * u.μm,
}
FULLSTACK_ICRS = ICRS(golden_spiral_grid(size=1000))


@pytest.mark.parametrize("fullstack_times", [Time("J2000.1"), Time("J2010")])
@pytest.mark.parametrize(
    "fullstack_locations",
    [
        EarthLocation(lat=lat * u.deg, lon=lon * u.deg, height=height * u.m)
        for lat, lon, height in [
            (0, 0, 0),
            (23, 0, 0),
            (-70, 0, 0),
            (0, 100, 0),
            (23, 0, 3000),
        ]
    ],
)
@pytest.mark.parametrize(
    "fullstack_obsconditions,min_alt,tol",
    [
        ({"pressure": 0 * u.bar}, -90 * u.deg, 5 * u.µas),
        ({"relative_humidity": 0 * u.one}, 10 * u.deg, 100 * u.mas),
        ({"temperature": 10 * u.deg_C}, 5 * u.deg, 750 * u.mas),
        ({"relative_humidity": 50 * u.percent}, 5 * u.deg, 750 * u.mas),
        ({"obswl": 21 * u.cm}, 5 * u.deg, 750 * u.mas),
    ],
)
def test_iau_fullstack(
    fullstack_times, fullstack_locations, fullstack_obsconditions, min_alt, tol
):
    """
    Test the full transform from ICRS <-> AltAz
    """
    obsconditions = DEFAULT_OBSCONDITIONS | fullstack_obsconditions

    # Transform to the altaz frame
    aacoo = FULLSTACK_ICRS.transform_to(
        AltAz(obstime=fullstack_times, location=fullstack_locations, **obsconditions)
    )

    # compare aacoo to the fiducial AltAz - should always be different
    fullstack_fiducial_altaz = FULLSTACK_ICRS.transform_to(ALTAZFRAME)
    assert np.all(np.abs(aacoo.alt - fullstack_fiducial_altaz.alt) > 50 * u.mas)
    assert np.all(np.abs(aacoo.az - fullstack_fiducial_altaz.az) > 50 * u.mas)

    # if the refraction correction is included, we *only* do the comparisons
    # where altitude is high enough.  The SOFA guides imply that below 5 deg is
    # where accuracy gets more problematic, and testing reveals that alt<~0
    # gives garbage round-tripping, and <10 can give ~1 arcsec uncertainty,
    # but if there is no refraction correction, we still check everything
    msk = aacoo.alt > min_alt

    # now make sure the full stack round-tripping works
    icrs2 = aacoo.transform_to(ICRS())

    assert_quantity_allclose(
        np.abs(FULLSTACK_ICRS.ra - icrs2.ra)[msk], 0 * u.μas, atol=tol, rtol=0
    )
    assert_quantity_allclose(
        np.abs(FULLSTACK_ICRS.dec - icrs2.dec)[msk], 0 * u.μas, atol=tol, rtol=0
    )

    # check that we're consistent with the ERFA alt/az result
    astrom, _ = erfa.apco13(
        *get_jd12(fullstack_times, "utc"),
        (iers_table := iers.earth_orientation_table.get()).ut1_utc(fullstack_times),
        (geodetic := fullstack_locations.geodetic).lon,
        geodetic.lat,
        geodetic.height,
        *iers_table.pm_xy(fullstack_times),
        obsconditions["pressure"],
        obsconditions["temperature"],
        obsconditions["relative_humidity"],
        obsconditions["obswl"],
    )
    erfa_az, erfa_zen, _, _, _ = erfa.atioq(
        *erfa.atciqz(FULLSTACK_ICRS.ra, FULLSTACK_ICRS.dec, astrom), astrom
    )
    assert_quantity_allclose(aacoo.alt, 90 * u.deg - erfa_zen, atol=1 * u.mas)
    assert_quantity_allclose(aacoo.az, erfa_az, atol=1 * u.mas)


def test_fiducial_roudtrip():
    """
    Test the full transform from ICRS <-> AltAz
    """
    aacoo = FULLSTACK_ICRS.transform_to(ALTAZFRAME)

    # make sure the round-tripping works
    icrs2 = aacoo.transform_to(ICRS())
    npt.assert_allclose(FULLSTACK_ICRS.ra.deg, icrs2.ra.deg)
    npt.assert_allclose(FULLSTACK_ICRS.dec.deg, icrs2.dec.deg)


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

    # check that these message(s) appear among any other warnings
    with (
        pytest.warns(erfa.core.ErfaWarning),
        pytest.warns(
            AstropyWarning,
            match="Tried to get polar motions for times after IERS data is valid.*",
        ),
        iers.conf.set_temp("auto_max_age", None),
    ):
        SkyCoord(1 * u.deg, 2 * u.deg).transform_to(AltAz(location=location, obstime=t))
