# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np
from numpy import testing as npt

from ... import units as u
from ...time import Time
from ..builtin_frames import ICRS, AltAz, CIRS
from .. import EarthLocation
from .. import Angle, SkyCoord
from ...tests.helper import pytest
from ... import erfa
from ...utils import iers
from .utils import randomly_sample_sphere


#  These fixtures are used in test_iau_fullstack
@pytest.fixture(scope="module")
def fullstack_icrs():
    ra, dec, _ = randomly_sample_sphere(1000)
    return ICRS(ra=ra, dec=dec)

@pytest.fixture(scope="module")
def fullstack_fiducial_altaz(fullstack_icrs):
    altazframe = AltAz(location=EarthLocation(lat=0*u.deg, lon=0*u.deg, height=0*u.m),
                       obstime=Time('J2000'))
    return fullstack_icrs.transform_to(altazframe)

@pytest.fixture(scope="module", params=['J2000.1', 'J2010'])
def fullstack_times(request):
    return Time(request.param)

@pytest.fixture(scope="module", params=[(0, 0, 0), (23, 0, 0), (-70, 0, 0), (0, 100, 0), (23, 0, 3000)])
def fullstack_locations(request):
    return EarthLocation(lat=request.param[0]*u.deg, lon=request.param[0]*u.deg,
                         height=request.param[0]*u.m)

@pytest.fixture(scope="module", params=[(0*u.bar, 0*u.deg_C, 0, 1*u.micron),
                                        (1*u.bar, 0*u.deg_C, 0, 1*u.micron),
                                        (1*u.bar, 10*u.deg_C, 0, 1*u.micron),
                                        (1*u.bar, 0*u.deg_C, .5, 1*u.micron),
                                        (1*u.bar, 0*u.deg_C, 0, 21*u.cm)])
def fullstack_obsconditions(request):
    return request.param

def _erfa_check(ira, idec, astrom):
    """
    This function does the same thing the astropy layer is supposed to do, but
    all in erfa
    """
    cra, cdec = erfa.atciq(ira, idec, 0, 0, 0, 0, astrom)
    az, zen, ha, odec, ora = erfa.atioq(cra, cdec, astrom)
    alt = np.pi/2-zen
    cra2, cdec2 = erfa.atoiq('A', az, zen, astrom)
    ira2, idec2 = erfa.aticq(cra2, cdec2, astrom)

    dct = locals()
    del dct['astrom']
    return dct

def test_iau_fullstack(fullstack_icrs,  fullstack_fiducial_altaz,
                       fullstack_times, fullstack_locations,
                       fullstack_obsconditions):
    """
    Test the full transform from ICRS <-> AltAz
    """

    #create the altaz frame
    altazframe = AltAz(obstime=fullstack_times, location=fullstack_locations,
                       pressure=fullstack_obsconditions[0],
                       temperature=fullstack_obsconditions[1],
                       relative_humidity=fullstack_obsconditions[2],
                       obswl=fullstack_obsconditions[3])

    aacoo = fullstack_icrs.transform_to(altazframe)

    #compare aacoo to the fiducial AltAz - should always be different
    assert np.all(np.abs(aacoo.alt - fullstack_fiducial_altaz.alt) > 50*u.milliarcsecond)
    assert np.all(np.abs(aacoo.az - fullstack_fiducial_altaz.az) > 50*u.milliarcsecond)

    #if the refraction correction is included, we *only* do the comparisons
    # where altitude >5 degrees.  The SOFA guides imply that below 5 is where
    # where accuracy gets more problematic, and testing reveals that alt<~0
    # gives garbage round-tripping, and <10 can give ~1 arcsec uncertainty
    if fullstack_obsconditions[0].value == 0:
        # but if there is no refraction correction, check everything
        msk = slice(None)
        tol = 5*u.microarcsecond
    else:
        msk = aacoo.alt > 5*u.deg
        # most of them aren't this bad, but some of those at low alt are offset
        # this much.  For alt > 10, this is always better than 100 masec
        tol = 750*u.milliarcsecond


    #now make sure the full stack round-tripping works
    icrs2 = aacoo.transform_to(ICRS)

    adras = np.abs(fullstack_icrs.ra - icrs2.ra)[msk]
    addecs = np.abs(fullstack_icrs.dec - icrs2.dec)[msk]
    assert np.all(adras < tol), 'largest RA change is {0} mas, > {1}'.format(np.max(adras.arcsec*1000), tol)
    assert np.all(addecs < tol), 'largest Dec change is {0} mas, > {1}'.format(np.max(addecs.arcsec*1000), tol)

    #check that we're consistent with the ERFA alt/az result
    xp, yp = u.Quantity(iers.IERS.open().pm_xy(fullstack_times.jd1, fullstack_times.jd2)).to(u.radian).value
    lon = fullstack_locations.geodetic[0].to(u.radian).value
    lat = fullstack_locations.geodetic[1].to(u.radian).value
    height = fullstack_locations.geodetic[2].to(u.m).value
    astrom, eo = erfa.apco13(fullstack_times.jd1, fullstack_times.jd2,
                             fullstack_times.delta_ut1_utc,
                             lon, lat, height,
                             xp, yp,
                             fullstack_obsconditions[0].to(u.hPa).value,
                             fullstack_obsconditions[1].to(u.deg_C).value,
                             fullstack_obsconditions[2],
                             fullstack_obsconditions[3].to(u.micron).value)
    erfadct = _erfa_check(fullstack_icrs.ra.rad, fullstack_icrs.dec.rad, astrom)
    npt.assert_allclose(erfadct['alt'], aacoo.alt.radian, atol=1e-7)
    npt.assert_allclose(erfadct['az'], aacoo.az.radian, atol=1e-7)

def test_fiducial_roudtrip(fullstack_icrs, fullstack_fiducial_altaz):
    """
    Test the full transform from ICRS <-> AltAz
    """
    aacoo = fullstack_icrs.transform_to(fullstack_fiducial_altaz)

    #make sure the round-tripping works
    icrs2 = aacoo.transform_to(ICRS)
    npt.assert_allclose(fullstack_icrs.ra.deg, icrs2.ra.deg)
    npt.assert_allclose(fullstack_icrs.dec.deg, icrs2.dec.deg)


def test_future_altaz(recwarn):
    """
    While this does test the full stack, it is mostly meant to check that a
    warning is raised when attempting to get to AltAz in the future (beyond
    IERS tables)
    """
    from ...utils.exceptions import AstropyWarning

    location = EarthLocation(lat=0*u.deg, lon=0*u.deg)
    t = Time('J2030')

    SkyCoord(1*u.deg, 2*u.deg).transform_to(AltAz(location=location, obstime=t))
    w1 = recwarn.pop(AstropyWarning)
    w2 = recwarn.pop(AstropyWarning)
    assert "Tried to get polar motions for times after IERS data is valid." in str(w1.message)
    assert "(some) times are outside of range covered by IERS table." in str(w2.message)


#<--------------- Below here are tests against "known good" examples ---------->
@pytest.mark.xfail
def test_against_hor2eq():
    """Check that Astropy gives consistent results with an IDL hor2eq example.

    See example input and output here:
    http://idlastro.gsfc.nasa.gov/ftp/pro/astro/hor2eq.pro

    Observatory position for `kpno` from here:
    http://idlastro.gsfc.nasa.gov/ftp/pro/astro/observatory.pro
    """
    obstime = Time('J2014.12')
    kpno = EarthLocation(lon=Angle('111d36.0m'), lat=Angle('31d57.8m'), height=2120.*u.m)

    #pressure=0 means no refraction correction, so need to use hor2eq with refract_=0
    aaframe = AltAz(obstime=obstime, location=kpno, pressure=0)
    aacoo = SkyCoord(alt=12.34*u.deg, az=56.78*u.deg, frame=aaframe)
    astropy_fk5 = aacoo.transform_to('fk5')
    hor2eq_fk5 = SkyCoord(39.208254*u.deg, 34.554463*u.deg, frame='fk5')

    distance = astropy_fk5.separation(hor2eq_fk5)
    assert distance < 1 * u.arcsec


@pytest.mark.xfail
def test_against_pyephem():
    """Check that Astropy gives consistent results with one PyEphem example.

    PyEphem: http://rhodesmill.org/pyephem/

    See example input and output here:
    https://gist.github.com/zonca/1672906
    https://github.com/phn/pytpm/issues/2#issuecomment-3698679
    """
    obstime = Time('2011-09-18 08:50:00')
    location = EarthLocation(lon=Angle('-109d24m53.1s'), lat=Angle('33d41m46.0s'),height=30000.*u.m)
    temperature = 15 * u.deg_C
    pressure = 1.010 * u.bar
    # relative_humidity = ?
    # obswl = ?
    aaframe = AltAz(obstime=obstime,location=location,
                    temperature=temperature, pressure=pressure)

    altaz = SkyCoord('-60.7665d 6.8927d', frame=aaframe)
    radec_actual = altaz.transform_to('icrs')
    radec_expected = SkyCoord('196.497518d -4.569323d', frame='icrs')
    distance = radec_actual.separation(radec_expected)

    # TODO: what assertion precision makes sense here?
    assert distance < 1 * u.arcsec


@pytest.mark.xfail
def test_against_jpl_horizons():
    """Check that Astropy gives consistent results with the JPL Horizons example.

    See example input and output here:
    http://ssd.jpl.nasa.gov/?horizons_tutorial
    """
    obstime = Time('1998-07-28 03:00')
    location = EarthLocation(lon=Angle('248.405300d'), lat=Angle('31.9585d'),height=2.06*u.km)
    temperature = 15 * u.deg_C # TODO: correct???
    pressure = 1.010 * u.bar # TODO: correct???
    # relative_humidity = ?
    # obswl = ?
    aaframe = AltAz(obstime=obstime,location=location,
                    temperature=temperature, pressure=pressure)

    altaz = SkyCoord('143.2970d 2.6223d', frame=aaframe)
    radec_actual = altaz.transform_to('icrs')
    radec_expected = SkyCoord('19h24m55.01s -40d56m28.9s', frame='icrs')
    distance = radec_actual.separation(radec_expected)

    # TODO: what assertion precision makes sense here?
    assert distance < 1 * u.arcsec
