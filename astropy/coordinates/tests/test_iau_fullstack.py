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
from ...tests.helper import pytest
from ... import erfa
from ...utils import iers, NumpyRNGContext


#  These fixtures are used in test_iau_fullstack
@pytest.fixture(scope="module")
def fullstack_icrs():
    N_TO_TEST = 500
    with NumpyRNGContext(12345):
        icrs = ICRS(ra=np.random.rand(N_TO_TEST)*360*u.deg,
                    dec=(np.random.rand(N_TO_TEST)*180-90)*u.deg)
    return icrs

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
    cra, cdec = erfa.atciq(ira, idec, 0, 0, 0, 0, astrom)
    az, zen, ha, odec, ora = erfa.atioq(cra, cdec, astrom)
    cra2, cdec2 = erfa.atoiq('A', az, zen, astrom)
    ira2, idec2 = erfa.aticq(cra2, cdec2, astrom)

    npt.assert_allclose(cra2, cra)
    npt.assert_allclose(cdec2, cdec)

    npt.assert_allclose(ira2, ira)
    npt.assert_allclose(idec2, idec)

    return cra, cdec, az, zen, np.pi/2-zen, cra2, cdec2, ira2, idec2



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
    assert not np.allclose(aacoo.alt.deg, fullstack_fiducial_altaz.alt.deg)
    assert not np.allclose(aacoo.az.deg, fullstack_fiducial_altaz.az.deg)

    #check that we're consistent with the ERFA result
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
    erfares = _erfa_check(fullstack_icrs.ra.rad, fullstack_icrs.dec.rad, astrom)


    #now make sure the round-tripping works

    #first check just the CIRS stages
    cirs1 = fullstack_icrs.transform_to(CIRS(obstime=fullstack_times))
    cirs2 = cirs1.transform_to(altazframe).transform_to(CIRS(obstime=fullstack_times))
    npt.assert_allclose(cirs1.ra.deg, cirs2.ra.deg)
    npt.assert_allclose(cirs1.dec.deg, cirs2.dec.deg)

    #now the full stack
    icrs2 = aacoo.transform_to(ICRS)
    icrs2 = cirs2.transform_to(ICRS)

    npt.assert_allclose(fullstack_icrs.ra.deg, icrs2.ra.deg)
    npt.assert_allclose(fullstack_icrs.dec.deg, icrs2.dec.deg)

def test_fiducial_roudtrip(fullstack_icrs, fullstack_fiducial_altaz):
    """
    Test the full transform from ICRS <-> AltAz
    """
    aacoo = fullstack_icrs.transform_to(fullstack_fiducial_altaz)

    #make sure the round-tripping works
    icrs2 = aacoo.transform_to(ICRS)
    npt.assert_allclose(fullstack_icrs.ra.deg, icrs2.ra.deg)
    npt.assert_allclose(fullstack_icrs.dec.deg, icrs2.dec.deg)

