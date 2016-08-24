# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Regression tests for coordinates-related bugs that don't have an obvious other
place to live
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np

from ... import units as u
from .. import (AltAz, EarthLocation, SkyCoord, get_sun, ICRS, GCRS,
                FK4, FK4NoETerms)
from ...time import Time

from ...tests.helper import pytest, assert_quantity_allclose
from .test_matching import HAS_SCIPY, OLDER_SCIPY


def test_regression_3920():
    """
    Issue: https://github.com/astropy/astropy/issues/3920
    """
    loc = EarthLocation.from_geodetic(0*u.deg, 0*u.deg, 0)
    time = Time('2010-1-1')

    aa = AltAz(location=loc, obstime=time)
    sc = SkyCoord(10*u.deg, 3*u.deg)
    assert sc.transform_to(aa).shape == tuple()
    # That part makes sense: the input is a scalar so the output is too

    sc2 = SkyCoord(10*u.deg, 3*u.deg, 1*u.AU)
    assert sc2.transform_to(aa).shape == tuple()
    # in 3920 that assert fails, because the shape is (1,)

    # check that the same behavior occurs even if transform is from low-level classes
    icoo = ICRS(sc.data)
    icoo2 = ICRS(sc2.data)
    assert icoo.transform_to(aa).shape == tuple()
    assert icoo2.transform_to(aa).shape == tuple()


def test_regression_3938():
    """
    Issue: https://github.com/astropy/astropy/issues/3938
    """
    # Set up list of targets - we don't use `from_name` here to avoid
    # remote_data requirements, but it does the same thing
    # vega = SkyCoord.from_name('Vega')
    vega = SkyCoord(279.23473479*u.deg, 38.78368896*u.deg)
    # capella = SkyCoord.from_name('Capella')
    capella = SkyCoord(79.17232794*u.deg, 45.99799147*u.deg)
    # sirius = SkyCoord.from_name('Sirius')
    sirius = SkyCoord(101.28715533*u.deg, -16.71611586*u.deg)
    targets = [vega, capella, sirius]

    # Feed list of targets into SkyCoord
    combined_coords = SkyCoord(targets)

    # Set up AltAz frame
    time = Time('2012-01-01 00:00:00')
    location = EarthLocation('10d', '45d', 0)
    aa = AltAz(location=location, obstime=time)

    combined_coords.transform_to(aa)
    # in 3938 the above yields ``UnitConversionError: '' (dimensionless) and 'pc' (length) are not convertible``


def test_regression_3998():
    """
    Issue: https://github.com/astropy/astropy/issues/3998
    """
    time = Time('2012-01-01 00:00:00')
    assert time.isscalar

    sun = get_sun(time)
    assert sun.isscalar
    # in 3998, the above yields False - `sun` is a length-1 vector

    assert sun.obstime is time


def test_regression_4033():
    """
    Issue: https://github.com/astropy/astropy/issues/4033
    """
    # alb = SkyCoord.from_name('Albireo')
    alb = SkyCoord(292.68033548*u.deg, 27.95968007*u.deg)
    alb_wdist = SkyCoord(alb, distance=133*u.pc)

    # de = SkyCoord.from_name('Deneb')
    de = SkyCoord(310.35797975*u.deg, 45.28033881*u.deg)
    de_wdist = SkyCoord(de, distance=802*u.pc)

    aa = AltAz(location=EarthLocation(lat=45*u.deg, lon=0*u.deg), obstime='2010-1-1')
    deaa = de.transform_to(aa)
    albaa = alb.transform_to(aa)
    alb_wdistaa = alb_wdist.transform_to(aa)
    de_wdistaa = de_wdist.transform_to(aa)

    # these work fine
    sepnod = deaa.separation(albaa)
    sepwd = deaa.separation(alb_wdistaa)
    assert_quantity_allclose(sepnod, 22.2862*u.deg, rtol=1e-6)
    assert_quantity_allclose(sepwd, 22.2862*u.deg, rtol=1e-6)
    # parallax should be present when distance added
    assert np.abs(sepnod - sepwd) > 1*u.marcsec

    # in 4033, the following fail with a recursion error
    assert_quantity_allclose(de_wdistaa.separation(alb_wdistaa), 22.2862*u.deg, rtol=1e-3)
    assert_quantity_allclose(alb_wdistaa.separation(deaa), 22.2862*u.deg, rtol=1e-3)


@pytest.mark.skipif(str('not HAS_SCIPY'))
@pytest.mark.skipif(str('OLDER_SCIPY'))
def test_regression_4082():
    """
    Issue: https://github.com/astropy/astropy/issues/4082
    """
    from .. import search_around_sky, search_around_3d
    cat = SkyCoord([10.076,10.00455], [18.54746, 18.54896], unit='deg')
    search_around_sky(cat[0:1], cat, seplimit=u.arcsec * 60, storekdtree=False)
    # in the issue, this raises a TypeError

    #also check 3d for good measure, although it's not really affected by this bug directly
    cat3d = SkyCoord([10.076,10.00455]*u.deg, [18.54746, 18.54896]*u.deg, distance=[0.1,1.5]*u.kpc)
    search_around_3d(cat3d[0:1], cat3d, 1*u.kpc, storekdtree=False)


def test_regression_4996():
    # this part is the actual regression test
    deltat = np.linspace(-12, 12, 1000)*u.hour
    times = Time('2012-7-13 00:00:00') + deltat
    suncoo = get_sun(times)
    assert suncoo.shape == (len(times),)

    # and this is an additional test to make sure more complex arrays work
    times2 = Time('2012-7-13 00:00:00') + deltat.reshape(10, 20, 5)
    suncoo2 = get_sun(times2)
    assert suncoo2.shape == times2.shape

    # this is intentially not allclose - they should be *exactly* the same
    assert np.all(suncoo.ra.ravel() == suncoo2.ra.ravel())


def test_regression_5209():
    "check that distances are not lost on SkyCoord init"
    time = Time('2015-01-01')
    # this regression test was written for a later astropy with get_moon, but we
    # can just mock it up
    #moon = get_moon(time)
    moon = GCRS(10*u.deg, 20*u.deg, distance=384000*u.km, obstime=time)
    new_coord = SkyCoord([moon])
    assert_quantity_allclose(new_coord[0].distance, moon.distance)


def test_regression_4293():
    """Really just an extra test on FK4 no e, after finding that the units
    were not always taken correctly.  This test is against explicitly doing
    the transformations on pp170 of Explanatory Supplement to the Astronomical
    Almanac (Seidelmann, 2005).

    See https://github.com/astropy/astropy/pull/4293#issuecomment-234973086
    """
    # Check all over sky, but avoiding poles (note that FK4 did not ignore
    # e terms within 10∘ of the poles...  see p170 of explan.supp.).
    ra, dec = np.meshgrid(np.arange(0, 359, 45), np.arange(-80, 81, 40))
    fk4 = FK4(ra.ravel() * u.deg, dec.ravel() * u.deg)

    Dc = -0.065838*u.arcsec
    Dd = +0.335299*u.arcsec
    # Dc * tan(obliquity), as given on p.170
    Dctano = -0.028553*u.arcsec

    fk4noe_dec = (fk4.dec - (Dd*np.cos(fk4.ra) -
                             Dc*np.sin(fk4.ra))*np.sin(fk4.dec) -
                  Dctano*np.cos(fk4.dec))
    fk4noe_ra = fk4.ra - (Dc*np.cos(fk4.ra) +
                          Dd*np.sin(fk4.ra)) / np.cos(fk4.dec)

    fk4noe = fk4.transform_to(FK4NoETerms)
    # Tolerance here just set to how well the coordinates match, which is much
    # better than the claimed accuracy of <1 mas for this first-order in
    # v_earth/c approximation.
    # Interestingly, if one divides by np.cos(fk4noe_dec) in the ra correction,
    # the match becomes good to 2 μas.
    assert_quantity_allclose(fk4noe.ra, fk4noe_ra, atol=11.*u.uas, rtol=0)
    assert_quantity_allclose(fk4noe.dec, fk4noe_dec, atol=3.*u.uas, rtol=0)
