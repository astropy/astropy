# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import pytest
import numpy as np
from ...tests.helper import quantity_allclose

from ... import units as u
from ... import constants
from ...time import Time
from ..builtin_frames import ICRS, AltAz, LSR, GCRS, Galactic, FK5
from ..baseframe import frame_transform_graph
from ..sites import get_builtin_sites
from .. import (TimeAttribute,
                FunctionTransformWithFiniteDifference, get_sun,
                CartesianRepresentation, SphericalRepresentation,
                CartesianDifferential, SphericalDifferential,
                DynamicMatrixTransform)

J2000 = Time('J2000')


@pytest.mark.parametrize("dt, symmetric", [(1*u.second, True),
                                           (1*u.year, True),
                                           (1*u.second, False),
                                           (1*u.year, False)])
def test_faux_lsr(dt, symmetric):
    class LSR2(LSR):
        obstime = TimeAttribute(default=J2000)

    @frame_transform_graph.transform(FunctionTransformWithFiniteDifference,
                                     ICRS, LSR2, finite_difference_dt=dt,
                                     symmetric_finite_difference=symmetric)
    def icrs_to_lsr(icrs_coo, lsr_frame):
        dt = lsr_frame.obstime - J2000
        offset = lsr_frame.v_bary * dt.to(u.second)
        return lsr_frame.realize_frame(icrs_coo.data.without_differentials() + offset)

    @frame_transform_graph.transform(FunctionTransformWithFiniteDifference,
                                     LSR2, ICRS, finite_difference_dt=dt,
                                     symmetric_finite_difference=symmetric)
    def lsr_to_icrs(lsr_coo, icrs_frame):
        dt = lsr_coo.obstime - J2000
        offset = lsr_coo.v_bary * dt.to(u.second)
        return icrs_frame.realize_frame(lsr_coo.data - offset)

    ic = ICRS(ra=12.3*u.deg, dec=45.6*u.deg, distance=7.8*u.au,
              pm_ra_cosdec=0*u.marcsec/u.yr, pm_dec=0*u.marcsec/u.yr,
              radial_velocity=0*u.km/u.s)
    lsrc = ic.transform_to(LSR2())

    assert quantity_allclose(ic.cartesian.xyz, lsrc.cartesian.xyz)

    idiff = ic.cartesian.differentials['s']
    ldiff = lsrc.cartesian.differentials['s']
    change = (ldiff.d_xyz - idiff.d_xyz).to(u.km/u.s)
    totchange = np.sum(change**2)**0.5
    assert quantity_allclose(totchange, np.sum(lsrc.v_bary.d_xyz**2)**0.5)

    ic2 = ICRS(ra=120.3*u.deg, dec=45.6*u.deg, distance=7.8*u.au,
              pm_ra_cosdec=0*u.marcsec/u.yr, pm_dec=10*u.marcsec/u.yr,
              radial_velocity=1000*u.km/u.s)
    lsrc2 = ic2.transform_to(LSR2())

    tot = np.sum(lsrc2.cartesian.differentials['s'].d_xyz**2)**0.5
    assert np.abs(tot.to('km/s') - 1000*u.km/u.s) < 20*u.km/u.s


def test_faux_fk5_galactic():

    from ..builtin_frames.galactic_transforms import fk5_to_gal, _gal_to_fk5

    class Galactic2(Galactic):
        pass

    dt = 1000*u.s

    @frame_transform_graph.transform(FunctionTransformWithFiniteDifference,
                                     FK5, Galactic2, finite_difference_dt=dt,
                                     symmetric_finite_difference=True,
                                     finite_difference_frameattr_name=None)
    def fk5_to_gal2(fk5_coo, gal_frame):
        trans = DynamicMatrixTransform(fk5_to_gal, FK5, Galactic2)
        return trans(fk5_coo, gal_frame)

    @frame_transform_graph.transform(FunctionTransformWithFiniteDifference,
                                     Galactic2, ICRS, finite_difference_dt=dt,
                                     symmetric_finite_difference=True,
                                     finite_difference_frameattr_name=None)
    def gal2_to_fk5(gal_coo, fk5_frame):
        trans = DynamicMatrixTransform(_gal_to_fk5, Galactic2, FK5)
        return trans(gal_coo, fk5_frame)

    c1 = FK5(ra=150*u.deg, dec=-17*u.deg, radial_velocity=83*u.km/u.s,
             pm_ra_cosdec=-41*u.mas/u.yr, pm_dec=16*u.mas/u.yr,
             distance=150*u.pc)
    c2 = c1.transform_to(Galactic2)
    c3 = c1.transform_to(Galactic)

    # compare the matrix and finite-difference calculations
    assert quantity_allclose(c2.pm_l_cosb, c3.pm_l_cosb, rtol=1e-4)
    assert quantity_allclose(c2.pm_b, c3.pm_b, rtol=1e-4)


def test_gcrs_diffs():
    time = Time('J2017')
    gf = GCRS(obstime=time)
    sung = get_sun(time)  # should have very little vhelio

    # qtr-year off sun location should be the direction of ~ maximal vhelio
    qtrsung = get_sun(time-.25*u.year)

    # now we use those essentially as directions where the velocities should
    # be either maximal or minimal - with or perpendiculat to Earh's orbit
    msungr = CartesianRepresentation(-sung.cartesian.xyz).represent_as(SphericalRepresentation)
    suni = ICRS(ra=msungr.lon, dec=msungr.lat, distance=100*u.au,
                pm_ra_cosdec=0*u.marcsec/u.yr, pm_dec=0*u.marcsec/u.yr,
                radial_velocity=0*u.km/u.s)
    qtrsuni = ICRS(ra=qtrsung.ra, dec=qtrsung.dec, distance=100*u.au,
                   pm_ra_cosdec=0*u.marcsec/u.yr, pm_dec=0*u.marcsec/u.yr,
                   radial_velocity=0*u.km/u.s)

    # Now we transform those parallel- and perpendicular-to Earth's orbit
    # directions to GCRS, which should shift the velocity to either include
    # the Earth's velocity vector, or not (for parallel and perpendicular,
    # respectively).
    sung = suni.transform_to(gf)
    qtrsung = qtrsuni.transform_to(gf)

    # should be high along the ecliptic-not-sun sun axis and
    # low along the sun axis
    assert np.abs(qtrsung.radial_velocity) > 30*u.km/u.s
    assert np.abs(qtrsung.radial_velocity) < 40*u.km/u.s
    assert np.abs(sung.radial_velocity) < 1*u.km/u.s

    suni2 = sung.transform_to(ICRS)
    assert np.all(np.abs(suni2.data.differentials['s'].d_xyz) < 3e-5*u.km/u.s)
    qtrisun2 = qtrsung.transform_to(ICRS)
    assert np.all(np.abs(qtrisun2.data.differentials['s'].d_xyz) < 3e-5*u.km/u.s)


def test_altaz_diffs():
    time = Time('J2015') + np.linspace(-1, 1, 1000)*u.day
    loc = get_builtin_sites()['greenwich']
    aa = AltAz(obstime=time, location=loc)

    icoo = ICRS(np.zeros_like(time)*u.deg, 10*u.deg, 100*u.au,
                pm_ra_cosdec=np.zeros_like(time)*u.marcsec/u.yr,
                pm_dec=0*u.marcsec/u.yr,
                radial_velocity=0*u.km/u.s)

    acoo = icoo.transform_to(aa)

    # Make sure the change in radial velocity over ~2 days isn't too much
    # more than the rotation speed of the Earth - some excess is expected
    # because the orbit also shifts the RV, but it should be pretty small
    # over this short a time.
    assert np.ptp(acoo.radial_velocity)/2 < (2*np.pi*constants.R_earth/u.day)*1.2  # MAGIC NUMBER

    cdiff = acoo.data.differentials['s'].represent_as(CartesianDifferential,
                                                    acoo.data)

    # The "total" velocity should be > c, because the *tangential* velocity
    # isn't a True velocity, but rather an induced velocity due to the Earth's
    # rotation at a distance of 100 AU
    assert np.all(np.sum(cdiff.d_xyz**2, axis=0)**0.5 > constants.c)


_xfail = pytest.mark.xfail


@pytest.mark.parametrize('distance', [1000*u.au,
                                      10*u.pc,
                                      10*u.kpc,
                                      100*u.kpc])
def test_numerical_limits(distance):
    """
    Tests the numerical stability of the default settings for the finite
    difference transformation calculation.  This is *known* to fail for at
    >~1kpc, but this may be improved in future versions.
    """

    if distance.unit == u.kpc:
        # pytest.mark.parametrize syntax changed in pytest 3.1 to handle
        # directly marking xfails, thus the workaround below to support
        # pytest <3.1 for the 2.0.x LTS
        pytest.xfail()

    time = Time('J2017') + np.linspace(-.5, .5, 100)*u.year

    icoo = ICRS(ra=0*u.deg, dec=10*u.deg, distance=distance,
                pm_ra_cosdec=0*u.marcsec/u.yr, pm_dec=0*u.marcsec/u.yr,
                radial_velocity=0*u.km/u.s)
    gcoo = icoo.transform_to(GCRS(obstime=time))
    rv = gcoo.radial_velocity.to('km/s')

    # if its a lot bigger than this - ~the maximal velocity shift along
    # the direction above with a small allowance for noise - finite-difference
    # rounding errors have ruined the calculation
    assert np.ptp(rv) < 65*u.km/u.s


def diff_info_plot(frame, time):
    """
    Useful for plotting a frame with multiple times. *Not* used in the testing
    suite per se, but extremely useful for interactive plotting of results from
    tests in this module.
    """
    from matplotlib import pyplot as plt

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(20, 12))
    ax1.plot_date(time.plot_date, frame.data.differentials['s'].d_xyz.to(u.km/u.s).T, fmt='-')
    ax1.legend(['x', 'y', 'z'])

    ax2.plot_date(time.plot_date, np.sum(frame.data.differentials['s'].d_xyz.to(u.km/u.s)**2, axis=0)**0.5, fmt='-')
    ax2.set_title('total')

    sd = frame.data.differentials['s'].represent_as(SphericalDifferential, frame.data)

    ax3.plot_date(time.plot_date, sd.d_distance.to(u.km/u.s), fmt='-')
    ax3.set_title('radial')

    ax4.plot_date(time.plot_date, sd.d_lat.to(u.marcsec/u.yr), fmt='-', label='lat')
    ax4.plot_date(time.plot_date, sd.d_lon.to(u.marcsec/u.yr), fmt='-', label='lon')

    return fig
