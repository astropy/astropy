# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import pytest
import numpy as np
from ...tests.helper import quantity_allclose

from ... import units as u
from ...time import Time
from ..builtin_frames import ICRS, AltAz, LSR, GCRS
from ..baseframe import frame_transform_graph
from .. import (EarthLocation, TimeFrameAttribute,
                FunctionTransformWithFiniteDifference, get_sun,
                CartesianRepresentation, SphericalRepresentation)

J2000 = Time('J2000')

@pytest.mark.parametrize("dt, symmetric", [(1*u.second, True),
                                           (1*u.year, True),
                                           (1*u.second, False),
                                           (1*u.year, False)])
def test_faux_lsr(dt, symmetric):
    class LSR2(LSR):
        obstime = TimeFrameAttribute(default=J2000)

    dt = 1*u.s
    @frame_transform_graph.transform(FunctionTransformWithFiniteDifference,
                                     ICRS, LSR2, finite_difference_dt=dt,
                                     symmetric_finite_difference=symmetric)
    def icrs_to_lsr(icrs_coo, lsr_frame):
        dt = lsr_frame.obstime - J2000
        offset = lsr_frame.v_bary.cartesian * dt.to(u.second)
        return  lsr_frame.realize_frame(icrs_coo.data.without_differentials() + offset)

    @frame_transform_graph.transform(FunctionTransformWithFiniteDifference,
                                     LSR2, ICRS, finite_difference_dt=dt,
                                     symmetric_finite_difference=symmetric)
    def lsr_to_icrs(lsr_coo, icrs_frame):
        dt = lsr_frame.obstime - J2000
        offset = lsr_frame.v_bary.cartesian * dt.to(u.second)
        return icrs_frame.realize_frame(lsr_coo.data - offset)


    ic = ICRS(ra=12.3*u.deg, dec=45.6*u.deg, distance=7.8*u.au,
              pm_ra=0*u.marcsec/u.yr, pm_dec=0*u.marcsec/u.yr,
              radial_velocity=0*u.km/u.s)
    lsrc = ic.transform_to(LSR2())

    assert quantity_allclose(ic.cartesian.xyz, lsrc.cartesian.xyz)

    idiff = ic.data.to_cartesian(True).differentials[0]
    ldiff = lsrc.data.to_cartesian(True).differentials[0]
    change = (ldiff.d_xyz - idiff.d_xyz).to(u.km/u.s)
    totchange = np.sum(change**2)**0.5
    assert quantity_allclose(totchange, lsrc.v_bary.distance)


    ic2 = ICRS(ra=120.3*u.deg, dec=45.6*u.deg, distance=7.8*u.au,
              pm_ra=0*u.marcsec/u.yr, pm_dec=10*u.marcsec/u.yr,
              radial_velocity=1000*u.km/u.s)
    lsrc2 = ic2.transform_to(LSR2())

    tot = np.sum(lsrc2.data.to_cartesian(True).differentials[0].d_xyz**2)**0.5
    assert np.abs(tot.to('km/s') - 1000*u.km/u.s) < 20*u.km/u.s

def test_gcrs_diffs():
    time = Time('J2017')
    gf = GCRS(obstime=time)
    sung = get_sun(time)  # should have very little vhelio

    # qtr-year off sun location should be the direction of ~ maximal vhelio
    qtrsung = get_sun(time-.25*u.year)

    #now we use those directions to
    msungr = CartesianRepresentation(-sung.cartesian.xyz).represent_as(SphericalRepresentation)
    suni = ICRS(ra=msungr.lon, dec=msungr.lat, distance=100*u.au,
                pm_ra=0*u.marcsec/u.yr, pm_dec=0*u.marcsec/u.yr,
                radial_velocity=0*u.km/u.s)
    qtrsuni = ICRS(ra=qtrsung.ra, dec=qtrsung.dec, distance=100*u.au,
                   pm_ra=0*u.marcsec/u.yr, pm_dec=0*u.marcsec/u.yr,
                   radial_velocity=0*u.km/u.s)
    sung = suni.transform_to(gf)
    qtrsung = qtrsuni.transform_to(gf)

    # should be high along the ecliptic-not-sun sun axis and
    # low along the sun axis
    assert np.abs(qtrsung.radial_velocity) > 30*u.km/u.s
    assert np.abs(qtrsung.radial_velocity) < 40*u.km/u.s
    assert np.abs(sung.radial_velocity) < 1*u.km/u.s

    suni2 = sung.transform_to(ICRS)
    assert np.all(suni2.data.differentials[0].d_xyz < 1e-5*u.km/u.s)
    qtrisun2 = qtrsung.transform_to(ICRS)
    assert np.all(qtrisun2.data.differentials[0].d_xyz < 1e-5*u.km/u.s)
