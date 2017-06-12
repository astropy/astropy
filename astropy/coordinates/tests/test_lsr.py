# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from copy import deepcopy
import numpy as np

from ...tests.helper import quantity_allclose
from ... import units as u
from ...extern import six

from ..builtin_frames import LSR, ICRS, Galactic
from ..representation import CartesianRepresentation

# TODO: this is a work in progress
def test_lsr_sanity():

    # random numbers, but zero velocity in ICRS frame
    icrs = ICRS(ra=15.1241*u.deg, dec=17.5143*u.deg, distance=150.12*u.pc,
                pm_ra=0*u.mas/u.yr, pm_dec=0*u.mas/u.yr,
                radial_velocity=0*u.km/u.s)
    lsr = icrs.transform_to(LSR)

    lsr_diff = lsr.data.differentials[0]
    cart_lsr_vel = lsr_diff.represent_as(CartesianRepresentation, base=lsr.data)
    lsr_vel = ICRS(cart_lsr_vel)
    gal_lsr = lsr_vel.transform_to(Galactic).cartesian.xyz
    assert quantity_allclose(gal_lsr.to(u.km/u.s, u.dimensionless_angles()),
                             lsr.v_bary.cartesian.xyz)

    # moving with LSR velocity
    lsr = LSR(ra=15.1241*u.deg, dec=17.5143*u.deg, distance=150.12*u.pc,
              pm_ra=0*u.mas/u.yr, pm_dec=0*u.mas/u.yr,
              radial_velocity=0*u.km/u.s)
    icrs = lsr.transform_to(ICRS)

    icrs_diff = icrs.data.differentials[0]
    cart_vel = icrs_diff.represent_as(CartesianRepresentation, base=icrs.data)
    vel = ICRS(cart_vel)
    gal_icrs = vel.transform_to(Galactic).cartesian.xyz
    assert quantity_allclose(gal_icrs.to(u.km/u.s, u.dimensionless_angles()),
                             -lsr.v_bary.cartesian.xyz)
