# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import pytest
import numpy as np
from ...tests.helper import quantity_allclose

from ... import units as u
from ...time import Time
from ..builtin_frames import ICRS, AltAz, LSR
from ..baseframe import frame_transform_graph
from .. import (EarthLocation, TimeFrameAttribute,
                FunctionTransformWithFiniteDifference)

J2000 = Time('J2000')

def test_faux_lsr():
    class LSR2(LSR):
        obstime = TimeFrameAttribute(default=J2000)

    @frame_transform_graph.transform(FunctionTransformWithFiniteDifference, ICRS, LSR2)
    def icrs_to_lsr(icrs_coo, lsr_frame):
        dt = lsr_frame.obstime - J2000
        offset = lsr_frame.v_bary.cartesian * dt
        return lsr_frame.realize_frame(icrs_coo.data.without_differentials() + offset)

        @frame_transform_graph.transform(FunctionTransformWithFiniteDifference, ICRS, LSR2)
        def lsr_to_icrs(lsr_coo, icrs_frame):
            dt = lsr_frame.obstime - J2000
            offset = lsr_frame.v_bary.cartesian * dt
            return icrs_frame.realize_frame(lsr_coo.data - offset)

    ic = ICRS(ra=0*u.deg, dec=0*u.deg, distance=10*u.kpc,
              pm_ra=0*u.marcsec/u.yr, pm_dec=0*u.marcsec/u.yr,
              radial_velocity=0*u.km/u.s)
    lsrc = ic.transform_to(LSR2())

    assert quantity_allclose(ic.cartesian.xyz, lsrc.cartesian.xyz)
    print(ic.data.differentials)
    print(lsrc.data.differentials)
    assert quantity_allclose(ic.data.to_cartesian(True).differentials[0],
                             lsrc.data.to_cartesian(True).differentials[0])
