# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Contains the transformation functions for getting to any system "lower" than
CIRS (various forms of earth-centric or observer-oriented systems).
"""
from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

import numpy as np

from ... import units as u
from ..baseframe import frame_transform_graph
from ..transformations import FunctionTransform
from ..representation import UnitSphericalRepresentation
from ... import erfa
from ...utils import iers

from .cirs import CIRS
from .itrs import ITRS
from .altaz import AltAz

_PIOVER2 = np.pi / 2.

@frame_transform_graph.transform(FunctionTransform, CIRS, AltAz)
def cirs_to_altaz(cirs_coo, altaz_frame):
    if cirs_coo.obstime != altaz_frame.obstime:
        # the only frame attribute for the current CIRS is the obstime, but this
        # would need to be updated if a future change allowed specifying an
        # Earth location algorithm or something
        cirs_coo = cirs_coo.transform_to(CIRS(obstime=altaz_frame.obstime))
    cirs_ra = cirs_coo.ra.to(u.radian).value
    cirs_dec = cirs_coo.dec.to(u.radian).value

    lon, lat, height = altaz_frame.location.geodetic  # assume EarthLocation

    #get the polar motion from the IERS table
    xp, yp = iers.IERS.open().pm_xy(cirs_coo.obstime.jd1, cirs_coo.obstime.jd2)

    #first set up the astrometry context for ICRS<->CIRS
    astrom = erfa.apio13(cirs_coo.obstime.jd1, cirs_coo.obstime.jd2,
                         cirs_coo.obstime.delta_ut1_utc[0],
                         lon.to(u.radian).value, lat.to(u.radian).value,
                         height.to(u.m).value,
                         xp.to(u.radian).value, yp.to(u.radian).value,
                         # all below are already in correct units
                         altaz_frame.pressure.value,
                         altaz_frame.temperature.value,
                         altaz_frame.relative_humidity,
                         altaz_frame.obswl.value)

    az, zen, ha, obs_dec, obs_ra = erfa.atioq(cirs_ra, cirs_dec, astrom)

    rep = UnitSphericalRepresentation(lat=u.Quantity(_PIOVER2 - zen, u.radian, copy=False),
                                      lon=u.Quantity(az, u.radian, copy=False),
                                      copy=False)
    return altaz_frame.realize_frame(rep)
