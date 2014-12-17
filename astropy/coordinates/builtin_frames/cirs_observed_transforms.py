# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Contains the transformation functions for getting to "observed" systems from CIRS.
Currently that just means AltAz.
"""
from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

import numpy as np

from ... import units as u
from ..baseframe import frame_transform_graph
from ..transformations import FunctionTransform
from ..representation import UnitSphericalRepresentation
from ... import erfa

from .cirs import CIRS
from .altaz import AltAz
from .utils import get_polar_motion, PIOVER2


@frame_transform_graph.transform(FunctionTransform, CIRS, AltAz)
def cirs_to_altaz(cirs_coo, altaz_frame):
    if np.all(cirs_coo.obstime != altaz_frame.obstime):
        # the only frame attribute for the current CIRS is the obstime, but this
        # would need to be updated if a future change allowed specifying an
        # Earth location algorithm or something
        cirs_coo = cirs_coo.transform_to(CIRS(obstime=altaz_frame.obstime))
    srepr = cirs_coo.represent_as(UnitSphericalRepresentation)
    cirs_ra = srepr.lon.to(u.radian).value
    cirs_dec = srepr.lat.to(u.radian).value

    lon, lat, height = altaz_frame.location.to_geodetic('WGS84')
    xp, yp = get_polar_motion(cirs_coo.obstime)

    #first set up the astrometry context for ICRS<->CIRS
    astrom = erfa.apio13(cirs_coo.obstime.jd1, cirs_coo.obstime.jd2,
                         cirs_coo.obstime.delta_ut1_utc[0],
                         lon.to(u.radian).value, lat.to(u.radian).value,
                         height.to(u.m).value,
                         xp, yp,  # polar motion
                         # all below are already in correct units because they are QuantityFrameAttribues
                         altaz_frame.pressure.value,
                         altaz_frame.temperature.value,
                         altaz_frame.relative_humidity,
                         altaz_frame.obswl.value)

    az, zen, ha, obs_dec, obs_ra = erfa.atioq(cirs_ra, cirs_dec, astrom)

    rep = UnitSphericalRepresentation(lat=u.Quantity(PIOVER2 - zen, u.radian, copy=False),
                                      lon=u.Quantity(az, u.radian, copy=False),
                                      copy=False)
    return altaz_frame.realize_frame(rep)


@frame_transform_graph.transform(FunctionTransform, AltAz, CIRS)
def altaz_to_cirs(altaz_coo, cirs_frame):
    srepr = altaz_coo.represent_as(UnitSphericalRepresentation)
    az = srepr.lon.to(u.radian).value
    zen = PIOVER2 - srepr.lat.to(u.radian).value

    lon, lat, height = altaz_coo.location.to_geodetic('WGS84')
    xp, yp = get_polar_motion(altaz_coo.obstime)

    #first set up the astrometry context for ICRS<->CIRS at the altaz_coo time
    astrom = erfa.apio13(altaz_coo.obstime.jd1, altaz_coo.obstime.jd2,
                         altaz_coo.obstime.delta_ut1_utc[0],
                         lon.to(u.radian).value, lat.to(u.radian).value,
                         height.to(u.m).value,
                         xp, yp,  # polar motion
                         # all below are already in correct units because they are QuantityFrameAttribues
                         altaz_coo.pressure.value,
                         altaz_coo.temperature.value,
                         altaz_coo.relative_humidity,
                         altaz_coo.obswl.value)

    # the 'A' indicates zen/az inputs
    cirs_ra, cirs_dec = erfa.atoiq('A', az, zen, astrom)

    #the final transform may be a no-op if the obstimes are the same
    return CIRS(ra=cirs_ra*u.radian, dec=cirs_dec*u.radian,
                obstime=altaz_coo.obstime).transform_to(cirs_frame)


@frame_transform_graph.transform(FunctionTransform, AltAz, AltAz)
def altaz_to_altaz(from_coo, to_frame):
    # for now we just implement this through CIRS to make sure we get everything
    # covered
    return from_coo.transform_to(CIRS(obstime=from_coo.obstime)).transform_to(to_frame)
