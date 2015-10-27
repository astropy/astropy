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
from ..representation import (SphericalRepresentation, CartesianRepresentation,
                              UnitSphericalRepresentation)
from ... import _erfa as erfa

from .cirs import CIRS
from .altaz import AltAz
from .utils import get_polar_motion, get_dut1utc, PIOVER2


@frame_transform_graph.transform(FunctionTransform, CIRS, AltAz)
def cirs_to_altaz(cirs_coo, altaz_frame):
    if np.any(cirs_coo.obstime != altaz_frame.obstime):
        # the only frame attribute for the current CIRS is the obstime, but this
        # would need to be updated if a future change allowed specifying an
        # Earth location algorithm or something
        cirs_coo = cirs_coo.transform_to(CIRS(obstime=altaz_frame.obstime))

    # we use the same obstime everywhere now that we know they're the same
    obstime = cirs_coo.obstime

    usrepr = cirs_coo.represent_as(UnitSphericalRepresentation)
    cirs_ra = usrepr.lon.to(u.radian).value
    cirs_dec = usrepr.lat.to(u.radian).value

    lon, lat, height = altaz_frame.location.to_geodetic('WGS84')
    xp, yp = get_polar_motion(obstime)

    #first set up the astrometry context for CIRS<->AltAz
    astrom = erfa.apio13(obstime.jd1, obstime.jd2,
                         get_dut1utc(obstime),
                         lon.to(u.radian).value, lat.to(u.radian).value,
                         height.to(u.m).value,
                         xp, yp,  # polar motion
                         # all below are already in correct units because they are QuantityFrameAttribues
                         altaz_frame.pressure.value,
                         altaz_frame.temperature.value,
                         altaz_frame.relative_humidity,
                         altaz_frame.obswl.value)

    az, zen, ha, obs_dec, obs_ra = erfa.atioq(cirs_ra, cirs_dec, astrom)

    dat = cirs_coo.data
    if dat.get_name() == 'unitspherical'  or dat.to_cartesian().x.unit == u.one:
        rep = UnitSphericalRepresentation(lat=u.Quantity(PIOVER2 - zen, u.radian, copy=False),
                                          lon=u.Quantity(az, u.radian, copy=False),
                                          copy=False)
    else:
        # now we get the distance as the cartesian distance from the earth
        # location to the coordinate location
        locitrs = altaz_frame.location.get_itrs(obstime)
        distance = locitrs.separation_3d(cirs_coo)
        rep = SphericalRepresentation(lat=u.Quantity(PIOVER2 - zen, u.radian, copy=False),
                                      lon=u.Quantity(az, u.radian, copy=False),
                                      distance=distance,
                                      copy=False)
    return altaz_frame.realize_frame(rep)


@frame_transform_graph.transform(FunctionTransform, AltAz, CIRS)
def altaz_to_cirs(altaz_coo, cirs_frame):
    usrepr = altaz_coo.represent_as(UnitSphericalRepresentation)
    az = usrepr.lon.to(u.radian).value
    zen = PIOVER2 - usrepr.lat.to(u.radian).value

    lon, lat, height = altaz_coo.location.to_geodetic('WGS84')
    xp, yp = get_polar_motion(altaz_coo.obstime)

    #first set up the astrometry context for ICRS<->CIRS at the altaz_coo time
    astrom = erfa.apio13(altaz_coo.obstime.jd1, altaz_coo.obstime.jd2,
                         get_dut1utc(altaz_coo.obstime),
                         lon.to(u.radian).value, lat.to(u.radian).value,
                         height.to(u.m).value,
                         xp, yp,  # polar motion
                         # all below are already in correct units because they are QuantityFrameAttribues
                         altaz_coo.pressure.value,
                         altaz_coo.temperature.value,
                         altaz_coo.relative_humidity,
                         altaz_coo.obswl.value)

    # the 'A' indicates zen/az inputs
    cirs_ra, cirs_dec = erfa.atoiq('A', az, zen, astrom)*u.radian

    if isinstance(altaz_coo.data, UnitSphericalRepresentation):
        distance = None
    else:
        locitrs = altaz_coo.location.get_itrs(altaz_coo.obstime)

        # To compute the distance in a way that is reversable with cirs_to_altaz
        # we use basic trigonometry.  The altaz_coo's distance ("d") is one leg
        # of a triangle, and the earth center -> EarthLocation distance ("r")
        # is a neighboring leg.  We can also easily get the angle between the
        # earth center->target (calculated above by apio13), and
        # earth center->EarthLocation vectors.  This is a Side-Side-Angle
        # situation, and the formula below is the trig formula to solve for
        # the remaining side

        ucirs = cirs_frame.realize_frame(UnitSphericalRepresentation(lon=cirs_ra, lat=cirs_dec))
        delta = ucirs.separation(locitrs)
        r = locitrs.spherical.distance
        d = altaz_coo.distance

        sindoverd = np.sin(delta) / d
        distance = np.sin(delta + np.arcsin(r*sindoverd))/sindoverd

    #the final transform may be a no-op if the obstimes are the same
    return CIRS(ra=cirs_ra, dec=cirs_dec, distance=distance,
                obstime=altaz_coo.obstime).transform_to(cirs_frame)


@frame_transform_graph.transform(FunctionTransform, AltAz, AltAz)
def altaz_to_altaz(from_coo, to_frame):
    # for now we just implement this through CIRS to make sure we get everything
    # covered
    return from_coo.transform_to(CIRS(obstime=from_coo.obstime)).transform_to(to_frame)
