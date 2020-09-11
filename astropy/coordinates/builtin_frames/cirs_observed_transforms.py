# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Contains the transformation functions for getting to "observed" systems from CIRS.
Currently that just means AltAz.
"""

import numpy as np
import erfa

from astropy import units as u
from astropy.coordinates.baseframe import frame_transform_graph
from astropy.coordinates.transformations import FunctionTransformWithFiniteDifference
from astropy.coordinates.representation import (SphericalRepresentation,
                                                UnitSphericalRepresentation)

from .cirs import CIRS
from .altaz import AltAz
from .utils import PIOVER2
from ..erfa_astrom import erfa_astrom


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference, CIRS, AltAz)
def cirs_to_altaz(cirs_coo, altaz_frame):
    if np.any(cirs_coo.obstime != altaz_frame.obstime):
        # the only frame attribute for the current CIRS is the obstime, but this
        # would need to be updated if a future change allowed specifying an
        # Earth location algorithm or something
        cirs_coo = cirs_coo.transform_to(CIRS(obstime=altaz_frame.obstime))

    # we use the same obstime everywhere now that we know they're the same
    obstime = cirs_coo.obstime

    # if the data are UnitSphericalRepresentation, we can skip the distance calculations
    is_unitspherical = (isinstance(cirs_coo.data, UnitSphericalRepresentation) or
                        cirs_coo.cartesian.x.unit == u.one)

    if is_unitspherical:
        usrepr = cirs_coo.represent_as(UnitSphericalRepresentation)
        cirs_ra = usrepr.lon.to_value(u.radian)
        cirs_dec = usrepr.lat.to_value(u.radian)
    else:
        # compute an "astrometric" ra/dec -i.e., the direction of the
        # displacement vector from the observer to the target in CIRS
        loccirs = altaz_frame.location.get_itrs(cirs_coo.obstime).transform_to(cirs_coo)
        diffrepr = (cirs_coo.cartesian - loccirs.cartesian).represent_as(UnitSphericalRepresentation)

        cirs_ra = diffrepr.lon.to_value(u.radian)
        cirs_dec = diffrepr.lat.to_value(u.radian)

    # first set up the astrometry context for CIRS<->AltAz
    astrom = erfa_astrom.get().apio13(altaz_frame)

    az, zen, _, _, _ = erfa.atioq(cirs_ra, cirs_dec, astrom)

    if is_unitspherical:
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


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference, AltAz, CIRS)
def altaz_to_cirs(altaz_coo, cirs_frame):
    usrepr = altaz_coo.represent_as(UnitSphericalRepresentation)
    az = usrepr.lon.to_value(u.radian)
    zen = PIOVER2 - usrepr.lat.to_value(u.radian)

    # first set up the astrometry context for ICRS<->CIRS at the altaz_coo time
    astrom = erfa_astrom.get().apio13(altaz_coo)

    # the 'A' indicates zen/az inputs
    cirs_ra, cirs_dec = erfa.atoiq('A', az, zen, astrom)*u.radian
    if isinstance(altaz_coo.data, UnitSphericalRepresentation) or altaz_coo.cartesian.x.unit == u.one:
        cirs_at_aa_time = CIRS(ra=cirs_ra, dec=cirs_dec, distance=None,
                               obstime=altaz_coo.obstime)
    else:
        # treat the output of atoiq as an "astrometric" RA/DEC, so to get the
        # actual RA/Dec from the observers vantage point, we have to reverse
        # the vector operation of cirs_to_altaz (see there for more detail)

        loccirs = altaz_coo.location.get_itrs(altaz_coo.obstime).transform_to(cirs_frame)

        astrometric_rep = SphericalRepresentation(lon=cirs_ra, lat=cirs_dec,
                                                  distance=altaz_coo.distance)
        newrepr = astrometric_rep + loccirs.cartesian
        cirs_at_aa_time = CIRS(newrepr, obstime=altaz_coo.obstime)

    # this final transform may be a no-op if the obstimes are the same
    return cirs_at_aa_time.transform_to(cirs_frame)


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference, AltAz, AltAz)
def altaz_to_altaz(from_coo, to_frame):
    # for now we just implement this through CIRS to make sure we get everything
    # covered
    return from_coo.transform_to(CIRS(obstime=from_coo.obstime)).transform_to(to_frame)
