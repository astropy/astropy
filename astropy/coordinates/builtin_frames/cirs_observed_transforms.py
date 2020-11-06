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

    obsgeoloc, obsgeovel = altaz_frame.location.get_gcrs_posvel(altaz_frame.obstime)
    if (np.any(cirs_coo.obstime != altaz_frame.obstime) or
            np.any(cirs_coo.obsgeoloc != obsgeoloc)):
        cirs_coo = cirs_coo.transform_to(CIRS(obstime=altaz_frame.obstime,
                                              obsgeoloc=obsgeoloc,
                                              obsgeovel=obsgeovel))

    # if the data are UnitSphericalRepresentation, we can skip the distance calculations
    is_unitspherical = (isinstance(cirs_coo.data, UnitSphericalRepresentation) or
                        cirs_coo.cartesian.x.unit == u.one)

    # We used to do "astrometric" corrections here, but these are no longer necesssary
    # CIRS has proper topocentric behaviour
    usrepr = cirs_coo.represent_as(UnitSphericalRepresentation)
    cirs_ra = usrepr.lon.to_value(u.radian)
    cirs_dec = usrepr.lat.to_value(u.radian)

    # first set up the astrometry context for CIRS<->AltAz
    astrom = erfa_astrom.get().apio(altaz_frame)
    az, zen, _, _, _ = erfa.atioq(cirs_ra, cirs_dec, astrom)

    if is_unitspherical:
        rep = UnitSphericalRepresentation(lat=u.Quantity(PIOVER2 - zen, u.radian, copy=False),
                                          lon=u.Quantity(az, u.radian, copy=False),
                                          copy=False)
    else:
        # since we've transformed to CIRS at the observatory location, just use CIRS distance
        rep = SphericalRepresentation(lat=u.Quantity(PIOVER2 - zen, u.radian, copy=False),
                                      lon=u.Quantity(az, u.radian, copy=False),
                                      distance=cirs_coo.distance,
                                      copy=False)
    return altaz_frame.realize_frame(rep)


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference, AltAz, CIRS)
def altaz_to_cirs(altaz_coo, cirs_frame):
    usrepr = altaz_coo.represent_as(UnitSphericalRepresentation)
    az = usrepr.lon.to_value(u.radian)
    zen = PIOVER2 - usrepr.lat.to_value(u.radian)

    # first set up the astrometry context for ICRS<->CIRS at the altaz_coo time
    astrom = erfa_astrom.get().apio(altaz_coo)

    # the 'A' indicates zen/az inputs
    cirs_ra, cirs_dec = erfa.atoiq('A', az, zen, astrom)*u.radian
    obsgeoloc, obsgeovel = altaz_coo.location.get_gcrs_posvel(altaz_coo.obstime)
    if isinstance(altaz_coo.data, UnitSphericalRepresentation) or altaz_coo.cartesian.x.unit == u.one:
        distance = None
    else:
        distance = altaz_coo.distance

    cirs_at_aa_time = CIRS(ra=cirs_ra, dec=cirs_dec, distance=distance,
                           obstime=altaz_coo.obstime,
                           obsgeoloc=obsgeoloc,
                           obsgeovel=obsgeovel)

    # this final transform may be a no-op if the obstimes and locations are the same
    return cirs_at_aa_time.transform_to(cirs_frame)


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference, AltAz, AltAz)
def altaz_to_altaz(from_coo, to_frame):
    # for now we just implement this through CIRS to make sure we get everything
    # covered
    return from_coo.transform_to(CIRS(obstime=from_coo.obstime)).transform_to(to_frame)
