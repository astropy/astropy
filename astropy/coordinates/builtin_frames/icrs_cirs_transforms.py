# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Contains the transofrmation functions for getting from ICRS to CIRS and anything
in between (currently that means GCRS)
"""
from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

import numpy as np

from ... import units as u
from ..baseframe import frame_transform_graph
from ..transformations import FunctionTransform
from ..representation import UnitSphericalRepresentation, SphericalRepresentation
from ... import erfa

from .icrs import ICRS
from .gcrs import GCRS
from .cirs import CIRS


#first the ICRS/CIRS related transforms

@frame_transform_graph.transform(FunctionTransform, ICRS, CIRS)
def icrs_to_cirs(icrs_coo, cirs_frame):
    #parallax in arcsec
    if isinstance(icrs_coo.data, UnitSphericalRepresentation):  # no distance
        px = 0
    else:
        px = 1 / icrs_coo.distance.to(u.parsec).value
    srepr = icrs_coo.represent_as(UnitSphericalRepresentation)
    i_ra = srepr.lon.to(u.radian).value
    i_dec = srepr.lat.to(u.radian).value

    #first set up the astrometry context for ICRS<->CIRS
    astrom, eo = erfa.apci13(cirs_frame.obstime.jd1, cirs_frame.obstime.jd2)

    # TODO: possibly switch to something that is like atciq, but skips the first
    # step, b/c that involves some  wasteful computations b/c pm=0  here
    cirs_ra, cirs_dec = erfa.atciq(i_ra, i_dec, 0, 0, px, 0, astrom)

    dat = icrs_coo.data
    if dat.get_name() == 'unitspherical'  or dat.to_cartesian().x.unit == u.one:
        rep = UnitSphericalRepresentation(lat=u.Quantity(cirs_dec, u.radian, copy=False),
                                          lon=u.Quantity(cirs_ra, u.radian, copy=False),
                                          copy=False)
    else:
        #compute the distance as just the cartesian sum from moving to the Earth
        #we have to do this because the ERFA functions throw away distance info
        distance = np.sum((astrom['eb']*u.au + icrs_coo.cartesian.xyz.T)**2, -1)**0.5
        rep = SphericalRepresentation(lat=u.Quantity(cirs_dec, u.radian, copy=False),
                                      lon=u.Quantity(cirs_ra, u.radian, copy=False),
                                      distance=distance, copy=False)

    return cirs_frame.realize_frame(rep)

@frame_transform_graph.transform(FunctionTransform, CIRS, ICRS)
def cirs_to_icrs(cirs_coo, icrs_frame):
    srepr = cirs_coo.represent_as(UnitSphericalRepresentation)
    cirs_ra = srepr.lon.to(u.radian).value
    cirs_dec = srepr.lat.to(u.radian).value

    #first set up the astrometry context for ICRS<->CIRS
    astrom, eo = erfa.apci13(cirs_coo.obstime.jd1, cirs_coo.obstime.jd2)

    icrs_ra, icrs_dec = erfa.aticq(cirs_ra, cirs_dec, astrom)

    dat = cirs_coo.data
    if dat.get_name() == 'unitspherical'  or dat.to_cartesian().x.unit == u.one:
        rep = UnitSphericalRepresentation(lat=u.Quantity(icrs_dec, u.radian, copy=False),
                                          lon=u.Quantity(icrs_ra, u.radian, copy=False),
                                          copy=False)
    else:
        #compute the distance as just the cartesian sum from moving to the SSB
        #we have to do this because the ERFA functions throw away distance info
        distance = np.sum((-astrom['eb']*u.au + cirs_coo.cartesian.xyz.T)**2, -1)**0.5
        rep = SphericalRepresentation(lat=u.Quantity(icrs_dec, u.radian, copy=False),
                                      lon=u.Quantity(icrs_ra, u.radian, copy=False),
                                      distance=distance, copy=False)
    return icrs_frame.realize_frame(rep)

@frame_transform_graph.transform(FunctionTransform, CIRS, CIRS)
def cirs_to_cirs(from_coo, to_frame):
    if np.all(from_coo.obstime == to_frame.obstime):
        return to_frame.realize_frame(from_coo.data)
    else:
        # the CIRS<-> CIRS transform actually goes through ICRS.  This has a
        # subtle implication that a point in CIRS is uniquely determined
        # by the corresponding astrometric ICRS coordinate *at its
        # current time*.  This has some subtle implications in terms of GR, but
        # is sort of glossed over in the current scheme because we are dropping
        # distances anyway.
        return from_coo.transform_to(ICRS).transform_to(to_frame)


# Now the GCRS-related transforms to/from ICRS

@frame_transform_graph.transform(FunctionTransform, ICRS, GCRS)
def icrs_to_gcrs(icrs_coo, gcrs_frame):
    #parallax in arcsec
    if isinstance(icrs_coo.data, UnitSphericalRepresentation):  # no distance
        px = 0
    else:
        px = 1 / icrs_coo.distance.to(u.parsec).value
    srepr = icrs_coo.represent_as(UnitSphericalRepresentation)
    i_ra = srepr.lon.to(u.radian).value
    i_dec = srepr.lat.to(u.radian).value

    #first set up the astrometry context for ICRS<->GCRS
    pv = np.array([gcrs_frame.obsgeoloc.value,
                   gcrs_frame.obsgeovel.value])
    astrom = erfa.apcs13(gcrs_frame.obstime.jd1, gcrs_frame.obstime.jd2, pv)

    # TODO: possibly switch to something that is like atciq, but skips the first
    # step, b/c that involves some  wasteful computations b/c pm=0  here
    gcrs_ra, gcrs_dec = erfa.atciq(i_ra, i_dec, 0, 0, px, 0, astrom)

    dat = icrs_coo.data
    if dat.get_name() == 'unitspherical'  or dat.to_cartesian().x.unit == u.one:
        rep = UnitSphericalRepresentation(lat=u.Quantity(gcrs_dec, u.radian, copy=False),
                                          lon=u.Quantity(gcrs_ra, u.radian, copy=False),
                                          copy=False)
    else:
        #compute the distance as just the cartesian sum from moving to the Earth
        #we have to do this because the ERFA functions throw away distance info
        distance = np.sum((astrom['eb']*u.au + icrs_coo.cartesian.xyz.T)**2, -1)**0.5
        rep = SphericalRepresentation(lat=u.Quantity(gcrs_dec, u.radian, copy=False),
                                      lon=u.Quantity(gcrs_ra, u.radian, copy=False),
                                      distance=distance, copy=False)
    return gcrs_frame.realize_frame(rep)


@frame_transform_graph.transform(FunctionTransform, GCRS, ICRS)
def gcrs_to_icrs(gcrs_coo, icrs_frame):
    srepr = gcrs_coo.represent_as(UnitSphericalRepresentation)
    cirs_ra = srepr.lon.to(u.radian).value
    cirs_dec = srepr.lat.to(u.radian).value

    #first set up the astrometry context for ICRS<->GCRS
    pv = np.array([gcrs_coo.obsgeoloc.value,
                   gcrs_coo.obsgeovel.value])
    astrom = erfa.apcs13(gcrs_coo.obstime.jd1, gcrs_coo.obstime.jd2, pv)

    icrs_ra, icrs_dec = erfa.aticq(cirs_ra, cirs_dec, astrom)

    dat = gcrs_coo.data
    if dat.get_name() == 'unitspherical'  or dat.to_cartesian().x.unit == u.one:
        rep = UnitSphericalRepresentation(lat=u.Quantity(icrs_dec, u.radian, copy=False),
                                          lon=u.Quantity(icrs_ra, u.radian, copy=False),
                                          copy=False)
    else:
        #compute the distance as just the cartesian sum from moving to the SSB
        #we have to do this because the ERFA functions throw away distance info
        distance = np.sum((-astrom['eb']*u.au + gcrs_coo.cartesian.xyz.T)**2, -1)**0.5
        rep = SphericalRepresentation(lat=u.Quantity(icrs_dec, u.radian, copy=False),
                                      lon=u.Quantity(icrs_ra, u.radian, copy=False),
                                      distance=distance, copy=False)
    return icrs_frame.realize_frame(rep)


@frame_transform_graph.transform(FunctionTransform, GCRS, GCRS)
def gcrs_to_gcrs(from_coo, to_frame):
    if np.all(from_coo.obstime == to_frame.obstime):
        return to_frame.realize_frame(from_coo.data)
    else:
        # like CIRS, we do this self-transform via ICRS
        return from_coo.transform_to(ICRS).transform_to(to_frame)
