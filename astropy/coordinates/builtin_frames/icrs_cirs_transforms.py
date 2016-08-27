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
from ..representation import (SphericalRepresentation, CartesianRepresentation,
                              UnitSphericalRepresentation)
from ... import _erfa as erfa

from .icrs import ICRS
from .gcrs import GCRS
from .cirs import CIRS
from .hcrs import HCRS
from .utils import get_jd12, aticq, atciqz


def prepare_earth_position_vel(time):
    """
    Get barycentric position and velocity, and heliocentric position of Earth

    Parameters
    -----------
    time : `~astropy.time.Time`
        time at which to calculate position and velocity of Earth

    Returns
    --------
    earth_pv : `np.ndarray`
        Barycentric position and velocity of Earth, in au and au/day
    earth_helio : `np.ndarray`
        Heliocentric position of Earth in au
    """
    # this goes here to avoid circular import errors
    from ..solar_system import get_body_barycentric, solar_system_ephemeris
    # get barycentric position and velocity of earth
    ephemeris = solar_system_ephemeris.get()
    earth_pv = get_body_barycentric('earth', time,
                                    ephemeris=ephemeris, get_velocity=True)

    # get heliocentric position of earth
    sun = get_body_barycentric('sun', time, ephemeris=ephemeris)
    earth_heliocentric = (earth_pv[0].xyz - sun.xyz).to(u.au)

    # prepare to pass to erfa
    earth_pv = np.array([earth_pv[0].xyz.to(u.au), earth_pv[1].xyz.to(u.au/u.d)])
    earth_pv = np.rollaxis(np.rollaxis(earth_pv, 0, earth_pv.ndim), 0, earth_pv.ndim)
    earth_heliocentric = np.rollaxis(earth_heliocentric, 0, earth_heliocentric.ndim)
    return earth_pv, earth_heliocentric


def get_cip(jd1, jd2):
    # classical NPB matrix, IAU 2006/2000A
    rpnb = erfa.pnm06a(jd1, jd2)
    # CIP X, Y coordinates from array
    x, y = erfa.bpn2xy(rpnb)
    # CIO locator, s
    s = erfa.s06(jd1, jd2, x, y)
    return x, y, s


# First the ICRS/CIRS related transforms
@frame_transform_graph.transform(FunctionTransform, ICRS, CIRS)
def icrs_to_cirs(icrs_coo, cirs_frame):
    # first set up the astrometry context for ICRS<->CIRS
    jd1, jd2 = get_jd12(cirs_frame.obstime, 'tdb')
    x, y, s = get_cip(jd1, jd2)
    earth_pv, earth_heliocentric = prepare_earth_position_vel(cirs_frame.obstime)
    astrom = erfa.apci(jd1, jd2, earth_pv, earth_heliocentric, x, y, s)

    if icrs_coo.data.get_name() == 'unitspherical' or icrs_coo.data.to_cartesian().x.unit == u.one:
        # if no distance, just do the infinite-distance/no parallax calculation
        usrepr = icrs_coo.represent_as(UnitSphericalRepresentation)
        i_ra = usrepr.lon.to(u.radian).value
        i_dec = usrepr.lat.to(u.radian).value
        cirs_ra, cirs_dec = atciqz(i_ra, i_dec, astrom)

        newrep = UnitSphericalRepresentation(lat=u.Quantity(cirs_dec, u.radian, copy=False),
                                             lon=u.Quantity(cirs_ra, u.radian, copy=False),
                                             copy=False)
    else:
        # When there is a distance,  we first offset for parallax to get the
        # astrometric coordinate direction and *then* run the ERFA transform for
        # no parallax/PM. This ensures reversibility and is more sensible for
        # inside solar system objects
        newxyz = icrs_coo.cartesian.xyz
        newxyz = np.rollaxis(newxyz, 0, newxyz.ndim) - astrom['eb'] * u.au
        # roll xyz back to the first axis
        newxyz = np.rollaxis(newxyz, -1, 0)
        newcart = CartesianRepresentation(newxyz)

        srepr = newcart.represent_as(SphericalRepresentation)
        i_ra = srepr.lon.to(u.radian).value
        i_dec = srepr.lat.to(u.radian).value
        cirs_ra, cirs_dec = atciqz(i_ra, i_dec, astrom)

        newrep = SphericalRepresentation(lat=u.Quantity(cirs_dec, u.radian, copy=False),
                                         lon=u.Quantity(cirs_ra, u.radian, copy=False),
                                         distance=srepr.distance, copy=False)

    return cirs_frame.realize_frame(newrep)


@frame_transform_graph.transform(FunctionTransform, CIRS, ICRS)
def cirs_to_icrs(cirs_coo, icrs_frame):
    srepr = cirs_coo.represent_as(UnitSphericalRepresentation)
    cirs_ra = srepr.lon.to(u.radian).value
    cirs_dec = srepr.lat.to(u.radian).value

    # set up the astrometry context for ICRS<->cirs and then convert to
    # astrometric coordinate direction
    jd1, jd2 = get_jd12(cirs_coo.obstime, 'tdb')
    x, y, s = get_cip(jd1, jd2)
    earth_pv, earth_heliocentric = prepare_earth_position_vel(cirs_coo.obstime)
    astrom = erfa.apci(jd1, jd2, earth_pv, earth_heliocentric, x, y, s)
    i_ra, i_dec = aticq(cirs_ra, cirs_dec, astrom)

    if cirs_coo.data.get_name() == 'unitspherical' or cirs_coo.data.to_cartesian().x.unit == u.one:
        # if no distance, just use the coordinate direction to yield the
        # infinite-distance/no parallax answer
        newrep = UnitSphericalRepresentation(lat=u.Quantity(i_dec, u.radian, copy=False),
                                             lon=u.Quantity(i_ra, u.radian, copy=False),
                                             copy=False)
    else:
        # When there is a distance, apply the parallax/offset to the SSB as the
        # last step - ensures round-tripping with the icrs_to_cirs transform

        # the distance in intermedrep is *not* a real distance as it does not
        # include the offset back to the SSB
        intermedrep = SphericalRepresentation(lat=u.Quantity(i_dec, u.radian, copy=False),
                                              lon=u.Quantity(i_ra, u.radian, copy=False),
                                              distance=cirs_coo.distance,
                                              copy=False)

        newxyz = intermedrep.to_cartesian().xyz
        # roll xyz to last axis and add the barycentre position
        newxyz = np.rollaxis(newxyz, 0, newxyz.ndim) + astrom['eb'] * u.au
        # roll xyz back to the first axis
        newxyz = np.rollaxis(newxyz, -1, 0)
        newrep = CartesianRepresentation(newxyz).represent_as(SphericalRepresentation)

    return icrs_frame.realize_frame(newrep)


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

    # first set up the astrometry context for ICRS<->GCRS
    pv = np.array([gcrs_frame.obsgeoloc.value,
                   gcrs_frame.obsgeovel.value])

    jd1, jd2 = get_jd12(gcrs_frame.obstime, 'tdb')
    earth_pv, earth_heliocentric = prepare_earth_position_vel(gcrs_frame.obstime)

    # get astrometry context
    astrom = erfa.apcs(jd1, jd2, pv, earth_pv, earth_heliocentric)

    if icrs_coo.data.get_name() == 'unitspherical' or icrs_coo.data.to_cartesian().x.unit == u.one:
        # if no distance, just do the infinite-distance/no parallax calculation
        usrepr = icrs_coo.represent_as(UnitSphericalRepresentation)
        i_ra = usrepr.lon.to(u.radian).value
        i_dec = usrepr.lat.to(u.radian).value
        gcrs_ra, gcrs_dec = atciqz(i_ra, i_dec, astrom)

        newrep = UnitSphericalRepresentation(lat=u.Quantity(gcrs_dec, u.radian, copy=False),
                                             lon=u.Quantity(gcrs_ra, u.radian, copy=False),
                                             copy=False)
    else:
        # When there is a distance,  we first offset for parallax to get the
        # BCRS coordinate direction and *then* run the ERFA transform for no
        # parallax/PM. This ensures reversibility and is more sensible for
        # inside solar system objects
        newxyz = icrs_coo.cartesian.xyz
        newxyz = np.rollaxis(newxyz, 0, newxyz.ndim) - astrom['eb'] * u.au
        newxyz = np.rollaxis(newxyz, -1, 0)
        newcart = CartesianRepresentation(newxyz)

        srepr = newcart.represent_as(SphericalRepresentation)
        i_ra = srepr.lon.to(u.radian).value
        i_dec = srepr.lat.to(u.radian).value
        gcrs_ra, gcrs_dec = atciqz(i_ra, i_dec, astrom)

        newrep = SphericalRepresentation(lat=u.Quantity(gcrs_dec, u.radian, copy=False),
                                         lon=u.Quantity(gcrs_ra, u.radian, copy=False),
                                         distance=srepr.distance, copy=False)

    return gcrs_frame.realize_frame(newrep)


@frame_transform_graph.transform(FunctionTransform, GCRS, ICRS)
def gcrs_to_icrs(gcrs_coo, icrs_frame):
    srepr = gcrs_coo.represent_as(UnitSphericalRepresentation)
    gcrs_ra = srepr.lon.to(u.radian).value
    gcrs_dec = srepr.lat.to(u.radian).value

    # set up the astrometry context for ICRS<->GCRS and then convert to BCRS
    # coordinate direction
    pv = np.array([gcrs_coo.obsgeoloc.value,
                   gcrs_coo.obsgeovel.value])
    jd1, jd2 = get_jd12(gcrs_coo.obstime, 'tdb')

    earth_pv, earth_heliocentric = prepare_earth_position_vel(gcrs_coo.obstime)
    astrom = erfa.apcs(jd1, jd2, pv, earth_pv, earth_heliocentric)

    i_ra, i_dec = aticq(gcrs_ra, gcrs_dec, astrom)

    if gcrs_coo.data.get_name() == 'unitspherical' or gcrs_coo.data.to_cartesian().x.unit == u.one:
        # if no distance, just use the coordinate direction to yield the
        # infinite-distance/no parallax answer
        newrep = UnitSphericalRepresentation(lat=u.Quantity(i_dec, u.radian, copy=False),
                                             lon=u.Quantity(i_ra, u.radian, copy=False),
                                             copy=False)
    else:
        # When there is a distance, apply the parallax/offset to the SSB as the
        # last step - ensures round-tripping with the icrs_to_gcrs transform

        # the distance in intermedrep is *not* a real distance as it does not
        # include the offset back to the SSB
        intermedrep = SphericalRepresentation(lat=u.Quantity(i_dec, u.radian, copy=False),
                                              lon=u.Quantity(i_ra, u.radian, copy=False),
                                              distance=gcrs_coo.distance,
                                              copy=False)

        newxyz = intermedrep.to_cartesian().xyz
        # roll xyz to last axis and add the heliocentre position
        newxyz = np.rollaxis(newxyz, 0, newxyz.ndim) + astrom['eb'] * u.au
        # roll xyz back to the first axis
        newxyz = np.rollaxis(newxyz, -1, 0)
        newrep = CartesianRepresentation(newxyz).represent_as(SphericalRepresentation)
    return icrs_frame.realize_frame(newrep)


@frame_transform_graph.transform(FunctionTransform, GCRS, GCRS)
def gcrs_to_gcrs(from_coo, to_frame):
    if (np.all(from_coo.obstime == to_frame.obstime)
        and np.all(from_coo.obsgeoloc == to_frame.obsgeoloc)):
        return to_frame.realize_frame(from_coo.data)
    else:
        # like CIRS, we do this self-transform via ICRS
        return from_coo.transform_to(ICRS).transform_to(to_frame)


@frame_transform_graph.transform(FunctionTransform, GCRS, HCRS)
def gcrs_to_hcrs(gcrs_coo, hcrs_frame):

    if np.any(gcrs_coo.obstime != hcrs_frame.obstime):
        # if they GCRS obstime and HCRS obstime are not the same, we first
        # have to move to a GCRS where they are.
        frameattrs = gcrs_coo.get_frame_attr_names()
        frameattrs['obstime'] = hcrs_frame.obstime
        gcrs_coo = gcrs_coo.transform_to(GCRS(**frameattrs))

    srepr = gcrs_coo.represent_as(UnitSphericalRepresentation)
    gcrs_ra = srepr.lon.to(u.radian).value
    gcrs_dec = srepr.lat.to(u.radian).value

    # set up the astrometry context for ICRS<->GCRS and then convert to ICRS
    # coordinate direction
    pv = np.array([gcrs_coo.obsgeoloc.value,
                   gcrs_coo.obsgeovel.value])

    jd1, jd2 = get_jd12(hcrs_frame.obstime, 'tdb')
    earth_pv, earth_heliocentric = prepare_earth_position_vel(gcrs_coo.obstime)
    astrom = erfa.apcs(jd1, jd2, pv, earth_pv, earth_heliocentric)

    i_ra, i_dec = aticq(gcrs_ra, gcrs_dec, astrom)

    # convert to Quantity objects
    i_ra = u.Quantity(i_ra, u.radian, copy=False)
    i_dec = u.Quantity(i_dec, u.radian, copy=False)
    if gcrs_coo.data.get_name() == 'unitspherical' or gcrs_coo.data.to_cartesian().x.unit == u.one:
        # if no distance, just use the coordinate direction to yield the
        # infinite-distance/no parallax answer
        newrep = UnitSphericalRepresentation(lat=i_dec, lon=i_ra, copy=False)
    else:
        # When there is a distance, apply the parallax/offset to the
        # Heliocentre as the last step to ensure round-tripping with the
        # hcrs_to_gcrs transform

        # Note that the distance in intermedrep is *not* a real distance as it
        # does not include the offset back to the Heliocentre
        intermedrep = SphericalRepresentation(lat=i_dec, lon=i_ra,
                                              distance=gcrs_coo.distance,
                                              copy=False)

        newxyz = intermedrep.to_cartesian().xyz
        # roll astrom['eh'] to the last axis and scale
        eh = np.rollaxis(astrom['eh'], 0, astrom['eh'].ndim) * astrom['em'] * u.au
        # roll astrom['eh'] back to the first axis
        eh = np.rollaxis(eh, -1, 0)

        # roll xyz to last axis and add the heliocentre position
        newxyz = np.rollaxis(newxyz, 0, newxyz.ndim) + eh
        # roll xyz back to the first axis
        newxyz = np.rollaxis(newxyz, -1, 0)
        newrep = CartesianRepresentation(newxyz).represent_as(SphericalRepresentation)

    return hcrs_frame.realize_frame(newrep)


@frame_transform_graph.transform(FunctionTransform, HCRS, ICRS)
def hcrs_to_icrs(hcrs_coo, icrs_frame):
    # this goes here to avoid circular import errors
    from ..solar_system import get_body_barycentric, solar_system_ephemeris

    # this is just an origin translation so without a distance it
    # does not change the data
    if hcrs_coo.data.__class__ == UnitSphericalRepresentation:
        return icrs_frame.realize_frame(hcrs_coo.data)

    ephemeris = solar_system_ephemeris.get()
    bary_sun_pos = get_body_barycentric('sun', hcrs_coo.obstime, ephemeris=ephemeris)
    sun_vector = bary_sun_pos.xyz
    new_vector = np.rollaxis(hcrs_coo.cartesian.xyz, -1, 0) + np.rollaxis(sun_vector, -1, 0)
    newrep = CartesianRepresentation(np.rollaxis(new_vector, 0, new_vector.ndim))
    return icrs_frame.realize_frame(newrep)


@frame_transform_graph.transform(FunctionTransform, ICRS, HCRS)
def icrs_to_hcrs(icrs_coo, hcrs_frame):
    # this goes here to avoid circular import errors
    from ..solar_system import get_body_barycentric, solar_system_ephemeris

    # this is just an origin translation so without a distance it
    # does not change the data
    if icrs_coo.data.__class__ == UnitSphericalRepresentation:
        return hcrs_frame.realize_frame(icrs_coo.data)

    ephemeris = solar_system_ephemeris.get()
    bary_sun_pos = get_body_barycentric('sun', hcrs_frame.obstime, ephemeris=ephemeris)
    sun_vector = bary_sun_pos.xyz
    new_vector = np.rollaxis(icrs_coo.cartesian.xyz, -1, 0) - np.rollaxis(sun_vector, -1, 0)
    newrep = CartesianRepresentation(np.rollaxis(new_vector, 0, new_vector.ndim))
    return hcrs_frame.realize_frame(newrep)


@frame_transform_graph.transform(FunctionTransform, HCRS, HCRS)
def hcrs_to_hcrs(from_coo, to_frame):
    if np.all(from_coo.obstime == to_frame.obstime):
        return to_frame.realize_frame(from_coo.data)
    else:
        # like CIRS, we do this self-transform via ICRS
        return from_coo.transform_to(ICRS).transform_to(to_frame)
