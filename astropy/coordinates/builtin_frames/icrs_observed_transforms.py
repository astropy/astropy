# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Contains the transformation functions for getting to "observed" systems from ICRS.
Currently that just means AltAz.
"""
import erfa

from astropy import units as u
from astropy.coordinates.builtin_frames.utils import atciqz, aticq
from astropy.coordinates.baseframe import frame_transform_graph
from astropy.coordinates.transformations import FunctionTransformWithFiniteDifference
from astropy.coordinates.representation import (SphericalRepresentation,
                                                CartesianRepresentation,
                                                UnitSphericalRepresentation)

from .icrs import ICRS
from .altaz import AltAz
from .utils import PIOVER2
from ..erfa_astrom import erfa_astrom


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference, ICRS, AltAz)
def icrs_to_altaz(icrs_coo, altaz_frame):
    # if the data are UnitSphericalRepresentation, we can skip the distance calculations
    is_unitspherical = (isinstance(icrs_coo.data, UnitSphericalRepresentation) or
                        icrs_coo.cartesian.x.unit == u.one)
    # first set up the astrometry context for ICRS<->AltAz
    astrom = erfa_astrom.get().apco(altaz_frame)

    # correct for parallax to find BCRS direction from observer (as in erfa.pmpx)
    if is_unitspherical:
        srepr = icrs_coo.spherical
    else:
        observer_icrs = CartesianRepresentation(astrom['eb'], unit=u.au, xyz_axis=-1, copy=False)
        srepr = (icrs_coo.cartesian - observer_icrs).represent_as(
            SphericalRepresentation)

    # convert to topocentric CIRS
    cirs_ra, cirs_dec = atciqz(srepr, astrom)

    # now perform AltAz conversion
    az, zen, ha, odec, ora = erfa.atioq(cirs_ra, cirs_dec, astrom)
    alt = PIOVER2-zen
    if is_unitspherical:
        aa_srepr = UnitSphericalRepresentation(az << u.radian, alt << u.radian, copy=False)
    else:
        aa_srepr = SphericalRepresentation(az << u.radian, alt << u.radian, srepr.distance, copy=False)
    return altaz_frame.realize_frame(aa_srepr)


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference, AltAz, ICRS)
def altaz_to_icrs(altaz_coo, icrs_frame):
    # if the data are UnitSphericalRepresentation, we can skip the distance calculations
    is_unitspherical = (isinstance(altaz_coo.data, UnitSphericalRepresentation) or
                        altaz_coo.cartesian.x.unit == u.one)

    usrepr = altaz_coo.represent_as(UnitSphericalRepresentation)
    az = usrepr.lon.to_value(u.radian)
    zen = PIOVER2 - usrepr.lat.to_value(u.radian)

    # first set up the astrometry context for ICRS<->CIRS at the altaz_coo time
    astrom = erfa_astrom.get().apco(altaz_coo)

    # Topocentric CIRS
    cirs_ra, cirs_dec = erfa.atoiq('A', az, zen, astrom) << u.radian
    if is_unitspherical:
        srepr = SphericalRepresentation(cirs_ra, cirs_dec, 1, copy=False)
    else:
        srepr = SphericalRepresentation(lon=cirs_ra, lat=cirs_dec,
                                        distance=altaz_coo.distance, copy=False)

    # BCRS (Astrometric) direction to source
    bcrs_ra, bcrs_dec = aticq(srepr, astrom) << u.radian

    # Correct for parallax to get ICRS representation
    if is_unitspherical:
        icrs_srepr = UnitSphericalRepresentation(bcrs_ra, bcrs_dec, copy=False)
    else:
        icrs_srepr = SphericalRepresentation(lon=bcrs_ra, lat=bcrs_dec,
                                             distance=altaz_coo.distance, copy=False)
        observer_icrs = CartesianRepresentation(astrom['eb'], unit=u.au, xyz_axis=-1, copy=False)
        newrepr = icrs_srepr.to_cartesian() + observer_icrs
        icrs_srepr = newrepr.represent_as(SphericalRepresentation)

    return icrs_frame.realize_frame(icrs_srepr)
