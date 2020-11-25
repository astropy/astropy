# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Contains the transformation functions for getting from ICRS/HCRS to CIRS and
anything in between (currently that means GCRS)
"""

import numpy as np

from astropy import units as u
from astropy.coordinates.baseframe import frame_transform_graph
from astropy.coordinates.transformations import (
    FunctionTransformWithFiniteDifference,
    AffineTransform,
)
from astropy.coordinates.representation import (
    SphericalRepresentation,
    CartesianRepresentation,
    UnitSphericalRepresentation,
    CartesianDifferential,
)

from .icrs import ICRS
from .gcrs import GCRS
from .cirs import CIRS
from .hcrs import HCRS
from .equatorial import TETE
from .utils import aticq, atciqz, get_offset_sun_from_barycenter

from ..erfa_astrom import erfa_astrom


# First the ICRS/CIRS related transforms
@frame_transform_graph.transform(FunctionTransformWithFiniteDifference, ICRS, CIRS)
def icrs_to_cirs(icrs_coo, cirs_frame):
    # first set up the astrometry context for ICRS<->CIRS
    astrom = erfa_astrom.get().apci(cirs_frame)

    if icrs_coo.data.get_name() == 'unitspherical' or icrs_coo.data.to_cartesian().x.unit == u.one:
        # if no distance, just do the infinite-distance/no parallax calculation
        srepr = icrs_coo.spherical
        cirs_ra, cirs_dec = atciqz(srepr.without_differentials(), astrom)

        newrep = UnitSphericalRepresentation(lat=u.Quantity(cirs_dec, u.radian, copy=False),
                                             lon=u.Quantity(cirs_ra, u.radian, copy=False),
                                             copy=False)
    else:
        # When there is a distance,  we first offset for parallax to get the
        # astrometric coordinate direction and *then* run the ERFA transform for
        # no parallax/PM. This ensures reversibility and is more sensible for
        # inside solar system objects
        astrom_eb = CartesianRepresentation(astrom['eb'], unit=u.au,
                                            xyz_axis=-1, copy=False)
        newcart = icrs_coo.cartesian - astrom_eb
        srepr = newcart.represent_as(SphericalRepresentation)
        cirs_ra, cirs_dec = atciqz(srepr.without_differentials(), astrom)

        newrep = SphericalRepresentation(lat=u.Quantity(cirs_dec, u.radian, copy=False),
                                         lon=u.Quantity(cirs_ra, u.radian, copy=False),
                                         distance=srepr.distance, copy=False)

    return cirs_frame.realize_frame(newrep)


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference, CIRS, ICRS)
def cirs_to_icrs(cirs_coo, icrs_frame):
    # set up the astrometry context for ICRS<->cirs and then convert to
    # astrometric coordinate direction
    astrom = erfa_astrom.get().apci(cirs_coo)
    srepr = cirs_coo.represent_as(SphericalRepresentation)
    i_ra, i_dec = aticq(srepr.without_differentials(), astrom)

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
                                              distance=srepr.distance,
                                              copy=False)

        astrom_eb = CartesianRepresentation(astrom['eb'], unit=u.au,
                                            xyz_axis=-1, copy=False)
        newrep = intermedrep + astrom_eb

    return icrs_frame.realize_frame(newrep)


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference, CIRS, CIRS)
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
        return from_coo.transform_to(ICRS()).transform_to(to_frame)


# Now the GCRS-related transforms to/from ICRS

@frame_transform_graph.transform(FunctionTransformWithFiniteDifference, ICRS, GCRS)
def icrs_to_gcrs(icrs_coo, gcrs_frame):
    # first set up the astrometry context for ICRS<->GCRS.
    astrom = erfa_astrom.get().apcs(gcrs_frame)

    if icrs_coo.data.get_name() == 'unitspherical' or icrs_coo.data.to_cartesian().x.unit == u.one:
        # if no distance, just do the infinite-distance/no parallax calculation
        srepr = icrs_coo.represent_as(SphericalRepresentation)
        gcrs_ra, gcrs_dec = atciqz(srepr.without_differentials(), astrom)

        newrep = UnitSphericalRepresentation(lat=u.Quantity(gcrs_dec, u.radian, copy=False),
                                             lon=u.Quantity(gcrs_ra, u.radian, copy=False),
                                             copy=False)
    else:
        # When there is a distance,  we first offset for parallax to get the
        # BCRS coordinate direction and *then* run the ERFA transform for no
        # parallax/PM. This ensures reversibility and is more sensible for
        # inside solar system objects
        astrom_eb = CartesianRepresentation(astrom['eb'], unit=u.au,
                                            xyz_axis=-1, copy=False)
        newcart = icrs_coo.cartesian - astrom_eb

        srepr = newcart.represent_as(SphericalRepresentation)
        gcrs_ra, gcrs_dec = atciqz(srepr.without_differentials(), astrom)

        newrep = SphericalRepresentation(lat=u.Quantity(gcrs_dec, u.radian, copy=False),
                                         lon=u.Quantity(gcrs_ra, u.radian, copy=False),
                                         distance=srepr.distance, copy=False)

    return gcrs_frame.realize_frame(newrep)


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference,
                                 GCRS, ICRS)
def gcrs_to_icrs(gcrs_coo, icrs_frame):
    # set up the astrometry context for ICRS<->GCRS and then convert to BCRS
    # coordinate direction
    astrom = erfa_astrom.get().apcs(gcrs_coo)

    srepr = gcrs_coo.represent_as(SphericalRepresentation)
    i_ra, i_dec = aticq(srepr.without_differentials(), astrom)

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
                                              distance=srepr.distance,
                                              copy=False)

        astrom_eb = CartesianRepresentation(astrom['eb'], unit=u.au,
                                            xyz_axis=-1, copy=False)
        newrep = intermedrep + astrom_eb

    return icrs_frame.realize_frame(newrep)


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference, GCRS, GCRS)
def gcrs_to_gcrs(from_coo, to_frame):
    if (np.all(from_coo.obstime == to_frame.obstime) and
            np.all(from_coo.obsgeoloc == to_frame.obsgeoloc)):
        return to_frame.realize_frame(from_coo.data)
    else:
        # like CIRS, we do this self-transform via ICRS
        return from_coo.transform_to(ICRS()).transform_to(to_frame)


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference, GCRS, HCRS)
def gcrs_to_hcrs(gcrs_coo, hcrs_frame):

    if np.any(gcrs_coo.obstime != hcrs_frame.obstime):
        # if they GCRS obstime and HCRS obstime are not the same, we first
        # have to move to a GCRS where they are.
        frameattrs = gcrs_coo.get_frame_attr_names()
        frameattrs['obstime'] = hcrs_frame.obstime
        gcrs_coo = gcrs_coo.transform_to(GCRS(**frameattrs))

    # set up the astrometry context for ICRS<->GCRS and then convert to ICRS
    # coordinate direction
    astrom = erfa_astrom.get().apcs(gcrs_coo)
    srepr = gcrs_coo.represent_as(SphericalRepresentation)
    i_ra, i_dec = aticq(srepr.without_differentials(), astrom)

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
                                              distance=srepr.distance,
                                              copy=False)

        # astrom['eh'] and astrom['em'] contain Sun to observer unit vector,
        # and distance, respectively. Shapes are (X) and (X,3), where (X) is the
        # shape resulting from broadcasting the shape of the times object
        # against the shape of the pv array.
        # broadcast em to eh and scale eh
        eh = astrom['eh'] * astrom['em'][..., np.newaxis]
        eh = CartesianRepresentation(eh, unit=u.au, xyz_axis=-1, copy=False)

        newrep = intermedrep.to_cartesian() + eh

    return hcrs_frame.realize_frame(newrep)


_NEED_ORIGIN_HINT = ("The input {0} coordinates do not have length units. This "
                     "probably means you created coordinates with lat/lon but "
                     "no distance.  Heliocentric<->ICRS transforms cannot "
                     "function in this case because there is an origin shift.")


@frame_transform_graph.transform(AffineTransform, HCRS, ICRS)
def hcrs_to_icrs(hcrs_coo, icrs_frame):
    # this is just an origin translation so without a distance it cannot go ahead
    if isinstance(hcrs_coo.data, UnitSphericalRepresentation):
        raise u.UnitsError(_NEED_ORIGIN_HINT.format(hcrs_coo.__class__.__name__))

    return None, get_offset_sun_from_barycenter(hcrs_coo.obstime,
                                                include_velocity=bool(hcrs_coo.data.differentials))


@frame_transform_graph.transform(AffineTransform, ICRS, HCRS)
def icrs_to_hcrs(icrs_coo, hcrs_frame):
    # this is just an origin translation so without a distance it cannot go ahead
    if isinstance(icrs_coo.data, UnitSphericalRepresentation):
        raise u.UnitsError(_NEED_ORIGIN_HINT.format(icrs_coo.__class__.__name__))

    return None, get_offset_sun_from_barycenter(hcrs_frame.obstime, reverse=True,
                                                include_velocity=bool(icrs_coo.data.differentials))


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference, HCRS, HCRS)
def hcrs_to_hcrs(from_coo, to_frame):
    if np.all(from_coo.obstime == to_frame.obstime):
        return to_frame.realize_frame(from_coo.data)
    else:
        # like CIRS, we do this self-transform via ICRS
        return from_coo.transform_to(ICRS()).transform_to(to_frame)


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference, TETE, TETE)
def tete_to_tete(from_coo, to_frame):
    if (np.all(from_coo.location == to_frame.location)
            and np.all(from_coo.obstime == to_frame.obstime)):
        return to_frame.realize_frame(from_coo.data)
    else:
        # We perform the self-transform through ICRS, since any change in obstime
        # will require handling the change in pos/vel of the Earth and it's impact
        # on aberration.
        return from_coo.transform_to(ICRS()).transform_to(to_frame)
