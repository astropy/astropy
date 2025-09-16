# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Contains the transformation functions for getting to "observed" systems from CIRS.
"""

import erfa
import numpy as np

from astropy import units as u
from astropy.coordinates.baseframe import frame_transform_graph
from astropy.coordinates.erfa_astrom import erfa_astrom
from astropy.coordinates.representation import (
    SphericalRepresentation,
    UnitSphericalRepresentation,
)
from astropy.coordinates.transformations import FunctionTransformWithFiniteDifference
from astropy.utils.compat import COPY_IF_NEEDED

from .altaz import AltAz
from .cirs import CIRS
from .hadec import HADec
from .utils import PIOVER2


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference, CIRS, AltAz)
@frame_transform_graph.transform(FunctionTransformWithFiniteDifference, CIRS, HADec)
def cirs_to_observed(cirs_coo, observed_frame):
    if np.any(observed_frame.location != cirs_coo.location) or np.any(
        cirs_coo.obstime != observed_frame.obstime
    ):
        cirs_coo = cirs_coo.transform_to(
            CIRS(obstime=observed_frame.obstime, location=observed_frame.location)
        )

    # if the data are UnitSphericalRepresentation, we can skip the distance calculations
    is_unitspherical = (
        isinstance(cirs_coo.data, UnitSphericalRepresentation)
        or cirs_coo.cartesian.x.unit == u.one
    )

    # We used to do "astrometric" corrections here, but these are no longer necessary
    # CIRS has proper topocentric behaviour
    usrepr = cirs_coo.represent_as(UnitSphericalRepresentation)
    cirs_ra = usrepr.lon.to_value(u.radian)
    cirs_dec = usrepr.lat.to_value(u.radian)
    # first set up the astrometry context for CIRS<->observed
    astrom = erfa_astrom.get().apio(observed_frame)

    if isinstance(observed_frame, AltAz):
        lon, zen, _, _, _ = erfa.atioq(cirs_ra, cirs_dec, astrom)
        lat = PIOVER2 - zen
    else:
        _, _, lon, lat, _ = erfa.atioq(cirs_ra, cirs_dec, astrom)

    if is_unitspherical:
        rep = UnitSphericalRepresentation(
            lat=u.Quantity(lat, u.radian, copy=COPY_IF_NEEDED),
            lon=u.Quantity(lon, u.radian, copy=COPY_IF_NEEDED),
            copy=False,
        )
    else:
        # since we've transformed to CIRS at the observatory location, just use CIRS distance
        rep = SphericalRepresentation(
            lat=u.Quantity(lat, u.radian, copy=COPY_IF_NEEDED),
            lon=u.Quantity(lon, u.radian, copy=COPY_IF_NEEDED),
            distance=cirs_coo.distance,
            copy=False,
        )
    return observed_frame.realize_frame(rep)


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference, AltAz, CIRS)
@frame_transform_graph.transform(FunctionTransformWithFiniteDifference, HADec, CIRS)
def observed_to_cirs(observed_coo, cirs_frame):
    usrepr = observed_coo.represent_as(UnitSphericalRepresentation)
    lon = usrepr.lon.to_value(u.radian)
    lat = usrepr.lat.to_value(u.radian)

    if isinstance(observed_coo, AltAz):
        # the 'A' indicates zen/az inputs
        coord_type = "A"
        lat = PIOVER2 - lat
    else:
        coord_type = "H"

    # first set up the astrometry context for ICRS<->CIRS at the observed_coo time
    astrom = erfa_astrom.get().apio(observed_coo)

    cirs_ra, cirs_dec = erfa.atoiq(coord_type, lon, lat, astrom) << u.radian
    if (
        isinstance(observed_coo.data, UnitSphericalRepresentation)
        or observed_coo.cartesian.x.unit == u.one
    ):
        distance = None
    else:
        distance = observed_coo.distance

    cirs_at_aa_time = CIRS(
        ra=cirs_ra,
        dec=cirs_dec,
        distance=distance,
        obstime=observed_coo.obstime,
        location=observed_coo.location,
    )

    # this final transform may be a no-op if the obstimes and locations are the same
    return cirs_at_aa_time.transform_to(cirs_frame)
