# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

import numpy as np

from ...extern.six.moves import range
from ... import units as u
from ..baseframe import frame_transform_graph
from ..attributes import TimeAttribute
from ..transformations import (FunctionTransformWithFiniteDifference,
                               FunctionTransform, DynamicMatrixTransform)
from ..representation import (CartesianRepresentation,
                              UnitSphericalRepresentation)
from .. import earth_orientation as earth

from .utils import EQUINOX_B1950
from .baseradec import _base_radec_docstring, BaseRADecFrame


class FK4(BaseRADecFrame):
    """
    A coordinate or frame in the FK4 system.

    Note that this is a barycentric version of FK4 - that is, the origin for
    this frame is the Solar System Barycenter, *not* the Earth geocenter.

    The frame attributes are listed under **Other Parameters**.

    {params}

    Other parameters
    ----------------
    equinox : `~astropy.time.Time`
        The equinox of this frame.
    obstime : `~astropy.time.Time`
        The time this frame was observed.  If ``None``, will be the same as
        ``equinox``.
    """

    equinox = TimeAttribute(default=EQUINOX_B1950)
    obstime = TimeAttribute(default=None, secondary_attribute='equinox')


FK4.__doc__ = FK4.__doc__.format(params=_base_radec_docstring)

# the "self" transform


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference, FK4, FK4)
def fk4_to_fk4(fk4coord1, fk4frame2):
    # deceptively complicated: need to transform to No E-terms FK4, precess, and
    # then come back, because precession is non-trivial with E-terms
    fnoe_w_eqx1 = fk4coord1.transform_to(FK4NoETerms(equinox=fk4coord1.equinox))
    fnoe_w_eqx2 = fnoe_w_eqx1.transform_to(FK4NoETerms(equinox=fk4frame2.equinox))
    return fnoe_w_eqx2.transform_to(fk4frame2)


class FK4NoETerms(BaseRADecFrame):
    """
    A coordinate or frame in the FK4 system, but with the E-terms of aberration
    removed.

    The frame attributes are listed under **Other Parameters**.

    {params}

    Other parameters
    ----------------
    equinox : `~astropy.time.Time`
        The equinox of this frame.
    obstime : `~astropy.time.Time`
        The time this frame was observed.  If ``None``, will be the same as
        ``equinox``.
    """

    equinox = TimeAttribute(default=EQUINOX_B1950)
    obstime = TimeAttribute(default=None, secondary_attribute='equinox')

    @staticmethod
    def _precession_matrix(oldequinox, newequinox):
        """
        Compute and return the precession matrix for FK4 using Newcomb's method.
        Used inside some of the transformation functions.

        Parameters
        ----------
        oldequinox : `~astropy.time.Time`
            The equinox to precess from.
        newequinox : `~astropy.time.Time`
            The equinox to precess to.

        Returns
        -------
        newcoord : array
            The precession matrix to transform to the new equinox
        """
        return earth._precession_matrix_besselian(oldequinox.byear, newequinox.byear)


FK4NoETerms.__doc__ = FK4NoETerms.__doc__.format(params=_base_radec_docstring)

# the "self" transform


@frame_transform_graph.transform(DynamicMatrixTransform, FK4NoETerms, FK4NoETerms)
def fk4noe_to_fk4noe(fk4necoord1, fk4neframe2):
    return fk4necoord1._precession_matrix(fk4necoord1.equinox, fk4neframe2.equinox)


# FK4-NO-E to/from FK4 ----------------------------->
# Unlike other frames, this module include *two* frame classes for FK4
# coordinates - one including the E-terms of aberration (FK4), and
# one not including them (FK4NoETerms). The following functions
# implement the transformation between these two.
def fk4_e_terms(equinox):
    """
    Return the e-terms of aberation vector

    Parameters
    ----------
    equinox : Time object
        The equinox for which to compute the e-terms
    """
    # Constant of aberration at J2000; from Explanatory Supplement to the
    # Astronomical Almanac (Seidelmann, 2005).
    k = 0.0056932  # in degrees (v_earth/c ~ 1e-4 rad ~ 0.0057 deg)
    k = np.radians(k)

    # Eccentricity of the Earth's orbit
    e = earth.eccentricity(equinox.jd)

    # Mean longitude of perigee of the solar orbit
    g = earth.mean_lon_of_perigee(equinox.jd)
    g = np.radians(g)

    # Obliquity of the ecliptic
    o = earth.obliquity(equinox.jd, algorithm=1980)
    o = np.radians(o)

    return e * k * np.sin(g), \
           -e * k * np.cos(g) * np.cos(o), \
           -e * k * np.cos(g) * np.sin(o)


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference, FK4, FK4NoETerms)
def fk4_to_fk4_no_e(fk4coord, fk4noeframe):
    # Extract cartesian vector
    rep = fk4coord.cartesian

    # Find distance (for re-normalization)
    d_orig = rep.norm()
    rep /= d_orig

    # Apply E-terms of aberration. Note that this depends on the equinox (not
    # the observing time/epoch) of the coordinates. See issue #1496 for a
    # discussion of this.
    eterms_a = CartesianRepresentation(
        u.Quantity(fk4_e_terms(fk4coord.equinox), u.dimensionless_unscaled,
                   copy=False), copy=False)
    rep = rep - eterms_a + eterms_a.dot(rep) * rep

    # Find new distance (for re-normalization)
    d_new = rep.norm()

    # Renormalize
    rep *= d_orig / d_new

    # now re-cast into an appropriate Representation, and precess if need be
    if isinstance(fk4coord.data, UnitSphericalRepresentation):
        rep = rep.represent_as(UnitSphericalRepresentation)

    # if no obstime was given in the new frame, use the old one for consistency
    newobstime = fk4coord._obstime if fk4noeframe._obstime is None else fk4noeframe._obstime

    fk4noe = FK4NoETerms(rep, equinox=fk4coord.equinox, obstime=newobstime)
    if fk4coord.equinox != fk4noeframe.equinox:
        # precession
        fk4noe = fk4noe.transform_to(fk4noeframe)
    return fk4noe


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference, FK4NoETerms, FK4)
def fk4_no_e_to_fk4(fk4noecoord, fk4frame):
    # first precess, if necessary
    if fk4noecoord.equinox != fk4frame.equinox:
        fk4noe_w_fk4equinox = FK4NoETerms(equinox=fk4frame.equinox,
                                          obstime=fk4noecoord.obstime)
        fk4noecoord = fk4noecoord.transform_to(fk4noe_w_fk4equinox)

    # Extract cartesian vector
    rep = fk4noecoord.cartesian

    # Find distance (for re-normalization)
    d_orig = rep.norm()
    rep /= d_orig

    # Apply E-terms of aberration. Note that this depends on the equinox (not
    # the observing time/epoch) of the coordinates. See issue #1496 for a
    # discussion of this.
    eterms_a = CartesianRepresentation(
        u.Quantity(fk4_e_terms(fk4noecoord.equinox), u.dimensionless_unscaled,
                   copy=False), copy=False)

    rep0 = rep.copy()
    for _ in range(10):
        rep = (eterms_a + rep0) / (1. + eterms_a.dot(rep))

    # Find new distance (for re-normalization)
    d_new = rep.norm()

    # Renormalize
    rep *= d_orig / d_new

    # now re-cast into an appropriate Representation, and precess if need be
    if isinstance(fk4noecoord.data, UnitSphericalRepresentation):
        rep = rep.represent_as(UnitSphericalRepresentation)

    return fk4frame.realize_frame(rep)
