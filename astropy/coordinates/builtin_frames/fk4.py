# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

import numpy as np

from ... import units as u
from ..representation import (SphericalRepresentation, CartesianRepresentation,
                              UnitSphericalRepresentation)
from ..baseframe import (BaseCoordinateFrame, frame_transform_graph,
                         TimeFrameAttribute, RepresentationMapping)
from ..transformations import FunctionTransform, DynamicMatrixTransform
from .. import earth_orientation as earth

from .fk5 import FK5
from .consts import EQUINOX_B1950, EQUINOX_J2000


class FK4(BaseCoordinateFrame):
    """
    A coordinate or frame in the FK4 system.

    Parameters
    ----------
    representation : `BaseRepresentation` or None
        A representation object or None to have no data (or use the other keywords)
    ra : `Angle`, optional, must be keyword
        The RA for this object (``dec`` must also be given and ``representation``
        must be None).
    dec : `Angle`, optional, must be keyword
        The Declination for this object (``ra`` must also be given and
        ``representation`` must be None).
    distance : :class:`~astropy.units.Quantity`, optional, must be keyword
        The Distance for this object along the line-of-sight.
        (``representation`` must be None).
    equinox : astropy.time.Time, optional, must be keyword
        The equinox of this frame.
    obstime : astropy.time.Time, optional, must be keyword
        The time this frame was observed.  If None, will be the same as
        ``equinox``.
    """
    frame_specific_representation_info = {
        'spherical': [RepresentationMapping('lon', 'ra'),
                      RepresentationMapping('lat', 'dec')]
    }
    frame_specific_representation_info['unitspherical'] = \
        frame_specific_representation_info['spherical']

    default_representation = SphericalRepresentation
    equinox = TimeFrameAttribute(default=EQUINOX_B1950)
    obstime = TimeFrameAttribute(default=None, secondary_attribute='equinox')


@frame_transform_graph.transform(FunctionTransform, FK4, FK4)
def fk4_to_fk4(fk4coord1, fk4frame2):
    # deceptively complicated: need to transform to No E-terms FK4, precess, and
    # then come back, because precession is non-trivial with E-terms
    fnoe_w_eqx1 = fk4coord1.transform_to(FK4NoETerms(equinox=fk4coord1.equinox))
    fnoe_w_eqx2 = fnoe_w_eqx1.transform_to(FK4NoETerms(equinox=fk4frame2.equinox))
    return fnoe_w_eqx2.transform_to(fk4frame2)


class FK4NoETerms(BaseCoordinateFrame):
    """
    A coordinate or frame in the FK4 system, but with the E-terms of aberration
    removed.

    Parameters
    ----------
    representation : `BaseRepresentation` or None
        A representation object or None to have no data (or use the other keywords)
    ra : `Angle`, optional, must be keyword
        The RA for this object (``dec`` must also be given and ``representation``
        must be None).
    dec : `Angle`, optional, must be keyword
        The Declination for this object (``ra`` must also be given and
        ``representation`` must be None).
    distance : :class:`~astropy.units.Quantity`, optional, must be keyword
        The Distance for this object along the line-of-sight.
        (``representation`` must be None).
    obstime : astropy.time.Time, optional, must be keyword
        The time this frame was observed.  If None, will be the same as
        ``equinox``.
    """
    frame_specific_representation_info = {
        'spherical': [RepresentationMapping('lon', 'ra'),
                      RepresentationMapping('lat', 'dec')]
    }
    frame_specific_representation_info['unitspherical'] = \
        frame_specific_representation_info['spherical']

    default_representation = SphericalRepresentation
    equinox = TimeFrameAttribute(default=EQUINOX_B1950)
    obstime = TimeFrameAttribute(default=None, secondary_attribute='equinox')

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

    @staticmethod
    def _fk4_B_matrix(obstime):
        """
        This is a correction term in the FK4 transformations because FK4 is a
        rotating system - see Murray 89 eqn 29
        """
        # Note this is *julian century*, not besselian
        T = (obstime.jyear - 1950.) / 100.
        return _B1950_TO_J2000_M + _FK4_CORR * T


@frame_transform_graph.transform(DynamicMatrixTransform, FK4NoETerms, FK4NoETerms)
def fk4noe_to_fk4noe(fk4necoord1, fk4neframe2):
    return fk4necoord1._precession_matrix(fk4necoord1.equinox, fk4neframe2.equinox)


# FK4-NO-E to/from FK4 ----------------------------->
# In the present framework, we include two coordinate classes for FK4
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
    # Constant of aberration at J2000
    k = 0.0056932

    # Eccentricity of the Earth's orbit
    e = earth.eccentricity(equinox.jd)
    e = np.radians(e)

    # Mean longitude of perigee of the solar orbit
    g = earth.mean_lon_of_perigee(equinox.jd)
    g = np.radians(g)

    # Obliquity of the ecliptic
    o = earth.obliquity(equinox.jd, algorithm=1980)
    o = np.radians(o)

    return e * k * np.sin(g), \
           -e * k * np.cos(g) * np.cos(o), \
           -e * k * np.cos(g) * np.sin(o)


@frame_transform_graph.transform(FunctionTransform, FK4, FK4NoETerms)
def fk4_to_fk4_no_e(fk4coord, fk4noeframe):
    # Extract cartesian vector
    c = fk4coord.cartesian.xyz
    r = np.asarray(c.reshape((3, c.size // 3)))

    # Find distance (for re-normalization)
    d_orig = np.sqrt(np.sum(r ** 2))

    # Apply E-terms of aberration. Note that this depends on the equinox (not
    # the observing time/epoch) of the coordinates. See issue #1496 for a
    # discussion of this.
    eterms_a = np.asarray(fk4_e_terms(fk4coord.equinox))
    r = r - eterms_a.reshape(3, 1) + np.dot(eterms_a, r) * r

    # Find new distance (for re-normalization)
    d_new = np.sqrt(np.sum(r ** 2))

    # Renormalize
    r = r * d_orig / d_new

    subshape = c.shape[1:]
    x = r[0].reshape(subshape)
    y = r[1].reshape(subshape)
    z = r[2].reshape(subshape)

    #now re-cast into an appropriate Representation, and precess if need be
    if isinstance(fk4coord.data, UnitSphericalRepresentation):
        representation = CartesianRepresentation(x=x*u.one, y=y*u.one, z=z*u.one)
        representation = representation.represent_as(UnitSphericalRepresentation)
    else:
        representation = CartesianRepresentation(x=x*c.unit, y=y*c.unit, z=z*c.unit)

    # if no obstime was given in the new frame, use the old one for consistency
    newobstime = fk4coord._obstime if fk4noeframe._obstime is None else fk4noeframe._obstime

    fk4noe = FK4NoETerms(representation, equinox=fk4coord.equinox,
                         obstime=newobstime)
    if fk4coord.equinox != fk4noeframe.equinox:
        #precession
        fk4noe = fk4noe.transform_to(fk4noeframe)
    return fk4noe


@frame_transform_graph.transform(FunctionTransform, FK4NoETerms, FK4)
def fk4_no_e_to_fk4(fk4noecoord, fk4frame):
    #first precess, if necessary
    if fk4noecoord.equinox != fk4frame.equinox:
        fk4noe_w_fk4equinox = FK4NoETerms(equinox=fk4frame.equinox,
                                          obstime=fk4noecoord.obstime)
        fk4noecoord = fk4noecoord.transform_to(fk4noe_w_fk4equinox)

    # Extract cartesian vector
    c = fk4noecoord.cartesian.xyz
    r = np.asarray(c.reshape((3, c.size // 3)))

    # Find distance (for re-normalization)
    d_orig = np.sqrt(np.sum(r ** 2))

    # Apply E-terms of aberration. Note that this depends on the equinox (not
    # the observing time/epoch) of the coordinates. See issue #1496 for a
    # discussion of this.
    eterms_a = np.asarray(fk4_e_terms(fk4noecoord.equinox))
    r0 = r.copy()
    for _ in range(10):
        r = (eterms_a.reshape(3, 1) + r0) / (1. + np.dot(eterms_a, r))

    # Find new distance (for re-normalization)
    d_new = np.sqrt(np.sum(r ** 2))

    # Renormalize
    r = r * d_orig / d_new

    subshape = c.shape[1:]
    x = r[0].reshape(subshape)
    y = r[1].reshape(subshape)
    z = r[2].reshape(subshape)

    #now re-cast into an appropriate Representation, and precess if need be
    if isinstance(fk4noecoord.data, UnitSphericalRepresentation):
        representation = CartesianRepresentation(x=x*u.one, y=y*u.one, z=z*u.one)
        representation = representation.represent_as(UnitSphericalRepresentation)
    else:
        representation = CartesianRepresentation(x=x*c.unit, y=y*c.unit, z=z*c.unit)

    return fk4frame.realize_frame(representation)


# FK5 to/from FK4 ------------------->
# B1950->J2000 matrix from Murray 1989 A&A 218,325 eqn 28
_B1950_TO_J2000_M = \
    np.mat([[0.9999256794956877, -0.0111814832204662, -0.0048590038153592],
            [0.0111814832391717,  0.9999374848933135, -0.0000271625947142],
            [0.0048590037723143, -0.0000271702937440,  0.9999881946023742]])

_FK4_CORR = \
    np.mat([[-0.0026455262, -1.1539918689, +2.1111346190],
            [+1.1540628161, -0.0129042997, +0.0236021478],
            [-2.1112979048, -0.0056024448, +0.0102587734]]) * 1.e-6


# This transformation can't be static because the observation date is needed.
@frame_transform_graph.transform(DynamicMatrixTransform, FK4NoETerms, FK5)
def fk4_no_e_to_fk5(fk4noecoord, fk5frame):
    # Correction terms for FK4 being a rotating system
    B = FK4NoETerms._fk4_B_matrix(fk4noecoord.obstime)

    # construct both precession matricies - if the equinoxes are B1950 and
    # J2000, these are just identity matricies
    pmat1 = fk4noecoord._precession_matrix(fk4noecoord.equinox, EQUINOX_B1950)
    pmat2 = fk5frame._precession_matrix(EQUINOX_J2000, fk5frame.equinox)

    return pmat2 * B * pmat1


# This transformation can't be static because the observation date is needed.
@frame_transform_graph.transform(DynamicMatrixTransform, FK5, FK4NoETerms)
def fk5_to_fk4_no_e(fk5coord, fk4noeframe):
    # Get transposed version of the rotating correction terms... so with the
    # transpose this takes us from FK5/J200 to FK4/B1950
    B = FK4NoETerms._fk4_B_matrix(fk4noeframe.obstime).T

    # construct both precession matricies - if the equinoxes are B1950 and
    # J2000, these are just identity matricies
    pmat1 = fk5coord._precession_matrix(fk5coord.equinox, EQUINOX_J2000)
    pmat2 = fk4noeframe._precession_matrix(EQUINOX_B1950, fk4noeframe.equinox)

    return pmat2 * B * pmat1
