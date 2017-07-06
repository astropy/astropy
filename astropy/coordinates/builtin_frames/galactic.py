# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

from ... import units as u
from ..angles import Angle
from .. import representation as r
from ..baseframe import BaseCoordinateFrame, RepresentationMapping

# these are needed for defining the NGP
from .fk5 import FK5
from .fk4 import FK4NoETerms


class Galactic(BaseCoordinateFrame):
    """
    A coordinate or frame in the Galactic coordinate system.

    This frame is used in a variety of Galactic contexts because it has as its
    x-y plane the plane of the Milky Way.  The positive x direction (i.e., the
    l=0, b=0 direction) points to the center of the Milky Way and the z-axis
    points toward the North Galactic Pole (following the IAU's 1958 definition
    [1]_). However, unlike the `~astropy.coordinates.Galactocentric` frame, the
    *origin* of this frame in 3D space is the solar system barycenter, not
    the center of the Milky Way.

    Parameters
    ----------
    representation : `BaseRepresentation` or None
        A representation object or None to have no data (or use the other keywords)

    l : `Angle`, optional, must be keyword
        The Galactic longitude for this object (``b`` must also be given and
        ``representation`` must be None).
    b : `Angle`, optional, must be keyword
        The Galactic latitude for this object (``l`` must also be given and
        ``representation`` must be None).
    distance : `~astropy.units.Quantity`, optional, must be keyword
        The Distance for this object along the line-of-sight.

    pm_l_cosb : :class:`~astropy.units.Quantity`, optional, must be keyword
        The proper motion in Galactic longitude (including the ``cos(b)`` term)
        for this object (``pm_b`` must also be given).
    pm_b : :class:`~astropy.units.Quantity`, optional, must be keyword
        The proper motion in Galactic latitude for this object (``pm_l_cosb``
        must also be given).
    radial_velocity : :class:`~astropy.units.Quantity`, optional, must be keyword
        The radial velocity of this object.

    copy : bool, optional
        If `True` (default), make copies of the input coordinate arrays.
        Can only be passed in as a keyword argument.

    differential_cls : `BaseDifferential`, dict, optional
        A differential class or dictionary of differential classes (currently
        only a velocity differential with key 's' is supported). This sets
        the expected input differential class, thereby changing the expected
        keyword arguments of the data passed in. For example, passing
        ``differential_cls=CartesianDifferential`` will make the classes
        expect velocity data with the argument names ``v_x, v_y, v_z``.

    Notes
    -----
    .. [1] Blaauw, A.; Gum, C. S.; Pawsey, J. L.; Westerhout, G. (1960), "The
       new I.A.U. system of galactic coordinates (1958 revision),"
       `MNRAS, Vol 121, pp.123 <http://adsabs.harvard.edu/abs/1960MNRAS.121..123B>`_.
    """

    frame_specific_representation_info = {
        r.SphericalRepresentation: [
            RepresentationMapping('lon', 'l'),
            RepresentationMapping('lat', 'b')
        ],
        r.CartesianRepresentation: [
            RepresentationMapping('x', 'u'),
            RepresentationMapping('y', 'v'),
            RepresentationMapping('z', 'w')
        ],
        r.CartesianDifferential: [
            RepresentationMapping('d_x', 'U', u.km/u.s),
            RepresentationMapping('d_y', 'V', u.km/u.s),
            RepresentationMapping('d_z', 'W', u.km/u.s)
        ],
        r.SphericalCosLatDifferential: [
            RepresentationMapping('d_lon_coslat', 'pm_l_cosb', u.mas/u.yr),
            RepresentationMapping('d_lat', 'pm_b', u.mas/u.yr),
            RepresentationMapping('d_distance', 'radial_velocity', u.km/u.s),
        ],
        r.SphericalDifferential: [
            RepresentationMapping('d_lon', 'pm_l', u.mas/u.yr),
            RepresentationMapping('d_lat', 'pm_b', u.mas/u.yr),
            RepresentationMapping('d_distance', 'radial_velocity', u.km/u.s),
        ]
    }
    frame_specific_representation_info[r.UnitSphericalRepresentation] = \
        frame_specific_representation_info[r.SphericalRepresentation]
    frame_specific_representation_info[r.UnitSphericalCosLatDifferential] = \
        frame_specific_representation_info[r.SphericalCosLatDifferential]
    frame_specific_representation_info[r.UnitSphericalDifferential] = \
        frame_specific_representation_info[r.SphericalDifferential]

    default_representation = r.SphericalRepresentation
    default_differential = r.SphericalCosLatDifferential

    # North galactic pole and zeropoint of l in FK4/FK5 coordinates. Needed for
    # transformations to/from FK4/5

    # These are from the IAU's definition of galactic coordinates
    _ngp_B1950 = FK4NoETerms(ra=192.25*u.degree, dec=27.4*u.degree)
    _lon0_B1950 = Angle(123, u.degree)

    # These are *not* from Reid & Brunthaler 2004 - instead, they were
    # derived by doing:
    #
    # >>> FK4NoETerms(ra=192.25*u.degree, dec=27.4*u.degree).transform_to(FK5)
    #
    # This gives better consistency with other codes than using the values
    # from Reid & Brunthaler 2004 and the best self-consistency between FK5
    # -> Galactic and FK5 -> FK4 -> Galactic. The lon0 angle was found by
    # optimizing the self-consistency.
    _ngp_J2000 = FK5(ra=192.8594812065348*u.degree, dec=27.12825118085622*u.degree)
    _lon0_J2000 = Angle(122.9319185680026, u.degree)
