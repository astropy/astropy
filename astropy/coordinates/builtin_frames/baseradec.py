# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

from ... import units as u
from .. import representation as r
from ..baseframe import BaseCoordinateFrame, RepresentationMapping

__all__ = ['BaseRADecFrame']

_base_radec_docstring = """Parameters
    ----------
    representation : `BaseRepresentation` or None
        A representation object or ``None`` to have no data (or use the other
        keywords below).

    ra : `Angle`, optional, must be keyword
        The RA for this object (``dec`` must also be given and ``representation``
        must be None).
    dec : `Angle`, optional, must be keyword
        The Declination for this object (``ra`` must also be given and
        ``representation`` must be None).
    distance : `~astropy.units.Quantity`, optional, must be keyword
        The Distance for this object along the line-of-sight.
        (``representation`` must be None).

    pm_ra_cosdec : :class:`~astropy.units.Quantity`, optional, must be keyword
        The proper motion in Right Ascension (including the ``cos(dec)`` factor)
        for this object (``pm_dec`` must also be given).
    pm_dec : :class:`~astropy.units.Quantity`, optional, must be keyword
        The proper motion in Declination for this object (``pm_ra_cosdec`` must
        also be given).
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
"""


class BaseRADecFrame(BaseCoordinateFrame):
    """
    A base class that defines default representation info for frames that
    represent longitude and latitude as Right Ascension and Declination
    following typical "equatorial" conventions.

    {params}
    """
    frame_specific_representation_info = {
        r.SphericalRepresentation: [
            RepresentationMapping('lon', 'ra'),
            RepresentationMapping('lat', 'dec')
        ],
        r.SphericalCosLatDifferential: [
            RepresentationMapping('d_lon_coslat', 'pm_ra_cosdec', u.mas/u.yr),
            RepresentationMapping('d_lat', 'pm_dec', u.mas/u.yr),
            RepresentationMapping('d_distance', 'radial_velocity', u.km/u.s)
        ],
        r.SphericalDifferential: [
            RepresentationMapping('d_lon', 'pm_ra', u.mas/u.yr),
            RepresentationMapping('d_lat', 'pm_dec', u.mas/u.yr),
            RepresentationMapping('d_distance', 'radial_velocity', u.km/u.s)
        ],
        r.CartesianDifferential: [
            RepresentationMapping('d_x', 'v_x', u.km/u.s),
            RepresentationMapping('d_y', 'v_y', u.km/u.s),
            RepresentationMapping('d_z', 'v_z', u.km/u.s)
        ],
    }
    frame_specific_representation_info[r.UnitSphericalRepresentation] = \
        frame_specific_representation_info[r.SphericalRepresentation]
    frame_specific_representation_info[r.UnitSphericalCosLatDifferential] = \
        frame_specific_representation_info[r.SphericalCosLatDifferential]
    frame_specific_representation_info[r.UnitSphericalDifferential] = \
        frame_specific_representation_info[r.SphericalDifferential]

    default_representation = r.SphericalRepresentation
    default_differential = r.SphericalCosLatDifferential


BaseRADecFrame.__doc__ = BaseRADecFrame.__doc__.format(
    params=_base_radec_docstring)
