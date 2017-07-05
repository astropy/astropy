# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

from ... import units as u
from .. import representation as r
from ..baseframe import BaseCoordinateFrame, RepresentationMapping
from .galactic import Galactic


class Supergalactic(BaseCoordinateFrame):
    """
    Supergalactic Coordinates
    (see Lahav et al. 2000, <http://adsabs.harvard.edu/abs/2000MNRAS.312..166L>,
    and references therein).

    Parameters
    ----------
    representation : `BaseRepresentation` or None
        A representation object or None to have no data (or use the other keywords)

    sgl : `Angle`, optional, must be keyword
        The supergalactic longitude for this object (``sgb`` must also be given and
        ``representation`` must be None).
    sgb : `Angle`, optional, must be keyword
        The supergalactic latitude for this object (``sgl`` must also be given and
        ``representation`` must be None).
    distance : `~astropy.units.Quantity`, optional, must be keyword
        The Distance for this object along the line-of-sight.

    pm_sgl_cossgb : :class:`~astropy.units.Quantity`, optional, must be keyword
        The proper motion in Right Ascension for this object (``pm_sgb`` must
        also be given).
    pm_sgb : :class:`~astropy.units.Quantity`, optional, must be keyword
        The proper motion in Declination for this object (``pm_sgl_cossgb`` must
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

    frame_specific_representation_info = {
        r.SphericalRepresentation: [
            RepresentationMapping('lon', 'sgl'),
            RepresentationMapping('lat', 'sgb')
        ],
        r.CartesianRepresentation: [
            RepresentationMapping('x', 'sgx'),
            RepresentationMapping('y', 'sgy'),
            RepresentationMapping('z', 'sgz')
        ],
        r.SphericalCosLatDifferential: [
            RepresentationMapping('d_lon_coslat', 'pm_sgl_cossgb', u.mas/u.yr),
            RepresentationMapping('d_lat', 'pm_sgb', u.mas/u.yr),
            RepresentationMapping('d_distance', 'radial_velocity', u.km/u.s),
        ],
        r.SphericalDifferential: [
            RepresentationMapping('d_lon', 'pm_sgl', u.mas/u.yr),
            RepresentationMapping('d_lat', 'pm_sgb', u.mas/u.yr),
            RepresentationMapping('d_distance', 'radial_velocity', u.km/u.s),
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

    # North supergalactic pole in Galactic coordinates.
    # Needed for transformations to/from Galactic coordinates.
    _nsgp_gal = Galactic(l=47.37*u.degree, b=+6.32*u.degree)
