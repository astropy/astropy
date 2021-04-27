# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

from astropy.utils.decorators import format_doc
from astropy.coordinates import representation as r
from astropy.coordinates.baseframe import BaseCoordinateFrame, RepresentationMapping, base_doc

__all__ = ['BaseRADecFrame']


doc_components = """
    ra : `~astropy.coordinates.Angle`, optional, keyword-only
        The RA for this object (``dec`` must also be given and ``representation``
        must be None).
    dec : `~astropy.coordinates.Angle`, optional, keyword-only
        The Declination for this object (``ra`` must also be given and
        ``representation`` must be None).
    distance : `~astropy.units.Quantity` ['length'], optional, keyword-only
        The Distance for this object along the line-of-sight.
        (``representation`` must be None).

    pm_ra_cosdec : `~astropy.units.Quantity` ['angular speed'], optional, keyword-only
        The proper motion in Right Ascension (including the ``cos(dec)`` factor)
        for this object (``pm_dec`` must also be given).
    pm_dec : `~astropy.units.Quantity` ['angular speed'], optional, keyword-only
        The proper motion in Declination for this object (``pm_ra_cosdec`` must
        also be given).
    radial_velocity : `~astropy.units.Quantity` ['speed'], optional, keyword-only
        The radial velocity of this object.
"""


@format_doc(base_doc, components=doc_components, footer="")
class BaseRADecFrame(BaseCoordinateFrame):
    """
    A base class that defines default representation info for frames that
    represent longitude and latitude as Right Ascension and Declination
    following typical "equatorial" conventions.
    """
    frame_specific_representation_info = {
        r.SphericalRepresentation: [
            RepresentationMapping('lon', 'ra'),
            RepresentationMapping('lat', 'dec')
        ]
    }

    default_representation = r.SphericalRepresentation
    default_differential = r.SphericalCosLatDifferential
