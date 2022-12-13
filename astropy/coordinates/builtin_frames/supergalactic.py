# Licensed under a 3-clause BSD style license - see LICENSE.rst

from astropy import units as u
from astropy.coordinates import representation as r
from astropy.coordinates.baseframe import (
    BaseCoordinateFrame,
    RepresentationMapping,
    base_doc,
)
from astropy.utils.decorators import format_doc

from .galactic import Galactic

__all__ = ["Supergalactic"]


doc_components = """
    sgl : `~astropy.coordinates.Angle`, optional, keyword-only
        The supergalactic longitude for this object (``sgb`` must also be given and
        ``representation`` must be None).
    sgb : `~astropy.coordinates.Angle`, optional, keyword-only
        The supergalactic latitude for this object (``sgl`` must also be given and
        ``representation`` must be None).
    distance : `~astropy.units.Quantity` ['speed'], optional, keyword-only
        The Distance for this object along the line-of-sight.

    pm_sgl_cossgb : `~astropy.units.Quantity` ['angular speed'], optional, keyword-only
        The proper motion in Right Ascension for this object (``pm_sgb`` must
        also be given).
    pm_sgb : `~astropy.units.Quantity` ['angular speed'], optional, keyword-only
        The proper motion in Declination for this object (``pm_sgl_cossgb`` must
        also be given).
    radial_velocity : `~astropy.units.Quantity` ['speed'], optional, keyword-only
        The radial velocity of this object.
"""


@format_doc(base_doc, components=doc_components, footer="")
class Supergalactic(BaseCoordinateFrame):
    """
    Supergalactic Coordinates
    (see Lahav et al. 2000, <https://ui.adsabs.harvard.edu/abs/2000MNRAS.312..166L>,
    and references therein).
    """

    frame_specific_representation_info = {
        r.SphericalRepresentation: [
            RepresentationMapping("lon", "sgl"),
            RepresentationMapping("lat", "sgb"),
        ],
        r.CartesianRepresentation: [
            RepresentationMapping("x", "sgx"),
            RepresentationMapping("y", "sgy"),
            RepresentationMapping("z", "sgz"),
        ],
        r.CartesianDifferential: [
            RepresentationMapping("d_x", "v_x", u.km / u.s),
            RepresentationMapping("d_y", "v_y", u.km / u.s),
            RepresentationMapping("d_z", "v_z", u.km / u.s),
        ],
    }

    default_representation = r.SphericalRepresentation
    default_differential = r.SphericalCosLatDifferential

    # North supergalactic pole in Galactic coordinates.
    # Needed for transformations to/from Galactic coordinates.
    _nsgp_gal = Galactic(l=47.37 * u.degree, b=+6.32 * u.degree)
