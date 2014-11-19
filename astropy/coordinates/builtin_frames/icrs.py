
from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

# Standard library
import inspect

# Dependencies
import numpy as np

# Project
from ...extern import six
from ...utils.compat.odict import OrderedDict
from ... import units as u
from ...time import Time
from ..angles import Angle
from ..representation import (SphericalRepresentation, CartesianRepresentation,
                             UnitSphericalRepresentation)
from ..baseframe import (BaseCoordinateFrame, frame_transform_graph, GenericFrame,
                        FrameAttribute, TimeFrameAttribute,
                        RepresentationMapping)
from ..transformations import FunctionTransform, DynamicMatrixTransform

class ICRS(BaseCoordinateFrame):
    """
    A coordinate or frame in the ICRS system.

    If you're looking for "J2000" coordinates, and aren't sure if you want to
    use this or `FK5`, you probably want to use ICRS. It's more well-defined as
    a catalog coordinate and is an inertial system, and is very close (within
    tens of milliarcseconds) to J2000 equatorial.

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
    distance : `~astropy.units.Quantity`, optional, must be keyword
        The Distance for this object along the line-of-sight.
        (``representation`` must be None).
    """

    frame_specific_representation_info = {
        'spherical': [RepresentationMapping('lon', 'ra'),
                      RepresentationMapping('lat', 'dec')]
    }
    frame_specific_representation_info['unitspherical'] = \
        frame_specific_representation_info['spherical']

    default_representation = SphericalRepresentation

    @staticmethod
    def _icrs_to_fk5_matrix():
        """
        B-matrix from USNO circular 179.  Used by the ICRS->FK5 transformation
        functions.
        """
        from ..angles import rotation_matrix

        eta0 = -19.9 / 3600000.
        xi0 = 9.1 / 3600000.
        da0 = -22.9 / 3600000.

        m1 = rotation_matrix(-eta0, 'x')
        m2 = rotation_matrix(xi0, 'y')
        m3 = rotation_matrix(da0, 'z')

        return m1 * m2 * m3

# define this because it only needs to be computed once
ICRS._ICRS_TO_FK5_J2000_MAT = ICRS._icrs_to_fk5_matrix()
