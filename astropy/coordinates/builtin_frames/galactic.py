
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



from .fk5 import FK5
from .fk4 import FK4NoETerms


# The UTC time scale is not properly defined prior to 1960, so Time('B1950',
# scale='utc') will emit a warning. Instead, we use Time('B1950', scale='tai')
# which is equivalent, but does not emit a warning.
_EQUINOX_J2000 = Time('J2000', scale='utc')
_EQUINOX_B1950 = Time('B1950', scale='tai')


class Galactic(BaseCoordinateFrame):
    """
    Galactic Coordinates.

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
    """

    frame_specific_representation_info = {
        'spherical': [RepresentationMapping('lon', 'l'),
                      RepresentationMapping('lat', 'b')],
        'cartesian': [RepresentationMapping('x', 'w'),
                      RepresentationMapping('y', 'u'),
                      RepresentationMapping('z', 'v')]
    }
    frame_specific_representation_info['unitspherical'] = \
        frame_specific_representation_info['spherical']

    default_representation = SphericalRepresentation

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


# Galactic to/from FK4/FK5 ----------------------->
# can't be static because the equinox is needed
@frame_transform_graph.transform(DynamicMatrixTransform, FK5, Galactic)
def fk5_to_gal(fk5coord, galframe):
    from ..angles import rotation_matrix

    #need precess to J2000 first
    pmat = fk5coord._precession_matrix(fk5coord.equinox, _EQUINOX_J2000)
    mat1 = rotation_matrix(180 - Galactic._lon0_J2000.degree, 'z')
    mat2 = rotation_matrix(90 - Galactic._ngp_J2000.dec.degree, 'y')
    mat3 = rotation_matrix(Galactic._ngp_J2000.ra.degree, 'z')

    return mat1 * mat2 * mat3 * pmat


@frame_transform_graph.transform(DynamicMatrixTransform, Galactic, FK5)
def _gal_to_fk5(galcoord, fk5frame):
    return fk5_to_gal(fk5frame, galcoord).T


@frame_transform_graph.transform(DynamicMatrixTransform, FK4NoETerms, Galactic)
def fk4_to_gal(fk4coords, galframe):
    from ..angles import rotation_matrix

    mat1 = rotation_matrix(180 - Galactic._lon0_B1950.degree, 'z')
    mat2 = rotation_matrix(90 - Galactic._ngp_B1950.dec.degree, 'y')
    mat3 = rotation_matrix(Galactic._ngp_B1950.ra.degree, 'z')
    matprec = fk4coords._precession_matrix(fk4coords.equinox, _EQUINOX_B1950)

    return mat1 * mat2 * mat3 * matprec


@frame_transform_graph.transform(DynamicMatrixTransform, Galactic, FK4NoETerms)
def gal_to_fk4(galcoords, fk4frame):
    return fk4_to_gal(fk4frame, galcoords).T
