# -*- coding: utf-8 -*-

import numpy as np

from .icrs import ICRS
from .altaz import AltAz
from ..transformations import DynamicMatrixTransform, FunctionTransform
from ..baseframe import CoordinateAttribute, frame_transform_graph
from ..angles import rotation_matrix


_astrometric_cache = {}


def make_astrometric_cls(framecls):
    """
    Create a new class that is the Astrometric frame for a specific class of
    origin frame. If such a class has already been created for this frame, the
    same class will be returned.

    Parameters
    ----------
    framecls : coordinate frame class (i.e., subclass of `~astropy.coordinates.BaseCoordinateFrame`)
        The class to create the Astrometric frame of.

    Notes
    -----
    This function is necessary because Astropy's frame transformations depend
    on connection between specific frame *classes*.  So each type of frame
    needs its own distinct astrometric frame class.  This function generates
    just that class, as well as ensuring that only one example of such a class
    actually gets created in any given python session.
    """

    if framecls in _astrometric_cache:
        return _astrometric_cache[framecls]

    class Astrometric(framecls):
        """
        A frame which is relative to some position on the sky. Useful for
        calculating offsets and dithers in the frame of the sky.

        Parameters
        ----------
        representation : `BaseRepresentation` or None
            A representation object or None to have no data (or use the other keywords)
        origin : `SkyCoord` or low-level coordinate object.
            the coordinate which specifiy the origin of this frame.

        """

        origin = CoordinateAttribute(default=None, frame=framecls)

    @frame_transform_graph.transform(FunctionTransform, Astrometric, Astrometric)
    def astrometric_to_astrometric(from_astrometric_coord, to_astrometric_frame):
        """Transform between two astrometric frames."""

        # If both frames have on-sky positions, then the transform should happen relative to both origins.
        if (from_astrometric_coord.origin is not None) and (to_astrometric_frame.origin is not None):
            return from_astrometric_coord.transform_to(framecls).transform_to(to_astrometric_frame)

        # Otherwise, the transform occurs just by setting the new origin.
        return to_astrometric_frame.realize_frame(from_astrometric_coord.cartesian)

    @frame_transform_graph.transform(DynamicMatrixTransform, framecls, Astrometric)
    def icrs_to_astrometric(reference_frame, astrometric_frame):
        """Convert an ICRS coordinate to an Astrometric frame."""

        # Define rotation matricies along the position angle vector, and
        # relative to the origin.
        origin = astrometric_frame.origin.spherical
        mat2 = rotation_matrix(-origin.lat, 'y')
        mat3 = rotation_matrix(origin.lon, 'z')
        R = mat2 * mat3
        return R

    @frame_transform_graph.transform(DynamicMatrixTransform, Astrometric, framecls)
    def astrometric_to_icrs(astrometric_coord, reference_frame):
        """Convert an Astrometric frame coordinate to an ICRS"""

        # use the forward transform, but just invert it
        R = icrs_to_astrometric(reference_frame, astrometric_coord)
        return R.T  # this is the inverse because R is a rotation matrix

    _astrometric_cache[framecls] = Astrometric
    return Astrometric

AstrometricICRS = make_astrometric_cls(ICRS)
AstrometricAltAz = make_astrometric_cls(AltAz)
