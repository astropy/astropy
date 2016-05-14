# -*- coding: utf-8 -*-

import numpy as np

from .icrs import ICRS
from .altaz import AltAz
from ..transformations import DynamicMatrixTransform, FunctionTransform
from ..baseframe import CoordinateAttribute, frame_transform_graph
from ..angles import rotation_matrix


_astrometric_cache = {}


def make_astrometric_cls(frame):
    """Create an Astrometric frame relative to some origin frame."""

    if frame in _astrometric_cache:
        return _astrometric_cache[frame]

    class Astrometric(frame):
        """
        A frame which is relative to some position on the sky. Useful for calculating offsets and dithers
        in the frame of the sky.

        Parameters
        ----------
        representation : `BaseRepresentation` or None
            A representation object or None to have no data (or use the other keywords)
        X : `Angle`, optional, must be keyword
        Y : `Angle`, optional, must be keyword
        origin : `SkyCoord`, optional
            the coordinates which specifiy the origin of this frame.

        """

        origin = CoordinateAttribute(default=None, frame=frame)

        def at_origin(self, origin):
            """This method returns a new frame, identical to this frame, except at a differnt origin.
            This can be used to apply e.g. a sequence of telescope offsets to different targets (different origins).

            Parameters
            ----------
            origin : `SkyCoord`, optional
                the coordinates which specifiy the origin of this frame.

            Returns
            -------
            frame : `Astrometric`
                A new `Astrometric` which is centered at `origin`.

            """
            attrs = {}
            for name, value in self.get_frame_attr_names().items():
                attrs[name] = getattr(self, name, value)
            attrs['origin'] = origin
            return self.__class__(**attrs)

    @frame_transform_graph.transform(FunctionTransform, Astrometric, Astrometric)
    def astrometric_to_astrometric(from_astrometric_coord, to_astrometric_frame):
        """Transform between two astrometric frames."""

        # If both frames have on-sky positions, then the transform should happen relative to both origins.
        if (from_astrometric_coord.origin is not None) and (to_astrometric_frame.origin is not None):
            return from_astrometric_coord.transform_to(frame).transform_to(to_astrometric_frame)

        # Otherwise, the transform occurs just by setting the new origin.
        return to_astrometric_frame.realize_frame(from_astrometric_coord.cartesian)

    @frame_transform_graph.transform(DynamicMatrixTransform, frame, Astrometric)
    def icrs_to_astrometric(reference_frame, astrometric_frame):
        """Convert an ICRS coordinate to an Astrometric frame."""

        # Define rotation matricies along the position angle vector, and
        # relative to the origin.
        origin = astrometric_frame.origin.spherical
        mat2 = rotation_matrix(-origin.lat, 'y')
        mat3 = rotation_matrix(origin.lon, 'z')
        R = mat2 * mat3
        return R

    @frame_transform_graph.transform(DynamicMatrixTransform, Astrometric, frame)
    def astrometric_to_icrs(astrometric_coord, reference_frame):
        """Convert an Astrometric frame coordinate to an ICRS"""

        # Define rotation matricies along the position angle vector, and
        # relative to the origin.
        origin = astrometric_coord.origin.spherical
        mat2 = rotation_matrix(-origin.lat, 'y')
        mat3 = rotation_matrix(origin.lon, 'z')
        R = mat2 * mat3
        Rinv = np.linalg.inv(R)
        return Rinv

    _astrometric_cache[frame] = Astrometric
    return Astrometric

AstrometricICRS = make_astrometric_cls(ICRS)
AstrometricAltAz = make_astrometric_cls(AltAz)
