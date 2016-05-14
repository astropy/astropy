# -*- coding: utf-8 -*-

import numpy as np

import astropy.units as u

from .icrs import ICRS
from ..sky_coordinate import SkyCoord
from ..transformations import DynamicMatrixTransform, FunctionTransform
from ..baseframe import FrameAttribute, frame_transform_graph
from ..angles import rotation_matrix


class CoordinateLocationAttribute(FrameAttribute):
    """A frame attribute which is a coordinates object."""

    def __init__(self, frame, default=None):
        self._frame = frame
        super(CoordinateLocationAttribute, self).__init__(default)

    def convert_input(self, value):
        """
        Checks that the input is a SkyCoord with the necessary units (or the
        special value ``None``).

        Parameters
        ----------
        value : object
            Input value to be converted.

        Returns
        -------
        out, converted : correctly-typed object, boolean
            Tuple consisting of the correctly-typed object and a boolean which
            indicates if conversion was actually performed.

        Raises
        ------
        ValueError
            If the input is not valid for this attribute.
        """
        if value is None:
            return None, False
        elif isinstance(value, self._frame):
            return value, False
        else:
            if not hasattr(value, 'transform_to'):
                raise ValueError('"{0}" was passed into an '
                                 'CoordinateLocationAttribute, but it does not have '
                                 '"transform_to" method'.format(value))
            transformedobj = value.transform_to(self._frame)
            if isinstance(transformedobj, SkyCoord):
                return transformedobj.frame, True
            return transformedobj, True


def astrometric(frame):
    """Define an Astrometric frame relative to some origin frame."""

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

        origin = CoordinateLocationAttribute(default=None, frame=frame)

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

    return Astrometric

AstrometricICRS = astrometric(ICRS)
