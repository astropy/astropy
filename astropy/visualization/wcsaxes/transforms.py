# Licensed under a 3-clause BSD style license - see LICENSE.rst


# Note: This file includes code derived from pywcsgrid2
#
# This file contains Matplotlib transformation objects (e.g. from pixel to world
# coordinates, but also world-to-world).

import abc

import numpy as np
from matplotlib.path import Path
from matplotlib.transforms import Transform

from astropy import units as u
from astropy.coordinates import (
    BaseCoordinateFrame,
    SkyCoord,
    UnitSphericalRepresentation,
    frame_transform_graph,
)

__all__ = [
    "CoordinateTransform",
    "CurvedTransform",
    "Pixel2WorldTransform",
    "World2PixelTransform",
]


class CurvedTransform(Transform, metaclass=abc.ABCMeta):
    """
    Abstract base class for non-affine curved transforms.
    """

    input_dims = 2
    output_dims = 2
    is_separable = False

    def transform_path(self, path):
        """
        Transform a Matplotlib Path.

        Parameters
        ----------
        path : :class:`~matplotlib.path.Path`
            The path to transform

        Returns
        -------
        path : :class:`~matplotlib.path.Path`
            The resulting path
        """
        return Path(self.transform(path.vertices), path.codes)

    transform_path_non_affine = transform_path

    def transform(self, input):
        raise NotImplementedError("")

    def inverted(self):
        raise NotImplementedError("")


class CoordinateTransform(CurvedTransform):
    has_inverse = True

    def __init__(self, input_system, output_system):
        super().__init__()
        self._input_system_name = input_system
        self._output_system_name = output_system

        if isinstance(self._input_system_name, str):
            frame_cls = frame_transform_graph.lookup_name(self._input_system_name)
            if frame_cls is None:
                raise ValueError(f"Frame {self._input_system_name} not found")
            else:
                self.input_system = frame_cls()
        elif isinstance(self._input_system_name, BaseCoordinateFrame):
            self.input_system = self._input_system_name
        else:
            raise TypeError(
                "input_system should be a WCS instance, string, or a coordinate frame"
                " instance"
            )

        if isinstance(self._output_system_name, str):
            frame_cls = frame_transform_graph.lookup_name(self._output_system_name)
            if frame_cls is None:
                raise ValueError(f"Frame {self._output_system_name} not found")
            else:
                self.output_system = frame_cls()

        elif isinstance(self._output_system_name, BaseCoordinateFrame):
            self.output_system = self._output_system_name
        else:
            raise TypeError(
                "output_system should be a WCS instance, string, or a coordinate frame"
                " instance"
            )

        if self.output_system == self.input_system:
            self.same_frames = True
        else:
            self.same_frames = False

    @property
    def same_frames(self):
        return self._same_frames

    @same_frames.setter
    def same_frames(self, same_frames):
        self._same_frames = same_frames

    def transform(self, input_coords):
        """
        Transform one set of coordinates to another.
        """
        if self.same_frames:
            return input_coords

        input_coords = input_coords * u.deg
        x_in, y_in = input_coords[:, 0], input_coords[:, 1]

        c_in = SkyCoord(
            UnitSphericalRepresentation(x_in, y_in), frame=self.input_system
        )

        # We often need to transform arrays that contain NaN values, and filtering
        # out the NaN values would have a performance hit, so instead we just pass
        # on all values and just ignore Numpy warnings
        with np.errstate(all="ignore"):
            c_out = c_in.transform_to(self.output_system)

        lon = c_out.spherical.lon.deg
        lat = c_out.spherical.lat.deg

        return np.concatenate((lon[:, np.newaxis], lat[:, np.newaxis]), axis=1)

    transform_non_affine = transform

    def inverted(self):
        """
        Return the inverse of the transform.
        """
        return CoordinateTransform(self._output_system_name, self._input_system_name)


class World2PixelTransform(CurvedTransform, metaclass=abc.ABCMeta):
    """
    Base transformation from world to pixel coordinates.
    """

    has_inverse = True
    frame_in = None

    @property
    @abc.abstractmethod
    def input_dims(self):
        """
        The number of input world dimensions.
        """

    @abc.abstractmethod
    def transform(self, world):
        """
        Transform world to pixel coordinates. You should pass in a NxM array
        where N is the number of points to transform, and M is the number of
        dimensions. This then returns the (x, y) pixel coordinates
        as a Nx2 array.
        """

    @abc.abstractmethod
    def inverted(self):
        """
        Return the inverse of the transform.
        """


class Pixel2WorldTransform(CurvedTransform, metaclass=abc.ABCMeta):
    """
    Base transformation from pixel to world coordinates.
    """

    has_inverse = True
    frame_out = None

    @property
    @abc.abstractmethod
    def output_dims(self):
        """
        The number of output world dimensions.
        """

    @abc.abstractmethod
    def transform(self, pixel):
        """
        Transform pixel to world coordinates. You should pass in a Nx2 array
        of (x, y) pixel coordinates to transform to world coordinates. This
        will then return an NxM array where M is the number of dimensions.
        """

    @abc.abstractmethod
    def inverted(self):
        """
        Return the inverse of the transform.
        """
