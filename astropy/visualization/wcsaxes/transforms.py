# Licensed under a 3-clause BSD style license - see LICENSE.rst


# Note: This file incldues code dervived from pywcsgrid2
#
# This file contains Matplotlib transformation objects (e.g. from pixel to world
# coordinates, but also world-to-world).

import abc

import numpy as np

from matplotlib.path import Path
from matplotlib.transforms import Transform

from astropy import units as u
from astropy.coordinates import (SkyCoord, frame_transform_graph,
                                 UnitSphericalRepresentation,
                                 BaseCoordinateFrame)

__all__ = ['CurvedTransform', 'CoordinateTransform']


class CurvedTransform(Transform, metaclass=abc.ABCMeta):
    """
    Abstract base class for non-affine curved transforms
    """

    input_dims = 2
    output_dims = 2
    is_separable = False

    def transform_path(self, path):
        """
        Transform a Matplotlib Path

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
            self.input_system = frame_transform_graph.lookup_name(self._input_system_name)
            if self.input_system is None:
                raise ValueError("Frame {0} not found".format(self._input_system_name))
        elif isinstance(self._input_system_name, BaseCoordinateFrame):
            self.input_system = self._input_system_name
        else:
            raise TypeError("input_system should be a WCS instance, string, or a coordinate frame instance")

        if isinstance(self._output_system_name, str):
            self.output_system = frame_transform_graph.lookup_name(self._output_system_name)
            if self.output_system is None:
                raise ValueError("Frame {0} not found".format(self._output_system_name))
        elif isinstance(self._output_system_name, BaseCoordinateFrame):
            self.output_system = self._output_system_name
        else:
            raise TypeError("output_system should be a WCS instance, string, or a coordinate frame instance")

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
        Transform one set of coordinates to another
        """
        if self.same_frames:
            return input_coords

        input_coords = input_coords*u.deg
        x_in, y_in = input_coords[:, 0], input_coords[:, 1]

        c_in = SkyCoord(UnitSphericalRepresentation(x_in, y_in),
                        frame=self.input_system)

        # We often need to transform arrays that contain NaN values, and filtering
        # out the NaN values would have a performance hit, so instead we just pass
        # on all values and just ignore Numpy warnings
        with np.errstate(all='ignore'):
            c_out = c_in.transform_to(self.output_system)

        lon = c_out.spherical.lon.deg
        lat = c_out.spherical.lat.deg

        return np.concatenate((lon[:, np.newaxis], lat[:, np.newaxis]), axis=1)

    transform_non_affine = transform

    def inverted(self):
        """
        Return the inverse of the transform
        """
        return CoordinateTransform(self._output_system_name, self._input_system_name)
