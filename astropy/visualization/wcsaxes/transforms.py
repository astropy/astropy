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
from astropy.wcs import WCS
from astropy.wcs.utils import wcs_to_celestial_frame
from astropy.coordinates import (SkyCoord, frame_transform_graph,
                            SphericalRepresentation,
                            UnitSphericalRepresentation,
                            BaseCoordinateFrame)


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


class WCSWorld2PixelTransform(CurvedTransform):
    """
    WCS transformation from world to pixel coordinates
    """

    has_inverse = True

    def __init__(self, wcs, slice=None):
        super().__init__()
        self.wcs = wcs
        if self.wcs.wcs.naxis > 2:
            if slice is None:
                raise ValueError("WCS has more than 2 dimensions, so ``slice`` should be set")
            elif len(slice) != self.wcs.wcs.naxis:
                raise ValueError("slice should have as many elements as WCS "
                                 "has dimensions (should be {0})".format(self.wcs.wcs.naxis))
            else:
                self.slice = slice
                self.x_index = slice.index('x')
                self.y_index = slice.index('y')
        else:
            self.slice = None

    def __eq__(self, other):
        return (isinstance(other, type(self)) and self.wcs == other.wcs
                and self.slice == other.slice)

    @property
    def input_dims(self):
        return self.wcs.wcs.naxis

    def transform(self, world):
        """
        Transform world to pixel coordinates. You should pass in a NxM array
        where N is the number of points to transform, and M is the number of
        dimensions in the WCS. This then returns the (x, y) pixel coordinates
        as a Nx2 array.
        """

        if world.shape[1] != self.wcs.wcs.naxis:
            raise ValueError("Second dimension of input values should match number of WCS coordinates")

        if world.shape[0] == 0:
            pixel = np.zeros((0, 2))
        else:
            pixel = self.wcs.wcs_world2pix(world, 1) - 1

        if self.slice is None:
            return pixel
        else:
            return pixel[:, (self.x_index, self.y_index)]

    transform_non_affine = transform

    def inverted(self):
        """
        Return the inverse of the transform
        """
        return WCSPixel2WorldTransform(self.wcs, slice=self.slice)


class WCSPixel2WorldTransform(CurvedTransform):
    """
    WCS transformation from pixel to world coordinates
    """

    has_inverse = True

    def __init__(self, wcs, slice=None):
        super().__init__()
        self.wcs = wcs
        self.slice = slice
        if self.slice is not None:
            self.x_index = slice.index('x')
            self.y_index = slice.index('y')

    def __eq__(self, other):
        return (isinstance(other, type(self)) and self.wcs == other.wcs
                and self.slice == other.slice)

    @property
    def output_dims(self):
        return self.wcs.wcs.naxis

    def transform(self, pixel):
        """
        Transform pixel to world coordinates. You should pass in a Nx2 array
        of (x, y) pixel coordinates to transform to world coordinates. This
        will then return an NxM array where M is the number of dimensions in
        the WCS
        """

        if self.slice is None:
            pixel_full = pixel.copy()
        else:
            pixel_full = []
            for index in self.slice:
                if index == 'x':
                    pixel_full.append(pixel[:, 0])
                elif index == 'y':
                    pixel_full.append(pixel[:, 1])
                else:
                    pixel_full.append(index)
            pixel_full = np.array(np.broadcast_arrays(*pixel_full)).transpose()

        pixel_full += 1

        if pixel_full.shape[0] == 0:
            world = np.zeros((0, 2))
        else:
            world = self.wcs.wcs_pix2world(pixel_full, 1)

        # At the moment, one has to manually check that the transformation
        # round-trips, otherwise it should be considered invalid.
        pixel_check = self.wcs.wcs_world2pix(world, 1)
        with np.errstate(invalid='ignore'):
            invalid = np.any(np.abs(pixel_check - pixel_full) > 1., axis=1)
        world[invalid] = np.nan

        return world

    transform_non_affine = transform

    def inverted(self):
        """
        Return the inverse of the transform
        """
        return WCSWorld2PixelTransform(self.wcs, slice=self.slice)


class CoordinateTransform(CurvedTransform):

    has_inverse = True

    def __init__(self, input_system, output_system):
        super().__init__()
        self._input_system_name = input_system
        self._output_system_name = output_system

        if isinstance(self._input_system_name, WCS):
            self.input_system = wcs_to_celestial_frame(self._input_system_name)
        elif isinstance(self._input_system_name, str):
            self.input_system = frame_transform_graph.lookup_name(self._input_system_name)
            if self.input_system is None:
                raise ValueError("Frame {0} not found".format(self._input_system_name))
        elif isinstance(self._input_system_name, BaseCoordinateFrame):
            self.input_system = self._input_system_name
        else:
            raise TypeError("input_system should be a WCS instance, string, or a coordinate frame instance")

        if isinstance(self._output_system_name, WCS):
            self.output_system = wcs_to_celestial_frame(self._output_system_name)
        elif isinstance(self._output_system_name, str):
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
