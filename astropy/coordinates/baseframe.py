# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Framework and base classes for coordinate frames/"low-level" coordinate classes.
"""
from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

# Standard library
import abc

# Dependencies
import numpy as np

# Project
from ..extern import six
from .. import units as u
from .transformations import TransformGraph

__all__ = ['BaseCoordinateFrame', 'spatial_transform_graph']

# the graph used for all transformations between frames
frame_transform_graph = TransformGraph()


@six.add_metaclass(abc.ABCMeta)
class BaseCoordinateFrame(object):
    """ docstring """

    preferred_representation = None
    preferred_attr_names = None
    frame_attr_names = ()

    def __init__(self, representation):

        # TODO: check that representation is a valid object

        self._cache = dict()

    @property
    def data(self):
        if self._data is None:
            raise AttributeError()
        raise NotImplementedError()

    def represent_as(self, new_representation):
        return self.data.represent_as(new_representation)

    def transform_to(self, new_frame):
        """
        Transform this object's coordinate data to a new frame.

        Parameters
        ----------
        new_frame : class
            The frame to transform this coordinate frame into.

        Returns
        -------
        transframe
            A new object with the coordinate data represented in the
            ``newframe`` system.

        Raises
        ------
        ValueError
            If there is no possible transformation route.
        """
        from .errors import ConvertError

        if new_frame is self.__class__:
            return copy.deepcopy(self)

        trans = frame_transform_graph.get_transform(self.__class__, new_frame)
        if trans is None:
            raise ConvertError('Cannot transform from {0} to '
                               '{1}'.format(self.__class__, new_frame))
        return trans(self)

    def is_transformable_to(self, new_frame):
        """
        Determines if this coordinate frame can be transformed to another
        given frame.

        Parameters
        ----------
        new_frame : class
            The proposed frame to transform to.

        Returns
        -------
        transformable : bool or str
            `True` if this can be transformed to ``new_frame``, `False` if
            not. The string 'same' if ``new_frame`` is the same system as this
            object (i.e. no transformation is needed).
        """
        if self.__class__ is new_frame:
            return 'same'
        else:
            trans = frame_transform_graph.get_transform(self.__class__,
                                                        new_frame)
            return trans is not None

    @property
    def cartesian(self):
        """
        Shorthand for a cartesian representation of the coordinates in this object.
        """
        # TODO: if representations are updated to use a full transform graph,
        #       the representation aliases should not be hard-coded like this
        cached_repr = self._cache.get('cartesian', None)
        if not cached_repr:
            self._cache['cartesian'] = self.represent_as(CartesianRepresentation)
        return self._cache['cartesian']

    @property
    def spherical(self):
        """
        Shorthand for a spherical representation of the coordinates in this object.
        """
        # TODO: if representations are updated to use a full transform graph,
        #       the representation aliases should not be hard-coded like this
        cached_repr = self._cache.get('spherical', None)
        if not cached_repr:
            self._cache['spherical'] = self.represent_as(SphericalRepresentation)
        return self._cache['spherical']
