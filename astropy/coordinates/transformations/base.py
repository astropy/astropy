# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Coordinate transformations base class.

This module contains the abstract base class for coordinate transformations. For
concrete implementations, see the submodules
:mod:`~astropy.coordinates.transformations.affine`,
:mod:`~astropy.coordinates.transformations.function`, and
:mod:`~astropy.coordinates.transformations.composite`.
"""

from abc import ABCMeta, abstractmethod

__all__ = ["CoordinateTransform"]


class CoordinateTransform(metaclass=ABCMeta):
    """
    An object that transforms a coordinate from one system to another.
    Subclasses must implement `__call__` with the provided signature.
    They should also call this superclass's ``__init__`` in their
    ``__init__``.

    Parameters
    ----------
    fromsys : `~astropy.coordinates.BaseCoordinateFrame` subclass
        The coordinate frame class to start from.
    tosys : `~astropy.coordinates.BaseCoordinateFrame` subclass
        The coordinate frame class to transform into.
    priority : float or int
        The priority if this transform when finding the shortest
        coordinate transform path - large numbers are lower priorities.
    register_graph : `~astropy.coordinates.TransformGraph` or None
        A graph to register this transformation with on creation, or
        `None` to leave it unregistered.
    """

    def __init__(self, fromsys, tosys, priority=1, register_graph=None):
        if not isinstance(fromsys, type):
            raise TypeError("fromsys must be a class")
        if not isinstance(tosys, type):
            raise TypeError("tosys must be a class")

        self.fromsys = fromsys
        self.tosys = tosys
        self.priority = float(priority)

        if register_graph:
            # this will do the type-checking when it adds to the graph
            self.register(register_graph)
        else:
            if not isinstance(fromsys, type) or not isinstance(tosys, type):
                raise TypeError("fromsys and tosys must be classes")

        self.overlapping_frame_attr_names = overlap = []
        if hasattr(fromsys, "frame_attributes") and hasattr(tosys, "frame_attributes"):
            # the if statement is there so that non-frame things might be usable
            # if it makes sense
            for from_nm in fromsys.frame_attributes:
                if from_nm in tosys.frame_attributes:
                    overlap.append(from_nm)

    def register(self, graph):
        """
        Add this transformation to the requested Transformation graph,
        replacing anything already connecting these two coordinates.

        Parameters
        ----------
        graph : `~astropy.coordinates.TransformGraph` object
            The graph to register this transformation with.
        """
        graph.add_transform(self.fromsys, self.tosys, self)

    def unregister(self, graph):
        """
        Remove this transformation from the requested transformation
        graph.

        Parameters
        ----------
        graph : a TransformGraph object
            The graph to unregister this transformation from.

        Raises
        ------
        ValueError
            If this is not currently in the transform graph.
        """
        graph.remove_transform(self.fromsys, self.tosys, self)

    @abstractmethod
    def __call__(self, fromcoord, toframe):
        """
        Does the actual coordinate transformation from the ``fromsys`` class to
        the ``tosys`` class.

        Parameters
        ----------
        fromcoord : `~astropy.coordinates.BaseCoordinateFrame` subclass instance
            An object of class matching ``fromsys`` that is to be transformed.
        toframe : object
            An object that has the attributes necessary to fully specify the
            frame.  That is, it must have attributes with names that match the
            keys of the dictionary ``tosys.frame_attributes``.
            Typically this is of class ``tosys``, but it *might* be
            some other class as long as it has the appropriate attributes.

        Returns
        -------
        tocoord : `~astropy.coordinates.BaseCoordinateFrame` subclass instance
            The new coordinate after the transform has been applied.
        """
