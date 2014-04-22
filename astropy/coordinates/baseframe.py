# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Framework and base classes for coordinate frames/"low-level" coordinate classes.
"""
from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

# Standard library
import copy

# Dependencies
import numpy as np

# Project
from ..extern import six
from .. import units as u
from .transformations import TransformGraph

__all__ = ['BaseCoordinateFrame', 'frame_transform_graph']


# the graph used for all transformations between frames
frame_transform_graph = TransformGraph()


class FrameMeta(type):
    def __new__(cls, name, parents, clsdct):
        if clsdct['preferred_representation']:
            # create properties for the preferred_attr_names
            for propnm, reprnm in six.iteritems(clsdct['preferred_attr_names']):
                clsdct[propnm] = property(FrameMeta.repr_getter_factory(reprnm))

            #update the class docstring with the preferred initialization
            #TODO: make this more-correct for the astropy doc style
            init_docstr = """
            The preferred representation is "{0}", allowing kwargs: {1}
            """.format(clsdct['preferred_representation'],
                       clsdct['preferred_attr_names'])
            clsdct['__doc__'] += init_docstr

        return super(FrameMeta, cls).__new__(cls, name, parents, clsdct)

    @staticmethod
    def repr_getter_factory(reprnm):
        def getter(self):
            rep = self.represent_as(self.preferred_representation)
            return getattr(rep, reprnm)
        return getter


@six.add_metaclass(FrameMeta)
class BaseCoordinateFrame(object):
    """ docstring """

    preferred_representation = None
    preferred_attr_names = {}  # maps preferred name to "real" name on repr obj
    frame_attr_names = {}  # maps attribute to default value

    def __init__(self, representation=None, **kwargs):
        from .representation import SphericalRepresentation
        from .representation import UnitSphericalRepresentation

        for fnm, fdefault in six.iteritems(self.frame_attr_names):
            if fnm in kwargs:
                setattr(self, fnm, kwargs.pop(fnm))
            else:
                setattr(self, fnm, fdefault)

        pref_rep = self.preferred_representation

        if pref_rep:
            pref_kwargs = {}
            for nmkw, nmrep in six.iteritems(self.preferred_attr_names):
                if nmkw in kwargs:
                    pref_kwargs[nmrep] = kwargs.pop(nmkw)

            if pref_kwargs:
                if representation:
                    msg = ('Cannot give both a representation object ({0}) and '
                           'arguments for the preferred representation ({1})')
                    raise ValueError(msg.format(representation, pref_kwargs))

                #special-case the Spherical->UnitSpherical if no `distance`
                #TODO: possibly generalize this somehow?

                if (pref_rep == SphericalRepresentation and
                    'distance' not in pref_kwargs):
                    representation = UnitSphericalRepresentation(**pref_kwargs)
                else:
                    representation = pref_rep(**pref_kwargs)

        if kwargs:
            raise TypeError('Coordinate frame got unexpected keywords: ' +
                            str(kwargs.keys()))

        self._data = representation

        if self._data:
            self._rep_cache = dict()
            self._rep_cache[representation.__class__.__name__] = representation

    @property
    def data(self):
        if self._data is None:
            raise AttributeError('The frame object "{0}" does not have '
                                 'associated data'.format(repr(self)))
        return self._data

    @property
    def has_data(self):
        return self._data is not None

    def realize_frame(self, representation):
        """
        Generates a new frame *with data* from a frame without data.

        Parameters
        ----------
        representation : BaseRepresentation
            The representation to use as the data for the new frame.

        Returns
        -------
        frameobj : same as this frame
            A new object with the same frame attributes as this one, but
            with the `representation` as the data.
        """
        frattrs = dict([(nm, getattr(self, nm))for nm in self.frame_attr_names])
        return self.__class__(representation, **frattrs)

    def represent_as(self, new_representation):
        cached_repr = self._rep_cache.get(new_representation.__name__, None)
        if not cached_repr:
            rep = self.data.represent_as(new_representation)
            self._rep_cache[new_representation.__name__] = rep
        return self._rep_cache[new_representation.__name__]

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

        trans = frame_transform_graph.get_transform(self.__class__, new_frame)
        if trans is None:
            if new_frame is self.__class__:
                # no special transform needed
                return copy.deepcopy(self)
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

    def __repr__(self):
        from .representation import SphericalRepresentation
        from .representation import UnitSphericalRepresentation

        content = ['<' + self.__class__.__name__]
        if self.has_data:
            content[-1] += ' coordinate:'
            if self.preferred_representation:
                #special-case the Spherical->UnitSpherical if no `distance`
                #TODO: possibly generalize this somehow?
                prefrep = self.represent_as(self.preferred_representation)
                nodistance = (prefrep.__class__ == SphericalRepresentation and
                             self.data.__class__ == UnitSphericalRepresentation)

                for prefnm, repnm in six.iteritems(self.preferred_attr_names):
                    if nodistance and repnm == 'distance':
                        continue
                    datastr = str(getattr(self._data, repnm))
                    content.append(prefnm + '=' + datastr + ',')
            else:
                content.append(repr(self.data) + ',')
        else:
            content[-1] += ' frame:'

        for nm in self.frame_attr_names:
            content.append(nm + '=' + str(getattr(self, nm)) + ',')

        if content[-1][-1] in (',', ':'):
            #remove trailing comma/colon
            content[-1] = content[-1][:-1]

        return ' '.join(content) + '>'

    def __getitem__(self, view):
        if self.has_data:
            return self.realize_frame(self.data[view])
        else:
            raise ValueError('Cannot index a frame with no data')

    @property
    def cartesian(self):
        """
        Shorthand for a cartesian representation of the coordinates in this object.
        """
        # TODO: if representations are updated to use a full transform graph,
        #       the representation aliases should not be hard-coded like this
        from .representation import CartesianRepresentation

        return self.represent_as(CartesianRepresentation)

    @property
    def spherical(self):
        """
        Shorthand for a spherical representation of the coordinates in this object.
        """
        # TODO: if representations are updated to use a full transform graph,
        #       the representation aliases should not be hard-coded like this
        from .representation import SphericalRepresentation

        return self.represent_as(SphericalRepresentation)
