# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Framework and base classes for coordinate frames/"low-level" coordinate classes.
"""
from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

# Standard library
import copy
import inspect

# Dependencies
import numpy as np

# Project
from ..extern import six
from .. import units as u
from .transformations import TransformGraph

__all__ = ['BaseCoordinateFrame', 'frame_transform_graph', 'GenericFrame']


# the graph used for all transformations between frames
frame_transform_graph = TransformGraph()


class FrameMeta(type):
    def __new__(cls, name, parents, clsdct):
        # somewhat hacky, but this is the best way to get the MRO according to
        # https://mail.python.org/pipermail/python-list/2002-December/167861.html
        mro = super(FrameMeta, cls).__new__(cls, name, parents, clsdct).__mro__
        parent_clsdcts = [c.__dict__ for c in mro]
        parent_clsdcts.insert(0, clsdct)

        #now look through the whole MRO for the relevant class attributes
        pref_repr = None
        pref_attrs = frame_attrs = {}
        found_pref_repr = found_pref_attrs = found_frame_attrs = False
        for clsdcti in parent_clsdcts:
            if not found_pref_repr and 'preferred_representation' in clsdcti:
                pref_repr = clsdcti['preferred_representation']
                found_pref_repr = True
            if not found_pref_attrs and 'preferred_attr_names' in clsdcti:
                pref_attrs = clsdcti['preferred_attr_names']
                found_pref_attrs = True
            if not found_frame_attrs and 'frame_attr_names' in clsdcti:
                frame_attrs = clsdcti['frame_attr_names']
                found_frame_attrs = True

            if found_pref_repr and found_pref_attrs and found_frame_attrs:
                break
        else:
            raise ValueError('Could not find the expected BaseCoordinateFrame '
                             'class attributes.  Are you mis-using FrameMeta?')


        if pref_repr:
            # create properties for the preferred_attr_names
            for propnm, reprnm in six.iteritems(pref_attrs):
                clsdct[propnm] = property(FrameMeta.repr_getter_factory(reprnm))

        #also make properties for the frame_attr_names to make them immutible
        #after creation
        for attrnm in frame_attrs:
            if attrnm in clsdct:
                if isinstance(clsdct[attrnm], property):
                    #means the property was defined already... so trust the
                    #subclasser and don't make a new one
                    continue
                else:
                    raise ValueError("A non-property exists that's *also* in frame_attr_names")
            clsdct[attrnm] = property(FrameMeta.frame_attr_factory(attrnm))

        return super(FrameMeta, cls).__new__(cls, name, parents, clsdct)

    @staticmethod
    def repr_getter_factory(reprnm):
        def getter(self):
            rep = self.represent_as(self.preferred_representation)
            return getattr(rep, reprnm)
        return getter

    @staticmethod
    def frame_attr_factory(attrnm):
        def getter(self):
            return getattr(self, '_'+ attrnm)
        return getter


@six.add_metaclass(FrameMeta)
class BaseCoordinateFrame(object):
    """
    The base class for coordinate frames.

    This class is intended to be subclassed to create instances of specific
    systems.  Subclasses can implement the following attributes.

    * `preferred_representation`
        A subclass of `~astropy.coordinates.BaseRepresentation` that will be
        treated as the "standard" representation of this frame, or None to have
        no special representation.

    * `preferred_attr_names`
        A dictionary mapping attribute names to be created on *this* class to
        names of attributes on the `preferred_representation`.  If
        `preferred_representation` is None, this does nothing.

    * `frame_attr_names`
        A dictionary with keys that are the additional attributes necessary to
        specify the frame, and values that are the default values of those
        attributes.
    """

    preferred_representation = None
    preferred_attr_names = {}  # maps preferred name to "real" name on repr obj
    frame_attr_names = {}  # maps attribute to default value

    def __init__(self, *args, **kwargs):
        from .representation import SphericalRepresentation, \
                                    UnitSphericalRepresentation, \
                                    BaseRepresentation

        representation = None  # if not set below, this is a frame with no data

        for fnm, fdefault in six.iteritems(self.frame_attr_names):

            # read-only properties for these attributes are made in the
            # metaclass  so we set the 'real' attrbiutes as the name prefaced
            # with an underscore

            if fnm in kwargs:
                setattr(self, '_' + fnm, kwargs.pop(fnm))
            else:
                setattr(self, '_' + fnm, fdefault)

        pref_rep = self.preferred_representation
        args = list(args)  # need to be able to pop them
        if (len(args) > 0) and (isinstance(args[0], BaseRepresentation) or
                                args[0] is None):
            representation = args.pop(0)
            if len(args) > 0:
                raise TypeError('Cannot create a frame with both a '
                                'representation and other positional arguments')

        elif pref_rep:
            pref_kwargs = {}
            for nmkw, nmrep in six.iteritems(self.preferred_attr_names):
                if len(args) > 0:
                    #first gather up positional args
                    pref_kwargs[nmrep] = args.pop(0)
                elif nmkw in kwargs:
                    pref_kwargs[nmrep] = kwargs.pop(nmkw)

            #special-case the Spherical->UnitSpherical if no `distance`
            #TODO: possibly generalize this somehow?

            if pref_kwargs:
                if (pref_rep == SphericalRepresentation and
                    'distance' not in pref_kwargs):
                    print('pan',self.preferred_attr_names)
                    representation = UnitSphericalRepresentation(**pref_kwargs)
                else:
                    print('pan2',self.preferred_attr_names)
                    representation = pref_rep(**pref_kwargs)

        if len(args) > 0:
            raise TypeError(self.__class__.__name__ + '.__init__ had {0} '
                            'remaining unprocessed arguments'.format(len(args)))
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
        new_frame : class or frame object
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

        if self._data is None:
            raise ValueError('Cannot transform a frame with no data')

        if inspect.isclass(new_frame):
            #means use the defaults for this class
            new_frame = new_frame()

        trans = frame_transform_graph.get_transform(self.__class__,
                                                    new_frame.__class__)
        if trans is None:
            if new_frame is self.__class__:
                # no special transform needed
                return copy.deepcopy(self)
            msg = 'Cannot transform from {0} to {1}'
            raise ConvertError(msg.format(self.__class__, new_frame.__class__))
        return trans(self, new_frame)

    def is_transformable_to(self, new_frame):
        """
        Determines if this coordinate frame can be transformed to another
        given frame.

        Parameters
        ----------
        new_frame : class or frame object
            The proposed frame to transform into.

        Returns
        -------
        transformable : bool or str
            `True` if this can be transformed to ``new_frame``, `False` if
            not, or the string 'same' if ``new_frame`` is the same system as
            this object but no transformation is defined.

        Notes
        -----
        A return value of 'same' means the transformation will work, but it will
        just give back a copy of this object.  The intended usage is::

            if coord.is_transformable_to(some_unknown_frame):
                coord2 = coord.transform_to(some_unknown_frame)

        This will work even if `some_unknown_frame`  turns out to be the same
        frame class as `coord`.  This is intended for cases where the frame is
        the same regardless of the frame attributes (e.g. ICRS), but be aware
        that it *might* also indicate that someone forgot to define the
        transformation between two objects of the same frame class but with
        different attributes.

        """

        new_frame_cls = new_frame if inspect.isclass(new_frame) else new_frame.__class__
        trans = frame_transform_graph.get_transform(self.__class__, new_frame_cls)

        if trans is None:
            if new_frame_cls is self.__class__:
                return 'same'
            else:
                return False
        else:
            return True

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
                    datastr = str(getattr(prefrep, repnm))
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

    def separation(self, other):
        """
        Computes on-sky separation between this coordinate and another.

        Parameters
        ----------
        other : `~astropy.coordinates.SphericalCoordinatesBase`
            The coordinate to get the separation to.

        Returns
        -------
        sep : `~astropy.coordinates.Angle`
            The on-sky separation between this and the ``other`` coordinate.

        Notes
        -----
        The separation is calculated using the Vincenty formula, which
        is stable at all locations, including poles and antipodes [1]_.

        .. [1] http://en.wikipedia.org/wiki/Great-circle_distance

        """
        from .representation import UnitSphericalRepresentation
        from .angle_utilities import angular_separation
        from .angles import Angle

        self_unit_sph = self.represent_as(UnitSphericalRepresentation)
        other_unit_sph = other.transform_to(self.__class__).represent_as(UnitSphericalRepresentation)

        # Get the separation as a Quantity, convert to Angle in degrees
        sep = angular_separation(self_unit_sph.lon, self_unit_sph.lat,
                                 other_unit_sph.lon, other_unit_sph.lat)
        return Angle(sep, unit=u.degree)

    def separation_3d(self, other):
        """
        Computes three dimensional separation between this coordinate
        and another.

        Parameters
        ----------
        other : `~astropy.coordinates.SphericalCoordinatesBase`
            The coordinate system to get the distance to.

        Returns
        -------
        sep : `~astropy.coordinates.Distance`
            The real-space distance between these two coordinates.

        Raises
        ------
        ValueError
            If this or the other coordinate do not have distances.
        """
        from .representation import UnitSphericalRepresentation
        from .distances import Distance

        if self.data.__class__ == UnitSphericalRepresentation:
            raise ValueError('This object does not have a distance; cannot '
                             'compute 3d separation.')

        # do this first just in case the conversion somehow creates a distance
        other_in_self_system = other.transform_to(self.__class__)

        if other_in_self_system.__class__ == UnitSphericalRepresentation:
            raise ValueError('The other object does not have a distance; '
                             'cannot compute 3d separation.')

        dx = self.cartesian.x - other_in_self_system.cartesian.x
        dy = self.cartesian.y - other_in_self_system.cartesian.y
        dz = self.cartesian.z - other_in_self_system.cartesian.z

        distval = (dx.value ** 2 + dy.value ** 2 + dz.value ** 2) ** 0.5
        return Distance(distval, dx.unit)

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


class GenericFrame(BaseCoordinateFrame):
    """
    A frame object that can't store data but can hold any arbitrary frame
    attributes. Mostly useful as a utility for the high-level class to store
    intermediate frame attributes.

    Parameters
    ----------
    frame_attrs : dict
        A dictionary of attributes to be used as the frame attributes for this
        frame.
    """

    def __init__(self, frame_attrs):
        super(GenericFrame, self).__init__(None)

        self.frame_attr_names = frame_attrs
        for attrnm, attrval in six.iteritems(frame_attrs):
            setattr(self, attrnm, attrval)

    def __setattr__(self, name, value):
        if name in self.frame_attr_names:
            raise AttributeError("can't set frame attribute '{0}'".format(name))
        else:
            super(GenericFrame, self).__setattr__(name, value)
