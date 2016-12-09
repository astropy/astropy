# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Framework and base classes for coordinate frames/"low-level" coordinate
classes.
"""

from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

# Standard library
import abc
import copy
import inspect
from collections import namedtuple, OrderedDict, defaultdict

# Dependencies
import numpy as np

# Project
from ..utils.compat.misc import override__dir__
from ..utils.compat.numpy import broadcast_to as np_broadcast_to
from ..utils.decorators import lazyproperty
from ..extern import six
from ..extern.six.moves import zip
from .. import units as u
from ..utils import (OrderedDescriptor, OrderedDescriptorContainer,
                     ShapedLikeNDArray, check_broadcast)
from .transformations import TransformGraph
from .representation import (BaseRepresentation, CartesianRepresentation,
                             SphericalRepresentation,
                             UnitSphericalRepresentation,
                             REPRESENTATION_CLASSES)


__all__ = ['BaseCoordinateFrame', 'frame_transform_graph', 'GenericFrame',
           'FrameAttribute', 'TimeFrameAttribute', 'QuantityFrameAttribute',
           'EarthLocationAttribute', 'RepresentationMapping',
           'CartesianRepresentationFrameAttribute', 'CoordinateAttribute']


# the graph used for all transformations between frames
frame_transform_graph = TransformGraph()


def _get_repr_cls(value):
    """
    Return a valid representation class from ``value`` or raise exception.
    """

    if value in REPRESENTATION_CLASSES:
        value = REPRESENTATION_CLASSES[value]
    try:
        # value might not be a class, so use try
        assert issubclass(value, BaseRepresentation)
    except (TypeError, AssertionError):
        raise ValueError(
            'Representation is {0!r} but must be a BaseRepresentation class '
            'or one of the string aliases {1}'.format(
                value, list(REPRESENTATION_CLASSES)))
    return value


# Need to subclass ABCMeta as well, so that this meta class can be combined
# with ShapedLikeNDArray below (which is an ABC); without it, one gets
# "TypeError: metaclass conflict: the metaclass of a derived class must be a
#  (non-strict) subclass of the metaclasses of all its bases"
class FrameMeta(OrderedDescriptorContainer, abc.ABCMeta):
    def __new__(mcls, name, bases, members):
        if 'default_representation' in members:
            default_repr = members.pop('default_representation')
            found_default_repr = True
        else:
            default_repr = None
            found_default_repr = False

        if 'frame_specific_representation_info' in members:
            repr_info = members.pop('frame_specific_representation_info')
            found_repr_info = True
        else:
            repr_info = None
            found_repr_info = False

        # somewhat hacky, but this is the best way to get the MRO according to
        # https://mail.python.org/pipermail/python-list/2002-December/167861.html
        tmp_cls = super(FrameMeta, mcls).__new__(mcls, name, bases, members)

        # now look through the whole MRO for the class attributes, raw for
        # frame_attr_names, and leading underscore for others
        for m in (c.__dict__ for c in tmp_cls.__mro__):
            if not found_default_repr and '_default_representation' in m:
                default_repr = m['_default_representation']
                found_default_repr = True
            if (not found_repr_info and
                    '_frame_specific_representation_info' in m):
                repr_info = m['_frame_specific_representation_info']
                found_repr_info = True

            if found_default_repr and found_repr_info:
                break
        else:
            raise ValueError(
                'Could not find all expected BaseCoordinateFrame class '
                'attributes.  Are you mis-using FrameMeta?')

        # Make read-only properties for the frame class attributes that should
        # be read-only to make them immutable after creation.
        # We copy attributes instead of linking to make sure there's no
        # accidental cross-talk between classes
        mcls.readonly_prop_factory(members, 'default_representation',
                                   default_repr)
        mcls.readonly_prop_factory(members,
                                   'frame_specific_representation_info',
                                   copy.deepcopy(repr_info))

        # now set the frame name as lower-case class name, if it isn't explicit
        if 'name' not in members:
            members['name'] = name.lower()

        return super(FrameMeta, mcls).__new__(mcls, name, bases, members)

    @staticmethod
    def readonly_prop_factory(members, attr, value):
        private_attr = '_' + attr

        def getter(self):
            return getattr(self, private_attr)

        members[private_attr] = value
        members[attr] = property(getter)


class FrameAttribute(OrderedDescriptor):
    """A non-mutable data descriptor to hold a frame attribute.

    This class must be used to define frame attributes (e.g. ``equinox`` or
    ``obstime``) that are included in a frame class definition.

    Examples
    --------
    The `~astropy.coordinates.FK4` class uses the following class attributes::

      class FK4(BaseCoordinateFrame):
          equinox = TimeFrameAttribute(default=_EQUINOX_B1950)
          obstime = TimeFrameAttribute(default=None,
                                       secondary_attribute='equinox')

    This means that ``equinox`` and ``obstime`` are available to be set as
    keyword arguments when creating an ``FK4`` class instance and are then
    accessible as instance attributes.  The instance value for the attribute
    must be stored in ``'_' + <attribute_name>`` by the frame ``__init__``
    method.

    Note in this example that ``equinox`` and ``obstime`` are time attributes
    and use the ``TimeAttributeFrame`` class.  This subclass overrides the
    ``convert_input`` method to validate and convert inputs into a ``Time``
    object.

    Parameters
    ----------
    default : object
        Default value for the attribute if not provided
    secondary_attribute : str
        Name of a secondary instance attribute which supplies the value if
        ``default is None`` and no value was supplied during initialization.
    """

    _class_attribute_ = 'frame_attributes'
    _name_attribute_ = 'name'
    name = '<unbound>'

    def __init__(self, default=None, secondary_attribute=''):
        self.default = default
        self.secondary_attribute = secondary_attribute
        super(FrameAttribute, self).__init__()

    def convert_input(self, value):
        """
        Validate the input ``value`` and convert to expected attribute class.

        The base method here does nothing, but subclasses can implement this
        as needed.  The method should catch any internal exceptions and raise
        ValueError with an informative message.

        The method returns the validated input along with a boolean that
        indicates whether the input value was actually converted.  If the input
        value was already the correct type then the ``converted`` return value
        should be ``False``.

        Parameters
        ----------
        value : object
            Input value to be converted.

        Returns
        -------
        output_value
            The ``value`` converted to the correct type (or just ``value`` if
            ``converted`` is False)
        converted : bool
            True if the conversion was actually performed, False otherwise.

        Raises
        ------
        ValueError
            If the input is not valid for this attribute.
        """
        return value, False

    def __get__(self, instance, frame_cls=None):
        if instance is None:
            out = self.default
        else:
            out = getattr(instance, '_' + self.name, self.default)
            if out is None:
                out = getattr(instance, self.secondary_attribute, self.default)

        out, converted = self.convert_input(out)
        if instance is not None:
            instance_shape = getattr(instance, 'shape', None)
            if instance_shape is not None and (getattr(out, 'size', 1) > 1 and
                                               out.shape != instance_shape):
                # If the shapes do not match, try broadcasting.
                try:
                    if isinstance(out, ShapedLikeNDArray):
                        out = out._apply(np_broadcast_to, shape=instance_shape,
                                         subok=True)
                    else:
                        out = np_broadcast_to(out, instance_shape, subok=True)
                except ValueError:
                    # raise more informative exception.
                    raise ValueError(
                        "attribute {0} should be scalar or have shape {1}, "
                        "but is has shape {2} and could not be broadcast."
                        .format(self.name, instance_shape, out.shape))

                converted = True

            if converted:
                setattr(instance, '_' + self.name, out)

        return out

    def __set__(self, instance, val):
        raise AttributeError('Cannot set frame attribute')


class TimeFrameAttribute(FrameAttribute):
    """
    Frame attribute descriptor for quantities that are Time objects.
    See the `~astropy.coordinates.FrameAttribute` API doc for further
    information.

    Parameters
    ----------
    default : object
        Default value for the attribute if not provided
    secondary_attribute : str
        Name of a secondary instance attribute which supplies the value if
        ``default is None`` and no value was supplied during initialization.
    """

    def convert_input(self, value):
        """
        Convert input value to a Time object and validate by running through
        the Time constructor.  Also check that the input was a scalar.

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

        from ..time import Time

        if value is None:
            return None, False

        if isinstance(value, Time):
            out = value
            converted = False
        else:
            try:
                out = Time(value)
            except Exception as err:
                raise ValueError(
                    'Invalid time input {0}={1!r}\n{2}'.format(self.name,
                                                               value, err))
            converted = True

        return out, converted


class CartesianRepresentationFrameAttribute(FrameAttribute):
    """
    A frame attribute that is a CartesianRepresentation with specified units.

    Parameters
    ----------
    default : object
        Default value for the attribute if not provided
    secondary_attribute : str
        Name of a secondary instance attribute which supplies the value if
        ``default is None`` and no value was supplied during initialization.
    unit : unit object or None
        Name of a unit that the input will be converted into. If None, no
        unit-checking or conversion is performed
    """
    def __init__(self, default=None, secondary_attribute='', unit=None):
        super(CartesianRepresentationFrameAttribute, self).__init__(default, secondary_attribute)
        self.unit = unit

    def convert_input(self, value):
        """
        Checks that the input is a CartesianRepresentation with the correct
        unit, or the special value ``[0, 0, 0]``.

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
        if (isinstance(value, list) and len(value) == 3 and
                all(v == 0 for v in value) and self.unit is not None):
            return CartesianRepresentation(np.zeros(3) * self.unit), True
        else:
            # is it a CartesianRepresentation with correct unit?
            if hasattr(value, 'xyz') and value.xyz.unit == self.unit:
                return value, False

            converted = True
            # if it's a CartesianRepresentation, get the xyz Quantity
            value = getattr(value, 'xyz', value)
            if not hasattr(value, 'unit'):
                raise TypeError('tried to set a CartesianRepresentationFrameAttribute with '
                                'something that does not have a unit.')

            value = value.to(self.unit)

            # now try and make a CartesianRepresentation.
            cartrep = CartesianRepresentation(value, copy=False)
            return cartrep, converted


class QuantityFrameAttribute(FrameAttribute):
    """
    A frame attribute that is a quantity with specified units and shape
    (optionally).

    Parameters
    ----------
    default : object
        Default value for the attribute if not provided
    secondary_attribute : str
        Name of a secondary instance attribute which supplies the value if
        ``default is None`` and no value was supplied during initialization.
    unit : unit object or None
        Name of a unit that the input will be converted into. If None, no
        unit-checking or conversion is performed
    shape : tuple or None
        If given, specifies the shape the attribute must be
    """
    def __init__(self, default=None, secondary_attribute='', unit=None, shape=None):
        super(QuantityFrameAttribute, self).__init__(default, secondary_attribute)
        self.unit = unit
        self.shape = shape

    def convert_input(self, value):
        """
        Checks that the input is a Quantity with the necessary units (or the
        special value ``0``).

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
        if np.all(value == 0) and self.unit is not None:
            return u.Quantity(np.zeros(self.shape), self.unit), True
        else:
            if not hasattr(value, 'unit'):
                raise TypeError('Tried to set a QuantityFrameAttribute with '
                                'something that does not have a unit.')
            oldvalue = value
            value = u.Quantity(oldvalue, self.unit, copy=False)
            if self.shape is not None and value.shape != self.shape:
                raise ValueError('The provided value has shape "{0}", but '
                                 'should have shape "{1}"'.format(value.shape,
                                                                  self.shape))
            converted = oldvalue is not value
            return value, converted


class EarthLocationAttribute(FrameAttribute):
    """
    A frame attribute that can act as a `~astropy.coordinates.EarthLocation`.
    It can be created as anything that can be transformed to the
    `~astropy.coordinates.ITRS` frame, but always presents as an `EarthLocation`
    when accessed after creation.

    Parameters
    ----------
    default : object
        Default value for the attribute if not provided
    secondary_attribute : str
        Name of a secondary instance attribute which supplies the value if
        ``default is None`` and no value was supplied during initialization.
    """

    def convert_input(self, value):
        """
        Checks that the input is a Quantity with the necessary units (or the
        special value ``0``).

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
        elif isinstance(value, EarthLocation):
            return value, False
        else:
            #we have to do the import here because of some tricky circular deps
            from .builtin_frames import ITRS

            if not hasattr(value, 'transform_to'):
                raise ValueError('"{0}" was passed into an '
                                 'EarthLocationAttribute, but it does not have '
                                 '"transform_to" method'.format(value))
            itrsobj = value.transform_to(ITRS)
            return itrsobj.earth_location, True


class CoordinateAttribute(FrameAttribute):
    """
    A frame attribute which is a coordinate object. It can be given as a
    low-level frame class *or* a `~astropy.coordinates.SkyCoord`, but will
    always be converted to the low-level frame class when accessed.

    Parameters
    ----------
    frame : a coordinate frame class
        The type of frame this attribute can be
    default : object
        Default value for the attribute if not provided
    secondary_attribute : str
        Name of a secondary instance attribute which supplies the value if
        ``default is None`` and no value was supplied during initialization.
    """
    def __init__(self, frame, default=None, secondary_attribute=''):
        self._frame = frame
        super(CoordinateAttribute, self).__init__(default, secondary_attribute)

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
                raise ValueError('"{0}" was passed into a '
                                 'CoordinateAttribute, but it does not have '
                                 '"transform_to" method'.format(value))
            transformedobj = value.transform_to(self._frame)
            if hasattr(transformedobj, 'frame'):
                transformedobj = transformedobj.frame
            return transformedobj, True


_RepresentationMappingBase = \
    namedtuple('RepresentationMapping',
               ('reprname', 'framename', 'defaultunit'))


class RepresentationMapping(_RepresentationMappingBase):
    """
    This `~collections.namedtuple` is used with the
    ``frame_specific_representation_info`` attribute to tell frames what
    attribute names (and default units) to use for a particular representation.
    ``reprname`` and ``framename`` should be strings, while ``defaultunit`` can
    be either an astropy unit, the string ``'recommended'`` (to use whatever
    the representation's ``recommended_units`` is), or None (to indicate that
    no unit mapping should be done).
    """

    def __new__(cls, reprname, framename, defaultunit='recommended'):
        # this trick just provides some defaults
        return super(RepresentationMapping, cls).__new__(cls, reprname,
                                                         framename,
                                                         defaultunit)

@six.add_metaclass(FrameMeta)
class BaseCoordinateFrame(ShapedLikeNDArray):
    """
    The base class for coordinate frames.

    This class is intended to be subclassed to create instances of specific
    systems.  Subclasses can implement the following attributes:

    * `default_representation`
        A subclass of `~astropy.coordinates.BaseRepresentation` that will be
        treated as the default representation of this frame.  This is the
        representation assumed by default when the frame is created.

    * `~astropy.coordinates.FrameAttribute` class attributes
       Frame attributes such as ``FK4.equinox`` or ``FK4.obstime`` are defined
       using a descriptor class.  See the narrative documentation or
       built-in classes code for details.

    * `frame_specific_representation_info`
        A dictionary mapping the name or class of a representation to a list of
        `~astropy.coordinates.RepresentationMapping` objects that tell what
        names and default units should be used on this frame for the components
        of that representation.

    Parameters
    ----------
    representation : `BaseRepresentation` or None
        A representation object or `None` to have no data (or use the other
        arguments)
    *args, **kwargs
        Coordinates, with names that depend on the subclass.
    copy : bool, optional
        If `True` (default), make copies of the input coordinate arrays.
        Can only be passed in as a keyword argument.
    """

    default_representation = None
    # specifies special names/units for representation attributes
    frame_specific_representation_info = {}

    _inherit_descriptors_ = (FrameAttribute,)

    frame_attributes = OrderedDict()
    # Default empty frame_attributes dict

    def __init__(self, *args, **kwargs):
        copy = kwargs.pop('copy', True)
        self._attr_names_with_defaults = []

        if 'representation' in kwargs:
            self.representation = kwargs.pop('representation')

        # if not set below, this is a frame with no data
        representation_data = None

        args = list(args)  # need to be able to pop them
        if (len(args) > 0) and (isinstance(args[0], BaseRepresentation) or
                                args[0] is None):
            representation_data = args.pop(0)
            if len(args) > 0:
                raise TypeError(
                    'Cannot create a frame with both a representation and '
                    'other positional arguments')

        elif self.representation:
            repr_kwargs = {}
            for nmkw, nmrep in self.representation_component_names.items():
                if len(args) > 0:
                    #first gather up positional args
                    repr_kwargs[nmrep] = args.pop(0)
                elif nmkw in kwargs:
                    repr_kwargs[nmrep] = kwargs.pop(nmkw)

            #special-case the Spherical->UnitSpherical if no `distance`
            #TODO: possibly generalize this somehow?

            if repr_kwargs:
                if repr_kwargs.get('distance', True) is None:
                    del repr_kwargs['distance']
                if (issubclass(self.representation, SphericalRepresentation) and
                        'distance' not in repr_kwargs):
                    representation = self.representation._unit_representation
                else:
                    representation = self.representation
                representation_data = representation(copy=copy, **repr_kwargs)

        if len(args) > 0:
            raise TypeError(
                '{0}.__init__ had {1} remaining unhandled arguments'.format(
                    self.__class__.__name__, len(args)))

        self._data = representation_data

        values = {}
        for fnm, fdefault in self.get_frame_attr_names().items():
            # Read-only frame attributes are defined as FrameAttribue
            # descriptors which are not settable, so set 'real' attributes as
            # the name prefaced with an underscore.

            if fnm in kwargs:
                value = kwargs.pop(fnm)
                setattr(self, '_' + fnm, value)
                # Validate attribute by getting it.  If the instance has data,
                # this also checks its shape is OK.  If not, we do it below.
                values[fnm] = getattr(self, fnm)
            else:
                setattr(self, '_' + fnm, fdefault)
                self._attr_names_with_defaults.append(fnm)

        if kwargs:
            raise TypeError(
                'Coordinate frame got unexpected keywords: {0}'.format(
                    list(kwargs)))

        # We do ``is None`` because self._data might evaluate to false for
        # empty arrays or data == 0
        if self._data is None:
            # No data: we still need to check that any non-scalar attributes
            # have consistent shapes. Collect them for all attributes with
            # size > 1 (which should be array-like and thus have a shape).
            shapes = {fnm: value.shape for fnm, value in values.items()
                      if getattr(value, 'size', 1) > 1}
            if shapes:
                if len(shapes) > 1:
                    try:
                        self._no_data_shape = check_broadcast(*shapes.values())
                    except ValueError:
                        raise ValueError(
                            "non-scalar attributes with inconsistent "
                            "shapes: {0}".format(shapes))

                    # Above, we checked that it is possible to broadcast all
                    # shapes.  By getting and thus validating the attributes,
                    # we verify that the attributes can in fact be broadcast.
                    for fnm in shapes:
                        getattr(self, fnm)
                else:
                    self._no_data_shape = shapes.popitem()[1]

            else:
                self._no_data_shape = ()
        else:
            # Set up representation cache.
            self.cache['representation'][self._data.__class__.__name__, False] = self._data

    @lazyproperty
    def cache(self):
        """
        Cache for this frame, a dict.  It stores anything that should be
        computed from the coordinate data (*not* from the frame attributes).
        This can be used in functions to store anything that might be
        expensive to compute but might be re-used by some other function.
        E.g.::

            if 'user_data' in myframe.cache:
                data = myframe.cache['user_data']
            else:
                myframe.cache['user_data'] = data = expensive_func(myframe.lat)

        If in-place modifications are made to the frame data, the cache should
        be cleared::

            myframe.cache.clear()

        """
        return defaultdict(dict)

    @property
    def data(self):
        """
        The coordinate data for this object.  If this frame has no data, an
        `ValueError` will be raised.  Use `has_data` to
        check if data is present on this frame object.
        """
        if self._data is None:
            raise ValueError('The frame object "{0}" does not have associated '
                             'data'.format(repr(self)))
        return self._data

    @property
    def has_data(self):
        """
        True if this frame has `data`, False otherwise.
        """
        return self._data is not None

    @property
    def shape(self):
        return self.data.shape if self.has_data else self._no_data_shape

    # We have to override the ShapedLikeNDArray definitions, since our shape
    # does not have to be that of the data.
    def __len__(self):
        return len(self.data)

    def __nonzero__(self):  # Py 2.x
        return self.has_data and self.size > 0

    def __bool__(self):  # Py 3.x
        return self.has_data and self.size > 0

    @property
    def size(self):
        return self.data.size

    @property
    def isscalar(self):
        return self.has_data and self.data.isscalar

    @classmethod
    def get_frame_attr_names(cls):
        return OrderedDict((name, getattr(cls, name))
                           for name in cls.frame_attributes)

    @property
    def representation(self):
        """
        The representation of the data in this frame, as a class that is
        subclassed from `~astropy.coordinates.BaseRepresentation`.  Can
        also be *set* using the string name of the representation.
        """
        if not hasattr(self, '_representation'):
            self._representation = self.default_representation
        return self._representation

    @representation.setter
    def representation(self, value):
        self._representation = _get_repr_cls(value)

    @classmethod
    def _get_representation_info(cls):
        # This exists as a class method only to support handling frame inputs
        # without units, which are deprecated and will be removed.  This can be
        # moved into the representation_info property at that time.

        repr_attrs = {}
        for repr_cls in REPRESENTATION_CLASSES.values():
            repr_attrs[repr_cls] = {'names': [], 'units': []}
            for c in repr_cls.attr_classes.keys():
                repr_attrs[repr_cls]['names'].append(c)
                rec_unit = repr_cls.recommended_units.get(c, None)
                repr_attrs[repr_cls]['units'].append(rec_unit)

        for repr_cls, mappings in cls._frame_specific_representation_info.items():
            # keys may be a class object or a name
            repr_cls = _get_repr_cls(repr_cls)

            # take the 'names' and 'units' tuples from repr_attrs,
            # and then use the RepresentationMapping objects
            # to update as needed for this frame.
            nms = repr_attrs[repr_cls]['names']
            uns = repr_attrs[repr_cls]['units']
            comptomap = dict([(m.reprname, m) for m in mappings])
            for i, c in enumerate(repr_cls.attr_classes.keys()):
                if c in comptomap:
                    mapp = comptomap[c]
                    nms[i] = mapp.framename
                    # need the isinstance because otherwise if it's a unit it
                    # will try to compare to the unit string representation
                    if not (isinstance(mapp.defaultunit, six.string_types) and
                            mapp.defaultunit == 'recommended'):
                        uns[i] = mapp.defaultunit
                        # else we just leave it as recommended_units says above
            # Convert to tuples so that this can't mess with frame internals
            repr_attrs[repr_cls]['names'] = tuple(nms)
            repr_attrs[repr_cls]['units'] = tuple(uns)

        return repr_attrs

    @property
    def representation_info(self):
        """
        A dictionary with the information of what attribute names for this frame
        apply to particular representations.
        """
        return self._get_representation_info()

    @property
    def representation_component_names(self):
        out = OrderedDict()
        if self.representation is None:
            return out
        data_names = self.representation.attr_classes.keys()
        repr_names = self.representation_info[self.representation]['names']
        for repr_name, data_name in zip(repr_names, data_names):
            out[repr_name] = data_name
        return out

    @property
    def representation_component_units(self):
        out = OrderedDict()
        if self.representation is None:
            return out
        repr_attrs = self.representation_info[self.representation]
        repr_names = repr_attrs['names']
        repr_units = repr_attrs['units']
        for repr_name, repr_unit in zip(repr_names, repr_units):
            if repr_unit:
                out[repr_name] = repr_unit
        return out

    def realize_frame(self, representation):
        """
        Generates a new frame *with new data* from another frame (which may or
        may not have data).

        Parameters
        ----------
        representation : BaseRepresentation
            The representation to use as the data for the new frame.

        Returns
        -------
        frameobj : same as this frame
            A new object with the same frame attributes as this one, but
            with the ``representation`` as the data.
        """
        frattrs = dict([(attr, getattr(self, attr))
                        for attr in self.get_frame_attr_names()
                        if attr not in self._attr_names_with_defaults])
        return self.__class__(representation, **frattrs)

    def represent_as(self, new_representation, in_frame_units=False):
        """
        Generate and return a new representation of this frame's `data`
        as a Representation object.

        Note: In order to make an in-place change of the representation
        of a Frame or SkyCoord object, set the ``representation``
        attribute of that object to the desired new representation.

        Parameters
        ----------
        new_representation : subclass of BaseRepresentation or string
            The type of representation to generate.  May be a *class*
            (not an instance), or the string name of the representation
            class.

        in_frame_units : bool
            Force the representation units to match the specified units
            particular to this frame

        Returns
        -------
        newrep : BaseRepresentation-derived object
            A new representation object of this frame's `data`.

        Raises
        ------
        AttributeError
            If this object had no `data`

        Examples
        --------
        >>> from astropy import units as u
        >>> from astropy.coordinates import SkyCoord, CartesianRepresentation
        >>> coord = SkyCoord(0*u.deg, 0*u.deg)
        >>> coord.represent_as(CartesianRepresentation)
        <CartesianRepresentation (x, y, z) [dimensionless]
                ( 1.,  0.,  0.)>

        >>> coord.representation = CartesianRepresentation
        >>> coord
        <SkyCoord (ICRS): (x, y, z) [dimensionless]
            ( 1.,  0.,  0.)>
        """
        new_representation = _get_repr_cls(new_representation)

        cached_repr = self.cache['representation'].get((new_representation.__name__,
                                                        in_frame_units))
        if not cached_repr:
            data = self.data.represent_as(new_representation)

            # If the new representation is known to this frame and has a defined
            # set of names and units, then use that.
            new_attrs = self.representation_info.get(new_representation)
            if new_attrs and in_frame_units:
                datakwargs = dict((comp, getattr(data, comp))
                                  for comp in data.components)
                for comp, new_attr_unit in zip(data.components, new_attrs['units']):
                    if new_attr_unit:
                        datakwargs[comp] = datakwargs[comp].to(new_attr_unit)
                data = data.__class__(copy=False, **datakwargs)

            self.cache['representation'][new_representation.__name__, in_frame_units] = data

        return self.cache['representation'][new_representation.__name__, in_frame_units]

    def transform_to(self, new_frame):
        """
        Transform this object's coordinate data to a new frame.

        Parameters
        ----------
        new_frame : class or frame object or SkyCoord object
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

        if hasattr(new_frame, '_sky_coord_frame'):
            # Input new_frame is not a frame instance or class and is most
            # likely a SkyCoord object.
            new_frame = new_frame._sky_coord_frame

        trans = frame_transform_graph.get_transform(self.__class__,
                                                    new_frame.__class__)
        if trans is None:
            if new_frame is self.__class__:
                # no special transform needed, but should update frame info
                return new_frame.realize_frame(self.data)
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

        This will work even if ``some_unknown_frame``  turns out to be the same
        frame class as ``coord``.  This is intended for cases where the frame
        is the same regardless of the frame attributes (e.g. ICRS), but be
        aware that it *might* also indicate that someone forgot to define the
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

    def is_frame_attr_default(self, attrnm):
        """
        Determine whether or not a frame attribute has its value because it's
        the default value, or because this frame was created with that value
        explicitly requested.

        Parameters
        ----------
        attrnm : str
            The name of the attribute to check.

        Returns
        -------
        isdefault : bool
            True if the attribute ``attrnm`` has its value by default, False if
            it was specified at creation of this frame.
        """
        return attrnm in self._attr_names_with_defaults

    def is_equivalent_frame(self, other):
        """
        Checks if this object is the same frame as the ``other`` object.

        To be the same frame, two objects must be the same frame class and have
        the same frame attributes.  Note that it does *not* matter what, if any,
        data either object has.

        Parameters
        ----------
        other : BaseCoordinateFrame
            the other frame to check

        Returns
        -------
        isequiv : bool
            True if the frames are the same, False if not.

        Raises
        ------
        TypeError
            If ``other`` isn't a `BaseCoordinateFrame` or subclass.
        """
        if self.__class__ == other.__class__:
            for frame_attr_name in self.get_frame_attr_names():
                if np.any(getattr(self, frame_attr_name) !=
                          getattr(other, frame_attr_name)):
                    return False
            return True
        elif not isinstance(other, BaseCoordinateFrame):
            raise TypeError("Tried to do is_equivalent_frame on something that "
                            "isn't a frame")
        else:
            return False

    def __repr__(self):
        frameattrs = self._frame_attrs_repr()
        data_repr = self._data_repr()

        if frameattrs:
            frameattrs = ' ({0})'.format(frameattrs)

        if data_repr:
            return '<{0} Coordinate{1}: {2}>'.format(self.__class__.__name__,
                                                     frameattrs, data_repr)
        else:
            return '<{0} Frame{1}>'.format(self.__class__.__name__,
                                            frameattrs)

    def _data_repr(self):
        """Returns a string representation of the coordinate data."""

        if not self.has_data:
            return ''

        if self.representation:
            if (issubclass(self.representation, SphericalRepresentation) and
                    isinstance(self.data, UnitSphericalRepresentation)):
                data = self.represent_as(self.data.__class__,
                                         in_frame_units=True)
            else:
                data = self.represent_as(self.representation,
                                         in_frame_units=True)

            data_repr = repr(data)
            for nmpref, nmrepr in self.representation_component_names.items():
                data_repr = data_repr.replace(nmrepr, nmpref)

        else:
            data = self.data
            data_repr = repr(self.data)

        if data_repr.startswith('<' + data.__class__.__name__):
            # remove both the leading "<" and the space after the name, as well
            # as the trailing ">"
            data_repr = data_repr[(len(data.__class__.__name__) + 2):-1]
        else:
            data_repr = 'Data:\n' + data_repr

        return data_repr

    def _frame_attrs_repr(self):
        """
        Returns a string representation of the frame's attributes, if any.
        """
        return ', '.join([attrnm + '=' + str(getattr(self, attrnm))
                          for attrnm in self.get_frame_attr_names()])

    def _apply(self, method, *args, **kwargs):
        """Create a new instance, applying a method to the underlying data.

        In typical usage, the method is any of the shape-changing methods for
        `~numpy.ndarray` (``reshape``, ``swapaxes``, etc.), as well as those
        picking particular elements (``__getitem__``, ``take``, etc.), which
        are all defined in `~astropy.utils.misc.ShapedLikeNDArray`. It will be
        applied to the underlying arrays in the representation (e.g., ``x``,
        ``y``, and ``z`` for `~astropy.coordinates.CartesianRepresentation`),
        as well as to any frame attributes that have a shape, with the results
        used to create a new instance.

        Internally, it is also used to apply functions to the above parts
        (in particular, `~numpy.broadcast_to`).

        Parameters
        ----------
        method : str or callable
            If str, it is the name of a method that is applied to the internal
            ``components``. If callable, the function is applied.
        args : tuple
            Any positional arguments for ``method``.
        kwargs : dict
            Any keyword arguments for ``method``.
        """
        def apply_method(value):
            if isinstance(value, ShapedLikeNDArray):
                return value._apply(method, *args, **kwargs)
            else:
                if callable(method):
                    return method(value, *args, **kwargs)
                else:
                    return getattr(value, method)(*args, **kwargs)


        data = apply_method(self.data) if self.has_data else None

        frattrs = {}
        for attr in self.get_frame_attr_names():
            if attr not in self._attr_names_with_defaults:
                value = getattr(self, attr)
                if getattr(value, 'size', 1) > 1:
                    value = apply_method(value)
                elif method == 'copy' or method == 'flatten':
                    # flatten should copy also for a single element array, but
                    # we cannot use it directly for array scalars, since it
                    # always returns a one-dimensional array. So, just copy.
                    value = copy.copy(value)

                frattrs[attr] = value

        out = self.__class__(data, **frattrs)
        out.representation = self.representation
        return out

    @override__dir__
    def __dir__(self):
        """
        Override the builtin `dir` behavior to include representation
        names.

        TODO: dynamic representation transforms (i.e. include cylindrical et al.).
        """
        dir_values = set(self.representation_component_names)

        return dir_values

    def __getattr__(self, attr):
        """
        Allow access to attributes defined in
        ``self.representation_component_names``.

        TODO: dynamic representation transforms (i.e. include cylindrical et
        al.).
        """

        # attr == '_representation' is likely from the hasattr() test in the
        # representation property which is used for
        # self.representation_component_names.
        #
        # Prevent infinite recursion here.
        if (attr == '_representation' or
                attr not in self.representation_component_names):
            raise AttributeError("'{0}' object has no attribute '{1}'"
                                 .format(self.__class__.__name__, attr))
        elif self._data is None:
            # The second clause of the ``or`` above means that
            # ``attr in self.representation_component_names``, otherwise we
            # would repeat that check here, because we only want to hit this
            # clause when the user tried to access one of the representation's
            # components

            self.data  # this raises the "no data" error by design - doing it
            # this way means we don't have to replicate the error message here

        else:
            rep = self.represent_as(self.representation, in_frame_units=True)
            val = getattr(rep, self.representation_component_names[attr])
            return val

    def __setattr__(self, attr, value):
        repr_attr_names = set()
        if hasattr(self, 'representation_info'):
            for representation_attr in self.representation_info.values():
                repr_attr_names.update(representation_attr['names'])
        if attr in repr_attr_names:
            raise AttributeError(
                'Cannot set any frame attribute {0}'.format(attr))
        else:
            super(BaseCoordinateFrame, self).__setattr__(attr, value)

    def separation(self, other):
        """
        Computes on-sky separation between this coordinate and another.

        Parameters
        ----------
        other : `~astropy.coordinates.BaseCoordinateFrame`
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
        from .angle_utilities import angular_separation
        from .angles import Angle

        self_unit_sph = self.represent_as(UnitSphericalRepresentation)
        other_transformed = other.transform_to(self)
        other_unit_sph = other_transformed.represent_as(UnitSphericalRepresentation)

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
        other : `~astropy.coordinates.BaseCoordinateFrame`
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

        from .distances import Distance

        if issubclass(self.data.__class__, UnitSphericalRepresentation):
            raise ValueError('This object does not have a distance; cannot '
                             'compute 3d separation.')

        # do this first just in case the conversion somehow creates a distance
        other_in_self_system = other.transform_to(self)

        if issubclass(other_in_self_system.__class__, UnitSphericalRepresentation):
            raise ValueError('The other object does not have a distance; '
                             'cannot compute 3d separation.')

        return Distance((self.cartesian -
                         other_in_self_system.cartesian).norm())

    @property
    def cartesian(self):
        """
        Shorthand for a cartesian representation of the coordinates in this
        object.
        """

        # TODO: if representations are updated to use a full transform graph,
        #       the representation aliases should not be hard-coded like this
        return self.represent_as(CartesianRepresentation, in_frame_units=True)

    @property
    def spherical(self):
        """
        Shorthand for a spherical representation of the coordinates in this object.
        """

        # TODO: if representations are updated to use a full transform graph,
        #       the representation aliases should not be hard-coded like this
        return self.represent_as(SphericalRepresentation, in_frame_units=True)


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

    name = None  # it's not a "real" frame so it doesn't have a name

    def __init__(self, frame_attrs):
        self.frame_attributes = OrderedDict()
        for name, default in frame_attrs.items():
            self.frame_attributes[name] = FrameAttribute(default)
            setattr(self, '_' + name, default)

        super(GenericFrame, self).__init__(None)

    def __getattr__(self, name):
        if '_' + name in self.__dict__:
            return getattr(self, '_' + name)
        else:
            raise AttributeError('no {0}'.format(name))

    def __setattr__(self, name, value):
        if name in self.get_frame_attr_names():
            raise AttributeError("can't set frame attribute '{0}'".format(name))
        else:
            super(GenericFrame, self).__setattr__(name, value)

# doing this import at the bottom prevents a circular import issue that is
# otherwise present due to EarthLocation needing to import ITRS
from .earth import EarthLocation
