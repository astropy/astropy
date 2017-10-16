# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Framework and base classes for coordinate frames/"low-level" coordinate
classes.
"""


# Standard library
import abc
import copy
import inspect
from collections import namedtuple, OrderedDict, defaultdict
import warnings

# Dependencies
import numpy as np

# Project
from ..utils.compat.misc import override__dir__
from ..utils.decorators import lazyproperty
from ..utils.exceptions import AstropyWarning
from .. import units as u
from ..utils import (OrderedDescriptorContainer, ShapedLikeNDArray,
                     check_broadcast)
from .transformations import TransformGraph
from . import representation as r

from .attributes import Attribute

# Import old names for Attributes so we don't break backwards-compatibility
# (some users rely on them being here, although that is not encouraged, as this
# is not the public API location -- see attributes.py).
from .attributes import (
    TimeFrameAttribute, QuantityFrameAttribute,
    EarthLocationAttribute, CoordinateAttribute,
    CartesianRepresentationFrameAttribute)  # pylint: disable=W0611


__all__ = ['BaseCoordinateFrame', 'frame_transform_graph',
           'GenericFrame', 'RepresentationMapping']


# the graph used for all transformations between frames
frame_transform_graph = TransformGraph()


def _get_repr_cls(value):
    """
    Return a valid representation class from ``value`` or raise exception.
    """

    if value in r.REPRESENTATION_CLASSES:
        value = r.REPRESENTATION_CLASSES[value]
    elif (not isinstance(value, type) or
          not issubclass(value, r.BaseRepresentation)):
        raise ValueError(
            'Representation is {0!r} but must be a BaseRepresentation class '
            'or one of the string aliases {1}'.format(
                value, list(r.REPRESENTATION_CLASSES)))
    return value


def _get_repr_classes(base, **differentials):
    """Get valid representation and differential classes.

    Parameters
    ----------
    base : str or `~astropy.coordinates.BaseRepresentation` subclass
        class for the representation of the base coordinates.  If a string,
        it is looked up among the known representation classes.
    **differentials : dict of str or `~astropy.coordinates.BaseDifferentials`
        Keys are like for normal differentials, i.e., 's' for a first
        derivative in time, etc.  If an item is set to `None`, it will be
        guessed from the base class.

    Returns
    -------
    repr_classes : dict of subclasses
        The base class is keyed by 'base'; the others by the keys of
        ``diffferentials``.
    """
    base = _get_repr_cls(base)
    repr_classes = {'base': base}

    for name, differential_cls in differentials.items():
        if differential_cls == 'base':
            # We don't want to fail for this case.
            differential_cls = r.DIFFERENTIAL_CLASSES.get(base.get_name(), None)

        elif differential_cls in r.DIFFERENTIAL_CLASSES:
            differential_cls = r.DIFFERENTIAL_CLASSES[differential_cls]

        elif (differential_cls is not None and
              (not isinstance(differential_cls, type) or
               not issubclass(differential_cls, r.BaseDifferential))):
            raise ValueError(
                'Differential is {0!r} but must be a BaseDifferential class '
                'or one of the string aliases {1}'.format(
                    differential_cls, list(r.DIFFERENTIAL_CLASSES)))
        repr_classes[name] = differential_cls
    return repr_classes


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

        if 'default_differential' in members:
            default_diff = members.pop('default_differential')
            found_default_diff = True
        else:
            default_diff = None
            found_default_diff = False

        if 'frame_specific_representation_info' in members:
            repr_info = members.pop('frame_specific_representation_info')
            found_repr_info = True
        else:
            repr_info = None
            found_repr_info = False

        # somewhat hacky, but this is the best way to get the MRO according to
        # https://mail.python.org/pipermail/python-list/2002-December/167861.html
        tmp_cls = super().__new__(mcls, name, bases, members)

        # now look through the whole MRO for the class attributes, raw for
        # frame_attr_names, and leading underscore for others
        for m in (c.__dict__ for c in tmp_cls.__mro__):
            if not found_default_repr and '_default_representation' in m:
                default_repr = m['_default_representation']
                found_default_repr = True

            if not found_default_diff and '_default_differential' in m:
                default_diff = m['_default_differential']
                found_default_diff = True

            if (not found_repr_info and
                    '_frame_specific_representation_info' in m):
                repr_info = m['_frame_specific_representation_info']
                found_repr_info = True

            if found_default_repr and found_default_diff and found_repr_info:
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
        mcls.readonly_prop_factory(members, 'default_differential',
                                   default_diff)
        mcls.readonly_prop_factory(members,
                                   'frame_specific_representation_info',
                                   copy.deepcopy(repr_info))

        # now set the frame name as lower-case class name, if it isn't explicit
        if 'name' not in members:
            members['name'] = name.lower()

        return super().__new__(mcls, name, bases, members)

    @staticmethod
    def readonly_prop_factory(members, attr, value):
        private_attr = '_' + attr

        def getter(self):
            return getattr(self, private_attr)

        members[private_attr] = value
        members[attr] = property(getter)


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
        return super().__new__(cls, reprname, framename, defaultunit)


class BaseCoordinateFrame(ShapedLikeNDArray, metaclass=FrameMeta):
    """
    The base class for coordinate frames.

    This class is intended to be subclassed to create instances of specific
    systems.  Subclasses can implement the following attributes:

    * `default_representation`
        A subclass of `~astropy.coordinates.BaseRepresentation` that will be
        treated as the default representation of this frame.  This is the
        representation assumed by default when the frame is created.

    * `default_differential`
        A subclass of `~astropy.coordinates.BaseDifferential` that will be
        treated as the default differential class of this frame.  This is the
        differential class assumed by default when the frame is created.

    * `~astropy.coordinates.Attribute` class attributes
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
    differential_cls : `BaseDifferential`, dict, optional
        A differential class or dictionary of differential classes (currently
        only a velocity differential with key 's' is supported). This sets
        the expected input differential class, thereby changing the expected
        keyword arguments of the data passed in. For example, passing
        ``differential_cls=CartesianDifferential`` will make the classes
        expect velocity data with the argument names ``v_x, v_y, v_z``.
    copy : bool, optional
        If `True` (default), make copies of the input coordinate arrays.
        Can only be passed in as a keyword argument.
    """

    default_representation = None
    default_differential = None

    # Specifies special names and units for representation and differential
    # attributes.
    frame_specific_representation_info = {}

    _inherit_descriptors_ = (Attribute,)

    frame_attributes = OrderedDict()
    # Default empty frame_attributes dict

    def __init__(self, *args, copy=True, representation=None,
                 differential_cls=None,**kwargs):
        self._attr_names_with_defaults = []

        # TODO: we should be able to deal with an instance, not just a
        # class or string for representation and differential_cls.

        if representation is not None or differential_cls is not None:

            if representation is None:
                representation = self.default_representation

            if (inspect.isclass(differential_cls) and
                    issubclass(differential_cls, r.BaseDifferential)):
                # TODO: assumes the differential class is for the velocity
                # differential
                differential_cls = {'s': differential_cls}

            elif differential_cls is None:
                differential_cls = {'s': 'base'} # see set_representation_cls()

            self.set_representation_cls(representation, **differential_cls)

        # if not set below, this is a frame with no data
        representation_data = None
        differential_data = None

        args = list(args)  # need to be able to pop them
        if (len(args) > 0) and (isinstance(args[0], r.BaseRepresentation) or
                                args[0] is None):
            representation_data = args.pop(0)
            if len(args) > 0:
                raise TypeError(
                    'Cannot create a frame with both a representation and '
                    'other positional arguments')

            if representation_data is not None:
                diffs = representation_data.differentials
                differential_data = diffs.get('s', None)
                if ((differential_data is None and len(diffs) > 0) or
                        (differential_data is not None and len(diffs) > 1)):
                    raise ValueError('Multiple differentials are associated '
                                     'with the representation object passed in '
                                     'to the frame initializer. Only a single '
                                     'velocity differential is supported. Got: '
                                     '{0}'.format(diffs))

        elif self.representation:
            representation_cls = self.representation
            # Get any representation data passed in to the frame initializer
            # using keyword or positional arguments for the component names
            repr_kwargs = {}
            for nmkw, nmrep in self.representation_component_names.items():
                if len(args) > 0:
                    # first gather up positional args
                    repr_kwargs[nmrep] = args.pop(0)
                elif nmkw in kwargs:
                    repr_kwargs[nmrep] = kwargs.pop(nmkw)

            # special-case the Spherical->UnitSpherical if no `distance`
            # TODO: possibly generalize this somehow?
            if repr_kwargs:
                if repr_kwargs.get('distance', True) is None:
                    del repr_kwargs['distance']

                if (issubclass(representation_cls, r.SphericalRepresentation)
                        and 'distance' not in repr_kwargs):
                    representation_cls = representation_cls._unit_representation

                representation_data = representation_cls(copy=copy, **repr_kwargs)

            # Now we handle the Differential data:
            # Get any differential data passed in to the frame initializer
            # using keyword or positional arguments for the component names
            differential_cls = self.get_representation_cls('s')
            diff_component_names = self.get_representation_component_names('s')
            diff_kwargs = {}
            for nmkw, nmrep in diff_component_names.items():
                if len(args) > 0:
                    # first gather up positional args
                    diff_kwargs[nmrep] = args.pop(0)
                elif nmkw in kwargs:
                    diff_kwargs[nmrep] = kwargs.pop(nmkw)

            if diff_kwargs:
                if (hasattr(differential_cls, '_unit_differential') and
                        'd_distance' not in diff_kwargs):
                    differential_cls = differential_cls._unit_differential

                elif len(diff_kwargs) == 1 and 'd_distance' in diff_kwargs:
                    differential_cls = r.RadialDifferential

                differential_data = differential_cls(copy=copy, **diff_kwargs)

        if len(args) > 0:
            raise TypeError(
                '{0}.__init__ had {1} remaining unhandled arguments'.format(
                    self.__class__.__name__, len(args)))

        if representation_data is None and differential_data is not None:
            raise ValueError("Cannot pass in differential component data "
                             "without positional (representation) data.")

        if differential_data:
            self._data = representation_data.with_differentials(
                {'s': differential_data})
        else:
            self._data = representation_data  # possibly None.

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
            # This makes the cache keys backwards-compatible, but also adds
            # support for having differentials attached to the frame data
            # representation object.
            if 's' in self._data.differentials:
                # TODO: assumes a velocity unit differential
                key = (self._data.__class__.__name__,
                       self._data.differentials['s'].__class__.__name__,
                       False)
            else:
                key = (self._data.__class__.__name__, False)

            # Set up representation cache.
            self.cache['representation'][key] = self._data

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
            raise ValueError('The frame object "{0!r}" does not have '
                             'associated data'.format(self))
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

    def __bool__(self):
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

    def get_representation_cls(self, which='base'):
        """The class used for part of this frame's data.

        Parameters
        ----------
        which : ('base', 's', `None`)
            The class of which part to return.  'base' means the class used to
            represent the coordinates; 's' the first derivative to time, i.e.,
            the class representing the proper motion and/or radial velocity.
            If `None`, return a dict with both.

        Returns
        -------
        representation : `~astropy.coordinates.BaseRepresentation` or `~astropy.coordinates.BaseDifferential`.
        """
        if not hasattr(self, '_representation'):
            self._representation = {'base': self.default_representation,
                                    's': self.default_differential}
        return self._representation[which] if which is not None else self._representation

    def set_representation_cls(self, base=None, s='base'):
        """Set representation and/or differential class for this frame's data.

        Parameters
        ----------
        base : str, `~astropy.coordinates.BaseRepresentation` subclass, optional
            The name or subclass to use to represent the coordinate data.
        s : `~astropy.coordinates.BaseDifferential` subclass, optional
            The differential subclass to use to represent any velocities,
            such as proper motion and radial velocity.  If equal to 'base',
            which is the default, it will be inferred from the representation.
            If `None`, the representation will drop any differentials.
        """
        if base is None:
            base = self._representation['base']
        self._representation = _get_repr_classes(base=base, s=s)

    representation = property(
        fget=get_representation_cls, fset=set_representation_cls,
        doc="""The representation class used for this frame's data.

        This will be a subclass from `~astropy.coordinates.BaseRepresentation`.
        Can also be *set* using the string name of the representation. If you
        wish to set an explicit differential class (rather than have it be
        inferred), use the ``set_represenation_cls`` method.
        """)

    @classmethod
    def _get_representation_info(cls):
        # This exists as a class method only to support handling frame inputs
        # without units, which are deprecated and will be removed.  This can be
        # moved into the representation_info property at that time.

        repr_attrs = {}
        for repr_diff_cls in (list(r.REPRESENTATION_CLASSES.values()) +
                              list(r.DIFFERENTIAL_CLASSES.values())):
            repr_attrs[repr_diff_cls] = {'names': [], 'units': []}
            for c in repr_diff_cls.attr_classes.keys():
                repr_attrs[repr_diff_cls]['names'].append(c)
                rec_unit = repr_diff_cls.recommended_units.get(c, None)
                repr_attrs[repr_diff_cls]['units'].append(rec_unit)

        for repr_diff_cls, mappings in cls._frame_specific_representation_info.items():

            if isinstance(repr_diff_cls, str):
                # TODO: this provides a layer of backwards compatibility in
                # case the key is a string, but now we want explicit classes.
                repr_diff_cls = _get_repr_cls(repr_diff_cls)

            # take the 'names' and 'units' tuples from repr_attrs,
            # and then use the RepresentationMapping objects
            # to update as needed for this frame.
            nms = repr_attrs[repr_diff_cls]['names']
            uns = repr_attrs[repr_diff_cls]['units']
            comptomap = dict([(m.reprname, m) for m in mappings])
            for i, c in enumerate(repr_diff_cls.attr_classes.keys()):
                if c in comptomap:
                    mapp = comptomap[c]
                    nms[i] = mapp.framename

                    # need the isinstance because otherwise if it's a unit it
                    # will try to compare to the unit string representation
                    if not (isinstance(mapp.defaultunit, str) and
                            mapp.defaultunit == 'recommended'):
                        uns[i] = mapp.defaultunit
                        # else we just leave it as recommended_units says above

            # Convert to tuples so that this can't mess with frame internals
            repr_attrs[repr_diff_cls]['names'] = tuple(nms)
            repr_attrs[repr_diff_cls]['units'] = tuple(uns)

        return repr_attrs

    @property
    def representation_info(self):
        """
        A dictionary with the information of what attribute names for this frame
        apply to particular representations.
        """
        return self._get_representation_info()

    def get_representation_component_names(self, which='base'):
        out = OrderedDict()
        repr_or_diff_cls = self.get_representation_cls(which)
        if repr_or_diff_cls is None:
            return out
        data_names = repr_or_diff_cls.attr_classes.keys()
        repr_names = self.representation_info[repr_or_diff_cls]['names']
        for repr_name, data_name in zip(repr_names, data_names):
            out[repr_name] = data_name
        return out

    def get_representation_component_units(self, which='base'):
        out = OrderedDict()
        repr_or_diff_cls = self.get_representation_cls(which)
        if repr_or_diff_cls is None:
            return out
        repr_attrs = self.representation_info[repr_or_diff_cls]
        repr_names = repr_attrs['names']
        repr_units = repr_attrs['units']
        for repr_name, repr_unit in zip(repr_names, repr_units):
            if repr_unit:
                out[repr_name] = repr_unit
        return out

    representation_component_names = property(get_representation_component_names)

    representation_component_units = property(get_representation_component_units)

    def replicate(self, copy=False, **kwargs):
        """
        Return a replica of the frame, optionally with new frame attributes.

        The replica is a new frame object that has the same data as this frame
        object and with frame attributes overriden if they are provided as extra
        keyword arguments to this method. If ``copy`` is set to `True` then a
        copy of the internal arrays will be made.  Otherwise the replica will
        use a reference to the original arrays when possible to save memory. The
        internal arrays are normally not changeable by the user so in most cases
        it should not be necessary to set ``copy`` to `True`.

        Parameters
        ----------
        copy : bool, optional
            If True, the resulting object is a copy of the data.  When False,
            references are used where  possible. This rule also applies to the
            frame attributes.

        Any additional keywords are treated as frame attributes to be set on the
        new frame object.

        Returns
        -------
        frameobj : same as this frame
            Replica of this object, but possibly with new frame attributes.
        """
        return self._apply('copy' if copy else 'replicate', **kwargs)

    def replicate_without_data(self, copy=False, **kwargs):
        """
        Return a replica without data, optionally with new frame attributes.

        The replica is a new frame object without data but with the same frame
        attributes as this object, except where overriden by extra keyword
        arguments to this method.  The ``copy`` keyword determines if the frame
        attributes are truly copied vs being references (which saves memory for
        cases where frame attributes are large).

        This method is essentially the converse of `realize_frame`.

        Parameters
        ----------
        copy : bool, optional
            If True, the resulting object has copies of the frame attributes.
            When False, references are used where  possible.

        Any additional keywords are treated as frame attributes to be set on the
        new frame object.

        Returns
        -------
        frameobj : same as this frame
            Replica of this object, but without data and possibly with new frame
            attributes.
        """
        kwargs['_framedata'] = None
        return self._apply('copy' if copy else 'replicate', **kwargs)

    def realize_frame(self, representation):
        """
        Generates a new frame *with new data* from another frame (which may or
        may not have data). Roughly speaking, the converse of
        `replicate_without_data`.

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
        # Here we pass representation_cls=None to _apply, since we do not want
        # to insist that the realized frame has the same representation as
        # self.  [Avoids breaking sunpy; see gh-6208]
        # TODO: should we expose this, so one has a choice?
        return self._apply('replicate', _framedata=representation,
                           representation_cls=None)

    def represent_as(self, base, s='base', in_frame_units=False):
        """
        Generate and return a new representation of this frame's `data`
        as a Representation object.

        Note: In order to make an in-place change of the representation
        of a Frame or SkyCoord object, set the ``representation``
        attribute of that object to the desired new representation, or
        use the ``set_representation_cls`` method to also set the differential.

        Parameters
        ----------
        base : subclass of BaseRepresentation or string
            The type of representation to generate.  Must be a *class*
            (not an instance), or the string name of the representation
            class.
        s : subclass of `~astropy.coordinates.BaseDifferential`, str, optional
            Class in which any velocities should be represented. Must be
            a *class* (not an instance), or the string name of the
            differential class.  If equal to 'base' (default), inferred from
            the base class.  If `None`, all velocity information is dropped.
        in_frame_units : bool, keyword only
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

        # For backwards compatibility (because in_frame_units used to be the
        # 2nd argument), we check to see if `new_differential` is a boolean. If
        # it is, we ignore the value of `new_differential` and warn about the
        # position change
        if isinstance(s, bool):
            warnings.warn("The argument position for `in_frame_units` in "
                          "`represent_as` has changed. Use as a keyword "
                          "argument if needed.", AstropyWarning)
            in_frame_units = s
            s = 'base'

        # In the future, we may want to support more differentials, in which
        # case one probably needs to define **kwargs above and use it here.
        # But for now, we only care about the velocity.
        repr_classes = _get_repr_classes(base=base, s=s)
        representation_cls = repr_classes['base']
        # We only keep velocity information
        if 's' in self.data.differentials:
            differential_cls = repr_classes['s']
        elif s is None or s == 'base':
            differential_cls = None
        else:
            raise TypeError('Frame data has no associated differentials '
                            '(i.e. the frame has no velocity data) - '
                            'represent_as() only accepts a new '
                            'representation.')

        if differential_cls:
            cache_key = (representation_cls.__name__,
                         differential_cls.__name__, in_frame_units)
        else:
            cache_key = (representation_cls.__name__, in_frame_units)

        cached_repr = self.cache['representation'].get(cache_key)
        if not cached_repr:
            if differential_cls:
                # TODO NOTE: only supports a single differential
                data = self.data.represent_as(representation_cls,
                                              differential_cls)
                diff = data.differentials['s']  # TODO: assumes velocity
            else:
                data = self.data.represent_as(representation_cls)

            # If the new representation is known to this frame and has a defined
            # set of names and units, then use that.
            new_attrs = self.representation_info.get(representation_cls)
            if new_attrs and in_frame_units:
                datakwargs = dict((comp, getattr(data, comp))
                                  for comp in data.components)
                for comp, new_attr_unit in zip(data.components, new_attrs['units']):
                    if new_attr_unit:
                        datakwargs[comp] = datakwargs[comp].to(new_attr_unit)
                data = data.__class__(copy=False, **datakwargs)

            if differential_cls:
                # the original differential
                data_diff = self.data.differentials['s']

                # If the new differential is known to this frame and has a
                # defined set of names and units, then use that.
                new_attrs = self.representation_info.get(differential_cls)
                if new_attrs and in_frame_units:
                    diffkwargs = dict((comp, getattr(diff, comp))
                                      for comp in diff.components)
                    for comp, new_attr_unit in zip(diff.components,
                                                   new_attrs['units']):
                        # Some special-casing to treat a situation where the
                        # input data has a UnitSphericalDifferential or a
                        # RadialDifferential. It is re-represented to the
                        # frame's differential class (which might be, e.g., a
                        # dimensional Differential), so we don't want to try to
                        # convert the empty component units
                        if (isinstance(data_diff,
                                       (r.UnitSphericalDifferential,
                                        r.UnitSphericalCosLatDifferential)) and
                                comp not in data_diff.__class__.attr_classes):
                            continue

                        elif (isinstance(data_diff, r.RadialDifferential) and
                              comp not in data_diff.__class__.attr_classes):
                            continue

                        if new_attr_unit and hasattr(diff, comp):
                            diffkwargs[comp] = diffkwargs[comp].to(new_attr_unit)

                    diff = diff.__class__(copy=False, **diffkwargs)

                    # Here we have to bypass using with_differentials() because
                    # it has a validation check. But because .representation and
                    # .differential_cls don't point to the original classes, if
                    # the input differential is a RadialDifferential, it usually
                    # gets turned into a SphericalCosLatDifferential (or
                    # whatever the default is) with strange units for the d_lon
                    # and d_lat attributes. This then causes the dictionary key
                    # check to fail (i.e. comparison against
                    # `diff._get_deriv_key()`)
                    data._differentials.update({'s': diff})
                    # data = data.with_differentials({'s': diff})

            self.cache['representation'][cache_key] = data

        return self.cache['representation'][cache_key]

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

        if (getattr(self.data, 'differentials', None) and
           hasattr(self, 'obstime') and hasattr(new_frame, 'obstime') and
           np.any(self.obstime != new_frame.obstime)):
            raise NotImplementedError('You cannot transform a frame that has '
                                      'velocities to another frame at a '
                                      'different obstime. If you think this '
                                      'should (or should not) be possible, '
                                      'please comment at https://github.com/astropy/astropy/issues/6280')

        if inspect.isclass(new_frame):
            # Use the default frame attributes for this class
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
            if (issubclass(self.representation, r.SphericalRepresentation) and
                    isinstance(self.data, r.UnitSphericalRepresentation)):
                rep_cls = self.data.__class__
            else:
                rep_cls = self.representation

            if 's' in self.data.differentials:
                dif_cls = self.get_representation_cls('s')
                dif_data = self.data.differentials['s']
                if isinstance(dif_data, (r.UnitSphericalDifferential,
                                         r.UnitSphericalCosLatDifferential,
                                         r.RadialDifferential)):
                    dif_cls = dif_data.__class__

            else:
                dif_cls = None

            data = self.represent_as(rep_cls, dif_cls, in_frame_units=True)

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

        if 's' in self.data.differentials:
            data_repr_spl = data_repr.split('\n')
            if 'has differentials' in data_repr_spl[-1]:
                diffrepr = repr(data.differentials['s']).split('\n')
                if diffrepr[0].startswith('<'):
                    diffrepr[0] = ' ' + ' '.join(diffrepr[0].split(' ')[1:])
                for frm_nm, rep_nm in self.get_representation_component_names('s').items():
                    diffrepr[0] = diffrepr[0].replace(rep_nm, frm_nm)
                if diffrepr[-1].endswith('>'):
                    diffrepr[-1] = diffrepr[-1][:-1]
                data_repr_spl[-1] = '\n'.join(diffrepr)

            data_repr = '\n'.join(data_repr_spl)

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
        if '_framedata' in kwargs:
            data = kwargs.pop('_framedata')
        else:
            data = self.data if self.has_data else None

        # This is to provide a slightly nicer error message if the user tries to
        # use frame_obj.representation instead of frame_obj.data to get the
        # underlying representation object [e.g., #2890]
        if inspect.isclass(data):
            raise TypeError('Class passed as data instead of a representation '
                            'instance. If you called frame.representation, this'
                            ' returns the representation class. frame.data '
                            'returns the instantiated object - you may want to '
                            ' use this instead.')

        # TODO: expose this trickery in docstring?
        representation_cls = kwargs.pop('representation_cls',
                                        self.representation)

        differential_cls = kwargs.pop('differential_cls',
                                      self.get_representation_cls('s'))

        def apply_method(value):
            if isinstance(value, ShapedLikeNDArray):
                if method == 'replicate' and not hasattr(value, method):
                    return value  # reference directly
                else:
                    return value._apply(method, *args, **kwargs)
            else:
                if callable(method):
                    return method(value, *args, **kwargs)
                else:
                    if method == 'replicate' and not hasattr(value, method):
                        return value  # reference directly
                    else:
                        return getattr(value, method)(*args, **kwargs)

        if data is not None:
            data = apply_method(data)

        # TODO: change to representation_cls in __init__ - gh-6219.
        frattrs = {'representation': representation_cls,
                   'differential_cls': differential_cls}
        for attr in self.get_frame_attr_names():
            if attr not in self._attr_names_with_defaults:
                if (method == 'copy' or method == 'replicate') and attr in kwargs:
                    value = kwargs[attr]
                else:
                    value = getattr(self, attr)
                if getattr(value, 'size', 1) > 1:
                    value = apply_method(value)
                elif method == 'copy' or method == 'flatten':
                    # flatten should copy also for a single element array, but
                    # we cannot use it directly for array scalars, since it
                    # always returns a one-dimensional array. So, just copy.
                    value = copy.copy(value)

                frattrs[attr] = value

        return self.__class__(data, **frattrs)

    @override__dir__
    def __dir__(self):
        """
        Override the builtin `dir` behavior to include representation
        names.

        TODO: dynamic representation transforms (i.e. include cylindrical et al.).
        """
        dir_values = set(self.representation_component_names)
        dir_values |= set(self.get_representation_component_names('s'))

        return dir_values

    def __getattr__(self, attr):
        """
        Allow access to attributes on the representation and differential as
        found via ``self.get_representation_component_names``.

        TODO: We should handle dynamic representation transforms here (e.g.,
        `.cylindrical`) instead of defining properties as below.
        """

        # attr == '_representation' is likely from the hasattr() test in the
        # representation property which is used for
        # self.representation_component_names.
        #
        # Prevent infinite recursion here.
        if attr.startswith('_'):
            return self.__getattribute__(attr)  # Raise AttributeError.

        repr_names = self.representation_component_names
        if attr in repr_names:
            if self._data is None:
                self.data  # this raises the "no data" error by design - doing it
                # this way means we don't have to replicate the error message here

            rep = self.represent_as(self.representation, in_frame_units=True)
            val = getattr(rep, repr_names[attr])
            return val

        diff_names = self.get_representation_component_names('s')
        if attr in diff_names:
            if self._data is None:
                self.data  # see above.
            # TODO: this doesn't work for the case when there is only
            # unitspherical information. The differential_cls gets set to the
            # default_differential, which expects full information, so the
            # units don't work out
            rep = self.represent_as(in_frame_units=True,
                                    **self.get_representation_cls(None))
            val = getattr(rep.differentials['s'], diff_names[attr])
            return val

        return self.__getattribute__(attr)  # Raise AttributeError.

    def __setattr__(self, attr, value):
        repr_attr_names = set()
        if hasattr(self, 'representation_info'):
            for representation_attr in self.representation_info.values():
                repr_attr_names.update(representation_attr['names'])

        if attr in repr_attr_names:
            raise AttributeError(
                'Cannot set any frame attribute {0}'.format(attr))
        else:
            super().__setattr__(attr, value)

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

        self_unit_sph = self.represent_as(r.UnitSphericalRepresentation)
        other_transformed = other.transform_to(self)
        other_unit_sph = other_transformed.represent_as(r.UnitSphericalRepresentation)

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

        if issubclass(self.data.__class__, r.UnitSphericalRepresentation):
            raise ValueError('This object does not have a distance; cannot '
                             'compute 3d separation.')

        # do this first just in case the conversion somehow creates a distance
        other_in_self_system = other.transform_to(self)

        if issubclass(other_in_self_system.__class__, r.UnitSphericalRepresentation):
            raise ValueError('The other object does not have a distance; '
                             'cannot compute 3d separation.')

        # drop the differentials to ensure they don't do anything odd in the
        # subtraction
        self_car = self.data.without_differentials().represent_as(r.CartesianRepresentation)
        other_car = other_in_self_system.data.without_differentials().represent_as(r.CartesianRepresentation)
        return Distance((self_car - other_car).norm())

    @property
    def cartesian(self):
        """
        Shorthand for a cartesian representation of the coordinates in this
        object.
        """

        # TODO: if representations are updated to use a full transform graph,
        #       the representation aliases should not be hard-coded like this
        return self.represent_as('cartesian', in_frame_units=True)

    @property
    def spherical(self):
        """
        Shorthand for a spherical representation of the coordinates in this
        object.
        """

        # TODO: if representations are updated to use a full transform graph,
        #       the representation aliases should not be hard-coded like this
        return self.represent_as('spherical', in_frame_units=True)

    @property
    def sphericalcoslat(self):
        """
        Shorthand for a spherical representation of the positional data and a
        `SphericalCosLatDifferential` for the velocity data in this object.
        """

        # TODO: if representations are updated to use a full transform graph,
        #       the representation aliases should not be hard-coded like this
        return self.represent_as('spherical', 'sphericalcoslat',
                                 in_frame_units=True)


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
            self.frame_attributes[name] = Attribute(default)
            setattr(self, '_' + name, default)

        super().__init__(None)

    def __getattr__(self, name):
        if '_' + name in self.__dict__:
            return getattr(self, '_' + name)
        else:
            raise AttributeError('no {0}'.format(name))

    def __setattr__(self, name, value):
        if name in self.get_frame_attr_names():
            raise AttributeError("can't set frame attribute '{0}'".format(name))
        else:
            super().__setattr__(name, value)
