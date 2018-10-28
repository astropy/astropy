"""
In this module, we define the coordinate representation classes, which are
used to represent low-level cartesian, spherical, cylindrical, and other
coordinates.
"""


import abc
import functools
import operator
from collections import OrderedDict
import inspect
import warnings

import numpy as np
import astropy.units as u

from .angles import Angle, Longitude, Latitude
from .distances import Distance
from .._erfa import ufunc as erfa_ufunc
from ..utils import ShapedLikeNDArray, classproperty

from ..utils import deprecated_attribute
from ..utils.exceptions import AstropyDeprecationWarning
from ..utils.misc import InheritDocstrings
from ..utils.compat import NUMPY_LT_1_14

__all__ = ["BaseRepresentationOrDifferential", "BaseRepresentation",
           "CartesianRepresentation", "SphericalRepresentation",
           "UnitSphericalRepresentation", "RadialRepresentation",
           "PhysicsSphericalRepresentation", "CylindricalRepresentation",
           "BaseDifferential", "CartesianDifferential",
           "BaseSphericalDifferential", "BaseSphericalCosLatDifferential",
           "SphericalDifferential", "SphericalCosLatDifferential",
           "UnitSphericalDifferential", "UnitSphericalCosLatDifferential",
           "RadialDifferential", "CylindricalDifferential",
           "PhysicsSphericalDifferential"]

# Module-level dict mapping representation string alias names to classes.
# This is populated by the metaclass init so all representation and differential
# classes get registered automatically.
REPRESENTATION_CLASSES = {}
DIFFERENTIAL_CLASSES = {}

# a hash for the content of the above two dicts, cached for speed.
_REPRDIFF_HASH = None
def get_reprdiff_cls_hash():
    """
    Returns a hash value that should be invariable if the
    `REPRESENTATION_CLASSES` and `DIFFERENTIAL_CLASSES` dictionaries have not
    changed.
    """
    global _REPRDIFF_HASH
    if _REPRDIFF_HASH is None:
        _REPRDIFF_HASH = (hash(tuple(REPRESENTATION_CLASSES.items())) +
                          hash(tuple(DIFFERENTIAL_CLASSES.items())) )
    return _REPRDIFF_HASH


def _invalidate_reprdiff_cls_hash():
    global _REPRDIFF_HASH
    _REPRDIFF_HASH = None



# recommended_units deprecation message; if the attribute is removed later,
# also remove its use in BaseFrame._get_representation_info.
_recommended_units_deprecation = """
The 'recommended_units' attribute is deprecated since 3.0 and may be removed
in a future version. Its main use, of representing angles in degrees in frames,
is now done automatically in frames. Further overrides are discouraged but can
be done using a frame's ``frame_specific_representation_info``.
"""


def _array2string(values, prefix=''):
    # Work around version differences for array2string.
    kwargs = {'separator': ', ', 'prefix': prefix}
    kwargs['formatter'] = {}
    if NUMPY_LT_1_14:  # in 1.14, style is no longer used (and deprecated)
        kwargs['style'] = repr

    return np.array2string(values, **kwargs)


def _combine_xyz(x, y, z, xyz_axis=0):
    """
    Combine components ``x``, ``y``, ``z`` into a single Quantity array.

    Parameters
    ----------
    x, y, z : `~astropy.units.Quantity`
        The individual x, y, and z components.
    xyz_axis : int, optional
        The axis in the final array along which the x, y, z components
        should be stored (default: 0).

    Returns
    -------
    xyz : `~astropy.units.Quantity`
        With dimension 3 along ``xyz_axis``, i.e., using the default of ``0``,
        the shape will be ``(3,) + x.shape``.
    """
    # Get x, y, z to the same units (this is very fast for identical units)
    # since np.stack cannot deal with quantity.
    cls = x.__class__
    unit = x.unit
    x = x.value
    y = y.to_value(unit)
    z = z.to_value(unit)

    xyz = np.stack([x, y, z], axis=xyz_axis)
    return cls(xyz, unit=unit, copy=False)


class BaseRepresentationOrDifferential(ShapedLikeNDArray):
    """3D coordinate representations and differentials.

    Parameters
    ----------
    comp1, comp2, comp3 : `~astropy.units.Quantity` or subclass
        The components of the 3D point or differential.  The names are the
        keys and the subclasses the values of the ``attr_classes`` attribute.
    copy : bool, optional
        If `True` (default), arrays will be copied rather than referenced.
    """

    # Ensure multiplication/division with ndarray or Quantity doesn't lead to
    # object arrays.
    __array_priority__ = 50000

    def __init__(self, *args, **kwargs):
        # make argument a list, so we can pop them off.
        args = list(args)
        components = self.components
        attrs = []
        for component in components:
            try:
                attrs.append(args.pop(0) if args else kwargs.pop(component))
            except KeyError:
                raise TypeError('__init__() missing 1 required positional '
                                'argument: {0!r}'.format(component))

        copy = args.pop(0) if args else kwargs.pop('copy', True)

        if args:
            raise TypeError('unexpected arguments: {0}'.format(args))

        if kwargs:
            for component in components:
                if component in kwargs:
                    raise TypeError("__init__() got multiple values for "
                                    "argument {0!r}".format(component))

            raise TypeError('unexpected keyword arguments: {0}'.format(kwargs))

        # Pass attributes through the required initializing classes.
        attrs = [self.attr_classes[component](attr, copy=copy)
                 for component, attr in zip(components, attrs)]
        try:
            attrs = np.broadcast_arrays(*attrs, subok=True)
        except ValueError:
            if len(components) <= 2:
                c_str = ' and '.join(components)
            else:
                c_str = ', '.join(components[:2]) + ', and ' + components[2]
            raise ValueError("Input parameters {0} cannot be broadcast"
                             .format(c_str))
        # Set private attributes for the attributes. (If not defined explicitly
        # on the class, the metaclass will define properties to access these.)
        for component, attr in zip(components, attrs):
            setattr(self, '_' + component, attr)

    @classmethod
    def get_name(cls):
        """Name of the representation or differential.

        In lower case, with any trailing 'representation' or 'differential'
        removed. (E.g., 'spherical' for
        `~astropy.coordinates.SphericalRepresentation` or
        `~astropy.coordinates.SphericalDifferential`.)
        """
        name = cls.__name__.lower()

        if name.endswith('representation'):
            name = name[:-14]
        elif name.endswith('differential'):
            name = name[:-12]

        return name

    # The two methods that any subclass has to define.
    @classmethod
    @abc.abstractmethod
    def from_cartesian(cls, other):
        """Create a representation of this class from a supplied Cartesian one.

        Parameters
        ----------
        other : `CartesianRepresentation`
            The representation to turn into this class

        Returns
        -------
        representation : object of this class
            A new representation of this class's type.
        """
        # Note: the above docstring gets overridden for differentials.
        raise NotImplementedError()

    @abc.abstractmethod
    def to_cartesian(self):
        """Convert the representation to its Cartesian form.

        Note that any differentials get dropped.

        Returns
        -------
        cartrepr : `CartesianRepresentation`
            The representation in Cartesian form.
        """
        # Note: the above docstring gets overridden for differentials.
        raise NotImplementedError()

    @property
    def components(self):
        """A tuple with the in-order names of the coordinate components."""
        return tuple(self.attr_classes)

    def _apply(self, method, *args, **kwargs):
        """Create a new representation or differential with ``method`` applied
        to the component data.

        In typical usage, the method is any of the shape-changing methods for
        `~numpy.ndarray` (``reshape``, ``swapaxes``, etc.), as well as those
        picking particular elements (``__getitem__``, ``take``, etc.), which
        are all defined in `~astropy.utils.misc.ShapedLikeNDArray`. It will be
        applied to the underlying arrays (e.g., ``x``, ``y``, and ``z`` for
        `~astropy.coordinates.CartesianRepresentation`), with the results used
        to create a new instance.

        Internally, it is also used to apply functions to the components
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
        if callable(method):
            apply_method = lambda array: method(array, *args, **kwargs)
        else:
            apply_method = operator.methodcaller(method, *args, **kwargs)

        new = super().__new__(self.__class__)
        for component in self.components:
            setattr(new, '_' + component,
                    apply_method(getattr(self, component)))
        return new

    @property
    def shape(self):
        """The shape of the instance and underlying arrays.

        Like `~numpy.ndarray.shape`, can be set to a new shape by assigning a
        tuple.  Note that if different instances share some but not all
        underlying data, setting the shape of one instance can make the other
        instance unusable.  Hence, it is strongly recommended to get new,
        reshaped instances with the ``reshape`` method.

        Raises
        ------
        AttributeError
            If the shape of any of the components cannot be changed without the
            arrays being copied.  For these cases, use the ``reshape`` method
            (which copies any arrays that cannot be reshaped in-place).
        """
        return getattr(self, self.components[0]).shape

    @shape.setter
    def shape(self, shape):
        # We keep track of arrays that were already reshaped since we may have
        # to return those to their original shape if a later shape-setting
        # fails. (This can happen since coordinates are broadcast together.)
        reshaped = []
        oldshape = self.shape
        for component in self.components:
            val = getattr(self, component)
            if val.size > 1:
                try:
                    val.shape = shape
                except AttributeError:
                    for val2 in reshaped:
                        val2.shape = oldshape
                    raise
                else:
                    reshaped.append(val)

    # Required to support multiplication and division, and defined by the base
    # representation and differential classes.
    @abc.abstractmethod
    def _scale_operation(self, op, *args):
        raise NotImplementedError()

    def __mul__(self, other):
        return self._scale_operation(operator.mul, other)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __truediv__(self, other):
        return self._scale_operation(operator.truediv, other)

    def __div__(self, other):  # pragma: py2
        return self._scale_operation(operator.truediv, other)

    def __neg__(self):
        return self._scale_operation(operator.neg)

    # Follow numpy convention and make an independent copy.
    def __pos__(self):
        return self.copy()

    # Required to support addition and subtraction, and defined by the base
    # representation and differential classes.
    @abc.abstractmethod
    def _combine_operation(self, op, other, reverse=False):
        raise NotImplementedError()

    def __add__(self, other):
        return self._combine_operation(operator.add, other)

    def __radd__(self, other):
        return self._combine_operation(operator.add, other, reverse=True)

    def __sub__(self, other):
        return self._combine_operation(operator.sub, other)

    def __rsub__(self, other):
        return self._combine_operation(operator.sub, other, reverse=True)

    # The following are used for repr and str
    @property
    def _values(self):
        """Turn the coordinates into a record array with the coordinate values.

        The record array fields will have the component names.
        """
        coo_items = [(c, getattr(self, c)) for c in self.components]
        result = np.empty(self.shape, [(c, coo.dtype) for c, coo in coo_items])
        for c, coo in coo_items:
            result[c] = coo.value
        return result

    @property
    def _units(self):
        """Return a dictionary with the units of the coordinate components."""
        return dict([(component, getattr(self, component).unit)
                     for component in self.components])

    @property
    def _unitstr(self):
        units_set = set(self._units.values())
        if len(units_set) == 1:
            unitstr = units_set.pop().to_string()
        else:
            unitstr = '({0})'.format(
                ', '.join([self._units[component].to_string()
                           for component in self.components]))
        return unitstr

    def __str__(self):
        return '{0} {1:s}'.format(_array2string(self._values), self._unitstr)

    def __repr__(self):
        prefixstr = '    '
        arrstr = _array2string(self._values, prefix=prefixstr)

        diffstr = ''
        if getattr(self, 'differentials', None):
            diffstr = '\n (has differentials w.r.t.: {0})'.format(
                ', '.join([repr(key) for key in self.differentials.keys()]))

        unitstr = ('in ' + self._unitstr) if self._unitstr else '[dimensionless]'
        return '<{0} ({1}) {2:s}\n{3}{4}{5}>'.format(
            self.__class__.__name__, ', '.join(self.components),
            unitstr, prefixstr, arrstr, diffstr)


def _make_getter(component):
    """Make an attribute getter for use in a property.

    Parameters
    ----------
    component : str
        The name of the component that should be accessed.  This assumes the
        actual value is stored in an attribute of that name prefixed by '_'.
    """
    # This has to be done in a function to ensure the reference to component
    # is not lost/redirected.
    component = '_' + component

    def get_component(self):
        return getattr(self, component)
    return get_component


# Need to also subclass ABCMeta rather than type, so that this meta class can
# be combined with a ShapedLikeNDArray subclass (which is an ABC).  Without it:
# "TypeError: metaclass conflict: the metaclass of a derived class must be a
#  (non-strict) subclass of the metaclasses of all its bases"
class MetaBaseRepresentation(InheritDocstrings, abc.ABCMeta):
    def __init__(cls, name, bases, dct):
        super().__init__(name, bases, dct)

        # Register representation name (except for BaseRepresentation)
        if cls.__name__ == 'BaseRepresentation':
            return

        if 'attr_classes' not in dct:
            raise NotImplementedError('Representations must have an '
                                      '"attr_classes" class attribute.')

        if 'recommended_units' in dct:
            warnings.warn(_recommended_units_deprecation,
                          AstropyDeprecationWarning)
            # Ensure we don't override the property that warns about the
            # deprecation, but that the value remains the same.
            dct.setdefault('_recommended_units', dct.pop('recommended_units'))

        repr_name = cls.get_name()

        if repr_name in REPRESENTATION_CLASSES:
            raise ValueError("Representation class {0} already defined"
                             .format(repr_name))

        REPRESENTATION_CLASSES[repr_name] = cls
        _invalidate_reprdiff_cls_hash()

        # define getters for any component that does not yet have one.
        for component in cls.attr_classes:
            if not hasattr(cls, component):
                setattr(cls, component,
                        property(_make_getter(component),
                                 doc=("The '{0}' component of the points(s)."
                                      .format(component))))


class BaseRepresentation(BaseRepresentationOrDifferential,
                         metaclass=MetaBaseRepresentation):
    """Base for representing a point in a 3D coordinate system.

    Parameters
    ----------
    comp1, comp2, comp3 : `~astropy.units.Quantity` or subclass
        The components of the 3D points.  The names are the keys and the
        subclasses the values of the ``attr_classes`` attribute.
    differentials : dict, `BaseDifferential`, optional
        Any differential classes that should be associated with this
        representation. The input must either be a single `BaseDifferential`
        subclass instance, or a dictionary with keys set to a string
        representation of the SI unit with which the differential (derivative)
        is taken. For example, for a velocity differential on a positional
        representation, the key would be ``'s'`` for seconds, indicating that
        the derivative is a time derivative.
    copy : bool, optional
        If `True` (default), arrays will be copied rather than referenced.

    Notes
    -----
    All representation classes should subclass this base representation class,
    and define an ``attr_classes`` attribute, an `~collections.OrderedDict`
    which maps component names to the class that creates them. They must also
    define a ``to_cartesian`` method and a ``from_cartesian`` class method. By
    default, transformations are done via the cartesian system, but classes
    that want to define a smarter transformation path can overload the
    ``represent_as`` method. If one wants to use an associated differential
    class, one should also define ``unit_vectors`` and ``scale_factors``
    methods (see those methods for details).
    """

    recommended_units = deprecated_attribute('recommended_units', since='3.0')
    _recommended_units = {}

    def __init__(self, *args, differentials=None, **kwargs):
        # Handle any differentials passed in.
        super().__init__(*args, **kwargs)
        self._differentials = self._validate_differentials(differentials)

    def _validate_differentials(self, differentials):
        """
        Validate that the provided differentials are appropriate for this
        representation and recast/reshape as necessary and then return.

        Note that this does *not* set the differentials on
        ``self._differentials``, but rather leaves that for the caller.
        """

        # Now handle the actual validation of any specified differential classes
        if differentials is None:
            differentials = dict()

        elif isinstance(differentials, BaseDifferential):
            # We can't handle auto-determining the key for this combo
            if (isinstance(differentials, RadialDifferential) and
                    isinstance(self, UnitSphericalRepresentation)):
                raise ValueError("To attach a RadialDifferential to a "
                                 "UnitSphericalRepresentation, you must supply "
                                 "a dictionary with an appropriate key.")

            key = differentials._get_deriv_key(self)
            differentials = {key: differentials}

        for key in differentials:
            try:
                diff = differentials[key]
            except TypeError:
                raise TypeError("'differentials' argument must be a "
                                "dictionary-like object")

            diff._check_base(self)

            if (isinstance(diff, RadialDifferential) and
                    isinstance(self, UnitSphericalRepresentation)):
                # We trust the passing of a key for a RadialDifferential
                # attached to a UnitSphericalRepresentation because it will not
                # have a paired component name (UnitSphericalRepresentation has
                # no .distance) to automatically determine the expected key
                pass

            else:
                expected_key = diff._get_deriv_key(self)
                if key != expected_key:
                    raise ValueError("For differential object '{0}', expected "
                                     "unit key = '{1}' but received key = '{2}'"
                                     .format(repr(diff), expected_key, key))

            # For now, we are very rigid: differentials must have the same shape
            # as the representation. This makes it easier to handle __getitem__
            # and any other shape-changing operations on representations that
            # have associated differentials
            if diff.shape != self.shape:
                # TODO: message of IncompatibleShapeError is not customizable,
                #       so use a valueerror instead?
                raise ValueError("Shape of differentials must be the same "
                                 "as the shape of the representation ({0} vs "
                                 "{1})".format(diff.shape, self.shape))

        return differentials

    def _raise_if_has_differentials(self, op_name):
        """
        Used to raise a consistent exception for any operation that is not
        supported when a representation has differentials attached.
        """
        if self.differentials:
            raise TypeError("Operation '{0}' is not supported when "
                            "differentials are attached to a {1}."
                            .format(op_name, self.__class__.__name__))

    @property
    def _compatible_differentials(self):
        return [DIFFERENTIAL_CLASSES[self.get_name()]]

    @property
    def differentials(self):
        """A dictionary of differential class instances.

        The keys of this dictionary must be a string representation of the SI
        unit with which the differential (derivative) is taken. For example, for
        a velocity differential on a positional representation, the key would be
        ``'s'`` for seconds, indicating that the derivative is a time
        derivative.
        """
        return self._differentials

    # We do not make unit_vectors and scale_factors abstract methods, since
    # they are only necessary if one also defines an associated Differential.
    # Also, doing so would break pre-differential representation subclasses.
    def unit_vectors(self):
        r"""Cartesian unit vectors in the direction of each component.

        Given unit vectors :math:`\hat{e}_c` and scale factors :math:`f_c`,
        a change in one component of :math:`\delta c` corresponds to a change
        in representation of :math:`\delta c \times f_c \times \hat{e}_c`.

        Returns
        -------
        unit_vectors : dict of `CartesianRepresentation`
            The keys are the component names.
        """
        raise NotImplementedError("{} has not implemented unit vectors"
                                  .format(type(self)))

    def scale_factors(self):
        r"""Scale factors for each component's direction.

        Given unit vectors :math:`\hat{e}_c` and scale factors :math:`f_c`,
        a change in one component of :math:`\delta c` corresponds to a change
        in representation of :math:`\delta c \times f_c \times \hat{e}_c`.

        Returns
        -------
        scale_factors : dict of `~astropy.units.Quantity`
            The keys are the component names.
        """
        raise NotImplementedError("{} has not implemented scale factors."
                                  .format(type(self)))

    def _re_represent_differentials(self, new_rep, differential_class):
        """Re-represent the differentials to the specified classes.

        This returns a new dictionary with the same keys but with the
        attached differentials converted to the new differential classes.
        """
        if differential_class is None:
            return dict()

        if not self.differentials and differential_class:
            raise ValueError("No differentials associated with this "
                             "representation!")

        elif (len(self.differentials) == 1 and
                inspect.isclass(differential_class) and
                issubclass(differential_class, BaseDifferential)):
            # TODO: is there a better way to do this?
            differential_class = {
                list(self.differentials.keys())[0]: differential_class
            }

        elif set(differential_class.keys()) != set(self.differentials.keys()):
            ValueError("Desired differential classes must be passed in "
                       "as a dictionary with keys equal to a string "
                       "representation of the unit of the derivative "
                       "for each differential stored with this "
                       "representation object ({0})"
                       .format(self.differentials))

        new_diffs = dict()
        for k in self.differentials:
            diff = self.differentials[k]
            try:
                new_diffs[k] = diff.represent_as(differential_class[k],
                                                 base=self)
            except Exception:
                if (differential_class[k] not in
                        new_rep._compatible_differentials):
                    raise TypeError("Desired differential class {0} is not "
                                    "compatible with the desired "
                                    "representation class {1}"
                                    .format(differential_class[k],
                                            new_rep.__class__))
                else:
                    raise

        return new_diffs

    def represent_as(self, other_class, differential_class=None):
        """Convert coordinates to another representation.

        If the instance is of the requested class, it is returned unmodified.
        By default, conversion is done via cartesian coordinates.

        Parameters
        ----------
        other_class : `~astropy.coordinates.BaseRepresentation` subclass
            The type of representation to turn the coordinates into.
        differential_class : dict of `~astropy.coordinates.BaseDifferential`, optional
            Classes in which the differentials should be represented.
            Can be a single class if only a single differential is attached,
            otherwise it should be a `dict` keyed by the same keys as the
            differentials.
        """
        if other_class is self.__class__ and not differential_class:
            return self.without_differentials()

        else:
            if isinstance(other_class, str):
                raise ValueError("Input to a representation's represent_as "
                                 "must be a class, not a string. For "
                                 "strings, use frame objects")

            if other_class is not self.__class__:
                # The default is to convert via cartesian coordinates
                new_rep = other_class.from_cartesian(self.to_cartesian())
            else:
                new_rep = self

            new_rep._differentials = self._re_represent_differentials(
                new_rep, differential_class)

            return new_rep

    def with_differentials(self, differentials):
        """
        Create a new representation with the same positions as this
        representation, but with these new differentials.

        Differential keys that already exist in this object's differential dict
        are overwritten.

        Parameters
        ----------
        differentials : Sequence of `~astropy.coordinates.BaseDifferential`
            The differentials for the new representation to have.

        Returns
        -------
        newrepr
            A copy of this representation, but with the ``differentials`` as
            its differentials.
        """
        if not differentials:
            return self

        args = [getattr(self, component) for component in self.components]

        # We shallow copy the differentials dictionary so we don't update the
        # current object's dictionary when adding new keys
        new_rep = self.__class__(*args, differentials=self.differentials.copy(),
                                 copy=False)
        new_rep._differentials.update(
            new_rep._validate_differentials(differentials))

        return new_rep

    def without_differentials(self):
        """Return a copy of the representation without attached differentials.

        Returns
        -------
        newrepr
            A shallow copy of this representation, without any differentials.
            If no differentials were present, no copy is made.
        """

        if not self._differentials:
            return self

        args = [getattr(self, component) for component in self.components]
        return self.__class__(*args, copy=False)

    @classmethod
    def from_representation(cls, representation):
        """Create a new instance of this representation from another one.

        Parameters
        ----------
        representation : `~astropy.coordinates.BaseRepresentation` instance
            The presentation that should be converted to this class.
        """
        return representation.represent_as(cls)

    def _apply(self, method, *args, **kwargs):
        """Create a new representation with ``method`` applied to the component
        data.

        This is not a simple inherit from ``BaseRepresentationOrDifferential``
        because we need to call ``._apply()`` on any associated differential
        classes.

        See docstring for `BaseRepresentationOrDifferential._apply`.

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
        rep = super()._apply(method, *args, **kwargs)

        rep._differentials = dict(
            [(k, diff._apply(method, *args, **kwargs))
             for k, diff in self._differentials.items()])
        return rep

    def _scale_operation(self, op, *args):
        """Scale all non-angular components, leaving angular ones unchanged.

        Parameters
        ----------
        op : `~operator` callable
            Operator to apply (e.g., `~operator.mul`, `~operator.neg`, etc.
        *args
            Any arguments required for the operator (typically, what is to
            be multiplied with, divided by).
        """

        self._raise_if_has_differentials(op.__name__)

        results = []
        for component, cls in self.attr_classes.items():
            value = getattr(self, component)
            if issubclass(cls, Angle):
                results.append(value)
            else:
                results.append(op(value, *args))

        # try/except catches anything that cannot initialize the class, such
        # as operations that returned NotImplemented or a representation
        # instead of a quantity (as would happen for, e.g., rep * rep).
        try:
            return self.__class__(*results)
        except Exception:
            return NotImplemented

    def _combine_operation(self, op, other, reverse=False):
        """Combine two representation.

        By default, operate on the cartesian representations of both.

        Parameters
        ----------
        op : `~operator` callable
            Operator to apply (e.g., `~operator.add`, `~operator.sub`, etc.
        other : `~astropy.coordinates.BaseRepresentation` instance
            The other representation.
        reverse : bool
            Whether the operands should be reversed (e.g., as we got here via
            ``self.__rsub__`` because ``self`` is a subclass of ``other``).
        """
        self._raise_if_has_differentials(op.__name__)

        result = self.to_cartesian()._combine_operation(op, other, reverse)
        if result is NotImplemented:
            return NotImplemented
        else:
            return self.from_cartesian(result)

    # We need to override this setter to support differentials
    @BaseRepresentationOrDifferential.shape.setter
    def shape(self, shape):
        orig_shape = self.shape

        # See: https://stackoverflow.com/questions/3336767/ for an example
        BaseRepresentationOrDifferential.shape.fset(self, shape)

        # also try to perform shape-setting on any associated differentials
        try:
            for k in self.differentials:
                self.differentials[k].shape = shape
        except Exception:
            BaseRepresentationOrDifferential.shape.fset(self, orig_shape)
            for k in self.differentials:
                self.differentials[k].shape = orig_shape

            raise

    def norm(self):
        """Vector norm.

        The norm is the standard Frobenius norm, i.e., the square root of the
        sum of the squares of all components with non-angular units.

        Note that any associated differentials will be dropped during this
        operation.

        Returns
        -------
        norm : `astropy.units.Quantity`
            Vector norm, with the same shape as the representation.
        """
        return np.sqrt(functools.reduce(
            operator.add, (getattr(self, component)**2
                           for component, cls in self.attr_classes.items()
                           if not issubclass(cls, Angle))))

    def mean(self, *args, **kwargs):
        """Vector mean.

        Averaging is done by converting the representation to cartesian, and
        taking the mean of the x, y, and z components. The result is converted
        back to the same representation as the input.

        Refer to `~numpy.mean` for full documentation of the arguments, noting
        that ``axis`` is the entry in the ``shape`` of the representation, and
        that the ``out`` argument cannot be used.

        Returns
        -------
        mean : representation
            Vector mean, in the same representation as that of the input.
        """
        self._raise_if_has_differentials('mean')
        return self.from_cartesian(self.to_cartesian().mean(*args, **kwargs))

    def sum(self, *args, **kwargs):
        """Vector sum.

        Adding is done by converting the representation to cartesian, and
        summing the x, y, and z components. The result is converted back to the
        same representation as the input.

        Refer to `~numpy.sum` for full documentation of the arguments, noting
        that ``axis`` is the entry in the ``shape`` of the representation, and
        that the ``out`` argument cannot be used.

        Returns
        -------
        sum : representation
            Vector sum, in the same representation as that of the input.
        """
        self._raise_if_has_differentials('sum')
        return self.from_cartesian(self.to_cartesian().sum(*args, **kwargs))

    def dot(self, other):
        """Dot product of two representations.

        The calculation is done by converting both ``self`` and ``other``
        to `~astropy.coordinates.CartesianRepresentation`.

        Note that any associated differentials will be dropped during this
        operation.

        Parameters
        ----------
        other : `~astropy.coordinates.BaseRepresentation`
            The representation to take the dot product with.

        Returns
        -------
        dot_product : `~astropy.units.Quantity`
            The sum of the product of the x, y, and z components of the
            cartesian representations of ``self`` and ``other``.
        """
        return self.to_cartesian().dot(other)

    def cross(self, other):
        """Vector cross product of two representations.

        The calculation is done by converting both ``self`` and ``other``
        to `~astropy.coordinates.CartesianRepresentation`, and converting the
        result back to the type of representation of ``self``.

        Parameters
        ----------
        other : representation
            The representation to take the cross product with.

        Returns
        -------
        cross_product : representation
            With vectors perpendicular to both ``self`` and ``other``, in the
            same type of representation as ``self``.
        """
        self._raise_if_has_differentials('cross')
        return self.from_cartesian(self.to_cartesian().cross(other))


class CartesianRepresentation(BaseRepresentation):
    """
    Representation of points in 3D cartesian coordinates.

    Parameters
    ----------
    x, y, z : `~astropy.units.Quantity` or array
        The x, y, and z coordinates of the point(s). If ``x``, ``y``, and ``z``
        have different shapes, they should be broadcastable. If not quantity,
        ``unit`` should be set.  If only ``x`` is given, it is assumed that it
        contains an array with the 3 coordinates stored along ``xyz_axis``.
    unit : `~astropy.units.Unit` or str
        If given, the coordinates will be converted to this unit (or taken to
        be in this unit if not given.
    xyz_axis : int, optional
        The axis along which the coordinates are stored when a single array is
        provided rather than distinct ``x``, ``y``, and ``z`` (default: 0).

    differentials : dict, `CartesianDifferential`, optional
        Any differential classes that should be associated with this
        representation. The input must either be a single
        `CartesianDifferential` instance, or a dictionary of
        `CartesianDifferential` s with keys set to a string representation of
        the SI unit with which the differential (derivative) is taken. For
        example, for a velocity differential on a positional representation, the
        key would be ``'s'`` for seconds, indicating that the derivative is a
        time derivative.

    copy : bool, optional
        If `True` (default), arrays will be copied rather than referenced.
    """

    attr_classes = OrderedDict([('x', u.Quantity),
                                ('y', u.Quantity),
                                ('z', u.Quantity)])

    _xyz = None

    def __init__(self, x, y=None, z=None, unit=None, xyz_axis=None,
                 differentials=None, copy=True):

        if y is None and z is None:
            if isinstance(x, np.ndarray) and x.dtype.kind not in 'OV':
                # Short-cut for 3-D array input.
                x = u.Quantity(x, unit, copy=copy, subok=True)
                # Keep a link to the array with all three coordinates
                # so that we can return it quickly if needed in get_xyz.
                self._xyz = x
                if xyz_axis:
                    x = np.moveaxis(x, xyz_axis, 0)
                    self._xyz_axis = xyz_axis
                else:
                    self._xyz_axis = 0

                self._x, self._y, self._z = x
                self._differentials = self._validate_differentials(differentials)
                return

            else:
                x, y, z = x

        if xyz_axis is not None:
            raise ValueError("xyz_axis should only be set if x, y, and z are "
                             "in a single array passed in through x, "
                             "i.e., y and z should not be not given.")

        if y is None or z is None:
            raise ValueError("x, y, and z are required to instantiate {0}"
                             .format(self.__class__.__name__))

        if unit is not None:
            x = u.Quantity(x, unit, copy=copy, subok=True)
            y = u.Quantity(y, unit, copy=copy, subok=True)
            z = u.Quantity(z, unit, copy=copy, subok=True)
            copy = False

        super().__init__(x, y, z, copy=copy, differentials=differentials)
        if not (self._x.unit.is_equivalent(self._y.unit) and
                self._x.unit.is_equivalent(self._z.unit)):
            raise u.UnitsError("x, y, and z should have matching physical types")

    def unit_vectors(self):
        l = np.broadcast_to(1.*u.one, self.shape, subok=True)
        o = np.broadcast_to(0.*u.one, self.shape, subok=True)
        return OrderedDict(
            (('x', CartesianRepresentation(l, o, o, copy=False)),
             ('y', CartesianRepresentation(o, l, o, copy=False)),
             ('z', CartesianRepresentation(o, o, l, copy=False))))

    def scale_factors(self):
        l = np.broadcast_to(1.*u.one, self.shape, subok=True)
        return OrderedDict((('x', l), ('y', l), ('z', l)))

    def get_xyz(self, xyz_axis=0):
        """Return a vector array of the x, y, and z coordinates.

        Parameters
        ----------
        xyz_axis : int, optional
            The axis in the final array along which the x, y, z components
            should be stored (default: 0).

        Returns
        -------
        xyz : `~astropy.units.Quantity`
            With dimension 3 along ``xyz_axis``.  Note that, if possible,
            this will be a view.
        """
        if self._xyz is not None:
            if self._xyz_axis == xyz_axis:
                return self._xyz
            else:
                return np.moveaxis(self._xyz, self._xyz_axis, xyz_axis)

        # Create combined array.  TO DO: keep it in _xyz for repeated use?
        # But then in-place changes have to cancel it. Likely best to
        # also update components.
        return _combine_xyz(self._x, self._y, self._z, xyz_axis=xyz_axis)

    xyz = property(get_xyz)

    @classmethod
    def from_cartesian(cls, other):
        return other

    def to_cartesian(self):
        return self

    def transform(self, matrix):
        """
        Transform the cartesian coordinates using a 3x3 matrix.

        This returns a new representation and does not modify the original one.
        Any differentials attached to this representation will also be
        transformed.

        Parameters
        ----------
        matrix : `~numpy.ndarray`
            A 3x3 transformation matrix, such as a rotation matrix.


        Examples
        --------

        We can start off by creating a cartesian representation object:

            >>> from astropy import units as u
            >>> from astropy.coordinates import CartesianRepresentation
            >>> rep = CartesianRepresentation([1, 2] * u.pc,
            ...                               [2, 3] * u.pc,
            ...                               [3, 4] * u.pc)

        We now create a rotation matrix around the z axis:

            >>> from astropy.coordinates.matrix_utilities import rotation_matrix
            >>> rotation = rotation_matrix(30 * u.deg, axis='z')

        Finally, we can apply this transformation:

            >>> rep_new = rep.transform(rotation)
            >>> rep_new.xyz  # doctest: +FLOAT_CMP
            <Quantity [[ 1.8660254 , 3.23205081],
                       [ 1.23205081, 1.59807621],
                       [ 3.        , 4.        ]] pc>
        """
        # erfa rxp: Multiply a p-vector by an r-matrix.
        p = erfa_ufunc.rxp(matrix, self.get_xyz(xyz_axis=-1))
        # Handle differentials attached to this representation
        if self.differentials:
            # TODO: speed this up going via d.d_xyz.
            new_diffs = dict(
                (k, d.from_cartesian(d.to_cartesian().transform(matrix)))
                for k, d in self.differentials.items())
        else:
            new_diffs = None

        return self.__class__(p, xyz_axis=-1, copy=False, differentials=new_diffs)

    def _combine_operation(self, op, other, reverse=False):
        self._raise_if_has_differentials(op.__name__)

        try:
            other_c = other.to_cartesian()
        except Exception:
            return NotImplemented

        first, second = ((self, other_c) if not reverse else
                         (other_c, self))
        return self.__class__(*(op(getattr(first, component),
                                   getattr(second, component))
                                for component in first.components))

    def norm(self):
        """Vector norm.

        The norm is the standard Frobenius norm, i.e., the square root of the
        sum of the squares of all components with non-angular units.

        Note that any associated differentials will be dropped during this
        operation.

        Returns
        -------
        norm : `astropy.units.Quantity`
            Vector norm, with the same shape as the representation.
        """
        # erfa pm: Modulus of p-vector.
        return erfa_ufunc.pm(self.get_xyz(xyz_axis=-1))

    def mean(self, *args, **kwargs):
        """Vector mean.

        Returns a new CartesianRepresentation instance with the means of the
        x, y, and z components.

        Refer to `~numpy.mean` for full documentation of the arguments, noting
        that ``axis`` is the entry in the ``shape`` of the representation, and
        that the ``out`` argument cannot be used.
        """
        self._raise_if_has_differentials('mean')
        return self._apply('mean', *args, **kwargs)

    def sum(self, *args, **kwargs):
        """Vector sum.

        Returns a new CartesianRepresentation instance with the sums of the
        x, y, and z components.

        Refer to `~numpy.sum` for full documentation of the arguments, noting
        that ``axis`` is the entry in the ``shape`` of the representation, and
        that the ``out`` argument cannot be used.
        """
        self._raise_if_has_differentials('sum')
        return self._apply('sum', *args, **kwargs)

    def dot(self, other):
        """Dot product of two representations.

        Note that any associated differentials will be dropped during this
        operation.

        Parameters
        ----------
        other : representation
            If not already cartesian, it is converted.

        Returns
        -------
        dot_product : `~astropy.units.Quantity`
            The sum of the product of the x, y, and z components of ``self``
            and ``other``.
        """
        try:
            other_c = other.to_cartesian()
        except Exception:
            raise TypeError("cannot only take dot product with another "
                            "representation, not a {0} instance."
                            .format(type(other)))
        # erfa pdp: p-vector inner (=scalar=dot) product.
        return erfa_ufunc.pdp(self.get_xyz(xyz_axis=-1),
                              other_c.get_xyz(xyz_axis=-1))

    def cross(self, other):
        """Cross product of two representations.

        Parameters
        ----------
        other : representation
            If not already cartesian, it is converted.

        Returns
        -------
        cross_product : `~astropy.coordinates.CartesianRepresentation`
            With vectors perpendicular to both ``self`` and ``other``.
        """
        self._raise_if_has_differentials('cross')
        try:
            other_c = other.to_cartesian()
        except Exception:
            raise TypeError("cannot only take cross product with another "
                            "representation, not a {0} instance."
                            .format(type(other)))
        # erfa pxp: p-vector outer (=vector=cross) product.
        sxo = erfa_ufunc.pxp(self.get_xyz(xyz_axis=-1),
                             other_c.get_xyz(xyz_axis=-1))
        return self.__class__(sxo, xyz_axis=-1)


class UnitSphericalRepresentation(BaseRepresentation):
    """
    Representation of points on a unit sphere.

    Parameters
    ----------
    lon, lat : `~astropy.units.Quantity` or str
        The longitude and latitude of the point(s), in angular units. The
        latitude should be between -90 and 90 degrees, and the longitude will
        be wrapped to an angle between 0 and 360 degrees. These can also be
        instances of `~astropy.coordinates.Angle`,
        `~astropy.coordinates.Longitude`, or `~astropy.coordinates.Latitude`.

    differentials : dict, `BaseDifferential`, optional
        Any differential classes that should be associated with this
        representation. The input must either be a single `BaseDifferential`
        instance (see `._compatible_differentials` for valid types), or a
        dictionary of of differential instances with keys set to a string
        representation of the SI unit with which the differential (derivative)
        is taken. For example, for a velocity differential on a positional
        representation, the key would be ``'s'`` for seconds, indicating that
        the derivative is a time derivative.

    copy : bool, optional
        If `True` (default), arrays will be copied rather than referenced.
    """

    attr_classes = OrderedDict([('lon', Longitude),
                                ('lat', Latitude)])

    @classproperty
    def _dimensional_representation(cls):
        return SphericalRepresentation

    def __init__(self, lon, lat, differentials=None, copy=True):
        super().__init__(lon, lat, differentials=differentials, copy=copy)

    @property
    def _compatible_differentials(self):
        return [UnitSphericalDifferential, UnitSphericalCosLatDifferential,
                SphericalDifferential, SphericalCosLatDifferential,
                RadialDifferential]

    # Could let the metaclass define these automatically, but good to have
    # a bit clearer docstrings.
    @property
    def lon(self):
        """
        The longitude of the point(s).
        """
        return self._lon

    @property
    def lat(self):
        """
        The latitude of the point(s).
        """
        return self._lat

    def unit_vectors(self):
        sinlon, coslon = np.sin(self.lon), np.cos(self.lon)
        sinlat, coslat = np.sin(self.lat), np.cos(self.lat)
        return OrderedDict(
            (('lon', CartesianRepresentation(-sinlon, coslon, 0., copy=False)),
             ('lat', CartesianRepresentation(-sinlat*coslon, -sinlat*sinlon,
                                             coslat, copy=False))))

    def scale_factors(self, omit_coslat=False):
        sf_lat = np.broadcast_to(1./u.radian, self.shape, subok=True)
        sf_lon = sf_lat if omit_coslat else np.cos(self.lat) / u.radian
        return OrderedDict((('lon', sf_lon),
                            ('lat', sf_lat)))

    def to_cartesian(self):
        """
        Converts spherical polar coordinates to 3D rectangular cartesian
        coordinates.
        """
        # NUMPY_LT_1_16 cannot create a vector automatically
        p = u.Quantity(np.empty(self.shape + (3,)), u.dimensionless_unscaled,
                       copy=False)
        # erfa s2c: Convert [unit]spherical coordinates to Cartesian.
        p = erfa_ufunc.s2c(self.lon, self.lat, p)
        return CartesianRepresentation(p, xyz_axis=-1, copy=False)

    @classmethod
    def from_cartesian(cls, cart):
        """
        Converts 3D rectangular cartesian coordinates to spherical polar
        coordinates.
        """
        p = cart.get_xyz(xyz_axis=-1)
        # erfa c2s: P-vector to [unit]spherical coordinates.
        return cls(*erfa_ufunc.c2s(p), copy=False)

    def represent_as(self, other_class, differential_class=None):
        # Take a short cut if the other class is a spherical representation

        # TODO: this could be optimized to shortcut even if a differential_class
        # is passed in, using the ._re_represent_differentials() method
        if inspect.isclass(other_class) and not differential_class:
            if issubclass(other_class, PhysicsSphericalRepresentation):
                return other_class(phi=self.lon, theta=90 * u.deg - self.lat, r=1.0,
                                   copy=False)
            elif issubclass(other_class, SphericalRepresentation):
                return other_class(lon=self.lon, lat=self.lat, distance=1.0,
                                   copy=False)

        return super().represent_as(other_class, differential_class)

    def __mul__(self, other):
        self._raise_if_has_differentials('multiplication')
        return self._dimensional_representation(lon=self.lon, lat=self.lat,
                                                distance=1. * other)

    def __truediv__(self, other):
        self._raise_if_has_differentials('division')
        return self._dimensional_representation(lon=self.lon, lat=self.lat,
                                                distance=1. / other)

    def __neg__(self):
        self._raise_if_has_differentials('negation')
        return self.__class__(self.lon + 180. * u.deg, -self.lat, copy=False)

    def norm(self):
        """Vector norm.

        The norm is the standard Frobenius norm, i.e., the square root of the
        sum of the squares of all components with non-angular units, which is
        always unity for vectors on the unit sphere.

        Returns
        -------
        norm : `~astropy.units.Quantity`
            Dimensionless ones, with the same shape as the representation.
        """
        return u.Quantity(np.ones(self.shape), u.dimensionless_unscaled,
                          copy=False)

    def _combine_operation(self, op, other, reverse=False):
        self._raise_if_has_differentials(op.__name__)

        result = self.to_cartesian()._combine_operation(op, other, reverse)
        if result is NotImplemented:
            return NotImplemented
        else:
            return self._dimensional_representation.from_cartesian(result)

    def mean(self, *args, **kwargs):
        """Vector mean.

        The representation is converted to cartesian, the means of the x, y,
        and z components are calculated, and the result is converted to a
        `~astropy.coordinates.SphericalRepresentation`.

        Refer to `~numpy.mean` for full documentation of the arguments, noting
        that ``axis`` is the entry in the ``shape`` of the representation, and
        that the ``out`` argument cannot be used.
        """
        self._raise_if_has_differentials('mean')
        return self._dimensional_representation.from_cartesian(
            self.to_cartesian().mean(*args, **kwargs))

    def sum(self, *args, **kwargs):
        """Vector sum.

        The representation is converted to cartesian, the sums of the x, y,
        and z components are calculated, and the result is converted to a
        `~astropy.coordinates.SphericalRepresentation`.

        Refer to `~numpy.sum` for full documentation of the arguments, noting
        that ``axis`` is the entry in the ``shape`` of the representation, and
        that the ``out`` argument cannot be used.
        """
        self._raise_if_has_differentials('sum')
        return self._dimensional_representation.from_cartesian(
            self.to_cartesian().sum(*args, **kwargs))

    def cross(self, other):
        """Cross product of two representations.

        The calculation is done by converting both ``self`` and ``other``
        to `~astropy.coordinates.CartesianRepresentation`, and converting the
        result back to `~astropy.coordinates.SphericalRepresentation`.

        Parameters
        ----------
        other : representation
            The representation to take the cross product with.

        Returns
        -------
        cross_product : `~astropy.coordinates.SphericalRepresentation`
            With vectors perpendicular to both ``self`` and ``other``.
        """
        self._raise_if_has_differentials('cross')
        return self._dimensional_representation.from_cartesian(
            self.to_cartesian().cross(other))


class RadialRepresentation(BaseRepresentation):
    """
    Representation of the distance of points from the origin.

    Note that this is mostly intended as an internal helper representation.
    It can do little else but being used as a scale in multiplication.

    Parameters
    ----------
    distance : `~astropy.units.Quantity`
        The distance of the point(s) from the origin.

    differentials : dict, `BaseDifferential`, optional
        Any differential classes that should be associated with this
        representation. The input must either be a single `BaseDifferential`
        instance (see `._compatible_differentials` for valid types), or a
        dictionary of of differential instances with keys set to a string
        representation of the SI unit with which the differential (derivative)
        is taken. For example, for a velocity differential on a positional
        representation, the key would be ``'s'`` for seconds, indicating that
        the derivative is a time derivative.

    copy : bool, optional
        If `True` (default), arrays will be copied rather than referenced.
    """

    attr_classes = OrderedDict([('distance', u.Quantity)])

    def __init__(self, distance, differentials=None, copy=True):
        super().__init__(distance, copy=copy, differentials=differentials)

    @property
    def distance(self):
        """
        The distance from the origin to the point(s).
        """
        return self._distance

    def unit_vectors(self):
        """Cartesian unit vectors are undefined for radial representation."""
        raise NotImplementedError('Cartesian unit vectors are undefined for '
                                  '{0} instances'.format(self.__class__))

    def scale_factors(self):
        l = np.broadcast_to(1.*u.one, self.shape, subok=True)
        return OrderedDict((('distance', l),))

    def to_cartesian(self):
        """Cannot convert radial representation to cartesian."""
        raise NotImplementedError('cannot convert {0} instance to cartesian.'
                                  .format(self.__class__))

    @classmethod
    def from_cartesian(cls, cart):
        """
        Converts 3D rectangular cartesian coordinates to radial coordinate.
        """
        return cls(distance=cart.norm(), copy=False)

    def _scale_operation(self, op, *args):
        self._raise_if_has_differentials(op.__name__)
        return op(self.distance, *args)

    def norm(self):
        """Vector norm.

        Just the distance itself.

        Returns
        -------
        norm : `~astropy.units.Quantity`
            Dimensionless ones, with the same shape as the representation.
        """
        return self.distance

    def _combine_operation(self, op, other, reverse=False):
        return NotImplemented


class SphericalRepresentation(BaseRepresentation):
    """
    Representation of points in 3D spherical coordinates.

    Parameters
    ----------
    lon, lat : `~astropy.units.Quantity`
        The longitude and latitude of the point(s), in angular units. The
        latitude should be between -90 and 90 degrees, and the longitude will
        be wrapped to an angle between 0 and 360 degrees. These can also be
        instances of `~astropy.coordinates.Angle`,
        `~astropy.coordinates.Longitude`, or `~astropy.coordinates.Latitude`.

    distance : `~astropy.units.Quantity`
        The distance to the point(s). If the distance is a length, it is
        passed to the :class:`~astropy.coordinates.Distance` class, otherwise
        it is passed to the :class:`~astropy.units.Quantity` class.

    differentials : dict, `BaseDifferential`, optional
        Any differential classes that should be associated with this
        representation. The input must either be a single `BaseDifferential`
        instance (see `._compatible_differentials` for valid types), or a
        dictionary of of differential instances with keys set to a string
        representation of the SI unit with which the differential (derivative)
        is taken. For example, for a velocity differential on a positional
        representation, the key would be ``'s'`` for seconds, indicating that
        the derivative is a time derivative.

    copy : bool, optional
        If `True` (default), arrays will be copied rather than referenced.
    """

    attr_classes = OrderedDict([('lon', Longitude),
                                ('lat', Latitude),
                                ('distance', u.Quantity)])
    _unit_representation = UnitSphericalRepresentation

    def __init__(self, lon, lat, distance, differentials=None, copy=True):
        super().__init__(lon, lat, distance, copy=copy,
                         differentials=differentials)
        if self._distance.unit.physical_type == 'length':
            try:
                self._distance = Distance(self._distance, copy=False)
            except ValueError as e:
                if e.args[0].startswith('Distance must be >= 0'):
                    raise ValueError("Distance must be >= 0. To allow negative "
                                     "distance values, you must explicitly pass"
                                     " in a `Distance` object with the the "
                                     "argument 'allow_negative=True'.")
                else:
                    raise

    @property
    def _compatible_differentials(self):
        return [UnitSphericalDifferential, UnitSphericalCosLatDifferential,
                SphericalDifferential, SphericalCosLatDifferential,
                RadialDifferential]

    @property
    def lon(self):
        """
        The longitude of the point(s).
        """
        return self._lon

    @property
    def lat(self):
        """
        The latitude of the point(s).
        """
        return self._lat

    @property
    def distance(self):
        """
        The distance from the origin to the point(s).
        """
        return self._distance

    def unit_vectors(self):
        sinlon, coslon = np.sin(self.lon), np.cos(self.lon)
        sinlat, coslat = np.sin(self.lat), np.cos(self.lat)
        return OrderedDict(
            (('lon', CartesianRepresentation(-sinlon, coslon, 0., copy=False)),
             ('lat', CartesianRepresentation(-sinlat*coslon, -sinlat*sinlon,
                                             coslat, copy=False)),
             ('distance', CartesianRepresentation(coslat*coslon, coslat*sinlon,
                                                  sinlat, copy=False))))

    def scale_factors(self, omit_coslat=False):
        sf_lat = self.distance / u.radian
        sf_lon = sf_lat if omit_coslat else sf_lat * np.cos(self.lat)
        sf_distance = np.broadcast_to(1.*u.one, self.shape, subok=True)
        return OrderedDict((('lon', sf_lon),
                            ('lat', sf_lat),
                            ('distance', sf_distance)))

    def represent_as(self, other_class, differential_class=None):
        # Take a short cut if the other class is a spherical representation

        # TODO: this could be optimized to shortcut even if a differential_class
        # is passed in, using the ._re_represent_differentials() method
        if inspect.isclass(other_class) and not differential_class:
            if issubclass(other_class, PhysicsSphericalRepresentation):
                return other_class(phi=self.lon, theta=90 * u.deg - self.lat,
                                   r=self.distance, copy=False)
            elif issubclass(other_class, UnitSphericalRepresentation):
                return other_class(lon=self.lon, lat=self.lat, copy=False)

        return super().represent_as(other_class, differential_class)

    def to_cartesian(self):
        """
        Converts spherical polar coordinates to 3D rectangular cartesian
        coordinates.
        """

        # We need to convert Distance to Quantity to allow negative values.
        if isinstance(self.distance, Distance):
            d = self.distance.view(u.Quantity)
        else:
            d = self.distance

        # NUMPY_LT_1_16 cannot create a vector automatically
        p = u.Quantity(np.empty(self.shape + (3,)), d.unit, copy=False)
        # erfa s2p: Convert spherical polar coordinates to p-vector.
        p = erfa_ufunc.s2p(self.lon, self.lat, d, p)

        return CartesianRepresentation(p, xyz_axis=-1, copy=False)

    @classmethod
    def from_cartesian(cls, cart):
        """
        Converts 3D rectangular cartesian coordinates to spherical polar
        coordinates.
        """
        p = cart.get_xyz(xyz_axis=-1)
        # erfa p2s: P-vector to spherical polar coordinates.
        return cls(*erfa_ufunc.p2s(p), copy=False)

    def norm(self):
        """Vector norm.

        The norm is the standard Frobenius norm, i.e., the square root of the
        sum of the squares of all components with non-angular units.  For
        spherical coordinates, this is just the absolute value of the distance.

        Returns
        -------
        norm : `astropy.units.Quantity`
            Vector norm, with the same shape as the representation.
        """
        return np.abs(self.distance)

    def __neg__(self):
        self._raise_if_has_differentials('negation')
        return self.__class__(self.lon + 180. * u.deg, -self.lat, self.distance,
                              copy=False)


class PhysicsSphericalRepresentation(BaseRepresentation):
    """
    Representation of points in 3D spherical coordinates (using the physics
    convention of using ``phi`` and ``theta`` for azimuth and inclination
    from the pole).

    Parameters
    ----------
    phi, theta : `~astropy.units.Quantity` or str
        The azimuth and inclination of the point(s), in angular units. The
        inclination should be between 0 and 180 degrees, and the azimuth will
        be wrapped to an angle between 0 and 360 degrees. These can also be
        instances of `~astropy.coordinates.Angle`.  If ``copy`` is False, `phi`
        will be changed inplace if it is not between 0 and 360 degrees.

    r : `~astropy.units.Quantity`
        The distance to the point(s). If the distance is a length, it is
        passed to the :class:`~astropy.coordinates.Distance` class, otherwise
        it is passed to the :class:`~astropy.units.Quantity` class.

    differentials : dict, `PhysicsSphericalDifferential`, optional
        Any differential classes that should be associated with this
        representation. The input must either be a single
        `PhysicsSphericalDifferential` instance, or a dictionary of of
        differential instances with keys set to a string representation of the
        SI unit with which the differential (derivative) is taken. For example,
        for a velocity differential on a positional representation, the key
        would be ``'s'`` for seconds, indicating that the derivative is a time
        derivative.

    copy : bool, optional
        If `True` (default), arrays will be copied rather than referenced.
    """

    attr_classes = OrderedDict([('phi', Angle),
                                ('theta', Angle),
                                ('r', u.Quantity)])

    def __init__(self, phi, theta, r, differentials=None, copy=True):
        super().__init__(phi, theta, r, copy=copy, differentials=differentials)

        # Wrap/validate phi/theta
        if copy:
            self._phi = self._phi.wrap_at(360 * u.deg)
        else:
            # necessary because the above version of `wrap_at` has to be a copy
            self._phi.wrap_at(360 * u.deg, inplace=True)

        if np.any(self._theta < 0.*u.deg) or np.any(self._theta > 180.*u.deg):
            raise ValueError('Inclination angle(s) must be within '
                             '0 deg <= angle <= 180 deg, '
                             'got {0}'.format(theta.to(u.degree)))

        if self._r.unit.physical_type == 'length':
            self._r = self._r.view(Distance)

    @property
    def phi(self):
        """
        The azimuth of the point(s).
        """
        return self._phi

    @property
    def theta(self):
        """
        The elevation of the point(s).
        """
        return self._theta

    @property
    def r(self):
        """
        The distance from the origin to the point(s).
        """
        return self._r

    def unit_vectors(self):
        sinphi, cosphi = np.sin(self.phi), np.cos(self.phi)
        sintheta, costheta = np.sin(self.theta), np.cos(self.theta)
        return OrderedDict(
            (('phi', CartesianRepresentation(-sinphi, cosphi, 0., copy=False)),
             ('theta', CartesianRepresentation(costheta*cosphi,
                                               costheta*sinphi,
                                               -sintheta, copy=False)),
             ('r', CartesianRepresentation(sintheta*cosphi, sintheta*sinphi,
                                           costheta, copy=False))))

    def scale_factors(self):
        r = self.r / u.radian
        sintheta = np.sin(self.theta)
        l = np.broadcast_to(1.*u.one, self.shape, subok=True)
        return OrderedDict((('phi', r * sintheta),
                            ('theta', r),
                            ('r', l)))

    def represent_as(self, other_class, differential_class=None):
        # Take a short cut if the other class is a spherical representation

        # TODO: this could be optimized to shortcut even if a differential_class
        # is passed in, using the ._re_represent_differentials() method
        if inspect.isclass(other_class) and not differential_class:
            if issubclass(other_class, SphericalRepresentation):
                return other_class(lon=self.phi, lat=90 * u.deg - self.theta,
                                   distance=self.r)
            elif issubclass(other_class, UnitSphericalRepresentation):
                return other_class(lon=self.phi, lat=90 * u.deg - self.theta)

        return super().represent_as(other_class, differential_class)

    def to_cartesian(self):
        """
        Converts spherical polar coordinates to 3D rectangular cartesian
        coordinates.
        """

        # We need to convert Distance to Quantity to allow negative values.
        if isinstance(self.r, Distance):
            d = self.r.view(u.Quantity)
        else:
            d = self.r

        x = d * np.sin(self.theta) * np.cos(self.phi)
        y = d * np.sin(self.theta) * np.sin(self.phi)
        z = d * np.cos(self.theta)

        return CartesianRepresentation(x=x, y=y, z=z, copy=False)

    @classmethod
    def from_cartesian(cls, cart):
        """
        Converts 3D rectangular cartesian coordinates to spherical polar
        coordinates.
        """

        s = np.hypot(cart.x, cart.y)
        r = np.hypot(s, cart.z)

        phi = np.arctan2(cart.y, cart.x)
        theta = np.arctan2(s, cart.z)

        return cls(phi=phi, theta=theta, r=r, copy=False)

    def norm(self):
        """Vector norm.

        The norm is the standard Frobenius norm, i.e., the square root of the
        sum of the squares of all components with non-angular units.  For
        spherical coordinates, this is just the absolute value of the radius.

        Returns
        -------
        norm : `astropy.units.Quantity`
            Vector norm, with the same shape as the representation.
        """
        return np.abs(self.r)


class CylindricalRepresentation(BaseRepresentation):
    """
    Representation of points in 3D cylindrical coordinates.

    Parameters
    ----------
    rho : `~astropy.units.Quantity`
        The distance from the z axis to the point(s).

    phi : `~astropy.units.Quantity` or str
        The azimuth of the point(s), in angular units, which will be wrapped
        to an angle between 0 and 360 degrees. This can also be instances of
        `~astropy.coordinates.Angle`,

    z : `~astropy.units.Quantity`
        The z coordinate(s) of the point(s)

    differentials : dict, `CylindricalDifferential`, optional
        Any differential classes that should be associated with this
        representation. The input must either be a single
        `CylindricalDifferential` instance, or a dictionary of of differential
        instances with keys set to a string representation of the SI unit with
        which the differential (derivative) is taken. For example, for a
        velocity differential on a positional representation, the key would be
        ``'s'`` for seconds, indicating that the derivative is a time
        derivative.

    copy : bool, optional
        If `True` (default), arrays will be copied rather than referenced.
    """

    attr_classes = OrderedDict([('rho', u.Quantity),
                                ('phi', Angle),
                                ('z', u.Quantity)])

    def __init__(self, rho, phi, z, differentials=None, copy=True):
        super().__init__(rho, phi, z, copy=copy, differentials=differentials)

        if not self._rho.unit.is_equivalent(self._z.unit):
            raise u.UnitsError("rho and z should have matching physical types")

    @property
    def rho(self):
        """
        The distance of the point(s) from the z-axis.
        """
        return self._rho

    @property
    def phi(self):
        """
        The azimuth of the point(s).
        """
        return self._phi

    @property
    def z(self):
        """
        The height of the point(s).
        """
        return self._z

    def unit_vectors(self):
        sinphi, cosphi = np.sin(self.phi), np.cos(self.phi)
        l = np.broadcast_to(1., self.shape)
        return OrderedDict(
            (('rho', CartesianRepresentation(cosphi, sinphi, 0, copy=False)),
             ('phi', CartesianRepresentation(-sinphi, cosphi, 0, copy=False)),
             ('z', CartesianRepresentation(0, 0, l, unit=u.one, copy=False))))

    def scale_factors(self):
        rho = self.rho / u.radian
        l = np.broadcast_to(1.*u.one, self.shape, subok=True)
        return OrderedDict((('rho', l),
                            ('phi', rho),
                            ('z', l)))

    @classmethod
    def from_cartesian(cls, cart):
        """
        Converts 3D rectangular cartesian coordinates to cylindrical polar
        coordinates.
        """

        rho = np.hypot(cart.x, cart.y)
        phi = np.arctan2(cart.y, cart.x)
        z = cart.z

        return cls(rho=rho, phi=phi, z=z, copy=False)

    def to_cartesian(self):
        """
        Converts cylindrical polar coordinates to 3D rectangular cartesian
        coordinates.
        """
        x = self.rho * np.cos(self.phi)
        y = self.rho * np.sin(self.phi)
        z = self.z

        return CartesianRepresentation(x=x, y=y, z=z, copy=False)


class MetaBaseDifferential(InheritDocstrings, abc.ABCMeta):
    """Set default ``attr_classes`` and component getters on a Differential.

    For these, the components are those of the base representation prefixed
    by 'd_', and the class is `~astropy.units.Quantity`.
    """
    def __init__(cls, name, bases, dct):
        super().__init__(name, bases, dct)

        # Don't do anything for base helper classes.
        if cls.__name__ in ('BaseDifferential', 'BaseSphericalDifferential',
                            'BaseSphericalCosLatDifferential'):
            return

        if 'base_representation' not in dct:
            raise NotImplementedError('Differential representations must have a'
                                      '"base_representation" class attribute.')

        # If not defined explicitly, create attr_classes.
        if not hasattr(cls, 'attr_classes'):
            base_attr_classes = cls.base_representation.attr_classes
            cls.attr_classes = OrderedDict([('d_' + c, u.Quantity)
                                            for c in base_attr_classes])

        if 'recommended_units' in dct:
            warnings.warn(_recommended_units_deprecation,
                          AstropyDeprecationWarning)
            # Ensure we don't override the property that warns about the
            # deprecation, but that the value remains the same.
            dct.setdefault('_recommended_units', dct.pop('recommended_units'))

        repr_name = cls.get_name()
        if repr_name in DIFFERENTIAL_CLASSES:
            raise ValueError("Differential class {0} already defined"
                             .format(repr_name))

        DIFFERENTIAL_CLASSES[repr_name] = cls
        _invalidate_reprdiff_cls_hash()

        # If not defined explicitly, create properties for the components.
        for component in cls.attr_classes:
            if not hasattr(cls, component):
                setattr(cls, component,
                        property(_make_getter(component),
                                 doc=("Component '{0}' of the Differential."
                                      .format(component))))


class BaseDifferential(BaseRepresentationOrDifferential,
                       metaclass=MetaBaseDifferential):
    r"""A base class representing differentials of representations.

    These represent differences or derivatives along each component.
    E.g., for physics spherical coordinates, these would be
    :math:`\delta r, \delta \theta, \delta \phi`.

    Parameters
    ----------
    d_comp1, d_comp2, d_comp3 : `~astropy.units.Quantity` or subclass
        The components of the 3D differentials.  The names are the keys and the
        subclasses the values of the ``attr_classes`` attribute.
    copy : bool, optional
        If `True` (default), arrays will be copied rather than referenced.

    Notes
    -----
    All differential representation classes should subclass this base class,
    and define an ``base_representation`` attribute with the class of the
    regular `~astropy.coordinates.BaseRepresentation` for which differential
    coordinates are provided. This will set up a default ``attr_classes``
    instance with names equal to the base component names prefixed by ``d_``,
    and all classes set to `~astropy.units.Quantity`, plus properties to access
    those, and a default ``__init__`` for initialization.
    """

    recommended_units = deprecated_attribute('recommended_units', since='3.0')
    _recommended_units = {}

    @classmethod
    def _check_base(cls, base):
        if cls not in base._compatible_differentials:
            raise TypeError("Differential class {0} is not compatible with the "
                            "base (representation) class {1}"
                            .format(cls, base.__class__))

    def _get_deriv_key(self, base):
        """Given a base (representation instance), determine the unit of the
        derivative by removing the representation unit from the component units
        of this differential.
        """

        # This check is just a last resort so we don't return a strange unit key
        # from accidentally passing in the wrong base.
        self._check_base(base)

        for name in base.components:
            comp = getattr(base, name)
            d_comp = getattr(self, 'd_{0}'.format(name), None)
            if d_comp is not None:
                d_unit = comp.unit / d_comp.unit

                # This is quite a bit faster than using to_system() or going
                # through Quantity()
                d_unit_si = d_unit.decompose(u.si.bases)
                d_unit_si._scale = 1 # remove the scale from the unit

                return str(d_unit_si)

        else:
            raise RuntimeError("Invalid representation-differential units! This"
                               " likely happened because either the "
                               "representation or the associated differential "
                               "have non-standard units. Check that the input "
                               "positional data have positional units, and the "
                               "input velocity data have velocity units, or "
                               "are both dimensionless.")

    @classmethod
    def _get_base_vectors(cls, base):
        """Get unit vectors and scale factors from base.

        Parameters
        ----------
        base : instance of ``self.base_representation``
            The points for which the unit vectors and scale factors should be
            retrieved.

        Returns
        -------
        unit_vectors : dict of `CartesianRepresentation`
            In the directions of the coordinates of base.
        scale_factors : dict of `~astropy.units.Quantity`
            Scale factors for each of the coordinates

        Raises
        ------
        TypeError : if the base is not of the correct type
        """
        cls._check_base(base)
        return base.unit_vectors(), base.scale_factors()

    def to_cartesian(self, base):
        """Convert the differential to 3D rectangular cartesian coordinates.

        Parameters
        ----------
        base : instance of ``self.base_representation``
             The points for which the differentials are to be converted: each of
             the components is multiplied by its unit vectors and scale factors.

        Returns
        -------
        This object as a `CartesianDifferential`
        """
        base_e, base_sf = self._get_base_vectors(base)
        return functools.reduce(
            operator.add, (getattr(self, d_c) * base_sf[c] * base_e[c]
                           for d_c, c in zip(self.components, base.components)))

    @classmethod
    def from_cartesian(cls, other, base):
        """Convert the differential from 3D rectangular cartesian coordinates to
        the desired class.

        Parameters
        ----------
        other :
            The object to convert into this differential.
        base : instance of ``self.base_representation``
             The points for which the differentials are to be converted: each of
             the components is multiplied by its unit vectors and scale factors.

        Returns
        -------
        A new differential object that is this class' type.
        """
        base_e, base_sf = cls._get_base_vectors(base)
        return cls(*(other.dot(e / base_sf[component])
                     for component, e in base_e.items()), copy=False)

    def represent_as(self, other_class, base):
        """Convert coordinates to another representation.

        If the instance is of the requested class, it is returned unmodified.
        By default, conversion is done via cartesian coordinates.

        Parameters
        ----------
        other_class : `~astropy.coordinates.BaseRepresentation` subclass
            The type of representation to turn the coordinates into.
        base : instance of ``self.base_representation``, optional
            Base relative to which the differentials are defined.  If the other
            class is a differential representation, the base will be converted
            to its ``base_representation``.
        """
        if other_class is self.__class__:
            return self

        # The default is to convert via cartesian coordinates.
        self_cartesian = self.to_cartesian(base)
        if issubclass(other_class, BaseDifferential):
            base = base.represent_as(other_class.base_representation)
            return other_class.from_cartesian(self_cartesian, base)
        else:
            return other_class.from_cartesian(self_cartesian)

    @classmethod
    def from_representation(cls, representation, base):
        """Create a new instance of this representation from another one.

        Parameters
        ----------
        representation : `~astropy.coordinates.BaseRepresentation` instance
            The presentation that should be converted to this class.
        base : instance of ``cls.base_representation``
            The base relative to which the differentials will be defined. If
            the representation is a differential itself, the base will be
            converted to its ``base_representation`` to help convert it.
        """
        if isinstance(representation, BaseDifferential):
            cartesian = representation.to_cartesian(
                base.represent_as(representation.base_representation))
        else:
            cartesian = representation.to_cartesian()

        return cls.from_cartesian(cartesian, base)

    def _scale_operation(self, op, *args):
        """Scale all components.

        Parameters
        ----------
        op : `~operator` callable
            Operator to apply (e.g., `~operator.mul`, `~operator.neg`, etc.
        *args
            Any arguments required for the operator (typically, what is to
            be multiplied with, divided by).
        """
        scaled_attrs = [op(getattr(self, c), *args) for c in self.components]
        return self.__class__(*scaled_attrs, copy=False)

    def _combine_operation(self, op, other, reverse=False):
        """Combine two differentials, or a differential with a representation.

        If ``other`` is of the same differential type as ``self``, the
        components will simply be combined.  If ``other`` is a representation,
        it will be used as a base for which to evaluate the differential,
        and the result is a new representation.

        Parameters
        ----------
        op : `~operator` callable
            Operator to apply (e.g., `~operator.add`, `~operator.sub`, etc.
        other : `~astropy.coordinates.BaseRepresentation` instance
            The other differential or representation.
        reverse : bool
            Whether the operands should be reversed (e.g., as we got here via
            ``self.__rsub__`` because ``self`` is a subclass of ``other``).
        """
        if isinstance(self, type(other)):
            first, second = (self, other) if not reverse else (other, self)
            return self.__class__(*[op(getattr(first, c), getattr(second, c))
                                    for c in self.components])
        else:
            try:
                self_cartesian = self.to_cartesian(other)
            except TypeError:
                return NotImplemented

            return other._combine_operation(op, self_cartesian, not reverse)

    def __sub__(self, other):
        # avoid "differential - representation".
        if isinstance(other, BaseRepresentation):
            return NotImplemented
        return super().__sub__(other)

    def norm(self, base=None):
        """Vector norm.

        The norm is the standard Frobenius norm, i.e., the square root of the
        sum of the squares of all components with non-angular units.

        Parameters
        ----------
        base : instance of ``self.base_representation``
            Base relative to which the differentials are defined. This is
            required to calculate the physical size of the differential for
            all but cartesian differentials.

        Returns
        -------
        norm : `astropy.units.Quantity`
            Vector norm, with the same shape as the representation.
        """
        return self.to_cartesian(base).norm()


class CartesianDifferential(BaseDifferential):
    """Differentials in of points in 3D cartesian coordinates.

    Parameters
    ----------
    d_x, d_y, d_z : `~astropy.units.Quantity` or array
        The x, y, and z coordinates of the differentials. If ``d_x``, ``d_y``,
        and ``d_z`` have different shapes, they should be broadcastable. If not
        quantities, ``unit`` should be set.  If only ``d_x`` is given, it is
        assumed that it contains an array with the 3 coordinates stored along
        ``xyz_axis``.
    unit : `~astropy.units.Unit` or str
        If given, the differentials will be converted to this unit (or taken to
        be in this unit if not given.
    xyz_axis : int, optional
        The axis along which the coordinates are stored when a single array is
        provided instead of distinct ``d_x``, ``d_y``, and ``d_z`` (default: 0).
    copy : bool, optional
        If `True` (default), arrays will be copied rather than referenced.
    """
    base_representation = CartesianRepresentation
    _d_xyz = None

    def __init__(self, d_x, d_y=None, d_z=None, unit=None, xyz_axis=None,
                 copy=True):

        if d_y is None and d_z is None:
            if isinstance(d_x, np.ndarray) and d_x.dtype.kind not in 'OV':
                # Short-cut for 3-D array input.
                d_x = u.Quantity(d_x, unit, copy=copy, subok=True)
                # Keep a link to the array with all three coordinates
                # so that we can return it quickly if needed in get_xyz.
                self._d_xyz = d_x
                if xyz_axis:
                    d_x = np.moveaxis(d_x, xyz_axis, 0)
                    self._xyz_axis = xyz_axis
                else:
                    self._xyz_axis = 0

                self._d_x, self._d_y, self._d_z = d_x
                return

            else:
                d_x, d_y, d_z = d_x

        if xyz_axis is not None:
            raise ValueError("xyz_axis should only be set if d_x, d_y, and d_z "
                             "are in a single array passed in through d_x, "
                             "i.e., d_y and d_z should not be not given.")

        if d_y is None or d_z is None:
            raise ValueError("d_x, d_y, and d_z are required to instantiate {0}"
                             .format(self.__class__.__name__))

        if unit is not None:
            d_x = u.Quantity(d_x, unit, copy=copy, subok=True)
            d_y = u.Quantity(d_y, unit, copy=copy, subok=True)
            d_z = u.Quantity(d_z, unit, copy=copy, subok=True)
            copy = False

        super().__init__(d_x, d_y, d_z, copy=copy)
        if not (self._d_x.unit.is_equivalent(self._d_y.unit) and
                self._d_x.unit.is_equivalent(self._d_z.unit)):
            raise u.UnitsError('d_x, d_y and d_z should have equivalent units.')

    def to_cartesian(self, base=None):
        return CartesianRepresentation(*[getattr(self, c) for c
                                         in self.components])

    @classmethod
    def from_cartesian(cls, other, base=None):
        return cls(*[getattr(other, c) for c in other.components])

    def get_d_xyz(self, xyz_axis=0):
        """Return a vector array of the x, y, and z coordinates.

        Parameters
        ----------
        xyz_axis : int, optional
            The axis in the final array along which the x, y, z components
            should be stored (default: 0).

        Returns
        -------
        d_xyz : `~astropy.units.Quantity`
            With dimension 3 along ``xyz_axis``.  Note that, if possible,
            this will be a view.
        """
        if self._d_xyz is not None:
            if self._xyz_axis == xyz_axis:
                return self._d_xyz
            else:
                return np.moveaxis(self._d_xyz, self._xyz_axis, xyz_axis)

        # Create combined array.  TO DO: keep it in _d_xyz for repeated use?
        # But then in-place changes have to cancel it. Likely best to
        # also update components.
        return _combine_xyz(self._d_x, self._d_y, self._d_z, xyz_axis=xyz_axis)

    d_xyz = property(get_d_xyz)


class BaseSphericalDifferential(BaseDifferential):
    def _d_lon_coslat(self, base):
        """Convert longitude differential d_lon to d_lon_coslat.

        Parameters
        ----------
        base : instance of ``cls.base_representation``
            The base from which the latitude will be taken.
        """
        self._check_base(base)
        return self.d_lon * np.cos(base.lat)

    @classmethod
    def _get_d_lon(cls, d_lon_coslat, base):
        """Convert longitude differential d_lon_coslat to d_lon.

        Parameters
        ----------
        d_lon_coslat : `~astropy.units.Quantity`
            Longitude differential that includes ``cos(lat)``.
        base : instance of ``cls.base_representation``
            The base from which the latitude will be taken.
        """
        cls._check_base(base)
        return d_lon_coslat / np.cos(base.lat)

    def _combine_operation(self, op, other, reverse=False):
        """Combine two differentials, or a differential with a representation.

        If ``other`` is of the same differential type as ``self``, the
        components will simply be combined.  If both are different parts of
        a `~astropy.coordinates.SphericalDifferential` (e.g., a
        `~astropy.coordinates.UnitSphericalDifferential` and a
        `~astropy.coordinates.RadialDifferential`), they will combined
        appropriately.

        If ``other`` is a representation, it will be used as a base for which
        to evaluate the differential, and the result is a new representation.

        Parameters
        ----------
        op : `~operator` callable
            Operator to apply (e.g., `~operator.add`, `~operator.sub`, etc.
        other : `~astropy.coordinates.BaseRepresentation` instance
            The other differential or representation.
        reverse : bool
            Whether the operands should be reversed (e.g., as we got here via
            ``self.__rsub__`` because ``self`` is a subclass of ``other``).
        """
        if (isinstance(other, BaseSphericalDifferential) and
                not isinstance(self, type(other)) or
                isinstance(other, RadialDifferential)):
            all_components = set(self.components) | set(other.components)
            first, second = (self, other) if not reverse else (other, self)
            result_args = {c: op(getattr(first, c, 0.), getattr(second, c, 0.))
                           for c in all_components}
            return SphericalDifferential(**result_args)

        return super()._combine_operation(op, other, reverse)


class UnitSphericalDifferential(BaseSphericalDifferential):
    """Differential(s) of points on a unit sphere.

    Parameters
    ----------
    d_lon, d_lat : `~astropy.units.Quantity`
        The longitude and latitude of the differentials.
    copy : bool, optional
        If `True` (default), arrays will be copied rather than referenced.
    """
    base_representation = UnitSphericalRepresentation

    @classproperty
    def _dimensional_differential(cls):
        return SphericalDifferential

    def __init__(self, d_lon, d_lat, copy=True):
        super().__init__(d_lon, d_lat, copy=copy)
        if not self._d_lon.unit.is_equivalent(self._d_lat.unit):
            raise u.UnitsError('d_lon and d_lat should have equivalent units.')

    def to_cartesian(self, base):
        if isinstance(base, SphericalRepresentation):
            scale = base.distance
        elif isinstance(base, PhysicsSphericalRepresentation):
            scale = base.r
        else:
            return super().to_cartesian(base)

        base = base.represent_as(UnitSphericalRepresentation)
        return scale * super().to_cartesian(base)

    def represent_as(self, other_class, base=None):
        # Only have enough information to represent other unit-spherical.
        if issubclass(other_class, UnitSphericalCosLatDifferential):
            return other_class(self._d_lon_coslat(base), self.d_lat)

        return super().represent_as(other_class, base)

    @classmethod
    def from_representation(cls, representation, base=None):
        # All spherical differentials can be done without going to Cartesian,
        # though CosLat needs base for the latitude.
        if isinstance(representation, SphericalDifferential):
            return cls(representation.d_lon, representation.d_lat)
        elif isinstance(representation, (SphericalCosLatDifferential,
                                         UnitSphericalCosLatDifferential)):
            d_lon = cls._get_d_lon(representation.d_lon_coslat, base)
            return cls(d_lon, representation.d_lat)
        elif isinstance(representation, PhysicsSphericalDifferential):
            return cls(representation.d_phi, -representation.d_theta)

        return super().from_representation(representation, base)


class SphericalDifferential(BaseSphericalDifferential):
    """Differential(s) of points in 3D spherical coordinates.

    Parameters
    ----------
    d_lon, d_lat : `~astropy.units.Quantity`
        The differential longitude and latitude.
    d_distance : `~astropy.units.Quantity`
        The differential distance.
    copy : bool, optional
        If `True` (default), arrays will be copied rather than referenced.
    """
    base_representation = SphericalRepresentation
    _unit_differential = UnitSphericalDifferential

    def __init__(self, d_lon, d_lat, d_distance, copy=True):
        super().__init__(d_lon, d_lat, d_distance, copy=copy)
        if not self._d_lon.unit.is_equivalent(self._d_lat.unit):
            raise u.UnitsError('d_lon and d_lat should have equivalent units.')

    def represent_as(self, other_class, base=None):
        # All spherical differentials can be done without going to Cartesian,
        # though CosLat needs base for the latitude.
        if issubclass(other_class, UnitSphericalDifferential):
            return other_class(self.d_lon, self.d_lat)
        elif issubclass(other_class, RadialDifferential):
            return other_class(self.d_distance)
        elif issubclass(other_class, SphericalCosLatDifferential):
            return other_class(self._d_lon_coslat(base), self.d_lat,
                               self.d_distance)
        elif issubclass(other_class, UnitSphericalCosLatDifferential):
            return other_class(self._d_lon_coslat(base), self.d_lat)
        elif issubclass(other_class, PhysicsSphericalDifferential):
            return other_class(self.d_lon, -self.d_lat, self.d_distance)
        else:
            return super().represent_as(other_class, base)

    @classmethod
    def from_representation(cls, representation, base=None):
        # Other spherical differentials can be done without going to Cartesian,
        # though CosLat needs base for the latitude.
        if isinstance(representation, SphericalCosLatDifferential):
            d_lon = cls._get_d_lon(representation.d_lon_coslat, base)
            return cls(d_lon, representation.d_lat, representation.d_distance)
        elif isinstance(representation, PhysicsSphericalDifferential):
            return cls(representation.d_phi, -representation.d_theta,
                       representation.d_r)

        return super().from_representation(representation, base)


class BaseSphericalCosLatDifferential(BaseDifferential):
    """Differentials from points on a spherical base representation.

    With cos(lat) assumed to be included in the longitude differential.
    """
    @classmethod
    def _get_base_vectors(cls, base):
        """Get unit vectors and scale factors from (unit)spherical base.

        Parameters
        ----------
        base : instance of ``self.base_representation``
            The points for which the unit vectors and scale factors should be
            retrieved.

        Returns
        -------
        unit_vectors : dict of `CartesianRepresentation`
            In the directions of the coordinates of base.
        scale_factors : dict of `~astropy.units.Quantity`
            Scale factors for each of the coordinates.  The scale factor for
            longitude does not include the cos(lat) factor.

        Raises
        ------
        TypeError : if the base is not of the correct type
        """
        cls._check_base(base)
        return base.unit_vectors(), base.scale_factors(omit_coslat=True)

    def _d_lon(self, base):
        """Convert longitude differential with cos(lat) to one without.

        Parameters
        ----------
        base : instance of ``cls.base_representation``
            The base from which the latitude will be taken.
        """
        self._check_base(base)
        return self.d_lon_coslat / np.cos(base.lat)

    @classmethod
    def _get_d_lon_coslat(cls, d_lon, base):
        """Convert longitude differential d_lon to d_lon_coslat.

        Parameters
        ----------
        d_lon : `~astropy.units.Quantity`
            Value of the longitude differential without ``cos(lat)``.
        base : instance of ``cls.base_representation``
            The base from which the latitude will be taken.
        """
        cls._check_base(base)
        return d_lon * np.cos(base.lat)

    def _combine_operation(self, op, other, reverse=False):
        """Combine two differentials, or a differential with a representation.

        If ``other`` is of the same differential type as ``self``, the
        components will simply be combined.  If both are different parts of
        a `~astropy.coordinates.SphericalDifferential` (e.g., a
        `~astropy.coordinates.UnitSphericalDifferential` and a
        `~astropy.coordinates.RadialDifferential`), they will combined
        appropriately.

        If ``other`` is a representation, it will be used as a base for which
        to evaluate the differential, and the result is a new representation.

        Parameters
        ----------
        op : `~operator` callable
            Operator to apply (e.g., `~operator.add`, `~operator.sub`, etc.
        other : `~astropy.coordinates.BaseRepresentation` instance
            The other differential or representation.
        reverse : bool
            Whether the operands should be reversed (e.g., as we got here via
            ``self.__rsub__`` because ``self`` is a subclass of ``other``).
        """
        if (isinstance(other, BaseSphericalCosLatDifferential) and
                not isinstance(self, type(other)) or
                isinstance(other, RadialDifferential)):
            all_components = set(self.components) | set(other.components)
            first, second = (self, other) if not reverse else (other, self)
            result_args = {c: op(getattr(first, c, 0.), getattr(second, c, 0.))
                           for c in all_components}
            return SphericalCosLatDifferential(**result_args)

        return super()._combine_operation(op, other, reverse)


class UnitSphericalCosLatDifferential(BaseSphericalCosLatDifferential):
    """Differential(s) of points on a unit sphere.

    Parameters
    ----------
    d_lon_coslat, d_lat : `~astropy.units.Quantity`
        The longitude and latitude of the differentials.
    copy : bool, optional
        If `True` (default), arrays will be copied rather than referenced.
    """
    base_representation = UnitSphericalRepresentation
    attr_classes = OrderedDict([('d_lon_coslat', u.Quantity),
                                ('d_lat', u.Quantity)])

    @classproperty
    def _dimensional_differential(cls):
        return SphericalCosLatDifferential

    def __init__(self, d_lon_coslat, d_lat, copy=True):
        super().__init__(d_lon_coslat, d_lat, copy=copy)
        if not self._d_lon_coslat.unit.is_equivalent(self._d_lat.unit):
            raise u.UnitsError('d_lon_coslat and d_lat should have equivalent '
                               'units.')

    def to_cartesian(self, base):
        if isinstance(base, SphericalRepresentation):
            scale = base.distance
        elif isinstance(base, PhysicsSphericalRepresentation):
            scale = base.r
        else:
            return super().to_cartesian(base)

        base = base.represent_as(UnitSphericalRepresentation)
        return scale * super().to_cartesian(base)

    def represent_as(self, other_class, base=None):
        # Only have enough information to represent other unit-spherical.
        if issubclass(other_class, UnitSphericalDifferential):
            return other_class(self._d_lon(base), self.d_lat)

        return super().represent_as(other_class, base)

    @classmethod
    def from_representation(cls, representation, base=None):
        # All spherical differentials can be done without going to Cartesian,
        # though w/o CosLat needs base for the latitude.
        if isinstance(representation, SphericalCosLatDifferential):
            return cls(representation.d_lon_coslat, representation.d_lat)
        elif isinstance(representation, (SphericalDifferential,
                                         UnitSphericalDifferential)):
            d_lon_coslat = cls._get_d_lon_coslat(representation.d_lon, base)
            return cls(d_lon_coslat, representation.d_lat)
        elif isinstance(representation, PhysicsSphericalDifferential):
            d_lon_coslat = cls._get_d_lon_coslat(representation.d_phi, base)
            return cls(d_lon_coslat, -representation.d_theta)

        return super().from_representation(representation, base)


class SphericalCosLatDifferential(BaseSphericalCosLatDifferential):
    """Differential(s) of points in 3D spherical coordinates.

    Parameters
    ----------
    d_lon_coslat, d_lat : `~astropy.units.Quantity`
        The differential longitude (with cos(lat) included) and latitude.
    d_distance : `~astropy.units.Quantity`
        The differential distance.
    copy : bool, optional
        If `True` (default), arrays will be copied rather than referenced.
    """
    base_representation = SphericalRepresentation
    _unit_differential = UnitSphericalCosLatDifferential
    attr_classes = OrderedDict([('d_lon_coslat', u.Quantity),
                                ('d_lat', u.Quantity),
                                ('d_distance', u.Quantity)])

    def __init__(self, d_lon_coslat, d_lat, d_distance, copy=True):
        super().__init__(d_lon_coslat, d_lat, d_distance, copy=copy)
        if not self._d_lon_coslat.unit.is_equivalent(self._d_lat.unit):
            raise u.UnitsError('d_lon_coslat and d_lat should have equivalent '
                               'units.')

    def represent_as(self, other_class, base=None):
        # All spherical differentials can be done without going to Cartesian,
        # though some need base for the latitude to remove cos(lat).
        if issubclass(other_class, UnitSphericalCosLatDifferential):
            return other_class(self.d_lon_coslat, self.d_lat)
        elif issubclass(other_class, RadialDifferential):
            return other_class(self.d_distance)
        elif issubclass(other_class, SphericalDifferential):
            return other_class(self._d_lon(base), self.d_lat, self.d_distance)
        elif issubclass(other_class, UnitSphericalDifferential):
            return other_class(self._d_lon(base), self.d_lat)
        elif issubclass(other_class, PhysicsSphericalDifferential):
            return other_class(self._d_lon(base), -self.d_lat, self.d_distance)

        return super().represent_as(other_class, base)

    @classmethod
    def from_representation(cls, representation, base=None):
        # Other spherical differentials can be done without going to Cartesian,
        # though we need base for the latitude to remove coslat.
        if isinstance(representation, SphericalDifferential):
            d_lon_coslat = cls._get_d_lon_coslat(representation.d_lon, base)
            return cls(d_lon_coslat, representation.d_lat,
                       representation.d_distance)
        elif isinstance(representation, PhysicsSphericalDifferential):
            d_lon_coslat = cls._get_d_lon_coslat(representation.d_phi, base)
            return cls(d_lon_coslat, -representation.d_theta,
                       representation.d_r)

        return super().from_representation(representation, base)


class RadialDifferential(BaseDifferential):
    """Differential(s) of radial distances.

    Parameters
    ----------
    d_distance : `~astropy.units.Quantity`
        The differential distance.
    copy : bool, optional
        If `True` (default), arrays will be copied rather than referenced.
    """
    base_representation = RadialRepresentation

    def to_cartesian(self, base):
        return self.d_distance * base.represent_as(
            UnitSphericalRepresentation).to_cartesian()

    @classmethod
    def from_cartesian(cls, other, base):
        return cls(other.dot(base.represent_as(UnitSphericalRepresentation)),
                   copy=False)

    @classmethod
    def from_representation(cls, representation, base=None):
        if isinstance(representation, (SphericalDifferential,
                                       SphericalCosLatDifferential)):
            return cls(representation.d_distance)
        elif isinstance(representation, PhysicsSphericalDifferential):
            return cls(representation.d_r)
        else:
            return super().from_representation(representation, base)

    def _combine_operation(self, op, other, reverse=False):
        if isinstance(other, self.base_representation):
            if reverse:
                first, second = other.distance, self.d_distance
            else:
                first, second = self.d_distance, other.distance
            return other.__class__(op(first, second), copy=False)
        elif isinstance(other, (BaseSphericalDifferential,
                                BaseSphericalCosLatDifferential)):
            all_components = set(self.components) | set(other.components)
            first, second = (self, other) if not reverse else (other, self)
            result_args = {c: op(getattr(first, c, 0.), getattr(second, c, 0.))
                           for c in all_components}
            return SphericalDifferential(**result_args)

        else:
            return super()._combine_operation(op, other, reverse)


class PhysicsSphericalDifferential(BaseDifferential):
    """Differential(s) of 3D spherical coordinates using physics convention.

    Parameters
    ----------
    d_phi, d_theta : `~astropy.units.Quantity`
        The differential azimuth and inclination.
    d_r : `~astropy.units.Quantity`
        The differential radial distance.
    copy : bool, optional
        If `True` (default), arrays will be copied rather than referenced.
    """
    base_representation = PhysicsSphericalRepresentation

    def __init__(self, d_phi, d_theta, d_r, copy=True):
        super().__init__(d_phi, d_theta, d_r, copy=copy)
        if not self._d_phi.unit.is_equivalent(self._d_theta.unit):
            raise u.UnitsError('d_phi and d_theta should have equivalent '
                               'units.')

    def represent_as(self, other_class, base=None):
        # All spherical differentials can be done without going to Cartesian,
        # though CosLat needs base for the latitude. For those, explicitly
        # do the equivalent of self._d_lon_coslat in SphericalDifferential.
        if issubclass(other_class, SphericalDifferential):
            return other_class(self.d_phi, -self.d_theta, self.d_r)
        elif issubclass(other_class, UnitSphericalDifferential):
            return other_class(self.d_phi, -self.d_theta)
        elif issubclass(other_class, SphericalCosLatDifferential):
            self._check_base(base)
            d_lon_coslat = self.d_phi * np.sin(base.theta)
            return other_class(d_lon_coslat, -self.d_theta, self.d_r)
        elif issubclass(other_class, UnitSphericalCosLatDifferential):
            self._check_base(base)
            d_lon_coslat = self.d_phi * np.sin(base.theta)
            return other_class(d_lon_coslat, -self.d_theta)
        elif issubclass(other_class, RadialDifferential):
            return other_class(self.d_r)

        return super().represent_as(other_class, base)

    @classmethod
    def from_representation(cls, representation, base=None):
        # Other spherical differentials can be done without going to Cartesian,
        # though we need base for the latitude to remove coslat. For that case,
        # do the equivalent of cls._d_lon in SphericalDifferential.
        if isinstance(representation, SphericalDifferential):
            return cls(representation.d_lon, -representation.d_lat,
                       representation.d_distance)
        elif isinstance(representation, SphericalCosLatDifferential):
            cls._check_base(base)
            d_phi = representation.d_lon_coslat / np.sin(base.theta)
            return cls(d_phi, -representation.d_lat, representation.d_distance)

        return super().from_representation(representation, base)


class CylindricalDifferential(BaseDifferential):
    """Differential(s) of points in cylindrical coordinates.

    Parameters
    ----------
    d_rho : `~astropy.units.Quantity`
        The differential cylindrical radius.
    d_phi : `~astropy.units.Quantity`
        The differential azimuth.
    d_z : `~astropy.units.Quantity`
        The differential height.
    copy : bool, optional
        If `True` (default), arrays will be copied rather than referenced.
    """
    base_representation = CylindricalRepresentation

    def __init__(self, d_rho, d_phi, d_z, copy=False):
        super().__init__(d_rho, d_phi, d_z, copy=copy)
        if not self._d_rho.unit.is_equivalent(self._d_z.unit):
            raise u.UnitsError("d_rho and d_z should have equivalent units.")
