"""
In this module, we define the coordinate representation classes, which are
used to represent low-level cartesian, spherical, cylindrical, and other
coordinates.
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import abc
import functools
import operator
from collections import OrderedDict

import numpy as np
import astropy.units as u

from .angles import Angle, Longitude, Latitude
from .distances import Distance
from ..extern import six
from ..utils import ShapedLikeNDArray, classproperty
from ..utils.misc import InheritDocstrings
from ..utils.compat import NUMPY_LT_1_12
from ..utils.compat.numpy import broadcast_arrays, broadcast_to
from ..utils.misc import isiterable, check_broadcast, IncompatibleShapeError

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

def _array2string(values, prefix=''):
    # Mimic numpy >=1.12 array2string, in which structured arrays are
    # typeset taking into account all printoptions.
    if NUMPY_LT_1_12:  # pragma: no cover
        # Mimic StructureFormat from numpy >=1.12 assuming float-only data.
        from numpy.core.arrayprint import FloatFormat
        opts = np.get_printoptions()
        format_functions = [FloatFormat(np.atleast_1d(values[component]).ravel(),
                                        precision=opts['precision'],
                                        suppress_small=opts['suppress'])
                            for component in values.dtype.names]

        def fmt(x):
            return '({})'.format(', '.join(format_function(field)
                                           for field, format_function in
                                           zip(x, format_functions)))
        # Before 1.12, structures arrays were set as "numpystr",
        # so that is the formmater we need to replace.
        formatter = {'numpystr': fmt}
    else:
        fmt = repr
        formatter = {}

    return np.array2string(values, formatter=formatter, style=fmt,
                           separator=', ', prefix=prefix)

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
    # Add new axis in x, y, z so one can concatenate them around it.
    # NOTE: just use np.stack once our minimum numpy version is 1.10.
    result_ndim = x.ndim + 1
    if not -result_ndim <= xyz_axis < result_ndim:
        raise IndexError('xyz_axis {0} out of bounds [-{1}, {1})'
                         .format(xyz_axis, result_ndim))

    if xyz_axis < 0:
        xyz_axis += result_ndim

    # Get x, y, z to the same units (this is very fast for identical units)
    # since np.concatenate cannot deal with quantity.
    cls = x.__class__
    y = cls(y, x.unit, copy=False)
    z = cls(z, x.unit, copy=False)

    sh = x.shape
    sh = sh[:xyz_axis] + (1,) + sh[xyz_axis:]
    xyz_value = np.concatenate([c.reshape(sh).value for c in (x, y, z)],
                               axis=xyz_axis)
    return cls(xyz_value, unit=x.unit, copy=False)


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

        if args:
            differentials = args.pop(0)
        else:
            differentials = kwargs.pop('differentials', tuple())

        if differentials is None:
            differentials = tuple()

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
            attrs = broadcast_arrays(*attrs, subok=True)
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

        # Now handle the actual validation of any specified differential classes
        if isinstance(differentials, BaseDifferential):
            differentials = (differentials, )

        for diff in differentials:
            if not isinstance(diff, BaseDifferential):
                raise TypeError("'differentials' must be an iterable of "
                                "BaseDifferential subclass instances. Got '{0}'"
                                .format(type(diff)))

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

        # store as a private name so users know this is not settable
        self._differentials = tuple(differentials)

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
    # Should be replaced by abstractclassmethod once we support only PY3
    @abc.abstractmethod
    def from_cartesian(self):
        raise NotImplementedError()

    @abc.abstractmethod
    def to_cartesian(self):
        raise NotImplementedError()

    @property
    def components(self):
        """A tuple with the in-order names of the coordinate components."""
        return tuple(self.attr_classes)

    @property
    def differentials(self):
        """A tuple of differential class instances"""
        return self._differentials

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

        return self.__class__(*[apply_method(getattr(self, component))
                                for component in self.components], copy=False)

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

        # also perform shape-setting on any associated differentials
        for diff in self.differentials:
            diff.shape = shape

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
        # The "str(c)" is needed for PY2; it can be removed for astropy 3.0.
        coo_items = [(str(c), getattr(self, c)) for c in self.components]
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


        unitstr = ('in ' + self._unitstr) if self._unitstr else '[dimensionless]'
        return '<{0} ({1}) {2:s}\n{3}{4}>'.format(
            self.__class__.__name__, ', '.join(self.components),
            unitstr, prefixstr, arrstr)


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
        super(MetaBaseRepresentation, cls).__init__(name, bases, dct)

        # Register representation name (except for BaseRepresentation)
        if cls.__name__ == 'BaseRepresentation':
            return

        if 'attr_classes' not in dct:
            raise NotImplementedError('Representations must have an '
                                      '"attr_classes" class attribute.')

        repr_name = cls.get_name()

        if repr_name in REPRESENTATION_CLASSES:
            raise ValueError("Representation class {0} already defined"
                             .format(repr_name))

        REPRESENTATION_CLASSES[repr_name] = cls

        # define getters for any component that does not yet have one.
        for component in cls.attr_classes:
            if not hasattr(cls, component):
                setattr(cls, component,
                        property(_make_getter(component),
                                 doc=("The '{0}' component of the points(s)."
                                      .format(component))))


@six.add_metaclass(MetaBaseRepresentation)
class BaseRepresentation(BaseRepresentationOrDifferential):
    """Base for representing a point in a 3D coordinate system.

    Parameters
    ----------
    comp1, comp2, comp3 : `~astropy.units.Quantity` or subclass
        The components of the 3D points.  The names are the keys and the
        subclasses the values of the ``attr_classes`` attribute.
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
    methods (see those methods for details). Finally, classes can also define a
    ``recommended_units`` dictionary, which maps component names to the units
    they are best presented to users in (this is used only in representations
    of coordinates, and may be overridden by frame classes).
    """

    recommended_units = {}  # subclasses can override

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

    def represent_as(self, other_class):
        """Convert coordinates to another representation.

        If the instance is of the requested class, it is returned unmodified.
        By default, conversion is done via cartesian coordinates.

        Parameters
        ----------
        other_class : `~astropy.coordinates.BaseRepresentation` subclass
            The type of representation to turn the coordinates into.
        """
        if other_class is self.__class__:
            return self
        else:
            # The default is to convert via cartesian coordinates
            if isiterable(other_class):
                if (self.differentials is None or
                        len(other_class) != 1+len(self.differentials)):
                    raise ValueError("If multiple classes are passed in to "
                                     "change the representations of attached "
                                     "differentials, the number of classes must "
                                     "be equal to 1 + the number of "
                                     "differentials. (got {0})"
                                     .format(len(other_class)))

                new_rep = other_class[0].from_cartesian(self.to_cartesian())

                new_diffs = []
                for i,diff in enumerate(self.differentials):
                    new_diffs.append(diff.represent_as(other_class[1+i],
                                                       base=new_rep))
                new_rep._differentials = tuple(new_diffs)

            else:
                new_rep = other_class.from_cartesian(self.to_cartesian())

            return new_rep

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
        rep = super(BaseRepresentation, self)._apply(method, *args, **kwargs)
        diffs = [diff._apply(method, *args, **kwargs)
                 for diff in self.differentials]
        rep._differentials = tuple(diffs)
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
        result = self.to_cartesian()._combine_operation(op, other, reverse)
        if result is NotImplemented:
            return NotImplemented
        else:
            return self.from_cartesian(result)

    def norm(self):
        """Vector norm.

        The norm is the standard Frobenius norm, i.e., the square root of the
        sum of the squares of all components with non-angular units.

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
        return self.from_cartesian(self.to_cartesian().sum(*args, **kwargs))

    def dot(self, other):
        """Dot product of two representations.

        The calculation is done by converting both ``self`` and ``other``
        to `~astropy.coordinates.CartesianRepresentation`.

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
    copy : bool, optional
        If `True` (default), arrays will be copied rather than referenced.
    """

    attr_classes = OrderedDict([('x', u.Quantity),
                                ('y', u.Quantity),
                                ('z', u.Quantity)])

    def __init__(self, x, y=None, z=None, unit=None, differentials=None,
                 xyz_axis=None, copy=True):

        if y is None and z is None:
            if xyz_axis is not None and xyz_axis != 0:
                x = np.rollaxis(x, xyz_axis, 0)
            x, y, z = x
        elif xyz_axis is not None:
            raise ValueError("xyz_axis should only be set if x, y, and z are "
                             "in a single array passed in through x, "
                             "i.e., y and z should not be not given.")
        elif (y is None and z is not None) or (y is not None and z is None):
            raise ValueError("x, y, and z are required to instantiate {0}"
                             .format(self.__class__.__name__))

        if unit is not None:
            x = u.Quantity(x, unit, copy=copy, subok=True)
            y = u.Quantity(y, unit, copy=copy, subok=True)
            z = u.Quantity(z, unit, copy=copy, subok=True)
            copy = False

        super(CartesianRepresentation, self).__init__(x, y, z, copy=copy,
                                                      differentials=differentials)
        if not (self._x.unit.physical_type ==
                self._y.unit.physical_type == self._z.unit.physical_type):
            raise u.UnitsError("x, y, and z should have matching physical types")

    def unit_vectors(self):
        l = broadcast_to(1.*u.one, self.shape, subok=True)
        o = broadcast_to(0.*u.one, self.shape, subok=True)
        return OrderedDict(
            (('x', CartesianRepresentation(l, o, o, copy=False)),
             ('y', CartesianRepresentation(o, l, o, copy=False)),
             ('z', CartesianRepresentation(o, o, l, copy=False))))

    def scale_factors(self):
        l = broadcast_to(1.*u.one, self.shape, subok=True)
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
            With dimension 3 along ``xyz_axis``.
        """
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
        # Avoid doing gratuitous np.array for things that look like arrays.
        try:
            matrix_shape = matrix.shape
        except AttributeError:
            matrix = np.array(matrix)
            matrix_shape = matrix.shape

        if matrix_shape[-2:] != (3, 3):
            raise ValueError("tried to do matrix multiplication with an array "
                             "that doesn't end in 3x3")

        # TODO: since this is likely to be a widely used function in coordinate
        # transforms, it should be optimized (for example in Cython).

        # Get xyz once since it's an expensive operation
        oldxyz = self.xyz
        # Note that neither dot nor einsum handles Quantity properly, so we use
        # the arrays and put the unit back in the end.
        if self.isscalar and not matrix_shape[:-2]:
            # a fast path for scalar coordinates.
            newxyz = matrix.dot(oldxyz.value)
        else:
            # Matrix multiply all pmat items and coordinates, broadcasting the
            # remaining dimensions.
            newxyz = np.einsum('...ij,j...->i...', matrix, oldxyz.value)

        newxyz = u.Quantity(newxyz, oldxyz.unit, copy=False)
        return self.__class__(*newxyz, copy=False)

    def _combine_operation(self, op, other, reverse=False):
        try:
            other_c = other.to_cartesian()
        except Exception:
            return NotImplemented

        first, second = ((self, other_c) if not reverse else
                         (other_c, self))
        return self.__class__(*(op(getattr(first, component),
                                   getattr(second, component))
                                for component in first.components))

    def mean(self, *args, **kwargs):
        """Vector mean.

        Returns a new CartesianRepresentation instance with the means of the
        x, y, and z components.

        Refer to `~numpy.mean` for full documentation of the arguments, noting
        that ``axis`` is the entry in the ``shape`` of the representation, and
        that the ``out`` argument cannot be used.
        """
        return self._apply('mean', *args, **kwargs)

    def sum(self, *args, **kwargs):
        """Vector sum.

        Returns a new CartesianRepresentation instance with the sums of the
        x, y, and z components.

        Refer to `~numpy.sum` for full documentation of the arguments, noting
        that ``axis`` is the entry in the ``shape`` of the representation, and
        that the ``out`` argument cannot be used.
        """
        return self._apply('sum', *args, **kwargs)

    def dot(self, other):
        """Dot product of two representations.

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
        return functools.reduce(operator.add,
                                (getattr(self, component) *
                                 getattr(other_c, component)
                                 for component in self.components))

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
        try:
            other_c = other.to_cartesian()
        except Exception:
            raise TypeError("cannot only take cross product with another "
                            "representation, not a {0} instance."
                            .format(type(other)))
        return self.__class__(self.y * other_c.z - self.z * other_c.y,
                              self.z * other_c.x - self.x * other_c.z,
                              self.x * other_c.y - self.y * other_c.x)


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

    copy : bool, optional
        If `True` (default), arrays will be copied rather than referenced.
    """

    attr_classes = OrderedDict([('lon', Longitude),
                                ('lat', Latitude)])
    recommended_units = {'lon': u.deg, 'lat': u.deg}

    @classproperty
    def _dimensional_representation(cls):
        return SphericalRepresentation

    def __init__(self, lon, lat, differentials=None, copy=True):
        super(UnitSphericalRepresentation,
              self).__init__(lon, lat, copy=copy, differentials=differentials)

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
        sf_lat = broadcast_to(1./u.radian, self.shape, subok=True)
        sf_lon  = sf_lat if omit_coslat else np.cos(self.lat) / u.radian
        return OrderedDict((('lon', sf_lon),
                            ('lat', sf_lat)))

    def to_cartesian(self):
        """
        Converts spherical polar coordinates to 3D rectangular cartesian
        coordinates.
        """
        x = np.cos(self.lat) * np.cos(self.lon)
        y = np.cos(self.lat) * np.sin(self.lon)
        z = np.sin(self.lat)

        return CartesianRepresentation(x=x, y=y, z=z, copy=False,
                                       differentials=self.differentials)

    @classmethod
    def from_cartesian(cls, cart):
        """
        Converts 3D rectangular cartesian coordinates to spherical polar
        coordinates.
        """

        s = np.hypot(cart.x, cart.y)

        lon = np.arctan2(cart.y, cart.x)
        lat = np.arctan2(cart.z, s)

        return cls(lon=lon, lat=lat, copy=False,
                   differentials=cart.differentials)

    def represent_as(self, other_class):
        # Take a short cut if the other class is a spherical representation
        if issubclass(other_class, PhysicsSphericalRepresentation):
            return other_class(phi=self.lon, theta=90 * u.deg - self.lat, r=1.0,
                               copy=False)
        elif issubclass(other_class, SphericalRepresentation):
            return other_class(lon=self.lon, lat=self.lat, distance=1.0,
                               copy=False)
        else:
            return super(UnitSphericalRepresentation,
                         self).represent_as(other_class)

    def __mul__(self, other):
        return self._dimensional_representation(lon=self.lon, lat=self.lat,
                                                distance=1. * other)

    def __truediv__(self, other):
        return self._dimensional_representation(lon=self.lon, lat=self.lat,
                                                distance=1. / other)

    def __neg__(self):
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

    copy : bool, optional
        If `True` (default), arrays will be copied rather than referenced.
    """

    attr_classes = OrderedDict([('distance', u.Quantity)])

    def __init__(self, distance, differentials=None, copy=True):
        super(RadialRepresentation, self).__init__(distance, copy=copy,
                                                   differentials=differentials)

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
        l = broadcast_to(1.*u.one, self.shape, subok=True)
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
        return cls(distance=cart.norm(), copy=False,
                   differentials=cart.differentials)

    def _scale_operation(self, op, *args):
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

    copy : bool, optional
        If `True` (default), arrays will be copied rather than referenced.
    """

    attr_classes = OrderedDict([('lon', Longitude),
                                ('lat', Latitude),
                                ('distance', u.Quantity)])
    recommended_units = {'lon': u.deg, 'lat': u.deg}
    _unit_representation = UnitSphericalRepresentation

    def __init__(self, lon, lat, distance, differentials=None, copy=True):
        super(SphericalRepresentation,
              self).__init__(lon, lat, distance, copy=copy,
                             differentials=differentials)
        if self._distance.unit.physical_type == 'length':
            self._distance = self._distance.view(Distance)

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
        sf_distance = broadcast_to(1.*u.one, self.shape, subok=True)
        return OrderedDict((('lon', sf_lon),
                            ('lat', sf_lat),
                            ('distance', sf_distance)))

    def represent_as(self, other_class):
        # Take a short cut if the other class is a spherical representation
        if issubclass(other_class, PhysicsSphericalRepresentation):
            return other_class(phi=self.lon, theta=90 * u.deg - self.lat,
                               r=self.distance, copy=False)
        elif issubclass(other_class, UnitSphericalRepresentation):
            return other_class(lon=self.lon, lat=self.lat, copy=False)
        else:
            return super(SphericalRepresentation,
                         self).represent_as(other_class)

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

        x = d * np.cos(self.lat) * np.cos(self.lon)
        y = d * np.cos(self.lat) * np.sin(self.lon)
        z = d * np.sin(self.lat)

        return CartesianRepresentation(x=x, y=y, z=z, copy=False,
                                       differentials=self.differentials)

    @classmethod
    def from_cartesian(cls, cart):
        """
        Converts 3D rectangular cartesian coordinates to spherical polar
        coordinates.
        """

        s = np.hypot(cart.x, cart.y)
        r = np.hypot(s, cart.z)

        lon = np.arctan2(cart.y, cart.x)
        lat = np.arctan2(cart.z, s)

        return cls(lon=lon, lat=lat, distance=r, copy=False,
                   differentials=cart.differentials)

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

    copy : bool, optional
        If `True` (default), arrays will be copied rather than referenced.
    """

    attr_classes = OrderedDict([('phi', Angle),
                                ('theta', Angle),
                                ('r', u.Quantity)])
    recommended_units = {'phi': u.deg, 'theta': u.deg}

    def __init__(self, phi, theta, r, differentials=None, copy=True):
        super(PhysicsSphericalRepresentation,
              self).__init__(phi, theta, r, copy=copy,
                             differentials=differentials)

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
        l = broadcast_to(1.*u.one, self.shape, subok=True)
        return OrderedDict((('phi', r * sintheta),
                            ('theta', r),
                            ('r', l)))

    def represent_as(self, other_class):
        # Take a short cut if the other class is a spherical representation
        if issubclass(other_class, SphericalRepresentation):
            return other_class(lon=self.phi, lat=90 * u.deg - self.theta,
                               distance=self.r)
        elif issubclass(other_class, UnitSphericalRepresentation):
            return other_class(lon=self.phi, lat=90 * u.deg - self.theta)
        else:
            return super(PhysicsSphericalRepresentation, self).represent_as(other_class)

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

        return CartesianRepresentation(x=x, y=y, z=z, copy=False,
                                       differentials=self.differentials)

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

        return cls(phi=phi, theta=theta, r=r, copy=False,
                   differentials=cart.differentials)

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

    copy : bool, optional
        If `True` (default), arrays will be copied rather than referenced.
    """

    attr_classes = OrderedDict([('rho', u.Quantity),
                                ('phi', Angle),
                                ('z', u.Quantity)])
    recommended_units = {'phi': u.deg}

    def __init__(self, rho, phi, z, differentials=None, copy=True):
        super(CylindricalRepresentation,
              self).__init__(rho, phi, z, copy=copy,
                             differentials=differentials)

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
        l = broadcast_to(1., self.shape)
        return OrderedDict(
            (('rho', CartesianRepresentation(cosphi, sinphi, 0, copy=False)),
             ('phi', CartesianRepresentation(-sinphi, cosphi, 0, copy=False)),
             ('z', CartesianRepresentation(0, 0, l, unit=u.one, copy=False))))

    def scale_factors(self):
        rho = self.rho / u.radian
        l = broadcast_to(1.*u.one, self.shape, subok=True)
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

        return cls(rho=rho, phi=phi, z=z, copy=False,
                   differentials=cart.differentials)

    def to_cartesian(self):
        """
        Converts cylindrical polar coordinates to 3D rectangular cartesian
        coordinates.
        """
        x = self.rho * np.cos(self.phi)
        y = self.rho * np.sin(self.phi)
        z = self.z

        return CartesianRepresentation(x=x, y=y, z=z, copy=False,
                                       differentials=self.differentials)


class MetaBaseDifferential(InheritDocstrings, abc.ABCMeta):
    """Set default ``attr_classes`` and component getters on a Differential.

    For these, the components are those of the base representation prefixed
    by 'd_', and the class is `~astropy.units.Quantity`.
    """
    def __init__(cls, name, bases, dct):
        super(MetaBaseDifferential, cls).__init__(name, bases, dct)

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

        repr_name = cls.get_name()
        if repr_name in DIFFERENTIAL_CLASSES:
            raise ValueError("Differential class {0} already defined"
                             .format(repr_name))

        DIFFERENTIAL_CLASSES[repr_name] = cls

        # If not defined explicitly, create properties for the components.
        for component in cls.attr_classes:
            if not hasattr(cls, component):
                setattr(cls, component,
                        property(_make_getter(component),
                                 doc=("Component '{0}' of the Differential."
                                      .format(component))))


@six.add_metaclass(MetaBaseDifferential)
class BaseDifferential(BaseRepresentationOrDifferential):
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

    @classmethod
    def _check_base(cls, base):
        if not isinstance(base, cls.base_representation):
            raise TypeError('need a base of the correct representation type, '
                            '{0}, not {1}.'.format(cls.base_representation,
                                                   type(base)))

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
        """
        base_e, base_sf = self._get_base_vectors(base)
        return functools.reduce(
            operator.add, (getattr(self, d_c) * base_sf[c] * base_e[c]
                           for d_c, c in zip(self.components, base.components)))

    @classmethod
    def from_cartesian(cls, other, base):
        """Interpret 3D rectangular cartesian coordinates as differentials.

        Parameters
        ----------
        base : instance of ``self.base_representation``
            Base relative to which the differentials are defined.
        """
        base_e, base_sf = cls._get_base_vectors(base)
        return cls(*(other.dot(e / base_sf[component])
                     for component, e in six.iteritems(base_e)), copy=False)

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
        return super(BaseDifferential, self).__sub__(other)

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
    def __init__(self, d_x, d_y=None, d_z=None, unit=None, xyz_axis=None,
                 copy=True):

        if d_y is None and d_z is None:
            if xyz_axis is not None and xyz_axis != 0:
                d_x = np.rollaxis(d_x, xyz_axis, 0)
            d_x, d_y, d_z = d_x
        elif xyz_axis is not None:
            raise ValueError("xyz_axis should only be set if d_x, d_y, and d_z "
                             "are in a single array passed in through d_x, "
                             "i.e., d_y and d_z should not be not given.")
        elif ((d_y is None and d_z is not None) or
              (d_y is not None and d_z is None)):
            raise ValueError("d_x, d_y, and d_z are required to instantiate {0}"
                             .format(self.__class__.__name__))

        if unit is not None:
            d_x = u.Quantity(d_x, unit, copy=copy, subok=True)
            d_y = u.Quantity(d_y, unit, copy=copy, subok=True)
            d_z = u.Quantity(d_z, unit, copy=copy, subok=True)
            copy = False

        super(CartesianDifferential, self).__init__(d_x, d_y, d_z, copy=copy)
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
        xyz : `~astropy.units.Quantity`
            With dimension 3 along ``xyz_axis``.
        """
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
            result_args = {c : op(getattr(first, c, 0.), getattr(second, c, 0.))
                           for c in all_components}
            return SphericalDifferential(**result_args)

        return super(BaseSphericalDifferential,
                     self)._combine_operation(op, other, reverse)


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

    def __init__(self, d_lon, d_lat, d_distance, copy=True):
        super(SphericalDifferential, self).__init__(d_lon, d_lat, d_distance,
                                                    copy=copy)
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
            return super(SphericalDifferential,
                         self).represent_as(other_class, base)

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

        return super(SphericalDifferential,
                     cls).from_representation(representation, base)


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

    def __init__(self, d_lon, d_lat, copy=True):
        super(UnitSphericalDifferential,
              self).__init__(d_lon, d_lat, copy=copy)
        if not self._d_lon.unit.is_equivalent(self._d_lat.unit):
            raise u.UnitsError('d_lon and d_lat should have equivalent units.')

    def to_cartesian(self, base):
        if isinstance(base, SphericalRepresentation):
            scale = base.distance
        elif isinstance(base, PhysicsSphericalRepresentation):
            scale = base.r
        else:
            return super(UnitSphericalDifferential, self).to_cartesian(base)

        base = base.represent_as(UnitSphericalRepresentation)
        return scale * super(UnitSphericalDifferential, self).to_cartesian(base)

    def represent_as(self, other_class, base=None):
        # Only have enough information to represent other unit-spherical.
        if issubclass(other_class, UnitSphericalCosLatDifferential):
            return other_class(self._d_lon_coslat(base), self.d_lat)

        return super(UnitSphericalDifferential,
                     self).represent_as(other_class, base)

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

        return super(UnitSphericalDifferential,
                     cls).from_representation(representation, base)


class BaseSphericalCosLatDifferential(BaseDifferential):
    """Differtials from points on a spherical base representation.

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
            result_args = {c : op(getattr(first, c, 0.), getattr(second, c, 0.))
                           for c in all_components}
            return SphericalCosLatDifferential(**result_args)

        return super(BaseSphericalCosLatDifferential,
                     self)._combine_operation(op, other, reverse)


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
    attr_classes = OrderedDict([('d_lon_coslat', u.Quantity),
                                ('d_lat', u.Quantity),
                                ('d_distance', u.Quantity)])

    def __init__(self, d_lon_coslat, d_lat, d_distance, copy=True):
        super(SphericalCosLatDifferential,
              self).__init__(d_lon_coslat, d_lat, d_distance, copy=copy)
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

        return super(SphericalCosLatDifferential,
                     self).represent_as(other_class, base)

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

        return super(SphericalCosLatDifferential,
                     cls).from_representation(representation, base)


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

    def __init__(self, d_lon_coslat, d_lat, copy=True):
        super(UnitSphericalCosLatDifferential,
              self).__init__(d_lon_coslat, d_lat, copy=copy)
        if not self._d_lon_coslat.unit.is_equivalent(self._d_lat.unit):
            raise u.UnitsError('d_lon_coslat and d_lat should have equivalent '
                               'units.')

    def to_cartesian(self, base):
        if isinstance(base, SphericalRepresentation):
            scale = base.distance
        elif isinstance(base, PhysicsSphericalRepresentation):
            scale = base.r
        else:
            return super(UnitSphericalCosLatDifferential,
                         self).to_cartesian(base)

        base = base.represent_as(UnitSphericalRepresentation)
        return scale * super(UnitSphericalCosLatDifferential,
                             self).to_cartesian(base)

    def represent_as(self, other_class, base=None):
        # Only have enough information to represent other unit-spherical.
        if issubclass(other_class, UnitSphericalDifferential):
            return other_class(self._d_lon(base), self.d_lat)

        return super(UnitSphericalCosLatDifferential,
                     self).represent_as(other_class, base)

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

        return super(UnitSphericalDifferential,
                     cls).from_representation(representation, base)


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
            return super(RadialDifferential,
                         cls).from_representation(representation, base)

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
            result_args = {c : op(getattr(first, c, 0.), getattr(second, c, 0.))
                           for c in all_components}
            return SphericalDifferential(**result_args)

        else:
            return super(RadialDifferential,
                         self)._combine_operation(op, other, reverse)


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
        super(PhysicsSphericalDifferential,
              self).__init__(d_phi, d_theta, d_r, copy=copy)
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

        return super(PhysicsSphericalDifferential,
                     self).represent_as(other_class, base)

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

        return super(PhysicsSphericalDifferential,
                     cls).from_representation(representation, base)


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
        super(CylindricalDifferential,
              self).__init__(d_rho, d_phi, d_z, copy=copy)
        if not self._d_rho.unit.is_equivalent(self._d_z.unit):
            raise u.UnitsError("d_rho and d_z should have equivalent units.")
