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
from ..utils import ShapedLikeNDArray, classproperty, sharedmethod
from ..utils.compat import NUMPY_LT_1_12
from ..utils.compat.numpy import broadcast_arrays, broadcast_to

# Note: Offset classes gets added to this at the end.
__all__ = ["BaseRepresentation", "CartesianRepresentation",
           "SphericalRepresentation", "UnitSphericalRepresentation",
           "PhysicsSphericalRepresentation", "CylindricalRepresentation"]

# Module-level dict mapping representation string alias names to class.
# This is populated by the metaclass init so all representation classes
# get registered automatically.
REPRESENTATION_CLASSES = {}

def _array2string(values, prefix=''):
    # Mimic numpy >=1.12 array2string, in which structured arrays are
    # typeset taking into account all printoptions.
    # TODO: in final numpy 1.12, the scalar case should work as well;
    # see https://github.com/numpy/numpy/issues/8172
    if NUMPY_LT_1_12:
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


class RepresentationBase(ShapedLikeNDArray):
    """Base for 3D coordinate representations and their offsets."""

    # Ensure multiplication/division with ndarray or Quantity doesn't lead to
    # object arrays.
    __array_priority__ = 50000

    # Should be replaced by abstractclassmethod once we support only Python 3
    @abc.abstractmethod
    def from_cartesian(self):
        raise NotImplementedError()

    @abc.abstractmethod
    def to_cartesian(self):
        raise NotImplementedError()

    @property
    def components(self):
        """A tuple with the in-order names of the coordinate components"""
        return tuple(self.attr_classes)

    def _apply(self, method, *args, **kwargs):
        """Create a new representation with ``method`` applied to the arrays.

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
        return self.__class__( *[apply_method(getattr(self, component))
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

    @abc.abstractmethod
    def _scale_operation(self, op, *args):
        raise NotImplementedError()

    def __mul__(self, other):
        return self._scale_operation(operator.mul, other)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __truediv__(self, other):
        return self._scale_operation(operator.truediv, other)

    def __div__(self, other):
        return self._scale_operation(operator.truediv, other)

    def __neg__(self):
        return self._scale_operation(operator.neg)

    # Follow numpy convention and make an independent copy.
    def __pos__(self):
        return self.copy()

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


# Need to subclass ABCMeta rather than type, so that this meta class can be
# combined with a ShapedLikeNDArray subclass (which is an ABC).  Without it:
# "TypeError: metaclass conflict: the metaclass of a derived class must be a
#  (non-strict) subclass of the metaclasses of all its bases"
class MetaBaseRepresentation(abc.ABCMeta):
    def __init__(cls, name, bases, dct):
        super(MetaBaseRepresentation, cls).__init__(name, bases, dct)

        # Register representation name (except for BaseRepresentation)
        if cls.__name__ == 'BaseRepresentation':
            return

        if name != 'BaseRepresentation' and 'attr_classes' not in dct:
            raise NotImplementedError('Representations must have an '
                                      '"attr_classes" class attribute.')

        repr_name = cls.get_name()

        if repr_name in REPRESENTATION_CLASSES:
            raise ValueError("Representation class {0} already defined".format(repr_name))

        REPRESENTATION_CLASSES[repr_name] = cls


@six.add_metaclass(MetaBaseRepresentation)
class BaseRepresentation(RepresentationBase):
    """Base for representing a point in a 3D coordinate system.

    Notes
    -----
    All representation classes should subclass this base representation
    class. All subclasses should then define a ``to_cartesian`` method and a
    ``from_cartesian`` class method. By default, transformations are done via
    the cartesian system, but classes that want to define a smarter
    transformation path can overload the ``represent_as`` method.
    Furthermore, all classes must define an ``attr_classes`` attribute, an
    `~collections.OrderedDict` which maps component names to the class that
    creates them.  They can also define a `recommended_units` dictionary, which
    maps component names to the units they are best presented to users in.  Note
    that frame classes may override this with their own preferred units.
    """

    recommended_units = {}  # subclasses can override

    def represent_as(self, other_class):
        if other_class == self.__class__:
            return self
        else:
            # The default is to convert via cartesian coordinates
            return other_class.from_cartesian(self.to_cartesian())

    @classmethod
    def from_representation(cls, representation):
        return representation.represent_as(cls)

    @classmethod
    def get_name(cls):
        name = cls.__name__.lower()
        if name.endswith('representation'):
            name = name[:-14]
        return name

    def _scale_operation(self, op, *args):
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
        contains an array with the 3 coordinates are stored along ``xyz_axis``.
    unit : `~astropy.units.Unit` or str
        If given, the coordinates will be converted to this unit (or taken to
        be in this unit if not given.
    xyz_axis : int, optional
        The axis along which the coordinates are stored when a single array is
        provided rather than distinct ``x``, ``y``, and ``z`` (default: 0).
    copy : bool, optional
        If True (default), arrays will be copied rather than referenced.
    """

    attr_classes = OrderedDict([('x', u.Quantity),
                                ('y', u.Quantity),
                                ('z', u.Quantity)])

    def __init__(self, x, y=None, z=None, unit=None, xyz_axis=None, copy=True):

        if unit is None and not hasattr(x, 'unit'):
            raise TypeError('x should have a unit unless an explicit unit '
                            'is passed in.')

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

        try:
            x = self.attr_classes['x'](x, unit=unit, copy=copy)
            y = self.attr_classes['y'](y, unit=unit, copy=copy)
            z = self.attr_classes['z'](z, unit=unit, copy=copy)
        except u.UnitsError:
            raise u.UnitsError('x, y, and z should have a unit consistent with '
                               '{0}'.format(unit))
        except Exception:
            raise TypeError('x, y, and z should be able to initialize ' +
                            ('a {0}'.format(self.attr_classes['x'].__name__))
                             if len(set(self.attr_classes.values)) == 1 else
                            ('{0}, {1}, and {2}, resp.'.format(
                                cls.__name__ for cls in
                                self.attr_classes.values())))

        if not (x.unit.physical_type ==
                y.unit.physical_type == z.unit.physical_type):
            raise u.UnitsError("x, y, and z should have matching physical types")

        try:
            x, y, z = broadcast_arrays(x, y, z, subok=True)
        except ValueError:
            raise ValueError("Input parameters x, y, and z cannot be broadcast")

        self._x = x
        self._y = y
        self._z = z

    @property
    def x(self):
        """
        The x component of the point(s).
        """
        return self._x

    @property
    def y(self):
        """
        The y component of the point(s).
        """
        return self._y

    @property
    def z(self):
        """
        The z component of the point(s).
        """
        return self._z

    def unit_vectors(self):
        """Cartesian unit vectors in the direction of each component.

        Returns
        -------
        unit_vectors : dict of `CartesianRepresentation`
            The keys are the component names.
        """
        l = broadcast_to(1.*u.one, self.shape, subok=True)
        o = broadcast_to(0.*u.one, self.shape, subok=True)
        return OrderedDict(
            (('x', CartesianRepresentation(l, o, o, copy=False)),
             ('y', CartesianRepresentation(o, l, o, copy=False)),
             ('z', CartesianRepresentation(o, o, l, copy=False))))

    def scale_factors(self):
        """Scale factors for each component's direction.

        Returns
        -------
        scale_factors : dict of `~astropy.units.Quantity`
            The keys are the component names.
        """
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
        # Add new axis in x, y, z so one can concatenate them around it.
        # NOTE: just use np.stack once our minimum numpy version is 1.10.
        result_ndim = self.ndim + 1
        if not -result_ndim <= xyz_axis < result_ndim:
            raise IndexError('xyz_axis {0} out of bounds [-{1}, {1})'
                             .format(xyz_axis, result_ndim))
        if xyz_axis < 0:
            xyz_axis += result_ndim

        # Get x, y, z to the same units (this is very fast for identical units)
        # since np.concatenate cannot deal with quantity.
        cls = self._x.__class__
        x = self._x
        y = cls(self._y, x.unit, copy=False)
        z = cls(self._z, x.unit, copy=False)

        sh = self.shape
        sh = sh[:xyz_axis] + (1,) + sh[xyz_axis:]
        xyz_value = np.concatenate([c.reshape(sh).value for c in (x, y, z)],
                                   axis=xyz_axis)
        return cls(xyz_value, unit=x.unit, copy=False)

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
        If True arrays will be copied rather than referenced.
    """

    attr_classes = OrderedDict([('lon', Longitude),
                                ('lat', Latitude)])
    recommended_units = {'lon': u.deg, 'lat': u.deg}

    @classproperty
    def _dimensional_representation(cls):
        return SphericalRepresentation

    def __init__(self, lon, lat, copy=True):

        if not isinstance(lon, u.Quantity) or isinstance(lon, Latitude):
            raise TypeError('lon should be a Quantity, Angle, or Longitude')

        if not isinstance(lat, u.Quantity) or isinstance(lat, Longitude):
            raise TypeError('lat should be a Quantity, Angle, or Latitude')
        # Let the Longitude and Latitude classes deal with e.g. parsing
        lon = self.attr_classes['lon'](lon, copy=copy)
        lat = self.attr_classes['lat'](lat, copy=copy)

        try:
            lon, lat = broadcast_arrays(lon, lat, subok=True)
        except ValueError:
            raise ValueError("Input parameters lon and lat cannot be broadcast")

        self._lon = lon
        self._lat = lat

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
        """Cartesian unit vectors in the direction of each component.

        Returns
        -------
        unit_vectors : dict of `CartesianRepresentation`
            The keys are the component names.
        """
        sinlon, coslon = np.sin(self.lon), np.cos(self.lon)
        sinlat, coslat = np.sin(self.lat), np.cos(self.lat)
        return OrderedDict(
            (('lon', CartesianRepresentation(-sinlon, coslon, 0., copy=False)),
             ('lat', CartesianRepresentation(-sinlat*coslon, -sinlat*sinlon,
                                             coslat, copy=False))))

    def scale_factors(self):
        """Scale factors for each component's direction.

        Returns
        -------
        scale_factors : dict of `~astropy.units.Quantity`
            The keys are the component names.
        """
        coslat = np.cos(self.lat)
        l = broadcast_to(1./u.radian, self.shape, subok=True)
        return OrderedDict((('lon', coslat / u.radian), ('lat', l)))

    def to_cartesian(self):
        """
        Converts spherical polar coordinates to 3D rectangular cartesian
        coordinates.
        """
        x = np.cos(self.lat) * np.cos(self.lon)
        y = np.cos(self.lat) * np.sin(self.lon)
        z = np.sin(self.lat)

        return CartesianRepresentation(x=x, y=y, z=z, copy=False)

    @classmethod
    def from_cartesian(cls, cart):
        """
        Converts 3D rectangular cartesian coordinates to spherical polar
        coordinates.
        """

        s = np.hypot(cart.x, cart.y)

        lon = np.arctan2(cart.y, cart.x)
        lat = np.arctan2(cart.z, s)

        return cls(lon=lon, lat=lat, copy=False)

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
        return self.__class__(self.lon + 180. * u.deg, -self.lat)

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
        If True arrays will be copied rather than referenced.
    """

    attr_classes = OrderedDict([('lon', Longitude),
                                ('lat', Latitude),
                                ('distance', u.Quantity)])
    recommended_units = {'lon': u.deg, 'lat': u.deg}
    _unit_representation = UnitSphericalRepresentation

    def __init__(self, lon, lat, distance, copy=True):

        if not isinstance(lon, u.Quantity) or isinstance(lon, Latitude):
            raise TypeError('lon should be a Quantity, Angle, or Longitude')

        if not isinstance(lat, u.Quantity) or isinstance(lat, Longitude):
            raise TypeError('lat should be a Quantity, Angle, or Latitude')

        # Let the Longitude and Latitude classes deal with e.g. parsing
        lon = self.attr_classes['lon'](lon, copy=copy)
        lat = self.attr_classes['lat'](lat, copy=copy)

        distance = self.attr_classes['distance'](distance, copy=copy)
        if distance.unit.physical_type == 'length':
            distance = distance.view(Distance)

        try:
            lon, lat, distance = broadcast_arrays(lon, lat, distance,
                                                  subok=True)
        except ValueError:
            raise ValueError("Input parameters lon, lat, and distance cannot be broadcast")

        self._lon = lon
        self._lat = lat
        self._distance = distance

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
        """Cartesian unit vectors in the direction of each component.

        Returns
        -------
        unit_vectors : dict of `CartesianRepresentation`
            The keys are the component names.
        """
        sinlon, coslon = np.sin(self.lon), np.cos(self.lon)
        sinlat, coslat = np.sin(self.lat), np.cos(self.lat)
        return OrderedDict(
            (('lon', CartesianRepresentation(-sinlon, coslon, 0., copy=False)),
             ('lat', CartesianRepresentation(-sinlat*coslon, -sinlat*sinlon,
                                             coslat, copy=False)),
             ('distance', CartesianRepresentation(coslat*coslon, coslat*sinlon,
                                                  sinlat, copy=False))))

    def scale_factors(self):
        """Scale factors for each component's direction.

        Returns
        -------
        scale_factors : dict of `~astropy.units.Quantity`
            The keys are the component names.
        """
        d = self.distance / u.radian
        coslat = np.cos(self.lat)
        l = broadcast_to(1.*u.one, self.shape, subok=True)
        return OrderedDict((('lon', d * coslat),
                            ('lat', d),
                            ('distance', l)))

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

        return CartesianRepresentation(x=x, y=y, z=z, copy=False)

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

        return cls(lon=lon, lat=lat, distance=r, copy=False)

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
        If True arrays will be copied rather than referenced.
    """

    attr_classes = OrderedDict([('phi', Angle),
                                ('theta', Angle),
                                ('r', u.Quantity)])
    recommended_units = {'phi': u.deg, 'theta': u.deg}

    def __init__(self, phi, theta, r, copy=True):

        if not isinstance(phi, u.Quantity) or isinstance(phi, Latitude):
            raise TypeError('phi should be a Quantity or Angle')

        if not isinstance(theta, u.Quantity) or isinstance(theta, Longitude):
            raise TypeError('phi should be a Quantity or Angle')

        # Let the Longitude and Latitude classes deal with e.g. parsing
        phi = self.attr_classes['phi'](phi, copy=copy)
        theta = self.attr_classes['theta'](theta, copy=copy)

        # Wrap/validate phi/theta
        if copy:
            phi = phi.wrap_at(360 * u.deg)
        else:
            # necessary because the above version of `wrap_at` has to be a copy
            phi.wrap_at(360 * u.deg, inplace=True)
        if np.any(theta.value < 0.) or np.any(theta.value > 180.):
            raise ValueError('Inclination angle(s) must be within 0 deg <= angle <= 180 deg, '
                             'got {0}'.format(theta.to(u.degree)))

        r = self.attr_classes['r'](r, copy=copy)
        if r.unit.physical_type == 'length':
            r = r.view(Distance)

        try:
            phi, theta, r = broadcast_arrays(phi, theta, r, subok=True)
        except ValueError:
            raise ValueError("Input parameters phi, theta, and r cannot be broadcast")

        self._phi = phi
        self._theta = theta
        self._distance = r

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
        return self._distance

    def unit_vectors(self):
        """Cartesian unit vectors in the direction of each component.

        Returns
        -------
        unit_vectors : dict of `CartesianRepresentation`
            The keys are the component names.
        """
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
        """Scale factors for each component's direction.

        Returns
        -------
        scale_factors : dict of `~astropy.units.Quantity`
            The keys are the component names.
        """
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

    phi : `~astropy.units.Quantity`
        The azimuth of the point(s), in angular units, which will be wrapped
        to an angle between 0 and 360 degrees. This can also be instances of
        `~astropy.coordinates.Angle`,

    z : `~astropy.units.Quantity`
        The z coordinate(s) of the point(s)

    copy : bool, optional
        If True arrays will be copied rather than referenced.
    """

    attr_classes = OrderedDict([('rho', u.Quantity),
                                ('phi', Angle),
                                ('z', u.Quantity)])
    recommended_units = {'phi': u.deg}

    def __init__(self, rho, phi, z, copy=True):

        if not isinstance(phi, u.Quantity) or isinstance(phi, Latitude):
            raise TypeError('phi should be a Quantity or Angle')

        rho = self.attr_classes['rho'](rho, copy=copy)
        phi = self.attr_classes['phi'](phi, copy=copy)
        z = self.attr_classes['z'](z, copy=copy)

        if not (rho.unit.physical_type == z.unit.physical_type):
            raise u.UnitsError("rho and z should have matching physical types")

        try:
            rho, phi, z = broadcast_arrays(rho, phi, z, subok=True)
        except ValueError:
            raise ValueError("Input parameters rho, phi, and z cannot be broadcast")

        self._rho = rho
        self._phi = phi
        self._z = z

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
        """Cartesian unit vectors in the direction of each component.

        Returns
        -------
        unit_vectors : dict of `CartesianRepresentation`
            The keys are the component names.
        """
        sinphi, cosphi = np.sin(self.phi), np.cos(self.phi)
        l = broadcast_to(1., self.shape)
        return OrderedDict(
            (('rho', CartesianRepresentation(cosphi, sinphi, 0, copy=False)),
             ('phi', CartesianRepresentation(-sinphi, cosphi, 0, copy=False)),
             ('z', CartesianRepresentation(0, 0, l, unit=u.one, copy=False))))

    def scale_factors(self):
        """Scale factors for each component's direction.

        Returns
        -------
        scale_factors : dict of `~astropy.units.Quantity`
            The keys are the component names.
        """
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


class BaseOffset(RepresentationBase):
    """Offsets from points on a base representation.

    Parameters
    ----------
    d_* : `~astropy.units.Quantity`
        The offsets in terms of the components of the base representation.
        They can either be given in order in ``args`` or by name in ``kwargs``,
        where the names are the base representation names prefixed by 'd_'.
        The units should be consistent with those of the components, but not
        necessarily equivalent (e.g., for a spherical representation, they can
        be rotation rates (angular per time) for longitude and latitude, and
        velocity for distance.
    copy : bool, optional
        If True arrays will be copied rather than referenced.
    """

    base_representation = None
    attr_classes = OrderedDict()

    def __init__(self, *args, **kwargs):
        self.attr_classes = OrderedDict()
        attrs = []
        names = []
        regular_units = []
        angular_units = []
        base_attr_classes = self.base_representation.attr_classes
        args = list(args)
        for component, base_attr_cls in base_attr_classes.items():
            name = 'd_' + component
            attr = args.pop(0) if args else kwargs.pop(name)
            attrs.append(attr)
            names.append(name)
            try:
                unit = attr.unit
            except AttributeError:
                raise TypeError('attributes should have units.')
            if issubclass(base_attr_cls, Angle):
                angular_units.append(unit)
            else:
                regular_units.append(unit)

        copy = args.pop(0) if args else kwargs.pop('copy', True)
        if args:
            raise TypeError('unexpected arguments: {0}'.format(args))

        if kwargs:
            for name in names:
                if name in kwargs:
                    raise TypeError("__init__() got multiple values for "
                                    "argument {0!r}".format(name))

            raise TypeError('unexpected keyword arguments: {0}'.format(kwargs))

        for units in regular_units, angular_units:
            for unit in units[1:]:
                if not unit.is_equivalent(units[0]):
                    raise u.UnitsError('inconsistent units: {0}.'
                                       .format(units))

        attrs = [u.Quantity(a, copy=copy, subok=True) for a in attrs]
        try:
            attrs = broadcast_arrays(*attrs, subok=True)
        except ValueError:
            raise ValueError("Input parameters cannot be broadcast")

        for name, attr in zip(names, attrs):
            self.attr_classes[name] = type(attr)
            setattr(self, '_' + name, attr)

    @classmethod
    def _check_base(cls, base):
        if not isinstance(base, cls.base_representation):
            raise TypeError('need a base of the correct representation type, '
                            '{0}, not {1}.'.format(cls.base_representation,
                                                   type(base)))

    def to_cartesian(self, base):
        """Convert the offset to 3D rectangular cartesian coordinates.

        Parameters
        ----------
        base : instance of ``self.base_representation``
             The points for which the offsets are to converted: each of the
             components is multiplied by its unit vectors and scale factors.
        """
        self._check_base(base)
        base_e = base.unit_vectors()
        base_sf = base.scale_factors()
        return functools.reduce(
            operator.add, (getattr(self, 'd_' + component) *
                           base_sf[component] * e
                           for component, e in six.iteritems(base_e)))

    @classmethod
    def from_cartesian(cls, other, base):
        """Interpret 3D rectangular cartesian coordinates as offsets.

        Parameters
        ----------
        base : instance of ``self.base_representation``
            Base relative to which the offsets are defined.
        """
        cls._check_base(base)
        base_e = base.unit_vectors()
        base_sf = base.scale_factors()
        return cls(*(other.dot(e / base_sf[component])
                     for component, e in six.iteritems(base_e)))

    def _scale_operation(self, op, *args):
        scaled_attrs = [op(getattr(self, c), *args) for c in self.components]
        return self.__class__(*scaled_attrs, copy=False)

    def _combine_operation(self, op, other, reverse=False):
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

    def norm(self, base=None):
        """Vector norm.

        The norm is the standard Frobenius norm, i.e., the square root of the
        sum of the squares of all components with non-angular units.

        Parameters
        ----------
        base : instance of ``self.base_representation``
            Base relative to which the offsets are defined, which is required
            to calculate the physical size of the offset.

        Returns
        -------
        norm : `astropy.units.Quantity`
            Vector norm, with the same shape as the representation.
        """
        return self.to_cartesian(base).norm()

    # ensure components behave as if they were properties.
    def __getattr__(self, attr):
        if attr in self.attr_classes:
            return getattr(self, '_' + attr)

        return self.__getattribute__(attr)

    def __setattr__(self, attr, value):
        if attr in self.attr_classes:
            raise AttributeError("can't set attribute.")
        return super(BaseOffset, self).__setattr__(attr, value)

    def __delattr__(self, attr):
        if attr in self.attr_classes:
            raise AttributeError("can't delete attribute.")
        return super(BaseOffset, self).__delattr__(attr)


class CartesianOffset(BaseOffset):
    def to_cartesian(self, base=None):
        return CartesianRepresentation(*[getattr(self, c) for c
                                         in self.components])

    @classmethod
    def from_cartesian(cls, other, base=None):
        return cls(*[getattr(other, c) for c in other.components])


for r, r_cls in REPRESENTATION_CLASSES.items():
    name = r_cls.__name__.replace('Representation', 'Offset')
    if name not in locals():
        locals()[name] = type(name, (BaseOffset,),
                              {'base_representation': r_cls})
    if name not in __all__:
        __all__.append(name)
