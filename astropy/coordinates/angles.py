# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module contains the fundamental classes used for representing
coordinates in astropy.
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import math

import numpy as np

from ..extern import six
from . import angle_utilities as util
from .. import units as u
from ..utils import deprecated


__all__ = ['Angle', 'Latitude', 'Longitude']


TWOPI = math.pi * 2.0  # no need to calculate this all the time


class Angle(u.Quantity):
    """
    One or more angular value(s) with units equivalent to radians or degrees.

    An angle can be specified either as an array, scalar, tuple (see
    below), string, `~astropy.units.Quantity` or another
    `~astropy.coordinates.Angle`.

    The input parser is flexible and supports a variety of formats::

      Angle('10.2345d')
      Angle(['10.2345d', '-20d'])
      Angle('1:2:30.43 degrees')
      Angle('1 2 0 hours')
      Angle(np.arange(1, 8), unit=u.deg)
      Angle(u'1°2′3″')
      Angle('1d2m3.4s')
      Angle('-1h2m3s')
      Angle((-1, 2, 3), unit=u.deg)  # (d, m, s)
      Angle(10.2345 * u.deg)
      Angle(Angle(10.2345 * u.deg))

    Parameters
    ----------
    angle : array, scalar, Quantity, Angle
        The angle value. If a tuple, will be interpreted as ``(h, m
        s)`` or ``(d, m, s)`` depending on `unit`. If a string, it
        will be interpreted following the rules described above.

        If `angle` is a sequence or array of strings, the resulting
        values will be in the given `unit`, or if None is provided,
        the unit will be taken from the first given value.

    unit : `~astropy.units.UnitBase`, str, optional
        The unit of the value specified for the angle.  This may be
        any string that `~astropy.units.Unit` understands, but it is
        better to give an actual unit object.  Must be an angular
        unit.

    dtype : ~numpy.dtype, optional
        See `~astropy.units.Quantity`.

    copy : bool, optional
        See `~astropy.units.Quantity`.

    Raises
    ------
    `~astropy.units.core.UnitsError`
        If a unit is not provided or it is not an angular unit.
    """
    _include_easy_conversion_members = True

    def __new__(cls, angle, unit=None, dtype=None, copy=True):
        unit = cls._convert_unit_to_angle_unit(unit)
        if (unit is not None and
            not unit.is_equivalent(u.radian)):
            raise u.UnitsError(
                "Given unit {0} is not convertible to an angle".format(
                    unit))

        if isinstance(angle, u.Quantity):
            # This includes Angle subclasses as well
            if unit is not None:
                angle = angle.to(unit).value
            else:
                unit = angle.unit
                unit = cls._convert_unit_to_angle_unit(unit)
                if not unit.is_equivalent(u.radian):
                    raise u.UnitsError(
                        "Given quantity {0} is not convertible to an angle".format(
                            angle))

                angle = angle.value

        angle = cls._tuple_to_float(angle, unit)

        try:
            angle = np.asarray(angle)
        except ValueError as e:
            raise TypeError(str(e))

        if angle.dtype.type in (np.string_, np.unicode_):
            # We need to modify this value from within
            # convert_string_to_angle, and the only way to do that
            # across Python 2.6 - 3.3 is to use this "store it in a
            # list" trick.
            determined_unit = [unit]

            def convert_string_to_angle(x):
                ang, new_unit = util.parse_angle(six.text_type(x), unit)
                if determined_unit[0] is None:
                    determined_unit[0] = new_unit
                if new_unit is not None:
                    return new_unit.to(
                        determined_unit[0], cls._tuple_to_float(ang, new_unit))
                else:
                    return cls._tuple_to_float(ang, determined_unit[0])

            convert_string_to_angle_ufunc = np.vectorize(
                convert_string_to_angle,
                otypes=[np.float_])
            angle = convert_string_to_angle_ufunc(angle)
            unit = determined_unit[0]

        elif angle.dtype.kind not in 'iuf':
            raise TypeError("Unsupported dtype '{0}'".format(angle.dtype))

        if unit is None:
            raise u.UnitsError("No unit was specified")

        self = super(Angle, cls).__new__(
            cls, angle, unit, dtype=dtype, copy=copy)

        return self

    @staticmethod
    def _tuple_to_float(angle, unit):
        """
        Converts an angle represented as a 3-tuple into a floating
        point number in the given unit.
        """
        if isinstance(angle, tuple):
            # TODO: Numpy array of tuples?
            if unit is u.hourangle:
                util.check_hms_ranges(*angle)
                angle = util.hms_to_hours(*angle)
            elif unit is u.degree:
                angle = util.dms_to_degrees(*angle)
            else:
                raise u.UnitsError(
                    "Can not parse '{0}' as unit '{1}'".format(
                        angle, unit))
        return angle

    @staticmethod
    def _convert_unit_to_angle_unit(unit):
        if unit is not None:
            unit = u.Unit(unit)

            if unit is u.hour:
                unit = u.hourangle
        return unit

    def __quantity_view__(self, obj, unit):
        unit = self._convert_unit_to_angle_unit(unit)
        if unit is not None and unit.is_equivalent(u.radian):
            result = obj.view(Angle)
            return result
        return super(Angle, self).__quantity_view__(
            obj, unit)

    def __quantity_instance__(self, val, unit, **kwargs):
        unit = self._convert_unit_to_angle_unit(unit)
        if unit is not None and unit.is_equivalent(u.radian):
            return self.__class__(val, unit, **kwargs)
        return super(Angle, self).__quantity_instance__(val, unit, **kwargs)

    def __array_wrap__(self, obj, context=None):
        obj = super(Angle, self).__array_wrap__(obj, context=context)

        if isinstance(obj, Angle):
            return Angle(obj.value, obj.unit)

        return obj

    def __add__(self, other):
        return super(Angle, self).__add__(other)

    def __sub__(self, other):
        return super(Angle, self).__sub__(other)

    def __mul__(self, other):
        if isinstance(other, type(self)):
            raise TypeError(
                "multiplication is not supported between two {0} "
                "objects".format(
                    type(self).__name__))
        return super(Angle, self).__mul__(other)

    def __div__(self, other):
        if isinstance(other, type(self)):
            raise TypeError(
                "division is not supported between two {0} objects".format(
                    type(self).__name__))
        return super(Angle, self).__div__(other)

    __truediv__ = __div__

    @property
    def hour(self):
        """
        The angle's value in hours (read-only property).
        """
        return self.hourangle

    @property
    def hms(self):
        """
        The angle's value in hours, as a ``(h, m, s)`` tuple
        (read-only property).
        """
        return util.hours_to_hms(self.hourangle)

    @property
    def dms(self):
        """
        The angle's value in degrees, as a ``(d, m, s)`` tuple
        (read-only property).
        """
        return util.degrees_to_dms(self.degree)

    def to_string(self, unit=None, decimal=False, sep='fromunit',
                  precision=None, alwayssign=False, pad=False,
                  fields=3, format=None):
        """ A string representation of the angle.

        Parameters
        ----------
        units : `~astropy.units.UnitBase`, optional
            Specifies the units.  Must be an angular unit.  If not
            provided, the unit used to initialize the angle will be
            used.

        decimal : bool, optional
            If True, a decimal respresentation will be used, otherwise
            the returned string will be in sexagesimal form.

        sep : str, optional
            The separator between numbers in a sexagesimal
            representation.  E.g., if it is ':', the result is
            "12:41:11.1241". Also accepts 2 or 3 separators. E.g.,
            ``sep='hms'`` would give the result "12h41m11.1241s", or
            sep='-:' would yield "11-21:17.124".  Alternatively, the
            special string 'fromunit' means 'dms' if the unit is
            degrees, or 'hms' if the unit is hours.

        precision : int, optional
            The level of decimal precision.  If `decimal` is True,
            this is the raw precision, otherwise it gives the
            precision of the last place of the sexagesimal
            representation (seconds).  If `None`, or not provided, the
            number of decimal places is determined by the value, and
            will be between 0-8 decimal places as required.

        alwayssign : bool, optional
            If `True`, include the sign no matter what.  If `False`,
            only include the sign if it is negative.

        pad : bool, optional
            If `True`, include leading zeros when needed to ensure a
            fixed number of characters for sexagesimal representation.

        fields : int, optional
            Specifies the number of fields to display when outputting
            sexagesimal notation.  For example:

                - fields == 1: `'5d'`
                - fields == 2: `'5d45m'`
                - fields == 3: `'5d45m32.5s'`

            By default, all fields are displayed.

        format : str, optional
            The format of the result.  If not provided, an unadorned
            string is returned.  Supported values are:

            - 'latex': Return a LaTeX-formatted string

            - 'unicode': Return a string containing non-ASCII unicode
              characters, such as the degree symbol

        Returns
        -------
        strrepr : str
            A string representation of the angle.

        """
        if unit is None:
            unit = self.unit
        unit = self._convert_unit_to_angle_unit(unit)

        separators = {
            None: {
                u.degree: 'dms',
                u.hourangle: 'hms'},
            'latex': {
                u.degree: [r'^\circ', r'{}^\prime', r'{}^{\prime\prime}'],
                u.hourangle: [r'^\mathrm{h}', r'^\mathrm{m}', r'^\mathrm{s}']},
            'unicode': {
                u.degree: '°′″',
                u.hourangle: 'ʰᵐˢ'}
            }

        if sep == 'fromunit':
            if format not in separators:
                raise ValueError("Unknown format '{0}'".format(format))
            seps = separators[format]
            if unit in seps:
                sep = seps[unit]

        # Create an iterator so we can format each element of what
        # might be an array.
        if unit is u.degree:
            if decimal:
                values = self.degree
                if precision is not None:
                    func = ("{0:0." + str(precision) + "f}").format
                else:
                    func = '{0:g}'.format
            else:
                if sep == 'fromunit':
                    sep = 'dms'
                values = self.degree
                func = lambda x: util.degrees_to_string(
                    x, precision=precision, sep=sep, pad=pad,
                    fields=fields)

        elif unit is u.hourangle:
            if decimal:
                values = self.hour
                if precision is not None:
                    func = ("{0:0." + str(precision) + "f}").format
                else:
                    func = '{0:g}'.format
            else:
                if sep == 'fromunit':
                    sep = 'hms'
                values = self.hour
                func = lambda x: util.hours_to_string(
                    x, precision=precision, sep=sep, pad=pad,
                    fields=fields)

        elif unit.is_equivalent(u.radian):
            if decimal:
                values = self.to(unit).value
                if precision is not None:
                    func = ("{0:1." + str(precision) + "f}").format
                else:
                    func = "{0:g}".format
            elif sep == 'fromunit':
                values = self.to(unit).value
                unit_string = unit.to_string(format=format)
                if format == 'latex':
                    unit_string = unit_string[1:-1]

                if precision is not None:
                    def plain_unit_format(val):
                        return ("{0:0." + str(precision) + "f}{1}").format(
                            val, unit_string)
                    func = plain_unit_format
                else:
                    def plain_unit_format(val):
                        return "{0:g}{1}".format(val, unit_string)
                    func = plain_unit_format
            else:
                raise ValueError(
                    "'{0}' can not be represented in sexagesimal "
                    "notation".format(
                        unit.name))

        else:
            raise u.UnitsError(
                "The unit value provided is not an angular unit.")

        def do_format(val):
            s = func(float(val))
            if alwayssign and not s.startswith('-'):
                s = '+' + s
            if format == 'latex':
                s = '${0}$'.format(s)
            return s

        format_ufunc = np.vectorize(do_format, otypes=[np.object])
        result = format_ufunc(values)
        if result.ndim == 0:
            result = result[()]
        return result

    def wrap_at(self, wrap_angle, inplace=False):
        """
        Wrap the Angle object at the given ``wrap_angle``.

        This method forces all the angle values to be within a contiguous 360 degree
        range so that ``wrap_angle - 360d <= angle < wrap_angle``.  By default a new
        Angle object is returned, but if the ``inplace`` argument is ``True`` then
        the Angle object is wrapped in place and nothing is returned.

        For instance::

          >>> from astropy.coordinates import Angle
          >>> import astropy.units as u
          >>> a = Angle([-20.0, 150.0, 350.0] * u.deg)

          >>> a.wrap_at(360 * u.deg).degree  # Wrap into range 0 to 360 degrees
          array([ 340.,  150.,  350.])

          >>> a.wrap_at('180d', inplace=True)  # Wrap into range -180 to 180 degrees
          >>> a.degree
          array([ -20.,  150.,  -10.])

        Parameters
        ----------
        wrap_angle : str, Angle, angular Quantity
            Specifies a single value for the wrap angle.  This can be any
            object that can initialize an Angle object, e.g. '180d', 180 * u.deg,
            or Angle(180, unit=u.deg).

        inplace : bool
            If ``True`` then wrap the object in place instead of returning a new Angle

        Returns
        -------
        out : Angle or None
            If ``inplace is False`` (default), return new Angle object with angles
            wrapped accordingly.  Otherwise wrap in place and return None.
        """
        wrap_angle = Angle(wrap_angle)  # Convert to an Angle
        wrapped = np.mod(self - wrap_angle, 360.0 * u.deg) - (360.0 * u.deg - wrap_angle)

        if inplace:
            self[()] = wrapped
        else:
            return wrapped

    def is_within_bounds(self, lower=None, upper=None):
        """
        Check if all angle(s) satisfy ``lower <= angle < upper``

        If ``lower`` is not specified (or ``None``) then no lower bounds check is
        performed.  Likewise ``upper`` can be left unspecified.  For example::

          >>> from astropy.coordinates import Angle
          >>> import astropy.units as u
          >>> a = Angle([-20, 150, 350] * u.deg)
          >>> a.is_within_bounds('0d', '360d')
          False
          >>> a.is_within_bounds(None, '360d')
          True
          >>> a.is_within_bounds(-30 * u.deg, None)
          True

        Parameters
        ----------
        lower : str, Angle, angular Quantity, None
            Specifies lower bound for checking.  This can be any object
            that can initialize an Angle object, e.g. '180d', 180 * u.deg,
            or Angle(180, unit=u.deg).
        upper : str, Angle, angular Quantity, None
            Specifies upper bound for checking.  This can be any object
            that can initialize an Angle object, e.g. '180d', 180 * u.deg,
            or Angle(180, unit=u.deg).

        Returns
        -------
        is_within_bounds : bool
            True if all angles satisfy ``lower <= angle < upper``
        """
        ok = True
        if lower is not None:
            ok &= np.all(Angle(lower) <= self)
        if ok and upper is not None:
            ok &= np.all(self < Angle(upper))
        return bool(ok)

    @deprecated("0.3", name="format", alternative="to_string")
    def format(self, unit=u.degree, decimal=False, sep='fromunit', precision=5,
               alwayssign=False, pad=False):
        return self.to_string(
            unit=unit, decimal=decimal, sep=sep, precision=precision,
            alwayssign=alwayssign, pad=pad)

    def __str__(self):
        return str(self.to_string())

    def _repr_latex_(self):
        return str(self.to_string(format='latex'))


class Latitude(Angle):
    """
    Latitude-like angle(s) which must be in the range -90 to +90 deg.

    A Latitude object is distinguished from a pure `Angle` by virtue
    of being constrained so that::

      -90.0 * u.deg <= angle(s) <= +90.0 * u.deg

    Any attempt to set a value outside that range will result in a `ValueError`.

    The input angle(s) can be specified either as an array, list, scalar, tuple (see
    below), string, :class:`~astropy.units.quantity.Quantity` or another `Angle`.

    The input parser is flexible and supports all of the input formats supported by `Angle`.

    Parameters
    ----------
    angle : array, list, scalar, Quantity, Angle
        The angle value(s). If a tuple, will be interpreted as ``(h, m
        s)`` or ``(d, m, s)`` depending on `unit`. If a string, it
        will be interpreted following the rules described for `Angle`.

        If `angle` is a sequence or array of strings, the resulting
        values will be in the given `unit`, or if None is provided,
        the unit will be taken from the first given value.

    unit : :class:`~astropy.units.core.UnitBase`, str, optional
        The unit of the value specified for the angle.  This may be
        any string that `~astropy.units.Unit` understands, but it is
        better to give an actual unit object.  Must be an angular
        unit.

    Raises
    ------
    `~astropy.units.core.UnitsError`
        If a unit is not provided or it is not an angular unit.
    """
    def __new__(cls, angle, unit=None, **kwargs):
        self = super(Latitude, cls).__new__(cls, angle, unit=unit, **kwargs)
        self._validate_angles()
        return self

    def _validate_angles(self):
        if np.any(self < -90.0 * u.deg) or np.any(self > 90.0 * u.deg):
            raise ValueError('Latitude angle(s) must be within -90 deg <= angle <= 90 deg, '
                             'got {0}'.format(self.degree))

    def __setitem__(self, item, value):
        super(Latitude, self).__setitem__(item, value)
        self._validate_angles()


class Longitude(Angle):
    """
    Longitude-like angle(s) which are wrapped within a contiguous 360 degree range.

    A ``Longitude`` object is distinguished from a pure `~astropy.coordinates.Angle` by virtue
    of a ``wrap_angle`` property.  The ``wrap_angle`` specifies that all angle values
    represented by the object will be in the range::

      wrap_angle - 360 * u.deg <= angle(s) < wrap_angle

    The default ``wrap_angle`` is 360 deg.  Setting ``wrap_angle=180 * u.deg`` would
    instead result in values between -180 and +180 deg.  Setting the ``wrap_angle``
    attribute of an existing ``Longitude`` object will result in re-wrapping the
    angle values in-place.

    The input angle(s) can be specified either as an array, list, scalar, tuple,
    string, :class:`~astropy.units.quantity.Quantity` or another `Angle`.

    The input parser is flexible and supports all of the input formats supported by `Angle`.

    Parameters
    ----------
    angle : array, list, scalar, Quantity, `Angle`
        The angle value(s). If a tuple, will be interpreted as ``(h, m
        s)`` or ``(d, m, s)`` depending on `unit`. If a string, it
        will be interpreted following the rules described for `Angle`.

        If `angle` is a sequence or array of strings, the resulting
        values will be in the given `unit`, or if None is provided,
        the unit will be taken from the first given value.

    unit : :class:`~astropy.units.core.UnitBase`, str, optional
        The unit of the value specified for the angle.  This may be
        any string that `~astropy.units.Unit` understands, but it is
        better to give an actual unit object.  Must be an angular
        unit.

    wrap_angle : `Angle` or equivalent
        Angle at which to wrap back to ``wrap_angle - 360 deg``.

    Raises
    ------
    `~astropy.units.core.UnitsError`
        If a unit is not provided or it is not an angular unit.
    """
    def __new__(cls, angle, unit=None, wrap_angle=360 * u.deg, **kwargs):
        self = super(Longitude, cls).__new__(cls, angle, unit=unit, **kwargs)
        self.wrap_angle = wrap_angle
        return self

    def __setitem__(self, item, value):
        super(Longitude, self).__setitem__(item, value)
        self._wrap_internal()

    def _wrap_internal(self):
        """
        Wrap the internal values in the Longitude object.  Using the `Angle`
        wrap_at() method causes recursion.
        """
        d360 = 360.0 * u.deg
        wrapped = np.mod(self - self.wrap_angle, d360) - (d360 - self.wrap_angle)
        super(Longitude, self).__setitem__((), wrapped)

    @property
    def wrap_angle(self):
        return self._wrap_angle

    @wrap_angle.setter
    def wrap_angle(self, value):
        self._wrap_angle = Angle(value)
        self._wrap_internal()

    def __quantity_view__(self, obj, unit):
        unit = self._convert_unit_to_angle_unit(unit)
        if unit is not None and unit.is_equivalent(u.radian):
            # by default, wrap_angle and equivalencies remain the same
            # TODO: generalize to some _things_to_copy once #1422, #1373 merged
            return obj.view(Longitude)
        return super(Angle, self).__quantity_view__(obj, unit)

    def __quantity_instance__(self, val, unit, **kwargs):
        unit = self._convert_unit_to_angle_unit(unit)
        if unit is not None and unit.is_equivalent(u.radian):
            # by default, wrap_angle remains the same
            if 'wrap_angle' not in kwargs:
                kwargs['wrap_angle'] = getattr(self, 'wrap_angle')
            return Longitude(val, unit, **kwargs)
        return super(Angle, self).__quantity_instance__(val, unit, **kwargs)

    def __getitem__(self, key):
        out = super(Longitude, self).__getitem__(key)
        out._wrap_angle = self._wrap_angle
        return out

    # deprecated; TODO: move to quantity later (once #1422, #1373 merged)
    def __getslice__(self, i, j):
        return self.__getitem__(slice(i, j))

#<----------------------------------Rotations--------------------------------->


def rotation_matrix(angle, axis='z', unit=None):
    """
    Generate a 3x3 cartesian rotation matrix in for rotation about
    a particular axis.

    Parameters
    ----------
    angle : convertible to Angle
        The amount of rotation this matrix should represent.

    axis : str or 3-sequence
        Either 'x','y', 'z', or a (x,y,z) specifying an axis to rotate
        about. If 'x','y', or 'z', the rotation sense is
        counterclockwise looking down the + axis (e.g. positive
        rotations obey left-hand-rule).

    unit : UnitBase, optional
        If `angle` does not have associated units, they are in this
        unit.  If neither are provided, it is assumed to be degrees.

    Returns
    -------
    rmat: `numpy.matrix`
        A unitary rotation matrix.
    """
    # TODO: This doesn't handle arrays of angles

    if unit is None:
        unit = u.degree

    angle = Angle(angle, unit=unit)

    if axis == 'z':
        s = np.sin(angle)
        c = np.cos(angle)
        return np.matrix(((c, s, 0),
                          (-s, c, 0),
                          (0, 0, 1)))
    elif axis == 'y':
        s = np.sin(angle)
        c = np.cos(angle)
        return np.matrix(((c, 0, -s),
                          (0, 1, 0),
                          (s, 0, c)))
    elif axis == 'x':
        s = np.sin(angle)
        c = np.cos(angle)
        return np.matrix(((1, 0, 0),
                          (0, c, s),
                          (0, -s, c)))
    else:
        x, y, z = axis
        w = np.cos(angle / 2)

        # normalize
        if w == 1:
            x = y = z = 0
        else:
            l = np.sqrt((x * x + y * y + z * z) / (1 - w * w))
            x /= l
            y /= l
            z /= l

        wsq = w * w
        xsq = x * x
        ysq = y * y
        zsq = z * z
        return np.matrix(((wsq + xsq - ysq - zsq, 2 * x * y - 2 * w * z, 2 * x * z + 2 * w * y),
                          (2 * x * y + 2 * w * z, wsq - xsq + ysq - zsq, 2 * y * z - 2 * w * x),
                          (2 * x * z - 2 * w * y, 2 * y * z + 2 * w * x, wsq - xsq - ysq + zsq)))


def angle_axis(matrix, unit=None):
    """
    Computes the angle of rotation and the rotation axis for a given rotation
    matrix.

    Parameters
    ----------
    matrix : array-like
        A 3 x 3 unitary rotation matrix.

    unit : UnitBase
        The output unit.  If `None`, the output unit is degrees.

    Returns
    -------
    angle : Angle
        The angle of rotation for this matrix.

    axis : array (length 3)
        The axis of rotation for this matrix.
    """
    # TODO: This doesn't handle arrays of angles

    m = np.asmatrix(matrix)
    if m.shape != (3, 3):
        raise ValueError('matrix is not 3x3')

    angle = np.acos((m[0, 0] + m[1, 1] + m[2, 2] - 1) / 2)
    denom = np.sqrt(2 * ((m[2, 1] - m[1, 2]) + (m[0, 2] - m[2, 0]) + (m[1, 0] - m[0, 1])))
    axis = np.array((m[2, 1] - m[1, 2], m[0, 2] - m[2, 0], m[1, 0] - m[0, 1])) / denom
    axis /= np.sqrt(np.sum(axis ** 2))

    angle = Angle(angle, u.radian)
    if unit is None:
        unit = u.degree
    return angle.to(unit), axis
