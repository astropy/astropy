# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module contains the classes and utility functions for distance and
cartesian coordinates.
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from abc import ABCMeta, abstractproperty

import numpy as np

from .. import units as u

__all__ = ['Distance', 'CartesianPoints', 'cartesian_to_spherical',
           'spherical_to_cartesian']


__doctest_requires__ = {'*': ['scipy.integrate']}


class Distance(u.Quantity):
    """
    A one-dimensional distance.

    This can be initialized in one of four ways:

    * A distance `value` (array or float) and a `unit`
    * A `~astropy.units.quantity.Quantity` object
    * A redshift and (optionally) a cosmology.
    * Providing a distance modulus

    Parameters
    ----------
    value : scalar or `~astropy.units.quantity.Quantity`
        The value of this distance
    unit : `~astropy.units.core.UnitBase`
        The units for this distance, *if* `value` is not a `Quantity`.
        Must have dimensions of distance.
    z : float
        A redshift for this distance.  It will be converted to a distance
        by computing the luminosity distance for this redshift given the
        cosmology specified by `cosmology`. Must be given as a keyword argument.
    cosmology : `~astropy.cosmology.Cosmology` or None
        A cosmology that will be used to compute the distance from `z`.
        If None, the current cosmology will be used (see
        `astropy.cosmology` for details).
    distmod : float or `~astropy.units.Quantity`
        The distance modulus for this distance.
    dtype : ~numpy.dtype, optional
        See `~astropy.units.Quantity`. Must be given as a keyword argument.
    copy : bool, optional
        See `~astropy.units.Quantity`. Must be given as a keyword argument.

    Raises
    ------
    astropy.units.core.UnitsError
        If the `unit` is not a distance.
    ValueError
        If `z` is provided with a `unit` or `cosmology` is provided when `z` is
        *not* given, or `value` is given as well as `z`

    Examples
    --------
    >>> from astropy import units as u
    >>> from astropy import cosmology
    >>> from astropy.cosmology import WMAP5, WMAP7
    >>> cosmology.set_current(WMAP7)
    >>> d1 = Distance(10, u.Mpc)
    >>> d2 = Distance(40, unit=u.au)
    >>> d3 = Distance(value=5, unit=u.kpc)
    >>> d4 = Distance(z=0.23)
    >>> d5 = Distance(z=0.23, cosmology=WMAP5)
    >>> d6 = Distance(distmod=24.47)
    """

    _include_easy_conversion_members = True

    def __new__(cls, value=None, unit=None, z=None, cosmology=None,
                distmod=None, dtype=None, copy=True):
        from ..cosmology import get_current

        if isinstance(value, u.Quantity):
            # This includes Distances as well
            if z is not None or distmod is not None:
                raise ValueError('`value` was given along with `z` or `distmod`'
                                 ' in Quantity constructor.')

            if unit is not None:
                value = value.to(unit).value
            else:
                unit = value.unit
                value = value.value
        elif value is None:
            if z is not None:
                if distmod is not None:
                    raise ValueError('both `z` and `distmod` given in Distance '
                                     'constructor')

                if cosmology is None:
                    cosmology = get_current()

                ld = cosmology.luminosity_distance(z)

                if unit is not None:
                    ld = ld.to(unit)
                value = ld.value
                unit = ld.unit

            elif distmod is not None:
                value = cls._distmod_to_pc(distmod)
                if unit is None:
                    # choose unit based on most reasonable of Mpc, kpc, or pc
                    if value > 1e6:
                        value = value / 1e6
                        unit = u.megaparsec
                    elif value > 1e3:
                        value = value / 1e3
                        unit = u.kiloparsec
                    else:
                        unit = u.parsec
                else:
                    value = u.Quantity(value, u.parsec).to(unit).value
            else:
                raise ValueError('none of `value`, `z`, or `distmod` were given'
                                 ' to Distance constructor')

                value = ld.value
                unit = ld.unit
        elif z is not None:  # and value is not None based on above
            raise ValueError('Both `z` and a `value` were provided in Distance '
                             'constructor')
        elif cosmology is not None:
            raise ValueError('A `cosmology` was given but `z` was not provided '
                             'in Distance constructor')
        elif unit is None:
            raise u.UnitsError('No unit was provided to Distance constructor')
        #"else" the baseline `value` + `unit` case

        unit = cls._convert_to_and_validate_distance_unit(unit)

        try:
            value = np.asarray(value)
        except ValueError as e:
            raise TypeError(str(e))

        if value.dtype.kind not in 'iuf':
            raise TypeError("Unsupported dtype '{0}'".format(value.dtype))

        return super(Distance, cls).__new__(cls, value, unit, dtype=dtype,
                                            copy=copy)

    def __quantity_view__(self, obj, unit):
        unit = self._convert_to_and_validate_distance_unit(unit)
        return super(Distance, self).__quantity_view__(obj, unit)

    def __quantity_instance__(self, val, unit, **kwargs):
        unit = self._convert_to_and_validate_distance_unit(unit)
        return super(Distance, self).__quantity_instance__(val, unit, **kwargs)

    @staticmethod
    def _convert_to_and_validate_distance_unit(unit):
        """
        raises astropy.units.UnitsError if not a distance unit
        """
        if unit is not None:
            unit = u.Unit(unit)

        if not unit.is_equivalent(u.kpc):
            raise u.UnitsError(six.u('Unit "{0}" is not a distance').format(unit))
        return unit

    @property
    def z(self):
        """Short for ``self.compute_z()``"""
        return self.compute_z()

    def compute_z(self, cosmology=None):
        """
        The redshift for this distance assuming its physical distance is
        a luminosity distance.

        Parameters
        ----------
        cosmology : `~astropy.cosmology.cosmology` or None
            The cosmology to assume for this calculation, or None to use the
            current cosmology.

        Returns
        -------
        z : float
            The redshift of this distance given the provided `cosmology`.
        """
        from ..cosmology import luminosity_distance
        from scipy import optimize

        # FIXME: array: need to make this calculation more vector-friendly

        f = lambda z, d, cos: (luminosity_distance(z, cos).value - d) ** 2
        return optimize.brent(f, (self.Mpc, cosmology))

    @property
    def distmod(self):
        """  The distance modulus of this distance as a Quantity """
        val = 5. * np.log10(self.to(u.pc).value) - 5.
        return u.Quantity(val, u.mag)

    @staticmethod
    def _distmod_to_pc(dm):
        return 10 ** ((dm + 5) / 5.)


class CartesianPoints(object):
    """
    A cartesian representation of a point in three-dimensional space.

    Attributes
    ----------
    x : number or array
        The first cartesian coordinate.
    y : number or array
        The second cartesian coordinate.
    z : number or array
        The third cartesian coordinate.
    unit : `~astropy.units.UnitBase` object or None
        The physical unit of the coordinate values.
    """

    def __init__(self, x, y, z, unit=None):
        self.x = x
        self.y = y
        self.z = z
        self.unit = unit

    def to_spherical(self):
        """
        Converts to the spherical representation of this point.

        Returns
        -------
        r : float or array
            The radial coordinate (in the same units as the inputs).
        lat : float or array
            The latitude in radians
        lon : float or array
            The longitude in radians

        """
        return cartesian_to_spherical(self.x, self.y, self.z)

    def __repr__(self):
        return '<CartesianPoints ({x}, {y}, {z}) {unit}>'.format(x=self.x,
                y=self.y, z=self.z, unit=self.unit)

    def __eq__(self, other):
        return (isinstance(other, CartesianPoints) and self.x == other.x and
                self.y == other.y and self.z == other.z and
                self.unit == other.unit)

    def __add__(self, other):
        if isinstance(other, CartesianPoints) or (hasattr(other, 'x') and
            hasattr(other, 'y') and hasattr(other, 'z') and
            hasattr(other, 'unit')):
            newx = self.x + other.unit.to(self.unit, other.x)
            newy = self.y + other.unit.to(self.unit, other.y)
            newz = self.z + other.unit.to(self.unit, other.z)
        else:
            msg = "unsupported operand type(s) for +: '{sel}' and '{other}'"
            raise TypeError(msg.format(type(self).__name__,
                                        type(other).__name__))
        return CartesianPoints(newx, newy, newz, self.unit)

    def __sub__(self, other):
        if isinstance(other, CartesianPoints) or (hasattr(other, 'x') and
            hasattr(other, 'y') and hasattr(other, 'z') and
            hasattr(other, 'unit')):
            newx = self.x - other.unit.to(self.unit, other.x)
            newy = self.y - other.unit.to(self.unit, other.y)
            newz = self.z - other.unit.to(self.unit, other.z)
        else:
            msg = "unsupported operand type(s) for -: '{sel}' and '{other}'"
            raise TypeError(msg.format(type(self).__name__,
                                        type(other).__name__))
        return CartesianPoints(newx, newy, newz, self.unit)

#<------------transformation-related utility functions----------------->


def cartesian_to_spherical(x, y, z):
    """
    Converts 3D rectangular cartesian coordinates to spherical polar
    coordinates.

    Note that the resulting angles are latitude/longitude or
    elevation/azimuthal form.  I.e., the origin is along the equator
    rather than at the north pole.

    .. note::
        This is a low-level function used internally in
        `astropy.coordinates`.  It is provided for users if they really
        want to use it, but it is recommended that you use the
        `astropy.coordinates` coordinate systems.

    Parameters
    ----------
    x : scalar or array-like
        The first cartesian coordinate.
    y : scalar or array-like
        The second cartesian coordinate.
    z : scalar or array-like
        The third cartesian coordinate.

    Returns
    -------
    r : float or array
        The radial coordinate (in the same units as the inputs).
    lat : float or array
        The latitude in radians
    lon : float or array
        The longitude in radians
    """
    import math

    xsq = x ** 2
    ysq = y ** 2
    zsq = z ** 2

    r = (xsq + ysq + zsq) ** 0.5
    s = (xsq + ysq) ** 0.5

    if np.isscalar(x) and np.isscalar(y) and np.isscalar(z):
        lon = math.atan2(y, x)
        lat = math.atan2(z, s)
    else:
        lon = np.arctan2(y, x)
        lat = np.arctan2(z, s)

    return r, lat, lon


def spherical_to_cartesian(r, lat, lon):
    """
    Converts spherical polar coordinates to rectangular cartesian
    coordinates.

    Note that the input angles should be in latitude/longitude or
    elevation/azimuthal form.  I.e., the origin is along the equator
    rather than at the north pole.

    .. note::
        This is a low-level function used internally in
        `astropy.coordinates`.  It is provided for users if they really
        want to use it, but it is recommended that you use the
        `astropy.coordinates` coordinate systems.

    Parameters
    ----------
    r : scalar or array-like
        The radial coordinate (in the same units as the inputs).
    lat : scalar or array-like
        The latitude in radians
    lon : scalar or array-like
        The longitude in radians

    Returns
    -------
    x : float or array
        The first cartesian coordinate.
    y : float or array
        The second cartesian coordinate.
    z : float or array
        The third cartesian coordinate.


    """
    import math

    if np.isscalar(r) and np.isscalar(lat) and np.isscalar(lon):
        x = r * math.cos(lat) * math.cos(lon)
        y = r * math.cos(lat) * math.sin(lon)
        z = r * math.sin(lat)
    else:
        x = r * np.cos(lat) * np.cos(lon)
        y = r * np.cos(lat) * np.sin(lon)
        z = r * np.sin(lat)

    return x, y, z
