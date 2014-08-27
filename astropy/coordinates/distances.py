# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module contains the classes and utility functions for distance and
cartesian coordinates.
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import math

import numpy as np

from .. import units as u
from ..utils import deprecated

__all__ = ['Distance', 'CartesianPoints', 'cartesian_to_spherical',
           'spherical_to_cartesian']


__doctest_requires__ = {'*': ['scipy.integrate']}


class Distance(u.Quantity):
    """
    A one-dimensional distance.

    This can be initialized in one of four ways:

    * A distance ``value`` (array or float) and a ``unit``
    * A `~astropy.units.Quantity` object
    * A redshift and (optionally) a cosmology.
    * Providing a distance modulus

    Parameters
    ----------
    value : scalar or `~astropy.units.Quantity`.
        The value of this distance.
    unit : `~astropy.units.UnitBase`
        The units for this distance, *if* ``value`` is not a
        `~astropy.units.Quantity`. Must have dimensions of distance.
    z : float
        A redshift for this distance.  It will be converted to a distance
        by computing the luminosity distance for this redshift given the
        cosmology specified by ``cosmology``. Must be given as a keyword
        argument.
    cosmology : ``Cosmology`` or `None`
        A cosmology that will be used to compute the distance from ``z``.
        If `None`, the current cosmology will be used (see
        `astropy.cosmology` for details).
    distmod : float or `~astropy.units.Quantity`
        The distance modulus for this distance. Note that if ``unit`` is not
        provided, a guess will be made at the unit between AU, pc, kpc, and Mpc.
    dtype : `~numpy.dtype`, optional
        See `~astropy.units.Quantity`.
    copy : bool, optional
        See `~astropy.units.Quantity`.
    order : {'C', 'F', 'A'}, optional
        See `~astropy.units.Quantity`.
    subok : bool, optional
        See `~astropy.units.Quantity`.
    ndmin : int, optional
        See `~astropy.units.Quantity`.
    allow_negative : bool, optional
        Whether to allow negative distances (which are possible is some
        cosmologies).  Default: ``False``.

    Raises
    ------
    `~astropy.units.UnitsError`
        If the ``unit`` is not a distance.
    ValueError
        If value specified is less than 0 and ``allow_negative=False``.

        If ``z`` is provided with a ``unit`` or ``cosmology`` is provided
        when ``z`` is *not* given, or ``value`` is given as well as ``z``.


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
    >>> d7 = Distance(Distance(10 * u.Mpc))
    """

    _include_easy_conversion_members = True

    def __new__(cls, value=None, unit=None, z=None, cosmology=None,
                distmod=None, dtype=None, copy=True, order=None,
                subok=False, ndmin=0, allow_negative=False):

        if z is not None:
            if value is not None or distmod is not None:
                raise ValueError('Should given only one of `value`, `z` '
                                 'or `distmod` in Distance constructor.')

            if cosmology is None:
                from ..cosmology import default_cosmology
                cosmology = default_cosmology.get()

            value = cosmology.luminosity_distance(z)
            # Continue on to take account of unit and other arguments
            # but a copy is already made, so no longer necessary
            copy = False

        else:
            if cosmology is not None:
                raise ValueError('A `cosmology` was given but `z` was not '
                                 'provided in Distance constructor')

            if distmod is not None:
                if value is not None:
                    raise ValueError('Should given only one of `value`, `z` '
                                     'or `distmod` in Distance constructor.')

                value = cls._distmod_to_pc(distmod)
                if unit is None:
                    # if the unit is not specified, guess based on the mean of
                    # the log of the distance
                    meanlogval = np.log10(value.value).mean()
                    if meanlogval > 6:
                        unit = u.Mpc
                    elif meanlogval > 3:
                        unit = u.kpc
                    elif meanlogval < -3: #~200 AU
                        unit = u.AU
                    else:
                        unit = u.pc

                # Continue on to take account of unit and other arguments
                # but a copy is already made, so no longer necessary
                copy = False

            elif value is None:
                raise ValueError('None of `value`, `z`, or `distmod` were '
                                 'given to Distance constructor')

        # now we have arguments like for a Quantity, so let it do the work
        distance = super(Distance, cls).__new__(
            cls, value, unit, dtype=dtype, copy=copy, order=order,
            subok=subok, ndmin=ndmin)

        if not distance.unit.is_equivalent(u.m):
            raise u.UnitsError('Unit "{0}" is not a length type'.format(unit))

        if not allow_negative and np.any(distance.value < 0):
            raise ValueError("Distance must be >= 0.  Use the argument "
                             "'allow_negative=True' to allow negative values.")

        return distance

    def __quantity_subclass__(self, unit):
        if unit.is_equivalent(u.m):
            return Distance, True
        else:
            return super(Distance, self).__quantity_subclass__(unit)[0], False

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
        cosmology : ``Cosmology`` or `None`
            The cosmology to assume for this calculation, or `None` to use the
            current cosmology (see `astropy.cosmology` for details).

        Returns
        -------
        z : float
            The redshift of this distance given the provided ``cosmology``.
        """
        from ..cosmology import luminosity_distance
        from scipy import optimize

        # FIXME: array: need to make this calculation more vector-friendly

        f = lambda z, d, cos: (luminosity_distance(z, cos).value - d) ** 2
        return optimize.brent(f, (self.Mpc, cosmology))

    @property
    def distmod(self):
        """The distance modulus as a `~astropy.units.Quantity`"""
        val = 5. * np.log10(self.to(u.pc).value) - 5.
        return u.Quantity(val, u.mag)

    @classmethod
    def _distmod_to_pc(cls, dm):
        dm = u.Quantity(dm, u.mag)
        return cls(10 ** ((dm.value + 5) / 5.), u.pc, copy=False)


@deprecated('v0.4', alternative='astropy.coordinates.CartesianRepresentation')
class CartesianPoints(u.Quantity):
    """
    A cartesian representation of a point in three-dimensional space.

    Parameters
    ----------
    x : `~astropy.units.Quantity` or array-like
        The first cartesian coordinate or a single array or
        `~astropy.units.Quantity` where the first dimension is length-3.
    y : `~astropy.units.Quantity` or array-like, optional
        The second cartesian coordinate.
    z : `~astropy.units.Quantity` or array-like, optional
        The third cartesian coordinate.
    unit : `~astropy.units.UnitBase` object or `None`
        The physical unit of the coordinate values. If ``x``, ``y``, or ``z``
        are quantities, they will be converted to this unit.
    dtype : `~numpy.dtype`, optional
        See `~astropy.units.Quantity`. Must be given as a keyword argument.
    copy : bool, optional
        See `~astropy.units.Quantity`. Must be given as a keyword argument.

    Raises
    ------
    UnitsError
        If the units on ``x``, ``y``, and ``z`` do not match or an invalid
        unit is given.
    ValueError
        If ``y`` and ``z`` don't match ``x``'s shape or ``x`` is not length-3
    TypeError
        If incompatible array types are passed into ``x``, ``y``, or ``z``

    """

    #this ensures that __array_wrap__ gets called for ufuncs even when
    #where a quantity is first, like ``3*u.m + c``
    __array_priority__ = 10001

    def __new__(cls, x, y=None, z=None, unit=None, dtype=None, copy=True):
        if y is None and z is None:
            if len(x) != 3:
                raise ValueError('Input to CartesianPoints is not length 3')

            qarr = x
            if unit is None and hasattr(qarr, 'unit'):
                unit = qarr.unit  # for when a Quantity is given
        elif y is not None and z is not None:
            if unit is None:
                #they must all match units or this fails
                for coo in (x, y, z):
                    if hasattr(coo, 'unit'):
                        if unit is not None and coo.unit != unit:
                            raise u.UnitsError('Units for `x`, `y`, and `z` do '
                                               'not match in CartesianPoints   ')
                        unit = coo.unit
                #if `unit`  is still None at this point, it means none were
                #Quantties, which is fine, because it means the user wanted
                #the unit to be None
            else:
                #convert any that are like a Quantity to the given unit
                if hasattr(x, 'to'):
                    x = x.to(unit)
                if hasattr(y, 'to'):
                    y = y.to(unit)
                if hasattr(z, 'to'):
                    z = z.to(unit)

            qarr = [np.asarray(coo) for coo in (x, y, z)]
            if not (qarr[0].shape == qarr[1].shape == qarr[2].shape):
                raise ValueError("Shapes for `x`, `y`, and `z` don't match in "
                                 "CartesianPoints")
                #let the unit be whatever it is
        else:
            raise TypeError('Must give all of `x`, `y`, and `z` or just array in '
                            'CartesianPoints')
        try:
            unit = _convert_to_and_validate_length_unit(unit, True)
        except TypeError as e:
            raise u.UnitsError(str(e))

        try:
            qarr = np.asarray(qarr)
        except ValueError as e:
            raise TypeError(str(e))

        if qarr.dtype.kind not in 'iuf':
            raise TypeError("Unsupported dtype '{0}'".format(qarr.dtype))

        return u.Quantity.__new__(cls, qarr, unit, dtype=dtype, copy=copy)

    def __quantity_subclass__(self, unit):
        if unit.is_equivalent(u.m):
            return CartesianPoints, True
        else:
            return u.Quantity.__quantity_subclass__(unit)[0], False

    def __array_wrap__(self, obj, context=None):
        #always convert to CartesianPoints because all operations that would
        #screw up the units are killed by _convert_to_and_validate_length_unit
        obj = u.Quantity.__array_wrap__(obj, context=context)

        #always prefer self's unit, if possible
        if obj.unit.is_equivalent(self.unit):
            return obj.to(self.unit)
        else:
            return obj

    @property
    def x(self):
        """
        The second cartesian coordinate as a `~astropy.units.Quantity`.
        """
        return self.view(u.Quantity)[0]

    @property
    def y(self):
        """
        The second cartesian coordinate as a `~astropy.units.Quantity`.
        """
        return self.view(u.Quantity)[1]

    @property
    def z(self):
        """
        The third cartesian coordinate as a `~astropy.units.Quantity`.
        """
        return self.view(u.Quantity)[2]

    def to_spherical(self):
        """
        Converts to the spherical representation of this point.

        Returns
        -------
        r : `~astropy.units.Quantity`
            The radial coordinate (in the same units as this `CartesianPoints`).
        lat : `~astropy.units.Quantity`
            The spherical coordinates latitude.
        lon : `~astropy.units.Quantity`
            The spherical coordinates longitude.

        """
        from .angles import Latitude, Longitude

        rarr, latarr, lonarr = cartesian_to_spherical(self.x, self.y, self.z)

        r = Distance(rarr, unit=self.unit)
        lat = Latitude(latarr, unit=u.radian)
        lon = Longitude(lonarr, unit=u.radian)

        return r, lat, lon


def _convert_to_and_validate_length_unit(unit, allow_dimensionless=False):
    """
    raises UnitsError if not a length unit
    """
    try:
        unit = u.Unit(unit)
        assert (unit.is_equivalent(u.kpc) or
                allow_dimensionless and unit == u.dimensionless_unscaled)
    except (TypeError, AssertionError):
        raise u.UnitsError('Unit "{0}" is not a length type'.format(unit))

    return unit

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

    if np.isscalar(r) and np.isscalar(lat) and np.isscalar(lon):
        x = r * math.cos(lat) * math.cos(lon)
        y = r * math.cos(lat) * math.sin(lon)
        z = r * math.sin(lat)
    else:
        x = r * np.cos(lat) * np.cos(lon)
        y = r * np.cos(lat) * np.sin(lon)
        z = r * np.sin(lat)

    return x, y, z
