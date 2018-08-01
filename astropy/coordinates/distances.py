# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module contains the classes and utility functions for distance and
cartesian coordinates.
"""

import numpy as np

from .. import units as u
from .angles import Angle

__all__ = ['Distance']


__doctest_requires__ = {'*': ['scipy.integrate']}


class Distance(u.SpecificTypeQuantity):
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
    parallax : `~astropy.units.Quantity` or `~astropy.coordinates.Angle`
        The parallax in angular units.
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
    >>> d8 = Distance(parallax=21.34*u.mas)
    """

    _equivalent_unit = u.m
    _include_easy_conversion_members = True

    def __new__(cls, value=None, unit=None, z=None, cosmology=None,
                distmod=None, parallax=None, dtype=None, copy=True, order=None,
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

            value_msg = ('Should given only one of `value`, `z`, `distmod`, or '
                         '`parallax` in Distance constructor.')
            n_not_none = np.sum([x is not None
                                 for x in [value, z, distmod, parallax]])
            if n_not_none > 1:
                raise ValueError(value_msg)

            if distmod is not None:
                value = cls._distmod_to_pc(distmod)
                if unit is None:
                    # if the unit is not specified, guess based on the mean of
                    # the log of the distance
                    meanlogval = np.log10(value.value).mean()
                    if meanlogval > 6:
                        unit = u.Mpc
                    elif meanlogval > 3:
                        unit = u.kpc
                    elif meanlogval < -3:  # ~200 AU
                        unit = u.AU
                    else:
                        unit = u.pc

                # Continue on to take account of unit and other arguments
                # but a copy is already made, so no longer necessary
                copy = False

            elif parallax is not None:
                value = parallax.to(u.pc, equivalencies=u.parallax()).value
                unit = u.pc

                # Continue on to take account of unit and other arguments
                # but a copy is already made, so no longer necessary
                copy = False

            elif value is None:
                raise ValueError('None of `value`, `z`, `distmod`, or '
                                 '`parallax` were given to Distance '
                                 'constructor')

        # now we have arguments like for a Quantity, so let it do the work
        distance = super().__new__(
            cls, value, unit, dtype=dtype, copy=copy, order=order,
            subok=subok, ndmin=ndmin)

        if not allow_negative and np.any(distance.value < 0):
            raise ValueError("Distance must be >= 0.  Use the argument "
                             "'allow_negative=True' to allow negative values.")

        return distance

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

        if cosmology is None:
            from ..cosmology import default_cosmology
            cosmology = default_cosmology.get()

        from ..cosmology import z_at_value
        return z_at_value(cosmology.luminosity_distance, self, ztol=1.e-10)

    @property
    def distmod(self):
        """The distance modulus as a `~astropy.units.Quantity`"""
        val = 5. * np.log10(self.to_value(u.pc)) - 5.
        return u.Quantity(val, u.mag, copy=False)

    @classmethod
    def _distmod_to_pc(cls, dm):
        dm = u.Quantity(dm, u.mag)
        return cls(10 ** ((dm.value + 5) / 5.), u.pc, copy=False)

    @property
    def parallax(self):
        """The parallax angle as an `~astropy.coordinates.Angle` object"""
        return Angle(self.to(u.milliarcsecond, u.parallax()))
