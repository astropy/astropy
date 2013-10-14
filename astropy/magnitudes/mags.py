# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""This module contains basic magnitude classes."""
from __future__ import division, print_function

# THIRD-PARTY
import numpy as np

# LOCAL
from .. import units as u


__all__ = ['Magnitude', 'STMAG', 'ABMAG']


class Magnitude(u.Quantity):
    """Class to handle generic magnitude system.

    If ``value`` given is a `~astropy.units.quantity.Quantity` but not
    in the unit of magnitude, the following conversion is done:

    .. math::

        mag = -2.5 \\times log_{10} flux

    Parameters
    ----------
    value : number, `~astropy.units.quantity.Quantity` object, or sequence of `~astropy.units.quantity.Quantity` objects.
        The numerical value of this quantity in the units given by
        unit.  If a `~astropy.units.quantity.Quantity` or sequence of them,
        creates a new `~astropy.units.quantity.Quantity` object, converting
        to magnitude as needed.

    dtype : `~numpy.dtype`, optional
        The ``dtype`` of the resulting Numpy array or scalar that will
        hold the value.  If not provided, is is determined
        automatically from the input value.

    equivalencies : list of equivalence pairs, optional
        A list of equivalence pairs. See :ref:`unit_equivalencies`.

    copy : bool, optional
        If `True` (default), then the value is copied.  Otherwise, a copy
        will only be made if :func:`__array__` returns a copy, if obj is a
        nested sequence, or if a copy is needed to satisfy ``dtype``.
        (The `False` option is intended mostly for internal use, to speed
        up initialization where it is known a copy has been made already.
        Use with care.)

    """
    _builtin_zpt = 0.0  # Used by subclass

    def __new__(cls, value, dtype=None, equivalencies=[], copy=True):
        from ..utils.misc import isiterable

        if isiterable(value):
            mags = [cls._from_flux(v) for v in value]
        else:
            mags = cls._from_flux(value)

        return super(Magnitude, cls).__new__(
            cls, mags, unit=u.mag, dtype=dtype, equivalencies=equivalencies,
            copy=copy)

    @classmethod
    def _from_flux(cls, value):
        """Convert to magnitude.
        No check on flux unit. Just matter if it is a magnitude or not.

        Parameters
        ----------
        value : number or `~astropy.units.quantity.Quantity`

        Returns
        -------
        mags : number

        """
        if isinstance(value, u.Quantity):
            if value.unit.decompose() != u.mag:
                mags = -2.5 * np.log10(value.value) + cls._builtin_zpt
            else:
                mags = value.value
        else:
            mags = value

        return mags

    def to_flux(self, zeropoint=0.0):
        """Convert magnitude to flux.

        .. math::

            flux = 10^{-0.4 \\; (mag - zeropoint)}

        Parameters
        ----------
        zeropoint : number or `~astropy.units.quantity.Quantity`
            Magnitude zeropoint to apply.

        Returns
        -------
        flux : number
            Flux values. No unit is provided, just the numbers.

        Raises
        ------
        astropy.units.core.UnitsError
            Zeropoint unit is not a magnitude.

        ValueError
            Zeropoint value is not scalar or a valid number.

        """
        if isinstance(zeropoint, u.Quantity):
            if zeropoint.unit.decompose() != u.mag:
                raise u.UnitsError(
                    'Zeropoint {0} is not a magnitude.'.format(zeropoint))
            zpt = zeropoint.value
        else:
            zpt = zeropoint

        if not np.isscalar(zpt) or not isinstance(zpt, (int, long, float)):
            raise ValueError('Zeropoint value {0} is invalid.'.format(zpt))

        return 10 ** (-0.4 * (self.value - zpt))


class STMAG(Magnitude):
    """This magnitude system is defined such that an object with
    constant flux distribution
    :math:`F_{\\lambda} = 3.63 \\times 10^{-9} \\; erg \\; cm^{-2} \\; s^{-1} \\; \\AA^{-1}`
    at all wavelengths will have zero color at all wavelengths.

    If ``value`` given is a `~astropy.units.quantity.Quantity` but not
    in the unit of magnitude, the following conversion is done:

    .. math::

        mag = -2.5 \\times log_{10} flux - 21.1

    Parameters
    ----------
    value, dtype, equivalencies, copy
        See `Magnitude` for more information.

    """
    _builtin_zpt = -21.1

    def to_flux(self):
        """Convert magnitude to flux.

        .. math::

            flux = 10^{-0.4 \\; (mag + 21.1)}

        Returns
        -------
        flux : number
            Flux values. No unit is provided, just the numbers.

        """
        return super(STMAG, self).to_flux(zeropoint=self._builtin_zpt)


class ABMAG(Magnitude):
    """This magnitude system is defined such that an object with
    constant flux distribution
    :math:`F_{\\nu} = 3.63 \\times 10^{-20} \\; erg \\; cm^{-2} \\; s^{-1} \\; Hz^{-1}`
    at all wavelengths will have zero color at all wavelengths.

    If ``value`` given is a `~astropy.units.quantity.Quantity` but not
    in the unit of magnitude, the following conversion is done:

    .. math::

        mag = -2.5 \\times log_{10} flux - 48.6

    Parameters
    ----------
    value, dtype, equivalencies, copy
        See `Magnitude` for more information.

    """
    _builtin_zpt = -48.6

    def to_flux(self):
        """Convert magnitude to flux.

        .. math::

            flux = 10^{-0.4 \\; (mag + 48.6)}

        Returns
        -------
        flux : number
            Flux values. No unit is provided, just the numbers.

        """
        return super(ABMAG, self).to_flux(zeropoint=self._builtin_zpt)
