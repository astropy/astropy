# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Functions related to blackbody radiation."""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

# STDLIB
import warnings

# THIRD-PARTY
import numpy as np

# LOCAL
from .. import constants as const
from .. import units as u
from ..utils.exceptions import AstropyUserWarning


__all__ = ['blackbody_nu', 'blackbody_lambda']

# Units
FNU = u.erg / (u.cm**2 * u.s * u.Hz)
FLAM = u.erg / (u.cm**2 * u.s * u.AA)


def blackbody_nu(in_x, temperature):
    """Calculate blackbody flux per steradian, :math:`B_{\\nu}(T)`.

    .. note::

        Use `numpy.errstate` to suppress Numpy warnings, if desired.

    .. warning::

        Output values might contain ``nan`` and ``inf``.

    Parameters
    ----------
    in_x : number, array_like, or `~astropy.units.Quantity`
        Frequency, wavelength, or wave number.
        If not a Quantity, it is assumed to be in Hz.

    temperature : number or `~astropy.units.Quantity`
        Blackbody temperature.
        If not a Quantity, it is assumed to be in Kelvin.

    Returns
    -------
    flux : `~astropy.units.Quantity`
        Blackbody monochromatic flux in
        :math:`erg \\; cm^{-2} s^{-1} Hz^{-1} sr^{-1}`.

    Raises
    ------
    ValueError
        Invalid temperature.

    ZeroDivisionError
        Wavelength is zero (when converting to frequency).

    """
    # Convert to units for calculations, also force double precision
    with u.add_enabled_equivalencies(u.spectral() + u.temperature()):
        freq = u.Quantity(in_x, u.Hz, dtype=np.float64)
        temp = u.Quantity(temperature, u.K, dtype=np.float64)

    # Check if input values are physically possible
    if temp < 0:
        raise ValueError('Invalid temperature {0}'.format(temp))
    if np.any(freq <= 0):  # pragma: no cover
        warnings.warn('Input contains invalid wavelength/frequency value(s)',
                      AstropyUserWarning)

    # Calculate blackbody flux
    bb_nu = (2.0 * const.h * freq ** 3 /
             (const.c ** 2 * np.expm1(const.h * freq / (const.k_B * temp))))
    flux = bb_nu.to(FNU, u.spectral_density(freq))

    return flux / u.sr  # Add per steradian to output flux unit


def blackbody_lambda(in_x, temperature):
    """Like :func:`blackbody_nu` but for :math:`B_{\\lambda}(T)`.

    Parameters
    ----------
    in_x : number, array_like, or `~astropy.units.Quantity`
        Frequency, wavelength, or wave number.
        If not a Quantity, it is assumed to be in Angstrom.

    temperature : number or `~astropy.units.Quantity`
        Blackbody temperature.
        If not a Quantity, it is assumed to be in Kelvin.

    Returns
    -------
    flux : `~astropy.units.Quantity`
        Blackbody monochromatic flux in
        :math:`erg \\; cm^{-2} s^{-1} \\AA^{-1} sr^{-1}`.

    """
    if getattr(in_x, 'unit', None) is None:
        in_x = u.Quantity(in_x, u.AA)

    bb_nu = blackbody_nu(in_x, temperature) * u.sr  # Remove sr for conversion
    flux = bb_nu.to(FLAM, u.spectral_density(in_x))

    return flux / u.sr  # Add per steradian to output flux unit
