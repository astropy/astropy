# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Functions related to blackbody radiation."""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

# THIRD-PARTY
import numpy as np

# LOCAL
from .. import constants as const
from .. import units as u


__all__ = ['planck']

# Underflow and overflow limits
_VERY_SMALL = 1e-4
_VERY_LARGE = 85.0

# Constants and units
HC = const.h.cgs * const.c.to('AA/s')
FNU = u.erg / (u.cm**2 * u.s * u.Hz)


def planck(in_x, temperature, flux_unit=FNU, silence_numpy=True):
    """Calculate blackbody flux per steradian.

    .. warning::

        Data points where overflow or underflow occurs will be set
        to zeroes.

    Parameters
    ----------
    in_x : number, array_like, or `~astropy.units.Quantity`
        Frequency, wavelength, or wave number.
        If not a Quantity, it is assumed to be in Angstrom.

    temp : number or `~astropy.units.Quantity`
        Blackbody temperature.
        If not a Quantity, it is assumed to be in Kelvin.

    flux_unit : `~astropy.units.Unit`
        Flux unit for the blackbody radiation.
        By default, it is :math:`erg \\; cm^{-2} s^{-1} Hz^{-1}`.
        Any given unit must be convertible from the default with
        :func:`~astropy.units.equivalencies.spectral_density`.

    silence_numpy : bool
        By default, this function temporarily silences Numpy warnings
        about overflow and underflow because they are already handled
        inside the function.

    Returns
    -------
    flux : `~astropy.units.Quantity`
        Blackbody monochromatic flux in given unit per steradian.

    """
    # Temporarily silence Numpy
    if silence_numpy:
        old_np_err_cfg = np.seterr(all='ignore')

    # Convert to Angstrom for calculations
    if not isinstance(in_x, u.Quantity):
        w = in_x * u.AA
    else:
        w = in_x.to(u.AA, u.spectral())

    # Convert to Kelvin for calculations
    if not isinstance(temperature, u.Quantity):
        temp = temperature * u.K
    else:
        temp = temperature.to(u.K, u.temperature())

    # Force double precision
    w = w.astype(np.float64)
    temp = temp.astype(np.float64)

    x = w * temp

    # Catch division by zero
    mask = x > 0
    x = np.where(mask, HC / (const.k_B.cgs * x), 0.0)

    # Catch overflow/underflow
    mask = (x >= _VERY_SMALL) & (x < _VERY_LARGE)
    factor = np.where(mask, 1.0 / np.expm1(x), 0.0)

    # Calculate blackbody flux
    freq = w.to(u.Hz, u.spectral())
    flux = 2.0 * const.h * factor * freq * freq * freq / const.c ** 2
    flux = flux.to(flux_unit, u.spectral_density(w))

    # Catch NaN
    flux[factor == 0] = 0.0

    # Add per steradian to output flux unit
    flux /= 1 * u.sr

    # Restore Numpy settings
    if silence_numpy:
        np.seterr(**old_np_err_cfg)

    return flux
