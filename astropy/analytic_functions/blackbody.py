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


# Some platform implementations of expm1() are buggy and Numpy uses
# them anyways--the bug is that on certain large inputs it returns
# NaN instead of INF like it should (it should only return NaN on a
# NaN input
# See https://github.com/astropy/astropy/issues/4171
_has_buggy_expm1 = np.isnan(np.expm1(1000))


def blackbody_nu(in_x, temperature):
    """Calculate blackbody flux per steradian, :math:`B_{\\nu}(T)`.

    .. note::

        Use `numpy.errstate` to suppress Numpy warnings, if desired.

    .. warning::

        Output values might contain ``nan`` and ``inf``.

    Parameters
    ----------
    in_x : number, array-like, or `~astropy.units.Quantity`
        Frequency, wavelength, or wave number.
        If not a Quantity, it is assumed to be in Hz.

    temperature : number, array-like, or `~astropy.units.Quantity`
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
    if np.any(temp < 0):
        raise ValueError('Temperature should be positive: {0}'.format(temp))
    if not np.all(np.isfinite(freq)) or np.any(freq <= 0):
        warnings.warn('Input contains invalid wavelength/frequency value(s)',
                      AstropyUserWarning)

    log_boltz = const.h * freq / (const.k_B * temp)
    boltzm1 = np.expm1(log_boltz)

    if _has_buggy_expm1:
        # Replace incorrect nan results with infs--any result of 'nan' is
        # incorrect unless the input (in log_boltz) happened to be nan to begin
        # with.  (As noted in #4393 ideally this would be replaced by a version
        # of expm1 that doesn't have this bug, rather than fixing incorrect
        # results after the fact...)
        boltzm1_nans = np.isnan(boltzm1)
        if np.any(boltzm1_nans):
            if boltzm1.isscalar and not np.isnan(log_boltz):
                boltzm1 = np.inf
            else:
                boltzm1[np.where(~np.isnan(log_boltz) & boltzm1_nans)] = np.inf

    # Calculate blackbody flux
    bb_nu = (2.0 * const.h * freq ** 3 / (const.c ** 2 * boltzm1))
    flux = bb_nu.to(FNU, u.spectral_density(freq))

    return flux / u.sr  # Add per steradian to output flux unit


def blackbody_lambda(in_x, temperature):
    """Like :func:`blackbody_nu` but for :math:`B_{\\lambda}(T)`.

    Parameters
    ----------
    in_x : number, array-like, or `~astropy.units.Quantity`
        Frequency, wavelength, or wave number.
        If not a Quantity, it is assumed to be in Angstrom.

    temperature : number, array-like, or `~astropy.units.Quantity`
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
