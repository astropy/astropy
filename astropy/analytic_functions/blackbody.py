# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Functions related to blackbody radiation."""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

# THIRD-PARTY
import numpy as np

# LOCAL
from .. import constants as const
from .. import units as u


__all__ = ['planck_lambda']


def planck_lambda(w, temp):
    """Calculate blackbody flux given wavelength or frequency
    and temperature.

    .. note:: No check for overflow.

    Parameters
    ----------
    w : `~astropy.units.Quantity`
        Frequency, wavelength, or wave number.

    temp : `~astropy.units.Quantity`
        Temperature.

    Returns
    -------
    flux : `~astropy.units.Quantity`
        Blackbody monochromatic flux in :math:`W m^{-2} \\mu m^{-1} sr^{-1}`

    """
    try:
        w = w.to(u.m, equivalencies=u.spectral())
    except u.UnitsError:
        raise ValueError(
            'First parameter must be in wavelength or frequency units')

    try:
        temp = temp.to(u.K, equivalencies=u.temperature())
    except u.UnitsError:
        raise ValueError(
            'Input temperature quantity must be in temperature units')

    flux = (2 * const.h * const.c**2 / w**5 /
            np.expm1(const.h*const.c/(w * const.k_B * temp)) / u.sr)
    flux = flux.to(u.Watt / (u.m * u.m * u.um * u.sr))

    return flux
