# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Functions related to blackbody radiation."""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

# LOCAL
from ..modeling import blackbody as _bb
from ..utils.decorators import deprecated


__all__ = ['blackbody_nu', 'blackbody_lambda']

# Units
FNU = _bb.FNU
FLAM = _bb.FLAM


@deprecated('2.0', alternative='astropy.modeling.blackbody.blackbody_nu')
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
    return _bb.blackbody_nu(in_x, temperature)


@deprecated('2.0', alternative='astropy.modeling.blackbody.blackbody_lambda')
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
        :math:`erg \\; cm^{-2} s^{-1} \\mathring{A}^{-1} sr^{-1}`.

    """
    return _bb.blackbody_lambda(in_x, temperature)
