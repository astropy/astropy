# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
A set of standard astronomical equivalencies.
"""
from __future__ import absolute_import, division, print_function, unicode_literals

from .._constants import si as _si
from . import si
from . import cgs

__all__ = ['spectral', 'spectral_density']


def spectral():
    """
    Returns a list of equivalence pairs that handle spectral
    wavelength, frequency, and energy equivalences.

    Allows conversions between wavelength units, frequency units and
    energy units as they relate to light.
    """

    return [
        (si.m, si.Hz, lambda x: _si.c / x),
        (si.m, si.J, lambda x: (_si.c * _si.h) / x),
        (si.Hz, si.J, lambda x: _si.h * x)
    ]


def spectral_density(sunit, sfactor):
    """
    Returns a list of equivalence pairs that handle spectral density
    with regard to wavelength and frequency.
    """
    c_Aps = _si.c * 10 ** 10

    flambda = cgs.erg / si.angstrom / si.cm ** 2 / si.s
    fnu = cgs.erg / si.Hz / si.cm ** 2 / si.s
    nufnu = cgs.erg / si.cm ** 2 / si.s

    def converter(x):
        return x * (sunit.to(si.AA, sfactor, spectral()) ** 2 / c_Aps)

    def iconverter(x):
        return x / (sunit.to(si.AA, sfactor, spectral()) ** 2 / c_Aps)

    def converter_fnu_nufnu(x):
        return x * sunit.to(si.Hz, sfactor, spectral())

    def iconverter_fnu_nufnu(x):
        return x / sunit.to(si.Hz, sfactor, spectral())

    return [
        (si.AA, fnu, converter, iconverter),
        (flambda, fnu, converter, iconverter),
        (si.AA, si.Hz, converter, iconverter),
        (flambda, si.Hz, converter, iconverter),
        (fnu, nufnu, converter_fnu_nufnu, iconverter_fnu_nufnu),
        ]
