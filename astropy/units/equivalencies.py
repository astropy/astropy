# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
A set of standard astronomical equivalencies.
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from ..constants import si as _si
from . import si
from . import cgs
from . import astrophys

__all__ = ['parallax', 'spectral', 'spectral_density']


def parallax():
    """
    Returns a list of equivalence pairs that handle the conversion
    between parallax angle and distance.
    """
    return [
        (si.arcsecond, astrophys.parsec, lambda x: 1. / x)
    ]


def spectral():
    """
    Returns a list of equivalence pairs that handle spectral
    wavelength, frequency, and energy equivalences.

    Allows conversions between wavelength units, frequency units and
    energy units as they relate to light.
    """

    return [
        (si.m, si.Hz, lambda x: _si.c.value / x),
        (si.m, si.J, lambda x: (_si.c.value * _si.h.value) / x),
        (si.Hz, si.J, lambda x: _si.h.value * x)
    ]


def spectral_density(sunit, sfactor):
    """
    Returns a list of equivalence pairs that handle spectral density
    with regard to wavelength and frequency.
    """
    c_Aps = _si.c.value * 10 ** 10

    fla = cgs.erg / si.angstrom / si.cm ** 2 / si.s
    fnu = cgs.erg / si.Hz / si.cm ** 2 / si.s
    nufnu = cgs.erg / si.cm ** 2 / si.s
    lafla = nufnu

    def converter(x):
        return x * (sunit.to(si.AA, sfactor, spectral()) ** 2 / c_Aps)

    def iconverter(x):
        return x / (sunit.to(si.AA, sfactor, spectral()) ** 2 / c_Aps)

    def converter_fnu_nufnu(x):
        return x * sunit.to(si.Hz, sfactor, spectral())

    def iconverter_fnu_nufnu(x):
        return x / sunit.to(si.Hz, sfactor, spectral())

    def converter_fla_lafla(x):
        return x * sunit.to(si.AA, sfactor, spectral())

    def iconverter_fla_lafla(x):
        return x / sunit.to(si.AA, sfactor, spectral())

    return [
        (si.AA, fnu, converter, iconverter),
        (fla, fnu, converter, iconverter),
        (si.AA, si.Hz, converter, iconverter),
        (fla, si.Hz, converter, iconverter),
        (fnu, nufnu, converter_fnu_nufnu, iconverter_fnu_nufnu),
        (fla, lafla, converter_fla_lafla, iconverter_fla_lafla),
    ]
