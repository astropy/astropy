# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
A set of standard astronomical equivalencies.
"""
from __future__ import absolute_import, division, print_function, unicode_literals

from .._constants import si as _si
from . import si
from . import cgs
from . import astrophys

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

    # Flux density
    fla = cgs.erg / si.angstrom / si.cm ** 2 / si.s
    fnu = cgs.erg / si.Hz / si.cm ** 2 / si.s
    nufnu = cgs.erg / si.cm ** 2 / si.s
    lafla = nufnu

    # Luminosity density
    lla = cgs.erg / si.angstrom / si.s
    lnu = cgs.erg / si.Hz / si.s
    nulnu = cgs.erg / si.s
    lalla = nulnu

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

    # Photons

    def converter_photons_to_ergs(x):
        return x * _si.h * sunit.to(si.Hz, sfactor, spectral())

    def converter_ergs_to_photons(x):
        return x / _si.h / sunit.to(si.Hz, sfactor, spectral())

    return [
        (fla, fnu, converter, iconverter),
        (fnu, nufnu, converter_fnu_nufnu, iconverter_fnu_nufnu),
        (fla, lafla, converter_fla_lafla, iconverter_fla_lafla),
        (lla, lnu, converter, iconverter),
        (lnu, nulnu, converter_fnu_nufnu, iconverter_fnu_nufnu),
        (lla, lalla, converter_fla_lafla, iconverter_fla_lafla),
        (astrophys.photon, cgs.erg, converter_photons_to_ergs, converter_ergs_to_photons),
        ]
