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

__all__ = ['parallax', 'spectral', 'spectral_density', 'doppler_radio',
           'doppler_optical', 'doppler_relativistic']

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

def doppler_radio(rest):
    """
    Return the equivalency pairs for the radio convention for velocity:

    http://www.gb.nrao.edu/~fghigo/gbtdoc/doppler.html
    Radio   :math:`V = c \frac{f_0 - f}{f_0}  ;  f(V) = f_0 ( 1 - V/c )`

    Parameters
    ----------
    rest : Quantity
        Any quantity supported by the standard spectral equivalencies
        (wavelength, energy, frequency)

    Example
    -------
    >>> CO_restfreq = 115.27120*u.GHz  # rest frequency of 12 CO 1-0 in GHz
    >>> radio_CO_equiv = doppler_radio(CO_restfreq)
    >>> measured_freq = 115.2832*u.GHz
    >>> doppler_radio = measured_freq.to(u.km/u.s, equivalencies=radio_CO_equiv)
    >>> print doppler_radio
    -31.2090920889 km / s
    """
    restfreq = rest.to(si.Hz, spectral())

    return [(si.Hz, si.km/si.s,
            lambda x: (restfreq.value-x) / restfreq.value * _si.c.to('km/s').value,
            lambda x: (1-x/_si.c.to('km/s').value) * restfreq)]

def doppler_optical(rest):
    """
    Return the equivalency pairs for the optical convention for velocity:

    http://www.gb.nrao.edu/~fghigo/gbtdoc/doppler.html
    Optical       :math:`V = c \frac{f_0 - f}{f  }  ;  f(V) = f_0 ( 1 + V/c )^{-1}`

    Parameters
    ----------
    rest : Quantity
        Any quantity supported by the standard spectral equivalencies
        (wavelength, energy, frequency)

    Example
    -------
    >>> CO_restfreq = 115.27120*u.GHz  # rest frequency of 12 CO 1-0 in GHz
    >>> optical_CO_equiv = doppler_optical(CO_restfreq)
    >>> measured_freq = 115.2832*u.GHz
    >>> doppler_optical = measured_freq.to(u.km/u.s, equivalencies=optical_CO_equiv)
    >>> print doppler_optical
    -31.205843488 km / s
    """
    restfreq = rest.to(si.Hz, spectral())

    return [(si.Hz, si.km/si.s,
            lambda x: (restfreq.value-x) / x * _si.c.to('km/s').value,
            lambda x: (1+x/_si.c.to('km/s').value)**(-1) * restfreq)]

def doppler_relativistic(rest):
    """
    Return the equivalency pairs for the relativistic convention for velocity:

    http://www.gb.nrao.edu/~fghigo/gbtdoc/doppler.html
    Relativistic  :math:`V = c \frac{f_0^2 - f^2}{f_0^2 + f^2} ;  f(V) = f_0 \frac{\left(1 - (V/c)^2\right)^{1/2}}{(1+V/c)}`

    Parameters
    ----------
    rest : Quantity
        Any quantity supported by the standard spectral equivalencies
        (wavelength, energy, frequency)

    Example
    -------
    >>> CO_restfreq = 115.27120*u.GHz  # rest frequency of 12 CO 1-0 in GHz
    >>> relativistic_CO_equiv = doppler_relativistic(CO_restfreq)
    >>> measured_freq = 115.2832*u.GHz
    >>> doppler_relativistic = measured_freq.to(u.km/u.s, equivalencies=relativistic_CO_equiv)
    >>> print doppler_relativistic
    -31.2074676194 km / s
    >>> measured_velocity = 1250 * u.km/u.s
    >>> relativistic_frequency = measured_freq.to(u.GHz, equivalencies=relativistic_CO_equiv)
    >>> print relativistic_frequency
    115.2832 GHz
    """
    restfreq = rest.to(si.Hz, spectral())

    return [(si.Hz, si.km/si.s,
            lambda x: (restfreq.value**2-x**2) / (restfreq.value**2+x**2) * _si.c.to('km/s').value,
            lambda x: (1-(x/_si.c.to('km/s').value)**2)**0.5 / (1+(x/_si.c.to('km/s').value)))]
