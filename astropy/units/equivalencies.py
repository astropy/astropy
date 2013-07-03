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
    r"""
    Return the equivalency pairs for the radio convention for velocity:

    Radio   :math:`V = c \frac{f_0 - f}{f_0}  ;  f(V) = f_0 ( 1 - V/c )`

    Parameters
    ----------
    rest : Quantity
        Any quantity supported by the standard spectral equivalencies
        (wavelength, energy, frequency)

    References
    ----------
    `http://www.gb.nrao.edu/~fghigo/gbtdoc/doppler.html`_

    Examples
    --------
    >>> import astropy.units as u
    >>> CO_restfreq = 115.27120*u.GHz  # rest frequency of 12 CO 1-0 in GHz
    >>> radio_CO_equiv = u.doppler_radio(CO_restfreq)
    >>> measured_freq = 115.2832*u.GHz
    >>> doppler_radio = measured_freq.to(u.km/u.s, equivalencies=radio_CO_equiv)
    >>> doppler_radio
    <Quantity -31.2090920889... km / s>
    """

    ckms = _si.c.to('km/s').value

    def to_vel_freq(x):
        restfreq = rest.to(si.Hz, equivalencies=spectral()).value
        return (restfreq-x) / (restfreq) * ckms

    def from_vel_freq(x):
        restfreq = rest.to(si.Hz, equivalencies=spectral()).value
        voverc = x/ckms
        return restfreq * (1-voverc)


    def to_vel_wav(x):
        restwav = rest.to(si.AA, spectral()).value
        return (x-restwav) / (x) * ckms

    def from_vel_wav(x):
        restwav = rest.to(si.AA, spectral()).value
        return restwav * ckms / (ckms-x)


    def to_vel_en(x):
        resten = rest.to(si.eV, equivalencies=spectral()).value
        return (resten-x) / (resten) * ckms

    def from_vel_en(x):
        resten = rest.to(si.eV, equivalencies=spectral()).value
        voverc = x/ckms
        return resten * (1-voverc)

    return [(si.Hz, si.km/si.s, to_vel_freq, from_vel_freq),
            (si.AA, si.km/si.s, to_vel_wav, from_vel_wav),
            (si.eV, si.km/si.s, to_vel_en, from_vel_en),
            ]


def doppler_optical(rest):
    r"""
    Return the equivalency pairs for the optical convention for velocity:

    Optical       :math:`V = c \frac{f_0 - f}{f  }  ;  f(V) = f_0 ( 1 + V/c )^{-1}`

    Parameters
    ----------
    rest : Quantity
        Any quantity supported by the standard spectral equivalencies
        (wavelength, energy, frequency)

    References
    ----------
    `http://www.gb.nrao.edu/~fghigo/gbtdoc/doppler.html`_

    Examples
    --------
    >>> import astropy.units as u
    >>> CO_restfreq = 115.27120*u.GHz  # rest frequency of 12 CO 1-0 in GHz
    >>> optical_CO_equiv = u.doppler_optical(CO_restfreq)
    >>> measured_freq = 115.2832*u.GHz
    >>> doppler_optical = measured_freq.to(u.km/u.s, equivalencies=optical_CO_equiv)
    >>> doppler_optical
    <Quantity -31.205843488... km / s>
    """

    ckms = _si.c.to('km/s').value

    def to_vel_freq(x):
        restfreq = rest.to(si.Hz, equivalencies=spectral()).value
        return ckms * (restfreq-x) / x

    def from_vel_freq(x):
        restfreq = rest.to(si.Hz, equivalencies=spectral()).value
        voverc = x/ckms
        return restfreq / (1+voverc)


    def to_vel_wav(x):
        restwav = rest.to(si.AA, spectral()).value
        return ckms * (x/restwav-1)

    def from_vel_wav(x):
        restwav = rest.to(si.AA, spectral()).value
        voverc = x/ckms
        return restwav * (1+voverc)


    def to_vel_en(x):
        resten = rest.to(si.eV, equivalencies=spectral()).value
        return ckms * (resten-x) / x

    def from_vel_en(x):
        resten = rest.to(si.eV, equivalencies=spectral()).value
        voverc = x/ckms
        return resten / (1+voverc)

    return [(si.Hz, si.km/si.s, to_vel_freq, from_vel_freq),
            (si.AA, si.km/si.s, to_vel_wav, from_vel_wav),
            (si.eV, si.km/si.s, to_vel_en, from_vel_en),
            ]


def doppler_relativistic(rest):
    r"""
    Return the equivalency pairs for the relativistic convention for velocity:

    Relativistic  :math:`V = c \frac{f_0^2 - f^2}{f_0^2 + f^2} ;  f(V) = f_0 \frac{\left(1 - (V/c)^2\right)^{1/2}}{(1+V/c)}`

    Parameters
    ----------
    rest : Quantity
        Any quantity supported by the standard spectral equivalencies
        (wavelength, energy, frequency)

    References
    ----------
    `http://www.gb.nrao.edu/~fghigo/gbtdoc/doppler.html`_

    Examples
    --------
    >>> import astropy.units as u
    >>> CO_restfreq = 115.27120*u.GHz  # rest frequency of 12 CO 1-0 in GHz
    >>> relativistic_CO_equiv = u.doppler_relativistic(CO_restfreq)
    >>> measured_freq = 115.2832*u.GHz
    >>> doppler_relativistic = measured_freq.to(u.km/u.s, equivalencies=relativistic_CO_equiv)
    >>> doppler_relativistic
    <Quantity -31.2074676194... km / s>
    >>> measured_velocity = 1250 * u.km/u.s
    >>> relativistic_frequency = measured_velocity.to(u.GHz, equivalencies=relativistic_CO_equiv)
    >>> relativistic_frequency
    <Quantity 114.79156867... GHz>
    >>> relativistic_wavelength = measured_velocity.to(u.mm, equivalencies=relativistic_CO_equiv)
    >>> relativistic_wavelength
    <Quantity 2.61162436818... mm>
    """

    ckms = _si.c.to('km/s').value

    def to_vel_freq(x):
        restfreq = rest.to(si.Hz, equivalencies=spectral()).value
        return (restfreq**2-x**2) / (restfreq**2+x**2) * ckms

    def from_vel_freq(x):
        restfreq = rest.to(si.Hz, equivalencies=spectral()).value
        voverc = x/ckms
        return restfreq * ((1-voverc) / (1+(voverc)))**0.5


    def to_vel_wav(x):
        restwav = rest.to(si.AA, spectral()).value
        return (x**2-restwav**2) / (restwav**2+x**2) * ckms

    def from_vel_wav(x):
        restwav = rest.to(si.AA, spectral()).value
        voverc = x/ckms
        return restwav * ((1+voverc) / (1-voverc))**0.5


    def to_vel_en(x):
        resten = rest.to(si.eV, spectral()).value
        return (resten**2-x**2) / (resten**2+x**2) * ckms

    def from_vel_en(x):
        resten = rest.to(si.eV, spectral()).value
        voverc = x/ckms
        return resten * ((1-voverc) / (1+(voverc)))**0.5

    return [(si.Hz, si.km/si.s, to_vel_freq, from_vel_freq),
            (si.AA, si.km/si.s, to_vel_wav, from_vel_wav),
            (si.eV, si.km/si.s, to_vel_en, from_vel_en),
            ]
