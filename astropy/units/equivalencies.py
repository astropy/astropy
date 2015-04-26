# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""A set of standard astronomical equivalencies."""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

# THIRD-PARTY
import numpy as np

# LOCAL
from ..constants import si as _si
from . import si
from . import cgs
from . import astrophys
from .function import units as function_units
from . import dimensionless_unscaled


__all__ = ['parallax', 'spectral', 'spectral_density', 'doppler_radio',
           'doppler_optical', 'doppler_relativistic', 'mass_energy',
           'brightness_temperature', 'dimensionless_angles',
           'logarithmic', 'temperature', 'temperature_energy']


def dimensionless_angles():
    """Allow angles to be equivalent to dimensionless (with 1 rad = 1 m/m = 1).

    It is special compared to other equivalency pairs in that it
    allows this independent of the power to which the angle is raised,
    and independent of whether it is part of a more complicated unit.
    """
    return [(si.radian, None)]


def logarithmic():
    """Allow logarithmic units to be converted to dimensionless fractions"""
    return [
        (dimensionless_unscaled, function_units.dex,
         lambda x: np.log10(x), lambda x: 10.**x)
    ]


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
    wavelength, wave number, frequency, and energy equivalences.

    Allows conversions between wavelength units, wave number units,
    frequency units, and energy units as they relate to light.

    There are two types of wave number:

        * spectroscopic - :math:`1 / \\lambda` (per meter)
        * angular - :math:`2 \\pi / \\lambda` (radian per meter)

    """
    hc = _si.h.value * _si.c.value
    two_pi = 2.0 * np.pi
    inv_m_spec = si.m ** -1
    inv_m_ang = si.radian / si.m

    return [
        (si.m, si.Hz, lambda x: _si.c.value / x),
        (si.m, si.J, lambda x: hc / x),
        (si.Hz, si.J, lambda x: _si.h.value * x, lambda x: x / _si.h.value),
        (si.m, inv_m_spec, lambda x: 1.0 / x),
        (si.Hz, inv_m_spec, lambda x: x / _si.c.value,
         lambda x: _si.c.value * x),
        (si.J, inv_m_spec, lambda x: x / hc, lambda x: hc * x),
        (inv_m_spec, inv_m_ang, lambda x: x * two_pi, lambda x: x / two_pi),
        (si.m, inv_m_ang, lambda x: two_pi / x),
        (si.Hz, inv_m_ang, lambda x: two_pi * x / _si.c.value,
         lambda x: _si.c.value * x / two_pi),
        (si.J, inv_m_ang, lambda x: x * two_pi / hc, lambda x: hc * x / two_pi)
    ]


def spectral_density(wav, factor=None):
    """
    Returns a list of equivalence pairs that handle spectral density
    with regard to wavelength and frequency.

    Parameters
    ----------
    wav : `~astropy.units.Quantity`
        `~astropy.units.Quantity` associated with values being converted
        (e.g., wavelength or frequency).

    Notes
    -----
    The ``factor`` argument is left for backward-compatibility with the syntax
    ``spectral_density(unit, factor)`` but users are encouraged to use
    ``spectral_density(factor * unit)`` instead.

    """
    from .core import UnitBase

    if isinstance(wav, UnitBase):
        if factor is None:
            raise ValueError(
                'If `wav` is specified as a unit, `factor` should be set')
        wav = factor * wav   # Convert to Quantity

    c_Aps = _si.c.to(si.AA / si.s).value  # Angstrom/s
    h_cgs = _si.h.cgs.value  # erg * s
    hc = c_Aps * h_cgs

    fla = cgs.erg / si.angstrom / si.cm ** 2 / si.s
    fnu = cgs.erg / si.Hz / si.cm ** 2 / si.s
    nufnu = cgs.erg / si.cm ** 2 / si.s
    lafla = nufnu
    photlam = astrophys.photon / (si.cm ** 2 * si.s * si.AA)
    photnu = astrophys.photon / (si.cm ** 2 * si.s * si.Hz)

    def converter(x):
        return x * (wav.to(si.AA, spectral()).value ** 2 / c_Aps)

    def iconverter(x):
        return x / (wav.to(si.AA, spectral()).value ** 2 / c_Aps)

    def converter_fnu_nufnu(x):
        return x * wav.to(si.Hz, spectral()).value

    def iconverter_fnu_nufnu(x):
        return x / wav.to(si.Hz, spectral()).value

    def converter_fla_lafla(x):
        return x * wav.to(si.AA, spectral()).value

    def iconverter_fla_lafla(x):
        return x / wav.to(si.AA, spectral()).value

    def converter_photlam_fla(x):
        return hc * x / wav.to(si.AA, spectral()).value

    def iconverter_photlam_fla(x):
        return x * wav.to(si.AA, spectral()).value / hc

    def converter_photlam_fnu(x):
        return h_cgs * x * wav.to(si.AA, spectral()).value

    def iconverter_photlam_fnu(x):
        return x / (wav.to(si.AA, spectral()).value * h_cgs)

    def converter_photlam_photnu(x):
        return x * wav.to(si.AA, spectral()).value ** 2 / c_Aps

    def iconverter_photlam_photnu(x):
        return c_Aps * x / wav.to(si.AA, spectral()).value ** 2

    converter_photnu_fnu = converter_photlam_fla

    iconverter_photnu_fnu = iconverter_photlam_fla

    def converter_photnu_fla(x):
        return x * hc * c_Aps / wav.to(si.AA, spectral()).value ** 3

    def iconverter_photnu_fla(x):
        return x * wav.to(si.AA, spectral()).value ** 3 / (hc * c_Aps)

    return [
        (fla, fnu, converter, iconverter),
        (fnu, nufnu, converter_fnu_nufnu, iconverter_fnu_nufnu),
        (fla, lafla, converter_fla_lafla, iconverter_fla_lafla),
        (photlam, fla, converter_photlam_fla, iconverter_photlam_fla),
        (photlam, fnu, converter_photlam_fnu, iconverter_photlam_fnu),
        (photlam, photnu, converter_photlam_photnu, iconverter_photlam_photnu),
        (photnu, fnu, converter_photnu_fnu, iconverter_photnu_fnu),
        (photnu, fla, converter_photnu_fla, iconverter_photnu_fla)
    ]


def doppler_radio(rest):
    r"""
    Return the equivalency pairs for the radio convention for velocity.

    The radio convention for the relation between velocity and frequency is:

    :math:`V = c \frac{f_0 - f}{f_0}  ;  f(V) = f_0 ( 1 - V/c )`

    Parameters
    ----------
    rest : `~astropy.units.Quantity`
        Any quantity supported by the standard spectral equivalencies
        (wavelength, energy, frequency, wave number).

    References
    ----------
    `NRAO site defining the conventions <http://www.gb.nrao.edu/~fghigo/gbtdoc/doppler.html>`_

    Examples
    --------
    >>> import astropy.units as u
    >>> CO_restfreq = 115.27120*u.GHz  # rest frequency of 12 CO 1-0 in GHz
    >>> radio_CO_equiv = u.doppler_radio(CO_restfreq)
    >>> measured_freq = 115.2832*u.GHz
    >>> radio_velocity = measured_freq.to(u.km/u.s, equivalencies=radio_CO_equiv)
    >>> radio_velocity  # doctest: +FLOAT_CMP
    <Quantity -31.209092088877583 km / s>
    """

    try:
        rest.to(u.GHz, u.spectral())
    except Exception as ex:
        raise u.UnitsError("The 'rest' value must be a spectral equivalent "
                           "(frequency, wavelength, or energy).")

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
    Return the equivalency pairs for the optical convention for velocity.

    The optical convention for the relation between velocity and frequency is:

    :math:`V = c \frac{f_0 - f}{f  }  ;  f(V) = f_0 ( 1 + V/c )^{-1}`

    Parameters
    ----------
    rest : `~astropy.units.Quantity`
        Any quantity supported by the standard spectral equivalencies
        (wavelength, energy, frequency, wave number).

    References
    ----------
    `NRAO site defining the conventions <http://www.gb.nrao.edu/~fghigo/gbtdoc/doppler.html>`_

    Examples
    --------
    >>> import astropy.units as u
    >>> CO_restfreq = 115.27120*u.GHz  # rest frequency of 12 CO 1-0 in GHz
    >>> optical_CO_equiv = u.doppler_optical(CO_restfreq)
    >>> measured_freq = 115.2832*u.GHz
    >>> optical_velocity = measured_freq.to(u.km/u.s, equivalencies=optical_CO_equiv)
    >>> optical_velocity  # doctest: +FLOAT_CMP
    <Quantity -31.20584348799674 km / s>
    """

    try:
        rest.to(u.GHz, u.spectral())
    except Exception as ex:
        raise u.UnitsError("The 'rest' value must be a spectral equivalent "
                           "(frequency, wavelength, or energy).")

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
    Return the equivalency pairs for the relativistic convention for velocity.

    The full relativistic convention for the relation between velocity and frequency is:

    :math:`V = c \frac{f_0^2 - f^2}{f_0^2 + f^2} ;  f(V) = f_0 \frac{\left(1 - (V/c)^2\right)^{1/2}}{(1+V/c)}`

    Parameters
    ----------
    rest : `~astropy.units.Quantity`
        Any quantity supported by the standard spectral equivalencies
        (wavelength, energy, frequency, wave number).

    References
    ----------
    `NRAO site defining the conventions <http://www.gb.nrao.edu/~fghigo/gbtdoc/doppler.html>`_

    Examples
    --------
    >>> import astropy.units as u
    >>> CO_restfreq = 115.27120*u.GHz  # rest frequency of 12 CO 1-0 in GHz
    >>> relativistic_CO_equiv = u.doppler_relativistic(CO_restfreq)
    >>> measured_freq = 115.2832*u.GHz
    >>> relativistic_velocity = measured_freq.to(u.km/u.s, equivalencies=relativistic_CO_equiv)
    >>> relativistic_velocity  # doctest: +FLOAT_CMP
    <Quantity -31.207467619351537 km / s>
    >>> measured_velocity = 1250 * u.km/u.s
    >>> relativistic_frequency = measured_velocity.to(u.GHz, equivalencies=relativistic_CO_equiv)
    >>> relativistic_frequency  # doctest: +FLOAT_CMP
    <Quantity 114.79156866993588 GHz>
    >>> relativistic_wavelength = measured_velocity.to(u.mm, equivalencies=relativistic_CO_equiv)
    >>> relativistic_wavelength  # doctest: +FLOAT_CMP
    <Quantity 2.6116243681798923 mm>
    """

    try:
        rest.to(u.GHz, u.spectral())
    except Exception as ex:
        raise u.UnitsError("The 'rest' value must be a spectral equivalent "
                           "(frequency, wavelength, or energy).")

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


def mass_energy():
    """
    Returns a list of equivalence pairs that handle the conversion
    between mass and energy.
    """

    return [(si.kg, si.J, lambda x: x * _si.c.value ** 2,
             lambda x: x / _si.c.value ** 2),
            (si.kg / si.m ** 2, si.J / si.m ** 2 ,
             lambda x: x * _si.c.value ** 2,
             lambda x: x / _si.c.value ** 2),
            (si.kg / si.m ** 3, si.J / si.m ** 3 ,
             lambda x: x * _si.c.value ** 2,
             lambda x: x / _si.c.value ** 2),
            (si.kg / si.s, si.J / si.s , lambda x: x * _si.c.value ** 2,
             lambda x: x / _si.c.value ** 2),
    ]

def brightness_temperature(beam_area, disp):
    """
    Defines the conversion between Jy/beam and "brightness temperature",
    :math:`T_B`, in Kelvins.  The brightness temperature is a unit very
    commonly used in radio astronomy.  See, e.g., "Tools of Radio Astronomy"
    (Wilson 2009) eqn 8.16 and eqn 8.19 (these pages are available on `google
    books
    <http://books.google.com/books?id=9KHw6R8rQEMC&pg=PA179&source=gbs_toc_r&cad=4#v=onepage&q&f=false>`__).

    :math:`T_B \equiv S_\\nu / \left(2 k \\nu^2 / c^2 \\right)`

    However, the beam area is essential for this computation: the brightness
    temperature is inversely proportional to the beam area

    Parameters
    ----------
    beam_area : Beam Area equivalent
        Beam area in angular units, i.e. steradian equivalent
    disp : `~astropy.units.Quantity` with spectral units
        The observed `spectral` equivalent `~astropy.units.Unit` (e.g.,
        frequency or wavelength)

    Examples
    --------
    Arecibo C-band beam::

        >>> import numpy as np
        >>> from astropy import units as u
        >>> beam_sigma = 50*u.arcsec
        >>> beam_area = 2*np.pi*(beam_sigma)**2
        >>> freq = 5*u.GHz
        >>> equiv = u.brightness_temperature(beam_area, freq)
        >>> u.Jy.to(u.K, equivalencies=equiv)  # doctest: +FLOAT_CMP
        3.526294429423223
        >>> (1*u.Jy).to(u.K, equivalencies=equiv)  # doctest: +FLOAT_CMP
        <Quantity 3.526294429423223 K>

    VLA synthetic beam::

        >>> bmaj = 15*u.arcsec
        >>> bmin = 15*u.arcsec
        >>> fwhm_to_sigma = 1./(8*np.log(2))**0.5
        >>> beam_area = 2.*np.pi*(bmaj*bmin/fwhm_to_sigma**2)
        >>> freq = 5*u.GHz
        >>> equiv = u.brightness_temperature(beam_area, freq)
        >>> u.Jy.to(u.K, equivalencies=equiv)  # doctest: +FLOAT_CMP
        7.065788175060084
    """
    beam = beam_area.to(si.sr).value
    nu = disp.to(si.GHz, spectral())

    def convert_Jy_to_K(x_jybm):
        factor = (2 * _si.k_B * si.K * nu**2 / _si.c**2).to(astrophys.Jy).value
        return (x_jybm / beam / factor)

    def convert_K_to_Jy(x_K):
        factor = (astrophys.Jy / (2 * _si.k_B * nu**2 / _si.c**2)).to(si.K).value
        return (x_K * beam / factor)

    return [(astrophys.Jy, si.K, convert_Jy_to_K, convert_K_to_Jy)]


def temperature():
    """Convert between Kelvin, Celsius, and Fahrenheit here because
    Unit and CompositeUnit cannot do addition or subtraction properly.
    """
    from .imperial import deg_F
    return [
        (si.K, si.deg_C, lambda x: x - 273.15, lambda x: x + 273.15),
        (si.deg_C, deg_F, lambda x: x * 1.8 + 32.0, lambda x: (x - 32.0) / 1.8),
        (si.K, deg_F, lambda x: (x - 273.15) * 1.8 + 32.0,
         lambda x: ((x - 32.0) / 1.8) + 273.15)]

def temperature_energy():
    """Convert between Kelvin and keV(eV) to an equivalent amount."""
    return [
        (si.K, si.eV, lambda x: x / (_si.e.value / _si.k_B),
         lambda x: x * (_si.e.value / _si.k_B))]
