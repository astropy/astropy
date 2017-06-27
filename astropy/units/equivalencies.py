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
from .core import UnitsError


__all__ = ['parallax', 'spectral', 'spectral_density', 'doppler_radio',
           'doppler_optical', 'doppler_relativistic', 'mass_energy',
           'brightness_temperature', 'dimensionless_angles',
           'logarithmic', 'temperature', 'temperature_energy', 'molar_mass_amu',
           'pixel_scale', 'plate_scale']


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
         np.log10, lambda x: 10.**x)
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

    c_Aps = _si.c.to_value(si.AA / si.s)  # Angstrom/s
    h_cgs = _si.h.cgs.value  # erg * s
    hc = c_Aps * h_cgs

    # flux density
    f_la = cgs.erg / si.angstrom / si.cm ** 2 / si.s
    f_nu = cgs.erg / si.Hz / si.cm ** 2 / si.s
    nu_f_nu = cgs.erg / si.cm ** 2 / si.s
    la_f_la = nu_f_nu
    phot_f_la = astrophys.photon / (si.cm ** 2 * si.s * si.AA)
    phot_f_nu = astrophys.photon / (si.cm ** 2 * si.s * si.Hz)

    # luminosity density
    L_nu = cgs.erg / si.s / si.Hz
    L_la = cgs.erg / si.s / si.angstrom
    nu_L_nu = cgs.erg / si.s
    la_L_la = nu_L_nu
    phot_L_la = astrophys.photon / (si.s * si.AA)
    phot_L_nu = astrophys.photon / (si.s * si.Hz)

    def converter(x):
        return x * (wav.to_value(si.AA, spectral()) ** 2 / c_Aps)

    def iconverter(x):
        return x / (wav.to_value(si.AA, spectral()) ** 2 / c_Aps)

    def converter_f_nu_to_nu_f_nu(x):
        return x * wav.to_value(si.Hz, spectral())

    def iconverter_f_nu_to_nu_f_nu(x):
        return x / wav.to_value(si.Hz, spectral())

    def converter_f_la_to_la_f_la(x):
        return x * wav.to_value(si.AA, spectral())

    def iconverter_f_la_to_la_f_la(x):
        return x / wav.to_value(si.AA, spectral())

    def converter_phot_f_la_to_f_la(x):
        return hc * x / wav.to_value(si.AA, spectral())

    def iconverter_phot_f_la_to_f_la(x):
        return x * wav.to_value(si.AA, spectral()) / hc

    def converter_phot_f_la_to_f_nu(x):
        return h_cgs * x * wav.to_value(si.AA, spectral())

    def iconverter_phot_f_la_to_f_nu(x):
        return x / (wav.to_value(si.AA, spectral()) * h_cgs)

    def converter_phot_f_la_phot_f_nu(x):
        return x * wav.to_value(si.AA, spectral()) ** 2 / c_Aps

    def iconverter_phot_f_la_phot_f_nu(x):
        return c_Aps * x / wav.to_value(si.AA, spectral()) ** 2

    converter_phot_f_nu_to_f_nu = converter_phot_f_la_to_f_la
    iconverter_phot_f_nu_to_f_nu = iconverter_phot_f_la_to_f_la

    def converter_phot_f_nu_to_f_la(x):
        return x * hc * c_Aps / wav.to_value(si.AA, spectral()) ** 3

    def iconverter_phot_f_nu_to_f_la(x):
        return x * wav.to_value(si.AA, spectral()) ** 3 / (hc * c_Aps)

    # for luminosity density
    converter_L_nu_to_nu_L_nu = converter_f_nu_to_nu_f_nu
    iconverter_L_nu_to_nu_L_nu = iconverter_f_nu_to_nu_f_nu
    converter_L_la_to_la_L_la = converter_f_la_to_la_f_la
    iconverter_L_la_to_la_L_la = iconverter_f_la_to_la_f_la

    converter_phot_L_la_to_L_la = converter_phot_f_la_to_f_la
    iconverter_phot_L_la_to_L_la = iconverter_phot_f_la_to_f_la
    converter_phot_L_la_to_L_nu = converter_phot_f_la_to_f_nu
    iconverter_phot_L_la_to_L_nu = iconverter_phot_f_la_to_f_nu
    converter_phot_L_la_phot_L_nu = converter_phot_f_la_phot_f_nu
    iconverter_phot_L_la_phot_L_nu = iconverter_phot_f_la_phot_f_nu
    converter_phot_L_nu_to_L_nu = converter_phot_f_nu_to_f_nu
    iconverter_phot_L_nu_to_L_nu = iconverter_phot_f_nu_to_f_nu
    converter_phot_L_nu_to_L_la = converter_phot_f_nu_to_f_la
    iconverter_phot_L_nu_to_L_la = iconverter_phot_f_nu_to_f_la

    return [
        # flux
        (f_la, f_nu, converter, iconverter),
        (f_nu, nu_f_nu, converter_f_nu_to_nu_f_nu, iconverter_f_nu_to_nu_f_nu),
        (f_la, la_f_la, converter_f_la_to_la_f_la, iconverter_f_la_to_la_f_la),
        (phot_f_la, f_la, converter_phot_f_la_to_f_la, iconverter_phot_f_la_to_f_la),
        (phot_f_la, f_nu, converter_phot_f_la_to_f_nu, iconverter_phot_f_la_to_f_nu),
        (phot_f_la, phot_f_nu, converter_phot_f_la_phot_f_nu, iconverter_phot_f_la_phot_f_nu),
        (phot_f_nu, f_nu, converter_phot_f_nu_to_f_nu, iconverter_phot_f_nu_to_f_nu),
        (phot_f_nu, f_la, converter_phot_f_nu_to_f_la, iconverter_phot_f_nu_to_f_la),
        # luminosity
        (L_la, L_nu, converter, iconverter),
        (L_nu, nu_L_nu, converter_L_nu_to_nu_L_nu, iconverter_L_nu_to_nu_L_nu),
        (L_la, la_L_la, converter_L_la_to_la_L_la, iconverter_L_la_to_la_L_la),
        (phot_L_la, L_la, converter_phot_L_la_to_L_la, iconverter_phot_L_la_to_L_la),
        (phot_L_la, L_nu, converter_phot_L_la_to_L_nu, iconverter_phot_L_la_to_L_nu),
        (phot_L_la, phot_L_nu, converter_phot_L_la_phot_L_nu, iconverter_phot_L_la_phot_L_nu),
        (phot_L_nu, L_nu, converter_phot_L_nu_to_L_nu, iconverter_phot_L_nu_to_L_nu),
        (phot_L_nu, L_la, converter_phot_L_nu_to_L_la, iconverter_phot_L_nu_to_L_la),
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

    assert_is_spectral_unit(rest)

    ckms = _si.c.to_value('km/s')

    def to_vel_freq(x):
        restfreq = rest.to_value(si.Hz, equivalencies=spectral())
        return (restfreq-x) / (restfreq) * ckms

    def from_vel_freq(x):
        restfreq = rest.to_value(si.Hz, equivalencies=spectral())
        voverc = x/ckms
        return restfreq * (1-voverc)

    def to_vel_wav(x):
        restwav = rest.to_value(si.AA, spectral())
        return (x-restwav) / (x) * ckms

    def from_vel_wav(x):
        restwav = rest.to_value(si.AA, spectral())
        return restwav * ckms / (ckms-x)

    def to_vel_en(x):
        resten = rest.to_value(si.eV, equivalencies=spectral())
        return (resten-x) / (resten) * ckms

    def from_vel_en(x):
        resten = rest.to_value(si.eV, equivalencies=spectral())
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

    assert_is_spectral_unit(rest)

    ckms = _si.c.to_value('km/s')

    def to_vel_freq(x):
        restfreq = rest.to_value(si.Hz, equivalencies=spectral())
        return ckms * (restfreq-x) / x

    def from_vel_freq(x):
        restfreq = rest.to_value(si.Hz, equivalencies=spectral())
        voverc = x/ckms
        return restfreq / (1+voverc)

    def to_vel_wav(x):
        restwav = rest.to_value(si.AA, spectral())
        return ckms * (x/restwav-1)

    def from_vel_wav(x):
        restwav = rest.to_value(si.AA, spectral())
        voverc = x/ckms
        return restwav * (1+voverc)

    def to_vel_en(x):
        resten = rest.to_value(si.eV, equivalencies=spectral())
        return ckms * (resten-x) / x

    def from_vel_en(x):
        resten = rest.to_value(si.eV, equivalencies=spectral())
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

    assert_is_spectral_unit(rest)

    ckms = _si.c.to_value('km/s')

    def to_vel_freq(x):
        restfreq = rest.to_value(si.Hz, equivalencies=spectral())
        return (restfreq**2-x**2) / (restfreq**2+x**2) * ckms

    def from_vel_freq(x):
        restfreq = rest.to_value(si.Hz, equivalencies=spectral())
        voverc = x/ckms
        return restfreq * ((1-voverc) / (1+(voverc)))**0.5

    def to_vel_wav(x):
        restwav = rest.to_value(si.AA, spectral())
        return (x**2-restwav**2) / (restwav**2+x**2) * ckms

    def from_vel_wav(x):
        restwav = rest.to_value(si.AA, spectral())
        voverc = x/ckms
        return restwav * ((1+voverc) / (1-voverc))**0.5

    def to_vel_en(x):
        resten = rest.to_value(si.eV, spectral())
        return (resten**2-x**2) / (resten**2+x**2) * ckms

    def from_vel_en(x):
        resten = rest.to_value(si.eV, spectral())
        voverc = x/ckms
        return resten * ((1-voverc) / (1+(voverc)))**0.5

    return [(si.Hz, si.km/si.s, to_vel_freq, from_vel_freq),
            (si.AA, si.km/si.s, to_vel_wav, from_vel_wav),
            (si.eV, si.km/si.s, to_vel_en, from_vel_en),
            ]


def molar_mass_amu():
    """
    Returns the equivalence between amu and molar mass.
    """
    return [
        (si.g/si.mol, astrophys.u)
    ]


def mass_energy():
    """
    Returns a list of equivalence pairs that handle the conversion
    between mass and energy.
    """

    return [(si.kg, si.J, lambda x: x * _si.c.value ** 2,
             lambda x: x / _si.c.value ** 2),
            (si.kg / si.m ** 2, si.J / si.m ** 2,
             lambda x: x * _si.c.value ** 2,
             lambda x: x / _si.c.value ** 2),
            (si.kg / si.m ** 3, si.J / si.m ** 3,
             lambda x: x * _si.c.value ** 2,
             lambda x: x / _si.c.value ** 2),
            (si.kg / si.s, si.J / si.s, lambda x: x * _si.c.value ** 2,
             lambda x: x / _si.c.value ** 2),
    ]


def brightness_temperature(beam_area, disp):
    r"""
    Defines the conversion between Jy/beam and "brightness temperature",
    :math:`T_B`, in Kelvins.  The brightness temperature is a unit very
    commonly used in radio astronomy.  See, e.g., "Tools of Radio Astronomy"
    (Wilson 2009) eqn 8.16 and eqn 8.19 (these pages are available on `google
    books
    <http://books.google.com/books?id=9KHw6R8rQEMC&pg=PA179&source=gbs_toc_r&cad=4#v=onepage&q&f=false>`__).

    :math:`T_B \equiv S_\nu / \left(2 k \nu^2 / c^2 \right)`

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
        >>> beam_area = 2.*np.pi*(bmaj*bmin*fwhm_to_sigma**2)
        >>> freq = 5*u.GHz
        >>> equiv = u.brightness_temperature(beam_area, freq)
        >>> u.Jy.to(u.K, equivalencies=equiv)  # doctest: +FLOAT_CMP
        217.2658703625732
    """
    beam = beam_area.to_value(si.sr)
    nu = disp.to(si.GHz, spectral())

    def convert_Jy_to_K(x_jybm):
        factor = (2 * _si.k_B * si.K * nu**2 / _si.c**2).to_value(astrophys.Jy)
        return (x_jybm / beam / factor)

    def convert_K_to_Jy(x_K):
        factor = (astrophys.Jy / (2 * _si.k_B * nu**2 / _si.c**2)).to_value(si.K)
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
        (si.K, si.eV, lambda x: x / (_si.e.value / _si.k_B.value),
         lambda x: x * (_si.e.value / _si.k_B.value))]


def assert_is_spectral_unit(value):
    try:
        value.to(si.Hz, spectral())
    except (AttributeError, UnitsError) as ex:
        raise UnitsError("The 'rest' value must be a spectral equivalent "
                         "(frequency, wavelength, or energy).")


def pixel_scale(pixscale):
    """
    Convert between pixel distances (in units of ``pix``) and angular units,
    given a particular ``pixscale``.

    Parameters
    ----------
    pixscale : `~astropy.units.Quantity`
        The pixel scale either in units of angle/pixel or pixel/angle.
    """
    if pixscale.unit.is_equivalent(si.arcsec/astrophys.pix):
        pixscale_val = pixscale.to_value(si.radian/astrophys.pix)
    elif pixscale.unit.is_equivalent(astrophys.pix/si.arcsec):
        pixscale_val = (1/pixscale).to_value(si.radian/astrophys.pix)
    else:
        raise UnitsError("The pixel scale must be in angle/pixel or "
                         "pixel/angle")

    return [(astrophys.pix, si.radian, lambda px: px*pixscale_val, lambda rad: rad/pixscale_val)]


def plate_scale(platescale):
    """
    Convert between lengths (to be interpreted as lengths in the focal plane)
    and angular units with a specified ``platescale``.

    Parameters
    ----------
    platescale : `~astropy.units.Quantity`
        The pixel scale either in units of distance/pixel or distance/angle.
    """
    if platescale.unit.is_equivalent(si.arcsec/si.m):
        platescale_val = platescale.to_value(si.radian/si.m)
    elif platescale.unit.is_equivalent(si.m/si.arcsec):
        platescale_val = (1/platescale).to_value(si.radian/si.m)
    else:
        raise UnitsError("The pixel scale must be in angle/distance or "
                         "distance/angle")

    return [(si.m, si.radian, lambda d: d*platescale_val, lambda rad: rad/platescale_val)]
