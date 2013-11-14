# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
A set of standard astronomical equivalencies.
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from math import pi
from ..constants import si as _si
from . import si
from . import cgs
from . import astrophys

__all__ = ['parallax', 'spectral', 'spectral_density', 'doppler_radio',
           'doppler_optical', 'doppler_relativistic', 'mass_energy',
           'brightness_temperature', 'dimensionless_angles', 
           'radio_lines_simple', 'radio_lines']


def dimensionless_angles():
    """Allow angles to be equivalent to dimensionless (with 1 rad = 1 m/m = 1).

    It is special compared to other equivalency pairs in that it
    allows this independent of the power to which the angle is raised,
    and indepedent of whether it is part of a more complicated unit.
    """
    return [(si.radian, None)]


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

    """
    hc = _si.h.value * _si.c.value
    inv_m = si.m ** -1
    return [
        (si.m, si.Hz, lambda x: _si.c.value / x),
        (si.m, si.J, lambda x: hc / x),
        (si.m, inv_m, lambda x: 1.0 / x),
        (si.Hz, si.J, lambda x: _si.h.value * x, lambda x: x / _si.h.value),
        (si.Hz, inv_m, lambda x: x / _si.c.value, lambda x: _si.c.value * x),
        (si.J, inv_m, lambda x: x / hc, lambda x: hc * x)
    ]


def spectral_density(wav, factor=None):
    """
    Returns a list of equivalence pairs that handle spectral density
    with regard to wavelength and frequency.

    Parameters
    ----------
    wav : Quantity
        Quantity associated with values being converted
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
                'If ``wav`` is specified as a unit, ``factor`` should be set')
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
        (si.AA, fnu, converter, iconverter),
        (fla, fnu, converter, iconverter),
        (si.AA, si.Hz, converter, iconverter),
        (fla, si.Hz, converter, iconverter),
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
    rest : Quantity
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
    >>> radio_velocity
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
    Return the equivalency pairs for the optical convention for velocity.

    The optical convention for the relation between velocity and frequency is:

    :math:`V = c \frac{f_0 - f}{f  }  ;  f(V) = f_0 ( 1 + V/c )^{-1}`

    Parameters
    ----------
    rest : Quantity
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
    >>> optical_velocity
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
    Return the equivalency pairs for the relativistic convention for velocity.

    The full relativistic convention for the relation between velocity and frequency is:

    :math:`V = c \frac{f_0^2 - f^2}{f_0^2 + f^2} ;  f(V) = f_0 \frac{\left(1 - (V/c)^2\right)^{1/2}}{(1+V/c)}`

    Parameters
    ----------
    rest : Quantity
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
    >>> relativistic_velocity
    <Quantity -31.207467619...
    >>> measured_velocity = 1250 * u.km/u.s
    >>> relativistic_frequency = measured_velocity.to(u.GHz, equivalencies=relativistic_CO_equiv)
    >>> relativistic_frequency
    <Quantity 114.7915686...
    >>> relativistic_wavelength = measured_velocity.to(u.mm, equivalencies=relativistic_CO_equiv)
    >>> relativistic_wavelength
    <Quantity 2.6116243681...
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
    "Antenna Gain" or "sensitivity" equivalency: Defines the conversion between
    Jy/beam and "brightness temperature", :math:`T_B`, in Kelvins.  This is a
    unit very commonly used in radio astronomy.  Typically, the gain refers to
    the conversion between corrected antenna temperature :math:`T_A^*` and flux
    density.  See, e.g., "Tools of Radio Astronomy" (Wilson 2009) eqn 8.16 and
    eqn 8.19 (these pages are available on `google books
    <http://books.google.com/books?id=9KHw6R8rQEMC&pg=PA179&source=gbs_toc_r&cad=4#v=onepage&q&f=false>`__).

    :math:`T_B \equiv S_\\nu / \left(2 k \\nu^2 / c^2 \\right)`

    However, the beam area is essential for this computation: the brighntess
    temperature is inversely proportional to the beam area

    Parameters
    ----------
    beam_area : Beam Area equivalent
        Beam area in angular units, i.e. steradian equivalent
    disp : `Quantity` with spectral units
        The observed `spectral` equivalent `Unit` (e.g., frequency or
        wavelength)

    Examples
    --------
    Arecibo C-band beam gain ~ 7 K/Jy::

        >>> import numpy as np
        >>> from astropy import units as u
        >>> beam_area = np.pi*(50*u.arcsec)**2
        >>> freq = 5*u.GHz
        >>> u.Jy.to(u.K, equivalencies=u.brightness_temperature(beam_area,freq))
        7.052588858...
        >>> (1*u.Jy).to(u.K, equivalencies=u.brightness_temperature(beam_area,freq))
        <Quantity 7.05258...

    VLA synthetic beam::

        >>> beam_area = np.pi*(15*u.arcsec)**2
        >>> freq = 5*u.GHz
        >>> u.Jy.to(u.K, equivalencies=u.brightness_temperature(beam_area,freq))
        78.36209843...
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

def radio_lines_simple(disp_obs):
    """ Returns equivalence pairs between Jy km/s, a line-strength
    unit often reported in radio astronomy, and SI units.

    Parameters
    ----------
    disp_obs : `Quantity` with spectral units
        The `spectral` equivalent `Unit` (e.g., frequency or
        wavelength) in the observer frame.
    """
    
    if not disp_obs.isscalar:
        raise ValueError("disp_obs must be scalar")
    fobs = disp_obs.to(si.Hz, equivalencies=spectral()).value
    c_si = _si.c.value # m / s
    return [(astrophys.Jy * si.km / si.s, si.W / si.m**2,
             lambda x: x * 1e-23 * fobs / c_si,
             lambda x: c_si * x * 1e23 / fobs)]
    

def radio_lines(disp_obs, disp_rest, lumdist):
    """
    Returns a list of equivalence pairs between observational
    and non-observational line strength units in astronomy.

    Observationally, line units are typically reported in Jy km/s.
    These are often converted to either solar luminosities or
    K km/s pc^2.  These equivalencies allow the user to go between these pairs.

    Parameters
    ----------
    disp_obs : `Quantity` with spectral units
        The `spectral` equivalent `Unit` (e.g., frequency or
        wavelength) in the observer frame.

    disp_rest : `Quantity` with spectral units
        The `spectral` equivalent `Unit` (e.g., frequency or
        wavelength) in the rest frame.

    lumdist : `Quantity` with units of distance
        The luminosity distance to the source
    """

    # The conversion from Jy km/s to W/m^2 is
    #     1 [Jy km/s] = 1e-26 * nu_obs / c_kms [W/m2]
    #  where 1e-26 is because of the definition of a jansky and c_kms
    #  is the speed of light in km/s

    # The conversion from W/m^2 to L_sun is
    #     1 [W/m^2] = 4 * pi * D_L^2 / Lsol [L_sun]
    # where D_L is the luminosity distance in m, and Lsol is the luminosity
    # of the sun in SI units

    # The conversion from K km/s pc^2 to W/m^2 is
    #     2 * 1000 * nu_obs**3 * k_B / (1+z) c^3 (D_A)**2
    # Where 2 nu^2/c^2 k_B comes from the R-J expression for a Blackbody
    # 1000 nu / c comes from km/s to frequency width
    # D_A is measured in pc because the expression is (L/D_A)**2 where
    #  L is the size of the emitting region, which is 1pc for K km/s pc^2

    # All the inter-conversoins follow from those.

    # Approximate versions (less precision on constants in front) 
    #  are given in Carilli and Walter, AARA, 2013, 51, section 2.4
    # Unfortunately, they don't derive the 'magic numbers' in front,
    # and I (@aconley) am not aware of anywhere where they are explained,
    # but they can be worked out by hand in terms of fundamental constants 
    # and unit conversions.

    speceqv = spectral() # Spectral equivalency
    if not disp_obs.isscalar:
        raise ValueError("disp_obs must be scalar")
    fobs = disp_obs.to(si.Hz, equivalencies=speceqv).value
    if not disp_rest.isscalar:
        raise ValueError("disp_rest must be scalar")
    frest = disp_rest.to(si.Hz, equivalencies=speceqv).value
    opz = frest / fobs

    dl = lumdist.to(astrophys.Mpc)
    if not dl.isscalar:
        raise ValueError("User provided luminosity distance must be scalar")
    if dl.value <= 0.0:
        raise ValueError("Invalid (non-positive) luminosity distance")

    # Set up quantities we will reuse below
    dl_mpc = dl.value # Mpc, convenient to keep around because of the pc^2
    dl_m = dl.to(si.m).value # m
    # SI 4pi D_L^2 / L_sun 
    fpidl2overLsun_si = 4.0 * pi * dl_m ** 2 / _si.L_sun.value 
    c_si = _si.c.value # m/s
    c2overkb = c_si ** 2 / _si.k_B.value # c^2 / Boltzmann in SI units
    c3overkb = c_si * c2overkb # c^3 / Boltzmann in SI units

    # This is L_sun * c**3 / (4 pi (1pc/1m)**2 * k_B) in SI units,
    # which shows up in the lsun to K km/s pc^2 conversion
    lsun_kkmspc2_const = _si.L_sun.value * c3overkb /\
                         (4.0 * pi * astrophys.parsec.to(si.m) ** 2)

    def jykms_to_wm2(x):
        return x * 1e-23 * fobs / c_si
    
    def wm2_to_jykms(x):
        return c_si * x * 1e23 / fobs

    def wm2_to_lsun(x):
        return x * fpidl2overLsun_si

    def lsun_to_wm2(x):
        return x / fpidl2overLsun_si

    def jykms_to_lsun(x):
        return x * 1e-23 * fobs * fpidl2overLsun_si / c_si

    def lsun_to_jykms(x):
        return 1e23 * x * c_si / (fobs * fpidl2overLsun_si)

    def wm2_to_kkmspc2(x):
        return 0.5e9 * x * c3overkb * dl_mpc**2 / (frest ** 3)
        
    def kkmspc2_to_wm2(x):
        return 2e-9 * x * frest ** 3 / (c3overkb * dl_mpc ** 2)

    def jykms_to_kkmspc2(x):
        return 0.5e-14 * x * c2overkb * dl_mpc**2 / (opz * frest ** 2)

    def kkmspc2_to_jykms(x):
        return 2e14 * x * opz * frest**2 / (c2overkb * dl_mpc**2)

    def lsun_to_kkmspc2(x):
        return x * 5e-4 * lsun_kkmspc2_const * frest ** (-3)

    def kkmspc2_to_lsun(x):
        return x * 2e3 * frest ** 3 / lsun_kkmspc2_const 
    
    wm2 = si.W / si.m**2
    jykms = astrophys.Jy * si.km / si.s
    kkmspc2 = si.K * si.km * astrophys.parsec**2 / si.s
    
    return [
        (wm2, jykms, wm2_to_jykms, jykms_to_wm2),
        (wm2, astrophys.solLum, wm2_to_lsun, lsun_to_wm2),
        (wm2, kkmspc2, wm2_to_kkmspc2, kkmspc2_to_wm2),
        (jykms, astrophys.solLum, jykms_to_lsun, lsun_to_jykms),
        (jykms, kkmspc2, jykms_to_kkmspc2, kkmspc2_to_jykms),
        (astrophys.solLum, kkmspc2, lsun_to_kkmspc2, kkmspc2_to_lsun)
    ]
