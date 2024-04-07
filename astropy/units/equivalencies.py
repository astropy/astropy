# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""A set of standard astronomical equivalencies."""

from collections import UserList

# THIRD-PARTY
import numpy as np

# LOCAL
from astropy.constants import si as _si
from astropy.utils.misc import isiterable

from . import astrophys, cgs, dimensionless_unscaled, misc, si
from .core import Unit, UnitsError
from .function import units as function_units

__all__ = [
    "parallax",
    "spectral",
    "spectral_density",
    "doppler_radio",
    "doppler_optical",
    "doppler_relativistic",
    "doppler_redshift",
    "mass_energy",
    "brightness_temperature",
    "thermodynamic_temperature",
    "beam_angular_area",
    "dimensionless_angles",
    "logarithmic",
    "temperature",
    "temperature_energy",
    "molar_mass_amu",
    "pixel_scale",
    "plate_scale",
    "Equivalency",
]


class Equivalency(UserList):
    """
    A container for a units equivalency.

    Attributes
    ----------
    name: `str`
        The name of the equivalency.
    kwargs: `dict`
        Any positional or keyword arguments used to make the equivalency.
    """

    def __init__(self, equiv_list, name="", kwargs=None):
        self.data = equiv_list
        self.name = [name]
        self.kwargs = [kwargs] if kwargs is not None else [{}]

    def __add__(self, other):
        if isinstance(other, Equivalency):
            new = super().__add__(other)
            new.name = self.name[:] + other.name
            new.kwargs = self.kwargs[:] + other.kwargs
            return new
        else:
            return self.data.__add__(other)

    def __eq__(self, other):
        return (
            isinstance(other, self.__class__)
            and self.name == other.name
            and self.kwargs == other.kwargs
        )


def dimensionless_angles():
    """Allow angles to be equivalent to dimensionless (with 1 rad = 1 m/m = 1).

    It is special compared to other equivalency pairs in that it
    allows this independent of the power to which the angle is raised,
    and independent of whether it is part of a more complicated unit.
    """
    return Equivalency([(si.radian, None)], "dimensionless_angles")


def logarithmic():
    """Allow logarithmic units to be converted to dimensionless fractions."""
    return Equivalency(
        [(dimensionless_unscaled, function_units.dex, np.log10, lambda x: 10.0**x)],
        "logarithmic",
    )


def parallax():
    """
    Returns a list of equivalence pairs that handle the conversion
    between parallax angle and distance.
    """

    def parallax_converter(x):
        x = np.asanyarray(x)
        d = 1 / x

        if isiterable(d):
            d[d < 0] = np.nan
            return d

        else:
            if d < 0:
                return np.array(np.nan)
            else:
                return d

    return Equivalency(
        [(si.arcsecond, astrophys.parsec, parallax_converter)], "parallax"
    )


def spectral():
    """
    Returns a list of equivalence pairs that handle spectral
    wavelength, wave number, frequency, and energy equivalencies.

    Allows conversions between wavelength units, wave number units,
    frequency units, and energy units as they relate to light.

    There are two types of wave number:

        * spectroscopic - :math:`1 / \\lambda` (per meter)
        * angular - :math:`2 \\pi / \\lambda` (radian per meter)

    """
    c = _si.c.value
    h = _si.h.value
    hc = h * c
    two_pi = 2.0 * np.pi
    inv_m_spec = si.m**-1
    inv_m_ang = si.radian / si.m

    return Equivalency(
        [
            (si.m, si.Hz, lambda x: c / x),
            (si.m, si.J, lambda x: hc / x),
            (si.Hz, si.J, lambda x: h * x, lambda x: x / h),
            (si.m, inv_m_spec, lambda x: 1.0 / x),
            (si.Hz, inv_m_spec, lambda x: x / c, lambda x: c * x),
            (si.J, inv_m_spec, lambda x: x / hc, lambda x: hc * x),
            (inv_m_spec, inv_m_ang, lambda x: x * two_pi, lambda x: x / two_pi),
            (si.m, inv_m_ang, lambda x: two_pi / x),
            (si.Hz, inv_m_ang, lambda x: two_pi * x / c, lambda x: c * x / two_pi),
            (si.J, inv_m_ang, lambda x: x * two_pi / hc, lambda x: hc * x / two_pi),
        ],
        "spectral",
    )


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
            raise ValueError("If `wav` is specified as a unit, `factor` should be set")
        wav = factor * wav  # Convert to Quantity
    c_Aps = _si.c.to_value(si.AA / si.s)  # Angstrom/s
    h_cgs = _si.h.cgs.value  # erg * s
    hc = c_Aps * h_cgs

    # flux density
    f_la = cgs.erg / si.angstrom / si.cm**2 / si.s
    f_nu = cgs.erg / si.Hz / si.cm**2 / si.s
    nu_f_nu = cgs.erg / si.cm**2 / si.s
    la_f_la = nu_f_nu
    phot_f_la = astrophys.photon / (si.cm**2 * si.s * si.AA)
    phot_f_nu = astrophys.photon / (si.cm**2 * si.s * si.Hz)
    la_phot_f_la = astrophys.photon / (si.cm**2 * si.s)

    # luminosity density
    L_nu = cgs.erg / si.s / si.Hz
    L_la = cgs.erg / si.s / si.angstrom
    nu_L_nu = cgs.erg / si.s
    la_L_la = nu_L_nu
    phot_L_la = astrophys.photon / (si.s * si.AA)
    phot_L_nu = astrophys.photon / (si.s * si.Hz)

    # surface brightness (flux equiv)
    S_la = cgs.erg / si.angstrom / si.cm**2 / si.s / si.sr
    S_nu = cgs.erg / si.Hz / si.cm**2 / si.s / si.sr
    nu_S_nu = cgs.erg / si.cm**2 / si.s / si.sr
    la_S_la = nu_S_nu
    phot_S_la = astrophys.photon / (si.cm**2 * si.s * si.AA * si.sr)
    phot_S_nu = astrophys.photon / (si.cm**2 * si.s * si.Hz * si.sr)

    # surface brightness (luminosity equiv)
    SL_nu = cgs.erg / si.s / si.Hz / si.sr
    SL_la = cgs.erg / si.s / si.angstrom / si.sr
    nu_SL_nu = cgs.erg / si.s / si.sr
    la_SL_la = nu_SL_nu
    phot_SL_la = astrophys.photon / (si.s * si.AA * si.sr)
    phot_SL_nu = astrophys.photon / (si.s * si.Hz * si.sr)

    def f_la_to_f_nu(x):
        return x * (wav.to_value(si.AA, spectral()) ** 2 / c_Aps)

    def f_la_from_f_nu(x):
        return x / (wav.to_value(si.AA, spectral()) ** 2 / c_Aps)

    def f_nu_to_nu_f_nu(x):
        return x * wav.to_value(si.Hz, spectral())

    def f_nu_from_nu_f_nu(x):
        return x / wav.to_value(si.Hz, spectral())

    def f_la_to_la_f_la(x):
        return x * wav.to_value(si.AA, spectral())

    def f_la_from_la_f_la(x):
        return x / wav.to_value(si.AA, spectral())

    def phot_f_la_to_f_la(x):
        return hc * x / wav.to_value(si.AA, spectral())

    def phot_f_la_from_f_la(x):
        return x * wav.to_value(si.AA, spectral()) / hc

    def phot_f_la_to_f_nu(x):
        return h_cgs * x * wav.to_value(si.AA, spectral())

    def phot_f_la_from_f_nu(x):
        return x / (wav.to_value(si.AA, spectral()) * h_cgs)

    def phot_f_la_to_phot_f_nu(x):
        return x * wav.to_value(si.AA, spectral()) ** 2 / c_Aps

    def phot_f_la_from_phot_f_nu(x):
        return c_Aps * x / wav.to_value(si.AA, spectral()) ** 2

    phot_f_nu_to_f_nu = phot_f_la_to_f_la
    phot_f_nu_from_f_nu = phot_f_la_from_f_la

    def phot_f_nu_to_f_la(x):
        return x * hc * c_Aps / wav.to_value(si.AA, spectral()) ** 3

    def phot_f_nu_from_f_la(x):
        return x * wav.to_value(si.AA, spectral()) ** 3 / (hc * c_Aps)

    # for luminosity density
    L_nu_to_nu_L_nu = f_nu_to_nu_f_nu
    L_nu_from_nu_L_nu = f_nu_from_nu_f_nu
    L_la_to_la_L_la = f_la_to_la_f_la
    L_la_from_la_L_la = f_la_from_la_f_la

    phot_L_la_to_L_la = phot_f_la_to_f_la
    phot_L_la_from_L_la = phot_f_la_from_f_la
    phot_L_la_to_L_nu = phot_f_la_to_f_nu
    phot_L_la_from_L_nu = phot_f_la_from_f_nu
    phot_L_la_to_phot_L_nu = phot_f_la_to_phot_f_nu
    phot_L_la_from_phot_L_nu = phot_f_la_from_phot_f_nu
    phot_L_nu_to_L_nu = phot_f_nu_to_f_nu
    phot_L_nu_from_L_nu = phot_f_nu_from_f_nu
    phot_L_nu_to_L_la = phot_f_nu_to_f_la
    phot_L_nu_from_L_la = phot_f_nu_from_f_la

    return Equivalency(
        [
            # flux
            (f_la, f_nu, f_la_to_f_nu, f_la_from_f_nu),
            (f_nu, nu_f_nu, f_nu_to_nu_f_nu, f_nu_from_nu_f_nu),
            (f_la, la_f_la, f_la_to_la_f_la, f_la_from_la_f_la),
            (phot_f_la, f_la, phot_f_la_to_f_la, phot_f_la_from_f_la),
            (phot_f_la, f_nu, phot_f_la_to_f_nu, phot_f_la_from_f_nu),
            (phot_f_la, phot_f_nu, phot_f_la_to_phot_f_nu, phot_f_la_from_phot_f_nu),
            (phot_f_nu, f_nu, phot_f_nu_to_f_nu, phot_f_nu_from_f_nu),
            (phot_f_nu, f_la, phot_f_nu_to_f_la, phot_f_nu_from_f_la),
            # integrated flux
            (la_phot_f_la, la_f_la, phot_f_la_to_f_la, phot_f_la_from_f_la),
            # luminosity
            (L_la, L_nu, f_la_to_f_nu, f_la_from_f_nu),
            (L_nu, nu_L_nu, L_nu_to_nu_L_nu, L_nu_from_nu_L_nu),
            (L_la, la_L_la, L_la_to_la_L_la, L_la_from_la_L_la),
            (phot_L_la, L_la, phot_L_la_to_L_la, phot_L_la_from_L_la),
            (phot_L_la, L_nu, phot_L_la_to_L_nu, phot_L_la_from_L_nu),
            (phot_L_la, phot_L_nu, phot_L_la_to_phot_L_nu, phot_L_la_from_phot_L_nu),
            (phot_L_nu, L_nu, phot_L_nu_to_L_nu, phot_L_nu_from_L_nu),
            (phot_L_nu, L_la, phot_L_nu_to_L_la, phot_L_nu_from_L_la),
            # surface brightness (flux equiv)
            (S_la, S_nu, f_la_to_f_nu, f_la_from_f_nu),
            (S_nu, nu_S_nu, f_nu_to_nu_f_nu, f_nu_from_nu_f_nu),
            (S_la, la_S_la, f_la_to_la_f_la, f_la_from_la_f_la),
            (phot_S_la, S_la, phot_f_la_to_f_la, phot_f_la_from_f_la),
            (phot_S_la, S_nu, phot_f_la_to_f_nu, phot_f_la_from_f_nu),
            (phot_S_la, phot_S_nu, phot_f_la_to_phot_f_nu, phot_f_la_from_phot_f_nu),
            (phot_S_nu, S_nu, phot_f_nu_to_f_nu, phot_f_nu_from_f_nu),
            (phot_S_nu, S_la, phot_f_nu_to_f_la, phot_f_nu_from_f_la),
            # surface brightness (luminosity equiv)
            (SL_la, SL_nu, f_la_to_f_nu, f_la_from_f_nu),
            (SL_nu, nu_SL_nu, L_nu_to_nu_L_nu, L_nu_from_nu_L_nu),
            (SL_la, la_SL_la, L_la_to_la_L_la, L_la_from_la_L_la),
            (phot_SL_la, SL_la, phot_L_la_to_L_la, phot_L_la_from_L_la),
            (phot_SL_la, SL_nu, phot_L_la_to_L_nu, phot_L_la_from_L_nu),
            (phot_SL_la, phot_SL_nu, phot_L_la_to_phot_L_nu, phot_L_la_from_phot_L_nu),
            (phot_SL_nu, SL_nu, phot_L_nu_to_L_nu, phot_L_nu_from_L_nu),
            (phot_SL_nu, SL_la, phot_L_nu_to_L_la, phot_L_nu_from_L_la),
        ],
        "spectral_density",
        {"wav": wav, "factor": factor},
    )


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
    `NRAO site defining the conventions <https://www.gb.nrao.edu/~fghigo/gbtdoc/doppler.html>`_

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

    ckms = _si.c.to_value("km/s")

    def to_vel_freq(x):
        restfreq = rest.to_value(si.Hz, equivalencies=spectral())
        return (restfreq - x) / (restfreq) * ckms

    def from_vel_freq(x):
        restfreq = rest.to_value(si.Hz, equivalencies=spectral())
        voverc = x / ckms
        return restfreq * (1 - voverc)

    def to_vel_wav(x):
        restwav = rest.to_value(si.AA, spectral())
        return (x - restwav) / (x) * ckms

    def from_vel_wav(x):
        restwav = rest.to_value(si.AA, spectral())
        return restwav * ckms / (ckms - x)

    def to_vel_en(x):
        resten = rest.to_value(si.eV, equivalencies=spectral())
        return (resten - x) / (resten) * ckms

    def from_vel_en(x):
        resten = rest.to_value(si.eV, equivalencies=spectral())
        voverc = x / ckms
        return resten * (1 - voverc)

    return Equivalency(
        [
            (si.Hz, si.km / si.s, to_vel_freq, from_vel_freq),
            (si.AA, si.km / si.s, to_vel_wav, from_vel_wav),
            (si.eV, si.km / si.s, to_vel_en, from_vel_en),
        ],
        "doppler_radio",
        {"rest": rest},
    )


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
    `NRAO site defining the conventions <https://www.gb.nrao.edu/~fghigo/gbtdoc/doppler.html>`_

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

    ckms = _si.c.to_value("km/s")

    def to_vel_freq(x):
        restfreq = rest.to_value(si.Hz, equivalencies=spectral())
        return ckms * (restfreq - x) / x

    def from_vel_freq(x):
        restfreq = rest.to_value(si.Hz, equivalencies=spectral())
        voverc = x / ckms
        return restfreq / (1 + voverc)

    def to_vel_wav(x):
        restwav = rest.to_value(si.AA, spectral())
        return ckms * (x / restwav - 1)

    def from_vel_wav(x):
        restwav = rest.to_value(si.AA, spectral())
        voverc = x / ckms
        return restwav * (1 + voverc)

    def to_vel_en(x):
        resten = rest.to_value(si.eV, equivalencies=spectral())
        return ckms * (resten - x) / x

    def from_vel_en(x):
        resten = rest.to_value(si.eV, equivalencies=spectral())
        voverc = x / ckms
        return resten / (1 + voverc)

    return Equivalency(
        [
            (si.Hz, si.km / si.s, to_vel_freq, from_vel_freq),
            (si.AA, si.km / si.s, to_vel_wav, from_vel_wav),
            (si.eV, si.km / si.s, to_vel_en, from_vel_en),
        ],
        "doppler_optical",
        {"rest": rest},
    )


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
    `NRAO site defining the conventions <https://www.gb.nrao.edu/~fghigo/gbtdoc/doppler.html>`_

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

    ckms = _si.c.to_value("km/s")

    def to_vel_freq(x):
        restfreq = rest.to_value(si.Hz, equivalencies=spectral())
        return (restfreq**2 - x**2) / (restfreq**2 + x**2) * ckms

    def from_vel_freq(x):
        restfreq = rest.to_value(si.Hz, equivalencies=spectral())
        voverc = x / ckms
        return restfreq * ((1 - voverc) / (1 + (voverc))) ** 0.5

    def to_vel_wav(x):
        restwav = rest.to_value(si.AA, spectral())
        return (x**2 - restwav**2) / (restwav**2 + x**2) * ckms

    def from_vel_wav(x):
        restwav = rest.to_value(si.AA, spectral())
        voverc = x / ckms
        return restwav * ((1 + voverc) / (1 - voverc)) ** 0.5

    def to_vel_en(x):
        resten = rest.to_value(si.eV, spectral())
        return (resten**2 - x**2) / (resten**2 + x**2) * ckms

    def from_vel_en(x):
        resten = rest.to_value(si.eV, spectral())
        voverc = x / ckms
        return resten * ((1 - voverc) / (1 + (voverc))) ** 0.5

    return Equivalency(
        [
            (si.Hz, si.km / si.s, to_vel_freq, from_vel_freq),
            (si.AA, si.km / si.s, to_vel_wav, from_vel_wav),
            (si.eV, si.km / si.s, to_vel_en, from_vel_en),
        ],
        "doppler_relativistic",
        {"rest": rest},
    )


def doppler_redshift():
    """
    Returns the equivalence between Doppler redshift (unitless) and radial velocity.

    .. note::

        This equivalency is not compatible with cosmological
        redshift in `astropy.cosmology.units`.

    """
    rv_unit = si.km / si.s
    C_KMS = _si.c.to_value(rv_unit)

    def convert_z_to_rv(z):
        zponesq = (1 + z) ** 2
        return C_KMS * (zponesq - 1) / (zponesq + 1)

    def convert_rv_to_z(rv):
        beta = rv / C_KMS
        return np.sqrt((1 + beta) / (1 - beta)) - 1

    return Equivalency(
        [(dimensionless_unscaled, rv_unit, convert_z_to_rv, convert_rv_to_z)],
        "doppler_redshift",
    )


def molar_mass_amu():
    """
    Returns the equivalence between amu and molar mass.
    """
    return Equivalency([(si.g / si.mol, misc.u)], "molar_mass_amu")


def mass_energy():
    """
    Returns a list of equivalence pairs that handle the conversion
    between mass and energy.
    """
    c2 = _si.c.value**2
    return Equivalency(
        [
            (si.kg, si.J, lambda x: x * c2, lambda x: x / c2),
            (si.kg / si.m**2, si.J / si.m**2, lambda x: x * c2, lambda x: x / c2),
            (si.kg / si.m**3, si.J / si.m**3, lambda x: x * c2, lambda x: x / c2),
            (si.kg / si.s, si.J / si.s, lambda x: x * c2, lambda x: x / c2),
        ],
        "mass_energy",
    )


def brightness_temperature(frequency, beam_area=None):
    r"""
    Defines the conversion between Jy/sr and "brightness temperature",
    :math:`T_B`, in Kelvins.  The brightness temperature is a unit very
    commonly used in radio astronomy.  See, e.g., "Tools of Radio Astronomy"
    (Wilson 2009) eqn 8.16 and eqn 8.19 (these pages are available on `google
    books
    <https://books.google.com/books?id=9KHw6R8rQEMC&pg=PA179&source=gbs_toc_r&cad=4#v=onepage&q&f=false>`__).

    :math:`T_B \equiv S_\nu / \left(2 k \nu^2 / c^2 \right)`

    If the input is in Jy/beam or Jy (assuming it came from a single beam), the
    beam area is essential for this computation: the brightness temperature is
    inversely proportional to the beam area.

    Parameters
    ----------
    frequency : `~astropy.units.Quantity`
        The observed ``spectral`` equivalent `~astropy.units.Unit` (e.g.,
        frequency or wavelength).  The variable is named 'frequency' because it
        is more commonly used in radio astronomy.
        BACKWARD COMPATIBILITY NOTE: previous versions of the brightness
        temperature equivalency used the keyword ``disp``, which is no longer
        supported.
    beam_area : `~astropy.units.Quantity` ['solid angle']
        Beam area in angular units, i.e. steradian equivalent

    Examples
    --------
    Arecibo C-band beam::

        >>> import numpy as np
        >>> from astropy import units as u
        >>> beam_sigma = 50*u.arcsec
        >>> beam_area = 2*np.pi*(beam_sigma)**2
        >>> freq = 5*u.GHz
        >>> equiv = u.brightness_temperature(freq)
        >>> (1*u.Jy/beam_area).to(u.K, equivalencies=equiv)  # doctest: +FLOAT_CMP
        <Quantity 3.526295144567176 K>

    VLA synthetic beam::

        >>> bmaj = 15*u.arcsec
        >>> bmin = 15*u.arcsec
        >>> fwhm_to_sigma = 1./(8*np.log(2))**0.5
        >>> beam_area = 2.*np.pi*(bmaj*bmin*fwhm_to_sigma**2)
        >>> freq = 5*u.GHz
        >>> equiv = u.brightness_temperature(freq)
        >>> (u.Jy/beam_area).to(u.K, equivalencies=equiv)  # doctest: +FLOAT_CMP
        <Quantity 217.2658703625732 K>

    Any generic surface brightness:

        >>> surf_brightness = 1e6*u.MJy/u.sr
        >>> surf_brightness.to(u.K, equivalencies=u.brightness_temperature(500*u.GHz)) # doctest: +FLOAT_CMP
        <Quantity 130.1931904778803 K>
    """
    nu = frequency.to(si.GHz, spectral())
    factor_Jy = (2 * _si.k_B * si.K * nu**2 / _si.c**2).to(astrophys.Jy).value
    factor_K = (astrophys.Jy / (2 * _si.k_B * nu**2 / _si.c**2)).to(si.K).value

    if beam_area is not None:
        beam = beam_area.to_value(si.sr)

        def convert_Jy_to_K(x_jybm):
            return x_jybm / beam / factor_Jy

        def convert_K_to_Jy(x_K):
            return x_K * beam / factor_K

        return Equivalency(
            [
                (astrophys.Jy, si.K, convert_Jy_to_K, convert_K_to_Jy),
                (astrophys.Jy / astrophys.beam, si.K, convert_Jy_to_K, convert_K_to_Jy),
            ],
            "brightness_temperature",
            {"frequency": frequency, "beam_area": beam_area},
        )
    else:

        def convert_JySr_to_K(x_jysr):
            return x_jysr / factor_Jy

        def convert_K_to_JySr(x_K):
            return x_K / factor_K  # multiplied by 1x for 1 steradian

        return Equivalency(
            [(astrophys.Jy / si.sr, si.K, convert_JySr_to_K, convert_K_to_JySr)],
            "brightness_temperature",
            {"frequency": frequency, "beam_area": beam_area},
        )


def beam_angular_area(beam_area):
    """
    Convert between the ``beam`` unit, which is commonly used to express the area
    of a radio telescope resolution element, and an area on the sky.
    This equivalency also supports direct conversion between ``Jy/beam`` and
    ``Jy/steradian`` units, since that is a common operation.

    Parameters
    ----------
    beam_area : unit-like
        The area of the beam in angular area units (e.g., steradians)
        Must have angular area equivalent units.
    """
    return Equivalency(
        [
            (astrophys.beam, Unit(beam_area)),
            (astrophys.beam**-1, Unit(beam_area) ** -1),
            (astrophys.Jy / astrophys.beam, astrophys.Jy / Unit(beam_area)),
        ],
        "beam_angular_area",
        {"beam_area": beam_area},
    )


def thermodynamic_temperature(frequency, T_cmb=None):
    r"""Defines the conversion between Jy/sr and "thermodynamic temperature",
    :math:`T_{CMB}`, in Kelvins.  The thermodynamic temperature is a unit very
    commonly used in cosmology. See eqn 8 in [1].

    :math:`K_{CMB} \equiv I_\nu / \left(2 k \nu^2 / c^2  f(\nu) \right)`

    with :math:`f(\nu) = \frac{ x^2 e^x}{(e^x - 1 )^2}`
    where :math:`x = h \nu / k T`

    Parameters
    ----------
    frequency : `~astropy.units.Quantity`
        The observed `spectral` equivalent `~astropy.units.Unit` (e.g.,
        frequency or wavelength). Must have spectral units.
    T_cmb :  `~astropy.units.Quantity` ['temperature'] or None
        The CMB temperature at z=0.  If `None`, the default cosmology will be
        used to get this temperature. Must have units of temperature.

    Notes
    -----
    For broad band receivers, this conversion do not hold
    as it highly depends on the frequency

    References
    ----------
    .. [1] Planck 2013 results. IX. HFI spectral response
       https://arxiv.org/abs/1303.5070

    Examples
    --------
    Planck HFI 143 GHz::

        >>> from astropy import units as u
        >>> from astropy.cosmology import Planck15
        >>> freq = 143 * u.GHz
        >>> equiv = u.thermodynamic_temperature(freq, Planck15.Tcmb0)
        >>> (1. * u.mK).to(u.MJy / u.sr, equivalencies=equiv)  # doctest: +FLOAT_CMP
        <Quantity 0.37993172 MJy / sr>

    """
    nu = frequency.to(si.GHz, spectral())

    if T_cmb is None:
        from astropy.cosmology import default_cosmology

        T_cmb = default_cosmology.get().Tcmb0

    def f(nu, T_cmb=T_cmb):
        x = _si.h * nu / _si.k_B / T_cmb
        return x**2 * np.exp(x) / np.expm1(x) ** 2

    def convert_Jy_to_K(x_jybm):
        factor = (f(nu) * 2 * _si.k_B * si.K * nu**2 / _si.c**2).to_value(astrophys.Jy)
        return x_jybm / factor

    def convert_K_to_Jy(x_K):
        factor = (astrophys.Jy / (f(nu) * 2 * _si.k_B * nu**2 / _si.c**2)).to_value(
            si.K
        )
        return x_K / factor

    return Equivalency(
        [(astrophys.Jy / si.sr, si.K, convert_Jy_to_K, convert_K_to_Jy)],
        "thermodynamic_temperature",
        {"frequency": frequency, "T_cmb": T_cmb},
    )


def temperature():
    """Convert between Kelvin, Celsius, Rankine and Fahrenheit here because
    Unit and CompositeUnit cannot do addition or subtraction properly.
    """
    from .imperial import deg_F as F
    from .imperial import deg_R as R

    K = si.K
    C = si.deg_C

    return Equivalency(
        [
            (K, C, lambda x: x - 273.15, lambda x: x + 273.15),
            (C, F, lambda x: x * 1.8 + 32.0, lambda x: (x - 32.0) / 1.8),
            (K, F, lambda x: x * 1.8 - 459.67, lambda x: (x + 459.67) / 1.8),
            (R, F, lambda x: x - 459.67, lambda x: x + 459.67),
            (R, C, lambda x: (x - 491.67) * (5 / 9), lambda x: x * 1.8 + 491.67),
            (R, K, lambda x: x * (5 / 9), lambda x: x * 1.8),
        ],
        "temperature",
    )


def temperature_energy():
    """Convert between Kelvin and keV(eV) to an equivalent amount."""
    e = _si.e.value
    k_B = _si.k_B.value
    return Equivalency(
        [(si.K, si.eV, lambda x: x / (e / k_B), lambda x: x * (e / k_B))],
        "temperature_energy",
    )


def assert_is_spectral_unit(value):
    try:
        value.to(si.Hz, spectral())
    except (AttributeError, UnitsError) as ex:
        raise UnitsError(
            "The 'rest' value must be a spectral equivalent "
            "(frequency, wavelength, or energy)."
        )


def pixel_scale(pixscale):
    """
    Convert between pixel distances (in units of ``pix``) and other units,
    given a particular ``pixscale``.

    Parameters
    ----------
    pixscale : `~astropy.units.Quantity`
        The pixel scale either in units of <unit>/pixel or pixel/<unit>.
    """
    decomposed = pixscale.unit.decompose()
    dimensions = dict(zip(decomposed.bases, decomposed.powers))
    pix_power = dimensions.get(misc.pix, 0)

    if pix_power == -1:
        physical_unit = Unit(pixscale * misc.pix)
    elif pix_power == 1:
        physical_unit = Unit(misc.pix / pixscale)
    else:
        raise UnitsError(
            "The pixel scale unit must have pixel dimensionality of 1 or -1."
        )

    return Equivalency(
        [(misc.pix, physical_unit)], "pixel_scale", {"pixscale": pixscale}
    )


def plate_scale(platescale):
    """
    Convert between lengths (to be interpreted as lengths in the focal plane)
    and angular units with a specified ``platescale``.

    Parameters
    ----------
    platescale : `~astropy.units.Quantity`
        The pixel scale either in units of distance/pixel or distance/angle.
    """
    if platescale.unit.is_equivalent(si.arcsec / si.m):
        platescale_val = platescale.to_value(si.radian / si.m)
    elif platescale.unit.is_equivalent(si.m / si.arcsec):
        platescale_val = (1 / platescale).to_value(si.radian / si.m)
    else:
        raise UnitsError("The pixel scale must be in angle/distance or distance/angle")

    return Equivalency(
        [(si.m, si.radian, lambda d: d * platescale_val, lambda a: a / platescale_val)],
        "plate_scale",
        {"platescale": platescale},
    )
