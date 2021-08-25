# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Cosmological units and equivalencies.
"""  # (newline needed for unit summary)

import astropy.units as u
from astropy.units.utils import generate_unit_summary as _generate_unit_summary

__all__ = ["littleh", "redshift",
           # equivalencies
           "dimensionless_redshift", "with_redshift", "with_H0"]

_ns = globals()


###############################################################################
# Cosmological Units

# This is not formally a unit, but is used in that way in many contexts, and
# an appropriate equivalency is only possible if it's treated as a unit.
redshift = u.def_unit(['redshift'], prefixes=False, namespace=_ns,
                      doc="Cosmological redshift.", format={'latex': r''})

# This is not formally a unit, but is used in that way in many contexts, and
# an appropriate equivalency is only possible if it's treated as a unit (see
# https://arxiv.org/pdf/1308.4150.pdf for more)
# Also note that h or h100 or h_100 would be a better name, but they either
# conflict or have numbers in them, which is disallowed
littleh = u.def_unit(['littleh'], namespace=_ns, prefixes=False,
                     doc='Reduced/"dimensionless" Hubble constant',
                     format={'latex': r'h_{100}'})


###############################################################################
# Equivalencies


def dimensionless_redshift():
    """Allow redshift to be 1-to-1 equivalent to dimensionless.

    It is special compared to other equivalency pairs in that it
    allows this independent of the power to which the angle is raised,
    and independent of whether it is part of a more complicated unit.
    """
    return u.Equivalency([(redshift, None)], "dimensionless_redshift")


def with_redshift(cosmology=None, *, Tcmb=True, atzkw=None):
    """Convert quantities between measures of cosmological distance.

    When this equivalency is enabled ALL redshifts are treated as cosmological.
    Care should be taken to not misinterpret a relativistic, gravitational, etc
    redshift as a cosmological one.

    Parameters
    ----------
    cosmology : `~astropy.cosmology.Cosmology`, str, or None, optional
        A cosmology realization or built-in cosmology's name (e.g. 'Planck18').
        If None, will use the default cosmology
        (controlled by :class:`~astropy.cosmology.default_cosmology`).
    Tcmb : bool (optional, keyword-only)
        Whether to create a CMB temperature <-> redshift equivalency.
    atzkw : dict or None (optional, keyword-only)
        keyword arguments for :func:`~astropy.cosmology.z_at_value`

    Returns
    -------
    `~astropy.units.equivalencies.Equivalency`
        With equivalencies between redshift and (temperature, distance).
    """
    from astropy.cosmology import default_cosmology, z_at_value

    cosmology = cosmology if cosmology is not None else default_cosmology.get()
    with default_cosmology.set(cosmology):  # if already cosmo, passes through
        cosmology = default_cosmology.get()
    atzkw = atzkw or {}  # None -> {}

    equivs = []  # will append as built

    # -----------
    # CMB Temperature <-> Redshift

    if Tcmb:

        def z_to_Tcmb(z):
            return cosmology.Tcmb(z)

        def Tcmb_to_z(T):
            return z_at_value(cosmology.Tcmb, T << u.K, **atzkw) << redshift

        equivs.append((redshift, u.K, z_to_Tcmb, Tcmb_to_z))

    # -----------

    return u.Equivalency(equivs, "with_redshift",
                         {'cosmology': cosmology, 'Tcmb': Tcmb})


def with_H0(H0=None):
    """
    Convert between quantities with little-h and the equivalent physical units.

    Parameters
    ----------
    H0 : None or `~astropy.units.Quantity` ['frequency']
        The value of the Hubble constant to assume. If a
        `~astropy.units.Quantity`, will assume the quantity *is* ``H0``. If
        `None` (default), use the ``H0`` attribute from
        :mod:`~astropy.cosmology.default_cosmology`.

    References
    ----------
    For an illuminating discussion on why you may or may not want to use
    little-h at all, see https://arxiv.org/pdf/1308.4150.pdf
    """
    if H0 is None:
        from .realizations import default_cosmology
        H0 = default_cosmology.get().H0

    h100_val_unit = u.Unit(100/(H0.to_value((u.km/u.s)/u.Mpc)) * littleh)

    return u.Equivalency([(h100_val_unit, None)], "with_H0", kwargs={"H0": H0})


# ===================================================================
# Enable the set of default equivalencies.
# If the cosmology package is imported, this is added to the list astropy-wide.

u.add_enabled_equivalencies(dimensionless_redshift())


# =============================================================================
# DOCSTRING

# This generates a docstring for this module that describes all of the
# standard units defined here.
if __doc__ is not None:
    __doc__ += _generate_unit_summary(_ns)
