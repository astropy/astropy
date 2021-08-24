# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Cosmological units and equivalencies.
"""  # (newline needed for unit summary)

import astropy.units as u
from astropy.units.utils import generate_unit_summary as _generate_unit_summary

__all__ = ["littleh", "with_H0"]

_ns = globals()

###############################################################################
# Cosmological Units

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


# =============================================================================
# DOCSTRING

# This generates a docstring for this module that describes all of the
# standard units defined here.
if __doc__ is not None:
    __doc__ += _generate_unit_summary(_ns)
