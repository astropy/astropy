# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Quantity helpers for the ERFA ufuncs."""


from astropy.units.core import UnitsError, UnitTypeError, dimensionless_unscaled
from . import UFUNC_HELPERS
from .helpers import get_converter, helper_invariant, helper_multiplication


erfa_ufuncs = ('s2c', 's2p', 'c2s', 'p2s', 'pm', 'pdp', 'pxp', 'rxp')


def helper_s2c(f, unit1, unit2):
    from astropy.units.si import radian
    try:
        return [get_converter(unit1, radian),
                get_converter(unit2, radian)], dimensionless_unscaled
    except UnitsError:
        raise UnitTypeError("Can only apply '{}' function to "
                            "quantities with angle units"
                            .format(f.__name__))


def helper_s2p(f, unit1, unit2, unit3):
    from astropy.units.si import radian
    try:
        return [get_converter(unit1, radian),
                get_converter(unit2, radian), None], unit3
    except UnitsError:
        raise UnitTypeError("Can only apply '{}' function to "
                            "quantities with angle units"
                            .format(f.__name__))


def helper_c2s(f, unit1):
    from astropy.units.si import radian
    return [None], (radian, radian)


def helper_p2s(f, unit1):
    from astropy.units.si import radian
    return [None], (radian, radian, unit1)


def get_erfa_helpers():
    from astropy._erfa import ufunc as erfa_ufunc

    ERFA_HELPERS = {}
    ERFA_HELPERS[erfa_ufunc.s2c] = helper_s2c
    ERFA_HELPERS[erfa_ufunc.s2p] = helper_s2p
    ERFA_HELPERS[erfa_ufunc.c2s] = helper_c2s
    ERFA_HELPERS[erfa_ufunc.p2s] = helper_p2s
    ERFA_HELPERS[erfa_ufunc.pm] = helper_invariant
    ERFA_HELPERS[erfa_ufunc.pdp] = helper_multiplication
    ERFA_HELPERS[erfa_ufunc.pxp] = helper_multiplication
    ERFA_HELPERS[erfa_ufunc.rxp] = helper_multiplication
    return ERFA_HELPERS


UFUNC_HELPERS.register_module('astropy._erfa.ufunc', erfa_ufuncs,
                              get_erfa_helpers)
