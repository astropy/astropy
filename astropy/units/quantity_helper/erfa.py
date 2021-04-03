# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Quantity helpers for the ERFA ufuncs."""
# Tests for these are in coordinates, not in units.

from erfa import ufunc as erfa_ufunc

from astropy.units.core import UnitsError, UnitTypeError, dimensionless_unscaled
from . import UFUNC_HELPERS
from .helpers import (get_converter, helper_invariant, helper_multiplication,
                      helper_twoarg_invariant, _d)
# TODO: cannot import StructuredUnit here; maybe try to rewrite so that is possible?

erfa_ufuncs = ('s2c', 's2p', 'c2s', 'p2s', 'pm', 'pdp', 'pxp', 'rxp',
               'cpv', 'p2pv', 'pv2p', 'pv2s', 'pvdpv', 'pvm', 'pvmpv', 'pvppv',
               'pvstar', 'pvtob', 'pvu', 'pvup', 'pvxpv', 'rxpv', 's2pv', 's2xpv',
               'starpv', 'sxpv', 'trxpv', 'gd2gc', 'gc2gd')


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


def helper_gc2gd(f, nounit, unit1):
    from astropy.units.si import m, radian
    if nounit is not None:
        raise UnitTypeError("ellipsoid cannot be a quantity.")
    try:
        return [None, get_converter(unit1, m)], (radian, radian, m, None)
    except UnitsError:
        raise UnitTypeError("Can only apply '{}' function to "
                            "quantities with length units"
                            .format(f.__name__))


def helper_gd2gc(f, nounit, unit1, unit2, unit3):
    from astropy.units.si import m, radian
    if nounit is not None:
        raise UnitTypeError("ellipsoid cannot be a quantity.")
    try:
        return [None,
                get_converter(unit1, radian),
                get_converter(unit2, radian),
                get_converter(unit3, m)], (m, None)
    except UnitsError:
        raise UnitTypeError("Can only apply '{}' function to lon, lat "
                            "with angle and height with length units"
                            .format(f.__name__))


def helper_p2pv(f, unit1):
    from astropy.units.structured import StructuredUnit
    from astropy.units.si import m, s
    if isinstance(unit1, StructuredUnit):
        raise UnitTypeError("p vector unit cannot be a structured unit.")
    return [None], StructuredUnit((unit1, m / s))


def helper_pv2p(f, unit1):
    try:
        return [None], unit1['p']
    except Exception as exc:
        raise UnitTypeError("pv vector should have a structured unit.") from exc


def helper_pv2s(f, unit1):
    from astropy.units.si import radian
    try:
        p_unit, v_unit = unit1['p'], unit1['v']
    except Exception as exc:
        raise UnitTypeError("pv vector should have a structured unit.") from exc
    else:
        ang_unit = radian * v_unit / p_unit
        return [None], (radian, radian, p_unit, ang_unit, ang_unit, v_unit)


def helper_s2pv(f, unit_theta, unit_phi, unit_r, unit_td, unit_pd, unit_rd):
    from astropy.units.si import radian
    from astropy.units.structured import StructuredUnit
    time_unit = unit_r / unit_rd
    return [get_converter(unit_theta, radian),
            get_converter(unit_phi, radian),
            None,
            get_converter(unit_td, radian / time_unit),
            get_converter(unit_pd, radian / time_unit),
            None], StructuredUnit((unit_r, unit_rd))


def helper_pv_multiplication(f, unit1, unit2):
    from astropy.units.structured import StructuredUnit
    try:
        result_unit = StructuredUnit((unit1['p'] * unit2['p'],
                                      unit1['v'] * unit2['p']))
    except Exception as exc:
        raise UnitTypeError("pv vectors should have a structured unit.") from exc
    converter = get_converter(unit2, StructuredUnit(
        (unit2['p'], unit1['v'] * unit2['p'] / unit1['p'])))
    return [None, converter], result_unit


def helper_pvm(f, unit1):
    try:
        return [None], (unit1['p'], unit1['v'])
    except Exception as exc:
        raise UnitTypeError("pv vector should have a structured unit.") from exc


def helper_pvstar(f, unit1):
    from astropy.units.structured import StructuredUnit
    from astropy.units.astrophys import AU
    from astropy.units.si import km, s, radian, day, year, arcsec

    return [get_converter(unit1, StructuredUnit((AU, AU/day)))], (
        radian, radian, radian / year, radian / year, arcsec, km / s, None)


def helper_starpv(f, unit_ra, unit_dec, unit_pmr, unit_pmd,
                  unit_px, unit_rv):
    from astropy.units.structured import StructuredUnit
    from astropy.units.si import km, s, day, year, radian, arcsec
    from astropy.units.astrophys import AU

    return [get_converter(unit_ra, radian),
            get_converter(unit_dec, radian),
            get_converter(unit_pmr, radian/year),
            get_converter(unit_pmd, radian/year),
            get_converter(unit_px, arcsec),
            get_converter(unit_rv, km/s)], (StructuredUnit((AU, AU/day)), None)


def helper_pvtob(f, unit_elong, unit_phi, unit_hm,
                 unit_xp, unit_yp, unit_sp, unit_theta):
    from astropy.units.structured import StructuredUnit
    from astropy.units.si import m, s, radian

    return [get_converter(unit_elong, radian),
            get_converter(unit_phi, radian),
            get_converter(unit_hm, m),
            get_converter(unit_xp, radian),
            get_converter(unit_yp, radian),
            get_converter(unit_sp, radian),
            get_converter(unit_theta, radian)], StructuredUnit((m, m/s))


def helper_pvu(f, unit_t, unit_pv):
    try:
        unit_p, unit_v = unit_pv['p'], unit_pv['v']
    except Exception as exc:
        raise UnitTypeError("pv vector should have a structured unit.") from exc
    return [get_converter(unit_t, unit_p/unit_v), None], unit_pv


def helper_pvup(f, unit_t, unit_pv):
    try:
        unit_p, unit_v = unit_pv['p'], unit_pv['v']
    except Exception as exc:
        raise UnitTypeError("pv vector should have a structured unit.") from exc
    return [get_converter(unit_t, unit_p/unit_v), None], unit_p


def helper_s2xpv(f, unit1, unit2, unit_pv):
    from astropy.units.structured import StructuredUnit
    try:
        return [None, None, None], StructuredUnit((_d(unit1) * unit_pv['p'],
                                                   _d(unit2) * unit_pv['v']))
    except Exception as exc:
        raise UnitTypeError("pv vector should have a structured unit.") from exc


def get_erfa_helpers():
    ERFA_HELPERS = {}
    ERFA_HELPERS[erfa_ufunc.s2c] = helper_s2c
    ERFA_HELPERS[erfa_ufunc.s2p] = helper_s2p
    ERFA_HELPERS[erfa_ufunc.c2s] = helper_c2s
    ERFA_HELPERS[erfa_ufunc.p2s] = helper_p2s
    ERFA_HELPERS[erfa_ufunc.pm] = helper_invariant
    ERFA_HELPERS[erfa_ufunc.cpv] = helper_invariant
    ERFA_HELPERS[erfa_ufunc.p2pv] = helper_p2pv
    ERFA_HELPERS[erfa_ufunc.pv2p] = helper_pv2p
    ERFA_HELPERS[erfa_ufunc.pv2s] = helper_pv2s
    ERFA_HELPERS[erfa_ufunc.pvdpv] = helper_pv_multiplication
    ERFA_HELPERS[erfa_ufunc.pvxpv] = helper_pv_multiplication
    ERFA_HELPERS[erfa_ufunc.pvm] = helper_pvm
    ERFA_HELPERS[erfa_ufunc.pvmpv] = helper_twoarg_invariant
    ERFA_HELPERS[erfa_ufunc.pvppv] = helper_twoarg_invariant
    ERFA_HELPERS[erfa_ufunc.pvstar] = helper_pvstar
    ERFA_HELPERS[erfa_ufunc.pvtob] = helper_pvtob
    ERFA_HELPERS[erfa_ufunc.pvu] = helper_pvu
    ERFA_HELPERS[erfa_ufunc.pvup] = helper_pvup
    ERFA_HELPERS[erfa_ufunc.pdp] = helper_multiplication
    ERFA_HELPERS[erfa_ufunc.pxp] = helper_multiplication
    ERFA_HELPERS[erfa_ufunc.rxp] = helper_multiplication
    ERFA_HELPERS[erfa_ufunc.rxpv] = helper_multiplication
    ERFA_HELPERS[erfa_ufunc.s2pv] = helper_s2pv
    ERFA_HELPERS[erfa_ufunc.s2xpv] = helper_s2xpv
    ERFA_HELPERS[erfa_ufunc.starpv] = helper_starpv
    ERFA_HELPERS[erfa_ufunc.sxpv] = helper_multiplication
    ERFA_HELPERS[erfa_ufunc.trxpv] = helper_multiplication
    ERFA_HELPERS[erfa_ufunc.gc2gd] = helper_gc2gd
    ERFA_HELPERS[erfa_ufunc.gd2gc] = helper_gd2gc
    return ERFA_HELPERS


UFUNC_HELPERS.register_module('erfa.ufunc', erfa_ufuncs,
                              get_erfa_helpers)
