# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Quantity helpers for the ERFA ufuncs."""
# Tests for these are in coordinates, not in units.

from erfa import dt_eraASTROM, dt_eraLDBODY, dt_pv
from erfa import ufunc as erfa_ufunc

from astropy.units.core import UnitsError, UnitTypeError, dimensionless_unscaled
from astropy.units.structured import StructuredUnit

from . import UFUNC_HELPERS
from .helpers import (
    _d,
    get_converter,
    helper_invariant,
    helper_multiplication,
    helper_twoarg_invariant,
)

erfa_ufuncs = (
    "s2c", "s2p", "c2s", "p2s", "pm", "pdp", "pxp", "rxp", "cpv", "p2pv", "pv2p",
    "pv2s", "pvdpv", "pvm", "pvmpv", "pvppv", "pvstar", "pvtob", "pvu", "pvup",
    "pvxpv", "rxpv", "s2pv", "s2xpv", "starpv", "sxpv", "trxpv", "gd2gc", "gd2gce",
    "gc2gd", "gc2gde", "ldn", "aper", "apio", "atciq", "atciqn", "atciqz", "aticq",
    "atioq", "atoiq",
)  # fmt: skip


def has_matching_structure(unit, dtype):
    dtype_fields = dtype.fields
    if dtype_fields:
        return (
            isinstance(unit, StructuredUnit)
            and len(unit) == len(dtype_fields)
            and all(
                has_matching_structure(u, df_v[0])
                for (u, df_v) in zip(unit.values(), dtype_fields.values())
            )
        )
    else:
        return not isinstance(unit, StructuredUnit)


def check_structured_unit(unit, dtype):
    if not has_matching_structure(unit, dtype):
        msg = {dt_pv: "pv", dt_eraLDBODY: "ldbody", dt_eraASTROM: "astrom"}.get(
            dtype, "function"
        )
        raise UnitTypeError(f"{msg} input needs unit matching dtype={dtype}.")


def helper_s2c(f, unit1, unit2):
    from astropy.units.si import radian

    try:
        return [
            get_converter(unit1, radian),
            get_converter(unit2, radian),
        ], dimensionless_unscaled
    except UnitsError:
        raise UnitTypeError(
            f"Can only apply '{f.__name__}' function to quantities with angle units"
        )


def helper_s2p(f, unit1, unit2, unit3):
    from astropy.units.si import radian

    try:
        return [get_converter(unit1, radian), get_converter(unit2, radian), None], unit3
    except UnitsError:
        raise UnitTypeError(
            f"Can only apply '{f.__name__}' function to quantities with angle units"
        )


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
        raise UnitTypeError(
            f"Can only apply '{f.__name__}' function to quantities with length units"
        )


def helper_gc2gde(f, unit_r, unit_flat, unit_xyz):
    from astropy.units.si import m, radian

    return [
        get_converter(unit_r, m),
        get_converter(_d(unit_flat), dimensionless_unscaled),
        get_converter(unit_xyz, m),
    ], (
        radian,
        radian,
        m,
        None,
    )


def helper_gd2gc(f, nounit, unit1, unit2, unit3):
    from astropy.units.si import m, radian

    if nounit is not None:
        raise UnitTypeError("ellipsoid cannot be a quantity.")
    try:
        return [
            None,
            get_converter(unit1, radian),
            get_converter(unit2, radian),
            get_converter(unit3, m),
        ], (m, None)
    except UnitsError:
        raise UnitTypeError(
            f"Can only apply '{f.__name__}' function to lon, lat "
            "with angle and height with length units"
        )


def helper_gd2gce(f, unit_r, unit_flat, unit_long, unit_lat, unit_h):
    from astropy.units.si import m, radian

    return [
        get_converter(unit_r, m),
        get_converter(_d(unit_flat), dimensionless_unscaled),
        get_converter(unit_long, radian),
        get_converter(unit_lat, radian),
        get_converter(unit_h, m),
    ], (m, None)


def helper_p2pv(f, unit1):
    from astropy.units.si import s

    if isinstance(unit1, StructuredUnit):
        raise UnitTypeError("p vector unit cannot be a structured unit.")
    return [None], StructuredUnit((unit1, unit1 / s))


def helper_pv2p(f, unit1):
    check_structured_unit(unit1, dt_pv)
    return [None], unit1[0]


def helper_pv2s(f, unit_pv):
    from astropy.units.si import radian

    check_structured_unit(unit_pv, dt_pv)
    ang_unit = radian * unit_pv[1] / unit_pv[0]
    return [None], (radian, radian, unit_pv[0], ang_unit, ang_unit, unit_pv[1])


def helper_s2pv(f, unit_theta, unit_phi, unit_r, unit_td, unit_pd, unit_rd):
    from astropy.units.si import radian

    time_unit = unit_r / unit_rd
    return [
        get_converter(unit_theta, radian),
        get_converter(unit_phi, radian),
        None,
        get_converter(unit_td, radian / time_unit),
        get_converter(unit_pd, radian / time_unit),
        None,
    ], StructuredUnit((unit_r, unit_rd))


def helper_pv_multiplication(f, unit1, unit2):
    check_structured_unit(unit1, dt_pv)
    check_structured_unit(unit2, dt_pv)
    result_unit = StructuredUnit((unit1[0] * unit2[0], unit1[1] * unit2[0]))
    converter = get_converter(
        unit2, StructuredUnit((unit2[0], unit1[1] * unit2[0] / unit1[0]))
    )
    return [None, converter], result_unit


def helper_pvm(f, unit1):
    check_structured_unit(unit1, dt_pv)
    return [None], (unit1[0], unit1[1])


def helper_pvstar(f, unit1):
    from astropy.units.astrophys import AU
    from astropy.units.si import arcsec, day, km, radian, s, year

    return [get_converter(unit1, StructuredUnit((AU, AU / day)))], (
        radian,
        radian,
        radian / year,
        radian / year,
        arcsec,
        km / s,
        None,
    )


def helper_starpv(f, unit_ra, unit_dec, unit_pmr, unit_pmd, unit_px, unit_rv):
    from astropy.units.astrophys import AU
    from astropy.units.si import arcsec, day, km, radian, s, year

    return [
        get_converter(unit_ra, radian),
        get_converter(unit_dec, radian),
        get_converter(unit_pmr, radian / year),
        get_converter(unit_pmd, radian / year),
        get_converter(unit_px, arcsec),
        get_converter(unit_rv, km / s),
    ], (StructuredUnit((AU, AU / day)), None)


def helper_pvtob(
    f, unit_elong, unit_phi, unit_hm, unit_xp, unit_yp, unit_sp, unit_theta
):
    from astropy.units.si import m, radian, s

    return [
        get_converter(unit_elong, radian),
        get_converter(unit_phi, radian),
        get_converter(unit_hm, m),
        get_converter(unit_xp, radian),
        get_converter(unit_yp, radian),
        get_converter(unit_sp, radian),
        get_converter(unit_theta, radian),
    ], StructuredUnit((m, m / s))


def helper_pvu(f, unit_t, unit_pv):
    check_structured_unit(unit_pv, dt_pv)
    return [get_converter(unit_t, unit_pv[0] / unit_pv[1]), None], unit_pv


def helper_pvup(f, unit_t, unit_pv):
    check_structured_unit(unit_pv, dt_pv)
    return [get_converter(unit_t, unit_pv[0] / unit_pv[1]), None], unit_pv[0]


def helper_s2xpv(f, unit1, unit2, unit_pv):
    check_structured_unit(unit_pv, dt_pv)
    return [None, None, None], StructuredUnit(
        (_d(unit1) * unit_pv[0], _d(unit2) * unit_pv[1])
    )


def ldbody_unit():
    from astropy.units.astrophys import AU, Msun
    from astropy.units.si import day, radian

    return StructuredUnit((Msun, radian, (AU, AU / day)), erfa_ufunc.dt_eraLDBODY)


def astrom_unit():
    from astropy.units.astrophys import AU
    from astropy.units.si import rad, year

    one = rel2c = dimensionless_unscaled

    return StructuredUnit(
        (
            year,
            AU,
            one,
            AU,
            rel2c,
            one,
            one,
            rad,
            rad,
            rad,
            rad,
            one,
            one,
            rel2c,
            rad,
            rad,
            rad,
        ),
        erfa_ufunc.dt_eraASTROM,
    )


def helper_ldn(f, unit_b, unit_ob, unit_sc):
    from astropy.units.astrophys import AU

    return [
        get_converter(unit_b, ldbody_unit()),
        get_converter(unit_ob, AU),
        get_converter(_d(unit_sc), dimensionless_unscaled),
    ], dimensionless_unscaled


def helper_aper(f, unit_theta, unit_astrom):
    check_structured_unit(unit_astrom, dt_eraASTROM)
    unit_along = unit_astrom[7]  # along

    if unit_astrom[14] is unit_along:  # eral
        result_unit = unit_astrom
    else:
        result_units = tuple(
            (unit_along if i == 14 else v) for i, v in enumerate(unit_astrom.values())
        )
        result_unit = unit_astrom.__class__(result_units, names=unit_astrom)
    return [get_converter(unit_theta, unit_along), None], result_unit


def helper_apio(
    f,
    unit_sp,
    unit_theta,
    unit_elong,
    unit_phi,
    unit_hm,
    unit_xp,
    unit_yp,
    unit_refa,
    unit_refb,
):
    from astropy.units.si import m, radian

    return [
        get_converter(unit_sp, radian),
        get_converter(unit_theta, radian),
        get_converter(unit_elong, radian),
        get_converter(unit_phi, radian),
        get_converter(unit_hm, m),
        get_converter(unit_xp, radian),
        get_converter(unit_xp, radian),
        get_converter(unit_xp, radian),
        get_converter(unit_xp, radian),
    ], astrom_unit()


def helper_atciq(f, unit_rc, unit_dc, unit_pr, unit_pd, unit_px, unit_rv, unit_astrom):
    from astropy.units.si import arcsec, km, radian, s, year

    return [
        get_converter(unit_rc, radian),
        get_converter(unit_dc, radian),
        get_converter(unit_pr, radian / year),
        get_converter(unit_pd, radian / year),
        get_converter(unit_px, arcsec),
        get_converter(unit_rv, km / s),
        get_converter(unit_astrom, astrom_unit()),
    ], (radian, radian)


def helper_atciqn(
    f, unit_rc, unit_dc, unit_pr, unit_pd, unit_px, unit_rv, unit_astrom, unit_b
):
    from astropy.units.si import arcsec, km, radian, s, year

    return [
        get_converter(unit_rc, radian),
        get_converter(unit_dc, radian),
        get_converter(unit_pr, radian / year),
        get_converter(unit_pd, radian / year),
        get_converter(unit_px, arcsec),
        get_converter(unit_rv, km / s),
        get_converter(unit_astrom, astrom_unit()),
        get_converter(unit_b, ldbody_unit()),
    ], (radian, radian)


def helper_atciqz_aticq(f, unit_rc, unit_dc, unit_astrom):
    from astropy.units.si import radian

    return [
        get_converter(unit_rc, radian),
        get_converter(unit_dc, radian),
        get_converter(unit_astrom, astrom_unit()),
    ], (radian, radian)


def helper_aticqn(f, unit_rc, unit_dc, unit_astrom, unit_b):
    from astropy.units.si import radian

    return [
        get_converter(unit_rc, radian),
        get_converter(unit_dc, radian),
        get_converter(unit_astrom, astrom_unit()),
        get_converter(unit_b, ldbody_unit()),
    ], (radian, radian)


def helper_atioq(f, unit_rc, unit_dc, unit_astrom):
    from astropy.units.si import radian

    return [
        get_converter(unit_rc, radian),
        get_converter(unit_dc, radian),
        get_converter(unit_astrom, astrom_unit()),
    ], (radian,) * 5


def helper_atoiq(f, unit_type, unit_ri, unit_di, unit_astrom):
    from astropy.units.si import radian

    if unit_type is not None:
        raise UnitTypeError("argument 'type' should not have a unit")

    return [
        None,
        get_converter(unit_ri, radian),
        get_converter(unit_di, radian),
        get_converter(unit_astrom, astrom_unit()),
    ], (radian, radian)


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
    ERFA_HELPERS[erfa_ufunc.gc2gde] = helper_gc2gde
    ERFA_HELPERS[erfa_ufunc.gd2gc] = helper_gd2gc
    ERFA_HELPERS[erfa_ufunc.gd2gce] = helper_gd2gce
    ERFA_HELPERS[erfa_ufunc.ldn] = helper_ldn
    ERFA_HELPERS[erfa_ufunc.aper] = helper_aper
    ERFA_HELPERS[erfa_ufunc.apio] = helper_apio
    ERFA_HELPERS[erfa_ufunc.atciq] = helper_atciq
    ERFA_HELPERS[erfa_ufunc.atciqn] = helper_atciqn
    ERFA_HELPERS[erfa_ufunc.atciqz] = helper_atciqz_aticq
    ERFA_HELPERS[erfa_ufunc.aticq] = helper_atciqz_aticq
    ERFA_HELPERS[erfa_ufunc.aticqn] = helper_aticqn
    ERFA_HELPERS[erfa_ufunc.atioq] = helper_atioq
    ERFA_HELPERS[erfa_ufunc.atoiq] = helper_atoiq
    return ERFA_HELPERS


UFUNC_HELPERS.register_module("erfa.ufunc", erfa_ufuncs, get_erfa_helpers)
