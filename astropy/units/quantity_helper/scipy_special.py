# UFUNCS FROM SCIPY.SPECIAL
# available ufuncs in this module are at
# https://docs.scipy.org/doc/scipy/reference/special.html

from ...utils import minversion
from ..core import UnitsError, UnitTypeError, dimensionless_unscaled
from . import UFUNC_HELPERS
from .helpers import (get_converter,
                      helper_dimensionless_to_dimensionless,
                      helper_cbrt,
                      helper_two_arg_dimensionless)


def helper_degree_to_dimensionless(f, unit):
    from ..si import degree
    try:
        return [get_converter(unit, degree)], dimensionless_unscaled
    except UnitsError:
        raise UnitTypeError("Can only apply '{0}' function to "
                            "quantities with angle units"
                            .format(f.__name__))


def helper_degree_minute_second_to_radian(f, unit1, unit2, unit3):
    from ..si import degree, arcmin, arcsec, radian
    try:
        return [get_converter(unit1, degree),
                get_converter(unit2, arcmin),
                get_converter(unit3, arcsec)], radian
    except UnitsError:
        raise UnitTypeError("Can only apply '{0}' function to "
                            "quantities with angle units"
                            .format(f.__name__))


try:
    import scipy
    import scipy.special as sps
except ImportError:
    pass
else:
    # ufuncs that require dimensionless input and give dimensionless output
    dimensionless_to_dimensionless_sps_ufuncs = [
        sps.erf, sps.gamma, sps.gammasgn,
        sps.psi, sps.rgamma, sps.erfc, sps.erfcx, sps.erfi, sps.wofz,
        sps.dawsn, sps.entr, sps.exprel, sps.expm1, sps.log1p, sps.exp2,
        sps.exp10, sps.j0, sps.j1, sps.y0, sps.y1, sps.i0, sps.i0e, sps.i1,
        sps.i1e, sps.k0, sps.k0e, sps.k1, sps.k1e, sps.itj0y0,
        sps.it2j0y0, sps.iti0k0, sps.it2i0k0]

    # TODO: Revert https://github.com/astropy/astropy/pull/7219 when astropy
    #       requires scipy>=0.18.
    # See https://github.com/astropy/astropy/issues/7159
    if minversion(scipy, "0.18"):
        dimensionless_to_dimensionless_sps_ufuncs.append(sps.loggamma)

    for ufunc in dimensionless_to_dimensionless_sps_ufuncs:
        UFUNC_HELPERS[ufunc] = helper_dimensionless_to_dimensionless

    # ufuncs that require input in degrees and give dimensionless output
    degree_to_dimensionless_sps_ufuncs = (
        sps.cosdg, sps.sindg, sps.tandg, sps.cotdg)
    for ufunc in degree_to_dimensionless_sps_ufuncs:
        UFUNC_HELPERS[ufunc] = helper_degree_to_dimensionless

    # ufuncs that require 2 dimensionless inputs and give dimensionless output.
    # note: sps.jv and sps.jn are aliases in some scipy versions, which will
    # cause the same key to be written twice, but since both are handled by the
    # same helper there is no harm done.
    two_arg_dimensionless_sps_ufuncs = (
        sps.jv, sps.jn, sps.jve, sps.yn, sps.yv, sps.yve, sps.kn, sps.kv,
        sps.kve, sps.iv, sps.ive, sps.hankel1, sps.hankel1e, sps.hankel2,
        sps.hankel2e)
    for ufunc in two_arg_dimensionless_sps_ufuncs:
        UFUNC_HELPERS[ufunc] = helper_two_arg_dimensionless

    # ufuncs handled as special cases
    UFUNC_HELPERS[sps.cbrt] = helper_cbrt
    UFUNC_HELPERS[sps.radian] = helper_degree_minute_second_to_radian
