# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Convenience functions for `astropy.cosmology`.
"""

import warnings
import numpy as np

from astropy.units import Quantity
from astropy.utils.exceptions import AstropyUserWarning
from .core import CosmologyError

__all__ = ['z_at_value']

__doctest_requires__ = {'*': ['scipy']}


def z_at_value(func, fval, zmin=1e-8, zmax=1000, ztol=1e-8, maxfun=500,
               method='Brent', bracket=None, verbose=False):
    """ Find the redshift ``z`` at which ``func(z) = fval``.

    This finds the redshift at which one of the cosmology functions or
    methods (for example Planck13.distmod) is equal to a known value.

    .. warning::
      Make sure you understand the behavior of the function that you
      are trying to invert! Depending on the cosmology, there may not
      be a unique solution. For example, in the standard Lambda CDM
      cosmology, there are two redshifts which give an angular
      diameter distance of 1500 Mpc, z ~ 0.7 and z ~ 3.8. To force
      ``z_at_value`` to find the solution you are interested in, use the
      ``zmin`` and ``zmax`` keywords to limit the search range (see the
      example below).

    Parameters
    ----------
    func : function or method
       A function that takes a redshift as input.

    fval : `~astropy.units.Quantity` instance
       The (scalar) value of ``func(z)`` to recover.

    zmin : float, optional
       The lower search limit for ``z``.  Beware of divergences
       in some cosmological functions, such as distance moduli,
       at z=0 (default 1e-8).

    zmax : float, optional
       The upper search limit for ``z`` (default 1000).

    ztol : float, optional
       The relative error in ``z`` acceptable for convergence.

    maxfun : int, optional
       The maximum number of function evaluations allowed in the
       optimization routine (default 500).

    method : str or callable, optional
       Type of solver to pass to the minimizer. The built-in options provided by
       :func:`~scipy.optimize.minimize_scalar` are 'Brent' (default), 'Golden' and
       'Bounded' with names case insensitive - see documentation there for details.
       It also accepts a custom solver by passing any user-provided callable object
       that meets the requirements listed therein under the Notes on
       "Custom minimizers" - or in more detail in :doc:`scipy:tutorial/optimize` -
       although their use is currently untested.

       .. versionadded:: 4.3

    bracket : sequence, optional
       For methods 'Brent' and 'Golden', ``bracket`` defines the bracketing
       interval and can either have three items (z1, z2, z3) so that z1 < z2 < z3
       and ``func(z2) < func(z1), func(z3)`` or two items z1 and z3 which are
       assumed to be a starting interval for a downhill bracket search.
       For non-monotone functions such as angular diameter distance this may be
       used to start the search on the desired side of the maximum, but see
       Examples below for usage notes.

       .. versionadded:: 4.3

    verbose : bool, optional
       Print diagnostic output from solver (default `False`).

       .. versionadded:: 4.3


    Returns
    -------
    z : float
      The redshift ``z`` satisfying ``zmin < z < zmax`` and ``func(z) =
      fval`` within ``ztol``.

    Notes
    -----
    This works for any arbitrary input cosmology, but is inefficient
    if you want to invert a large number of values for the same
    cosmology. In this case, it is faster to instead generate an array
    of values at many closely-spaced redshifts that cover the relevant
    redshift range, and then use interpolation to find the redshift at
    each value you are interested in. For example, to efficiently find
    the redshifts corresponding to 10^6 values of the distance modulus
    in a Planck13 cosmology, you could do the following:

    >>> import astropy.units as u
    >>> from astropy.cosmology import Planck13, z_at_value

    Generate 10^6 distance moduli between 24 and 44 for which we
    want to find the corresponding redshifts:

    >>> Dvals = (24 + np.random.rand(1000000) * 20) * u.mag

    Make a grid of distance moduli covering the redshift range we
    need using 50 equally log-spaced values between zmin and
    zmax. We use log spacing to adequately sample the steep part of
    the curve at low distance moduli:

    >>> zmin = z_at_value(Planck13.distmod, Dvals.min())
    >>> zmax = z_at_value(Planck13.distmod, Dvals.max())
    >>> zgrid = np.logspace(np.log10(zmin), np.log10(zmax), 50)
    >>> Dgrid = Planck13.distmod(zgrid)

    Finally interpolate to find the redshift at each distance modulus:

    >>> zvals = np.interp(Dvals.value, Dgrid.value, zgrid)

    Examples
    --------
    >>> import astropy.units as u
    >>> from astropy.cosmology import Planck13, Planck18, z_at_value

    The age and lookback time are monotonic with redshift, and so a
    unique solution can be found:

    >>> z_at_value(Planck13.age, 2 * u.Gyr)                # doctest: +FLOAT_CMP
    3.19812268

    The angular diameter is not monotonic however, and there are two
    redshifts that give a value of 1500 Mpc. You can use the zmin and
    zmax keywords to find the one you are interested in:

    >>> z_at_value(Planck18.angular_diameter_distance,
    ...            1500 * u.Mpc, zmax=1.5)                 # doctest: +FLOAT_CMP
    0.68044452
    >>> z_at_value(Planck18.angular_diameter_distance,
    ...            1500 * u.Mpc, zmin=2.5)                 # doctest: +FLOAT_CMP
    3.7823268

    Alternatively the ``bracket`` option may be used to initialize the
    function solver on a desired region, but one should be aware that this
    does not guarantee it will remain close to this starting bracket.
    For the example of angular diameter distance, which has a maximum near
    a redshift of 1.6 in this cosmology, defining a bracket on either side
    of this maximum will often return a solution on the same side:

    >>> z_at_value(Planck18.angular_diameter_distance,
    ...            1500 * u.Mpc, bracket=(1.0, 1.2))       # doctest: +FLOAT_CMP +IGNORE_WARNINGS
    0.68044452

    But this is not ascertained especially if the bracket is chosen too wide
    and/or too close to the turning point:

    >>> z_at_value(Planck18.angular_diameter_distance,
    ...            1500 * u.Mpc, bracket=(0.1, 1.5))       # doctest: +SKIP
    3.7823268                                              # doctest: +SKIP

    Likewise, even for the same minimizer and same starting conditions different
    results can be found depending on architecture or library versions:

    >>> z_at_value(Planck18.angular_diameter_distance,
    ...            1500 * u.Mpc, bracket=(2.0, 2.5))       # doctest: +SKIP
    3.7823268                                              # doctest: +SKIP

    >>> z_at_value(Planck18.angular_diameter_distance,
    ...            1500 * u.Mpc, bracket=(2.0, 2.5))       # doctest: +SKIP
    0.68044452                                             # doctest: +SKIP

    It is therefore generally safer to use the 3-parameter variant to ensure
    the solution stays within the bracketing limits:

    >>> z_at_value(Planck18.angular_diameter_distance,
    ...            1500 * u.Mpc, bracket=(0.1, 1.0, 1.5))  # doctest: +FLOAT_CMP
    0.68044452

    Also note that the luminosity distance and distance modulus (two
    other commonly inverted quantities) are monotonic in flat and open
    universes, but not in closed universes.
    """
    from scipy.optimize import minimize_scalar

    opt = {'maxiter': maxfun}
    # Assume custom methods support the same options as default; otherwise user will see warnings.
    if str(method).lower() == 'bounded':
        opt['xatol'] = ztol
        if bracket is not None:
            warnings.warn(f"Option 'bracket' is ignored by method {method}.")
            bracket = None
    else:
        opt['xtol'] = ztol

    # fval falling inside the interval of bracketing function values does not guarantee
    # it has a unique solution, but for Standard Cosmological quantities normally should
    # (being monotonic or having a single extremum).
    # In these cases keep solver from returning solutions outside of bracket.
    fval_zmin, fval_zmax = func((zmin, zmax))
    nobracket = False
    if np.sign(fval - fval_zmin) != np.sign(fval_zmax - fval):
        if bracket is None:
            nobracket = True
        else:
            fval_brac = func(bracket)
            if np.sign(fval - fval_brac[0]) != np.sign(fval_brac[-1] - fval):
                nobracket = True
            else:
                zmin, zmax = bracket[0], bracket[-1]
                fval_zmin, fval_zmax = fval_brac[[0, -1]]
    if nobracket:
        warnings.warn(f"fval is not bracketed by func(zmin)={fval_zmin} and func(zmax)="
                      f"{fval_zmax}. This means either there is no solution, or that there is "
                      "more than one solution between zmin and zmax satisfying fval = func(z).",
                      AstropyUserWarning)

    if isinstance(fval_zmin, Quantity):
        val = fval.to_value(fval_zmin.unit)
    else:
        val = fval

    # 'Brent' and 'Golden' ignore `bounds`, force solution inside zlim
    def f(z):
        if z > zmax:
            return 1.e300 * (1.0 + z - zmax)
        elif z < zmin:
            return 1.e300 * (1.0 + zmin - z)
        elif isinstance(fval_zmin, Quantity):
            return abs(func(z).value - val)
        else:
            return abs(func(z) - val)

    res = minimize_scalar(f, method=method, bounds=(zmin, zmax), bracket=bracket, options=opt)

    # Scipy docs state that `OptimizeResult` always has 'status' and 'message' attributes,
    # but only `_minimize_scalar_bounded()` seems to have really implemented them.
    if not res.success:
        warnings.warn(f"Solver returned {res.get('status')}: {res.get('message', 'Unsuccessful')}\n"
                      f"Precision {res.fun} reached after {res.nfev} function calls.",
                      AstropyUserWarning)

    if verbose:
        print(res)

    if np.allclose(res.x, zmax):
        raise CosmologyError(f"Best guess z={res.x} is very close to the upper z limit {zmax}.\n"
                             "Try re-running with a different zmax.")
    elif np.allclose(res.x, zmin):
        raise CosmologyError(f"Best guess z={res.x} is very close to the lower z limit {zmin}.\n"
                             "Try re-running with a different zmin.")
    return res.x
