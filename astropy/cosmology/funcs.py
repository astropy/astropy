# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Convenience functions for `astropy.cosmology`.
"""

import warnings
import numpy as np

try:
    from scipy.interpolate import CubicSpline
except ImportError:
    HAS_SCIPY = False
else:
    HAS_SCIPY = True

from .core import CosmologyError
from astropy.units import Quantity

__all__ = ['z_at_value']

__doctest_requires__ = {'*': ['scipy']}

_interpolation_choices = ['cubic', 'linear']


def _z_at_array(func, fvals, zmin, zmax, nbins=10000,
                logspace=True, interpolation=None):
    """Helper function to interpolate (func, z) over a grid for array input"""
    if logspace:
        zgrid = np.logspace(np.log10(zmin), np.log10(zmax), nbins)
    else:
        zgrid = np.linspace(zmin, zmax, nbins)

    fgrid = func(zgrid)

    fvals_val = fvals.value
    fgrid_val = fgrid.to_value(fvals.unit)

    if interpolation is None:
        if HAS_SCIPY:
            interpolation = 'cubic'
        else:
            warnings.warn((
                "SciPy not found, so falling back to linear numpy interpolation"
                " scheme for array z_at_value, instead of the normal cubic"
                " spline."))
            interpolation = 'linear'

    if interpolation not in _interpolation_choices:
        raise ValueError(
                f"{interpolation} is an invalid interpolation method."
                f" Please use an option from {_interpolation_choices}."
            )

    if interpolation == 'cubic':
        interpolator = CubicSpline(fgrid_val, zgrid)
        zvals = interpolator(fvals_val)
    else:
        zvals = np.interp(fvals_val, fgrid_val, zgrid)

    return zvals


def z_at_value(func, fval, zmin=1e-8, zmax=1000, ztol=1e-8, maxfun=500,
               nbins=10000, logspace=True, interpolation=None):
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

    .. warning::
      The results of this function will likely differ depending on whether
      you input a single value for ``fval``, or an array of values.
      In the single value case, this function will attempt to find
      something close to the exact value using ``fminbound`` (proximity
      depending on ``ztol``), but in the array case, this function
      will use a cubic spline interpolation scheme, which will give
      approximate answers which increase in accuracy as ``nbins``
      is increased.

    Parameters
    ----------
    func : callable
       A function that takes a redshift as input.
    fval : astropy.Quantity instance or array_like
       The value of ``func(z)``.
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
    nbins : float, optional
        If passing an array of ``fval``, this specifies the number
        of gridpoints to use for interpolation.
    logspace: bool, optional
        If passing an array of ``fval``, this specifies whether
        to create the gridpoints in logarithmic space
    interpolation: None or string, optional
        To force linear interpolation, use ``'linear'``. Otherwise,
        the choice depends on whether scipy is installed (cubic) or not
        (linear).

    Returns
    -------
    z : float or array_like
      The redshift ``z`` satisfying ``zmin < z < zmax`` and ``func(z) =
      fval`` within ``ztol``.

    Examples
    --------
    >>> import astropy.units as u
    >>> from astropy.cosmology import Planck13, z_at_value

    The age and lookback time are monotonic with redshift, and so a
    unique solution can be found:

    >>> z_at_value(Planck13.age, 2 * u.Gyr)  # doctest: +FLOAT_CMP
    3.19812268

    The angular diameter is not monotonic however, and there are two
    redshifts that give a value of 1500 Mpc. Use the zmin and zmax keywords
    to find the one you're interested in:

    >>> z_at_value(Planck13.angular_diameter_distance,
    ...            1500 * u.Mpc, zmax=1.5)  # doctest: +FLOAT_CMP
    0.6812769577
    >>> z_at_value(Planck13.angular_diameter_distance,
    ...            1500 * u.Mpc, zmin=2.5)  # doctest: +FLOAT_CMP
    3.7914913242

    Also note that the luminosity distance and distance modulus (two
    other commonly inverted quantities) are monotonic in flat and open
    universes, but not in closed universes.
    """
    is_array = (
           (isinstance(fval, Quantity) and not fval.isscalar) or
           (not isinstance(fval, Quantity) and isinstance(fval, np.ndarray))
           )
    if is_array:
        return _z_at_array(func, fval,
                           zmin=zmin, zmax=zmax,
                           nbins=nbins, logspace=logspace,
                           interpolation=interpolation)

    from scipy.optimize import fminbound

    fval_zmin = func(zmin)
    fval_zmax = func(zmax)
    if np.sign(fval - fval_zmin) != np.sign(fval_zmax - fval):
        warnings.warn("""\
fval is not bracketed by func(zmin) and func(zmax). This means either
there is no solution, or that there is more than one solution between
zmin and zmax satisfying fval = func(z).""")

    if isinstance(fval_zmin, Quantity):
        val = fval.to_value(fval_zmin.unit)
        f = lambda z: abs(func(z).value - val)
    else:
        f = lambda z: abs(func(z) - fval)

    zbest, resval, ierr, ncall = fminbound(f, zmin, zmax, maxfun=maxfun,
                                           full_output=1, xtol=ztol)

    if ierr != 0:
        warnings.warn(f'Maximum number of function calls ({ncall}) reached')

    if np.allclose(zbest, zmax):
        raise CosmologyError("Best guess z is very close the upper z limit.\n"
                             "Try re-running with a different zmax.")
    elif np.allclose(zbest, zmin):
        raise CosmologyError("Best guess z is very close the lower z limit.\n"
                             "Try re-running with a different zmin.")

    return zbest
