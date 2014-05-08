# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Convenience functions for `astropy.cosmology`.
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import warnings
import numpy as np

from .core import default_cosmology as _default_cosmology
from .core import CosmologyError
from ..units import Quantity
from ..utils import deprecated

__all__ = ['H', 'angular_diameter_distance', 'arcsec_per_kpc_comoving',
           'arcsec_per_kpc_proper', 'comoving_distance', 'critical_density',
           'distmod', 'kpc_comoving_per_arcmin', 'kpc_proper_per_arcmin',
           'lookback_time', 'luminosity_distance', 'scale_factor',
           'z_at_value']

__doctest_requires__ = {'*': ['scipy.integrate']}


def z_at_value(func, fval, zmin=0, zmax=1000, ztol=1e-5, maxfun=500):
    """ Find the redshift ``z`` at which ``func(z) = fval``.

    This finds the redshift at which one of the cosmology functions or
    methods (for example Planck13.distmod) is equal to a known value.

    .. warning::
      Make sure you understand the behaviour of the function that you
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
    fval : astropy.Quantity instance
       The value of ``func(z)``.
    zmin : float, optional
       The lower search limit for ``z`` (default 0).
    zmax : float, optional
       The upper search limit for ``z`` (default 1000).
    ztol : float, optional
       The relative error in ``z`` acceptable for convergence.
    maxfun : int, optional
       The maximum number of function evaluations allowed in the
       optimization routine (default 500).

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
    each value you're interested in. For example, to efficiently find
    the redshifts corresponding to 10^6 values of the distance modulus
    in a Planck13 cosmology, you could do the following:

    >>> import astropy.units as u
    >>> from astropy.cosmology import Planck13, z_at_value

    Generate 10^6 distance moduli between 23 and 43 for which we
    want to find the corresponding redshifts:

    >>> Dvals = (23 + np.random.rand(1e6) * 20) * u.mag

    Make a grid of distance moduli covering the redshift range we
    need using 50 equally log-spaced values between zmin and
    zmax. We use log spacing to adequately sample the steep part of
    the curve at low distance moduli:

    >>> zmin = z_at_value(Planck13.distmod, Dvals.min())
    >>> zmax = z_at_value(Planck13.distmod, Dvals.max())
    >>> zgrid = np.logspace(zmin, zmax)
    >>> Dgrid = Planck13.distmod(zgrid)

    Finally interpolate to find the redshift at each distance modulus:

    >>> zvals = np.interp(Dvals.value, zgrid, Dgrid.value)

    Examples
    --------
    >>> import astropy.units as u
    >>> from astropy.cosmology import Planck13, z_at_value

    The age and lookback time are monotonic with redshift, and so a
    unique solution can be found:

    >>> z_at_value(Planck13.age, 2 * u.Gyr)
    3.1981191749374629

    The angular diameter is not monotonic however, and there are two
    redshifts that give a value of 1500 Mpc. Use the zmin and zmax keywords
    to find the one you're interested in:

    >>> z_at_value(Planck13.angular_diameter_distance, 1500 * u.Mpc, zmax=1.5)
    0.68127769625288614
    >>> z_at_value(Planck13.angular_diameter_distance, 1500 * u.Mpc, zmin=2.5)
    3.7914918534022011

    Also note that the luminosity distance and distance modulus (two
    other commonly inverted quantities) are monotonic in flat and open
    universes, but not in closed universes.
    """
    from scipy.optimize import fminbound

    fval_zmin = func(zmin)
    fval_zmax = func(zmax)
    if np.sign(fval - fval_zmin) != np.sign(fval_zmax - fval):
        warnings.warn("""\
fval is not bracketed by func(zmin) and func(zmax). This means either
there is no solution, or that there is more than one solution between
zmin and zmax satisfying fval = func(z).""")

    if isinstance(fval_zmin, Quantity):
        unit = fval_zmin.unit
        val = fval.to(unit).value
        f = lambda z: abs(func(z).value - val)
    else:
        f = lambda z: abs(func(z) - fval)

    zbest, resval, ierr, ncall = fminbound(f, zmin, zmax, maxfun=maxfun,
                                           full_output=1)

    if ierr != 0:
        warnings.warn('Maximum number of function calls ({}) reached'.format(
            ncall))

    if np.allclose(zbest, zmax):
        raise CosmologyError("Best guess z is very close the upper z limit.\n"
                             "Try re-running with a different zmax.")
    elif np.allclose(zbest, zmin):
        raise CosmologyError("Best guess z is very close the lower z limit.\n"
                             "Try re-running with a different zmin.")

    return zbest


@deprecated(since='0.4', alternative='<Cosmology object>.kpc_comoving_per_arcmin')
def kpc_comoving_per_arcmin(z, cosmo=None):
    """ Separation in transverse comoving kpc corresponding to an
    arcminute at redshift ``z``.

    Parameters
    ----------
    z : array_like
      Input redshifts.

    Returns
    -------
    d : `~astropy.units.Quantity`
      The distance in comoving kpc corresponding to an arcmin at each
      input redshift.
    """
    if cosmo is None:
        cosmo = _default_cosmology.get()
    return cosmo.kpc_comoving_per_arcmin(z)


@deprecated(since='0.4', alternative='<Cosmology object>.kpc_proper_per_arcmin')
def kpc_proper_per_arcmin(z, cosmo=None):
    """ Separation in transverse proper kpc corresponding to an
    arcminute at redshift ``z``.

    Parameters
    ----------
    z : array_like
      Input redshifts.

    Returns
    -------
    d : `~astropy.units.Quantity`
      The distance in proper kpc corresponding to an arcmin at each
      input redshift.
    """
    if cosmo is None:
        cosmo = _default_cosmology.get()
    return cosmo.kpc_proper_per_arcmin(z)


@deprecated(since='0.4', alternative='<Cosmology object>.arcsec_per_kpc_comoving')
def arcsec_per_kpc_comoving(z, cosmo=None):
    """ Angular separation in arcsec corresponding to a comoving kpc
    at redshift ``z``.

    Parameters
    ----------
    z : array_like
      Input redshifts.

    Returns
    -------
    theta : `~astropy.units.Quantity`
      The angular separation in arcsec corresponding to a comoving kpc
      at each input redshift.
    """
    if cosmo is None:
        cosmo = _default_cosmology.get()
    return cosmo.arcsec_per_kpc_comoving(z)


@deprecated(since='0.4', alternative='<Cosmology object>.arcsec_per_kpc_proper')
def arcsec_per_kpc_proper(z, cosmo=None):
    """ Angular separation in arcsec corresponding to a proper kpc at
    redshift ``z``.

    Parameters
    ----------
    z : array_like
      Input redshifts.

    Returns
    -------
    theta : `~astropy.units.Quantity`
      The angular separation in arcsec corresponding to a proper kpc
      at each input redshift.
    """
    if cosmo is None:
        cosmo = _default_cosmology.get()
    return cosmo.arcsec_per_kpc_proper(z)


@deprecated(since='0.4', alternative='<Cosmology object>.distmod')
def distmod(z, cosmo=None):
    """ Distance modulus at redshift ``z``.

    The distance modulus is defined as the (apparent magnitude -
    absolute magnitude) for an object at redshift ``z``.

    Parameters
    ----------
    z : array_like
      Input redshifts.

    Returns
    -------
    distmod : `~astropy.units.Quantity`
      Distance modulus at each input redshift.

    See Also
    --------
    z_at_value : Find the redshift corresponding to a distance modulus.
    """
    if cosmo is None:
        cosmo = _default_cosmology.get()
    return cosmo.distmod(z)


@deprecated(since='0.4', alternative='<Cosmology object>.H')
def H(z, cosmo=None):
    """ Hubble parameter (km/s/Mpc) at redshift ``z``.

    Parameters
    ----------
    z : array_like
      Input redshifts.

    Returns
    -------
    H : `~astropy.units.Quantity`
      Hubble parameter at each input redshift.
    """
    if cosmo is None:
        cosmo = _default_cosmology.get()
    return cosmo.H(z)


@deprecated(since='0.4', alternative='<Cosmology object>.scale_factor')
def scale_factor(z, cosmo=None):
    """ Scale factor at redshift ``z``.

    The scale factor is defined as ``a = 1 / (1 + z)``.

    Parameters
    ----------
    z : array_like
      Input redshifts.

    Returns
    -------
    scalefac : ndarray, or float if input scalar
      Scale factor at each input redshift.
    """
    if cosmo is None:
        cosmo = _default_cosmology.get()
    return cosmo.scale_factor(z)


@deprecated(since='0.4', alternative='<Cosmology object>.critical_density')
def critical_density(z, cosmo=None):
    """ Critical density in grams per cubic cm at redshift ``z``.

    Parameters
    ----------
    z : array_like
      Input redshifts.

    Returns
    -------
    critdens : `~astropy.units.Quantity`
      Critical density at each input redshift.
    """
    if cosmo is None:
        cosmo = _default_cosmology.get()
    return cosmo.critical_density(z)


@deprecated(since='0.4', alternative='<Cosmology object>.lookback_time')
def lookback_time(z, cosmo=None):
    """ Lookback time in Gyr to redshift ``z``.

    The lookback time is the difference between the age of the
    Universe now and the age at redshift ``z``.

    Parameters
    ----------
    z : array_like
      Input redshifts.

    Returns
    -------
    t : `~astropy.units.Quantity`
      Lookback time at each input redshift.

    See Also
    --------
    z_at_value : Find the redshift corresponding to a lookback time.
    """
    if cosmo is None:
        cosmo = _default_cosmology.get()
    return cosmo.lookback_time(z)


@deprecated(since='0.4', alternative='<Cosmology object>.comoving_distance')
def comoving_distance(z, cosmo=None):
    """ Comoving distance in Mpc at redshift ``z``.

    The comoving distance along the line-of-sight between two objects
    remains constant with time for objects in the Hubble flow.

    Parameters
    ----------
    z : array_like
      Input redshifts.

    Returns
    -------
    codist : `~astropy.units.Quantity`
      Comoving distance at each input redshift.
    """
    if cosmo is None:
        cosmo = _default_cosmology.get()
    return cosmo.comoving_distance(z)


@deprecated(since='0.4', alternative='<Cosmology object>.angular_diameter_distance')
def angular_diameter_distance(z, cosmo=None):
    """ Angular diameter distance in Mpc at a given redshift.

    This gives the proper (sometimes called 'physical') transverse
    distance corresponding to an angle of 1 radian for an object at
    redshift ``z``.

    Parameters
    ----------
    z : array_like
      Input redshifts.

    Returns
    -------
    angdist : `~astropy.units.Quantity`
      Angular diameter distance at each input redshift.
    """
    if cosmo is None:
        cosmo = _default_cosmology.get()
    return cosmo.angular_diameter_distance(z)


@deprecated(since='0.4', alternative='<Cosmology object>.luminosity_distance')
def luminosity_distance(z, cosmo=None):
    """ Luminosity distance in Mpc at redshift ``z``.

    This is the distance to use when converting between the bolometric
    flux from an object at redshift ``z`` and its bolometric luminosity.

    Parameters
    ----------
    z : array_like
      Input redshifts.

    Returns
    -------
    lumdist : `~astropy.units.Quantity`
      Luminosity distance at each input redshift.

    See Also
    --------
    z_at_value : Find the redshift corresponding to a luminosity distance.
    """
    if cosmo is None:
        cosmo = _default_cosmology.get()
    return cosmo.luminosity_distance(z)
