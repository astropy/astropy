# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Convenience functions for `astropy.cosmology`.
"""
import warnings
import numpy as np
from .core import get_current as _get_current
from .core import CosmologyError
from ..units import Quantity


def z_at_value(func, fval, zmin=0, zmax=1e4, ztol=1e-5):
    """ Find the redshift `z` at which `func(z) = fval`.

    This function finds the redshift at which one of the cosmology
    functions (for example Planck13.distmod) is equal to a known
    value. WARNING: Make sure you understand the behaviour of the
    function that you are trying to invert! Depending on the
    cosmology, there may not be a unique solution. For example, in the
    standard Lambda CDM cosmology, there are two redshifts which give
    an angular diameter distance of 1500 Mpc, z ~ 0.7 and z ~ 3.8. To
    force `z_at_value` to find the solution you are interested in, use
    the `zmin` and `zmax` keywords to limit the search range (see the
    example below).

    Parameters
    ----------
    func : function or method
       A function that takes a redshift as input.
    fval : astropy.Quantity instance
       The value of `func(z)`.
    zmin : float
       The lower search limit for `z` (default 0).
    zmax : float
       The upper search limit for `z` (default 10,000).
    ztol : float
       The relative error in `z` acceptable for convergence.

    Returns
    -------
    z : float
      The redshift `z` satisfying `zmin < z < zmax` and `func(z) =
      fval` within `ztol`.

    Notes
    -----
    This function works for any arbitrary input cosmology, but is slow
    and inefficient if you want to invert a very large number of
    values for the same cosmology. In this case, you may want to
    generate an array of function values at many closely-spaced
    redshifts that cover the redshift range you're interested in, and
    then use interpolation to find the redshifts.

    Examples
    --------
    >>> import astropy.units as u
    >>> from astropy.cosmology import Planck13, z_at_value
    >>> z_at_value(Planck13.age, 2 * u.Gyr)
    3.1981191749374629
    >>> z_at_value(Planck13.angular_diameter_distance, 1500 * u.Mpc, zmax=1.5)
    0.68127769625288614
    >>> z_at_value(Planck13.angular_diameter_distance, 1500 * u.Mpc, zmin=2.5)
    3.7914918534022011
    """
    from scipy.optimize import fminbound

    testval = func(zmin)
    if isinstance(testval, Quantity):
        unit = testval.unit
        val = fval.to(unit).value
        f = lambda z: abs(func(z).value - val)
    else:
        f = lambda z: abs(func(z) - fval)

    zbest, resval, ierr, ncall = fminbound(f, zmin, zmax, full_output=1)

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


def age(z, cosmo=None):
    """ Age of the universe in Gyr at redshift `z`.

    Parameters
    ----------
    z : array_like
      Input redshifts.

    Returns
    -------
    t : astropy.units.Quantity
      The age of the universe in Gyr at each input redshift.
        """
    if cosmo is None:
        cosmo = _get_current()
    return cosmo.age(z)


def comoving_volume(z, cosmo=None):
    """ Comoving volume in cubic Mpc at redshift `z`.

    This is the volume of the universe encompassed by redshifts
    less than `z`. For the case of omega_k = 0 it is a sphere of
    radius `comoving_distance(z)` but it is less intuitive if
    omega_k is not 0.

    Parameters
    ----------
    z : array_like
      Input redshifts.

    Returns
    -------
    V : astropy.units.Quantity
      Comoving volume in :math:`Mpc^3` at each input redshift.
    """
    if cosmo is None:
        cosmo = _get_current()
    return cosmo.comoving_volume(z)


def kpc_comoving_per_arcmin(z, cosmo=None):
    """ Separation in transverse comoving kpc corresponding to an
    arcminute at redshift `z`.

    Parameters
    ----------
    z : array_like
      Input redshifts.

    Returns
    -------
    d : astropy.units.Quantity
      The distance in comoving kpc corresponding to an arcmin at each
      input redshift.
    """
    if cosmo is None:
        cosmo = _get_current()
    return cosmo.kpc_comoving_per_arcmin(z)


def kpc_proper_per_arcmin(z, cosmo=None):
    """ Separation in transverse proper kpc corresponding to an
    arcminute at redshift `z`.

    Parameters
    ----------
    z : array_like
      Input redshifts.

    Returns
    -------
    d : astropy.units.Quantity
      The distance in proper kpc corresponding to an arcmin at each
      input redshift.
    """
    if cosmo is None:
        cosmo = _get_current()
    return cosmo.kpc_proper_per_arcmin(z)


def arcsec_per_kpc_comoving(z, cosmo=None):
    """ Angular separation in arcsec corresponding to a comoving kpc
    at redshift `z`.

    Parameters
    ----------
    z : array_like
      Input redshifts.

    Returns
    -------
    theta : astropy.units.Quantity
      The angular separation in arcsec corresponding to a comoving kpc
      at each input redshift.
    """
    if cosmo is None:
        cosmo = _get_current()
    return cosmo.arcsec_per_kpc_comoving(z)


def arcsec_per_kpc_proper(z, cosmo=None):
    """ Angular separation in arcsec corresponding to a proper kpc at
    redshift `z`.

    Parameters
    ----------
    z : array_like
      Input redshifts.

    Returns
    -------
    theta : astropy.units.Quantity
      The angular separation in arcsec corresponding to a proper kpc
      at each input redshift.
    """
    if cosmo is None:
        cosmo = _get_current()
    return cosmo.arcsec_per_kpc_proper(z)


def distmod(z, cosmo=None):
    """ Distance modulus at redshift `z`.

    The distance modulus is defined as the (apparent magnitude -
    absolute magnitude) for an object at redshift `z`.

    Parameters
    ----------
    z : array_like
      Input redshifts.

    Returns
    -------
    distmod : astropy.units.Quantity
      Distance modulus at each input redshift.
    """
    if cosmo is None:
        cosmo = _get_current()
    return cosmo.distmod(z)


def H(z, cosmo=None):
    """ Hubble parameter (km/s/Mpc) at redshift `z`.

    Parameters
    ----------
    z : array_like
      Input redshifts.

    Returns
    -------
    H : astropy.units.Quantity
      Hubble parameter at each input redshift.
    """
    if cosmo is None:
        cosmo = _get_current()
    return cosmo.H(z)


def scale_factor(z, cosmo=None):
    """ Scale factor at redshift `z`.

    The scale factor is defined as `a = 1 / (1 + z)`.

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
        cosmo = _get_current()
    return cosmo.scale_factor(z)


def critical_density(z, cosmo=None):
    """ Critical density in grams per cubic cm at redshift `z`.

    Parameters
    ----------
    z : array_like
      Input redshifts.

    Returns
    -------
    critdens : astropy.units.Quantity
      Critical density at each input redshift.
    """
    if cosmo is None:
        cosmo = _get_current()
    return cosmo.critical_density(z)


def lookback_time(z, cosmo=None):
    """ Lookback time in Gyr to redshift `z`.

    The lookback time is the difference between the age of the
    Universe now and the age at redshift `z`.

    Parameters
    ----------
    z : array_like
      Input redshifts.

    Returns
    -------
    t : astropy.units.Quantity
      Lookback time at each input redshift.
    """
    if cosmo is None:
        cosmo = _get_current()
    return cosmo.lookback_time(z)


def comoving_distance(z, cosmo=None):
    """ Comoving distance in Mpc at redshift `z`.

    The comoving distance along the line-of-sight between two objects
    remains constant with time for objects in the Hubble flow.

    Parameters
    ----------
    z : array_like
      Input redshifts.

    Returns
    -------
    codist : astropy.units.Quantity
      Comoving distance at each input redshift.
    """
    if cosmo is None:
        cosmo = _get_current()
    return cosmo.comoving_distance(z)


def angular_diameter_distance(z, cosmo=None):
    """ Angular diameter distance in Mpc at a given redshift.

    This gives the proper (sometimes called 'physical') transverse
    distance corresponding to an angle of 1 radian for an object at
    redshift `z`.

    Parameters
    ----------
    z : array_like
      Input redshifts.

    Returns
    -------
    angdist : astropy.units.Quantity
      Angular diameter distance at each input redshift.
    """
    if cosmo is None:
        cosmo = _get_current()
    return cosmo.angular_diameter_distance(z)


def luminosity_distance(z, cosmo=None):
    """ Luminosity distance in Mpc at redshift `z`.

    This is the distance to use when converting between the bolometric
    flux from an object at redshift `z` and its bolometric luminosity.

    Parameters
    ----------
    z : array_like
      Input redshifts.

    Returns
    -------
    lumdist : astropy.units.Quantity
      Angular diameter distance at each input redshift.
    """
    if cosmo is None:
        cosmo = _get_current()
    return cosmo.luminosity_distance(z)
