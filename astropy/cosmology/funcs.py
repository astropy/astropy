# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Convenience functions for `astropy.cosmology`.
"""
from .core import get_current as _get_current


def kpc_comoving_per_arcmin(z, cosmo=None):
    """ Separation in transverse comoving kpc corresponding to an
    arcminute at redshift `z`.

    Parameters
    ----------
    z : array_like
      Input redshifts.

    Returns
    -------
    d : ndarray, or float if input scalar
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
    d : ndarray, or float if input scalar
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
    theta : ndarray, or float if input scalar
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
    theta : ndarray, or float if input scalar
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
    distmod : ndarray, or float if input scalar
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
    H : ndarray, or float if input scalar
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
    distmod : ndarray, or float if input scalar
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
    critdens : ndarray, or float if input scalar
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
    t : ndarray, or float if input scalar
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
    codist : ndarray, or float if input scalar
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
    angdist : ndarray, or float if input scalar
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
    lumdist : ndarray, or float if input scalar
      Angular diameter distance at each input redshift.
    """
    if cosmo is None:
        cosmo = _get_current()
    return cosmo.luminosity_distance(z)
