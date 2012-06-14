# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Convenience functions for `astropy.cosmology`.
"""
from .core import get_current as _get_current
from math import pi as _pi

_arcsec_in_radians = 1 / 3600. * _pi / 180
_arcmin_in_radians = 1 / 60. * _pi / 180


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
    return cosmo.comoving_transverse_distance(z) * 1.e3 * _arcmin_in_radians


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
    return cosmo.angular_diameter_distance(z) * 1.e3 * _arcmin_in_radians


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
    return 1 / (cosmo.comoving_transverse_distance(z) *
                1.e3 * _arcsec_in_radians)


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
    return 1 / (cosmo.angular_diameter_distance(z) * 1.e3 * _arcsec_in_radians)


def distmod(z, cosmo=None):
    """ Distance modulus at redshift `z`.

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


def H(z=0, cosmo=None):
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


def scale_factor(z=0, cosmo=None):
    """ Scale factor at redshift `z`.

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


def critical_density(z=0, cosmo=None):
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


def lookback_time(z=0, cosmo=None):
    """ Lookback time in Gyr to redshift `z`.

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


def comoving_distance(z=0, cosmo=None):
    """ Comoving distance in Mpc at redshift `z`.

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


def angular_diameter_distance(z=0, cosmo=None):
    """ Angular diameter distance in Mpc at a given redshift.

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


def luminosity_distance(z=0, cosmo=None):
    """ Angular diameter distance in Mpc at a given redshift.

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
