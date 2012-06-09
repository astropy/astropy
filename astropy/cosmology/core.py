# Licensed under a 3-clause BSD style license - see LICENSE.rst
import sys
import warnings
from math import sqrt, pi

import numpy as np

from ..constants.cgs import pc, G, c
from ..config import ConfigurationItem
from ..utils.misc import isiterable

import parameters

# Originally authored by Andrew Becker (becker@astro.washington.edu),
# and modified by Neil Crighton (neilcrighton@gmail.com) and Roban
# Kramer (robanhk@gmail.com).

# Many of these adapted from astro-ph/9905116

__all__ = ("FLRWCosmology kpc_comoving_per_arcmin kpc_proper_per_arcmin "
           "arcsec_per_kpc_comoving arcsec_per_kpc_proper distmod "
           "get_current set_current WMAP5 WMAP7").split()

# Constants

# speed of light in km/s
c_kms = c * 1e-5

# Mpc in cm
Mpc = 1e6 * pc

# Mpc in km
Mpc_km = 1e-5 * Mpc

# Gyr in seconds
Gyr = 1e9 * 365.25 * 24 * 60 * 60

arcsec_in_radians = 1 / 3600. * pi / 180
arcmin_in_radians = 1 / 60. * pi / 180

DEFAULT_COSMOLOGY = ConfigurationItem(
    'default_cosmology', 'no_default',
    'The default cosmology to use. Note this is only read on import, '
    'changing this value at runtime has no effect.')


class CosmologyError(Exception):
    pass

class Cosmology(object):
    """ Placeholder for when a more general Cosmology class is
    implemented. """
    pass

class FLRWCosmology(Cosmology):
    """ A class describing an isotropic and homogeneous
    (Friedmann-Lemaitre-Robertson-Walker) cosmology.

    Attributes
    ----------
    H0 : float
      Hubble parameter at z=0 in km/s/Mpc
    Om : float
      Omega matter; matter density / critical density at z=0
    Ol : float
      Omega lambda; dark energy density / critical density at z=0
    Ok : float
      Omega_k, the curvature density at z=0. Defined as 1 - Om - Ol
    h : float
      Dimensionless Hubble parameter (H0 = 100*h km/s/Mpc).
      Often used to quote cosmological values independently of H0.
    hubble_time : float
      Hubble time in Gyr.
    hubble_distance : float
      Hubble distance in Mpc.
    critical_density0 : float
      Critical density in g cm^-3 at z=0.

    Notes
    -----
    Note the energy density from radiation, omega_r, is ignored (valid
    for redshifts < ~10).

    Examples
    --------
    >>> from astro.cosmology import FLRWCosmology
    >>> cosmo = FLRWCosmology(H0=70, Om=0.3, Ol=0.7)

    The comoving distance in Mpc at redshift z:

    >>> dc = cosmo.comoving_distance(z)
    """
    def __init__(self, H0, Om, Ol, name='FLRWCosmology'):

        # all densities are in units of the critical density
        self.Om = float(Om)
        self.Ol = float(Ol)
        Ok = 1 - self.Om - self.Ol
        if abs(Ok) < 1e-5:
            Ok = 0
        self.Ok = Ok
        self.name = name

        # Hubble parameter at z=0, km/s/Mpc
        self.H0 = float(H0)
        # H0 in s^-1
        H0_s = self.H0 / Mpc_km
        # 100 km/s/Mpc * h = H0 (so h is dimensionless)
        self.h = self.H0 / 100.
        # Hubble time in Gyr
        self.hubble_time = 1. / H0_s / Gyr
        # Hubble distance in Mpc
        self.hubble_distance = c_kms / self.H0

        # critical density at z=0 (grams per cubic cm)
        self.critical_density0 = 3. * H0_s**2 / (8. * pi * G)

    def __repr__(self):
        s = "%s(H0=%.3g, Om=%.3g, Ol=%.3g, Ok=%.3g)" % (
            self.name, self.H0, self.Om, self.Ol, self.Ok)
        return s

    def _efunc(self, z):
        """ Function used to calculate the hubble parameter as a
        function of redshift. Eqn 14 from Hogg."""
        if isiterable(z):
            z = np.asarray(z)
        zp1 = 1. + z
        return np.sqrt(self.Om*zp1**3 + self.Ok*zp1**2 + self.Ol)

    def _inv_efunc(self, z):
        """ Integrand of the comoving distance.
        """
        zp1 = 1. + z
        return 1. / np.sqrt(self.Om*zp1**3 + self.Ok*zp1**2 + self.Ol)

    def _tfunc(self, z):
        """ Integrand of the lookback time.

        Eqn 30 from Hogg."""
        zp1 = 1. + z
        return 1. / (zp1*np.sqrt(self.Om*zp1**3 + self.Ok*zp1**2 + self.Ol))

    def _xfunc(self, z):
        """ Integrand of the absorption distance.

        See Hogg 1999 section 11.
        """
        zp1 = 1. + z
        return zp1**2 / np.sqrt(self.Om*zp1**3 + self.Ok*zp1**2 + self.Ol)

    def H(self, z):
        """ Hubble parameter (km/s/Mpc) at redshift `z`.

        Parameters
        ----------
        z : array_like
          Input redshifts.

        Returns
        -------
        H : ndarray, or float if input scalar
          Hubble parameter in km/s/Mpc at each input redshift.
        """
        return self.H0 * self._efunc(z)

    def scale_factor(self, z):
        """ Scale factor at redshift `z`.

        The scale factor is defined as `a = 1 / (1 + z)`.

        Parameters
        ----------
        z : array_like
          Input redshifts.

        Returns
        -------
        a : ndarray, or float if input scalar
          Scale factor at each input redshift.
        """

        if isiterable(z):
            z = np.asarray(z)

        return 1. / (1. + z)

    def lookback_time(self, z):
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
          Lookback time in Gyr to each input redshift.
        """
        from scipy.integrate import quad
        if not isiterable(z):
            return self.hubble_time * quad(self._tfunc, 0, z)[0]

        out = np.array([quad(self._tfunc, 0, redshift)[0] for redshift in z])
        return self.hubble_time * np.array(out)

    def age(self, z):
        """ Age of the universe in Gyr at redshift `z`.

        Parameters
        ----------
        z : array_like
          Input redshifts.

        Returns
        -------
        t : ndarray, or float if input scalar
          The age of the universe in Gyr at each input redshift.
        """
        from scipy.integrate import quad
        if not isiterable(z):
            return self.hubble_time * quad(self._tfunc, z, np.inf)[0]

        out = [quad(self._tfunc, redshift, np.inf)[0] for redshift in z]
        return self.hubble_time * np.array(out)

    def critical_density(self, z):
        """ Critical density in grams per cubic cm at redshift `z`.

        Parameters
        ----------
        z : array_like
          Input redshifts.

        Returns
        -------
        rho : ndarray, or float if input scalar
          Critical density in g/cm^3 at each input redshift.
        """
        return self.critical_density0 * (self._efunc(z))**2

    def comoving_distance(self, z):
        """ Comoving line-of-sight distance in Mpc at a given
        redshift.

        The comoving distance along the line-of-sight between two
        objects remains constant with time for objects in the Hubble
        flow.

        Parameters
        ----------
        z : array_like
          Input redshifts.

        Returns
        -------
        d : ndarray, or float if input scalar
          Comoving distance in Mpc to each input redshift.
        """
        from scipy.integrate import quad
        if not isiterable(z):
            return self.hubble_distance * quad(self._inv_efunc, 0, z)[0]

        out = [quad(self._inv_efunc, 0, redshift)[0] for redshift in z]
        return self.hubble_distance * np.array(out)

    def comoving_transverse_distance(self, z):
        """ Comoving transverse distance in Mpc at a given redshift.

        This value is the transverse comoving distance at redshift `z`
        corresponding to an angular separation of 1 radian. This is
        the same as the comoving distance if omega_k is zero (as in
        the current concordance lambda CDM model).

        Parameters
        ----------
        z : array_like
          Input redshifts.

        Returns
        -------
        d : ndarray, or float if input scalar
          Comoving transverse distance in Mpc at each input redshift.

        Notes
        -----
        This quantity also called the 'proper motion distance' in some
        texts.
        """
        Ok = self.Ok
        dc = self.comoving_distance(z)
        if Ok == 0:
            return dc
        sqrtOk = sqrt(abs(Ok))
        dh = self.hubble_distance
        if Ok > 0:
            return dh / sqrtOk * np.sinh(sqrtOk * dc / dh)
        else:
            return dh / sqrtOk * np.sin(sqrtOk * dc / dh)

    def angular_diameter_distance(self, z):
        """ Angular diameter distance in Mpc at a given redshift.

        This gives the proper (sometimes called 'physical') transverse
        distance corresponding to an angle of 1 radian for an object
        at redshift `z`.

        Weinberg, 1972, pp 421-424; Weedman, 1986, pp 65-67; Peebles,
        1993, pp 325-327.

        Parameters
        ----------
        z : array_like
          Input redshifts.

        Returns
        -------
        d : ndarray, or float if input scalar
          Angular diameter distance in Mpc at each input redshift.
        """
        if isiterable(z):
            z = np.asarray(z)

        return self.comoving_transverse_distance(z) / (1. + z)

    def luminosity_distance(self, z):
        """ Luminosity distance in Mpc at redshift `z`.

        This is the distance to use when converting between the
        bolometric flux from an object at redshift `z` and its
        bolometric luminosity.

        Parameters
        ----------
        z : array_like
          Input redshifts.

        Returns
        -------
        d : ndarray, or float if input scalar
          Luminosity distance in Mpc at each input redshift.

        References
        ----------
        Weinberg, 1972, pp 420-424; Weedman, 1986, pp 60-62.
        """
        if isiterable(z):
            z = np.asarray(z)

        return (1. + z) * self.comoving_transverse_distance(z)

    def angular_diameter_distance_z1z2(self, z1, z2):
        """ Angular diameter distance between objects at 2 redshifts.
        Useful for gravitational lensing.

        Parameters
        ----------
        z1, z2 : array_like, shape (N,)
          Input redshifts. z2 must be large than z1.

        Returns
        -------
        d : ndarray, shape (N,) or float if input scalar
          The angular diameter distance between each input redshift
          pair.

        Raises
        ------
        CosmologyError
          If omega_k is < 0.

        Notes
        -----
        This method only works for flat or open curvature
        (omega_k >= 0).
        """
        # does not work for negative curvature
        Ok = self.Ok
        if Ok < 0:
            raise CosmologyError('Ok must be > 0 to use this method.')

        outscalar = False
        if not isiterable(z1) and not isiterable(z2):
            outscalar = True

        z1 = np.atleast_1d(z1)
        z2 = np.atleast_1d(z2)

        if z1.size != z2.size:
            raise ValueError('z1 and z2 must be the same size.')

        if (z1 > z2).any():
            raise ValueError('z2 must greater than z1')

        # z1 < z2
        if (z2 < z1).any():
            z1, z2 = z2, z1

        dm1 = self.comoving_transverse_distance(z1)
        dm2 = self.comoving_transverse_distance(z2)
        dh_2 = self.hubble_distance**2

        out = 1. / (1. + z2) * (dm2*np.sqrt(1. + Ok*dm1**2 / dh_2) -
                                dm1*np.sqrt(1. + Ok*dm2**2 / dh_2))

        if outscalar:
            return out[0]

        return out

    def absorption_distance(self, z):
        """ Absorption distance at redshift `z`.

        This is used to calculate the number of objects with some
        cross section of absorption and number density intersecting a
        sightline along per unit redshift path.

        Parameters
        ----------
        z : array_like
          Input redshifts.

        Returns
        -------
        d : ndarray, or float if input scalar
          Absorption distance (dimensionless) at each input redshift.

        References
        ----------
        Hogg 1999 Section 11. (astro-ph/9905116)
        Bahcall, John N. and Peebles, P.J.E. 1969, ApJ, 156L, 7B
        """
        from scipy.integrate import quad
        return quad(self._xfunc, 0, z)[0]

    def distmod(self, z):
        """ Distance modulus at redshift `z`.

        The distance modulus is defined as the (apparent magnitude -
        absolute magnitude for an object at redshift `z`).

        Parameters
        ----------
        z : array_like
          Input redshifts.

        Returns
        -------
        distmod : ndarray, or float if input scalar
          Distance modulus at each input redshift.
        """
        # Remember that the luminosity distance is in Mpc
        return 5. * np.log10(self.luminosity_distance(z) * 1.e5)

    def comoving_volume(self, z):
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
        V : ndarray, or float if input scalar
          Comoving volume in Mpc^3 at each input redshift.
        """
        Ok = self.Ok
        if Ok == 0:
            return 4. / 3. * pi * self.comoving_distance(z)**3

        dh = self.hubble_distance
        dm = self.comoving_transverse_distance(z)
        term1 = 4. * pi * dh**3 / (2. * Ok)
        term2 = dm / dh * sqrt(1 + Ok * (dm / dh)**2)
        term3 = sqrt(abs(Ok)) * dm / dh

        if Ok > 0:
            return term1 * (term2 - 1. / sqrt(abs(Ok)) * np.arcsinh(term3))
        else:
            return term1 * (term2 - 1. / sqrt(abs(Ok)) * np.arcsin(term3))


# Pre-defined cosmologies. This loops over the parameter sets in the
# parameters module and creates a FLRWCosmology instance with the same
# name as the parameter set in the current module's namespace.

for key in parameters.available:
    par = getattr(parameters, key)
    cosmo = FLRWCosmology(par['H0'], par['Om'], par['Ol'], name=key)
    cosmo.__doc__ = "%s cosmology\n\n(from %s)" % (key, par['reference'])
    setattr(sys.modules[__name__], key, cosmo)

# don't leave these variables floating around in the namespace
del key, par, cosmo

#########################################################################
# The variable below contains the current cosmology used by the
# convenience functions below and by other astropy functions if no
# cosmology is explicitly given. It can be set with set_current() and
# should be accessed using get_current().
#########################################################################


def get_cosmology_from_string(arg):
    """ Return a cosmology instance from a string.
    """
    if arg == 'no_default':
        cosmo = None
    else:
        try:
            cosmo = getattr(sys.modules[__name__], arg)
        except AttributeError:
            s = "Unknown cosmology '%s'. Valid cosmologies:\n%s" % (
                    arg, parameters.available)
            raise ValueError(s)
    return cosmo

_current = get_cosmology_from_string(DEFAULT_COSMOLOGY())


def get_current():
    """ Get the current cosmology.

    If no current has been set, the WMAP7 comology is returned and a
    warning is given.

    Returns
    -------
    cosmo : `Cosmology` instance

    See Also
    --------
    set_current : sets the current cosmology
    """
    if _current is None:
        warnings.warn('No default cosmology has been specified, '
                      'using 7-year WMAP.')
        return WMAP7

    return _current


def set_current(arg):
    """ Set the current cosmology.

    Call this with an empty string ('') to get a list of the strings
    that map to available pre-defined cosmologies.

    Parameters
    ----------
    arg : str or `Cosmology` instance
      The cosmology to use.

    Notes
    -----
    **Warning:** `set_current` is the only way to change the current
    cosmology at runtime! The current cosmology can also be read from
    an option in the astropy configuration file when astropy.cosmology
    is first imported. However, any subsequent changes to the
    cosmology configuration option using
    `astropy.config.ConfigurationItem.set()` at run-time will not
    update the current cosmology.

    See Also
    --------
    get_current : returns the currently-set cosmology
    """
    global _current
    if isinstance(arg, basestring):
        _current = get_cosmology_from_string(arg)
    elif isinstance(arg, Cosmology):
        _current = arg
    else:
        raise ValueError(
            "Argument must be a string or cosmology instance. Valid strings:"
            "\n%s" % parameters.available)


# convenience functions
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
        cosmo = get_current()
    return cosmo.comoving_transverse_distance(z) * 1.e3 * arcmin_in_radians


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
        cosmo = get_current()
    return cosmo.angular_diameter_distance(z) * 1.e3 * arcmin_in_radians


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
        cosmo = get_current()
    return 1 / (cosmo.comoving_transverse_distance(z) *
                1.e3 * arcsec_in_radians)


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
        cosmo = get_current()
    return 1 / (cosmo.angular_diameter_distance(z) * 1.e3 * arcsec_in_radians)


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
        cosmo = get_current()
    return cosmo.distmod(z)
