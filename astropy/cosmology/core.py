# Licensed under a 3-clause BSD style license - see LICENSE.rst
import sys
import warnings
from math import sqrt, pi, exp, log
from abc import ABCMeta, abstractmethod

import numpy as np

from ..constants.cgs import pc, G, c
from ..config import ConfigurationItem
from ..utils.misc import isiterable

import parameters

# Originally authored by Andrew Becker (becker@astro.washington.edu),
# and modified by Neil Crighton (neilcrighton@gmail.com) and Roban
# Kramer (robanhk@gmail.com).

# Many of these adapted from Hogg 1999, astro-ph/9905116
# and Linder 2003, PRL 90, 91301

__all__ = ["FLRWCosmology","LambdaCDMCosmology","FlatLambdaCDMCosmology",
           "wCDMCosmology","w0waCDMCosmology","wpwaCDMCosmology",
           "w0wzCDMCosmology","get_current","set_current","WMAP5","WMAP7"]

# Constants

# speed of light in km/s
c_kms = c * 1e-5

# Mpc in cm
Mpc = 1e6 * pc

# Mpc in km
Mpc_km = 1e-5 * Mpc

# Gyr in seconds; note these are Julian years, which are defined
#  to be exactly 365.25 days of 86400 seconds each.
Gyr = 1e9 * 365.25 * 24 * 60 * 60


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

    This is an abstract base class which must have at least
    get_w overwridden to be instantiated.

    Examples of available subclasses are LambdaCDMCosmology and wCDMCosmology.

    Attributes
    ----------
    H0 : float
      Hubble parameter at z=0 in km/s/Mpc
    Om : float
      Omega matter; matter density / critical density at z=0
    Ode : float
      Omega dark energy; dark energy density / critical density at z=0
    Ok : float
      Omega_k, the curvature density at z=0. Defined as 1 - Om - Ode
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
    for redshifts < ~10).  In many cases it will be useful to also
    overload other functions such as _efunc and _inv_efunc to
    allow more efficient computation.  
    """
    __metaclass__ = ABCMeta

    def __init__(self, H0, Om, Ode, name='FLRWCosmology'):

        # all densities are in units of the critical density
        self._Om = float(Om)
        self._Ode = float(Ode)
        self._Ok = 1.0 - self._Om - self._Ode
        self.name = name

        # Hubble parameter at z=0, km/s/Mpc
        self._H0 = float(H0)
        # H0 in s^-1
        H0_s = self.H0 / Mpc_km
        # 100 km/s/Mpc * h = H0 (so h is dimensionless)
        self.h = self.H0 / 100.
        # Hubble time in Gyr
        self.hubble_time = 1. / H0_s / Gyr
        # Hubble distance in Mpc
        self.hubble_distance = c_kms / self._H0

        # critical density at z=0 (grams per cubic cm)
        self.critical_density0 = 3. * H0_s**2 / (8. * pi * G)

    def __repr__(self):
        return "%s(H0=%.3g, Om=%.3g, Ode=%.3g, Ok=%.3g)" % \
            (self.name, self._H0, self._Om, self._Ode, self._Ok)

    #Set up a set of properties for H0, Om, Ode, Ok so that the
    # cosmology can't get into an inconsistent state (e.g., Om+Ode+Ok != 1).

    @property
    def H0(self):
        return self._H0

    @H0.setter
    def H0(self, value):
        self._H0 = float(value)
        self.h = self._H0 / 100.
        H0_s = self._H0 / Mpc_km
        self.hubble_time = 1. / H0_s / Gyr
        self.hubble_distance = c_kms / self._H0
        self.critical_density0 = 3. * H0_s**2 / (8. * pi * G)

    @property
    def Om(self):
        return self._Om

    @Om.setter
    def Om(self, value):
        self._Om = float(value)
        self._Ok = 1.0 - self._Om - self._Ode

    @property
    def Ode(self):
        return self._Ode

    @Ode.setter
    def Ode(self, value):
        self._Ode = float(value)
        self._Ok = 1.0 - self._Om - self._Ode

    @property
    def Ok(self):
        return self._Ok

    @Ok.setter
    def Ok(self, value):
        # This one is tricky because there is no unique way to update
        # Om/Ode when Ok changes -- we arbitrarily decide to update
        # Ode when the user changes Ok by hand
        self._Ok = float(value)
        self._Ode = 1.0 - self._Ok - self._Om

    @abstractmethod
    def get_w(self, z):
        """ Return the dark energy equation of state at redshift z.

        This must be overridden by subclasses.
        """
        raise NotImplementedError("get_w is not implemented")

    def _w_integral(self, ln1pz):
        """ Internal convenience function for w(z) integral"""
        
        #See Linder 2003, PRL 90, 91301 eq (5)
        #Assumes scalar input, since this should only be called
        # inside an integral

        z = exp(ln1pz)-1.0
        return 1.0 + self.get_w(z)

    def _efunc(self, z):
        """ Function used to calculate H(z), the Hubble parameter.

        Parameters
        ----------
        z : array_like
          Input redshifts.

        Returns
        -------
        A value E such that H(z) = H_0 * E

        For most popular forms of w(z), this should be overridden
        to allow more efficient computation.
        """

        # This allows for an arbitrary w(z) following eq (5) of
        # Linder 2003, PRL 90, 91301.  However, it is quite inefficient
        # because it involves a numerical integration.  Most forms of w(z)
        # in common use are chosen to make this integral analytic, and
        # subclasses should probably take advantage of that if possible
        # by overloading this function.  
        #
        # In particular, the integral that needs to be evaluated is
        #   I = \int_{0}^{\log(1+z)} d \log(1+z') [ 1 + w(z') ]
        # The code below does this numerically, but if there is an 
        # analytic expression for this integral, it is probably worth
        # your while to make that explicit.  In that case, _efunc should 
        # return sqrt( (1+z)^3 Om + (1+z) Ok + Ode * exp(3 * I) )
        # and _inv_efunc one over that expression.
        #
        # See w0waCDMCosmology for an example.

        from scipy.integrate import quad

        if isiterable(z):
            z = np.asarray(z)
        Om, Ode, Ok = self.Om, self.Ode, self.Ok
        zp1 = 1.0 + z

        if not isiterable(z) :
            ival = quad(self._w_integral,0,log(1+z))[0]
        else:
            ival = np.array([quad(self._w_integral,0,log(1+redshift))[0]
                             for redshift in z])

        return np.sqrt(zp1**2 * (Om * zp1 + Ok) + Ode * np.exp(3. * ival))

    def _inv_efunc(self, z):
        """Inverse of _efunc"""

        # For efficiency, it is probably a good idea to also overload
        # this for specific dark energy models, as is the case for _efunc
        return 1.0 / self._efunc(z)

    def _tfunc(self, z):
        """ Integrand of the lookback time.

        Eqn 30 from Hogg 1999.
        """

        if isiterable(z):
            zp1 = 1.0 + np.asarray(z)
        else:
            zp1 = 1. + z

        return 1.0 / (zp1 * self._efunc(z))

    def _xfunc(self, z):
        """ Integrand of the absorption distance.
        
        See Hogg 1999 section 11.
        """

        if isiterable(z):
            zp1 = 1.0 + np.asarray(z)
        else:
            zp1 = 1. + z
        return zp1**2 / self._efunc(z)

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
        return self._H0 * self._efunc(z)

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
            raise CosmologyError('Ok must be >= 0 to use this method.')

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
        sightline per unit redshift path.

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
        if not isiterable(z):
            return quad(self._xfunc, 0, z)[0]

        out = [quad(self._xfunc, 0, redshift)[0] for redshift in z]
        return np.array(out)

    def distmod(self, z):
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


class FlatLambdaCDMCosmology(FLRWCosmology):
    """Flat FLRW cosmology with a cosmological constant.

    Examples
    --------
    >>> from astro.cosmology import FlatLambdaCDMCosmology
    >>> cosmo = FlatLambdaCDMCosmology(H0=70, Om=0.3)

    The comoving distance in Mpc at redshift z:

    >>> dc = cosmo.comoving_distance(z)
    """

    def __init__(self, H0, Om, name='FlatLambdaCDMCosmology'):
        FLRWCosmology.__init__(self, H0, Om, 1.0-Ode, name=name)


    def __repr__(self):
        return "%s(H0=%.3g, Om=%.3g, Ode=%.3g)" % \
            (self.name, self._H0, self._Om, self._Ode)

    def get_w(self,z):
        """Returns dark energy equation of state at redshift `z`.

        Parameters
        ----------
        z : array_like
          Input redshifts.

        Returns
        -------
        w : ndarray, or float if input scalar
          Dark energy equation of state, P(z)/rho(z).
        """

        return -1.0*np.ones_like(z)
    
    @property
    def Om(self):
        return self._Om

    @Om.setter
    def Om(self, value):
        self._Om = float(value)
        self._Ode = 1.0 - self._Om

    @property
    def Ode(self):
        return self._Ode

    @Ode.setter
    def Ode(self, value):
        self._Ode = float(value)
        self._Om = 1.0 - self._Ode

    @property
    def Ok(self):
        return 0.0

    @Ok.setter
    def Ok(self):
        raise CosmologyError("Can't set Ok for FlatLambdaCDMCosmology")

    def _efunc(self, z):
        """ Function used to calculate H(z), the Hubble parameter.

        Parameters
        ----------
        z : array_like
          Input redshifts.

        Returns
        -------
        A value E such that H(z) = H_0 * E
        """

        if isiterable(z):
            z = np.asarray(z)
        Om, Ode = self.Om, self.Ode
        zp1 = 1.0 + z
        return np.sqrt(zp1**3 * Om + Ode)

    def _inv_efunc(self, z):
        """ Function used to calculate 1.0/H(z)

        Parameters
        ----------
        z : array_like
          Input redshifts.

        Returns
        -------
        A value E such that H(z) = H_0 / E
        """
        
        if isiterable(z):
            z = np.asarray(z)
        Om, Ode = self.Om, self.Ode
        zp1 = 1.0 + z
        return 1.0 / np.sqrt(zp1**3 * Om + Ode)


FLRWCosmology.register(FlatLambdaCDMCosmology)


class LambdaCDMCosmology(FLRWCosmology):
    """FLRW cosmology with a cosmological constant and curvature.

    Examples
    --------
    >>> from astro.cosmology import LambdaCDMCosmology
    >>> cosmo = LambdaCDMCosmology(H0=70, Om=0.3, Ode=0.7)

    The comoving distance in Mpc at redshift z:

    >>> dc = cosmo.comoving_distance(z)
    """

    def __init__(self, H0, Om, Ode, name='LambdaCDMCosmology'):
        FLRWCosmology.__init__(self, H0, Om, Ode, name=name)

    def get_w(self, z):
        """Returns dark energy equation of state at redshift `z`.

        Parameters
        ----------
        z : array_like
          Input redshifts.

        Returns
        -------
        w : ndarray, or float if input scalar
          Dark energy equation of state, P(z)/rho(z).
        """

        return -1.0*np.ones_like(z)
    
    def _efunc(self, z):
        """ Function used to calculate H(z), the Hubble parameter.

        Parameters
        ----------
        z : array_like
          Input redshifts.

        Returns
        -------
        A value E such that H(z) = H_0 * E
        """

        if isiterable(z):
            z = np.asarray(z)
        Om, Ode, Ok = self.Om, self.Ode, self.Ok
        zp1 = 1.0 + z
        return np.sqrt(zp1**2 * (Om * zp1 + Ok) + Ode)

    def _inv_efunc(self, z):
        """ Function used to calculate 1.0/H(z)

        Parameters
        ----------
        z : array_like
          Input redshifts.

        Returns
        -------
        A value E such that H(z) = H_0 / E
        """
        
        if isiterable(z):
            z = np.asarray(z)
        Om, Ode, Ok = self.Om, self.Ode, self.Ok
        zp1 = 1.0 + z
        return 1.0 / np.sqrt(zp1**2 * (Om * zp1 + Ok) + Ode)


FLRWCosmology.register(LambdaCDMCosmology)


class wCDMCosmology(FLRWCosmology):
    """FLRW cosmology with a constant dark energy equation of state
    and curvature.

    Examples
    --------
    >>> from astro.cosmology import wCDMCosmology
    >>> cosmo = wCDMCosmology(H0=70, Om=0.3, Ode=0.7, w=-0.9)

    The comoving distance in Mpc at redshift z:

    >>> dc = cosmo.comoving_distance(z)
    """

    def __init__(self, H0, Om, Ode, w=-1., name='wCDMCosmology'):
        FLRWCosmology.__init__(self, H0, Om, Ode, name=name)
        self.w = float(w)

    def __repr__(self):
        return "%s(H0=%.3g, Om=%.3g, Ode=%.3g, w=%.3g)" % \
            (self.name, self._H0, self._Om, self._Ode, self.w)

    def get_w(self, z):
        """Returns dark energy equation of state at redshift `z`.

        Parameters
        ----------
        z : array_like
          Input redshifts.

        Returns
        -------
        w : ndarray, or float if input scalar
          Dark energy equation of state, P(z)/rho(z).
        """

        return self.w*np.ones_like(z)
    
    def _efunc(self, z):
        """ Function used to calculate H(z), the Hubble parameter.

        Parameters
        ----------
        z : array_like
          Input redshifts.

        Returns
        -------
        A value E such that H(z) = H_0 * E
        """
        
        if isiterable(z):
            z = np.asarray(z)
        Om, Ode, Ok, w = self.Om, self.Ode, self.Ok, self.w
        zp1 = 1.0 + z
        return np.sqrt(zp1**2 * (Om * zp1 + Ok) + Ode * zp1**(3*(1+w)))

    def _inv_efunc(self, z):
        """ Function used to calculate 1.0/H(z)

        Parameters
        ----------
        z : array_like
          Input redshifts.

        Returns
        -------
        A value E such that H(z) = H_0 / E
        """
        
        if isiterable(z):
            z = np.asarray(z)
        Om, Ode, Ok, w = self.Om, self.Ode, self.Ok, self.w
        zp1 = 1.0 + z
        return 1.0 / np.sqrt(zp1**2 * (Om * zp1 + Ok) + Ode * zp1**(3*(1+w)))

FLRWCosmology.register(wCDMCosmology)

class w0waCDMCosmology(FLRWCosmology):
    """FLRW cosmology with a CPL dark energy equation of state and curvature.

    The equation for the dark energy equation of state uses the
    CPL form as described in Chevallier & Polarski Int. J. Mod. Phys.
    D10, 213 (2001) and Linder PRL 90, 91301 (2003):
    w(z) = w_0 + w_a (1-a) = w_0 + w_a z / (1+z)


    Examples
    --------
    >>> from astro.cosmology import w0waCDMCosmology
    >>> cosmo = w0waCDMCosmology(H0=70, Om=0.3, Ode=0.7, w0=-0.9, wa=0.2)

    The comoving distance in Mpc at redshift z:

    >>> dc = cosmo.comoving_distance(z)
    """

    def __init__(self, H0, Om, Ode, w0=-1., wa=0., name='w0waCDMCosmology'):
        FLRWCosmology.__init__(self, H0, Om, Ode, name=name)
        self.w0 = float(w0)
        self.wa = float(wa)

    def __repr__(self):
        return "%s(H0=%.3g, Om=%.3g, Ode=%.3g, w0=%.3g, wa=%.3g)" % \
            (self.name, self._H0, self._Om, self._Ode, self.w0, self.wa)

    def get_w(self, z):
        """Returns dark energy equation of state at redshift `z`.

        Parameters
        ----------
        z : array_like
          Input redshifts.

        Returns
        -------
        w : ndarray, or float if input scalar
          Dark energy equation of state, P(z)/rho(z).
        """

        if isiterable(z):
            z = np.asarray(z)

        return self.w0 + self.wa * z / (1.0 + z)

    def _efunc(self, z):
        """ Function used to calculate H(z), the Hubble parameter.
        
        Parameters
        ----------
        z : array_like
          Input redshifts.

        Returns
        -------
        A value E such that H(z) = H_0 * E
        """
        
        if isiterable(z):
            z = np.asarray(z)
        Om, Ode, Ok, w0, wa = self.Om, self.Ode, self.Ok, self.w0, self.wa
        zp1 = 1.0 + z

        if abs(wa) < 1e-5 :
            if abs(w0+1) < 1e-5 :
                #Cosmological constant, or at least close enough
                return np.sqrt(zp1**2 * (Om * zp1 + Ok) + Ode)
            else :
                #Dark energy constant, but not cosmological constant
                return np.sqrt(zp1**2 * (Om * zp1 + Ok) +
                               Ode * zp1**(3 * (1 + w0)))
        else :
            #General form from Linder 2003, PRL 90, 91301 in the discussion
            #after eq (7)
            return np.sqrt(zp1**2 * (Om * zp1 + Ok) +
                           Ode * zp1**(3 * (1 + w0 + wa)) *
                           exp(-3 * wa * z / zp1))
                               
    def _inv_efunc(self, z):
        """ Function used to calculate 1.0/H(z)

        Parameters
        ----------
        z : array_like
          Input redshifts.

        Returns
        -------
        A value E such that H(z) = H_0 / E
        """

        #For efficiency, don't just call _efunc for this one
        #See comments for _efunc for explanation of these formulae
        if isiterable(z):
            z = np.asarray(z)
        Om, Ode, Ok, w0, wa = self.Om, self.Ode, self.Ok, self.w0, self.wa
        zp1 = 1.0 + z

        if abs(wa) < 1e-5 :
            if abs(w0+1) < 1e-5 :
                return 1.0/np.sqrt(zp1**2 * (Om * zp1 + Ok) + Ode)
            else :
                return 1.0/np.sqrt(zp1**2 * (Om * zp1 + Ok) +
                                   Ode * zp1**(3 * (1 + w0)))
        else :
            return 1.0/np.sqrt(zp1**2 * (Om * zp1 + Ok) +
                               Ode * zp1**(3 * (1 + w0 + wa)) *
                               exp(-3 * wa * z / zp1))


FLRWCosmology.register(w0waCDMCosmology)

class wpwaCDMCosmology(FLRWCosmology):
    """FLRW cosmology with a CPL dark energy equation of state, a pivot
    redshift, and curvature.

    The equation for the dark energy equation of state uses the
    CPL form as described in Chevallier & Polarski Int. J. Mod. Phys.
    D10, 213 (2001) and Linder PRL 90, 91301 (2003), but modified
    to have a pivot redshift as in the findings of the Dark Energy
    Task Force (Albrecht et al. arXiv:0901.0721 (2009))
    w(a) = w_p + w_a (a_p - a) = w_p + w_a( 1/(1+zp) - 1/(1+z) )

    Examples
    --------
    >>> from astro.cosmology import wpwaCDMCosmology
    >>> cosmo = wpwaCDMCosmology(H0=70,Om=0.3,Ode=0.7,wp=-0.9,wa=0.2,zp=0.4)

    The comoving distance in Mpc at redshift z:

    >>> dc = cosmo.comoving_distance(z)
    """

    def __init__(self, H0, Om, Ode, wp=-1., wa=0., zp=0, 
                 name='wpwaCDMCosmology'):
        FLRWCosmology.__init__(self, H0, Om, Ode, name=name)
        self.wp = float(wp)
        self.wa = float(wa)
        self.zp = float(zp)

    def __repr__(self):
        return "%s(H0=%.3g, Om=%.3g, Ode=%.3g, wp=%.3g, wa=%.3g, zp=%.3g)" % \
            (self.name, self._H0, self._Om, self._Ode, self.wp, self.wa,
             self.zp)

    def get_w(self, z):
        """Returns dark energy equation of state at redshift `z`.

        Parameters
        ----------
        z : array_like
          Input redshifts.

        Returns
        -------
        w : ndarray, or float if input scalar
          Dark energy equation of state, P(z)/rho(z).
        """

        if isiterable(z):
            z = np.asarray(z)

        apiv = 1.0 / (1.0 + self.zp)
        return self.wp + self.wa * (apiv - 1.0 / (1. + z))

    def _efunc(self, z):
        """ Function used to calculate H(z), the Hubble parameter.
        
        Parameters
        ----------
        z : array_like
          Input redshifts.

        Returns
        -------
        A value E such that H(z) = H_0 * E
        """
        
        if isiterable(z):
            z = np.asarray(z)
        Om, Ode, Ok, wp, wa = self.Om, self.Ode, self.Ok, self.wp, self.wa
        zp1 = 1.0 + z

        if abs(wa) < 1e-5 :
            if abs(wp+1) < 1e-5 :
                #Cosmological constant, or at least close enough
                return np.sqrt(zp1**2 * (Om * zp1 + Ok) + Ode)
            else :
                #Dark energy constant, but not cosmological constant
                return np.sqrt(zp1**2 * (Om * zp1 + Ok) +
                               Ode * zp1**(3 * (1 + wp)))
        else :
            #General form from Linder 2003, PRL 90, 91301 in the discussion
            #after eq (7), but modified for w(z) = w_p + w_a (a_p - a)
            apiv = 1.0 / (1.0 + self.zp)
            return np.sqrt(zp1**2 * (Om * zp1 + Ok) +
                           Ode * zp1**(3 * (1 + wp + apiv*wa)) *
                           exp(-3 * wa * z / zp1))
                               
    def _inv_efunc(self, z):
        """ Function used to calculate 1.0/H(z)

        Parameters
        ----------
        z : array_like
          Input redshifts.

        Returns
        -------
        A value E such that H(z) = H_0 / E
        """

        #For efficiency, don't just call _efunc for this one
        #See comments for _efunc for explanation of these formulae
        if isiterable(z):
            z = np.asarray(z)
        Om, Ode, Ok, wp, wa = self.Om, self.Ode, self.Ok, self.wp, self.wa
        zp1 = 1.0 + z

        if abs(wa) < 1e-5 :
            if abs(wp+1) < 1e-5 :
                return 1.0/np.sqrt(zp1**2 * (Om * zp1 + Ok) + Ode)
            else :
                return 1.0/np.sqrt(zp1**2 * (Om * zp1 + Ok) +
                                   Ode * zp1**(3 * (1 + wp)))
        else :
            apiv = 1.0 / (1.0 + self.zp)
            return 1.0/np.sqrt(zp1**2 * (Om * zp1 + Ok) +
                               Ode * zp1**(3 * (1 + wp + apiv*wa)) *
                               exp(-3 * wa * z / zp1))


FLRWCosmology.register(wpwaCDMCosmology)


class w0wzCDMCosmology(FLRWCosmology):
    """FLRW cosmology with a variable dark energy equation of state
    and curvature.

    The equation for the dark energy equation of state uses the
    simple form:
    w(z) = w_0 + w_z * z

    This form is not recommended for z > 1.


    Examples
    --------
    >>> from astro.cosmology import wawzCDMCosmology
    >>> cosmo = wawzCDMCosmology(H0=70, Om=0.3, Ode=0.7, w0=-0.9, wz=0.2)

    The comoving distance in Mpc at redshift z:

    >>> dc = cosmo.comoving_distance(z)
    """

    def __init__(self, H0, Om, Ode, w0=-1., wz=0., name='w0wzCDMCosmology'):
        FLRWCosmology.__init__(self, H0, Om, Ode, name=name)
        self.w0 = float(w0)
        self.wz = float(wz)

    def __repr__(self):
        return "%s(H0=%.3g, Om=%.3g, Ode=%.3g, w0=%.3g, wz=%.3g)" % \
            (self.name, self._H0, self._Om, self._Ode, self.w0, self.wz)

    def get_w(self, z):
        """Returns dark energy equation of state at redshift `z`.

        Parameters
        ----------
        z : array_like
          Input redshifts.

        Returns
        -------
        w : ndarray, or float if input scalar
          Dark energy equation of state, P(z)/rho(z).
        """

        if not isiterable(z):
            return self.w0 + self.wz * z
        else:
            z = np.asarray(z)
            return self.wz * z + self.w0


    def _efunc(self, z):
        """ Function used to calculate H(z), the Hubble parameter.
        
        Parameters
        ----------
        z : array_like
          Input redshifts.

        Returns
        -------
        A value E such that H(z) = H_0 * E
        """
        
        if isiterable(z):
            z = np.asarray(z)
        Om, Ode, Ok, w0, wz = self.Om, self.Ode, self.Ok, self.w0, self.wz
        zp1 = 1.0 + z

        if abs(wz) < 1e-5 :
            if abs(w0+1) < 1e-5 :
                #Cosmological constant, or at least close enough
                return np.sqrt(zp1**2 * (Om * zp1 + Ok) + Ode)
            else :
                #Dark energy constant, but not cosmological constant
                return np.sqrt(zp1**2 * (Om * zp1 + Ok) +
                               Ode * zp1**(3 * (1 + w0)))
        else :
            #General form from Linder 2003, PRL 90, 91301 in the discussion
            #after eq (5)
            return np.sqrt(zp1**2 * (Om * zp1 + Ok) +
                           Ode * zp1**(3 * (1 + w0 - wz)) *
                           exp(-3 * wz * z))
                               
    def _inv_efunc(self, z):
        """ Function used to calculate 1.0/H(z)

        Parameters
        ----------
        z : array_like
          Input redshifts.

        Returns
        -------
        A value E such that H(z) = H_0 / E
        """

        #For efficiency, don't just call _efunc for this one
        #See comments for _efunc for explanation of these formulae
        if isiterable(z):
            z = np.asarray(z)
        Om, Ode, Ok, w0, wz = self.Om, self.Ode, self.Ok, self.w0, self.wz
        zp1 = 1.0 + z

        if abs(wz) < 1e-5 :
            if abs(w0+1) < 1e-5 :
                return 1.0/np.sqrt(zp1**2 * (Om * zp1 + Ok) + Ode)
            else :
                return 1.0/np.sqrt(zp1**2 * (Om * zp1 + Ok) +
                                   Ode * zp1**(3 * (1 + w0)))
        else :
            return 1.0/np.sqrt(zp1**2 * (Om * zp1 + Ok) +
                               Ode * zp1**(3 * (1 + w0 - wz)) *
                               exp(-3 * wz * z))


FLRWCosmology.register(w0wzCDMCosmology)


# Pre-defined cosmologies. This loops over the parameter sets in the
# parameters module and creates a LambdaCDMCosmology instance with the same
# name as the parameter set in the current module's namespace.
# Note this assumes all the cosmologies in parameters are LambdaCDM,
# which is true at least as of this writing.        

for key in parameters.available:
    par = getattr(parameters, key)
    cosmo = LambdaCDMCosmology(par['H0'], par['Om'], par['Ode'], name=key)
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


def set_current(cosmo):
    """ Set the current cosmology.

    Call this with an empty string ('') to get a list of the strings
    that map to available pre-defined cosmologies.

    .. warning::
        `set_current` is the only way to change the current cosmology at
        runtime! The current cosmology can also be read from an option
        in the astropy configuration file when astropy.cosmology is first
        imported. However, any subsequent changes to the cosmology
        configuration option using `ConfigurationItem.set
        <astropy.config.configuration.ConfigurationItem.set>` at run-time
        will not update the current cosmology.

    Parameters
    ----------
    cosmo : str or `Cosmology` instance
      The cosmology to use.



    See Also
    --------
    get_current : returns the currently-set cosmology
    """
    global _current
    if isinstance(cosmo, basestring):
        _current = get_cosmology_from_string(cosmo)
    elif isinstance(cosmo, Cosmology):
        _current = cosmo
    else:
        raise ValueError(
            "Argument must be a string or cosmology instance. Valid strings:"
            "\n%s" % parameters.available)
