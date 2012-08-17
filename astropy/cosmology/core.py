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

__all__ = ["FLRW", "LambdaCDM", "wCDM", "w0waCDM", "wpwaCDM", 
           "w0wzCDM","get_current", "set_current", "WMAP5", "WMAP7"]

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

class FLRW(Cosmology):
    """ A class describing an isotropic and homogeneous
    (Friedmann-Lemaitre-Robertson-Walker) cosmology.

    This is an abstract base class -- you can't instantiate
    examples of this class, but must work with one of its
    subclasses such as `LambdaCDM` or `wCDM`.

    Attributes
    ----------
    H0 : float
      Hubble parameter at z=0 in km/s/Mpc
    Om0 : float
      Omega matter; matter density / critical density at z=0
    Ode0 : float
      Omega dark energy; dark energy density / critical density at z=0
    Ok0 : float
      Omega_k, the curvature density at z=0. Defined as 1 - Om0 - Ode0
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
    Class instances are static -- you can't change the values
    of the parameters.  That is, all of the attributes above are
    read only.

    The energy density from radiation, omega_r, is ignored (valid
    for redshifts < ~10).  
    """
    __metaclass__ = ABCMeta

    def __init__(self, H0, Om0, Ode0, name='FLRW'):

        # all densities are in units of the critical density
        self._Om0 = float(Om0)
        self._Ode0 = float(Ode0)
        self._Ok0 = 1.0 - self._Om0 - self._Ode0
        self.name = name

        # Hubble parameter at z=0, km/s/Mpc
        self._H0 = float(H0)
        # H0 in s^-1
        H0_s = self._H0 / Mpc_km
        # 100 km/s/Mpc * h = H0 (so h is dimensionless)
        self._h = self._H0 / 100.
        # Hubble time in Gyr
        self._hubble_time = 1. / H0_s / Gyr
        # Hubble distance in Mpc
        self._hubble_distance = c_kms / self._H0

        # critical density at z=0 (grams per cubic cm)
        self._critical_density0 = 3. * H0_s**2 / (8. * pi * G)

    def __repr__(self):
        return "%s(H0=%.3g, Om0=%.3g, Ode0=%.3g, Ok0=%.3g)" % \
            (self.name, self._H0, self._Om0, self._Ode0, self._Ok0)

    #Set up a set of properties for H0, Om, Ode, Ok for user access.
    #Note that we don't let these be set (so, obj.Om = value fails)

    @property
    def H0(self):
        return self._H0

    @property
    def Om0(self):
        return self._Om0

    @property
    def Ode0(self):
        return self._Ode0

    @property
    def Ok0(self):
        return self._Ok0

    @property
    def h(self):
        return self._h

    @property
    def hubble_time(self):
        return self._hubble_time

    @property
    def hubble_distance(self):
        return self._hubble_distance

    @property
    def critical_density0(self):
        return self._critical_density0

    @abstractmethod
    def w(self, z):
        """ Return the dark energy equation of state at redshift `z`.

        The dark energy equation of state is Pressure/density
        in units where c=1.

        This must be overridden by subclasses.
        """
        raise NotImplementedError("w(z) is not implemented")

    def Om(self, z):
        """ Return the density parameter for non-relativistic matter 
        at redshift `z`.

        Parameters
        ----------
        z : array_like
          Input redshifts.

        Returns
        -------
          The density of non-relativistic matter relative to the critical
          density at each redshift.
        """

        if isiterable(z):
            z = np.asarray(z)
        return self._Om0 * (1. + z)**3 * self.inv_efunc(z)**2

    def Ok(self, z):
        """ Return the equivalent density parameter for curvature
        at redshift `z`.

        Parameters
        ----------
        z : array_like
          Input redshifts.

        Returns
        -------
          The equivalent density parameter for curvature at each redshift.
        """

        if self._Ok0 == 0:
            #Common enough case to be worth checking
            return np.zeros_like(z)

        if isiterable(z):
            z = np.asarray(z)
        return self._Ok0 * (1. + z)**2 * self.inv_efunc(z)**2

    def Ode(self, z):
        """ Return the density parameter for non-relativistic matter 
        at redshift `z`.

        Parameters
        ----------
        z : array_like
          Input redshifts.

        Returns
        -------
          The density of non-relativistic matter relative to the critical
          density at each redshift.
        """

        if self._Ode0 == 0:
            return np.zeros_like(z)

        return self._Ode0 * self.de_density_scale(z) * self.inv_efunc(z)**2


    def _w_integrand(self, ln1pz):
        """ Internal convenience function for w(z) integral"""
        
        #See Linder 2003, PRL 90, 91301 eq (5)
        #Assumes scalar input, since this should only be called
        # inside an integral

        z = exp(ln1pz)-1.0
        return 1.0 + self.w(z)

    def de_density_scale(self, z):
        """ Evaluates the redshift dependence of the dark energy density.
        
        Parameters
        ----------
        z : array_like
          Input redshifts.

        Returns
        -------
        A value I such that :math:`\rho(z) = \rho_0 I`
        
        It will generally helpful for subclasses to overload this method if
        the integral can be done analytically for the particular dark
        energy equation of state that they implement.
        """

        # This allows for an arbitrary w(z) following eq (5) of
        # Linder 2003, PRL 90, 91301.  The equation that has to be
        # integrated is
        #
        #   I = \exp \left( 3 \int_{a}^1 \frac{ da^{\prime} }{ a^{\prime} }
        #      \left[ 1 + w\left( a^{\prime} \right) \right] \right)
        #
        # The code here does this numerically.  However, most popular
        # forms of w(z) are designed to make this integral analytic,
        # so it is probably a good idea for subclasses to overload this 
        # method if an analytic form is available.
        # 
        # The integral we actually use (the one given in Linder)
        # is rewritten in terms of z, but it's the same thing.

        from scipy.integrate import quad

        if isiterable(z):
            z = np.asarray(z)
            ival = np.array([quad(self._w_integrand,0,log(1+redshift))[0]
                             for redshift in z])
            return np.exp(3 * ival)
        else:
            ival = quad(self._w_integrand,0,log(1+z))[0]
            return exp(3 * ival)

    def efunc(self, z):
        """ Function used to calculate H(z), the Hubble parameter.

        Parameters
        ----------
        z : array_like
          Input redshifts.

        Returns
        -------
        A value E such that :math:`H(z) = H_0 E`

        It is not necessary to override this method, but if de_density_scale
        takes a particularly simple form, it may be advantageous to.
        """

        if isiterable(z):
            z = np.asarray(z)
        Om0, Ode0, Ok0 = self._Om0, self._Ode0, self._Ok0
        zp1 = 1.0 + z

        return np.sqrt(zp1**2 * (Om0 * zp1 + Ok0) +
                       Ode0 * self.de_density_scale(z))

    def inv_efunc(self, z):
        """Inverse of efunc"""

        #Avoid the function overhead by repeating code
        if isiterable(z):
            z = np.asarray(z)
        Om0, Ode0, Ok0 = self._Om0, self._Ode0, self._Ok0
        zp1 = 1.0 + z

        return 1.0/np.sqrt(zp1**2 * (Om0 * zp1 + Ok0) +
                           Ode0 * self.de_density_scale(z))

    def _tfunc(self, z):
        """ Integrand of the lookback time.

        Eqn 30 from Hogg 1999.
        """

        if isiterable(z):
            zp1 = 1.0 + np.asarray(z)
        else:
            zp1 = 1. + z

        return 1.0 / (zp1 * self.efunc(z))

    def _xfunc(self, z):
        """ Integrand of the absorption distance.
        
        See Hogg 1999 section 11.
        """

        if isiterable(z):
            zp1 = 1.0 + np.asarray(z)
        else:
            zp1 = 1. + z
        return zp1**2 / self.efunc(z)

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

        return self._H0 * self.efunc(z)

    def scale_factor(self, z):
        """ Scale factor at redshift `z`.

        The scale factor is defined as :math:`a = 1 / (1 + z)`.

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
            return self._hubble_time * quad(self._tfunc, 0, z)[0]

        out = np.array([quad(self._tfunc, 0, redshift)[0] for redshift in z])
        return self._hubble_time * np.array(out)

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
            return self._hubble_time * quad(self._tfunc, z, np.inf)[0]

        out = [quad(self._tfunc, redshift, np.inf)[0] for redshift in z]
        return self._hubble_time * np.array(out)

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

        return self._critical_density0 * (self.efunc(z))**2

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
            return self._hubble_distance * quad(self.inv_efunc, 0, z)[0]

        out = [quad(self.inv_efunc, 0, redshift)[0] for redshift in z]
        return self._hubble_distance * np.array(out)

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

        Ok0 = self._Ok0
        dc = self.comoving_distance(z)
        if Ok0 == 0:
            return dc
        sqrtOk0 = sqrt(abs(Ok0))
        dh = self._hubble_distance
        if Ok0 > 0:
            return dh / sqrtOk0 * np.sinh(sqrtOk0 * dc / dh)
        else:
            return dh / sqrtOk0 * np.sin(sqrtOk0 * dc / dh)

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
        Ok0 = self._Ok0
        if Ok0 < 0:
            raise CosmologyError('Ok0 must be >= 0 to use this method.')

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
        dh_2 = self._hubble_distance**2

        out = 1. / (1. + z2) * (dm2*np.sqrt(1. + Ok0*dm1**2 / dh_2) -
                                dm1*np.sqrt(1. + Ok0*dm2**2 / dh_2))

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

        Ok0 = self._Ok0
        if Ok0 == 0:
            return 4. / 3. * pi * self.comoving_distance(z)**3

        dh = self._hubble_distance
        dm = self.comoving_transverse_distance(z)
        term1 = 4. * pi * dh**3 / (2. * Ok0)
        term2 = dm / dh * sqrt(1 + Ok0 * (dm / dh)**2)
        term3 = sqrt(abs(Ok0)) * dm / dh

        if Ok0 > 0:
            return term1 * (term2 - 1. / sqrt(abs(Ok0)) * np.arcsinh(term3))
        else:
            return term1 * (term2 - 1. / sqrt(abs(Ok0)) * np.arcsin(term3))


class LambdaCDM(FLRW):
    """FLRW cosmology with a cosmological constant and curvature.

    This has no additional attributes beyond those of FLRW.

    Examples
    --------
    >>> from astro.cosmology import LambdaCDM
    >>> cosmo = LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)

    The comoving distance in Mpc at redshift z:

    >>> dc = cosmo.comoving_distance(z)
    """

    def __init__(self, H0, Om0, Ode0, name='LambdaCDM'):
        FLRW.__init__(self, H0, Om0, Ode0, name=name)

    def w(self, z):
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
    
    def de_density_scale(self, z):
        """ Density evolution factor for dark energy"""
        return np.ones_like(z)

    def efunc(self, z):
        """ Function used to calculate H(z), the Hubble parameter.

        Parameters
        ----------
        z : array_like
          Input redshifts.

        Returns
        -------
        A value E such that :math:`H(z) = H_0 E`
        """

        if isiterable(z):
            z = np.asarray(z)

        #We override this because it takes a particularly simple
        # form for a cosmological constant
        Om0, Ode0, Ok0 = self._Om0, self._Ode0, self._Ok0
        zp1 = 1.0 + z

        return np.sqrt(zp1**2 * (Om0 * zp1 + Ok0) + Ode0)

    def inv_efunc(self, z):
        """ Function used to calculate 1.0/H(z)

        Parameters
        ----------
        z : array_like
          Input redshifts.

        Returns
        -------
        A value E such that :math:`H(z) = H_0 / E`
        """
        
        if isiterable(z):
            z = np.asarray(z)
        Om0, Ode0, Ok0 = self._Om0, self._Ode0, self._Ok0
        zp1 = 1.0 + z

        return 1.0 / np.sqrt(zp1**2 * (Om0 * zp1 + Ok0) + Ode0)


class wCDM(FLRW):
    """FLRW cosmology with a constant dark energy equation of state
    and curvature.

    This has one additional attribute beyond those of FLRW.

    Attributes
    ----------
    w0 : float
      Dark energy equation of state (P/rho). -1 is a
      cosmological constant.

    Examples
    --------
    >>> from astro.cosmology import wCDM
    >>> cosmo = wCDM(H0=70, Om0=0.3, Ode0=0.7, w0=-0.9)

    The comoving distance in Mpc at redshift z:

    >>> dc = cosmo.comoving_distance(z)
    """

    def __init__(self, H0, Om0, Ode0, w0=-1., name='wCDM'):
        FLRW.__init__(self, H0, Om0, Ode0, name=name)
        self._w0 = float(w0)

    def __repr__(self):
        return "%s(H0=%.3g, Om0=%.3g, Ode0=%.3g, Ok0=%.3g, w0=%.3g)" % \
            (self.name, self._H0, self._Om0, 
             self._Ode0, self._Ok0, self._w0)

    @property
    def w0(self):
        return self._w0

    def w(self, z):
        """Returns dark energy equation of state at redshift `z`.

        Parameters
        ----------
        z : array_like
          Input redshifts.

        Returns
        -------
        w : ndarray, or float if input scalar
          Dark energy equation of state, :math:`w(z)=P(z)/rho(z)`
        """

        return self._w0*np.ones_like(z)
    
    def de_density_scale(self, z):
        """ Density evolution factor for dark energy"""

        if isiterable(z):
            z = np.asarray(z)
        return (1.0 + z)**(3 * (1 + self._w0))

    def efunc(self, z):
        """ Function used to calculate H(z), the Hubble parameter.

        Parameters
        ----------
        z : array_like
          Input redshifts.

        Returns
        -------
        A value E such that :math:`H(z) = H_0 E`
        """
        
        if isiterable(z):
            z = np.asarray(z)
        Om0, Ode0, Ok0, w0 = self._Om0, self._Ode, self._Ok0, self._w0
        zp1 = 1.0 + z
        return np.sqrt(zp1**2 * (Om0 * zp1 + Ok0) +
                       Ode0 * zp1**(3.0 * (1 + w0)))

    def inv_efunc(self, z):
        """ Function used to calculate 1.0/H(z)

        Parameters
        ----------
        z : array_like
          Input redshifts.

        Returns
        -------
        A value E such that :math:`H(z) = H_0 / E`
        """
        
        if isiterable(z):
            z = np.asarray(z)
        Om0, Ode0, Ok0, w0 = self._Om0, self._Ode0, self._Ok0, self._w0
        zp1 = 1.0 + z
        return 1.0 / np.sqrt(zp1**2 * (Om0 * zp1 + Ok0) + 
                             Ode0 * zp1**(3 * (1 + w0)))


class w0waCDM(FLRW):
    """FLRW cosmology with a CPL dark energy equation of state and curvature.

    The equation for the dark energy equation of state uses the
    CPL form as described in Chevallier & Polarski Int. J. Mod. Phys.
    D10, 213 (2001) and Linder PRL 90, 91301 (2003):
    :math:`w(z) = w_0 + w_a (1-a) = w_0 + w_a z / (1+z)`

    This has two additional attributes beyond those of FLRW.

    Attributes
    ----------
    w0 : float
      Dark energy equation of state (P/rho) at current epoch.
    wa : float
      Negative derivative of the dark energy equation of state with
      respect to the scale factor.

    Examples
    --------
    >>> from astro.cosmology import w0waCDM
    >>> cosmo = w0waCDM(H0=70, Om0=0.3, Ode0=0.7, w0=-0.9, wa=0.2)

    The comoving distance in Mpc at redshift z:

    >>> dc = cosmo.comoving_distance(z)
    """

    def __init__(self, H0, Om0, Ode0, w0=-1., wa=0., name='w0waCDM'):
        FLRW.__init__(self, H0, Om0, Ode0, name=name)
        self._w0 = float(w0)
        self._wa = float(wa)

    def __repr__(self):
        return "%s(H0=%.3g, Om0=%.3g, Ode0=%.3g, Ok0=%.3g, w0=%.3g, wa=%.3g)" %\
            (self.name, self._H0, self._Om0, self._Ode0, self._Ok0,
             self._w0, self._wa)

    @property
    def w0(self):
        return self._w0

    @property
    def wa(self):
        return self._wa

    def w(self, z):
        """Returns dark energy equation of state at redshift `z`.

        Parameters
        ----------
        z : array_like
          Input redshifts.

        Returns
        -------
        w : ndarray, or float if input scalar
          Dark energy equation of state, :math:`w(z) = P(z)/rho(z)`
        """

        if isiterable(z):
            z = np.asarray(z)

        return self._w0 + self._wa * z / (1.0 + z)

    def de_density_scale(self, z):
        """ Density evolution factor for dark energy"""
        if isiterable(z):
            z = np.asarray(z)
        zp1 = 1.0 + z
        return zp1**(3 * (1 + self._w0 + self._wa)) * \
            exp(-3 * self._wa * z / zp1)


class wpwaCDM(FLRW):
    """FLRW cosmology with a CPL dark energy equation of state, a pivot
    redshift, and curvature.

    The equation for the dark energy equation of state uses the
    CPL form as described in Chevallier & Polarski Int. J. Mod. Phys.
    D10, 213 (2001) and Linder PRL 90, 91301 (2003), but modified
    to have a pivot redshift as in the findings of the Dark Energy
    Task Force (Albrecht et al. arXiv:0901.0721 (2009)):
    :math:`w(a) = w_p + w_a (a_p - a) = w_p + w_a( 1/(1+zp) - 1/(1+z) )`

    This has three additional attributes beyond those of FLRW.

    Attributes
    ----------
    wp : float
      Dark energy equation of state (P/rho) at the pivot redshift
    wa : float
      Negative derivative of the dark energy equation of state with
      respect to the scale factor.
    zp : float
      Pivot redshift


    Examples
    --------
    >>> from astro.cosmology import wpwaCDM
    >>> cosmo = wpwaCDM(H0=70,Om0=0.3,Ode0=0.7,wp=-0.9,wa=0.2,zp=0.4)

    The comoving distance in Mpc at redshift z:

    >>> dc = cosmo.comoving_distance(z)
    """

    def __init__(self, H0, Om0, Ode0, wp=-1., wa=0., zp=0, 
                 name='wpwaCDM'):
        FLRW.__init__(self, H0, Om0, Ode0, name=name)
        self._wp = float(wp)
        self._wa = float(wa)
        self._zp = float(zp)

    def __repr__(self):
        str = "%s(H0=%.3g, Om0=%.3g, Ode0=%.3g, Ok0=%.3g, wp=%.3g, "+\
            "wa=%.3g, zp=%.3g)"
        return str % (self.name, self._H0, self._Om0, self._Ode0, 
                      self._Ok0, self._wp, self._wa, self._zp)

    @property
    def wp(self):
        return self._wp

    @property
    def wa(self):
        return self._wa

    @property
    def zp(self):
        return self._zp

    def w(self, z):
        """Returns dark energy equation of state at redshift `z`.

        Parameters
        ----------
        z : array_like
          Input redshifts.

        Returns
        -------
        w : ndarray, or float if input scalar
          Dark energy equation of state, :math:`w(z) = P(z)/rho(z)`
        """

        if isiterable(z):
            z = np.asarray(z)

        apiv = 1.0 / (1.0 + self._zp)
        return self._wp + self._wa * (apiv - 1.0 / (1. + z))

    def de_density_scale(self, z):
        """ Density evolution factor for dark energy"""

        if isiterable(z):
            z = np.asarray(z)
        zp1 = 1.0 + z
        apiv = 1.0 / (1.0 + self._zp)
        return zp1**(3 * (1 + self._wp + apiv * self._wa)) * \
            exp(-3 * self._wa * z / zp1)


class w0wzCDM(FLRW):
    """FLRW cosmology with a variable dark energy equation of state
    and curvature.

    The equation for the dark energy equation of state uses the
    simple form: :math:`w(z) = w_0 + w_z z`

    This form is not recommended for z > 1.

    This has two additional attributes beyond those of FLRW.

    Attributes
    ----------
    w0 : float
      Dark energy equation of state (P/rho) at current epoch.
    wz : float
      Derivative of the dark energy equation of state with respect to z.

    Examples
    --------
    >>> from astro.cosmology import wawzCDM
    >>> cosmo = wawzCDM(H0=70, Om0=0.3, Ode0=0.7, w0=-0.9, wz=0.2)

    The comoving distance in Mpc at redshift z:

    >>> dc = cosmo.comoving_distance(z)
    """

    def __init__(self, H0, Om0, Ode0, w0=-1., wz=0., name='w0wzCDM'):
        FLRW.__init__(self, H0, Om0, Ode0, name=name)
        self._w0 = float(w0)
        self._wz = float(wz)

    def __repr__(self):
        return "%s(H0=%.3g, Om0=%.3g, Ode0=%.3g, w0=%.3g, wz=%.3g)" % \
            (self.name, self._H0, self._Om0, 
             self._Ode0, self._w0, self._wz)

    @property
    def w0(self):
        return self._w0

    @property
    def wz(self):
        return self._wz

    def w(self, z):
        """Returns dark energy equation of state at redshift `z`.

        Parameters
        ----------
        z : array_like
          Input redshifts.

        Returns
        -------
        w : ndarray, or float if input scalar
          Dark energy equation of state, :math:w(z) = P(z)/rho(z)`
        """

        if isiterable(z):
            z = np.asarray(z)

        return self._w0 + self._wz * z

    def de_density_scale(self, z):
        """ Density evolution factor for dark energy"""

        if isiterable(z):
            z = np.asarray(z)
        zp1 = 1.0 + z
        return zp1**(3 * (1 + self._w0 - self._wz)) * exp(-3 * self._wz * z)

# Pre-defined cosmologies. This loops over the parameter sets in the
# parameters module and creates a LambdaCDM instance with the same
# name as the parameter set in the current module's namespace.
# Note this assumes all the cosmologies in parameters are LambdaCDM,
# which is true at least as of this writing.        

for key in parameters.available:
    par = getattr(parameters, key)
    cosmo = LambdaCDM(par['H0'], par['Om0'], par['Ode0'], name=key)
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
