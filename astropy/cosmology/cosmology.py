import sys
from math import sqrt, log10, sin, sinh, asin, asinh, pi

import numpy as np
from scipy import integrate

from ..constants.cgs import pc, G, c
from parameters import WMAP5, WMAP7

# Originally authored by Andrew Becker (becker@astro.washington.edu),
# and modified by Neil Crighton (neilcrighton@gmail.com) and Roban
# Kramer (robanhk@gmail.com).

# Many of these adapted from
# astro-ph/9905116

# Constants

# speed of light in km/s
c_kms  = c * 1e-5

# Mpc in cm
Mpc = 1e6 * pc

# Mpc in km
Mpc_km = 1e-5 * Mpc

# Gyr in seconds
Gyr = 1e9 * 365.25 * 24 * 60 * 60

arcsec_in_radians = 1 / 3600. * pi / 180
arcmin_in_radians = 1 / 60. * pi / 180

class Cosmology(object):
    """
    Creating a new Cosmology instance
    ---------------------------------

    Create a new cosmology object with these three keywords:

    H0:  Hubble parameter at z=0 in km/s/Mpc
    Om:  Omega matter; matter density / critical density at z=0
    Ol:  Omega lambda; dark energy density / critical density at z=0

    or with a single name, either 'wmap7' or 'wmap5'.

    Variables derived from H0, Om, and Ol
    -------------------------------------
    
    Ok:  Omega_k, the curvature density at z=0. Defined as 1 - Om - Ol
    h:  Dimensionless Hubble parameter (H0 = 100*h km/s/Mpc).
        Often used to quote cosmological values independent of H0.
    hubble_time:    Hubble time in Gyr
    hubble_distance:   Hubble distance in Mpc
    critical_density0:  Critical density in g cm^-3 at z=0
    
    Methods
    -------

    H(z):
      Hubble parameter at redshift z (km/s/Mpc)
    critical_density(z):
      The critical density such that the universe is flat at redshift z
    scale_factor(z):
      Scale factor at redshift z
    lookback_time(z):
      Lookback time to redshift z (Gyr)
    age(z):
      Age of the universe at redshift z (Gyr)
    comoving_distance(z):
      Line of sight comoving distance to z (Mpc)
    comoving_transverse_distance(z):
      Transverse comoving distance to z (Mpc)
    angular_diameter_distance(z):
      Angular diameter distance (Mpc)
    luminosity_distance(z):
      Luminosity distance (Mpc)
    angular_diameter_distance_z1z2(z1, z2):
      Angular diameter distance between objects at z1 and z2
    absorption_distance(z):
      Absorption distance corresponding to redshift z
    distmod(z):
      Distance modulus
    comoving_volume(z):
      Comoving volume at redshift z

    Note the energy density from radiation, omega_r, is ignored (valid
    for redshifts < ~10).

    Examples
    --------
    from astro.cosmology import Cosmology
    cosmo = Cosmology(H0=70, Om=0.3, Ol=0.7)

    # get comoving distance in Mpc at redshift z
    dc = cosmo.comoving_distance(z)
    """
    def __init__(self, name=None, H0=None, Om=None, Ol=None):

        if name in (None, 'wmap7'):
            self.Om = WMAP7['Om']
            self.Ol = WMAP7['Ol']
            self.H0 = WMAP7['H0']
        elif name == 'wmap5':
            self.Om = WMAP5['Om']
            self.Ol = WMAP5['Ol']
            self.H0 = WMAP5['H0']
        else:
            s = "Cosmology must be one of ('wmap5', 'wmap7')"
            raise ValueError(s)

        # all densities are in units of the critical density
        if Om is not None:
            self.Om = float(Om)
        if Ol is not None:
            self.Ol = float(Ol)
        Ok = 1 - self.Om - self.Ol
        if abs(Ok) < 1e-5:
            Ok = 0
        self.Ok = Ok

        # Hubble parameter at z=0, km/s/Mpc
        if H0 is not None:
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
        s = "Cosmology(H0=%.3g, Om=%.3g, Ol=%.3g, Ok=%.3g)" % (
            self.H0, self.Om, self.Ol, self.Ok)
        return s

    def _efunc(self, z):
        """ Function for integration. Eqn 14 from Hogg."""
        zp1 = 1. + z
        return sqrt(self.Om*zp1**3 + self.Ok*zp1**2 + self.Ol)

    def _inv_efunc(self, z):
        """ integrand of the comoving distance.
        """
        zp1 = 1. + z
        return 1. / sqrt(self.Om*zp1**3 + self.Ok*zp1**2 + self.Ol)

    def _tfunc(self, z):
        """ Integrand of the lookback time.

        Eqn 30 from Hogg."""
        zp1 = 1. + z
        return 1. / (zp1*sqrt(self.Om*zp1**3 + self.Ok*zp1**2 + self.Ol))

    def _xfunc(self, z):
        """ Integrand of the absorption distance.

        See Hogg 1999 section 11.
        """
        zp1 = 1. + z
        return zp1**2 / sqrt(self.Om*zp1**3 + self.Ok*zp1**2 + self.Ol)

    def H(self, z):
        """ Hubble parameter (km/s/Mpc) at redshift `z`. """
        return self.H0 * self._efunc(z)

    def scale_factor(self, z):
        """ Scale factor at redshift `z`.

        The scale factor `a = 1 / (1 + z)`. """
        return 1. / (1. + z)

    def lookback_time(self, z):
        """ Lookback time in Gyr to a given redshift.

        The lookback time is the difference between the age of the
        Universe now and the age at redshift `z`.
        """
        return self.hubble_time * integrate.quad(self._tfunc, 0, z)[0]

    def age(self, z):
        """ Age of the universe in Gyr at a given redshift. """
        return self.hubble_time * integrate.quad(self._tfunc, z, np.inf)[0]

    def critical_density(self, z):
        """ Critical density in grams per cubic cm at redshift `z`.
        """
        return self.critical_density0 * (self._efunc(z))**2

    def comoving_distance(self, z):
        """ Comoving line-of-sight distance in Mpc at a given
        redshift.

        The comoving distance along the line-of-sight between two
        objects remains constant with time for objects in the Hubble
        flow.
        """
        return self.hubble_distance * integrate.quad(self._inv_efunc, 0, z)[0]

    def comoving_transverse_distance(self, z):
        """ Comoving transverse distance in Mpc at a given redshift.

        This value is the transverse comoving distance at redshift `z`
        corresponding to an angular separation of 1 radian.

        This is the same as the comoving distance if omega_k is zero
        (as in the current concordance lambda CDM model).
        """
        Ok = self.Ok
        dc  = self.comoving_distance(z)
        if Ok == 0:
            return dc
        sqrtOk = sqrt(abs(Ok))
        dh  = self.hubble_distance
        if Ok > 0:
            return dh / sqrtOk * sinh(sqrtOk * dc / dh)
        else:
            return dh / sqrtOk * sin(sqrtOk * dc / dh)

    def angular_diameter_distance(self, z):
        """ Angular diameter distance in Mpc at a given redshift.

        This gives the proper (sometimes called 'physical') transverse
        distance corresponding to an angle of 1 radian for an object
        at redshift `z`.

        Weinberg, 1972, pp 421-424; Weedman, 1986, pp 65-67; Peebles,
        1993, pp 325-327.
        """
        return self.comoving_transverse_distance(z) / (1. + z)

    def luminosity_distance(self, z):
        """ Luminosity distance in Mpc at redshift `z`.

        This is the distance to use when converting between the
        bolometric flux from an object at redshift `z` and its
        bolometric luminosity.

        Weinberg, 1972, pp 420-424; Weedman, 1986, pp 60-62.
        """
        return (1. + z) * self.comoving_transverse_distance(z)

    def angular_diameter_distance_z1z2(self, z1, z2):
        """ Angular diameter distance between objects at 2 redshifts.
        Useful for gravitational lensing.

        This method only works for flat or open curvature
        (omega_k >= 0).
        """
        # does not work for negative curvature
        Ok  = self.Ok
        assert(Ok) >= 0

        # z1 < z2
        if (z2 < z1):
            z1, z2 = z2, z1

        dm1 = self.comoving_transverse_distance(z1)
        dm2 = self.comoving_transverse_distance(z2)
        dh_2  = self.hubble_distance**2

        return 1. / (1. + z2) * (dm2*sqrt(1. + Ok*dm1**2 / dh_2) -
                                 dm1*sqrt(1. + Ok*dm2**2 / dh_2))

    def absorption_distance(self, z):
        """ Absorption distance at redshift `z`.

        This quantity is dimensionless (despite being called a
        distance)."""
        return integrate.quad(self._xfunc, 0, z)[0]

    def distmod(self, z):
        """ Distance modulus at redshift `z`.

        The distance modulus is defined as the (apparent magnitude -
        absolute magnitude for an object at redshift `z`)."""
        # Remember that the luminosity distance is in Mpc
        return 5. * log10(self.luminosity_distance(z) * 1.e5)

    def comoving_volume(self, z):
        """ Comoving volume in cubic Mpc at redshift `z`.

        This is the volume of the universe encompassed by redshifts
        less than `z`. For the case of omega_k = 0 it is a sphere of
        radius `comoving_distance(z)` but it is less intuitive if
        omega_k is not 0.
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
            return term1 * (term2 - 1. / sqrt(abs(Ok)) * asinh(term3))
        else:
            return term1 * (term2 - 1. / sqrt(abs(Ok)) * asin(term3))

# convenience functions
def kpc_comoving_per_arcmin(z, cosmo=None):
    """ Separation in transverse comoving kpc corresponding to an
    arcminute at redshift `z`.
    """
    if cosmo is None:
        cosmo = Cosmology()
    return cosmo.comoving_transverse_distance(z) * 1.e3 * arcmin_in_radians

def kpc_proper_per_arcmin(z, cosmo=None):
    """ Separation in transverse proper kpc corresponding to an
    arcminute at redshift `z`.
    """
    if cosmo is None:
        cosmo = Cosmology()
    return cosmo.angular_diameter_distance(z) * 1.e3 * arcmin_in_radians

def arcsec_per_kpc_comoving(z, cosmo=None):
    """ Angular separation in arcsec corresponding to a comoving
    transverse separation in kpc at redshift `z`.
    """
    if cosmo is None:
        cosmo = Cosmology()
    return 1 / (cosmo.comoving_transverse_distance(z) *
                 1.e3 * arcsec_in_radians)

def arcsec_per_kpc_proper(z, cosmo=None):
    """ Angular separation in arcsec corresponding to a proper
    transverse separation in kpc at redshift `z`.
    """
    if cosmo is None:
        cosmo = Cosmology()
    return 1 / (cosmo.angular_diameter_distance(z) * 1.e3 * arcsec_in_radians)

def distmod(z, cosmo=None):
    """ Distance modulus at redshift `z`.
    """
    if cosmo is None:
        cosmo = Cosmology()
    return cosmo.distmod(z)

def to_xyz(ra, dec, r, deg2rad=np.pi/180.):
    """ Convert a ra, dec and comoving distance to 3d comoving
    coordinates.

    Convert a vector pointing from the origin with comoving length
    r in the direction ra,dec (in degrees) to comoving 3-d
    coordinates x,y,z.

    Assumes universe is flat.

    ra = 0 corresponds to the positive x-z half plane. dec = 0
    corresponds to the whole x-y plane. If this is confusing, draw a
    picture :)

    To align the axes such that a vector r with direction (ra0, dec0)
    points along the positive x-axis, use ra-ra0 and dec-dec0 instead
    of ra and dec. In this case increasing dec corresponds to
    increasing z, and increasing ra (from ra=0) corresponds to
    increasing y.
    """
    # radians
    ra1 = deg2rad * ra
    dec1 = deg2rad * dec

    cos_dec1 = np.cos(dec1)
    X = r * (cos_dec1 * np.cos(ra1))
    Y = r * (cos_dec1 * np.sin(ra1))
    Z = r * np.sin(dec1)

    return X,Y,Z
