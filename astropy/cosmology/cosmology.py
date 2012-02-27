import sys
import numpy as np

try:
    from scipy import integrate
except ImportError:
    print ('WARNING: No Scipy found, so most cosmology functions will '
           'not work.')
from math import sqrt,log10,sin,sinh,pi

# originally by Andrew Becker; becker@astro.washington.edu, some
# changes by nhmc

# Many of these adapted from
# astro-ph/9905116

# CONSTANTS
c_kms   = 299792.458        # speed of light in km/s (exact)
pc      = 3.08567782e16     # 1 parsec in m (from wikipedia, which
                            # cites P. Kenneth Seidelmann,
                            # Ed. (1992). Explanatory Supplement to
                            # the Astronomical Almanac. Sausalito, CA:
                            # University Science Books. p. 716 and
                            # s.v. parsec in Glossary.)
Gyr = 1e9 * 365.25 * 24 * 60 * 60  # seconds per Gyr 

arcsec_in_radians = 1 / 3600. * pi / 180
arcmin_in_radians = 1 / 60. * pi / 180

class Cosmology(object):
    """
    Instantiate with optional keyword arguments (defaults shown).

    Om = 0.27:   Omega matter
    Ol = 0.73:   Omega lambda
    w  = -1.0:   pressure/density dark energy ratio
    H0 = 73:     Present day hubble constant in km/s/Mpc

    Derived values
    
    Ok:    Curvature density  (assumes Omega_total=1)
    h:     Dimensionless Hubble parameter (H0 = 100*h km/s/Mpc)
           Used to give H0-independent distances
    th:    Hubble time in Gyr
    dh:    Hubble distance in Mpc

    Methods

    Hz(z):      Hubble constant at redshift z (km/s/Mpc)
    a(z):       scale factor at redshift z
    tl(z):      Lookback time for redshift z (Gyr)
    dc(z):      Line of sight comoving distance at z (Mpc)
    dm(z):      Transverse comoving distance at z (Mpc)
    da(z):      Angular diameter distance (Mpc)
    da2(z1, z2): Angular diameter distance between objects at z1 and z2
    dl(z):      Luminosity distance (Mpc)
    distmod(z): Distance modulus
    X(z):       Absorption distance corresponding to redshift z

    Constants
    
    c_kms  = 2.9979E5:        Speed of light (km/s)
    pc = 3.08567782E16:       Parsec (metres)

    
    Note that the integration routine used, scipy.integrate.quad, will
    start giving nonsense for very large redshifts (>100?).

    Note also the energy density from radiation, Omega_r, is ignored
    (valid for redshifts < ~10).

    Examples
    --------

    from astro.cosmology import Cosmology
    cosmo = Cosmology(Om=0.3, Ol=0.7, H0=70)

    # get comoving distance at redshift z
    dc = cosmo.dc(z)

    """
    # Note all the distance units are the same as dh (Mpc)

    def __init__(self, Om=0.27, Ol=0.73, w=-1., H0=73.):
        # The defaults above imply Ok = 0 if we assume Omega total = 1
        # (which we do here).

        # all densities are in units of the critical density
        self.Om = Om                   # matter density
        self.Ol = Ol                   # lambda density
        self.Ok = 1. - Om - Ol         # curvature density
        self.w = w                     # pressure / density for dark energy
        self.H0 = H0                   # present day Hubble constant, km/s/Mpc
        self.h = H0 / 100.             # 100 km/s/Mpc * h = H0
        self.th = 1. / H0 * pc * 1.e3 / Gyr # Hubble time in s
        self.dh = c_kms / H0           # Hubble distance in Mpc

    def __repr__(self):
        return "Cosmology(Om=%(Om)s, Ol=%(Ol)s, Ok=%(Ok)s, H0=%(H0)s, w=%(w)s)" % self.__dict__

    def efunc(self, z):
        """ Function for integration. Eqn 14 from Hogg, adding 'w' to
        Omega lambda."""
        zp1 = 1. + z
        temp = 3. * (1. + self.w)
        return sqrt(self.Om*zp1**3 + self.Ok*zp1**2 + self.Ol*zp1**temp)

    def inv_efunc(self, z):
        """ Function for integration. Eqn 14 from Hogg, adding 'w' to
        Omega lambda."""
        zp1 = 1. + z
        temp = 3. * (1. + self.w)
        return 1. / sqrt(self.Om*zp1**3 + self.Ok*zp1**2
                         + self.Ol*zp1**temp)

    def tfunc(self, z):
        """ Function for integration to find lookback time. Eqn 30
        from Hogg, adding 'w' to Omega lambda."""
        zp1 = 1. + z
        temp = 3.*(1. + self.w)
        return 1. / (zp1*sqrt(self.Om*zp1**3 + self.Ok*zp1**2
                              + self.Ol*zp1**temp))

    def xfunc(self, z):
        """ Function for integration to find the absorption distance.
        """
        zp1 = 1. + z
        temp = 3.*(1. + self.w)
        return zp1**2 / sqrt(self.Om*zp1**3 + self.Ok*zp1**2
                             + self.Ol*zp1**temp)

    def Hz(self, z):
        """ Returns the Hubble constant at redshift z. """
        return self.H0 * self.efunc(z)

    def a(self, z):
        """ Returns the scale factor at z,   1 / (1+z). """
        return 1. / (1.+z)

    def tl(self, z):
        """ Returns the lookback time in Gyr; the difference
        between the age of the Universe now and the age at z."""
        return self.th * integrate.quad(self.tfunc, 0, z)[0]

    def dc(self, z):
        """ Returns the comoving distance to an object at redshift z
        in metres. Remains constant with time if the objects are in
        the Hubble flow."""
        return self.dh * integrate.quad(self.inv_efunc, 0, z)[0]

    def dm(self, z):
        """ Returns the transverse comoving distance in metres.
        dm*dtheta gives the transverse comoving distance between two
        objects separated by an angle dtheta radians at redshift
        z. Note dm = dc if Ok is zero (as in current lambda CDM
        model)."""
        Ok = self.Ok
        dc  = self.dc(z)
        if Ok == 0:
            return dc
        sqrtOk = sqrt(abs(Ok))
        dh  = self.dh
        if Ok > 0:
            return dh / sqrtOk * sinh(sqrtOk * dc / dh)
        else:
            return dh / sqrtOk * sin(sqrtOk * dc / dh)

    def da(self, z):
        """ Angular diameter distance in metres.  Ratio of an object's
        physical transverse size to its angular size in radians."""
        return self.dm(z) / (1. + z)

    def da2(self, z1, z2):
        """ Angular diameter distance in metres between objects at 2
        redshifts.  Useful for gravitational lensing."""
        # does not work for negative curvature
        assert(self.Ok) >= 0

        # z1 < z2
        if (z2 < z1):
            z1, z2 = z2, z1

        dm1 = self.dm(z1)
        dm2 = self.dm(z2)
        Ok  = self.Ok
        dh_2  = self.dh * self.dh

        return 1. / (1.+z2) * (dm2*sqrt(1. + Ok*dm1**2 / dh_2) -
                               dm1*sqrt(1. + Ok*dm2**2 / dh_2))

    def dl(self, z):
        """ Returns the luminosity distance in metres.  (Relationship
        between bolometric flux and bolometric luminosity.)"""
        return (1.+z) * self.dm(z)

    def distmod(self, z):
        """ Returns the distance modulus (apparent magnitude -
        absolute magnitude for an object at redshift z)."""
        # Remember that dl is in Mpc
        return 5. * log10(self.dl(z) * 1.e5)

    def X(self, z2, z1=0):
        """ Return the absorption distance from redshift
        z2 to z1. Dimensionless (despite being called a distance) """
        return integrate.quad(self.xfunc, z1, z2)[0]

# convenience functions
def kpc_per_arcmin(z, cosmo=None):
    """ Find the number of transverse comoving kpc corresponding to an
    arcminute at the given redshift. 
    """
    if cosmo is None:
        cosmo = Cosmology(Om=0.27, Ol=0.73, w=-1., H0=70.)
    return cosmo.dm(z) * 1.e3 * arcmin_in_radians

def kpc_proper_per_arcmin(z, cosmo=None):
    """ Find the number of transverse proper kpc corresponding to an
    arcminute at the given redshift.

    proper kpc = comoving kpc * (scale factor)

    where scale factor = 1 / (1 + redshift)
    """
    if cosmo is None:
        cosmo = Cosmology(Om=0.27, Ol=0.73, w=-1., H0=70.)
    return cosmo.dm(z) * 1.e3 * arcmin_in_radians / (1 + z)

def arcsec_per_kpc(z, cosmo=None):
    """ Find the number of arcsec corresponding to a
    comoving transverse separation in kpc at the given redshift.
    """
    if cosmo is None:
        cosmo = Cosmology(Om=0.27, Ol=0.73, w=-1., H0=70.)
    return 1 / (cosmo.dm(z) * 1.e3 * arcsec_in_radians)

def arcsec_per_kpc_proper(z, cosmo=None):
    """ Find the number of arcsec corresponding to a proper transverse
    separation in kpc at the given redshift.

    proper kpc = comoving kpc * (scale factor)

    where scale factor = 1 / (1 + redshift)
    """
    if cosmo is None:
        cosmo = Cosmology(Om=0.27, Ol=0.73, w=-1., H0=70.)
    return (1 + z) / (cosmo.dm(z) * 1.e3 * arcsec_in_radians)


def distmod(z, cosmo=None):
    """Find the distance modulus at the given redshift.
    """
    if cosmo is None:
        cosmo = Cosmology(Om=0.27, Ol=0.73, w=-1., H0=70.)
    return cosmo.distmod(z)

def to_xyz(ra, dec, r, deg2rad=np.pi/180.):
    """ Convert a comoving distance, ra and dec into 3d comoving
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

