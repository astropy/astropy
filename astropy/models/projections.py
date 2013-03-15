# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Implements sky projections defined in WCS Paper II [2]_

All angles are set and reported in deg but internally the code works and 
keeps all angles in radians. For this to work, the mechanism of setting Model's 
properties is bypassed by passing an empty parlist to Model.__init__. This also
has the effect of not creating Projection.parameters. Projection.parnames is 
created within the Projection class. 
For an example see AZP 

References
----------
.. [2] Calabretta, M.R., Greisen, E.W., 2002, A&A, 395, 1077 (Paper II)
              
"""
from __future__ import division, print_function
from .models import Model
from . import parameters
import numpy as np

projcodes = ['TAN', 'AZP', 'SZP', 'STG', 'SIN', 'ARC', 'ZPN', 'ZEA', 'AIR', 
                    'CYP', 'CEA', 'MER']

__all__ = ['Pix2Sky_AZP', 'Sky2Pix_AZP', 'Pix2Sky_CAR', 'Sky2Pix_CAR',
           'Pix2Sky_CEA', 'Sky2Pix_CEA', 'Pix2Sky_COP', 'Sky2Pix_COP',
           'Pix2Sky_CYP', 'Sky2Pix_CYP', 'Pix2Sky_MER', 'Sky2Pix_MER',
           'Pix2Sky_SIN', 'Sky2Pix_SIN', 'Pix2Sky_STG', 'Sky2Pix_STG',
           'Pix2Sky_TAN', 'Sky2Pix_TAN']

class Projection(Model):
    """
    Base class for all sky projections
    """
    def __init__(self, parnames):        
        """
        Parameters
        ----------
        parnames : list of strings
            parameter names
        """
        super(Projection, self).__init__(parnames, ndim=2, outdim=2)
        self._pdim = 1
        # the radius of the projection sphere, by which x,y are scaled
        # not sure if it's necessary to make this a user parameter
        self.r0 = 180/np.pi
        
    @property
    def ndim(self):
        return self._ndim
    
    @property
    def outdim(self):
        return self._outdim
    
    @property
    def pdim(self):
        return self._pdim
    
class Zenithal(Projection):
    """
    Base class for all Zenithal projections
    """
    def __init__(self, parnames):
        """
        Parameters
        ----------
        parnames : list of strings
            parameter names
        """
        self.phi0 = 0.
        self.theta0 = 90.
        super(Zenithal, self).__init__(parnames)
        self.has_inverse = True

    def _compute_rtheta(self):
        # Subclasses must implement this method
        raise NotImplementedError("Subclasses should implement this")
    
    def __call__(self, x, y):
        raise NotImplementedError("Subclasses should implement this")
        
class Pix2Sky_AZP(Zenithal):
    """
    AZP : Zenital perspective projection - pixel to sky
    
    """
    parnames = ['mu', 'gamma']
    def __init__(self, mu=0., gamma=0.):
        """
        Parameters
        --------------
        mu : float
            distance from point of projection to center of sphere
            in spherical radii, default is 0.
        gamma : float
            look angle in deg, default is 0.
        """
        if mu == -1:
            raise ValueError("AZP projection is not defined for mu=-1")
        # units : mu - in spherical radii, gamma - in deg
        self._mu = parameters._Parameter('mu', mu, self, 1) 
        self._gamma = parameters._Parameter('gamma', np.deg2rad(gamma), self, 1)
        # bypass initialization of self.parameters
        self.parcheck = {'mu': self.check_mu}
        super(Pix2Sky_AZP, self).__init__(parnames=self.parnames)
    
    
    def check_mu(self, val):
        if   val == -1:
            raise ValueError("AZP projection is not defined for mu=-1")
    
    def inverse(self):
        return Sky2Pix_AZP(self.mu[0], self.gamma[0])
    
    def _compute_rtheta(self, x, y):
        return np.sqrt(x**2 + y**2 * (np.cos(self._gamma))**2)
        
    def __call__(self, x, y):
        x = np.asarray(x) + 0.
        y = np.asarray(y) + 0.
        phi = np.rad2deg(np.arctan2(x / np.cos(self._gamma), -y))
        r = self._compute_rtheta(x, y)
        pho = r / (self.r0 * (self.mu + 1) + y * np.sin(self._gamma))
        psi = np.arctan2( 1, pho)
        omega = np.arcsin((pho * self.mu) / (np.sqrt(pho**2 + 1)))
        theta1 = np.rad2deg(psi - omega)
        theta2 = np.rad2deg(psi + omega) + 180
        if np.abs(self.mu) < 1:
            if theta1 < 90 and theta1 > -90:
                theta = theta1
            else:
                theta = theta2
        else:
            #theta1dif = 90 - theta1
            #theta2dif = 90 - theta2
            if theta1 < theta2:
                theta = theta1
            else:
                theta = theta2
        return phi, theta
 
class Sky2Pix_AZP(Zenithal):
    """
    AZP : Zenital perspective projection - sky to pixel
    
    """
    parnames = ['mu', 'gamma']
    def __init__(self, mu=0., gamma=0.):
        """
        Parameters
        --------------
        mu : float
            distance from point of projection to center of sphere
            in spherical radii, default is 0.
        gamma : float
            look angle in deg, default is 0.

        """
        if mu == -1:
            raise ValueError("AZP projections is not defined for mu=-1")
        self._mu = parameters._Parameter('mu', mu, self, 1)
        self._gamma = parameters._Parameter('gamma', np.deg2rad(gamma), self, 1)
        self.parcheck = {'mu': self.check_mu}
        super(Sky2Pix_AZP, self).__init__(parnames=self.parnames)       
    
    def check_mu(self, val):
        if   val == -1:
            raise ValueError("AZP projection is not defined for mu=-1")

    def inverse(self):
        return Pix2Sky_AZP(self.mu[0], self.gamma[0])
    
    def _compute_rtheta(self, phi, theta):
        rtheta = (self.r0* (self.mu + 1) * 
                  np.cos(theta)) / ((self.mu + np.sin(theta)) +
                                    np.cos(theta) * np.cos(phi) *
                                    np.tan(self._gamma))
        return rtheta
   
    def __call__(self, phi, theta):
        phi = np.deg2rad(np.asarray(phi) + 0.)
        theta = np.deg2rad(np.asarray(theta) + 0.)
        r = self._compute_rtheta(phi, theta)
        x = r * np.sin(phi)
        y = (-r * np.cos(phi)) / np.cos(self._gamma)
        return x, y
    
class Pix2Sky_TAN(Zenithal):
    """
    TAN : Gnomonic projection - pixel to sky
    
    """
    def __init__(self):
        super(Pix2Sky_TAN, self).__init__(parnames=[])
    
    def _compute_rtheta(self, x, y):
        return np.sqrt(x**2 + y**2)
    
    def inverse(self):
        return Sky2Pix_TAN() 
        
    def __call__(self, x, y):
        x = np.asarray(x) + 0.
        y = np.asarray(y) + 0.
        phi = np.rad2deg(np.arctan2(x, -y))
        rtheta = self._compute_rtheta(x, y)
        theta = np.rad2deg(np.arctan2(self.r0, rtheta))
        return phi, theta

class Sky2Pix_TAN(Zenithal):
    """
    TAN : Gnomonic Projection - sky to pixel
    
    """
    def __init__(self):
        super(Sky2Pix_TAN, self).__init__(parnames=[])
            
    def _compute_rtheta(self, theta):
        return 1/ np.tan(theta)
    
    def inverse(self):
        return Pix2Sky_TAN() 
        
    def __call__(self, phi, theta):
        phi = np.deg2rad(np.asarray(phi) + 0.)
        theta = np.deg2rad(np.asarray(theta) + 0.)
        rtheta = self._compute_rtheta(theta)
        x = np.rad2deg(rtheta * np.sin(phi))
        y = - np.rad2deg(rtheta * np.cos(phi))
        return x, y
        
class Pix2Sky_STG(Zenithal):
    """
    STG : Stereographic Projection - pixel to sky
    
    """
    def __init__(self):
        super(Pix2Sky_STG, self).__init__(parnames=[])
        
    def _compute_rtheta(self, x, y):
        return np.sqrt(x**2 + y**2)
    
    def inverse(self):
        return Sky2Pix_STG() 
        
    def __call__(self, x, y):
        x = np.asarray(x) + 0.
        y = np.asarray(y) + 0.
        phi = np.rad2deg(np.arctan2(x, -y))
        rtheta = self._compute_rtheta(x, y)
        theta = 90 - np.rad2deg(2 * np.arctan(rtheta/(2*self.r0)))
        return phi, theta

class Sky2Pix_STG(Zenithal):
    """
    STG : Stereographic Projection - sky to pixel

    """
    def __init__(self):
        super(Sky2Pix_STG, self).__init__(parnames=[])
            
    def inverse(self):
        return Pix2Sky_STG()
    
    def _compute_rtheta(self, theta):
        return (self.r0 *2 * np.cos(theta)) / (1 + np.sin(theta))
    
    def __call__(self, phi, theta):
        phi = np.deg2rad(np.asarray(phi) + 0.)
        theta = np.deg2rad(np.asarray(theta) + 0.)
        rtheta = self._compute_rtheta(theta)
        x = rtheta * np.sin(phi)
        y = - rtheta * np.cos(phi)
        return x, y
       
class Pix2Sky_SIN(Zenithal):
    """
    SIN : Slant orthographic projection - pixel to sky
    
    """
    def __init__(self):
        super(Pix2Sky_SIN, self).__init__([])
    
    def inverse(self):
        return Sky2Pix_SIN()
    
    def _compute_rtheta(self, x, y):
        return np.sqrt(x**2 + y**2)
    
    def __call__(self, x, y):
        x = np.asarray(x) + 0.
        y = np.asarray(y) + 0.
        phi = np.rad2deg(np.arctan2(x, -y))
        rtheta = self._compute_rtheta(x, y)
        theta = np.rad2deg(np.arccos(rtheta / self.r0))
        phi = np.rad2deg(np.arctan2(x, -y))
        return phi, theta
    
class Sky2Pix_SIN(Zenithal):
    """
    SIN : Slant othographic projection - sky to pixel
    
    """
    def __init__(self):
        super(Sky2Pix_SIN, self).__init__([])
        
    def inverse(self):
        return Pix2Sky_SIN()
    
    def _compute_rtheta(self, theta):
        return self.r0 * np.cos(theta)

    def __call__(self, phi, theta):
        phi = np.deg2rad(np.asarray(phi) + 0.)
        theta = np.deg2rad(np.asarray(theta) + 0.)
        rtheta = self._compute_rtheta(theta)
        x = rtheta * np.sin(phi)
        y = - rtheta * np.cos(phi)
        return x, y
    
class Cylindrical(Projection):
    """
    Base class for Cylindrical projections
    
    """
    ##TODO: define parnames
    def __init__(self, parnames):
        self.phi0 = 0.0
        self.theta0 = 0.0
        super(Cylindrical, self).__init__(parnames)
        self.has_inverse = True
        
    def inverse(self):
        raise NotImplementedError()
    
    def __call__(self, x, y):
        raise NotImplementedError()
    
class Pix2Sky_CYP(Cylindrical):
    """
    CYP : Cylindrical perspective - pixel to sky
    
    """
    parnames = ['mu', 'lam']
    def __init__(self, mu, lam):
        if mu == -lam:
            ValueError("CYP projection is not defined for mu=-lambda")
        self._mu = parameters._Parameter('mu', mu, self, 1)
        self._lam = parameters._Parameter('lam', lam, self, 1)
        super(Pix2Sky_CYP, self).__init__(parnames=self.parnames)
        
    def check_mu(self, val):
        if   val == -1:
            raise ValueError("CYP projection is not defined for mu=-lambda")
        
    def inverse(self):
        return Sky2Pix_CYP(self.mu[0], self.lam[0])
    
    def __call__(self, x, y):
        x = np.asarray(x) + 0.
        y = np.asarray(y) + 0.
        phi = np.rad2deg(x / self.lam)
        eta = y / (self.r0 * (self.mu + self.lam))
        theta = np.arctan2(eta, 1) + np.arcsin(eta * self.mu /
                                               (np.sqrt(eta**2 + 1)))
        return phi, np.rad2deg(theta)
    
class Sky2Pix_CYP(Cylindrical):
    """
    CYP : Cylindrical Perspective - sky to pixel
    
    """
    parnames = ['mu', 'lam']
    def __init__(self, mu, lam):
        if mu == -lam:
            ValueError("CYP projection is not defined for mu=-lambda")
        self._mu = parameters._Parameter('mu', mu, self, 1)
        self._lam = parameters._Parameter('lam', lam, self, 1)
        super(Sky2Pix_CYP, self).__init__(parnames=self.parnames)
        
    def check_mu(self, val):
        if   val == -1:
            raise ValueError("CYP projection is not defined for mu=-lambda")    
   
    def inverse(self):
        return Pix2Sky_CYP(self.mu[0], self.lam[0])
    
    def __call__(self, phi, theta):
        phi = np.asarray(np.deg2rad(phi)) +0.
        theta = np.asarray(np.deg2rad(theta)) +0.
        x = self._lam * phi
        y = self.r0 * ((self.mu + self.lam) / (self.mu+np.cos(theta))) * \
                                                np.sin(theta)
        return x, y
    
class Pix2Sky_CEA(Cylindrical):
    """
    CEA : Cylindrical equal area projection - pixel to sky
    
    """
    parnames = ['lam']
    def __init__(self, lam=1):
        self._lam = parameters._Parameter('lam', lam, self, 1)
        super(Pix2Sky_CEA, self).__init__(parnames=self.parnames)
        
    def inverse(self):
        return Sky2Pix_CEA(self.lam[0])
    
    def __call__(self, x, y):
        phi = np.rad2deg(x)
        lam = np.asarray(lam)
        theta = np.rad2deg(np.arcsin(self.r0 * lam * y))
        return phi, theta
    
class Sky2Pix_CEA(Cylindrical):
    """
    CEA: Cylindrical equal area projection - sky to pixel
    
    """
    parnames = ['lam']
    def __init__(self, lam=1):
        self._lam = parameters._Parameter('lam', lam, self, 1)
        super(Sky2Pix_CEA, self).__init__(parnames=self.parnames)
        
    def inverse(self):
        return Pix2Sky_CEA(self.lam[0])
    
    def __call__(self, phi, theta):
        phi = np.asarray(np.deg2rad(phi)) +0.
        theta = np.asarray(np.deg2rad(theta)) +0.
        x = phi
        y = self.r0 * np.sin(theta) / self._lam
        return x, y
    
class Pix2Sky_CAR(Cylindrical):
    """
    CAR: Plate carree projection - pixel to sky
    
    """
    def __init__(self):
        super(Pix2Sky_CAR, self).__init__([])
        
    def inverse(self):
        return Sky2Pix_CAR()
    
    def __call__(self, x, y):
        phi = np.rad2deg(x)
        theta = np.rad2deg(y)
        return phi, theta
    
class Sky2Pix_CAR(Cylindrical):
    """
    CAR: Plate carree projection - sky to pixel
    
    """
    def __init__(self):
        super(Sky2Pix_CAR, self).__init__([])

    def inverse(self):
        return Pix2Sky_CAR()
    
    def __call__(self, phi, theta):
        x = np.asarray(np.deg2rad(phi)) +0.
        y = np.asarray(np.deg2rad(theta)) +0.
        return x, y
    
class Pix2Sky_MER(Cylindrical):
    """
    MER: Mercator - pixel to sky
    
    """
    def __init__(self):
        super(Pix2Sky_MER, self).__init__([])

    def inverse(self):
        return Sky2Pix_MER()
    
    def __call__(self, x, y):
        phi = np.rad2deg(np.asarray(x)+0.)
        theta = np.rad2deg(2*np.arctan(np.e**(y*self.r0))) - 90.
        return phi, theta
    
class Sky2Pix_MER(Cylindrical):
    """
    MER: Mercator - sky to pixel
    
    """
    def __init__(self):
        super(Sky2Pix_MER, self).__init__([])
     
    def inverse(self):
        return Pix2Sky_MER()
    
    def __call__(self, phi, theta):
        x = np.asarray(np.deg2rad(phi)) +0.
        theta = np.asarray(np.deg2rad(theta)) +0.
        y = self.r0*np.log(np.tan((np.pi/2 + theta)/2))
        return x, y
    
class Conic(Projection):
    """
    Base class for conic projections
    
    """
    def __init__(self, thetaA, eta, parnames):
        """
        Parameters
        ----------
        thetaA : float
            angle in deg
        eta : float
        
        parnames : list of strings
            parameter names
                  
        """
        self._phi0 = 0.0
        theta1 = thetaA - eta
        theta2 = thetaA + eta
        if theta1 >= 90 or theta1 <= -90:
            raise ValueError("Conic projection with thetaA={0} and eta={1} is not",
                                            " defined".format(thetaA, eta))
        if theta2 >= 90 or theta2 <= -90:
            raise ValueError("Conic projection with thetaA={0} and eta={1} is not",
                                            " defined".format(thetaA, eta))
        self._theta1 = np.deg2rad(theta1)
        self._theta2 = np.deg2rad(theta2)
        self._thetaA = np.deg2rad(thetaA)
        self._eta = np.deg2rad(eta)
        super(Conic, self).__init__(parnames)
        self.has_inverse = True
        
    def inverse(self):
        raise NotImplementedError()
    
    def _compute_rtheta(self):
        # Subclasses must implement this method
        raise NotImplementedError()
    
    def __call__(self, x, y):
        raise NotImplementedError()
    
class Pix2Sky_COP(Conic):
    """
    COP: Conic perspective - pixel to sky
    
    """
    def __init__(self, thetaA, eta):
        super(Pix2Sky_COP, self).__init__(thetaA, eta, parnames=[])
        
    def _compute_rtheta(self, x, y, Y0):
        return np.sign(self._thetaA) * np.sqrt(x**2 + (Y0 - y)**2)
    
    def inverse(self):
        return Sky2Pix_COP(self._thetaA, self._eta)
    
    def __call__(self, x, y):
        Y0 = self.r0 * np.cos(self._eta) / np.tan(self._thetaA)
        rtheta = self._compute_rtheta(x, y, Y0)
        phi = np.arctan2(x / rtheta, (Y0 - y) / rtheta) / np.sin(self._thetaA)
        theta = self._thetaA + np.arctan(1/np.tan(self._thetaA) - rtheta /
                    (self.r0 * np.cos(self._eta)))
        return np.rad2deg(phi), np.rad2deg(theta)
    
class Sky2Pix_COP(Conic):
    """
    COP: Conic perspective - sky to pixel
    
    """
    def __init__(self, thetaA, eta):
        super(Sky2Pix_COP, self).__init__(thetaA, eta, parnames=[])
        
    def _compute_rtheta(self, theta):
        return self.r0 * np.cos(self._eta) * \
                    (1/np.tan(self._thetaA) - np.tan(theta - self._thetaA))
    
    def inverse(self):
        return Pix2Sky_COP(self._thetaA, self._eta)
    
    def __call__(self, phi, theta):
        phi = np.deg2rad(phi)
        theta = np.deg2rad(theta)
        rtheta = self._compute_rtheta(phi, theta)
        Y0 = self.r0 * np.cos(self._eta) / np.tan(self._thetaA)
        C = np.sin(self._thetaA)
        x = rtheta * np.sin(C)
        y = -rtheta * np.cos(C) + Y0
        return x, y