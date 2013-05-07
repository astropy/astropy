# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Implements spherical rotations, defined in WCS Paper II [1]_

RotateNative2Celestial and RotateCelestial2Native follow the convention 
in WCS paper II to rotate to/from a native sphere and the celestial sphere.

The user interface uses angles in deg but internally radians are used. 
This is managed through properties. To bypass the use of Model's properties, 
an empty param_names list is passed to Model.__init__ and the properties are 
defined in the rotations classes.
      
References
----------
.. [1] Calabretta, M.R., Greisen, E.W., 2002, A&A, 395, 1077 (Paper II)

"""
from __future__ import division, print_function
import math
import numbers
import numpy as np
from .core import *
from .parameters import Parameter
from .utils import InputParameterError

__all__ = ['RotateCelestial2Native', 'RotateNative2Celestial',
           'MatrixRotation2D']

class RotateNative2Celestial(Model):
    """
    Transformation from Native to Celestial Spherical Coordinates.
    
    Defines a ZXZ rotation.
        
    Parameters
    ----------
    phi, theta, psi : float
        Euler angles in deg
    """
    param_names = ['phi', 'theta', 'psi']
    def __init__(self, phi, theta, psi):
        self._phi = Parameter('phi', np.deg2rad(phi), self, 1)
        self._theta = Parameter('theta', np.deg2rad(theta), self, 1)
        self._psi = Parameter('psi', np.deg2rad(psi), self, 1)
        super(RotateNative2Celestial, self).__init__(param_names=[], n_inputs=2, n_outputs=2)
        self.has_inverse = True
            
    @property
    def phi(self):
        return np.rad2deg(self._phi)
    
    @phi.setter
    def phi(self, val):
        self._phi = Parameter('phi', np.deg2rad(val), self)
        
    @property
    def theta(self):
        return np.rad2deg(self._theta)
    
    @theta.setter
    def theta(self, val):
        self._theta = Parameter('theta', np.deg2rad(val), self)

    @property
    def psi(self):
        return np.rad2deg(self._psi)
    
    @psi.setter
    def psi(self, val):
        self._psi = Parameter('psi', np.deg2rad(val), self)
                
    def inverse(self, phi, theta, psi):
        return RotateCelestial2Native(self.phi, self.theta, self.psi)
    
    def __call__(self, nphi, ntheta):
        nphi = np.deg2rad(nphi)
        ntheta = np.deg2rad(ntheta)
        calpha = np.rad2deg(self._phi + np.arctan2(-np.cos(ntheta) * 
                                                   np.sin(nphi - self._psi),
                        np.sin(ntheta) * np.cos(self._theta) - np.cos(ntheta)
                        * np.sin(self._theta) * np.cos(nphi - self._psi)))
        cdelta = np.rad2deg(np.arcsin(np.sin(ntheta) * np.sin(self._theta) + 
                                      np.cos(ntheta) * np.cos(self._theta) *
                                      np.cos(nphi-self._psi)))
        ind = calpha < 0
        calpha[ind] += 360
        return calpha, cdelta
    
class RotateCelestial2Native(Model):
    """
    Transformation from Celestial to Native to Spherical Coordinates.
    
    Defines a ZXZ rotation.
        
    Parameters
    ----------
    phi, theta, psi : float
        Euler angles in deg
    """
    param_names = ['phi', 'theta', 'psi']
    def __init__(self, phi, theta, psi):
        self._phi = Parameter('phi', np.deg2rad(phi), self, 1)
        self._theta = Parameter('theta', np.deg2rad(theta), self, 1)
        self._psi = Parameter('psi', np.deg2rad(psi), self, 1)
        super(RotateCelestial2Native, self).__init__(param_names=[], n_inputs=2, n_outputs=2)
        
        self.has_inverse = True
    
    @property
    def phi(self):
        return np.rad2deg(self._phi)
    
    @phi.setter
    def phi(self, val):
        self._phi = Parameter('phi', np.deg2rad(val), self)
        
    @property
    def theta(self):
        return np.rad2deg(self._theta)
    
    @theta.setter
    def theta(self, val):
        self._theta = Parameter('theta', np.deg2rad(val), self)

    @property
    def psi(self):
        return np.rad2deg(self._psi)
    
    @psi.setter
    def psi(self, val):
        self._psi = Parameter('psi', np.deg2rad(val), self)

    def inverse(self, phi, theta, psi):
        return RotateNative2Celestial(self.phi, self.theta, self.psi)
    
    def __call__(self, calpha, cdelta):
        calpha = np.deg2rad(calpha)
        cdelta = np.deg2rad(cdelta)
        nphi = np.rad2deg(self._psi + np.arctan2(-np.cos(cdelta) * 
                                                 np.sin(calpha - self._phi),
                    np.sin(cdelta) * np.cos(self._theta) - np.cos(cdelta) *
                    np.sin(self._theta) * np.cos(calpha-self._phi)))
        ntheta = np.rad2deg(np.arcsin(np.sin(cdelta) * np.sin(self._theta) +
                                      np.cos(cdelta) * np.cos(self._theta) *
                                      np.cos(calpha-self._phi)))
        ind = nphi > 180
        nphi[ind] -= 360
        return nphi, ntheta
    
class MatrixRotation2D(Model):
    """
    Perform a clockwise 2D matrix rotation given either an angle or a 
    rotation matrix.
    
    If both rotmat and angle are given, angle will be ignored.
    
    Parameters
    ----------
    rotmat : ndarray
        rotation matrix
    angle : float
        angle of rotation in deg
    """
    def __init__(self, rotmat=None, angle=None):
        if rotmat is None and angle is None:
            raise InputParameterError("Expected at least one argument - " 
                                           "a rotation matrix or an angle")
        if rotmat is not None:
            self._validate_rotmat(rotmat)
            self._rotmat = Parameter('rotmat', np.asarray(rotmat) + 0.,
                                      self, 1)
            super(MatrixRotation2D, self).__init__(param_names=['rotmat'], n_inputs=1,
                                                                                n_outputs=1, param_dim=1)
        else:
            self._validate_angle(angle)
            self._angle = Parameter('angle', np.deg2rad(angle), self, 1)
            super(MatrixRotation2D, self).__init__(param_names=[], n_inputs=1,
                                                                            n_outputs=1,param_dim=1)
            self.param_names = ['angle']
            self._rotmat = Parameter('rotmat', self._compute_matrix(angle),
                                      self, 1)
        self._n_inputs = self._rotmat[0].shape[0]
        self._n_outputs = self.n_inputs
        self._parcheck = {'rotmat': self._validate_rotmat,
                                     'angle': self._validate_angle}
    
    @property
    def angle(self):
        return Parameter('angle', np.rad2deg(self._angle), self, param_dim=1)
    
    @angle.setter
    def angle(self, val):
        self._angle = Parameter('angle', np.deg2rad(val), self, param_dim=1)
        
    def _validate_rotmat(self, rotmat):
        assert rotmat.n_inputs == 2, "Expected rotation matrix to be a 2D array"
                    
    def _validate_angle(self, angle):
        a = np.asarray(angle)
        assert isinstance(a, numbers.Number), \
                    "Expected angle to be a number"
        assert a.ndim == 0, \
                    "Expected angle to be a number"
                        
    def _compute_matrix(self, angle):
        return np.array([[math.cos(angle), math.sin(angle)], \
                        [-math.sin(angle), math.cos(angle)]], \
                                        dtype=np.float64)
    
    def inverse(self):
        nrot = np.linalg.inv(self._rotmat[0])
        return MatrixRotation2D(rotmat=nrot)
                                        
    def __call__(self, x, y):
        """
        Parameters
        ----------
        x, y : 1D array or list
              x and y coordinates
        """
        x = np.asarray(x)
        y = np.asarray(y)
        assert x.shape == y.shape
        shape = x.shape
        inarr = np.array([x.flatten(), y.flatten()], dtype=np.float64)
        assert inarr.shape[0] == 2 and inarr.ndim == 2, \
                    "Incompatible shape in MatrixRotation"
        result = np.dot(self._rotmat[0], inarr)
        x, y = result[0], result[1]
        if x.shape != shape:
            x.shape = shape
            y.shape = shape
        return x, y
    