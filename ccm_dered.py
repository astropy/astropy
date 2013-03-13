# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import print_function, division
import numpy as np

def ccm_dered(wave, flux, EBV, R_V=3.1, version='ccm89'):
    """Deredden a flux array using the CCM parameterization.
    
    Parameters
    ----------
        wave: np.array
            wavelength in Angstroms
        flux: np.array
        EBV: float
            E(B-V) extinction
        R_V: float, optional
            defaults to standard Milky Way average of 3.1
        version: string, optional
            'ccm89' is the default Cardelli, Clayton, & Mathis (1989), but does
                include the O'Donnell (1994) parameters to match IDL astrolib.
            'gcc09' is Gordon, Cartledge, & Clayton (2009).
    
    Returns
    ----------
        dered_flux: np.array
            dereddened flux vector
        
    Notes
    ----------
    Cardelli, Clayton, & Mathis (1989) parameterization is used.
    In the optical range, the updated parameters of O'Donnell (1994) are used 
        instead of CCM.
    Function valid between 910 A and 3.3 microns, although note the original CCM
        values were derived using only >1250 A data.
    Gordon, Cartledge, & Clayton (2009) has updated UV coefficients.

    References
    ----------
    Cardelli, J. A., Clayton, G. C., & Mathis, J. S. 1989, ApJ, 345, 245
    Gordon, K. D., Cartledge, S., & Clayton, G. C. 2009, ApJ, 705, 1320
    O'Donnell, J. E. 1994, ApJ, 422, 158O

    """
    
    version = version.lower()
    if version not in ['ccm89','gcc09']:
        raise ValueError('ccm_dered: version must be ccm89 or gcc09')
    
    x = 1e4 / wave      # inverse microns
    
    if any(x < 0.3) or any(x > 11):
        raise ValueError('ccm_dered valid only for wavelengths from 910 A to '+
            '3.3 microns')
    
    a = np.zeros(x.size)
    b = np.zeros(x.size)
    
    # NIR
    valid = np.where((0.3 <= x) & (x < 1.1))
    a[valid] =  0.574 * x[valid]**1.61
    b[valid] = -0.527 * x[valid]**1.61
    
    # optical, using O'Donnell (1994) values
    valid = np.where((1.1 <= x) & (x < 3.3))
    y = x[valid] - 1.82
    coef_a = np.array([-0.505, 1.647, -0.827, -1.718, 1.137, 0.701, -0.609,
        0.104, 1.])
    coef_b = np.array([3.347, -10.805, 5.491, 11.102, -7.985, -3.989, 2.908,
        1.952, 0.])
    function_a = np.poly1d(coef_a)
    function_b = np.poly1d(coef_b)
    a[valid] = function_a(y)
    b[valid] = function_b(y)
    
    # UV
    valid = np.where((3.3 <= x) & (x < 8))
    y = x[valid]
    F_a = np.zeros(valid[0].size)
    F_b = np.zeros(valid[0].size)
    select = np.where(y >= 5.9)
    yselect = y[select] - 5.9
    F_a[select] = -0.04473 * yselect**2 - 0.009779 * yselect**3
    F_b[select] = 0.2130 * yselect**2 + 0.1207 * yselect**3
    a[valid] = 1.752 - 0.316*y - (0.104 / ((y-4.67)**2 + 0.341)) + F_a
    b[valid] = -3.090 + 1.825*y + (1.206 / ((y-4.62)**2 + 0.263)) + F_b
    
    # far-UV CCM89 extrapolation
    valid = np.where((8 <= x) & (x < 11))
    y = x[valid] - 8.
    coef_a = np.array([-0.070, 0.137, -0.628, -1.073])
    coef_b = np.array([0.374, -0.420, 4.257, 13.670])  
    function_a = np.poly1d(coef_a)
    function_b = np.poly1d(coef_b)
    a[valid] = function_a(y)
    b[valid] = function_b(y)
    
    # Overwrite UV with GCC09 version if applicable. Not an extrapolation.
    if version == 'gcc09':
        valid = np.where((3.3 <= x) & (x < 11))
        y = x[valid]
        F_a = np.zeros(valid[0].size)
        F_b = np.zeros(valid[0].size)
        select = np.where((5.9 <= y) & (y < 8))
        yselect = y[select] - 5.9
        F_a[select] = -0.110 * yselect**2 - 0.0099 * yselect**3
        F_b[select] = 0.537 * yselect**2 + 0.0530 * yselect**3
        a[valid] = 1.896 - 0.372*y - (0.0108 / ((y-4.57)**2 + 0.0422)) + F_a
        b[valid] = -3.503 + 2.057*y + (0.718 / ((y-4.59)**2 + 0.0530)) + F_b
    
    
    A_V = EBV * R_V
    A_lambda = A_V * (a + b/R_V)
    dered_flux = 10**(0.4 * A_lambda)
    
    return dered_flux