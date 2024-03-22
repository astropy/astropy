# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module contains standard functions for earth orientation, such as
precession and nutation.

This module is (currently) not intended to be part of the public API, but
is instead primarily for internal use in `coordinates`
"""

import erfa
import numpy as np

from astropy.time import Time

from .builtin_frames.utils import get_jd12
from .matrix_utilities import matrix_transpose, rotation_matrix

jd1950 = Time("B1950").jd
jd2000 = Time("J2000").jd


def eccentricity(jd):
    """
    Eccentricity of the Earth's orbit at the requested Julian Date.

    Parameters
    ----------
    jd : scalar or array-like
        Julian date at which to compute the eccentricity

    Returns
    -------
    eccentricity : scalar or array
        The eccentricity (or array of eccentricities)

    References
    ----------
    * Explanatory Supplement to the Astronomical Almanac: P. Kenneth
      Seidelmann (ed), University Science Books (1992).
    """
    T = (jd - jd1950) / 36525.0

    p = (-0.000000126, -0.00004193, 0.01673011)

    return np.polyval(p, T)


def mean_lon_of_perigee(jd):
    """
    Computes the mean longitude of perigee of the Earth's orbit at the
    requested Julian Date.

    Parameters
    ----------
    jd : scalar or array-like
        Julian date at which to compute the mean longitude of perigee

    Returns
    -------
    mean_lon_of_perigee : scalar or array
        Mean longitude of perigee in degrees (or array of mean longitudes)

    References
    ----------
    * Explanatory Supplement to the Astronomical Almanac: P. Kenneth
      Seidelmann (ed), University Science Books (1992).
    """
    T = (jd - jd1950) / 36525.0

    p = (0.012, 1.65, 6190.67, 1015489.951)

    return np.polyval(p, T) / 3600.0


def obliquity(jd, algorithm=2006):
    """
    Computes the obliquity of the Earth at the requested Julian Date.

    Parameters
    ----------
    jd : scalar or array-like
        Julian date (TT) at which to compute the obliquity
    algorithm : int
        Year of algorithm based on IAU adoption. Can be 2006, 2000 or 1980.
        The IAU 2006 algorithm is based on Hilton et al. 2006.
        The IAU 1980 algorithm is based on the Explanatory Supplement to the
        Astronomical Almanac (1992).
        The IAU 2000 algorithm starts with the IAU 1980 algorithm and applies a
        precession-rate correction from the IAU 2000 precession model.

    Returns
    -------
    obliquity : scalar or array
        Mean obliquity in degrees (or array of obliquities)

    References
    ----------
    * Hilton, J. et al., 2006, Celest.Mech.Dyn.Astron. 94, 351
    * Capitaine, N., et al., 2003, Astron.Astrophys. 400, 1145-1154
    * Explanatory Supplement to the Astronomical Almanac: P. Kenneth
      Seidelmann (ed), University Science Books (1992).
    """
    if algorithm == 2006:
        return np.rad2deg(erfa.obl06(jd, 0))
    elif algorithm == 2000:
        return np.rad2deg(erfa.obl80(jd, 0) + erfa.pr00(jd, 0)[1])
    elif algorithm == 1980:
        return np.rad2deg(erfa.obl80(jd, 0))
    else:
        raise ValueError("invalid algorithm year for computing obliquity")


def precession_matrix_Capitaine(fromepoch, toepoch):
    """
    Computes the precession matrix from one Julian epoch to another, per IAU 2006.

    Parameters
    ----------
    fromepoch : `~astropy.time.Time`
        The epoch to precess from.
    toepoch : `~astropy.time.Time`
        The epoch to precess to.

    Returns
    -------
    pmatrix : 3x3 array
        Precession matrix to get from ``fromepoch`` to ``toepoch``

    References
    ----------
    Hilton, J. et al., 2006, Celest.Mech.Dyn.Astron. 94, 351
    """
    # Multiply the two precession matrices (without frame bias) through J2000.0
    fromepoch_to_J2000 = matrix_transpose(erfa.bp06(*get_jd12(fromepoch, "tt"))[1])
    J2000_to_toepoch = erfa.bp06(*get_jd12(toepoch, "tt"))[1]
    return J2000_to_toepoch @ fromepoch_to_J2000


def _precession_matrix_besselian(epoch1, epoch2):
    """
    Computes the precession matrix from one Besselian epoch to another using
    Newcomb's method.

    ``epoch1`` and ``epoch2`` are in Besselian year numbers.
    """
    # tropical years
    t1 = (epoch1 - 1850.0) / 1000.0
    t2 = (epoch2 - 1850.0) / 1000.0
    dt = t2 - t1

    zeta1 = 23035.545 + t1 * 139.720 + 0.060 * t1 * t1
    zeta2 = 30.240 - 0.27 * t1
    zeta3 = 17.995
    pzeta = (zeta3, zeta2, zeta1, 0)
    zeta = np.polyval(pzeta, dt) / 3600

    z1 = 23035.545 + t1 * 139.720 + 0.060 * t1 * t1
    z2 = 109.480 + 0.39 * t1
    z3 = 18.325
    pz = (z3, z2, z1, 0)
    z = np.polyval(pz, dt) / 3600

    theta1 = 20051.12 - 85.29 * t1 - 0.37 * t1 * t1
    theta2 = -42.65 - 0.37 * t1
    theta3 = -41.8
    ptheta = (theta3, theta2, theta1, 0)
    theta = np.polyval(ptheta, dt) / 3600

    return (
        rotation_matrix(-z, "z")
        @ rotation_matrix(theta, "y")
        @ rotation_matrix(-zeta, "z")
    )


def nutation_components2000B(jd):
    """
    Computes nutation components following the IAU 2000B specification.

    Parameters
    ----------
    jd : scalar
        Julian date (TT) at which to compute the nutation components

    Returns
    -------
    eps : float
        epsilon in radians
    dpsi : float
        dpsi in radians
    deps : float
        depsilon in raidans
    """
    dpsi, deps, epsa, _, _, _, _, _ = erfa.pn00b(jd, 0)
    return epsa, dpsi, deps


def nutation_matrix(epoch):
    """
    Nutation matrix generated from nutation components, IAU 2000B model.

    Matrix converts from mean coordinate to true coordinate as
    r_true = M * r_mean

    Parameters
    ----------
    epoch : `~astropy.time.Time`
        The epoch at which to compute the nutation matrix

    Returns
    -------
    nmatrix : 3x3 array
        Nutation matrix for the specified epoch

    References
    ----------
    * Explanatory Supplement to the Astronomical Almanac: P. Kenneth
      Seidelmann (ed), University Science Books (1992).
    """
    # TODO: implement higher precision 2006/2000A model if requested/needed
    return erfa.num00b(*get_jd12(epoch, "tt"))
