# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module contains standard functions for earth orientation, such as
precession and nutation.

This module is (currently) not intended to be part of the public API, but
is instead primarily for internal use in `coordinates`
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np

from ..time import Time
from .. import units as u
from .matrix_utilities import rotation_matrix, matrix_product, matrix_transpose


jd1950 = Time('B1950', scale='tai').jd
jd2000 = Time('J2000', scale='utc').jd
_asecperrad = u.radian.to(u.arcsec)


def eccentricity(jd):
    """
    Eccentricity of the Earth's orbit at the requested Julian Date.

    Parameters
    ----------
    jd : scalar or array-like
        Julian date at which to compute the eccentricity

    returns
    -------
    eccentricity : scalar or array
        The eccentricity (or array of eccentricities)

    References
    ----------
    * Explanatory Supplement to the Astronomical Almanac: P. Kenneth
      Seidelmann (ed), University Science Books (1992).
    """
    T = (jd - jd1950) / 36525.0

    p = (-0.000000126, - 0.00004193, 0.01673011)

    return np.polyval(p, T)


def mean_lon_of_perigee(jd):
    """
    Computes the mean longitude of perigee of the Earth's orbit at the
    requested Julian Date.

    Parameters
    ----------
    jd : scalar or array-like
        Julian date at which to compute the mean longitude of perigee

    returns
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

    return np.polyval(p, T) / 3600.


def obliquity(jd, algorithm=2006):
    """
    Computes the obliquity of the Earth at the requested Julian Date.

    Parameters
    ----------
    jd : scalar or array-like
        Julian date at which to compute the obliquity
    algorithm : int
        Year of algorithm based on IAU adoption. Can be 2006, 2000 or 1980. The
        2006 algorithm is mentioned in Circular 179, but the canonical reference
        for the IAU adoption is apparently Hilton et al. 06 is composed of the
        1980 algorithm with a precession-rate correction due to the 2000
        precession models, and a description of the 1980 algorithm can be found
        in the Explanatory Supplement to the Astronomical Almanac.

    returns
    -------
    obliquity : scalar or array
        Mean obliquity in degrees (or array of obliquities)

    References
    ----------
    * Hilton, J. et al., 2006, Celest.Mech.Dyn.Astron. 94, 351. 2000
    * USNO Circular 179
    * Explanatory Supplement to the Astronomical Almanac: P. Kenneth
      Seidelmann (ed), University Science Books (1992).
    """
    T = (jd - jd2000) / 36525.0

    if algorithm == 2006:
        p = (-0.0000000434, -0.000000576, 0.00200340, -0.0001831, -46.836769, 84381.406)
        corr = 0
    elif algorithm == 2000:
        p = (0.001813, -0.00059, -46.8150, 84381.448)
        corr = -0.02524 * T
    elif algorithm == 1980:
        p = (0.001813, -0.00059, -46.8150, 84381.448)
        corr = 0
    else:
        raise ValueError('invalid algorithm year for computing obliquity')

    return (np.polyval(p, T) + corr) / 3600.


# TODO: replace this with SOFA equivalent
def precession_matrix_Capitaine(fromepoch, toepoch):
    """
    Computes the precession matrix from one Julian epoch to another.
    The exact method is based on Capitaine et al. 2003, which should
    match the IAU 2006 standard.

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
    USNO Circular 179
    """
    mat_fromto2000 = matrix_transpose(
        _precess_from_J2000_Capitaine(fromepoch.jyear))
    mat_2000toto = _precess_from_J2000_Capitaine(toepoch.jyear)

    return np.dot(mat_2000toto, mat_fromto2000)


def _precess_from_J2000_Capitaine(epoch):
    """
    Computes the precession matrix from J2000 to the given Julian Epoch.
    Expression from from Capitaine et al. 2003 as expressed in the USNO
    Circular 179.  This should match the IAU 2006 standard from SOFA.

    Parameters
    ----------
    epoch : scalar
        The epoch as a Julian year number (e.g. J2000 is 2000.0)

    """
    T = (epoch - 2000.0) / 100.0
    # from USNO circular
    pzeta = (-0.0000003173, -0.000005971, 0.01801828, 0.2988499, 2306.083227, 2.650545)
    pz = (-0.0000002904, -0.000028596, 0.01826837, 1.0927348, 2306.077181, -2.650545)
    ptheta = (-0.0000001274, -0.000007089, -0.04182264, -0.4294934, 2004.191903, 0)
    zeta = np.polyval(pzeta, T) / 3600.0
    z = np.polyval(pz, T) / 3600.0
    theta = np.polyval(ptheta, T) / 3600.0

    return matrix_product(rotation_matrix(-z, 'z'),
                          rotation_matrix(theta, 'y'),
                          rotation_matrix(-zeta, 'z'))


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

    return matrix_product(rotation_matrix(-z, 'z'),
                          rotation_matrix(theta, 'y'),
                          rotation_matrix(-zeta, 'z'))


def _load_nutation_data(datastr, seriestype):
    """
    Loads nutation series from data stored in string form.

    Seriestype can be 'lunisolar' or 'planetary'
    """

    if seriestype == 'lunisolar':
        dtypes = [('nl', int),
                  ('nlp', int),
                  ('nF', int),
                  ('nD', int),
                  ('nOm', int),
                  ('ps', float),
                  ('pst', float),
                  ('pc', float),
                  ('ec', float),
                  ('ect', float),
                  ('es', float)]
    elif seriestype == 'planetary':
        dtypes = [('nl', int),
                  ('nF', int),
                  ('nD', int),
                  ('nOm', int),
                  ('nme', int),
                  ('nve', int),
                  ('nea', int),
                  ('nma', int),
                  ('nju', int),
                  ('nsa', int),
                  ('nur', int),
                  ('nne', int),
                  ('npa', int),
                  ('sp', int),
                  ('cp', int),
                  ('se', int),
                  ('ce', int)]
    else:
        raise ValueError('requested invalid nutation series type')

    lines = [l for l in datastr.split('\n')
             if not l.startswith('#') if not l.strip() == '']

    lists = [[] for _ in dtypes]
    for l in lines:
        for i, e in enumerate(l.split(' ')):
            lists[i].append(dtypes[i][1](e))
    return np.rec.fromarrays(lists, names=[e[0] for e in dtypes])


_nut_data_00b = """
#l lprime F D Omega longitude_sin longitude_sin*t longitude_cos obliquity_cos obliquity_cos*t,obliquity_sin

0 0 0 0 1 -172064161.0 -174666.0 33386.0 92052331.0 9086.0 15377.0
0 0 2 -2 2 -13170906.0 -1675.0 -13696.0 5730336.0 -3015.0 -4587.0
0 0 2 0 2 -2276413.0 -234.0 2796.0 978459.0 -485.0 1374.0
0 0 0 0 2 2074554.0 207.0 -698.0 -897492.0 470.0 -291.0
0 1 0 0 0 1475877.0 -3633.0 11817.0 73871.0 -184.0 -1924.0
0 1 2 -2 2 -516821.0 1226.0 -524.0 224386.0 -677.0 -174.0
1 0 0 0 0 711159.0 73.0 -872.0 -6750.0 0.0 358.0
0 0 2 0 1 -387298.0 -367.0 380.0 200728.0 18.0 318.0
1 0 2 0 2 -301461.0 -36.0 816.0 129025.0 -63.0 367.0
0 -1 2 -2 2 215829.0 -494.0 111.0 -95929.0 299.0 132.0
0 0 2 -2 1 128227.0 137.0 181.0 -68982.0 -9.0 39.0
-1 0 2 0 2 123457.0 11.0 19.0 -53311.0 32.0 -4.0
-1 0 0 2 0 156994.0 10.0 -168.0 -1235.0 0.0 82.0
1 0 0 0 1 63110.0 63.0 27.0 -33228.0 0.0 -9.0
-1 0 0 0 1 -57976.0 -63.0 -189.0 31429.0 0.0 -75.0
-1 0 2 2 2 -59641.0 -11.0 149.0 25543.0 -11.0 66.0
1 0 2 0 1 -51613.0 -42.0 129.0 26366.0 0.0 78.0
-2 0 2 0 1 45893.0 50.0 31.0 -24236.0 -10.0 20.0
0 0 0 2 0 63384.0 11.0 -150.0 -1220.0 0.0 29.0
0 0 2 2 2 -38571.0 -1.0 158.0 16452.0 -11.0 68.0
0 -2 2 -2 2 32481.0 0.0 0.0 -13870.0 0.0 0.0
-2 0 0 2 0 -47722.0 0.0 -18.0 477.0 0.0 -25.0
2 0 2 0 2 -31046.0 -1.0 131.0 13238.0 -11.0 59.0
1 0 2 -2 2 28593.0 0.0 -1.0 -12338.0 10.0 -3.0
-1 0 2 0 1 20441.0 21.0 10.0 -10758.0 0.0 -3.0
2 0 0 0 0 29243.0 0.0 -74.0 -609.0 0.0 13.0
0 0 2 0 0 25887.0 0.0 -66.0 -550.0 0.0 11.0
0 1 0 0 1 -14053.0 -25.0 79.0 8551.0 -2.0 -45.0
-1 0 0 2 1 15164.0 10.0 11.0 -8001.0 0.0 -1.0
0 2 2 -2 2 -15794.0 72.0 -16.0 6850.0 -42.0 -5.0
0 0 -2 2 0 21783.0 0.0 13.0 -167.0 0.0 13.0
1 0 0 -2 1 -12873.0 -10.0 -37.0 6953.0 0.0 -14.0
0 -1 0 0 1 -12654.0 11.0 63.0 6415.0 0.0 26.0
-1 0 2 2 1 -10204.0 0.0 25.0 5222.0 0.0 15.0
0 2 0 0 0 16707.0 -85.0 -10.0 168.0 -1.0 10.0
1 0 2 2 2 -7691.0 0.0 44.0 3268.0 0.0 19.0
-2 0 2 0 0 -11024.0 0.0 -14.0 104.0 0.0 2.0
0 1 2 0 2 7566.0 -21.0 -11.0 -3250.0 0.0 -5.0
0 0 2 2 1 -6637.0 -11.0 25.0 3353.0 0.0 14.0
0 -1 2 0 2 -7141.0 21.0 8.0 3070.0 0.0 4.0
0 0 0 2 1 -6302.0 -11.0 2.0 3272.0 0.0 4.0
1 0 2 -2 1 5800.0 10.0 2.0 -3045.0 0.0 -1.0
2 0 2 -2 2 6443.0 0.0 -7.0 -2768.0 0.0 -4.0
-2 0 0 2 1 -5774.0 -11.0 -15.0 3041.0 0.0 -5.0
2 0 2 0 1 -5350.0 0.0 21.0 2695.0 0.0 12.0
0 -1 2 -2 1 -4752.0 -11.0 -3.0 2719.0 0.0 -3.0
0 0 0 -2 1 -4940.0 -11.0 -21.0 2720.0 0.0 -9.0
-1 -1 0 2 0 7350.0 0.0 -8.0 -51.0 0.0 4.0
2 0 0 -2 1 4065.0 0.0 6.0 -2206.0 0.0 1.0
1 0 0 2 0 6579.0 0.0 -24.0 -199.0 0.0 2.0
0 1 2 -2 1 3579.0 0.0 5.0 -1900.0 0.0 1.0
1 -1 0 0 0 4725.0 0.0 -6.0 -41.0 0.0 3.0
-2 0 2 0 2 -3075.0 0.0 -2.0 1313.0 0.0 -1.0
3 0 2 0 2 -2904.0 0.0 15.0 1233.0 0.0 7.0
0 -1 0 2 0 4348.0 0.0 -10.0 -81.0 0.0 2.0
1 -1 2 0 2 -2878.0 0.0 8.0 1232.0 0.0 4.0
0 0 0 1 0 -4230.0 0.0 5.0 -20.0 0.0 -2.0
-1 -1 2 2 2 -2819.0 0.0 7.0 1207.0 0.0 3.0
-1 0 2 0 0 -4056.0 0.0 5.0 40.0 0.0 -2.0
0 -1 2 2 2 -2647.0 0.0 11.0 1129.0 0.0 5.0
-2 0 0 0 1 -2294.0 0.0 -10.0 1266.0 0.0 -4.0
1 1 2 0 2 2481.0 0.0 -7.0 -1062.0 0.0 -3.0
2 0 0 0 1 2179.0 0.0 -2.0 -1129.0 0.0 -2.0
-1 1 0 1 0 3276.0 0.0 1.0 -9.0 0.0 0.0
1 1 0 0 0 -3389.0 0.0 5.0 35.0 0.0 -2.0
1 0 2 0 0 3339.0 0.0 -13.0 -107.0 0.0 1.0
-1 0 2 -2 1 -1987.0 0.0 -6.0 1073.0 0.0 -2.0
1 0 0 0 2 -1981.0 0.0 0.0 854.0 0.0 0.0
-1 0 0 1 0 4026.0 0.0 -353.0 -553.0 0.0 -139.0
0 0 2 1 2 1660.0 0.0 -5.0 -710.0 0.0 -2.0
-1 0 2 4 2 -1521.0 0.0 9.0 647.0 0.0 4.0
-1 1 0 1 1 1314.0 0.0 0.0 -700.0 0.0 0.0
0 -2 2 -2 1 -1283.0 0.0 0.0 672.0 0.0 0.0
1 0 2 2 1 -1331.0 0.0 8.0 663.0 0.0 4.0
-2 0 2 2 2 1383.0 0.0 -2.0 -594.0 0.0 -2.0
-1 0 0 0 2 1405.0 0.0 4.0 -610.0 0.0 2.0
1 1 2 -2 2 1290.0 0.0 0.0 -556.0 0.0 0.0
"""[1:-1]
_nut_data_00b = _load_nutation_data(_nut_data_00b, 'lunisolar')

# TODO: replace w/SOFA equivalent


def nutation_components2000B(jd):
    """
    Computes nutation components following the IAU 2000B specification

    Parameters
    ----------
    jd : scalar
        epoch at which to compute the nutation components as a JD

    Returns
    -------
    eps : float
        epsilon in radians
    dpsi : float
        dpsi in radians
    deps : float
        depsilon in raidans
    """
    epsa = np.radians(obliquity(jd, 2000))
    t = (jd - jd2000) / 36525

    # Fundamental (Delaunay) arguments from Simon et al. (1994) via SOFA
    # Mean anomaly of moon
    el = ((485868.249036 + 1717915923.2178 * t) % 1296000) / _asecperrad
    # Mean anomaly of sun
    elp = ((1287104.79305 + 129596581.0481 * t) % 1296000) / _asecperrad
    # Mean argument of the latitude of Moon
    F = ((335779.526232 + 1739527262.8478 * t) % 1296000) / _asecperrad
    # Mean elongation of the Moon from Sun
    D = ((1072260.70369 + 1602961601.2090 * t) % 1296000) / _asecperrad
    # Mean longitude of the ascending node of Moon
    Om = ((450160.398036 + -6962890.5431 * t) % 1296000) / _asecperrad

    # compute nutation series using array loaded from data directory
    dat = _nut_data_00b
    arg = dat.nl * el + dat.nlp * elp + dat.nF * F + dat.nD * D + dat.nOm * Om
    sarg = np.sin(arg)
    carg = np.cos(arg)

    p1u_asecperrad = _asecperrad * 1e7  # 0.1 microasrcsecperrad
    dpsils = np.sum((dat.ps + dat.pst * t) * sarg + dat.pc * carg) / p1u_asecperrad
    depsls = np.sum((dat.ec + dat.ect * t) * carg + dat.es * sarg) / p1u_asecperrad
    # fixed offset in place of planetary tersm
    m_asecperrad = _asecperrad * 1e3  # milliarcsec per rad
    dpsipl = -0.135 / m_asecperrad
    depspl = 0.388 / m_asecperrad

    return epsa, dpsils + dpsipl, depsls + depspl  # all in radians


def nutation_matrix(epoch):
    """
    Nutation matrix generated from nutation components.

    Matrix converts from mean coordinate to true coordinate as
    r_true = M * r_mean
    """
    # TODO: implement higher precision 2006/2000A model if requested/needed
    epsa, dpsi, deps = nutation_components2000B(epoch.jd)  # all in radians

    return matrix_product(rotation_matrix(-(epsa + deps), 'x', False),
                          rotation_matrix(-dpsi, 'z', False),
                          rotation_matrix(epsa, 'x', False))
