# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module contains convenience functions implementing some of the
algorithms contained within Jean Meeus, 'Astronomical Algorithms',
second edition, 1998, Willmann-Bell.
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import numpy as np
from numpy.polynomial.polynomial import polyval

from .. import units as u, coordinates as coord
from .. import _erfa as erfa
from . import GCRS, SkyCoord, GeocentricTrueEcliptic, ICRS
from .builtin_frames.utils import get_jd12, cartrepr_from_matmul


# Meeus 1998: table 47.A
#   D   M   M'  F   l    r
_MOON_L_R = (
     (0, 0,  1,  0,  6288774,    -20905355),
     (2, 0,  -1, 0,  1274027,    -3699111),
     (2, 0,  0,  0,  658314,     -2955968),
     (0, 0,  2,  0,  213618,     -569925),
     (0, 1,  0,  0,  -185116,    48888),
     (0, 0,  0,  2,  -114332,    -3149),
     (2, 0,  -2, 0,  58793,      246158),
     (2, -1, -1, 0,  57066,      -152138),
     (2, 0,  1,  0,  53322,      -170733),
     (2, -1, 0,  0,  45758,      -204586),
     (0, 1,  -1, 0,  -40923,     -129620),
     (1, 0,  0,  0,  -34720,     108743),
     (0, 1,  1,  0,  -30383,     104755),
     (2, 0,  0,  -2, 15327,      10321),
     (0, 0,  1,  2,  -12528,     0),
     (0, 0,  1,  -2, 10980,      79661),
     (4, 0,  -1, 0,  10675,      -34782),
     (0, 0,  3,  0,  10034,      -23210),
     (4, 0,  -2, 0,  8548,       -21636),
     (2, 1,  -1, 0,  -7888,      24208),
     (2, 1,  0,  0,  -6766,      30824),
     (1, 0,  -1, 0,  -5163,      -8379),
     (1, 1,  0,  0,  4987,       -16675),
     (2, -1, 1,  0,  4036,       -12831),
     (2, 0,  2,  0,  3994,       -10445),
     (4, 0,  0,  0,  3861,       -11650),
     (2, 0,  -3, 0,  3665,       14403),
     (0, 1,  -2, 0,  -2689,      -7003),
     (2, 0,  -1, 2,  -2602,      0),
     (2, -1, -2, 0,  2390,       10056),
     (1, 0,  1,  0,  -2348,      6322),
     (2, -2, 0,  0,  2236,       -9884),
     (0, 1,  2,  0,  -2120,      5751),
     (0, 2,  0,  0,  -2069,      0),
     (2, -2, -1, 0,  2048,       -4950),
     (2, 0,  1,  -2, -1773,      4130),
     (2, 0,  0,  2,  -1595,      0),
     (4, -1, -1, 0,  1215,       -3958),
     (0, 0,  2,  2,  -1110,      0),
     (3, 0,  -1, 0,  -892,       3258),
     (2, 1,  1,  0,  -810,       2616),
     (4, -1, -2, 0,  759,        -1897),
     (0, 2,  -1, 0,  -713,       -2117),
     (2, 2,  -1, 0,  -700,       2354),
     (2, 1,  -2, 0,  691,        0),
     (2, -1, 0,  -2, 596,        0),
     (4, 0,  1,  0,  549,        -1423),
     (0, 0,  4,  0,  537,        -1117),
     (4, -1, 0,  0,  520,        -1571),
     (1, 0,  -2, 0,  -487,       -1739),
     (2, 1,  0,  -2, -399,       0),
     (0, 0,  2,  -2, -381,       -4421),
     (1, 1,  1,  0,  351,        0),
     (3, 0,  -2, 0,  -340,       0),
     (4, 0,  -3, 0,  330,        0),
     (2, -1, 2,  0,  327,        0),
     (0, 2,  1,  0,  -323,       1165),
     (1, 1,  -1, 0,  299,        0),
     (2, 0,  3,  0,  294,        0),
     (2, 0,  -1, -2, 0,          8752)
)

# Meeus 1998: table 47.B
#   D   M   M'  F   b
_MOON_B = (
     (0, 0,  0,  1,  5128122),
     (0, 0,  1,  1,  280602),
     (0, 0,  1,  -1, 277693),
     (2, 0,  0,  -1, 173237),
     (2, 0,  -1, 1,  55413),
     (2, 0,  -1, -1, 46271),
     (2, 0,  0,  1,  32573),
     (0, 0,  2,  1,  17198),
     (2, 0,  1,  -1, 9266),
     (0, 0,  2,  -1, 8822),
     (2, -1, 0,  -1, 8216),
     (2, 0,  -2, -1, 4324),
     (2, 0,  1,  1,  4200),
     (2, 1,  0,  -1, -3359),
     (2, -1, -1, 1,  2463),
     (2, -1, 0,  1,  2211),
     (2, -1, -1, -1, 2065),
     (0, 1,  -1, -1, -1870),
     (4, 0,  -1, -1, 1828),
     (0, 1,  0,  1,  -1794),
     (0, 0,  0,  3,  -1749),
     (0, 1,  -1, 1,  -1565),
     (1, 0,  0,  1,  -1491),
     (0, 1,  1,  1,  -1475),
     (0, 1,  1,  -1, -1410),
     (0, 1,  0,  -1, -1344),
     (1, 0,  0,  -1, -1335),
     (0, 0,  3,  1,  1107),
     (4, 0,  0,  -1, 1021),
     (4, 0,  -1, 1,  833),
#second column
     (0, 0,  1,  -3, 777),
     (4, 0,  -2, 1,  671),
     (2, 0,  0,  -3, 607),
     (2, 0,  2,  -1, 596),
     (2, -1, 1,  -1, 491),
     (2, 0,  -2, 1,  -451),
     (0, 0,  3,  -1, 439),
     (2, 0,  2,  1,  422),
     (2, 0,  -3, -1, 421),
     (2, 1,  -1, 1,  -366),
     (2, 1,  0,  1,  -351),
     (4, 0,  0,  1,  331),
     (2, -1, 1,  1,  315),
     (2, -2, 0,  -1, 302),
     (0, 0,  1,  3,  -283),
     (2, 1,  1,  -1, -229),
     (1, 1,  0,  -1, 223),
     (1, 1,  0,  1,  223),
     (0, 1,  -2, -1, -220),
     (2, 1,  -1, -1, -220),
     (1, 0,  1,  1,  -185),
     (2, -1, -2, -1, 181),
     (0, 1,  2,  1,  -177),
     (4, 0,  -2, -1, 176),
     (4, -1, -1, -1, 166),
     (1, 0,  1,  -1, -164),
     (4, 0,  1,  -1, 132),
     (1, 0,  -1, -1, -119),
     (4, -1, 0,  -1, 115),
     (2, -2, 0,  1,  107)
)

"""Coefficients of polynomials for various terms"""
_coLc = (2.18316448e+02,  4.81267881e+05, -1.57860000e-03, 
         1.85583502e-06, -1.53388349e-08)
_coD = (2.97850192e+02,  4.45267111e+05, -1.88190000e-03,
        1.83194472e-06, -8.84447000e-09)
_coM = (3.57529109e+02,  3.59990503e+04, -1.53600000e-04,
        4.08329931e-08)
_coMc = (1.34963396e+02,  4.77198868e+05,  8.74140000e-03,
         1.43474081e-05, -6.79717238e-08)
_coF = (9.32720950e+01,  4.83202018e+05, -3.65390000e-03,
        -2.83607487e-07,  1.15833246e-09)
_coA1 = (119.75, 131.849)
_coA2 = (53.09, 479264.290)
_coA3 = (313.45 , 481266.484)
_coE = (1.0, -0.002516, -0.0000074)

def calcMoon(t, location=None):
    """
    Lunar position model ELP2000-82 of (Chapront-Touze' and Chapront, 1983, 124, 50)

    This is the simplified version of Jean Meeus, Astronomical Algorithms,
    second edition, 1998, Willmann-Bell. Meeus claims approximate accuracy of 10" 
    in longitude and 4" in latitude, with no specified time range.
    
    Tests against JPL ephemerides show accuracy of 10 arcseconds and 50 km over the
    date range CE 1950-2050.
    
    Parameters
    -----------
    time : `~astropy.time.Time`
        Time of observation.
    location : `~astropy.coordinates.EarthLocation`, optional
        Location of observer on the Earth.  If not given, will be taken from
        ``time`` (if not present, a geocentric observer will be assumed).
    
    Returns
    --------
    skycoord : `~astropy.coordinates.SkyCoord`
        GCRS Coordinate for the body    
    """    
    if location is None:
        location = t.location
            
    # number of centuries since J2000.0.
    # This should strictly speaking be in Ephemeris Time, but TDB or TT
    # will introduce error smaller than intrinsic accuracy of algorithm.
    T = (t.tdb.jyear-2000.0)/100.
    
    # constants that are needed for all calculations
    Lc = u.deg*polyval(T, _coLc)
    D = u.deg*polyval(T, _coD)
    M = u.deg*polyval(T, _coM)
    Mc = u.deg*polyval(T, _coMc)
    F = u.deg*polyval(T, _coF)

    # convert to angle objects and wrap at 360 degrees
    Lc = coord.Angle(Lc).wrap_at(360*u.deg)
    D = coord.Angle(D).wrap_at(360*u.deg)
    M = coord.Angle(M).wrap_at(360*u.deg)
    Mc = coord.Angle(Mc).wrap_at(360*u.deg)
    F = coord.Angle(F).wrap_at(360*u.deg)

    A1 = coord.Angle(u.deg*polyval(T, _coA1)).wrap_at(360*u.deg)
    A2 = coord.Angle(u.deg*polyval(T, _coA2)).wrap_at(360*u.deg)
    A3 = coord.Angle(u.deg*polyval(T, _coA3)).wrap_at(360*u.deg)
    E = polyval(T, _coE)

    suml = sumr = 0.0
    lrnum = len(_MOON_L_R)
    for i in range(lrnum):
        DNum = _MOON_L_R[i][0]
        MNum = _MOON_L_R[i][1]
        McNum = _MOON_L_R[i][2]
        FNum = _MOON_L_R[i][3]

        corr = 1
        if MNum == 1 or MNum == -1:
            corr = E
        if MNum == 2 or MNum == -2:
            corr = E*E
 
        suml += _MOON_L_R[i][4]*corr*np.sin(D*DNum+M*MNum+Mc*McNum+F*FNum)
        sumr += _MOON_L_R[i][5]*corr*np.cos(D*DNum+M*MNum+Mc*McNum+F*FNum)

    sumb = 0.0
    bnum = len(_MOON_B)
    for i in range(bnum):
        DNum = _MOON_B[i][0]
        MNum = _MOON_B[i][1]
        McNum = _MOON_B[i][2]
        FNum = _MOON_B[i][3]

        corr = 1
        if MNum == 1 or MNum == -1:
            corr = E
        if MNum == 2 or MNum == -2:
            corr = E*E
 
        sumb += _MOON_B[i][4]*corr*np.sin(D*DNum+M*MNum+Mc*McNum+F*FNum)

    suml += (3958*np.sin(A1)+1962*np.sin(Lc-F)+318*np.sin(A2))
    sumb += (-2235*np.sin(Lc)+382*np.sin(A3)+175*np.sin(A1-F)+
             175*np.sin(A1+F)+127*np.sin(Lc-Mc)-115*np.sin(Lc+Mc))

    # ensure units
    suml = suml*u.deg
    sumb = sumb*u.deg
    
    # nutation of longitude
    jd1, jd2 = get_jd12(t, 'tt')
    nut, _ = erfa.nut06a(jd1, jd2)
    nut = nut*u.rad
    
    # calculate ecliptic coordinates
    lon = Lc + suml/1000000 + nut
    lat = sumb/1000000
    dist = (385000.56+sumr/1000)*u.km
    
    # Meeus algorithm gives GeocentricTrueEcliptic coordinates
    ecliptic_coo = GeocentricTrueEcliptic(lon, lat, distance=dist,
                                          equinox=t)
    if location is None:
        return SkyCoord(ecliptic_coo.transform_to(GCRS(obstime=t)))
        
    loc, vel = location.get_gcrs_posvel(t)
    """
    I suspect a bug in the coordinate transformations between GeocentricTrueEcliptic 
    and GCRS.  This version, with an intermediate transformation to ICRS, works.
    """
    icrs = ecliptic_coo.transform_to(ICRS())
    return SkyCoord(icrs.transform_to(GCRS(obstime=t, 
                                           obsgeoloc=loc,
                                           obsgeovel=vel)))

    