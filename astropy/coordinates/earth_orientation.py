# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module contains standard functions for earth orientation, such as
precession and nutation.

This module is (currently) not intended to be part of the public API, but
is instead primarily for internal use in `coordinates`
"""
import numpy as np

from ..time import Time
from .. import units as u

jd2000 = Time('J2000', scale='utc').jd
_asecperrad = u.radian.to(u.arcsec)


def obliquity(jd, algorithm=2006):
    """
    Computes the obliquity of the Earth at the requested Julian Date.

    Parameters
    ----------
    jd : scalar or array-like
        julian date at which to compute obliquity
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
    T = (jd-jd2000)/36525.0

    if algorithm==2006:
        p = (-0.0000000434,-0.000000576,0.00200340,-0.0001831,-46.836769,84381.406)
        corr = 0
    elif algorithm==2000:
        p = (0.001813,-0.00059,-46.8150,84381.448)
        corr = -0.02524*T
    elif algorithm==1980:
        p = (0.001813,-0.00059,-46.8150,84381.448)
        corr = 0
    else:
        raise ValueError('invalid algorithm year for computing obliquity')

    return (np.polyval(p,T)+corr)/3600.


#TODO: replace this with SOFA equivalent
def precession_matrix_Capitaine(fromepoch, toepoch):
        """
        Parameters
        ----------
        fromepoch : `~astropy.time.Time`
            The epoch to precess from.
        toepoch : `~astropy.time.Time`
            The epoch to precess to.

        Returns
        -------
        pmatrix : 3x3 array
            Precession matrix to get from `fromepoch` to `toepoch`

        References
        ----------
        USNO Circular 179
        """
        mat_fromto2000 = _precession_matrix_J2000_Capitaine(fromepoch.jyear)
        mat_2000toto = _precession_matrix_J2000_Capitaine(fromepoch.jyear).T

        return np.dot(mat_fromto2000, mat_2000toto)


def _precession_matrix_J2000_Capitaine(epoch):
    """
    Computes the precession matrix from J2000 to the given Julian Epoch.
    Expression from from Capitaine et al. 2003 as expressed in the USNO
    Circular 179.  This should match the IAU 2006 standard from SOFA.
    """
    from .angles import rotation_matrix

    T = (epoch-2000.0)/100.0
    #from USNO circular
    pzeta = (-0.0000003173,-0.000005971,0.01801828,0.2988499,2306.083227,2.650545)
    pz = (-0.0000002904,-0.000028596,0.01826837,1.0927348,2306.077181,-2.650545)
    ptheta = (-0.0000001274,-0.000007089,-0.04182264,-0.4294934,2004.191903,0)
    zeta = np.polyval(pzeta,T)/3600.0
    z = np.polyval(pz,T)/3600.0
    theta = np.polyval(ptheta,T)/3600.0

    return rotation_matrix(-z,'z') *\
           rotation_matrix(theta,'y') *\
           rotation_matrix(-zeta,'z')




def _load_nutation_data(datafn, seriestype):
    """
    Loads nutation series from saved data files.

    Seriestype can be 'lunisolar' or 'planetary'
    """
    from os.path import join

    from ..config.data import get_data_contents

    if seriestype == 'lunisolar':
        dtypes = [('nl',int),
                  ('nlp',int),
                  ('nF',int),
                  ('nD',int),
                  ('nOm',int),
                  ('ps',float),
                  ('pst',float),
                  ('pc',float),
                  ('ec',float),
                  ('ect',float),
                  ('es',float)]
    elif seriestype == 'planetary':
        dtypes = [('nl',int),
                  ('nF',int),
                  ('nD',int),
                  ('nOm',int),
                  ('nme',int),
                  ('nve',int),
                  ('nea',int),
                  ('nma',int),
                  ('nju',int),
                  ('nsa',int),
                  ('nur',int),
                  ('nne',int),
                  ('npa',int),
                  ('sp',int),
                  ('cp',int),
                  ('se',int),
                  ('ce',int)]
    else:
        raise ValueError('requested invalid nutation series type')

    lines = [l for l in get_data_contents(join('data', datafn)).split('\n') if not l.startswith('#') if not l.strip()=='']

    lists = [[] for n in dtypes]
    for l in lines:
        for i,e in enumerate(l.split(' ')):
            lists[i].append(dtypes[i][1](e))
    return np.rec.fromarrays(lists,names=[e[0] for e in dtypes])


_nut_data_00b = _load_nutation_data('iau00b_nutation.tab','lunisolar')
#TODO: replace w/SOFA equivalent
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

    #Fundamental (Delaunay) arguments from Simon et al. (1994) via SOFA
    #Mean anomaly of moon
    el = ((485868.249036 + 1717915923.2178*t)%1296000)/_asecperrad
    #Mean anomaly of sun
    elp = ((1287104.79305 + 129596581.0481*t)%1296000)/_asecperrad
    #Mean argument of the latitude of Moon
    F = ((335779.526232 + 1739527262.8478*t)%1296000)/_asecperrad
    #Mean elongation of the Moon from Sun
    D = ((1072260.70369 + 1602961601.2090*t)%1296000)/_asecperrad
    #Mean longitude of the ascending node of Moon
    Om = ((450160.398036 + -6962890.5431*t)%1296000)/_asecperrad

    #compute nutation series using array loaded from data directory
    dat = _nut_data_00b
    arg = dat.nl*el + dat.nlp*elp + dat.nF*F + dat.nD*D + dat.nOm*Om
    sarg = np.sin(arg)
    carg = np.cos(arg)

    p1u_asecperrad = _asecperrad*1e7 #0.1 microasrcsecperrad
    dpsils = np.sum((dat.ps + dat.pst*t)*sarg + dat.pc*carg)/p1u_asecperrad
    depsls = np.sum((dat.ec + dat.ect*t)*carg + dat.es*sarg)/p1u_asecperrad
    #fixed offset in place of planetary tersm
    m_asecperrad = _asecperrad*1e3 #milliarcsec per rad
    dpsipl = -0.135/m_asecperrad
    depspl =  0.388/m_asecperrad

    return epsa,dpsils+dpsipl,depsls+depspl #all in radians

def nutation_matrix(epoch):
    """
    Nutation matrix generated from nutation components.

    Matrix converts from mean coordinate to true coordinate as
    r_true = M * r_mean
    """
    from .angles import rotation_matrix

    #TODO: implement higher precision 2006/2000A model if requested/needed
    epsa,dpsi,deps = nutation_components2000B(epoch.jd) #all in radians

    return rotation_matrix(-(epsa + deps),'x',False) *\
           rotation_matrix(-dpsi,'z',False) *\
           rotation_matrix(epsa,'x',False)
