# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module contains standard functions for earth orientation, such as
precession and nutation.

This module is (currently) not intended to be part of the public API, but
is instead primarily for internal use in `coordinates`
"""

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
    from ..obstools import jd2000

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
def precession_matrix_J2000_Capitaine(epoch):
        """
        Computes the precession matrix from J2000 to the given Julian Epoch.
        Expression from from Capitaine et al. 2003 as expressed in the USNO
        Circular 179.  This should match the IAU 2006 standard from SOFA.

        Parameters
        ----------
            epoch : int
                The epoch at which to compute the precession matrix.

        Returns
        -------
        pmatrix : 3x3 array
            Precession matrix at `epoch`
        """
        from ..utils import rotation_matrix

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


def load_nutation_data(datafn, seriestype):
    """
    Loads nutation series from saved data files.

    Seriestype can be 'lunisolar' or 'planetary'
    """
    from ..utils.data import get_package_data

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

    lines = [l for l in get_package_data(datafn).split('\n') if not l.startswith('#') if not l.strip()=='']

    lists = [[] for n in dtypes]
    for l in lines:
        for i,e in enumerate(l.split(' ')):
            lists[i].append(dtypes[i][1](e))
    return np.rec.fromarrays(lists,names=[e[0] for e in dtypes])

#TODO: replace w/SOFA equivalent
_nut_data_00a_ls = load_nutation_data('iau00a_nutation_ls.tab','lunisolar')
_nut_data_00a_pl = load_nutation_data('iau00a_nutation_pl.tab','planetary')
def nutation_components20062000A(epoch):
    """
    :returns: eps,dpsi,deps in radians
    """
    from ..obstools import epoch_to_jd
    from .funcs import obliquity

    epsa = obliquity(epoch_to_jd(epoch),2006)

    raise NotImplementedError('2006/2000A nutation model not implemented')

    return epsa,dpsi,deps



_nut_data_00b = load_nutation_data('iau00b_nutation.tab','lunisolar')
def nutation_components2000B(intime,asepoch=True):
    """
    :param intime: time to compute the nutation components as a JD or epoch
    :type intime: scalar
    :param asepoch: if True, `intime` is interpreted as an epoch, otherwise JD
    :type asepoch: bool

    :returns: eps,dpsi,deps in radians
    """
    from ..constants import asecperrad
    from ..obstools import epoch_to_jd,jd2000
    from .funcs import obliquity

    if asepoch:
        jd = epoch_to_jd(intime)
    else:
        jd = intime
    epsa = np.radians(obliquity(jd,2000))
    t = (jd-jd2000)/36525

    #Fundamental (Delaunay) arguments from Simon et al. (1994) via SOFA
    #Mean anomaly of moon
    el = ((485868.249036 + 1717915923.2178*t)%1296000)/asecperrad
    #Mean anomaly of sun
    elp = ((1287104.79305 + 129596581.0481*t)%1296000)/asecperrad
    #Mean argument of the latitude of Moon
    F = ((335779.526232 + 1739527262.8478*t)%1296000)/asecperrad
    #Mean elongation of the Moon from Sun
    D = ((1072260.70369 + 1602961601.2090*t)%1296000)/asecperrad
    #Mean longitude of the ascending node of Moon
    Om = ((450160.398036 + -6962890.5431*t)%1296000)/asecperrad

    #compute nutation series using array loaded from data directory
    dat = _nut_data_00b
    arg = dat.nl*el + dat.nlp*elp + dat.nF*F + dat.nD*D + dat.nOm*Om
    sarg = np.sin(arg)
    carg = np.cos(arg)

    p1uasecperrad = asecperrad*1e7 #0.1 microasrcsecperrad
    dpsils = np.sum((dat.ps + dat.pst*t)*sarg + dat.pc*carg)/p1uasecperrad
    depsls = np.sum((dat.ec + dat.ect*t)*carg + dat.es*sarg)/p1uasecperrad
    #fixed offset in place of planetary tersm
    masecperrad = asecperrad*1e3 #milliarcsec per rad
    dpsipl = -0.135/masecperrad
    depspl =  0.388/masecperrad

    return epsa,dpsils+dpsipl,depsls+depspl #all in radians

def nutation_matrix(epoch):
    """
    Nutation matrix generated from nutation components.

    Matrix converts from mean coordinate to true coordinate as
    r_true = M * r_mean
    """
    from ..utils import rotation_matrix

    #TODO: implement higher precision 2006/2000A model if requested/needed
    epsa,dpsi,deps = nutation_components2000B(epoch) #all in radians

    return rotation_matrix(-(epsa + deps),'x',False) *\
           rotation_matrix(-dpsi,'z',False) *\
           rotation_matrix(epsa,'x',False)
