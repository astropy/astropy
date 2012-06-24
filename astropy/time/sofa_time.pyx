import numpy as np
cimport numpy as np
import cython

ctypedef np.double_t DOUBLE_T

cdef extern from "../../cextern/sofa/sofa.h": 
    double iauEpb(double dj1, double dj2)
    double iauEpj(double dj1, double dj2)
    void iauEpj2jd(double epj, double *djm0, double *djm)
    void iauEpb2jd(double epb, double *djm0, double *djm)
    int iauCal2jd(int iy, int im, int id, double *djm0, double *djm)
    int iauJd2cal(double dj1, double dj2,int *iy, int *im, int *id, double *fd)
    int iauJdcalf(int ndp, double dj1, double dj2, int iymdf[4])

    # Internal JD to/from datetime format converters
    int iauD2dtf(char *scale, int ndp, double d1, double d2, int *iy, int *im, int *id, int ihmsf[4])
    int iauDtf2d(char *scale, int iy, int im, int id, int ihr, int imn, double sec, double *d1, double *d2)

    # Time scale helper routines
    double iauDtdb(double date1, double date2, double ut, double elong, double u, double v)
    int iauDat(int iy, int im, int id, double fd, double *deltat)

    # Time scale conversion routines
    int iauTaitt(double tai1, double tai2, double *tt1, double *tt2)
    int iauTttai(double tt1, double tt2, double *tai1, double *tai2)
    int iauTaiutc(double tai1, double tai2, double *utc1, double *utc2)
    int iauUtctai(double utc1, double utc2, double *tai1, double *tai2)
    int iauTcbtdb(double tcb1, double tcb2, double *tdb1, double *tdb2)
    int iauTdbtcb(double tdb1, double tdb2, double *tcb1, double *tcb2)
    int iauTcgtt(double tcg1, double tcg2, double *tt1, double *tt2)
    int iauTttcg(double tt1, double tt2, double *tcg1, double *tcg2)

    int iauTaiut1(double tai1, double tai2, double dta, double *ut11, double *ut12)
    int iauUt1tai(double ut11, double ut12, double dta, double *tai1, double *tai2)
    int iauTtut1(double tt1, double tt2, double dt, double *ut11, double *ut12)
    int iauUt1tt(double ut11, double ut12, double dt, double *tt1, double *tt2)
    int iauTdbtt(double tdb1, double tdb2, double dtr, double *tt1, double *tt2)
    int iauTttdb(double tt1, double tt2, double dtr, double *tdb1, double *tdb2)
    int iauUt1utc(double ut11, double ut12, double dut1, double *utc1, double *utc2)
    int iauUtcut1(double utc1, double utc2, double dut1, double *ut11, double *ut12)

    # Geodetic
    int iauAf2a(char s, int ideg, int iamin, double asec, double *rad)
    int iauGd2gc(int n, double elong, double phi, double height, double xyz[3])


@cython.wraparound(False)
@cython.boundscheck(False)
def cal2jd( 
    np.ndarray[int, ndim=1] iy,
    np.ndarray[int, ndim=1] im,
    np.ndarray[int, ndim=1] id,
    np.ndarray[double, ndim=1] djm0,
    np.ndarray[double, ndim=1] djm):
    cdef unsigned int i
    cdef unsigned n = iy.shape[0]
    for i in range(n):
        ret = iauCal2jd( iy[i], im[i], id[i], &djm0[i], &djm[i])
        if ret != 0:
            raise ValueError('Fail: {}'.format(ret))
    return

@cython.wraparound(False)
@cython.boundscheck(False)
def cal2jd( 
    np.ndarray[int, ndim=1] iy,
    np.ndarray[int, ndim=1] im,
    np.ndarray[int, ndim=1] id,
    np.ndarray[double, ndim=1] djm0,
    np.ndarray[double, ndim=1] djm):
    """
    """
    cdef unsigned int i
    cdef unsigned n = iy.shape[0]
    for i in range(n):
        ret = iauCal2jd( iy[i], im[i], id[i], &djm0[i], &djm[i])
        if ret != 0:
            raise ValueError('Fail: {}'.format(ret))
    return


@cython.wraparound(False)
@cython.boundscheck(False)
def d_tai_utc(np.ndarray[int, ndim=1] iy,
              np.ndarray[int, ndim=1] im,
              np.ndarray[int, ndim=1] id,
              np.ndarray[double, ndim=1] fd):
    """
    int iauDat(int iy, int im, int id, double fd, double *deltat)
    For a given UTC date, calculate delta(AT) = TAI-UTC.
    """
    cdef int i
    cdef int n = iy.shape[0]
    assert (iy.shape[0] == im.shape[0] == id.shape[0] == fd.shape[0])

    cdef np.ndarray[double, ndim=1] out = np.empty(n, dtype=np.double)

    for i in range(n):
        ret = iauDat(iy[i], im[i], id[i], fd[i],
                     &out[i])
        if ret != 0:
            raise ValueError('Fail: {}'.format(ret))
    return out


@cython.wraparound(False)
@cython.boundscheck(False)
def jd_dtf(scale, ndp,
              np.ndarray[double, ndim=1] d1,
              np.ndarray[double, ndim=1] d2):
    """
        int iauD2dtf(const char *scale, int ndp, double d1, double d2,
                 int *iy, int *im, int *id, int ihmsf[4])
    """
    cdef int i
    cdef int n = d1.shape[0]

    cdef np.ndarray[int, ndim=1] iy = np.empty(n, dtype=np.intc)
    cdef np.ndarray[int, ndim=1] im = np.empty(n, dtype=np.intc)
    cdef np.ndarray[int, ndim=1] id = np.empty(n, dtype=np.intc)
    cdef np.ndarray[int, ndim=2] ihmsf = np.empty((n, 4), dtype=np.intc)

    for i in range(n):
        ret = iauD2dtf(scale, ndp, d1[i], d2[i],
                     &iy[i], &im[i], &id[i], &ihmsf[i, 0])

        if ret != 0:
            raise ValueError('Fail: {}'.format(ret))
    return iy, im, id, ihmsf

@cython.wraparound(False)
@cython.boundscheck(False)
def dtf_jd(scale, 
              np.ndarray[int, ndim=1] iy,
              np.ndarray[int, ndim=1] im,
              np.ndarray[int, ndim=1] id,
              np.ndarray[int, ndim=1] ihr,
              np.ndarray[int, ndim=1] imn,
              np.ndarray[double, ndim=1] sec):
    """
    int iauDtf2d(char *scale, int iy, int im, int id, int ihr, int imn, double sec, double *d1, double *d2)
    """
    cdef int i
    cdef int n = iy.shape[0]
    assert (iy.shape[0] == im.shape[0] == id.shape[0] ==
            ihr.shape[0] == imn.shape[0] == sec.shape[0])

    cdef np.ndarray[double, ndim=1] out1 = np.empty(n, dtype=np.double)
    cdef np.ndarray[double, ndim=1] out2 = np.empty(n, dtype=np.double)

    for i in range(n):
        ret = iauDtf2d(scale, iy[i], im[i], id[i], ihr[i], imn[i], sec[i],
                       &out1[i], &out2[i])
        if ret != 0:
            raise ValueError('Fail: {}'.format(ret))
    return out1, out2


@cython.wraparound(False)
@cython.boundscheck(False)
def tai_tt( 
    np.ndarray[double, ndim=1] in1,
    np.ndarray[double, ndim=1] in2):
    """
    int iauTaitt(double tai1, double tai2, double *tt1, double *tt2)
    """
    assert in1.shape[0] == in2.shape[0]
    cdef unsigned n = in1.shape[0]
    cdef unsigned int i
    cdef np.ndarray[double, ndim=1] out1 = np.empty(n, dtype=np.double)
    cdef np.ndarray[double, ndim=1] out2 = np.empty(n, dtype=np.double)

    for i in range(n):
        ret = iauTaitt(in1[i], in2[i], &out1[i], &out2[i])
        if ret != 0:
            raise ValueError('Error code {} in iauTaitt'.format(ret))
    return out1, out2


@cython.wraparound(False)
@cython.boundscheck(False)
def tcb_tdb( 
    np.ndarray[double, ndim=1] in1,
    np.ndarray[double, ndim=1] in2):
    """
    int iauTcbtdb(double tcb1, double tcb2, double *tdb1, double *tdb2)
    """

    assert in1.shape[0] == in2.shape[0]
    cdef unsigned n = in1.shape[0]
    cdef unsigned int i
    cdef np.ndarray[double, ndim=1] out1 = np.empty(n, dtype=np.double)
    cdef np.ndarray[double, ndim=1] out2 = np.empty(n, dtype=np.double)

    for i in range(n):
        ret = iauTcbtdb(in1[i], in2[i], &out1[i], &out2[i])
        if ret != 0:
            raise ValueError('Error code {} in iauTcbtdb'.format(ret))
    return out1, out2


@cython.wraparound(False)
@cython.boundscheck(False)
def tcg_tt( 
    np.ndarray[double, ndim=1] in1,
    np.ndarray[double, ndim=1] in2):
    """
   int iauTcgtt(double tcg1, double tcg2, double *tt1, double *tt2)
    """

    assert in1.shape[0] == in2.shape[0]
    cdef unsigned n = in1.shape[0]
    cdef unsigned int i
    cdef np.ndarray[double, ndim=1] out1 = np.empty(n, dtype=np.double)
    cdef np.ndarray[double, ndim=1] out2 = np.empty(n, dtype=np.double)

    for i in range(n):
        ret = iauTcgtt(in1[i], in2[i], &out1[i], &out2[i])
        if ret != 0:
            raise ValueError('Error code {} in iauTcgtt'.format(ret))
    return out1, out2


@cython.wraparound(False)
@cython.boundscheck(False)
def tdb_tcb( 
    np.ndarray[double, ndim=1] in1,
    np.ndarray[double, ndim=1] in2):
    """
   int iauTdbtcb(double tdb1, double tdb2, double *tcb1, double *tcb2)
    """

    assert in1.shape[0] == in2.shape[0]
    cdef unsigned n = in1.shape[0]
    cdef unsigned int i
    cdef np.ndarray[double, ndim=1] out1 = np.empty(n, dtype=np.double)
    cdef np.ndarray[double, ndim=1] out2 = np.empty(n, dtype=np.double)

    for i in range(n):
        ret = iauTdbtcb(in1[i], in2[i], &out1[i], &out2[i])
        if ret != 0:
            raise ValueError('Error code {} in iauTdbtcb'.format(ret))
    return out1, out2


@cython.wraparound(False)
@cython.boundscheck(False)
def tt_tai( 
    np.ndarray[double, ndim=1] in1,
    np.ndarray[double, ndim=1] in2):
    """
   int iauTttai(double tt1, double tt2, double *tai1, double *tai2)
    """

    assert in1.shape[0] == in2.shape[0]
    cdef unsigned n = in1.shape[0]
    cdef unsigned int i
    cdef np.ndarray[double, ndim=1] out1 = np.empty(n, dtype=np.double)
    cdef np.ndarray[double, ndim=1] out2 = np.empty(n, dtype=np.double)

    for i in range(n):
        ret = iauTttai(in1[i], in2[i], &out1[i], &out2[i])
        if ret != 0:
            raise ValueError('Error code {} in iauTttai'.format(ret))
    return out1, out2


@cython.wraparound(False)
@cython.boundscheck(False)
def tt_tcg( 
    np.ndarray[double, ndim=1] in1,
    np.ndarray[double, ndim=1] in2):
    """
    int iauTttcg(double tt1, double tt2, double *tcg1, double *tcg2)
    """

    assert in1.shape[0] == in2.shape[0]
    cdef unsigned n = in1.shape[0]
    cdef unsigned int i
    cdef np.ndarray[double, ndim=1] out1 = np.empty(n, dtype=np.double)
    cdef np.ndarray[double, ndim=1] out2 = np.empty(n, dtype=np.double)

    for i in range(n):
        ret = iauTttcg(in1[i], in2[i], &out1[i], &out2[i])
        if ret != 0:
            raise ValueError('Error code {} in iauTttcg'.format(ret))
    return out1, out2


@cython.wraparound(False)
@cython.boundscheck(False)
def utc_tai( 
    np.ndarray[double, ndim=1] in1,
    np.ndarray[double, ndim=1] in2):
    """
    int iauUtctai(double utc1, double utc2, double *tai1, double *tai2)
    """
    assert in1.shape[0] == in2.shape[0]
    cdef unsigned n = in1.shape[0]
    cdef unsigned int i
    cdef np.ndarray[double, ndim=1] out1 = np.empty(n, dtype=np.double)
    cdef np.ndarray[double, ndim=1] out2 = np.empty(n, dtype=np.double)

    for i in range(n):
        ret = iauUtctai(in1[i], in2[i], &out1[i], &out2[i])
        if ret != 0:
            raise ValueError('Error code {} in iauUtctai'.format(ret))
    return out1, out2


@cython.wraparound(False)
@cython.boundscheck(False)
def tai_utc( 
    np.ndarray[double, ndim=1] in1,
    np.ndarray[double, ndim=1] in2):
    """
    int iauTaiutc(double tai1, double tai2, double *utc1, double *utc2)
    """
    assert in1.shape[0] == in2.shape[0]
    cdef unsigned n = in1.shape[0]
    cdef unsigned int i
    cdef np.ndarray[double, ndim=1] out1 = np.empty(n, dtype=np.double)
    cdef np.ndarray[double, ndim=1] out2 = np.empty(n, dtype=np.double)

    for i in range(n):
        ret = iauTaiutc(in1[i], in2[i], &out1[i], &out2[i])
        if ret != 0:
            raise ValueError('Error code {} in iauTaiutc'.format(ret))
    return out1, out2


@cython.wraparound(False)
@cython.boundscheck(False)
def tai_ut1( 
    np.ndarray[double, ndim=1] in1,
    np.ndarray[double, ndim=1] in2,
    np.ndarray[double, ndim=1] dt):
    """
    int iauTaiut1(double tai1, double tai2, double dta, double *ut11, double *ut12)
    """
    assert in1.shape[0] == in2.shape[0] == dt.shape[0]
    cdef unsigned n = in1.shape[0]
    cdef unsigned int i
    cdef np.ndarray[double, ndim=1] out1 = np.empty(n, dtype=np.double)
    cdef np.ndarray[double, ndim=1] out2 = np.empty(n, dtype=np.double)

    for i in range(n):
        ret = iauTaiut1(in1[i], in2[i], dt[i], &out1[i], &out2[i])
        if ret != 0:
            raise ValueError('Error code {} in iauTaiut1'.format(ret))
    return out1, out2


@cython.wraparound(False)
@cython.boundscheck(False)
def ut1_tai( 
    np.ndarray[double, ndim=1] in1,
    np.ndarray[double, ndim=1] in2,
    np.ndarray[double, ndim=1] dt):
    """
    int iauUt1tai(double ut11, double ut12, double dta, double *tai1, double *tai2)
    """
    assert in1.shape[0] == in2.shape[0] == dt.shape[0]
    cdef unsigned n = in1.shape[0]
    cdef unsigned int i
    cdef np.ndarray[double, ndim=1] out1 = np.empty(n, dtype=np.double)
    cdef np.ndarray[double, ndim=1] out2 = np.empty(n, dtype=np.double)

    for i in range(n):
        ret = iauUt1tai(in1[i], in2[i], dt[i], &out1[i], &out2[i])
        if ret != 0:
            raise ValueError('Error code {} in iauUt1tai'.format(ret))
    return out1, out2


@cython.wraparound(False)
@cython.boundscheck(False)
def tt_ut1( 
    np.ndarray[double, ndim=1] in1,
    np.ndarray[double, ndim=1] in2,
    np.ndarray[double, ndim=1] dt):
    """
    int iauTtut1(double tt1, double tt2, double dt, double *ut11, double *ut12)
    """
    assert in1.shape[0] == in2.shape[0] == dt.shape[0]
    cdef unsigned n = in1.shape[0]
    cdef unsigned int i
    cdef np.ndarray[double, ndim=1] out1 = np.empty(n, dtype=np.double)
    cdef np.ndarray[double, ndim=1] out2 = np.empty(n, dtype=np.double)

    for i in range(n):
        ret = iauTtut1(in1[i], in2[i], dt[i], &out1[i], &out2[i])
        if ret != 0:
            raise ValueError('Error code {} in iauTtut1'.format(ret))
    return out1, out2


@cython.wraparound(False)
@cython.boundscheck(False)
def ut1_tt( 
    np.ndarray[double, ndim=1] in1,
    np.ndarray[double, ndim=1] in2,
    np.ndarray[double, ndim=1] dt):
    """
    int iauUt1tt(double ut11, double ut12, double dt, double *tt1, double *tt2)
    """
    assert in1.shape[0] == in2.shape[0] == dt.shape[0]
    cdef unsigned n = in1.shape[0]
    cdef unsigned int i
    cdef np.ndarray[double, ndim=1] out1 = np.empty(n, dtype=np.double)
    cdef np.ndarray[double, ndim=1] out2 = np.empty(n, dtype=np.double)

    for i in range(n):
        ret = iauUt1tt(in1[i], in2[i], dt[i], &out1[i], &out2[i])
        if ret != 0:
            raise ValueError('Error code {} in iauUt1tt'.format(ret))
    return out1, out2


@cython.wraparound(False)
@cython.boundscheck(False)
def tdb_tt( 
    np.ndarray[double, ndim=1] in1,
    np.ndarray[double, ndim=1] in2,
    np.ndarray[double, ndim=1] dt):
    """
    int iauTdbtt(double tdb1, double tdb2, double dtr, double *tt1, double *tt2)
    """
    assert in1.shape[0] == in2.shape[0] == dt.shape[0]
    cdef unsigned n = in1.shape[0]
    cdef unsigned int i
    cdef np.ndarray[double, ndim=1] out1 = np.empty(n, dtype=np.double)
    cdef np.ndarray[double, ndim=1] out2 = np.empty(n, dtype=np.double)

    for i in range(n):
        ret = iauTdbtt(in1[i], in2[i], dt[i], &out1[i], &out2[i])
        if ret != 0:
            raise ValueError('Error code {} in iauTdbtt'.format(ret))
    return out1, out2


@cython.wraparound(False)
@cython.boundscheck(False)
def tt_tdb( 
    np.ndarray[double, ndim=1] in1,
    np.ndarray[double, ndim=1] in2,
    np.ndarray[double, ndim=1] dt):
    """
    int iauTttdb(double tt1, double tt2, double dtr, double *tdb1, double *tdb2)
    """
    assert in1.shape[0] == in2.shape[0] == dt.shape[0]
    cdef unsigned n = in1.shape[0]
    cdef unsigned int i
    cdef np.ndarray[double, ndim=1] out1 = np.empty(n, dtype=np.double)
    cdef np.ndarray[double, ndim=1] out2 = np.empty(n, dtype=np.double)

    for i in range(n):
        ret = iauTttdb(in1[i], in2[i], dt[i], &out1[i], &out2[i])
        if ret != 0:
            raise ValueError('Error code {} in iauTttdb'.format(ret))
    return out1, out2


@cython.wraparound(False)
@cython.boundscheck(False)
def ut1_utc( 
    np.ndarray[double, ndim=1] in1,
    np.ndarray[double, ndim=1] in2,
    np.ndarray[double, ndim=1] dt):
    """
    int iauUt1utc(double ut11, double ut12, double dut1, double *utc1, double *utc2)
    """
    assert in1.shape[0] == in2.shape[0] == dt.shape[0]
    cdef unsigned n = in1.shape[0]
    cdef unsigned int i
    cdef np.ndarray[double, ndim=1] out1 = np.empty(n, dtype=np.double)
    cdef np.ndarray[double, ndim=1] out2 = np.empty(n, dtype=np.double)

    for i in range(n):
        ret = iauUt1utc(in1[i], in2[i], dt[i], &out1[i], &out2[i])
        if ret != 0:
            raise ValueError('Error code {} in iauUt1utc'.format(ret))
    return out1, out2


@cython.wraparound(False)
@cython.boundscheck(False)
def utc_ut1( 
    np.ndarray[double, ndim=1] in1,
    np.ndarray[double, ndim=1] in2,
    np.ndarray[double, ndim=1] dt):
    """
    int iauUtcut1(double utc1, double utc2, double dut1, double *ut11, double *ut12)
    """
    assert in1.shape[0] == in2.shape[0] == dt.shape[0]
    cdef unsigned n = in1.shape[0]
    cdef unsigned int i
    cdef np.ndarray[double, ndim=1] out1 = np.empty(n, dtype=np.double)
    cdef np.ndarray[double, ndim=1] out2 = np.empty(n, dtype=np.double)

    for i in range(n):
        ret = iauUtcut1(in1[i], in2[i], dt[i], &out1[i], &out2[i])
        if ret != 0:
            raise ValueError('Error code {} in iauUtcut1'.format(ret))
    return out1, out2


@cython.wraparound(False)
@cython.boundscheck(False)
def d_tdb_tt(np.ndarray[double, ndim=1] in1,
             np.ndarray[double, ndim=1] in2,
             np.ndarray[double, ndim=1] ut,
             elong, u, v):
    """
    compute DTR = TDB-TT
    double iauDtdb(double date1, double date2, double ut,
                   double elong, double u, double v)
    """
    assert in1.shape[0] == in2.shape[0] == ut.shape[0]
    cdef unsigned n = in1.shape[0]
    cdef unsigned int i
    cdef np.ndarray[double, ndim=1] out = np.empty(n, dtype=np.double)
    cdef double c_elong = elong
    cdef double c_u = u
    cdef double c_v = v

    for i in range(n):
        out[i] = iauDtdb(in1[i], in2[i], ut[i], c_elong, c_u, c_v)
    return out


def iau_af2a(sign, ideg, iamin, asec):
    """
    Wrap
    int iauAf2a(char s, int ideg, int iamin, double asec, double *rad)
    """
    cdef double rad
    s = ord(sign)
    ret = iauAf2a(s, ideg, iamin, asec, &rad)
    if ret != 0:
        raise ValueError('Error code {}'.format(ret))

    return rad

def iau_gd2gc(n, elong, phi, height):
    """
    Wrap
    int iauGd2gc(int n, double elong, double phi, double height, double xyz[3])
    """
    cdef np.ndarray[double, ndim=1] xyz = np.empty(3, dtype=np.double)
    ret = iauGd2gc(n, elong, phi, height, &xyz[0])
    if ret != 0:
        raise ValueError('Error code {}'.format(ret))

    return xyz


@cython.wraparound(False)
@cython.boundscheck(False)
def jd_julian_epoch(np.ndarray[double, ndim=1] jd1,
                    np.ndarray[double, ndim=1] jd2):
    """ Wrap double iauEpj(double dj1, double dj2)
    **  Julian Date to Julian Epoch.
    """
    assert jd1.shape[0] == jd2.shape[0]
    cdef unsigned n = jd1.shape[0]
    cdef unsigned int i
    cdef np.ndarray[double, ndim=1] epd = np.empty(n, dtype=np.double)

    for i in range(n):
        epd[i] = iauEpj(jd1[i], jd2[i])
    return epd


@cython.wraparound(False)
@cython.boundscheck(False)
def julian_epoch_jd(np.ndarray[double, ndim=1] epd):
    """ Wrap void iauEpj2jd(double epj, double *djm0, double *djm)
    **  Julian Epoch to Julian Date.
    """
    cdef unsigned n = epd.shape[0]
    cdef unsigned int i
    cdef np.ndarray[double, ndim=1] jd1 = np.empty(n, dtype=np.double)
    cdef np.ndarray[double, ndim=1] jd2 = np.empty(n, dtype=np.double)

    for i in range(n):
        iauEpj2jd(epd[i], &jd1[i], &jd2[i])
    return jd1, jd2


@cython.wraparound(False)
@cython.boundscheck(False)
def jd_besselian_epoch(np.ndarray[double, ndim=1] jd1,
                       np.ndarray[double, ndim=1] jd2):
    """ Wrap double iauEpb(double dj1, double dj2)
    **  Julian Date to Besselian Epoch.
    """
    assert jd1.shape[0] == jd2.shape[0]
    cdef unsigned n = jd1.shape[0]
    cdef unsigned int i
    cdef np.ndarray[double, ndim=1] epd = np.empty(n, dtype=np.double)

    for i in range(n):
        epd[i] = iauEpb(jd1[i], jd2[i])
    return epd


@cython.wraparound(False)
@cython.boundscheck(False)
def besselian_epoch_jd(np.ndarray[double, ndim=1] epd):
    """ Wrap void iauEpb2jd(double epj, double *djm0, double *djm)
    **  Besselian Epoch to Julian Date.
    """
    cdef unsigned n = epd.shape[0]
    cdef unsigned int i
    cdef np.ndarray[double, ndim=1] jd1 = np.empty(n, dtype=np.double)
    cdef np.ndarray[double, ndim=1] jd2 = np.empty(n, dtype=np.double)

    for i in range(n):
        iauEpb2jd(epd[i], &jd1[i], &jd2[i])
    return jd1, jd2
