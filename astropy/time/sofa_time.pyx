import warnings

import numpy as np
cimport numpy as np
import cython

ctypedef np.double_t DOUBLE_T

cdef extern from "sofa.h":
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

DUBIOUS = 'dubious year for UTC (before 1960.0 or 5 years ' \
          'beyond last known leap second)'

def check_return(ret, func_name, warns={}, errors={}):
    """Check the return value from an iau routine"""
    if ret in warns:
        warnings.warn('{0}: {1}'.format(func_name, warns[ret]))
    elif ret in errors:
        raise ValueError('{0}: {1}'.format(func_name, errors[ret]))
    elif ret != 0:
        raise ValueError('Unexpected return code {0} from {1}'
                         .format(repr(ret), func_name))


@cython.wraparound(False)
@cython.boundscheck(False)
def cal2jd( 
    np.ndarray[int, ndim=1] iy,
    np.ndarray[int, ndim=1] im,
    np.ndarray[int, ndim=1] id,
    np.ndarray[double, ndim=1] djm0,
    np.ndarray[double, ndim=1] djm):
    """
    int iauCal2jd(int iy, int im, int id, double *djm0, double *djm)
    Calendar date to high-precision JD.

    **  Given:
    **     iy,im,id  int     year, month, day in Gregorian calendar (Note 1)
    **
    **  Returned:
    **     djm0      double  MJD zero-point: always 2400000.5
    **     djm       double  Modified Julian Date for 0 hrs
    **
    **  Returned (function value):
    **               int     status:
    **                           0 = OK
    **                          -1 = bad year   (Note 3: JD not computed)
    **                          -2 = bad month  (JD not computed)
    **                          -3 = bad day    (JD computed)
    **
    **  Notes:
    **
    **  1) The algorithm used is valid from -4800 March 1, but this
    **     implementation rejects dates before -4799 January 1.
    **
    **  2) The Julian Date is returned in two pieces, in the usual SOFA
    **     manner, which is designed to preserve time resolution.  The
    **     Julian Date is available as a single number by adding djm0 and
    **     djm.
    **
    **  3) In early eras the conversion is from the "Proleptic Gregorian
    **     Calendar";  no account is taken of the date(s) of adoption of
    **     the Gregorian Calendar, nor is the AD/BC numbering convention
    **     observed.
    """
    cdef unsigned int i
    cdef unsigned n = iy.shape[0]
    warns = {-3: 'Bad input day (JD still computed)'}
    errs = {-1: 'Bad input year',
             -2: 'Bad input month'}
    for i in range(n):
        ret = iauCal2jd( iy[i], im[i], id[i], &djm0[i], &djm[i])
        check_return(ret, 'iauCal2jd', warns, errs)
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

    **  Given:
    **     iy     int      UTC:  year (Notes 1 and 2)
    **     im     int            month (Note 2)
    **     id     int            day (Notes 2 and 3)
    **     fd     double         fraction of day (Note 4)
    **
    **  Returned:
    **     deltat double   TAI minus UTC, seconds
    **
    **  Returned (function value):
    **            int      status (Note 5):
    **                       1 = dubious year (Note 1)
    **                       0 = OK
    **                      -1 = bad year
    **                      -2 = bad month
    **                      -3 = bad day (Note 3)
    **                      -4 = bad fraction (Note 4)
    **
    **  Notes:
    **
    **  1) UTC began at 1960 January 1.0 (JD 2436934.5) and it is improper
    **     to call the function with an earlier date.  If this is attempted,
    **     zero is returned together with a warning status.
    **
    **     Because leap seconds cannot, in principle, be predicted in
    **     advance, a reliable check for dates beyond the valid range is
    **     impossible.  To guard against gross errors, a year five or more
    **     after the release year of the present function (see parameter
    **     IYV) is considered dubious.  In this case a warning status is
    **     returned but the result is computed in the normal way.
    **
    **     For both too-early and too-late years, the warning status is
    **     j=+1.  This is distinct from the error status j=-1, which
    **     signifies a year so early that JD could not be computed.
    **
    **  2) If the specified date is for a day which ends with a leap second,
    **     the UTC-TAI value returned is for the period leading up to the
    **     leap second.  If the date is for a day which begins as a leap
    **     second ends, the UTC-TAI returned is for the period following the
    **     leap second.
    **
    **  3) The day number must be in the normal calendar range, for example
    **     1 through 30 for April.  The "almanac" convention of allowing
    **     such dates as January 0 and December 32 is not supported in this
    **     function, in order to avoid confusion near leap seconds.
    **
    **  4) The fraction of day is used only for dates before the
    **     introduction of leap seconds, the first of which occurred at the
    **     end of 1971.  It is tested for validity (zero to less than 1 is
    **     the valid range) even if not used;  if invalid, zero is used and
    **     status j=-4 is returned.  For many applications, setting fd to
    **     zero is acceptable;  the resulting error is always less than 3 ms
    **     (and occurs only pre-1972).
    **
    **  5) The status value returned in the case where there are multiple
    **     errors refers to the first error detected.  For example, if the
    **     month and day are 13 and 32 respectively, j=-2 (bad month)
    **     will be returned.
    **
    **  6) In cases where a valid result is not available, zero is returned.

    **                       1 = dubious year (Note 1)
    **                       0 = OK
    **                      -1 = bad year
    **                      -2 = bad month
    **                      -3 = bad day (Note 3)
    **                      -4 = bad fraction (Note 4)

    """
    cdef int i
    cdef int n = iy.shape[0]
    assert (iy.shape[0] == im.shape[0] == id.shape[0] == fd.shape[0])

    cdef np.ndarray[double, ndim=1] out = np.empty(n, dtype=np.double)

    warns = {1: DUBIOUS}
    errs = {-1: 'bad year',
             -2: 'bad month (must be 1 to 12)',
             -3: 'bad day (must be within normal calendar date for a month)',
             -4: 'bad fraction of day'}
    for i in range(n):
        ret = iauDat(iy[i], im[i], id[i], fd[i],
                     &out[i])
        check_return(ret, 'iauDat', warns, errs)

    return out


@cython.wraparound(False)
@cython.boundscheck(False)
def jd_dtf(scale, ndp,
              np.ndarray[double, ndim=1] d1,
              np.ndarray[double, ndim=1] d2):
    """
    int iauD2dtf(const char *scale, int ndp, double d1, double d2,
             int *iy, int *im, int *id, int ihmsf[4])

    **  Given:
    **     scale     char[]  time scale ID (Note 1)
    **     ndp       int     resolution (Note 2)
    **     d1,d2     double  time as a 2-part Julian Date (Notes 3,4)
    **
    **  Returned:
    **     iy,im,id  int     year, month, day in Gregorian calendar (Note 5)
    **     ihmsf     int[4]  hours, minutes, seconds, fraction (Note 1)
    **
    **  Returned (function value):
    **               int     status: +1 = dubious year (Note 5)
    **                                0 = OK
    **                               -1 = unacceptable date (Note 6)
    **
    **  Notes:
    **
    **  1) scale identifies the time scale.  Only the value "UTC" (in upper
    **     case) is significant, and enables handling of leap seconds (see
    **     Note 4).
    **
    **  2) ndp is the number of decimal places in the seconds field, and can
    **     have negative as well as positive values, such as:
    **
    **     ndp         resolution
    **     -4            1 00 00
    **     -3            0 10 00
    **     -2            0 01 00
    **     -1            0 00 10
    **      0            0 00 01
    **      1            0 00 00.1
    **      2            0 00 00.01
    **      3            0 00 00.001
    **
    **     The limits are platform dependent, but a safe range is -5 to +9.
    **
    **  3) d1+d2 is Julian Date, apportioned in any convenient way between
    **     the two arguments, for example where d1 is the Julian Day Number
    **     and d2 is the fraction of a day.  In the case of UTC, where the
    **     use of JD is problematical, special conventions apply:  see the
    **     next note.
    **
    **  4) JD cannot unambiguously represent UTC during a leap second unless
    **     special measures are taken.  The SOFA internal convention is that
    **     the quasi-JD day represents UTC days whether the length is 86399,
    **     86400 or 86401 SI seconds.
    **
    **  5) The warning status "dubious year" flags UTCs that predate the
    **     introduction of the time scale and that are too far in the future
    **     to be trusted.  See iauDat for further details.
    **
    **  6) For calendar conventions and limitations, see iauCal2jd.
    """
    cdef int i
    cdef int n = d1.shape[0]

    cdef np.ndarray[int, ndim=1] iy = np.empty(n, dtype=np.intc)
    cdef np.ndarray[int, ndim=1] im = np.empty(n, dtype=np.intc)
    cdef np.ndarray[int, ndim=1] id = np.empty(n, dtype=np.intc)
    cdef np.ndarray[int, ndim=2] ihmsf = np.empty((n, 4), dtype=np.intc)

    warns = {1: DUBIOUS}
    errs = {-1: 'unacceptable date'}

    for i in range(n):
        ret = iauD2dtf(scale, ndp, d1[i], d2[i],
                     &iy[i], &im[i], &id[i], &ihmsf[i, 0])
        check_return(ret, 'iauD2dtf', warns, errs)

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

    **  Given:
    **     scale     char[]  time scale ID (Note 1)
    **     iy,im,id  int     year, month, day in Gregorian calendar (Note 2)
    **     ihr,imn   int     hour, minute
    **     sec       double  seconds
    **
    **  Returned:
    **     d1,d2     double  2-part Julian Date (Notes 3,4)
    **
    **  Returned (function value):
    **               int     status: +3 = both of next two
    **                               +2 = time is after end of day (Note 5)
    **                               +1 = dubious year (Note 6)
    **                                0 = OK
    **                               -1 = bad year
    **                               -2 = bad month
    **                               -3 = bad day
    **                               -4 = bad hour
    **                               -5 = bad minute
    **                               -6 = bad second (<0)
    **
    **  Notes:
    **
    **  1) scale identifies the time scale.  Only the value "UTC" (in upper
    **     case) is significant, and enables handling of leap seconds (see
    **     Note 4).
    **
    **  2) For calendar conventions and limitations, see iauCal2jd.
    **
    **  3) The sum of the results, d1+d2, is Julian Date, where normally d1
    **     is the Julian Day Number and d2 is the fraction of a day.  In the
    **     case of UTC, where the use of JD is problematical, special
    **     conventions apply:  see the next note.
    **
    **  4) JD cannot unambiguously represent UTC during a leap second unless
    **     special measures are taken.  The SOFA internal convention is that
    **     the quasi-JD day represents UTC days whether the length is 86399,
    **     86400 or 86401 SI seconds.
    **
    **  5) The warning status "time is after end of day" usually means that
    **     the sec argument is greater than 60.0.  However, in a day ending
    **     in a leap second the limit changes to 61.0 (or 59.0 in the case
    **     of a negative leap second).
    **
    **  6) The warning status "dubious year" flags UTCs that predate the
    **     introduction of the time scale and that are too far in the future
    **     to be trusted.  See iauDat for further details.
    **
    **  7) Only in the case of continuous and regular time scales (TAI, TT,
    **     TCG, TCB and TDB) is the result d1+d2 a Julian Date, strictly
    **     speaking.  In the other cases (UT1 and UTC) the result must be
    **     used with circumspection;  in particular the difference between
    **     two such results cannot be interpreted as a precise time
    **     interval.

    """
    cdef int i
    cdef int n = iy.shape[0]
    assert (iy.shape[0] == im.shape[0] == id.shape[0] ==
            ihr.shape[0] == imn.shape[0] == sec.shape[0])

    cdef np.ndarray[double, ndim=1] out1 = np.empty(n, dtype=np.double)
    cdef np.ndarray[double, ndim=1] out2 = np.empty(n, dtype=np.double)

    warns = {3: 'time is after end of day and ' + DUBIOUS,
             2: 'time is after end of day',
             1: DUBIOUS}
    errs = {-1: 'bad year',
            -2: 'bad month',
            -3: 'bad day',
            -4: 'bad hour',
            -5: 'bad minute',
            -6: 'bad second (< 0)'}

    for i in range(n):
        ret = iauDtf2d(scale, iy[i], im[i], id[i], ihr[i], imn[i], sec[i],
                       &out1[i], &out2[i])
        check_return(ret, 'iauDtf2d', warns, errs)

    return out1, out2


@cython.wraparound(False)
@cython.boundscheck(False)
def tai_tt( 
    np.ndarray[double, ndim=1] in1,
    np.ndarray[double, ndim=1] in2):
    """
    int iauTaitt(double tai1, double tai2, double *tt1, double *tt2)

    **  Given:
    **     tai1,tai2  double    TAI as a 2-part Julian Date
    **
    **  Returned:
    **     tt1,tt2    double    TT as a 2-part Julian Date
    **
    **  Returned (function value):
    **                int       status:  0 = OK
    **
    **  Note:
    **
    **     tai1+tai2 is Julian Date, apportioned in any convenient way
    **     between the two arguments, for example where tai1 is the Julian
    **     Day Number and tai2 is the fraction of a day.  The returned
    **     tt1,tt2 follow suit.
    """
    assert in1.shape[0] == in2.shape[0]
    cdef unsigned n = in1.shape[0]
    cdef unsigned int i
    cdef np.ndarray[double, ndim=1] out1 = np.empty(n, dtype=np.double)
    cdef np.ndarray[double, ndim=1] out2 = np.empty(n, dtype=np.double)

    for i in range(n):
        ret = iauTaitt(in1[i], in2[i], &out1[i], &out2[i])
        check_return(ret, 'iauTaitt')

    return out1, out2


@cython.wraparound(False)
@cython.boundscheck(False)
def tcb_tdb( 
    np.ndarray[double, ndim=1] in1,
    np.ndarray[double, ndim=1] in2):
    """
    int iauTcbtdb(double tcb1, double tcb2, double *tdb1, double *tdb2)

    **  Given:
    **     tcb1,tcb2  double    TCB as a 2-part Julian Date
    **
    **  Returned:
    **     tdb1,tdb2  double    TDB as a 2-part Julian Date
    **
    **  Returned (function value):
    **                int       status:  0 = OK
    **
    **  Notes:
    **
    **  1) tcb1+tcb2 is Julian Date, apportioned in any convenient way
    **     between the two arguments, for example where tcb1 is the Julian
    **     Day Number and tcb2 is the fraction of a day.  The returned
    **     tdb1,tdb2 follow suit.
    **
    **  2) The 2006 IAU General Assembly introduced a conventional linear
    **     transformation between TDB and TCB.  This transformation
    **     compensates for the drift between TCB and terrestrial time TT,
    **     and keeps TDB approximately centered on TT.  Because the
    **     relationship between TT and TCB depends on the adopted solar
    **     system ephemeris, the degree of alignment between TDB and TT over
    **     long intervals will vary according to which ephemeris is used.
    **     Former definitions of TDB attempted to avoid this problem by
    **     stipulating that TDB and TT should differ only by periodic
    **     effects.  This is a good description of the nature of the
    **     relationship but eluded precise mathematical formulation.  The
    **     conventional linear relationship adopted in 2006 sidestepped
    **     these difficulties whilst delivering a TDB that in practice was
    **     consistent with values before that date.
    **
    **  3) TDB is essentially the same as Teph, the time argument for the
    **     JPL solar system ephemerides.
    """

    assert in1.shape[0] == in2.shape[0]
    cdef unsigned n = in1.shape[0]
    cdef unsigned int i
    cdef np.ndarray[double, ndim=1] out1 = np.empty(n, dtype=np.double)
    cdef np.ndarray[double, ndim=1] out2 = np.empty(n, dtype=np.double)

    for i in range(n):
        ret = iauTcbtdb(in1[i], in2[i], &out1[i], &out2[i])
        check_return(ret, 'iauTcbtdb')

    return out1, out2


@cython.wraparound(False)
@cython.boundscheck(False)
def tcg_tt( 
    np.ndarray[double, ndim=1] in1,
    np.ndarray[double, ndim=1] in2):
    """
   int iauTcgtt(double tcg1, double tcg2, double *tt1, double *tt2)

   **  Given:
   **     tcg1,tcg2  double    TCG as a 2-part Julian Date
   **
   **  Returned:
   **     tt1,tt2    double    TT as a 2-part Julian Date
   **
   **  Returned (function value):
   **                int       status:  0 = OK
    """

    assert in1.shape[0] == in2.shape[0]
    cdef unsigned n = in1.shape[0]
    cdef unsigned int i
    cdef np.ndarray[double, ndim=1] out1 = np.empty(n, dtype=np.double)
    cdef np.ndarray[double, ndim=1] out2 = np.empty(n, dtype=np.double)

    for i in range(n):
        ret = iauTcgtt(in1[i], in2[i], &out1[i], &out2[i])
        check_return(ret, 'iauTcgtt')

    return out1, out2


@cython.wraparound(False)
@cython.boundscheck(False)
def tdb_tcb( 
    np.ndarray[double, ndim=1] in1,
    np.ndarray[double, ndim=1] in2):
    """
    int iauTdbtcb(double tdb1, double tdb2, double *tcb1, double *tcb2)

    **  Given:
    **     tdb1,tdb2  double    TDB as a 2-part Julian Date
    **
    **  Returned:
    **     tcb1,tcb2  double    TCB as a 2-part Julian Date
    **
    **  Returned (function value):
    **                int       status:  0 = OK
    """
    assert in1.shape[0] == in2.shape[0]
    cdef unsigned n = in1.shape[0]
    cdef unsigned int i
    cdef np.ndarray[double, ndim=1] out1 = np.empty(n, dtype=np.double)
    cdef np.ndarray[double, ndim=1] out2 = np.empty(n, dtype=np.double)

    for i in range(n):
        ret = iauTdbtcb(in1[i], in2[i], &out1[i], &out2[i])
        check_return(ret, 'iauTdbtcb')

    return out1, out2


@cython.wraparound(False)
@cython.boundscheck(False)
def tt_tai( 
    np.ndarray[double, ndim=1] in1,
    np.ndarray[double, ndim=1] in2):
    """
    int iauTttai(double tt1, double tt2, double *tai1, double *tai2)

    **  Given:
    **     tt1,tt2    double    TT as a 2-part Julian Date
    **
    **  Returned:
    **     tai1,tai2  double    TAI as a 2-part Julian Date
    **
    **  Returned (function value):
    **                int       status:  0 = OK
    """
    assert in1.shape[0] == in2.shape[0]
    cdef unsigned n = in1.shape[0]
    cdef unsigned int i
    cdef np.ndarray[double, ndim=1] out1 = np.empty(n, dtype=np.double)
    cdef np.ndarray[double, ndim=1] out2 = np.empty(n, dtype=np.double)

    for i in range(n):
        ret = iauTttai(in1[i], in2[i], &out1[i], &out2[i])
        check_return(ret, 'iauTttai')

    return out1, out2


@cython.wraparound(False)
@cython.boundscheck(False)
def tt_tcg( 
    np.ndarray[double, ndim=1] in1,
    np.ndarray[double, ndim=1] in2):
    """
    int iauTttcg(double tt1, double tt2, double *tcg1, double *tcg2)
    **  Given:
    **     tt1,tt2    double    TT as a 2-part Julian Date
    **
    **  Returned:
    **     tcg1,tcg2  double    TCG as a 2-part Julian Date
    **
    **  Returned (function value):
    **                int       status:  0 = OK
    """
    assert in1.shape[0] == in2.shape[0]
    cdef unsigned n = in1.shape[0]
    cdef unsigned int i
    cdef np.ndarray[double, ndim=1] out1 = np.empty(n, dtype=np.double)
    cdef np.ndarray[double, ndim=1] out2 = np.empty(n, dtype=np.double)

    for i in range(n):
        ret = iauTttcg(in1[i], in2[i], &out1[i], &out2[i])
        check_return(ret, 'iauTttcg')

    return out1, out2


@cython.wraparound(False)
@cython.boundscheck(False)
def utc_tai( 
    np.ndarray[double, ndim=1] in1,
    np.ndarray[double, ndim=1] in2):
    """
    int iauUtctai(double utc1, double utc2, double *tai1, double *tai2)

    **  Given:
    **     utc1,utc2  double   UTC as a 2-part quasi Julian Date (Notes 1-4)
    **
    **  Returned:
    **     tai1,tai2  double   TAI as a 2-part Julian Date (Note 5)
    **
    **  Returned (function value):
    **                int      status: +1 = dubious year (Note 3)
    **                                  0 = OK
    **                                 -1 = unacceptable date
    **
    **  Notes:
    **
    **  1) utc1+utc2 is quasi Julian Date (see Note 2), apportioned in any
    **     convenient way between the two arguments, for example where utc1
    **     is the Julian Day Number and utc2 is the fraction of a day.
    **
    **  2) JD cannot unambiguously represent UTC during a leap second unless
    **     special measures are taken.  The convention in the present
    **     function is that the JD day represents UTC days whether the
    **     length is 86399, 86400 or 86401 SI seconds.
    **
    **  3) The warning status "dubious year" flags UTCs that predate the
    **     introduction of the time scale and that are too far in the future
    **     to be trusted.  See iauDat  for further details.
    **
    **  4) The function iauDtf2d converts from calendar date and time of day
    **     into 2-part Julian Date, and in the case of UTC implements the
    **     leap-second-ambiguity convention described above.
    **
    **  5) The returned TAI1,TAI2 are such that their sum is the TAI Julian
    **     Date.
    **
    """
    assert in1.shape[0] == in2.shape[0]
    cdef unsigned n = in1.shape[0]
    cdef unsigned int i
    cdef np.ndarray[double, ndim=1] out1 = np.empty(n, dtype=np.double)
    cdef np.ndarray[double, ndim=1] out2 = np.empty(n, dtype=np.double)

    warns = {1: DUBIOUS}
    errs = {-1: 'unacceptable date'}

    for i in range(n):
        ret = iauUtctai(in1[i], in2[i], &out1[i], &out2[i])
        check_return(ret, 'iauUtctai', warns, errs)

    return out1, out2


@cython.wraparound(False)
@cython.boundscheck(False)
def tai_utc( 
    np.ndarray[double, ndim=1] in1,
    np.ndarray[double, ndim=1] in2):
    """
    int iauTaiutc(double tai1, double tai2, double *utc1, double *utc2)

    **  Given:
    **     tai1,tai2  double   TAI as a 2-part Julian Date (Note 1)
    **
    **  Returned:
    **     utc1,utc2  double   UTC as a 2-part quasi Julian Date (Notes 1-3)
    **
    **  Returned (function value):
    **                int      status: +1 = dubious year (Note 4)
    **                                  0 = OK
    **                                 -1 = unacceptable date
    **
    **  Notes:
    **
    **  1) tai1+tai2 is Julian Date, apportioned in any convenient way
    **     between the two arguments, for example where tai1 is the Julian
    **     Day Number and tai2 is the fraction of a day.  The returned utc1
    **     and utc2 form an analogous pair, except that a special convention
    **     is used, to deal with the problem of leap seconds - see the next
    **     note.
    **
    **  2) JD cannot unambiguously represent UTC during a leap second unless
    **     special measures are taken.  The convention in the present
    **     function is that the JD day represents UTC days whether the
    **     length is 86399, 86400 or 86401 SI seconds.
    **
    **  3) The function iauD2dtf can be used to transform the UTC quasi-JD
    **     into calendar date and clock time, including UTC leap second
    **     handling.
    **
    **  4) The warning status "dubious year" flags UTCs that predate the
    **     introduction of the time scale and that are too far in the future
    **     to be trusted.  See iauDat for further details.
        """
    assert in1.shape[0] == in2.shape[0]
    cdef unsigned n = in1.shape[0]
    cdef unsigned int i
    cdef np.ndarray[double, ndim=1] out1 = np.empty(n, dtype=np.double)
    cdef np.ndarray[double, ndim=1] out2 = np.empty(n, dtype=np.double)

    warns = {1: DUBIOUS}
    errs = {-1: 'unacceptable date'}

    for i in range(n):
        ret = iauTaiutc(in1[i], in2[i], &out1[i], &out2[i])
        check_return(ret, 'iauTaiutc', warns, errs)

    return out1, out2


@cython.wraparound(False)
@cython.boundscheck(False)
def tai_ut1( 
    np.ndarray[double, ndim=1] in1,
    np.ndarray[double, ndim=1] in2,
    np.ndarray[double, ndim=1] dt):
    """
    int iauTaiut1(double tai1, double tai2, double dta, double *ut11, double *ut12)

    **  Given:
    **     tai1,tai2  double    TAI as a 2-part Julian Date
    **     dta        double    UT1-TAI in seconds
    **
    **  Returned:
    **     ut11,ut12  double    UT1 as a 2-part Julian Date
    **
    **  Returned (function value):
    **                int       status:  0 = OK
    **
    **  Notes:
    **
    **  1) tai1+tai2 is Julian Date, apportioned in any convenient way
    **     between the two arguments, for example where tai1 is the Julian
    **     Day Number and tai2 is the fraction of a day.  The returned
    **     UT11,UT12 follow suit.
    **
    **  2) The argument dta, i.e. UT1-TAI, is an observed quantity, and is
    **     available from IERS tabulations.
    """
    assert in1.shape[0] == in2.shape[0] == dt.shape[0]
    cdef unsigned n = in1.shape[0]
    cdef unsigned int i
    cdef np.ndarray[double, ndim=1] out1 = np.empty(n, dtype=np.double)
    cdef np.ndarray[double, ndim=1] out2 = np.empty(n, dtype=np.double)

    for i in range(n):
        ret = iauTaiut1(in1[i], in2[i], dt[i], &out1[i], &out2[i])
        check_return(ret, 'iauTaiut1')

    return out1, out2


@cython.wraparound(False)
@cython.boundscheck(False)
def ut1_tai( 
    np.ndarray[double, ndim=1] in1,
    np.ndarray[double, ndim=1] in2,
    np.ndarray[double, ndim=1] dt):
    """
    int iauUt1tai(double ut11, double ut12, double dta, double *tai1, double *tai2)

    **  Given:
    **     ut11,ut12  double    UT1 as a 2-part Julian Date
    **     dta        double    UT1-TAI in seconds
    **
    **  Returned:
    **     tai1,tai2  double    TAI as a 2-part Julian Date
    **
    **  Returned (function value):
    **                int       status:  0 = OK
    **
    **  Notes:
    **
    **  1) ut11+ut12 is Julian Date, apportioned in any convenient way
    **     between the two arguments, for example where ut11 is the Julian
    **     Day Number and ut12 is the fraction of a day.  The returned
    **     tai1,tai2 follow suit.
    **
    **  2) The argument dta, i.e. UT1-TAI, is an observed quantity, and is
    **     available from IERS tabulations.
    """
    assert in1.shape[0] == in2.shape[0] == dt.shape[0]
    cdef unsigned n = in1.shape[0]
    cdef unsigned int i
    cdef np.ndarray[double, ndim=1] out1 = np.empty(n, dtype=np.double)
    cdef np.ndarray[double, ndim=1] out2 = np.empty(n, dtype=np.double)

    for i in range(n):
        ret = iauUt1tai(in1[i], in2[i], dt[i], &out1[i], &out2[i])
        check_return(ret, 'iauUt1tai')

    return out1, out2


@cython.wraparound(False)
@cython.boundscheck(False)
def tt_ut1( 
    np.ndarray[double, ndim=1] in1,
    np.ndarray[double, ndim=1] in2,
    np.ndarray[double, ndim=1] dt):
    """
    int iauTtut1(double tt1, double tt2, double dt, double *ut11, double *ut12)

    **  Given:
    **     tt1,tt2    double    TT as a 2-part Julian Date
    **     dt         double    TT-UT1 in seconds
    **
    **  Returned:
    **     ut11,ut12  double    UT1 as a 2-part Julian Date
    **
    **  Returned (function value):
    **                int       status:  0 = OK
    **
    **  Notes:
    **
    **  1) tt1+tt2 is Julian Date, apportioned in any convenient way between
    **     the two arguments, for example where tt1 is the Julian Day Number
    **     and tt2 is the fraction of a day.  The returned ut11,ut12 follow
    **     suit.
    **
    **  2) The argument dt is classical Delta T.
    """
    assert in1.shape[0] == in2.shape[0] == dt.shape[0]
    cdef unsigned n = in1.shape[0]
    cdef unsigned int i
    cdef np.ndarray[double, ndim=1] out1 = np.empty(n, dtype=np.double)
    cdef np.ndarray[double, ndim=1] out2 = np.empty(n, dtype=np.double)

    for i in range(n):
        ret = iauTtut1(in1[i], in2[i], dt[i], &out1[i], &out2[i])
        check_return(ret, 'iauTtut1')

    return out1, out2


@cython.wraparound(False)
@cython.boundscheck(False)
def ut1_tt( 
    np.ndarray[double, ndim=1] in1,
    np.ndarray[double, ndim=1] in2,
    np.ndarray[double, ndim=1] dt):
    """
    int iauUt1tt(double ut11, double ut12, double dt, double *tt1, double *tt2)

    **  Given:
    **     ut11,ut12  double    UT1 as a 2-part Julian Date
    **     dt         double    TT-UT1 in seconds
    **
    **  Returned:
    **     tt1,tt2    double    TT as a 2-part Julian Date
    **
    **  Returned (function value):
    **                int       status:  0 = OK
    **
    **  Notes:
    **
    **  1) ut11+ut12 is Julian Date, apportioned in any convenient way
    **     between the two arguments, for example where ut11 is the Julian
    **     Day Number and ut12 is the fraction of a day.  The returned
    **     tt1,tt2 follow suit.
    **
    **  2) The argument dt is classical Delta T.
    """
    assert in1.shape[0] == in2.shape[0] == dt.shape[0]
    cdef unsigned n = in1.shape[0]
    cdef unsigned int i
    cdef np.ndarray[double, ndim=1] out1 = np.empty(n, dtype=np.double)
    cdef np.ndarray[double, ndim=1] out2 = np.empty(n, dtype=np.double)

    for i in range(n):
        ret = iauUt1tt(in1[i], in2[i], dt[i], &out1[i], &out2[i])
        check_return(ret, 'iauUt1tt')

    return out1, out2


@cython.wraparound(False)
@cython.boundscheck(False)
def tdb_tt( 
    np.ndarray[double, ndim=1] in1,
    np.ndarray[double, ndim=1] in2,
    np.ndarray[double, ndim=1] dt):
    """
    int iauTdbtt(double tdb1, double tdb2, double dtr, double *tt1, double *tt2)

    **  Given:
    **     tdb1,tdb2  double    TDB as a 2-part Julian Date
    **     dtr        double    TDB-TT in seconds
    **
    **  Returned:
    **     tt1,tt2    double    TT as a 2-part Julian Date
    **
    **  Returned (function value):
    **                int       status:  0 = OK
    **
    **  Notes:
    **
    **  1) tdb1+tdb2 is Julian Date, apportioned in any convenient way
    **     between the two arguments, for example where tdb1 is the Julian
    **     Day Number and tdb2 is the fraction of a day.  The returned
    **     tt1,tt2 follow suit.
    **
    **  2) The argument dtr represents the quasi-periodic component of the
    **     GR transformation between TT and TCB.  It is dependent upon the
    **     adopted solar-system ephemeris, and can be obtained by numerical
    **     integration, by interrogating a precomputed time ephemeris or by
    **     evaluating a model such as that implemented in the SOFA function
    **     iauDtdb.   The quantity is dominated by an annual term of 1.7 ms
    **     amplitude.
    **
    **  3) TDB is essentially the same as Teph, the time argument for the
    **     JPL solar system ephemerides.
    """
    assert in1.shape[0] == in2.shape[0] == dt.shape[0]
    cdef unsigned n = in1.shape[0]
    cdef unsigned int i
    cdef np.ndarray[double, ndim=1] out1 = np.empty(n, dtype=np.double)
    cdef np.ndarray[double, ndim=1] out2 = np.empty(n, dtype=np.double)

    for i in range(n):
        ret = iauTdbtt(in1[i], in2[i], dt[i], &out1[i], &out2[i])
        check_return(ret, 'iauTdbtt')

    return out1, out2


@cython.wraparound(False)
@cython.boundscheck(False)
def tt_tdb( 
    np.ndarray[double, ndim=1] in1,
    np.ndarray[double, ndim=1] in2,
    np.ndarray[double, ndim=1] dt):
    """
    int iauTttdb(double tt1, double tt2, double dtr, double *tdb1, double *tdb2)

    **  Given:
    **     tt1,tt2    double    TT as a 2-part Julian Date
    **     dtr        double    TDB-TT in seconds
    **
    **  Returned:
    **     tdb1,tdb2  double    TDB as a 2-part Julian Date
    **
    **  Returned (function value):
    **                int       status:  0 = OK
    **
    **  Notes:
    **
    **  1) tt1+tt2 is Julian Date, apportioned in any convenient way between
    **     the two arguments, for example where tt1 is the Julian Day Number
    **     and tt2 is the fraction of a day.  The returned tdb1,tdb2 follow
    **     suit.
    **
    **  2) The argument dtr represents the quasi-periodic component of the
    **     GR transformation between TT and TCB.  It is dependent upon the
    **     adopted solar-system ephemeris, and can be obtained by numerical
    **     integration, by interrogating a precomputed time ephemeris or by
    **     evaluating a model such as that implemented in the SOFA function
    **     iauDtdb.   The quantity is dominated by an annual term of 1.7 ms
    **     amplitude.
    **
    **  3) TDB is essentially the same as Teph, the time argument for the JPL
    **     solar system ephemerides.
    """
    assert in1.shape[0] == in2.shape[0] == dt.shape[0]
    cdef unsigned n = in1.shape[0]
    cdef unsigned int i
    cdef np.ndarray[double, ndim=1] out1 = np.empty(n, dtype=np.double)
    cdef np.ndarray[double, ndim=1] out2 = np.empty(n, dtype=np.double)

    for i in range(n):
        ret = iauTttdb(in1[i], in2[i], dt[i], &out1[i], &out2[i])
        check_return(ret, 'iauTttdb')

    return out1, out2


@cython.wraparound(False)
@cython.boundscheck(False)
def ut1_utc( 
    np.ndarray[double, ndim=1] in1,
    np.ndarray[double, ndim=1] in2,
    np.ndarray[double, ndim=1] dt):
    """
    int iauUt1utc(double ut11, double ut12, double dut1, double *utc1, double *utc2)

    **  Given:
    **     ut11,ut12  double   UT1 as a 2-part Julian Date (Note 1)
    **     dut1       double   Delta UT1: UT1-UTC in seconds (Note 2)
    **
    **  Returned:
    **     utc1,utc2  double   UTC as a 2-part quasi Julian Date (Notes 3,4)
    **
    **  Returned (function value):
    **                int      status: +1 = dubious year (Note 5)
    **                                  0 = OK
    **                                 -1 = unacceptable date
    **
    **  Notes:
    **
    **  1) ut11+ut12 is Julian Date, apportioned in any convenient way
    **     between the two arguments, for example where ut11 is the Julian
    **     Day Number and ut12 is the fraction of a day.  The returned utc1
    **     and utc2 form an analogous pair, except that a special convention
    **     is used, to deal with the problem of leap seconds - see Note 3.
    **
    **  2) Delta UT1 can be obtained from tabulations provided by the
    **     International Earth Rotation and Reference Systems Service.  The
    **     value changes abruptly by 1s at a leap second;  however, close to
    **     a leap second the algorithm used here is tolerant of the "wrong"
    **     choice of value being made.
    **
    **  3) JD cannot unambiguously represent UTC during a leap second unless
    **     special measures are taken.  The convention in the present
    **     function is that the returned quasi JD day UTC1+UTC2 represents
    **     UTC days whether the length is 86399, 86400 or 86401 SI seconds.
    **
    **  4) The function iauD2dtf can be used to transform the UTC quasi-JD
    **     into calendar date and clock time, including UTC leap second
    **     handling.
    **
    **  5) The warning status "dubious year" flags UTCs that predate the
    **     introduction of the time scale and that are too far in the future
    **     to be trusted.  See iauDat for further details.
    """
    assert in1.shape[0] == in2.shape[0] == dt.shape[0]
    cdef unsigned n = in1.shape[0]
    cdef unsigned int i
    cdef np.ndarray[double, ndim=1] out1 = np.empty(n, dtype=np.double)
    cdef np.ndarray[double, ndim=1] out2 = np.empty(n, dtype=np.double)

    warns = {1: DUBIOUS}
    errs = {-1: 'unacceptable date'}

    for i in range(n):
        ret = iauUt1utc(in1[i], in2[i], dt[i], &out1[i], &out2[i])
        check_return(ret, 'iauUt1utc', warns, errs)

    return out1, out2


@cython.wraparound(False)
@cython.boundscheck(False)
def utc_ut1( 
    np.ndarray[double, ndim=1] in1,
    np.ndarray[double, ndim=1] in2,
    np.ndarray[double, ndim=1] dt):
    """
    int iauUtcut1(double utc1, double utc2, double dut1, double *ut11, double *ut12)

    **  Given:
    **     utc1,utc2  double   UTC as a 2-part quasi Julian Date (Notes 1-4)
    **     dut1       double   Delta UT1 = UT1-UTC in seconds (Note 5)
    **
    **  Returned:
    **     ut11,ut12  double   UT1 as a 2-part Julian Date (Note 6)
    **
    **  Returned (function value):
    **                int      status: +1 = dubious year (Note 7)
    **                                  0 = OK
    **                                 -1 = unacceptable date
    **
    **  Notes:
    **
    **  1) utc1+utc2 is quasi Julian Date (see Note 2), apportioned in any
    **     convenient way between the two arguments, for example where utc1
    **     is the Julian Day Number and utc2 is the fraction of a day.
    **
    **  2) JD cannot unambiguously represent UTC during a leap second unless
    **     special measures are taken.  The convention in the present
    **     function is that the JD day represents UTC days whether the
    **     length is 86399, 86400 or 86401 SI seconds.
    **
    **  3) The warning status "dubious year" flags UTCs that predate the
    **     introduction of the time scale and that are too far in the future
    **     to be trusted.  See iauDat  for further details.
    **
    **  4) The function iauDtf2d  converts from calendar date and time of
    **     day into 2-part Julian Date, and in the case of UTC implements
    **     the leap-second-ambiguity convention described above.
    **
    **  5) Delta UT1 can be obtained from tabulations provided by the
    **     International Earth Rotation and Reference Systems Service.  It
    **     It is the caller's responsibility to supply a DUT argument
    **     containing the UT1-UTC value that matches the given UTC.
    **
    **  6) The returned ut11,ut12 are such that their sum is the UT1 Julian
    **     Date.
    **
    **  7) The warning status "dubious year" flags UTCs that predate the
    **     introduction of the time scale and that are too far in the future
    **     to be trusted.  See iauDat for further details.
    """
    assert in1.shape[0] == in2.shape[0] == dt.shape[0]
    cdef unsigned n = in1.shape[0]
    cdef unsigned int i
    cdef np.ndarray[double, ndim=1] out1 = np.empty(n, dtype=np.double)
    cdef np.ndarray[double, ndim=1] out2 = np.empty(n, dtype=np.double)

    warns = {1: DUBIOUS}
    errs = {-1: 'unacceptable date'}

    for i in range(n):
        ret = iauUtcut1(in1[i], in2[i], dt[i], &out1[i], &out2[i])
        check_return(ret, 'iauUtcut1', warns, errs)

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

    **  Given:
    **     date1,date2   double  date, TDB (Notes 1-3)
    **     ut            double  universal time (UT1, fraction of one day)
    **     elong         double  longitude (east positive, radians)
    **     u             double  distance from Earth spin axis (km)
    **     v             double  distance north of equatorial plane (km)
    **
    **  Returned (function value):
    **                   double  TDB-TT (seconds)
    **
    **  Notes:
    **
    **  1) The date date1+date2 is a Julian Date, apportioned in any
    **     convenient way between the two arguments.  For example,
    **     JD(TT)=2450123.7 could be expressed in any of these ways,
    **     among others:
    **
    **            date1          date2
    **
    **         2450123.7           0.0       (JD method)
    **         2451545.0       -1421.3       (J2000 method)
    **         2400000.5       50123.2       (MJD method)
    **         2450123.5           0.2       (date & time method)
    **
    **     The JD method is the most natural and convenient to use in
    **     cases where the loss of several decimal digits of resolution
    **     is acceptable.  The J2000 method is best matched to the way
    **     the argument is handled internally and will deliver the
    **     optimum resolution.  The MJD method and the date & time methods
    **     are both good compromises between resolution and convenience.
    **
    **     Although the date is, formally, barycentric dynamical time (TDB),
    **     the terrestrial dynamical time (TT) can be used with no practical
    **     effect on the accuracy of the prediction.
    **
    **  2) TT can be regarded as a coordinate time that is realized as an
    **     offset of 32.184s from International Atomic Time, TAI.  TT is a
    **     specific linear transformation of geocentric coordinate time TCG,
    **     which is the time scale for the Geocentric Celestial Reference
    **     System, GCRS.
    **
    **  3) TDB is a coordinate time, and is a specific linear transformation
    **     of barycentric coordinate time TCB, which is the time scale for
    **     the Barycentric Celestial Reference System, BCRS.
    **
    **  4) The difference TCG-TCB depends on the masses and positions of the
    **     bodies of the solar system and the velocity of the Earth.  It is
    **     dominated by a rate difference, the residual being of a periodic
    **     character.  The latter, which is modeled by the present function,
    **     comprises a main (annual) sinusoidal term of amplitude
    **     approximately 0.00166 seconds, plus planetary terms up to about
    **     20 microseconds, and lunar and diurnal terms up to 2 microseconds.
    **     These effects come from the changing transverse Doppler effect
    **     and gravitational red-shift as the observer (on the Earth's
    **     surface) experiences variations in speed (with respect to the
    **     BCRS) and gravitational potential.
    **
    **  5) TDB can be regarded as the same as TCB but with a rate adjustment
    **     to keep it close to TT, which is convenient for many applications.
    **     The history of successive attempts to define TDB is set out in
    **     Resolution 3 adopted by the IAU General Assembly in 2006, which
    **     defines a fixed TDB(TCB) transformation that is consistent with
    **     contemporary solar-system ephemerides.  Future ephemerides will
    **     imply slightly changed transformations between TCG and TCB, which
    **     could introduce a linear drift between TDB and TT;  however, any
    **     such drift is unlikely to exceed 1 nanosecond per century.
    **
    **  6) The geocentric TDB-TT model used in the present function is that of
    **     Fairhead & Bretagnon (1990), in its full form.  It was originally
    **     supplied by Fairhead (private communications with P.T.Wallace,
    **     1990) as a Fortran subroutine.  The present C function contains an
    **     adaptation of the Fairhead code.  The numerical results are
    **     essentially unaffected by the changes, the differences with
    **     respect to the Fairhead & Bretagnon original being at the 1e-20 s
    **     level.
    **
    **     The topocentric part of the model is from Moyer (1981) and
    **     Murray (1983), with fundamental arguments adapted from
    **     Simon et al. 1994.  It is an approximation to the expression
    **     ( v / c ) . ( r / c ), where v is the barycentric velocity of
    **     the Earth, r is the geocentric position of the observer and
    **     c is the speed of light.
    **
    **     By supplying zeroes for u and v, the topocentric part of the
    **     model can be nullified, and the function will return the Fairhead
    **     & Bretagnon result alone.
    **
    **  7) During the interval 1950-2050, the absolute accuracy is better
    **     than +/- 3 nanoseconds relative to time ephemerides obtained by
    **     direct numerical integrations based on the JPL DE405 solar system
    **     ephemeris.
    **
    **  8) It must be stressed that the present function is merely a model,
    **     and that numerical integration of solar-system ephemerides is the
    **     definitive method for predicting the relationship between TCG and
    **     TCB and hence between TT and TDB.
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
    int iauAf2a(char s, int ideg, int iamin, double asec, double *rad)

    **  Given:
    **     s         char    sign:  '-' = negative, otherwise positive
    **     ideg      int     degrees
    **     iamin     int     arcminutes
    **     asec      double  arcseconds
    **
    **  Returned:
    **     rad       double  angle in radians
    **
    **  Returned (function value):
    **               int     status:  0 = OK
    **                                1 = ideg outside range 0-359
    **                                2 = iamin outside range 0-59
    **                                3 = asec outside range 0-59.999...
    **
    **  Notes:
    **
    **  1)  The result is computed even if any of the range checks fail.
    **
    **  2)  Negative ideg, iamin and/or asec produce a warning status, but
    **      the absolute value is used in the conversion.
    **
    **  3)  If there are multiple errors, the status value reflects only the
    **      first, the smallest taking precedence.
    """
    cdef double rad
    s = ord(sign)

    warns = {1: 'ideg outside range 0-359',
             2: 'iamin outside range 0-59',
             3: 'asec outside range 0-59.999...'}

    ret = iauAf2a(s, ideg, iamin, asec, &rad)
    check_return(ret, 'iauAf2a', warns)

    return rad

def iau_gd2gc(n, elong, phi, height):
    """
    Wrap
    int iauGd2gc(int n, double elong, double phi, double height, double xyz[3])

    **  Given:
    **     n       int        ellipsoid identifier (Note 1)
    **     elong   double     longitude (radians, east +ve)
    **     phi     double     latitude (geodetic, radians, Note 3)
    **     height  double     height above ellipsoid (geodetic, Notes 2,3)
    **
    **  Returned:
    **     xyz     double[3]  geocentric vector (Note 2)
    **
    **  Returned (function value):
    **             int        status:  0 = OK
    **                                -1 = illegal identifier (Note 3)
    **                                -2 = illegal case (Note 3)
    **
    **  Notes:
    **
    **  1) The identifier n is a number that specifies the choice of
    **     reference ellipsoid.  The following are supported:
    **
    **        n    ellipsoid
    **
    **        1     WGS84
    **        2     GRS80
    **        3     WGS72
    **
    **     The n value has no significance outside the SOFA software.  For
    **     convenience, symbols WGS84 etc. are defined in sofam.h.
    **
    **  2) The height (height, given) and the geocentric vector (xyz,
    **     returned) are in meters.
    **
    **  3) No validation is performed on the arguments elong, phi and
    **     height.  An error status -1 means that the identifier n is
    **     illegal.  An error status -2 protects against cases that would
    **     lead to arithmetic exceptions.  In all error cases, xyz is set
    **     to zeros.
    **
    **  4) The inverse transformation is performed in the function iauGc2gd.
    """
    cdef np.ndarray[double, ndim=1] xyz = np.empty(3, dtype=np.double)

    errs = {-1: 'illegal identifier',
             -2: 'illegal case'}

    ret = iauGd2gc(n, elong, phi, height, &xyz[0])
    check_return(ret, 'iauGd2gc', errors=errs)

    return xyz


@cython.wraparound(False)
@cython.boundscheck(False)
def jd_julian_epoch(np.ndarray[double, ndim=1] jd1,
                    np.ndarray[double, ndim=1] jd2):
    """ Wrap double iauEpj(double dj1, double dj2)
    **  Julian Date to Julian Epoch.

    **  Given:
    **     dj1,dj2    double     Julian Date (see note)
    **
    **  Returned (function value):
    **                double     Julian Epoch
    **
    **  Note:
    **
    **     The Julian Date is supplied in two pieces, in the usual SOFA
    **     manner, which is designed to preserve time resolution.  The
    **     Julian Date is available as a single number by adding dj1 and
    **     dj2.  The maximum resolution is achieved if dj1 is 2451545D0
    **     (J2000.0).
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
    **  Given:
    **     epj      double    Julian Epoch (e.g. 1996.8D0)
    **
    **  Returned:
    **     djm0     double    MJD zero-point: always 2400000.5
    **     djm      double    Modified Julian Date
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

    **  Given:
    **     dj1,dj2    double     Julian Date (see note)
    **
    **  Returned (function value):
    **                double     Besselian Epoch.
    **
    **  Note:
    **
    **     The Julian Date is supplied in two pieces, in the usual SOFA
    **     manner, which is designed to preserve time resolution.  The
    **     Julian Date is available as a single number by adding dj1 and
    **     dj2.  The maximum resolution is achieved if dj1 is 2451545D0
    **     (J2000.0).
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

    **  Given:
    **     epb      double    Besselian Epoch (e.g. 1957.3D0)
    **
    **  Returned:
    **     djm0     double    MJD zero-point: always 2400000.5
    **     djm      double    Modified Julian Date
    **
    **  Note:
    **
    **     The Julian Date is returned in two pieces, in the usual SOFA
    **     manner, which is designed to preserve time resolution.  The
    **     Julian Date is available as a single number by adding djm0 and
    **     djm.
    """
    cdef unsigned n = epd.shape[0]
    cdef unsigned int i
    cdef np.ndarray[double, ndim=1] jd1 = np.empty(n, dtype=np.double)
    cdef np.ndarray[double, ndim=1] jd2 = np.empty(n, dtype=np.double)

    for i in range(n):
        iauEpb2jd(epd[i], &jd1[i], &jd2[i])
    return jd1, jd2
