import warnings

from ..utils.exceptions import AstropyUserWarning

import numpy as np
cimport numpy as np
import cython

ctypedef np.double_t DOUBLE_T

cdef extern from "erfa.h":
    double eraEpb(double dj1, double dj2)
    double eraEpj(double dj1, double dj2)
    void eraEpj2jd(double epj, double *djm0, double *djm)
    void eraEpb2jd(double epb, double *djm0, double *djm)
    int eraCal2jd(int iy, int im, int id, double *djm0, double *djm)
    int eraJd2cal(double dj1, double dj2,int *iy, int *im, int *id, double *fd)
    int eraJdcalf(int ndp, double dj1, double dj2, int iymdf[4])

    # Internal JD to/from datetime format converters
    int eraD2dtf(char *scale, int ndp, double d1, double d2, int *iy, int *im, int *id, int ihmsf[4])
    int eraDtf2d(char *scale, int iy, int im, int id, int ihr, int imn, double sec, double *d1, double *d2)

    # Time scale helper routines
    double eraDtdb(double date1, double date2, double ut, double elong, double u, double v)
    int eraDat(int iy, int im, int id, double fd, double *deltat)

    # Time scale conversion routines
    int eraTaitt(double tai1, double tai2, double *tt1, double *tt2)
    int eraTttai(double tt1, double tt2, double *tai1, double *tai2)
    int eraTaiutc(double tai1, double tai2, double *utc1, double *utc2)
    int eraUtctai(double utc1, double utc2, double *tai1, double *tai2)
    int eraTcbtdb(double tcb1, double tcb2, double *tdb1, double *tdb2)
    int eraTdbtcb(double tdb1, double tdb2, double *tcb1, double *tcb2)
    int eraTcgtt(double tcg1, double tcg2, double *tt1, double *tt2)
    int eraTttcg(double tt1, double tt2, double *tcg1, double *tcg2)

    int eraTaiut1(double tai1, double tai2, double dta, double *ut11, double *ut12)
    int eraUt1tai(double ut11, double ut12, double dta, double *tai1, double *tai2)
    int eraTtut1(double tt1, double tt2, double dt, double *ut11, double *ut12)
    int eraUt1tt(double ut11, double ut12, double dt, double *tt1, double *tt2)
    int eraTdbtt(double tdb1, double tdb2, double dtr, double *tt1, double *tt2)
    int eraTttdb(double tt1, double tt2, double dtr, double *tdb1, double *tdb2)
    int eraUt1utc(double ut11, double ut12, double dut1, double *utc1, double *utc2)
    int eraUtcut1(double utc1, double utc2, double dut1, double *ut11, double *ut12)

    # Geodetic
    int eraAf2a(char s, int ideg, int iamin, double asec, double *rad)
    int eraGd2gc(int n, double elong, double phi, double height, double xyz[3])
    int eraGc2gd(int n, double xyz[3], double *elong, double *phi, double *height )

    # Sidereal time
    double eraGmst06(double uta, double utb, double tta, double ttb)
    double eraGmst00(double uta, double utb, double tta, double ttb)
    double eraGmst82(double dj1, double dj2)
    double eraGst00a(double uta, double utb, double tta, double ttb)
    double eraGst00b(double uta, double utb)
    double eraGst06a(double uta, double utb, double tta, double ttb)
    double eraGst94(double uta, double utb)

DUBIOUS = 'dubious year for UTC (before 1960.0 or 5 years ' \
          'beyond last known leap second)'

def check_return(ret, func_name, warns={}, errors={}):
    """Check the return value from an era routine"""
    if ret in warns:
        warnings.warn('{0}: {1}'.format(func_name, warns[ret]), AstropyUserWarning)
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
    int eraCal2jd(int iy, int im, int id, double *djm0, double *djm)
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
    **  2) The Julian Date is returned in two pieces, in the usual ERFA
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
        ret = eraCal2jd( iy[i], im[i], id[i], &djm0[i], &djm[i])
        check_return(ret, 'eraCal2jd', warns, errs)
    return

@cython.wraparound(False)
@cython.boundscheck(False)
def d_tai_utc(np.ndarray[int, ndim=1] iy,
              np.ndarray[int, ndim=1] im,
              np.ndarray[int, ndim=1] id,
              np.ndarray[double, ndim=1] fd):
    """
    int eraDat(int iy, int im, int id, double fd, double *deltat)
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
        ret = eraDat(iy[i], im[i], id[i], fd[i],
                     &out[i])
        check_return(ret, 'eraDat', warns, errs)

    return out


@cython.wraparound(False)
@cython.boundscheck(False)
def jd_dtf(scale, ndp,
              np.ndarray[double, ndim=1] d1,
              np.ndarray[double, ndim=1] d2):
    """
    int eraD2dtf(const char *scale, int ndp, double d1, double d2,
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
    **     special measures are taken.  The ERFA internal convention is that
    **     the quasi-JD day represents UTC days whether the length is 86399,
    **     86400 or 86401 SI seconds.
    **
    **  5) The warning status "dubious year" flags UTCs that predate the
    **     introduction of the time scale and that are too far in the future
    **     to be trusted.  See eraDat for further details.
    **
    **  6) For calendar conventions and limitations, see eraCal2jd.
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
        ret = eraD2dtf(scale, ndp, d1[i], d2[i],
                     &iy[i], &im[i], &id[i], &ihmsf[i, 0])
        check_return(ret, 'eraD2dtf', warns, errs)

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
    int eraDtf2d(char *scale, int iy, int im, int id, int ihr, int imn, double sec, double *d1, double *d2)

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
    **  2) For calendar conventions and limitations, see eraCal2jd.
    **
    **  3) The sum of the results, d1+d2, is Julian Date, where normally d1
    **     is the Julian Day Number and d2 is the fraction of a day.  In the
    **     case of UTC, where the use of JD is problematical, special
    **     conventions apply:  see the next note.
    **
    **  4) JD cannot unambiguously represent UTC during a leap second unless
    **     special measures are taken.  The ERFA internal convention is that
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
    **     to be trusted.  See eraDat for further details.
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
        ret = eraDtf2d(scale, iy[i], im[i], id[i], ihr[i], imn[i], sec[i],
                       &out1[i], &out2[i])
        check_return(ret, 'eraDtf2d', warns, errs)

    return out1, out2


@cython.wraparound(False)
@cython.boundscheck(False)
def tai_tt( 
    np.ndarray[double, ndim=1] in1,
    np.ndarray[double, ndim=1] in2):
    """
    int eraTaitt(double tai1, double tai2, double *tt1, double *tt2)

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
        ret = eraTaitt(in1[i], in2[i], &out1[i], &out2[i])
        check_return(ret, 'eraTaitt')

    return out1, out2


@cython.wraparound(False)
@cython.boundscheck(False)
def tcb_tdb( 
    np.ndarray[double, ndim=1] in1,
    np.ndarray[double, ndim=1] in2):
    """
    int eraTcbtdb(double tcb1, double tcb2, double *tdb1, double *tdb2)

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
        ret = eraTcbtdb(in1[i], in2[i], &out1[i], &out2[i])
        check_return(ret, 'eraTcbtdb')

    return out1, out2


@cython.wraparound(False)
@cython.boundscheck(False)
def tcg_tt( 
    np.ndarray[double, ndim=1] in1,
    np.ndarray[double, ndim=1] in2):
    """
   int eraTcgtt(double tcg1, double tcg2, double *tt1, double *tt2)

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
        ret = eraTcgtt(in1[i], in2[i], &out1[i], &out2[i])
        check_return(ret, 'eraTcgtt')

    return out1, out2


@cython.wraparound(False)
@cython.boundscheck(False)
def tdb_tcb( 
    np.ndarray[double, ndim=1] in1,
    np.ndarray[double, ndim=1] in2):
    """
    int eraTdbtcb(double tdb1, double tdb2, double *tcb1, double *tcb2)

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
        ret = eraTdbtcb(in1[i], in2[i], &out1[i], &out2[i])
        check_return(ret, 'eraTdbtcb')

    return out1, out2


@cython.wraparound(False)
@cython.boundscheck(False)
def tt_tai( 
    np.ndarray[double, ndim=1] in1,
    np.ndarray[double, ndim=1] in2):
    """
    int eraTttai(double tt1, double tt2, double *tai1, double *tai2)

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
        ret = eraTttai(in1[i], in2[i], &out1[i], &out2[i])
        check_return(ret, 'eraTttai')

    return out1, out2


@cython.wraparound(False)
@cython.boundscheck(False)
def tt_tcg( 
    np.ndarray[double, ndim=1] in1,
    np.ndarray[double, ndim=1] in2):
    """
    int eraTttcg(double tt1, double tt2, double *tcg1, double *tcg2)
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
        ret = eraTttcg(in1[i], in2[i], &out1[i], &out2[i])
        check_return(ret, 'eraTttcg')

    return out1, out2


@cython.wraparound(False)
@cython.boundscheck(False)
def utc_tai( 
    np.ndarray[double, ndim=1] in1,
    np.ndarray[double, ndim=1] in2):
    """
    int eraUtctai(double utc1, double utc2, double *tai1, double *tai2)

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
    **     to be trusted.  See eraDat  for further details.
    **
    **  4) The function eraDtf2d converts from calendar date and time of day
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
        ret = eraUtctai(in1[i], in2[i], &out1[i], &out2[i])
        check_return(ret, 'eraUtctai', warns, errs)

    return out1, out2


@cython.wraparound(False)
@cython.boundscheck(False)
def tai_utc( 
    np.ndarray[double, ndim=1] in1,
    np.ndarray[double, ndim=1] in2):
    """
    int eraTaiutc(double tai1, double tai2, double *utc1, double *utc2)

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
    **  3) The function eraD2dtf can be used to transform the UTC quasi-JD
    **     into calendar date and clock time, including UTC leap second
    **     handling.
    **
    **  4) The warning status "dubious year" flags UTCs that predate the
    **     introduction of the time scale and that are too far in the future
    **     to be trusted.  See eraDat for further details.
        """
    assert in1.shape[0] == in2.shape[0]
    cdef unsigned n = in1.shape[0]
    cdef unsigned int i
    cdef np.ndarray[double, ndim=1] out1 = np.empty(n, dtype=np.double)
    cdef np.ndarray[double, ndim=1] out2 = np.empty(n, dtype=np.double)

    warns = {1: DUBIOUS}
    errs = {-1: 'unacceptable date'}

    for i in range(n):
        ret = eraTaiutc(in1[i], in2[i], &out1[i], &out2[i])
        check_return(ret, 'eraTaiutc', warns, errs)

    return out1, out2


@cython.wraparound(False)
@cython.boundscheck(False)
def tai_ut1( 
    np.ndarray[double, ndim=1] in1,
    np.ndarray[double, ndim=1] in2,
    np.ndarray[double, ndim=1] dt):
    """
    int eraTaiut1(double tai1, double tai2, double dta, double *ut11, double *ut12)

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
        ret = eraTaiut1(in1[i], in2[i], dt[i], &out1[i], &out2[i])
        check_return(ret, 'eraTaiut1')

    return out1, out2


@cython.wraparound(False)
@cython.boundscheck(False)
def ut1_tai( 
    np.ndarray[double, ndim=1] in1,
    np.ndarray[double, ndim=1] in2,
    np.ndarray[double, ndim=1] dt):
    """
    int eraUt1tai(double ut11, double ut12, double dta, double *tai1, double *tai2)

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
        ret = eraUt1tai(in1[i], in2[i], dt[i], &out1[i], &out2[i])
        check_return(ret, 'eraUt1tai')

    return out1, out2


@cython.wraparound(False)
@cython.boundscheck(False)
def tt_ut1( 
    np.ndarray[double, ndim=1] in1,
    np.ndarray[double, ndim=1] in2,
    np.ndarray[double, ndim=1] dt):
    """
    int eraTtut1(double tt1, double tt2, double dt, double *ut11, double *ut12)

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
        ret = eraTtut1(in1[i], in2[i], dt[i], &out1[i], &out2[i])
        check_return(ret, 'eraTtut1')

    return out1, out2


@cython.wraparound(False)
@cython.boundscheck(False)
def ut1_tt( 
    np.ndarray[double, ndim=1] in1,
    np.ndarray[double, ndim=1] in2,
    np.ndarray[double, ndim=1] dt):
    """
    int eraUt1tt(double ut11, double ut12, double dt, double *tt1, double *tt2)

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
        ret = eraUt1tt(in1[i], in2[i], dt[i], &out1[i], &out2[i])
        check_return(ret, 'eraUt1tt')

    return out1, out2


@cython.wraparound(False)
@cython.boundscheck(False)
def tdb_tt( 
    np.ndarray[double, ndim=1] in1,
    np.ndarray[double, ndim=1] in2,
    np.ndarray[double, ndim=1] dt):
    """
    int eraTdbtt(double tdb1, double tdb2, double dtr, double *tt1, double *tt2)

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
    **     evaluating a model such as that implemented in the ERFA function
    **     eraDtdb.   The quantity is dominated by an annual term of 1.7 ms
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
        ret = eraTdbtt(in1[i], in2[i], dt[i], &out1[i], &out2[i])
        check_return(ret, 'eraTdbtt')

    return out1, out2


@cython.wraparound(False)
@cython.boundscheck(False)
def tt_tdb( 
    np.ndarray[double, ndim=1] in1,
    np.ndarray[double, ndim=1] in2,
    np.ndarray[double, ndim=1] dt):
    """
    int eraTttdb(double tt1, double tt2, double dtr, double *tdb1, double *tdb2)

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
    **     evaluating a model such as that implemented in the ERFA function
    **     eraDtdb.   The quantity is dominated by an annual term of 1.7 ms
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
        ret = eraTttdb(in1[i], in2[i], dt[i], &out1[i], &out2[i])
        check_return(ret, 'eraTttdb')

    return out1, out2


@cython.wraparound(False)
@cython.boundscheck(False)
def ut1_utc( 
    np.ndarray[double, ndim=1] in1,
    np.ndarray[double, ndim=1] in2,
    np.ndarray[double, ndim=1] dt):
    """
    int eraUt1utc(double ut11, double ut12, double dut1, double *utc1, double *utc2)

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
    **  4) The function eraD2dtf can be used to transform the UTC quasi-JD
    **     into calendar date and clock time, including UTC leap second
    **     handling.
    **
    **  5) The warning status "dubious year" flags UTCs that predate the
    **     introduction of the time scale and that are too far in the future
    **     to be trusted.  See eraDat for further details.
    """
    assert in1.shape[0] == in2.shape[0] == dt.shape[0]
    cdef unsigned n = in1.shape[0]
    cdef unsigned int i
    cdef np.ndarray[double, ndim=1] out1 = np.empty(n, dtype=np.double)
    cdef np.ndarray[double, ndim=1] out2 = np.empty(n, dtype=np.double)

    warns = {1: DUBIOUS}
    errs = {-1: 'unacceptable date'}

    for i in range(n):
        ret = eraUt1utc(in1[i], in2[i], dt[i], &out1[i], &out2[i])
        check_return(ret, 'eraUt1utc', warns, errs)

    return out1, out2


@cython.wraparound(False)
@cython.boundscheck(False)
def utc_ut1( 
    np.ndarray[double, ndim=1] in1,
    np.ndarray[double, ndim=1] in2,
    np.ndarray[double, ndim=1] dt):
    """
    int eraUtcut1(double utc1, double utc2, double dut1, double *ut11, double *ut12)

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
    **     to be trusted.  See eraDat  for further details.
    **
    **  4) The function eraDtf2d  converts from calendar date and time of
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
    **     to be trusted.  See eraDat for further details.
    """
    assert in1.shape[0] == in2.shape[0] == dt.shape[0]
    cdef unsigned n = in1.shape[0]
    cdef unsigned int i
    cdef np.ndarray[double, ndim=1] out1 = np.empty(n, dtype=np.double)
    cdef np.ndarray[double, ndim=1] out2 = np.empty(n, dtype=np.double)

    warns = {1: DUBIOUS}
    errs = {-1: 'unacceptable date'}

    for i in range(n):
        ret = eraUtcut1(in1[i], in2[i], dt[i], &out1[i], &out2[i])
        check_return(ret, 'eraUtcut1', warns, errs)

    return out1, out2


@cython.wraparound(False)
@cython.boundscheck(False)
def d_tdb_tt(np.ndarray[double, ndim=1] in1,
             np.ndarray[double, ndim=1] in2,
             np.ndarray[double, ndim=1] ut,
             np.ndarray[double, ndim=1] elong,
	     np.ndarray[double, ndim=1] u,
	     np.ndarray[double, ndim=1] v):
    """
    compute DTR = TDB-TT
    double eraDtdb(double date1, double date2, double ut,
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
    assert elong.shape[0] == u.shape[0] == v.shape[0]
    assert elong.shape[0] == 1 or elong.shape[0] == in1.shape[0]
    cdef unsigned n = in1.shape[0]
    cdef unsigned int i, j
    cdef np.ndarray[double, ndim=1] out = np.empty(n, dtype=np.double)

    for i in range(n):
        j = min(i, elong.shape[0]-1)
        out[i] = eraDtdb(in1[i], in2[i], ut[i], elong[j], u[j], v[j])
    return out


def era_af2a(sign, ideg, iamin, asec):
    """
    int eraAf2a(char s, int ideg, int iamin, double asec, double *rad)

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

    ret = eraAf2a(s, ideg, iamin, asec, &rad)
    check_return(ret, 'eraAf2a', warns)

    return rad

def era_gd2gc(n, elong, phi, height):
    """
    Wrap
    int eraGd2gc(int n, double elong, double phi, double height, double xyz[3])

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
    **     The n value has no significance outside the ERFA software.  For
    **     convenience, symbols WGS84 etc. are defined in erfam.h.
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
    **  4) The inverse transformation is performed in the function eraGc2gd.
    """
    assert elong.shape[0] == phi.shape[0] == height.shape[0]
    cdef unsigned int i
    cdef unsigned int nitems = elong.shape[0]
    cdef np.ndarray[double, ndim=1] xyz = np.empty(3, dtype=np.double)
    cdef np.ndarray[double, ndim=2] out = np.empty((nitems, 3), dtype=np.double)

    errs = {-1: 'illegal identifier',
             -2: 'illegal case'}

    for i in range(nitems):
        ret = eraGd2gc(n, elong[i], phi[i], height[i], &xyz[0])
        check_return(ret, 'eraGd2gc', errors=errs)
        out[i] = xyz

    return out

@cython.wraparound(False)
@cython.boundscheck(False)
def era_gc2gd(n, xyz):
    """
    Wrap
    int eraGc2gd(int n, double xyz[3], double *elong, double *phi, double *height )

    **  Given:
    **     n       int        ellipsoid identifier (Note 1)
    **     xyz     double[3]  geocentric vector (Note 2)
    **
    **  Returned:
    **     elong   double     longitude (radians, east +ve)
    **     phi     double     latitude (geodetic, radians, Note 3)
    **     height  double     height above ellipsoid (geodetic, Notes 2,3)
    **
    **  Returned (function value):
    **            int         status:  0 = OK
    **                                -1 = illegal identifier (Note 3)
    **                                -2 = internal error (Note 3)
    **
    **  Notes:
    **
    **  1) The identifier n is a number that specifies the choice of
    **     reference ellipsoid.  The following are supported:
    **
    **        n    ellipsoid
    **
    **        1     ERFA_WGS84
    **        2     ERFA_GRS80
    **        3     ERFA_WGS72
    **
    **     The n value has no significance outside the ERFA software.  For
    **     convenience, symbols ERFA_WGS84 etc. are defined in erfam.h.
    **
    **  2) The geocentric vector (xyz, given) and height (height, returned)
    **     are in meters.
    **
    **  3) An error status -1 means that the identifier n is illegal.  An
    **     error status -2 is theoretically impossible.  In all error cases,
    **     phi and height are both set to -1e9.
    **
    **  4) The inverse transformation is performed in the function eraGd2gc.
    **
    **  Called:
    **     eraEform     Earth reference ellipsoids
    **     eraGc2gde    geocentric to geodetic transformation, general
    **
    **  Copyright (C) 2013, NumFOCUS Foundation.
    **  Derived, with permission, from the SOFA library.  See notes at end of file.
    """
    assert xyz.shape[1] == 3
    cdef unsigned int i
    cdef unsigned int nitems = xyz.shape[0]
    cdef np.ndarray[double, ndim=1] elong = np.empty(nitems, dtype=np.double)
    cdef np.ndarray[double, ndim=1] phi = np.empty(nitems, dtype=np.double)
    cdef np.ndarray[double, ndim=1] height = np.empty(nitems, dtype=np.double)
    cdef double xyz_item[3]

    errs = {-1: 'illegal identifier',
            -2: 'illegal case'}
    for i in range(nitems):
        # ensure xyz are in a contiguous array
        for j in range(3):
            xyz_item[j] = xyz[i, j]
        ret = eraGc2gd(n, xyz_item, &elong[i], &phi[i], &height[i])
        check_return(ret, 'eraGd2gc', errors=errs)

    return elong, phi, height

@cython.wraparound(False)
@cython.boundscheck(False)
def jd_julian_epoch(np.ndarray[double, ndim=1] jd1,
                    np.ndarray[double, ndim=1] jd2):
    """ Wrap double eraEpj(double dj1, double dj2)
    **  Julian Date to Julian Epoch.

    **  Given:
    **     dj1,dj2    double     Julian Date (see note)
    **
    **  Returned (function value):
    **                double     Julian Epoch
    **
    **  Note:
    **
    **     The Julian Date is supplied in two pieces, in the usual ERFA
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
        epd[i] = eraEpj(jd1[i], jd2[i])
    return epd


@cython.wraparound(False)
@cython.boundscheck(False)
def julian_epoch_jd(np.ndarray[double, ndim=1] epd):
    """ Wrap void eraEpj2jd(double epj, double *djm0, double *djm)
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
        eraEpj2jd(epd[i], &jd1[i], &jd2[i])
    return jd1, jd2


@cython.wraparound(False)
@cython.boundscheck(False)
def jd_besselian_epoch(np.ndarray[double, ndim=1] jd1,
                       np.ndarray[double, ndim=1] jd2):
    """ Wrap double eraEpb(double dj1, double dj2)
    **  Julian Date to Besselian Epoch.

    **  Given:
    **     dj1,dj2    double     Julian Date (see note)
    **
    **  Returned (function value):
    **                double     Besselian Epoch.
    **
    **  Note:
    **
    **     The Julian Date is supplied in two pieces, in the usual ERFA
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
        epd[i] = eraEpb(jd1[i], jd2[i])
    return epd


@cython.wraparound(False)
@cython.boundscheck(False)
def besselian_epoch_jd(np.ndarray[double, ndim=1] epd):
    """ Wrap void eraEpb2jd(double epj, double *djm0, double *djm)
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
    **     The Julian Date is returned in two pieces, in the usual ERFA
    **     manner, which is designed to preserve time resolution.  The
    **     Julian Date is available as a single number by adding djm0 and
    **     djm.
    """
    cdef unsigned n = epd.shape[0]
    cdef unsigned int i
    cdef np.ndarray[double, ndim=1] jd1 = np.empty(n, dtype=np.double)
    cdef np.ndarray[double, ndim=1] jd2 = np.empty(n, dtype=np.double)

    for i in range(n):
        eraEpb2jd(epd[i], &jd1[i], &jd2[i])
    return jd1, jd2

@cython.wraparound(False)
@cython.boundscheck(False)
def gmst00(np.ndarray[double, ndim=1] ut11,
           np.ndarray[double, ndim=1] ut12,
           np.ndarray[double, ndim=1] tt1,
           np.ndarray[double, ndim=1] tt2):
    """Wrap double eraGmst00(double uta, double utb, double tta, double ttb)
    **  Greenwich mean sidereal time (model consistent with IAU 2000
    **  resolutions).
    **
    **  Given:
    **     uta,utb    double    UT1 as a 2-part Julian Date (Notes 1,2)
    **     tta,ttb    double    TT as a 2-part Julian Date (Notes 1,2)
    **
    **  Returned (function value):
    **                double    Greenwich mean sidereal time (radians)
    **
    **  Notes:
    **
    **  1) The UT1 and TT dates uta+utb and tta+ttb respectively, are both
    **     Julian Dates, apportioned in any convenient way between the
    **     argument pairs.  For example, JD=2450123.7 could be expressed in
    **     any of these ways, among others:
    **
    **            Part A         Part B
    **
    **         2450123.7           0.0       (JD method)
    **         2451545.0       -1421.3       (J2000 method)
    **         2400000.5       50123.2       (MJD method)
    **         2450123.5           0.2       (date & time method)
    **
    **     The JD method is the most natural and convenient to use in
    **     cases where the loss of several decimal digits of resolution
    **     is acceptable (in the case of UT;  the TT is not at all critical
    **     in this respect).  The J2000 and MJD methods are good compromises
    **     between resolution and convenience.  For UT, the date & time
    **     method is best matched to the algorithm that is used by the Earth
    **     Rotation Angle function, called internally:  maximum precision is
    **     delivered when the uta argument is for 0hrs UT1 on the day in
    **     question and the utb argument lies in the range 0 to 1, or vice
    **     versa.
    **
    **  2) Both UT1 and TT are required, UT1 to predict the Earth rotation
    **     and TT to predict the effects of precession.  If UT1 is used for
    **     both purposes, errors of order 100 microarcseconds result.
    **
    **  3) This GMST is compatible with the IAU 2000 resolutions and must be
    **     used only in conjunction with other IAU 2000 compatible
    **     components such as precession-nutation and equation of the
    **     equinoxes.
    **
    **  4) The result is returned in the range 0 to 2pi.
    **
    **  5) The algorithm is from Capitaine et al. (2003) and IERS
    **     Conventions 2003.
    **
    **  Called:
    **     eraEra00     Earth rotation angle, IAU 2000
    **     eraAnp       normalize angle into range 0 to 2pi
    **
    **  References:
    **
    **     Capitaine, N., Wallace, P.T. and McCarthy, D.D., "Expressions to
    **     implement the IAU 2000 definition of UT1", Astronomy &
    **     Astrophysics, 406, 1135-1149 (2003)
    **
    **     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
    **     IERS Technical Note No. 32, BKG (2004)
    **
    **  Copyright (C) 2013, NumFOCUS Foundation.
    **  Derived, with permission, from the SOFA library.  See notes at end of file.
    """
    assert ut11.shape[0] == ut12.shape[0] == tt1.shape[0] == tt2.shape[0]
    cdef unsigned n = ut11.shape[0]
    cdef unsigned int i
    cdef np.ndarray[double, ndim=1] gmst = np.empty(n, dtype=np.double)

    for i in range(n):
        gmst[i] = eraGmst00(ut11[i], ut12[i], tt1[i], tt2[i])

    return gmst

@cython.wraparound(False)
@cython.boundscheck(False)
def gmst06(np.ndarray[double, ndim=1] ut11,
           np.ndarray[double, ndim=1] ut12,
           np.ndarray[double, ndim=1] tt1,
           np.ndarray[double, ndim=1] tt2):
    """Wrap double eraGmst06(double uta, double utb, double tta, double ttb)
    **  Greenwich mean sidereal time (consistent with IAU 2006 precession).
    **
    **  Given:
    **     uta,utb    double    UT1 as a 2-part Julian Date (Notes 1,2)
    **     tta,ttb    double    TT as a 2-part Julian Date (Notes 1,2)
    **
    **  Returned (function value):
    **                double    Greenwich mean sidereal time (radians)
    **
    **  Notes:
    **
    **  1) The UT1 and TT dates uta+utb and tta+ttb respectively, are both
    **     Julian Dates, apportioned in any convenient way between the
    **     argument pairs.  For example, JD=2450123.7 could be expressed in
    **     any of these ways, among others:
    **
    **            Part A        Part B
    **
    **         2450123.7           0.0       (JD method)
    **         2451545.0       -1421.3       (J2000 method)
    **         2400000.5       50123.2       (MJD method)
    **         2450123.5           0.2       (date & time method)
    **
    **     The JD method is the most natural and convenient to use in
    **     cases where the loss of several decimal digits of resolution
    **     is acceptable (in the case of UT;  the TT is not at all critical
    **     in this respect).  The J2000 and MJD methods are good compromises
    **     between resolution and convenience.  For UT, the date & time
    **     method is best matched to the algorithm that is used by the Earth
    **     rotation angle function, called internally:  maximum precision is
    **     delivered when the uta argument is for 0hrs UT1 on the day in
    **     question and the utb argument lies in the range 0 to 1, or vice
    **     versa.
    **
    **  2) Both UT1 and TT are required, UT1 to predict the Earth rotation
    **     and TT to predict the effects of precession.  If UT1 is used for
    **     both purposes, errors of order 100 microarcseconds result.
    **
    **  3) This GMST is compatible with the IAU 2006 precession and must not
    **     be used with other precession models.
    **
    **  4) The result is returned in the range 0 to 2pi.
    **
    **  Called:
    **     eraEra00     Earth rotation angle, IAU 2000
    **     eraAnp       normalize angle into range 0 to 2pi
    **
    **  Reference:
    **
    **     Capitaine, N., Wallace, P.T. & Chapront, J., 2005,
    **     Astron.Astrophys. 432, 355
    """
    assert ut11.shape[0] == ut12.shape[0] == tt1.shape[0] == tt2.shape[0]
    cdef unsigned n = ut11.shape[0]
    cdef unsigned int i
    cdef np.ndarray[double, ndim=1] gmst = np.empty(n, dtype=np.double)

    for i in range(n):
        gmst[i] = eraGmst06(ut11[i], ut12[i], tt1[i], tt2[i])

    return gmst

@cython.wraparound(False)
@cython.boundscheck(False)
def gmst82(np.ndarray[double, ndim=1] ut11,
           np.ndarray[double, ndim=1] ut12):
    """Wrap double double eraGmst82(double dj1, double dj2)
    **  Universal Time to Greenwich mean sidereal time (IAU 1982 model).
    **
    **  Given:
    **     dj1,dj2    double    UT1 Julian Date (see note)
    **
    **  Returned (function value):
    **                double    Greenwich mean sidereal time (radians)
    **
    **  Notes:
    **
    **  1) The UT1 date dj1+dj2 is a Julian Date, apportioned in any
    **     convenient way between the arguments dj1 and dj2.  For example,
    **     JD(UT1)=2450123.7 could be expressed in any of these ways,
    **     among others:
    **
    **             dj1            dj2
    **
    **         2450123.7D0        0D0        (JD method)
    **          2451545D0      -1421.3D0     (J2000 method)
    **         2400000.5D0     50123.2D0     (MJD method)
    **         2450123.5D0       0.2D0       (date & time method)
    **
    **     The JD method is the most natural and convenient to use in
    **     cases where the loss of several decimal digits of resolution
    **     is acceptable.  The J2000 and MJD methods are good compromises
    **     between resolution and convenience.  The date & time method is
    **     best matched to the algorithm used:  maximum accuracy (or, at
    **     least, minimum noise) is delivered when the dj1 argument is for
    **     0hrs UT1 on the day in question and the dj2 argument lies in the
    **     range 0 to 1, or vice versa.
    **
    **  2) The algorithm is based on the IAU 1982 expression.  This is
    **     always described as giving the GMST at 0 hours UT1.  In fact, it
    **     gives the difference between the GMST and the UT, the steady
    **     4-minutes-per-day drawing-ahead of ST with respect to UT.  When
    **     whole days are ignored, the expression happens to equal the GMST
    **     at 0 hours UT1 each day.
    **
    **  3) In this function, the entire UT1 (the sum of the two arguments
    **     dj1 and dj2) is used directly as the argument for the standard
    **     formula, the constant term of which is adjusted by 12 hours to
    **     take account of the noon phasing of Julian Date.  The UT1 is then
    **     added, but omitting whole days to conserve accuracy.
    **
    **  Called:
    **     eraAnp       normalize angle into range 0 to 2pi
    **
    **  References:
    **
    **     Transactions of the International Astronomical Union,
    **     XVIII B, 67 (1983).
    **
    **     Aoki et al., Astron. Astrophys. 105, 359-361 (1982).
    """
    assert ut11.shape[0] == ut12.shape[0]
    cdef unsigned n = ut11.shape[0]
    cdef unsigned int i
    cdef np.ndarray[double, ndim=1] gmst = np.empty(n, dtype=np.double)

    for i in range(n):
        gmst[i] = eraGmst82(ut11[i], ut12[i]
)
    return gmst

@cython.wraparound(False)
@cython.boundscheck(False)
def gst00a(np.ndarray[double, ndim=1] ut11,
           np.ndarray[double, ndim=1] ut12,
           np.ndarray[double, ndim=1] tt1,
           np.ndarray[double, ndim=1] tt2):
    """Wrap double eraGst00a(double uta, double utb, double tta, double ttb)
    **  Greenwich apparent sidereal time (consistent with IAU 2000
    **  resolutions).
    **
    **  Given:
    **     uta,utb    double    UT1 as a 2-part Julian Date (Notes 1,2)
    **     tta,ttb    double    TT as a 2-part Julian Date (Notes 1,2)
    **
    **  Returned (function value):
    **                double    Greenwich apparent sidereal time (radians)
    **
    **  Notes:
    **
    **  1) The UT1 and TT dates uta+utb and tta+ttb respectively, are both
    **     Julian Dates, apportioned in any convenient way between the
    **     argument pairs.  For example, JD=2450123.7 could be expressed in
    **     any of these ways, among others:
    **
    **            Part A        Part B
    **
    **         2450123.7           0.0       (JD method)
    **         2451545.0       -1421.3       (J2000 method)
    **         2400000.5       50123.2       (MJD method)
    **         2450123.5           0.2       (date & time method)
    **
    **     The JD method is the most natural and convenient to use in
    **     cases where the loss of several decimal digits of resolution
    **     is acceptable (in the case of UT;  the TT is not at all critical
    **     in this respect).  The J2000 and MJD methods are good compromises
    **     between resolution and convenience.  For UT, the date & time
    **     method is best matched to the algorithm that is used by the Earth
    **     Rotation Angle function, called internally:  maximum precision is
    **     delivered when the uta argument is for 0hrs UT1 on the day in
    **     question and the utb argument lies in the range 0 to 1, or vice
    **     versa.
    **
    **  2) Both UT1 and TT are required, UT1 to predict the Earth rotation
    **     and TT to predict the effects of precession-nutation.  If UT1 is
    **     used for both purposes, errors of order 100 microarcseconds
    **     result.
    **
    **  3) This GAST is compatible with the IAU 2000 resolutions and must be
    **     used only in conjunction with other IAU 2000 compatible
    **     components such as precession-nutation.
    **
    **  4) The result is returned in the range 0 to 2pi.
    **
    **  5) The algorithm is from Capitaine et al. (2003) and IERS
    **     Conventions 2003.
    **
    **  Called:
    **     eraGmst00    Greenwich mean sidereal time, IAU 2000
    **     eraEe00a     equation of the equinoxes, IAU 2000A
    **     eraAnp       normalize angle into range 0 to 2pi
    **
    **  References:
    **
    **     Capitaine, N., Wallace, P.T. and McCarthy, D.D., "Expressions to
    **     implement the IAU 2000 definition of UT1", Astronomy &
    **     Astrophysics, 406, 1135-1149 (2003)
    **
    **     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
    **     IERS Technical Note No. 32, BKG (2004)
    """
    assert ut11.shape[0] == ut12.shape[0] == tt1.shape[0] == tt2.shape[0]
    cdef unsigned n = ut11.shape[0]
    cdef unsigned int i
    cdef np.ndarray[double, ndim=1] gst = np.empty(n, dtype=np.double)

    for i in range(n):
        gst[i] = eraGst00a(ut11[i], ut12[i], tt1[i], tt2[i])

    return gst

@cython.wraparound(False)
@cython.boundscheck(False)
def gst00b(np.ndarray[double, ndim=1] ut11,
            np.ndarray[double, ndim=1] ut12):
    """Wrap double eraGst00b(double uta, double utb)
    **  Greenwich apparent sidereal time (consistent with IAU 2000
    **  resolutions but using the truncated nutation model IAU 2000B).
    **
    **  Given:
    **     uta,utb    double    UT1 as a 2-part Julian Date (Notes 1,2)
    **
    **  Returned (function value):
    **                double    Greenwich apparent sidereal time (radians)
    **
    **  Notes:
    **
    **  1) The UT1 date uta+utb is a Julian Date, apportioned in any
    **     convenient way between the argument pair.  For example,
    **     JD=2450123.7 could be expressed in any of these ways, among
    **     others:
    **
    **             uta            utb
    **
    **         2450123.7           0.0       (JD method)
    **         2451545.0       -1421.3       (J2000 method)
    **         2400000.5       50123.2       (MJD method)
    **         2450123.5           0.2       (date & time method)
    **
    **     The JD method is the most natural and convenient to use in cases
    **     where the loss of several decimal digits of resolution is
    **     acceptable.  The J2000 and MJD methods are good compromises
    **     between resolution and convenience.  For UT, the date & time
    **     method is best matched to the algorithm that is used by the Earth
    **     Rotation Angle function, called internally:  maximum precision is
    **     delivered when the uta argument is for 0hrs UT1 on the day in
    **     question and the utb argument lies in the range 0 to 1, or vice
    **     versa.
    **
    **  2) The result is compatible with the IAU 2000 resolutions, except
    **     that accuracy has been compromised for the sake of speed and
    **     convenience in two respects:
    **
    **     . UT is used instead of TDB (or TT) to compute the precession
    **       component of GMST and the equation of the equinoxes.  This
    **       results in errors of order 0.1 mas at present.
    **
    **     . The IAU 2000B abridged nutation model (McCarthy & Luzum, 2001)
    **       is used, introducing errors of up to 1 mas.
    **
    **  3) This GAST is compatible with the IAU 2000 resolutions and must be
    **     used only in conjunction with other IAU 2000 compatible
    **     components such as precession-nutation.
    **
    **  4) The result is returned in the range 0 to 2pi.
    **
    **  5) The algorithm is from Capitaine et al. (2003) and IERS
    **     Conventions 2003.
    **
    **  Called:
    **     eraGmst00    Greenwich mean sidereal time, IAU 2000
    **     eraEe00b     equation of the equinoxes, IAU 2000B
    **     eraAnp       normalize angle into range 0 to 2pi
    **
    **  References:
    **
    **     Capitaine, N., Wallace, P.T. and McCarthy, D.D., "Expressions to
    **     implement the IAU 2000 definition of UT1", Astronomy &
    **     Astrophysics, 406, 1135-1149 (2003)
    **
    **     McCarthy, D.D. & Luzum, B.J., "An abridged model of the
    **     precession-nutation of the celestial pole", Celestial Mechanics &
    **     Dynamical Astronomy, 85, 37-49 (2003)
    **
    **     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
    **     IERS Technical Note No. 32, BKG (2004)
    """
    assert ut11.shape[0] == ut12.shape[0]
    cdef unsigned n = ut11.shape[0]
    cdef unsigned int i
    cdef np.ndarray[double, ndim=1] gst = np.empty(n, dtype=np.double)

    for i in range(n):
        gst[i] = eraGst00b(ut11[i], ut12[i]
)
    return gst

@cython.wraparound(False)
@cython.boundscheck(False)
def gst06a(np.ndarray[double, ndim=1] ut11,
           np.ndarray[double, ndim=1] ut12,
           np.ndarray[double, ndim=1] tt1,
           np.ndarray[double, ndim=1] tt2):
    """Wrap double eraGst06a(double uta, double utb, double tta, double ttb)
    **  Greenwich apparent sidereal time (consistent with IAU 2000 and 2006
    **  resolutions).
    **
    **  Given:
    **     uta,utb    double    UT1 as a 2-part Julian Date (Notes 1,2)
    **     tta,ttb    double    TT as a 2-part Julian Date (Notes 1,2)
    **
    **  Returned (function value):
    **                double    Greenwich apparent sidereal time (radians)
    **
    **  Notes:
    **
    **  1) The UT1 and TT dates uta+utb and tta+ttb respectively, are both
    **     Julian Dates, apportioned in any convenient way between the
    **     argument pairs.  For example, JD=2450123.7 could be expressed in
    **     any of these ways, among others:
    **
    **            Part A        Part B
    **
    **         2450123.7           0.0       (JD method)
    **         2451545.0       -1421.3       (J2000 method)
    **         2400000.5       50123.2       (MJD method)
    **         2450123.5           0.2       (date & time method)
    **
    **     The JD method is the most natural and convenient to use in
    **     cases where the loss of several decimal digits of resolution
    **     is acceptable (in the case of UT;  the TT is not at all critical
    **     in this respect).  The J2000 and MJD methods are good compromises
    **     between resolution and convenience.  For UT, the date & time
    **     method is best matched to the algorithm that is used by the Earth
    **     rotation angle function, called internally:  maximum precision is
    **     delivered when the uta argument is for 0hrs UT1 on the day in
    **     question and the utb argument lies in the range 0 to 1, or vice
    **     versa.
    **
    **  2) Both UT1 and TT are required, UT1 to predict the Earth rotation
    **     and TT to predict the effects of precession-nutation.  If UT1 is
    **     used for both purposes, errors of order 100 microarcseconds
    **     result.
    **
    **  3) This GAST is compatible with the IAU 2000/2006 resolutions and
    **     must be used only in conjunction with IAU 2006 precession and
    **     IAU 2000A nutation.
    **
    **  4) The result is returned in the range 0 to 2pi.
    **
    **  Called:
    **     eraPnm06a    classical NPB matrix, IAU 2006/2000A
    **     eraGst06     Greenwich apparent ST, IAU 2006, given NPB matrix
    **
    **  Reference:
    **
    **     Wallace, P.T. & Capitaine, N., 2006, Astron.Astrophys. 459, 981
    """
    assert ut11.shape[0] == ut12.shape[0] == tt1.shape[0] == tt2.shape[0]
    cdef unsigned n = ut11.shape[0]
    cdef unsigned int i
    cdef np.ndarray[double, ndim=1] gst = np.empty(n, dtype=np.double)

    for i in range(n):
        gst[i] = eraGst06a(ut11[i], ut12[i], tt1[i], tt2[i])

    return gst

@cython.wraparound(False)
@cython.boundscheck(False)
def gst94(np.ndarray[double, ndim=1] ut11,
          np.ndarray[double, ndim=1] ut12):
    """Wrap double eraGst94(double uta, double utb)
    **  Greenwich apparent sidereal time (consistent with IAU 1982/94
    **  resolutions).
    **
    **  Given:
    **     uta,utb    double    UT1 as a 2-part Julian Date (Notes 1,2)
    **
    **  Returned (function value):
    **                double    Greenwich apparent sidereal time (radians)
    **
    **  Notes:
    **
    **  1) The UT1 date uta+utb is a Julian Date, apportioned in any
    **     convenient way between the argument pair.  For example,
    **     JD=2450123.7 could be expressed in any of these ways, among
    **     others:
    **
    **             uta            utb
    **
    **         2450123.7           0.0       (JD method)
    **         2451545.0       -1421.3       (J2000 method)
    **         2400000.5       50123.2       (MJD method)
    **         2450123.5           0.2       (date & time method)
    **
    **     The JD method is the most natural and convenient to use in cases
    **     where the loss of several decimal digits of resolution is
    **     acceptable.  The J2000 and MJD methods are good compromises
    **     between resolution and convenience.  For UT, the date & time
    **     method is best matched to the algorithm that is used by the Earth
    **     Rotation Angle function, called internally:  maximum precision is
    **     delivered when the uta argument is for 0hrs UT1 on the day in
    **     question and the utb argument lies in the range 0 to 1, or vice
    **     versa.
    **
    **  2) The result is compatible with the IAU 1982 and 1994 resolutions,
    **     except that accuracy has been compromised for the sake of
    **     convenience in that UT is used instead of TDB (or TT) to compute
    **     the equation of the equinoxes.
    **
    **  3) This GAST must be used only in conjunction with contemporaneous
    **     IAU standards such as 1976 precession, 1980 obliquity and 1982
    **     nutation.  It is not compatible with the IAU 2000 resolutions.
    **
    **  4) The result is returned in the range 0 to 2pi.
    **
    **  Called:
    **     eraGmst82    Greenwich mean sidereal time, IAU 1982
    **     eraEqeq94    equation of the equinoxes, IAU 1994
    **     eraAnp       normalize angle into range 0 to 2pi
    **
    **  References:
    **
    **     Explanatory Supplement to the Astronomical Almanac,
    **     P. Kenneth Seidelmann (ed), University Science Books (1992)
    **
    **     IAU Resolution C7, Recommendation 3 (1994)
    """
    assert ut11.shape[0] == ut12.shape[0]
    cdef unsigned n = ut11.shape[0]
    cdef unsigned int i
    cdef np.ndarray[double, ndim=1] gst = np.empty(n, dtype=np.double)

    for i in range(n):
        gst[i] = eraGst94(ut11[i], ut12[i])

    return gst
