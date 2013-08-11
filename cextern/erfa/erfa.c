#include "erfa.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

void eraA2af(int ndp, double angle, char *sign, int idmsf[4])
/*
**  - - - - - - - -
**   e r a A 2 a f
**  - - - - - - - -
**
**  Decompose radians into degrees, arcminutes, arcseconds, fraction.
**
**  Given:
**     ndp     int     resolution (Note 1)
**     angle   double  angle in radians
**
**  Returned:
**     sign    char    '+' or '-'
**     idmsf   int[4]  degrees, arcminutes, arcseconds, fraction
**
**  Called:
**     eraD2tf      decompose days to hms
**
**  Notes:
**
**  1) The argument ndp is interpreted as follows:
**
**     ndp         resolution
**      :      ...0000 00 00
**     -7         1000 00 00
**     -6          100 00 00
**     -5           10 00 00
**     -4            1 00 00
**     -3            0 10 00
**     -2            0 01 00
**     -1            0 00 10
**      0            0 00 01
**      1            0 00 00.1
**      2            0 00 00.01
**      3            0 00 00.001
**      :            0 00 00.000...
**
**  2) The largest positive useful value for ndp is determined by the
**     size of angle, the format of doubles on the target platform, and
**     the risk of overflowing idmsf[3].  On a typical platform, for
**     angle up to 2pi, the available floating-point precision might
**     correspond to ndp=12.  However, the practical limit is typically
**     ndp=9, set by the capacity of a 32-bit int, or ndp=4 if int is
**     only 16 bits.
**
**  3) The absolute value of angle may exceed 2pi.  In cases where it
**     does not, it is up to the caller to test for and handle the
**     case where angle is very nearly 2pi and rounds up to 360 degrees,
**     by testing for idmsf[0]=360 and setting idmsf[0-3] to zero.
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
/* Hours to degrees * radians to turns */
   const double F = 15.0 / D2PI;


/* Scale then use days to h,m,s function. */
   eraD2tf(ndp, angle*F, sign, idmsf);

   return;

}

void eraA2tf(int ndp, double angle, char *sign, int ihmsf[4])
/*
**  - - - - - - - -
**   e r a A 2 t f
**  - - - - - - - -
**
**  Decompose radians into hours, minutes, seconds, fraction.
**
**  Given:
**     ndp     int     resolution (Note 1)
**     angle   double  angle in radians
**
**  Returned:
**     sign    char    '+' or '-'
**     ihmsf   int[4]  hours, minutes, seconds, fraction
**
**  Called:
**     eraD2tf      decompose days to hms
**
**  Notes:
**
**  1) The argument ndp is interpreted as follows:
**
**     ndp         resolution
**      :      ...0000 00 00
**     -7         1000 00 00
**     -6          100 00 00
**     -5           10 00 00
**     -4            1 00 00
**     -3            0 10 00
**     -2            0 01 00
**     -1            0 00 10
**      0            0 00 01
**      1            0 00 00.1
**      2            0 00 00.01
**      3            0 00 00.001
**      :            0 00 00.000...
**
**  2) The largest positive useful value for ndp is determined by the
**     size of angle, the format of doubles on the target platform, and
**     the risk of overflowing ihmsf[3].  On a typical platform, for
**     angle up to 2pi, the available floating-point precision might
**     correspond to ndp=12.  However, the practical limit is typically
**     ndp=9, set by the capacity of a 32-bit int, or ndp=4 if int is
**     only 16 bits.
**
**  3) The absolute value of angle may exceed 2pi.  In cases where it
**     does not, it is up to the caller to test for and handle the
**     case where angle is very nearly 2pi and rounds up to 24 hours,
**     by testing for ihmsf[0]=24 and setting ihmsf(0-3) to zero.
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
/* Scale then use days to h,m,s function. */
   eraD2tf(ndp, angle/D2PI, sign, ihmsf);

   return;

}

int eraAf2a(char s, int ideg, int iamin, double asec, double *rad)
/*
**  - - - - - - - -
**   e r a A f 2 a
**  - - - - - - - -
**
**  Convert degrees, arcminutes, arcseconds to radians.
**
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
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{

/* Compute the interval. */
   *rad  = ( s == '-' ? -1.0 : 1.0 ) *
           ( 60.0 * ( 60.0 * ( (double) abs(ideg) ) +
                             ( (double) abs(iamin) ) ) +
                                        fabs(asec) ) * DAS2R;

/* Validate arguments and return status. */
   if ( ideg < 0 || ideg > 359 ) return 1;
   if ( iamin < 0 || iamin > 59 ) return 2;
   if ( asec < 0.0 || asec >= 60.0 ) return 3;
   return 0;

}

double eraAnp(double a)
/*
**  - - - - - - -
**   e r a A n p
**  - - - - - - -
**
**  Normalize angle into the range 0 <= a < 2pi.
**
**  Given:
**     a        double     angle (radians)
**
**  Returned (function value):
**              double     angle in range 0-2pi
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double w;


   w = fmod(a, D2PI);
   if (w < 0) w += D2PI;

   return w;

}

double eraAnpm(double a)
/*
**  - - - - - - - -
**   e r a A n p m
**  - - - - - - - -
**
**  Normalize angle into the range -pi <= a < +pi.
**
**  Given:
**     a        double     angle (radians)
**
**  Returned (function value):
**              double     angle in range +/-pi
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double w;


   w = fmod(a, D2PI);
   if (fabs(w) >= DPI) w -= dsign(D2PI, a);

   return w;

}

void eraBi00(double *dpsibi, double *depsbi, double *dra)
/*
**  - - - - - - - -
**   e r a B i 0 0
**  - - - - - - - -
**
**  Frame bias components of IAU 2000 precession-nutation models (part
**  of MHB2000 with additions).
**
**  Returned:
**     dpsibi,depsbi  double  longitude and obliquity corrections
**     dra            double  the ICRS RA of the J2000.0 mean equinox
**
**  Notes:
**
**  1) The frame bias corrections in longitude and obliquity (radians)
**     are required in order to correct for the offset between the GCRS
**     pole and the mean J2000.0 pole.  They define, with respect to the
**     GCRS frame, a J2000.0 mean pole that is consistent with the rest
**     of the IAU 2000A precession-nutation model.
**
**  2) In addition to the displacement of the pole, the complete
**     description of the frame bias requires also an offset in right
**     ascension.  This is not part of the IAU 2000A model, and is from
**     Chapront et al. (2002).  It is returned in radians.
**
**  3) This is a supplemented implementation of one aspect of the IAU
**     2000A nutation model, formally adopted by the IAU General
**     Assembly in 2000, namely MHB2000 (Mathews et al. 2002).
**
**  References:
**
**     Chapront, J., Chapront-Touze, M. & Francou, G., Astron.
**     Astrophys., 387, 700, 2002.
**
**     Mathews, P.M., Herring, T.A., Buffet, B.A., "Modeling of nutation
**     and precession   New nutation series for nonrigid Earth and
**     insights into the Earth's interior", J.Geophys.Res., 107, B4,
**     2002.  The MHB2000 code itself was obtained on 9th September 2002
**     from ftp://maia.usno.navy.mil/conv2000/chapter5/IAU2000A.
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
/* The frame bias corrections in longitude and obliquity */
   const double DPBIAS = -0.041775  * DAS2R,
                DEBIAS = -0.0068192 * DAS2R;

/* The ICRS RA of the J2000.0 equinox (Chapront et al., 2002) */
   const double DRA0 = -0.0146 * DAS2R;


/* Return the results (which are fixed). */
   *dpsibi = DPBIAS;
   *depsbi = DEBIAS;
   *dra = DRA0;

   return;

}

void eraBp00(double date1, double date2,
             double rb[3][3], double rp[3][3], double rbp[3][3])
/*
**  - - - - - - - -
**   e r a B p 0 0
**  - - - - - - - -
**
**  Frame bias and precession, IAU 2000.
**
**  Given:
**     date1,date2  double         TT as a 2-part Julian Date (Note 1)
**
**  Returned:
**     rb           double[3][3]   frame bias matrix (Note 2)
**     rp           double[3][3]   precession matrix (Note 3)
**     rbp          double[3][3]   bias-precession matrix (Note 4)
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
**     convenient way between the two arguments.  For example,
**     JD(TT)=2450123.7 could be expressed in any of these ways,
**     among others:
**
**             date1         date2
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
**  2) The matrix rb transforms vectors from GCRS to mean J2000.0 by
**     applying frame bias.
**
**  3) The matrix rp transforms vectors from J2000.0 mean equator and
**     equinox to mean equator and equinox of date by applying
**     precession.
**
**  4) The matrix rbp transforms vectors from GCRS to mean equator and
**     equinox of date by applying frame bias then precession.  It is
**     the product rp x rb.
**
**  5) It is permissible to re-use the same array in the returned
**     arguments.  The arrays are filled in the order given.
**
**  Called:
**     eraBi00      frame bias components, IAU 2000
**     eraPr00      IAU 2000 precession adjustments
**     eraIr        initialize r-matrix to identity
**     eraRx        rotate around X-axis
**     eraRy        rotate around Y-axis
**     eraRz        rotate around Z-axis
**     eraCr        copy r-matrix
**     eraRxr       product of two r-matrices
**
**  Reference:
**     "Expressions for the Celestial Intermediate Pole and Celestial
**     Ephemeris Origin consistent with the IAU 2000A precession-
**     nutation model", Astron.Astrophys. 400, 1145-1154 (2003)
**
**     n.b. The celestial ephemeris origin (CEO) was renamed "celestial
**          intermediate origin" (CIO) by IAU 2006 Resolution 2.
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
/* J2000.0 obliquity (Lieske et al. 1977) */
   const double EPS0 = 84381.448 * DAS2R;

   double t, dpsibi, depsbi;
   double dra0, psia77, oma77, chia, dpsipr, depspr, psia, oma,
          rbw[3][3];


/* Interval between fundamental epoch J2000.0 and current date (JC). */
   t = ((date1 - DJ00) + date2) / DJC;

/* Frame bias. */
   eraBi00(&dpsibi, &depsbi, &dra0);

/* Precession angles (Lieske et al. 1977) */
   psia77 = (5038.7784 + (-1.07259 + (-0.001147) * t) * t) * t * DAS2R;
   oma77  =       EPS0 + ((0.05127 + (-0.007726) * t) * t) * t * DAS2R;
   chia   = (  10.5526 + (-2.38064 + (-0.001125) * t) * t) * t * DAS2R;

/* Apply IAU 2000 precession corrections. */
   eraPr00(date1, date2, &dpsipr,  &depspr);
   psia = psia77 + dpsipr;
   oma  = oma77  + depspr;

/* Frame bias matrix: GCRS to J2000.0. */
   eraIr(rbw);
   eraRz(dra0, rbw);
   eraRy(dpsibi * sin(EPS0), rbw);
   eraRx(-depsbi, rbw);
   eraCr(rbw, rb);

/* Precession matrix: J2000.0 to mean of date. */
   eraIr(rp);
   eraRx(EPS0,  rp);
   eraRz(-psia, rp);
   eraRx(-oma,  rp);
   eraRz(chia,  rp);

/* Bias-precession matrix: GCRS to mean of date. */
   eraRxr(rp, rbw, rbp);

   return;

}

void eraBp06(double date1, double date2,
             double rb[3][3], double rp[3][3], double rbp[3][3])
/*
**  - - - - - - - -
**   e r a B p 0 6
**  - - - - - - - -
**
**  Frame bias and precession, IAU 2006.
**
**  Given:
**     date1,date2  double         TT as a 2-part Julian Date (Note 1)
**
**  Returned:
**     rb           double[3][3]   frame bias matrix (Note 2)
**     rp           double[3][3]   precession matrix (Note 3)
**     rbp          double[3][3]   bias-precession matrix (Note 4)
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
**     convenient way between the two arguments.  For example,
**     JD(TT)=2450123.7 could be expressed in any of these ways,
**     among others:
**
**             date1         date2
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
**  2) The matrix rb transforms vectors from GCRS to mean J2000.0 by
**     applying frame bias.
**
**  3) The matrix rp transforms vectors from mean J2000.0 to mean of
**     date by applying precession.
**
**  4) The matrix rbp transforms vectors from GCRS to mean of date by
**     applying frame bias then precession.  It is the product rp x rb.
**
**  Called:
**     eraPfw06     bias-precession F-W angles, IAU 2006
**     eraFw2m      F-W angles to r-matrix
**     eraPmat06    PB matrix, IAU 2006
**     eraTr        transpose r-matrix
**     eraRxr       product of two r-matrices
**
**  References:
**
**     Capitaine, N. & Wallace, P.T., 2006, Astron.Astrophys. 450, 855
**
**     Wallace, P.T. & Capitaine, N., 2006, Astron.Astrophys. 459, 981
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double gamb, phib, psib, epsa, rbt[3][3];


/* B matrix. */
   eraPfw06(DJM0, DJM00, &gamb, &phib, &psib, &epsa);
   eraFw2m(gamb, phib, psib, epsa, rb);

/* PxB matrix. */
   eraPmat06(date1, date2, rbp);

/* P matrix. */
   eraTr(rb, rbt);
   eraRxr(rbp, rbt, rp);

   return;

}

void eraBpn2xy(double rbpn[3][3], double *x, double *y)
/*
**  - - - - - - - - - -
**   e r a B p n 2 x y
**  - - - - - - - - - -
**
**  Extract from the bias-precession-nutation matrix the X,Y coordinates
**  of the Celestial Intermediate Pole.
**
**  Given:
**     rbpn      double[3][3]  celestial-to-true matrix (Note 1)
**
**  Returned:
**     x,y       double        Celestial Intermediate Pole (Note 2)
**
**  Notes:
**
**  1) The matrix rbpn transforms vectors from GCRS to true equator (and
**     CIO or equinox) of date, and therefore the Celestial Intermediate
**     Pole unit vector is the bottom row of the matrix.
**
**  2) The arguments x,y are components of the Celestial Intermediate
**     Pole unit vector in the Geocentric Celestial Reference System.
**
**  Reference:
**
**     "Expressions for the Celestial Intermediate Pole and Celestial
**     Ephemeris Origin consistent with the IAU 2000A precession-
**     nutation model", Astron.Astrophys. 400, 1145-1154
**     (2003)
**
**     n.b. The celestial ephemeris origin (CEO) was renamed "celestial
**          intermediate origin" (CIO) by IAU 2006 Resolution 2.
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
/* Extract the X,Y coordinates. */
   *x = rbpn[2][0];
   *y = rbpn[2][1];

   return;

}

void eraC2i00a(double date1, double date2, double rc2i[3][3])
/*
**  - - - - - - - - - -
**   e r a C 2 i 0 0 a
**  - - - - - - - - - -
**
**  Form the celestial-to-intermediate matrix for a given date using the
**  IAU 2000A precession-nutation model.
**
**  Given:
**     date1,date2 double       TT as a 2-part Julian Date (Note 1)
**
**  Returned:
**     rc2i        double[3][3] celestial-to-intermediate matrix (Note 2)
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
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
**  2) The matrix rc2i is the first stage in the transformation from
**     celestial to terrestrial coordinates:
**
**        [TRS]  =  RPOM * R_3(ERA) * rc2i * [CRS]
**
**               =  rc2t * [CRS]
**
**     where [CRS] is a vector in the Geocentric Celestial Reference
**     System and [TRS] is a vector in the International Terrestrial
**     Reference System (see IERS Conventions 2003), ERA is the Earth
**     Rotation Angle and RPOM is the polar motion matrix.
**
**  3) A faster, but slightly less accurate result (about 1 mas), can be
**     obtained by using instead the eraC2i00b function.
**
**  Called:
**     eraPnm00a    classical NPB matrix, IAU 2000A
**     eraC2ibpn    celestial-to-intermediate matrix, given NPB matrix
**
**  References:
**
**     "Expressions for the Celestial Intermediate Pole and Celestial
**     Ephemeris Origin consistent with the IAU 2000A precession-
**     nutation model", Astron.Astrophys. 400, 1145-1154
**     (2003)
**
**     n.b. The celestial ephemeris origin (CEO) was renamed "celestial
**          intermediate origin" (CIO) by IAU 2006 Resolution 2.
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double rbpn[3][3];


/* Obtain the celestial-to-true matrix (IAU 2000A). */
   eraPnm00a(date1, date2, rbpn);

/* Form the celestial-to-intermediate matrix. */
   eraC2ibpn(date1, date2, rbpn, rc2i);

   return;

}

void eraC2i00b(double date1, double date2, double rc2i[3][3])
/*
**  - - - - - - - - - -
**   e r a C 2 i 0 0 b
**  - - - - - - - - - -
**
**  Form the celestial-to-intermediate matrix for a given date using the
**  IAU 2000B precession-nutation model.
**
**  Given:
**     date1,date2 double       TT as a 2-part Julian Date (Note 1)
**
**  Returned:
**     rc2i        double[3][3] celestial-to-intermediate matrix (Note 2)
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
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
**  2) The matrix rc2i is the first stage in the transformation from
**     celestial to terrestrial coordinates:
**
**        [TRS]  =  RPOM * R_3(ERA) * rc2i * [CRS]
**
**               =  rc2t * [CRS]
**
**     where [CRS] is a vector in the Geocentric Celestial Reference
**     System and [TRS] is a vector in the International Terrestrial
**     Reference System (see IERS Conventions 2003), ERA is the Earth
**     Rotation Angle and RPOM is the polar motion matrix.
**
**  3) The present function is faster, but slightly less accurate (about
**     1 mas), than the eraC2i00a function.
**
**  Called:
**     eraPnm00b    classical NPB matrix, IAU 2000B
**     eraC2ibpn    celestial-to-intermediate matrix, given NPB matrix
**
**  References:
**
**     "Expressions for the Celestial Intermediate Pole and Celestial
**     Ephemeris Origin consistent with the IAU 2000A precession-
**     nutation model", Astron.Astrophys. 400, 1145-1154
**     (2003)
**
**     n.b. The celestial ephemeris origin (CEO) was renamed "celestial
**          intermediate origin" (CIO) by IAU 2006 Resolution 2.
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double rbpn[3][3];


/* Obtain the celestial-to-true matrix (IAU 2000B). */
   eraPnm00b(date1, date2, rbpn);

/* Form the celestial-to-intermediate matrix. */
   eraC2ibpn(date1, date2, rbpn, rc2i);

   return;

}

void eraC2i06a(double date1, double date2, double rc2i[3][3])
/*
**  - - - - - - - - - -
**   e r a _ c 2 i 0 6 a
**  - - - - - - - - - -
**
**  Form the celestial-to-intermediate matrix for a given date using the
**  IAU 2006 precession and IAU 2000A nutation models.
**
**  Given:
**     date1,date2 double       TT as a 2-part Julian Date (Note 1)
**
**  Returned:
**     rc2i        double[3][3] celestial-to-intermediate matrix (Note 2)
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
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
**  2) The matrix rc2i is the first stage in the transformation from
**     celestial to terrestrial coordinates:
**
**        [TRS]  =  RPOM * R_3(ERA) * rc2i * [CRS]
**
**               =  RC2T * [CRS]
**
**     where [CRS] is a vector in the Geocentric Celestial Reference
**     System and [TRS] is a vector in the International Terrestrial
**     Reference System (see IERS Conventions 2003), ERA is the Earth
**     Rotation Angle and RPOM is the polar motion matrix.
**
**  Called:
**     eraPnm06a    classical NPB matrix, IAU 2006/2000A
**     eraBpn2xy    extract CIP X,Y coordinates from NPB matrix
**     eraS06       the CIO locator s, Given X,Y, IAU 2006
**     eraC2ixys    celestial-to-intermediate matrix, Given X,Y and s
**
**  References:
**
**     McCarthy, D. D., Petit, G. (eds.), 2004, IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double rbpn[3][3], x, y, s;


/* Obtain the celestial-to-true matrix (IAU 2006/2000A). */
   eraPnm06a(date1, date2, rbpn);

/* Extract the X,Y coordinates. */
   eraBpn2xy(rbpn, &x, &y);

/* Obtain the CIO locator. */
   s = eraS06(date1, date2, x, y);

/* Form the celestial-to-intermediate matrix. */
   eraC2ixys(x, y, s, rc2i);

   return;

}

void eraC2ibpn(double date1, double date2, double rbpn[3][3],
               double rc2i[3][3])
/*
**  - - - - - - - - - -
**   e r a C 2 i b p n
**  - - - - - - - - - -
**
**  Form the celestial-to-intermediate matrix for a given date given
**  the bias-precession-nutation matrix.  IAU 2000.
**
**  Given:
**     date1,date2 double       TT as a 2-part Julian Date (Note 1)
**     rbpn        double[3][3] celestial-to-true matrix (Note 2)
**
**  Returned:
**     rc2i        double[3][3] celestial-to-intermediate matrix (Note 3)
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
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
**  2) The matrix rbpn transforms vectors from GCRS to true equator (and
**     CIO or equinox) of date.  Only the CIP (bottom row) is used.
**
**  3) The matrix rc2i is the first stage in the transformation from
**     celestial to terrestrial coordinates:
**
**        [TRS] = RPOM * R_3(ERA) * rc2i * [CRS]
**
**              = RC2T * [CRS]
**
**     where [CRS] is a vector in the Geocentric Celestial Reference
**     System and [TRS] is a vector in the International Terrestrial
**     Reference System (see IERS Conventions 2003), ERA is the Earth
**     Rotation Angle and RPOM is the polar motion matrix.
**
**  4) Although its name does not include "00", This function is in fact
**     specific to the IAU 2000 models.
**
**  Called:
**     eraBpn2xy    extract CIP X,Y coordinates from NPB matrix
**     eraC2ixy     celestial-to-intermediate matrix, given X,Y
**
**  References:
**     "Expressions for the Celestial Intermediate Pole and Celestial
**     Ephemeris Origin consistent with the IAU 2000A precession-
**     nutation model", Astron.Astrophys. 400, 1145-1154 (2003)
**
**     n.b. The celestial ephemeris origin (CEO) was renamed "celestial
**          intermediate origin" (CIO) by IAU 2006 Resolution 2.
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double x, y;


/* Extract the X,Y coordinates. */
   eraBpn2xy(rbpn, &x, &y);

/* Form the celestial-to-intermediate matrix (n.b. IAU 2000 specific). */
   eraC2ixy(date1, date2, x, y, rc2i);

   return;

}

void eraC2ixy(double date1, double date2, double x, double y,
              double rc2i[3][3])
/*
**  - - - - - - - - -
**   e r a C 2 i x y
**  - - - - - - - - -
**
**  Form the celestial to intermediate-frame-of-date matrix for a given
**  date when the CIP X,Y coordinates are known.  IAU 2000.
**
**  Given:
**     date1,date2 double       TT as a 2-part Julian Date (Note 1)
**     x,y         double       Celestial Intermediate Pole (Note 2)
**
**  Returned:
**     rc2i        double[3][3] celestial-to-intermediate matrix (Note 3)
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
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
**  2) The Celestial Intermediate Pole coordinates are the x,y components
**     of the unit vector in the Geocentric Celestial Reference System.
**
**  3) The matrix rc2i is the first stage in the transformation from
**     celestial to terrestrial coordinates:
**
**        [TRS] = RPOM * R_3(ERA) * rc2i * [CRS]
**
**              = RC2T * [CRS]
**
**     where [CRS] is a vector in the Geocentric Celestial Reference
**     System and [TRS] is a vector in the International Terrestrial
**     Reference System (see IERS Conventions 2003), ERA is the Earth
**     Rotation Angle and RPOM is the polar motion matrix.
**
**  4) Although its name does not include "00", This function is in fact
**     specific to the IAU 2000 models.
**
**  Called:
**     eraC2ixys    celestial-to-intermediate matrix, given X,Y and s
**     eraS00       the CIO locator s, given X,Y, IAU 2000A
**
**  Reference:
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/

{
/* Compute s and then the matrix. */
   eraC2ixys(x, y, eraS00(date1, date2, x, y), rc2i);

   return;

}

void eraC2ixys(double x, double y, double s, double rc2i[3][3])
/*
**  - - - - - - - - - -
**   e r a C 2 i x y s
**  - - - - - - - - - -
**
**  Form the celestial to intermediate-frame-of-date matrix given the CIP
**  X,Y and the CIO locator s.
**
**  Given:
**     x,y      double         Celestial Intermediate Pole (Note 1)
**     s        double         the CIO locator s (Note 2)
**
**  Returned:
**     rc2i     double[3][3]   celestial-to-intermediate matrix (Note 3)
**
**  Notes:
**
**  1) The Celestial Intermediate Pole coordinates are the x,y
**     components of the unit vector in the Geocentric Celestial
**     Reference System.
**
**  2) The CIO locator s (in radians) positions the Celestial
**     Intermediate Origin on the equator of the CIP.
**
**  3) The matrix rc2i is the first stage in the transformation from
**     celestial to terrestrial coordinates:
**
**        [TRS] = RPOM * R_3(ERA) * rc2i * [CRS]
**
**              = RC2T * [CRS]
**
**     where [CRS] is a vector in the Geocentric Celestial Reference
**     System and [TRS] is a vector in the International Terrestrial
**     Reference System (see IERS Conventions 2003), ERA is the Earth
**     Rotation Angle and RPOM is the polar motion matrix.
**
**  Called:
**     eraIr        initialize r-matrix to identity
**     eraRz        rotate around Z-axis
**     eraRy        rotate around Y-axis
**
**  Reference:
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double r2, e, d;


/* Obtain the spherical angles E and d. */
   r2 = x*x + y*y;
   e = (r2 != 0.0) ? atan2(y, x) : 0.0;
   d = atan(sqrt(r2 / (1.0 - r2)));

/* Form the matrix. */
   eraIr(rc2i);
   eraRz(e, rc2i);
   eraRy(d, rc2i);
   eraRz(-(e+s), rc2i);

   return;

}

void eraC2s(double p[3], double *theta, double *phi)
/*
**  - - - - - - -
**   e r a C 2 s
**  - - - - - - -
**
**  P-vector to spherical coordinates.
**
**  Given:
**     p      double[3]    p-vector
**
**  Returned:
**     theta  double       longitude angle (radians)
**     phi    double       latitude angle (radians)
**
**  Notes:
**
**  1) The vector p can have any magnitude; only its direction is used.
**
**  2) If p is null, zero theta and phi are returned.
**
**  3) At either pole, zero theta is returned.
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double x, y, z, d2;


   x  = p[0];
   y  = p[1];
   z  = p[2];
   d2 = x*x + y*y;

   *theta = (d2 == 0.0) ? 0.0 : atan2(y, x);
   *phi = (z == 0.0) ? 0.0 : atan2(z, sqrt(d2));

   return;

}

void eraC2t00a(double tta, double ttb, double uta, double utb,
               double xp, double yp, double rc2t[3][3])
/*
**  - - - - - - - - - -
**   e r a C 2 t 0 0 a
**  - - - - - - - - - -
**
**  Form the celestial to terrestrial matrix given the date, the UT1 and
**  the polar motion, using the IAU 2000A nutation model.
**
**  Given:
**     tta,ttb  double         TT as a 2-part Julian Date (Note 1)
**     uta,utb  double         UT1 as a 2-part Julian Date (Note 1)
**     xp,yp    double         coordinates of the pole (radians, Note 2)
**
**  Returned:
**     rc2t     double[3][3]   celestial-to-terrestrial matrix (Note 3)
**
**  Notes:
**
**  1) The TT and UT1 dates tta+ttb and uta+utb are Julian Dates,
**     apportioned in any convenient way between the arguments uta and
**     utb.  For example, JD(UT1)=2450123.7 could be expressed in any of
**     these ways, among others:
**
**             uta            utb
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in
**     cases where the loss of several decimal digits of resolution is
**     acceptable.  The J2000 and MJD methods are good compromises
**     between resolution and convenience.  In the case of uta,utb, the
**     date & time method is best matched to the Earth rotation angle
**     algorithm used:  maximum precision is delivered when the uta
**     argument is for 0hrs UT1 on the day in question and the utb
**     argument lies in the range 0 to 1, or vice versa.
**
**  2) The arguments xp and yp are the coordinates (in radians) of the
**     Celestial Intermediate Pole with respect to the International
**     Terrestrial Reference System (see IERS Conventions 2003),
**     measured along the meridians to 0 and 90 deg west respectively.
**
**  3) The matrix rc2t transforms from celestial to terrestrial
**     coordinates:
**
**        [TRS] = RPOM * R_3(ERA) * RC2I * [CRS]
**
**              = rc2t * [CRS]
**
**     where [CRS] is a vector in the Geocentric Celestial Reference
**     System and [TRS] is a vector in the International Terrestrial
**     Reference System (see IERS Conventions 2003), RC2I is the
**     celestial-to-intermediate matrix, ERA is the Earth rotation
**     angle and RPOM is the polar motion matrix.
**
**  4) A faster, but slightly less accurate result (about 1 mas), can
**     be obtained by using instead the eraC2t00b function.
**
**  Called:
**     eraC2i00a    celestial-to-intermediate matrix, IAU 2000A
**     eraEra00     Earth rotation angle, IAU 2000
**     eraSp00      the TIO locator s', IERS 2000
**     eraPom00     polar motion matrix
**     eraC2tcio    form CIO-based celestial-to-terrestrial matrix
**
**  Reference:
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double rc2i[3][3], era, sp, rpom[3][3];


/* Form the celestial-to-intermediate matrix for this TT (IAU 2000A). */
   eraC2i00a(tta, ttb, rc2i );

/* Predict the Earth rotation angle for this UT1. */
   era = eraEra00(uta, utb);

/* Estimate s'. */
   sp = eraSp00(tta, ttb);

/* Form the polar motion matrix. */
   eraPom00(xp, yp, sp, rpom);

/* Combine to form the celestial-to-terrestrial matrix. */
   eraC2tcio(rc2i, era, rpom, rc2t);

   return;

}

void eraC2t00b(double tta, double ttb, double uta, double utb,
               double xp, double yp, double rc2t[3][3])
/*
**  - - - - - - - - - -
**   e r a C 2 t 0 0 b
**  - - - - - - - - - -
**
**  Form the celestial to terrestrial matrix given the date, the UT1 and
**  the polar motion, using the IAU 2000B nutation model.
**
**  Given:
**     tta,ttb  double         TT as a 2-part Julian Date (Note 1)
**     uta,utb  double         UT1 as a 2-part Julian Date (Note 1)
**     xp,yp    double         coordinates of the pole (radians, Note 2)
**
**  Returned:
**     rc2t     double[3][3]   celestial-to-terrestrial matrix (Note 3)
**
**  Notes:
**
**  1) The TT and UT1 dates tta+ttb and uta+utb are Julian Dates,
**     apportioned in any convenient way between the arguments uta and
**     utb.  For example, JD(UT1)=2450123.7 could be expressed in any of
**     these ways, among others:
**
**             uta            utb
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in
**     cases where the loss of several decimal digits of resolution is
**     acceptable.  The J2000 and MJD methods are good compromises
**     between resolution and convenience.  In the case of uta,utb, the
**     date & time method is best matched to the Earth rotation angle
**     algorithm used:  maximum precision is delivered when the uta
**     argument is for 0hrs UT1 on the day in question and the utb
**     argument lies in the range 0 to 1, or vice versa.
**
**  2) The arguments xp and yp are the coordinates (in radians) of the
**     Celestial Intermediate Pole with respect to the International
**     Terrestrial Reference System (see IERS Conventions 2003),
**     measured along the meridians to 0 and 90 deg west respectively.
**
**  3) The matrix rc2t transforms from celestial to terrestrial
**     coordinates:
**
**        [TRS] = RPOM * R_3(ERA) * RC2I * [CRS]
**
**              = rc2t * [CRS]
**
**     where [CRS] is a vector in the Geocentric Celestial Reference
**     System and [TRS] is a vector in the International Terrestrial
**     Reference System (see IERS Conventions 2003), RC2I is the
**     celestial-to-intermediate matrix, ERA is the Earth rotation
**     angle and RPOM is the polar motion matrix.
**
**  4) The present function is faster, but slightly less accurate (about
**     1 mas), than the eraC2t00a function.
**
**  Called:
**     eraC2i00b    celestial-to-intermediate matrix, IAU 2000B
**     eraEra00     Earth rotation angle, IAU 2000
**     eraPom00     polar motion matrix
**     eraC2tcio    form CIO-based celestial-to-terrestrial matrix
**
**  Reference:
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double rc2i[3][3], era, rpom[3][3];


/* Form the celestial-to-intermediate matrix for this TT (IAU 2000B). */
   eraC2i00b(tta, ttb, rc2i);

/* Predict the Earth rotation angle for this UT1. */
   era = eraEra00(uta, utb);

/* Form the polar motion matrix (neglecting s'). */
   eraPom00(xp, yp, 0.0, rpom);

/* Combine to form the celestial-to-terrestrial matrix. */
   eraC2tcio(rc2i, era, rpom, rc2t);

   return;

}

void eraC2t06a(double tta, double ttb, double uta, double utb,
               double xp, double yp, double rc2t[3][3])
/*
**  - - - - - - - - - -
**   e r a C 2 t 0 6 a
**  - - - - - - - - - -
**
**  Form the celestial to terrestrial matrix given the date, the UT1 and
**  the polar motion, using the IAU 2006 precession and IAU 2000A
**  nutation models.
**
**  Given:
**     tta,ttb  double         TT as a 2-part Julian Date (Note 1)
**     uta,utb  double         UT1 as a 2-part Julian Date (Note 1)
**     xp,yp    double         coordinates of the pole (radians, Note 2)
**
**  Returned:
**     rc2t     double[3][3]   celestial-to-terrestrial matrix (Note 3)
**
**  Notes:
**
**  1) The TT and UT1 dates tta+ttb and uta+utb are Julian Dates,
**     apportioned in any convenient way between the arguments uta and
**     utb.  For example, JD(UT1)=2450123.7 could be expressed in any of
**     these ways, among others:
**
**             uta            utb
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in
**     cases where the loss of several decimal digits of resolution is
**     acceptable.  The J2000 and MJD methods are good compromises
**     between resolution and convenience.  In the case of uta,utb, the
**     date & time method is best matched to the Earth rotation angle
**     algorithm used:  maximum precision is delivered when the uta
**     argument is for 0hrs UT1 on the day in question and the utb
**     argument lies in the range 0 to 1, or vice versa.
**
**  2) The arguments xp and yp are the coordinates (in radians) of the
**     Celestial Intermediate Pole with respect to the International
**     Terrestrial Reference System (see IERS Conventions 2003),
**     measured along the meridians to 0 and 90 deg west respectively.
**
**  3) The matrix rc2t transforms from celestial to terrestrial
**     coordinates:
**
**        [TRS] = RPOM * R_3(ERA) * RC2I * [CRS]
**
**              = rc2t * [CRS]
**
**     where [CRS] is a vector in the Geocentric Celestial Reference
**     System and [TRS] is a vector in the International Terrestrial
**     Reference System (see IERS Conventions 2003), RC2I is the
**     celestial-to-intermediate matrix, ERA is the Earth rotation
**     angle and RPOM is the polar motion matrix.
**
**  Called:
**     eraC2i06a    celestial-to-intermediate matrix, IAU 2006/2000A
**     eraEra00     Earth rotation angle, IAU 2000
**     eraSp00      the TIO locator s', IERS 2000
**     eraPom00     polar motion matrix
**     eraC2tcio    form CIO-based celestial-to-terrestrial matrix
**
**  Reference:
**
**     McCarthy, D. D., Petit, G. (eds.), 2004, IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double rc2i[3][3], era, sp, rpom[3][3];


/* Form the celestial-to-intermediate matrix for this TT. */
   eraC2i06a(tta, ttb, rc2i);

/* Predict the Earth rotation angle for this UT1. */
   era = eraEra00(uta, utb);

/* Estimate s'. */
   sp = eraSp00(tta, ttb);

/* Form the polar motion matrix. */
   eraPom00(xp, yp, sp, rpom);

/* Combine to form the celestial-to-terrestrial matrix. */
   eraC2tcio(rc2i, era, rpom, rc2t);

   return;

}

void eraC2tcio(double rc2i[3][3], double era, double rpom[3][3],
               double rc2t[3][3])
/*
**  - - - - - - - - - -
**   e r a C 2 t c i o
**  - - - - - - - - - -
**
**  Assemble the celestial to terrestrial matrix from CIO-based
**  components (the celestial-to-intermediate matrix, the Earth Rotation
**  Angle and the polar motion matrix).
**
**  Given:
**     rc2i     double[3][3]    celestial-to-intermediate matrix
**     era      double          Earth rotation angle
**     rpom     double[3][3]    polar-motion matrix
**
**  Returned:
**     rc2t     double[3][3]    celestial-to-terrestrial matrix
**
**  Notes:
**
**  1) This function constructs the rotation matrix that transforms
**     vectors in the celestial system into vectors in the terrestrial
**     system.  It does so starting from precomputed components, namely
**     the matrix which rotates from celestial coordinates to the
**     intermediate frame, the Earth rotation angle and the polar motion
**     matrix.  One use of the present function is when generating a
**     series of celestial-to-terrestrial matrices where only the Earth
**     Rotation Angle changes, avoiding the considerable overhead of
**     recomputing the precession-nutation more often than necessary to
**     achieve given accuracy objectives.
**
**  2) The relationship between the arguments is as follows:
**
**        [TRS] = RPOM * R_3(ERA) * rc2i * [CRS]
**
**              = rc2t * [CRS]
**
**     where [CRS] is a vector in the Geocentric Celestial Reference
**     System and [TRS] is a vector in the International Terrestrial
**     Reference System (see IERS Conventions 2003).
**
**  Called:
**     eraCr        copy r-matrix
**     eraRz        rotate around Z-axis
**     eraRxr       product of two r-matrices
**
**  Reference:
**
**     McCarthy, D. D., Petit, G. (eds.), 2004, IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double r[3][3];


/* Construct the matrix. */
   eraCr(rc2i, r);
   eraRz(era, r);
   eraRxr(rpom, r, rc2t);

   return;

}

void eraC2teqx(double rbpn[3][3], double gst, double rpom[3][3],
               double rc2t[3][3])
/*
**  - - - - - - - - - -
**   e r a C 2 t e q x
**  - - - - - - - - - -
**
**  Assemble the celestial to terrestrial matrix from equinox-based
**  components (the celestial-to-true matrix, the Greenwich Apparent
**  Sidereal Time and the polar motion matrix).
**
**  Given:
**     rbpn     double[3][3]    celestial-to-true matrix
**     gst      double          Greenwich (apparent) Sidereal Time
**     rpom     double[3][3]    polar-motion matrix
**
**  Returned:
**     rc2t     double[3][3]    celestial-to-terrestrial matrix (Note 2)
**
**  Notes:
**
**  1) This function constructs the rotation matrix that transforms
**     vectors in the celestial system into vectors in the terrestrial
**     system.  It does so starting from precomputed components, namely
**     the matrix which rotates from celestial coordinates to the
**     true equator and equinox of date, the Greenwich Apparent Sidereal
**     Time and the polar motion matrix.  One use of the present function
**     is when generating a series of celestial-to-terrestrial matrices
**     where only the Sidereal Time changes, avoiding the considerable
**     overhead of recomputing the precession-nutation more often than
**     necessary to achieve given accuracy objectives.
**
**  2) The relationship between the arguments is as follows:
**
**        [TRS] = rpom * R_3(gst) * rbpn * [CRS]
**
**              = rc2t * [CRS]
**
**     where [CRS] is a vector in the Geocentric Celestial Reference
**     System and [TRS] is a vector in the International Terrestrial
**     Reference System (see IERS Conventions 2003).
**
**  Called:
**     eraCr        copy r-matrix
**     eraRz        rotate around Z-axis
**     eraRxr       product of two r-matrices
**
**  Reference:
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double r[3][3];


/* Construct the matrix. */
   eraCr(rbpn, r);
   eraRz(gst, r);
   eraRxr(rpom, r, rc2t);

   return;

}

void eraC2tpe(double tta, double ttb, double uta, double utb,
              double dpsi, double deps, double xp, double yp,
              double rc2t[3][3])
/*
**  - - - - - - - - -
**   e r a C 2 t p e
**  - - - - - - - - -
**
**  Form the celestial to terrestrial matrix given the date, the UT1,
**  the nutation and the polar motion.  IAU 2000.
**
**  Given:
**     tta,ttb    double        TT as a 2-part Julian Date (Note 1)
**     uta,utb    double        UT1 as a 2-part Julian Date (Note 1)
**     dpsi,deps  double        nutation (Note 2)
**     xp,yp      double        coordinates of the pole (radians, Note 3)
**
**  Returned:
**     rc2t       double[3][3]  celestial-to-terrestrial matrix (Note 4)
**
**  Notes:
**
**  1) The TT and UT1 dates tta+ttb and uta+utb are Julian Dates,
**     apportioned in any convenient way between the arguments uta and
**     utb.  For example, JD(UT1)=2450123.7 could be expressed in any of
**     these ways, among others:
**
**             uta            utb
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in
**     cases where the loss of several decimal digits of resolution is
**     acceptable.  The J2000 and MJD methods are good compromises
**     between resolution and convenience.  In the case of uta,utb, the
**     date & time method is best matched to the Earth rotation angle
**     algorithm used:  maximum precision is delivered when the uta
**     argument is for 0hrs UT1 on the day in question and the utb
**     argument lies in the range 0 to 1, or vice versa.
**
**  2) The caller is responsible for providing the nutation components;
**     they are in longitude and obliquity, in radians and are with
**     respect to the equinox and ecliptic of date.  For high-accuracy
**     applications, free core nutation should be included as well as
**     any other relevant corrections to the position of the CIP.
**
**  3) The arguments xp and yp are the coordinates (in radians) of the
**     Celestial Intermediate Pole with respect to the International
**     Terrestrial Reference System (see IERS Conventions 2003),
**     measured along the meridians to 0 and 90 deg west respectively.
**
**  4) The matrix rc2t transforms from celestial to terrestrial
**     coordinates:
**
**        [TRS] = RPOM * R_3(GST) * RBPN * [CRS]
**
**              = rc2t * [CRS]
**
**     where [CRS] is a vector in the Geocentric Celestial Reference
**     System and [TRS] is a vector in the International Terrestrial
**     Reference System (see IERS Conventions 2003), RBPN is the
**     bias-precession-nutation matrix, GST is the Greenwich (apparent)
**     Sidereal Time and RPOM is the polar motion matrix.
**
**  5) Although its name does not include "00", This function is in fact
**     specific to the IAU 2000 models.
**
**  Called:
**     eraPn00      bias/precession/nutation results, IAU 2000
**     eraGmst00    Greenwich mean sidereal time, IAU 2000
**     eraSp00      the TIO locator s', IERS 2000
**     eraEe00      equation of the equinoxes, IAU 2000
**     eraPom00     polar motion matrix
**     eraC2teqx    form equinox-based celestial-to-terrestrial matrix
**
**  Reference:
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double epsa, rb[3][3], rp[3][3], rbp[3][3], rn[3][3],
          rbpn[3][3], gmst, ee, sp, rpom[3][3];


/* Form the celestial-to-true matrix for this TT. */
   eraPn00(tta, ttb, dpsi, deps, &epsa, rb, rp, rbp, rn, rbpn);

/* Predict the Greenwich Mean Sidereal Time for this UT1 and TT. */
   gmst = eraGmst00(uta, utb, tta, ttb);

/* Predict the equation of the equinoxes given TT and nutation. */
   ee = eraEe00(tta, ttb, epsa, dpsi);

/* Estimate s'. */
   sp = eraSp00(tta, ttb);

/* Form the polar motion matrix. */
   eraPom00(xp, yp, sp, rpom);

/* Combine to form the celestial-to-terrestrial matrix. */
   eraC2teqx(rbpn, gmst + ee, rpom, rc2t);

   return;

}

void eraC2txy(double tta, double ttb, double uta, double utb,
              double x, double y, double xp, double yp,
              double rc2t[3][3])
/*
**  - - - - - - - - -
**   e r a C 2 t x y
**  - - - - - - - - -
**
**  Form the celestial to terrestrial matrix given the date, the UT1,
**  the CIP coordinates and the polar motion.  IAU 2000.
**
**  Given:
**     tta,ttb  double         TT as a 2-part Julian Date (Note 1)
**     uta,utb  double         UT1 as a 2-part Julian Date (Note 1)
**     x,y      double         Celestial Intermediate Pole (Note 2)
**     xp,yp    double         coordinates of the pole (radians, Note 3)
**
**  Returned:
**     rc2t     double[3][3]   celestial-to-terrestrial matrix (Note 4)
**
**  Notes:
**
**  1) The TT and UT1 dates tta+ttb and uta+utb are Julian Dates,
**     apportioned in any convenient way between the arguments uta and
**     utb.  For example, JD(UT1)=2450123.7 could be expressed in any o
**     these ways, among others:
**
**             uta            utb
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in
**     cases where the loss of several decimal digits of resolution is
**     acceptable.  The J2000 and MJD methods are good compromises
**     between resolution and convenience.  In the case of uta,utb, the
**     date & time method is best matched to the Earth rotation angle
**     algorithm used:  maximum precision is delivered when the uta
**     argument is for 0hrs UT1 on the day in question and the utb
**     argument lies in the range 0 to 1, or vice versa.
**
**  2) The Celestial Intermediate Pole coordinates are the x,y
**     components of the unit vector in the Geocentric Celestial
**     Reference System.
**
**  3) The arguments xp and yp are the coordinates (in radians) of the
**     Celestial Intermediate Pole with respect to the International
**     Terrestrial Reference System (see IERS Conventions 2003),
**     measured along the meridians to 0 and 90 deg west respectively.
**
**  4) The matrix rc2t transforms from celestial to terrestrial
**     coordinates:
**
**        [TRS] = RPOM * R_3(ERA) * RC2I * [CRS]
**
**              = rc2t * [CRS]
**
**     where [CRS] is a vector in the Geocentric Celestial Reference
**     System and [TRS] is a vector in the International Terrestrial
**     Reference System (see IERS Conventions 2003), ERA is the Earth
**     Rotation Angle and RPOM is the polar motion matrix.
**
**  5) Although its name does not include "00", This function is in fact
**     specific to the IAU 2000 models.
**
**  Called:
**     eraC2ixy     celestial-to-intermediate matrix, given X,Y
**     eraEra00     Earth rotation angle, IAU 2000
**     eraSp00      the TIO locator s', IERS 2000
**     eraPom00     polar motion matrix
**     eraC2tcio    form CIO-based celestial-to-terrestrial matrix
**
** Reference:
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double rc2i[3][3], era, sp, rpom[3][3];


/* Form the celestial-to-intermediate matrix for this TT. */
   eraC2ixy(tta, ttb, x, y, rc2i);

/* Predict the Earth rotation angle for this UT1. */
   era = eraEra00(uta, utb);

/* Estimate s'. */
   sp = eraSp00(tta, ttb);

/* Form the polar motion matrix. */
   eraPom00(xp, yp, sp, rpom);

/* Combine to form the celestial-to-terrestrial matrix. */
   eraC2tcio(rc2i, era, rpom, rc2t);

   return;

}

int eraCal2jd(int iy, int im, int id, double *djm0, double *djm)
/*
**  - - - - - - - - - -
**   e r a C a l 2 j d
**  - - - - - - - - - -
**
**  Gregorian Calendar to Julian Date.
**
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
**
**  Reference:
**
**     Explanatory Supplement to the Astronomical Almanac,
**     P. Kenneth Seidelmann (ed), University Science Books (1992),
**     Section 12.92 (p604).
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   int j, ly, my;
   long iypmy;

/* Earliest year allowed (4800BC) */
   const int IYMIN = -4799;

/* Month lengths in days */
   static const int mtab[]
                     = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};


/* Preset status. */
   j = 0;

/* Validate year and month. */
   if (iy < IYMIN) return -1;
   if (im < 1 || im > 12) return -2;

/* If February in a leap year, 1, otherwise 0. */
   ly = ((im == 2) && !(iy%4) && (iy%100 || !(iy%400)));

/* Validate day, taking into account leap years. */
   if ( (id < 1) || (id > (mtab[im-1] + ly))) j = -3;

/* Return result. */
   my = (im - 14) / 12;
   iypmy = (long) (iy + my);
   *djm0 = 2400000.5;
   *djm = (double)((1461L * (iypmy + 4800L)) / 4L
                 + (367L * (long) (im - 2 - 12 * my)) / 12L
                 - (3L * ((iypmy + 4900L) / 100L)) / 4L
                 + (long) id - 2432076L);

/* Return status. */
   return j;

}

void eraCp(double p[3], double c[3])
/*
**  - - - - - -
**   e r a C p
**  - - - - - -
**
**  Copy a p-vector.
**
**  Given:
**     p        double[3]     p-vector to be copied
**
**  Returned:
**     c        double[3]     copy
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   c[0] = p[0];
   c[1] = p[1];
   c[2] = p[2];

   return;

}

void eraCpv(double pv[2][3], double c[2][3])
/*
**  - - - - - - -
**   e r a C p v
**  - - - - - - -
**
**  Copy a position/velocity vector.
**
**  Given:
**     pv     double[2][3]    position/velocity vector to be copied
**
**  Returned:
**     c      double[2][3]    copy
**
**  Called:
**     eraCp        copy p-vector
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   eraCp(pv[0], c[0]);
   eraCp(pv[1], c[1]);

   return;

}

void eraCr(double r[3][3], double c[3][3])
/*
**  - - - - - -
**   e r a C r
**  - - - - - -
**
**  Copy an r-matrix.
**
**  Given:
**     r        double[3][3]    r-matrix to be copied
**
**  Returned:
**   char[]     double[3][3]    copy
**
**  Called:
**     eraCp        copy p-vector
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   eraCp(r[0], c[0]);
   eraCp(r[1], c[1]);
   eraCp(r[2], c[2]);

   return;

}

int eraD2dtf(const char *scale, int ndp, double d1, double d2,
             int *iy, int *im, int *id, int ihmsf[4])
/*
**  - - - - - - - - -
**   e r a D 2 d t f
**  - - - - - - - - -
**
**  Format for output a 2-part Julian Date (or in the case of UTC a
**  quasi-JD form that includes special provision for leap seconds).
**
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
**
**  Called:
**     eraJd2cal    JD to Gregorian calendar
**     eraD2tf      decompose days to hms
**     eraDat       delta(AT) = TAI-UTC
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   int leap;
   char s;
   int iy1, im1, id1, js, iy2, im2, id2, ihmsf1[4], i;
   double a1, b1, fd, dat1, w, dat2, ddt;


/* The two-part JD. */
   a1 = d1;
   b1 = d2;

/* Provisional calendar date. */
   js = eraJd2cal(a1, b1, &iy1, &im1, &id1, &fd);
   if ( js ) return js < 0 ? -1 : js;

/* Is this a leap second day? */
   leap = 0;
   if ( ! strcmp(scale,"UTC") ) {

   /* TAI-UTC today. */
      js = eraDat(iy1, im1, id1, fd, &dat1);
      if ( js < 0 ) return -1;

   /* TAI-UTC tomorrow (at noon, to avoid rounding effects). */
      js = eraJd2cal(a1+1.5, b1-fd, &iy2, &im2, &id2, &w);
      js = eraDat(iy2, im2, id2, 0.0, &dat2);
      if ( js < 0 ) return -1;

   /* The change in TAI-UTC (seconds). */
      ddt = dat2 - dat1;

   /* If leap second day, scale the fraction of a day into SI. */
      leap = fabs(ddt) > 0.5;
      if (leap) fd += fd * ddt/DAYSEC;
   }

/* Provisional time of day. */
   eraD2tf ( ndp, fd, &s, ihmsf1 );

/* Is this a leap second day? */
   if ( ! leap ) {

   /* No.  Has the time rounded up to 24h? */
      if ( ihmsf1[0] > 23 ) {

      /* Yes.  We will need tomorrow's calendar date. */
         js = eraJd2cal(a1+1.5, b1-fd, &iy2, &im2, &id2, &w);

      /* Use 0h tomorrow. */
         iy1 = iy2;
         im1 = im2;
         id1 = id2;
         for ( i = 0; i < 4; i++ ) {
            ihmsf1[i] = 0;
         }
      }
   } else {

   /* This is a leap second day.  Has the time reached or passed 24h? */
      if ( ihmsf1[0] > 23 ) {

      /* Yes.  Use 23 59 60... */
         ihmsf1[0] = 23;
         ihmsf1[1] = 59;
         ihmsf1[2] = 60;
      }
   }

/* Results. */
   *iy = iy1;
   *im = im1;
   *id = id1;
   for ( i = 0; i < 4; i++ ) {
      ihmsf[i] = ihmsf1[i];
   }

/* Status. */
   return js < 0 ? -1 : js;

}

void eraD2tf(int ndp, double days, char *sign, int ihmsf[4])
/*
**  - - - - - - - -
**   e r a D 2 t f
**  - - - - - - - -
**
**  Decompose days to hours, minutes, seconds, fraction.
**
**  Given:
**     ndp     int     resolution (Note 1)
**     days    double  interval in days
**
**  Returned:
**     sign    char    '+' or '-'
**     ihmsf   int[4]  hours, minutes, seconds, fraction
**
**  Notes:
**
**  1) The argument ndp is interpreted as follows:
**
**     ndp         resolution
**      :      ...0000 00 00
**     -7         1000 00 00
**     -6          100 00 00
**     -5           10 00 00
**     -4            1 00 00
**     -3            0 10 00
**     -2            0 01 00
**     -1            0 00 10
**      0            0 00 01
**      1            0 00 00.1
**      2            0 00 00.01
**      3            0 00 00.001
**      :            0 00 00.000...
**
**  2) The largest positive useful value for ndp is determined by the
**     size of days, the format of double on the target platform, and
**     the risk of overflowing ihmsf[3].  On a typical platform, for
**     days up to 1.0, the available floating-point precision might
**     correspond to ndp=12.  However, the practical limit is typically
**     ndp=9, set by the capacity of a 32-bit int, or ndp=4 if int is
**     only 16 bits.
**
**  3) The absolute value of days may exceed 1.0.  In cases where it
**     does not, it is up to the caller to test for and handle the
**     case where days is very nearly 1.0 and rounds up to 24 hours,
**     by testing for ihmsf[0]=24 and setting ihmsf[0-3] to zero.
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   int nrs, n;
   double rs, rm, rh, a, w, ah, am, as, af;


/* Handle sign. */
   *sign = (char) ( ( days >= 0.0 ) ? '+' : '-' );

/* Interval in seconds. */
   a = DAYSEC * fabs(days);

/* Pre-round if resolution coarser than 1s (then pretend ndp=1). */
   if (ndp < 0) {
      nrs = 1;
      for (n = 1; n <= -ndp; n++) {
          nrs *= (n == 2 || n == 4) ? 6 : 10;
      }
      rs = (double) nrs;
      w = a / rs;
      a = rs * dnint(w);
   }

/* Express the unit of each field in resolution units. */
   nrs = 1;
   for (n = 1; n <= ndp; n++) {
      nrs *= 10;
   }
   rs = (double) nrs;
   rm = rs * 60.0;
   rh = rm * 60.0;

/* Round the interval and express in resolution units. */
   a = dnint(rs * a);

/* Break into fields. */
   ah = a / rh;
   ah = dint(ah);
   a -= ah * rh;
   am = a / rm;
   am = dint(am);
   a -= am * rm;
   as = a / rs;
   as = dint(as);
   af = a - as * rs;

/* Return results. */
   ihmsf[0] = (int) ah;
   ihmsf[1] = (int) am;
   ihmsf[2] = (int) as;
   ihmsf[3] = (int) af;

   return;

}

int eraDat(int iy, int im, int id, double fd, double *deltat )
/*
**  - - - - - - -
**   e r a D a t
**  - - - - - - -
**
**  For a given UTC date, calculate delta(AT) = TAI-UTC.
**
**     :------------------------------------------:
**     :                                          :
**     :                 IMPORTANT                :
**     :                                          :
**     :  A new version of this function must be  :
**     :  produced whenever a new leap second is  :
**     :  announced.  There are four items to     :
**     :  change on each such occasion:           :
**     :                                          :
**     :  1) A new line must be added to the set  :
**     :     of statements that initialize the    :
**     :     array "changes".                     :
**     :                                          :
**     :  2) The parameter IYV must be set to     :
**     :     the current year.                    :
**     :                                          :
**     :  3) The "Latest leap second" comment     :
**     :     below must be set to the new leap    :
**     :     second date.                         :
**     :                                          :
**     :  4) The "This revision" comment, later,  :
**     :     must be set to the current date.     :
**     :                                          :
**     :  Change (2) must also be carried out     :
**     :  whenever the function is re-issued,     :
**     :  even if no leap seconds have been       :
**     :  added.                                  :
**     :                                          :
**     :  Latest leap second:  2012 June 30       :
**     :                                          :
**     :__________________________________________:
**
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
**     end of 1971.  It is tested for validity (0 to 1 is the valid
**     range) even if not used;  if invalid, zero is used and status
**     j=-4 is returned.  For many applications, setting fd to zero is
**     acceptable;  the resulting error is always less than 3 ms (and
**     occurs only pre-1972).
**
**  5) The status value returned in the case where there are multiple
**     errors refers to the first error detected.  For example, if the
**     month and day are 13 and 32 respectively, j=-2 (bad month)
**     will be returned.
**
**  6) In cases where a valid result is not available, zero is returned.
**
**  References:
**
**  1) For dates from 1961 January 1 onwards, the expressions from the
**     file ftp://maia.usno.navy.mil/ser7/tai-utc.dat are used.
**
**  2) The 5ms timestep at 1961 January 1 is taken from 2.58.1 (p87) of
**     the 1992 Explanatory Supplement.
**
**  Called:
**     eraCal2jd    Gregorian calendar to Julian Day number
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
/* Release year for this version of eraDat */
#define IYV (2012)

/* Reference dates (MJD) and drift rates (s/day), pre leap seconds */
   static const double drift[][2] = {
      { 37300.0, 0.0012960 },
      { 37300.0, 0.0012960 },
      { 37300.0, 0.0012960 },
      { 37665.0, 0.0011232 },
      { 37665.0, 0.0011232 },
      { 38761.0, 0.0012960 },
      { 38761.0, 0.0012960 },
      { 38761.0, 0.0012960 },
      { 38761.0, 0.0012960 },
      { 38761.0, 0.0012960 },
      { 38761.0, 0.0012960 },
      { 38761.0, 0.0012960 },
      { 39126.0, 0.0025920 },
      { 39126.0, 0.0025920 }
   };

/* Number of Delta(AT) expressions before leap seconds were introduced */
#define NERA1 ((int) (sizeof drift / sizeof (double) / 2))

/* Dates and Delta(AT)s */
   static const struct {
      int iyear, month;
      double delat;
   } changes[] = {
      { 1960,  1,  1.4178180 },
      { 1961,  1,  1.4228180 },
      { 1961,  8,  1.3728180 },
      { 1962,  1,  1.8458580 },
      { 1963, 11,  1.9458580 },
      { 1964,  1,  3.2401300 },
      { 1964,  4,  3.3401300 },
      { 1964,  9,  3.4401300 },
      { 1965,  1,  3.5401300 },
      { 1965,  3,  3.6401300 },
      { 1965,  7,  3.7401300 },
      { 1965,  9,  3.8401300 },
      { 1966,  1,  4.3131700 },
      { 1968,  2,  4.2131700 },
      { 1972,  1, 10.0       },
      { 1972,  7, 11.0       },
      { 1973,  1, 12.0       },
      { 1974,  1, 13.0       },
      { 1975,  1, 14.0       },
      { 1976,  1, 15.0       },
      { 1977,  1, 16.0       },
      { 1978,  1, 17.0       },
      { 1979,  1, 18.0       },
      { 1980,  1, 19.0       },
      { 1981,  7, 20.0       },
      { 1982,  7, 21.0       },
      { 1983,  7, 22.0       },
      { 1985,  7, 23.0       },
      { 1988,  1, 24.0       },
      { 1990,  1, 25.0       },
      { 1991,  1, 26.0       },
      { 1992,  7, 27.0       },
      { 1993,  7, 28.0       },
      { 1994,  7, 29.0       },
      { 1996,  1, 30.0       },
      { 1997,  7, 31.0       },
      { 1999,  1, 32.0       },
      { 2006,  1, 33.0       },
      { 2009,  1, 34.0       },
      { 2012,  7, 35.0       }
   };

/* Number of Delta(AT) changes */
   const int NDAT = sizeof changes / sizeof changes[0];

/* Miscellaneous local variables */
   int j, i, m;
   double da, djm0, djm;


/* Initialize the result to zero. */
   *deltat = da = 0.0;

/* If invalid fraction of a day, set error status and give up. */
   if (fd < 0.0 || fd > 1.0) return -4;

/* Convert the date into an MJD. */
   j = eraCal2jd(iy, im, id, &djm0, &djm);

/* If invalid year, month, or day, give up. */
   if (j < 0) return j;

/* If pre-UTC year, set warning status and give up. */
   if (iy < changes[0].iyear) return 1;

/* If suspiciously late year, set warning status but proceed. */
   if (iy > IYV + 5) j = 1;

/* Combine year and month to form a date-ordered integer... */
   m = 12*iy + im;

/* ...and use it to find the preceding table entry. */
   for (i = NDAT-1; i >=0; i--) {
      if (m >= (12 * changes[i].iyear + changes[i].month)) break;
   }

/* Get the Delta(AT). */
   da = changes[i].delat;

/* If pre-1972, adjust for drift. */
   if (i < NERA1) da += (djm + fd - drift[i][0]) * drift[i][1];

/* Return the Delta(AT) value. */
   *deltat = da;

/* Return the status. */
   return j;

}

double eraDtdb(double date1, double date2,
               double ut, double elong, double u, double v)
/*
**  - - - - - - - -
**   e r a D t d b
**  - - - - - - - -
**
**  An approximation to TDB-TT, the difference between barycentric
**  dynamical time and terrestrial time, for an observer on the Earth.
**
**  The different time scales - proper, coordinate and realized - are
**  related to each other:
**
**            TAI             <-  physically realized
**             :
**          offset            <-  observed (nominally +32.184s)
**             :
**            TT              <-  terrestrial time
**             :
**    rate adjustment (L_G)   <-  definition of TT
**             :
**            TCG             <-  time scale for GCRS
**             :
**      "periodic" terms      <-  eraDtdb  is an implementation
**             :
**    rate adjustment (L_C)   <-  function of solar-system ephemeris
**             :
**            TCB             <-  time scale for BCRS
**             :
**    rate adjustment (-L_B)  <-  definition of TDB
**             :
**            TDB             <-  TCB scaled to track TT
**             :
**      "periodic" terms      <-  -eraDtdb is an approximation
**             :
**            TT              <-  terrestrial time
**
**  Adopted values for the various constants can be found in the IERS
**  Conventions (McCarthy & Petit 2003).
**
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
**
**  References:
**
**     Fairhead, L., & Bretagnon, P., Astron.Astrophys., 229, 240-247
**     (1990).
**
**     IAU 2006 Resolution 3.
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**     Moyer, T.D., Cel.Mech., 23, 33 (1981).
**
**     Murray, C.A., Vectorial Astrometry, Adam Hilger (1983).
**
**     Seidelmann, P.K. et al., Explanatory Supplement to the
**     Astronomical Almanac, Chapter 2, University Science Books (1992).
**
**     Simon, J.L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
**     Francou, G. & Laskar, J., Astron.Astrophys., 282, 663-683 (1994).
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double t, tsol, w, elsun, emsun, d, elj, els, wt, w0, w1, w2, w3, w4,
          wf, wj;
   int j;

/*
** =====================
** Fairhead et al. model
** =====================
**
** 787 sets of three coefficients.
**
** Each set is
**    amplitude (microseconds)
**      frequency (radians per Julian millennium since J2000.0)
**      phase (radians)
**
** Sets   1-474 are the T**0 terms
**  "   475-679  "   "  T**1
**  "   680-764  "   "  T**2
**  "   765-784  "   "  T**3
**  "   785-787  "   "  T**4
*/

   static const double fairhd[787][3] = {
   /* 1, 10 */
      { 1656.674564e-6,     6283.075849991,  6.240054195 },
      {   22.417471e-6,     5753.384884897,  4.296977442 },
      {   13.839792e-6,    12566.151699983,  6.196904410 },
      {    4.770086e-6,      529.690965095,  0.444401603 },
      {    4.676740e-6,     6069.776754553,  4.021195093 },
      {    2.256707e-6,      213.299095438,  5.543113262 },
      {    1.694205e-6,      -3.523118349,   5.025132748 },
      {    1.554905e-6,    77713.771467920,  5.198467090 },
      {    1.276839e-6,     7860.419392439,  5.988822341 },
      {    1.193379e-6,     5223.693919802,  3.649823730 },
   /* 11, 20 */
      {    1.115322e-6,     3930.209696220,  1.422745069 },
      {    0.794185e-6,    11506.769769794,  2.322313077 },
      {    0.447061e-6,       26.298319800,  3.615796498 },
      {    0.435206e-6,     -398.149003408,  4.349338347 },
      {    0.600309e-6,     1577.343542448,  2.678271909 },
      {    0.496817e-6,     6208.294251424,  5.696701824 },
      {    0.486306e-6,     5884.926846583,  0.520007179 },
      {    0.432392e-6,       74.781598567,  2.435898309 },
      {    0.468597e-6,     6244.942814354,  5.866398759 },
      {    0.375510e-6,     5507.553238667,  4.103476804 },
   /* 21, 30 */
      {    0.243085e-6,     -775.522611324,  3.651837925 },
      {    0.173435e-6,    18849.227549974,  6.153743485 },
      {    0.230685e-6,     5856.477659115,  4.773852582 },
      {    0.203747e-6,    12036.460734888,  4.333987818 },
      {    0.143935e-6,     -796.298006816,  5.957517795 },
      {    0.159080e-6,    10977.078804699,  1.890075226 },
      {    0.119979e-6,       38.133035638,  4.551585768 },
      {    0.118971e-6,     5486.777843175,  1.914547226 },
      {    0.116120e-6,     1059.381930189,  0.873504123 },
      {    0.137927e-6,    11790.629088659,  1.135934669 },
   /* 31, 40 */
      {    0.098358e-6,     2544.314419883,  0.092793886 },
      {    0.101868e-6,    -5573.142801634,  5.984503847 },
      {    0.080164e-6,      206.185548437,  2.095377709 },
      {    0.079645e-6,     4694.002954708,  2.949233637 },
      {    0.062617e-6,       20.775395492,  2.654394814 },
      {    0.075019e-6,     2942.463423292,  4.980931759 },
      {    0.064397e-6,     5746.271337896,  1.280308748 },
      {    0.063814e-6,     5760.498431898,  4.167901731 },
      {    0.048042e-6,     2146.165416475,  1.495846011 },
      {    0.048373e-6,      155.420399434,  2.251573730 },
   /* 41, 50 */
      {    0.058844e-6,      426.598190876,  4.839650148 },
      {    0.046551e-6,       -0.980321068,  0.921573539 },
      {    0.054139e-6,    17260.154654690,  3.411091093 },
      {    0.042411e-6,     6275.962302991,  2.869567043 },
      {    0.040184e-6,       -7.113547001,  3.565975565 },
      {    0.036564e-6,     5088.628839767,  3.324679049 },
      {    0.040759e-6,    12352.852604545,  3.981496998 },
      {    0.036507e-6,      801.820931124,  6.248866009 },
      {    0.036955e-6,     3154.687084896,  5.071801441 },
      {    0.042732e-6,      632.783739313,  5.720622217 },
   /* 51, 60 */
      {    0.042560e-6,   161000.685737473,  1.270837679 },
      {    0.040480e-6,    15720.838784878,  2.546610123 },
      {    0.028244e-6,    -6286.598968340,  5.069663519 },
      {    0.033477e-6,     6062.663207553,  4.144987272 },
      {    0.034867e-6,      522.577418094,  5.210064075 },
      {    0.032438e-6,     6076.890301554,  0.749317412 },
      {    0.030215e-6,     7084.896781115,  3.389610345 },
      {    0.029247e-6,   -71430.695617928,  4.183178762 },
      {    0.033529e-6,     9437.762934887,  2.404714239 },
      {    0.032423e-6,     8827.390269875,  5.541473556 },
   /* 61, 70 */
      {    0.027567e-6,     6279.552731642,  5.040846034 },
      {    0.029862e-6,    12139.553509107,  1.770181024 },
      {    0.022509e-6,    10447.387839604,  1.460726241 },
      {    0.020937e-6,     8429.241266467,  0.652303414 },
      {    0.020322e-6,      419.484643875,  3.735430632 },
      {    0.024816e-6,    -1194.447010225,  1.087136918 },
      {    0.025196e-6,     1748.016413067,  2.901883301 },
      {    0.021691e-6,    14143.495242431,  5.952658009 },
      {    0.017673e-6,     6812.766815086,  3.186129845 },
      {    0.022567e-6,     6133.512652857,  3.307984806 },
   /* 71, 80 */
      {    0.016155e-6,    10213.285546211,  1.331103168 },
      {    0.014751e-6,     1349.867409659,  4.308933301 },
      {    0.015949e-6,     -220.412642439,  4.005298270 },
      {    0.015974e-6,    -2352.866153772,  6.145309371 },
      {    0.014223e-6,    17789.845619785,  2.104551349 },
      {    0.017806e-6,       73.297125859,  3.475975097 },
      {    0.013671e-6,     -536.804512095,  5.971672571 },
      {    0.011942e-6,     8031.092263058,  2.053414715 },
      {    0.014318e-6,    16730.463689596,  3.016058075 },
      {    0.012462e-6,      103.092774219,  1.737438797 },
   /* 81, 90 */
      {    0.010962e-6,        3.590428652,  2.196567739 },
      {    0.015078e-6,    19651.048481098,  3.969480770 },
      {    0.010396e-6,      951.718406251,  5.717799605 },
      {    0.011707e-6,    -4705.732307544,  2.654125618 },
      {    0.010453e-6,     5863.591206116,  1.913704550 },
      {    0.012420e-6,     4690.479836359,  4.734090399 },
      {    0.011847e-6,     5643.178563677,  5.489005403 },
      {    0.008610e-6,     3340.612426700,  3.661698944 },
      {    0.011622e-6,     5120.601145584,  4.863931876 },
      {    0.010825e-6,      553.569402842,  0.842715011 },
   /* 91, 100 */
      {    0.008666e-6,     -135.065080035,  3.293406547 },
      {    0.009963e-6,      149.563197135,  4.870690598 },
      {    0.009858e-6,     6309.374169791,  1.061816410 },
      {    0.007959e-6,      316.391869657,  2.465042647 },
      {    0.010099e-6,      283.859318865,  1.942176992 },
      {    0.007147e-6,     -242.728603974,  3.661486981 },
      {    0.007505e-6,     5230.807466803,  4.920937029 },
      {    0.008323e-6,    11769.853693166,  1.229392026 },
      {    0.007490e-6,    -6256.777530192,  3.658444681 },
      {    0.009370e-6,   149854.400134205,  0.673880395 },
   /* 101, 110 */
      {    0.007117e-6,       38.027672636,  5.294249518 },
      {    0.007857e-6,    12168.002696575,  0.525733528 },
      {    0.007019e-6,     6206.809778716,  0.837688810 },
      {    0.006056e-6,      955.599741609,  4.194535082 },
      {    0.008107e-6,    13367.972631107,  3.793235253 },
      {    0.006731e-6,     5650.292110678,  5.639906583 },
      {    0.007332e-6,       36.648562930,  0.114858677 },
      {    0.006366e-6,     4164.311989613,  2.262081818 },
      {    0.006858e-6,     5216.580372801,  0.642063318 },
      {    0.006919e-6,     6681.224853400,  6.018501522 },
   /* 111, 120 */
      {    0.006826e-6,     7632.943259650,  3.458654112 },
      {    0.005308e-6,    -1592.596013633,  2.500382359 },
      {    0.005096e-6,    11371.704689758,  2.547107806 },
      {    0.004841e-6,     5333.900241022,  0.437078094 },
      {    0.005582e-6,     5966.683980335,  2.246174308 },
      {    0.006304e-6,    11926.254413669,  2.512929171 },
      {    0.006603e-6,    23581.258177318,  5.393136889 },
      {    0.005123e-6,       -1.484472708,  2.999641028 },
      {    0.004648e-6,     1589.072895284,  1.275847090 },
      {    0.005119e-6,     6438.496249426,  1.486539246 },
   /* 121, 130 */
      {    0.004521e-6,     4292.330832950,  6.140635794 },
      {    0.005680e-6,    23013.539539587,  4.557814849 },
      {    0.005488e-6,       -3.455808046,  0.090675389 },
      {    0.004193e-6,     7234.794256242,  4.869091389 },
      {    0.003742e-6,     7238.675591600,  4.691976180 },
      {    0.004148e-6,     -110.206321219,  3.016173439 },
      {    0.004553e-6,    11499.656222793,  5.554998314 },
      {    0.004892e-6,     5436.993015240,  1.475415597 },
      {    0.004044e-6,     4732.030627343,  1.398784824 },
      {    0.004164e-6,    12491.370101415,  5.650931916 },
   /* 131, 140 */
      {    0.004349e-6,    11513.883316794,  2.181745369 },
      {    0.003919e-6,    12528.018664345,  5.823319737 },
      {    0.003129e-6,     6836.645252834,  0.003844094 },
      {    0.004080e-6,    -7058.598461315,  3.690360123 },
      {    0.003270e-6,       76.266071276,  1.517189902 },
      {    0.002954e-6,     6283.143160294,  4.447203799 },
      {    0.002872e-6,       28.449187468,  1.158692983 },
      {    0.002881e-6,      735.876513532,  0.349250250 },
      {    0.003279e-6,     5849.364112115,  4.893384368 },
      {    0.003625e-6,     6209.778724132,  1.473760578 },
   /* 141, 150 */
      {    0.003074e-6,      949.175608970,  5.185878737 },
      {    0.002775e-6,     9917.696874510,  1.030026325 },
      {    0.002646e-6,    10973.555686350,  3.918259169 },
      {    0.002575e-6,    25132.303399966,  6.109659023 },
      {    0.003500e-6,      263.083923373,  1.892100742 },
      {    0.002740e-6,    18319.536584880,  4.320519510 },
      {    0.002464e-6,      202.253395174,  4.698203059 },
      {    0.002409e-6,        2.542797281,  5.325009315 },
      {    0.003354e-6,   -90955.551694697,  1.942656623 },
      {    0.002296e-6,     6496.374945429,  5.061810696 },
   /* 151, 160 */
      {    0.003002e-6,     6172.869528772,  2.797822767 },
      {    0.003202e-6,    27511.467873537,  0.531673101 },
      {    0.002954e-6,    -6283.008539689,  4.533471191 },
      {    0.002353e-6,      639.897286314,  3.734548088 },
      {    0.002401e-6,    16200.772724501,  2.605547070 },
      {    0.003053e-6,   233141.314403759,  3.029030662 },
      {    0.003024e-6,    83286.914269554,  2.355556099 },
      {    0.002863e-6,    17298.182327326,  5.240963796 },
      {    0.002103e-6,    -7079.373856808,  5.756641637 },
      {    0.002303e-6,    83996.847317911,  2.013686814 },
   /* 161, 170 */
      {    0.002303e-6,    18073.704938650,  1.089100410 },
      {    0.002381e-6,       63.735898303,  0.759188178 },
      {    0.002493e-6,     6386.168624210,  0.645026535 },
      {    0.002366e-6,        3.932153263,  6.215885448 },
      {    0.002169e-6,    11015.106477335,  4.845297676 },
      {    0.002397e-6,     6243.458341645,  3.809290043 },
      {    0.002183e-6,     1162.474704408,  6.179611691 },
      {    0.002353e-6,     6246.427287062,  4.781719760 },
      {    0.002199e-6,     -245.831646229,  5.956152284 },
      {    0.001729e-6,     3894.181829542,  1.264976635 },
   /* 171, 180 */
      {    0.001896e-6,    -3128.388765096,  4.914231596 },
      {    0.002085e-6,       35.164090221,  1.405158503 },
      {    0.002024e-6,    14712.317116458,  2.752035928 },
      {    0.001737e-6,     6290.189396992,  5.280820144 },
      {    0.002229e-6,      491.557929457,  1.571007057 },
      {    0.001602e-6,    14314.168113050,  4.203664806 },
      {    0.002186e-6,      454.909366527,  1.402101526 },
      {    0.001897e-6,    22483.848574493,  4.167932508 },
      {    0.001825e-6,    -3738.761430108,  0.545828785 },
      {    0.001894e-6,     1052.268383188,  5.817167450 },
   /* 181, 190 */
      {    0.001421e-6,       20.355319399,  2.419886601 },
      {    0.001408e-6,    10984.192351700,  2.732084787 },
      {    0.001847e-6,    10873.986030480,  2.903477885 },
      {    0.001391e-6,    -8635.942003763,  0.593891500 },
      {    0.001388e-6,       -7.046236698,  1.166145902 },
      {    0.001810e-6,   -88860.057071188,  0.487355242 },
      {    0.001288e-6,    -1990.745017041,  3.913022880 },
      {    0.001297e-6,    23543.230504682,  3.063805171 },
      {    0.001335e-6,     -266.607041722,  3.995764039 },
      {    0.001376e-6,    10969.965257698,  5.152914309 },
   /* 191, 200 */
      {    0.001745e-6,   244287.600007027,  3.626395673 },
      {    0.001649e-6,    31441.677569757,  1.952049260 },
      {    0.001416e-6,     9225.539273283,  4.996408389 },
      {    0.001238e-6,     4804.209275927,  5.503379738 },
      {    0.001472e-6,     4590.910180489,  4.164913291 },
      {    0.001169e-6,     6040.347246017,  5.841719038 },
      {    0.001039e-6,     5540.085789459,  2.769753519 },
      {    0.001004e-6,     -170.672870619,  0.755008103 },
      {    0.001284e-6,    10575.406682942,  5.306538209 },
      {    0.001278e-6,       71.812653151,  4.713486491 },
   /* 201, 210 */
      {    0.001321e-6,    18209.330263660,  2.624866359 },
      {    0.001297e-6,    21228.392023546,  0.382603541 },
      {    0.000954e-6,     6282.095528923,  0.882213514 },
      {    0.001145e-6,     6058.731054289,  1.169483931 },
      {    0.000979e-6,     5547.199336460,  5.448375984 },
      {    0.000987e-6,    -6262.300454499,  2.656486959 },
      {    0.001070e-6,  -154717.609887482,  1.827624012 },
      {    0.000991e-6,     4701.116501708,  4.387001801 },
      {    0.001155e-6,      -14.227094002,  3.042700750 },
      {    0.001176e-6,      277.034993741,  3.335519004 },
   /* 211, 220 */
      {    0.000890e-6,    13916.019109642,  5.601498297 },
      {    0.000884e-6,    -1551.045222648,  1.088831705 },
      {    0.000876e-6,     5017.508371365,  3.969902609 },
      {    0.000806e-6,    15110.466119866,  5.142876744 },
      {    0.000773e-6,    -4136.910433516,  0.022067765 },
      {    0.001077e-6,      175.166059800,  1.844913056 },
      {    0.000954e-6,    -6284.056171060,  0.968480906 },
      {    0.000737e-6,     5326.786694021,  4.923831588 },
      {    0.000845e-6,     -433.711737877,  4.749245231 },
      {    0.000819e-6,     8662.240323563,  5.991247817 },
   /* 221, 230 */
      {    0.000852e-6,      199.072001436,  2.189604979 },
      {    0.000723e-6,    17256.631536341,  6.068719637 },
      {    0.000940e-6,     6037.244203762,  6.197428148 },
      {    0.000885e-6,    11712.955318231,  3.280414875 },
      {    0.000706e-6,    12559.038152982,  2.824848947 },
      {    0.000732e-6,     2379.164473572,  2.501813417 },
      {    0.000764e-6,    -6127.655450557,  2.236346329 },
      {    0.000908e-6,      131.541961686,  2.521257490 },
      {    0.000907e-6,    35371.887265976,  3.370195967 },
      {    0.000673e-6,     1066.495477190,  3.876512374 },
   /* 231, 240 */
      {    0.000814e-6,    17654.780539750,  4.627122566 },
      {    0.000630e-6,       36.027866677,  0.156368499 },
      {    0.000798e-6,      515.463871093,  5.151962502 },
      {    0.000798e-6,      148.078724426,  5.909225055 },
      {    0.000806e-6,      309.278322656,  6.054064447 },
      {    0.000607e-6,      -39.617508346,  2.839021623 },
      {    0.000601e-6,      412.371096874,  3.984225404 },
      {    0.000646e-6,    11403.676995575,  3.852959484 },
      {    0.000704e-6,    13521.751441591,  2.300991267 },
      {    0.000603e-6,   -65147.619767937,  4.140083146 },
   /* 241, 250 */
      {    0.000609e-6,    10177.257679534,  0.437122327 },
      {    0.000631e-6,     5767.611978898,  4.026532329 },
      {    0.000576e-6,    11087.285125918,  4.760293101 },
      {    0.000674e-6,    14945.316173554,  6.270510511 },
      {    0.000726e-6,     5429.879468239,  6.039606892 },
      {    0.000710e-6,    28766.924424484,  5.672617711 },
      {    0.000647e-6,    11856.218651625,  3.397132627 },
      {    0.000678e-6,    -5481.254918868,  6.249666675 },
      {    0.000618e-6,    22003.914634870,  2.466427018 },
      {    0.000738e-6,     6134.997125565,  2.242668890 },
   /* 251, 260 */
      {    0.000660e-6,      625.670192312,  5.864091907 },
      {    0.000694e-6,     3496.032826134,  2.668309141 },
      {    0.000531e-6,     6489.261398429,  1.681888780 },
      {    0.000611e-6,  -143571.324284214,  2.424978312 },
      {    0.000575e-6,    12043.574281889,  4.216492400 },
      {    0.000553e-6,    12416.588502848,  4.772158039 },
      {    0.000689e-6,     4686.889407707,  6.224271088 },
      {    0.000495e-6,     7342.457780181,  3.817285811 },
      {    0.000567e-6,     3634.621024518,  1.649264690 },
      {    0.000515e-6,    18635.928454536,  3.945345892 },
   /* 261, 270 */
      {    0.000486e-6,     -323.505416657,  4.061673868 },
      {    0.000662e-6,    25158.601719765,  1.794058369 },
      {    0.000509e-6,      846.082834751,  3.053874588 },
      {    0.000472e-6,   -12569.674818332,  5.112133338 },
      {    0.000461e-6,     6179.983075773,  0.513669325 },
      {    0.000641e-6,    83467.156352816,  3.210727723 },
      {    0.000520e-6,    10344.295065386,  2.445597761 },
      {    0.000493e-6,    18422.629359098,  1.676939306 },
      {    0.000478e-6,     1265.567478626,  5.487314569 },
      {    0.000472e-6,      -18.159247265,  1.999707589 },
   /* 271, 280 */
      {    0.000559e-6,    11190.377900137,  5.783236356 },
      {    0.000494e-6,     9623.688276691,  3.022645053 },
      {    0.000463e-6,     5739.157790895,  1.411223013 },
      {    0.000432e-6,    16858.482532933,  1.179256434 },
      {    0.000574e-6,    72140.628666286,  1.758191830 },
      {    0.000484e-6,    17267.268201691,  3.290589143 },
      {    0.000550e-6,     4907.302050146,  0.864024298 },
      {    0.000399e-6,       14.977853527,  2.094441910 },
      {    0.000491e-6,      224.344795702,  0.878372791 },
      {    0.000432e-6,    20426.571092422,  6.003829241 },
   /* 281, 290 */
      {    0.000481e-6,     5749.452731634,  4.309591964 },
      {    0.000480e-6,     5757.317038160,  1.142348571 },
      {    0.000485e-6,     6702.560493867,  0.210580917 },
      {    0.000426e-6,     6055.549660552,  4.274476529 },
      {    0.000480e-6,     5959.570433334,  5.031351030 },
      {    0.000466e-6,    12562.628581634,  4.959581597 },
      {    0.000520e-6,    39302.096962196,  4.788002889 },
      {    0.000458e-6,    12132.439962106,  1.880103788 },
      {    0.000470e-6,    12029.347187887,  1.405611197 },
      {    0.000416e-6,    -7477.522860216,  1.082356330 },
   /* 291, 300 */
      {    0.000449e-6,    11609.862544012,  4.179989585 },
      {    0.000465e-6,    17253.041107690,  0.353496295 },
      {    0.000362e-6,    -4535.059436924,  1.583849576 },
      {    0.000383e-6,    21954.157609398,  3.747376371 },
      {    0.000389e-6,       17.252277143,  1.395753179 },
      {    0.000331e-6,    18052.929543158,  0.566790582 },
      {    0.000430e-6,    13517.870106233,  0.685827538 },
      {    0.000368e-6,    -5756.908003246,  0.731374317 },
      {    0.000330e-6,    10557.594160824,  3.710043680 },
      {    0.000332e-6,    20199.094959633,  1.652901407 },
   /* 301, 310 */
      {    0.000384e-6,    11933.367960670,  5.827781531 },
      {    0.000387e-6,    10454.501386605,  2.541182564 },
      {    0.000325e-6,    15671.081759407,  2.178850542 },
      {    0.000318e-6,      138.517496871,  2.253253037 },
      {    0.000305e-6,     9388.005909415,  0.578340206 },
      {    0.000352e-6,     5749.861766548,  3.000297967 },
      {    0.000311e-6,     6915.859589305,  1.693574249 },
      {    0.000297e-6,    24072.921469776,  1.997249392 },
      {    0.000363e-6,     -640.877607382,  5.071820966 },
      {    0.000323e-6,    12592.450019783,  1.072262823 },
   /* 311, 320 */
      {    0.000341e-6,    12146.667056108,  4.700657997 },
      {    0.000290e-6,     9779.108676125,  1.812320441 },
      {    0.000342e-6,     6132.028180148,  4.322238614 },
      {    0.000329e-6,     6268.848755990,  3.033827743 },
      {    0.000374e-6,    17996.031168222,  3.388716544 },
      {    0.000285e-6,     -533.214083444,  4.687313233 },
      {    0.000338e-6,     6065.844601290,  0.877776108 },
      {    0.000276e-6,       24.298513841,  0.770299429 },
      {    0.000336e-6,    -2388.894020449,  5.353796034 },
      {    0.000290e-6,     3097.883822726,  4.075291557 },
   /* 321, 330 */
      {    0.000318e-6,      709.933048357,  5.941207518 },
      {    0.000271e-6,    13095.842665077,  3.208912203 },
      {    0.000331e-6,     6073.708907816,  4.007881169 },
      {    0.000292e-6,      742.990060533,  2.714333592 },
      {    0.000362e-6,    29088.811415985,  3.215977013 },
      {    0.000280e-6,    12359.966151546,  0.710872502 },
      {    0.000267e-6,    10440.274292604,  4.730108488 },
      {    0.000262e-6,      838.969287750,  1.327720272 },
      {    0.000250e-6,    16496.361396202,  0.898769761 },
      {    0.000325e-6,    20597.243963041,  0.180044365 },
   /* 331, 340 */
      {    0.000268e-6,     6148.010769956,  5.152666276 },
      {    0.000284e-6,     5636.065016677,  5.655385808 },
      {    0.000301e-6,     6080.822454817,  2.135396205 },
      {    0.000294e-6,     -377.373607916,  3.708784168 },
      {    0.000236e-6,     2118.763860378,  1.733578756 },
      {    0.000234e-6,     5867.523359379,  5.575209112 },
      {    0.000268e-6,  -226858.238553767,  0.069432392 },
      {    0.000265e-6,   167283.761587465,  4.369302826 },
      {    0.000280e-6,    28237.233459389,  5.304829118 },
      {    0.000292e-6,    12345.739057544,  4.096094132 },
   /* 341, 350 */
      {    0.000223e-6,    19800.945956225,  3.069327406 },
      {    0.000301e-6,    43232.306658416,  6.205311188 },
      {    0.000264e-6,    18875.525869774,  1.417263408 },
      {    0.000304e-6,    -1823.175188677,  3.409035232 },
      {    0.000301e-6,      109.945688789,  0.510922054 },
      {    0.000260e-6,      813.550283960,  2.389438934 },
      {    0.000299e-6,   316428.228673312,  5.384595078 },
      {    0.000211e-6,     5756.566278634,  3.789392838 },
      {    0.000209e-6,     5750.203491159,  1.661943545 },
      {    0.000240e-6,    12489.885628707,  5.684549045 },
   /* 351, 360 */
      {    0.000216e-6,     6303.851245484,  3.862942261 },
      {    0.000203e-6,     1581.959348283,  5.549853589 },
      {    0.000200e-6,     5642.198242609,  1.016115785 },
      {    0.000197e-6,      -70.849445304,  4.690702525 },
      {    0.000227e-6,     6287.008003254,  2.911891613 },
      {    0.000197e-6,      533.623118358,  1.048982898 },
      {    0.000205e-6,    -6279.485421340,  1.829362730 },
      {    0.000209e-6,   -10988.808157535,  2.636140084 },
      {    0.000208e-6,     -227.526189440,  4.127883842 },
      {    0.000191e-6,      415.552490612,  4.401165650 },
   /* 361, 370 */
      {    0.000190e-6,    29296.615389579,  4.175658539 },
      {    0.000264e-6,    66567.485864652,  4.601102551 },
      {    0.000256e-6,    -3646.350377354,  0.506364778 },
      {    0.000188e-6,    13119.721102825,  2.032195842 },
      {    0.000185e-6,     -209.366942175,  4.694756586 },
      {    0.000198e-6,    25934.124331089,  3.832703118 },
      {    0.000195e-6,     4061.219215394,  3.308463427 },
      {    0.000234e-6,     5113.487598583,  1.716090661 },
      {    0.000188e-6,     1478.866574064,  5.686865780 },
      {    0.000222e-6,    11823.161639450,  1.942386641 },
   /* 371, 380 */
      {    0.000181e-6,    10770.893256262,  1.999482059 },
      {    0.000171e-6,     6546.159773364,  1.182807992 },
      {    0.000206e-6,       70.328180442,  5.934076062 },
      {    0.000169e-6,    20995.392966449,  2.169080622 },
      {    0.000191e-6,    10660.686935042,  5.405515999 },
      {    0.000228e-6,    33019.021112205,  4.656985514 },
      {    0.000184e-6,    -4933.208440333,  3.327476868 },
      {    0.000220e-6,     -135.625325010,  1.765430262 },
      {    0.000166e-6,    23141.558382925,  3.454132746 },
      {    0.000191e-6,     6144.558353121,  5.020393445 },
   /* 381, 390 */
      {    0.000180e-6,     6084.003848555,  0.602182191 },
      {    0.000163e-6,    17782.732072784,  4.960593133 },
      {    0.000225e-6,    16460.333529525,  2.596451817 },
      {    0.000222e-6,     5905.702242076,  3.731990323 },
      {    0.000204e-6,      227.476132789,  5.636192701 },
      {    0.000159e-6,    16737.577236597,  3.600691544 },
      {    0.000200e-6,     6805.653268085,  0.868220961 },
      {    0.000187e-6,    11919.140866668,  2.629456641 },
      {    0.000161e-6,      127.471796607,  2.862574720 },
      {    0.000205e-6,     6286.666278643,  1.742882331 },
   /* 391, 400 */
      {    0.000189e-6,      153.778810485,  4.812372643 },
      {    0.000168e-6,    16723.350142595,  0.027860588 },
      {    0.000149e-6,    11720.068865232,  0.659721876 },
      {    0.000189e-6,     5237.921013804,  5.245313000 },
      {    0.000143e-6,     6709.674040867,  4.317625647 },
      {    0.000146e-6,     4487.817406270,  4.815297007 },
      {    0.000144e-6,     -664.756045130,  5.381366880 },
      {    0.000175e-6,     5127.714692584,  4.728443327 },
      {    0.000162e-6,     6254.626662524,  1.435132069 },
      {    0.000187e-6,    47162.516354635,  1.354371923 },
   /* 401, 410 */
      {    0.000146e-6,    11080.171578918,  3.369695406 },
      {    0.000180e-6,     -348.924420448,  2.490902145 },
      {    0.000148e-6,      151.047669843,  3.799109588 },
      {    0.000157e-6,     6197.248551160,  1.284375887 },
      {    0.000167e-6,      146.594251718,  0.759969109 },
      {    0.000133e-6,    -5331.357443741,  5.409701889 },
      {    0.000154e-6,       95.979227218,  3.366890614 },
      {    0.000148e-6,    -6418.140930027,  3.384104996 },
      {    0.000128e-6,    -6525.804453965,  3.803419985 },
      {    0.000130e-6,    11293.470674356,  0.939039445 },
   /* 411, 420 */
      {    0.000152e-6,    -5729.506447149,  0.734117523 },
      {    0.000138e-6,      210.117701700,  2.564216078 },
      {    0.000123e-6,     6066.595360816,  4.517099537 },
      {    0.000140e-6,    18451.078546566,  0.642049130 },
      {    0.000126e-6,    11300.584221356,  3.485280663 },
      {    0.000119e-6,    10027.903195729,  3.217431161 },
      {    0.000151e-6,     4274.518310832,  4.404359108 },
      {    0.000117e-6,     6072.958148291,  0.366324650 },
      {    0.000165e-6,    -7668.637425143,  4.298212528 },
      {    0.000117e-6,    -6245.048177356,  5.379518958 },
   /* 421, 430 */
      {    0.000130e-6,    -5888.449964932,  4.527681115 },
      {    0.000121e-6,     -543.918059096,  6.109429504 },
      {    0.000162e-6,     9683.594581116,  5.720092446 },
      {    0.000141e-6,     6219.339951688,  0.679068671 },
      {    0.000118e-6,    22743.409379516,  4.881123092 },
      {    0.000129e-6,     1692.165669502,  0.351407289 },
      {    0.000126e-6,     5657.405657679,  5.146592349 },
      {    0.000114e-6,      728.762966531,  0.520791814 },
      {    0.000120e-6,       52.596639600,  0.948516300 },
      {    0.000115e-6,       65.220371012,  3.504914846 },
   /* 431, 440 */
      {    0.000126e-6,     5881.403728234,  5.577502482 },
      {    0.000158e-6,   163096.180360983,  2.957128968 },
      {    0.000134e-6,    12341.806904281,  2.598576764 },
      {    0.000151e-6,    16627.370915377,  3.985702050 },
      {    0.000109e-6,     1368.660252845,  0.014730471 },
      {    0.000131e-6,     6211.263196841,  0.085077024 },
      {    0.000146e-6,     5792.741760812,  0.708426604 },
      {    0.000146e-6,      -77.750543984,  3.121576600 },
      {    0.000107e-6,     5341.013788022,  0.288231904 },
      {    0.000138e-6,     6281.591377283,  2.797450317 },
   /* 441, 450 */
      {    0.000113e-6,    -6277.552925684,  2.788904128 },
      {    0.000115e-6,     -525.758811831,  5.895222200 },
      {    0.000138e-6,     6016.468808270,  6.096188999 },
      {    0.000139e-6,    23539.707386333,  2.028195445 },
      {    0.000146e-6,    -4176.041342449,  4.660008502 },
      {    0.000107e-6,    16062.184526117,  4.066520001 },
      {    0.000142e-6,    83783.548222473,  2.936315115 },
      {    0.000128e-6,     9380.959672717,  3.223844306 },
      {    0.000135e-6,     6205.325306007,  1.638054048 },
      {    0.000101e-6,     2699.734819318,  5.481603249 },
   /* 451, 460 */
      {    0.000104e-6,     -568.821874027,  2.205734493 },
      {    0.000103e-6,     6321.103522627,  2.440421099 },
      {    0.000119e-6,     6321.208885629,  2.547496264 },
      {    0.000138e-6,     1975.492545856,  2.314608466 },
      {    0.000121e-6,      137.033024162,  4.539108237 },
      {    0.000123e-6,    19402.796952817,  4.538074405 },
      {    0.000119e-6,    22805.735565994,  2.869040566 },
      {    0.000133e-6,    64471.991241142,  6.056405489 },
      {    0.000129e-6,      -85.827298831,  2.540635083 },
      {    0.000131e-6,    13613.804277336,  4.005732868 },
   /* 461, 470 */
      {    0.000104e-6,     9814.604100291,  1.959967212 },
      {    0.000112e-6,    16097.679950283,  3.589026260 },
      {    0.000123e-6,     2107.034507542,  1.728627253 },
      {    0.000121e-6,    36949.230808424,  6.072332087 },
      {    0.000108e-6,   -12539.853380183,  3.716133846 },
      {    0.000113e-6,    -7875.671863624,  2.725771122 },
      {    0.000109e-6,     4171.425536614,  4.033338079 },
      {    0.000101e-6,     6247.911759770,  3.441347021 },
      {    0.000113e-6,     7330.728427345,  0.656372122 },
      {    0.000113e-6,    51092.726050855,  2.791483066 },
   /* 471, 480 */
      {    0.000106e-6,     5621.842923210,  1.815323326 },
      {    0.000101e-6,      111.430161497,  5.711033677 },
      {    0.000103e-6,      909.818733055,  2.812745443 },
      {    0.000101e-6,     1790.642637886,  1.965746028 },

   /* T */
      {  102.156724e-6,     6283.075849991,  4.249032005 },
      {    1.706807e-6,    12566.151699983,  4.205904248 },
      {    0.269668e-6,      213.299095438,  3.400290479 },
      {    0.265919e-6,      529.690965095,  5.836047367 },
      {    0.210568e-6,       -3.523118349,  6.262738348 },
      {    0.077996e-6,     5223.693919802,  4.670344204 },
   /* 481, 490 */
      {    0.054764e-6,     1577.343542448,  4.534800170 },
      {    0.059146e-6,       26.298319800,  1.083044735 },
      {    0.034420e-6,     -398.149003408,  5.980077351 },
      {    0.032088e-6,    18849.227549974,  4.162913471 },
      {    0.033595e-6,     5507.553238667,  5.980162321 },
      {    0.029198e-6,     5856.477659115,  0.623811863 },
      {    0.027764e-6,      155.420399434,  3.745318113 },
      {    0.025190e-6,     5746.271337896,  2.980330535 },
      {    0.022997e-6,     -796.298006816,  1.174411803 },
      {    0.024976e-6,     5760.498431898,  2.467913690 },
   /* 491, 500 */
      {    0.021774e-6,      206.185548437,  3.854787540 },
      {    0.017925e-6,     -775.522611324,  1.092065955 },
      {    0.013794e-6,      426.598190876,  2.699831988 },
      {    0.013276e-6,     6062.663207553,  5.845801920 },
      {    0.011774e-6,    12036.460734888,  2.292832062 },
      {    0.012869e-6,     6076.890301554,  5.333425680 },
      {    0.012152e-6,     1059.381930189,  6.222874454 },
      {    0.011081e-6,       -7.113547001,  5.154724984 },
      {    0.010143e-6,     4694.002954708,  4.044013795 },
      {    0.009357e-6,     5486.777843175,  3.416081409 },
   /* 501, 510 */
      {    0.010084e-6,      522.577418094,  0.749320262 },
      {    0.008587e-6,    10977.078804699,  2.777152598 },
      {    0.008628e-6,     6275.962302991,  4.562060226 },
      {    0.008158e-6,     -220.412642439,  5.806891533 },
      {    0.007746e-6,     2544.314419883,  1.603197066 },
      {    0.007670e-6,     2146.165416475,  3.000200440 },
      {    0.007098e-6,       74.781598567,  0.443725817 },
      {    0.006180e-6,     -536.804512095,  1.302642751 },
      {    0.005818e-6,     5088.628839767,  4.827723531 },
      {    0.004945e-6,    -6286.598968340,  0.268305170 },
   /* 511, 520 */
      {    0.004774e-6,     1349.867409659,  5.808636673 },
      {    0.004687e-6,     -242.728603974,  5.154890570 },
      {    0.006089e-6,     1748.016413067,  4.403765209 },
      {    0.005975e-6,    -1194.447010225,  2.583472591 },
      {    0.004229e-6,      951.718406251,  0.931172179 },
      {    0.005264e-6,      553.569402842,  2.336107252 },
      {    0.003049e-6,     5643.178563677,  1.362634430 },
      {    0.002974e-6,     6812.766815086,  1.583012668 },
      {    0.003403e-6,    -2352.866153772,  2.552189886 },
      {    0.003030e-6,      419.484643875,  5.286473844 },
   /* 521, 530 */
      {    0.003210e-6,       -7.046236698,  1.863796539 },
      {    0.003058e-6,     9437.762934887,  4.226420633 },
      {    0.002589e-6,    12352.852604545,  1.991935820 },
      {    0.002927e-6,     5216.580372801,  2.319951253 },
      {    0.002425e-6,     5230.807466803,  3.084752833 },
      {    0.002656e-6,     3154.687084896,  2.487447866 },
      {    0.002445e-6,    10447.387839604,  2.347139160 },
      {    0.002990e-6,     4690.479836359,  6.235872050 },
      {    0.002890e-6,     5863.591206116,  0.095197563 },
      {    0.002498e-6,     6438.496249426,  2.994779800 },
   /* 531, 540 */
      {    0.001889e-6,     8031.092263058,  3.569003717 },
      {    0.002567e-6,      801.820931124,  3.425611498 },
      {    0.001803e-6,   -71430.695617928,  2.192295512 },
      {    0.001782e-6,        3.932153263,  5.180433689 },
      {    0.001694e-6,    -4705.732307544,  4.641779174 },
      {    0.001704e-6,    -1592.596013633,  3.997097652 },
      {    0.001735e-6,     5849.364112115,  0.417558428 },
      {    0.001643e-6,     8429.241266467,  2.180619584 },
      {    0.001680e-6,       38.133035638,  4.164529426 },
      {    0.002045e-6,     7084.896781115,  0.526323854 },
   /* 541, 550 */
      {    0.001458e-6,     4292.330832950,  1.356098141 },
      {    0.001437e-6,       20.355319399,  3.895439360 },
      {    0.001738e-6,     6279.552731642,  0.087484036 },
      {    0.001367e-6,    14143.495242431,  3.987576591 },
      {    0.001344e-6,     7234.794256242,  0.090454338 },
      {    0.001438e-6,    11499.656222793,  0.974387904 },
      {    0.001257e-6,     6836.645252834,  1.509069366 },
      {    0.001358e-6,    11513.883316794,  0.495572260 },
      {    0.001628e-6,     7632.943259650,  4.968445721 },
      {    0.001169e-6,      103.092774219,  2.838496795 },
   /* 551, 560 */
      {    0.001162e-6,     4164.311989613,  3.408387778 },
      {    0.001092e-6,     6069.776754553,  3.617942651 },
      {    0.001008e-6,    17789.845619785,  0.286350174 },
      {    0.001008e-6,      639.897286314,  1.610762073 },
      {    0.000918e-6,    10213.285546211,  5.532798067 },
      {    0.001011e-6,    -6256.777530192,  0.661826484 },
      {    0.000753e-6,    16730.463689596,  3.905030235 },
      {    0.000737e-6,    11926.254413669,  4.641956361 },
      {    0.000694e-6,     3340.612426700,  2.111120332 },
      {    0.000701e-6,     3894.181829542,  2.760823491 },
   /* 561, 570 */
      {    0.000689e-6,     -135.065080035,  4.768800780 },
      {    0.000700e-6,    13367.972631107,  5.760439898 },
      {    0.000664e-6,     6040.347246017,  1.051215840 },
      {    0.000654e-6,     5650.292110678,  4.911332503 },
      {    0.000788e-6,     6681.224853400,  4.699648011 },
      {    0.000628e-6,     5333.900241022,  5.024608847 },
      {    0.000755e-6,     -110.206321219,  4.370971253 },
      {    0.000628e-6,     6290.189396992,  3.660478857 },
      {    0.000635e-6,    25132.303399966,  4.121051532 },
      {    0.000534e-6,     5966.683980335,  1.173284524 },
   /* 571, 580 */
      {    0.000543e-6,     -433.711737877,  0.345585464 },
      {    0.000517e-6,    -1990.745017041,  5.414571768 },
      {    0.000504e-6,     5767.611978898,  2.328281115 },
      {    0.000485e-6,     5753.384884897,  1.685874771 },
      {    0.000463e-6,     7860.419392439,  5.297703006 },
      {    0.000604e-6,      515.463871093,  0.591998446 },
      {    0.000443e-6,    12168.002696575,  4.830881244 },
      {    0.000570e-6,      199.072001436,  3.899190272 },
      {    0.000465e-6,    10969.965257698,  0.476681802 },
      {    0.000424e-6,    -7079.373856808,  1.112242763 },
   /* 581, 590 */
      {    0.000427e-6,      735.876513532,  1.994214480 },
      {    0.000478e-6,    -6127.655450557,  3.778025483 },
      {    0.000414e-6,    10973.555686350,  5.441088327 },
      {    0.000512e-6,     1589.072895284,  0.107123853 },
      {    0.000378e-6,    10984.192351700,  0.915087231 },
      {    0.000402e-6,    11371.704689758,  4.107281715 },
      {    0.000453e-6,     9917.696874510,  1.917490952 },
      {    0.000395e-6,      149.563197135,  2.763124165 },
      {    0.000371e-6,     5739.157790895,  3.112111866 },
      {    0.000350e-6,    11790.629088659,  0.440639857 },
   /* 591, 600 */
      {    0.000356e-6,     6133.512652857,  5.444568842 },
      {    0.000344e-6,      412.371096874,  5.676832684 },
      {    0.000383e-6,      955.599741609,  5.559734846 },
      {    0.000333e-6,     6496.374945429,  0.261537984 },
      {    0.000340e-6,     6055.549660552,  5.975534987 },
      {    0.000334e-6,     1066.495477190,  2.335063907 },
      {    0.000399e-6,    11506.769769794,  5.321230910 },
      {    0.000314e-6,    18319.536584880,  2.313312404 },
      {    0.000424e-6,     1052.268383188,  1.211961766 },
      {    0.000307e-6,       63.735898303,  3.169551388 },
   /* 601, 610 */
      {    0.000329e-6,       29.821438149,  6.106912080 },
      {    0.000357e-6,     6309.374169791,  4.223760346 },
      {    0.000312e-6,    -3738.761430108,  2.180556645 },
      {    0.000301e-6,      309.278322656,  1.499984572 },
      {    0.000268e-6,    12043.574281889,  2.447520648 },
      {    0.000257e-6,    12491.370101415,  3.662331761 },
      {    0.000290e-6,      625.670192312,  1.272834584 },
      {    0.000256e-6,     5429.879468239,  1.913426912 },
      {    0.000339e-6,     3496.032826134,  4.165930011 },
      {    0.000283e-6,     3930.209696220,  4.325565754 },
   /* 611, 620 */
      {    0.000241e-6,    12528.018664345,  3.832324536 },
      {    0.000304e-6,     4686.889407707,  1.612348468 },
      {    0.000259e-6,    16200.772724501,  3.470173146 },
      {    0.000238e-6,    12139.553509107,  1.147977842 },
      {    0.000236e-6,     6172.869528772,  3.776271728 },
      {    0.000296e-6,    -7058.598461315,  0.460368852 },
      {    0.000306e-6,    10575.406682942,  0.554749016 },
      {    0.000251e-6,    17298.182327326,  0.834332510 },
      {    0.000290e-6,     4732.030627343,  4.759564091 },
      {    0.000261e-6,     5884.926846583,  0.298259862 },
   /* 621, 630 */
      {    0.000249e-6,     5547.199336460,  3.749366406 },
      {    0.000213e-6,    11712.955318231,  5.415666119 },
      {    0.000223e-6,     4701.116501708,  2.703203558 },
      {    0.000268e-6,     -640.877607382,  0.283670793 },
      {    0.000209e-6,     5636.065016677,  1.238477199 },
      {    0.000193e-6,    10177.257679534,  1.943251340 },
      {    0.000182e-6,     6283.143160294,  2.456157599 },
      {    0.000184e-6,     -227.526189440,  5.888038582 },
      {    0.000182e-6,    -6283.008539689,  0.241332086 },
      {    0.000228e-6,    -6284.056171060,  2.657323816 },
   /* 631, 640 */
      {    0.000166e-6,     7238.675591600,  5.930629110 },
      {    0.000167e-6,     3097.883822726,  5.570955333 },
      {    0.000159e-6,     -323.505416657,  5.786670700 },
      {    0.000154e-6,    -4136.910433516,  1.517805532 },
      {    0.000176e-6,    12029.347187887,  3.139266834 },
      {    0.000167e-6,    12132.439962106,  3.556352289 },
      {    0.000153e-6,      202.253395174,  1.463313961 },
      {    0.000157e-6,    17267.268201691,  1.586837396 },
      {    0.000142e-6,    83996.847317911,  0.022670115 },
      {    0.000152e-6,    17260.154654690,  0.708528947 },
   /* 641, 650 */
      {    0.000144e-6,     6084.003848555,  5.187075177 },
      {    0.000135e-6,     5756.566278634,  1.993229262 },
      {    0.000134e-6,     5750.203491159,  3.457197134 },
      {    0.000144e-6,     5326.786694021,  6.066193291 },
      {    0.000160e-6,    11015.106477335,  1.710431974 },
      {    0.000133e-6,     3634.621024518,  2.836451652 },
      {    0.000134e-6,    18073.704938650,  5.453106665 },
      {    0.000134e-6,     1162.474704408,  5.326898811 },
      {    0.000128e-6,     5642.198242609,  2.511652591 },
      {    0.000160e-6,      632.783739313,  5.628785365 },
   /* 651, 660 */
      {    0.000132e-6,    13916.019109642,  0.819294053 },
      {    0.000122e-6,    14314.168113050,  5.677408071 },
      {    0.000125e-6,    12359.966151546,  5.251984735 },
      {    0.000121e-6,     5749.452731634,  2.210924603 },
      {    0.000136e-6,     -245.831646229,  1.646502367 },
      {    0.000120e-6,     5757.317038160,  3.240883049 },
      {    0.000134e-6,    12146.667056108,  3.059480037 },
      {    0.000137e-6,     6206.809778716,  1.867105418 },
      {    0.000141e-6,    17253.041107690,  2.069217456 },
      {    0.000129e-6,    -7477.522860216,  2.781469314 },
   /* 661, 670 */
      {    0.000116e-6,     5540.085789459,  4.281176991 },
      {    0.000116e-6,     9779.108676125,  3.320925381 },
      {    0.000129e-6,     5237.921013804,  3.497704076 },
      {    0.000113e-6,     5959.570433334,  0.983210840 },
      {    0.000122e-6,     6282.095528923,  2.674938860 },
      {    0.000140e-6,      -11.045700264,  4.957936982 },
      {    0.000108e-6,    23543.230504682,  1.390113589 },
      {    0.000106e-6,   -12569.674818332,  0.429631317 },
      {    0.000110e-6,     -266.607041722,  5.501340197 },
      {    0.000115e-6,    12559.038152982,  4.691456618 },
   /* 671, 680 */
      {    0.000134e-6,    -2388.894020449,  0.577313584 },
      {    0.000109e-6,    10440.274292604,  6.218148717 },
      {    0.000102e-6,     -543.918059096,  1.477842615 },
      {    0.000108e-6,    21228.392023546,  2.237753948 },
      {    0.000101e-6,    -4535.059436924,  3.100492232 },
      {    0.000103e-6,       76.266071276,  5.594294322 },
      {    0.000104e-6,      949.175608970,  5.674287810 },
      {    0.000101e-6,    13517.870106233,  2.196632348 },
      {    0.000100e-6,    11933.367960670,  4.056084160 },

   /* T^2 */
      {    4.322990e-6,     6283.075849991,  2.642893748 },
   /* 681, 690 */
      {    0.406495e-6,        0.000000000,  4.712388980 },
      {    0.122605e-6,    12566.151699983,  2.438140634 },
      {    0.019476e-6,      213.299095438,  1.642186981 },
      {    0.016916e-6,      529.690965095,  4.510959344 },
      {    0.013374e-6,       -3.523118349,  1.502210314 },
      {    0.008042e-6,       26.298319800,  0.478549024 },
      {    0.007824e-6,      155.420399434,  5.254710405 },
      {    0.004894e-6,     5746.271337896,  4.683210850 },
      {    0.004875e-6,     5760.498431898,  0.759507698 },
      {    0.004416e-6,     5223.693919802,  6.028853166 },
   /* 691, 700 */
      {    0.004088e-6,       -7.113547001,  0.060926389 },
      {    0.004433e-6,    77713.771467920,  3.627734103 },
      {    0.003277e-6,    18849.227549974,  2.327912542 },
      {    0.002703e-6,     6062.663207553,  1.271941729 },
      {    0.003435e-6,     -775.522611324,  0.747446224 },
      {    0.002618e-6,     6076.890301554,  3.633715689 },
      {    0.003146e-6,      206.185548437,  5.647874613 },
      {    0.002544e-6,     1577.343542448,  6.232904270 },
      {    0.002218e-6,     -220.412642439,  1.309509946 },
      {    0.002197e-6,     5856.477659115,  2.407212349 },
   /* 701, 710 */
      {    0.002897e-6,     5753.384884897,  5.863842246 },
      {    0.001766e-6,      426.598190876,  0.754113147 },
      {    0.001738e-6,     -796.298006816,  2.714942671 },
      {    0.001695e-6,      522.577418094,  2.629369842 },
      {    0.001584e-6,     5507.553238667,  1.341138229 },
      {    0.001503e-6,     -242.728603974,  0.377699736 },
      {    0.001552e-6,     -536.804512095,  2.904684667 },
      {    0.001370e-6,     -398.149003408,  1.265599125 },
      {    0.001889e-6,    -5573.142801634,  4.413514859 },
      {    0.001722e-6,     6069.776754553,  2.445966339 },
   /* 711, 720 */
      {    0.001124e-6,     1059.381930189,  5.041799657 },
      {    0.001258e-6,      553.569402842,  3.849557278 },
      {    0.000831e-6,      951.718406251,  2.471094709 },
      {    0.000767e-6,     4694.002954708,  5.363125422 },
      {    0.000756e-6,     1349.867409659,  1.046195744 },
      {    0.000775e-6,      -11.045700264,  0.245548001 },
      {    0.000597e-6,     2146.165416475,  4.543268798 },
      {    0.000568e-6,     5216.580372801,  4.178853144 },
      {    0.000711e-6,     1748.016413067,  5.934271972 },
      {    0.000499e-6,    12036.460734888,  0.624434410 },
   /* 721, 730 */
      {    0.000671e-6,    -1194.447010225,  4.136047594 },
      {    0.000488e-6,     5849.364112115,  2.209679987 },
      {    0.000621e-6,     6438.496249426,  4.518860804 },
      {    0.000495e-6,    -6286.598968340,  1.868201275 },
      {    0.000456e-6,     5230.807466803,  1.271231591 },
      {    0.000451e-6,     5088.628839767,  0.084060889 },
      {    0.000435e-6,     5643.178563677,  3.324456609 },
      {    0.000387e-6,    10977.078804699,  4.052488477 },
      {    0.000547e-6,   161000.685737473,  2.841633844 },
      {    0.000522e-6,     3154.687084896,  2.171979966 },
   /* 731, 740 */
      {    0.000375e-6,     5486.777843175,  4.983027306 },
      {    0.000421e-6,     5863.591206116,  4.546432249 },
      {    0.000439e-6,     7084.896781115,  0.522967921 },
      {    0.000309e-6,     2544.314419883,  3.172606705 },
      {    0.000347e-6,     4690.479836359,  1.479586566 },
      {    0.000317e-6,      801.820931124,  3.553088096 },
      {    0.000262e-6,      419.484643875,  0.606635550 },
      {    0.000248e-6,     6836.645252834,  3.014082064 },
      {    0.000245e-6,    -1592.596013633,  5.519526220 },
      {    0.000225e-6,     4292.330832950,  2.877956536 },
   /* 741, 750 */
      {    0.000214e-6,     7234.794256242,  1.605227587 },
      {    0.000205e-6,     5767.611978898,  0.625804796 },
      {    0.000180e-6,    10447.387839604,  3.499954526 },
      {    0.000229e-6,      199.072001436,  5.632304604 },
      {    0.000214e-6,      639.897286314,  5.960227667 },
      {    0.000175e-6,     -433.711737877,  2.162417992 },
      {    0.000209e-6,      515.463871093,  2.322150893 },
      {    0.000173e-6,     6040.347246017,  2.556183691 },
      {    0.000184e-6,     6309.374169791,  4.732296790 },
      {    0.000227e-6,   149854.400134205,  5.385812217 },
   /* 751, 760 */
      {    0.000154e-6,     8031.092263058,  5.120720920 },
      {    0.000151e-6,     5739.157790895,  4.815000443 },
      {    0.000197e-6,     7632.943259650,  0.222827271 },
      {    0.000197e-6,       74.781598567,  3.910456770 },
      {    0.000138e-6,     6055.549660552,  1.397484253 },
      {    0.000149e-6,    -6127.655450557,  5.333727496 },
      {    0.000137e-6,     3894.181829542,  4.281749907 },
      {    0.000135e-6,     9437.762934887,  5.979971885 },
      {    0.000139e-6,    -2352.866153772,  4.715630782 },
      {    0.000142e-6,     6812.766815086,  0.513330157 },
   /* 761, 770 */
      {    0.000120e-6,    -4705.732307544,  0.194160689 },
      {    0.000131e-6,   -71430.695617928,  0.000379226 },
      {    0.000124e-6,     6279.552731642,  2.122264908 },
      {    0.000108e-6,    -6256.777530192,  0.883445696 },

   /* T^3 */
      {    0.143388e-6,     6283.075849991,  1.131453581 },
      {    0.006671e-6,    12566.151699983,  0.775148887 },
      {    0.001480e-6,      155.420399434,  0.480016880 },
      {    0.000934e-6,      213.299095438,  6.144453084 },
      {    0.000795e-6,      529.690965095,  2.941595619 },
      {    0.000673e-6,     5746.271337896,  0.120415406 },
   /* 771, 780 */
      {    0.000672e-6,     5760.498431898,  5.317009738 },
      {    0.000389e-6,     -220.412642439,  3.090323467 },
      {    0.000373e-6,     6062.663207553,  3.003551964 },
      {    0.000360e-6,     6076.890301554,  1.918913041 },
      {    0.000316e-6,      -21.340641002,  5.545798121 },
      {    0.000315e-6,     -242.728603974,  1.884932563 },
      {    0.000278e-6,      206.185548437,  1.266254859 },
      {    0.000238e-6,     -536.804512095,  4.532664830 },
      {    0.000185e-6,      522.577418094,  4.578313856 },
      {    0.000245e-6,    18849.227549974,  0.587467082 },
   /* 781, 787 */
      {    0.000180e-6,      426.598190876,  5.151178553 },
      {    0.000200e-6,      553.569402842,  5.355983739 },
      {    0.000141e-6,     5223.693919802,  1.336556009 },
      {    0.000104e-6,     5856.477659115,  4.239842759 },

   /* T^4 */
      {    0.003826e-6,     6283.075849991,  5.705257275 },
      {    0.000303e-6,    12566.151699983,  5.407132842 },
      {    0.000209e-6,      155.420399434,  1.989815753 }
   };


/* Time since J2000.0 in Julian millennia. */
   t = ((date1 - DJ00) + date2) / DJM;

/* ================= */
/* Topocentric terms */
/* ================= */

/* Convert UT to local solar time in radians. */
   tsol = fmod(ut, 1.0) * D2PI + elong;

/* FUNDAMENTAL ARGUMENTS:  Simon et al. 1994. */

/* Combine time argument (millennia) with deg/arcsec factor. */
   w = t / 3600.0;

/* Sun Mean Longitude. */
   elsun = fmod(280.46645683 + 1296027711.03429 * w, 360.0) * DD2R;

/* Sun Mean Anomaly. */
   emsun = fmod(357.52910918 + 1295965810.481 * w, 360.0) * DD2R;

/* Mean Elongation of Moon from Sun. */
   d = fmod(297.85019547 + 16029616012.090 * w, 360.0) * DD2R;

/* Mean Longitude of Jupiter. */
   elj = fmod(34.35151874 + 109306899.89453 * w, 360.0) * DD2R;

/* Mean Longitude of Saturn. */
   els = fmod(50.07744430 + 44046398.47038 * w, 360.0) * DD2R;

/* TOPOCENTRIC TERMS:  Moyer 1981 and Murray 1983. */
   wt =   +  0.00029e-10 * u * sin(tsol + elsun - els)
          +  0.00100e-10 * u * sin(tsol - 2.0 * emsun)
          +  0.00133e-10 * u * sin(tsol - d)
          +  0.00133e-10 * u * sin(tsol + elsun - elj)
          -  0.00229e-10 * u * sin(tsol + 2.0 * elsun + emsun)
          -  0.02200e-10 * v * cos(elsun + emsun)
          +  0.05312e-10 * u * sin(tsol - emsun)
          -  0.13677e-10 * u * sin(tsol + 2.0 * elsun)
          -  1.31840e-10 * v * cos(elsun)
          +  3.17679e-10 * u * sin(tsol);

/* ===================== */
/* Fairhead et al. model */
/* ===================== */

/* T**0 */
   w0 = 0;
   for (j = 473; j >= 0; j--) {
      w0 += fairhd[j][0] * sin(fairhd[j][1] * t + fairhd[j][2]);
   }

/* T**1 */
   w1 = 0;
   for (j = 678; j >= 474; j--) {
      w1 += fairhd[j][0] * sin(fairhd[j][1] * t + fairhd[j][2]);
   }

/* T**2 */
   w2 = 0;
   for (j = 763; j >= 679; j--) {
      w2 += fairhd[j][0] * sin(fairhd[j][1] * t + fairhd[j][2]);
   }

/* T**3 */
   w3 = 0;
   for (j = 783; j >= 764; j--) {
      w3 += fairhd[j][0] * sin(fairhd[j][1] * t + fairhd[j][2]);
   }

/* T**4 */
   w4 = 0;
   for (j = 786; j >= 784; j--) {
      w4 += fairhd[j][0] * sin(fairhd[j][1] * t + fairhd[j][2]);
   }

/* Multiply by powers of T and combine. */
   wf = t * (t * (t * (t * w4 + w3) + w2) + w1) + w0;

/* Adjustments to use JPL planetary masses instead of IAU. */
   wj =   0.00065e-6 * sin(6069.776754 * t + 4.021194) +
          0.00033e-6 * sin( 213.299095 * t + 5.543132) +
        (-0.00196e-6 * sin(6208.294251 * t + 5.696701)) +
        (-0.00173e-6 * sin(  74.781599 * t + 2.435900)) +
          0.03638e-6 * t * t;

/* ============ */
/* Final result */
/* ============ */

/* TDB-TT in seconds. */
   w = wt + wf + wj;

   return w;

}

int eraDtf2d(const char *scale, int iy, int im, int id,
             int ihr, int imn, double sec, double *d1, double *d2)
/*
**  - - - - - - - - -
**   e r a D t f 2 d
**  - - - - - - - - -
**
**  Encode date and time fields into 2-part Julian Date (or in the case
**  of UTC a quasi-JD form that includes special provision for leap
**  seconds).
**
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
**
**  Called:
**     eraCal2jd    Gregorian calendar to JD
**     eraDat       delta(AT) = TAI-UTC
**     eraJd2cal    JD to Gregorian calendar
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   int js, iy2, im2, id2;
   double dj, w, day, seclim, dat1, dat2, ddt, time;


/* Today's Julian Day Number. */
   js = eraCal2jd(iy, im, id, &dj, &w);
   if ( js ) return js;
   dj += w;

/* Day length and final minute length in seconds (provisional). */
   day = DAYSEC;
   seclim = 60;

/* Deal with the UTC leap second case. */
   if ( ! strcmp(scale,"UTC") ) {

   /* TAI-UTC today. */
      js = eraDat(iy, im, id, 0.0, &dat1);
      if ( js < 0 ) return js;

   /* TAI-UTC tomorrow. */
      js = eraJd2cal ( dj, 1.0, &iy2, &im2, &id2, &w);
      if ( js ) return js;
      js = eraDat(iy2, im2, id2, 0.0, &dat2);
      if ( js < 0 ) return js;

   /* The change in TAI-UTC (seconds). */
      ddt = dat2 - dat1;

   /* If leap second day, correct the day and final minute lengths. */
      if ( fabs(ddt) > 0.5 ) {
         day += ddt;
         if ( ihr == 23 && imn == 59 ) seclim += ddt;
      }
   }

/* Validate the time. */
   if ( ihr >= 0 && ihr <= 23 ) {
      if ( imn >= 0 && imn <= 59 ) {
         if ( sec >= 0 ) {
            if ( sec >= seclim ) {
               js += 2;
            }
         } else {
            js = -6;
         }
      } else {
         js = -5;
      }
   } else {
      js = -4;
   }
   if ( js < 0 ) return js;

/* The time in days. */
   time  = ( 60.0 * ( (double) ( 60 * ihr + imn ) ) + sec ) / day;

/* Return the date and time. */
   *d1 = dj;
   *d2 = time;

/* Status. */
   return js;

}

double eraEe00(double date1, double date2, double epsa, double dpsi)
/*
**  - - - - - - - -
**   e r a E e 0 0
**  - - - - - - - -
**
**  The equation of the equinoxes, compatible with IAU 2000 resolutions,
**  given the nutation in longitude and the mean obliquity.
**
**  Given:
**     date1,date2  double    TT as a 2-part Julian Date (Note 1)
**     epsa         double    mean obliquity (Note 2)
**     dpsi         double    nutation in longitude (Note 3)
**
**  Returned (function value):
**                  double    equation of the equinoxes (Note 4)
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
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
**  2) The obliquity, in radians, is mean of date.
**
**  3) The result, which is in radians, operates in the following sense:
**
**        Greenwich apparent ST = GMST + equation of the equinoxes
**
**  4) The result is compatible with the IAU 2000 resolutions.  For
**     further details, see IERS Conventions 2003 and Capitaine et al.
**     (2002).
**
**  Called:
**     eraEect00    equation of the equinoxes complementary terms
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
*/
{
   double ee;


/* Equation of the equinoxes. */
   ee = dpsi * cos(epsa) + eraEect00(date1, date2);

   return ee;

}

double eraEe00a(double date1, double date2)
/*
**  - - - - - - - - -
**   e r a E e 0 0 a
**  - - - - - - - - -
**
**  Equation of the equinoxes, compatible with IAU 2000 resolutions.
**
**  Given:
**     date1,date2  double    TT as a 2-part Julian Date (Note 1)
**
**  Returned (function value):
**                  double    equation of the equinoxes (Note 2)
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
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
**  2) The result, which is in radians, operates in the following sense:
**
**        Greenwich apparent ST = GMST + equation of the equinoxes
**
**  3) The result is compatible with the IAU 2000 resolutions.  For
**     further details, see IERS Conventions 2003 and Capitaine et al.
**     (2002).
**
**  Called:
**     eraPr00      IAU 2000 precession adjustments
**     eraObl80     mean obliquity, IAU 1980
**     eraNut00a    nutation, IAU 2000A
**     eraEe00      equation of the equinoxes, IAU 2000
**
**  References:
**
**     Capitaine, N., Wallace, P.T. and McCarthy, D.D., "Expressions to
**     implement the IAU 2000 definition of UT1", Astronomy &
**     Astrophysics, 406, 1135-1149 (2003).
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004).
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double dpsipr, depspr, epsa, dpsi, deps, ee;


/* IAU 2000 precession-rate adjustments. */
   eraPr00(date1, date2, &dpsipr, &depspr);

/* Mean obliquity, consistent with IAU 2000 precession-nutation. */
   epsa = eraObl80(date1, date2) + depspr;

/* Nutation in longitude. */
   eraNut00a(date1, date2, &dpsi, &deps);

/* Equation of the equinoxes. */
   ee = eraEe00(date1, date2, epsa, dpsi);

   return ee;

}

double eraEe00b(double date1, double date2)
/*
**  - - - - - - - - -
**   e r a E e 0 0 b
**  - - - - - - - - -
**
**  Equation of the equinoxes, compatible with IAU 2000 resolutions but
**  using the truncated nutation model IAU 2000B.
**
**  Given:
**     date1,date2  double    TT as a 2-part Julian Date (Note 1)
**
**  Returned (function value):
**                  double    equation of the equinoxes (Note 2)
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
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
**  2) The result, which is in radians, operates in the following sense:
**
**        Greenwich apparent ST = GMST + equation of the equinoxes
**
**  3) The result is compatible with the IAU 2000 resolutions except
**     that accuracy has been compromised for the sake of speed.  For
**     further details, see McCarthy & Luzum (2001), IERS Conventions
**     2003 and Capitaine et al. (2003).
**
**  Called:
**     eraPr00      IAU 2000 precession adjustments
**     eraObl80     mean obliquity, IAU 1980
**     eraNut00b    nutation, IAU 2000B
**     eraEe00      equation of the equinoxes, IAU 2000
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
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double dpsipr, depspr, epsa, dpsi, deps, ee;


/* IAU 2000 precession-rate adjustments. */
   eraPr00(date1, date2, &dpsipr, &depspr);

/* Mean obliquity, consistent with IAU 2000 precession-nutation. */
   epsa = eraObl80(date1, date2) + depspr;

/* Nutation in longitude. */
   eraNut00b(date1, date2, &dpsi, &deps);

/* Equation of the equinoxes. */
   ee = eraEe00(date1, date2, epsa, dpsi);

   return ee;

}

double eraEe06a(double date1, double date2)
/*
**  - - - - - - - - -
**   e r a E e 0 6 a
**  - - - - - - - - -
**
**  Equation of the equinoxes, compatible with IAU 2000 resolutions and
**  IAU 2006/2000A precession-nutation.
**
**  Given:
**     date1,date2  double    TT as a 2-part Julian Date (Note 1)
**
**  Returned (function value):
**                  double    equation of the equinoxes (Note 2)
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
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
**  2) The result, which is in radians, operates in the following sense:
**
**        Greenwich apparent ST = GMST + equation of the equinoxes
**
**  Called:
**     eraAnpm      normalize angle into range +/- pi
**     eraGst06a    Greenwich apparent sidereal time, IAU 2006/2000A
**     eraGmst06    Greenwich mean sidereal time, IAU 2006
**
**  Reference:
**
**     McCarthy, D. D., Petit, G. (eds.), 2004, IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double gst06a, gmst06, ee;


/* Apparent and mean sidereal times. */
   gst06a = eraGst06a(0.0, 0.0, date1, date2);
   gmst06 = eraGmst06(0.0, 0.0, date1, date2);

/* Equation of the equinoxes. */
   ee  = eraAnpm(gst06a - gmst06);

   return ee;

}

double eraEect00(double date1, double date2)
/*
**  - - - - - - - - - -
**   e r a E e c t 0 0
**  - - - - - - - - - -
**
**  Equation of the equinoxes complementary terms, consistent with
**  IAU 2000 resolutions.
**
**  Given:
**     date1,date2  double   TT as a 2-part Julian Date (Note 1)
**
**  Returned (function value):
**                  double   complementary terms (Note 2)
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
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
**  2) The "complementary terms" are part of the equation of the
**     equinoxes (EE), classically the difference between apparent and
**     mean Sidereal Time:
**
**        GAST = GMST + EE
**
**     with:
**
**        EE = dpsi * cos(eps)
**
**     where dpsi is the nutation in longitude and eps is the obliquity
**     of date.  However, if the rotation of the Earth were constant in
**     an inertial frame the classical formulation would lead to
**     apparent irregularities in the UT1 timescale traceable to side-
**     effects of precession-nutation.  In order to eliminate these
**     effects from UT1, "complementary terms" were introduced in 1994
**     (IAU, 1994) and took effect from 1997 (Capitaine and Gontier,
**     1993):
**
**        GAST = GMST + CT + EE
**
**     By convention, the complementary terms are included as part of
**     the equation of the equinoxes rather than as part of the mean
**     Sidereal Time.  This slightly compromises the "geometrical"
**     interpretation of mean sidereal time but is otherwise
**     inconsequential.
**
**     The present function computes CT in the above expression,
**     compatible with IAU 2000 resolutions (Capitaine et al., 2002, and
**     IERS Conventions 2003).
**
**  Called:
**     eraFal03     mean anomaly of the Moon
**     eraFalp03    mean anomaly of the Sun
**     eraFaf03     mean argument of the latitude of the Moon
**     eraFad03     mean elongation of the Moon from the Sun
**     eraFaom03    mean longitude of the Moon's ascending node
**     eraFave03    mean longitude of Venus
**     eraFae03     mean longitude of Earth
**     eraFapa03    general accumulated precession in longitude
**
**  References:
**
**     Capitaine, N. & Gontier, A.-M., Astron. Astrophys., 275,
**     645-650 (1993)
**
**     Capitaine, N., Wallace, P.T. and McCarthy, D.D., "Expressions to
**     implement the IAU 2000 definition of UT1", Astronomy &
**     Astrophysics, 406, 1135-1149 (2003)
**
**     IAU Resolution C7, Recommendation 3 (1994)
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
/* Time since J2000.0, in Julian centuries */
   double t;

/* Miscellaneous */
   int i, j;
   double a, s0, s1;

/* Fundamental arguments */
   double fa[14];

/* Returned value. */
   double eect;

/* ----------------------------------------- */
/* The series for the EE complementary terms */
/* ----------------------------------------- */

   typedef struct {
      int nfa[8];      /* coefficients of l,l',F,D,Om,LVe,LE,pA */
      double s, c;     /* sine and cosine coefficients */
   } TERM;

/* Terms of order t^0 */
   static const TERM e0[] = {

   /* 1-10 */
      {{ 0,  0,  0,  0,  1,  0,  0,  0}, 2640.96e-6, -0.39e-6 },
      {{ 0,  0,  0,  0,  2,  0,  0,  0},   63.52e-6, -0.02e-6 },
      {{ 0,  0,  2, -2,  3,  0,  0,  0},   11.75e-6,  0.01e-6 },
      {{ 0,  0,  2, -2,  1,  0,  0,  0},   11.21e-6,  0.01e-6 },
      {{ 0,  0,  2, -2,  2,  0,  0,  0},   -4.55e-6,  0.00e-6 },
      {{ 0,  0,  2,  0,  3,  0,  0,  0},    2.02e-6,  0.00e-6 },
      {{ 0,  0,  2,  0,  1,  0,  0,  0},    1.98e-6,  0.00e-6 },
      {{ 0,  0,  0,  0,  3,  0,  0,  0},   -1.72e-6,  0.00e-6 },
      {{ 0,  1,  0,  0,  1,  0,  0,  0},   -1.41e-6, -0.01e-6 },
      {{ 0,  1,  0,  0, -1,  0,  0,  0},   -1.26e-6, -0.01e-6 },

   /* 11-20 */
      {{ 1,  0,  0,  0, -1,  0,  0,  0},   -0.63e-6,  0.00e-6 },
      {{ 1,  0,  0,  0,  1,  0,  0,  0},   -0.63e-6,  0.00e-6 },
      {{ 0,  1,  2, -2,  3,  0,  0,  0},    0.46e-6,  0.00e-6 },
      {{ 0,  1,  2, -2,  1,  0,  0,  0},    0.45e-6,  0.00e-6 },
      {{ 0,  0,  4, -4,  4,  0,  0,  0},    0.36e-6,  0.00e-6 },
      {{ 0,  0,  1, -1,  1, -8, 12,  0},   -0.24e-6, -0.12e-6 },
      {{ 0,  0,  2,  0,  0,  0,  0,  0},    0.32e-6,  0.00e-6 },
      {{ 0,  0,  2,  0,  2,  0,  0,  0},    0.28e-6,  0.00e-6 },
      {{ 1,  0,  2,  0,  3,  0,  0,  0},    0.27e-6,  0.00e-6 },
      {{ 1,  0,  2,  0,  1,  0,  0,  0},    0.26e-6,  0.00e-6 },

   /* 21-30 */
      {{ 0,  0,  2, -2,  0,  0,  0,  0},   -0.21e-6,  0.00e-6 },
      {{ 0,  1, -2,  2, -3,  0,  0,  0},    0.19e-6,  0.00e-6 },
      {{ 0,  1, -2,  2, -1,  0,  0,  0},    0.18e-6,  0.00e-6 },
      {{ 0,  0,  0,  0,  0,  8,-13, -1},   -0.10e-6,  0.05e-6 },
      {{ 0,  0,  0,  2,  0,  0,  0,  0},    0.15e-6,  0.00e-6 },
      {{ 2,  0, -2,  0, -1,  0,  0,  0},   -0.14e-6,  0.00e-6 },
      {{ 1,  0,  0, -2,  1,  0,  0,  0},    0.14e-6,  0.00e-6 },
      {{ 0,  1,  2, -2,  2,  0,  0,  0},   -0.14e-6,  0.00e-6 },
      {{ 1,  0,  0, -2, -1,  0,  0,  0},    0.14e-6,  0.00e-6 },
      {{ 0,  0,  4, -2,  4,  0,  0,  0},    0.13e-6,  0.00e-6 },

   /* 31-33 */
      {{ 0,  0,  2, -2,  4,  0,  0,  0},   -0.11e-6,  0.00e-6 },
      {{ 1,  0, -2,  0, -3,  0,  0,  0},    0.11e-6,  0.00e-6 },
      {{ 1,  0, -2,  0, -1,  0,  0,  0},    0.11e-6,  0.00e-6 }
   };

/* Terms of order t^1 */
   static const TERM e1[] = {
      {{ 0,  0,  0,  0,  1,  0,  0,  0},    -0.87e-6,  0.00e-6 }
   };

/* Number of terms in the series */
   const int NE0 = (int) (sizeof e0 / sizeof (TERM));
   const int NE1 = (int) (sizeof e1 / sizeof (TERM));

/*--------------------------------------------------------------------*/

/* Interval between fundamental epoch J2000.0 and current date (JC). */
   t = ((date1 - DJ00) + date2) / DJC;

/* Fundamental Arguments (from IERS Conventions 2003) */

/* Mean anomaly of the Moon. */
   fa[0] = eraFal03(t);

/* Mean anomaly of the Sun. */
   fa[1] = eraFalp03(t);

/* Mean longitude of the Moon minus that of the ascending node. */
   fa[2] = eraFaf03(t);

/* Mean elongation of the Moon from the Sun. */
   fa[3] = eraFad03(t);

/* Mean longitude of the ascending node of the Moon. */
   fa[4] = eraFaom03(t);

/* Mean longitude of Venus. */
   fa[5] = eraFave03(t);

/* Mean longitude of Earth. */
   fa[6] = eraFae03(t);

/* General precession in longitude. */
   fa[7] = eraFapa03(t);

/* Evaluate the EE complementary terms. */
   s0 = 0.0;
   s1 = 0.0;

   for (i = NE0-1; i >= 0; i--) {
      a = 0.0;
      for (j = 0; j < 8; j++) {
         a += (double)(e0[i].nfa[j]) * fa[j];
      }
      s0 += e0[i].s * sin(a) + e0[i].c * cos(a);
   }

   for (i = NE1-1; i >= 0; i--) {
      a = 0.0;
      for (j = 0; j < 8; j++) {
         a += (double)(e1[i].nfa[j]) * fa[j];
      }
      s1 += e1[i].s * sin(a) + e1[i].c * cos(a);
   }

   eect = (s0 + s1 * t ) * DAS2R;

   return eect;

}

int eraEform ( int n, double *a, double *f )
/*
**  - - - - - - - - -
**   e r a E f o r m
**  - - - - - - - - -
**
**  Earth reference ellipsoids.
**
**  Given:
**     n    int         ellipsoid identifier (Note 1)
**
**  Returned:
**     a    double      equatorial radius (meters, Note 2)
**     f    double      flattening (Note 2)
**
**  Returned (function value):
**          int         status:  0 = OK
**                              -1 = illegal identifier (Note 3)
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
**  2) The ellipsoid parameters are returned in the form of equatorial
**     radius in meters (a) and flattening (f).  The latter is a number
**     around 0.00335, i.e. around 1/298.
**
**  3) For the case where an unsupported n value is supplied, zero a and
**     f are returned, as well as error status.
**
**  References:
**
**     Department of Defense World Geodetic System 1984, National
**     Imagery and Mapping Agency Technical Report 8350.2, Third
**     Edition, p3-2.
**
**     Moritz, H., Bull. Geodesique 66-2, 187 (1992).
**
**     The Department of Defense World Geodetic System 1972, World
**     Geodetic System Committee, May 1974.
**
**     Explanatory Supplement to the Astronomical Almanac,
**     P. Kenneth Seidelmann (ed), University Science Books (1992),
**     p220.
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{

/* Look up a and f for the specified reference ellipsoid. */
   switch ( n ) {

   case WGS84:
      *a = 6378137.0;
      *f = 1.0 / 298.257223563;
      break;

   case GRS80:
      *a = 6378137.0;
      *f = 1.0 / 298.257222101;
      break;

   case WGS72:
      *a = 6378135.0;
      *f = 1.0 / 298.26;
      break;

   default:

   /* Invalid identifier. */
      *a = 0.0;
      *f = 0.0;
      return -1;

   }

/* OK status. */
   return 0;

}

double eraEo06a(double date1, double date2)
/*
**  - - - - - - - - -
**   e r a E o 0 6 a
**  - - - - - - - - -
**
**  Equation of the origins, IAU 2006 precession and IAU 2000A nutation.
**
**  Given:
**     date1,date2  double    TT as a 2-part Julian Date (Note 1)
**
**  Returned (function value):
**                  double    equation of the origins in radians
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
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
**  2) The equation of the origins is the distance between the true
**     equinox and the celestial intermediate origin and, equivalently,
**     the difference between Earth rotation angle and Greenwich
**     apparent sidereal time (ERA-GST).  It comprises the precession
**     (since J2000.0) in right ascension plus the equation of the
**     equinoxes (including the small correction terms).
**
**  Called:
**     eraPnm06a    classical NPB matrix, IAU 2006/2000A
**     eraBpn2xy    extract CIP X,Y coordinates from NPB matrix
**     eraS06       the CIO locator s, given X,Y, IAU 2006
**     eraEors      equation of the origins, Given NPB matrix and s
**
**  References:
**
**     Capitaine, N. & Wallace, P.T., 2006, Astron.Astrophys. 450, 855
**
**     Wallace, P.T. & Capitaine, N., 2006, Astron.Astrophys. 459, 981
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double r[3][3], x, y, s, eo;


/* Classical nutation x precession x bias matrix. */
   eraPnm06a(date1, date2, r);

/* Extract CIP coordinates. */
   eraBpn2xy(r, &x, &y);

/* The CIO locator, s. */
   s = eraS06(date1, date2, x, y);

/* Solve for the EO. */
   eo = eraEors(r, s);

   return eo;

}

double eraEors(double rnpb[3][3], double s)
/*
**  - - - - - - - -
**   e r a E o r s
**  - - - - - - - -
**
**  Equation of the origins, given the classical NPB matrix and the
**  quantity s.
**
**  Given:
**     rnpb  double[3][3]  classical nutation x precession x bias matrix
**     s     double        the quantity s (the CIO locator)
**
**  Returned (function value):
**           double        the equation of the origins in radians.
**
**  Notes:
**
**  1)  The equation of the origins is the distance between the true
**      equinox and the celestial intermediate origin and, equivalently,
**      the difference between Earth rotation angle and Greenwich
**      apparent sidereal time (ERA-GST).  It comprises the precession
**      (since J2000.0) in right ascension plus the equation of the
**      equinoxes (including the small correction terms).
**
**  2)  The algorithm is from Wallace & Capitaine (2006).
**
** References:
**
**     Capitaine, N. & Wallace, P.T., 2006, Astron.Astrophys. 450, 855
**
**     Wallace, P. & Capitaine, N., 2006, Astron.Astrophys. 459, 981
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double x, ax, xs, ys, zs, p, q, eo;


/* Evaluate Wallace & Capitaine (2006) expression (16). */
   x = rnpb[2][0];
   ax = x / (1.0 + rnpb[2][2]);
   xs = 1.0 - ax * x;
   ys = -ax * rnpb[2][1];
   zs = -x;
   p = rnpb[0][0] * xs + rnpb[0][1] * ys + rnpb[0][2] * zs;
   q = rnpb[1][0] * xs + rnpb[1][1] * ys + rnpb[1][2] * zs;
   eo = ((p != 0) || (q != 0)) ? s - atan2(q, p) : s;

   return eo;

}

double eraEpb(double dj1, double dj2)
/*
**  - - - - - - -
**   e r a E p b
**  - - - - - - -
**
**  Julian Date to Besselian Epoch.
**
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
**
**  Reference:
**
**     Lieske,J.H., 1979. Astron.Astrophys.,73,282.
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
/* J2000.0 minus B1900.0 (2415019.81352) in Julian days */
   static const double D1900 = 36524.68648;

   double epb;


   epb = 1900.0 + ((dj1 - DJ00) + (dj2 + D1900)) / DTY;

   return epb;

}

void eraEpb2jd(double epb, double *djm0, double *djm)
/*
**  - - - - - - - - - -
**   e r a E p b 2 j d
**  - - - - - - - - - -
**
**  Besselian Epoch to Julian Date.
**
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
**
**  Reference:
**
**     Lieske, J.H., 1979, Astron.Astrophys. 73, 282.
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   *djm0 = 2400000.5;
   *djm  =   15019.81352 + (epb - 1900.0) * DTY;

   return;

}

double eraEpj(double dj1, double dj2)
/*
**  - - - - - - -
**   e r a E p j
**  - - - - - - -
**
**  Julian Date to Julian Epoch.
**
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
**
**  Reference:
**
**     Lieske, J.H., 1979, Astron.Astrophys. 73, 282.
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double epj;


   epj = 2000.0 + ((dj1 - DJ00) + dj2) / DJY;

   return epj;

}

void eraEpj2jd(double epj, double *djm0, double *djm)
/*
**  - - - - - - - - - -
**   e r a E p j 2 j d
**  - - - - - - - - - -
**
**  Julian Epoch to Julian Date.
**
**  Given:
**     epj      double    Julian Epoch (e.g. 1996.8D0)
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
**
**  Reference:
**
**     Lieske, J.H., 1979, Astron.Astrophys. 73, 282.
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   *djm0 = 2400000.5;
   *djm  =   51544.5 + (epj - 2000.0) * 365.25;

   return;

}

int eraEpv00(double date1, double date2,
             double pvh[2][3], double pvb[2][3])
/*
**  - - - - - - - - -
**   e r a E p v 0 0
**  - - - - - - - - -
**
**  Earth position and velocity, heliocentric and barycentric, with
**  respect to the Barycentric Celestial Reference System.
**
**  Given:
**     date1,date2  double        TDB date (Note 1)
**
**  Returned:
**     pvh          double[2][3]  heliocentric Earth position/velocity
**     pvb          double[2][3]  barycentric Earth position/velocity
**
**  Returned (function value):
**                  int           status: 0 = OK
**                                       +1 = warning: date outside
**                                            the range 1900-2100 AD
**
**  Notes:
**
**  1) The TDB date date1+date2 is a Julian Date, apportioned in any
**     convenient way between the two arguments.  For example,
**     JD(TDB)=2450123.7 could be expressed in any of these ways, among
**     others:
**
**            date1          date2
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in cases
**     where the loss of several decimal digits of resolution is
**     acceptable.  The J2000 method is best matched to the way the
**     argument is handled internally and will deliver the optimum
**     resolution.  The MJD method and the date & time methods are both
**     good compromises between resolution and convenience.  However,
**     the accuracy of the result is more likely to be limited by the
**     algorithm itself than the way the date has been expressed.
**
**     n.b. TT can be used instead of TDB in most applications.
**
**  2) On return, the arrays pvh and pvb contain the following:
**
**        pvh[0][0]  x       }
**        pvh[0][1]  y       } heliocentric position, AU
**        pvh[0][2]  z       }
**
**        pvh[1][0]  xdot    }
**        pvh[1][1]  ydot    } heliocentric velocity, AU/d
**        pvh[1][2]  zdot    }
**
**        pvb[0][0]  x       }
**        pvb[0][1]  y       } barycentric position, AU
**        pvb[0][2]  z       }
**
**        pvb[1][0]  xdot    }
**        pvb[1][1]  ydot    } barycentric velocity, AU/d
**        pvb[1][2]  zdot    }
**
**     The vectors are with respect to the Barycentric Celestial
**     Reference System.  The time unit is one day in TDB.
**
**  3) The function is a SIMPLIFIED SOLUTION from the planetary theory
**     VSOP2000 (X. Moisson, P. Bretagnon, 2001, Celes. Mechanics &
**     Dyn. Astron., 80, 3/4, 205-213) and is an adaptation of original
**     Fortran code supplied by P. Bretagnon (private comm., 2000).
**
**  4) Comparisons over the time span 1900-2100 with this simplified
**     solution and the JPL DE405 ephemeris give the following results:
**
**                                RMS    max
**           Heliocentric:
**              position error    3.7   11.2   km
**              velocity error    1.4    5.0   mm/s
**
**           Barycentric:
**              position error    4.6   13.4   km
**              velocity error    1.4    4.9   mm/s
**
**     Comparisons with the JPL DE406 ephemeris show that by 1800 and
**     2200 the position errors are approximately double their 1900-2100
**     size.  By 1500 and 2500 the deterioration is a factor of 10 and
**     by 1000 and 3000 a factor of 60.  The velocity accuracy falls off
**     at about half that rate.
**
**  5) It is permissible to use the same array for pvh and pvb, which
**     will receive the barycentric values.
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
/*
** Matrix elements for orienting the analytical model to DE405.
**
** The corresponding Euler angles are:
**
**                       d  '  "
**   1st rotation    -  23 26 21.4091 about the x-axis  (obliquity)
**   2nd rotation    +         0.0475 about the z-axis  (RA offset)
**
** These were obtained empirically, by comparisons with DE405 over
** 1900-2100.
*/
   static const double am12 =  0.000000211284,
                       am13 = -0.000000091603,
                       am21 = -0.000000230286,
                       am22 =  0.917482137087,
                       am23 = -0.397776982902,
                       am32 =  0.397776982902,
                       am33 =  0.917482137087;

/*
** ----------------------
** Ephemeris Coefficients
** ----------------------
**
** The ephemeris consists of harmonic terms for predicting (i) the Sun
** to Earth vector and (ii) the Solar-System-barycenter to Sun vector
** respectively.  The coefficients are stored in arrays which, although
** 1-demensional, contain groups of three.  Each triplet of
** coefficients is the amplitude, phase and frequency for one term in
** the model, and each array contains the number of terms called for by
** the model.
**
** There are eighteen such arrays, named as follows:
**
**     array         model      power of T      component
**
**      e0x      Sun-to-Earth        0              x
**      e0y      Sun-to-Earth        0              y
**      e0z      Sun-to-Earth        0              z
**
**      e1x      Sun-to-Earth        1              x
**      e1y      Sun-to-Earth        1              y
**      e1z      Sun-to-Earth        1              z
**
**      e2x      Sun-to-Earth        2              x
**      e2y      Sun-to-Earth        2              y
**      e2z      Sun-to-Earth        2              z
**
**      s0x      SSB-to-Sun          0              x
**      s0y      SSB-to-Sun          0              y
**      s0z      SSB-to-Sun          0              z
**
**      s1x      SSB-to-Sun          1              x
**      s1y      SSB-to-Sun          1              y
**      s1z      SSB-to-Sun          1              z
**
**      s2x      SSB-to-Sun          2              x
**      s2y      SSB-to-Sun          2              y
**      s2z      SSB-to-Sun          2              z
*/

/* Sun-to-Earth, T^0, X */
   static const double e0x[] = {
      0.9998292878132e+00, 0.1753485171504e+01, 0.6283075850446e+01,
      0.8352579567414e-02, 0.1710344404582e+01, 0.1256615170089e+02,
      0.5611445335148e-02, 0.0000000000000e+00, 0.0000000000000e+00,
      0.1046664295572e-03, 0.1667225416770e+01, 0.1884922755134e+02,
      0.3110842534677e-04, 0.6687513390251e+00, 0.8399684731857e+02,
      0.2552413503550e-04, 0.5830637358413e+00, 0.5296909721118e+00,
      0.2137207845781e-04, 0.1092330954011e+01, 0.1577343543434e+01,
      0.1680240182951e-04, 0.4955366134987e+00, 0.6279552690824e+01,
      0.1679012370795e-04, 0.6153014091901e+01, 0.6286599010068e+01,
      0.1445526946777e-04, 0.3472744100492e+01, 0.2352866153506e+01,

      0.1091038246184e-04, 0.3689845786119e+01, 0.5223693906222e+01,
      0.9344399733932e-05, 0.6073934645672e+01, 0.1203646072878e+02,
      0.8993182910652e-05, 0.3175705249069e+01, 0.1021328554739e+02,
      0.5665546034116e-05, 0.2152484672246e+01, 0.1059381944224e+01,
      0.6844146703035e-05, 0.1306964099750e+01, 0.5753384878334e+01,
      0.7346610905565e-05, 0.4354980070466e+01, 0.3981490189893e+00,
      0.6815396474414e-05, 0.2218229211267e+01, 0.4705732307012e+01,
      0.6112787253053e-05, 0.5384788425458e+01, 0.6812766822558e+01,
      0.4518120711239e-05, 0.6087604012291e+01, 0.5884926831456e+01,
      0.4521963430706e-05, 0.1279424524906e+01, 0.6256777527156e+01,

      0.4497426764085e-05, 0.5369129144266e+01, 0.6309374173736e+01,
      0.4062190566959e-05, 0.5436473303367e+00, 0.6681224869435e+01,
      0.5412193480192e-05, 0.7867838528395e+00, 0.7755226100720e+00,
      0.5469839049386e-05, 0.1461440311134e+01, 0.1414349524433e+02,
      0.5205264083477e-05, 0.4432944696116e+01, 0.7860419393880e+01,
      0.2149759935455e-05, 0.4502237496846e+01, 0.1150676975667e+02,
      0.2279109618501e-05, 0.1239441308815e+01, 0.7058598460518e+01,
      0.2259282939683e-05, 0.3272430985331e+01, 0.4694002934110e+01,
      0.2558950271319e-05, 0.2265471086404e+01, 0.1216800268190e+02,
      0.2561581447555e-05, 0.1454740653245e+01, 0.7099330490126e+00,

      0.1781441115440e-05, 0.2962068630206e+01, 0.7962980379786e+00,
      0.1612005874644e-05, 0.1473255041006e+01, 0.5486777812467e+01,
      0.1818630667105e-05, 0.3743903293447e+00, 0.6283008715021e+01,
      0.1818601377529e-05, 0.6274174354554e+01, 0.6283142985870e+01,
      0.1554475925257e-05, 0.1624110906816e+01, 0.2513230340178e+02,
      0.2090948029241e-05, 0.5852052276256e+01, 0.1179062909082e+02,
      0.2000176345460e-05, 0.4072093298513e+01, 0.1778984560711e+02,
      0.1289535917759e-05, 0.5217019331069e+01, 0.7079373888424e+01,
      0.1281135307881e-05, 0.4802054538934e+01, 0.3738761453707e+01,
      0.1518229005692e-05, 0.8691914742502e+00, 0.2132990797783e+00,

      0.9450128579027e-06, 0.4601859529950e+01, 0.1097707878456e+02,
      0.7781119494996e-06, 0.1844352816694e+01, 0.8827390247185e+01,
      0.7733407759912e-06, 0.3582790154750e+01, 0.5507553240374e+01,
      0.7350644318120e-06, 0.2695277788230e+01, 0.1589072916335e+01,
      0.6535928827023e-06, 0.3651327986142e+01, 0.1176985366291e+02,
      0.6324624183656e-06, 0.2241302375862e+01, 0.6262300422539e+01,
      0.6298565300557e-06, 0.4407122406081e+01, 0.6303851278352e+01,
      0.8587037089179e-06, 0.3024307223119e+01, 0.1672837615881e+03,
      0.8299954491035e-06, 0.6192539428237e+01, 0.3340612434717e+01,
      0.6311263503401e-06, 0.2014758795416e+01, 0.7113454667900e-02,

      0.6005646745452e-06, 0.3399500503397e+01, 0.4136910472696e+01,
      0.7917715109929e-06, 0.2493386877837e+01, 0.6069776770667e+01,
      0.7556958099685e-06, 0.4159491740143e+01, 0.6496374930224e+01,
      0.6773228244949e-06, 0.4034162934230e+01, 0.9437762937313e+01,
      0.5370708577847e-06, 0.1562219163734e+01, 0.1194447056968e+01,
      0.5710804266203e-06, 0.2662730803386e+01, 0.6282095334605e+01,
      0.5709824583726e-06, 0.3985828430833e+01, 0.6284056366286e+01,
      0.5143950896447e-06, 0.1308144688689e+01, 0.6290189305114e+01,
      0.5088010604546e-06, 0.5352817214804e+01, 0.6275962395778e+01,
      0.4960369085172e-06, 0.2644267922349e+01, 0.6127655567643e+01,

      0.4803137891183e-06, 0.4008844192080e+01, 0.6438496133249e+01,
      0.5731747768225e-06, 0.3794550174597e+01, 0.3154687086868e+01,
      0.4735947960579e-06, 0.6107118308982e+01, 0.3128388763578e+01,
      0.4808348796625e-06, 0.4771458618163e+01, 0.8018209333619e+00,
      0.4115073743137e-06, 0.3327111335159e+01, 0.8429241228195e+01,
      0.5230575889287e-06, 0.5305708551694e+01, 0.1336797263425e+02,
      0.5133977889215e-06, 0.5784230738814e+01, 0.1235285262111e+02,
      0.5065815825327e-06, 0.2052064793679e+01, 0.1185621865188e+02,
      0.4339831593868e-06, 0.3644994195830e+01, 0.1726015463500e+02,
      0.3952928638953e-06, 0.4930376436758e+01, 0.5481254917084e+01,

      0.4898498111942e-06, 0.4542084219731e+00, 0.9225539266174e+01,
      0.4757490209328e-06, 0.3161126388878e+01, 0.5856477690889e+01,
      0.4727701669749e-06, 0.6214993845446e+00, 0.2544314396739e+01,
      0.3800966681863e-06, 0.3040132339297e+01, 0.4265981595566e+00,
      0.3257301077939e-06, 0.8064977360087e+00, 0.3930209696940e+01,
      0.3255810528674e-06, 0.1974147981034e+01, 0.2146165377750e+01,
      0.3252029748187e-06, 0.2845924913135e+01, 0.4164311961999e+01,
      0.3255505635308e-06, 0.3017900824120e+01, 0.5088628793478e+01,
      0.2801345211990e-06, 0.6109717793179e+01, 0.1256967486051e+02,
      0.3688987740970e-06, 0.2911550235289e+01, 0.1807370494127e+02,

      0.2475153429458e-06, 0.2179146025856e+01, 0.2629832328990e-01,
      0.3033457749150e-06, 0.1994161050744e+01, 0.4535059491685e+01,
      0.2186743763110e-06, 0.5125687237936e+01, 0.1137170464392e+02,
      0.2764777032774e-06, 0.4822646860252e+00, 0.1256262854127e+02,
      0.2199028768592e-06, 0.4637633293831e+01, 0.1255903824622e+02,
      0.2046482824760e-06, 0.1467038733093e+01, 0.7084896783808e+01,
      0.2611209147507e-06, 0.3044718783485e+00, 0.7143069561767e+02,
      0.2286079656818e-06, 0.4764220356805e+01, 0.8031092209206e+01,
      0.1855071202587e-06, 0.3383637774428e+01, 0.1748016358760e+01,
      0.2324669506784e-06, 0.6189088449251e+01, 0.1831953657923e+02,

      0.1709528015688e-06, 0.5874966729774e+00, 0.4933208510675e+01,
      0.2168156875828e-06, 0.4302994009132e+01, 0.1044738781244e+02,
      0.2106675556535e-06, 0.3800475419891e+01, 0.7477522907414e+01,
      0.1430213830465e-06, 0.1294660846502e+01, 0.2942463415728e+01,
      0.1388396901944e-06, 0.4594797202114e+01, 0.8635942003952e+01,
      0.1922258844190e-06, 0.4943044543591e+00, 0.1729818233119e+02,
      0.1888460058292e-06, 0.2426943912028e+01, 0.1561374759853e+03,
      0.1789449386107e-06, 0.1582973303499e+00, 0.1592596075957e+01,
      0.1360803685374e-06, 0.5197240440504e+01, 0.1309584267300e+02,
      0.1504038014709e-06, 0.3120360916217e+01, 0.1649636139783e+02,

      0.1382769533389e-06, 0.6164702888205e+01, 0.7632943190217e+01,
      0.1438059769079e-06, 0.1437423770979e+01, 0.2042657109477e+02,
      0.1326303260037e-06, 0.3609688799679e+01, 0.1213955354133e+02,
      0.1159244950540e-06, 0.5463018167225e+01, 0.5331357529664e+01,
      0.1433118149136e-06, 0.6028909912097e+01, 0.7342457794669e+01,
      0.1234623148594e-06, 0.3109645574997e+01, 0.6279485555400e+01,
      0.1233949875344e-06, 0.3539359332866e+01, 0.6286666145492e+01,
      0.9927196061299e-07, 0.1259321569772e+01, 0.7234794171227e+01,
      0.1242302191316e-06, 0.1065949392609e+01, 0.1511046609763e+02,
      0.1098402195201e-06, 0.2192508743837e+01, 0.1098880815746e+02,

      0.1158191395315e-06, 0.4054411278650e+01, 0.5729506548653e+01,
      0.9048475596241e-07, 0.5429764748518e+01, 0.9623688285163e+01,
      0.8889853269023e-07, 0.5046586206575e+01, 0.6148010737701e+01,
      0.1048694242164e-06, 0.2628858030806e+01, 0.6836645152238e+01,
      0.1112308378646e-06, 0.4177292719907e+01, 0.1572083878776e+02,
      0.8631729709901e-07, 0.1601345232557e+01, 0.6418140963190e+01,
      0.8527816951664e-07, 0.2463888997513e+01, 0.1471231707864e+02,
      0.7892139456991e-07, 0.3154022088718e+01, 0.2118763888447e+01,
      0.1051782905236e-06, 0.4795035816088e+01, 0.1349867339771e+01,
      0.1048219943164e-06, 0.2952983395230e+01, 0.5999216516294e+01,

      0.7435760775143e-07, 0.5420547991464e+01, 0.6040347114260e+01,
      0.9869574106949e-07, 0.3695646753667e+01, 0.6566935184597e+01,
      0.9156886364226e-07, 0.3922675306609e+01, 0.5643178611111e+01,
      0.7006834356188e-07, 0.1233968624861e+01, 0.6525804586632e+01,
      0.9806170182601e-07, 0.1919542280684e+01, 0.2122839202813e+02,
      0.9052289673607e-07, 0.4615902724369e+01, 0.4690479774488e+01,
      0.7554200867893e-07, 0.1236863719072e+01, 0.1253985337760e+02,
      0.8215741286498e-07, 0.3286800101559e+00, 0.1097355562493e+02,
      0.7185178575397e-07, 0.5880942158367e+01, 0.6245048154254e+01,
      0.7130726476180e-07, 0.7674871987661e+00, 0.6321103546637e+01,

      0.6650894461162e-07, 0.6987129150116e+00, 0.5327476111629e+01,
      0.7396888823688e-07, 0.3576824794443e+01, 0.5368044267797e+00,
      0.7420588884775e-07, 0.5033615245369e+01, 0.2354323048545e+02,
      0.6141181642908e-07, 0.9449927045673e+00, 0.1296430071988e+02,
      0.6373557924058e-07, 0.6206342280341e+01, 0.9517183207817e+00,
      0.6359474329261e-07, 0.5036079095757e+01, 0.1990745094947e+01,
      0.5740173582646e-07, 0.6105106371350e+01, 0.9555997388169e+00,
      0.7019864084602e-07, 0.7237747359018e+00, 0.5225775174439e+00,
      0.6398054487042e-07, 0.3976367969666e+01, 0.2407292145756e+02,
      0.7797092650498e-07, 0.4305423910623e+01, 0.2200391463820e+02,

      0.6466760000900e-07, 0.3500136825200e+01, 0.5230807360890e+01,
      0.7529417043890e-07, 0.3514779246100e+01, 0.1842262939178e+02,
      0.6924571140892e-07, 0.2743457928679e+01, 0.1554202828031e+00,
      0.6220798650222e-07, 0.2242598118209e+01, 0.1845107853235e+02,
      0.5870209391853e-07, 0.2332832707527e+01, 0.6398972393349e+00,
      0.6263953473888e-07, 0.2191105358956e+01, 0.6277552955062e+01,
      0.6257781390012e-07, 0.4457559396698e+01, 0.6288598745829e+01,
      0.5697304945123e-07, 0.3499234761404e+01, 0.1551045220144e+01,
      0.6335438746791e-07, 0.6441691079251e+00, 0.5216580451554e+01,
      0.6377258441152e-07, 0.2252599151092e+01, 0.5650292065779e+01,

      0.6484841818165e-07, 0.1992812417646e+01, 0.1030928125552e+00,
      0.4735551485250e-07, 0.3744672082942e+01, 0.1431416805965e+02,
      0.4628595996170e-07, 0.1334226211745e+01, 0.5535693017924e+00,
      0.6258152336933e-07, 0.4395836159154e+01, 0.2608790314060e+02,
      0.6196171366594e-07, 0.2587043007997e+01, 0.8467247584405e+02,
      0.6159556952126e-07, 0.4782499769128e+01, 0.2394243902548e+03,
      0.4987741172394e-07, 0.7312257619924e+00, 0.7771377146812e+02,
      0.5459280703142e-07, 0.3001376372532e+01, 0.6179983037890e+01,
      0.4863461189999e-07, 0.3767222128541e+01, 0.9027992316901e+02,
      0.5349912093158e-07, 0.3663594450273e+01, 0.6386168663001e+01,

      0.5673725607806e-07, 0.4331187919049e+01, 0.6915859635113e+01,
      0.4745485060512e-07, 0.5816195745518e+01, 0.6282970628506e+01,
      0.4745379005326e-07, 0.8323672435672e+00, 0.6283181072386e+01,
      0.4049002796321e-07, 0.3785023976293e+01, 0.6254626709878e+01,
      0.4247084014515e-07, 0.2378220728783e+01, 0.7875671926403e+01,
      0.4026912363055e-07, 0.2864103423269e+01, 0.6311524991013e+01,
      0.4062935011774e-07, 0.2415408595975e+01, 0.3634620989887e+01,
      0.5347771048509e-07, 0.3343479309801e+01, 0.2515860172507e+02,
      0.4829494136505e-07, 0.2821742398262e+01, 0.5760498333002e+01,
      0.4342554404599e-07, 0.5624662458712e+01, 0.7238675589263e+01,

      0.4021599184361e-07, 0.5557250275009e+00, 0.1101510648075e+02,
      0.4104900474558e-07, 0.3296691780005e+01, 0.6709674010002e+01,
      0.4376532905131e-07, 0.3814443999443e+01, 0.6805653367890e+01,
      0.3314590480650e-07, 0.3560229189250e+01, 0.1259245002418e+02,
      0.3232421839643e-07, 0.5185389180568e+01, 0.1066495398892e+01,
      0.3541176318876e-07, 0.3921381909679e+01, 0.9917696840332e+01,
      0.3689831242681e-07, 0.4190658955386e+01, 0.1192625446156e+02,
      0.3890605376774e-07, 0.5546023371097e+01, 0.7478166569050e-01,
      0.3038559339780e-07, 0.6231032794494e+01, 0.1256621883632e+02,
      0.3137083969782e-07, 0.6207063419190e+01, 0.4292330755499e+01,

      0.4024004081854e-07, 0.1195257375713e+01, 0.1334167431096e+02,
      0.3300234879283e-07, 0.1804694240998e+01, 0.1057540660594e+02,
      0.3635399155575e-07, 0.5597811343500e+01, 0.6208294184755e+01,
      0.3032668691356e-07, 0.3191059366530e+01, 0.1805292951336e+02,
      0.2809652069058e-07, 0.4094348032570e+01, 0.3523159621801e-02,
      0.3696955383823e-07, 0.5219282738794e+01, 0.5966683958112e+01,
      0.3562894142503e-07, 0.1037247544554e+01, 0.6357857516136e+01,
      0.3510598524148e-07, 0.1430020816116e+01, 0.6599467742779e+01,
      0.3617736142953e-07, 0.3002911403677e+01, 0.6019991944201e+01,
      0.2624524910730e-07, 0.2437046757292e+01, 0.6702560555334e+01,

      0.2535824204490e-07, 0.1581594689647e+01, 0.3141537925223e+02,
      0.3519787226257e-07, 0.5379863121521e+01, 0.2505706758577e+03,
      0.2578406709982e-07, 0.4904222639329e+01, 0.1673046366289e+02,
      0.3423887981473e-07, 0.3646448997315e+01, 0.6546159756691e+01,
      0.2776083886467e-07, 0.3307829300144e+01, 0.1272157198369e+02,
      0.3379592818379e-07, 0.1747541251125e+01, 0.1494531617769e+02,
      0.3050255426284e-07, 0.1784689432607e-01, 0.4732030630302e+01,
      0.2652378350236e-07, 0.4420055276260e+01, 0.5863591145557e+01,
      0.2374498173768e-07, 0.3629773929208e+01, 0.2388894113936e+01,
      0.2716451255140e-07, 0.3079623706780e+01, 0.1202934727411e+02,

      0.3038583699229e-07, 0.3312487903507e+00, 0.1256608456547e+02,
      0.2220681228760e-07, 0.5265520401774e+01, 0.1336244973887e+02,
      0.3044156540912e-07, 0.4766664081250e+01, 0.2908881142201e+02,
      0.2731859923561e-07, 0.5069146530691e+01, 0.1391601904066e+02,
      0.2285603018171e-07, 0.5954935112271e+01, 0.6076890225335e+01,
      0.2025006454555e-07, 0.4061789589267e+01, 0.4701116388778e+01,
      0.2012597519804e-07, 0.2485047705241e+01, 0.6262720680387e+01,
      0.2003406962258e-07, 0.4163779209320e+01, 0.6303431020504e+01,
      0.2207863441371e-07, 0.6923839133828e+00, 0.6489261475556e+01,
      0.2481374305624e-07, 0.5944173595676e+01, 0.1204357418345e+02,

      0.2130923288870e-07, 0.4641013671967e+01, 0.5746271423666e+01,
      0.2446370543391e-07, 0.6125796518757e+01, 0.1495633313810e+00,
      0.1932492759052e-07, 0.2234572324504e+00, 0.1352175143971e+02,
      0.2600122568049e-07, 0.4281012405440e+01, 0.4590910121555e+01,
      0.2431754047488e-07, 0.1429943874870e+00, 0.1162474756779e+01,
      0.1875902869209e-07, 0.9781803816948e+00, 0.6279194432410e+01,
      0.1874381139426e-07, 0.5670368130173e+01, 0.6286957268481e+01,
      0.2156696047173e-07, 0.2008985006833e+01, 0.1813929450232e+02,
      0.1965076182484e-07, 0.2566186202453e+00, 0.4686889479442e+01,
      0.2334816372359e-07, 0.4408121891493e+01, 0.1002183730415e+02,

      0.1869937408802e-07, 0.5272745038656e+01, 0.2427287361862e+00,
      0.2436236460883e-07, 0.4407720479029e+01, 0.9514313292143e+02,
      0.1761365216611e-07, 0.1943892315074e+00, 0.1351787002167e+02,
      0.2156289480503e-07, 0.1418570924545e+01, 0.6037244212485e+01,
      0.2164748979255e-07, 0.4724603439430e+01, 0.2301353951334e+02,
      0.2222286670853e-07, 0.2400266874598e+01, 0.1266924451345e+02,
      0.2070901414929e-07, 0.5230348028732e+01, 0.6528907488406e+01,
      0.1792745177020e-07, 0.2099190328945e+01, 0.6819880277225e+01,
      0.1841802068445e-07, 0.3467527844848e+00, 0.6514761976723e+02,
      0.1578401631718e-07, 0.7098642356340e+00, 0.2077542790660e-01,

      0.1561690152531e-07, 0.5943349620372e+01, 0.6272439236156e+01,
      0.1558591045463e-07, 0.7040653478980e+00, 0.6293712464735e+01,
      0.1737356469576e-07, 0.4487064760345e+01, 0.1765478049437e+02,
      0.1434755619991e-07, 0.2993391570995e+01, 0.1102062672231e+00,
      0.1482187806654e-07, 0.2278049198251e+01, 0.1052268489556e+01,
      0.1424812827089e-07, 0.1682114725827e+01, 0.1311972100268e+02,
      0.1380282448623e-07, 0.3262668602579e+01, 0.1017725758696e+02,
      0.1811481244566e-07, 0.3187771221777e+01, 0.1887552587463e+02,
      0.1504446185696e-07, 0.5650162308647e+01, 0.7626583626240e-01,
      0.1740776154137e-07, 0.5487068607507e+01, 0.1965104848470e+02,

      0.1374339536251e-07, 0.5745688172201e+01, 0.6016468784579e+01,
      0.1761377477704e-07, 0.5748060203659e+01, 0.2593412433514e+02,
      0.1535138225795e-07, 0.6226848505790e+01, 0.9411464614024e+01,
      0.1788140543676e-07, 0.6189318878563e+01, 0.3301902111895e+02,
      0.1375002807996e-07, 0.5371812884394e+01, 0.6327837846670e+00,
      0.1242115758632e-07, 0.1471687569712e+01, 0.3894181736510e+01,
      0.1450977333938e-07, 0.4143836662127e+01, 0.1277945078067e+02,
      0.1297579575023e-07, 0.9003477661957e+00, 0.6549682916313e+01,
      0.1462667934821e-07, 0.5760505536428e+01, 0.1863592847156e+02,
      0.1381774374799e-07, 0.1085471729463e+01, 0.2379164476796e+01,

      0.1682333169307e-07, 0.5409870870133e+01, 0.1620077269078e+02,
      0.1190812918837e-07, 0.1397205174601e+01, 0.1149965630200e+02,
      0.1221434762106e-07, 0.9001804809095e+00, 0.1257326515556e+02,
      0.1549934644860e-07, 0.4262528275544e+01, 0.1820933031200e+02,
      0.1252138953050e-07, 0.1411642012027e+01, 0.6993008899458e+01,
      0.1237078905387e-07, 0.2844472403615e+01, 0.2435678079171e+02,
      0.1446953389615e-07, 0.5295835522223e+01, 0.3813291813120e-01,
      0.1388446457170e-07, 0.4969428135497e+01, 0.2458316379602e+00,
      0.1019339179228e-07, 0.2491369561806e+01, 0.6112403035119e+01,
      0.1258880815343e-07, 0.4679426248976e+01, 0.5429879531333e+01,

      0.1297768238261e-07, 0.1074509953328e+01, 0.1249137003520e+02,
      0.9913505718094e-08, 0.4735097918224e+01, 0.6247047890016e+01,
      0.9830453155969e-08, 0.4158649187338e+01, 0.6453748665772e+01,
      0.1192615865309e-07, 0.3438208613699e+01, 0.6290122169689e+01,
      0.9835874798277e-08, 0.1913300781229e+01, 0.6319103810876e+01,
      0.9639087569277e-08, 0.9487683644125e+00, 0.8273820945392e+01,
      0.1175716107001e-07, 0.3228141664287e+01, 0.6276029531202e+01,
      0.1018926508678e-07, 0.2216607854300e+01, 0.1254537627298e+02,
      0.9500087869225e-08, 0.2625116459733e+01, 0.1256517118505e+02,
      0.9664192916575e-08, 0.5860562449214e+01, 0.6259197520765e+01,

      0.9612858712203e-08, 0.7885682917381e+00, 0.6306954180126e+01,
      0.1117645675413e-07, 0.3932148831189e+01, 0.1779695906178e+02,
      0.1158864052160e-07, 0.9995605521691e+00, 0.1778273215245e+02,
      0.9021043467028e-08, 0.5263769742673e+01, 0.6172869583223e+01,
      0.8836134773563e-08, 0.1496843220365e+01, 0.1692165728891e+01,
      0.1045872200691e-07, 0.7009039517214e+00, 0.2204125344462e+00,
      0.1211463487798e-07, 0.4041544938511e+01, 0.8257698122054e+02,
      0.8541990804094e-08, 0.1447586692316e+01, 0.6393282117669e+01,
      0.1038720703636e-07, 0.4594249718112e+00, 0.1550861511662e+02,
      0.1126722351445e-07, 0.3925550579036e+01, 0.2061856251104e+00,

      0.8697373859631e-08, 0.4411341856037e+01, 0.9491756770005e+00,
      0.8869380028441e-08, 0.2402659724813e+01, 0.3903911373650e+01,
      0.9247014693258e-08, 0.1401579743423e+01, 0.6267823317922e+01,
      0.9205062930950e-08, 0.5245978000814e+01, 0.6298328382969e+01,
      0.8000745038049e-08, 0.3590803356945e+01, 0.2648454860559e+01,
      0.9168973650819e-08, 0.2470150501679e+01, 0.1498544001348e+03,
      0.1075444949238e-07, 0.1328606161230e+01, 0.3694923081589e+02,
      0.7817298525817e-08, 0.6162256225998e+01, 0.4804209201333e+01,
      0.9541469226356e-08, 0.3942568967039e+01, 0.1256713221673e+02,
      0.9821910122027e-08, 0.2360246287233e+00, 0.1140367694411e+02,

      0.9897822023777e-08, 0.4619805634280e+01, 0.2280573557157e+02,
      0.7737289283765e-08, 0.3784727847451e+01, 0.7834121070590e+01,
      0.9260204034710e-08, 0.2223352487601e+01, 0.2787043132925e+01,
      0.7320252888486e-08, 0.1288694636874e+01, 0.6282655592598e+01,
      0.7319785780946e-08, 0.5359869567774e+01, 0.6283496108294e+01,
      0.7147219933778e-08, 0.5516616675856e+01, 0.1725663147538e+02,
      0.7946502829878e-08, 0.2630459984567e+01, 0.1241073141809e+02,
      0.9001711808932e-08, 0.2849815827227e+01, 0.6281591679874e+01,
      0.8994041507257e-08, 0.3795244450750e+01, 0.6284560021018e+01,
      0.8298582787358e-08, 0.5236413127363e+00, 0.1241658836951e+02,

      0.8526596520710e-08, 0.4794605424426e+01, 0.1098419223922e+02,
      0.8209822103197e-08, 0.1578752370328e+01, 0.1096996532989e+02,
      0.6357049861094e-08, 0.5708926113761e+01, 0.1596186371003e+01,
      0.7370473179049e-08, 0.3842402530241e+01, 0.4061219149443e+01,
      0.7232154664726e-08, 0.3067548981535e+01, 0.1610006857377e+03,
      0.6328765494903e-08, 0.1313930030069e+01, 0.1193336791622e+02,
      0.8030064908595e-08, 0.3488500408886e+01, 0.8460828644453e+00,
      0.6275464259232e-08, 0.1532061626198e+01, 0.8531963191132e+00,
      0.7051897446325e-08, 0.3285859929993e+01, 0.5849364236221e+01,
      0.6161593705428e-08, 0.1477341999464e+01, 0.5573142801433e+01,

      0.7754683957278e-08, 0.1586118663096e+01, 0.8662240327241e+01,
      0.5889928990701e-08, 0.1304887868803e+01, 0.1232342296471e+02,
      0.5705756047075e-08, 0.4555333589350e+01, 0.1258692712880e+02,
      0.5964178808332e-08, 0.3001762842062e+01, 0.5333900173445e+01,
      0.6712446027467e-08, 0.4886780007595e+01, 0.1171295538178e+02,
      0.5941809275464e-08, 0.4701509603824e+01, 0.9779108567966e+01,
      0.5466993627395e-08, 0.4588357817278e+01, 0.1884211409667e+02,
      0.6340512090980e-08, 0.1164543038893e+01, 0.5217580628120e+02,
      0.6325505710045e-08, 0.3919171259645e+01, 0.1041998632314e+02,
      0.6164789509685e-08, 0.2143828253542e+01, 0.6151533897323e+01,

      0.5263330812430e-08, 0.6066564434241e+01, 0.1885275071096e+02,
      0.5597087780221e-08, 0.2926316429472e+01, 0.4337116142245e+00,
      0.5396556236817e-08, 0.3244303591505e+01, 0.6286362197481e+01,
      0.5396615148223e-08, 0.3404304703662e+01, 0.6279789503410e+01,
      0.7091832443341e-08, 0.8532377803192e+00, 0.4907302013889e+01,
      0.6572352589782e-08, 0.4901966774419e+01, 0.1176433076753e+02,
      0.5960236060795e-08, 0.1874672315797e+01, 0.1422690933580e-01,
      0.5125480043511e-08, 0.3735726064334e+01, 0.1245594543367e+02,
      0.5928241866410e-08, 0.4502033899935e+01, 0.6414617803568e+01,
      0.5249600357424e-08, 0.4372334799878e+01, 0.1151388321134e+02,

      0.6059171276087e-08, 0.2581617302908e+01, 0.6062663316000e+01,
      0.5295235081662e-08, 0.2974811513158e+01, 0.3496032717521e+01,
      0.5820561875933e-08, 0.1796073748244e+00, 0.2838593341516e+00,
      0.4754696606440e-08, 0.1981998136973e+01, 0.3104930017775e+01,
      0.6385053548955e-08, 0.2559174171605e+00, 0.6133512519065e+01,
      0.6589828273941e-08, 0.2750967106776e+01, 0.4087944051283e+02,
      0.5383376567189e-08, 0.6325947523578e+00, 0.2248384854122e+02,
      0.5928941683538e-08, 0.1672304519067e+01, 0.1581959461667e+01,
      0.4816060709794e-08, 0.3512566172575e+01, 0.9388005868221e+01,
      0.6003381586512e-08, 0.5610932219189e+01, 0.5326786718777e+01,

      0.5504225393105e-08, 0.4037501131256e+01, 0.6503488384892e+01,
      0.5353772620129e-08, 0.6122774968240e+01, 0.1735668374386e+03,
      0.5786253768544e-08, 0.5527984999515e+01, 0.1350651127443e+00,
      0.5065706702002e-08, 0.9980765573624e+00, 0.1248988586463e+02,
      0.5972838885276e-08, 0.6044489493203e+01, 0.2673594526851e+02,
      0.5323585877961e-08, 0.3924265998147e+01, 0.4171425416666e+01,
      0.5210772682858e-08, 0.6220111376901e+01, 0.2460261242967e+02,
      0.4726549040535e-08, 0.3716043206862e+01, 0.7232251527446e+01,
      0.6029425105059e-08, 0.8548704071116e+00, 0.3227113045244e+03,
      0.4481542826513e-08, 0.1426925072829e+01, 0.5547199253223e+01,

      0.5836024505068e-08, 0.7135651752625e-01, 0.7285056171570e+02,
      0.4137046613272e-08, 0.5330767643283e+01, 0.1087398597200e+02,
      0.5171977473924e-08, 0.4494262335353e+00, 0.1884570439172e+02,
      0.5694429833732e-08, 0.2952369582215e+01, 0.9723862754494e+02,
      0.4009158925298e-08, 0.3500003416535e+01, 0.6244942932314e+01,
      0.4784939596873e-08, 0.6196709413181e+01, 0.2929661536378e+02,
      0.3983725022610e-08, 0.5103690031897e+01, 0.4274518229222e+01,
      0.3870535232462e-08, 0.3187569587401e+01, 0.6321208768577e+01,
      0.5140501213951e-08, 0.1668924357457e+01, 0.1232032006293e+02,
      0.3849034819355e-08, 0.4445722510309e+01, 0.1726726808967e+02,

      0.4002383075060e-08, 0.5226224152423e+01, 0.7018952447668e+01,
      0.3890719543549e-08, 0.4371166550274e+01, 0.1491901785440e+02,
      0.4887084607881e-08, 0.5973556689693e+01, 0.1478866649112e+01,
      0.3739939287592e-08, 0.2089084714600e+01, 0.6922973089781e+01,
      0.5031925918209e-08, 0.4658371936827e+01, 0.1715706182245e+02,
      0.4387748764954e-08, 0.4825580552819e+01, 0.2331413144044e+03,
      0.4147398098865e-08, 0.3739003524998e+01, 0.1376059875786e+02,
      0.3719089993586e-08, 0.1148941386536e+01, 0.6297302759782e+01,
      0.3934238461056e-08, 0.1559893008343e+01, 0.7872148766781e+01,
      0.3672471375622e-08, 0.5516145383612e+01, 0.6268848941110e+01,

      0.3768911277583e-08, 0.6116053700563e+01, 0.4157198507331e+01,
      0.4033388417295e-08, 0.5076821746017e+01, 0.1567108171867e+02,
      0.3764194617832e-08, 0.8164676232075e+00, 0.3185192151914e+01,
      0.4840628226284e-08, 0.1360479453671e+01, 0.1252801878276e+02,
      0.4949443923785e-08, 0.2725622229926e+01, 0.1617106187867e+03,
      0.4117393089971e-08, 0.6054459628492e+00, 0.5642198095270e+01,
      0.3925754020428e-08, 0.8570462135210e+00, 0.2139354194808e+02,
      0.3630551757923e-08, 0.3552067338279e+01, 0.6294805223347e+01,
      0.3627274802357e-08, 0.3096565085313e+01, 0.6271346477544e+01,
      0.3806143885093e-08, 0.6367751709777e+00, 0.1725304118033e+02,

      0.4433254641565e-08, 0.4848461503937e+01, 0.7445550607224e+01,
      0.3712319846576e-08, 0.1331950643655e+01, 0.4194847048887e+00,
      0.3849847534783e-08, 0.4958368297746e+00, 0.9562891316684e+00,
      0.3483955430165e-08, 0.2237215515707e+01, 0.1161697602389e+02,
      0.3961912730982e-08, 0.3332402188575e+01, 0.2277943724828e+02,
      0.3419978244481e-08, 0.5785600576016e+01, 0.1362553364512e+02,
      0.3329417758177e-08, 0.9812676559709e-01, 0.1685848245639e+02,
      0.4207206893193e-08, 0.9494780468236e+00, 0.2986433403208e+02,
      0.3268548976410e-08, 0.1739332095686e+00, 0.5749861718712e+01,
      0.3321880082685e-08, 0.1423354800666e+01, 0.6279143387820e+01,

      0.4503173010852e-08, 0.2314972675293e+00, 0.1385561574497e+01,
      0.4316599090954e-08, 0.1012646782616e+00, 0.4176041334900e+01,
      0.3283493323850e-08, 0.5233306881265e+01, 0.6287008313071e+01,
      0.3164033542343e-08, 0.4005597257511e+01, 0.2099539292909e+02,
      0.4159720956725e-08, 0.5365676242020e+01, 0.5905702259363e+01,
      0.3565176892217e-08, 0.4284440620612e+01, 0.3932462625300e-02,
      0.3514440950221e-08, 0.4270562636575e+01, 0.7335344340001e+01,
      0.3540596871909e-08, 0.5953553201060e+01, 0.1234573916645e+02,
      0.2960769905118e-08, 0.1115180417718e+01, 0.2670964694522e+02,
      0.2962213739684e-08, 0.3863811918186e+01, 0.6408777551755e+00,

      0.3883556700251e-08, 0.1268617928302e+01, 0.6660449441528e+01,
      0.2919225516346e-08, 0.4908605223265e+01, 0.1375773836557e+01,
      0.3115158863370e-08, 0.3744519976885e+01, 0.3802769619140e-01,
      0.4099438144212e-08, 0.4173244670532e+01, 0.4480965020977e+02,
      0.2899531858964e-08, 0.5910601428850e+01, 0.2059724391010e+02,
      0.3289733429855e-08, 0.2488050078239e+01, 0.1081813534213e+02,
      0.3933075612875e-08, 0.1122363652883e+01, 0.3773735910827e+00,
      0.3021403764467e-08, 0.4951973724904e+01, 0.2982630633589e+02,
      0.2798598949757e-08, 0.5117057845513e+01, 0.1937891852345e+02,
      0.3397421302707e-08, 0.6104159180476e+01, 0.6923953605621e+01,

      0.3720398002179e-08, 0.1184933429829e+01, 0.3066615496545e+02,
      0.3598484186267e-08, 0.3505282086105e+01, 0.6147450479709e+01,
      0.3694594027310e-08, 0.2286651088141e+01, 0.2636725487657e+01,
      0.2680444152969e-08, 0.1871816775482e+00, 0.6816289982179e+01,
      0.3497574865641e-08, 0.3143251755431e+01, 0.6418701221183e+01,
      0.3130274129494e-08, 0.2462167316018e+01, 0.1235996607578e+02,
      0.3241119069551e-08, 0.4256374004686e+01, 0.1652265972112e+02,
      0.2601960842061e-08, 0.4970362941425e+01, 0.1045450126711e+02,
      0.2690601527504e-08, 0.2372657824898e+01, 0.3163918923335e+00,
      0.2908688152664e-08, 0.4232652627721e+01, 0.2828699048865e+02,

      0.3120456131875e-08, 0.3925747001137e+00, 0.2195415756911e+02,
      0.3148855423384e-08, 0.3093478330445e+01, 0.1172006883645e+02,
      0.3051044261017e-08, 0.5560948248212e+01, 0.6055599646783e+01,
      0.2826006876660e-08, 0.5072790310072e+01, 0.5120601093667e+01,
      0.3100034191711e-08, 0.4998530231096e+01, 0.1799603123222e+02,
      0.2398771640101e-08, 0.2561739802176e+01, 0.6255674361143e+01,
      0.2384002842728e-08, 0.4087420284111e+01, 0.6310477339748e+01,
      0.2842146517568e-08, 0.2515048217955e+01, 0.5469525544182e+01,
      0.2847674371340e-08, 0.5235326497443e+01, 0.1034429499989e+02,
      0.2903722140764e-08, 0.1088200795797e+01, 0.6510552054109e+01,

      0.3187610710605e-08, 0.4710624424816e+01, 0.1693792562116e+03,
      0.3048869992813e-08, 0.2857975896445e+00, 0.8390110365991e+01,
      0.2860216950984e-08, 0.2241619020815e+01, 0.2243449970715e+00,
      0.2701117683113e-08, 0.6651573305272e-01, 0.6129297044991e+01,
      0.2509891590152e-08, 0.1285135324585e+01, 0.1044027435778e+02,
      0.2623200252223e-08, 0.2981229834530e+00, 0.6436854655901e+01,
      0.2622541669202e-08, 0.6122470726189e+01, 0.9380959548977e+01,
      0.2818435667099e-08, 0.4251087148947e+01, 0.5934151399930e+01,
      0.2365196797465e-08, 0.3465070460790e+01, 0.2470570524223e+02,
      0.2358704646143e-08, 0.5791603815350e+01, 0.8671969964381e+01,

      0.2388299481390e-08, 0.4142483772941e+01, 0.7096626156709e+01,
      0.1996041217224e-08, 0.2101901889496e+01, 0.1727188400790e+02,
      0.2687593060336e-08, 0.1526689456959e+01, 0.7075506709219e+02,
      0.2618913670810e-08, 0.2397684236095e+01, 0.6632000300961e+01,
      0.2571523050364e-08, 0.5751929456787e+00, 0.6206810014183e+01,
      0.2582135006946e-08, 0.5595464352926e+01, 0.4873985990671e+02,
      0.2372530190361e-08, 0.5092689490655e+01, 0.1590676413561e+02,
      0.2357178484712e-08, 0.4444363527851e+01, 0.3097883698531e+01,
      0.2451590394723e-08, 0.3108251687661e+01, 0.6612329252343e+00,
      0.2370045949608e-08, 0.2608133861079e+01, 0.3459636466239e+02,

      0.2268997267358e-08, 0.3639717753384e+01, 0.2844914056730e-01,
      0.1731432137906e-08, 0.1741898445707e+00, 0.2019909489111e+02,
      0.1629869741622e-08, 0.3902225646724e+01, 0.3035599730800e+02,
      0.2206215801974e-08, 0.4971131250731e+01, 0.6281667977667e+01,
      0.2205469554680e-08, 0.1677462357110e+01, 0.6284483723224e+01,
      0.2148792362509e-08, 0.4236259604006e+01, 0.1980482729015e+02,
      0.1873733657847e-08, 0.5926814998687e+01, 0.2876692439167e+02,
      0.2026573758959e-08, 0.4349643351962e+01, 0.2449240616245e+02,
      0.1807770325110e-08, 0.5700940482701e+01, 0.2045286941806e+02,
      0.1881174408581e-08, 0.6601286363430e+00, 0.2358125818164e+02,

      0.1368023671690e-08, 0.2211098592752e+01, 0.2473415438279e+02,
      0.1720017916280e-08, 0.4942488551129e+01, 0.1679593901136e+03,
      0.1702427665131e-08, 0.1452233856386e+01, 0.3338575901272e+03,
      0.1414032510054e-08, 0.5525357721439e+01, 0.1624205518357e+03,
      0.1652626045364e-08, 0.4108794283624e+01, 0.8956999012000e+02,
      0.1642957769686e-08, 0.7344335209984e+00, 0.5267006960365e+02,
      0.1614952403624e-08, 0.3541213951363e+01, 0.3332657872986e+02,
      0.1535988291188e-08, 0.4031094072151e+01, 0.3852657435933e+02,
      0.1593193738177e-08, 0.4185136203609e+01, 0.2282781046519e+03,
      0.1074569126382e-08, 0.1720485636868e+01, 0.8397383534231e+02,

      0.1074408214509e-08, 0.2758613420318e+01, 0.8401985929482e+02,
      0.9700199670465e-09, 0.4216686842097e+01, 0.7826370942180e+02,
      0.1258433517061e-08, 0.2575068876639e+00, 0.3115650189215e+03,
      0.1240303229539e-08, 0.4800844956756e+00, 0.1784300471910e+03,
      0.9018345948127e-09, 0.3896756361552e+00, 0.5886454391678e+02,
      0.1135301432805e-08, 0.3700805023550e+00, 0.7842370451713e+02,
      0.9215887951370e-09, 0.4364579276638e+01, 0.1014262087719e+03,
      0.1055401054147e-08, 0.2156564222111e+01, 0.5660027930059e+02,
      0.1008725979831e-08, 0.5454015785234e+01, 0.4245678405627e+02,
      0.7217398104321e-09, 0.1597772562175e+01, 0.2457074661053e+03,

      0.6912033134447e-09, 0.5824090621461e+01, 0.1679936946371e+03,
      0.6833881523549e-09, 0.3578778482835e+01, 0.6053048899753e+02,
      0.4887304205142e-09, 0.3724362812423e+01, 0.9656299901946e+02,
      0.5173709754788e-09, 0.5422427507933e+01, 0.2442876000072e+03,
      0.4671353097145e-09, 0.2396106924439e+01, 0.1435713242844e+03,
      0.5652608439480e-09, 0.2804028838685e+01, 0.8365903305582e+02,
      0.5604061331253e-09, 0.1638816006247e+01, 0.8433466158131e+02,
      0.4712723365400e-09, 0.8979003224474e+00, 0.3164282286739e+03,
      0.4909967465112e-09, 0.3210426725516e+01, 0.4059982187939e+03,
      0.4771358267658e-09, 0.5308027211629e+01, 0.1805255418145e+03,

      0.3943451445989e-09, 0.2195145341074e+01, 0.2568537517081e+03,
      0.3952109120244e-09, 0.5081189491586e+01, 0.2449975330562e+03,
      0.3788134594789e-09, 0.4345171264441e+01, 0.1568131045107e+03,
      0.3738330190479e-09, 0.2613062847997e+01, 0.3948519331910e+03,
      0.3099866678136e-09, 0.2846760817689e+01, 0.1547176098872e+03,
      0.2002962716768e-09, 0.4921360989412e+01, 0.2268582385539e+03,
      0.2198291338754e-09, 0.1130360117454e+00, 0.1658638954901e+03,
      0.1491958330784e-09, 0.4228195232278e+01, 0.2219950288015e+03,
      0.1475384076173e-09, 0.3005721811604e+00, 0.3052819430710e+03,
      0.1661626624624e-09, 0.7830125621203e+00, 0.2526661704812e+03,

      0.9015823460025e-10, 0.3807792942715e+01, 0.4171445043968e+03 };

/* Sun-to-Earth, T^0, Y */
   static const double e0y[] = {
      0.9998921098898e+00, 0.1826583913846e+00, 0.6283075850446e+01,
     -0.2442700893735e-01, 0.0000000000000e+00, 0.0000000000000e+00,
      0.8352929742915e-02, 0.1395277998680e+00, 0.1256615170089e+02,
      0.1046697300177e-03, 0.9641423109763e-01, 0.1884922755134e+02,
      0.3110841876663e-04, 0.5381140401712e+01, 0.8399684731857e+02,
      0.2570269094593e-04, 0.5301016407128e+01, 0.5296909721118e+00,
      0.2147389623610e-04, 0.2662510869850e+01, 0.1577343543434e+01,
      0.1680344384050e-04, 0.5207904119704e+01, 0.6279552690824e+01,
      0.1679117312193e-04, 0.4582187486968e+01, 0.6286599010068e+01,
      0.1440512068440e-04, 0.1900688517726e+01, 0.2352866153506e+01,

      0.1135139664999e-04, 0.5273108538556e+01, 0.5223693906222e+01,
      0.9345482571018e-05, 0.4503047687738e+01, 0.1203646072878e+02,
      0.9007418719568e-05, 0.1605621059637e+01, 0.1021328554739e+02,
      0.5671536712314e-05, 0.5812849070861e+00, 0.1059381944224e+01,
      0.7451401861666e-05, 0.2807346794836e+01, 0.3981490189893e+00,
      0.6393470057114e-05, 0.6029224133855e+01, 0.5753384878334e+01,
      0.6814275881697e-05, 0.6472990145974e+00, 0.4705732307012e+01,
      0.6113705628887e-05, 0.3813843419700e+01, 0.6812766822558e+01,
      0.4503851367273e-05, 0.4527804370996e+01, 0.5884926831456e+01,
      0.4522249141926e-05, 0.5991783029224e+01, 0.6256777527156e+01,

      0.4501794307018e-05, 0.3798703844397e+01, 0.6309374173736e+01,
      0.5514927480180e-05, 0.3961257833388e+01, 0.5507553240374e+01,
      0.4062862799995e-05, 0.5256247296369e+01, 0.6681224869435e+01,
      0.5414900429712e-05, 0.5499032014097e+01, 0.7755226100720e+00,
      0.5463153987424e-05, 0.6173092454097e+01, 0.1414349524433e+02,
      0.5071611859329e-05, 0.2870244247651e+01, 0.7860419393880e+01,
      0.2195112094455e-05, 0.2952338617201e+01, 0.1150676975667e+02,
      0.2279139233919e-05, 0.5951775132933e+01, 0.7058598460518e+01,
      0.2278386100876e-05, 0.4845456398785e+01, 0.4694002934110e+01,
      0.2559088003308e-05, 0.6945321117311e+00, 0.1216800268190e+02,

      0.2561079286856e-05, 0.6167224608301e+01, 0.7099330490126e+00,
      0.1792755796387e-05, 0.1400122509632e+01, 0.7962980379786e+00,
      0.1818715656502e-05, 0.4703347611830e+01, 0.6283142985870e+01,
      0.1818744924791e-05, 0.5086748900237e+01, 0.6283008715021e+01,
      0.1554518791390e-05, 0.5331008042713e-01, 0.2513230340178e+02,
      0.2063265737239e-05, 0.4283680484178e+01, 0.1179062909082e+02,
      0.1497613520041e-05, 0.6074207826073e+01, 0.5486777812467e+01,
      0.2000617940427e-05, 0.2501426281450e+01, 0.1778984560711e+02,
      0.1289731195580e-05, 0.3646340599536e+01, 0.7079373888424e+01,
      0.1282657998934e-05, 0.3232864804902e+01, 0.3738761453707e+01,

      0.1528915968658e-05, 0.5581433416669e+01, 0.2132990797783e+00,
      0.1187304098432e-05, 0.5453576453694e+01, 0.9437762937313e+01,
      0.7842782928118e-06, 0.2823953922273e+00, 0.8827390247185e+01,
      0.7352892280868e-06, 0.1124369580175e+01, 0.1589072916335e+01,
      0.6570189360797e-06, 0.2089154042840e+01, 0.1176985366291e+02,
      0.6324967590410e-06, 0.6704855581230e+00, 0.6262300422539e+01,
      0.6298289872283e-06, 0.2836414855840e+01, 0.6303851278352e+01,
      0.6476686465855e-06, 0.4852433866467e+00, 0.7113454667900e-02,
      0.8587034651234e-06, 0.1453511005668e+01, 0.1672837615881e+03,
      0.8068948788113e-06, 0.9224087798609e+00, 0.6069776770667e+01,

      0.8353786011661e-06, 0.4631707184895e+01, 0.3340612434717e+01,
      0.6009324532132e-06, 0.1829498827726e+01, 0.4136910472696e+01,
      0.7558158559566e-06, 0.2588596800317e+01, 0.6496374930224e+01,
      0.5809279504503e-06, 0.5516818853476e+00, 0.1097707878456e+02,
      0.5374131950254e-06, 0.6275674734960e+01, 0.1194447056968e+01,
      0.5711160507326e-06, 0.1091905956872e+01, 0.6282095334605e+01,
      0.5710183170746e-06, 0.2415001635090e+01, 0.6284056366286e+01,
      0.5144373590610e-06, 0.6020336443438e+01, 0.6290189305114e+01,
      0.5103108927267e-06, 0.3775634564605e+01, 0.6275962395778e+01,
      0.4960654697891e-06, 0.1073450946756e+01, 0.6127655567643e+01,

      0.4786385689280e-06, 0.2431178012310e+01, 0.6438496133249e+01,
      0.6109911263665e-06, 0.5343356157914e+01, 0.3154687086868e+01,
      0.4839898944024e-06, 0.5830833594047e-01, 0.8018209333619e+00,
      0.4734822623919e-06, 0.4536080134821e+01, 0.3128388763578e+01,
      0.4834741473290e-06, 0.2585090489754e+00, 0.7084896783808e+01,
      0.5134858581156e-06, 0.4213317172603e+01, 0.1235285262111e+02,
      0.5064004264978e-06, 0.4814418806478e+00, 0.1185621865188e+02,
      0.3753476772761e-06, 0.1599953399788e+01, 0.8429241228195e+01,
      0.4935264014283e-06, 0.2157417556873e+01, 0.2544314396739e+01,
      0.3950929600897e-06, 0.3359394184254e+01, 0.5481254917084e+01,

      0.4895849789777e-06, 0.5165704376558e+01, 0.9225539266174e+01,
      0.4215241688886e-06, 0.2065368800993e+01, 0.1726015463500e+02,
      0.3796773731132e-06, 0.1468606346612e+01, 0.4265981595566e+00,
      0.3114178142515e-06, 0.3615638079474e+01, 0.2146165377750e+01,
      0.3260664220838e-06, 0.4417134922435e+01, 0.4164311961999e+01,
      0.3976996123008e-06, 0.4700866883004e+01, 0.5856477690889e+01,
      0.2801459672924e-06, 0.4538902060922e+01, 0.1256967486051e+02,
      0.3638931868861e-06, 0.1334197991475e+01, 0.1807370494127e+02,
      0.2487013269476e-06, 0.3749275558275e+01, 0.2629832328990e-01,
      0.3034165481994e-06, 0.4236622030873e+00, 0.4535059491685e+01,

      0.2676278825586e-06, 0.5970848007811e+01, 0.3930209696940e+01,
      0.2764903818918e-06, 0.5194636754501e+01, 0.1256262854127e+02,
      0.2485149930507e-06, 0.1002434207846e+01, 0.5088628793478e+01,
      0.2199305540941e-06, 0.3066773098403e+01, 0.1255903824622e+02,
      0.2571106500435e-06, 0.7588312459063e+00, 0.1336797263425e+02,
      0.2049751817158e-06, 0.3444977434856e+01, 0.1137170464392e+02,
      0.2599707296297e-06, 0.1873128542205e+01, 0.7143069561767e+02,
      0.1785018072217e-06, 0.5015891306615e+01, 0.1748016358760e+01,
      0.2324833891115e-06, 0.4618271239730e+01, 0.1831953657923e+02,
      0.1709711119545e-06, 0.5300003455669e+01, 0.4933208510675e+01,

      0.2107159351716e-06, 0.2229819815115e+01, 0.7477522907414e+01,
      0.1750333080295e-06, 0.6161485880008e+01, 0.1044738781244e+02,
      0.2000598210339e-06, 0.2967357299999e+01, 0.8031092209206e+01,
      0.1380920248681e-06, 0.3027007923917e+01, 0.8635942003952e+01,
      0.1412460470299e-06, 0.6037597163798e+01, 0.2942463415728e+01,
      0.1888459803001e-06, 0.8561476243374e+00, 0.1561374759853e+03,
      0.1788370542585e-06, 0.4869736290209e+01, 0.1592596075957e+01,
      0.1360893296167e-06, 0.3626411886436e+01, 0.1309584267300e+02,
      0.1506846530160e-06, 0.1550975377427e+01, 0.1649636139783e+02,
      0.1800913376176e-06, 0.2075826033190e+01, 0.1729818233119e+02,

      0.1436261390649e-06, 0.6148876420255e+01, 0.2042657109477e+02,
      0.1220227114151e-06, 0.4382583879906e+01, 0.7632943190217e+01,
      0.1337883603592e-06, 0.2036644327361e+01, 0.1213955354133e+02,
      0.1159326650738e-06, 0.3892276994687e+01, 0.5331357529664e+01,
      0.1352853128569e-06, 0.1447950649744e+01, 0.1673046366289e+02,
      0.1433408296083e-06, 0.4457854692961e+01, 0.7342457794669e+01,
      0.1234701666518e-06, 0.1538818147151e+01, 0.6279485555400e+01,
      0.1234027192007e-06, 0.1968523220760e+01, 0.6286666145492e+01,
      0.1244024091797e-06, 0.5779803499985e+01, 0.1511046609763e+02,
      0.1097934945516e-06, 0.6210975221388e+00, 0.1098880815746e+02,

      0.1254611329856e-06, 0.2591963807998e+01, 0.1572083878776e+02,
      0.1158247286784e-06, 0.2483612812670e+01, 0.5729506548653e+01,
      0.9039078252960e-07, 0.3857554579796e+01, 0.9623688285163e+01,
      0.9108024978836e-07, 0.5826368512984e+01, 0.7234794171227e+01,
      0.8887068108436e-07, 0.3475694573987e+01, 0.6148010737701e+01,
      0.8632374035438e-07, 0.3059070488983e-01, 0.6418140963190e+01,
      0.7893186992967e-07, 0.1583194837728e+01, 0.2118763888447e+01,
      0.8297650201172e-07, 0.8519770534637e+00, 0.1471231707864e+02,
      0.1019759578988e-06, 0.1319598738732e+00, 0.1349867339771e+01,
      0.1010037696236e-06, 0.9937860115618e+00, 0.6836645152238e+01,

      0.1047727548266e-06, 0.1382138405399e+01, 0.5999216516294e+01,
      0.7351993881086e-07, 0.3833397851735e+01, 0.6040347114260e+01,
      0.9868771092341e-07, 0.2124913814390e+01, 0.6566935184597e+01,
      0.7007321959390e-07, 0.5946305343763e+01, 0.6525804586632e+01,
      0.6861411679709e-07, 0.4574654977089e+01, 0.7238675589263e+01,
      0.7554519809614e-07, 0.5949232686844e+01, 0.1253985337760e+02,
      0.9541880448335e-07, 0.3495242990564e+01, 0.2122839202813e+02,
      0.7185606722155e-07, 0.4310113471661e+01, 0.6245048154254e+01,
      0.7131360871710e-07, 0.5480309323650e+01, 0.6321103546637e+01,
      0.6651142021039e-07, 0.5411097713654e+01, 0.5327476111629e+01,

      0.8538618213667e-07, 0.1827849973951e+01, 0.1101510648075e+02,
      0.8634954288044e-07, 0.5443584943349e+01, 0.5643178611111e+01,
      0.7449415051484e-07, 0.2011535459060e+01, 0.5368044267797e+00,
      0.7421047599169e-07, 0.3464562529249e+01, 0.2354323048545e+02,
      0.6140694354424e-07, 0.5657556228815e+01, 0.1296430071988e+02,
      0.6353525143033e-07, 0.3463816593821e+01, 0.1990745094947e+01,
      0.6221964013447e-07, 0.1532259498697e+01, 0.9517183207817e+00,
      0.5852480257244e-07, 0.1375396598875e+01, 0.9555997388169e+00,
      0.6398637498911e-07, 0.2405645801972e+01, 0.2407292145756e+02,
      0.7039744069878e-07, 0.5397541799027e+01, 0.5225775174439e+00,

      0.6977997694382e-07, 0.4762347105419e+01, 0.1097355562493e+02,
      0.7460629558396e-07, 0.2711944692164e+01, 0.2200391463820e+02,
      0.5376577536101e-07, 0.2352980430239e+01, 0.1431416805965e+02,
      0.7530607893556e-07, 0.1943940180699e+01, 0.1842262939178e+02,
      0.6822928971605e-07, 0.4337651846959e+01, 0.1554202828031e+00,
      0.6220772380094e-07, 0.6716871369278e+00, 0.1845107853235e+02,
      0.6586950799043e-07, 0.2229714460505e+01, 0.5216580451554e+01,
      0.5873800565771e-07, 0.7627013920580e+00, 0.6398972393349e+00,
      0.6264346929745e-07, 0.6202785478961e+00, 0.6277552955062e+01,
      0.6257929115669e-07, 0.2886775596668e+01, 0.6288598745829e+01,

      0.5343536033409e-07, 0.1977241012051e+01, 0.4690479774488e+01,
      0.5587849781714e-07, 0.1922923484825e+01, 0.1551045220144e+01,
      0.6905100845603e-07, 0.3570757164631e+01, 0.1030928125552e+00,
      0.6178957066649e-07, 0.5197558947765e+01, 0.5230807360890e+01,
      0.6187270224331e-07, 0.8193497368922e+00, 0.5650292065779e+01,
      0.5385664291426e-07, 0.5406336665586e+01, 0.7771377146812e+02,
      0.6329363917926e-07, 0.2837760654536e+01, 0.2608790314060e+02,
      0.4546018761604e-07, 0.2933580297050e+01, 0.5535693017924e+00,
      0.6196091049375e-07, 0.4157871494377e+01, 0.8467247584405e+02,
      0.6159555108218e-07, 0.3211703561703e+01, 0.2394243902548e+03,

      0.4995340539317e-07, 0.1459098102922e+01, 0.4732030630302e+01,
      0.5457031243572e-07, 0.1430457676136e+01, 0.6179983037890e+01,
      0.4863461418397e-07, 0.2196425916730e+01, 0.9027992316901e+02,
      0.5342947626870e-07, 0.2086612890268e+01, 0.6386168663001e+01,
      0.5674296648439e-07, 0.2760204966535e+01, 0.6915859635113e+01,
      0.4745783120161e-07, 0.4245368971862e+01, 0.6282970628506e+01,
      0.4745676961198e-07, 0.5544725787016e+01, 0.6283181072386e+01,
      0.4049796869973e-07, 0.2213984363586e+01, 0.6254626709878e+01,
      0.4248333596940e-07, 0.8075781952896e+00, 0.7875671926403e+01,
      0.4027178070205e-07, 0.1293268540378e+01, 0.6311524991013e+01,

      0.4066543943476e-07, 0.3986141175804e+01, 0.3634620989887e+01,
      0.4858863787880e-07, 0.1276112738231e+01, 0.5760498333002e+01,
      0.5277398263530e-07, 0.4916111741527e+01, 0.2515860172507e+02,
      0.4105635656559e-07, 0.1725805864426e+01, 0.6709674010002e+01,
      0.4376781925772e-07, 0.2243642442106e+01, 0.6805653367890e+01,
      0.3235827894693e-07, 0.3614135118271e+01, 0.1066495398892e+01,
      0.3073244740308e-07, 0.2460873393460e+01, 0.5863591145557e+01,
      0.3088609271373e-07, 0.5678431771790e+01, 0.9917696840332e+01,
      0.3393022279836e-07, 0.3814017477291e+01, 0.1391601904066e+02,
      0.3038686508802e-07, 0.4660216229171e+01, 0.1256621883632e+02,

      0.4019677752497e-07, 0.5906906243735e+01, 0.1334167431096e+02,
      0.3288834998232e-07, 0.9536146445882e+00, 0.1620077269078e+02,
      0.3889973794631e-07, 0.3942205097644e+01, 0.7478166569050e-01,
      0.3050438987141e-07, 0.1624810271286e+01, 0.1805292951336e+02,
      0.3601142564638e-07, 0.4030467142575e+01, 0.6208294184755e+01,
      0.3689015557141e-07, 0.3648878818694e+01, 0.5966683958112e+01,
      0.3563471893565e-07, 0.5749584017096e+01, 0.6357857516136e+01,
      0.2776183170667e-07, 0.2630124187070e+01, 0.3523159621801e-02,
      0.2922350530341e-07, 0.1790346403629e+01, 0.1272157198369e+02,
      0.3511076917302e-07, 0.6142198301611e+01, 0.6599467742779e+01,

      0.3619351007632e-07, 0.1432421386492e+01, 0.6019991944201e+01,
      0.2561254711098e-07, 0.2302822475792e+01, 0.1259245002418e+02,
      0.2626903942920e-07, 0.8660470994571e+00, 0.6702560555334e+01,
      0.2550187397083e-07, 0.6069721995383e+01, 0.1057540660594e+02,
      0.2535873526138e-07, 0.1079020331795e-01, 0.3141537925223e+02,
      0.3519786153847e-07, 0.3809066902283e+01, 0.2505706758577e+03,
      0.3424651492873e-07, 0.2075435114417e+01, 0.6546159756691e+01,
      0.2372676630861e-07, 0.2057803120154e+01, 0.2388894113936e+01,
      0.2710980779541e-07, 0.1510068488010e+01, 0.1202934727411e+02,
      0.3038710889704e-07, 0.5043617528901e+01, 0.1256608456547e+02,

      0.2220364130585e-07, 0.3694793218205e+01, 0.1336244973887e+02,
      0.3025880825460e-07, 0.5450618999049e-01, 0.2908881142201e+02,
      0.2784493486864e-07, 0.3381164084502e+01, 0.1494531617769e+02,
      0.2294414142438e-07, 0.4382309025210e+01, 0.6076890225335e+01,
      0.2012723294724e-07, 0.9142212256518e+00, 0.6262720680387e+01,
      0.2036357831958e-07, 0.5676172293154e+01, 0.4701116388778e+01,
      0.2003474823288e-07, 0.2592767977625e+01, 0.6303431020504e+01,
      0.2207144900109e-07, 0.5404976271180e+01, 0.6489261475556e+01,
      0.2481664905135e-07, 0.4373284587027e+01, 0.1204357418345e+02,
      0.2674949182295e-07, 0.5859182188482e+01, 0.4590910121555e+01,

      0.2450554720322e-07, 0.4555381557451e+01, 0.1495633313810e+00,
      0.2601975986457e-07, 0.3933165584959e+01, 0.1965104848470e+02,
      0.2199860022848e-07, 0.5227977189087e+01, 0.1351787002167e+02,
      0.2448121172316e-07, 0.4858060353949e+01, 0.1162474756779e+01,
      0.1876014864049e-07, 0.5690546553605e+01, 0.6279194432410e+01,
      0.1874513219396e-07, 0.4099539297446e+01, 0.6286957268481e+01,
      0.2156380842559e-07, 0.4382594769913e+00, 0.1813929450232e+02,
      0.1981691240061e-07, 0.1829784152444e+01, 0.4686889479442e+01,
      0.2329992648539e-07, 0.2836254278973e+01, 0.1002183730415e+02,
      0.1765184135302e-07, 0.2803494925833e+01, 0.4292330755499e+01,

      0.2436368366085e-07, 0.2836897959677e+01, 0.9514313292143e+02,
      0.2164089203889e-07, 0.6127522446024e+01, 0.6037244212485e+01,
      0.1847755034221e-07, 0.3683163635008e+01, 0.2427287361862e+00,
      0.1674798769966e-07, 0.3316993867246e+00, 0.1311972100268e+02,
      0.2222542124356e-07, 0.8294097805480e+00, 0.1266924451345e+02,
      0.2071074505925e-07, 0.3659492220261e+01, 0.6528907488406e+01,
      0.1608224471835e-07, 0.4774492067182e+01, 0.1352175143971e+02,
      0.1857583439071e-07, 0.2873120597682e+01, 0.8662240327241e+01,
      0.1793018836159e-07, 0.5282441177929e+00, 0.6819880277225e+01,
      0.1575391221692e-07, 0.1320789654258e+01, 0.1102062672231e+00,

      0.1840132009557e-07, 0.1917110916256e+01, 0.6514761976723e+02,
      0.1760917288281e-07, 0.2972635937132e+01, 0.5746271423666e+01,
      0.1561779518516e-07, 0.4372569261981e+01, 0.6272439236156e+01,
      0.1558687885205e-07, 0.5416424926425e+01, 0.6293712464735e+01,
      0.1951359382579e-07, 0.3094448898752e+01, 0.2301353951334e+02,
      0.1569144275614e-07, 0.2802103689808e+01, 0.1765478049437e+02,
      0.1479130389462e-07, 0.2136435020467e+01, 0.2077542790660e-01,
      0.1467828510764e-07, 0.7072627435674e+00, 0.1052268489556e+01,
      0.1627627337440e-07, 0.3947607143237e+01, 0.6327837846670e+00,
      0.1503498479758e-07, 0.4079248909190e+01, 0.7626583626240e-01,

      0.1297967708237e-07, 0.6269637122840e+01, 0.1149965630200e+02,
      0.1374416896634e-07, 0.4175657970702e+01, 0.6016468784579e+01,
      0.1783812325219e-07, 0.1476540547560e+01, 0.3301902111895e+02,
      0.1525884228756e-07, 0.4653477715241e+01, 0.9411464614024e+01,
      0.1451067396763e-07, 0.2573001128225e+01, 0.1277945078067e+02,
      0.1297713111950e-07, 0.5612799618771e+01, 0.6549682916313e+01,
      0.1462784012820e-07, 0.4189661623870e+01, 0.1863592847156e+02,
      0.1384185980007e-07, 0.2656915472196e+01, 0.2379164476796e+01,
      0.1221497599801e-07, 0.5612515760138e+01, 0.1257326515556e+02,
      0.1560574525896e-07, 0.4783414317919e+01, 0.1887552587463e+02,

      0.1544598372036e-07, 0.2694431138063e+01, 0.1820933031200e+02,
      0.1531678928696e-07, 0.4105103489666e+01, 0.2593412433514e+02,
      0.1349321503795e-07, 0.3082437194015e+00, 0.5120601093667e+01,
      0.1252030290917e-07, 0.6124072334087e+01, 0.6993008899458e+01,
      0.1459243816687e-07, 0.3733103981697e+01, 0.3813291813120e-01,
      0.1226103625262e-07, 0.1267127706817e+01, 0.2435678079171e+02,
      0.1019449641504e-07, 0.4367790112269e+01, 0.1725663147538e+02,
      0.1380789433607e-07, 0.3387201768700e+01, 0.2458316379602e+00,
      0.1019453421658e-07, 0.9204143073737e+00, 0.6112403035119e+01,
      0.1297929434405e-07, 0.5786874896426e+01, 0.1249137003520e+02,

      0.9912677786097e-08, 0.3164232870746e+01, 0.6247047890016e+01,
      0.9829386098599e-08, 0.2586762413351e+01, 0.6453748665772e+01,
      0.1226807746104e-07, 0.6239068436607e+01, 0.5429879531333e+01,
      0.1192691755997e-07, 0.1867380051424e+01, 0.6290122169689e+01,
      0.9836499227081e-08, 0.3424716293727e+00, 0.6319103810876e+01,
      0.9642862564285e-08, 0.5661372990657e+01, 0.8273820945392e+01,
      0.1165184404862e-07, 0.5768367239093e+01, 0.1778273215245e+02,
      0.1175794418818e-07, 0.1657351222943e+01, 0.6276029531202e+01,
      0.1018948635601e-07, 0.6458292350865e+00, 0.1254537627298e+02,
      0.9500383606676e-08, 0.1054306140741e+01, 0.1256517118505e+02,

      0.1227512202906e-07, 0.2505278379114e+01, 0.2248384854122e+02,
      0.9664792009993e-08, 0.4289737277000e+01, 0.6259197520765e+01,
      0.9613285666331e-08, 0.5500597673141e+01, 0.6306954180126e+01,
      0.1117906736211e-07, 0.2361405953468e+01, 0.1779695906178e+02,
      0.9611378640782e-08, 0.2851310576269e+01, 0.2061856251104e+00,
      0.8845354852370e-08, 0.6208777705343e+01, 0.1692165728891e+01,
      0.1054046966600e-07, 0.5413091423934e+01, 0.2204125344462e+00,
      0.1215539124483e-07, 0.5613969479755e+01, 0.8257698122054e+02,
      0.9932460955209e-08, 0.1106124877015e+01, 0.1017725758696e+02,
      0.8785804715043e-08, 0.2869224476477e+01, 0.9491756770005e+00,

      0.8538084097562e-08, 0.6159640899344e+01, 0.6393282117669e+01,
      0.8648994369529e-08, 0.1374901198784e+01, 0.4804209201333e+01,
      0.1039063219067e-07, 0.5171080641327e+01, 0.1550861511662e+02,
      0.8867983926439e-08, 0.8317320304902e+00, 0.3903911373650e+01,
      0.8327495955244e-08, 0.3605591969180e+01, 0.6172869583223e+01,
      0.9243088356133e-08, 0.6114299196843e+01, 0.6267823317922e+01,
      0.9205657357835e-08, 0.3675153683737e+01, 0.6298328382969e+01,
      0.1033269714606e-07, 0.3313328813024e+01, 0.5573142801433e+01,
      0.8001706275552e-08, 0.2019980960053e+01, 0.2648454860559e+01,
      0.9171858254191e-08, 0.8992015524177e+00, 0.1498544001348e+03,

      0.1075327150242e-07, 0.2898669963648e+01, 0.3694923081589e+02,
      0.9884866689828e-08, 0.4946715904478e+01, 0.1140367694411e+02,
      0.9541835576677e-08, 0.2371787888469e+01, 0.1256713221673e+02,
      0.7739903376237e-08, 0.2213775190612e+01, 0.7834121070590e+01,
      0.7311962684106e-08, 0.3429378787739e+01, 0.1192625446156e+02,
      0.9724904869624e-08, 0.6195878564404e+01, 0.2280573557157e+02,
      0.9251628983612e-08, 0.6511509527390e+00, 0.2787043132925e+01,
      0.7320763787842e-08, 0.6001083639421e+01, 0.6282655592598e+01,
      0.7320296650962e-08, 0.3789073265087e+01, 0.6283496108294e+01,
      0.7947032271039e-08, 0.1059659582204e+01, 0.1241073141809e+02,

      0.9005277053115e-08, 0.1280315624361e+01, 0.6281591679874e+01,
      0.8995601652048e-08, 0.2224439106766e+01, 0.6284560021018e+01,
      0.8288040568796e-08, 0.5234914433867e+01, 0.1241658836951e+02,
      0.6359381347255e-08, 0.4137989441490e+01, 0.1596186371003e+01,
      0.8699572228626e-08, 0.1758411009497e+01, 0.6133512519065e+01,
      0.6456797542736e-08, 0.5919285089994e+01, 0.1685848245639e+02,
      0.7424573475452e-08, 0.5414616938827e+01, 0.4061219149443e+01,
      0.7235671196168e-08, 0.1496516557134e+01, 0.1610006857377e+03,
      0.8104015182733e-08, 0.1919918242764e+01, 0.8460828644453e+00,
      0.8098576535937e-08, 0.3819615855458e+01, 0.3894181736510e+01,

      0.6275292346625e-08, 0.6244264115141e+01, 0.8531963191132e+00,
      0.6052432989112e-08, 0.5037731872610e+00, 0.1567108171867e+02,
      0.5705651535817e-08, 0.2984557271995e+01, 0.1258692712880e+02,
      0.5789650115138e-08, 0.6087038140697e+01, 0.1193336791622e+02,
      0.5512132153377e-08, 0.5855668994076e+01, 0.1232342296471e+02,
      0.7388890819102e-08, 0.2443128574740e+01, 0.4907302013889e+01,
      0.5467593991798e-08, 0.3017561234194e+01, 0.1884211409667e+02,
      0.6388519802999e-08, 0.5887386712935e+01, 0.5217580628120e+02,
      0.6106777149944e-08, 0.3483461059895e+00, 0.1422690933580e-01,
      0.7383420275489e-08, 0.5417387056707e+01, 0.2358125818164e+02,

      0.5505208141738e-08, 0.2848193644783e+01, 0.1151388321134e+02,
      0.6310757462877e-08, 0.2349882520828e+01, 0.1041998632314e+02,
      0.6166904929691e-08, 0.5728575944077e+00, 0.6151533897323e+01,
      0.5263442042754e-08, 0.4495796125937e+01, 0.1885275071096e+02,
      0.5591828082629e-08, 0.1355441967677e+01, 0.4337116142245e+00,
      0.5397051680497e-08, 0.1673422864307e+01, 0.6286362197481e+01,
      0.5396992745159e-08, 0.1833502206373e+01, 0.6279789503410e+01,
      0.6572913000726e-08, 0.3331122065824e+01, 0.1176433076753e+02,
      0.5123421866413e-08, 0.2165327142679e+01, 0.1245594543367e+02,
      0.5930495725999e-08, 0.2931146089284e+01, 0.6414617803568e+01,

      0.6431797403933e-08, 0.4134407994088e+01, 0.1350651127443e+00,
      0.5003182207604e-08, 0.3805420303749e+01, 0.1096996532989e+02,
      0.5587731032504e-08, 0.1082469260599e+01, 0.6062663316000e+01,
      0.5935263407816e-08, 0.8384333678401e+00, 0.5326786718777e+01,
      0.4756019827760e-08, 0.3552588749309e+01, 0.3104930017775e+01,
      0.6599951172637e-08, 0.4320826409528e+01, 0.4087944051283e+02,
      0.5902606868464e-08, 0.4811879454445e+01, 0.5849364236221e+01,
      0.5921147809031e-08, 0.9942628922396e-01, 0.1581959461667e+01,
      0.5505382581266e-08, 0.2466557607764e+01, 0.6503488384892e+01,
      0.5353771071862e-08, 0.4551978748683e+01, 0.1735668374386e+03,

      0.5063282210946e-08, 0.5710812312425e+01, 0.1248988586463e+02,
      0.5926120403383e-08, 0.1333998428358e+01, 0.2673594526851e+02,
      0.5211016176149e-08, 0.4649315360760e+01, 0.2460261242967e+02,
      0.5347075084894e-08, 0.5512754081205e+01, 0.4171425416666e+01,
      0.4872609773574e-08, 0.1308025299938e+01, 0.5333900173445e+01,
      0.4727711321420e-08, 0.2144908368062e+01, 0.7232251527446e+01,
      0.6029426018652e-08, 0.5567259412084e+01, 0.3227113045244e+03,
      0.4321485284369e-08, 0.5230667156451e+01, 0.9388005868221e+01,
      0.4476406760553e-08, 0.6134081115303e+01, 0.5547199253223e+01,
      0.5835268277420e-08, 0.4783808492071e+01, 0.7285056171570e+02,

      0.5172183602748e-08, 0.5161817911099e+01, 0.1884570439172e+02,
      0.5693571465184e-08, 0.1381646203111e+01, 0.9723862754494e+02,
      0.4060634965349e-08, 0.3876705259495e+00, 0.4274518229222e+01,
      0.3967398770473e-08, 0.5029491776223e+01, 0.3496032717521e+01,
      0.3943754005255e-08, 0.1923162955490e+01, 0.6244942932314e+01,
      0.4781323427824e-08, 0.4633332586423e+01, 0.2929661536378e+02,
      0.3871483781204e-08, 0.1616650009743e+01, 0.6321208768577e+01,
      0.5141741733997e-08, 0.9817316704659e-01, 0.1232032006293e+02,
      0.4002385978497e-08, 0.3656161212139e+01, 0.7018952447668e+01,
      0.4901092604097e-08, 0.4404098713092e+01, 0.1478866649112e+01,

      0.3740932630345e-08, 0.5181188732639e+00, 0.6922973089781e+01,
      0.4387283718538e-08, 0.3254859566869e+01, 0.2331413144044e+03,
      0.5019197802033e-08, 0.3086773224677e+01, 0.1715706182245e+02,
      0.3834931695175e-08, 0.2797882673542e+01, 0.1491901785440e+02,
      0.3760413942497e-08, 0.2892676280217e+01, 0.1726726808967e+02,
      0.3719717204628e-08, 0.5861046025739e+01, 0.6297302759782e+01,
      0.4145623530149e-08, 0.2168239627033e+01, 0.1376059875786e+02,
      0.3932788425380e-08, 0.6271811124181e+01, 0.7872148766781e+01,
      0.3686377476857e-08, 0.3936853151404e+01, 0.6268848941110e+01,
      0.3779077950339e-08, 0.1404148734043e+01, 0.4157198507331e+01,

      0.4091334550598e-08, 0.2452436180854e+01, 0.9779108567966e+01,
      0.3926694536146e-08, 0.6102292739040e+01, 0.1098419223922e+02,
      0.4841000253289e-08, 0.6072760457276e+01, 0.1252801878276e+02,
      0.4949340130240e-08, 0.1154832815171e+01, 0.1617106187867e+03,
      0.3761557737360e-08, 0.5527545321897e+01, 0.3185192151914e+01,
      0.3647396268188e-08, 0.1525035688629e+01, 0.6271346477544e+01,
      0.3932405074189e-08, 0.5570681040569e+01, 0.2139354194808e+02,
      0.3631322501141e-08, 0.1981240601160e+01, 0.6294805223347e+01,
      0.4130007425139e-08, 0.2050060880201e+01, 0.2195415756911e+02,
      0.4433905965176e-08, 0.3277477970321e+01, 0.7445550607224e+01,

      0.3851814176947e-08, 0.5210690074886e+01, 0.9562891316684e+00,
      0.3485807052785e-08, 0.6653274904611e+00, 0.1161697602389e+02,
      0.3979772816991e-08, 0.1767941436148e+01, 0.2277943724828e+02,
      0.3402607460500e-08, 0.3421746306465e+01, 0.1087398597200e+02,
      0.4049993000926e-08, 0.1127144787547e+01, 0.3163918923335e+00,
      0.3420511182382e-08, 0.4214794779161e+01, 0.1362553364512e+02,
      0.3640772365012e-08, 0.5324905497687e+01, 0.1725304118033e+02,
      0.3323037987501e-08, 0.6135761838271e+01, 0.6279143387820e+01,
      0.4503141663637e-08, 0.1802305450666e+01, 0.1385561574497e+01,
      0.4314560055588e-08, 0.4812299731574e+01, 0.4176041334900e+01,

      0.3294226949110e-08, 0.3657547059723e+01, 0.6287008313071e+01,
      0.3215657197281e-08, 0.4866676894425e+01, 0.5749861718712e+01,
      0.4129362656266e-08, 0.3809342558906e+01, 0.5905702259363e+01,
      0.3137762976388e-08, 0.2494635174443e+01, 0.2099539292909e+02,
      0.3514010952384e-08, 0.2699961831678e+01, 0.7335344340001e+01,
      0.3327607571530e-08, 0.3318457714816e+01, 0.5436992986000e+01,
      0.3541066946675e-08, 0.4382703582466e+01, 0.1234573916645e+02,
      0.3216179847052e-08, 0.5271066317054e+01, 0.3802769619140e-01,
      0.2959045059570e-08, 0.5819591585302e+01, 0.2670964694522e+02,
      0.3884040326665e-08, 0.5980934960428e+01, 0.6660449441528e+01,

      0.2922027539886e-08, 0.3337290282483e+01, 0.1375773836557e+01,
      0.4110846382042e-08, 0.5742978187327e+01, 0.4480965020977e+02,
      0.2934508411032e-08, 0.2278075804200e+01, 0.6408777551755e+00,
      0.3966896193000e-08, 0.5835747858477e+01, 0.3773735910827e+00,
      0.3286695827610e-08, 0.5838898193902e+01, 0.3932462625300e-02,
      0.3720643094196e-08, 0.1122212337858e+01, 0.1646033343740e+02,
      0.3285508906174e-08, 0.9182250996416e+00, 0.1081813534213e+02,
      0.3753880575973e-08, 0.5174761973266e+01, 0.5642198095270e+01,
      0.3022129385587e-08, 0.3381611020639e+01, 0.2982630633589e+02,
      0.2798569205621e-08, 0.3546193723922e+01, 0.1937891852345e+02,

      0.3397872070505e-08, 0.4533203197934e+01, 0.6923953605621e+01,
      0.3708099772977e-08, 0.2756168198616e+01, 0.3066615496545e+02,
      0.3599283541510e-08, 0.1934395469918e+01, 0.6147450479709e+01,
      0.3688702753059e-08, 0.7149920971109e+00, 0.2636725487657e+01,
      0.2681084724003e-08, 0.4899819493154e+01, 0.6816289982179e+01,
      0.3495993460759e-08, 0.1572418915115e+01, 0.6418701221183e+01,
      0.3130770324995e-08, 0.8912190180489e+00, 0.1235996607578e+02,
      0.2744353821941e-08, 0.3800821940055e+01, 0.2059724391010e+02,
      0.2842732906341e-08, 0.2644717440029e+01, 0.2828699048865e+02,
      0.3046882682154e-08, 0.3987793020179e+01, 0.6055599646783e+01,

      0.2399072455143e-08, 0.9908826440764e+00, 0.6255674361143e+01,
      0.2384306274204e-08, 0.2516149752220e+01, 0.6310477339748e+01,
      0.2977324500559e-08, 0.5849195642118e+01, 0.1652265972112e+02,
      0.3062835258972e-08, 0.1681660100162e+01, 0.1172006883645e+02,
      0.3109682589231e-08, 0.5804143987737e+00, 0.2751146787858e+02,
      0.2903920355299e-08, 0.5800768280123e+01, 0.6510552054109e+01,
      0.2823221989212e-08, 0.9241118370216e+00, 0.5469525544182e+01,
      0.3187949696649e-08, 0.3139776445735e+01, 0.1693792562116e+03,
      0.2922559771655e-08, 0.3549440782984e+01, 0.2630839062450e+00,
      0.2436302066603e-08, 0.4735540696319e+01, 0.3946258593675e+00,

      0.3049473043606e-08, 0.4998289124561e+01, 0.8390110365991e+01,
      0.2863682575784e-08, 0.6709515671102e+00, 0.2243449970715e+00,
      0.2641750517966e-08, 0.5410978257284e+01, 0.2986433403208e+02,
      0.2704093466243e-08, 0.4778317207821e+01, 0.6129297044991e+01,
      0.2445522177011e-08, 0.6009020662222e+01, 0.1171295538178e+02,
      0.2623608810230e-08, 0.5010449777147e+01, 0.6436854655901e+01,
      0.2079259704053e-08, 0.5980943768809e+01, 0.2019909489111e+02,
      0.2820225596771e-08, 0.2679965110468e+01, 0.5934151399930e+01,
      0.2365221950927e-08, 0.1894231148810e+01, 0.2470570524223e+02,
      0.2359682077149e-08, 0.4220752950780e+01, 0.8671969964381e+01,

      0.2387577137206e-08, 0.2571783940617e+01, 0.7096626156709e+01,
      0.1982102089816e-08, 0.5169765997119e+00, 0.1727188400790e+02,
      0.2687502389925e-08, 0.6239078264579e+01, 0.7075506709219e+02,
      0.2207751669135e-08, 0.2031184412677e+01, 0.4377611041777e+01,
      0.2618370214274e-08, 0.8266079985979e+00, 0.6632000300961e+01,
      0.2591951887361e-08, 0.8819350522008e+00, 0.4873985990671e+02,
      0.2375055656248e-08, 0.3520944177789e+01, 0.1590676413561e+02,
      0.2472019978911e-08, 0.1551431908671e+01, 0.6612329252343e+00,
      0.2368157127199e-08, 0.4178610147412e+01, 0.3459636466239e+02,
      0.1764846605693e-08, 0.1506764000157e+01, 0.1980094587212e+02,

      0.2291769608798e-08, 0.2118250611782e+01, 0.2844914056730e-01,
      0.2209997316943e-08, 0.3363255261678e+01, 0.2666070658668e+00,
      0.2292699097923e-08, 0.4200423956460e+00, 0.1484170571900e-02,
      0.1629683015329e-08, 0.2331362582487e+01, 0.3035599730800e+02,
      0.2206492862426e-08, 0.3400274026992e+01, 0.6281667977667e+01,
      0.2205746568257e-08, 0.1066051230724e+00, 0.6284483723224e+01,
      0.2026310767991e-08, 0.2779066487979e+01, 0.2449240616245e+02,
      0.1762977622163e-08, 0.9951450691840e+00, 0.2045286941806e+02,
      0.1368535049606e-08, 0.6402447365817e+00, 0.2473415438279e+02,
      0.1720598775450e-08, 0.2303524214705e+00, 0.1679593901136e+03,

      0.1702429015449e-08, 0.6164622655048e+01, 0.3338575901272e+03,
      0.1414033197685e-08, 0.3954561185580e+01, 0.1624205518357e+03,
      0.1573768958043e-08, 0.2028286308984e+01, 0.3144167757552e+02,
      0.1650705184447e-08, 0.2304040666128e+01, 0.5267006960365e+02,
      0.1651087618855e-08, 0.2538461057280e+01, 0.8956999012000e+02,
      0.1616409518983e-08, 0.5111054348152e+01, 0.3332657872986e+02,
      0.1537175173581e-08, 0.5601130666603e+01, 0.3852657435933e+02,
      0.1593191980553e-08, 0.2614340453411e+01, 0.2282781046519e+03,
      0.1499480170643e-08, 0.3624721577264e+01, 0.2823723341956e+02,
      0.1493807843235e-08, 0.4214569879008e+01, 0.2876692439167e+02,

      0.1074571199328e-08, 0.1496911744704e+00, 0.8397383534231e+02,
      0.1074406983417e-08, 0.1187817671922e+01, 0.8401985929482e+02,
      0.9757576855851e-09, 0.2655703035858e+01, 0.7826370942180e+02,
      0.1258432887565e-08, 0.4969896184844e+01, 0.3115650189215e+03,
      0.1240336343282e-08, 0.5192460776926e+01, 0.1784300471910e+03,
      0.9016107005164e-09, 0.1960356923057e+01, 0.5886454391678e+02,
      0.1135392360918e-08, 0.5082427809068e+01, 0.7842370451713e+02,
      0.9216046089565e-09, 0.2793775037273e+01, 0.1014262087719e+03,
      0.1061276615030e-08, 0.3726144311409e+01, 0.5660027930059e+02,
      0.1010110596263e-08, 0.7404080708937e+00, 0.4245678405627e+02,

      0.7217424756199e-09, 0.2697449980577e-01, 0.2457074661053e+03,
      0.6912003846756e-09, 0.4253296276335e+01, 0.1679936946371e+03,
      0.6871814664847e-09, 0.5148072412354e+01, 0.6053048899753e+02,
      0.4887158016343e-09, 0.2153581148294e+01, 0.9656299901946e+02,
      0.5161802866314e-09, 0.3852750634351e+01, 0.2442876000072e+03,
      0.5652599559057e-09, 0.1233233356270e+01, 0.8365903305582e+02,
      0.4710812608586e-09, 0.5610486976767e+01, 0.3164282286739e+03,
      0.4909977500324e-09, 0.1639629524123e+01, 0.4059982187939e+03,
      0.4772641839378e-09, 0.3737100368583e+01, 0.1805255418145e+03,
      0.4487562567153e-09, 0.1158417054478e+00, 0.8433466158131e+02,

      0.3943441230497e-09, 0.6243502862796e+00, 0.2568537517081e+03,
      0.3952236913598e-09, 0.3510377382385e+01, 0.2449975330562e+03,
      0.3788898363417e-09, 0.5916128302299e+01, 0.1568131045107e+03,
      0.3738329328831e-09, 0.1042266763456e+01, 0.3948519331910e+03,
      0.2451199165151e-09, 0.1166788435700e+01, 0.1435713242844e+03,
      0.2436734402904e-09, 0.3254726114901e+01, 0.2268582385539e+03,
      0.2213605274325e-09, 0.1687210598530e+01, 0.1658638954901e+03,
      0.1491521204829e-09, 0.2657541786794e+01, 0.2219950288015e+03,
      0.1474995329744e-09, 0.5013089805819e+01, 0.3052819430710e+03,
      0.1661939475656e-09, 0.5495315428418e+01, 0.2526661704812e+03,

      0.9015946748003e-10, 0.2236989966505e+01, 0.4171445043968e+03 };

/* Sun-to-Earth, T^0, Z */
   static const double e0z[] = {
      0.2796207639075e-05, 0.3198701560209e+01, 0.8433466158131e+02,
      0.1016042198142e-05, 0.5422360395913e+01, 0.5507553240374e+01,
      0.8044305033647e-06, 0.3880222866652e+01, 0.5223693906222e+01,
      0.4385347909274e-06, 0.3704369937468e+01, 0.2352866153506e+01,
      0.3186156414906e-06, 0.3999639363235e+01, 0.1577343543434e+01,
      0.2272412285792e-06, 0.3984738315952e+01, 0.1047747311755e+01,
      0.1645620103007e-06, 0.3565412516841e+01, 0.5856477690889e+01,
      0.1815836921166e-06, 0.4984507059020e+01, 0.6283075850446e+01,
      0.1447461676364e-06, 0.3702753570108e+01, 0.9437762937313e+01,
      0.1430760876382e-06, 0.3409658712357e+01, 0.1021328554739e+02,

      0.1120445753226e-06, 0.4829561570246e+01, 0.1414349524433e+02,
      0.1090232840797e-06, 0.2080729178066e+01, 0.6812766822558e+01,
      0.9715727346551e-07, 0.3476295881948e+01, 0.4694002934110e+01,
      0.1036267136217e-06, 0.4056639536648e+01, 0.7109288135493e+02,
      0.8752665271340e-07, 0.4448159519911e+01, 0.5753384878334e+01,
      0.8331864956004e-07, 0.4991704044208e+01, 0.7084896783808e+01,
      0.6901658670245e-07, 0.4325358994219e+01, 0.6275962395778e+01,
      0.9144536848998e-07, 0.1141826375363e+01, 0.6620890113188e+01,
      0.7205085037435e-07, 0.3624344170143e+01, 0.5296909721118e+00,
      0.7697874654176e-07, 0.5554257458998e+01, 0.1676215758509e+03,

      0.5197545738384e-07, 0.6251760961735e+01, 0.1807370494127e+02,
      0.5031345378608e-07, 0.2497341091913e+01, 0.4705732307012e+01,
      0.4527110205840e-07, 0.2335079920992e+01, 0.6309374173736e+01,
      0.4753355798089e-07, 0.7094148987474e+00, 0.5884926831456e+01,
      0.4296951977516e-07, 0.1101916352091e+01, 0.6681224869435e+01,
      0.3855341568387e-07, 0.1825495405486e+01, 0.5486777812467e+01,
      0.5253930970990e-07, 0.4424740687208e+01, 0.7860419393880e+01,
      0.4024630496471e-07, 0.5120498157053e+01, 0.1336797263425e+02,
      0.4061069791453e-07, 0.6029771435451e+01, 0.3930209696940e+01,
      0.3797883804205e-07, 0.4435193600836e+00, 0.3154687086868e+01,

      0.2933033225587e-07, 0.5124157356507e+01, 0.1059381944224e+01,
      0.3503000930426e-07, 0.5421830162065e+01, 0.6069776770667e+01,
      0.3670096214050e-07, 0.4582101667297e+01, 0.1219403291462e+02,
      0.2905609437008e-07, 0.1926566420072e+01, 0.1097707878456e+02,
      0.2466827821713e-07, 0.6090174539834e+00, 0.6496374930224e+01,
      0.2691647295332e-07, 0.1393432595077e+01, 0.2200391463820e+02,
      0.2150554667946e-07, 0.4308671715951e+01, 0.5643178611111e+01,
      0.2237481922680e-07, 0.8133968269414e+00, 0.8635942003952e+01,
      0.1817741038157e-07, 0.3755205127454e+01, 0.3340612434717e+01,
      0.2227820762132e-07, 0.2759558596664e+01, 0.1203646072878e+02,

      0.1944713772307e-07, 0.5699645869121e+01, 0.1179062909082e+02,
      0.1527340520662e-07, 0.1986749091746e+01, 0.3981490189893e+00,
      0.1577282574914e-07, 0.3205017217983e+01, 0.5088628793478e+01,
      0.1424738825424e-07, 0.6256747903666e+01, 0.2544314396739e+01,
      0.1616563121701e-07, 0.2601671259394e+00, 0.1729818233119e+02,
      0.1401210391692e-07, 0.4686939173506e+01, 0.7058598460518e+01,
      0.1488726974214e-07, 0.2815862451372e+01, 0.2593412433514e+02,
      0.1692626442388e-07, 0.4956894109797e+01, 0.1564752902480e+03,
      0.1123571582910e-07, 0.2381192697696e+01, 0.3738761453707e+01,
      0.9903308606317e-08, 0.4294851657684e+01, 0.9225539266174e+01,

      0.9174533187191e-08, 0.3075171510642e+01, 0.4164311961999e+01,
      0.8645985631457e-08, 0.5477534821633e+00, 0.8429241228195e+01,
     -0.1085876492688e-07, 0.0000000000000e+00, 0.0000000000000e+00,
      0.9264309077815e-08, 0.5968571670097e+01, 0.7079373888424e+01,
      0.8243116984954e-08, 0.1489098777643e+01, 0.1044738781244e+02,
      0.8268102113708e-08, 0.3512977691983e+01, 0.1150676975667e+02,
      0.9043613988227e-08, 0.1290704408221e+00, 0.1101510648075e+02,
      0.7432912038789e-08, 0.1991086893337e+01, 0.2608790314060e+02,
      0.8586233727285e-08, 0.4238357924414e+01, 0.2986433403208e+02,
      0.7612230060131e-08, 0.2911090150166e+01, 0.4732030630302e+01,

      0.7097787751408e-08, 0.1908938392390e+01, 0.8031092209206e+01,
      0.7640237040175e-08, 0.6129219000168e+00, 0.7962980379786e+00,
      0.7070445688081e-08, 0.1380417036651e+01, 0.2146165377750e+01,
      0.7690770957702e-08, 0.1680504249084e+01, 0.2122839202813e+02,
      0.8051292542594e-08, 0.5127423484511e+01, 0.2942463415728e+01,
      0.5902709104515e-08, 0.2020274190917e+01, 0.7755226100720e+00,
      0.5134567496462e-08, 0.2606778676418e+01, 0.1256615170089e+02,
      0.5525802046102e-08, 0.1613011769663e+01, 0.8018209333619e+00,
      0.5880724784221e-08, 0.4604483417236e+01, 0.4690479774488e+01,
      0.5211699081370e-08, 0.5718964114193e+01, 0.8827390247185e+01,

      0.4891849573562e-08, 0.3689658932196e+01, 0.2132990797783e+00,
      0.5150246069997e-08, 0.4099769855122e+01, 0.6480980550449e+02,
      0.5102434319633e-08, 0.5660834602509e+01, 0.3379454372902e+02,
      0.5083405254252e-08, 0.9842221218974e+00, 0.4136910472696e+01,
      0.4206562585682e-08, 0.1341363634163e+00, 0.3128388763578e+01,
      0.4663249683579e-08, 0.8130132735866e+00, 0.5216580451554e+01,
      0.4099474416530e-08, 0.5791497770644e+01, 0.4265981595566e+00,
      0.4628251220767e-08, 0.1249802769331e+01, 0.1572083878776e+02,
      0.5024068728142e-08, 0.4795684802743e+01, 0.6290189305114e+01,
      0.5120234327758e-08, 0.3810420387208e+01, 0.5230807360890e+01,

      0.5524029815280e-08, 0.1029264714351e+01, 0.2397622045175e+03,
      0.4757415718860e-08, 0.3528044781779e+01, 0.1649636139783e+02,
      0.3915786131127e-08, 0.5593889282646e+01, 0.1589072916335e+01,
      0.4869053149991e-08, 0.3299636454433e+01, 0.7632943190217e+01,
      0.3649365703729e-08, 0.1286049002584e+01, 0.6206810014183e+01,
      0.3992493949002e-08, 0.3100307589464e+01, 0.2515860172507e+02,
      0.3320247477418e-08, 0.6212683940807e+01, 0.1216800268190e+02,
      0.3287123739696e-08, 0.4699118445928e+01, 0.7234794171227e+01,
      0.3472776811103e-08, 0.2630507142004e+01, 0.7342457794669e+01,
      0.3423253294767e-08, 0.2946432844305e+01, 0.9623688285163e+01,

      0.3896173898244e-08, 0.1224834179264e+01, 0.6438496133249e+01,
      0.3388455337924e-08, 0.1543807616351e+01, 0.1494531617769e+02,
      0.3062704716523e-08, 0.1191777572310e+01, 0.8662240327241e+01,
      0.3270075600400e-08, 0.5483498767737e+01, 0.1194447056968e+01,
      0.3101209215259e-08, 0.8000833804348e+00, 0.3772475342596e+02,
      0.2780883347311e-08, 0.4077980721888e+00, 0.5863591145557e+01,
      0.2903605931824e-08, 0.2617490302147e+01, 0.1965104848470e+02,
      0.2682014743119e-08, 0.2634703158290e+01, 0.7238675589263e+01,
      0.2534360108492e-08, 0.6102446114873e+01, 0.6836645152238e+01,
      0.2392564882509e-08, 0.3681820208691e+01, 0.5849364236221e+01,

      0.2656667254856e-08, 0.6216045388886e+01, 0.6133512519065e+01,
      0.2331242096773e-08, 0.5864949777744e+01, 0.4535059491685e+01,
      0.2287898363668e-08, 0.4566628532802e+01, 0.7477522907414e+01,
      0.2336944521306e-08, 0.2442722126930e+01, 0.1137170464392e+02,
      0.3156632236269e-08, 0.1626628050682e+01, 0.2509084901204e+03,
      0.2982612402766e-08, 0.2803604512609e+01, 0.1748016358760e+01,
      0.2774031674807e-08, 0.4654002897158e+01, 0.8223916695780e+02,
      0.2295236548638e-08, 0.4326518333253e+01, 0.3378142627421e+00,
      0.2190714699873e-08, 0.4519614578328e+01, 0.2908881142201e+02,
      0.2191495845045e-08, 0.3012626912549e+01, 0.1673046366289e+02,

      0.2492901628386e-08, 0.1290101424052e+00, 0.1543797956245e+03,
      0.1993778064319e-08, 0.3864046799414e+01, 0.1778984560711e+02,
      0.1898146479022e-08, 0.5053777235891e+01, 0.2042657109477e+02,
      0.1918280127634e-08, 0.2222470192548e+01, 0.4165496312290e+02,
      0.1916351061607e-08, 0.8719067257774e+00, 0.7737595720538e+02,
      0.1834720181466e-08, 0.4031491098040e+01, 0.2358125818164e+02,
      0.1249201523806e-08, 0.5938379466835e+01, 0.3301902111895e+02,
      0.1477304050539e-08, 0.6544722606797e+00, 0.9548094718417e+02,
      0.1264316431249e-08, 0.2059072853236e+01, 0.8399684731857e+02,
      0.1203526495039e-08, 0.3644813532605e+01, 0.4558517281984e+02,

      0.9221681059831e-09, 0.3241815055602e+01, 0.7805158573086e+02,
      0.7849278367646e-09, 0.5043812342457e+01, 0.5217580628120e+02,
      0.7983392077387e-09, 0.5000024502753e+01, 0.1501922143975e+03,
      0.7925395431654e-09, 0.1398734871821e-01, 0.9061773743175e+02,
      0.7640473285886e-09, 0.5067111723130e+01, 0.4951538251678e+02,
      0.5398937754482e-09, 0.5597382200075e+01, 0.1613385000004e+03,
      0.5626247550193e-09, 0.2601338209422e+01, 0.7318837597844e+02,
      0.5525197197855e-09, 0.5814832109256e+01, 0.1432335100216e+03,
      0.5407629837898e-09, 0.3384820609076e+01, 0.3230491187871e+03,
      0.3856739119801e-09, 0.1072391840473e+01, 0.2334791286671e+03,

      0.3856425239987e-09, 0.2369540393327e+01, 0.1739046517013e+03,
      0.4350867755983e-09, 0.5255575751082e+01, 0.1620484330494e+03,
      0.3844113924996e-09, 0.5482356246182e+01, 0.9757644180768e+02,
      0.2854869155431e-09, 0.9573634763143e+00, 0.1697170704744e+03,
      0.1719227671416e-09, 0.1887203025202e+01, 0.2265204242912e+03,
      0.1527846879755e-09, 0.3982183931157e+01, 0.3341954043900e+03,
      0.1128229264847e-09, 0.2787457156298e+01, 0.3119028331842e+03 };

/* Sun-to-Earth, T^1, X */
   static const double e1x[] = {
      0.1234046326004e-05, 0.0000000000000e+00, 0.0000000000000e+00,
      0.5150068824701e-06, 0.6002664557501e+01, 0.1256615170089e+02,
      0.1290743923245e-07, 0.5959437664199e+01, 0.1884922755134e+02,
      0.1068615564952e-07, 0.2015529654209e+01, 0.6283075850446e+01,
      0.2079619142538e-08, 0.1732960531432e+01, 0.6279552690824e+01,
      0.2078009243969e-08, 0.4915604476996e+01, 0.6286599010068e+01,
      0.6206330058856e-09, 0.3616457953824e+00, 0.4705732307012e+01,
      0.5989335313746e-09, 0.3802607304474e+01, 0.6256777527156e+01,
      0.5958495663840e-09, 0.2845866560031e+01, 0.6309374173736e+01,
      0.4866923261539e-09, 0.5213203771824e+01, 0.7755226100720e+00,

      0.4267785823142e-09, 0.4368189727818e+00, 0.1059381944224e+01,
      0.4610675141648e-09, 0.1837249181372e-01, 0.7860419393880e+01,
      0.3626989993973e-09, 0.2161590545326e+01, 0.5753384878334e+01,
      0.3563071194389e-09, 0.1452631954746e+01, 0.5884926831456e+01,
      0.3557015642807e-09, 0.4470593393054e+01, 0.6812766822558e+01,
      0.3210412089122e-09, 0.5195926078314e+01, 0.6681224869435e+01,
      0.2875473577986e-09, 0.5916256610193e+01, 0.2513230340178e+02,
      0.2842913681629e-09, 0.1149902426047e+01, 0.6127655567643e+01,
      0.2751248215916e-09, 0.5502088574662e+01, 0.6438496133249e+01,
      0.2481432881127e-09, 0.2921989846637e+01, 0.5486777812467e+01,

      0.2059885976560e-09, 0.3718070376585e+01, 0.7079373888424e+01,
      0.2015522342591e-09, 0.5979395259740e+01, 0.6290189305114e+01,
      0.1995364084253e-09, 0.6772087985494e+00, 0.6275962395778e+01,
      0.1957436436943e-09, 0.2899210654665e+01, 0.5507553240374e+01,
      0.1651609818948e-09, 0.6228206482192e+01, 0.1150676975667e+02,
      0.1822980550699e-09, 0.1469348746179e+01, 0.1179062909082e+02,
      0.1675223159760e-09, 0.3813910555688e+01, 0.7058598460518e+01,
      0.1706491764745e-09, 0.3004380506684e+00, 0.7113454667900e-02,
      0.1392952362615e-09, 0.1440393973406e+01, 0.7962980379786e+00,
      0.1209868266342e-09, 0.4150425791727e+01, 0.4694002934110e+01,

      0.1009827202611e-09, 0.3290040429843e+01, 0.3738761453707e+01,
      0.1047261388602e-09, 0.4229590090227e+01, 0.6282095334605e+01,
      0.1047006652004e-09, 0.2418967680575e+01, 0.6284056366286e+01,
      0.9609993143095e-10, 0.4627943659201e+01, 0.6069776770667e+01,
      0.9590900593873e-10, 0.1894393939924e+01, 0.4136910472696e+01,
      0.9146249188071e-10, 0.2010647519562e+01, 0.6496374930224e+01,
      0.8545274480290e-10, 0.5529846956226e-01, 0.1194447056968e+01,
      0.8224377881194e-10, 0.1254304102174e+01, 0.1589072916335e+01,
      0.6183529510410e-10, 0.3360862168815e+01, 0.8827390247185e+01,
      0.6259255147141e-10, 0.4755628243179e+01, 0.8429241228195e+01,

      0.5539291694151e-10, 0.5371746955142e+01, 0.4933208510675e+01,
      0.7328259466314e-10, 0.4927699613906e+00, 0.4535059491685e+01,
      0.6017835843560e-10, 0.5776682001734e-01, 0.1255903824622e+02,
      0.7079827775243e-10, 0.4395059432251e+01, 0.5088628793478e+01,
      0.5170358878213e-10, 0.5154062619954e+01, 0.1176985366291e+02,
      0.4872301838682e-10, 0.6289611648973e+00, 0.6040347114260e+01,
      0.5249869411058e-10, 0.5617272046949e+01, 0.3154687086868e+01,
      0.4716172354411e-10, 0.3965901800877e+01, 0.5331357529664e+01,
      0.4871214940964e-10, 0.4627507050093e+01, 0.1256967486051e+02,
      0.4598076850751e-10, 0.6023631226459e+01, 0.6525804586632e+01,

      0.4562196089485e-10, 0.4138562084068e+01, 0.3930209696940e+01,
      0.4325493872224e-10, 0.1330845906564e+01, 0.7632943190217e+01,
      0.5673781176748e-10, 0.2558752615657e+01, 0.5729506548653e+01,
      0.3961436642503e-10, 0.2728071734630e+01, 0.7234794171227e+01,
      0.5101868209058e-10, 0.4113444965144e+01, 0.6836645152238e+01,
      0.5257043167676e-10, 0.6195089830590e+01, 0.8031092209206e+01,
      0.5076613989393e-10, 0.2305124132918e+01, 0.7477522907414e+01,
      0.3342169352778e-10, 0.5415998155071e+01, 0.1097707878456e+02,
      0.3545881983591e-10, 0.3727160564574e+01, 0.4164311961999e+01,
      0.3364063738599e-10, 0.2901121049204e+00, 0.1137170464392e+02,

      0.3357039670776e-10, 0.1652229354331e+01, 0.5223693906222e+01,
      0.4307412268687e-10, 0.4938909587445e+01, 0.1592596075957e+01,
      0.3405769115435e-10, 0.2408890766511e+01, 0.3128388763578e+01,
      0.3001926198480e-10, 0.4862239006386e+01, 0.1748016358760e+01,
      0.2778264787325e-10, 0.5241168661353e+01, 0.7342457794669e+01,
      0.2676159480666e-10, 0.3423593942199e+01, 0.2146165377750e+01,
      0.2954273399939e-10, 0.1881721265406e+01, 0.5368044267797e+00,
      0.3309362888795e-10, 0.1931525677349e+01, 0.8018209333619e+00,
      0.2810283608438e-10, 0.2414659495050e+01, 0.5225775174439e+00,
      0.3378045637764e-10, 0.4238019163430e+01, 0.1554202828031e+00,

      0.2558134979840e-10, 0.1828225235805e+01, 0.5230807360890e+01,
      0.2273755578447e-10, 0.5858184283998e+01, 0.7084896783808e+01,
      0.2294176037690e-10, 0.4514589779057e+01, 0.1726015463500e+02,
      0.2533506099435e-10, 0.2355717851551e+01, 0.5216580451554e+01,
      0.2716685375812e-10, 0.2221003625100e+01, 0.8635942003952e+01,
      0.2419043435198e-10, 0.5955704951635e+01, 0.4690479774488e+01,
      0.2521232544812e-10, 0.1395676848521e+01, 0.5481254917084e+01,
      0.2630195021491e-10, 0.5727468918743e+01, 0.2629832328990e-01,
      0.2548395840944e-10, 0.2628351859400e-03, 0.1349867339771e+01 };

/* Sun-to-Earth, T^1, Y */
   static const double e1y[] = {
      0.9304690546528e-06, 0.0000000000000e+00, 0.0000000000000e+00,
      0.5150715570663e-06, 0.4431807116294e+01, 0.1256615170089e+02,
      0.1290825411056e-07, 0.4388610039678e+01, 0.1884922755134e+02,
      0.4645466665386e-08, 0.5827263376034e+01, 0.6283075850446e+01,
      0.2079625310718e-08, 0.1621698662282e+00, 0.6279552690824e+01,
      0.2078189850907e-08, 0.3344713435140e+01, 0.6286599010068e+01,
      0.6207190138027e-09, 0.5074049319576e+01, 0.4705732307012e+01,
      0.5989826532569e-09, 0.2231842216620e+01, 0.6256777527156e+01,
      0.5961360812618e-09, 0.1274975769045e+01, 0.6309374173736e+01,
      0.4874165471016e-09, 0.3642277426779e+01, 0.7755226100720e+00,

      0.4283834034360e-09, 0.5148765510106e+01, 0.1059381944224e+01,
      0.4652389287529e-09, 0.4715794792175e+01, 0.7860419393880e+01,
      0.3751707476401e-09, 0.6617207370325e+00, 0.5753384878334e+01,
      0.3559998806198e-09, 0.6155548875404e+01, 0.5884926831456e+01,
      0.3558447558857e-09, 0.2898827297664e+01, 0.6812766822558e+01,
      0.3211116927106e-09, 0.3625813502509e+01, 0.6681224869435e+01,
      0.2875609914672e-09, 0.4345435813134e+01, 0.2513230340178e+02,
      0.2843109704069e-09, 0.5862263940038e+01, 0.6127655567643e+01,
      0.2744676468427e-09, 0.3926419475089e+01, 0.6438496133249e+01,
      0.2481285237789e-09, 0.1351976572828e+01, 0.5486777812467e+01,

      0.2060338481033e-09, 0.2147556998591e+01, 0.7079373888424e+01,
      0.2015822358331e-09, 0.4408358972216e+01, 0.6290189305114e+01,
      0.2001195944195e-09, 0.5385829822531e+01, 0.6275962395778e+01,
      0.1953667642377e-09, 0.1304933746120e+01, 0.5507553240374e+01,
      0.1839744078713e-09, 0.6173567228835e+01, 0.1179062909082e+02,
      0.1643334294845e-09, 0.4635942997523e+01, 0.1150676975667e+02,
      0.1768051018652e-09, 0.5086283558874e+01, 0.7113454667900e-02,
      0.1674874205489e-09, 0.2243332137241e+01, 0.7058598460518e+01,
      0.1421445397609e-09, 0.6186899771515e+01, 0.7962980379786e+00,
      0.1255163958267e-09, 0.5730238465658e+01, 0.4694002934110e+01,

      0.1013945281961e-09, 0.1726055228402e+01, 0.3738761453707e+01,
      0.1047294335852e-09, 0.2658801228129e+01, 0.6282095334605e+01,
      0.1047103879392e-09, 0.8481047835035e+00, 0.6284056366286e+01,
      0.9530343962826e-10, 0.3079267149859e+01, 0.6069776770667e+01,
      0.9604637611690e-10, 0.3258679792918e+00, 0.4136910472696e+01,
      0.9153518537177e-10, 0.4398599886584e+00, 0.6496374930224e+01,
      0.8562458214922e-10, 0.4772686794145e+01, 0.1194447056968e+01,
      0.8232525360654e-10, 0.5966220721679e+01, 0.1589072916335e+01,
      0.6150223411438e-10, 0.1780985591923e+01, 0.8827390247185e+01,
      0.6272087858000e-10, 0.3184305429012e+01, 0.8429241228195e+01,

      0.5540476311040e-10, 0.3801260595433e+01, 0.4933208510675e+01,
      0.7331901699361e-10, 0.5205948591865e+01, 0.4535059491685e+01,
      0.6018528702791e-10, 0.4770139083623e+01, 0.1255903824622e+02,
      0.5150530724804e-10, 0.3574796899585e+01, 0.1176985366291e+02,
      0.6471933741811e-10, 0.2679787266521e+01, 0.5088628793478e+01,
      0.5317460644174e-10, 0.9528763345494e+00, 0.3154687086868e+01,
      0.4832187748783e-10, 0.5329322498232e+01, 0.6040347114260e+01,
      0.4716763555110e-10, 0.2395235316466e+01, 0.5331357529664e+01,
      0.4871509139861e-10, 0.3056663648823e+01, 0.1256967486051e+02,
      0.4598417696768e-10, 0.4452762609019e+01, 0.6525804586632e+01,

      0.5674189533175e-10, 0.9879680872193e+00, 0.5729506548653e+01,
      0.4073560328195e-10, 0.5939127696986e+01, 0.7632943190217e+01,
      0.5040994945359e-10, 0.4549875824510e+01, 0.8031092209206e+01,
      0.5078185134679e-10, 0.7346659893982e+00, 0.7477522907414e+01,
      0.3769343537061e-10, 0.1071317188367e+01, 0.7234794171227e+01,
      0.4980331365299e-10, 0.2500345341784e+01, 0.6836645152238e+01,
      0.3458236594757e-10, 0.3825159450711e+01, 0.1097707878456e+02,
      0.3578859493602e-10, 0.5299664791549e+01, 0.4164311961999e+01,
      0.3370504646419e-10, 0.5002316301593e+01, 0.1137170464392e+02,
      0.3299873338428e-10, 0.2526123275282e+01, 0.3930209696940e+01,

      0.4304917318409e-10, 0.3368078557132e+01, 0.1592596075957e+01,
      0.3402418753455e-10, 0.8385495425800e+00, 0.3128388763578e+01,
      0.2778460572146e-10, 0.3669905203240e+01, 0.7342457794669e+01,
      0.2782710128902e-10, 0.2691664812170e+00, 0.1748016358760e+01,
      0.2711725179646e-10, 0.4707487217718e+01, 0.5296909721118e+00,
      0.2981760946340e-10, 0.3190260867816e+00, 0.5368044267797e+00,
      0.2811672977772e-10, 0.3196532315372e+01, 0.7084896783808e+01,
      0.2863454474467e-10, 0.2263240324780e+00, 0.5223693906222e+01,
      0.3333464634051e-10, 0.3498451685065e+01, 0.8018209333619e+00,
      0.3312991747609e-10, 0.5839154477412e+01, 0.1554202828031e+00,

      0.2813255564006e-10, 0.8268044346621e+00, 0.5225775174439e+00,
      0.2665098083966e-10, 0.3934021725360e+01, 0.5216580451554e+01,
      0.2349795705216e-10, 0.5197620913779e+01, 0.2146165377750e+01,
      0.2330352293961e-10, 0.2984999231807e+01, 0.1726015463500e+02,
      0.2728001683419e-10, 0.6521679638544e+00, 0.8635942003952e+01,
      0.2484061007669e-10, 0.3468955561097e+01, 0.5230807360890e+01,
      0.2646328768427e-10, 0.1013724533516e+01, 0.2629832328990e-01,
      0.2518630264831e-10, 0.6108081057122e+01, 0.5481254917084e+01,
      0.2421901455384e-10, 0.1651097776260e+01, 0.1349867339771e+01,
      0.6348533267831e-11, 0.3220226560321e+01, 0.8433466158131e+02 };

/* Sun-to-Earth, T^1, Z */
   static const double e1z[] = {
      0.2278290449966e-05, 0.3413716033863e+01, 0.6283075850446e+01,
      0.5429458209830e-07, 0.0000000000000e+00, 0.0000000000000e+00,
      0.1903240492525e-07, 0.3370592358297e+01, 0.1256615170089e+02,
      0.2385409276743e-09, 0.3327914718416e+01, 0.1884922755134e+02,
      0.8676928342573e-10, 0.1824006811264e+01, 0.5223693906222e+01,
      0.7765442593544e-10, 0.3888564279247e+01, 0.5507553240374e+01,
      0.7066158332715e-10, 0.5194267231944e+01, 0.2352866153506e+01,
      0.7092175288657e-10, 0.2333246960021e+01, 0.8399684731857e+02,
      0.5357582213535e-10, 0.2224031176619e+01, 0.5296909721118e+00,
      0.3828035865021e-10, 0.2156710933584e+01, 0.6279552690824e+01,

      0.3824857220427e-10, 0.1529755219915e+01, 0.6286599010068e+01,
      0.3286995181628e-10, 0.4879512900483e+01, 0.1021328554739e+02 };

/* Sun-to-Earth, T^2, X */
   static const double e2x[] = {
     -0.4143818297913e-10, 0.0000000000000e+00, 0.0000000000000e+00,
      0.2171497694435e-10, 0.4398225628264e+01, 0.1256615170089e+02,
      0.9845398442516e-11, 0.2079720838384e+00, 0.6283075850446e+01,
      0.9256833552682e-12, 0.4191264694361e+01, 0.1884922755134e+02,
      0.1022049384115e-12, 0.5381133195658e+01, 0.8399684731857e+02 };

/* Sun-to-Earth, T^2, Y */
   static const double e2y[] = {
      0.5063375872532e-10, 0.0000000000000e+00, 0.0000000000000e+00,
      0.2173815785980e-10, 0.2827805833053e+01, 0.1256615170089e+02,
      0.1010231999920e-10, 0.4634612377133e+01, 0.6283075850446e+01,
      0.9259745317636e-12, 0.2620612076189e+01, 0.1884922755134e+02,
      0.1022202095812e-12, 0.3809562326066e+01, 0.8399684731857e+02 };

/* Sun-to-Earth, T^2, Z */
   static const double e2z[] = {
      0.9722666114891e-10, 0.5152219582658e+01, 0.6283075850446e+01,
     -0.3494819171909e-11, 0.0000000000000e+00, 0.0000000000000e+00,
      0.6713034376076e-12, 0.6440188750495e+00, 0.1256615170089e+02 };

/* SSB-to-Sun, T^0, X */
   static const double s0x[] = {
      0.4956757536410e-02, 0.3741073751789e+01, 0.5296909721118e+00,
      0.2718490072522e-02, 0.4016011511425e+01, 0.2132990797783e+00,
      0.1546493974344e-02, 0.2170528330642e+01, 0.3813291813120e-01,
      0.8366855276341e-03, 0.2339614075294e+01, 0.7478166569050e-01,
      0.2936777942117e-03, 0.0000000000000e+00, 0.0000000000000e+00,
      0.1201317439469e-03, 0.4090736353305e+01, 0.1059381944224e+01,
      0.7578550887230e-04, 0.3241518088140e+01, 0.4265981595566e+00,
      0.1941787367773e-04, 0.1012202064330e+01, 0.2061856251104e+00,
      0.1889227765991e-04, 0.3892520416440e+01, 0.2204125344462e+00,
      0.1937896968613e-04, 0.4797779441161e+01, 0.1495633313810e+00,

      0.1434506110873e-04, 0.3868960697933e+01, 0.5225775174439e+00,
      0.1406659911580e-04, 0.4759766557397e+00, 0.5368044267797e+00,
      0.1179022300202e-04, 0.7774961520598e+00, 0.7626583626240e-01,
      0.8085864460959e-05, 0.3254654471465e+01, 0.3664874755930e-01,
      0.7622752967615e-05, 0.4227633103489e+01, 0.3961708870310e-01,
      0.6209171139066e-05, 0.2791828325711e+00, 0.7329749511860e-01,
      0.4366435633970e-05, 0.4440454875925e+01, 0.1589072916335e+01,
      0.3792124889348e-05, 0.5156393842356e+01, 0.7113454667900e-02,
      0.3154548963402e-05, 0.6157005730093e+01, 0.4194847048887e+00,
      0.3088359882942e-05, 0.2494567553163e+01, 0.6398972393349e+00,

      0.2788440902136e-05, 0.4934318747989e+01, 0.1102062672231e+00,
      0.3039928456376e-05, 0.4895077702640e+01, 0.6283075850446e+01,
      0.2272258457679e-05, 0.5278394064764e+01, 0.1030928125552e+00,
      0.2162007057957e-05, 0.5802978019099e+01, 0.3163918923335e+00,
      0.1767632855737e-05, 0.3415346595193e-01, 0.1021328554739e+02,
      0.1349413459362e-05, 0.2001643230755e+01, 0.1484170571900e-02,
      0.1170141900476e-05, 0.2424750491620e+01, 0.6327837846670e+00,
      0.1054355266820e-05, 0.3123311487576e+01, 0.4337116142245e+00,
      0.9800822461610e-06, 0.3026258088130e+01, 0.1052268489556e+01,
      0.1091203749931e-05, 0.3157811670347e+01, 0.1162474756779e+01,

      0.6960236715913e-06, 0.8219570542313e+00, 0.1066495398892e+01,
      0.5689257296909e-06, 0.1323052375236e+01, 0.9491756770005e+00,
      0.6613172135802e-06, 0.2765348881598e+00, 0.8460828644453e+00,
      0.6277702517571e-06, 0.5794064466382e+01, 0.1480791608091e+00,
      0.6304884066699e-06, 0.7323555380787e+00, 0.2243449970715e+00,
      0.4897850467382e-06, 0.3062464235399e+01, 0.3340612434717e+01,
      0.3759148598786e-06, 0.4588290469664e+01, 0.3516457698740e-01,
      0.3110520548195e-06, 0.1374299536572e+01, 0.6373574839730e-01,
      0.3064708359780e-06, 0.4222267485047e+01, 0.1104591729320e-01,
      0.2856347168241e-06, 0.3714202944973e+01, 0.1510475019529e+00,

      0.2840945514288e-06, 0.2847972875882e+01, 0.4110125927500e-01,
      0.2378951599405e-06, 0.3762072563388e+01, 0.2275259891141e+00,
      0.2714229481417e-06, 0.1036049980031e+01, 0.2535050500000e-01,
      0.2323551717307e-06, 0.4682388599076e+00, 0.8582758298370e-01,
      0.1881790512219e-06, 0.4790565425418e+01, 0.2118763888447e+01,
      0.2261353968371e-06, 0.1669144912212e+01, 0.7181332454670e-01,
      0.2214546389848e-06, 0.3937717281614e+01, 0.2968341143800e-02,
      0.2184915594933e-06, 0.1129169845099e+00, 0.7775000683430e-01,
      0.2000164937936e-06, 0.4030009638488e+01, 0.2093666171530e+00,
      0.1966105136719e-06, 0.8745955786834e+00, 0.2172315424036e+00,

      0.1904742332624e-06, 0.5919743598964e+01, 0.2022531624851e+00,
      0.1657399705031e-06, 0.2549141484884e+01, 0.7358765972222e+00,
      0.1574070533987e-06, 0.5277533020230e+01, 0.7429900518901e+00,
      0.1832261651039e-06, 0.3064688127777e+01, 0.3235053470014e+00,
      0.1733615346569e-06, 0.3011432799094e+01, 0.1385174140878e+00,
      0.1549124014496e-06, 0.4005569132359e+01, 0.5154640627760e+00,
      0.1637044713838e-06, 0.1831375966632e+01, 0.8531963191132e+00,
      0.1123420082383e-06, 0.1180270407578e+01, 0.1990721704425e+00,
      0.1083754165740e-06, 0.3414101320863e+00, 0.5439178814476e+00,
      0.1156638012655e-06, 0.6130479452594e+00, 0.5257585094865e+00,

      0.1142548785134e-06, 0.3724761948846e+01, 0.5336234347371e+00,
      0.7921463895965e-07, 0.2435425589361e+01, 0.1478866649112e+01,
      0.7428600285231e-07, 0.3542144398753e+01, 0.2164800718209e+00,
      0.8323211246747e-07, 0.3525058072354e+01, 0.1692165728891e+01,
      0.7257595116312e-07, 0.1364299431982e+01, 0.2101180877357e+00,
      0.7111185833236e-07, 0.2460478875808e+01, 0.4155522422634e+00,
      0.6868090383716e-07, 0.4397327670704e+01, 0.1173197218910e+00,
      0.7226419974175e-07, 0.4042647308905e+01, 0.1265567569334e+01,
      0.6955642383177e-07, 0.2865047906085e+01, 0.9562891316684e+00,
      0.7492139296331e-07, 0.5014278994215e+01, 0.1422690933580e-01,

      0.6598363128857e-07, 0.2376730020492e+01, 0.6470106940028e+00,
      0.7381147293385e-07, 0.3272990384244e+01, 0.1581959461667e+01,
      0.6402909624032e-07, 0.5302290955138e+01, 0.9597935788730e-01,
      0.6237454263857e-07, 0.5444144425332e+01, 0.7084920306520e-01,
      0.5241198544016e-07, 0.4215359579205e+01, 0.5265099800692e+00,
      0.5144463853918e-07, 0.1218916689916e+00, 0.5328719641544e+00,
      0.5868164772299e-07, 0.2369402002213e+01, 0.7871412831580e-01,
      0.6233195669151e-07, 0.1254922242403e+01, 0.2608790314060e+02,
      0.6068463791422e-07, 0.5679713760431e+01, 0.1114304132498e+00,
      0.4359361135065e-07, 0.6097219641646e+00, 0.1375773836557e+01,

      0.4686510366826e-07, 0.4786231041431e+01, 0.1143987543936e+00,
      0.3758977287225e-07, 0.1167368068139e+01, 0.1596186371003e+01,
      0.4282051974778e-07, 0.1519471064319e+01, 0.2770348281756e+00,
      0.5153765386113e-07, 0.1860532322984e+01, 0.2228608264996e+00,
      0.4575129387188e-07, 0.7632857887158e+00, 0.1465949902372e+00,
      0.3326844933286e-07, 0.1298219485285e+01, 0.5070101000000e-01,
      0.3748617450984e-07, 0.1046510321062e+01, 0.4903339079539e+00,
      0.2816756661499e-07, 0.3434522346190e+01, 0.2991266627620e+00,
      0.3412750405039e-07, 0.2523766270318e+01, 0.3518164938661e+00,
      0.2655796761776e-07, 0.2904422260194e+01, 0.6256703299991e+00,

      0.2963597929458e-07, 0.5923900431149e+00, 0.1099462426779e+00,
      0.2539523734781e-07, 0.4851947722567e+01, 0.1256615170089e+02,
      0.2283087914139e-07, 0.3400498595496e+01, 0.6681224869435e+01,
      0.2321309799331e-07, 0.5789099148673e+01, 0.3368040641550e-01,
      0.2549657649750e-07, 0.3991856479792e-01, 0.1169588211447e+01,
      0.2290462303977e-07, 0.2788567577052e+01, 0.1045155034888e+01,
      0.1945398522914e-07, 0.3290896998176e+01, 0.1155361302111e+01,
      0.1849171512638e-07, 0.2698060129367e+01, 0.4452511715700e-02,
      0.1647199834254e-07, 0.3016735644085e+01, 0.4408250688924e+00,
      0.1529530765273e-07, 0.5573043116178e+01, 0.6521991896920e-01,

      0.1433199339978e-07, 0.1481192356147e+01, 0.9420622223326e+00,
      0.1729134193602e-07, 0.1422817538933e+01, 0.2108507877249e+00,
      0.1716463931346e-07, 0.3469468901855e+01, 0.2157473718317e+00,
      0.1391206061378e-07, 0.6122436220547e+01, 0.4123712502208e+00,
      0.1404746661924e-07, 0.1647765641936e+01, 0.4258542984690e-01,
      0.1410452399455e-07, 0.5989729161964e+01, 0.2258291676434e+00,
      0.1089828772168e-07, 0.2833705509371e+01, 0.4226656969313e+00,
      0.1047374564948e-07, 0.5090690007331e+00, 0.3092784376656e+00,
      0.1358279126532e-07, 0.5128990262836e+01, 0.7923417740620e-01,
      0.1020456476148e-07, 0.9632772880808e+00, 0.1456308687557e+00,

      0.1033428735328e-07, 0.3223779318418e+01, 0.1795258541446e+01,
      0.1412435841540e-07, 0.2410271572721e+01, 0.1525316725248e+00,
      0.9722759371574e-08, 0.2333531395690e+01, 0.8434341241180e-01,
      0.9657334084704e-08, 0.6199270974168e+01, 0.1272681024002e+01,
      0.1083641148690e-07, 0.2864222292929e+01, 0.7032915397480e-01,
      0.1067318403838e-07, 0.5833458866568e+00, 0.2123349582968e+00,
      0.1062366201976e-07, 0.4307753989494e+01, 0.2142632012598e+00,
      0.1236364149266e-07, 0.2873917870593e+01, 0.1847279083684e+00,
      0.1092759489593e-07, 0.2959887266733e+01, 0.1370332435159e+00,
      0.8912069362899e-08, 0.5141213702562e+01, 0.2648454860559e+01,

      0.9656467707970e-08, 0.4532182462323e+01, 0.4376440768498e+00,
      0.8098386150135e-08, 0.2268906338379e+01, 0.2880807454688e+00,
      0.7857714675000e-08, 0.4055544260745e+01, 0.2037373330570e+00,
      0.7288455940646e-08, 0.5357901655142e+01, 0.1129145838217e+00,
      0.9450595950552e-08, 0.4264926963939e+01, 0.5272426800584e+00,
      0.9381718247537e-08, 0.7489366976576e-01, 0.5321392641652e+00,
      0.7079052646038e-08, 0.1923311052874e+01, 0.6288513220417e+00,
      0.9259004415344e-08, 0.2970256853438e+01, 0.1606092486742e+00,
      0.8259801499742e-08, 0.3327056314697e+01, 0.8389694097774e+00,
      0.6476334355779e-08, 0.2954925505727e+01, 0.2008557621224e+01,

      0.5984021492007e-08, 0.9138753105829e+00, 0.2042657109477e+02,
      0.5989546863181e-08, 0.3244464082031e+01, 0.2111650433779e+01,
      0.6233108606023e-08, 0.4995232638403e+00, 0.4305306221819e+00,
      0.6877299149965e-08, 0.2834987233449e+01, 0.9561746721300e-02,
      0.8311234227190e-08, 0.2202951835758e+01, 0.3801276407308e+00,
      0.6599472832414e-08, 0.4478581462618e+01, 0.1063314406849e+01,
      0.6160491096549e-08, 0.5145858696411e+01, 0.1368660381889e+01,
      0.6164772043891e-08, 0.3762976697911e+00, 0.4234171675140e+00,
      0.6363248684450e-08, 0.3162246718685e+01, 0.1253008786510e-01,
      0.6448587520999e-08, 0.3442693302119e+01, 0.5287268506303e+00,

      0.6431662283977e-08, 0.8977549136606e+00, 0.5306550935933e+00,
      0.6351223158474e-08, 0.4306447410369e+01, 0.5217580628120e+02,
      0.5476721393451e-08, 0.3888529177855e+01, 0.2221856701002e+01,
      0.5341772572619e-08, 0.2655560662512e+01, 0.7466759693650e-01,
      0.5337055758302e-08, 0.5164990735946e+01, 0.7489573444450e-01,
      0.5373120816787e-08, 0.6041214553456e+01, 0.1274714967946e+00,
      0.5392351705426e-08, 0.9177763485932e+00, 0.1055449481598e+01,
      0.6688495850205e-08, 0.3089608126937e+01, 0.2213766559277e+00,
      0.5072003660362e-08, 0.4311316541553e+01, 0.2132517061319e+00,
      0.5070726650455e-08, 0.5790675464444e+00, 0.2133464534247e+00,

      0.5658012950032e-08, 0.2703945510675e+01, 0.7287631425543e+00,
      0.4835509924854e-08, 0.2975422976065e+01, 0.7160067364790e-01,
      0.6479821978012e-08, 0.1324168733114e+01, 0.2209183458640e-01,
      0.6230636494980e-08, 0.2860103632836e+01, 0.3306188016693e+00,
      0.4649239516213e-08, 0.4832259763403e+01, 0.7796265773310e-01,
      0.6487325792700e-08, 0.2726165825042e+01, 0.3884652414254e+00,
      0.4682823682770e-08, 0.6966602455408e+00, 0.1073608853559e+01,
      0.5704230804976e-08, 0.5669634104606e+01, 0.8731175355560e-01,
      0.6125413585489e-08, 0.1513386538915e+01, 0.7605151500000e-01,
      0.6035825038187e-08, 0.1983509168227e+01, 0.9846002785331e+00,

      0.4331123462303e-08, 0.2782892992807e+01, 0.4297791515992e+00,
      0.4681107685143e-08, 0.5337232886836e+01, 0.2127790306879e+00,
      0.4669105829655e-08, 0.5837133792160e+01, 0.2138191288687e+00,
      0.5138823602365e-08, 0.3080560200507e+01, 0.7233337363710e-01,
      0.4615856664534e-08, 0.1661747897471e+01, 0.8603097737811e+00,
      0.4496916702197e-08, 0.2112508027068e+01, 0.7381754420900e-01,
      0.4278479042945e-08, 0.5716528462627e+01, 0.7574578717200e-01,
      0.3840525503932e-08, 0.6424172726492e+00, 0.3407705765729e+00,
      0.4866636509685e-08, 0.4919244697715e+01, 0.7722995774390e-01,
      0.3526100639296e-08, 0.2550821052734e+01, 0.6225157782540e-01,

      0.3939558488075e-08, 0.3939331491710e+01, 0.5268983110410e-01,
      0.4041268772576e-08, 0.2275337571218e+01, 0.3503323232942e+00,
      0.3948761842853e-08, 0.1999324200790e+01, 0.1451108196653e+00,
      0.3258394550029e-08, 0.9121001378200e+00, 0.5296435984654e+00,
      0.3257897048761e-08, 0.3428428660869e+01, 0.5297383457582e+00,
      0.3842559031298e-08, 0.6132927720035e+01, 0.9098186128426e+00,
      0.3109920095448e-08, 0.7693650193003e+00, 0.3932462625300e-02,
      0.3132237775119e-08, 0.3621293854908e+01, 0.2346394437820e+00,
      0.3942189421510e-08, 0.4841863659733e+01, 0.3180992042600e-02,
      0.3796972285340e-08, 0.1814174994268e+01, 0.1862120789403e+00,

      0.3995640233688e-08, 0.1386990406091e+01, 0.4549093064213e+00,
      0.2875013727414e-08, 0.9178318587177e+00, 0.1905464808669e+01,
      0.3073719932844e-08, 0.2688923811835e+01, 0.3628624111593e+00,
      0.2731016580075e-08, 0.1188259127584e+01, 0.2131850110243e+00,
      0.2729549896546e-08, 0.3702160634273e+01, 0.2134131485323e+00,
      0.3339372892449e-08, 0.7199163960331e+00, 0.2007689919132e+00,
      0.2898833764204e-08, 0.1916709364999e+01, 0.5291709230214e+00,
      0.2894536549362e-08, 0.2424043195547e+01, 0.5302110212022e+00,
      0.3096872473843e-08, 0.4445894977497e+01, 0.2976424921901e+00,
      0.2635672326810e-08, 0.3814366984117e+01, 0.1485980103780e+01,

      0.3649302697001e-08, 0.2924200596084e+01, 0.6044726378023e+00,
      0.3127954585895e-08, 0.1842251648327e+01, 0.1084620721060e+00,
      0.2616040173947e-08, 0.4155841921984e+01, 0.1258454114666e+01,
      0.2597395859860e-08, 0.1158045978874e+00, 0.2103781122809e+00,
      0.2593286172210e-08, 0.4771850408691e+01, 0.2162200472757e+00,
      0.2481823585747e-08, 0.4608842558889e+00, 0.1062562936266e+01,
      0.2742219550725e-08, 0.1538781127028e+01, 0.5651155736444e+00,
      0.3199558469610e-08, 0.3226647822878e+00, 0.7036329877322e+00,
      0.2666088542957e-08, 0.1967991731219e+00, 0.1400015846597e+00,
      0.2397067430580e-08, 0.3707036669873e+01, 0.2125476091956e+00,

      0.2376570772738e-08, 0.1182086628042e+01, 0.2140505503610e+00,
      0.2547228007887e-08, 0.4906256820629e+01, 0.1534957940063e+00,
      0.2265575594114e-08, 0.3414949866857e+01, 0.2235935264888e+00,
      0.2464381430585e-08, 0.4599122275378e+01, 0.2091065926078e+00,
      0.2433408527044e-08, 0.2830751145445e+00, 0.2174915669488e+00,
      0.2443605509076e-08, 0.4212046432538e+01, 0.1739420156204e+00,
      0.2319779262465e-08, 0.9881978408630e+00, 0.7530171478090e-01,
      0.2284622835465e-08, 0.5565347331588e+00, 0.7426161660010e-01,
      0.2467268750783e-08, 0.5655708150766e+00, 0.2526561439362e+00,
      0.2808513492782e-08, 0.1418405053408e+01, 0.5636314030725e+00,

      0.2329528932532e-08, 0.4069557545675e+01, 0.1056200952181e+01,
      0.9698639532817e-09, 0.1074134313634e+01, 0.7826370942180e+02 };

/* SSB-to-Sun, T^0, Y */
   static const double s0y[] = {
      0.4955392320126e-02, 0.2170467313679e+01, 0.5296909721118e+00,
      0.2722325167392e-02, 0.2444433682196e+01, 0.2132990797783e+00,
      0.1546579925346e-02, 0.5992779281546e+00, 0.3813291813120e-01,
      0.8363140252966e-03, 0.7687356310801e+00, 0.7478166569050e-01,
      0.3385792683603e-03, 0.0000000000000e+00, 0.0000000000000e+00,
      0.1201192221613e-03, 0.2520035601514e+01, 0.1059381944224e+01,
      0.7587125720554e-04, 0.1669954006449e+01, 0.4265981595566e+00,
      0.1964155361250e-04, 0.5707743963343e+01, 0.2061856251104e+00,
      0.1891900364909e-04, 0.2320960679937e+01, 0.2204125344462e+00,
      0.1937373433356e-04, 0.3226940689555e+01, 0.1495633313810e+00,

      0.1437139941351e-04, 0.2301626908096e+01, 0.5225775174439e+00,
      0.1406267683099e-04, 0.5188579265542e+01, 0.5368044267797e+00,
      0.1178703080346e-04, 0.5489483248476e+01, 0.7626583626240e-01,
      0.8079835186041e-05, 0.1683751835264e+01, 0.3664874755930e-01,
      0.7623253594652e-05, 0.2656400462961e+01, 0.3961708870310e-01,
      0.6248667483971e-05, 0.4992775362055e+01, 0.7329749511860e-01,
      0.4366353695038e-05, 0.2869706279678e+01, 0.1589072916335e+01,
      0.3829101568895e-05, 0.3572131359950e+01, 0.7113454667900e-02,
      0.3175733773908e-05, 0.4535372530045e+01, 0.4194847048887e+00,
      0.3092437902159e-05, 0.9230153317909e+00, 0.6398972393349e+00,

      0.2874168812154e-05, 0.3363143761101e+01, 0.1102062672231e+00,
      0.3040119321826e-05, 0.3324250895675e+01, 0.6283075850446e+01,
      0.2699723308006e-05, 0.2917882441928e+00, 0.1030928125552e+00,
      0.2134832683534e-05, 0.4220997202487e+01, 0.3163918923335e+00,
      0.1770412139433e-05, 0.4747318496462e+01, 0.1021328554739e+02,
      0.1377264209373e-05, 0.4305058462401e+00, 0.1484170571900e-02,
      0.1127814538960e-05, 0.8538177240740e+00, 0.6327837846670e+00,
      0.1055608090130e-05, 0.1551800742580e+01, 0.4337116142245e+00,
      0.9802673861420e-06, 0.1459646735377e+01, 0.1052268489556e+01,
      0.1090329461951e-05, 0.1587351228711e+01, 0.1162474756779e+01,

      0.6959590025090e-06, 0.5534442628766e+01, 0.1066495398892e+01,
      0.5664914529542e-06, 0.6030673003297e+01, 0.9491756770005e+00,
      0.6607787763599e-06, 0.4989507233927e+01, 0.8460828644453e+00,
      0.6269725742838e-06, 0.4222951804572e+01, 0.1480791608091e+00,
      0.6301889697863e-06, 0.5444316669126e+01, 0.2243449970715e+00,
      0.4891042662861e-06, 0.1490552839784e+01, 0.3340612434717e+01,
      0.3457083123290e-06, 0.3030475486049e+01, 0.3516457698740e-01,
      0.3032559967314e-06, 0.2652038793632e+01, 0.1104591729320e-01,
      0.2841133988903e-06, 0.1276744786829e+01, 0.4110125927500e-01,
      0.2855564444432e-06, 0.2143368674733e+01, 0.1510475019529e+00,

      0.2765157135038e-06, 0.5444186109077e+01, 0.6373574839730e-01,
      0.2382312465034e-06, 0.2190521137593e+01, 0.2275259891141e+00,
      0.2808060365077e-06, 0.5735195064841e+01, 0.2535050500000e-01,
      0.2332175234405e-06, 0.9481985524859e-01, 0.7181332454670e-01,
      0.2322488199659e-06, 0.5180499361533e+01, 0.8582758298370e-01,
      0.1881850258423e-06, 0.3219788273885e+01, 0.2118763888447e+01,
      0.2196111392808e-06, 0.2366941159761e+01, 0.2968341143800e-02,
      0.2183810335519e-06, 0.4825445110915e+01, 0.7775000683430e-01,
      0.2002733093326e-06, 0.2457148995307e+01, 0.2093666171530e+00,
      0.1967111767229e-06, 0.5586291545459e+01, 0.2172315424036e+00,

      0.1568473250543e-06, 0.3708003123320e+01, 0.7429900518901e+00,
      0.1852528314300e-06, 0.4310638151560e+01, 0.2022531624851e+00,
      0.1832111226447e-06, 0.1494665322656e+01, 0.3235053470014e+00,
      0.1746805502310e-06, 0.1451378500784e+01, 0.1385174140878e+00,
      0.1555730966650e-06, 0.1068040418198e+01, 0.7358765972222e+00,
      0.1554883462559e-06, 0.2442579035461e+01, 0.5154640627760e+00,
      0.1638380568746e-06, 0.2597913420625e+00, 0.8531963191132e+00,
      0.1159938593640e-06, 0.5834512021280e+01, 0.1990721704425e+00,
      0.1083427965695e-06, 0.5054033177950e+01, 0.5439178814476e+00,
      0.1156480369431e-06, 0.5325677432457e+01, 0.5257585094865e+00,

      0.1141308860095e-06, 0.2153403923857e+01, 0.5336234347371e+00,
      0.7913146470946e-07, 0.8642846847027e+00, 0.1478866649112e+01,
      0.7439752463733e-07, 0.1970628496213e+01, 0.2164800718209e+00,
      0.7280277104079e-07, 0.6073307250609e+01, 0.2101180877357e+00,
      0.8319567719136e-07, 0.1954371928334e+01, 0.1692165728891e+01,
      0.7137705549290e-07, 0.8904989440909e+00, 0.4155522422634e+00,
      0.6900825396225e-07, 0.2825717714977e+01, 0.1173197218910e+00,
      0.7245757216635e-07, 0.2481677513331e+01, 0.1265567569334e+01,
      0.6961165696255e-07, 0.1292955312978e+01, 0.9562891316684e+00,
      0.7571804456890e-07, 0.3427517575069e+01, 0.1422690933580e-01,

      0.6605425721904e-07, 0.8052192701492e+00, 0.6470106940028e+00,
      0.7375477357248e-07, 0.1705076390088e+01, 0.1581959461667e+01,
      0.7041664951470e-07, 0.4848356967891e+00, 0.9597935788730e-01,
      0.6322199535763e-07, 0.3878069473909e+01, 0.7084920306520e-01,
      0.5244380279191e-07, 0.2645560544125e+01, 0.5265099800692e+00,
      0.5143125704988e-07, 0.4834486101370e+01, 0.5328719641544e+00,
      0.5871866319373e-07, 0.7981472548900e+00, 0.7871412831580e-01,
      0.6300822573871e-07, 0.5979398788281e+01, 0.2608790314060e+02,
      0.6062154271548e-07, 0.4108655402756e+01, 0.1114304132498e+00,
      0.4361912339976e-07, 0.5322624319280e+01, 0.1375773836557e+01,

      0.4417005920067e-07, 0.6240817359284e+01, 0.2770348281756e+00,
      0.4686806749936e-07, 0.3214977301156e+01, 0.1143987543936e+00,
      0.3758892132305e-07, 0.5879809634765e+01, 0.1596186371003e+01,
      0.5151351332319e-07, 0.2893377688007e+00, 0.2228608264996e+00,
      0.4554683578572e-07, 0.5475427144122e+01, 0.1465949902372e+00,
      0.3442381385338e-07, 0.5992034796640e+01, 0.5070101000000e-01,
      0.2831093954933e-07, 0.5367350273914e+01, 0.3092784376656e+00,
      0.3756267090084e-07, 0.5758171285420e+01, 0.4903339079539e+00,
      0.2816374679892e-07, 0.1863718700923e+01, 0.2991266627620e+00,
      0.3419307025569e-07, 0.9524347534130e+00, 0.3518164938661e+00,

      0.2904250494239e-07, 0.5304471615602e+01, 0.1099462426779e+00,
      0.2471734511206e-07, 0.1297069793530e+01, 0.6256703299991e+00,
      0.2539620831872e-07, 0.3281126083375e+01, 0.1256615170089e+02,
      0.2281017868007e-07, 0.1829122133165e+01, 0.6681224869435e+01,
      0.2275319473335e-07, 0.5797198160181e+01, 0.3932462625300e-02,
      0.2547755368442e-07, 0.4752697708330e+01, 0.1169588211447e+01,
      0.2285979669317e-07, 0.1223205292886e+01, 0.1045155034888e+01,
      0.1913386560994e-07, 0.1757532993389e+01, 0.1155361302111e+01,
      0.1809020525147e-07, 0.4246116108791e+01, 0.3368040641550e-01,
      0.1649213300201e-07, 0.1445162890627e+01, 0.4408250688924e+00,

      0.1834972793932e-07, 0.1126917567225e+01, 0.4452511715700e-02,
      0.1439550648138e-07, 0.6160756834764e+01, 0.9420622223326e+00,
      0.1487645457041e-07, 0.4358761931792e+01, 0.4123712502208e+00,
      0.1731729516660e-07, 0.6134456753344e+01, 0.2108507877249e+00,
      0.1717747163567e-07, 0.1898186084455e+01, 0.2157473718317e+00,
      0.1418190430374e-07, 0.4180286741266e+01, 0.6521991896920e-01,
      0.1404844134873e-07, 0.7654053565412e-01, 0.4258542984690e-01,
      0.1409842846538e-07, 0.4418612420312e+01, 0.2258291676434e+00,
      0.1090948346291e-07, 0.1260615686131e+01, 0.4226656969313e+00,
      0.1357577323612e-07, 0.3558248818690e+01, 0.7923417740620e-01,

      0.1018154061960e-07, 0.5676087241256e+01, 0.1456308687557e+00,
      0.1412073972109e-07, 0.8394392632422e+00, 0.1525316725248e+00,
      0.1030938326496e-07, 0.1653593274064e+01, 0.1795258541446e+01,
      0.1180081567104e-07, 0.1285802592036e+01, 0.7032915397480e-01,
      0.9708510575650e-08, 0.7631889488106e+00, 0.8434341241180e-01,
      0.9637689663447e-08, 0.4630642649176e+01, 0.1272681024002e+01,
      0.1068910429389e-07, 0.5294934032165e+01, 0.2123349582968e+00,
      0.1063716179336e-07, 0.2736266800832e+01, 0.2142632012598e+00,
      0.1234858713814e-07, 0.1302891146570e+01, 0.1847279083684e+00,
      0.8912631189738e-08, 0.3570415993621e+01, 0.2648454860559e+01,

      0.1036378285534e-07, 0.4236693440949e+01, 0.1370332435159e+00,
      0.9667798501561e-08, 0.2960768892398e+01, 0.4376440768498e+00,
      0.8108314201902e-08, 0.6987781646841e+00, 0.2880807454688e+00,
      0.7648364324628e-08, 0.2499017863863e+01, 0.2037373330570e+00,
      0.7286136828406e-08, 0.3787426951665e+01, 0.1129145838217e+00,
      0.9448237743913e-08, 0.2694354332983e+01, 0.5272426800584e+00,
      0.9374276106428e-08, 0.4787121277064e+01, 0.5321392641652e+00,
      0.7100226287462e-08, 0.3530238792101e+00, 0.6288513220417e+00,
      0.9253056659571e-08, 0.1399478925664e+01, 0.1606092486742e+00,
      0.6636432145504e-08, 0.3479575438447e+01, 0.1368660381889e+01,

      0.6469975312932e-08, 0.1383669964800e+01, 0.2008557621224e+01,
      0.7335849729765e-08, 0.1243698166898e+01, 0.9561746721300e-02,
      0.8743421205855e-08, 0.3776164289301e+01, 0.3801276407308e+00,
      0.5993635744494e-08, 0.5627122113596e+01, 0.2042657109477e+02,
      0.5981008479693e-08, 0.1674336636752e+01, 0.2111650433779e+01,
      0.6188535145838e-08, 0.5214925208672e+01, 0.4305306221819e+00,
      0.6596074017566e-08, 0.2907653268124e+01, 0.1063314406849e+01,
      0.6630815126226e-08, 0.2127643669658e+01, 0.8389694097774e+00,
      0.6156772830040e-08, 0.5082160803295e+01, 0.4234171675140e+00,
      0.6446960563014e-08, 0.1872100916905e+01, 0.5287268506303e+00,

      0.6429324424668e-08, 0.5610276103577e+01, 0.5306550935933e+00,
      0.6302232396465e-08, 0.1592152049607e+01, 0.1253008786510e-01,
      0.6399244436159e-08, 0.2746214421532e+01, 0.5217580628120e+02,
      0.5474965172558e-08, 0.2317666374383e+01, 0.2221856701002e+01,
      0.5339293190692e-08, 0.1084724961156e+01, 0.7466759693650e-01,
      0.5334733683389e-08, 0.3594106067745e+01, 0.7489573444450e-01,
      0.5392665782110e-08, 0.5630254365606e+01, 0.1055449481598e+01,
      0.6682075673789e-08, 0.1518480041732e+01, 0.2213766559277e+00,
      0.5079130495960e-08, 0.2739765115711e+01, 0.2132517061319e+00,
      0.5077759793261e-08, 0.5290711290094e+01, 0.2133464534247e+00,

      0.4832037368310e-08, 0.1404473217200e+01, 0.7160067364790e-01,
      0.6463279674802e-08, 0.6038381695210e+01, 0.2209183458640e-01,
      0.6240592771560e-08, 0.1290170653666e+01, 0.3306188016693e+00,
      0.4672013521493e-08, 0.3261895939677e+01, 0.7796265773310e-01,
      0.6500650750348e-08, 0.1154522312095e+01, 0.3884652414254e+00,
      0.6344161389053e-08, 0.6206111545062e+01, 0.7605151500000e-01,
      0.4682518370646e-08, 0.5409118796685e+01, 0.1073608853559e+01,
      0.5329460015591e-08, 0.1202985784864e+01, 0.7287631425543e+00,
      0.5701588675898e-08, 0.4098715257064e+01, 0.8731175355560e-01,
      0.6030690867211e-08, 0.4132033218460e+00, 0.9846002785331e+00,

      0.4336256312655e-08, 0.1211415991827e+01, 0.4297791515992e+00,
      0.4688498808975e-08, 0.3765479072409e+01, 0.2127790306879e+00,
      0.4675578609335e-08, 0.4265540037226e+01, 0.2138191288687e+00,
      0.4225578112158e-08, 0.5237566010676e+01, 0.3407705765729e+00,
      0.5139422230028e-08, 0.1507173079513e+01, 0.7233337363710e-01,
      0.4619995093571e-08, 0.9023957449848e-01, 0.8603097737811e+00,
      0.4494776255461e-08, 0.5414930552139e+00, 0.7381754420900e-01,
      0.4274026276788e-08, 0.4145735303659e+01, 0.7574578717200e-01,
      0.5018141789353e-08, 0.3344408829055e+01, 0.3180992042600e-02,
      0.4866163952181e-08, 0.3348534657607e+01, 0.7722995774390e-01,

      0.4111986020501e-08, 0.4198823597220e+00, 0.1451108196653e+00,
      0.3356142784950e-08, 0.5609144747180e+01, 0.1274714967946e+00,
      0.4070575554551e-08, 0.7028411059224e+00, 0.3503323232942e+00,
      0.3257451857278e-08, 0.5624697983086e+01, 0.5296435984654e+00,
      0.3256973703026e-08, 0.1857842076707e+01, 0.5297383457582e+00,
      0.3830771508640e-08, 0.4562887279931e+01, 0.9098186128426e+00,
      0.3725024005962e-08, 0.2358058692652e+00, 0.1084620721060e+00,
      0.3136763921756e-08, 0.2049731526845e+01, 0.2346394437820e+00,
      0.3795147256194e-08, 0.2432356296933e+00, 0.1862120789403e+00,
      0.2877342229911e-08, 0.5631101279387e+01, 0.1905464808669e+01,

      0.3076931798805e-08, 0.1117615737392e+01, 0.3628624111593e+00,
      0.2734765945273e-08, 0.5899826516955e+01, 0.2131850110243e+00,
      0.2733405296885e-08, 0.2130562964070e+01, 0.2134131485323e+00,
      0.2898552353410e-08, 0.3462387048225e+00, 0.5291709230214e+00,
      0.2893736103681e-08, 0.8534352781543e+00, 0.5302110212022e+00,
      0.3095717734137e-08, 0.2875061429041e+01, 0.2976424921901e+00,
      0.2636190425832e-08, 0.2242512846659e+01, 0.1485980103780e+01,
      0.3645512095537e-08, 0.1354016903958e+01, 0.6044726378023e+00,
      0.2808173547723e-08, 0.6705114365631e-01, 0.6225157782540e-01,
      0.2625012866888e-08, 0.4775705748482e+01, 0.5268983110410e-01,

      0.2572233995651e-08, 0.2638924216139e+01, 0.1258454114666e+01,
      0.2604238824792e-08, 0.4826358927373e+01, 0.2103781122809e+00,
      0.2596886385239e-08, 0.3200388483118e+01, 0.2162200472757e+00,
      0.3228057304264e-08, 0.5384848409563e+01, 0.2007689919132e+00,
      0.2481601798252e-08, 0.5173373487744e+01, 0.1062562936266e+01,
      0.2745977498864e-08, 0.6250966149853e+01, 0.5651155736444e+00,
      0.2669878833811e-08, 0.4906001352499e+01, 0.1400015846597e+00,
      0.3203986611711e-08, 0.5034333010005e+01, 0.7036329877322e+00,
      0.3354961227212e-08, 0.6108262423137e+01, 0.4549093064213e+00,
      0.2400407324558e-08, 0.2135399294955e+01, 0.2125476091956e+00,

      0.2379905859802e-08, 0.5893721933961e+01, 0.2140505503610e+00,
      0.2550844302187e-08, 0.3331940762063e+01, 0.1534957940063e+00,
      0.2268824211001e-08, 0.1843418461035e+01, 0.2235935264888e+00,
      0.2464700891204e-08, 0.3029548547230e+01, 0.2091065926078e+00,
      0.2436814726024e-08, 0.4994717970364e+01, 0.2174915669488e+00,
      0.2443623894745e-08, 0.2645102591375e+01, 0.1739420156204e+00,
      0.2318701783838e-08, 0.5700547397897e+01, 0.7530171478090e-01,
      0.2284448700256e-08, 0.5268898905872e+01, 0.7426161660010e-01,
      0.2468848123510e-08, 0.5276280575078e+01, 0.2526561439362e+00,
      0.2814052350303e-08, 0.6130168623475e+01, 0.5636314030725e+00,

      0.2243662755220e-08, 0.6631692457995e+00, 0.8886590321940e-01,
      0.2330795855941e-08, 0.2499435487702e+01, 0.1056200952181e+01,
      0.9757679038404e-09, 0.5796846023126e+01, 0.7826370942180e+02 };

/* SSB-to-Sun, T^0, Z */
   static const double s0z[] = {
      0.1181255122986e-03, 0.4607918989164e+00, 0.2132990797783e+00,
      0.1127777651095e-03, 0.4169146331296e+00, 0.5296909721118e+00,
      0.4777754401806e-04, 0.4582657007130e+01, 0.3813291813120e-01,
      0.1129354285772e-04, 0.5758735142480e+01, 0.7478166569050e-01,
     -0.1149543637123e-04, 0.0000000000000e+00, 0.0000000000000e+00,
      0.3298730512306e-05, 0.5978801994625e+01, 0.4265981595566e+00,
      0.2733376706079e-05, 0.7665413691040e+00, 0.1059381944224e+01,
      0.9426389657270e-06, 0.3710201265838e+01, 0.2061856251104e+00,
      0.8187517749552e-06, 0.3390675605802e+00, 0.2204125344462e+00,
      0.4080447871819e-06, 0.4552296640088e+00, 0.5225775174439e+00,

      0.3169973017028e-06, 0.3445455899321e+01, 0.5368044267797e+00,
      0.2438098615549e-06, 0.5664675150648e+01, 0.3664874755930e-01,
      0.2601897517235e-06, 0.1931894095697e+01, 0.1495633313810e+00,
      0.2314558080079e-06, 0.3666319115574e+00, 0.3961708870310e-01,
      0.1962549548002e-06, 0.3167411699020e+01, 0.7626583626240e-01,
      0.2180518287925e-06, 0.1544420746580e+01, 0.7113454667900e-02,
      0.1451382442868e-06, 0.1583756740070e+01, 0.1102062672231e+00,
      0.1358439007389e-06, 0.5239941758280e+01, 0.6398972393349e+00,
      0.1050585898028e-06, 0.2266958352859e+01, 0.3163918923335e+00,
      0.1050029870186e-06, 0.2711495250354e+01, 0.4194847048887e+00,

      0.9934920679800e-07, 0.1116208151396e+01, 0.1589072916335e+01,
      0.1048395331560e-06, 0.3408619600206e+01, 0.1021328554739e+02,
      0.8370147196668e-07, 0.3810459401087e+01, 0.2535050500000e-01,
      0.7989856510998e-07, 0.3769910473647e+01, 0.7329749511860e-01,
      0.5441221655233e-07, 0.2416994903374e+01, 0.1030928125552e+00,
      0.4610812906784e-07, 0.5858503336994e+01, 0.4337116142245e+00,
      0.3923022803444e-07, 0.3354170010125e+00, 0.1484170571900e-02,
      0.2610725582128e-07, 0.5410600646324e+01, 0.6327837846670e+00,
      0.2455279767721e-07, 0.6120216681403e+01, 0.1162474756779e+01,
      0.2375530706525e-07, 0.6055443426143e+01, 0.1052268489556e+01,

      0.1782967577553e-07, 0.3146108708004e+01, 0.8460828644453e+00,
      0.1581687095238e-07, 0.6255496089819e+00, 0.3340612434717e+01,
      0.1594657672461e-07, 0.3782604300261e+01, 0.1066495398892e+01,
      0.1563448615040e-07, 0.1997775733196e+01, 0.2022531624851e+00,
      0.1463624258525e-07, 0.1736316792088e+00, 0.3516457698740e-01,
      0.1331585056673e-07, 0.4331941830747e+01, 0.9491756770005e+00,
      0.1130634557637e-07, 0.6152017751825e+01, 0.2968341143800e-02,
      0.1028949607145e-07, 0.2101792614637e+00, 0.2275259891141e+00,
      0.1024074971618e-07, 0.4071833211074e+01, 0.5070101000000e-01,
      0.8826956060303e-08, 0.4861633688145e+00, 0.2093666171530e+00,

      0.8572230171541e-08, 0.5268190724302e+01, 0.4110125927500e-01,
      0.7649332643544e-08, 0.5134543417106e+01, 0.2608790314060e+02,
      0.8581673291033e-08, 0.2920218146681e+01, 0.1480791608091e+00,
      0.8430589300938e-08, 0.3604576619108e+01, 0.2172315424036e+00,
      0.7776165501012e-08, 0.3772942249792e+01, 0.6373574839730e-01,
      0.8311070234408e-08, 0.6200412329888e+01, 0.3235053470014e+00,
      0.6927365212582e-08, 0.4543353113437e+01, 0.8531963191132e+00,
      0.6791574208598e-08, 0.2882188406238e+01, 0.7181332454670e-01,
      0.5593100811839e-08, 0.1776646892780e+01, 0.7429900518901e+00,
      0.4553381853021e-08, 0.3949617611240e+01, 0.7775000683430e-01,

      0.5758000450068e-08, 0.3859251775075e+01, 0.1990721704425e+00,
      0.4281283457133e-08, 0.1466294631206e+01, 0.2118763888447e+01,
      0.4206935661097e-08, 0.5421776011706e+01, 0.1104591729320e-01,
      0.4213751641837e-08, 0.3412048993322e+01, 0.2243449970715e+00,
      0.5310506239878e-08, 0.5421641370995e+00, 0.5154640627760e+00,
      0.3827450341320e-08, 0.8887314524995e+00, 0.1510475019529e+00,
      0.4292435241187e-08, 0.1405043757194e+01, 0.1422690933580e-01,
      0.3189780702289e-08, 0.1060049293445e+01, 0.1173197218910e+00,
      0.3226611928069e-08, 0.6270858897442e+01, 0.2164800718209e+00,
      0.2893897608830e-08, 0.5117563223301e+01, 0.6470106940028e+00,

      0.3239852024578e-08, 0.4079092237983e+01, 0.2101180877357e+00,
      0.2956892222200e-08, 0.1594917021704e+01, 0.3092784376656e+00,
      0.2980177912437e-08, 0.5258787667564e+01, 0.4155522422634e+00,
      0.3163725690776e-08, 0.3854589225479e+01, 0.8582758298370e-01,
      0.2662262399118e-08, 0.3561326430187e+01, 0.5257585094865e+00,
      0.2766689135729e-08, 0.3180732086830e+00, 0.1385174140878e+00,
      0.2411600278464e-08, 0.3324798335058e+01, 0.5439178814476e+00,
      0.2483527695131e-08, 0.4169069291947e+00, 0.5336234347371e+00,
      0.7788777276590e-09, 0.1900569908215e+01, 0.5217580628120e+02 };

/* SSB-to-Sun, T^1, X */
   static const double s1x[] = {
     -0.1296310361520e-07, 0.0000000000000e+00, 0.0000000000000e+00,
      0.8975769009438e-08, 0.1128891609250e+01, 0.4265981595566e+00,
      0.7771113441307e-08, 0.2706039877077e+01, 0.2061856251104e+00,
      0.7538303866642e-08, 0.2191281289498e+01, 0.2204125344462e+00,
      0.6061384579336e-08, 0.3248167319958e+01, 0.1059381944224e+01,
      0.5726994235594e-08, 0.5569981398610e+01, 0.5225775174439e+00,
      0.5616492836424e-08, 0.5057386614909e+01, 0.5368044267797e+00,
      0.1010881584769e-08, 0.3473577116095e+01, 0.7113454667900e-02,
      0.7259606157626e-09, 0.3651858593665e+00, 0.6398972393349e+00,
      0.8755095026935e-09, 0.1662835408338e+01, 0.4194847048887e+00,

      0.5370491182812e-09, 0.1327673878077e+01, 0.4337116142245e+00,
      0.5743773887665e-09, 0.4250200846687e+01, 0.2132990797783e+00,
      0.4408103140300e-09, 0.3598752574277e+01, 0.1589072916335e+01,
      0.3101892374445e-09, 0.4887822983319e+01, 0.1052268489556e+01,
      0.3209453713578e-09, 0.9702272295114e+00, 0.5296909721118e+00,
      0.3017228286064e-09, 0.5484462275949e+01, 0.1066495398892e+01,
      0.3200700038601e-09, 0.2846613338643e+01, 0.1495633313810e+00,
      0.2137637279911e-09, 0.5692163292729e+00, 0.3163918923335e+00,
      0.1899686386727e-09, 0.2061077157189e+01, 0.2275259891141e+00,
      0.1401994545308e-09, 0.4177771136967e+01, 0.1102062672231e+00,

      0.1578057810499e-09, 0.5782460597335e+01, 0.7626583626240e-01,
      0.1237713253351e-09, 0.5705900866881e+01, 0.5154640627760e+00,
      0.1313076837395e-09, 0.5163438179576e+01, 0.3664874755930e-01,
      0.1184963304860e-09, 0.3054804427242e+01, 0.6327837846670e+00,
      0.1238130878565e-09, 0.2317292575962e+01, 0.3961708870310e-01,
      0.1015959527736e-09, 0.2194643645526e+01, 0.7329749511860e-01,
      0.9017954423714e-10, 0.2868603545435e+01, 0.1990721704425e+00,
      0.8668024955603e-10, 0.4923849675082e+01, 0.5439178814476e+00,
      0.7756083930103e-10, 0.3014334135200e+01, 0.9491756770005e+00,
      0.7536503401741e-10, 0.2704886279769e+01, 0.1030928125552e+00,

      0.5483308679332e-10, 0.6010983673799e+01, 0.8531963191132e+00,
      0.5184339620428e-10, 0.1952704573291e+01, 0.2093666171530e+00,
      0.5108658712030e-10, 0.2958575786649e+01, 0.2172315424036e+00,
      0.5019424524650e-10, 0.1736317621318e+01, 0.2164800718209e+00,
      0.4909312625978e-10, 0.3167216416257e+01, 0.2101180877357e+00,
      0.4456638901107e-10, 0.7697579923471e+00, 0.3235053470014e+00,
      0.4227030350925e-10, 0.3490910137928e+01, 0.6373574839730e-01,
      0.4095456040093e-10, 0.5178888984491e+00, 0.6470106940028e+00,
      0.4990537041422e-10, 0.3323887668974e+01, 0.1422690933580e-01,
      0.4321170010845e-10, 0.4288484987118e+01, 0.7358765972222e+00,

      0.3544072091802e-10, 0.6021051579251e+01, 0.5265099800692e+00,
      0.3480198638687e-10, 0.4600027054714e+01, 0.5328719641544e+00,
      0.3440287244435e-10, 0.4349525970742e+01, 0.8582758298370e-01,
      0.3330628322713e-10, 0.2347391505082e+01, 0.1104591729320e-01,
      0.2973060707184e-10, 0.4789409286400e+01, 0.5257585094865e+00,
      0.2932606766089e-10, 0.5831693799927e+01, 0.5336234347371e+00,
      0.2876972310953e-10, 0.2692638514771e+01, 0.1173197218910e+00,
      0.2827488278556e-10, 0.2056052487960e+01, 0.2022531624851e+00,
      0.2515028239756e-10, 0.7411863262449e+00, 0.9597935788730e-01,
      0.2853033744415e-10, 0.3948481024894e+01, 0.2118763888447e+01 };

/* SSB-to-Sun, T^1, Y */
   static const double s1y[] = {
      0.8989047573576e-08, 0.5840593672122e+01, 0.4265981595566e+00,
      0.7815938401048e-08, 0.1129664707133e+01, 0.2061856251104e+00,
      0.7550926713280e-08, 0.6196589104845e+00, 0.2204125344462e+00,
      0.6056556925895e-08, 0.1677494667846e+01, 0.1059381944224e+01,
      0.5734142698204e-08, 0.4000920852962e+01, 0.5225775174439e+00,
      0.5614341822459e-08, 0.3486722577328e+01, 0.5368044267797e+00,
      0.1028678147656e-08, 0.1877141024787e+01, 0.7113454667900e-02,
      0.7270792075266e-09, 0.5077167301739e+01, 0.6398972393349e+00,
      0.8734141726040e-09, 0.9069550282609e-01, 0.4194847048887e+00,
      0.5377371402113e-09, 0.6039381844671e+01, 0.4337116142245e+00,

      0.4729719431571e-09, 0.2153086311760e+01, 0.2132990797783e+00,
      0.4458052820973e-09, 0.5059830025565e+01, 0.5296909721118e+00,
      0.4406855467908e-09, 0.2027971692630e+01, 0.1589072916335e+01,
      0.3101659310977e-09, 0.3317677981860e+01, 0.1052268489556e+01,
      0.3016749232545e-09, 0.3913703482532e+01, 0.1066495398892e+01,
      0.3198541352656e-09, 0.1275513098525e+01, 0.1495633313810e+00,
      0.2142065389871e-09, 0.5301351614597e+01, 0.3163918923335e+00,
      0.1902615247592e-09, 0.4894943352736e+00, 0.2275259891141e+00,
      0.1613410990871e-09, 0.2449891130437e+01, 0.1102062672231e+00,
      0.1576992165097e-09, 0.4211421447633e+01, 0.7626583626240e-01,

      0.1241637259894e-09, 0.4140803368133e+01, 0.5154640627760e+00,
      0.1313974830355e-09, 0.3591920305503e+01, 0.3664874755930e-01,
      0.1181697118258e-09, 0.1506314382788e+01, 0.6327837846670e+00,
      0.1238239742779e-09, 0.7461405378404e+00, 0.3961708870310e-01,
      0.1010107068241e-09, 0.6271010795475e+00, 0.7329749511860e-01,
      0.9226316616509e-10, 0.1259158839583e+01, 0.1990721704425e+00,
      0.8664946419555e-10, 0.3353244696934e+01, 0.5439178814476e+00,
      0.7757230468978e-10, 0.1447677295196e+01, 0.9491756770005e+00,
      0.7693168628139e-10, 0.1120509896721e+01, 0.1030928125552e+00,
      0.5487897454612e-10, 0.4439380426795e+01, 0.8531963191132e+00,

      0.5196118677218e-10, 0.3788856619137e+00, 0.2093666171530e+00,
      0.5110853339935e-10, 0.1386879372016e+01, 0.2172315424036e+00,
      0.5027804534813e-10, 0.1647881805466e+00, 0.2164800718209e+00,
      0.4922485922674e-10, 0.1594315079862e+01, 0.2101180877357e+00,
      0.6155599524400e-10, 0.0000000000000e+00, 0.0000000000000e+00,
      0.4447147832161e-10, 0.5480720918976e+01, 0.3235053470014e+00,
      0.4144691276422e-10, 0.1931371033660e+01, 0.6373574839730e-01,
      0.4099950625452e-10, 0.5229611294335e+01, 0.6470106940028e+00,
      0.5060541682953e-10, 0.1731112486298e+01, 0.1422690933580e-01,
      0.4293615946300e-10, 0.2714571038925e+01, 0.7358765972222e+00,

      0.3545659845763e-10, 0.4451041444634e+01, 0.5265099800692e+00,
      0.3479112041196e-10, 0.3029385448081e+01, 0.5328719641544e+00,
      0.3438516493570e-10, 0.2778507143731e+01, 0.8582758298370e-01,
      0.3297341285033e-10, 0.7898709807584e+00, 0.1104591729320e-01,
      0.2972585818015e-10, 0.3218785316973e+01, 0.5257585094865e+00,
      0.2931707295017e-10, 0.4260731012098e+01, 0.5336234347371e+00,
      0.2897198149403e-10, 0.1120753978101e+01, 0.1173197218910e+00,
      0.2832293240878e-10, 0.4597682717827e+00, 0.2022531624851e+00,
      0.2864348326612e-10, 0.2169939928448e+01, 0.9597935788730e-01,
      0.2852714675471e-10, 0.2377659870578e+01, 0.2118763888447e+01 };

/* SSB-to-Sun, T^1, Z */
   static const double s1z[] = {
      0.5444220475678e-08, 0.1803825509310e+01, 0.2132990797783e+00,
      0.3883412695596e-08, 0.4668616389392e+01, 0.5296909721118e+00,
      0.1334341434551e-08, 0.0000000000000e+00, 0.0000000000000e+00,
      0.3730001266883e-09, 0.5401405918943e+01, 0.2061856251104e+00,
      0.2894929197956e-09, 0.4932415609852e+01, 0.2204125344462e+00,
      0.2857950357701e-09, 0.3154625362131e+01, 0.7478166569050e-01,
      0.2499226432292e-09, 0.3657486128988e+01, 0.4265981595566e+00,
      0.1937705443593e-09, 0.5740434679002e+01, 0.1059381944224e+01,
      0.1374894396320e-09, 0.1712857366891e+01, 0.5368044267797e+00,
      0.1217248678408e-09, 0.2312090870932e+01, 0.5225775174439e+00,

      0.7961052740870e-10, 0.5283368554163e+01, 0.3813291813120e-01,
      0.4979225949689e-10, 0.4298290471860e+01, 0.4194847048887e+00,
      0.4388552286597e-10, 0.6145515047406e+01, 0.7113454667900e-02,
      0.2586835212560e-10, 0.3019448001809e+01, 0.6398972393349e+00 };

/* SSB-to-Sun, T^2, X */
   static const double s2x[] = {
      0.1603551636587e-11, 0.4404109410481e+01, 0.2061856251104e+00,
      0.1556935889384e-11, 0.4818040873603e+00, 0.2204125344462e+00,
      0.1182594414915e-11, 0.9935762734472e+00, 0.5225775174439e+00,
      0.1158794583180e-11, 0.3353180966450e+01, 0.5368044267797e+00,
      0.9597358943932e-12, 0.5567045358298e+01, 0.2132990797783e+00,
      0.6511516579605e-12, 0.5630872420788e+01, 0.4265981595566e+00,
      0.7419792747688e-12, 0.2156188581957e+01, 0.5296909721118e+00,
      0.3951972655848e-12, 0.1981022541805e+01, 0.1059381944224e+01,
      0.4478223877045e-12, 0.0000000000000e+00, 0.0000000000000e+00 };

/* SSB-to-Sun, T^2, Y */
   static const double s2y[] = {
      0.1609114495091e-11, 0.2831096993481e+01, 0.2061856251104e+00,
      0.1560330784946e-11, 0.5193058213906e+01, 0.2204125344462e+00,
      0.1183535479202e-11, 0.5707003443890e+01, 0.5225775174439e+00,
      0.1158183066182e-11, 0.1782400404928e+01, 0.5368044267797e+00,
      0.1032868027407e-11, 0.4036925452011e+01, 0.2132990797783e+00,
      0.6540142847741e-12, 0.4058241056717e+01, 0.4265981595566e+00,
      0.7305236491596e-12, 0.6175401942957e+00, 0.5296909721118e+00,
     -0.5580725052968e-12, 0.0000000000000e+00, 0.0000000000000e+00,
      0.3946122651015e-12, 0.4108265279171e+00, 0.1059381944224e+01 };

/* SSB-to-Sun, T^2, Z */
   static const double s2z[] = {
      0.3749920358054e-12, 0.3230285558668e+01, 0.2132990797783e+00,
      0.2735037220939e-12, 0.6154322683046e+01, 0.5296909721118e+00 };

/* Pointers to coefficient arrays, in x,y,z sets */
   static const double *ce0[] = { e0x, e0y, e0z },
                       *ce1[] = { e1x, e1y, e1z },
                       *ce2[] = { e2x, e2y, e2z },
                       *cs0[] = { s0x, s0y, s0z },
                       *cs1[] = { s1x, s1y, s1z },
                       *cs2[] = { s2x, s2y, s2z };
   const double *coeffs;

/* Numbers of terms for each component of the model, in x,y,z sets */
   static const int ne0[3] = {(int)(sizeof e0x / sizeof (double) / 3),
                              (int)(sizeof e0y / sizeof (double) / 3),
                              (int)(sizeof e0z / sizeof (double) / 3) },
                    ne1[3] = {(int)(sizeof e1x / sizeof (double) / 3),
                              (int)(sizeof e1y / sizeof (double) / 3),
                              (int)(sizeof e1z / sizeof (double) / 3) },
                    ne2[3] = {(int)(sizeof e2x / sizeof (double) / 3),
                              (int)(sizeof e2y / sizeof (double) / 3),
                              (int)(sizeof e2z / sizeof (double) / 3) },
                    ns0[3] = {(int)(sizeof s0x / sizeof (double) / 3),
                              (int)(sizeof s0y / sizeof (double) / 3),
                              (int)(sizeof s0z / sizeof (double) / 3) },
                    ns1[3] = {(int)(sizeof s1x / sizeof (double) / 3),
                              (int)(sizeof s1y / sizeof (double) / 3),
                              (int)(sizeof s1z / sizeof (double) / 3) },
                    ns2[3] = {(int)(sizeof s2x / sizeof (double) / 3),
                              (int)(sizeof s2y / sizeof (double) / 3),
                              (int)(sizeof s2z / sizeof (double) / 3) };
   int nterms;

/* Miscellaneous */
   int jstat, i, j;
   double t, t2, xyz, xyzd, a, b, c, ct, p, cp,
          ph[3], vh[3], pb[3], vb[3], x, y, z;

/*--------------------------------------------------------------------*/

/* Time since reference epoch, Julian years. */
   t = ((date1 - DJ00) + date2) / DJY;
   t2 = t*t;

/* Set status. */
   jstat = fabs(t) <= 100.0 ? 0 : 1;

/* X then Y then Z. */
   for (i = 0; i < 3; i++) {

   /* Initialize position and velocity component. */
      xyz = 0.0;
      xyzd = 0.0;

   /* ------------------------------------------------ */
   /* Obtain component of Sun to Earth ecliptic vector */
   /* ------------------------------------------------ */

   /* Sun to Earth, T^0 terms. */
      coeffs = ce0[i];
      nterms = ne0[i];
      for (j = 0; j < nterms; j++) {
         a = *coeffs++;
         b = *coeffs++;
         c = *coeffs++;
         p = b + c*t;
         xyz  += a*cos(p);
         xyzd -= a*c*sin(p);
      }

   /* Sun to Earth, T^1 terms. */
      coeffs = ce1[i];
      nterms = ne1[i];
      for (j = 0; j < nterms; j++) {
         a = *coeffs++;
         b = *coeffs++;
         c = *coeffs++;
         ct = c*t;
         p = b + ct;
         cp = cos(p);
         xyz  += a*t*cp;
         xyzd += a*( cp - ct*sin(p) );
      }

   /* Sun to Earth, T^2 terms. */
      coeffs = ce2[i];
      nterms = ne2[i];
      for (j = 0; j < nterms; j++) {
         a = *coeffs++;
         b = *coeffs++;
         c = *coeffs++;
         ct = c*t;
         p = b + ct;
         cp = cos(p);
         xyz  += a*t2*cp;
         xyzd += a*t*( 2.0*cp - ct*sin(p) );
      }

   /* Heliocentric Earth position and velocity component. */
      ph[i] = xyz;
      vh[i] = xyzd / DJY;

   /* ------------------------------------------------ */
   /* Obtain component of SSB to Earth ecliptic vector */
   /* ------------------------------------------------ */

   /* SSB to Sun, T^0 terms. */
      coeffs = cs0[i];
      nterms = ns0[i];
      for (j = 0; j < nterms; j++) {
         a = *coeffs++;
         b = *coeffs++;
         c = *coeffs++;
         p = b + c*t;
         xyz  += a*cos(p);
         xyzd -= a*c*sin(p);
      }

   /* SSB to Sun, T^1 terms. */
      coeffs = cs1[i];
      nterms = ns1[i];
      for (j = 0; j < nterms; j++) {
         a = *coeffs++;
         b = *coeffs++;
         c = *coeffs++;
         ct = c*t;
         p = b + ct;
         cp = cos(p);
         xyz  += a*t*cp;
         xyzd += a*(cp - ct*sin(p));
      }

   /* SSB to Sun, T^2 terms. */
      coeffs = cs2[i];
      nterms = ns2[i];
      for (j = 0; j < nterms; j++) {
         a = *coeffs++;
         b = *coeffs++;
         c = *coeffs++;
         ct = c*t;
         p = b + ct;
         cp = cos(p);
         xyz  += a*t2*cp;
         xyzd += a*t*(2.0*cp - ct*sin(p));
     }

   /* Barycentric Earth position and velocity component. */
     pb[i] = xyz;
     vb[i] = xyzd / DJY;

   /* Next Cartesian component. */
   }

/* Rotate from ecliptic to BCRS coordinates. */

   x = ph[0];
   y = ph[1];
   z = ph[2];
   pvh[0][0] =      x + am12*y + am13*z;
   pvh[0][1] = am21*x + am22*y + am23*z;
   pvh[0][2] =          am32*y + am33*z;

   x = vh[0];
   y = vh[1];
   z = vh[2];
   pvh[1][0] =      x + am12*y + am13*z;
   pvh[1][1] = am21*x + am22*y + am23*z;
   pvh[1][2] =          am32*y + am33*z;

   x = pb[0];
   y = pb[1];
   z = pb[2];
   pvb[0][0] =      x + am12*y + am13*z;
   pvb[0][1] = am21*x + am22*y + am23*z;
   pvb[0][2] =          am32*y + am33*z;

   x = vb[0];
   y = vb[1];
   z = vb[2];
   pvb[1][0] =      x + am12*y + am13*z;
   pvb[1][1] = am21*x + am22*y + am23*z;
   pvb[1][2] =          am32*y + am33*z;

/* Return the status. */
   return jstat;

}

double eraEqeq94(double date1, double date2)
/*
**  - - - - - - - - - -
**   e r a E q e q 9 4
**  - - - - - - - - - -
**
**  Equation of the equinoxes, IAU 1994 model.
**
**  Given:
**     date1,date2   double     TDB date (Note 1)
**
**  Returned (function value):
**                   double     equation of the equinoxes (Note 2)
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
**  2) The result, which is in radians, operates in the following sense:
**
**        Greenwich apparent ST = GMST + equation of the equinoxes
**
**  Called:
**     eraNut80     nutation, IAU 1980
**     eraObl80     mean obliquity, IAU 1980
**
**  References:
**
**     IAU Resolution C7, Recommendation 3 (1994).
**
**     Capitaine, N. & Gontier, A.-M., 1993, Astron. Astrophys., 275,
**     645-650.
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double t,  om,  dpsi,  deps,  eps0, ee;


/* Interval between fundamental epoch J2000.0 and given date (JC). */
   t = ((date1 - DJ00) + date2) / DJC;

/* Longitude of the mean ascending node of the lunar orbit on the */
/* ecliptic, measured from the mean equinox of date. */
   om = eraAnpm((450160.280 + (-482890.539
           + (7.455 + 0.008 * t) * t) * t) * DAS2R
           + fmod(-5.0 * t, 1.0) * D2PI);

/* Nutation components and mean obliquity. */
   eraNut80(date1, date2, &dpsi, &deps);
   eps0 = eraObl80(date1, date2);

/* Equation of the equinoxes. */
   ee = dpsi*cos(eps0) + DAS2R*(0.00264*sin(om) + 0.000063*sin(om+om));

   return ee;

}

double eraEra00(double dj1, double dj2)
/*
**  - - - - - - - - -
**   e r a E r a 0 0
**  - - - - - - - - -
**
**  Earth rotation angle (IAU 2000 model).
**
**  Given:
**     dj1,dj2   double    UT1 as a 2-part Julian Date (see note)
**
**  Returned (function value):
**               double    Earth rotation angle (radians), range 0-2pi
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
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in
**     cases where the loss of several decimal digits of resolution
**     is acceptable.  The J2000 and MJD methods are good compromises
**     between resolution and convenience.  The date & time method is
**     best matched to the algorithm used:  maximum precision is
**     delivered when the dj1 argument is for 0hrs UT1 on the day in
**     question and the dj2 argument lies in the range 0 to 1, or vice
**     versa.
**
**  2) The algorithm is adapted from Expression 22 of Capitaine et al.
**     2000.  The time argument has been expressed in days directly,
**     and, to retain precision, integer contributions have been
**     eliminated.  The same formulation is given in IERS Conventions
**     (2003), Chap. 5, Eq. 14.
**
**  Called:
**     eraAnp       normalize angle into range 0 to 2pi
**
**  References:
**
**     Capitaine N., Guinot B. and McCarthy D.D, 2000, Astron.
**     Astrophys., 355, 398-405.
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double d1, d2, t, f, theta;


/* Days since fundamental epoch. */
   if (dj1 < dj2) {
      d1 = dj1;
      d2 = dj2;
   } else {
      d1 = dj2;
      d2 = dj1;
   }
   t = d1 + (d2- DJ00);

/* Fractional part of T (days). */
   f = fmod(d1, 1.0) + fmod(d2, 1.0);

/* Earth rotation angle at this UT1. */
   theta = eraAnp(D2PI * (f + 0.7790572732640
                            + 0.00273781191135448 * t));

   return theta;

}

double eraFad03(double t)
/*
**  - - - - - - - - -
**   e r a F a d 0 3
**  - - - - - - - - -
**
**  Fundamental argument, IERS Conventions (2003):
**  mean elongation of the Moon from the Sun.
**
**  Given:
**     t     double    TDB, Julian centuries since J2000.0 (Note 1)
**
**  Returned (function value):
**           double    D, radians (Note 2)
**
**  Notes:
**
**  1) Though t is strictly TDB, it is usually more convenient to use
**     TT, which makes no significant difference.
**
**  2) The expression used is as adopted in IERS Conventions (2003) and
**     is from Simon et al. (1994).
**
**  References:
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
**     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double a;


/* Mean elongation of the Moon from the Sun (IERS Conventions 2003). */
   a = fmod(          1072260.703692 +
             t * ( 1602961601.2090 +
             t * (        - 6.3706 +
             t * (          0.006593 +
             t * (        - 0.00003169 ) ) ) ), TURNAS ) * DAS2R;

   return a;

}

double eraFae03(double t)
/*
**  - - - - - - - - -
**   e r a F a e 0 3
**  - - - - - - - - -
**
**  Fundamental argument, IERS Conventions (2003):
**  mean longitude of Earth.
**
**  Given:
**     t     double    TDB, Julian centuries since J2000.0 (Note 1)
**
**  Returned (function value):
**           double    mean longitude of Earth, radians (Note 2)
**
**  Notes:
**
**  1) Though t is strictly TDB, it is usually more convenient to use
**     TT, which makes no significant difference.
**
**  2) The expression used is as adopted in IERS Conventions (2003) and
**     comes from Souchay et al. (1999) after Simon et al. (1994).
**
**  References:
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
**     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
**
**     Souchay, J., Loysel, B., Kinoshita, H., Folgueira, M. 1999,
**     Astron.Astrophys.Supp.Ser. 135, 111
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double a;


/* Mean longitude of Earth (IERS Conventions 2003). */
   a = fmod(1.753470314 + 628.3075849991 * t, D2PI);

   return a;

}

double eraFaf03(double t)
/*
**  - - - - - - - - -
**   e r a F a f 0 3
**  - - - - - - - - -
**
**  Fundamental argument, IERS Conventions (2003):
**  mean longitude of the Moon minus mean longitude of the ascending
**  node.
**
**  Given:
**     t     double    TDB, Julian centuries since J2000.0 (Note 1)
**
**  Returned (function value):
**           double    F, radians (Note 2)
**
**  Notes:
**
**  1) Though t is strictly TDB, it is usually more convenient to use
**     TT, which makes no significant difference.
**
**  2) The expression used is as adopted in IERS Conventions (2003) and
**     is from Simon et al. (1994).
**
**  References:
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
**     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double a;


/* Mean longitude of the Moon minus that of the ascending node */
/* (IERS Conventions 2003).                                    */
   a = fmod(           335779.526232 +
             t * ( 1739527262.8478 +
             t * (       - 12.7512 +
             t * (        - 0.001037 +
             t * (          0.00000417 ) ) ) ), TURNAS ) * DAS2R;

   return a;


}

double eraFaju03(double t)
/*
**  - - - - - - - - - -
**   e r a F a j u 0 3
**  - - - - - - - - - -
**
**  Fundamental argument, IERS Conventions (2003):
**  mean longitude of Jupiter.
**
**  Given:
**     t     double    TDB, Julian centuries since J2000.0 (Note 1)
**
**  Returned (function value):
**           double    mean longitude of Jupiter, radians (Note 2)
**
**  Notes:
**
**  1) Though t is strictly TDB, it is usually more convenient to use
**     TT, which makes no significant difference.
**
**  2) The expression used is as adopted in IERS Conventions (2003) and
**     comes from Souchay et al. (1999) after Simon et al. (1994).
**
**  References:
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
**     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
**
**     Souchay, J., Loysel, B., Kinoshita, H., Folgueira, M. 1999,
**     Astron.Astrophys.Supp.Ser. 135, 111
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double a;


/* Mean longitude of Jupiter (IERS Conventions 2003). */
   a = fmod(0.599546497 + 52.9690962641 * t, D2PI);

   return a;

}

double eraFal03(double t)
/*
**  - - - - - - - - -
**   e r a F a l 0 3
**  - - - - - - - - -
**
**  Fundamental argument, IERS Conventions (2003):
**  mean anomaly of the Moon.
**
**  Given:
**     t     double    TDB, Julian centuries since J2000.0 (Note 1)
**
**  Returned (function value):
**           double    l, radians (Note 2)
**
**  Notes:
**
**  1) Though t is strictly TDB, it is usually more convenient to use
**     TT, which makes no significant difference.
**
**  2) The expression used is as adopted in IERS Conventions (2003) and
**     is from Simon et al. (1994).
**
**  References:
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
**     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double a;


/* Mean anomaly of the Moon (IERS Conventions 2003). */
   a = fmod(           485868.249036  +
             t * ( 1717915923.2178 +
             t * (         31.8792 +
             t * (          0.051635 +
             t * (        - 0.00024470 ) ) ) ), TURNAS ) * DAS2R;

   return a;

}

double eraFalp03(double t)
/*
**  - - - - - - - - - -
**   e r a F a l p 0 3
**  - - - - - - - - - -
**
**  Fundamental argument, IERS Conventions (2003):
**  mean anomaly of the Sun.
**
**  Given:
**     t     double    TDB, Julian centuries since J2000.0 (Note 1)
**
**  Returned (function value):
**           double    l', radians (Note 2)
**
**  Notes:
**
**  1) Though t is strictly TDB, it is usually more convenient to use
**     TT, which makes no significant difference.
**
**  2) The expression used is as adopted in IERS Conventions (2003) and
**     is from Simon et al. (1994).
**
**  References:
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
**     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double a;


/* Mean anomaly of the Sun (IERS Conventions 2003). */
   a = fmod(         1287104.793048 +
             t * ( 129596581.0481 +
             t * (       - 0.5532 +
             t * (         0.000136 +
             t * (       - 0.00001149 ) ) ) ), TURNAS ) * DAS2R;

   return a;

}

double eraFama03(double t)
/*
**  - - - - - - - - - -
**   e r a F a m a 0 3
**  - - - - - - - - - -
**
**  Fundamental argument, IERS Conventions (2003):
**  mean longitude of Mars.
**
**  Given:
**     t     double    TDB, Julian centuries since J2000.0 (Note 1)
**
**  Returned (function value):
**           double    mean longitude of Mars, radians (Note 2)
**
**  Notes:
**
**  1) Though t is strictly TDB, it is usually more convenient to use
**     TT, which makes no significant difference.
**
**  2) The expression used is as adopted in IERS Conventions (2003) and
**     comes from Souchay et al. (1999) after Simon et al. (1994).
**
**  References:
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
**     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
**
**     Souchay, J., Loysel, B., Kinoshita, H., Folgueira, M. 1999,
**     Astron.Astrophys.Supp.Ser. 135, 111
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double a;


/* Mean longitude of Mars (IERS Conventions 2003). */
   a = fmod(6.203480913 + 334.0612426700 * t, D2PI);

   return a;

}

double eraFame03(double t)
/*
**  - - - - - - - - - -
**   e r a F a m e 0 3
**  - - - - - - - - - -
**
**  Fundamental argument, IERS Conventions (2003):
**  mean longitude of Mercury.
**
**  Given:
**     t     double    TDB, Julian centuries since J2000.0 (Note 1)
**
**  Returned (function value):
**           double    mean longitude of Mercury, radians (Note 2)
**
**  Notes:
**
**  1) Though t is strictly TDB, it is usually more convenient to use
**     TT, which makes no significant difference.
**
**  2) The expression used is as adopted in IERS Conventions (2003) and
**     comes from Souchay et al. (1999) after Simon et al. (1994).
**
**  References:
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
**     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
**
**     Souchay, J., Loysel, B., Kinoshita, H., Folgueira, M. 1999,
**     Astron.Astrophys.Supp.Ser. 135, 111
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double a;


/* Mean longitude of Mercury (IERS Conventions 2003). */
   a = fmod(4.402608842 + 2608.7903141574 * t, D2PI);

   return a;

}

double eraFane03(double t)
/*
**  - - - - - - - - - -
**   e r a F a n e 0 3
**  - - - - - - - - - -
**
**  Fundamental argument, IERS Conventions (2003):
**  mean longitude of Neptune.
**
**  Given:
**     t     double    TDB, Julian centuries since J2000.0 (Note 1)
**
**  Returned (function value):
**           double    mean longitude of Neptune, radians (Note 2)
**
**  Notes:
**
**  1) Though t is strictly TDB, it is usually more convenient to use
**     TT, which makes no significant difference.
**
**  2) The expression used is as adopted in IERS Conventions (2003) and
**     is adapted from Simon et al. (1994).
**
**  References:
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
**     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double a;


/* Mean longitude of Neptune (IERS Conventions 2003). */
   a = fmod(5.311886287 + 3.8133035638 * t, D2PI);

   return a;

}

double eraFaom03(double t)
/*
**  - - - - - - - - - -
**   e r a F a o m 0 3
**  - - - - - - - - - -
**
**  Fundamental argument, IERS Conventions (2003):
**  mean longitude of the Moon's ascending node.
**
**  Given:
**     t     double    TDB, Julian centuries since J2000.0 (Note 1)
**
**  Returned (function value):
**           double    Omega, radians (Note 2)
**
**  Notes:
**
**  1) Though t is strictly TDB, it is usually more convenient to use
**     TT, which makes no significant difference.
**
**  2) The expression used is as adopted in IERS Conventions (2003) and
**     is from Simon et al. (1994).
**
**  References:
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
**     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double a;


/* Mean longitude of the Moon's ascending node */
/* (IERS Conventions 2003).                    */
   a = fmod(          450160.398036 +
             t * ( - 6962890.5431 +
             t * (         7.4722 +
             t * (         0.007702 +
             t * (       - 0.00005939 ) ) ) ), TURNAS ) * DAS2R;

   return a;

}

double eraFapa03(double t)
/*
**  - - - - - - - - - -
**   e r a F a p a 0 3
**  - - - - - - - - - -
**
**  Fundamental argument, IERS Conventions (2003):
**  general accumulated precession in longitude.
**
**  Given:
**     t     double    TDB, Julian centuries since J2000.0 (Note 1)
**
**  Returned (function value):
**           double    general precession in longitude, radians (Note 2)
**
**  Notes:
**
**  1) Though t is strictly TDB, it is usually more convenient to use
**     TT, which makes no significant difference.
**
**  2) The expression used is as adopted in IERS Conventions (2003).  It
**     is taken from Kinoshita & Souchay (1990) and comes originally
**     from Lieske et al. (1977).
**
**  References:
**
**     Kinoshita, H. and Souchay J. 1990, Celest.Mech. and Dyn.Astron.
**     48, 187
**
**     Lieske, J.H., Lederle, T., Fricke, W. & Morando, B. 1977,
**     Astron.Astrophys. 58, 1-16
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double a;


/* General accumulated precession in longitude. */
   a = (0.024381750 + 0.00000538691 * t) * t;

   return a;

}

double eraFasa03(double t)
/*
**  - - - - - - - - - -
**   e r a F a s a 0 3
**  - - - - - - - - - -
**
**  Fundamental argument, IERS Conventions (2003):
**  mean longitude of Saturn.
**
**  Given:
**     t     double    TDB, Julian centuries since J2000.0 (Note 1)
**
**  Returned (function value):
**           double    mean longitude of Saturn, radians (Note 2)
**
**  Notes:
**
**  1) Though t is strictly TDB, it is usually more convenient to use
**     TT, which makes no significant difference.
**
**  2) The expression used is as adopted in IERS Conventions (2003) and
**     comes from Souchay et al. (1999) after Simon et al. (1994).
**
**  References:
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
**     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
**
**     Souchay, J., Loysel, B., Kinoshita, H., Folgueira, M. 1999,
**     Astron.Astrophys.Supp.Ser. 135, 111
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double a;


/* Mean longitude of Saturn (IERS Conventions 2003). */
   a = fmod(0.874016757 + 21.3299104960 * t, D2PI);

   return a;

}

double eraFaur03(double t)
/*
**  - - - - - - - - - -
**   e r a F a u r 0 3
**  - - - - - - - - - -
**
**  Fundamental argument, IERS Conventions (2003):
**  mean longitude of Uranus.
**
**  Given:
**     t     double    TDB, Julian centuries since J2000.0 (Note 1)
**
**  Returned  (function value):
**           double    mean longitude of Uranus, radians (Note 2)
**
**  Notes:
**
**  1) Though t is strictly TDB, it is usually more convenient to use
**     TT, which makes no significant difference.
**
**  2) The expression used is as adopted in IERS Conventions (2003) and
**     is adapted from Simon et al. (1994).
**
**  References:
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
**     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double a;


/* Mean longitude of Uranus (IERS Conventions 2003). */
   a = fmod(5.481293872 + 7.4781598567 * t, D2PI);

   return a;

}

double eraFave03(double t)
/*
**  - - - - - - - - - -
**   e r a F a v e 0 3
**  - - - - - - - - - -
**
**  Fundamental argument, IERS Conventions (2003):
**  mean longitude of Venus.
**
**  Given:
**     t     double    TDB, Julian centuries since J2000.0 (Note 1)
**
**  Returned (function value):
**           double    mean longitude of Venus, radians (Note 2)
**
**  Notes:
**
**  1) Though t is strictly TDB, it is usually more convenient to use
**     TT, which makes no significant difference.
**
**  2) The expression used is as adopted in IERS Conventions (2003) and
**     comes from Souchay et al. (1999) after Simon et al. (1994).
**
**  References:
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
**     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
**
**     Souchay, J., Loysel, B., Kinoshita, H., Folgueira, M. 1999,
**     Astron.Astrophys.Supp.Ser. 135, 111
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double a;


/* Mean longitude of Venus (IERS Conventions 2003). */
   a = fmod(3.176146697 + 1021.3285546211 * t, D2PI);

   return a;

}

void eraFk52h(double r5, double d5,
              double dr5, double dd5, double px5, double rv5,
              double *rh, double *dh,
              double *drh, double *ddh, double *pxh, double *rvh)
/*
**  - - - - - - - - -
**   e r a F k 5 2 h
**  - - - - - - - - -
**
**  Transform FK5 (J2000.0) star data into the Hipparcos system.
**
**  Given (all FK5, equinox J2000.0, epoch J2000.0):
**     r5      double    RA (radians)
**     d5      double    Dec (radians)
**     dr5     double    proper motion in RA (dRA/dt, rad/Jyear)
**     dd5     double    proper motion in Dec (dDec/dt, rad/Jyear)
**     px5     double    parallax (arcsec)
**     rv5     double    radial velocity (km/s, positive = receding)
**
**  Returned (all Hipparcos, epoch J2000.0):
**     rh      double    RA (radians)
**     dh      double    Dec (radians)
**     drh     double    proper motion in RA (dRA/dt, rad/Jyear)
**     ddh     double    proper motion in Dec (dDec/dt, rad/Jyear)
**     pxh     double    parallax (arcsec)
**     rvh     double    radial velocity (km/s, positive = receding)
**
**  Notes:
**
**  1) This function transforms FK5 star positions and proper motions
**     into the system of the Hipparcos catalog.
**
**  2) The proper motions in RA are dRA/dt rather than
**     cos(Dec)*dRA/dt, and are per year rather than per century.
**
**  3) The FK5 to Hipparcos transformation is modeled as a pure
**     rotation and spin;  zonal errors in the FK5 catalog are not
**     taken into account.
**
**  4) See also eraH2fk5, eraFk5hz, eraHfk5z.
**
**  Called:
**     eraStarpv    star catalog data to space motion pv-vector
**     eraFk5hip    FK5 to Hipparcos rotation and spin
**     eraRxp       product of r-matrix and p-vector
**     eraPxp       vector product of two p-vectors
**     eraPpp       p-vector plus p-vector
**     eraPvstar    space motion pv-vector to star catalog data
**
**  Reference:
**
**     F.Mignard & M.Froeschle, Astron. Astrophys. 354, 732-739 (2000).
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   int i;
   double pv5[2][3], r5h[3][3], s5h[3], wxp[3], vv[3], pvh[2][3];


/* FK5 barycentric position/velocity pv-vector (normalized). */
   eraStarpv(r5, d5, dr5, dd5, px5, rv5, pv5);

/* FK5 to Hipparcos orientation matrix and spin vector. */
   eraFk5hip(r5h, s5h);

/* Make spin units per day instead of per year. */
   for ( i = 0; i < 3; s5h[i++] /= 365.25 );

/* Orient the FK5 position into the Hipparcos system. */
   eraRxp(r5h, pv5[0], pvh[0]);

/* Apply spin to the position giving an extra space motion component. */
   eraPxp(pv5[0], s5h, wxp);

/* Add this component to the FK5 space motion. */
   eraPpp(wxp, pv5[1], vv);

/* Orient the FK5 space motion into the Hipparcos system. */
   eraRxp(r5h, vv, pvh[1]);

/* Hipparcos pv-vector to spherical. */
   eraPvstar(pvh, rh, dh, drh, ddh, pxh, rvh);

   return;

}

void eraFk5hip(double r5h[3][3], double s5h[3])
/*
**  - - - - - - - - - -
**   e r a F k 5 h i p
**  - - - - - - - - - -
**
**  FK5 to Hipparcos rotation and spin.
**
**  Returned:
**     r5h   double[3][3]  r-matrix: FK5 rotation wrt Hipparcos (Note 2)
**     s5h   double[3]     r-vector: FK5 spin wrt Hipparcos (Note 3)
**
**  Notes:
**
**  1) This function models the FK5 to Hipparcos transformation as a
**     pure rotation and spin;  zonal errors in the FK5 catalogue are
**     not taken into account.
**
**  2) The r-matrix r5h operates in the sense:
**
**           P_Hipparcos = r5h x P_FK5
**
**     where P_FK5 is a p-vector in the FK5 frame, and P_Hipparcos is
**     the equivalent Hipparcos p-vector.
**
**  3) The r-vector s5h represents the time derivative of the FK5 to
**     Hipparcos rotation.  The units are radians per year (Julian,
**     TDB).
**
**  Called:
**     eraRv2m      r-vector to r-matrix
**
**  Reference:
**
**     F.Mignard & M.Froeschle, Astron. Astrophys. 354, 732-739 (2000).
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double v[3];

/* FK5 wrt Hipparcos orientation and spin (radians, radians/year) */
   double epx, epy, epz;
   double omx, omy, omz;


   epx = -19.9e-3 * DAS2R;
   epy =  -9.1e-3 * DAS2R;
   epz =  22.9e-3 * DAS2R;

   omx = -0.30e-3 * DAS2R;
   omy =  0.60e-3 * DAS2R;
   omz =  0.70e-3 * DAS2R;

/* FK5 to Hipparcos orientation expressed as an r-vector. */
   v[0] = epx;
   v[1] = epy;
   v[2] = epz;

/* Re-express as an r-matrix. */
   eraRv2m(v, r5h);

/* Hipparcos wrt FK5 spin expressed as an r-vector. */
   s5h[0] = omx;
   s5h[1] = omy;
   s5h[2] = omz;

   return;

}

void eraFk5hz(double r5, double d5, double date1, double date2,
              double *rh, double *dh)
/*
**  - - - - - - - - -
**   e r a F k 5 h z
**  - - - - - - - - -
**
**  Transform an FK5 (J2000.0) star position into the system of the
**  Hipparcos catalogue, assuming zero Hipparcos proper motion.
**
**  Given:
**     r5           double   FK5 RA (radians), equinox J2000.0, at date
**     d5           double   FK5 Dec (radians), equinox J2000.0, at date
**     date1,date2  double   TDB date (Notes 1,2)
**
**  Returned:
**     rh           double   Hipparcos RA (radians)
**     dh           double   Hipparcos Dec (radians)
**
**  Notes:
**
**  1) This function converts a star position from the FK5 system to
**     the Hipparcos system, in such a way that the Hipparcos proper
**     motion is zero.  Because such a star has, in general, a non-zero
**     proper motion in the FK5 system, the function requires the date
**     at which the position in the FK5 system was determined.
**
**  2) The TT date date1+date2 is a Julian Date, apportioned in any
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
**  3) The FK5 to Hipparcos transformation is modeled as a pure
**     rotation and spin;  zonal errors in the FK5 catalogue are not
**     taken into account.
**
**  4) The position returned by this function is in the Hipparcos
**     reference system but at date date1+date2.
**
**  5) See also eraFk52h, eraH2fk5, eraHfk5z.
**
**  Called:
**     eraS2c       spherical coordinates to unit vector
**     eraFk5hip    FK5 to Hipparcos rotation and spin
**     eraSxp       multiply p-vector by scalar
**     eraRv2m      r-vector to r-matrix
**     eraTrxp      product of transpose of r-matrix and p-vector
**     eraPxp       vector product of two p-vectors
**     eraC2s       p-vector to spherical
**     eraAnp       normalize angle into range 0 to 2pi
**
**  Reference:
**
**     F.Mignard & M.Froeschle, 2000, Astron.Astrophys. 354, 732-739.
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double t, p5e[3], r5h[3][3], s5h[3], vst[3], rst[3][3], p5[3],
          ph[3], w;


/* Interval from given date to fundamental epoch J2000.0 (JY). */
   t = - ((date1 - DJ00) + date2) / DJY;

/* FK5 barycentric position vector. */
   eraS2c(r5, d5, p5e);

/* FK5 to Hipparcos orientation matrix and spin vector. */
   eraFk5hip(r5h, s5h);

/* Accumulated Hipparcos wrt FK5 spin over that interval. */
   eraSxp(t, s5h, vst);

/* Express the accumulated spin as a rotation matrix. */
   eraRv2m(vst, rst);

/* Derotate the vector's FK5 axes back to date. */
   eraTrxp(rst, p5e, p5);

/* Rotate the vector into the Hipparcos system. */
   eraRxp(r5h, p5, ph);

/* Hipparcos vector to spherical. */
   eraC2s(ph, &w, dh);
   *rh = eraAnp(w);

   return;

}

void eraFw2m(double gamb, double phib, double psi, double eps,
             double r[3][3])
/*
**  - - - - - - - -
**   e r a F w 2 m
**  - - - - - - - -
**
**  Form rotation matrix given the Fukushima-Williams angles.
**
**  Given:
**     gamb     double         F-W angle gamma_bar (radians)
**     phib     double         F-W angle phi_bar (radians)
**     psi      double         F-W angle psi (radians)
**     eps      double         F-W angle epsilon (radians)
**
**  Returned:
**     r        double[3][3]   rotation matrix
**
**  Notes:
**
**  1) Naming the following points:
**
**           e = J2000.0 ecliptic pole,
**           p = GCRS pole,
**           E = ecliptic pole of date,
**     and   P = CIP,
**
**     the four Fukushima-Williams angles are as follows:
**
**        gamb = gamma = epE
**        phib = phi = pE
**        psi = psi = pEP
**        eps = epsilon = EP
**
**  2) The matrix representing the combined effects of frame bias,
**     precession and nutation is:
**
**        NxPxB = R_1(-eps).R_3(-psi).R_1(phib).R_3(gamb)
**
**  3) Three different matrices can be constructed, depending on the
**     supplied angles:
**
**     o  To obtain the nutation x precession x frame bias matrix,
**        generate the four precession angles, generate the nutation
**        components and add them to the psi_bar and epsilon_A angles,
**        and call the present function.
**
**     o  To obtain the precession x frame bias matrix, generate the
**        four precession angles and call the present function.
**
**     o  To obtain the frame bias matrix, generate the four precession
**        angles for date J2000.0 and call the present function.
**
**     The nutation-only and precession-only matrices can if necessary
**     be obtained by combining these three appropriately.
**
**  Called:
**     eraIr        initialize r-matrix to identity
**     eraRz        rotate around Z-axis
**     eraRx        rotate around X-axis
**
**  Reference:
**
**     Hilton, J. et al., 2006, Celest.Mech.Dyn.Astron. 94, 351
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
/* Construct the matrix. */
   eraIr(r);
   eraRz(gamb, r);
   eraRx(phib, r);
   eraRz(-psi, r);
   eraRx(-eps, r);

   return;

}

void eraFw2xy(double gamb, double phib, double psi, double eps,
              double *x, double *y)
/*
**  - - - - - - - - -
**   e r a F w 2 x y
**  - - - - - - - - -
**
**  CIP X,Y given Fukushima-Williams bias-precession-nutation angles.
**
**  Given:
**     gamb     double    F-W angle gamma_bar (radians)
**     phib     double    F-W angle phi_bar (radians)
**     psi      double    F-W angle psi (radians)
**     eps      double    F-W angle epsilon (radians)
**
**  Returned:
**     x,y      double    CIP X,Y ("radians")
**
**  Notes:
**
**  1) Naming the following points:
**
**           e = J2000.0 ecliptic pole,
**           p = GCRS pole
**           E = ecliptic pole of date,
**     and   P = CIP,
**
**     the four Fukushima-Williams angles are as follows:
**
**        gamb = gamma = epE
**        phib = phi = pE
**        psi = psi = pEP
**        eps = epsilon = EP
**
**  2) The matrix representing the combined effects of frame bias,
**     precession and nutation is:
**
**        NxPxB = R_1(-epsA).R_3(-psi).R_1(phib).R_3(gamb)
**
**     X,Y are elements (3,1) and (3,2) of the matrix.
**
**  Called:
**     eraFw2m      F-W angles to r-matrix
**     eraBpn2xy    extract CIP X,Y coordinates from NPB matrix
**
**  Reference:
**
**     Hilton, J. et al., 2006, Celest.Mech.Dyn.Astron. 94, 351
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double r[3][3];


/* Form NxPxB matrix. */
   eraFw2m(gamb, phib, psi, eps, r);

/* Extract CIP X,Y. */
   eraBpn2xy(r, x, y);

   return;

}

int eraGc2gd ( int n, double xyz[3],
               double *elong, double *phi, double *height )
/*
**  - - - - - - - - -
**   e r a G c 2 g d
**  - - - - - - - - -
**
**  Transform geocentric coordinates to geodetic using the specified
**  reference ellipsoid.
**
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
**        1     WGS84
**        2     GRS80
**        3     WGS72
**
**     The n value has no significance outside the ERFA software.  For
**     convenience, symbols WGS84 etc. are defined in erfam.h.
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
*/
{
   int j;
   double a, f;


/* Obtain reference ellipsoid parameters. */
   j = eraEform ( n, &a, &f );

/* If OK, transform x,y,z to longitude, geodetic latitude, height. */
   if ( j == 0 ) {
      j = eraGc2gde ( a, f, xyz, elong, phi, height );
      if ( j < 0 ) j = -2;
   }

/* Deal with any errors. */
   if ( j < 0 ) {
      *phi = -1e9;
      *height = -1e9;
   }

/* Return the status. */
   return j;

}

int eraGc2gde ( double a, double f, double xyz[3],
                double *elong, double *phi, double *height )
/*
**  - - - - - - - - - -
**   e r a G c 2 g d e
**  - - - - - - - - - -
**
**  Transform geocentric coordinates to geodetic for a reference
**  ellipsoid of specified form.
**
**  Given:
**     a       double     equatorial radius (Notes 2,4)
**     f       double     flattening (Note 3)
**     xyz     double[3]  geocentric vector (Note 4)
**
**  Returned:
**     elong   double     longitude (radians, east +ve)
**     phi     double     latitude (geodetic, radians)
**     height  double     height above ellipsoid (geodetic, Note 4)
**
**  Returned (function value):
**             int        status:  0 = OK
**                                -1 = illegal f
**                                -2 = illegal a
**
**  Notes:
**
**  1) This function is based on the GCONV2H Fortran subroutine by
**     Toshio Fukushima (see reference).
**
**  2) The equatorial radius, a, can be in any units, but meters is
**     the conventional choice.
**
**  3) The flattening, f, is (for the Earth) a value around 0.00335,
**     i.e. around 1/298.
**
**  4) The equatorial radius, a, and the geocentric vector, xyz,
**     must be given in the same units, and determine the units of
**     the returned height, height.
**
**  5) If an error occurs (status < 0), elong, phi and height are
**     unchanged.
**
**  6) The inverse transformation is performed in the function
**     eraGd2gce.
**
**  7) The transformation for a standard ellipsoid (such as WGS84) can
**     more conveniently be performed by calling eraGc2gd, which uses a
**     numerical code to identify the required A and F values.
**
**  Reference:
**
**     Fukushima, T., "Transformation from Cartesian to geodetic
**     coordinates accelerated by Halley's method", J.Geodesy (2006)
**     79: 689-693
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double aeps2, e2, e4t, ec2, ec, b, x, y, z, p2, absz, p, s0, pn, zc,
                 c0, c02, c03, s02, s03, a02, a0, a03, d0, f0, b0, s1,
                 cc, s12, cc2;


/* ------------- */
/* Preliminaries */
/* ------------- */

/* Validate ellipsoid parameters. */
   if ( f < 0.0 || f >= 1.0 ) return -1;
   if ( a <= 0.0 ) return -2;

/* Functions of ellipsoid parameters (with further validation of f). */
   aeps2 = a*a * 1e-32;
   e2 = (2.0 - f) * f;
   e4t = e2*e2 * 1.5;
   ec2 = 1.0 - e2;
   if ( ec2 <= 0.0 ) return -1;
   ec = sqrt(ec2);
   b = a * ec;

/* Cartesian components. */
   x = xyz[0];
   y = xyz[1];
   z = xyz[2];

/* Distance from polar axis squared. */
   p2 = x*x + y*y;

/* Longitude. */
   *elong = p2 != 0.0 ? atan2(y, x) : 0.0;

/* Unsigned z-coordinate. */
   absz = fabs(z);

/* Proceed unless polar case. */
   if ( p2 > aeps2 ) {

   /* Distance from polar axis. */
      p = sqrt(p2);

   /* Normalization. */
      s0 = absz / a;
      pn = p / a;
      zc = ec * s0;

   /* Prepare Newton correction factors. */
      c0 = ec * pn;
      c02 = c0 * c0;
      c03 = c02 * c0;
      s02 = s0 * s0;
      s03 = s02 * s0;
      a02 = c02 + s02;
      a0 = sqrt(a02);
      a03 = a02 * a0;
      d0 = zc*a03 + e2*s03;
      f0 = pn*a03 - e2*c03;

   /* Prepare Halley correction factor. */
      b0 = e4t * s02 * c02 * pn * (a0 - ec);
      s1 = d0*f0 - b0*s0;
      cc = ec * (f0*f0 - b0*c0);

   /* Evaluate latitude and height. */
      *phi = atan(s1/cc);
      s12 = s1 * s1;
      cc2 = cc * cc;
      *height = (p*cc + absz*s1 - a * sqrt(ec2*s12 + cc2)) /
                                                        sqrt(s12 + cc2);
   } else {

   /* Exception: pole. */
      *phi = DPI / 2.0;
      *height = absz - b;
   }

/* Restore sign of latitude. */
   if ( z < 0 ) *phi = -*phi;

/* OK status. */
   return 0;

}

int eraGd2gc ( int n, double elong, double phi, double height,
               double xyz[3] )
/*
**  - - - - - - - - -
**   e r a G d 2 g c
**  - - - - - - - - -
**
**  Transform geodetic coordinates to geocentric using the specified
**  reference ellipsoid.
**
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
**
**  Called:
**     eraEform     Earth reference ellipsoids
**     eraGd2gce    geodetic to geocentric transformation, general
**     eraZp        zero p-vector
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   int j;
   double a, f;


/* Obtain reference ellipsoid parameters. */
   j = eraEform ( n, &a, &f );

/* If OK, transform longitude, geodetic latitude, height to x,y,z. */
   if ( j == 0 ) {
      j = eraGd2gce ( a, f, elong, phi, height, xyz );
      if ( j != 0 ) j = -2;
   }

/* Deal with any errors. */
   if ( j != 0 ) eraZp ( xyz );

/* Return the status. */
   return j;

}

int eraGd2gce ( double a, double f, double elong, double phi,
                double height, double xyz[3] )
/*
**  - - - - - - - - - -
**   e r a G d 2 g c e
**  - - - - - - - - - -
**
**  Transform geodetic coordinates to geocentric for a reference
**  ellipsoid of specified form.
**
**  Given:
**     a       double     equatorial radius (Notes 1,4)
**     f       double     flattening (Notes 2,4)
**     elong   double     longitude (radians, east +ve)
**     phi     double     latitude (geodetic, radians, Note 4)
**     height  double     height above ellipsoid (geodetic, Notes 3,4)
**
**  Returned:
**     xyz     double[3]  geocentric vector (Note 3)
**
**  Returned (function value):
**             int        status:  0 = OK
**                                -1 = illegal case (Note 4)
**  Notes:
**
**  1) The equatorial radius, a, can be in any units, but meters is
**     the conventional choice.
**
**  2) The flattening, f, is (for the Earth) a value around 0.00335,
**     i.e. around 1/298.
**
**  3) The equatorial radius, a, and the height, height, must be
**     given in the same units, and determine the units of the
**     returned geocentric vector, xyz.
**
**  4) No validation is performed on individual arguments.  The error
**     status -1 protects against (unrealistic) cases that would lead
**     to arithmetic exceptions.  If an error occurs, xyz is unchanged.
**
**  5) The inverse transformation is performed in the function
**     eraGc2gde.
**
**  6) The transformation for a standard ellipsoid (such as WGS84) can
**     more conveniently be performed by calling eraGd2gc,  which uses a
**     numerical code to identify the required a and f values.
**
**  References:
**
**     Green, R.M., Spherical Astronomy, Cambridge University Press,
**     (1985) Section 4.5, p96.
**
**     Explanatory Supplement to the Astronomical Almanac,
**     P. Kenneth Seidelmann (ed), University Science Books (1992),
**     Section 4.22, p202.
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double sp, cp, w, d, ac, as, r;


/* Functions of geodetic latitude. */
   sp = sin(phi);
   cp = cos(phi);
   w = 1.0 - f;
   w = w * w;
   d = cp*cp + w*sp*sp;
   if ( d <= 0.0 ) return -1;
   ac = a / sqrt(d);
   as = w * ac;

/* Geocentric vector. */
   r = (ac + height) * cp;
   xyz[0] = r * cos(elong);
   xyz[1] = r * sin(elong);
   xyz[2] = (as + height) * sp;

/* Success. */
   return 0;

}

double eraGmst00(double uta, double utb, double tta, double ttb)
/*
**  - - - - - - - - - -
**   e r a G m s t 0 0
**  - - - - - - - - - -
**
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
*/
{
   double t, gmst;


/* TT Julian centuries since J2000.0. */
   t = ((tta - DJ00) + ttb) / DJC;

/* Greenwich Mean Sidereal Time, IAU 2000. */
   gmst = eraAnp(eraEra00(uta, utb) +
                   (     0.014506   +
                   (  4612.15739966 +
                   (     1.39667721 +
                   (    -0.00009344 +
                   (     0.00001882 )
          * t) * t) * t) * t) * DAS2R);

   return gmst;

}

double eraGmst06(double uta, double utb, double tta, double ttb)
/*
**  - - - - - - - - - -
**   e r a G m s t 0 6
**  - - - - - - - - - -
**
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
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double t, gmst;


/* TT Julian centuries since J2000.0. */
   t = ((tta - DJ00) + ttb) / DJC;

/* Greenwich mean sidereal time, IAU 2006. */
   gmst = eraAnp(eraEra00(uta, utb) +
                  (    0.014506     +
                  (  4612.156534    +
                  (     1.3915817   +
                  (    -0.00000044  +
                  (    -0.000029956 +
                  (    -0.0000000368 )
          * t) * t) * t) * t) * t) * DAS2R);

   return gmst;

}

double eraGmst82(double dj1, double dj2)
/*
**  - - - - - - - - - -
**   e r a G m s t 8 2
**  - - - - - - - - - -
**
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
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
/* Coefficients of IAU 1982 GMST-UT1 model */
   double A = 24110.54841  -  DAYSEC / 2.0;
   double B = 8640184.812866;
   double C = 0.093104;
   double D =  -6.2e-6;

/* Note: the first constant, A, has to be adjusted by 12 hours */
/* because the UT1 is supplied as a Julian date, which begins  */
/* at noon.                                                    */

   double d1, d2, t, f, gmst;


/* Julian centuries since fundamental epoch. */
   if (dj1 < dj2) {
      d1 = dj1;
      d2 = dj2;
   } else {
      d1 = dj2;
      d2 = dj1;
   }
   t = (d1 + (d2 - DJ00)) / DJC;

/* Fractional part of JD(UT1), in seconds. */
   f = DAYSEC * (fmod(d1, 1.0) + fmod(d2, 1.0));

/* GMST at this UT1. */
   gmst = eraAnp(DS2R * ((A + (B + (C + D * t) * t) * t) + f));

   return gmst;

}

double eraGst00a(double uta, double utb, double tta, double ttb)
/*
**  - - - - - - - - - -
**   e r a G s t 0 0 a
**  - - - - - - - - - -
**
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
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double gmst00, ee00a, gst;


   gmst00 = eraGmst00(uta, utb, tta, ttb);
   ee00a = eraEe00a(tta, ttb);
   gst = eraAnp(gmst00 + ee00a);

   return gst;

}

double eraGst00b(double uta, double utb)
/*
**  - - - - - - - - - -
**   e r a G s t 0 0 b
**  - - - - - - - - - -
**
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
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double gmst00, ee00b, gst;


   gmst00 = eraGmst00(uta, utb, uta, utb);
   ee00b = eraEe00b(uta, utb);
   gst = eraAnp(gmst00 + ee00b);

   return gst;

}

double eraGst06(double uta, double utb, double tta, double ttb,
                double rnpb[3][3])
/*
**  - - - - - - - - -
**   e r a G s t 0 6
**  - - - - - - - - -
**
**  Greenwich apparent sidereal time, IAU 2006, given the NPB matrix.
**
**  Given:
**     uta,utb  double        UT1 as a 2-part Julian Date (Notes 1,2)
**     tta,ttb  double        TT as a 2-part Julian Date (Notes 1,2)
**     rnpb     double[3][3]  nutation x precession x bias matrix
**
**  Returned (function value):
**              double        Greenwich apparent sidereal time (radians)
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
**  3) Although the function uses the IAU 2006 series for s+XY/2, it is
**     otherwise independent of the precession-nutation model and can in
**     practice be used with any equinox-based NPB matrix.
**
**  4) The result is returned in the range 0 to 2pi.
**
**  Called:
**     eraBpn2xy    extract CIP X,Y coordinates from NPB matrix
**     eraS06       the CIO locator s, given X,Y, IAU 2006
**     eraAnp       normalize angle into range 0 to 2pi
**     eraEra00     Earth rotation angle, IAU 2000
**     eraEors      equation of the origins, given NPB matrix and s
**
**  Reference:
**
**     Wallace, P.T. & Capitaine, N., 2006, Astron.Astrophys. 459, 981
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double x, y, s, era, eors, gst;


/* Extract CIP coordinates. */
   eraBpn2xy(rnpb, &x, &y);

/* The CIO locator, s. */
   s = eraS06(tta, ttb, x, y);

/* Greenwich apparent sidereal time. */
   era = eraEra00(uta, utb);
   eors = eraEors(rnpb, s);
   gst = eraAnp(era - eors);

   return gst;

}

double eraGst06a(double uta, double utb, double tta, double ttb)
/*
**  - - - - - - - - - -
**   e r a G s t 0 6 a
**  - - - - - - - - - -
**
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
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double rnpb[3][3], gst;


/* Classical nutation x precession x bias matrix, IAU 2000A. */
   eraPnm06a(tta, ttb, rnpb);

/* Greenwich apparent sidereal time. */
   gst = eraGst06(uta, utb, tta, ttb, rnpb);

   return gst;

}

double eraGst94(double uta, double utb)
/*
**  - - - - - - - - -
**   e r a G s t 9 4
**  - - - - - - - - -
**
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
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double gmst82, eqeq94, gst;


   gmst82 = eraGmst82(uta, utb);
   eqeq94 = eraEqeq94(uta, utb);
   gst = eraAnp(gmst82  + eqeq94);

   return gst;

}

void eraH2fk5(double rh, double dh,
              double drh, double ddh, double pxh, double rvh,
              double *r5, double *d5,
              double *dr5, double *dd5, double *px5, double *rv5)
/*
**  - - - - - - - - -
**   e r a H 2 f k 5
**  - - - - - - - - -
**
**  Transform Hipparcos star data into the FK5 (J2000.0) system.
**
**  Given (all Hipparcos, epoch J2000.0):
**     rh      double    RA (radians)
**     dh      double    Dec (radians)
**     drh     double    proper motion in RA (dRA/dt, rad/Jyear)
**     ddh     double    proper motion in Dec (dDec/dt, rad/Jyear)
**     pxh     double    parallax (arcsec)
**     rvh     double    radial velocity (km/s, positive = receding)
**
**  Returned (all FK5, equinox J2000.0, epoch J2000.0):
**     r5      double    RA (radians)
**     d5      double    Dec (radians)
**     dr5     double    proper motion in RA (dRA/dt, rad/Jyear)
**     dd5     double    proper motion in Dec (dDec/dt, rad/Jyear)
**     px5     double    parallax (arcsec)
**     rv5     double    radial velocity (km/s, positive = receding)
**
**  Notes:
**
**  1) This function transforms Hipparcos star positions and proper
**     motions into FK5 J2000.0.
**
**  2) The proper motions in RA are dRA/dt rather than
**     cos(Dec)*dRA/dt, and are per year rather than per century.
**
**  3) The FK5 to Hipparcos transformation is modeled as a pure
**     rotation and spin;  zonal errors in the FK5 catalog are not
**     taken into account.
**
**  4) See also eraFk52h, eraFk5hz, eraHfk5z.
**
**  Called:
**     eraStarpv    star catalog data to space motion pv-vector
**     eraFk5hip    FK5 to Hipparcos rotation and spin
**     eraRv2m      r-vector to r-matrix
**     eraRxp       product of r-matrix and p-vector
**     eraTrxp      product of transpose of r-matrix and p-vector
**     eraPxp       vector product of two p-vectors
**     eraPmp       p-vector minus p-vector
**     eraPvstar    space motion pv-vector to star catalog data
**
**  Reference:
**
**     F.Mignard & M.Froeschle, Astron. Astrophys. 354, 732-739 (2000).
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   int i;
   double pvh[2][3], r5h[3][3], s5h[3], sh[3], wxp[3], vv[3], pv5[2][3];


/* Hipparcos barycentric position/velocity pv-vector (normalized). */
   eraStarpv(rh, dh, drh, ddh, pxh, rvh, pvh);

/* FK5 to Hipparcos orientation matrix and spin vector. */
   eraFk5hip(r5h, s5h);

/* Make spin units per day instead of per year. */
   for ( i = 0; i < 3; s5h[i++] /= 365.25 );

/* Orient the spin into the Hipparcos system. */
   eraRxp(r5h, s5h, sh);

/* De-orient the Hipparcos position into the FK5 system. */
   eraTrxp(r5h, pvh[0], pv5[0]);

/* Apply spin to the position giving an extra space motion component. */
   eraPxp(pvh[0], sh, wxp);

/* Subtract this component from the Hipparcos space motion. */
   eraPmp(pvh[1], wxp, vv);

/* De-orient the Hipparcos space motion into the FK5 system. */
   eraTrxp(r5h, vv, pv5[1]);

/* FK5 pv-vector to spherical. */
   eraPvstar(pv5, r5, d5, dr5, dd5, px5, rv5);

   return;

}

void eraHfk5z(double rh, double dh, double date1, double date2,
              double *r5, double *d5, double *dr5, double *dd5)
/*
**  - - - - - - - - -
**   e r a H f k 5 z
**  - - - - - - - - -
**
**  Transform a Hipparcos star position into FK5 J2000.0, assuming
**  zero Hipparcos proper motion.
**
**  Given:
**     rh            double    Hipparcos RA (radians)
**     dh            double    Hipparcos Dec (radians)
**     date1,date2   double    TDB date (Note 1)
**
**  Returned (all FK5, equinox J2000.0, date date1+date2):
**     r5            double    RA (radians)
**     d5            double    Dec (radians)
**     dr5           double    FK5 RA proper motion (rad/year, Note 4)
**     dd5           double    Dec proper motion (rad/year, Note 4)
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
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
**  2) The proper motion in RA is dRA/dt rather than cos(Dec)*dRA/dt.
**
**  3) The FK5 to Hipparcos transformation is modeled as a pure rotation
**     and spin;  zonal errors in the FK5 catalogue are not taken into
**     account.
**
**  4) It was the intention that Hipparcos should be a close
**     approximation to an inertial frame, so that distant objects have
**     zero proper motion;  such objects have (in general) non-zero
**     proper motion in FK5, and this function returns those fictitious
**     proper motions.
**
**  5) The position returned by this function is in the FK5 J2000.0
**     reference system but at date date1+date2.
**
**  6) See also eraFk52h, eraH2fk5, eraFk5zhz.
**
**  Called:
**     eraS2c       spherical coordinates to unit vector
**     eraFk5hip    FK5 to Hipparcos rotation and spin
**     eraRxp       product of r-matrix and p-vector
**     eraSxp       multiply p-vector by scalar
**     eraRxr       product of two r-matrices
**     eraTrxp      product of transpose of r-matrix and p-vector
**     eraPxp       vector product of two p-vectors
**     eraPv2s      pv-vector to spherical
**     eraAnp       normalize angle into range 0 to 2pi
**
**  Reference:
**
**     F.Mignard & M.Froeschle, 2000, Astron.Astrophys. 354, 732-739.
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double t, ph[3], r5h[3][3], s5h[3], sh[3], vst[3],
   rst[3][3], r5ht[3][3], pv5e[2][3], vv[3],
   w, r, v;


/* Time interval from fundamental epoch J2000.0 to given date (JY). */
   t = ((date1 - DJ00) + date2) / DJY;

/* Hipparcos barycentric position vector (normalized). */
   eraS2c(rh, dh, ph);

/* FK5 to Hipparcos orientation matrix and spin vector. */
   eraFk5hip(r5h, s5h);

/* Rotate the spin into the Hipparcos system. */
   eraRxp(r5h, s5h, sh);

/* Accumulated Hipparcos wrt FK5 spin over that interval. */
   eraSxp(t, s5h, vst);

/* Express the accumulated spin as a rotation matrix. */
   eraRv2m(vst, rst);

/* Rotation matrix:  accumulated spin, then FK5 to Hipparcos. */
   eraRxr(r5h, rst, r5ht);

/* De-orient & de-spin the Hipparcos position into FK5 J2000.0. */
   eraTrxp(r5ht, ph, pv5e[0]);

/* Apply spin to the position giving a space motion. */
   eraPxp(sh, ph, vv);

/* De-orient & de-spin the Hipparcos space motion into FK5 J2000.0. */
   eraTrxp(r5ht, vv, pv5e[1]);

/* FK5 position/velocity pv-vector to spherical. */
   eraPv2s(pv5e, &w, d5, &r, dr5, dd5, &v);
   *r5 = eraAnp(w);

   return;

}

void eraIr(double r[3][3])
/*
**  - - - - - -
**   e r a I r
**  - - - - - -
**
**  Initialize an r-matrix to the identity matrix.
**
**  Returned:
**     r       double[3][3]    r-matrix
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   r[0][0] = 1.0;
   r[0][1] = 0.0;
   r[0][2] = 0.0;
   r[1][0] = 0.0;
   r[1][1] = 1.0;
   r[1][2] = 0.0;
   r[2][0] = 0.0;
   r[2][1] = 0.0;
   r[2][2] = 1.0;

   return;

}

int eraJd2cal(double dj1, double dj2,
              int *iy, int *im, int *id, double *fd)
/*
**  - - - - - - - - - -
**   e r a J d 2 c a l
**  - - - - - - - - - -
**
**  Julian Date to Gregorian year, month, day, and fraction of a day.
**
**  Given:
**     dj1,dj2   double   Julian Date (Notes 1, 2)
**
**  Returned (arguments):
**     iy        int      year
**     im        int      month
**     id        int      day
**     fd        double   fraction of day
**
**  Returned (function value):
**               int      status:
**                           0 = OK
**                          -1 = unacceptable date (Note 3)
**
**  Notes:
**
**  1) The earliest valid date is -68569.5 (-4900 March 1).  The
**     largest value accepted is 10^9.
**
**  2) The Julian Date is apportioned in any convenient way between
**     the arguments dj1 and dj2.  For example, JD=2450123.7 could
**     be expressed in any of these ways, among others:
**
**            dj1             dj2
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**  3) In early eras the conversion is from the "proleptic Gregorian
**     calendar";  no account is taken of the date(s) of adoption of
**     the Gregorian calendar, nor is the AD/BC numbering convention
**     observed.
**
**  Reference:
**
**     Explanatory Supplement to the Astronomical Almanac,
**     P. Kenneth Seidelmann (ed), University Science Books (1992),
**     Section 12.92 (p604).
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
/* Minimum and maximum allowed JD */
   static const double djmin = -68569.5;
   static const double djmax = 1e9;

   long jd, l, n, i, k;
   double dj, d1, d2, f1, f2, f, d;


/* Verify date is acceptable. */
   dj = dj1 + dj2;
   if (dj < djmin || dj > djmax) return -1;

/* Copy the date, big then small, and re-align to midnight. */
   if (dj1 >= dj2) {
      d1 = dj1;
      d2 = dj2;
   } else {
      d1 = dj2;
      d2 = dj1;
   }
   d2 -= 0.5;

/* Separate day and fraction. */
   f1 = fmod(d1, 1.0);
   f2 = fmod(d2, 1.0);
   f = fmod(f1 + f2, 1.0);
   if (f < 0.0) f += 1.0;
   d = floor(d1 - f1) + floor(d2 - f2) + floor(f1 + f2 - f);
   jd = (long) floor(d) + 1L;

/* Express day in Gregorian calendar. */
   l = jd + 68569L;
   n = (4L * l) / 146097L;
   l -= (146097L * n + 3L) / 4L;
   i = (4000L * (l + 1L)) / 1461001L;
   l -= (1461L * i) / 4L - 31L;
   k = (80L * l) / 2447L;
   *id = (int) (l - (2447L * k) / 80L);
   l = k / 11L;
   *im = (int) (k + 2L - 12L * l);
   *iy = (int) (100L * (n - 49L) + i + l);
   *fd = f;

   return 0;

}

int eraJdcalf(int ndp, double dj1, double dj2, int iymdf[4])
/*
**  - - - - - - - - - -
**   e r a J d c a l f
**  - - - - - - - - - -
**
**  Julian Date to Gregorian Calendar, expressed in a form convenient
**  for formatting messages:  rounded to a specified precision.
**
**  Given:
**     ndp       int      number of decimal places of days in fraction
**     dj1,dj2   double   dj1+dj2 = Julian Date (Note 1)
**
**  Returned:
**     iymdf     int[4]   year, month, day, fraction in Gregorian
**                        calendar
**
**  Returned (function value):
**               int      status:
**                          -1 = date out of range
**                           0 = OK
**                          +1 = NDP not 0-9 (interpreted as 0)
**
**  Notes:
**
**  1) The Julian Date is apportioned in any convenient way between
**     the arguments dj1 and dj2.  For example, JD=2450123.7 could
**     be expressed in any of these ways, among others:
**
**             dj1            dj2
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**  2) In early eras the conversion is from the "Proleptic Gregorian
**     Calendar";  no account is taken of the date(s) of adoption of
**     the Gregorian Calendar, nor is the AD/BC numbering convention
**     observed.
**
**  3) Refer to the function eraJd2cal.
**
**  4) NDP should be 4 or less if internal overflows are to be
**     avoided on machines which use 16-bit integers.
**
**  Called:
**     eraJd2cal    JD to Gregorian calendar
**
**  Reference:
**
**     Explanatory Supplement to the Astronomical Almanac,
**     P. Kenneth Seidelmann (ed), University Science Books (1992),
**     Section 12.92 (p604).
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   int j, js;
   double denom, d1, d2, f1, f2, f;


/* Denominator of fraction (e.g. 100 for 2 decimal places). */
   if ((ndp >= 0) && (ndp <= 9)) {
      j = 0;
      denom = pow(10.0, ndp);
   } else {
      j = 1;
      denom = 1.0;
   }

/* Copy the date, big then small, and realign to midnight. */
   if (dj1 >= dj2) {
      d1 = dj1;
      d2 = dj2;
   } else {
      d1 = dj2;
      d2 = dj1;
   }
   d2 -= 0.5;

/* Separate days and fractions. */
   f1 = fmod(d1, 1.0);
   f2 = fmod(d2, 1.0);
   d1 = floor(d1 - f1);
   d2 = floor(d2 - f2);

/* Round the total fraction to the specified number of places. */
   f = floor((f1+f2)*denom + 0.5) / denom;

/* Re-assemble the rounded date and re-align to noon. */
   d2 += f + 0.5;

/* Convert to Gregorian calendar. */
   js = eraJd2cal(d1, d2, &iymdf[0], &iymdf[1], &iymdf[2], &f);
   if (js == 0) {
      iymdf[3] = (int) (f * denom);
   } else {
      j = js;
   }

/* Return the status. */
   return j;

}

void eraNum00a(double date1, double date2, double rmatn[3][3])
/*
**  - - - - - - - - - -
**   e r a N u m 0 0 a
**  - - - - - - - - - -
**
**  Form the matrix of nutation for a given date, IAU 2000A model.
**
**  Given:
**     date1,date2  double          TT as a 2-part Julian Date (Note 1)
**
**  Returned:
**     rmatn        double[3][3]    nutation matrix
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
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
**  2) The matrix operates in the sense V(true) = rmatn * V(mean), where
**     the p-vector V(true) is with respect to the true equatorial triad
**     of date and the p-vector V(mean) is with respect to the mean
**     equatorial triad of date.
**
**  3) A faster, but slightly less accurate result (about 1 mas), can be
**     obtained by using instead the eraNum00b function.
**
**  Called:
**     eraPn00a     bias/precession/nutation, IAU 2000A
**
**  Reference:
**
**     Explanatory Supplement to the Astronomical Almanac,
**     P. Kenneth Seidelmann (ed), University Science Books (1992),
**     Section 3.222-3 (p114).
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double dpsi, deps, epsa, rb[3][3], rp[3][3], rbp[3][3], rbpn[3][3];


/* Obtain the required matrix (discarding other results). */
   eraPn00a(date1, date2,
            &dpsi, &deps, &epsa, rb, rp, rbp, rmatn, rbpn);

   return;

}

void eraNum00b(double date1, double date2, double rmatn[3][3])
/*
**  - - - - - - - - - -
**   e r a N u m 0 0 b
**  - - - - - - - - - -
**
**  Form the matrix of nutation for a given date, IAU 2000B model.
**
**  Given:
**     date1,date2  double         TT as a 2-part Julian Date (Note 1)
**
**  Returned:
**     rmatn        double[3][3]   nutation matrix
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
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
**  2) The matrix operates in the sense V(true) = rmatn * V(mean), where
**     the p-vector V(true) is with respect to the true equatorial triad
**     of date and the p-vector V(mean) is with respect to the mean
**     equatorial triad of date.
**
**  3) The present function is faster, but slightly less accurate (about
**     1 mas), than the eraNum00a function.
**
**  Called:
**     eraPn00b     bias/precession/nutation, IAU 2000B
**
**  Reference:
**
**     Explanatory Supplement to the Astronomical Almanac,
**     P. Kenneth Seidelmann (ed), University Science Books (1992),
**     Section 3.222-3 (p114).
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double dpsi, deps, epsa, rb[3][3], rp[3][3], rbp[3][3], rbpn[3][3];


/* Obtain the required matrix (discarding other results). */
   eraPn00b(date1, date2,
            &dpsi, &deps, &epsa, rb, rp, rbp, rmatn, rbpn);

   return;

}

void eraNum06a(double date1, double date2, double rmatn[3][3])
/*
**  - - - - - - - - - -
**   e r a N u m 0 6 a
**  - - - - - - - - - -
**
**  Form the matrix of nutation for a given date, IAU 2006/2000A model.
**
**  Given:
**     date1,date2   double          TT as a 2-part Julian Date (Note 1)
**
**  Returned:
**     rmatn         double[3][3]    nutation matrix
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
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
**  2) The matrix operates in the sense V(true) = rmatn * V(mean), where
**     the p-vector V(true) is with respect to the true equatorial triad
**     of date and the p-vector V(mean) is with respect to the mean
**     equatorial triad of date.
**
**  Called:
**     eraObl06     mean obliquity, IAU 2006
**     eraNut06a    nutation, IAU 2006/2000A
**     eraNumat     form nutation matrix
**
**  Reference:
**
**     Explanatory Supplement to the Astronomical Almanac,
**     P. Kenneth Seidelmann (ed), University Science Books (1992),
**     Section 3.222-3 (p114).
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double eps, dp, de;


/* Mean obliquity. */
   eps = eraObl06(date1, date2);

/* Nutation components. */
   eraNut06a(date1, date2, &dp, &de);

/* Nutation matrix. */
   eraNumat(eps, dp, de, rmatn);

   return;

}

void eraNumat(double epsa, double dpsi, double deps, double rmatn[3][3])
/*
**  - - - - - - - - -
**   e r a N u m a t
**  - - - - - - - - -
**
**  Form the matrix of nutation.
**
**  Given:
**     epsa        double         mean obliquity of date (Note 1)
**     dpsi,deps   double         nutation (Note 2)
**
**  Returned:
**     rmatn       double[3][3]   nutation matrix (Note 3)
**
**  Notes:
**
**
**  1) The supplied mean obliquity epsa, must be consistent with the
**     precession-nutation models from which dpsi and deps were obtained.
**
**  2) The caller is responsible for providing the nutation components;
**     they are in longitude and obliquity, in radians and are with
**     respect to the equinox and ecliptic of date.
**
**  3) The matrix operates in the sense V(true) = rmatn * V(mean),
**     where the p-vector V(true) is with respect to the true
**     equatorial triad of date and the p-vector V(mean) is with
**     respect to the mean equatorial triad of date.
**
**  Called:
**     eraIr        initialize r-matrix to identity
**     eraRx        rotate around X-axis
**     eraRz        rotate around Z-axis
**
**  Reference:
**
**     Explanatory Supplement to the Astronomical Almanac,
**     P. Kenneth Seidelmann (ed), University Science Books (1992),
**     Section 3.222-3 (p114).
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
/* Build the rotation matrix. */
   eraIr(rmatn);
   eraRx(epsa, rmatn);
   eraRz(-dpsi, rmatn);
   eraRx(-(epsa + deps), rmatn);

   return;

}

void eraNut00a(double date1, double date2, double *dpsi, double *deps)
/*
**  - - - - - - - - - -
**   e r a N u t 0 0 a
**  - - - - - - - - - -
**
**  Nutation, IAU 2000A model (MHB2000 luni-solar and planetary nutation
**  with free core nutation omitted).
**
**  Given:
**     date1,date2   double   TT as a 2-part Julian Date (Note 1)
**
**  Returned:
**     dpsi,deps     double   nutation, luni-solar + planetary (Note 2)
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
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
**  2) The nutation components in longitude and obliquity are in radians
**     and with respect to the equinox and ecliptic of date.  The
**     obliquity at J2000.0 is assumed to be the Lieske et al. (1977)
**     value of 84381.448 arcsec.
**
**     Both the luni-solar and planetary nutations are included.  The
**     latter are due to direct planetary nutations and the
**     perturbations of the lunar and terrestrial orbits.
**
**  3) The function computes the MHB2000 nutation series with the
**     associated corrections for planetary nutations.  It is an
**     implementation of the nutation part of the IAU 2000A precession-
**     nutation model, formally adopted by the IAU General Assembly in
**     2000, namely MHB2000 (Mathews et al. 2002), but with the free
**     core nutation (FCN - see Note 4) omitted.
**
**  4) The full MHB2000 model also contains contributions to the
**     nutations in longitude and obliquity due to the free-excitation
**     of the free-core-nutation during the period 1979-2000.  These FCN
**     terms, which are time-dependent and unpredictable, are NOT
**     included in the present function and, if required, must be
**     independently computed.  With the FCN corrections included, the
**     present function delivers a pole which is at current epochs
**     accurate to a few hundred microarcseconds.  The omission of FCN
**     introduces further errors of about that size.
**
**  5) The present function provides classical nutation.  The MHB2000
**     algorithm, from which it is adapted, deals also with (i) the
**     offsets between the GCRS and mean poles and (ii) the adjustments
**     in longitude and obliquity due to the changed precession rates.
**     These additional functions, namely frame bias and precession
**     adjustments, are supported by the ERFA functions eraBi00  and
**     eraPr00.
**
**  6) The MHB2000 algorithm also provides "total" nutations, comprising
**     the arithmetic sum of the frame bias, precession adjustments,
**     luni-solar nutation and planetary nutation.  These total
**     nutations can be used in combination with an existing IAU 1976
**     precession implementation, such as eraPmat76,  to deliver GCRS-
**     to-true predictions of sub-mas accuracy at current dates.
**     However, there are three shortcomings in the MHB2000 model that
**     must be taken into account if more accurate or definitive results
**     are required (see Wallace 2002):
**
**       (i) The MHB2000 total nutations are simply arithmetic sums,
**           yet in reality the various components are successive Euler
**           rotations.  This slight lack of rigor leads to cross terms
**           that exceed 1 mas after a century.  The rigorous procedure
**           is to form the GCRS-to-true rotation matrix by applying the
**           bias, precession and nutation in that order.
**
**      (ii) Although the precession adjustments are stated to be with
**           respect to Lieske et al. (1977), the MHB2000 model does
**           not specify which set of Euler angles are to be used and
**           how the adjustments are to be applied.  The most literal
**           and straightforward procedure is to adopt the 4-rotation
**           epsilon_0, psi_A, omega_A, xi_A option, and to add DPSIPR
**           to psi_A and DEPSPR to both omega_A and eps_A.
**
**     (iii) The MHB2000 model predates the determination by Chapront
**           et al. (2002) of a 14.6 mas displacement between the
**           J2000.0 mean equinox and the origin of the ICRS frame.  It
**           should, however, be noted that neglecting this displacement
**           when calculating star coordinates does not lead to a
**           14.6 mas change in right ascension, only a small second-
**           order distortion in the pattern of the precession-nutation
**           effect.
**
**     For these reasons, the ERFA functions do not generate the "total
**     nutations" directly, though they can of course easily be
**     generated by calling eraBi00, eraPr00 and the present function
**     and adding the results.
**
**  7) The MHB2000 model contains 41 instances where the same frequency
**     appears multiple times, of which 38 are duplicates and three are
**     triplicates.  To keep the present code close to the original MHB
**     algorithm, this small inefficiency has not been corrected.
**
**  Called:
**     eraFal03     mean anomaly of the Moon
**     eraFaf03     mean argument of the latitude of the Moon
**     eraFaom03    mean longitude of the Moon's ascending node
**     eraFame03    mean longitude of Mercury
**     eraFave03    mean longitude of Venus
**     eraFae03     mean longitude of Earth
**     eraFama03    mean longitude of Mars
**     eraFaju03    mean longitude of Jupiter
**     eraFasa03    mean longitude of Saturn
**     eraFaur03    mean longitude of Uranus
**     eraFapa03    general accumulated precession in longitude
**
**  References:
**
**     Chapront, J., Chapront-Touze, M. & Francou, G. 2002,
**     Astron.Astrophys. 387, 700
**
**     Lieske, J.H., Lederle, T., Fricke, W. & Morando, B. 1977,
**     Astron.Astrophys. 58, 1-16
**
**     Mathews, P.M., Herring, T.A., Buffet, B.A. 2002, J.Geophys.Res.
**     107, B4.  The MHB_2000 code itself was obtained on 9th September
**     2002 from ftp//maia.usno.navy.mil/conv2000/chapter5/IAU2000A.
**
**     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
**     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
**
**     Souchay, J., Loysel, B., Kinoshita, H., Folgueira, M. 1999,
**     Astron.Astrophys.Supp.Ser. 135, 111
**
**     Wallace, P.T., "Software for Implementing the IAU 2000
**     Resolutions", in IERS Workshop 5.1 (2002)
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   int i;
   double t, el, elp, f, d, om, arg, dp, de, sarg, carg,
          al, af, ad, aom, alme, alve, alea, alma,
          alju, alsa, alur, alne, apa, dpsils, depsls,
          dpsipl, depspl;

/* Units of 0.1 microarcsecond to radians */
   const double U2R = DAS2R / 1e7;

/* ------------------------- */
/* Luni-Solar nutation model */
/* ------------------------- */

/* The units for the sine and cosine coefficients are */
/* 0.1 microarcsecond and the same per Julian century */

   static const struct {
      int nl,nlp,nf,nd,nom; /* coefficients of l,l',F,D,Om */
      double sp,spt,cp;     /* longitude sin, t*sin, cos coefficients */
      double ce,cet,se;     /* obliquity cos, t*cos, sin coefficients */
   } xls[] = {

   /* 1- 10 */
      { 0, 0, 0, 0, 1,
         -172064161.0, -174666.0, 33386.0, 92052331.0, 9086.0, 15377.0},
      { 0, 0, 2,-2, 2,
           -13170906.0, -1675.0, -13696.0, 5730336.0, -3015.0, -4587.0},
      { 0, 0, 2, 0, 2,-2276413.0,-234.0,2796.0,978459.0,-485.0, 1374.0},
      { 0, 0, 0, 0, 2,2074554.0, 207.0, -698.0,-897492.0,470.0, -291.0},
      { 0, 1, 0, 0, 0,1475877.0,-3633.0,11817.0,73871.0,-184.0,-1924.0},
      { 0, 1, 2,-2, 2,-516821.0,1226.0, -524.0,224386.0,-677.0, -174.0},
      { 1, 0, 0, 0, 0, 711159.0,  73.0, -872.0,  -6750.0,  0.0,  358.0},
      { 0, 0, 2, 0, 1,-387298.0,-367.0,  380.0, 200728.0, 18.0,  318.0},
      { 1, 0, 2, 0, 2,-301461.0, -36.0,  816.0, 129025.0,-63.0,  367.0},
      { 0,-1, 2,-2, 2, 215829.0,-494.0,  111.0, -95929.0,299.0,  132.0},

   /* 11-20 */
      { 0, 0, 2,-2, 1, 128227.0, 137.0,  181.0, -68982.0, -9.0,   39.0},
      {-1, 0, 2, 0, 2, 123457.0,  11.0,   19.0, -53311.0, 32.0,   -4.0},
      {-1, 0, 0, 2, 0, 156994.0,  10.0, -168.0,  -1235.0,  0.0,   82.0},
      { 1, 0, 0, 0, 1,  63110.0,  63.0,   27.0, -33228.0,  0.0,   -9.0},
      {-1, 0, 0, 0, 1, -57976.0, -63.0, -189.0,  31429.0,  0.0,  -75.0},
      {-1, 0, 2, 2, 2, -59641.0, -11.0,  149.0,  25543.0,-11.0,   66.0},
      { 1, 0, 2, 0, 1, -51613.0, -42.0,  129.0,  26366.0,  0.0,   78.0},
      {-2, 0, 2, 0, 1,  45893.0,  50.0,   31.0, -24236.0,-10.0,   20.0},
      { 0, 0, 0, 2, 0,  63384.0,  11.0, -150.0,  -1220.0,  0.0,   29.0},
      { 0, 0, 2, 2, 2, -38571.0,  -1.0,  158.0,  16452.0,-11.0,   68.0},

   /* 21-30 */
      { 0,-2, 2,-2, 2,  32481.0,   0.0,    0.0, -13870.0,  0.0,    0.0},
      {-2, 0, 0, 2, 0, -47722.0,   0.0,  -18.0,    477.0,  0.0,  -25.0},
      { 2, 0, 2, 0, 2, -31046.0,  -1.0,  131.0,  13238.0,-11.0,   59.0},
      { 1, 0, 2,-2, 2,  28593.0,   0.0,   -1.0, -12338.0, 10.0,   -3.0},
      {-1, 0, 2, 0, 1,  20441.0,  21.0,   10.0, -10758.0,  0.0,   -3.0},
      { 2, 0, 0, 0, 0,  29243.0,   0.0,  -74.0,   -609.0,  0.0,   13.0},
      { 0, 0, 2, 0, 0,  25887.0,   0.0,  -66.0,   -550.0,  0.0,   11.0},
      { 0, 1, 0, 0, 1, -14053.0, -25.0,   79.0,   8551.0, -2.0,  -45.0},
      {-1, 0, 0, 2, 1,  15164.0,  10.0,   11.0,  -8001.0,  0.0,   -1.0},
      { 0, 2, 2,-2, 2, -15794.0,  72.0,  -16.0,   6850.0,-42.0,   -5.0},

   /* 31-40 */
      { 0, 0,-2, 2, 0,  21783.0,   0.0,   13.0,   -167.0,  0.0,   13.0},
      { 1, 0, 0,-2, 1, -12873.0, -10.0,  -37.0,   6953.0,  0.0,  -14.0},
      { 0,-1, 0, 0, 1, -12654.0,  11.0,   63.0,   6415.0,  0.0,   26.0},
      {-1, 0, 2, 2, 1, -10204.0,   0.0,   25.0,   5222.0,  0.0,   15.0},
      { 0, 2, 0, 0, 0,  16707.0, -85.0,  -10.0,    168.0, -1.0,   10.0},
      { 1, 0, 2, 2, 2,  -7691.0,   0.0,   44.0,   3268.0,  0.0,   19.0},
      {-2, 0, 2, 0, 0, -11024.0,   0.0,  -14.0,    104.0,  0.0,    2.0},
      { 0, 1, 2, 0, 2,   7566.0, -21.0,  -11.0,  -3250.0,  0.0,   -5.0},
      { 0, 0, 2, 2, 1,  -6637.0, -11.0,   25.0,   3353.0,  0.0,   14.0},
      { 0,-1, 2, 0, 2,  -7141.0,  21.0,    8.0,   3070.0,  0.0,    4.0},

   /* 41-50 */
      { 0, 0, 0, 2, 1,  -6302.0, -11.0,    2.0,   3272.0,  0.0,    4.0},
      { 1, 0, 2,-2, 1,   5800.0,  10.0,    2.0,  -3045.0,  0.0,   -1.0},
      { 2, 0, 2,-2, 2,   6443.0,   0.0,   -7.0,  -2768.0,  0.0,   -4.0},
      {-2, 0, 0, 2, 1,  -5774.0, -11.0,  -15.0,   3041.0,  0.0,   -5.0},
      { 2, 0, 2, 0, 1,  -5350.0,   0.0,   21.0,   2695.0,  0.0,   12.0},
      { 0,-1, 2,-2, 1,  -4752.0, -11.0,   -3.0,   2719.0,  0.0,   -3.0},
      { 0, 0, 0,-2, 1,  -4940.0, -11.0,  -21.0,   2720.0,  0.0,   -9.0},
      {-1,-1, 0, 2, 0,   7350.0,   0.0,   -8.0,    -51.0,  0.0,    4.0},
      { 2, 0, 0,-2, 1,   4065.0,   0.0,    6.0,  -2206.0,  0.0,    1.0},
      { 1, 0, 0, 2, 0,   6579.0,   0.0,  -24.0,   -199.0,  0.0,    2.0},

   /* 51-60 */
      { 0, 1, 2,-2, 1,   3579.0,   0.0,    5.0,  -1900.0,  0.0,    1.0},
      { 1,-1, 0, 0, 0,   4725.0,   0.0,   -6.0,    -41.0,  0.0,    3.0},
      {-2, 0, 2, 0, 2,  -3075.0,   0.0,   -2.0,   1313.0,  0.0,   -1.0},
      { 3, 0, 2, 0, 2,  -2904.0,   0.0,   15.0,   1233.0,  0.0,    7.0},
      { 0,-1, 0, 2, 0,   4348.0,   0.0,  -10.0,    -81.0,  0.0,    2.0},
      { 1,-1, 2, 0, 2,  -2878.0,   0.0,    8.0,   1232.0,  0.0,    4.0},
      { 0, 0, 0, 1, 0,  -4230.0,   0.0,    5.0,    -20.0,  0.0,   -2.0},
      {-1,-1, 2, 2, 2,  -2819.0,   0.0,    7.0,   1207.0,  0.0,    3.0},
      {-1, 0, 2, 0, 0,  -4056.0,   0.0,    5.0,     40.0,  0.0,   -2.0},
      { 0,-1, 2, 2, 2,  -2647.0,   0.0,   11.0,   1129.0,  0.0,    5.0},

   /* 61-70 */
      {-2, 0, 0, 0, 1,  -2294.0,   0.0,  -10.0,   1266.0,  0.0,   -4.0},
      { 1, 1, 2, 0, 2,   2481.0,   0.0,   -7.0,  -1062.0,  0.0,   -3.0},
      { 2, 0, 0, 0, 1,   2179.0,   0.0,   -2.0,  -1129.0,  0.0,   -2.0},
      {-1, 1, 0, 1, 0,   3276.0,   0.0,    1.0,     -9.0,  0.0,    0.0},
      { 1, 1, 0, 0, 0,  -3389.0,   0.0,    5.0,     35.0,  0.0,   -2.0},
      { 1, 0, 2, 0, 0,   3339.0,   0.0,  -13.0,   -107.0,  0.0,    1.0},
      {-1, 0, 2,-2, 1,  -1987.0,   0.0,   -6.0,   1073.0,  0.0,   -2.0},
      { 1, 0, 0, 0, 2,  -1981.0,   0.0,    0.0,    854.0,  0.0,    0.0},
      {-1, 0, 0, 1, 0,   4026.0,   0.0, -353.0,   -553.0,  0.0, -139.0},
      { 0, 0, 2, 1, 2,   1660.0,   0.0,   -5.0,   -710.0,  0.0,   -2.0},

   /* 71-80 */
      {-1, 0, 2, 4, 2,  -1521.0,   0.0,    9.0,    647.0,  0.0,    4.0},
      {-1, 1, 0, 1, 1,   1314.0,   0.0,    0.0,   -700.0,  0.0,    0.0},
      { 0,-2, 2,-2, 1,  -1283.0,   0.0,    0.0,    672.0,  0.0,    0.0},
      { 1, 0, 2, 2, 1,  -1331.0,   0.0,    8.0,    663.0,  0.0,    4.0},
      {-2, 0, 2, 2, 2,   1383.0,   0.0,   -2.0,   -594.0,  0.0,   -2.0},
      {-1, 0, 0, 0, 2,   1405.0,   0.0,    4.0,   -610.0,  0.0,    2.0},
      { 1, 1, 2,-2, 2,   1290.0,   0.0,    0.0,   -556.0,  0.0,    0.0},
      {-2, 0, 2, 4, 2,  -1214.0,   0.0,    5.0,    518.0,  0.0,    2.0},
      {-1, 0, 4, 0, 2,   1146.0,   0.0,   -3.0,   -490.0,  0.0,   -1.0},
      { 2, 0, 2,-2, 1,   1019.0,   0.0,   -1.0,   -527.0,  0.0,   -1.0},

   /* 81-90 */
      { 2, 0, 2, 2, 2,  -1100.0,   0.0,    9.0,    465.0,  0.0,    4.0},
      { 1, 0, 0, 2, 1,   -970.0,   0.0,    2.0,    496.0,  0.0,    1.0},
      { 3, 0, 0, 0, 0,   1575.0,   0.0,   -6.0,    -50.0,  0.0,    0.0},
      { 3, 0, 2,-2, 2,    934.0,   0.0,   -3.0,   -399.0,  0.0,   -1.0},
      { 0, 0, 4,-2, 2,    922.0,   0.0,   -1.0,   -395.0,  0.0,   -1.0},
      { 0, 1, 2, 0, 1,    815.0,   0.0,   -1.0,   -422.0,  0.0,   -1.0},
      { 0, 0,-2, 2, 1,    834.0,   0.0,    2.0,   -440.0,  0.0,    1.0},
      { 0, 0, 2,-2, 3,   1248.0,   0.0,    0.0,   -170.0,  0.0,    1.0},
      {-1, 0, 0, 4, 0,   1338.0,   0.0,   -5.0,    -39.0,  0.0,    0.0},
      { 2, 0,-2, 0, 1,    716.0,   0.0,   -2.0,   -389.0,  0.0,   -1.0},

   /* 91-100 */
      {-2, 0, 0, 4, 0,   1282.0,   0.0,   -3.0,    -23.0,  0.0,    1.0},
      {-1,-1, 0, 2, 1,    742.0,   0.0,    1.0,   -391.0,  0.0,    0.0},
      {-1, 0, 0, 1, 1,   1020.0,   0.0,  -25.0,   -495.0,  0.0,  -10.0},
      { 0, 1, 0, 0, 2,    715.0,   0.0,   -4.0,   -326.0,  0.0,    2.0},
      { 0, 0,-2, 0, 1,   -666.0,   0.0,   -3.0,    369.0,  0.0,   -1.0},
      { 0,-1, 2, 0, 1,   -667.0,   0.0,    1.0,    346.0,  0.0,    1.0},
      { 0, 0, 2,-1, 2,   -704.0,   0.0,    0.0,    304.0,  0.0,    0.0},
      { 0, 0, 2, 4, 2,   -694.0,   0.0,    5.0,    294.0,  0.0,    2.0},
      {-2,-1, 0, 2, 0,  -1014.0,   0.0,   -1.0,      4.0,  0.0,   -1.0},
      { 1, 1, 0,-2, 1,   -585.0,   0.0,   -2.0,    316.0,  0.0,   -1.0},

   /* 101-110 */
      {-1, 1, 0, 2, 0,   -949.0,   0.0,    1.0,      8.0,  0.0,   -1.0},
      {-1, 1, 0, 1, 2,   -595.0,   0.0,    0.0,    258.0,  0.0,    0.0},
      { 1,-1, 0, 0, 1,    528.0,   0.0,    0.0,   -279.0,  0.0,    0.0},
      { 1,-1, 2, 2, 2,   -590.0,   0.0,    4.0,    252.0,  0.0,    2.0},
      {-1, 1, 2, 2, 2,    570.0,   0.0,   -2.0,   -244.0,  0.0,   -1.0},
      { 3, 0, 2, 0, 1,   -502.0,   0.0,    3.0,    250.0,  0.0,    2.0},
      { 0, 1,-2, 2, 0,   -875.0,   0.0,    1.0,     29.0,  0.0,    0.0},
      {-1, 0, 0,-2, 1,   -492.0,   0.0,   -3.0,    275.0,  0.0,   -1.0},
      { 0, 1, 2, 2, 2,    535.0,   0.0,   -2.0,   -228.0,  0.0,   -1.0},
      {-1,-1, 2, 2, 1,   -467.0,   0.0,    1.0,    240.0,  0.0,    1.0},

   /* 111-120 */
      { 0,-1, 0, 0, 2,    591.0,   0.0,    0.0,   -253.0,  0.0,    0.0},
      { 1, 0, 2,-4, 1,   -453.0,   0.0,   -1.0,    244.0,  0.0,   -1.0},
      {-1, 0,-2, 2, 0,    766.0,   0.0,    1.0,      9.0,  0.0,    0.0},
      { 0,-1, 2, 2, 1,   -446.0,   0.0,    2.0,    225.0,  0.0,    1.0},
      { 2,-1, 2, 0, 2,   -488.0,   0.0,    2.0,    207.0,  0.0,    1.0},
      { 0, 0, 0, 2, 2,   -468.0,   0.0,    0.0,    201.0,  0.0,    0.0},
      { 1,-1, 2, 0, 1,   -421.0,   0.0,    1.0,    216.0,  0.0,    1.0},
      {-1, 1, 2, 0, 2,    463.0,   0.0,    0.0,   -200.0,  0.0,    0.0},
      { 0, 1, 0, 2, 0,   -673.0,   0.0,    2.0,     14.0,  0.0,    0.0},
      { 0,-1,-2, 2, 0,    658.0,   0.0,    0.0,     -2.0,  0.0,    0.0},

   /* 121-130 */
      { 0, 3, 2,-2, 2,   -438.0,   0.0,    0.0,    188.0,  0.0,    0.0},
      { 0, 0, 0, 1, 1,   -390.0,   0.0,    0.0,    205.0,  0.0,    0.0},
      {-1, 0, 2, 2, 0,    639.0, -11.0,   -2.0,    -19.0,  0.0,    0.0},
      { 2, 1, 2, 0, 2,    412.0,   0.0,   -2.0,   -176.0,  0.0,   -1.0},
      { 1, 1, 0, 0, 1,   -361.0,   0.0,    0.0,    189.0,  0.0,    0.0},
      { 1, 1, 2, 0, 1,    360.0,   0.0,   -1.0,   -185.0,  0.0,   -1.0},
      { 2, 0, 0, 2, 0,    588.0,   0.0,   -3.0,    -24.0,  0.0,    0.0},
      { 1, 0,-2, 2, 0,   -578.0,   0.0,    1.0,      5.0,  0.0,    0.0},
      {-1, 0, 0, 2, 2,   -396.0,   0.0,    0.0,    171.0,  0.0,    0.0},
      { 0, 1, 0, 1, 0,    565.0,   0.0,   -1.0,     -6.0,  0.0,    0.0},

   /* 131-140 */
      { 0, 1, 0,-2, 1,   -335.0,   0.0,   -1.0,    184.0,  0.0,   -1.0},
      {-1, 0, 2,-2, 2,    357.0,   0.0,    1.0,   -154.0,  0.0,    0.0},
      { 0, 0, 0,-1, 1,    321.0,   0.0,    1.0,   -174.0,  0.0,    0.0},
      {-1, 1, 0, 0, 1,   -301.0,   0.0,   -1.0,    162.0,  0.0,    0.0},
      { 1, 0, 2,-1, 2,   -334.0,   0.0,    0.0,    144.0,  0.0,    0.0},
      { 1,-1, 0, 2, 0,    493.0,   0.0,   -2.0,    -15.0,  0.0,    0.0},
      { 0, 0, 0, 4, 0,    494.0,   0.0,   -2.0,    -19.0,  0.0,    0.0},
      { 1, 0, 2, 1, 2,    337.0,   0.0,   -1.0,   -143.0,  0.0,   -1.0},
      { 0, 0, 2, 1, 1,    280.0,   0.0,   -1.0,   -144.0,  0.0,    0.0},
      { 1, 0, 0,-2, 2,    309.0,   0.0,    1.0,   -134.0,  0.0,    0.0},

   /* 141-150 */
      {-1, 0, 2, 4, 1,   -263.0,   0.0,    2.0,    131.0,  0.0,    1.0},
      { 1, 0,-2, 0, 1,    253.0,   0.0,    1.0,   -138.0,  0.0,    0.0},
      { 1, 1, 2,-2, 1,    245.0,   0.0,    0.0,   -128.0,  0.0,    0.0},
      { 0, 0, 2, 2, 0,    416.0,   0.0,   -2.0,    -17.0,  0.0,    0.0},
      {-1, 0, 2,-1, 1,   -229.0,   0.0,    0.0,    128.0,  0.0,    0.0},
      {-2, 0, 2, 2, 1,    231.0,   0.0,    0.0,   -120.0,  0.0,    0.0},
      { 4, 0, 2, 0, 2,   -259.0,   0.0,    2.0,    109.0,  0.0,    1.0},
      { 2,-1, 0, 0, 0,    375.0,   0.0,   -1.0,     -8.0,  0.0,    0.0},
      { 2, 1, 2,-2, 2,    252.0,   0.0,    0.0,   -108.0,  0.0,    0.0},
      { 0, 1, 2, 1, 2,   -245.0,   0.0,    1.0,    104.0,  0.0,    0.0},

   /* 151-160 */
      { 1, 0, 4,-2, 2,    243.0,   0.0,   -1.0,   -104.0,  0.0,    0.0},
      {-1,-1, 0, 0, 1,    208.0,   0.0,    1.0,   -112.0,  0.0,    0.0},
      { 0, 1, 0, 2, 1,    199.0,   0.0,    0.0,   -102.0,  0.0,    0.0},
      {-2, 0, 2, 4, 1,   -208.0,   0.0,    1.0,    105.0,  0.0,    0.0},
      { 2, 0, 2, 0, 0,    335.0,   0.0,   -2.0,    -14.0,  0.0,    0.0},
      { 1, 0, 0, 1, 0,   -325.0,   0.0,    1.0,      7.0,  0.0,    0.0},
      {-1, 0, 0, 4, 1,   -187.0,   0.0,    0.0,     96.0,  0.0,    0.0},
      {-1, 0, 4, 0, 1,    197.0,   0.0,   -1.0,   -100.0,  0.0,    0.0},
      { 2, 0, 2, 2, 1,   -192.0,   0.0,    2.0,     94.0,  0.0,    1.0},
      { 0, 0, 2,-3, 2,   -188.0,   0.0,    0.0,     83.0,  0.0,    0.0},

   /* 161-170 */
      {-1,-2, 0, 2, 0,    276.0,   0.0,    0.0,     -2.0,  0.0,    0.0},
      { 2, 1, 0, 0, 0,   -286.0,   0.0,    1.0,      6.0,  0.0,    0.0},
      { 0, 0, 4, 0, 2,    186.0,   0.0,   -1.0,    -79.0,  0.0,    0.0},
      { 0, 0, 0, 0, 3,   -219.0,   0.0,    0.0,     43.0,  0.0,    0.0},
      { 0, 3, 0, 0, 0,    276.0,   0.0,    0.0,      2.0,  0.0,    0.0},
      { 0, 0, 2,-4, 1,   -153.0,   0.0,   -1.0,     84.0,  0.0,    0.0},
      { 0,-1, 0, 2, 1,   -156.0,   0.0,    0.0,     81.0,  0.0,    0.0},
      { 0, 0, 0, 4, 1,   -154.0,   0.0,    1.0,     78.0,  0.0,    0.0},
      {-1,-1, 2, 4, 2,   -174.0,   0.0,    1.0,     75.0,  0.0,    0.0},
      { 1, 0, 2, 4, 2,   -163.0,   0.0,    2.0,     69.0,  0.0,    1.0},

   /* 171-180 */
      {-2, 2, 0, 2, 0,   -228.0,   0.0,    0.0,      1.0,  0.0,    0.0},
      {-2,-1, 2, 0, 1,     91.0,   0.0,   -4.0,    -54.0,  0.0,   -2.0},
      {-2, 0, 0, 2, 2,    175.0,   0.0,    0.0,    -75.0,  0.0,    0.0},
      {-1,-1, 2, 0, 2,   -159.0,   0.0,    0.0,     69.0,  0.0,    0.0},
      { 0, 0, 4,-2, 1,    141.0,   0.0,    0.0,    -72.0,  0.0,    0.0},
      { 3, 0, 2,-2, 1,    147.0,   0.0,    0.0,    -75.0,  0.0,    0.0},
      {-2,-1, 0, 2, 1,   -132.0,   0.0,    0.0,     69.0,  0.0,    0.0},
      { 1, 0, 0,-1, 1,    159.0,   0.0,  -28.0,    -54.0,  0.0,   11.0},
      { 0,-2, 0, 2, 0,    213.0,   0.0,    0.0,     -4.0,  0.0,    0.0},
      {-2, 0, 0, 4, 1,    123.0,   0.0,    0.0,    -64.0,  0.0,    0.0},

   /* 181-190 */
      {-3, 0, 0, 0, 1,   -118.0,   0.0,   -1.0,     66.0,  0.0,    0.0},
      { 1, 1, 2, 2, 2,    144.0,   0.0,   -1.0,    -61.0,  0.0,    0.0},
      { 0, 0, 2, 4, 1,   -121.0,   0.0,    1.0,     60.0,  0.0,    0.0},
      { 3, 0, 2, 2, 2,   -134.0,   0.0,    1.0,     56.0,  0.0,    1.0},
      {-1, 1, 2,-2, 1,   -105.0,   0.0,    0.0,     57.0,  0.0,    0.0},
      { 2, 0, 0,-4, 1,   -102.0,   0.0,    0.0,     56.0,  0.0,    0.0},
      { 0, 0, 0,-2, 2,    120.0,   0.0,    0.0,    -52.0,  0.0,    0.0},
      { 2, 0, 2,-4, 1,    101.0,   0.0,    0.0,    -54.0,  0.0,    0.0},
      {-1, 1, 0, 2, 1,   -113.0,   0.0,    0.0,     59.0,  0.0,    0.0},
      { 0, 0, 2,-1, 1,   -106.0,   0.0,    0.0,     61.0,  0.0,    0.0},

   /* 191-200 */
      { 0,-2, 2, 2, 2,   -129.0,   0.0,    1.0,     55.0,  0.0,    0.0},
      { 2, 0, 0, 2, 1,   -114.0,   0.0,    0.0,     57.0,  0.0,    0.0},
      { 4, 0, 2,-2, 2,    113.0,   0.0,   -1.0,    -49.0,  0.0,    0.0},
      { 2, 0, 0,-2, 2,   -102.0,   0.0,    0.0,     44.0,  0.0,    0.0},
      { 0, 2, 0, 0, 1,    -94.0,   0.0,    0.0,     51.0,  0.0,    0.0},
      { 1, 0, 0,-4, 1,   -100.0,   0.0,   -1.0,     56.0,  0.0,    0.0},
      { 0, 2, 2,-2, 1,     87.0,   0.0,    0.0,    -47.0,  0.0,    0.0},
      {-3, 0, 0, 4, 0,    161.0,   0.0,    0.0,     -1.0,  0.0,    0.0},
      {-1, 1, 2, 0, 1,     96.0,   0.0,    0.0,    -50.0,  0.0,    0.0},
      {-1,-1, 0, 4, 0,    151.0,   0.0,   -1.0,     -5.0,  0.0,    0.0},

   /* 201-210 */
      {-1,-2, 2, 2, 2,   -104.0,   0.0,    0.0,     44.0,  0.0,    0.0},
      {-2,-1, 2, 4, 2,   -110.0,   0.0,    0.0,     48.0,  0.0,    0.0},
      { 1,-1, 2, 2, 1,   -100.0,   0.0,    1.0,     50.0,  0.0,    0.0},
      {-2, 1, 0, 2, 0,     92.0,   0.0,   -5.0,     12.0,  0.0,   -2.0},
      {-2, 1, 2, 0, 1,     82.0,   0.0,    0.0,    -45.0,  0.0,    0.0},
      { 2, 1, 0,-2, 1,     82.0,   0.0,    0.0,    -45.0,  0.0,    0.0},
      {-3, 0, 2, 0, 1,    -78.0,   0.0,    0.0,     41.0,  0.0,    0.0},
      {-2, 0, 2,-2, 1,    -77.0,   0.0,    0.0,     43.0,  0.0,    0.0},
      {-1, 1, 0, 2, 2,      2.0,   0.0,    0.0,     54.0,  0.0,    0.0},
      { 0,-1, 2,-1, 2,     94.0,   0.0,    0.0,    -40.0,  0.0,    0.0},

   /* 211-220 */
      {-1, 0, 4,-2, 2,    -93.0,   0.0,    0.0,     40.0,  0.0,    0.0},
      { 0,-2, 2, 0, 2,    -83.0,   0.0,   10.0,     40.0,  0.0,   -2.0},
      {-1, 0, 2, 1, 2,     83.0,   0.0,    0.0,    -36.0,  0.0,    0.0},
      { 2, 0, 0, 0, 2,    -91.0,   0.0,    0.0,     39.0,  0.0,    0.0},
      { 0, 0, 2, 0, 3,    128.0,   0.0,    0.0,     -1.0,  0.0,    0.0},
      {-2, 0, 4, 0, 2,    -79.0,   0.0,    0.0,     34.0,  0.0,    0.0},
      {-1, 0,-2, 0, 1,    -83.0,   0.0,    0.0,     47.0,  0.0,    0.0},
      {-1, 1, 2, 2, 1,     84.0,   0.0,    0.0,    -44.0,  0.0,    0.0},
      { 3, 0, 0, 0, 1,     83.0,   0.0,    0.0,    -43.0,  0.0,    0.0},
      {-1, 0, 2, 3, 2,     91.0,   0.0,    0.0,    -39.0,  0.0,    0.0},

   /* 221-230 */
      { 2,-1, 2, 0, 1,    -77.0,   0.0,    0.0,     39.0,  0.0,    0.0},
      { 0, 1, 2, 2, 1,     84.0,   0.0,    0.0,    -43.0,  0.0,    0.0},
      { 0,-1, 2, 4, 2,    -92.0,   0.0,    1.0,     39.0,  0.0,    0.0},
      { 2,-1, 2, 2, 2,    -92.0,   0.0,    1.0,     39.0,  0.0,    0.0},
      { 0, 2,-2, 2, 0,    -94.0,   0.0,    0.0,      0.0,  0.0,    0.0},
      {-1,-1, 2,-1, 1,     68.0,   0.0,    0.0,    -36.0,  0.0,    0.0},
      { 0,-2, 0, 0, 1,    -61.0,   0.0,    0.0,     32.0,  0.0,    0.0},
      { 1, 0, 2,-4, 2,     71.0,   0.0,    0.0,    -31.0,  0.0,    0.0},
      { 1,-1, 0,-2, 1,     62.0,   0.0,    0.0,    -34.0,  0.0,    0.0},
      {-1,-1, 2, 0, 1,    -63.0,   0.0,    0.0,     33.0,  0.0,    0.0},

   /* 231-240 */
      { 1,-1, 2,-2, 2,    -73.0,   0.0,    0.0,     32.0,  0.0,    0.0},
      {-2,-1, 0, 4, 0,    115.0,   0.0,    0.0,     -2.0,  0.0,    0.0},
      {-1, 0, 0, 3, 0,   -103.0,   0.0,    0.0,      2.0,  0.0,    0.0},
      {-2,-1, 2, 2, 2,     63.0,   0.0,    0.0,    -28.0,  0.0,    0.0},
      { 0, 2, 2, 0, 2,     74.0,   0.0,    0.0,    -32.0,  0.0,    0.0},
      { 1, 1, 0, 2, 0,   -103.0,   0.0,   -3.0,      3.0,  0.0,   -1.0},
      { 2, 0, 2,-1, 2,    -69.0,   0.0,    0.0,     30.0,  0.0,    0.0},
      { 1, 0, 2, 1, 1,     57.0,   0.0,    0.0,    -29.0,  0.0,    0.0},
      { 4, 0, 0, 0, 0,     94.0,   0.0,    0.0,     -4.0,  0.0,    0.0},
      { 2, 1, 2, 0, 1,     64.0,   0.0,    0.0,    -33.0,  0.0,    0.0},

   /* 241-250 */
      { 3,-1, 2, 0, 2,    -63.0,   0.0,    0.0,     26.0,  0.0,    0.0},
      {-2, 2, 0, 2, 1,    -38.0,   0.0,    0.0,     20.0,  0.0,    0.0},
      { 1, 0, 2,-3, 1,    -43.0,   0.0,    0.0,     24.0,  0.0,    0.0},
      { 1, 1, 2,-4, 1,    -45.0,   0.0,    0.0,     23.0,  0.0,    0.0},
      {-1,-1, 2,-2, 1,     47.0,   0.0,    0.0,    -24.0,  0.0,    0.0},
      { 0,-1, 0,-1, 1,    -48.0,   0.0,    0.0,     25.0,  0.0,    0.0},
      { 0,-1, 0,-2, 1,     45.0,   0.0,    0.0,    -26.0,  0.0,    0.0},
      {-2, 0, 0, 0, 2,     56.0,   0.0,    0.0,    -25.0,  0.0,    0.0},
      {-2, 0,-2, 2, 0,     88.0,   0.0,    0.0,      2.0,  0.0,    0.0},
      {-1, 0,-2, 4, 0,    -75.0,   0.0,    0.0,      0.0,  0.0,    0.0},

   /* 251-260 */
      { 1,-2, 0, 0, 0,     85.0,   0.0,    0.0,      0.0,  0.0,    0.0},
      { 0, 1, 0, 1, 1,     49.0,   0.0,    0.0,    -26.0,  0.0,    0.0},
      {-1, 2, 0, 2, 0,    -74.0,   0.0,   -3.0,     -1.0,  0.0,   -1.0},
      { 1,-1, 2,-2, 1,    -39.0,   0.0,    0.0,     21.0,  0.0,    0.0},
      { 1, 2, 2,-2, 2,     45.0,   0.0,    0.0,    -20.0,  0.0,    0.0},
      { 2,-1, 2,-2, 2,     51.0,   0.0,    0.0,    -22.0,  0.0,    0.0},
      { 1, 0, 2,-1, 1,    -40.0,   0.0,    0.0,     21.0,  0.0,    0.0},
      { 2, 1, 2,-2, 1,     41.0,   0.0,    0.0,    -21.0,  0.0,    0.0},
      {-2, 0, 0,-2, 1,    -42.0,   0.0,    0.0,     24.0,  0.0,    0.0},
      { 1,-2, 2, 0, 2,    -51.0,   0.0,    0.0,     22.0,  0.0,    0.0},

   /* 261-270 */
      { 0, 1, 2, 1, 1,    -42.0,   0.0,    0.0,     22.0,  0.0,    0.0},
      { 1, 0, 4,-2, 1,     39.0,   0.0,    0.0,    -21.0,  0.0,    0.0},
      {-2, 0, 4, 2, 2,     46.0,   0.0,    0.0,    -18.0,  0.0,    0.0},
      { 1, 1, 2, 1, 2,    -53.0,   0.0,    0.0,     22.0,  0.0,    0.0},
      { 1, 0, 0, 4, 0,     82.0,   0.0,    0.0,     -4.0,  0.0,    0.0},
      { 1, 0, 2, 2, 0,     81.0,   0.0,   -1.0,     -4.0,  0.0,    0.0},
      { 2, 0, 2, 1, 2,     47.0,   0.0,    0.0,    -19.0,  0.0,    0.0},
      { 3, 1, 2, 0, 2,     53.0,   0.0,    0.0,    -23.0,  0.0,    0.0},
      { 4, 0, 2, 0, 1,    -45.0,   0.0,    0.0,     22.0,  0.0,    0.0},
      {-2,-1, 2, 0, 0,    -44.0,   0.0,    0.0,     -2.0,  0.0,    0.0},

   /* 271-280 */
      { 0, 1,-2, 2, 1,    -33.0,   0.0,    0.0,     16.0,  0.0,    0.0},
      { 1, 0,-2, 1, 0,    -61.0,   0.0,    0.0,      1.0,  0.0,    0.0},
      { 0,-1,-2, 2, 1,     28.0,   0.0,    0.0,    -15.0,  0.0,    0.0},
      { 2,-1, 0,-2, 1,    -38.0,   0.0,    0.0,     19.0,  0.0,    0.0},
      {-1, 0, 2,-1, 2,    -33.0,   0.0,    0.0,     21.0,  0.0,    0.0},
      { 1, 0, 2,-3, 2,    -60.0,   0.0,    0.0,      0.0,  0.0,    0.0},
      { 0, 1, 2,-2, 3,     48.0,   0.0,    0.0,    -10.0,  0.0,    0.0},
      { 0, 0, 2,-3, 1,     27.0,   0.0,    0.0,    -14.0,  0.0,    0.0},
      {-1, 0,-2, 2, 1,     38.0,   0.0,    0.0,    -20.0,  0.0,    0.0},
      { 0, 0, 2,-4, 2,     31.0,   0.0,    0.0,    -13.0,  0.0,    0.0},

   /* 281-290 */
      {-2, 1, 0, 0, 1,    -29.0,   0.0,    0.0,     15.0,  0.0,    0.0},
      {-1, 0, 0,-1, 1,     28.0,   0.0,    0.0,    -15.0,  0.0,    0.0},
      { 2, 0, 2,-4, 2,    -32.0,   0.0,    0.0,     15.0,  0.0,    0.0},
      { 0, 0, 4,-4, 4,     45.0,   0.0,    0.0,     -8.0,  0.0,    0.0},
      { 0, 0, 4,-4, 2,    -44.0,   0.0,    0.0,     19.0,  0.0,    0.0},
      {-1,-2, 0, 2, 1,     28.0,   0.0,    0.0,    -15.0,  0.0,    0.0},
      {-2, 0, 0, 3, 0,    -51.0,   0.0,    0.0,      0.0,  0.0,    0.0},
      { 1, 0,-2, 2, 1,    -36.0,   0.0,    0.0,     20.0,  0.0,    0.0},
      {-3, 0, 2, 2, 2,     44.0,   0.0,    0.0,    -19.0,  0.0,    0.0},
      {-3, 0, 2, 2, 1,     26.0,   0.0,    0.0,    -14.0,  0.0,    0.0},

   /* 291-300 */
      {-2, 0, 2, 2, 0,    -60.0,   0.0,    0.0,      2.0,  0.0,    0.0},
      { 2,-1, 0, 0, 1,     35.0,   0.0,    0.0,    -18.0,  0.0,    0.0},
      {-2, 1, 2, 2, 2,    -27.0,   0.0,    0.0,     11.0,  0.0,    0.0},
      { 1, 1, 0, 1, 0,     47.0,   0.0,    0.0,     -1.0,  0.0,    0.0},
      { 0, 1, 4,-2, 2,     36.0,   0.0,    0.0,    -15.0,  0.0,    0.0},
      {-1, 1, 0,-2, 1,    -36.0,   0.0,    0.0,     20.0,  0.0,    0.0},
      { 0, 0, 0,-4, 1,    -35.0,   0.0,    0.0,     19.0,  0.0,    0.0},
      { 1,-1, 0, 2, 1,    -37.0,   0.0,    0.0,     19.0,  0.0,    0.0},
      { 1, 1, 0, 2, 1,     32.0,   0.0,    0.0,    -16.0,  0.0,    0.0},
      {-1, 2, 2, 2, 2,     35.0,   0.0,    0.0,    -14.0,  0.0,    0.0},

   /* 301-310 */
      { 3, 1, 2,-2, 2,     32.0,   0.0,    0.0,    -13.0,  0.0,    0.0},
      { 0,-1, 0, 4, 0,     65.0,   0.0,    0.0,     -2.0,  0.0,    0.0},
      { 2,-1, 0, 2, 0,     47.0,   0.0,    0.0,     -1.0,  0.0,    0.0},
      { 0, 0, 4, 0, 1,     32.0,   0.0,    0.0,    -16.0,  0.0,    0.0},
      { 2, 0, 4,-2, 2,     37.0,   0.0,    0.0,    -16.0,  0.0,    0.0},
      {-1,-1, 2, 4, 1,    -30.0,   0.0,    0.0,     15.0,  0.0,    0.0},
      { 1, 0, 0, 4, 1,    -32.0,   0.0,    0.0,     16.0,  0.0,    0.0},
      { 1,-2, 2, 2, 2,    -31.0,   0.0,    0.0,     13.0,  0.0,    0.0},
      { 0, 0, 2, 3, 2,     37.0,   0.0,    0.0,    -16.0,  0.0,    0.0},
      {-1, 1, 2, 4, 2,     31.0,   0.0,    0.0,    -13.0,  0.0,    0.0},

   /* 311-320 */
      { 3, 0, 0, 2, 0,     49.0,   0.0,    0.0,     -2.0,  0.0,    0.0},
      {-1, 0, 4, 2, 2,     32.0,   0.0,    0.0,    -13.0,  0.0,    0.0},
      { 1, 1, 2, 2, 1,     23.0,   0.0,    0.0,    -12.0,  0.0,    0.0},
      {-2, 0, 2, 6, 2,    -43.0,   0.0,    0.0,     18.0,  0.0,    0.0},
      { 2, 1, 2, 2, 2,     26.0,   0.0,    0.0,    -11.0,  0.0,    0.0},
      {-1, 0, 2, 6, 2,    -32.0,   0.0,    0.0,     14.0,  0.0,    0.0},
      { 1, 0, 2, 4, 1,    -29.0,   0.0,    0.0,     14.0,  0.0,    0.0},
      { 2, 0, 2, 4, 2,    -27.0,   0.0,    0.0,     12.0,  0.0,    0.0},
      { 1, 1,-2, 1, 0,     30.0,   0.0,    0.0,      0.0,  0.0,    0.0},
      {-3, 1, 2, 1, 2,    -11.0,   0.0,    0.0,      5.0,  0.0,    0.0},

   /* 321-330 */
      { 2, 0,-2, 0, 2,    -21.0,   0.0,    0.0,     10.0,  0.0,    0.0},
      {-1, 0, 0, 1, 2,    -34.0,   0.0,    0.0,     15.0,  0.0,    0.0},
      {-4, 0, 2, 2, 1,    -10.0,   0.0,    0.0,      6.0,  0.0,    0.0},
      {-1,-1, 0, 1, 0,    -36.0,   0.0,    0.0,      0.0,  0.0,    0.0},
      { 0, 0,-2, 2, 2,     -9.0,   0.0,    0.0,      4.0,  0.0,    0.0},
      { 1, 0, 0,-1, 2,    -12.0,   0.0,    0.0,      5.0,  0.0,    0.0},
      { 0,-1, 2,-2, 3,    -21.0,   0.0,    0.0,      5.0,  0.0,    0.0},
      {-2, 1, 2, 0, 0,    -29.0,   0.0,    0.0,     -1.0,  0.0,    0.0},
      { 0, 0, 2,-2, 4,    -15.0,   0.0,    0.0,      3.0,  0.0,    0.0},
      {-2,-2, 0, 2, 0,    -20.0,   0.0,    0.0,      0.0,  0.0,    0.0},

   /* 331-340 */
      {-2, 0,-2, 4, 0,     28.0,   0.0,    0.0,      0.0,  0.0,   -2.0},
      { 0,-2,-2, 2, 0,     17.0,   0.0,    0.0,      0.0,  0.0,    0.0},
      { 1, 2, 0,-2, 1,    -22.0,   0.0,    0.0,     12.0,  0.0,    0.0},
      { 3, 0, 0,-4, 1,    -14.0,   0.0,    0.0,      7.0,  0.0,    0.0},
      {-1, 1, 2,-2, 2,     24.0,   0.0,    0.0,    -11.0,  0.0,    0.0},
      { 1,-1, 2,-4, 1,     11.0,   0.0,    0.0,     -6.0,  0.0,    0.0},
      { 1, 1, 0,-2, 2,     14.0,   0.0,    0.0,     -6.0,  0.0,    0.0},
      {-3, 0, 2, 0, 0,     24.0,   0.0,    0.0,      0.0,  0.0,    0.0},
      {-3, 0, 2, 0, 2,     18.0,   0.0,    0.0,     -8.0,  0.0,    0.0},
      {-2, 0, 0, 1, 0,    -38.0,   0.0,    0.0,      0.0,  0.0,    0.0},

   /* 341-350 */
      { 0, 0,-2, 1, 0,    -31.0,   0.0,    0.0,      0.0,  0.0,    0.0},
      {-3, 0, 0, 2, 1,    -16.0,   0.0,    0.0,      8.0,  0.0,    0.0},
      {-1,-1,-2, 2, 0,     29.0,   0.0,    0.0,      0.0,  0.0,    0.0},
      { 0, 1, 2,-4, 1,    -18.0,   0.0,    0.0,     10.0,  0.0,    0.0},
      { 2, 1, 0,-4, 1,    -10.0,   0.0,    0.0,      5.0,  0.0,    0.0},
      { 0, 2, 0,-2, 1,    -17.0,   0.0,    0.0,     10.0,  0.0,    0.0},
      { 1, 0, 0,-3, 1,      9.0,   0.0,    0.0,     -4.0,  0.0,    0.0},
      {-2, 0, 2,-2, 2,     16.0,   0.0,    0.0,     -6.0,  0.0,    0.0},
      {-2,-1, 0, 0, 1,     22.0,   0.0,    0.0,    -12.0,  0.0,    0.0},
      {-4, 0, 0, 2, 0,     20.0,   0.0,    0.0,      0.0,  0.0,    0.0},

   /* 351-360 */
      { 1, 1, 0,-4, 1,    -13.0,   0.0,    0.0,      6.0,  0.0,    0.0},
      {-1, 0, 2,-4, 1,    -17.0,   0.0,    0.0,      9.0,  0.0,    0.0},
      { 0, 0, 4,-4, 1,    -14.0,   0.0,    0.0,      8.0,  0.0,    0.0},
      { 0, 3, 2,-2, 2,      0.0,   0.0,    0.0,     -7.0,  0.0,    0.0},
      {-3,-1, 0, 4, 0,     14.0,   0.0,    0.0,      0.0,  0.0,    0.0},
      {-3, 0, 0, 4, 1,     19.0,   0.0,    0.0,    -10.0,  0.0,    0.0},
      { 1,-1,-2, 2, 0,    -34.0,   0.0,    0.0,      0.0,  0.0,    0.0},
      {-1,-1, 0, 2, 2,    -20.0,   0.0,    0.0,      8.0,  0.0,    0.0},
      { 1,-2, 0, 0, 1,      9.0,   0.0,    0.0,     -5.0,  0.0,    0.0},
      { 1,-1, 0, 0, 2,    -18.0,   0.0,    0.0,      7.0,  0.0,    0.0},

   /* 361-370 */
      { 0, 0, 0, 1, 2,     13.0,   0.0,    0.0,     -6.0,  0.0,    0.0},
      {-1,-1, 2, 0, 0,     17.0,   0.0,    0.0,      0.0,  0.0,    0.0},
      { 1,-2, 2,-2, 2,    -12.0,   0.0,    0.0,      5.0,  0.0,    0.0},
      { 0,-1, 2,-1, 1,     15.0,   0.0,    0.0,     -8.0,  0.0,    0.0},
      {-1, 0, 2, 0, 3,    -11.0,   0.0,    0.0,      3.0,  0.0,    0.0},
      { 1, 1, 0, 0, 2,     13.0,   0.0,    0.0,     -5.0,  0.0,    0.0},
      {-1, 1, 2, 0, 0,    -18.0,   0.0,    0.0,      0.0,  0.0,    0.0},
      { 1, 2, 0, 0, 0,    -35.0,   0.0,    0.0,      0.0,  0.0,    0.0},
      {-1, 2, 2, 0, 2,      9.0,   0.0,    0.0,     -4.0,  0.0,    0.0},
      {-1, 0, 4,-2, 1,    -19.0,   0.0,    0.0,     10.0,  0.0,    0.0},

   /* 371-380 */
      { 3, 0, 2,-4, 2,    -26.0,   0.0,    0.0,     11.0,  0.0,    0.0},
      { 1, 2, 2,-2, 1,      8.0,   0.0,    0.0,     -4.0,  0.0,    0.0},
      { 1, 0, 4,-4, 2,    -10.0,   0.0,    0.0,      4.0,  0.0,    0.0},
      {-2,-1, 0, 4, 1,     10.0,   0.0,    0.0,     -6.0,  0.0,    0.0},
      { 0,-1, 0, 2, 2,    -21.0,   0.0,    0.0,      9.0,  0.0,    0.0},
      {-2, 1, 0, 4, 0,    -15.0,   0.0,    0.0,      0.0,  0.0,    0.0},
      {-2,-1, 2, 2, 1,      9.0,   0.0,    0.0,     -5.0,  0.0,    0.0},
      { 2, 0,-2, 2, 0,    -29.0,   0.0,    0.0,      0.0,  0.0,    0.0},
      { 1, 0, 0, 1, 1,    -19.0,   0.0,    0.0,     10.0,  0.0,    0.0},
      { 0, 1, 0, 2, 2,     12.0,   0.0,    0.0,     -5.0,  0.0,    0.0},

   /* 381-390 */
      { 1,-1, 2,-1, 2,     22.0,   0.0,    0.0,     -9.0,  0.0,    0.0},
      {-2, 0, 4, 0, 1,    -10.0,   0.0,    0.0,      5.0,  0.0,    0.0},
      { 2, 1, 0, 0, 1,    -20.0,   0.0,    0.0,     11.0,  0.0,    0.0},
      { 0, 1, 2, 0, 0,    -20.0,   0.0,    0.0,      0.0,  0.0,    0.0},
      { 0,-1, 4,-2, 2,    -17.0,   0.0,    0.0,      7.0,  0.0,    0.0},
      { 0, 0, 4,-2, 4,     15.0,   0.0,    0.0,     -3.0,  0.0,    0.0},
      { 0, 2, 2, 0, 1,      8.0,   0.0,    0.0,     -4.0,  0.0,    0.0},
      {-3, 0, 0, 6, 0,     14.0,   0.0,    0.0,      0.0,  0.0,    0.0},
      {-1,-1, 0, 4, 1,    -12.0,   0.0,    0.0,      6.0,  0.0,    0.0},
      { 1,-2, 0, 2, 0,     25.0,   0.0,    0.0,      0.0,  0.0,    0.0},

   /* 391-400 */
      {-1, 0, 0, 4, 2,    -13.0,   0.0,    0.0,      6.0,  0.0,    0.0},
      {-1,-2, 2, 2, 1,    -14.0,   0.0,    0.0,      8.0,  0.0,    0.0},
      {-1, 0, 0,-2, 2,     13.0,   0.0,    0.0,     -5.0,  0.0,    0.0},
      { 1, 0,-2,-2, 1,    -17.0,   0.0,    0.0,      9.0,  0.0,    0.0},
      { 0, 0,-2,-2, 1,    -12.0,   0.0,    0.0,      6.0,  0.0,    0.0},
      {-2, 0,-2, 0, 1,    -10.0,   0.0,    0.0,      5.0,  0.0,    0.0},
      { 0, 0, 0, 3, 1,     10.0,   0.0,    0.0,     -6.0,  0.0,    0.0},
      { 0, 0, 0, 3, 0,    -15.0,   0.0,    0.0,      0.0,  0.0,    0.0},
      {-1, 1, 0, 4, 0,    -22.0,   0.0,    0.0,      0.0,  0.0,    0.0},
      {-1,-1, 2, 2, 0,     28.0,   0.0,    0.0,     -1.0,  0.0,    0.0},

   /* 401-410 */
      {-2, 0, 2, 3, 2,     15.0,   0.0,    0.0,     -7.0,  0.0,    0.0},
      { 1, 0, 0, 2, 2,     23.0,   0.0,    0.0,    -10.0,  0.0,    0.0},
      { 0,-1, 2, 1, 2,     12.0,   0.0,    0.0,     -5.0,  0.0,    0.0},
      { 3,-1, 0, 0, 0,     29.0,   0.0,    0.0,     -1.0,  0.0,    0.0},
      { 2, 0, 0, 1, 0,    -25.0,   0.0,    0.0,      1.0,  0.0,    0.0},
      { 1,-1, 2, 0, 0,     22.0,   0.0,    0.0,      0.0,  0.0,    0.0},
      { 0, 0, 2, 1, 0,    -18.0,   0.0,    0.0,      0.0,  0.0,    0.0},
      { 1, 0, 2, 0, 3,     15.0,   0.0,    0.0,      3.0,  0.0,    0.0},
      { 3, 1, 0, 0, 0,    -23.0,   0.0,    0.0,      0.0,  0.0,    0.0},
      { 3,-1, 2,-2, 2,     12.0,   0.0,    0.0,     -5.0,  0.0,    0.0},

   /* 411-420 */
      { 2, 0, 2,-1, 1,     -8.0,   0.0,    0.0,      4.0,  0.0,    0.0},
      { 1, 1, 2, 0, 0,    -19.0,   0.0,    0.0,      0.0,  0.0,    0.0},
      { 0, 0, 4,-1, 2,    -10.0,   0.0,    0.0,      4.0,  0.0,    0.0},
      { 1, 2, 2, 0, 2,     21.0,   0.0,    0.0,     -9.0,  0.0,    0.0},
      {-2, 0, 0, 6, 0,     23.0,   0.0,    0.0,     -1.0,  0.0,    0.0},
      { 0,-1, 0, 4, 1,    -16.0,   0.0,    0.0,      8.0,  0.0,    0.0},
      {-2,-1, 2, 4, 1,    -19.0,   0.0,    0.0,      9.0,  0.0,    0.0},
      { 0,-2, 2, 2, 1,    -22.0,   0.0,    0.0,     10.0,  0.0,    0.0},
      { 0,-1, 2, 2, 0,     27.0,   0.0,    0.0,     -1.0,  0.0,    0.0},
      {-1, 0, 2, 3, 1,     16.0,   0.0,    0.0,     -8.0,  0.0,    0.0},

   /* 421-430 */
      {-2, 1, 2, 4, 2,     19.0,   0.0,    0.0,     -8.0,  0.0,    0.0},
      { 2, 0, 0, 2, 2,      9.0,   0.0,    0.0,     -4.0,  0.0,    0.0},
      { 2,-2, 2, 0, 2,     -9.0,   0.0,    0.0,      4.0,  0.0,    0.0},
      {-1, 1, 2, 3, 2,     -9.0,   0.0,    0.0,      4.0,  0.0,    0.0},
      { 3, 0, 2,-1, 2,     -8.0,   0.0,    0.0,      4.0,  0.0,    0.0},
      { 4, 0, 2,-2, 1,     18.0,   0.0,    0.0,     -9.0,  0.0,    0.0},
      {-1, 0, 0, 6, 0,     16.0,   0.0,    0.0,     -1.0,  0.0,    0.0},
      {-1,-2, 2, 4, 2,    -10.0,   0.0,    0.0,      4.0,  0.0,    0.0},
      {-3, 0, 2, 6, 2,    -23.0,   0.0,    0.0,      9.0,  0.0,    0.0},
      {-1, 0, 2, 4, 0,     16.0,   0.0,    0.0,     -1.0,  0.0,    0.0},

   /* 431-440 */
      { 3, 0, 0, 2, 1,    -12.0,   0.0,    0.0,      6.0,  0.0,    0.0},
      { 3,-1, 2, 0, 1,     -8.0,   0.0,    0.0,      4.0,  0.0,    0.0},
      { 3, 0, 2, 0, 0,     30.0,   0.0,    0.0,     -2.0,  0.0,    0.0},
      { 1, 0, 4, 0, 2,     24.0,   0.0,    0.0,    -10.0,  0.0,    0.0},
      { 5, 0, 2,-2, 2,     10.0,   0.0,    0.0,     -4.0,  0.0,    0.0},
      { 0,-1, 2, 4, 1,    -16.0,   0.0,    0.0,      7.0,  0.0,    0.0},
      { 2,-1, 2, 2, 1,    -16.0,   0.0,    0.0,      7.0,  0.0,    0.0},
      { 0, 1, 2, 4, 2,     17.0,   0.0,    0.0,     -7.0,  0.0,    0.0},
      { 1,-1, 2, 4, 2,    -24.0,   0.0,    0.0,     10.0,  0.0,    0.0},
      { 3,-1, 2, 2, 2,    -12.0,   0.0,    0.0,      5.0,  0.0,    0.0},

   /* 441-450 */
      { 3, 0, 2, 2, 1,    -24.0,   0.0,    0.0,     11.0,  0.0,    0.0},
      { 5, 0, 2, 0, 2,    -23.0,   0.0,    0.0,      9.0,  0.0,    0.0},
      { 0, 0, 2, 6, 2,    -13.0,   0.0,    0.0,      5.0,  0.0,    0.0},
      { 4, 0, 2, 2, 2,    -15.0,   0.0,    0.0,      7.0,  0.0,    0.0},
      { 0,-1, 1,-1, 1,      0.0,   0.0,-1988.0,      0.0,  0.0,-1679.0},
      {-1, 0, 1, 0, 3,      0.0,   0.0,  -63.0,      0.0,  0.0,  -27.0},
      { 0,-2, 2,-2, 3,     -4.0,   0.0,    0.0,      0.0,  0.0,    0.0},
      { 1, 0,-1, 0, 1,      0.0,   0.0,    5.0,      0.0,  0.0,    4.0},
      { 2,-2, 0,-2, 1,      5.0,   0.0,    0.0,     -3.0,  0.0,    0.0},
      {-1, 0, 1, 0, 2,      0.0,   0.0,  364.0,      0.0,  0.0,  176.0},

   /* 451-460 */
      {-1, 0, 1, 0, 1,      0.0,   0.0,-1044.0,      0.0,  0.0, -891.0},
      {-1,-1, 2,-1, 2,     -3.0,   0.0,    0.0,      1.0,  0.0,    0.0},
      {-2, 2, 0, 2, 2,      4.0,   0.0,    0.0,     -2.0,  0.0,    0.0},
      {-1, 0, 1, 0, 0,      0.0,   0.0,  330.0,      0.0,  0.0,    0.0},
      {-4, 1, 2, 2, 2,      5.0,   0.0,    0.0,     -2.0,  0.0,    0.0},
      {-3, 0, 2, 1, 1,      3.0,   0.0,    0.0,     -2.0,  0.0,    0.0},
      {-2,-1, 2, 0, 2,     -3.0,   0.0,    0.0,      1.0,  0.0,    0.0},
      { 1, 0,-2, 1, 1,     -5.0,   0.0,    0.0,      2.0,  0.0,    0.0},
      { 2,-1,-2, 0, 1,      3.0,   0.0,    0.0,     -1.0,  0.0,    0.0},
      {-4, 0, 2, 2, 0,      3.0,   0.0,    0.0,      0.0,  0.0,    0.0},

   /* 461-470 */
      {-3, 1, 0, 3, 0,      3.0,   0.0,    0.0,      0.0,  0.0,    0.0},
      {-1, 0,-1, 2, 0,      0.0,   0.0,    5.0,      0.0,  0.0,    0.0},
      { 0,-2, 0, 0, 2,      0.0,   0.0,    0.0,      1.0,  0.0,    0.0},
      { 0,-2, 0, 0, 2,      4.0,   0.0,    0.0,     -2.0,  0.0,    0.0},
      {-3, 0, 0, 3, 0,      6.0,   0.0,    0.0,      0.0,  0.0,    0.0},
      {-2,-1, 0, 2, 2,      5.0,   0.0,    0.0,     -2.0,  0.0,    0.0},
      {-1, 0,-2, 3, 0,     -7.0,   0.0,    0.0,      0.0,  0.0,    0.0},
      {-4, 0, 0, 4, 0,    -12.0,   0.0,    0.0,      0.0,  0.0,    0.0},
      { 2, 1,-2, 0, 1,      5.0,   0.0,    0.0,     -3.0,  0.0,    0.0},
      { 2,-1, 0,-2, 2,      3.0,   0.0,    0.0,     -1.0,  0.0,    0.0},

   /* 471-480 */
      { 0, 0, 1,-1, 0,     -5.0,   0.0,    0.0,      0.0,  0.0,    0.0},
      {-1, 2, 0, 1, 0,      3.0,   0.0,    0.0,      0.0,  0.0,    0.0},
      {-2, 1, 2, 0, 2,     -7.0,   0.0,    0.0,      3.0,  0.0,    0.0},
      { 1, 1, 0,-1, 1,      7.0,   0.0,    0.0,     -4.0,  0.0,    0.0},
      { 1, 0, 1,-2, 1,      0.0,   0.0,  -12.0,      0.0,  0.0,  -10.0},
      { 0, 2, 0, 0, 2,      4.0,   0.0,    0.0,     -2.0,  0.0,    0.0},
      { 1,-1, 2,-3, 1,      3.0,   0.0,    0.0,     -2.0,  0.0,    0.0},
      {-1, 1, 2,-1, 1,     -3.0,   0.0,    0.0,      2.0,  0.0,    0.0},
      {-2, 0, 4,-2, 2,     -7.0,   0.0,    0.0,      3.0,  0.0,    0.0},
      {-2, 0, 4,-2, 1,     -4.0,   0.0,    0.0,      2.0,  0.0,    0.0},

   /* 481-490 */
      {-2,-2, 0, 2, 1,     -3.0,   0.0,    0.0,      1.0,  0.0,    0.0},
      {-2, 0,-2, 4, 0,      0.0,   0.0,    0.0,      0.0,  0.0,    0.0},
      { 1, 2, 2,-4, 1,     -3.0,   0.0,    0.0,      1.0,  0.0,    0.0},
      { 1, 1, 2,-4, 2,      7.0,   0.0,    0.0,     -3.0,  0.0,    0.0},
      {-1, 2, 2,-2, 1,     -4.0,   0.0,    0.0,      2.0,  0.0,    0.0},
      { 2, 0, 0,-3, 1,      4.0,   0.0,    0.0,     -2.0,  0.0,    0.0},
      {-1, 2, 0, 0, 1,     -5.0,   0.0,    0.0,      3.0,  0.0,    0.0},
      { 0, 0, 0,-2, 0,      5.0,   0.0,    0.0,      0.0,  0.0,    0.0},
      {-1,-1, 2,-2, 2,     -5.0,   0.0,    0.0,      2.0,  0.0,    0.0},
      {-1, 1, 0, 0, 2,      5.0,   0.0,    0.0,     -2.0,  0.0,    0.0},

   /* 491-500 */
      { 0, 0, 0,-1, 2,     -8.0,   0.0,    0.0,      3.0,  0.0,    0.0},
      {-2, 1, 0, 1, 0,      9.0,   0.0,    0.0,      0.0,  0.0,    0.0},
      { 1,-2, 0,-2, 1,      6.0,   0.0,    0.0,     -3.0,  0.0,    0.0},
      { 1, 0,-2, 0, 2,     -5.0,   0.0,    0.0,      2.0,  0.0,    0.0},
      {-3, 1, 0, 2, 0,      3.0,   0.0,    0.0,      0.0,  0.0,    0.0},
      {-1, 1,-2, 2, 0,     -7.0,   0.0,    0.0,      0.0,  0.0,    0.0},
      {-1,-1, 0, 0, 2,     -3.0,   0.0,    0.0,      1.0,  0.0,    0.0},
      {-3, 0, 0, 2, 0,      5.0,   0.0,    0.0,      0.0,  0.0,    0.0},
      {-3,-1, 0, 2, 0,      3.0,   0.0,    0.0,      0.0,  0.0,    0.0},
      { 2, 0, 2,-6, 1,     -3.0,   0.0,    0.0,      2.0,  0.0,    0.0},

   /* 501-510 */
      { 0, 1, 2,-4, 2,      4.0,   0.0,    0.0,     -2.0,  0.0,    0.0},
      { 2, 0, 0,-4, 2,      3.0,   0.0,    0.0,     -1.0,  0.0,    0.0},
      {-2, 1, 2,-2, 1,     -5.0,   0.0,    0.0,      2.0,  0.0,    0.0},
      { 0,-1, 2,-4, 1,      4.0,   0.0,    0.0,     -2.0,  0.0,    0.0},
      { 0, 1, 0,-2, 2,      9.0,   0.0,    0.0,     -3.0,  0.0,    0.0},
      {-1, 0, 0,-2, 0,      4.0,   0.0,    0.0,      0.0,  0.0,    0.0},
      { 2, 0,-2,-2, 1,      4.0,   0.0,    0.0,     -2.0,  0.0,    0.0},
      {-4, 0, 2, 0, 1,     -3.0,   0.0,    0.0,      2.0,  0.0,    0.0},
      {-1,-1, 0,-1, 1,     -4.0,   0.0,    0.0,      2.0,  0.0,    0.0},
      { 0, 0,-2, 0, 2,      9.0,   0.0,    0.0,     -3.0,  0.0,    0.0},

   /* 511-520 */
      {-3, 0, 0, 1, 0,     -4.0,   0.0,    0.0,      0.0,  0.0,    0.0},
      {-1, 0,-2, 1, 0,     -4.0,   0.0,    0.0,      0.0,  0.0,    0.0},
      {-2, 0,-2, 2, 1,      3.0,   0.0,    0.0,     -2.0,  0.0,    0.0},
      { 0, 0,-4, 2, 0,      8.0,   0.0,    0.0,      0.0,  0.0,    0.0},
      {-2,-1,-2, 2, 0,      3.0,   0.0,    0.0,      0.0,  0.0,    0.0},
      { 1, 0, 2,-6, 1,     -3.0,   0.0,    0.0,      2.0,  0.0,    0.0},
      {-1, 0, 2,-4, 2,      3.0,   0.0,    0.0,     -1.0,  0.0,    0.0},
      { 1, 0, 0,-4, 2,      3.0,   0.0,    0.0,     -1.0,  0.0,    0.0},
      { 2, 1, 2,-4, 2,     -3.0,   0.0,    0.0,      1.0,  0.0,    0.0},
      { 2, 1, 2,-4, 1,      6.0,   0.0,    0.0,     -3.0,  0.0,    0.0},

   /* 521-530 */
      { 0, 1, 4,-4, 4,      3.0,   0.0,    0.0,      0.0,  0.0,    0.0},
      { 0, 1, 4,-4, 2,     -3.0,   0.0,    0.0,      1.0,  0.0,    0.0},
      {-1,-1,-2, 4, 0,     -7.0,   0.0,    0.0,      0.0,  0.0,    0.0},
      {-1,-3, 0, 2, 0,      9.0,   0.0,    0.0,      0.0,  0.0,    0.0},
      {-1, 0,-2, 4, 1,     -3.0,   0.0,    0.0,      2.0,  0.0,    0.0},
      {-2,-1, 0, 3, 0,     -3.0,   0.0,    0.0,      0.0,  0.0,    0.0},
      { 0, 0,-2, 3, 0,     -4.0,   0.0,    0.0,      0.0,  0.0,    0.0},
      {-2, 0, 0, 3, 1,     -5.0,   0.0,    0.0,      3.0,  0.0,    0.0},
      { 0,-1, 0, 1, 0,    -13.0,   0.0,    0.0,      0.0,  0.0,    0.0},
      {-3, 0, 2, 2, 0,     -7.0,   0.0,    0.0,      0.0,  0.0,    0.0},

   /* 531-540 */
      { 1, 1,-2, 2, 0,     10.0,   0.0,    0.0,      0.0,  0.0,    0.0},
      {-1, 1, 0, 2, 2,      3.0,   0.0,    0.0,     -1.0,  0.0,    0.0},
      { 1,-2, 2,-2, 1,     10.0,   0.0,   13.0,      6.0,  0.0,   -5.0},
      { 0, 0, 1, 0, 2,      0.0,   0.0,   30.0,      0.0,  0.0,   14.0},
      { 0, 0, 1, 0, 1,      0.0,   0.0, -162.0,      0.0,  0.0, -138.0},
      { 0, 0, 1, 0, 0,      0.0,   0.0,   75.0,      0.0,  0.0,    0.0},
      {-1, 2, 0, 2, 1,     -7.0,   0.0,    0.0,      4.0,  0.0,    0.0},
      { 0, 0, 2, 0, 2,     -4.0,   0.0,    0.0,      2.0,  0.0,    0.0},
      {-2, 0, 2, 0, 2,      4.0,   0.0,    0.0,     -2.0,  0.0,    0.0},
      { 2, 0, 0,-1, 1,      5.0,   0.0,    0.0,     -2.0,  0.0,    0.0},

   /* 541-550 */
      { 3, 0, 0,-2, 1,      5.0,   0.0,    0.0,     -3.0,  0.0,    0.0},
      { 1, 0, 2,-2, 3,     -3.0,   0.0,    0.0,      0.0,  0.0,    0.0},
      { 1, 2, 0, 0, 1,     -3.0,   0.0,    0.0,      2.0,  0.0,    0.0},
      { 2, 0, 2,-3, 2,     -4.0,   0.0,    0.0,      2.0,  0.0,    0.0},
      {-1, 1, 4,-2, 2,     -5.0,   0.0,    0.0,      2.0,  0.0,    0.0},
      {-2,-2, 0, 4, 0,      6.0,   0.0,    0.0,      0.0,  0.0,    0.0},
      { 0,-3, 0, 2, 0,      9.0,   0.0,    0.0,      0.0,  0.0,    0.0},
      { 0, 0,-2, 4, 0,      5.0,   0.0,    0.0,      0.0,  0.0,    0.0},
      {-1,-1, 0, 3, 0,     -7.0,   0.0,    0.0,      0.0,  0.0,    0.0},
      {-2, 0, 0, 4, 2,     -3.0,   0.0,    0.0,      1.0,  0.0,    0.0},

   /* 551-560 */
      {-1, 0, 0, 3, 1,     -4.0,   0.0,    0.0,      2.0,  0.0,    0.0},
      { 2,-2, 0, 0, 0,      7.0,   0.0,    0.0,      0.0,  0.0,    0.0},
      { 1,-1, 0, 1, 0,     -4.0,   0.0,    0.0,      0.0,  0.0,    0.0},
      {-1, 0, 0, 2, 0,      4.0,   0.0,    0.0,      0.0,  0.0,    0.0},
      { 0,-2, 2, 0, 1,     -6.0,   0.0,   -3.0,      3.0,  0.0,    1.0},
      {-1, 0, 1, 2, 1,      0.0,   0.0,   -3.0,      0.0,  0.0,   -2.0},
      {-1, 1, 0, 3, 0,     11.0,   0.0,    0.0,      0.0,  0.0,    0.0},
      {-1,-1, 2, 1, 2,      3.0,   0.0,    0.0,     -1.0,  0.0,    0.0},
      { 0,-1, 2, 0, 0,     11.0,   0.0,    0.0,      0.0,  0.0,    0.0},
      {-2, 1, 2, 2, 1,     -3.0,   0.0,    0.0,      2.0,  0.0,    0.0},

   /* 561-570 */
      { 2,-2, 2,-2, 2,     -1.0,   0.0,    3.0,      3.0,  0.0,   -1.0},
      { 1, 1, 0, 1, 1,      4.0,   0.0,    0.0,     -2.0,  0.0,    0.0},
      { 1, 0, 1, 0, 1,      0.0,   0.0,  -13.0,      0.0,  0.0,  -11.0},
      { 1, 0, 1, 0, 0,      3.0,   0.0,    6.0,      0.0,  0.0,    0.0},
      { 0, 2, 0, 2, 0,     -7.0,   0.0,    0.0,      0.0,  0.0,    0.0},
      { 2,-1, 2,-2, 1,      5.0,   0.0,    0.0,     -3.0,  0.0,    0.0},
      { 0,-1, 4,-2, 1,     -3.0,   0.0,    0.0,      1.0,  0.0,    0.0},
      { 0, 0, 4,-2, 3,      3.0,   0.0,    0.0,      0.0,  0.0,    0.0},
      { 0, 1, 4,-2, 1,      5.0,   0.0,    0.0,     -3.0,  0.0,    0.0},
      { 4, 0, 2,-4, 2,     -7.0,   0.0,    0.0,      3.0,  0.0,    0.0},

   /* 571-580 */
      { 2, 2, 2,-2, 2,      8.0,   0.0,    0.0,     -3.0,  0.0,    0.0},
      { 2, 0, 4,-4, 2,     -4.0,   0.0,    0.0,      2.0,  0.0,    0.0},
      {-1,-2, 0, 4, 0,     11.0,   0.0,    0.0,      0.0,  0.0,    0.0},
      {-1,-3, 2, 2, 2,     -3.0,   0.0,    0.0,      1.0,  0.0,    0.0},
      {-3, 0, 2, 4, 2,      3.0,   0.0,    0.0,     -1.0,  0.0,    0.0},
      {-3, 0, 2,-2, 1,     -4.0,   0.0,    0.0,      2.0,  0.0,    0.0},
      {-1,-1, 0,-2, 1,      8.0,   0.0,    0.0,     -4.0,  0.0,    0.0},
      {-3, 0, 0, 0, 2,      3.0,   0.0,    0.0,     -1.0,  0.0,    0.0},
      {-3, 0,-2, 2, 0,     11.0,   0.0,    0.0,      0.0,  0.0,    0.0},
      { 0, 1, 0,-4, 1,     -6.0,   0.0,    0.0,      3.0,  0.0,    0.0},

   /* 581-590 */
      {-2, 1, 0,-2, 1,     -4.0,   0.0,    0.0,      2.0,  0.0,    0.0},
      {-4, 0, 0, 0, 1,     -8.0,   0.0,    0.0,      4.0,  0.0,    0.0},
      {-1, 0, 0,-4, 1,     -7.0,   0.0,    0.0,      3.0,  0.0,    0.0},
      {-3, 0, 0,-2, 1,     -4.0,   0.0,    0.0,      2.0,  0.0,    0.0},
      { 0, 0, 0, 3, 2,      3.0,   0.0,    0.0,     -1.0,  0.0,    0.0},
      {-1, 1, 0, 4, 1,      6.0,   0.0,    0.0,     -3.0,  0.0,    0.0},
      { 1,-2, 2, 0, 1,     -6.0,   0.0,    0.0,      3.0,  0.0,    0.0},
      { 0, 1, 0, 3, 0,      6.0,   0.0,    0.0,      0.0,  0.0,    0.0},
      {-1, 0, 2, 2, 3,      6.0,   0.0,    0.0,     -1.0,  0.0,    0.0},
      { 0, 0, 2, 2, 2,      5.0,   0.0,    0.0,     -2.0,  0.0,    0.0},

   /* 591-600 */
      {-2, 0, 2, 2, 2,     -5.0,   0.0,    0.0,      2.0,  0.0,    0.0},
      {-1, 1, 2, 2, 0,     -4.0,   0.0,    0.0,      0.0,  0.0,    0.0},
      { 3, 0, 0, 0, 2,     -4.0,   0.0,    0.0,      2.0,  0.0,    0.0},
      { 2, 1, 0, 1, 0,      4.0,   0.0,    0.0,      0.0,  0.0,    0.0},
      { 2,-1, 2,-1, 2,      6.0,   0.0,    0.0,     -3.0,  0.0,    0.0},
      { 0, 0, 2, 0, 1,     -4.0,   0.0,    0.0,      2.0,  0.0,    0.0},
      { 0, 0, 3, 0, 3,      0.0,   0.0,  -26.0,      0.0,  0.0,  -11.0},
      { 0, 0, 3, 0, 2,      0.0,   0.0,  -10.0,      0.0,  0.0,   -5.0},
      {-1, 2, 2, 2, 1,      5.0,   0.0,    0.0,     -3.0,  0.0,    0.0},
      {-1, 0, 4, 0, 0,    -13.0,   0.0,    0.0,      0.0,  0.0,    0.0},

   /* 601-610 */
      { 1, 2, 2, 0, 1,      3.0,   0.0,    0.0,     -2.0,  0.0,    0.0},
      { 3, 1, 2,-2, 1,      4.0,   0.0,    0.0,     -2.0,  0.0,    0.0},
      { 1, 1, 4,-2, 2,      7.0,   0.0,    0.0,     -3.0,  0.0,    0.0},
      {-2,-1, 0, 6, 0,      4.0,   0.0,    0.0,      0.0,  0.0,    0.0},
      { 0,-2, 0, 4, 0,      5.0,   0.0,    0.0,      0.0,  0.0,    0.0},
      {-2, 0, 0, 6, 1,     -3.0,   0.0,    0.0,      2.0,  0.0,    0.0},
      {-2,-2, 2, 4, 2,     -6.0,   0.0,    0.0,      2.0,  0.0,    0.0},
      { 0,-3, 2, 2, 2,     -5.0,   0.0,    0.0,      2.0,  0.0,    0.0},
      { 0, 0, 0, 4, 2,     -7.0,   0.0,    0.0,      3.0,  0.0,    0.0},
      {-1,-1, 2, 3, 2,      5.0,   0.0,    0.0,     -2.0,  0.0,    0.0},

   /* 611-620 */
      {-2, 0, 2, 4, 0,     13.0,   0.0,    0.0,      0.0,  0.0,    0.0},
      { 2,-1, 0, 2, 1,     -4.0,   0.0,    0.0,      2.0,  0.0,    0.0},
      { 1, 0, 0, 3, 0,     -3.0,   0.0,    0.0,      0.0,  0.0,    0.0},
      { 0, 1, 0, 4, 1,      5.0,   0.0,    0.0,     -2.0,  0.0,    0.0},
      { 0, 1, 0, 4, 0,    -11.0,   0.0,    0.0,      0.0,  0.0,    0.0},
      { 1,-1, 2, 1, 2,      5.0,   0.0,    0.0,     -2.0,  0.0,    0.0},
      { 0, 0, 2, 2, 3,      4.0,   0.0,    0.0,      0.0,  0.0,    0.0},
      { 1, 0, 2, 2, 2,      4.0,   0.0,    0.0,     -2.0,  0.0,    0.0},
      {-1, 0, 2, 2, 2,     -4.0,   0.0,    0.0,      2.0,  0.0,    0.0},
      {-2, 0, 4, 2, 1,      6.0,   0.0,    0.0,     -3.0,  0.0,    0.0},

   /* 621-630 */
      { 2, 1, 0, 2, 1,      3.0,   0.0,    0.0,     -2.0,  0.0,    0.0},
      { 2, 1, 0, 2, 0,    -12.0,   0.0,    0.0,      0.0,  0.0,    0.0},
      { 2,-1, 2, 0, 0,      4.0,   0.0,    0.0,      0.0,  0.0,    0.0},
      { 1, 0, 2, 1, 0,     -3.0,   0.0,    0.0,      0.0,  0.0,    0.0},
      { 0, 1, 2, 2, 0,     -4.0,   0.0,    0.0,      0.0,  0.0,    0.0},
      { 2, 0, 2, 0, 3,      3.0,   0.0,    0.0,      0.0,  0.0,    0.0},
      { 3, 0, 2, 0, 2,      3.0,   0.0,    0.0,     -1.0,  0.0,    0.0},
      { 1, 0, 2, 0, 2,     -3.0,   0.0,    0.0,      1.0,  0.0,    0.0},
      { 1, 0, 3, 0, 3,      0.0,   0.0,   -5.0,      0.0,  0.0,   -2.0},
      { 1, 1, 2, 1, 1,     -7.0,   0.0,    0.0,      4.0,  0.0,    0.0},

   /* 631-640 */
      { 0, 2, 2, 2, 2,      6.0,   0.0,    0.0,     -3.0,  0.0,    0.0},
      { 2, 1, 2, 0, 0,     -3.0,   0.0,    0.0,      0.0,  0.0,    0.0},
      { 2, 0, 4,-2, 1,      5.0,   0.0,    0.0,     -3.0,  0.0,    0.0},
      { 4, 1, 2,-2, 2,      3.0,   0.0,    0.0,     -1.0,  0.0,    0.0},
      {-1,-1, 0, 6, 0,      3.0,   0.0,    0.0,      0.0,  0.0,    0.0},
      {-3,-1, 2, 6, 2,     -3.0,   0.0,    0.0,      1.0,  0.0,    0.0},
      {-1, 0, 0, 6, 1,     -5.0,   0.0,    0.0,      3.0,  0.0,    0.0},
      {-3, 0, 2, 6, 1,     -3.0,   0.0,    0.0,      2.0,  0.0,    0.0},
      { 1,-1, 0, 4, 1,     -3.0,   0.0,    0.0,      2.0,  0.0,    0.0},
      { 1,-1, 0, 4, 0,     12.0,   0.0,    0.0,      0.0,  0.0,    0.0},

   /* 641-650 */
      {-2, 0, 2, 5, 2,      3.0,   0.0,    0.0,     -1.0,  0.0,    0.0},
      { 1,-2, 2, 2, 1,     -4.0,   0.0,    0.0,      2.0,  0.0,    0.0},
      { 3,-1, 0, 2, 0,      4.0,   0.0,    0.0,      0.0,  0.0,    0.0},
      { 1,-1, 2, 2, 0,      6.0,   0.0,    0.0,      0.0,  0.0,    0.0},
      { 0, 0, 2, 3, 1,      5.0,   0.0,    0.0,     -3.0,  0.0,    0.0},
      {-1, 1, 2, 4, 1,      4.0,   0.0,    0.0,     -2.0,  0.0,    0.0},
      { 0, 1, 2, 3, 2,     -6.0,   0.0,    0.0,      3.0,  0.0,    0.0},
      {-1, 0, 4, 2, 1,      4.0,   0.0,    0.0,     -2.0,  0.0,    0.0},
      { 2, 0, 2, 1, 1,      6.0,   0.0,    0.0,     -3.0,  0.0,    0.0},
      { 5, 0, 0, 0, 0,      6.0,   0.0,    0.0,      0.0,  0.0,    0.0},

   /* 651-660 */
      { 2, 1, 2, 1, 2,     -6.0,   0.0,    0.0,      3.0,  0.0,    0.0},
      { 1, 0, 4, 0, 1,      3.0,   0.0,    0.0,     -2.0,  0.0,    0.0},
      { 3, 1, 2, 0, 1,      7.0,   0.0,    0.0,     -4.0,  0.0,    0.0},
      { 3, 0, 4,-2, 2,      4.0,   0.0,    0.0,     -2.0,  0.0,    0.0},
      {-2,-1, 2, 6, 2,     -5.0,   0.0,    0.0,      2.0,  0.0,    0.0},
      { 0, 0, 0, 6, 0,      5.0,   0.0,    0.0,      0.0,  0.0,    0.0},
      { 0,-2, 2, 4, 2,     -6.0,   0.0,    0.0,      3.0,  0.0,    0.0},
      {-2, 0, 2, 6, 1,     -6.0,   0.0,    0.0,      3.0,  0.0,    0.0},
      { 2, 0, 0, 4, 1,     -4.0,   0.0,    0.0,      2.0,  0.0,    0.0},
      { 2, 0, 0, 4, 0,     10.0,   0.0,    0.0,      0.0,  0.0,    0.0},

   /* 661-670 */
      { 2,-2, 2, 2, 2,     -4.0,   0.0,    0.0,      2.0,  0.0,    0.0},
      { 0, 0, 2, 4, 0,      7.0,   0.0,    0.0,      0.0,  0.0,    0.0},
      { 1, 0, 2, 3, 2,      7.0,   0.0,    0.0,     -3.0,  0.0,    0.0},
      { 4, 0, 0, 2, 0,      4.0,   0.0,    0.0,      0.0,  0.0,    0.0},
      { 2, 0, 2, 2, 0,     11.0,   0.0,    0.0,      0.0,  0.0,    0.0},
      { 0, 0, 4, 2, 2,      5.0,   0.0,    0.0,     -2.0,  0.0,    0.0},
      { 4,-1, 2, 0, 2,     -6.0,   0.0,    0.0,      2.0,  0.0,    0.0},
      { 3, 0, 2, 1, 2,      4.0,   0.0,    0.0,     -2.0,  0.0,    0.0},
      { 2, 1, 2, 2, 1,      3.0,   0.0,    0.0,     -2.0,  0.0,    0.0},
      { 4, 1, 2, 0, 2,      5.0,   0.0,    0.0,     -2.0,  0.0,    0.0},

   /* 671-678 */
      {-1,-1, 2, 6, 2,     -4.0,   0.0,    0.0,      2.0,  0.0,    0.0},
      {-1, 0, 2, 6, 1,     -4.0,   0.0,    0.0,      2.0,  0.0,    0.0},
      { 1,-1, 2, 4, 1,     -3.0,   0.0,    0.0,      2.0,  0.0,    0.0},
      { 1, 1, 2, 4, 2,      4.0,   0.0,    0.0,     -2.0,  0.0,    0.0},
      { 3, 1, 2, 2, 2,      3.0,   0.0,    0.0,     -1.0,  0.0,    0.0},
      { 5, 0, 2, 0, 1,     -3.0,   0.0,    0.0,      1.0,  0.0,    0.0},
      { 2,-1, 2, 4, 2,     -3.0,   0.0,    0.0,      1.0,  0.0,    0.0},
      { 2, 0, 2, 4, 1,     -3.0,   0.0,    0.0,      2.0,  0.0,    0.0}
   };

/* Number of terms in the luni-solar nutation model */
   const int NLS = (int) (sizeof xls / sizeof xls[0]);

/* ------------------------ */
/* Planetary nutation model */
/* ------------------------ */

/* The units for the sine and cosine coefficients are */
/* 0.1 microarcsecond                                 */

   static const struct {
      int nl,               /* coefficients of l, F, D and Omega */
          nf,
          nd,
          nom,
          nme,              /* coefficients of planetary longitudes */
          nve,
          nea,
          nma,
          nju,
          nsa,
          nur,
          nne,
          npa;              /* coefficient of general precession */
      int sp,cp;            /* longitude sin, cos coefficients */
      int se,ce;            /* obliquity sin, cos coefficients */
   } xpl[] = {

   /* 1-10 */
      { 0, 0, 0, 0, 0,  0,  8,-16, 4, 5, 0, 0, 0, 1440,   0,    0,   0},
      { 0, 0, 0, 0, 0,  0, -8, 16,-4,-5, 0, 0, 2,   56,-117,  -42, -40},
      { 0, 0, 0, 0, 0,  0,  8,-16, 4, 5, 0, 0, 2,  125, -43,    0, -54},
      { 0, 0, 0, 0, 0,  0,  0,  0, 0, 0,-1, 2, 2,    0,   5,    0,   0},
      { 0, 0, 0, 0, 0,  0, -4,  8,-1,-5, 0, 0, 2,    3,  -7,   -3,   0},
      { 0, 0, 0, 0, 0,  0,  4, -8, 3, 0, 0, 0, 1,    3,   0,    0,  -2},
      { 0, 1,-1, 1, 0,  0,  3, -8, 3, 0, 0, 0, 0, -114,   0,    0,  61},
      {-1, 0, 0, 0, 0, 10, -3,  0, 0, 0, 0, 0, 0, -219,  89,    0,   0},
      { 0, 0, 0, 0, 0,  0,  0,  0,-2, 6,-3, 0, 2,   -3,   0,    0,   0},
      { 0, 0, 0, 0, 0,  0,  4, -8, 3, 0, 0, 0, 0, -462,1604,    0,   0},

   /* 11-20 */
      { 0, 1,-1, 1, 0,  0, -5,  8,-3, 0, 0, 0, 0,   99,   0,    0, -53},
      { 0, 0, 0, 0, 0,  0, -4,  8,-3, 0, 0, 0, 1,   -3,   0,    0,   2},
      { 0, 0, 0, 0, 0,  0,  4, -8, 1, 5, 0, 0, 2,    0,   6,    2,   0},
      { 0, 0, 0, 0, 0, -5,  6,  4, 0, 0, 0, 0, 2,    3,   0,    0,   0},
      { 0, 0, 0, 0, 0,  0,  0,  0, 2,-5, 0, 0, 2,  -12,   0,    0,   0},
      { 0, 0, 0, 0, 0,  0,  0,  0, 2,-5, 0, 0, 1,   14,-218,  117,   8},
      { 0, 1,-1, 1, 0,  0, -1,  0, 2,-5, 0, 0, 0,   31,-481, -257, -17},
      { 0, 0, 0, 0, 0,  0,  0,  0, 2,-5, 0, 0, 0, -491, 128,    0,   0},
      { 0, 1,-1, 1, 0,  0, -1,  0,-2, 5, 0, 0, 0,-3084,5123, 2735,1647},
      { 0, 0, 0, 0, 0,  0,  0,  0,-2, 5, 0, 0, 1,-1444,2409,-1286,-771},

   /* 21-30 */
      { 0, 0, 0, 0, 0,  0,  0,  0,-2, 5, 0, 0, 2,   11, -24,  -11,  -9},
      { 2,-1,-1, 0, 0,  0,  3, -7, 0, 0, 0, 0, 0,   26,  -9,    0,   0},
      { 1, 0,-2, 0, 0, 19,-21,  3, 0, 0, 0, 0, 0,  103, -60,    0,   0},
      { 0, 1,-1, 1, 0,  2, -4,  0,-3, 0, 0, 0, 0,    0, -13,   -7,   0},
      { 1, 0,-1, 1, 0,  0, -1,  0, 2, 0, 0, 0, 0,  -26, -29,  -16,  14},
      { 0, 1,-1, 1, 0,  0, -1,  0,-4,10, 0, 0, 0,    9, -27,  -14,  -5},
      {-2, 0, 2, 1, 0,  0,  2,  0, 0,-5, 0, 0, 0,   12,   0,    0,  -6},
      { 0, 0, 0, 0, 0,  3, -7,  4, 0, 0, 0, 0, 0,   -7,   0,    0,   0},
      { 0,-1, 1, 0, 0,  0,  1,  0, 1,-1, 0, 0, 0,    0,  24,    0,   0},
      {-2, 0, 2, 1, 0,  0,  2,  0,-2, 0, 0, 0, 0,  284,   0,    0,-151},

   /* 31-40 */
      {-1, 0, 0, 0, 0, 18,-16,  0, 0, 0, 0, 0, 0,  226, 101,    0,   0},
      {-2, 1, 1, 2, 0,  0,  1,  0,-2, 0, 0, 0, 0,    0,  -8,   -2,   0},
      {-1, 1,-1, 1, 0, 18,-17,  0, 0, 0, 0, 0, 0,    0,  -6,   -3,   0},
      {-1, 0, 1, 1, 0,  0,  2, -2, 0, 0, 0, 0, 0,    5,   0,    0,  -3},
      { 0, 0, 0, 0, 0, -8, 13,  0, 0, 0, 0, 0, 2,  -41, 175,   76,  17},
      { 0, 2,-2, 2, 0, -8, 11,  0, 0, 0, 0, 0, 0,    0,  15,    6,   0},
      { 0, 0, 0, 0, 0, -8, 13,  0, 0, 0, 0, 0, 1,  425, 212, -133, 269},
      { 0, 1,-1, 1, 0, -8, 12,  0, 0, 0, 0, 0, 0, 1200, 598,  319,-641},
      { 0, 0, 0, 0, 0,  8,-13,  0, 0, 0, 0, 0, 0,  235, 334,    0,   0},
      { 0, 1,-1, 1, 0,  8,-14,  0, 0, 0, 0, 0, 0,   11, -12,   -7,  -6},

   /* 41-50 */
      { 0, 0, 0, 0, 0,  8,-13,  0, 0, 0, 0, 0, 1,    5,  -6,    3,   3},
      {-2, 0, 2, 1, 0,  0,  2,  0,-4, 5, 0, 0, 0,   -5,   0,    0,   3},
      {-2, 0, 2, 2, 0,  3, -3,  0, 0, 0, 0, 0, 0,    6,   0,    0,  -3},
      {-2, 0, 2, 0, 0,  0,  2,  0,-3, 1, 0, 0, 0,   15,   0,    0,   0},
      { 0, 0, 0, 1, 0,  3, -5,  0, 2, 0, 0, 0, 0,   13,   0,    0,  -7},
      {-2, 0, 2, 0, 0,  0,  2,  0,-4, 3, 0, 0, 0,   -6,  -9,    0,   0},
      { 0,-1, 1, 0, 0,  0,  0,  2, 0, 0, 0, 0, 0,  266, -78,    0,   0},
      { 0, 0, 0, 1, 0,  0, -1,  2, 0, 0, 0, 0, 0, -460,-435, -232, 246},
      { 0, 1,-1, 2, 0,  0, -2,  2, 0, 0, 0, 0, 0,    0,  15,    7,   0},
      {-1, 1, 0, 1, 0,  3, -5,  0, 0, 0, 0, 0, 0,   -3,   0,    0,   2},

   /* 51-60 */
      {-1, 0, 1, 0, 0,  3, -4,  0, 0, 0, 0, 0, 0,    0, 131,    0,   0},
      {-2, 0, 2, 0, 0,  0,  2,  0,-2,-2, 0, 0, 0,    4,   0,    0,   0},
      {-2, 2, 0, 2, 0,  0, -5,  9, 0, 0, 0, 0, 0,    0,   3,    0,   0},
      { 0, 1,-1, 1, 0,  0, -1,  0, 0, 0,-1, 0, 0,    0,   4,    2,   0},
      { 0, 0, 0, 0, 0,  0,  0,  0, 0, 0, 1, 0, 0,    0,   3,    0,   0},
      { 0, 1,-1, 1, 0,  0, -1,  0, 0, 0, 0, 2, 0,  -17, -19,  -10,   9},
      { 0, 0, 0, 0, 0,  0,  0,  0, 0, 0, 0, 2, 1,   -9, -11,    6,  -5},
      { 0, 0, 0, 0, 0,  0,  0,  0, 0, 0, 0, 2, 2,   -6,   0,    0,   3},
      {-1, 0, 1, 0, 0,  0,  3, -4, 0, 0, 0, 0, 0,  -16,   8,    0,   0},
      { 0,-1, 1, 0, 0,  0,  1,  0, 0, 2, 0, 0, 0,    0,   3,    0,   0},

   /* 61-70 */
      { 0, 1,-1, 2, 0,  0, -1,  0, 0, 2, 0, 0, 0,   11,  24,   11,  -5},
      { 0, 0, 0, 1, 0,  0, -9, 17, 0, 0, 0, 0, 0,   -3,  -4,   -2,   1},
      { 0, 0, 0, 2, 0, -3,  5,  0, 0, 0, 0, 0, 0,    3,   0,    0,  -1},
      { 0, 1,-1, 1, 0,  0, -1,  0,-1, 2, 0, 0, 0,    0,  -8,   -4,   0},
      { 0, 0, 0, 0, 0,  0,  0,  0, 1,-2, 0, 0, 0,    0,   3,    0,   0},
      { 1, 0,-2, 0, 0, 17,-16,  0,-2, 0, 0, 0, 0,    0,   5,    0,   0},
      { 0, 1,-1, 1, 0,  0, -1,  0, 1,-3, 0, 0, 0,    0,   3,    2,   0},
      {-2, 0, 2, 1, 0,  0,  5, -6, 0, 0, 0, 0, 0,   -6,   4,    2,   3},
      { 0,-2, 2, 0, 0,  0,  9,-13, 0, 0, 0, 0, 0,   -3,  -5,    0,   0},
      { 0, 1,-1, 2, 0,  0, -1,  0, 0, 1, 0, 0, 0,   -5,   0,    0,   2},

   /* 71-80 */
      { 0, 0, 0, 1, 0,  0,  0,  0, 0, 1, 0, 0, 0,    4,  24,   13,  -2},
      { 0,-1, 1, 0, 0,  0,  1,  0, 0, 1, 0, 0, 0,  -42,  20,    0,   0},
      { 0,-2, 2, 0, 0,  5, -6,  0, 0, 0, 0, 0, 0,  -10, 233,    0,   0},
      { 0,-1, 1, 1, 0,  5, -7,  0, 0, 0, 0, 0, 0,   -3,   0,    0,   1},
      {-2, 0, 2, 0, 0,  6, -8,  0, 0, 0, 0, 0, 0,   78, -18,    0,   0},
      { 2, 1,-3, 1, 0, -6,  7,  0, 0, 0, 0, 0, 0,    0,   3,    1,   0},
      { 0, 0, 0, 2, 0,  0,  0,  0, 1, 0, 0, 0, 0,    0,  -3,   -1,   0},
      { 0,-1, 1, 1, 0,  0,  1,  0, 1, 0, 0, 0, 0,    0,  -4,   -2,   1},
      { 0, 1,-1, 1, 0,  0, -1,  0, 0, 0, 2, 0, 0,    0,  -8,   -4,  -1},
      { 0, 0, 0, 0, 0,  0,  0,  0, 0, 0, 2, 0, 1,    0,  -5,    3,   0},

   /* 81-90 */
      { 0, 0, 0, 0, 0,  0,  0,  0, 0, 0, 2, 0, 2,   -7,   0,    0,   3},
      { 0, 0, 0, 0, 0,  0, -8, 15, 0, 0, 0, 0, 2,  -14,   8,    3,   6},
      { 0, 0, 0, 0, 0,  0, -8, 15, 0, 0, 0, 0, 1,    0,   8,   -4,   0},
      { 0, 1,-1, 1, 0,  0, -9, 15, 0, 0, 0, 0, 0,    0,  19,   10,   0},
      { 0, 0, 0, 0, 0,  0,  8,-15, 0, 0, 0, 0, 0,   45, -22,    0,   0},
      { 1,-1,-1, 0, 0,  0,  8,-15, 0, 0, 0, 0, 0,   -3,   0,    0,   0},
      { 2, 0,-2, 0, 0,  2, -5,  0, 0, 0, 0, 0, 0,    0,  -3,    0,   0},
      {-2, 0, 2, 0, 0,  0,  2,  0,-5, 5, 0, 0, 0,    0,   3,    0,   0},
      { 2, 0,-2, 1, 0,  0, -6,  8, 0, 0, 0, 0, 0,    3,   5,    3,  -2},
      { 2, 0,-2, 1, 0,  0, -2,  0, 3, 0, 0, 0, 0,   89, -16,   -9, -48},

   /* 91-100 */
      {-2, 1, 1, 0, 0,  0,  1,  0,-3, 0, 0, 0, 0,    0,   3,    0,   0},
      {-2, 1, 1, 1, 0,  0,  1,  0,-3, 0, 0, 0, 0,   -3,   7,    4,   2},
      {-2, 0, 2, 0, 0,  0,  2,  0,-3, 0, 0, 0, 0, -349, -62,    0,   0},
      {-2, 0, 2, 0, 0,  0,  6, -8, 0, 0, 0, 0, 0,  -15,  22,    0,   0},
      {-2, 0, 2, 0, 0,  0,  2,  0,-1,-5, 0, 0, 0,   -3,   0,    0,   0},
      {-1, 0, 1, 0, 0,  0,  1,  0,-1, 0, 0, 0, 0,  -53,   0,    0,   0},
      {-1, 1, 1, 1, 0,-20, 20,  0, 0, 0, 0, 0, 0,    5,   0,    0,  -3},
      { 1, 0,-2, 0, 0, 20,-21,  0, 0, 0, 0, 0, 0,    0,  -8,    0,   0},
      { 0, 0, 0, 1, 0,  0,  8,-15, 0, 0, 0, 0, 0,   15,  -7,   -4,  -8},
      { 0, 2,-2, 1, 0,  0,-10, 15, 0, 0, 0, 0, 0,   -3,   0,    0,   1},

   /* 101-110 */
      { 0,-1, 1, 0, 0,  0,  1,  0, 1, 0, 0, 0, 0,  -21, -78,    0,   0},
      { 0, 0, 0, 1, 0,  0,  0,  0, 1, 0, 0, 0, 0,   20, -70,  -37, -11},
      { 0, 1,-1, 2, 0,  0, -1,  0, 1, 0, 0, 0, 0,    0,   6,    3,   0},
      { 0, 1,-1, 1, 0,  0, -1,  0,-2, 4, 0, 0, 0,    5,   3,    2,  -2},
      { 2, 0,-2, 1, 0, -6,  8,  0, 0, 0, 0, 0, 0,  -17,  -4,   -2,   9},
      { 0,-2, 2, 1, 0,  5, -6,  0, 0, 0, 0, 0, 0,    0,   6,    3,   0},
      { 0, 0, 0, 0, 0,  0,  0,  0, 0,-1, 0, 0, 1,   32,  15,   -8,  17},
      { 0, 1,-1, 1, 0,  0, -1,  0, 0,-1, 0, 0, 0,  174,  84,   45, -93},
      { 0, 0, 0, 0, 0,  0,  0,  0, 0, 1, 0, 0, 0,   11,  56,    0,   0},
      { 0, 1,-1, 1, 0,  0, -1,  0, 0, 1, 0, 0, 0,  -66, -12,   -6,  35},

   /* 111-120 */
      { 0, 0, 0, 0, 0,  0,  0,  0, 0, 1, 0, 0, 1,   47,   8,    4, -25},
      { 0, 0, 0, 0, 0,  0,  0,  0, 0, 1, 0, 0, 2,    0,   8,    4,   0},
      { 0, 2,-2, 1, 0,  0, -9, 13, 0, 0, 0, 0, 0,   10, -22,  -12,  -5},
      { 0, 0, 0, 1, 0,  0,  7,-13, 0, 0, 0, 0, 0,   -3,   0,    0,   2},
      {-2, 0, 2, 0, 0,  0,  5, -6, 0, 0, 0, 0, 0,  -24,  12,    0,   0},
      { 0, 0, 0, 0, 0,  0,  9,-17, 0, 0, 0, 0, 0,    5,  -6,    0,   0},
      { 0, 0, 0, 0, 0,  0, -9, 17, 0, 0, 0, 0, 2,    3,   0,    0,  -2},
      { 1, 0,-1, 1, 0,  0, -3,  4, 0, 0, 0, 0, 0,    4,   3,    1,  -2},
      { 1, 0,-1, 1, 0, -3,  4,  0, 0, 0, 0, 0, 0,    0,  29,   15,   0},
      { 0, 0, 0, 2, 0,  0, -1,  2, 0, 0, 0, 0, 0,   -5,  -4,   -2,   2},

   /* 121-130 */
      { 0,-1, 1, 1, 0,  0,  0,  2, 0, 0, 0, 0, 0,    8,  -3,   -1,  -5},
      { 0,-2, 2, 0, 1,  0, -2,  0, 0, 0, 0, 0, 0,    0,  -3,    0,   0},
      { 0, 0, 0, 0, 0,  3, -5,  0, 2, 0, 0, 0, 0,   10,   0,    0,   0},
      {-2, 0, 2, 1, 0,  0,  2,  0,-3, 1, 0, 0, 0,    3,   0,    0,  -2},
      {-2, 0, 2, 1, 0,  3, -3,  0, 0, 0, 0, 0, 0,   -5,   0,    0,   3},
      { 0, 0, 0, 1, 0,  8,-13,  0, 0, 0, 0, 0, 0,   46,  66,   35, -25},
      { 0,-1, 1, 0, 0,  8,-12,  0, 0, 0, 0, 0, 0,  -14,   7,    0,   0},
      { 0, 2,-2, 1, 0, -8, 11,  0, 0, 0, 0, 0, 0,    0,   3,    2,   0},
      {-1, 0, 1, 0, 0,  0,  2, -2, 0, 0, 0, 0, 0,   -5,   0,    0,   0},
      {-1, 0, 0, 1, 0, 18,-16,  0, 0, 0, 0, 0, 0,  -68, -34,  -18,  36},

   /* 131-140 */
      { 0, 1,-1, 1, 0,  0, -1,  0,-1, 1, 0, 0, 0,    0,  14,    7,   0},
      { 0, 0, 0, 1, 0,  3, -7,  4, 0, 0, 0, 0, 0,   10,  -6,   -3,  -5},
      {-2, 1, 1, 1, 0,  0, -3,  7, 0, 0, 0, 0, 0,   -5,  -4,   -2,   3},
      { 0, 1,-1, 2, 0,  0, -1,  0,-2, 5, 0, 0, 0,   -3,   5,    2,   1},
      { 0, 0, 0, 1, 0,  0,  0,  0,-2, 5, 0, 0, 0,   76,  17,    9, -41},
      { 0, 0, 0, 1, 0,  0, -4,  8,-3, 0, 0, 0, 0,   84, 298,  159, -45},
      { 1, 0, 0, 1, 0,-10,  3,  0, 0, 0, 0, 0, 0,    3,   0,    0,  -1},
      { 0, 2,-2, 1, 0,  0, -2,  0, 0, 0, 0, 0, 0,   -3,   0,    0,   2},
      {-1, 0, 0, 1, 0, 10, -3,  0, 0, 0, 0, 0, 0,   -3,   0,    0,   1},
      { 0, 0, 0, 1, 0,  0,  4, -8, 3, 0, 0, 0, 0,  -82, 292,  156,  44},

   /* 141-150 */
      { 0, 0, 0, 1, 0,  0,  0,  0, 2,-5, 0, 0, 0,  -73,  17,    9,  39},
      { 0,-1, 1, 0, 0,  0,  1,  0, 2,-5, 0, 0, 0,   -9, -16,    0,   0},
      { 2,-1,-1, 1, 0,  0,  3, -7, 0, 0, 0, 0, 0,    3,   0,   -1,  -2},
      {-2, 0, 2, 0, 0,  0,  2,  0, 0,-5, 0, 0, 0,   -3,   0,    0,   0},
      { 0, 0, 0, 1, 0, -3,  7, -4, 0, 0, 0, 0, 0,   -9,  -5,   -3,   5},
      {-2, 0, 2, 0, 0,  0,  2,  0,-2, 0, 0, 0, 0, -439,   0,    0,   0},
      { 1, 0, 0, 1, 0,-18, 16,  0, 0, 0, 0, 0, 0,   57, -28,  -15, -30},
      {-2, 1, 1, 1, 0,  0,  1,  0,-2, 0, 0, 0, 0,    0,  -6,   -3,   0},
      { 0, 1,-1, 2, 0, -8, 12,  0, 0, 0, 0, 0, 0,   -4,   0,    0,   2},
      { 0, 0, 0, 1, 0, -8, 13,  0, 0, 0, 0, 0, 0,  -40,  57,   30,  21},

   /* 151-160 */
      { 0, 0, 0, 0, 0,  0,  1, -2, 0, 0, 0, 0, 1,   23,   7,    3, -13},
      { 0, 1,-1, 1, 0,  0,  0, -2, 0, 0, 0, 0, 0,  273,  80,   43,-146},
      { 0, 0, 0, 0, 0,  0,  1, -2, 0, 0, 0, 0, 0, -449, 430,    0,   0},
      { 0, 1,-1, 1, 0,  0, -2,  2, 0, 0, 0, 0, 0,   -8, -47,  -25,   4},
      { 0, 0, 0, 0, 0,  0, -1,  2, 0, 0, 0, 0, 1,    6,  47,   25,  -3},
      {-1, 0, 1, 1, 0,  3, -4,  0, 0, 0, 0, 0, 0,    0,  23,   13,   0},
      {-1, 0, 1, 1, 0,  0,  3, -4, 0, 0, 0, 0, 0,   -3,   0,    0,   2},
      { 0, 1,-1, 1, 0,  0, -1,  0, 0,-2, 0, 0, 0,    3,  -4,   -2,  -2},
      { 0, 1,-1, 1, 0,  0, -1,  0, 0, 2, 0, 0, 0,  -48,-110,  -59,  26},
      { 0, 0, 0, 0, 0,  0,  0,  0, 0, 2, 0, 0, 1,   51, 114,   61, -27},

   /* 161-170 */
      { 0, 0, 0, 0, 0,  0,  0,  0, 0, 2, 0, 0, 2, -133,   0,    0,  57},
      { 0, 1,-1, 0, 0,  3, -6,  0, 0, 0, 0, 0, 0,    0,   4,    0,   0},
      { 0, 0, 0, 1, 0, -3,  5,  0, 0, 0, 0, 0, 0,  -21,  -6,   -3,  11},
      { 0, 1,-1, 2, 0, -3,  4,  0, 0, 0, 0, 0, 0,    0,  -3,   -1,   0},
      { 0, 0, 0, 1, 0,  0, -2,  4, 0, 0, 0, 0, 0,  -11, -21,  -11,   6},
      { 0, 2,-2, 1, 0, -5,  6,  0, 0, 0, 0, 0, 0,  -18,-436, -233,   9},
      { 0,-1, 1, 0, 0,  5, -7,  0, 0, 0, 0, 0, 0,   35,  -7,    0,   0},
      { 0, 0, 0, 1, 0,  5, -8,  0, 0, 0, 0, 0, 0,    0,   5,    3,   0},
      {-2, 0, 2, 1, 0,  6, -8,  0, 0, 0, 0, 0, 0,   11,  -3,   -1,  -6},
      { 0, 0, 0, 1, 0,  0, -8, 15, 0, 0, 0, 0, 0,   -5,  -3,   -1,   3},

   /* 171-180 */
      {-2, 0, 2, 1, 0,  0,  2,  0,-3, 0, 0, 0, 0,  -53,  -9,   -5,  28},
      {-2, 0, 2, 1, 0,  0,  6, -8, 0, 0, 0, 0, 0,    0,   3,    2,   1},
      { 1, 0,-1, 1, 0,  0, -1,  0, 1, 0, 0, 0, 0,    4,   0,    0,  -2},
      { 0, 0, 0, 0, 0,  0,  0,  0, 3,-5, 0, 0, 0,    0,  -4,    0,   0},
      { 0, 1,-1, 1, 0,  0, -1,  0,-1, 0, 0, 0, 0,  -50, 194,  103,  27},
      { 0, 0, 0, 0, 0,  0,  0,  0,-1, 0, 0, 0, 1,  -13,  52,   28,   7},
      { 0, 0, 0, 0, 0,  0,  0,  0, 1, 0, 0, 0, 0,  -91, 248,    0,   0},
      { 0, 0, 0, 0, 0,  0,  0,  0, 1, 0, 0, 0, 1,    6,  49,   26,  -3},
      { 0, 1,-1, 1, 0,  0, -1,  0, 1, 0, 0, 0, 0,   -6, -47,  -25,   3},
      { 0, 0, 0, 0, 0,  0,  0,  0, 1, 0, 0, 0, 1,    0,   5,    3,   0},

   /* 181-190 */
      { 0, 0, 0, 0, 0,  0,  0,  0, 1, 0, 0, 0, 2,   52,  23,   10, -23},
      { 0, 1,-1, 2, 0,  0, -1,  0, 0,-1, 0, 0, 0,   -3,   0,    0,   1},
      { 0, 0, 0, 1, 0,  0,  0,  0, 0,-1, 0, 0, 0,    0,   5,    3,   0},
      { 0,-1, 1, 0, 0,  0,  1,  0, 0,-1, 0, 0, 0,   -4,   0,    0,   0},
      { 0, 0, 0, 0, 0,  0, -7, 13, 0, 0, 0, 0, 2,   -4,   8,    3,   2},
      { 0, 0, 0, 0, 0,  0,  7,-13, 0, 0, 0, 0, 0,   10,   0,    0,   0},
      { 2, 0,-2, 1, 0,  0, -5,  6, 0, 0, 0, 0, 0,    3,   0,    0,  -2},
      { 0, 2,-2, 1, 0,  0, -8, 11, 0, 0, 0, 0, 0,    0,   8,    4,   0},
      { 0, 2,-2, 1,-1,  0,  2,  0, 0, 0, 0, 0, 0,    0,   8,    4,   1},
      {-2, 0, 2, 0, 0,  0,  4, -4, 0, 0, 0, 0, 0,   -4,   0,    0,   0},

   /* 191-200 */
      { 0, 0, 0, 0, 0,  0,  0,  0, 2,-2, 0, 0, 0,   -4,   0,    0,   0},
      { 0, 1,-1, 1, 0,  0, -1,  0, 0, 3, 0, 0, 0,   -8,   4,    2,   4},
      { 0, 0, 0, 0, 0,  0,  0,  0, 0, 3, 0, 0, 1,    8,  -4,   -2,  -4},
      { 0, 0, 0, 0, 0,  0,  0,  0, 0, 3, 0, 0, 2,    0,  15,    7,   0},
      {-2, 0, 2, 0, 0,  3, -3,  0, 0, 0, 0, 0, 0, -138,   0,    0,   0},
      { 0, 0, 0, 2, 0,  0, -4,  8,-3, 0, 0, 0, 0,    0,  -7,   -3,   0},
      { 0, 0, 0, 2, 0,  0,  4, -8, 3, 0, 0, 0, 0,    0,  -7,   -3,   0},
      { 2, 0,-2, 1, 0,  0, -2,  0, 2, 0, 0, 0, 0,   54,   0,    0, -29},
      { 0, 1,-1, 2, 0,  0, -1,  0, 2, 0, 0, 0, 0,    0,  10,    4,   0},
      { 0, 1,-1, 2, 0,  0,  0, -2, 0, 0, 0, 0, 0,   -7,   0,    0,   3},

   /* 201-210 */
      { 0, 0, 0, 1, 0,  0,  1, -2, 0, 0, 0, 0, 0,  -37,  35,   19,  20},
      { 0,-1, 1, 0, 0,  0,  2, -2, 0, 0, 0, 0, 0,    0,   4,    0,   0},
      { 0,-1, 1, 0, 0,  0,  1,  0, 0,-2, 0, 0, 0,   -4,   9,    0,   0},
      { 0, 2,-2, 1, 0,  0, -2,  0, 0, 2, 0, 0, 0,    8,   0,    0,  -4},
      { 0, 1,-1, 1, 0,  3, -6,  0, 0, 0, 0, 0, 0,   -9, -14,   -8,   5},
      { 0, 0, 0, 0, 0,  3, -5,  0, 0, 0, 0, 0, 1,   -3,  -9,   -5,   3},
      { 0, 0, 0, 0, 0,  3, -5,  0, 0, 0, 0, 0, 0, -145,  47,    0,   0},
      { 0, 1,-1, 1, 0, -3,  4,  0, 0, 0, 0, 0, 0,  -10,  40,   21,   5},
      { 0, 0, 0, 0, 0, -3,  5,  0, 0, 0, 0, 0, 1,   11, -49,  -26,  -7},
      { 0, 0, 0, 0, 0, -3,  5,  0, 0, 0, 0, 0, 2,-2150,   0,    0, 932},

   /* 211-220 */
      { 0, 2,-2, 2, 0, -3,  3,  0, 0, 0, 0, 0, 0,  -12,   0,    0,   5},
      { 0, 0, 0, 0, 0, -3,  5,  0, 0, 0, 0, 0, 2,   85,   0,    0, -37},
      { 0, 0, 0, 0, 0,  0,  2, -4, 0, 0, 0, 0, 1,    4,   0,    0,  -2},
      { 0, 1,-1, 1, 0,  0,  1, -4, 0, 0, 0, 0, 0,    3,   0,    0,  -2},
      { 0, 0, 0, 0, 0,  0,  2, -4, 0, 0, 0, 0, 0,  -86, 153,    0,   0},
      { 0, 0, 0, 0, 0,  0, -2,  4, 0, 0, 0, 0, 1,   -6,   9,    5,   3},
      { 0, 1,-1, 1, 0,  0, -3,  4, 0, 0, 0, 0, 0,    9, -13,   -7,  -5},
      { 0, 0, 0, 0, 0,  0, -2,  4, 0, 0, 0, 0, 1,   -8,  12,    6,   4},
      { 0, 0, 0, 0, 0,  0, -2,  4, 0, 0, 0, 0, 2,  -51,   0,    0,  22},
      { 0, 0, 0, 0, 0, -5,  8,  0, 0, 0, 0, 0, 2,  -11,-268, -116,   5},

   /* 221-230 */
      { 0, 2,-2, 2, 0, -5,  6,  0, 0, 0, 0, 0, 0,    0,  12,    5,   0},
      { 0, 0, 0, 0, 0, -5,  8,  0, 0, 0, 0, 0, 2,    0,   7,    3,   0},
      { 0, 0, 0, 0, 0, -5,  8,  0, 0, 0, 0, 0, 1,   31,   6,    3, -17},
      { 0, 1,-1, 1, 0, -5,  7,  0, 0, 0, 0, 0, 0,  140,  27,   14, -75},
      { 0, 0, 0, 0, 0, -5,  8,  0, 0, 0, 0, 0, 1,   57,  11,    6, -30},
      { 0, 0, 0, 0, 0,  5, -8,  0, 0, 0, 0, 0, 0,  -14, -39,    0,   0},
      { 0, 1,-1, 2, 0,  0, -1,  0,-1, 0, 0, 0, 0,    0,  -6,   -2,   0},
      { 0, 0, 0, 1, 0,  0,  0,  0,-1, 0, 0, 0, 0,    4,  15,    8,  -2},
      { 0,-1, 1, 0, 0,  0,  1,  0,-1, 0, 0, 0, 0,    0,   4,    0,   0},
      { 0, 2,-2, 1, 0,  0, -2,  0, 1, 0, 0, 0, 0,   -3,   0,    0,   1},

   /* 231-240 */
      { 0, 0, 0, 0, 0,  0, -6, 11, 0, 0, 0, 0, 2,    0,  11,    5,   0},
      { 0, 0, 0, 0, 0,  0,  6,-11, 0, 0, 0, 0, 0,    9,   6,    0,   0},
      { 0, 0, 0, 0,-1,  0,  4,  0, 0, 0, 0, 0, 2,   -4,  10,    4,   2},
      { 0, 0, 0, 0, 1,  0, -4,  0, 0, 0, 0, 0, 0,    5,   3,    0,   0},
      { 2, 0,-2, 1, 0, -3,  3,  0, 0, 0, 0, 0, 0,   16,   0,    0,  -9},
      {-2, 0, 2, 0, 0,  0,  2,  0, 0,-2, 0, 0, 0,   -3,   0,    0,   0},
      { 0, 2,-2, 1, 0,  0, -7,  9, 0, 0, 0, 0, 0,    0,   3,    2,  -1},
      { 0, 0, 0, 0, 0,  0,  0,  0, 4,-5, 0, 0, 2,    7,   0,    0,  -3},
      { 0, 0, 0, 0, 0,  0,  0,  0, 2, 0, 0, 0, 0,  -25,  22,    0,   0},
      { 0, 0, 0, 0, 0,  0,  0,  0, 2, 0, 0, 0, 1,   42, 223,  119, -22},

   /* 241-250 */
      { 0, 1,-1, 1, 0,  0, -1,  0, 2, 0, 0, 0, 0,  -27,-143,  -77,  14},
      { 0, 0, 0, 0, 0,  0,  0,  0, 2, 0, 0, 0, 1,    9,  49,   26,  -5},
      { 0, 0, 0, 0, 0,  0,  0,  0, 2, 0, 0, 0, 2,-1166,   0,    0, 505},
      { 0, 2,-2, 2, 0,  0, -2,  0, 2, 0, 0, 0, 0,   -5,   0,    0,   2},
      { 0, 0, 0, 0, 0,  0,  0,  0, 0, 5, 0, 0, 2,   -6,   0,    0,   3},
      { 0, 0, 0, 1, 0,  3, -5,  0, 0, 0, 0, 0, 0,   -8,   0,    1,   4},
      { 0,-1, 1, 0, 0,  3, -4,  0, 0, 0, 0, 0, 0,    0,  -4,    0,   0},
      { 0, 2,-2, 1, 0, -3,  3,  0, 0, 0, 0, 0, 0,  117,   0,    0, -63},
      { 0, 0, 0, 1, 0,  0,  2, -4, 0, 0, 0, 0, 0,   -4,   8,    4,   2},
      { 0, 2,-2, 1, 0,  0, -4,  4, 0, 0, 0, 0, 0,    3,   0,    0,  -2},

   /* 251-260 */
      { 0, 1,-1, 2, 0, -5,  7,  0, 0, 0, 0, 0, 0,   -5,   0,    0,   2},
      { 0, 0, 0, 0, 0,  0,  3, -6, 0, 0, 0, 0, 0,    0,  31,    0,   0},
      { 0, 0, 0, 0, 0,  0, -3,  6, 0, 0, 0, 0, 1,   -5,   0,    1,   3},
      { 0, 1,-1, 1, 0,  0, -4,  6, 0, 0, 0, 0, 0,    4,   0,    0,  -2},
      { 0, 0, 0, 0, 0,  0, -3,  6, 0, 0, 0, 0, 1,   -4,   0,    0,   2},
      { 0, 0, 0, 0, 0,  0, -3,  6, 0, 0, 0, 0, 2,  -24, -13,   -6,  10},
      { 0,-1, 1, 0, 0,  2, -2,  0, 0, 0, 0, 0, 0,    3,   0,    0,   0},
      { 0, 0, 0, 1, 0,  2, -3,  0, 0, 0, 0, 0, 0,    0, -32,  -17,   0},
      { 0, 0, 0, 0, 0,  0, -5,  9, 0, 0, 0, 0, 2,    8,  12,    5,  -3},
      { 0, 0, 0, 0, 0,  0, -5,  9, 0, 0, 0, 0, 1,    3,   0,    0,  -1},

   /* 261-270 */
      { 0, 0, 0, 0, 0,  0,  5, -9, 0, 0, 0, 0, 0,    7,  13,    0,   0},
      { 0,-1, 1, 0, 0,  0,  1,  0,-2, 0, 0, 0, 0,   -3,  16,    0,   0},
      { 0, 2,-2, 1, 0,  0, -2,  0, 2, 0, 0, 0, 0,   50,   0,    0, -27},
      {-2, 1, 1, 1, 0,  0,  1,  0, 0, 0, 0, 0, 0,    0,  -5,   -3,   0},
      { 0,-2, 2, 0, 0,  3, -3,  0, 0, 0, 0, 0, 0,   13,   0,    0,   0},
      { 0, 0, 0, 0, 0, -6, 10,  0, 0, 0, 0, 0, 1,    0,   5,    3,   1},
      { 0, 0, 0, 0, 0, -6, 10,  0, 0, 0, 0, 0, 2,   24,   5,    2, -11},
      { 0, 0, 0, 0, 0, -2,  3,  0, 0, 0, 0, 0, 2,    5, -11,   -5,  -2},
      { 0, 0, 0, 0, 0, -2,  3,  0, 0, 0, 0, 0, 1,   30,  -3,   -2, -16},
      { 0, 1,-1, 1, 0, -2,  2,  0, 0, 0, 0, 0, 0,   18,   0,    0,  -9},

   /* 271-280 */
      { 0, 0, 0, 0, 0,  2, -3,  0, 0, 0, 0, 0, 0,    8, 614,    0,   0},
      { 0, 0, 0, 0, 0,  2, -3,  0, 0, 0, 0, 0, 1,    3,  -3,   -1,  -2},
      { 0, 0, 0, 0, 0,  0,  0,  0, 3, 0, 0, 0, 1,    6,  17,    9,  -3},
      { 0, 1,-1, 1, 0,  0, -1,  0, 3, 0, 0, 0, 0,   -3,  -9,   -5,   2},
      { 0, 0, 0, 0, 0,  0,  0,  0, 3, 0, 0, 0, 1,    0,   6,    3,  -1},
      { 0, 0, 0, 0, 0,  0,  0,  0, 3, 0, 0, 0, 2, -127,  21,    9,  55},
      { 0, 0, 0, 0, 0,  0,  4, -8, 0, 0, 0, 0, 0,    3,   5,    0,   0},
      { 0, 0, 0, 0, 0,  0, -4,  8, 0, 0, 0, 0, 2,   -6, -10,   -4,   3},
      { 0,-2, 2, 0, 0,  0,  2,  0,-2, 0, 0, 0, 0,    5,   0,    0,   0},
      { 0, 0, 0, 0, 0,  0, -4,  7, 0, 0, 0, 0, 2,   16,   9,    4,  -7},

   /* 281-290 */
      { 0, 0, 0, 0, 0,  0, -4,  7, 0, 0, 0, 0, 1,    3,   0,    0,  -2},
      { 0, 0, 0, 0, 0,  0,  4, -7, 0, 0, 0, 0, 0,    0,  22,    0,   0},
      { 0, 0, 0, 1, 0, -2,  3,  0, 0, 0, 0, 0, 0,    0,  19,   10,   0},
      { 0, 2,-2, 1, 0,  0, -2,  0, 3, 0, 0, 0, 0,    7,   0,    0,  -4},
      { 0, 0, 0, 0, 0,  0, -5, 10, 0, 0, 0, 0, 2,    0,  -5,   -2,   0},
      { 0, 0, 0, 1, 0, -1,  2,  0, 0, 0, 0, 0, 0,    0,   3,    1,   0},
      { 0, 0, 0, 0, 0,  0,  0,  0, 4, 0, 0, 0, 2,   -9,   3,    1,   4},
      { 0, 0, 0, 0, 0,  0, -3,  5, 0, 0, 0, 0, 2,   17,   0,    0,  -7},
      { 0, 0, 0, 0, 0,  0, -3,  5, 0, 0, 0, 0, 1,    0,  -3,   -2,  -1},
      { 0, 0, 0, 0, 0,  0,  3, -5, 0, 0, 0, 0, 0,  -20,  34,    0,   0},

   /* 291-300 */
      { 0, 0, 0, 0, 0,  1, -2,  0, 0, 0, 0, 0, 1,  -10,   0,    1,   5},
      { 0, 1,-1, 1, 0,  1, -3,  0, 0, 0, 0, 0, 0,   -4,   0,    0,   2},
      { 0, 0, 0, 0, 0,  1, -2,  0, 0, 0, 0, 0, 0,   22, -87,    0,   0},
      { 0, 0, 0, 0, 0, -1,  2,  0, 0, 0, 0, 0, 1,   -4,   0,    0,   2},
      { 0, 0, 0, 0, 0, -1,  2,  0, 0, 0, 0, 0, 2,   -3,  -6,   -2,   1},
      { 0, 0, 0, 0, 0, -7, 11,  0, 0, 0, 0, 0, 2,  -16,  -3,   -1,   7},
      { 0, 0, 0, 0, 0, -7, 11,  0, 0, 0, 0, 0, 1,    0,  -3,   -2,   0},
      { 0,-2, 2, 0, 0,  4, -4,  0, 0, 0, 0, 0, 0,    4,   0,    0,   0},
      { 0, 0, 0, 0, 0,  0,  2, -3, 0, 0, 0, 0, 0,  -68,  39,    0,   0},
      { 0, 2,-2, 1, 0, -4,  4,  0, 0, 0, 0, 0, 0,   27,   0,    0, -14},

   /* 301-310 */
      { 0,-1, 1, 0, 0,  4, -5,  0, 0, 0, 0, 0, 0,    0,  -4,    0,   0},
      { 0, 0, 0, 0, 0,  0,  1, -1, 0, 0, 0, 0, 0,  -25,   0,    0,   0},
      { 0, 0, 0, 0, 0, -4,  7,  0, 0, 0, 0, 0, 1,  -12,  -3,   -2,   6},
      { 0, 1,-1, 1, 0, -4,  6,  0, 0, 0, 0, 0, 0,    3,   0,    0,  -1},
      { 0, 0, 0, 0, 0, -4,  7,  0, 0, 0, 0, 0, 2,    3,  66,   29,  -1},
      { 0, 0, 0, 0, 0, -4,  6,  0, 0, 0, 0, 0, 2,  490,   0,    0,-213},
      { 0, 0, 0, 0, 0, -4,  6,  0, 0, 0, 0, 0, 1,  -22,  93,   49,  12},
      { 0, 1,-1, 1, 0, -4,  5,  0, 0, 0, 0, 0, 0,   -7,  28,   15,   4},
      { 0, 0, 0, 0, 0, -4,  6,  0, 0, 0, 0, 0, 1,   -3,  13,    7,   2},
      { 0, 0, 0, 0, 0,  4, -6,  0, 0, 0, 0, 0, 0,  -46,  14,    0,   0},

   /* 311-320 */
      {-2, 0, 2, 0, 0,  2, -2,  0, 0, 0, 0, 0, 0,   -5,   0,    0,   0},
      { 0, 0, 0, 0, 0,  0,  0,  1, 0, 0, 0, 0, 0,    2,   1,    0,   0},
      { 0,-1, 1, 0, 0,  1,  0,  0, 0, 0, 0, 0, 0,    0,  -3,    0,   0},
      { 0, 0, 0, 1, 0,  1, -1,  0, 0, 0, 0, 0, 0,  -28,   0,    0,  15},
      { 0, 0, 0, 0, 0,  0, -1,  0, 5, 0, 0, 0, 2,    5,   0,    0,  -2},
      { 0, 0, 0, 0, 0,  0,  1, -3, 0, 0, 0, 0, 0,    0,   3,    0,   0},
      { 0, 0, 0, 0, 0,  0, -1,  3, 0, 0, 0, 0, 2,  -11,   0,    0,   5},
      { 0, 0, 0, 0, 0,  0, -7, 12, 0, 0, 0, 0, 2,    0,   3,    1,   0},
      { 0, 0, 0, 0, 0, -1,  1,  0, 0, 0, 0, 0, 2,   -3,   0,    0,   1},
      { 0, 0, 0, 0, 0, -1,  1,  0, 0, 0, 0, 0, 1,   25, 106,   57, -13},

   /* 321-330 */
      { 0, 1,-1, 1, 0, -1,  0,  0, 0, 0, 0, 0, 0,    5,  21,   11,  -3},
      { 0, 0, 0, 0, 0,  1, -1,  0, 0, 0, 0, 0, 0, 1485,   0,    0,   0},
      { 0, 0, 0, 0, 0,  1, -1,  0, 0, 0, 0, 0, 1,   -7, -32,  -17,   4},
      { 0, 1,-1, 1, 0,  1, -2,  0, 0, 0, 0, 0, 0,    0,   5,    3,   0},
      { 0, 0, 0, 0, 0,  0, -2,  5, 0, 0, 0, 0, 2,   -6,  -3,   -2,   3},
      { 0, 0, 0, 0, 0,  0, -1,  0, 4, 0, 0, 0, 2,   30,  -6,   -2, -13},
      { 0, 0, 0, 0, 0,  0,  1,  0,-4, 0, 0, 0, 0,   -4,   4,    0,   0},
      { 0, 0, 0, 1, 0, -1,  1,  0, 0, 0, 0, 0, 0,  -19,   0,    0,  10},
      { 0, 0, 0, 0, 0,  0, -6, 10, 0, 0, 0, 0, 2,    0,   4,    2,  -1},
      { 0, 0, 0, 0, 0,  0, -6, 10, 0, 0, 0, 0, 0,    0,   3,    0,   0},

   /* 331-340 */
      { 0, 2,-2, 1, 0,  0, -3,  0, 3, 0, 0, 0, 0,    4,   0,    0,  -2},
      { 0, 0, 0, 0, 0,  0, -3,  7, 0, 0, 0, 0, 2,    0,  -3,   -1,   0},
      {-2, 0, 2, 0, 0,  4, -4,  0, 0, 0, 0, 0, 0,   -3,   0,    0,   0},
      { 0, 0, 0, 0, 0,  0, -5,  8, 0, 0, 0, 0, 2,    5,   3,    1,  -2},
      { 0, 0, 0, 0, 0,  0,  5, -8, 0, 0, 0, 0, 0,    0,  11,    0,   0},
      { 0, 0, 0, 0, 0,  0, -1,  0, 3, 0, 0, 0, 2,  118,   0,    0, -52},
      { 0, 0, 0, 0, 0,  0, -1,  0, 3, 0, 0, 0, 1,    0,  -5,   -3,   0},
      { 0, 0, 0, 0, 0,  0,  1,  0,-3, 0, 0, 0, 0,  -28,  36,    0,   0},
      { 0, 0, 0, 0, 0,  2, -4,  0, 0, 0, 0, 0, 0,    5,  -5,    0,   0},
      { 0, 0, 0, 0, 0, -2,  4,  0, 0, 0, 0, 0, 1,   14, -59,  -31,  -8},

   /* 341-350 */
      { 0, 1,-1, 1, 0, -2,  3,  0, 0, 0, 0, 0, 0,    0,   9,    5,   1},
      { 0, 0, 0, 0, 0, -2,  4,  0, 0, 0, 0, 0, 2, -458,   0,    0, 198},
      { 0, 0, 0, 0, 0, -6,  9,  0, 0, 0, 0, 0, 2,    0, -45,  -20,   0},
      { 0, 0, 0, 0, 0, -6,  9,  0, 0, 0, 0, 0, 1,    9,   0,    0,  -5},
      { 0, 0, 0, 0, 0,  6, -9,  0, 0, 0, 0, 0, 0,    0,  -3,    0,   0},
      { 0, 0, 0, 1, 0,  0,  1,  0,-2, 0, 0, 0, 0,    0,  -4,   -2,  -1},
      { 0, 2,-2, 1, 0, -2,  2,  0, 0, 0, 0, 0, 0,   11,   0,    0,  -6},
      { 0, 0, 0, 0, 0,  0, -4,  6, 0, 0, 0, 0, 2,    6,   0,    0,  -2},
      { 0, 0, 0, 0, 0,  0,  4, -6, 0, 0, 0, 0, 0,  -16,  23,    0,   0},
      { 0, 0, 0, 1, 0,  3, -4,  0, 0, 0, 0, 0, 0,    0,  -4,   -2,   0},

   /* 351-360 */
      { 0, 0, 0, 0, 0,  0, -1,  0, 2, 0, 0, 0, 2,   -5,   0,    0,   2},
      { 0, 0, 0, 0, 0,  0,  1,  0,-2, 0, 0, 0, 0, -166, 269,    0,   0},
      { 0, 0, 0, 1, 0,  0,  1,  0,-1, 0, 0, 0, 0,   15,   0,    0,  -8},
      { 0, 0, 0, 0, 0, -5,  9,  0, 0, 0, 0, 0, 2,   10,   0,    0,  -4},
      { 0, 0, 0, 0, 0,  0,  3, -4, 0, 0, 0, 0, 0,  -78,  45,    0,   0},
      { 0, 0, 0, 0, 0, -3,  4,  0, 0, 0, 0, 0, 2,    0,  -5,   -2,   0},
      { 0, 0, 0, 0, 0, -3,  4,  0, 0, 0, 0, 0, 1,    7,   0,    0,  -4},
      { 0, 0, 0, 0, 0,  3, -4,  0, 0, 0, 0, 0, 0,   -5, 328,    0,   0},
      { 0, 0, 0, 0, 0,  3, -4,  0, 0, 0, 0, 0, 1,    3,   0,    0,  -2},
      { 0, 0, 0, 1, 0,  0,  2, -2, 0, 0, 0, 0, 0,    5,   0,    0,  -2},

   /* 361-370 */
      { 0, 0, 0, 1, 0,  0, -1,  0, 2, 0, 0, 0, 0,    0,   3,    1,   0},
      { 0, 0, 0, 0, 0,  0,  1,  0, 0,-3, 0, 0, 0,   -3,   0,    0,   0},
      { 0, 0, 0, 0, 0,  0,  1,  0, 1,-5, 0, 0, 0,   -3,   0,    0,   0},
      { 0, 0, 0, 0, 0,  0, -1,  0, 1, 0, 0, 0, 1,    0,  -4,   -2,   0},
      { 0, 0, 0, 0, 0,  0,  1,  0,-1, 0, 0, 0, 0,-1223, -26,    0,   0},
      { 0, 0, 0, 0, 0,  0,  1,  0,-1, 0, 0, 0, 1,    0,   7,    3,   0},
      { 0, 0, 0, 0, 0,  0,  1,  0,-3, 5, 0, 0, 0,    3,   0,    0,   0},
      { 0, 0, 0, 1, 0, -3,  4,  0, 0, 0, 0, 0, 0,    0,   3,    2,   0},
      { 0, 0, 0, 0, 0,  0,  1,  0, 0,-2, 0, 0, 0,   -6,  20,    0,   0},
      { 0, 0, 0, 0, 0,  0,  2, -2, 0, 0, 0, 0, 0, -368,   0,    0,   0},

   /* 371-380 */
      { 0, 0, 0, 0, 0,  0,  1,  0, 0,-1, 0, 0, 0,  -75,   0,    0,   0},
      { 0, 0, 0, 1, 0,  0, -1,  0, 1, 0, 0, 0, 0,   11,   0,    0,  -6},
      { 0, 0, 0, 1, 0,  0, -2,  2, 0, 0, 0, 0, 0,    3,   0,    0,  -2},
      { 0, 0, 0, 0, 0, -8, 14,  0, 0, 0, 0, 0, 2,   -3,   0,    0,   1},
      { 0, 0, 0, 0, 0,  0,  1,  0, 2,-5, 0, 0, 0,  -13, -30,    0,   0},
      { 0, 0, 0, 0, 0,  0,  5, -8, 3, 0, 0, 0, 0,   21,   3,    0,   0},
      { 0, 0, 0, 0, 0,  0,  5, -8, 3, 0, 0, 0, 2,   -3,   0,    0,   1},
      { 0, 0, 0, 0, 0,  0, -1,  0, 0, 0, 0, 0, 1,   -4,   0,    0,   2},
      { 0, 0, 0, 0, 0,  0,  1,  0, 0, 0, 0, 0, 0,    8, -27,    0,   0},
      { 0, 0, 0, 0, 0,  0,  3, -8, 3, 0, 0, 0, 0,  -19, -11,    0,   0},

   /* 381-390 */
      { 0, 0, 0, 0, 0,  0, -3,  8,-3, 0, 0, 0, 2,   -4,   0,    0,   2},
      { 0, 0, 0, 0, 0,  0,  1,  0,-2, 5, 0, 0, 2,    0,   5,    2,   0},
      { 0, 0, 0, 0, 0, -8, 12,  0, 0, 0, 0, 0, 2,   -6,   0,    0,   2},
      { 0, 0, 0, 0, 0, -8, 12,  0, 0, 0, 0, 0, 0,   -8,   0,    0,   0},
      { 0, 0, 0, 0, 0,  0,  1,  0, 1,-2, 0, 0, 0,   -1,   0,    0,   0},
      { 0, 0, 0, 0, 0,  0,  1,  0, 0, 1, 0, 0, 2,  -14,   0,    0,   6},
      { 0, 0, 0, 0, 0,  0,  0,  2, 0, 0, 0, 0, 0,    6,   0,    0,   0},
      { 0, 0, 0, 0, 0,  0,  0,  2, 0, 0, 0, 0, 2,  -74,   0,    0,  32},
      { 0, 0, 0, 0, 0,  0,  1,  0, 0, 2, 0, 0, 2,    0,  -3,   -1,   0},
      { 0, 2,-2, 1, 0, -5,  5,  0, 0, 0, 0, 0, 0,    4,   0,    0,  -2},

   /* 391-400 */
      { 0, 0, 0, 0, 0,  0,  1,  0, 1, 0, 0, 0, 0,    8,  11,    0,   0},
      { 0, 0, 0, 0, 0,  0,  1,  0, 1, 0, 0, 0, 1,    0,   3,    2,   0},
      { 0, 0, 0, 0, 0,  0,  1,  0, 1, 0, 0, 0, 2, -262,   0,    0, 114},
      { 0, 0, 0, 0, 0,  3, -6,  0, 0, 0, 0, 0, 0,    0,  -4,    0,   0},
      { 0, 0, 0, 0, 0, -3,  6,  0, 0, 0, 0, 0, 1,   -7,   0,    0,   4},
      { 0, 0, 0, 0, 0, -3,  6,  0, 0, 0, 0, 0, 2,    0, -27,  -12,   0},
      { 0, 0, 0, 0, 0,  0, -1,  4, 0, 0, 0, 0, 2,  -19,  -8,   -4,   8},
      { 0, 0, 0, 0, 0, -5,  7,  0, 0, 0, 0, 0, 2,  202,   0,    0, -87},
      { 0, 0, 0, 0, 0, -5,  7,  0, 0, 0, 0, 0, 1,   -8,  35,   19,   5},
      { 0, 1,-1, 1, 0, -5,  6,  0, 0, 0, 0, 0, 0,    0,   4,    2,   0},

   /* 401-410 */
      { 0, 0, 0, 0, 0,  5, -7,  0, 0, 0, 0, 0, 0,   16,  -5,    0,   0},
      { 0, 2,-2, 1, 0,  0, -1,  0, 1, 0, 0, 0, 0,    5,   0,    0,  -3},
      { 0, 0, 0, 0, 0,  0, -1,  0, 1, 0, 0, 0, 0,    0,  -3,    0,   0},
      { 0, 0, 0, 0,-1,  0,  3,  0, 0, 0, 0, 0, 2,    1,   0,    0,   0},
      { 0, 0, 0, 0, 0,  0,  1,  0, 2, 0, 0, 0, 2,  -35, -48,  -21,  15},
      { 0, 0, 0, 0, 0,  0, -2,  6, 0, 0, 0, 0, 2,   -3,  -5,   -2,   1},
      { 0, 0, 0, 1, 0,  2, -2,  0, 0, 0, 0, 0, 0,    6,   0,    0,  -3},
      { 0, 0, 0, 0, 0,  0, -6,  9, 0, 0, 0, 0, 2,    3,   0,    0,  -1},
      { 0, 0, 0, 0, 0,  0,  6, -9, 0, 0, 0, 0, 0,    0,  -5,    0,   0},
      { 0, 0, 0, 0, 0, -2,  2,  0, 0, 0, 0, 0, 1,   12,  55,   29,  -6},

   /* 411-420 */
      { 0, 1,-1, 1, 0, -2,  1,  0, 0, 0, 0, 0, 0,    0,   5,    3,   0},
      { 0, 0, 0, 0, 0,  2, -2,  0, 0, 0, 0, 0, 0, -598,   0,    0,   0},
      { 0, 0, 0, 0, 0,  2, -2,  0, 0, 0, 0, 0, 1,   -3, -13,   -7,   1},
      { 0, 0, 0, 0, 0,  0,  1,  0, 3, 0, 0, 0, 2,   -5,  -7,   -3,   2},
      { 0, 0, 0, 0, 0,  0, -5,  7, 0, 0, 0, 0, 2,    3,   0,    0,  -1},
      { 0, 0, 0, 0, 0,  0,  5, -7, 0, 0, 0, 0, 0,    5,  -7,    0,   0},
      { 0, 0, 0, 1, 0, -2,  2,  0, 0, 0, 0, 0, 0,    4,   0,    0,  -2},
      { 0, 0, 0, 0, 0,  0,  4, -5, 0, 0, 0, 0, 0,   16,  -6,    0,   0},
      { 0, 0, 0, 0, 0,  1, -3,  0, 0, 0, 0, 0, 0,    8,  -3,    0,   0},
      { 0, 0, 0, 0, 0, -1,  3,  0, 0, 0, 0, 0, 1,    8, -31,  -16,  -4},

   /* 421-430 */
      { 0, 1,-1, 1, 0, -1,  2,  0, 0, 0, 0, 0, 0,    0,   3,    1,   0},
      { 0, 0, 0, 0, 0, -1,  3,  0, 0, 0, 0, 0, 2,  113,   0,    0, -49},
      { 0, 0, 0, 0, 0, -7, 10,  0, 0, 0, 0, 0, 2,    0, -24,  -10,   0},
      { 0, 0, 0, 0, 0, -7, 10,  0, 0, 0, 0, 0, 1,    4,   0,    0,  -2},
      { 0, 0, 0, 0, 0,  0,  3, -3, 0, 0, 0, 0, 0,   27,   0,    0,   0},
      { 0, 0, 0, 0, 0, -4,  8,  0, 0, 0, 0, 0, 2,   -3,   0,    0,   1},
      { 0, 0, 0, 0, 0, -4,  5,  0, 0, 0, 0, 0, 2,    0,  -4,   -2,   0},
      { 0, 0, 0, 0, 0, -4,  5,  0, 0, 0, 0, 0, 1,    5,   0,    0,  -2},
      { 0, 0, 0, 0, 0,  4, -5,  0, 0, 0, 0, 0, 0,    0,  -3,    0,   0},
      { 0, 0, 0, 0, 0,  0,  1,  1, 0, 0, 0, 0, 2,  -13,   0,    0,   6},

   /* 431-440 */
      { 0, 0, 0, 0, 0,  0, -2,  0, 5, 0, 0, 0, 2,    5,   0,    0,  -2},
      { 0, 0, 0, 0, 0,  0,  0,  3, 0, 0, 0, 0, 2,  -18, -10,   -4,   8},
      { 0, 0, 0, 0, 0,  1,  0,  0, 0, 0, 0, 0, 0,   -4, -28,    0,   0},
      { 0, 0, 0, 0, 0,  1,  0,  0, 0, 0, 0, 0, 2,   -5,   6,    3,   2},
      { 0, 0, 0, 0, 0, -9, 13,  0, 0, 0, 0, 0, 2,   -3,   0,    0,   1},
      { 0, 0, 0, 0, 0,  0, -1,  5, 0, 0, 0, 0, 2,   -5,  -9,   -4,   2},
      { 0, 0, 0, 0, 0,  0, -2,  0, 4, 0, 0, 0, 2,   17,   0,    0,  -7},
      { 0, 0, 0, 0, 0,  0,  2,  0,-4, 0, 0, 0, 0,   11,   4,    0,   0},
      { 0, 0, 0, 0, 0,  0, -2,  7, 0, 0, 0, 0, 2,    0,  -6,   -2,   0},
      { 0, 0, 0, 0, 0,  0,  2,  0,-3, 0, 0, 0, 0,   83,  15,    0,   0},

   /* 441-450 */
      { 0, 0, 0, 0, 0, -2,  5,  0, 0, 0, 0, 0, 1,   -4,   0,    0,   2},
      { 0, 0, 0, 0, 0, -2,  5,  0, 0, 0, 0, 0, 2,    0,-114,  -49,   0},
      { 0, 0, 0, 0, 0, -6,  8,  0, 0, 0, 0, 0, 2,  117,   0,    0, -51},
      { 0, 0, 0, 0, 0, -6,  8,  0, 0, 0, 0, 0, 1,   -5,  19,   10,   2},
      { 0, 0, 0, 0, 0,  6, -8,  0, 0, 0, 0, 0, 0,   -3,   0,    0,   0},
      { 0, 0, 0, 1, 0,  0,  2,  0,-2, 0, 0, 0, 0,   -3,   0,    0,   2},
      { 0, 0, 0, 0, 0,  0, -3,  9, 0, 0, 0, 0, 2,    0,  -3,   -1,   0},
      { 0, 0, 0, 0, 0,  0,  5, -6, 0, 0, 0, 0, 0,    3,   0,    0,   0},
      { 0, 0, 0, 0, 0,  0,  5, -6, 0, 0, 0, 0, 2,    0,  -6,   -2,   0},
      { 0, 0, 0, 0, 0,  0,  2,  0,-2, 0, 0, 0, 0,  393,   3,    0,   0},

   /* 451-460 */
      { 0, 0, 0, 0, 0,  0,  2,  0,-2, 0, 0, 0, 1,   -4,  21,   11,   2},
      { 0, 0, 0, 0, 0,  0,  2,  0,-2, 0, 0, 0, 2,   -6,   0,   -1,   3},
      { 0, 0, 0, 0, 0, -5, 10,  0, 0, 0, 0, 0, 2,   -3,   8,    4,   1},
      { 0, 0, 0, 0, 0,  0,  4, -4, 0, 0, 0, 0, 0,    8,   0,    0,   0},
      { 0, 0, 0, 0, 0,  0,  4, -4, 0, 0, 0, 0, 2,   18, -29,  -13,  -8},
      { 0, 0, 0, 0, 0, -3,  3,  0, 0, 0, 0, 0, 1,    8,  34,   18,  -4},
      { 0, 0, 0, 0, 0,  3, -3,  0, 0, 0, 0, 0, 0,   89,   0,    0,   0},
      { 0, 0, 0, 0, 0,  3, -3,  0, 0, 0, 0, 0, 1,    3,  12,    6,  -1},
      { 0, 0, 0, 0, 0,  3, -3,  0, 0, 0, 0, 0, 2,   54, -15,   -7, -24},
      { 0, 0, 0, 0, 0,  0,  2,  0, 0,-3, 0, 0, 0,    0,   3,    0,   0},

   /* 461-470 */
      { 0, 0, 0, 0, 0,  0, -5, 13, 0, 0, 0, 0, 2,    3,   0,    0,  -1},
      { 0, 0, 0, 0, 0,  0,  2,  0,-1, 0, 0, 0, 0,    0,  35,    0,   0},
      { 0, 0, 0, 0, 0,  0,  2,  0,-1, 0, 0, 0, 2, -154, -30,  -13,  67},
      { 0, 0, 0, 0, 0,  0,  2,  0, 0,-2, 0, 0, 0,   15,   0,    0,   0},
      { 0, 0, 0, 0, 0,  0,  2,  0, 0,-2, 0, 0, 1,    0,   4,    2,   0},
      { 0, 0, 0, 0, 0,  0,  3, -2, 0, 0, 0, 0, 0,    0,   9,    0,   0},
      { 0, 0, 0, 0, 0,  0,  3, -2, 0, 0, 0, 0, 2,   80, -71,  -31, -35},
      { 0, 0, 0, 0, 0,  0,  2,  0, 0,-1, 0, 0, 2,    0, -20,   -9,   0},
      { 0, 0, 0, 0, 0,  0, -6, 15, 0, 0, 0, 0, 2,   11,   5,    2,  -5},
      { 0, 0, 0, 0, 0, -8, 15,  0, 0, 0, 0, 0, 2,   61, -96,  -42, -27},

   /* 471-480 */
      { 0, 0, 0, 0, 0, -3,  9, -4, 0, 0, 0, 0, 2,   14,   9,    4,  -6},
      { 0, 0, 0, 0, 0,  0,  2,  0, 2,-5, 0, 0, 2,  -11,  -6,   -3,   5},
      { 0, 0, 0, 0, 0,  0, -2,  8,-1,-5, 0, 0, 2,    0,  -3,   -1,   0},
      { 0, 0, 0, 0, 0,  0,  6, -8, 3, 0, 0, 0, 2,  123,-415, -180, -53},
      { 0, 0, 0, 0, 0,  0,  2,  0, 0, 0, 0, 0, 0,    0,   0,    0, -35},
      { 0, 0, 0, 0, 0,  0,  2,  0, 0, 0, 0, 0, 0,   -5,   0,    0,   0},
      { 0, 0, 0, 0, 0,  0,  2,  0, 0, 0, 0, 0, 1,    7, -32,  -17,  -4},
      { 0, 1,-1, 1, 0,  0,  1,  0, 0, 0, 0, 0, 0,    0,  -9,   -5,   0},
      { 0, 0, 0, 0, 0,  0,  2,  0, 0, 0, 0, 0, 1,    0,  -4,    2,   0},
      { 0, 0, 0, 0, 0,  0,  2,  0, 0, 0, 0, 0, 2,  -89,   0,    0,  38},

   /* 481-490 */
      { 0, 0, 0, 0, 0,  0, -6, 16,-4,-5, 0, 0, 2,    0, -86,  -19,  -6},
      { 0, 0, 0, 0, 0,  0, -2,  8,-3, 0, 0, 0, 2,    0,   0,  -19,   6},
      { 0, 0, 0, 0, 0,  0, -2,  8,-3, 0, 0, 0, 2, -123,-416, -180,  53},
      { 0, 0, 0, 0, 0,  0,  6, -8, 1, 5, 0, 0, 2,    0,  -3,   -1,   0},
      { 0, 0, 0, 0, 0,  0,  2,  0,-2, 5, 0, 0, 2,   12,  -6,   -3,  -5},
      { 0, 0, 0, 0, 0,  3, -5,  4, 0, 0, 0, 0, 2,  -13,   9,    4,   6},
      { 0, 0, 0, 0, 0, -8, 11,  0, 0, 0, 0, 0, 2,    0, -15,   -7,   0},
      { 0, 0, 0, 0, 0, -8, 11,  0, 0, 0, 0, 0, 1,    3,   0,    0,  -1},
      { 0, 0, 0, 0, 0, -8, 11,  0, 0, 0, 0, 0, 2,  -62, -97,  -42,  27},
      { 0, 0, 0, 0, 0,  0, 11,  0, 0, 0, 0, 0, 2,  -11,   5,    2,   5},

   /* 491-500 */
      { 0, 0, 0, 0, 0,  0,  2,  0, 0, 1, 0, 0, 2,    0, -19,   -8,   0},
      { 0, 0, 0, 0, 0,  3, -3,  0, 2, 0, 0, 0, 2,   -3,   0,    0,   1},
      { 0, 2,-2, 1, 0,  0,  4, -8, 3, 0, 0, 0, 0,    0,   4,    2,   0},
      { 0, 1,-1, 0, 0,  0,  1,  0, 0, 0, 0, 0, 0,    0,   3,    0,   0},
      { 0, 2,-2, 1, 0,  0, -4,  8,-3, 0, 0, 0, 0,    0,   4,    2,   0},
      { 0, 0, 0, 0, 0,  0,  1,  2, 0, 0, 0, 0, 2,  -85, -70,  -31,  37},
      { 0, 0, 0, 0, 0,  0,  2,  0, 1, 0, 0, 0, 2,  163, -12,   -5, -72},
      { 0, 0, 0, 0, 0, -3,  7,  0, 0, 0, 0, 0, 2,  -63, -16,   -7,  28},
      { 0, 0, 0, 0, 0,  0,  0,  4, 0, 0, 0, 0, 2,  -21, -32,  -14,   9},
      { 0, 0, 0, 0, 0, -5,  6,  0, 0, 0, 0, 0, 2,    0,  -3,   -1,   0},

   /* 501-510 */
      { 0, 0, 0, 0, 0, -5,  6,  0, 0, 0, 0, 0, 1,    3,   0,    0,  -2},
      { 0, 0, 0, 0, 0,  5, -6,  0, 0, 0, 0, 0, 0,    0,   8,    0,   0},
      { 0, 0, 0, 0, 0,  5, -6,  0, 0, 0, 0, 0, 2,    3,  10,    4,  -1},
      { 0, 0, 0, 0, 0,  0,  2,  0, 2, 0, 0, 0, 2,    3,   0,    0,  -1},
      { 0, 0, 0, 0, 0,  0, -1,  6, 0, 0, 0, 0, 2,    0,  -7,   -3,   0},
      { 0, 0, 0, 0, 0,  0,  7, -9, 0, 0, 0, 0, 2,    0,  -4,   -2,   0},
      { 0, 0, 0, 0, 0,  2, -1,  0, 0, 0, 0, 0, 0,    6,  19,    0,   0},
      { 0, 0, 0, 0, 0,  2, -1,  0, 0, 0, 0, 0, 2,    5,-173,  -75,  -2},
      { 0, 0, 0, 0, 0,  0,  6, -7, 0, 0, 0, 0, 2,    0,  -7,   -3,   0},
      { 0, 0, 0, 0, 0,  0,  5, -5, 0, 0, 0, 0, 2,    7, -12,   -5,  -3},

   /* 511-520 */
      { 0, 0, 0, 0, 0, -1,  4,  0, 0, 0, 0, 0, 1,   -3,   0,    0,   2},
      { 0, 0, 0, 0, 0, -1,  4,  0, 0, 0, 0, 0, 2,    3,  -4,   -2,  -1},
      { 0, 0, 0, 0, 0, -7,  9,  0, 0, 0, 0, 0, 2,   74,   0,    0, -32},
      { 0, 0, 0, 0, 0, -7,  9,  0, 0, 0, 0, 0, 1,   -3,  12,    6,   2},
      { 0, 0, 0, 0, 0,  0,  4, -3, 0, 0, 0, 0, 2,   26, -14,   -6, -11},
      { 0, 0, 0, 0, 0,  0,  3, -1, 0, 0, 0, 0, 2,   19,   0,    0,  -8},
      { 0, 0, 0, 0, 0, -4,  4,  0, 0, 0, 0, 0, 1,    6,  24,   13,  -3},
      { 0, 0, 0, 0, 0,  4, -4,  0, 0, 0, 0, 0, 0,   83,   0,    0,   0},
      { 0, 0, 0, 0, 0,  4, -4,  0, 0, 0, 0, 0, 1,    0, -10,   -5,   0},
      { 0, 0, 0, 0, 0,  4, -4,  0, 0, 0, 0, 0, 2,   11,  -3,   -1,  -5},

   /* 521-530 */
      { 0, 0, 0, 0, 0,  0,  2,  1, 0, 0, 0, 0, 2,    3,   0,    1,  -1},
      { 0, 0, 0, 0, 0,  0, -3,  0, 5, 0, 0, 0, 2,    3,   0,    0,  -1},
      { 0, 0, 0, 0, 0,  1,  1,  0, 0, 0, 0, 0, 0,   -4,   0,    0,   0},
      { 0, 0, 0, 0, 0,  1,  1,  0, 0, 0, 0, 0, 1,    5, -23,  -12,  -3},
      { 0, 0, 0, 0, 0,  1,  1,  0, 0, 0, 0, 0, 2, -339,   0,    0, 147},
      { 0, 0, 0, 0, 0, -9, 12,  0, 0, 0, 0, 0, 2,    0, -10,   -5,   0},
      { 0, 0, 0, 0, 0,  0,  3,  0,-4, 0, 0, 0, 0,    5,   0,    0,   0},
      { 0, 2,-2, 1, 0,  1, -1,  0, 0, 0, 0, 0, 0,    3,   0,    0,  -1},
      { 0, 0, 0, 0, 0,  0,  7, -8, 0, 0, 0, 0, 2,    0,  -4,   -2,   0},
      { 0, 0, 0, 0, 0,  0,  3,  0,-3, 0, 0, 0, 0,   18,  -3,    0,   0},

   /* 531-540 */
      { 0, 0, 0, 0, 0,  0,  3,  0,-3, 0, 0, 0, 2,    9, -11,   -5,  -4},
      { 0, 0, 0, 0, 0, -2,  6,  0, 0, 0, 0, 0, 2,   -8,   0,    0,   4},
      { 0, 0, 0, 0, 0, -6,  7,  0, 0, 0, 0, 0, 1,    3,   0,    0,  -1},
      { 0, 0, 0, 0, 0,  6, -7,  0, 0, 0, 0, 0, 0,    0,   9,    0,   0},
      { 0, 0, 0, 0, 0,  0,  6, -6, 0, 0, 0, 0, 2,    6,  -9,   -4,  -2},
      { 0, 0, 0, 0, 0,  0,  3,  0,-2, 0, 0, 0, 0,   -4, -12,    0,   0},
      { 0, 0, 0, 0, 0,  0,  3,  0,-2, 0, 0, 0, 2,   67, -91,  -39, -29},
      { 0, 0, 0, 0, 0,  0,  5, -4, 0, 0, 0, 0, 2,   30, -18,   -8, -13},
      { 0, 0, 0, 0, 0,  3, -2,  0, 0, 0, 0, 0, 0,    0,   0,    0,   0},
      { 0, 0, 0, 0, 0,  3, -2,  0, 0, 0, 0, 0, 2,    0,-114,  -50,   0},

   /* 541-550 */
      { 0, 0, 0, 0, 0,  0,  3,  0,-1, 0, 0, 0, 2,    0,   0,    0,  23},
      { 0, 0, 0, 0, 0,  0,  3,  0,-1, 0, 0, 0, 2,  517,  16,    7,-224},
      { 0, 0, 0, 0, 0,  0,  3,  0, 0,-2, 0, 0, 2,    0,  -7,   -3,   0},
      { 0, 0, 0, 0, 0,  0,  4, -2, 0, 0, 0, 0, 2,  143,  -3,   -1, -62},
      { 0, 0, 0, 0, 0,  0,  3,  0, 0,-1, 0, 0, 2,   29,   0,    0, -13},
      { 0, 2,-2, 1, 0,  0,  1,  0,-1, 0, 0, 0, 0,   -4,   0,    0,   2},
      { 0, 0, 0, 0, 0, -8, 16,  0, 0, 0, 0, 0, 2,   -6,   0,    0,   3},
      { 0, 0, 0, 0, 0,  0,  3,  0, 2,-5, 0, 0, 2,    5,  12,    5,  -2},
      { 0, 0, 0, 0, 0,  0,  7, -8, 3, 0, 0, 0, 2,  -25,   0,    0,  11},
      { 0, 0, 0, 0, 0,  0, -5, 16,-4,-5, 0, 0, 2,   -3,   0,    0,   1},

   /* 551-560 */
      { 0, 0, 0, 0, 0,  0,  3,  0, 0, 0, 0, 0, 2,    0,   4,    2,   0},
      { 0, 0, 0, 0, 0,  0, -1,  8,-3, 0, 0, 0, 2,  -22,  12,    5,  10},
      { 0, 0, 0, 0, 0, -8, 10,  0, 0, 0, 0, 0, 2,   50,   0,    0, -22},
      { 0, 0, 0, 0, 0, -8, 10,  0, 0, 0, 0, 0, 1,    0,   7,    4,   0},
      { 0, 0, 0, 0, 0, -8, 10,  0, 0, 0, 0, 0, 2,    0,   3,    1,   0},
      { 0, 0, 0, 0, 0,  0,  2,  2, 0, 0, 0, 0, 2,   -4,   4,    2,   2},
      { 0, 0, 0, 0, 0,  0,  3,  0, 1, 0, 0, 0, 2,   -5, -11,   -5,   2},
      { 0, 0, 0, 0, 0, -3,  8,  0, 0, 0, 0, 0, 2,    0,   4,    2,   0},
      { 0, 0, 0, 0, 0, -5,  5,  0, 0, 0, 0, 0, 1,    4,  17,    9,  -2},
      { 0, 0, 0, 0, 0,  5, -5,  0, 0, 0, 0, 0, 0,   59,   0,    0,   0},

   /* 561-570 */
      { 0, 0, 0, 0, 0,  5, -5,  0, 0, 0, 0, 0, 1,    0,  -4,   -2,   0},
      { 0, 0, 0, 0, 0,  5, -5,  0, 0, 0, 0, 0, 2,   -8,   0,    0,   4},
      { 0, 0, 0, 0, 0,  2,  0,  0, 0, 0, 0, 0, 0,   -3,   0,    0,   0},
      { 0, 0, 0, 0, 0,  2,  0,  0, 0, 0, 0, 0, 1,    4, -15,   -8,  -2},
      { 0, 0, 0, 0, 0,  2,  0,  0, 0, 0, 0, 0, 2,  370,  -8,    0,-160},
      { 0, 0, 0, 0, 0,  0,  7, -7, 0, 0, 0, 0, 2,    0,   0,   -3,   0},
      { 0, 0, 0, 0, 0,  0,  7, -7, 0, 0, 0, 0, 2,    0,   3,    1,   0},
      { 0, 0, 0, 0, 0,  0,  6, -5, 0, 0, 0, 0, 2,   -6,   3,    1,   3},
      { 0, 0, 0, 0, 0,  7, -8,  0, 0, 0, 0, 0, 0,    0,   6,    0,   0},
      { 0, 0, 0, 0, 0,  0,  5, -3, 0, 0, 0, 0, 2,  -10,   0,    0,   4},

   /* 571-580 */
      { 0, 0, 0, 0, 0,  4, -3,  0, 0, 0, 0, 0, 2,    0,   9,    4,   0},
      { 0, 0, 0, 0, 0,  1,  2,  0, 0, 0, 0, 0, 2,    4,  17,    7,  -2},
      { 0, 0, 0, 0, 0, -9, 11,  0, 0, 0, 0, 0, 2,   34,   0,    0, -15},
      { 0, 0, 0, 0, 0, -9, 11,  0, 0, 0, 0, 0, 1,    0,   5,    3,   0},
      { 0, 0, 0, 0, 0,  0,  4,  0,-4, 0, 0, 0, 2,   -5,   0,    0,   2},
      { 0, 0, 0, 0, 0,  0,  4,  0,-3, 0, 0, 0, 2,  -37,  -7,   -3,  16},
      { 0, 0, 0, 0, 0, -6,  6,  0, 0, 0, 0, 0, 1,    3,  13,    7,  -2},
      { 0, 0, 0, 0, 0,  6, -6,  0, 0, 0, 0, 0, 0,   40,   0,    0,   0},
      { 0, 0, 0, 0, 0,  6, -6,  0, 0, 0, 0, 0, 1,    0,  -3,   -2,   0},
      { 0, 0, 0, 0, 0,  0,  4,  0,-2, 0, 0, 0, 2, -184,  -3,   -1,  80},

   /* 581-590 */
      { 0, 0, 0, 0, 0,  0,  6, -4, 0, 0, 0, 0, 2,   -3,   0,    0,   1},
      { 0, 0, 0, 0, 0,  3, -1,  0, 0, 0, 0, 0, 0,   -3,   0,    0,   0},
      { 0, 0, 0, 0, 0,  3, -1,  0, 0, 0, 0, 0, 1,    0, -10,   -6,  -1},
      { 0, 0, 0, 0, 0,  3, -1,  0, 0, 0, 0, 0, 2,   31,  -6,    0, -13},
      { 0, 0, 0, 0, 0,  0,  4,  0,-1, 0, 0, 0, 2,   -3, -32,  -14,   1},
      { 0, 0, 0, 0, 0,  0,  4,  0, 0,-2, 0, 0, 2,   -7,   0,    0,   3},
      { 0, 0, 0, 0, 0,  0,  5, -2, 0, 0, 0, 0, 2,    0,  -8,   -4,   0},
      { 0, 0, 0, 0, 0,  0,  4,  0, 0, 0, 0, 0, 0,    3,  -4,    0,   0},
      { 0, 0, 0, 0, 0,  8, -9,  0, 0, 0, 0, 0, 0,    0,   4,    0,   0},
      { 0, 0, 0, 0, 0,  5, -4,  0, 0, 0, 0, 0, 2,    0,   3,    1,   0},

   /* 591-600 */
      { 0, 0, 0, 0, 0,  2,  1,  0, 0, 0, 0, 0, 2,   19, -23,  -10,   2},
      { 0, 0, 0, 0, 0,  2,  1,  0, 0, 0, 0, 0, 1,    0,   0,    0, -10},
      { 0, 0, 0, 0, 0,  2,  1,  0, 0, 0, 0, 0, 1,    0,   3,    2,   0},
      { 0, 0, 0, 0, 0, -7,  7,  0, 0, 0, 0, 0, 1,    0,   9,    5,  -1},
      { 0, 0, 0, 0, 0,  7, -7,  0, 0, 0, 0, 0, 0,   28,   0,    0,   0},
      { 0, 0, 0, 0, 0,  4, -2,  0, 0, 0, 0, 0, 1,    0,  -7,   -4,   0},
      { 0, 0, 0, 0, 0,  4, -2,  0, 0, 0, 0, 0, 2,    8,  -4,    0,  -4},
      { 0, 0, 0, 0, 0,  4, -2,  0, 0, 0, 0, 0, 0,    0,   0,   -2,   0},
      { 0, 0, 0, 0, 0,  4, -2,  0, 0, 0, 0, 0, 0,    0,   3,    0,   0},
      { 0, 0, 0, 0, 0,  0,  5,  0,-4, 0, 0, 0, 2,   -3,   0,    0,   1},

   /* 601-610 */
      { 0, 0, 0, 0, 0,  0,  5,  0,-3, 0, 0, 0, 2,   -9,   0,    1,   4},
      { 0, 0, 0, 0, 0,  0,  5,  0,-2, 0, 0, 0, 2,    3,  12,    5,  -1},
      { 0, 0, 0, 0, 0,  3,  0,  0, 0, 0, 0, 0, 2,   17,  -3,   -1,   0},
      { 0, 0, 0, 0, 0, -8,  8,  0, 0, 0, 0, 0, 1,    0,   7,    4,   0},
      { 0, 0, 0, 0, 0,  8, -8,  0, 0, 0, 0, 0, 0,   19,   0,    0,   0},
      { 0, 0, 0, 0, 0,  5, -3,  0, 0, 0, 0, 0, 1,    0,  -5,   -3,   0},
      { 0, 0, 0, 0, 0,  5, -3,  0, 0, 0, 0, 0, 2,   14,  -3,    0,  -1},
      { 0, 0, 0, 0, 0, -9,  9,  0, 0, 0, 0, 0, 1,    0,   0,   -1,   0},
      { 0, 0, 0, 0, 0, -9,  9,  0, 0, 0, 0, 0, 1,    0,   0,    0,  -5},
      { 0, 0, 0, 0, 0, -9,  9,  0, 0, 0, 0, 0, 1,    0,   5,    3,   0},

   /* 611-620 */
      { 0, 0, 0, 0, 0,  9, -9,  0, 0, 0, 0, 0, 0,   13,   0,    0,   0},
      { 0, 0, 0, 0, 0,  6, -4,  0, 0, 0, 0, 0, 1,    0,  -3,   -2,   0},
      { 0, 0, 0, 0, 0,  0,  6,  0, 0, 0, 0, 0, 2,    2,   9,    4,   3},
      { 0, 0, 0, 0, 0,  0,  6,  0, 0, 0, 0, 0, 0,    0,   0,    0,  -4},
      { 0, 0, 0, 0, 0,  0,  6,  0, 0, 0, 0, 0, 0,    8,   0,    0,   0},
      { 0, 0, 0, 0, 0,  0,  6,  0, 0, 0, 0, 0, 1,    0,   4,    2,   0},
      { 0, 0, 0, 0, 0,  0,  6,  0, 0, 0, 0, 0, 2,    6,   0,    0,  -3},
      { 0, 0, 0, 0, 0,  0,  6,  0, 0, 0, 0, 0, 0,    6,   0,    0,   0},
      { 0, 0, 0, 0, 0,  0,  6,  0, 0, 0, 0, 0, 1,    0,   3,    1,   0},
      { 0, 0, 0, 0, 0,  0,  6,  0, 0, 0, 0, 0, 2,    5,   0,    0,  -2},

   /* 621-630 */
      { 0, 0, 0, 0, 0,  0,  0,  0, 0, 0, 0, 0, 2,    3,   0,    0,  -1},
      { 1, 0,-2, 0, 0,  0,  2,  0,-2, 0, 0, 0, 0,   -3,   0,    0,   0},
      { 1, 0,-2, 0, 0,  2, -2,  0, 0, 0, 0, 0, 0,    6,   0,    0,   0},
      { 1, 0,-2, 0, 0,  0,  1,  0,-1, 0, 0, 0, 0,    7,   0,    0,   0},
      { 1, 0,-2, 0, 0,  1, -1,  0, 0, 0, 0, 0, 0,   -4,   0,    0,   0},
      {-1, 0, 0, 0, 0,  3, -3,  0, 0, 0, 0, 0, 0,    4,   0,    0,   0},
      {-1, 0, 0, 0, 0,  0,  2,  0,-2, 0, 0, 0, 0,    6,   0,    0,   0},
      {-1, 0, 2, 0, 0,  0,  4, -8, 3, 0, 0, 0, 0,    0,  -4,    0,   0},
      { 1, 0,-2, 0, 0,  0,  4, -8, 3, 0, 0, 0, 0,    0,  -4,    0,   0},
      {-2, 0, 2, 0, 0,  0,  4, -8, 3, 0, 0, 0, 0,    5,   0,    0,   0},

   /* 631-640 */
      {-1, 0, 0, 0, 0,  0,  2,  0,-3, 0, 0, 0, 0,   -3,   0,    0,   0},
      {-1, 0, 0, 0, 0,  0,  1,  0,-1, 0, 0, 0, 0,    4,   0,    0,   0},
      {-1, 0, 0, 0, 0,  1, -1,  0, 0, 0, 0, 0, 0,   -5,   0,    0,   0},
      {-1, 0, 2, 0, 0,  2, -2,  0, 0, 0, 0, 0, 0,    4,   0,    0,   0},
      { 1,-1, 1, 0, 0,  0,  1,  0, 0, 0, 0, 0, 0,    0,   3,    0,   0},
      {-1, 0, 2, 0, 0,  0,  2,  0,-3, 0, 0, 0, 0,   13,   0,    0,   0},
      {-2, 0, 0, 0, 0,  0,  2,  0,-3, 0, 0, 0, 0,   21,  11,    0,   0},
      { 1, 0, 0, 0, 0,  0,  4, -8, 3, 0, 0, 0, 0,    0,  -5,    0,   0},
      {-1, 1,-1, 1, 0,  0, -1,  0, 0, 0, 0, 0, 0,    0,  -5,   -2,   0},
      { 1, 1,-1, 1, 0,  0, -1,  0, 0, 0, 0, 0, 0,    0,   5,    3,   0},

   /* 641-650 */
      {-1, 0, 0, 0, 0,  0,  4, -8, 3, 0, 0, 0, 0,    0,  -5,    0,   0},
      {-1, 0, 2, 1, 0,  0,  2,  0,-2, 0, 0, 0, 0,   -3,   0,    0,   2},
      { 0, 0, 0, 0, 0,  0,  2,  0,-2, 0, 0, 0, 0,   20,  10,    0,   0},
      {-1, 0, 2, 0, 0,  0,  2,  0,-2, 0, 0, 0, 0,  -34,   0,    0,   0},
      {-1, 0, 2, 0, 0,  3, -3,  0, 0, 0, 0, 0, 0,  -19,   0,    0,   0},
      { 1, 0,-2, 1, 0,  0, -2,  0, 2, 0, 0, 0, 0,    3,   0,    0,  -2},
      { 1, 2,-2, 2, 0, -3,  3,  0, 0, 0, 0, 0, 0,   -3,   0,    0,   1},
      { 1, 2,-2, 2, 0,  0, -2,  0, 2, 0, 0, 0, 0,   -6,   0,    0,   3},
      { 1, 0, 0, 0, 0,  1, -1,  0, 0, 0, 0, 0, 0,   -4,   0,    0,   0},
      { 1, 0, 0, 0, 0,  0,  1,  0,-1, 0, 0, 0, 0,    3,   0,    0,   0},

   /* 651-660 */
      { 0, 0,-2, 0, 0,  2, -2,  0, 0, 0, 0, 0, 0,    3,   0,    0,   0},
      { 0, 0,-2, 0, 0,  0,  1,  0,-1, 0, 0, 0, 0,    4,   0,    0,   0},
      { 0, 2, 0, 2, 0, -2,  2,  0, 0, 0, 0, 0, 0,    3,   0,    0,  -1},
      { 0, 2, 0, 2, 0,  0, -1,  0, 1, 0, 0, 0, 0,    6,   0,    0,  -3},
      { 0, 2, 0, 2, 0, -1,  1,  0, 0, 0, 0, 0, 0,   -8,   0,    0,   3},
      { 0, 2, 0, 2, 0, -2,  3,  0, 0, 0, 0, 0, 0,    0,   3,    1,   0},
      { 0, 0, 2, 0, 0,  0,  2,  0,-2, 0, 0, 0, 0,   -3,   0,    0,   0},
      { 0, 1, 1, 2, 0,  0,  1,  0, 0, 0, 0, 0, 0,    0,  -3,   -2,   0},
      { 1, 2, 0, 2, 0,  0,  1,  0, 0, 0, 0, 0, 0,  126, -63,  -27, -55},
      {-1, 2, 0, 2, 0, 10, -3,  0, 0, 0, 0, 0, 0,   -5,   0,    1,   2},

   /* 661-670 */
      { 0, 1, 1, 1, 0,  0,  1,  0, 0, 0, 0, 0, 0,   -3,  28,   15,   2},
      { 1, 2, 0, 2, 0,  0,  1,  0, 0, 0, 0, 0, 0,    5,   0,    1,  -2},
      { 0, 2, 0, 2, 0,  0,  4, -8, 3, 0, 0, 0, 0,    0,   9,    4,   1},
      { 0, 2, 0, 2, 0,  0, -4,  8,-3, 0, 0, 0, 0,    0,   9,    4,  -1},
      {-1, 2, 0, 2, 0,  0, -4,  8,-3, 0, 0, 0, 0, -126, -63,  -27,  55},
      { 2, 2,-2, 2, 0,  0, -2,  0, 3, 0, 0, 0, 0,    3,   0,    0,  -1},
      { 1, 2, 0, 1, 0,  0, -2,  0, 3, 0, 0, 0, 0,   21, -11,   -6, -11},
      { 0, 1, 1, 0, 0,  0,  1,  0, 0, 0, 0, 0, 0,    0,  -4,    0,   0},
      {-1, 2, 0, 1, 0,  0,  1,  0, 0, 0, 0, 0, 0,  -21, -11,   -6,  11},
      {-2, 2, 2, 2, 0,  0,  2,  0,-2, 0, 0, 0, 0,   -3,   0,    0,   1},

   /* 671-680 */
      { 0, 2, 0, 2, 0,  2, -3,  0, 0, 0, 0, 0, 0,    0,   3,    1,   0},
      { 0, 2, 0, 2, 0,  1, -1,  0, 0, 0, 0, 0, 0,    8,   0,    0,  -4},
      { 0, 2, 0, 2, 0,  0,  1,  0,-1, 0, 0, 0, 0,   -6,   0,    0,   3},
      { 0, 2, 0, 2, 0,  2, -2,  0, 0, 0, 0, 0, 0,   -3,   0,    0,   1},
      {-1, 2, 2, 2, 0,  0, -1,  0, 1, 0, 0, 0, 0,    3,   0,    0,  -1},
      { 1, 2, 0, 2, 0, -1,  1,  0, 0, 0, 0, 0, 0,   -3,   0,    0,   1},
      {-1, 2, 2, 2, 0,  0,  2,  0,-3, 0, 0, 0, 0,   -5,   0,    0,   2},
      { 2, 2, 0, 2, 0,  0,  2,  0,-3, 0, 0, 0, 0,   24, -12,   -5, -11},
      { 1, 2, 0, 2, 0,  0, -4,  8,-3, 0, 0, 0, 0,    0,   3,    1,   0},
      { 1, 2, 0, 2, 0,  0,  4, -8, 3, 0, 0, 0, 0,    0,   3,    1,   0},

   /* 681-687 */
      { 1, 1, 1, 1, 0,  0,  1,  0, 0, 0, 0, 0, 0,    0,   3,    2,   0},
      { 0, 2, 0, 2, 0,  0,  1,  0, 0, 0, 0, 0, 0,  -24, -12,   -5,  10},
      { 2, 2, 0, 1, 0,  0,  1,  0, 0, 0, 0, 0, 0,    4,   0,   -1,  -2},
      {-1, 2, 2, 2, 0,  0,  2,  0,-2, 0, 0, 0, 0,   13,   0,    0,  -6},
      {-1, 2, 2, 2, 0,  3, -3,  0, 0, 0, 0, 0, 0,    7,   0,    0,  -3},
      { 1, 2, 0, 2, 0,  1, -1,  0, 0, 0, 0, 0, 0,    3,   0,    0,  -1},
      { 0, 2, 2, 2, 0,  0,  2,  0,-2, 0, 0, 0, 0,    3,   0,    0,  -1}
   };

/* Number of terms in the planetary nutation model */
   const int NPL = (int) (sizeof xpl / sizeof xpl[0]);

/*--------------------------------------------------------------------*/

/* Interval between fundamental date J2000.0 and given date (JC). */
   t = ((date1 - DJ00) + date2) / DJC;

/* ------------------- */
/* LUNI-SOLAR NUTATION */
/* ------------------- */

/* Fundamental (Delaunay) arguments */

/* Mean anomaly of the Moon (IERS 2003). */
   el = eraFal03(t);

/* Mean anomaly of the Sun (MHB2000). */
   elp = fmod(1287104.79305  +
            t * (129596581.0481  +
            t * (-0.5532  +
            t * (0.000136  +
            t * (-0.00001149)))), TURNAS) * DAS2R;

/* Mean longitude of the Moon minus that of the ascending node */
/* (IERS 2003. */
   f = eraFaf03(t);

/* Mean elongation of the Moon from the Sun (MHB2000). */
   d = fmod(1072260.70369  +
          t * (1602961601.2090  +
          t * (-6.3706  +
          t * (0.006593  +
          t * (-0.00003169)))), TURNAS) * DAS2R;

/* Mean longitude of the ascending node of the Moon (IERS 2003). */
   om = eraFaom03(t);

/* Initialize the nutation values. */
   dp = 0.0;
   de = 0.0;

/* Summation of luni-solar nutation series (in reverse order). */
   for (i = NLS-1; i >= 0; i--) {

   /* Argument and functions. */
      arg = fmod((double)xls[i].nl  * el +
                 (double)xls[i].nlp * elp +
                 (double)xls[i].nf  * f +
                 (double)xls[i].nd  * d +
                 (double)xls[i].nom * om, D2PI);
      sarg = sin(arg);
      carg = cos(arg);

   /* Term. */
      dp += (xls[i].sp + xls[i].spt * t) * sarg + xls[i].cp * carg;
      de += (xls[i].ce + xls[i].cet * t) * carg + xls[i].se * sarg;
   }

/* Convert from 0.1 microarcsec units to radians. */
   dpsils = dp * U2R;
   depsls = de * U2R;

/* ------------------ */
/* PLANETARY NUTATION */
/* ------------------ */

/* n.b.  The MHB2000 code computes the luni-solar and planetary nutation */
/* in different functions, using slightly different Delaunay */
/* arguments in the two cases.  This behaviour is faithfully */
/* reproduced here.  Use of the IERS 2003 expressions for both */
/* cases leads to negligible changes, well below */
/* 0.1 microarcsecond. */

/* Mean anomaly of the Moon (MHB2000). */
   al = fmod(2.35555598 + 8328.6914269554 * t, D2PI);

/* Mean longitude of the Moon minus that of the ascending node */
/*(MHB2000). */
   af = fmod(1.627905234 + 8433.466158131 * t, D2PI);

/* Mean elongation of the Moon from the Sun (MHB2000). */
   ad = fmod(5.198466741 + 7771.3771468121 * t, D2PI);

/* Mean longitude of the ascending node of the Moon (MHB2000). */
   aom = fmod(2.18243920 - 33.757045 * t, D2PI);

/* General accumulated precession in longitude (IERS 2003). */
   apa = eraFapa03(t);

/* Planetary longitudes, Mercury through Uranus (IERS 2003). */
   alme = eraFame03(t);
   alve = eraFave03(t);
   alea = eraFae03(t);
   alma = eraFama03(t);
   alju = eraFaju03(t);
   alsa = eraFasa03(t);
   alur = eraFaur03(t);

/* Neptune longitude (MHB2000). */
   alne = fmod(5.321159000 + 3.8127774000 * t, D2PI);

/* Initialize the nutation values. */
   dp = 0.0;
   de = 0.0;

/* Summation of planetary nutation series (in reverse order). */
   for (i = NPL-1; i >= 0; i--) {

   /* Argument and functions. */
      arg = fmod((double)xpl[i].nl  * al   +
                 (double)xpl[i].nf  * af   +
                 (double)xpl[i].nd  * ad   +
                 (double)xpl[i].nom * aom  +
                 (double)xpl[i].nme * alme +
                 (double)xpl[i].nve * alve +
                 (double)xpl[i].nea * alea +
                 (double)xpl[i].nma * alma +
                 (double)xpl[i].nju * alju +
                 (double)xpl[i].nsa * alsa +
                 (double)xpl[i].nur * alur +
                 (double)xpl[i].nne * alne +
                 (double)xpl[i].npa * apa, D2PI);
      sarg = sin(arg);
      carg = cos(arg);

   /* Term. */
      dp += (double)xpl[i].sp * sarg + (double)xpl[i].cp * carg;
      de += (double)xpl[i].se * sarg + (double)xpl[i].ce * carg;

   }

/* Convert from 0.1 microarcsec units to radians. */
   dpsipl = dp * U2R;
   depspl = de * U2R;

/* ------- */
/* RESULTS */
/* ------- */

/* Add luni-solar and planetary components. */
   *dpsi = dpsils + dpsipl;
   *deps = depsls + depspl;

   return;

}

void eraNut00b(double date1, double date2, double *dpsi, double *deps)
/*
**  - - - - - - - - - -
**   e r a N u t 0 0 b
**  - - - - - - - - - -
**
**  Nutation, IAU 2000B model.
**
**  Given:
**     date1,date2   double    TT as a 2-part Julian Date (Note 1)
**
**  Returned:
**     dpsi,deps     double    nutation, luni-solar + planetary (Note 2)
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
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
**  2) The nutation components in longitude and obliquity are in radians
**     and with respect to the equinox and ecliptic of date.  The
**     obliquity at J2000.0 is assumed to be the Lieske et al. (1977)
**     value of 84381.448 arcsec.  (The errors that result from using
**     this function with the IAU 2006 value of 84381.406 arcsec can be
**     neglected.)
**
**     The nutation model consists only of luni-solar terms, but
**     includes also a fixed offset which compensates for certain long-
**     period planetary terms (Note 7).
**
**  3) This function is an implementation of the IAU 2000B abridged
**     nutation model formally adopted by the IAU General Assembly in
**     2000.  The function computes the MHB_2000_SHORT luni-solar
**     nutation series (Luzum 2001), but without the associated
**     corrections for the precession rate adjustments and the offset
**     between the GCRS and J2000.0 mean poles.
**
**  4) The full IAU 2000A (MHB2000) nutation model contains nearly 1400
**     terms.  The IAU 2000B model (McCarthy & Luzum 2003) contains only
**     77 terms, plus additional simplifications, yet still delivers
**     results of 1 mas accuracy at present epochs.  This combination of
**     accuracy and size makes the IAU 2000B abridged nutation model
**     suitable for most practical applications.
**
**     The function delivers a pole accurate to 1 mas from 1900 to 2100
**     (usually better than 1 mas, very occasionally just outside
**     1 mas).  The full IAU 2000A model, which is implemented in the
**     function eraNut00a (q.v.), delivers considerably greater accuracy
**     at current dates;  however, to realize this improved accuracy,
**     corrections for the essentially unpredictable free-core-nutation
**     (FCN) must also be included.
**
**  5) The present function provides classical nutation.  The
**     MHB_2000_SHORT algorithm, from which it is adapted, deals also
**     with (i) the offsets between the GCRS and mean poles and (ii) the
**     adjustments in longitude and obliquity due to the changed
**     precession rates.  These additional functions, namely frame bias
**     and precession adjustments, are supported by the ERFA functions
**     eraBi00  and eraPr00.
**
**  6) The MHB_2000_SHORT algorithm also provides "total" nutations,
**     comprising the arithmetic sum of the frame bias, precession
**     adjustments, and nutation (luni-solar + planetary).  These total
**     nutations can be used in combination with an existing IAU 1976
**     precession implementation, such as eraPmat76,  to deliver GCRS-
**     to-true predictions of mas accuracy at current epochs.  However,
**     for symmetry with the eraNut00a  function (q.v. for the reasons),
**     the ERFA functions do not generate the "total nutations"
**     directly.  Should they be required, they could of course easily
**     be generated by calling eraBi00, eraPr00 and the present function
**     and adding the results.
**
**  7) The IAU 2000B model includes "planetary bias" terms that are
**     fixed in size but compensate for long-period nutations.  The
**     amplitudes quoted in McCarthy & Luzum (2003), namely
**     Dpsi = -1.5835 mas and Depsilon = +1.6339 mas, are optimized for
**     the "total nutations" method described in Note 6.  The Luzum
**     (2001) values used in this ERFA implementation, namely -0.135 mas
**     and +0.388 mas, are optimized for the "rigorous" method, where
**     frame bias, precession and nutation are applied separately and in
**     that order.  During the interval 1995-2050, the ERFA
**     implementation delivers a maximum error of 1.001 mas (not
**     including FCN).
**
**  References:
**
**     Lieske, J.H., Lederle, T., Fricke, W., Morando, B., "Expressions
**     for the precession quantities based upon the IAU /1976/ system of
**     astronomical constants", Astron.Astrophys. 58, 1-2, 1-16. (1977)
**
**     Luzum, B., private communication, 2001 (Fortran code
**     MHB_2000_SHORT)
**
**     McCarthy, D.D. & Luzum, B.J., "An abridged model of the
**     precession-nutation of the celestial pole", Cel.Mech.Dyn.Astron.
**     85, 37-49 (2003)
**
**     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
**     Francou, G., Laskar, J., Astron.Astrophys. 282, 663-683 (1994)
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double t, el, elp, f, d, om, arg, dp, de, sarg, carg,
          dpsils, depsls, dpsipl, depspl;
   int i;

/* Units of 0.1 microarcsecond to radians */
   static const double U2R = DAS2R / 1e7;

/* ---------------------------------------- */
/* Fixed offsets in lieu of planetary terms */
/* ---------------------------------------- */

   static const double DPPLAN = -0.135 * DMAS2R;
   static const double DEPLAN =  0.388 * DMAS2R;

/* --------------------------------------------------- */
/* Luni-solar nutation: argument and term coefficients */
/* --------------------------------------------------- */

/* The units for the sine and cosine coefficients are */
/* 0.1 microarcsec and the same per Julian century    */

   static const struct {
      int nl,nlp,nf,nd,nom; /* coefficients of l,l',F,D,Om */
      double ps,pst,pc;     /* longitude sin, t*sin, cos coefficients */
      double ec,ect,es;     /* obliquity cos, t*cos, sin coefficients */

   } x[] = {

   /* 1-10 */
      { 0, 0, 0, 0,1,
         -172064161.0, -174666.0, 33386.0, 92052331.0, 9086.0, 15377.0},
      { 0, 0, 2,-2,2,
           -13170906.0, -1675.0, -13696.0, 5730336.0, -3015.0, -4587.0},
      { 0, 0, 2, 0,2,-2276413.0,-234.0, 2796.0, 978459.0,-485.0,1374.0},
      { 0, 0, 0, 0,2,2074554.0,  207.0, -698.0,-897492.0, 470.0,-291.0},
      { 0, 1, 0, 0,0,1475877.0,-3633.0,11817.0, 73871.0,-184.0,-1924.0},
      { 0, 1, 2,-2,2,-516821.0, 1226.0, -524.0, 224386.0,-677.0,-174.0},
      { 1, 0, 0, 0,0, 711159.0,   73.0, -872.0,  -6750.0,   0.0, 358.0},
      { 0, 0, 2, 0,1,-387298.0, -367.0,  380.0, 200728.0,  18.0, 318.0},
      { 1, 0, 2, 0,2,-301461.0,  -36.0,  816.0, 129025.0, -63.0, 367.0},
      { 0,-1, 2,-2,2, 215829.0, -494.0,  111.0, -95929.0, 299.0, 132.0},

   /* 11-20 */
      { 0, 0, 2,-2,1, 128227.0,  137.0,  181.0, -68982.0,  -9.0,  39.0},
      {-1, 0, 2, 0,2, 123457.0,   11.0,   19.0, -53311.0,  32.0,  -4.0},
      {-1, 0, 0, 2,0, 156994.0,   10.0, -168.0,  -1235.0,   0.0,  82.0},
      { 1, 0, 0, 0,1,  63110.0,   63.0,   27.0, -33228.0,   0.0,  -9.0},
      {-1, 0, 0, 0,1, -57976.0,  -63.0, -189.0,  31429.0,   0.0, -75.0},
      {-1, 0, 2, 2,2, -59641.0,  -11.0,  149.0,  25543.0, -11.0,  66.0},
      { 1, 0, 2, 0,1, -51613.0,  -42.0,  129.0,  26366.0,   0.0,  78.0},
      {-2, 0, 2, 0,1,  45893.0,   50.0,   31.0, -24236.0, -10.0,  20.0},
      { 0, 0, 0, 2,0,  63384.0,   11.0, -150.0,  -1220.0,   0.0,  29.0},
      { 0, 0, 2, 2,2, -38571.0,   -1.0,  158.0,  16452.0, -11.0,  68.0},

   /* 21-30 */
      { 0,-2, 2,-2,2,  32481.0,    0.0,    0.0, -13870.0,   0.0,   0.0},
      {-2, 0, 0, 2,0, -47722.0,    0.0,  -18.0,    477.0,   0.0, -25.0},
      { 2, 0, 2, 0,2, -31046.0,   -1.0,  131.0,  13238.0, -11.0,  59.0},
      { 1, 0, 2,-2,2,  28593.0,    0.0,   -1.0, -12338.0,  10.0,  -3.0},
      {-1, 0, 2, 0,1,  20441.0,   21.0,   10.0, -10758.0,   0.0,  -3.0},
      { 2, 0, 0, 0,0,  29243.0,    0.0,  -74.0,   -609.0,   0.0,  13.0},
      { 0, 0, 2, 0,0,  25887.0,    0.0,  -66.0,   -550.0,   0.0,  11.0},
      { 0, 1, 0, 0,1, -14053.0,  -25.0,   79.0,   8551.0,  -2.0, -45.0},
      {-1, 0, 0, 2,1,  15164.0,   10.0,   11.0,  -8001.0,   0.0,  -1.0},
      { 0, 2, 2,-2,2, -15794.0,   72.0,  -16.0,   6850.0, -42.0,  -5.0},

   /* 31-40 */
      { 0, 0,-2, 2,0,  21783.0,    0.0,   13.0,   -167.0,   0.0,  13.0},
      { 1, 0, 0,-2,1, -12873.0,  -10.0,  -37.0,   6953.0,   0.0, -14.0},
      { 0,-1, 0, 0,1, -12654.0,   11.0,   63.0,   6415.0,   0.0,  26.0},
      {-1, 0, 2, 2,1, -10204.0,    0.0,   25.0,   5222.0,   0.0,  15.0},
      { 0, 2, 0, 0,0,  16707.0,  -85.0,  -10.0,    168.0,  -1.0,  10.0},
      { 1, 0, 2, 2,2,  -7691.0,    0.0,   44.0,   3268.0,   0.0,  19.0},
      {-2, 0, 2, 0,0, -11024.0,    0.0,  -14.0,    104.0,   0.0,   2.0},
      { 0, 1, 2, 0,2,   7566.0,  -21.0,  -11.0,  -3250.0,   0.0,  -5.0},
      { 0, 0, 2, 2,1,  -6637.0,  -11.0,   25.0,   3353.0,   0.0,  14.0},
      { 0,-1, 2, 0,2,  -7141.0,   21.0,    8.0,   3070.0,   0.0,   4.0},

   /* 41-50 */
      { 0, 0, 0, 2,1,  -6302.0,  -11.0,    2.0,   3272.0,   0.0,   4.0},
      { 1, 0, 2,-2,1,   5800.0,   10.0,    2.0,  -3045.0,   0.0,  -1.0},
      { 2, 0, 2,-2,2,   6443.0,    0.0,   -7.0,  -2768.0,   0.0,  -4.0},
      {-2, 0, 0, 2,1,  -5774.0,  -11.0,  -15.0,   3041.0,   0.0,  -5.0},
      { 2, 0, 2, 0,1,  -5350.0,    0.0,   21.0,   2695.0,   0.0,  12.0},
      { 0,-1, 2,-2,1,  -4752.0,  -11.0,   -3.0,   2719.0,   0.0,  -3.0},
      { 0, 0, 0,-2,1,  -4940.0,  -11.0,  -21.0,   2720.0,   0.0,  -9.0},
      {-1,-1, 0, 2,0,   7350.0,    0.0,   -8.0,    -51.0,   0.0,   4.0},
      { 2, 0, 0,-2,1,   4065.0,    0.0,    6.0,  -2206.0,   0.0,   1.0},
      { 1, 0, 0, 2,0,   6579.0,    0.0,  -24.0,   -199.0,   0.0,   2.0},

   /* 51-60 */
      { 0, 1, 2,-2,1,   3579.0,    0.0,    5.0,  -1900.0,   0.0,   1.0},
      { 1,-1, 0, 0,0,   4725.0,    0.0,   -6.0,    -41.0,   0.0,   3.0},
      {-2, 0, 2, 0,2,  -3075.0,    0.0,   -2.0,   1313.0,   0.0,  -1.0},
      { 3, 0, 2, 0,2,  -2904.0,    0.0,   15.0,   1233.0,   0.0,   7.0},
      { 0,-1, 0, 2,0,   4348.0,    0.0,  -10.0,    -81.0,   0.0,   2.0},
      { 1,-1, 2, 0,2,  -2878.0,    0.0,    8.0,   1232.0,   0.0,   4.0},
      { 0, 0, 0, 1,0,  -4230.0,    0.0,    5.0,    -20.0,   0.0,  -2.0},
      {-1,-1, 2, 2,2,  -2819.0,    0.0,    7.0,   1207.0,   0.0,   3.0},
      {-1, 0, 2, 0,0,  -4056.0,    0.0,    5.0,     40.0,   0.0,  -2.0},
      { 0,-1, 2, 2,2,  -2647.0,    0.0,   11.0,   1129.0,   0.0,   5.0},

   /* 61-70 */
      {-2, 0, 0, 0,1,  -2294.0,    0.0,  -10.0,   1266.0,   0.0,  -4.0},
      { 1, 1, 2, 0,2,   2481.0,    0.0,   -7.0,  -1062.0,   0.0,  -3.0},
      { 2, 0, 0, 0,1,   2179.0,    0.0,   -2.0,  -1129.0,   0.0,  -2.0},
      {-1, 1, 0, 1,0,   3276.0,    0.0,    1.0,     -9.0,   0.0,   0.0},
      { 1, 1, 0, 0,0,  -3389.0,    0.0,    5.0,     35.0,   0.0,  -2.0},
      { 1, 0, 2, 0,0,   3339.0,    0.0,  -13.0,   -107.0,   0.0,   1.0},
      {-1, 0, 2,-2,1,  -1987.0,    0.0,   -6.0,   1073.0,   0.0,  -2.0},
      { 1, 0, 0, 0,2,  -1981.0,    0.0,    0.0,    854.0,   0.0,   0.0},
      {-1, 0, 0, 1,0,   4026.0,    0.0, -353.0,   -553.0,   0.0,-139.0},
      { 0, 0, 2, 1,2,   1660.0,    0.0,   -5.0,   -710.0,   0.0,  -2.0},

   /* 71-77 */
      {-1, 0, 2, 4,2,  -1521.0,    0.0,    9.0,    647.0,   0.0,   4.0},
      {-1, 1, 0, 1,1,   1314.0,    0.0,    0.0,   -700.0,   0.0,   0.0},
      { 0,-2, 2,-2,1,  -1283.0,    0.0,    0.0,    672.0,   0.0,   0.0},
      { 1, 0, 2, 2,1,  -1331.0,    0.0,    8.0,    663.0,   0.0,   4.0},
      {-2, 0, 2, 2,2,   1383.0,    0.0,   -2.0,   -594.0,   0.0,  -2.0},
      {-1, 0, 0, 0,2,   1405.0,    0.0,    4.0,   -610.0,   0.0,   2.0},
      { 1, 1, 2,-2,2,   1290.0,    0.0,    0.0,   -556.0,   0.0,   0.0}
   };

/* Number of terms in the series */
   const int NLS = (int) (sizeof x / sizeof x[0]);

/*--------------------------------------------------------------------*/

/* Interval between fundamental epoch J2000.0 and given date (JC). */
   t = ((date1 - DJ00) + date2) / DJC;

/* --------------------*/
/* LUNI-SOLAR NUTATION */
/* --------------------*/

/* Fundamental (Delaunay) arguments from Simon et al. (1994) */

/* Mean anomaly of the Moon. */
   el = fmod(485868.249036 + (1717915923.2178) * t, TURNAS) * DAS2R;

/* Mean anomaly of the Sun. */
   elp = fmod(1287104.79305 + (129596581.0481) * t, TURNAS) * DAS2R;

/* Mean argument of the latitude of the Moon. */
   f = fmod(335779.526232 + (1739527262.8478) * t, TURNAS) * DAS2R;

/* Mean elongation of the Moon from the Sun. */
   d = fmod(1072260.70369 + (1602961601.2090) * t, TURNAS) * DAS2R;

/* Mean longitude of the ascending node of the Moon. */
   om = fmod(450160.398036 + (-6962890.5431) * t, TURNAS) * DAS2R;

/* Initialize the nutation values. */
   dp = 0.0;
   de = 0.0;

/* Summation of luni-solar nutation series (smallest terms first). */
   for (i = NLS-1; i >= 0; i--) {

   /* Argument and functions. */
      arg = fmod( (double)x[i].nl  * el  +
                  (double)x[i].nlp * elp +
                  (double)x[i].nf  * f   +
                  (double)x[i].nd  * d   +
                  (double)x[i].nom * om, D2PI  );
      sarg = sin(arg);
      carg = cos(arg);

   /* Term. */
      dp += (x[i].ps + x[i].pst * t) * sarg + x[i].pc * carg;
      de += (x[i].ec + x[i].ect * t) * carg + x[i].es * sarg;
   }

/* Convert from 0.1 microarcsec units to radians. */
   dpsils = dp * U2R;
   depsls = de * U2R;

/* ------------------------------*/
/* IN LIEU OF PLANETARY NUTATION */
/* ------------------------------*/

/* Fixed offset to correct for missing terms in truncated series. */
   dpsipl = DPPLAN;
   depspl = DEPLAN;

/* --------*/
/* RESULTS */
/* --------*/

/* Add luni-solar and planetary components. */
   *dpsi = dpsils + dpsipl;
   *deps = depsls + depspl;

   return;

}

void eraNut06a(double date1, double date2, double *dpsi, double *deps)
/*
**  - - - - - - - - - -
**   e r a N u t 0 6 a
**  - - - - - - - - - -
**
**  IAU 2000A nutation with adjustments to match the IAU 2006
**  precession.
**
**  Given:
**     date1,date2   double   TT as a 2-part Julian Date (Note 1)
**
**  Returned:
**     dpsi,deps     double   nutation, luni-solar + planetary (Note 2)
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
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
**  2) The nutation components in longitude and obliquity are in radians
**     and with respect to the mean equinox and ecliptic of date,
**     IAU 2006 precession model (Hilton et al. 2006, Capitaine et al.
**     2005).
**
**  3) The function first computes the IAU 2000A nutation, then applies
**     adjustments for (i) the consequences of the change in obliquity
**     from the IAU 1980 ecliptic to the IAU 2006 ecliptic and (ii) the
**     secular variation in the Earth's dynamical form factor J2.
**
**  4) The present function provides classical nutation, complementing
**     the IAU 2000 frame bias and IAU 2006 precession.  It delivers a
**     pole which is at current epochs accurate to a few tens of
**     microarcseconds, apart from the free core nutation.
**
**  Called:
**     eraNut00a    nutation, IAU 2000A
**
**  References:
**
**     Chapront, J., Chapront-Touze, M. & Francou, G. 2002,
**     Astron.Astrophys. 387, 700
**
**     Lieske, J.H., Lederle, T., Fricke, W. & Morando, B. 1977,
**     Astron.Astrophys. 58, 1-16
**
**     Mathews, P.M., Herring, T.A., Buffet, B.A. 2002, J.Geophys.Res.
**     107, B4.  The MHB_2000 code itself was obtained on 9th September
**     2002 from ftp//maia.usno.navy.mil/conv2000/chapter5/IAU2000A.
**
**     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
**     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
**
**     Souchay, J., Loysel, B., Kinoshita, H., Folgueira, M. 1999,
**     Astron.Astrophys.Supp.Ser. 135, 111
**
**     Wallace, P.T., "Software for Implementing the IAU 2000
**     Resolutions", in IERS Workshop 5.1 (2002)
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double t, fj2, dp, de;


/* Interval between fundamental date J2000.0 and given date (JC). */
   t = ((date1 - DJ00) + date2) / DJC;

/* Factor correcting for secular variation of J2. */
   fj2 = -2.7774e-6 * t;

/* Obtain IAU 2000A nutation. */
   eraNut00a(date1, date2, &dp, &de);

/* Apply P03 adjustments (Wallace & Capitaine, 2006, Eqs.5). */
   *dpsi = dp + dp * (0.4697e-6 + fj2);
   *deps = de + de * fj2;

   return;

}

void eraNut80(double date1, double date2, double *dpsi, double *deps)
/*
**  - - - - - - - - -
**   e r a N u t 8 0
**  - - - - - - - - -
**
**  Nutation, IAU 1980 model.
**
**  Given:
**     date1,date2   double    TT as a 2-part Julian Date (Note 1)
**
**  Returned:
**     dpsi          double    nutation in longitude (radians)
**     deps          double    nutation in obliquity (radians)
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
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
**  2) The nutation components are with respect to the ecliptic of
**     date.
**
**  Called:
**     eraAnpm      normalize angle into range +/- pi
**
**  Reference:
**
**     Explanatory Supplement to the Astronomical Almanac,
**     P. Kenneth Seidelmann (ed), University Science Books (1992),
**     Section 3.222 (p111).
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double t, el, elp, f, d, om, dp, de, arg, s, c;
   int j;

/* Units of 0.1 milliarcsecond to radians */
   const double U2R = DAS2R / 1e4;

/* ------------------------------------------------ */
/* Table of multiples of arguments and coefficients */
/* ------------------------------------------------ */

/* The units for the sine and cosine coefficients are 0.1 mas and */
/* the same per Julian century */

   static const struct {
      int nl,nlp,nf,nd,nom; /* coefficients of l,l',F,D,Om */
      double sp,spt;        /* longitude sine, 1 and t coefficients */
      double ce,cet;        /* obliquity cosine, 1 and t coefficients */
   } x[] = {

   /* 1-10 */
      {  0,  0,  0,  0,  1, -171996.0, -174.2,  92025.0,    8.9 },
      {  0,  0,  0,  0,  2,    2062.0,    0.2,   -895.0,    0.5 },
      { -2,  0,  2,  0,  1,      46.0,    0.0,    -24.0,    0.0 },
      {  2,  0, -2,  0,  0,      11.0,    0.0,      0.0,    0.0 },
      { -2,  0,  2,  0,  2,      -3.0,    0.0,      1.0,    0.0 },
      {  1, -1,  0, -1,  0,      -3.0,    0.0,      0.0,    0.0 },
      {  0, -2,  2, -2,  1,      -2.0,    0.0,      1.0,    0.0 },
      {  2,  0, -2,  0,  1,       1.0,    0.0,      0.0,    0.0 },
      {  0,  0,  2, -2,  2,  -13187.0,   -1.6,   5736.0,   -3.1 },
      {  0,  1,  0,  0,  0,    1426.0,   -3.4,     54.0,   -0.1 },

   /* 11-20 */
      {  0,  1,  2, -2,  2,    -517.0,    1.2,    224.0,   -0.6 },
      {  0, -1,  2, -2,  2,     217.0,   -0.5,    -95.0,    0.3 },
      {  0,  0,  2, -2,  1,     129.0,    0.1,    -70.0,    0.0 },
      {  2,  0,  0, -2,  0,      48.0,    0.0,      1.0,    0.0 },
      {  0,  0,  2, -2,  0,     -22.0,    0.0,      0.0,    0.0 },
      {  0,  2,  0,  0,  0,      17.0,   -0.1,      0.0,    0.0 },
      {  0,  1,  0,  0,  1,     -15.0,    0.0,      9.0,    0.0 },
      {  0,  2,  2, -2,  2,     -16.0,    0.1,      7.0,    0.0 },
      {  0, -1,  0,  0,  1,     -12.0,    0.0,      6.0,    0.0 },
      { -2,  0,  0,  2,  1,      -6.0,    0.0,      3.0,    0.0 },

   /* 21-30 */
      {  0, -1,  2, -2,  1,      -5.0,    0.0,      3.0,    0.0 },
      {  2,  0,  0, -2,  1,       4.0,    0.0,     -2.0,    0.0 },
      {  0,  1,  2, -2,  1,       4.0,    0.0,     -2.0,    0.0 },
      {  1,  0,  0, -1,  0,      -4.0,    0.0,      0.0,    0.0 },
      {  2,  1,  0, -2,  0,       1.0,    0.0,      0.0,    0.0 },
      {  0,  0, -2,  2,  1,       1.0,    0.0,      0.0,    0.0 },
      {  0,  1, -2,  2,  0,      -1.0,    0.0,      0.0,    0.0 },
      {  0,  1,  0,  0,  2,       1.0,    0.0,      0.0,    0.0 },
      { -1,  0,  0,  1,  1,       1.0,    0.0,      0.0,    0.0 },
      {  0,  1,  2, -2,  0,      -1.0,    0.0,      0.0,    0.0 },

   /* 31-40 */
      {  0,  0,  2,  0,  2,   -2274.0,   -0.2,    977.0,   -0.5 },
      {  1,  0,  0,  0,  0,     712.0,    0.1,     -7.0,    0.0 },
      {  0,  0,  2,  0,  1,    -386.0,   -0.4,    200.0,    0.0 },
      {  1,  0,  2,  0,  2,    -301.0,    0.0,    129.0,   -0.1 },
      {  1,  0,  0, -2,  0,    -158.0,    0.0,     -1.0,    0.0 },
      { -1,  0,  2,  0,  2,     123.0,    0.0,    -53.0,    0.0 },
      {  0,  0,  0,  2,  0,      63.0,    0.0,     -2.0,    0.0 },
      {  1,  0,  0,  0,  1,      63.0,    0.1,    -33.0,    0.0 },
      { -1,  0,  0,  0,  1,     -58.0,   -0.1,     32.0,    0.0 },
      { -1,  0,  2,  2,  2,     -59.0,    0.0,     26.0,    0.0 },

   /* 41-50 */
      {  1,  0,  2,  0,  1,     -51.0,    0.0,     27.0,    0.0 },
      {  0,  0,  2,  2,  2,     -38.0,    0.0,     16.0,    0.0 },
      {  2,  0,  0,  0,  0,      29.0,    0.0,     -1.0,    0.0 },
      {  1,  0,  2, -2,  2,      29.0,    0.0,    -12.0,    0.0 },
      {  2,  0,  2,  0,  2,     -31.0,    0.0,     13.0,    0.0 },
      {  0,  0,  2,  0,  0,      26.0,    0.0,     -1.0,    0.0 },
      { -1,  0,  2,  0,  1,      21.0,    0.0,    -10.0,    0.0 },
      { -1,  0,  0,  2,  1,      16.0,    0.0,     -8.0,    0.0 },
      {  1,  0,  0, -2,  1,     -13.0,    0.0,      7.0,    0.0 },
      { -1,  0,  2,  2,  1,     -10.0,    0.0,      5.0,    0.0 },

   /* 51-60 */
      {  1,  1,  0, -2,  0,      -7.0,    0.0,      0.0,    0.0 },
      {  0,  1,  2,  0,  2,       7.0,    0.0,     -3.0,    0.0 },
      {  0, -1,  2,  0,  2,      -7.0,    0.0,      3.0,    0.0 },
      {  1,  0,  2,  2,  2,      -8.0,    0.0,      3.0,    0.0 },
      {  1,  0,  0,  2,  0,       6.0,    0.0,      0.0,    0.0 },
      {  2,  0,  2, -2,  2,       6.0,    0.0,     -3.0,    0.0 },
      {  0,  0,  0,  2,  1,      -6.0,    0.0,      3.0,    0.0 },
      {  0,  0,  2,  2,  1,      -7.0,    0.0,      3.0,    0.0 },
      {  1,  0,  2, -2,  1,       6.0,    0.0,     -3.0,    0.0 },
      {  0,  0,  0, -2,  1,      -5.0,    0.0,      3.0,    0.0 },

   /* 61-70 */
      {  1, -1,  0,  0,  0,       5.0,    0.0,      0.0,    0.0 },
      {  2,  0,  2,  0,  1,      -5.0,    0.0,      3.0,    0.0 },
      {  0,  1,  0, -2,  0,      -4.0,    0.0,      0.0,    0.0 },
      {  1,  0, -2,  0,  0,       4.0,    0.0,      0.0,    0.0 },
      {  0,  0,  0,  1,  0,      -4.0,    0.0,      0.0,    0.0 },
      {  1,  1,  0,  0,  0,      -3.0,    0.0,      0.0,    0.0 },
      {  1,  0,  2,  0,  0,       3.0,    0.0,      0.0,    0.0 },
      {  1, -1,  2,  0,  2,      -3.0,    0.0,      1.0,    0.0 },
      { -1, -1,  2,  2,  2,      -3.0,    0.0,      1.0,    0.0 },
      { -2,  0,  0,  0,  1,      -2.0,    0.0,      1.0,    0.0 },

   /* 71-80 */
      {  3,  0,  2,  0,  2,      -3.0,    0.0,      1.0,    0.0 },
      {  0, -1,  2,  2,  2,      -3.0,    0.0,      1.0,    0.0 },
      {  1,  1,  2,  0,  2,       2.0,    0.0,     -1.0,    0.0 },
      { -1,  0,  2, -2,  1,      -2.0,    0.0,      1.0,    0.0 },
      {  2,  0,  0,  0,  1,       2.0,    0.0,     -1.0,    0.0 },
      {  1,  0,  0,  0,  2,      -2.0,    0.0,      1.0,    0.0 },
      {  3,  0,  0,  0,  0,       2.0,    0.0,      0.0,    0.0 },
      {  0,  0,  2,  1,  2,       2.0,    0.0,     -1.0,    0.0 },
      { -1,  0,  0,  0,  2,       1.0,    0.0,     -1.0,    0.0 },
      {  1,  0,  0, -4,  0,      -1.0,    0.0,      0.0,    0.0 },

   /* 81-90 */
      { -2,  0,  2,  2,  2,       1.0,    0.0,     -1.0,    0.0 },
      { -1,  0,  2,  4,  2,      -2.0,    0.0,      1.0,    0.0 },
      {  2,  0,  0, -4,  0,      -1.0,    0.0,      0.0,    0.0 },
      {  1,  1,  2, -2,  2,       1.0,    0.0,     -1.0,    0.0 },
      {  1,  0,  2,  2,  1,      -1.0,    0.0,      1.0,    0.0 },
      { -2,  0,  2,  4,  2,      -1.0,    0.0,      1.0,    0.0 },
      { -1,  0,  4,  0,  2,       1.0,    0.0,      0.0,    0.0 },
      {  1, -1,  0, -2,  0,       1.0,    0.0,      0.0,    0.0 },
      {  2,  0,  2, -2,  1,       1.0,    0.0,     -1.0,    0.0 },
      {  2,  0,  2,  2,  2,      -1.0,    0.0,      0.0,    0.0 },

   /* 91-100 */
      {  1,  0,  0,  2,  1,      -1.0,    0.0,      0.0,    0.0 },
      {  0,  0,  4, -2,  2,       1.0,    0.0,      0.0,    0.0 },
      {  3,  0,  2, -2,  2,       1.0,    0.0,      0.0,    0.0 },
      {  1,  0,  2, -2,  0,      -1.0,    0.0,      0.0,    0.0 },
      {  0,  1,  2,  0,  1,       1.0,    0.0,      0.0,    0.0 },
      { -1, -1,  0,  2,  1,       1.0,    0.0,      0.0,    0.0 },
      {  0,  0, -2,  0,  1,      -1.0,    0.0,      0.0,    0.0 },
      {  0,  0,  2, -1,  2,      -1.0,    0.0,      0.0,    0.0 },
      {  0,  1,  0,  2,  0,      -1.0,    0.0,      0.0,    0.0 },
      {  1,  0, -2, -2,  0,      -1.0,    0.0,      0.0,    0.0 },

   /* 101-106 */
      {  0, -1,  2,  0,  1,      -1.0,    0.0,      0.0,    0.0 },
      {  1,  1,  0, -2,  1,      -1.0,    0.0,      0.0,    0.0 },
      {  1,  0, -2,  2,  0,      -1.0,    0.0,      0.0,    0.0 },
      {  2,  0,  0,  2,  0,       1.0,    0.0,      0.0,    0.0 },
      {  0,  0,  2,  4,  2,      -1.0,    0.0,      0.0,    0.0 },
      {  0,  1,  0,  1,  0,       1.0,    0.0,      0.0,    0.0 }
   };

/* Number of terms in the series */
   const int NT = (int) (sizeof x / sizeof x[0]);

/*--------------------------------------------------------------------*/

/* Interval between fundamental epoch J2000.0 and given date (JC). */
   t = ((date1 - DJ00) + date2) / DJC;

/* --------------------- */
/* Fundamental arguments */
/* --------------------- */

/* Mean longitude of Moon minus mean longitude of Moon's perigee. */
   el = eraAnpm(
        (485866.733 + (715922.633 + (31.310 + 0.064 * t) * t) * t)
        * DAS2R + fmod(1325.0 * t, 1.0) * D2PI);

/* Mean longitude of Sun minus mean longitude of Sun's perigee. */
   elp = eraAnpm(
         (1287099.804 + (1292581.224 + (-0.577 - 0.012 * t) * t) * t)
         * DAS2R + fmod(99.0 * t, 1.0) * D2PI);

/* Mean longitude of Moon minus mean longitude of Moon's node. */
   f = eraAnpm(
       (335778.877 + (295263.137 + (-13.257 + 0.011 * t) * t) * t)
       * DAS2R + fmod(1342.0 * t, 1.0) * D2PI);

/* Mean elongation of Moon from Sun. */
   d = eraAnpm(
       (1072261.307 + (1105601.328 + (-6.891 + 0.019 * t) * t) * t)
       * DAS2R + fmod(1236.0 * t, 1.0) * D2PI);

/* Longitude of the mean ascending node of the lunar orbit on the */
/* ecliptic, measured from the mean equinox of date. */
   om = eraAnpm(
        (450160.280 + (-482890.539 + (7.455 + 0.008 * t) * t) * t)
        * DAS2R + fmod(-5.0 * t, 1.0) * D2PI);

/* --------------- */
/* Nutation series */
/* --------------- */

/* Initialize nutation components. */
   dp = 0.0;
   de = 0.0;

/* Sum the nutation terms, ending with the biggest. */
   for (j = NT-1; j >= 0; j--) {

   /* Form argument for current term. */
      arg = (double)x[j].nl  * el
          + (double)x[j].nlp * elp
          + (double)x[j].nf  * f
          + (double)x[j].nd  * d
          + (double)x[j].nom * om;

   /* Accumulate current nutation term. */
      s = x[j].sp + x[j].spt * t;
      c = x[j].ce + x[j].cet * t;
      if (s != 0.0) dp += s * sin(arg);
      if (c != 0.0) de += c * cos(arg);
   }

/* Convert results from 0.1 mas units to radians. */
   *dpsi = dp * U2R;
   *deps = de * U2R;

   return;

}

void eraNutm80(double date1, double date2, double rmatn[3][3])
/*
**  - - - - - - - - - -
**   e r a N u t m 8 0
**  - - - - - - - - - -
**
**  Form the matrix of nutation for a given date, IAU 1980 model.
**
**  Given:
**     date1,date2    double          TDB date (Note 1)
**
**  Returned:
**     rmatn          double[3][3]    nutation matrix
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
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
**  2) The matrix operates in the sense V(true) = rmatn * V(mean),
**     where the p-vector V(true) is with respect to the true
**     equatorial triad of date and the p-vector V(mean) is with
**     respect to the mean equatorial triad of date.
**
**  Called:
**     eraNut80     nutation, IAU 1980
**     eraObl80     mean obliquity, IAU 1980
**     eraNumat     form nutation matrix
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double dpsi, deps, epsa;


/* Nutation components and mean obliquity. */
   eraNut80(date1, date2, &dpsi, &deps);
   epsa = eraObl80(date1, date2);

/* Build the rotation matrix. */
   eraNumat(epsa, dpsi, deps, rmatn);

   return;

}

double eraObl06(double date1, double date2)
/*
**  - - - - - - - - -
**   e r a O b l 0 6
**  - - - - - - - - -
**
**  Mean obliquity of the ecliptic, IAU 2006 precession model.
**
**  Given:
**     date1,date2  double   TT as a 2-part Julian Date (Note 1)
**
**  Returned (function value):
**                  double   obliquity of the ecliptic (radians, Note 2)
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
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
**  2) The result is the angle between the ecliptic and mean equator of
**     date date1+date2.
**
**  Reference:
**
**     Hilton, J. et al., 2006, Celest.Mech.Dyn.Astron. 94, 351
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double t, eps0;


/* Interval between fundamental date J2000.0 and given date (JC). */
   t = ((date1 - DJ00) + date2) / DJC;

/* Mean obliquity. */
   eps0 = (84381.406     +
          (-46.836769    +
          ( -0.0001831   +
          (  0.00200340  +
          ( -0.000000576 +
          ( -0.0000000434) * t) * t) * t) * t) * t) * DAS2R;

   return eps0;

}

double eraObl80(double date1, double date2)
/*
**  - - - - - - - - -
**   e r a O b l 8 0
**  - - - - - - - - -
**
**  Mean obliquity of the ecliptic, IAU 1980 model.
**
**  Given:
**     date1,date2   double    TT as a 2-part Julian Date (Note 1)
**
**  Returned (function value):
**                   double    obliquity of the ecliptic (radians, Note 2)
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
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
**  2) The result is the angle between the ecliptic and mean equator of
**     date date1+date2.
**
**  Reference:
**
**     Explanatory Supplement to the Astronomical Almanac,
**     P. Kenneth Seidelmann (ed), University Science Books (1992),
**     Expression 3.222-1 (p114).
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double t, eps0;


/* Interval between fundamental epoch J2000.0 and given date (JC). */
   t = ((date1 - DJ00) + date2) / DJC;

/* Mean obliquity of date. */
   eps0 = DAS2R * (84381.448  +
                  (-46.8150   +
                  (-0.00059   +
                  ( 0.001813) * t) * t) * t);

   return eps0;

}

void eraP06e(double date1, double date2,
             double *eps0, double *psia, double *oma, double *bpa,
             double *bqa, double *pia, double *bpia,
             double *epsa, double *chia, double *za, double *zetaa,
             double *thetaa, double *pa,
             double *gam, double *phi, double *psi)
/*
**  - - - - - - - -
**   e r a P 0 6 e
**  - - - - - - - -
**
**  Precession angles, IAU 2006, equinox based.
**
**  Given:
**     date1,date2   double   TT as a 2-part Julian Date (Note 1)
**
**  Returned (see Note 2):
**     eps0          double   epsilon_0
**     psia          double   psi_A
**     oma           double   omega_A
**     bpa           double   P_A
**     bqa           double   Q_A
**     pia           double   pi_A
**     bpia          double   Pi_A
**     epsa          double   obliquity epsilon_A
**     chia          double   chi_A
**     za            double   z_A
**     zetaa         double   zeta_A
**     thetaa        double   theta_A
**     pa            double   p_A
**     gam           double   F-W angle gamma_J2000
**     phi           double   F-W angle phi_J2000
**     psi           double   F-W angle psi_J2000
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
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
**  2) This function returns the set of equinox based angles for the
**     Capitaine et al. "P03" precession theory, adopted by the IAU in
**     2006.  The angles are set out in Table 1 of Hilton et al. (2006):
**
**     eps0   epsilon_0   obliquity at J2000.0
**     psia   psi_A       luni-solar precession
**     oma    omega_A     inclination of equator wrt J2000.0 ecliptic
**     bpa    P_A         ecliptic pole x, J2000.0 ecliptic triad
**     bqa    Q_A         ecliptic pole -y, J2000.0 ecliptic triad
**     pia    pi_A        angle between moving and J2000.0 ecliptics
**     bpia   Pi_A        longitude of ascending node of the ecliptic
**     epsa   epsilon_A   obliquity of the ecliptic
**     chia   chi_A       planetary precession
**     za     z_A         equatorial precession: -3rd 323 Euler angle
**     zetaa  zeta_A      equatorial precession: -1st 323 Euler angle
**     thetaa theta_A     equatorial precession: 2nd 323 Euler angle
**     pa     p_A         general precession
**     gam    gamma_J2000 J2000.0 RA difference of ecliptic poles
**     phi    phi_J2000   J2000.0 codeclination of ecliptic pole
**     psi    psi_J2000   longitude difference of equator poles, J2000.0
**
**     The returned values are all radians.
**
**  3) Hilton et al. (2006) Table 1 also contains angles that depend on
**     models distinct from the P03 precession theory itself, namely the
**     IAU 2000A frame bias and nutation.  The quoted polynomials are
**     used in other ERFA functions:
**
**     . eraXy06  contains the polynomial parts of the X and Y series.
**
**     . eraS06  contains the polynomial part of the s+XY/2 series.
**
**     . eraPfw06  implements the series for the Fukushima-Williams
**       angles that are with respect to the GCRS pole (i.e. the variants
**       that include frame bias).
**
**  4) The IAU resolution stipulated that the choice of parameterization
**     was left to the user, and so an IAU compliant precession
**     implementation can be constructed using various combinations of
**     the angles returned by the present function.
**
**  5) The parameterization used by ERFA is the version of the Fukushima-
**     Williams angles that refers directly to the GCRS pole.  These
**     angles may be calculated by calling the function eraPfw06.  ERFA
**     also supports the direct computation of the CIP GCRS X,Y by
**     series, available by calling eraXy06.
**
**  6) The agreement between the different parameterizations is at the
**     1 microarcsecond level in the present era.
**
**  7) When constructing a precession formulation that refers to the GCRS
**     pole rather than the dynamical pole, it may (depending on the
**     choice of angles) be necessary to introduce the frame bias
**     explicitly.
**
**  8) It is permissible to re-use the same variable in the returned
**     arguments.  The quantities are stored in the stated order.
**
**  Reference:
**
**     Hilton, J. et al., 2006, Celest.Mech.Dyn.Astron. 94, 351
**
**  Called:
**     eraObl06     mean obliquity, IAU 2006
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double t;


/* Interval between fundamental date J2000.0 and given date (JC). */
   t = ((date1 - DJ00) + date2) / DJC;

/* Obliquity at J2000.0. */

   *eps0 = 84381.406 * DAS2R;

/* Luni-solar precession. */

   *psia = ( 5038.481507     +
           (   -1.0790069    +
           (   -0.00114045   +
           (    0.000132851  +
           (   -0.0000000951 )
           * t) * t) * t) * t) * t * DAS2R;

/* Inclination of mean equator with respect to the J2000.0 ecliptic. */

   *oma = *eps0 + ( -0.025754     +
                  (  0.0512623    +
                  ( -0.00772503   +
                  ( -0.000000467  +
                  (  0.0000003337 )
                  * t) * t) * t) * t) * t * DAS2R;

/* Ecliptic pole x, J2000.0 ecliptic triad. */

   *bpa = (  4.199094     +
          (  0.1939873    +
          ( -0.00022466   +
          ( -0.000000912  +
          (  0.0000000120 )
          * t) * t) * t) * t) * t * DAS2R;

/* Ecliptic pole -y, J2000.0 ecliptic triad. */

   *bqa = ( -46.811015     +
          (   0.0510283    +
          (   0.00052413   +
          (  -0.000000646  +
          (  -0.0000000172 )
          * t) * t) * t) * t) * t * DAS2R;

/* Angle between moving and J2000.0 ecliptics. */

   *pia = ( 46.998973     +
          ( -0.0334926    +
          ( -0.00012559   +
          (  0.000000113  +
          ( -0.0000000022 )
          * t) * t) * t) * t) * t * DAS2R;

/* Longitude of ascending node of the moving ecliptic. */

   *bpia = ( 629546.7936      +
           (   -867.95758     +
           (      0.157992    +
           (     -0.0005371   +
           (     -0.00004797  +
           (      0.000000072 )
           * t) * t) * t) * t) * t) * DAS2R;

/* Mean obliquity of the ecliptic. */

   *epsa = eraObl06(date1, date2);

/* Planetary precession. */

   *chia = ( 10.556403     +
           ( -2.3814292    +
           ( -0.00121197   +
           (  0.000170663  +
           ( -0.0000000560 )
           * t) * t) * t) * t) * t * DAS2R;

/* Equatorial precession: minus the third of the 323 Euler angles. */

   *za = (   -2.650545     +
         ( 2306.077181     +
         (    1.0927348    +
         (    0.01826837   +
         (   -0.000028596  +
         (   -0.0000002904 )
         * t) * t) * t) * t) * t) * DAS2R;

/* Equatorial precession: minus the first of the 323 Euler angles. */

   *zetaa = (    2.650545     +
            ( 2306.083227     +
            (    0.2988499    +
            (    0.01801828   +
            (   -0.000005971  +
            (   -0.0000003173 )
            * t) * t) * t) * t) * t) * DAS2R;

/* Equatorial precession: second of the 323 Euler angles. */

   *thetaa = ( 2004.191903     +
             (   -0.4294934    +
             (   -0.04182264   +
             (   -0.000007089  +
             (   -0.0000001274 )
             * t) * t) * t) * t) * t * DAS2R;

/* General precession. */

   *pa = ( 5028.796195     +
         (    1.1054348    +
         (    0.00007964   +
         (   -0.000023857  +
         (    0.0000000383 )
         * t) * t) * t) * t) * t * DAS2R;

/* Fukushima-Williams angles for precession. */

   *gam = ( 10.556403     +
          (  0.4932044    +
          ( -0.00031238   +
          ( -0.000002788  +
          (  0.0000000260 )
          * t) * t) * t) * t) * t * DAS2R;

   *phi = *eps0 + ( -46.811015     +
                  (   0.0511269    +
                  (   0.00053289   +
                  (  -0.000000440  +
                  (  -0.0000000176 )
                  * t) * t) * t) * t) * t * DAS2R;

   *psi = ( 5038.481507     +
          (    1.5584176    +
          (   -0.00018522   +
          (   -0.000026452  +
          (   -0.0000000148 )
          * t) * t) * t) * t) * t * DAS2R;

   return;

}

void eraP2pv(double p[3], double pv[2][3])
/*
**  - - - - - - - -
**   e r a P 2 p v
**  - - - - - - - -
**
**  Extend a p-vector to a pv-vector by appending a zero velocity.
**
**  Given:
**     p        double[3]       p-vector
**
**  Returned:
**     pv       double[2][3]    pv-vector
**
**  Called:
**     eraCp        copy p-vector
**     eraZp        zero p-vector
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   eraCp(p, pv[0]);
   eraZp(pv[1]);

   return;

}

void eraP2s(double p[3], double *theta, double *phi, double *r)
/*
**  - - - - - - -
**   e r a P 2 s
**  - - - - - - -
**
**  P-vector to spherical polar coordinates.
**
**  Given:
**     p        double[3]    p-vector
**
**  Returned:
**     theta    double       longitude angle (radians)
**     phi      double       latitude angle (radians)
**     r        double       radial distance
**
**  Notes:
**
**  1) If P is null, zero theta, phi and r are returned.
**
**  2) At either pole, zero theta is returned.
**
**  Called:
**     eraC2s       p-vector to spherical
**     eraPm        modulus of p-vector
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   eraC2s(p, theta, phi);
   *r = eraPm(p);

   return;

}

double eraPap(double a[3], double b[3])
/*
**  - - - - - - -
**   e r a P a p
**  - - - - - - -
**
**  Position-angle from two p-vectors.
**
**  Given:
**     a      double[3]  direction of reference point
**     b      double[3]  direction of point whose PA is required
**
**  Returned (function value):
**            double     position angle of b with respect to a (radians)
**
**  Notes:
**
**  1) The result is the position angle, in radians, of direction b with
**     respect to direction a.  It is in the range -pi to +pi.  The
**     sense is such that if b is a small distance "north" of a the
**     position angle is approximately zero, and if b is a small
**     distance "east" of a the position angle is approximately +pi/2.
**
**  2) The vectors a and b need not be of unit length.
**
**  3) Zero is returned if the two directions are the same or if either
**     vector is null.
**
**  4) If vector a is at a pole, the result is ill-defined.
**
**  Called:
**     eraPn        decompose p-vector into modulus and direction
**     eraPm        modulus of p-vector
**     eraPxp       vector product of two p-vectors
**     eraPmp       p-vector minus p-vector
**     eraPdp       scalar product of two p-vectors
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double am, au[3], bm, st, ct, xa, ya, za, eta[3], xi[3], a2b[3], pa;


/* Modulus and direction of the a vector. */
   eraPn(a, &am, au);

/* Modulus of the b vector. */
   bm = eraPm(b);

/* Deal with the case of a null vector. */
   if ((am == 0.0) || (bm == 0.0)) {
      st = 0.0;
      ct = 1.0;
   } else {

   /* The "north" axis tangential from a (arbitrary length). */
      xa = a[0];
      ya = a[1];
      za = a[2];
      eta[0] = -xa * za;
      eta[1] = -ya * za;
      eta[2] =  xa*xa + ya*ya;

   /* The "east" axis tangential from a (same length). */
      eraPxp(eta, au, xi);

   /* The vector from a to b. */
      eraPmp(b, a, a2b);

   /* Resolve into components along the north and east axes. */
      st = eraPdp(a2b, xi);
      ct = eraPdp(a2b, eta);

   /* Deal with degenerate cases. */
      if ((st == 0.0) && (ct == 0.0)) ct = 1.0;
   }

/* Position angle. */
   pa = atan2(st, ct);

   return pa;

}

double eraPas(double al, double ap, double bl, double bp)
/*
**  - - - - - - -
**   e r a P a s
**  - - - - - - -
**
**  Position-angle from spherical coordinates.
**
**  Given:
**     al     double     longitude of point A (e.g. RA) in radians
**     ap     double     latitude of point A (e.g. Dec) in radians
**     bl     double     longitude of point B
**     bp     double     latitude of point B
**
**  Returned (function value):
**            double     position angle of B with respect to A
**
**  Notes:
**
**  1) The result is the bearing (position angle), in radians, of point
**     B with respect to point A.  It is in the range -pi to +pi.  The
**     sense is such that if B is a small distance "east" of point A,
**     the bearing is approximately +pi/2.
**
**  2) Zero is returned if the two points are coincident.
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double dl, x, y, pa;


   dl = bl - al;
   y = sin(dl) * cos(bp);
   x = sin(bp) * cos(ap) - cos(bp) * sin(ap) * cos(dl);
   pa = ((x != 0.0) || (y != 0.0)) ? atan2(y, x) : 0.0;

   return pa;

}

void eraPb06(double date1, double date2,
             double *bzeta, double *bz, double *btheta)
/*
**  - - - - - - - -
**   e r a P b 0 6
**  - - - - - - - -
**
**  This function forms three Euler angles which implement general
**  precession from epoch J2000.0, using the IAU 2006 model.  Frame
**  bias (the offset between ICRS and mean J2000.0) is included.
**
**  Given:
**     date1,date2  double   TT as a 2-part Julian Date (Note 1)
**
**  Returned:
**     bzeta        double   1st rotation: radians cw around z
**     bz           double   3rd rotation: radians cw around z
**     btheta       double   2nd rotation: radians ccw around y
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
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
**  2) The traditional accumulated precession angles zeta_A, z_A,
**     theta_A cannot be obtained in the usual way, namely through
**     polynomial expressions, because of the frame bias.  The latter
**     means that two of the angles undergo rapid changes near this
**     date.  They are instead the results of decomposing the
**     precession-bias matrix obtained by using the Fukushima-Williams
**     method, which does not suffer from the problem.  The
**     decomposition returns values which can be used in the
**     conventional formulation and which include frame bias.
**
**  3) The three angles are returned in the conventional order, which
**     is not the same as the order of the corresponding Euler
**     rotations.  The precession-bias matrix is
**     R_3(-z) x R_2(+theta) x R_3(-zeta).
**
**  4) Should zeta_A, z_A, theta_A angles be required that do not
**     contain frame bias, they are available by calling the ERFA
**     function eraP06e.
**
**  Called:
**     eraPmat06    PB matrix, IAU 2006
**     eraRz        rotate around Z-axis
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double r[3][3], r31, r32;


/* Precession matrix via Fukushima-Williams angles. */
   eraPmat06(date1, date2, r);

/* Solve for z. */
   *bz = atan2(r[1][2], r[0][2]);

/* Remove it from the matrix. */
   eraRz(*bz, r);

/* Solve for the remaining two angles. */
   *bzeta = atan2 (r[1][0], r[1][1]);
   r31 = r[2][0];
   r32 = r[2][1];
   *btheta = atan2(-dsign(sqrt(r31 * r31 + r32 * r32), r[0][2]),
                   r[2][2]);

   return;

}

double eraPdp(double a[3], double b[3])
/*
**  - - - - - - -
**   e r a P d p
**  - - - - - - -
**
**  p-vector inner (=scalar=dot) product.
**
**  Given:
**     a      double[3]     first p-vector
**     b      double[3]     second p-vector
**
**  Returned (function value):
**            double        a . b
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double w;


   w  = a[0] * b[0]
      + a[1] * b[1]
      + a[2] * b[2];

   return w;

}

void eraPfw06(double date1, double date2,
              double *gamb, double *phib, double *psib, double *epsa)
/*
**  - - - - - - - - -
**   e r a P f w 0 6
**  - - - - - - - - -
**
**  Precession angles, IAU 2006 (Fukushima-Williams 4-angle formulation).
**
**  Given:
**     date1,date2  double   TT as a 2-part Julian Date (Note 1)
**
**  Returned:
**     gamb         double   F-W angle gamma_bar (radians)
**     phib         double   F-W angle phi_bar (radians)
**     psib         double   F-W angle psi_bar (radians)
**     epsa         double   F-W angle epsilon_A (radians)
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
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
**  2) Naming the following points:
**
**           e = J2000.0 ecliptic pole,
**           p = GCRS pole,
**           E = mean ecliptic pole of date,
**     and   P = mean pole of date,
**
**     the four Fukushima-Williams angles are as follows:
**
**        gamb = gamma_bar = epE
**        phib = phi_bar = pE
**        psib = psi_bar = pEP
**        epsa = epsilon_A = EP
**
**  3) The matrix representing the combined effects of frame bias and
**     precession is:
**
**        PxB = R_1(-epsa).R_3(-psib).R_1(phib).R_3(gamb)
**
**  4) The matrix representing the combined effects of frame bias,
**     precession and nutation is simply:
**
**        NxPxB = R_1(-epsa-dE).R_3(-psib-dP).R_1(phib).R_3(gamb)
**
**     where dP and dE are the nutation components with respect to the
**     ecliptic of date.
**
**  Reference:
**
**     Hilton, J. et al., 2006, Celest.Mech.Dyn.Astron. 94, 351
**
**  Called:
**     eraObl06     mean obliquity, IAU 2006
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double t;


/* Interval between fundamental date J2000.0 and given date (JC). */
   t = ((date1 - DJ00) + date2) / DJC;

/* P03 bias+precession angles. */
   *gamb = (    -0.052928     +
           (    10.556378     +
           (     0.4932044    +
           (    -0.00031238   +
           (    -0.000002788  +
           (     0.0000000260 )
           * t) * t) * t) * t) * t) * DAS2R;
   *phib = ( 84381.412819     +
           (   -46.811016     +
           (     0.0511268    +
           (     0.00053289   +
           (    -0.000000440  +
           (    -0.0000000176 )
           * t) * t) * t) * t) * t) * DAS2R;
   *psib = (    -0.041775     +
           (  5038.481484     +
           (     1.5584175    +
           (    -0.00018522   +
           (    -0.000026452  +
           (    -0.0000000148 )
           * t) * t) * t) * t) * t) * DAS2R;
   *epsa =  eraObl06(date1, date2);

   return;

}

int eraPlan94(double date1, double date2, int np, double pv[2][3])
/*
**  - - - - - - - - - -
**   e r a P l a n 9 4
**  - - - - - - - - - -
**
**  Approximate heliocentric position and velocity of a nominated major
**  planet:  Mercury, Venus, EMB, Mars, Jupiter, Saturn, Uranus or
**  Neptune (but not the Earth itself).
**
**  Given:
**     date1  double       TDB date part A (Note 1)
**     date2  double       TDB date part B (Note 1)
**     np     int          planet (1=Mercury, 2=Venus, 3=EMB, 4=Mars,
**                             5=Jupiter, 6=Saturn, 7=Uranus, 8=Neptune)
**
**  Returned (argument):
**     pv     double[2][3] planet p,v (heliocentric, J2000.0, AU,AU/d)
**
**  Returned (function value):
**            int          status: -1 = illegal NP (outside 1-8)
**                                  0 = OK
**                                 +1 = warning: year outside 1000-3000
**                                 +2 = warning: failed to converge
**
**  Notes:
**
**  1) The date date1+date2 is in the TDB time scale (in practice TT can
**     be used) and is a Julian Date, apportioned in any convenient way
**     between the two arguments.  For example, JD(TDB)=2450123.7 could
**     be expressed in any of these ways, among others:
**
**            date1          date2
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in cases
**     where the loss of several decimal digits of resolution is
**     acceptable.  The J2000 method is best matched to the way the
**     argument is handled internally and will deliver the optimum
**     resolution.  The MJD method and the date & time methods are both
**     good compromises between resolution and convenience.  The limited
**     accuracy of the present algorithm is such that any of the methods
**     is satisfactory.
**
**  2) If an np value outside the range 1-8 is supplied, an error status
**     (function value -1) is returned and the pv vector set to zeroes.
**
**  3) For np=3 the result is for the Earth-Moon Barycenter.  To obtain
**     the heliocentric position and velocity of the Earth, use instead
**     the ERFA function eraEpv00.
**
**  4) On successful return, the array pv contains the following:
**
**        pv[0][0]   x      }
**        pv[0][1]   y      } heliocentric position, AU
**        pv[0][2]   z      }
**
**        pv[1][0]   xdot   }
**        pv[1][1]   ydot   } heliocentric velocity, AU/d
**        pv[1][2]   zdot   }
**
**     The reference frame is equatorial and is with respect to the
**     mean equator and equinox of epoch J2000.0.
**
**  5) The algorithm is due to J.L. Simon, P. Bretagnon, J. Chapront,
**     M. Chapront-Touze, G. Francou and J. Laskar (Bureau des
**     Longitudes, Paris, France).  From comparisons with JPL
**     ephemeris DE102, they quote the following maximum errors
**     over the interval 1800-2050:
**
**                     L (arcsec)    B (arcsec)      R (km)
**
**        Mercury          4             1             300
**        Venus            5             1             800
**        EMB              6             1            1000
**        Mars            17             1            7700
**        Jupiter         71             5           76000
**        Saturn          81            13          267000
**        Uranus          86             7          712000
**        Neptune         11             1          253000
**
**     Over the interval 1000-3000, they report that the accuracy is no
**     worse than 1.5 times that over 1800-2050.  Outside 1000-3000 the
**     accuracy declines.
**
**     Comparisons of the present function with the JPL DE200 ephemeris
**     give the following RMS errors over the interval 1960-2025:
**
**                      position (km)     velocity (m/s)
**
**        Mercury            334               0.437
**        Venus             1060               0.855
**        EMB               2010               0.815
**        Mars              7690               1.98
**        Jupiter          71700               7.70
**        Saturn          199000              19.4
**        Uranus          564000              16.4
**        Neptune         158000              14.4
**
**     Comparisons against DE200 over the interval 1800-2100 gave the
**     following maximum absolute differences.  (The results using
**     DE406 were essentially the same.)
**
**                   L (arcsec)   B (arcsec)     R (km)   Rdot (m/s)
**
**        Mercury        7            1            500       0.7
**        Venus          7            1           1100       0.9
**        EMB            9            1           1300       1.0
**        Mars          26            1           9000       2.5
**        Jupiter       78            6          82000       8.2
**        Saturn        87           14         263000      24.6
**        Uranus        86            7         661000      27.4
**        Neptune       11            2         248000      21.4
**
**  6) The present ERFA re-implementation of the original Simon et al.
**     Fortran code differs from the original in the following respects:
**
**       *  C instead of Fortran.
**
**       *  The date is supplied in two parts.
**
**       *  The result is returned only in equatorial Cartesian form;
**          the ecliptic longitude, latitude and radius vector are not
**          returned.
**
**       *  The result is in the J2000.0 equatorial frame, not ecliptic.
**
**       *  More is done in-line: there are fewer calls to subroutines.
**
**       *  Different error/warning status values are used.
**
**       *  A different Kepler's-equation-solver is used (avoiding
**          use of double precision complex).
**
**       *  Polynomials in t are nested to minimize rounding errors.
**
**       *  Explicit double constants are used to avoid mixed-mode
**          expressions.
**
**     None of the above changes affects the result significantly.
**
**  7) The returned status indicates the most serious condition
**     encountered during execution of the function.  Illegal np is
**     considered the most serious, overriding failure to converge,
**     which in turn takes precedence over the remote date warning.
**
**  Called:
**     eraAnp       normalize angle into range 0 to 2pi
**
**  Reference:  Simon, J.L, Bretagnon, P., Chapront, J.,
**              Chapront-Touze, M., Francou, G., and Laskar, J.,
**              Astron. Astrophys. 282, 663 (1994).
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
/* Gaussian constant */
   static const double GK = 0.017202098950;

/* Sin and cos of J2000.0 mean obliquity (IAU 1976) */
   static const double SINEPS = 0.3977771559319137;
   static const double COSEPS = 0.9174820620691818;

/* Maximum number of iterations allowed to solve Kepler's equation */
   static const int KMAX = 10;

   int jstat, i, k;
   double t, da, dl, de, dp, di, dom, dmu, arga, argl, am,
          ae, dae, ae2, at, r, v, si2, xq, xp, tl, xsw,
          xcw, xm2, xf, ci2, xms, xmc, xpxq2, x, y, z;

/* Planetary inverse masses */
   static const double amas[] = { 6023600.0,       /* Mercury */
                                   408523.5,       /* Venus   */
                                   328900.5,       /* EMB     */
                                  3098710.0,       /* Mars    */
                                     1047.355,     /* Jupiter */
                                     3498.5,       /* Saturn  */
                                    22869.0,       /* Uranus  */
                                    19314.0 };     /* Neptune */

/*
** Tables giving the mean Keplerian elements, limited to t^2 terms:
**
**   a       semi-major axis (AU)
**   dlm     mean longitude (degree and arcsecond)
**   e       eccentricity
**   pi      longitude of the perihelion (degree and arcsecond)
**   dinc    inclination (degree and arcsecond)
**   omega   longitude of the ascending node (degree and arcsecond)
*/

   static const double a[][3] = {
       {  0.3870983098,           0.0,     0.0 },  /* Mercury */
       {  0.7233298200,           0.0,     0.0 },  /* Venus   */
       {  1.0000010178,           0.0,     0.0 },  /* EMB     */
       {  1.5236793419,         3e-10,     0.0 },  /* Mars    */
       {  5.2026032092,     19132e-10, -39e-10 },  /* Jupiter */
       {  9.5549091915, -0.0000213896, 444e-10 },  /* Saturn  */
       { 19.2184460618,     -3716e-10, 979e-10 },  /* Uranus  */
       { 30.1103868694,    -16635e-10, 686e-10 }   /* Neptune */
   };

   static const double dlm[][3] = {
       { 252.25090552, 5381016286.88982,  -1.92789 },
       { 181.97980085, 2106641364.33548,   0.59381 },
       { 100.46645683, 1295977422.83429,  -2.04411 },
       { 355.43299958,  689050774.93988,   0.94264 },
       {  34.35151874,  109256603.77991, -30.60378 },
       {  50.07744430,   43996098.55732,  75.61614 },
       { 314.05500511,   15424811.93933,  -1.75083 },
       { 304.34866548,    7865503.20744,   0.21103 }
   };

   static const double e[][3] = {
       { 0.2056317526,  0.0002040653,    -28349e-10 },
       { 0.0067719164, -0.0004776521,     98127e-10 },
       { 0.0167086342, -0.0004203654, -0.0000126734 },
       { 0.0934006477,  0.0009048438,    -80641e-10 },
       { 0.0484979255,  0.0016322542, -0.0000471366 },
       { 0.0555481426, -0.0034664062, -0.0000643639 },
       { 0.0463812221, -0.0002729293,  0.0000078913 },
       { 0.0094557470,  0.0000603263,           0.0 }
   };

   static const double pi[][3] = {
       {  77.45611904,  5719.11590,   -4.83016 },
       { 131.56370300,   175.48640, -498.48184 },
       { 102.93734808, 11612.35290,   53.27577 },
       { 336.06023395, 15980.45908,  -62.32800 },
       {  14.33120687,  7758.75163,  259.95938 },
       {  93.05723748, 20395.49439,  190.25952 },
       { 173.00529106,  3215.56238,  -34.09288 },
       {  48.12027554,  1050.71912,   27.39717 }
   };

   static const double dinc[][3] = {
       { 7.00498625, -214.25629,   0.28977 },
       { 3.39466189,  -30.84437, -11.67836 },
       {        0.0,  469.97289,  -3.35053 },
       { 1.84972648, -293.31722,  -8.11830 },
       { 1.30326698,  -71.55890,  11.95297 },
       { 2.48887878,   91.85195, -17.66225 },
       { 0.77319689,  -60.72723,   1.25759 },
       { 1.76995259,    8.12333,   0.08135 }
   };

   static const double omega[][3] = {
       {  48.33089304,  -4515.21727,  -31.79892 },
       {  76.67992019, -10008.48154,  -51.32614 },
       { 174.87317577,  -8679.27034,   15.34191 },
       {  49.55809321, -10620.90088, -230.57416 },
       { 100.46440702,   6362.03561,  326.52178 },
       { 113.66550252,  -9240.19942,  -66.23743 },
       {  74.00595701,   2669.15033,  145.93964 },
       { 131.78405702,   -221.94322,   -0.78728 }
   };

/* Tables for trigonometric terms to be added to the mean elements of */
/* the semi-major axes */

   static const double kp[][9] = {
    {   69613, 75645, 88306, 59899, 15746, 71087, 142173,  3086,    0 },
    {   21863, 32794, 26934, 10931, 26250, 43725,  53867, 28939,    0 },
    {   16002, 21863, 32004, 10931, 14529, 16368,  15318, 32794,    0 },
    {    6345,  7818, 15636,  7077,  8184, 14163,   1107,  4872,    0 },
    {    1760,  1454,  1167,   880,   287,  2640,     19,  2047, 1454 },
    {     574,     0,   880,   287,    19,  1760,   1167,   306,  574 },
    {     204,     0,   177,  1265,     4,   385,    200,   208,  204 },
    {       0,   102,   106,     4,    98,  1367,    487,   204,    0 }
   };

   static const double ca[][9] = {
    {       4,    -13,    11,   -9,    -9,   -3,     -1,     4,     0 },
    {    -156,     59,   -42,    6,    19,  -20,    -10,   -12,     0 },
    {      64,   -152,    62,   -8,    32,  -41,     19,   -11,     0 },
    {     124,    621,  -145,  208,    54,  -57,     30,    15,     0 },
    {  -23437,  -2634,  6601, 6259, -1507,-1821,   2620, -2115, -1489 },
    {   62911,-119919, 79336,17814,-24241,12068,   8306, -4893,  8902 },
    {  389061,-262125,-44088, 8387,-22976,-2093,   -615, -9720,  6633 },
    { -412235,-157046,-31430,37817, -9740,  -13,  -7449,  9644,     0 }
   };

   static const double sa[][9] = {
    {     -29,    -1,     9,     6,    -6,     5,     4,     0,     0 },
    {     -48,  -125,   -26,   -37,    18,   -13,   -20,    -2,     0 },
    {    -150,   -46,    68,    54,    14,    24,   -28,    22,     0 },
    {    -621,   532,  -694,   -20,   192,   -94,    71,   -73,     0 },
    {  -14614,-19828, -5869,  1881, -4372, -2255,   782,   930,   913 },
    {  139737,     0, 24667, 51123, -5102,  7429, -4095, -1976, -9566 },
    { -138081,     0, 37205,-49039,-41901,-33872,-27037,-12474, 18797 },
    {       0, 28492,133236, 69654, 52322,-49577,-26430, -3593,     0 }
   };

/* Tables giving the trigonometric terms to be added to the mean */
/* elements of the mean longitudes */

   static const double kq[][10] = {
    {   3086,15746,69613,59899,75645,88306, 12661,  2658,    0,     0 },
    {  21863,32794,10931,   73, 4387,26934,  1473,  2157,    0,     0 },
    {     10,16002,21863,10931, 1473,32004,  4387,    73,    0,     0 },
    {     10, 6345, 7818, 1107,15636, 7077,  8184,   532,   10,     0 },
    {     19, 1760, 1454,  287, 1167,  880,   574,  2640,   19,  1454 },
    {     19,  574,  287,  306, 1760,   12,    31,    38,   19,   574 },
    {      4,  204,  177,    8,   31,  200,  1265,   102,    4,   204 },
    {      4,  102,  106,    8,   98, 1367,   487,   204,    4,   102 }
   };

   static const double cl[][10] = {
    {      21,   -95, -157,   41,   -5,   42,  23,  30,      0,     0 },
    {    -160,  -313, -235,   60,  -74,  -76, -27,  34,      0,     0 },
    {    -325,  -322,  -79,  232,  -52,   97,  55, -41,      0,     0 },
    {    2268,  -979,  802,  602, -668,  -33, 345, 201,    -55,     0 },
    {    7610, -4997,-7689,-5841,-2617, 1115,-748,-607,   6074,   354 },
    {  -18549, 30125,20012, -730,  824,   23,1289,-352, -14767, -2062 },
    { -135245,-14594, 4197,-4030,-5630,-2898,2540,-306,   2939,  1986 },
    {   89948,  2103, 8963, 2695, 3682, 1648, 866,-154,  -1963,  -283 }
   };

   static const double sl[][10] = {
    {   -342,   136,  -23,   62,   66,  -52, -33,    17,     0,     0 },
    {    524,  -149,  -35,  117,  151,  122, -71,   -62,     0,     0 },
    {   -105,  -137,  258,   35, -116,  -88,-112,   -80,     0,     0 },
    {    854,  -205, -936, -240,  140, -341, -97,  -232,   536,     0 },
    { -56980,  8016, 1012, 1448,-3024,-3710, 318,   503,  3767,   577 },
    { 138606,-13478,-4964, 1441,-1319,-1482, 427,  1236, -9167, -1918 },
    {  71234,-41116, 5334,-4935,-1848,   66, 434, -1748,  3780,  -701 },
    { -47645, 11647, 2166, 3194,  679,    0,-244,  -419, -2531,    48 }
   };

/*--------------------------------------------------------------------*/

/* Validate the planet number. */
   if ((np < 1) || (np > 8)) {
      jstat = -1;

   /* Reset the result in case of failure. */
      for (k = 0; k < 2; k++) {
         for (i = 0; i < 3; i++) {
            pv[k][i] = 0.0;
         }
      }

   } else {

   /* Decrement the planet number to start at zero. */
      np--;

   /* Time: Julian millennia since J2000.0. */
      t = ((date1 - DJ00) + date2) / DJM;

   /* OK status unless remote date. */
      jstat = fabs(t) <= 1.0 ? 0 : 1;

   /* Compute the mean elements. */
      da = a[np][0] +
          (a[np][1] +
           a[np][2] * t) * t;
      dl = (3600.0 * dlm[np][0] +
                    (dlm[np][1] +
                     dlm[np][2] * t) * t) * DAS2R;
      de = e[np][0] +
         ( e[np][1] +
           e[np][2] * t) * t;
      dp = eraAnpm((3600.0 * pi[np][0] +
                            (pi[np][1] +
                             pi[np][2] * t) * t) * DAS2R);
      di = (3600.0 * dinc[np][0] +
                    (dinc[np][1] +
                     dinc[np][2] * t) * t) * DAS2R;
      dom = eraAnpm((3600.0 * omega[np][0] +
                             (omega[np][1] +
                              omega[np][2] * t) * t) * DAS2R);

   /* Apply the trigonometric terms. */
      dmu = 0.35953620 * t;
      for (k = 0; k < 8; k++) {
         arga = kp[np][k] * dmu;
         argl = kq[np][k] * dmu;
         da += (ca[np][k] * cos(arga) +
                sa[np][k] * sin(arga)) * 1e-7;
         dl += (cl[np][k] * cos(argl) +
                sl[np][k] * sin(argl)) * 1e-7;
      }
      arga = kp[np][8] * dmu;
      da += t * (ca[np][8] * cos(arga) +
                 sa[np][8] * sin(arga)) * 1e-7;
      for (k = 8; k < 10; k++) {
         argl = kq[np][k] * dmu;
         dl += t * (cl[np][k] * cos(argl) +
                    sl[np][k] * sin(argl)) * 1e-7;
      }
      dl = fmod(dl, D2PI);

   /* Iterative soln. of Kepler's equation to get eccentric anomaly. */
      am = dl - dp;
      ae = am + de * sin(am);
      k = 0;
      dae = 1.0;
      while (k < KMAX && fabs(dae) > 1e-12) {
         dae = (am - ae + de * sin(ae)) / (1.0 - de * cos(ae));
         ae += dae;
         k++;
         if (k == KMAX-1) jstat = 2;
      }

   /* True anomaly. */
      ae2 = ae / 2.0;
      at = 2.0 * atan2(sqrt((1.0 + de) / (1.0 - de)) * sin(ae2),
                                                       cos(ae2));

   /* Distance (AU) and speed (radians per day). */
      r = da * (1.0 - de * cos(ae));
      v = GK * sqrt((1.0 + 1.0 / amas[np]) / (da * da * da));

      si2 = sin(di / 2.0);
      xq = si2 * cos(dom);
      xp = si2 * sin(dom);
      tl = at + dp;
      xsw = sin(tl);
      xcw = cos(tl);
      xm2 = 2.0 * (xp * xcw - xq * xsw);
      xf = da / sqrt(1  -  de * de);
      ci2 = cos(di / 2.0);
      xms = (de * sin(dp) + xsw) * xf;
      xmc = (de * cos(dp) + xcw) * xf;
      xpxq2 = 2 * xp * xq;

   /* Position (J2000.0 ecliptic x,y,z in AU). */
      x = r * (xcw - xm2 * xp);
      y = r * (xsw + xm2 * xq);
      z = r * (-xm2 * ci2);

   /* Rotate to equatorial. */
      pv[0][0] = x;
      pv[0][1] = y * COSEPS - z * SINEPS;
      pv[0][2] = y * SINEPS + z * COSEPS;

   /* Velocity (J2000.0 ecliptic xdot,ydot,zdot in AU/d). */
      x = v * (( -1.0 + 2.0 * xp * xp) * xms + xpxq2 * xmc);
      y = v * ((  1.0 - 2.0 * xq * xq) * xmc - xpxq2 * xms);
      z = v * (2.0 * ci2 * (xp * xms + xq * xmc));

   /* Rotate to equatorial. */
      pv[1][0] = x;
      pv[1][1] = y * COSEPS - z * SINEPS;
      pv[1][2] = y * SINEPS + z * COSEPS;

   }

/* Return the status. */
   return jstat;

}

double eraPm(double p[3])
/*
**  - - - - - -
**   e r a P m
**  - - - - - -
**
**  Modulus of p-vector.
**
**  Given:
**     p      double[3]     p-vector
**
**  Returned (function value):
**            double        modulus
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double w;


   w  = sqrt( p[0] * p[0]
            + p[1] * p[1]
            + p[2] * p[2] );

   return w;

}

void eraPmat00(double date1, double date2, double rbp[3][3])
/*
**  - - - - - - - - - -
**   e r a P m a t 0 0
**  - - - - - - - - - -
**
**  Precession matrix (including frame bias) from GCRS to a specified
**  date, IAU 2000 model.
**
**  Given:
**     date1,date2  double          TT as a 2-part Julian Date (Note 1)
**
**  Returned:
**     rbp          double[3][3]    bias-precession matrix (Note 2)
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
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
**  2) The matrix operates in the sense V(date) = rbp * V(GCRS), where
**     the p-vector V(GCRS) is with respect to the Geocentric Celestial
**     Reference System (IAU, 2000) and the p-vector V(date) is with
**     respect to the mean equatorial triad of the given date.
**
**  Called:
**     eraBp00      frame bias and precession matrices, IAU 2000
**
**  Reference:
**
**     IAU: Trans. International Astronomical Union, Vol. XXIVB;  Proc.
**     24th General Assembly, Manchester, UK.  Resolutions B1.3, B1.6.
**     (2000)
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double rb[3][3], rp[3][3];


/* Obtain the required matrix (discarding others). */
   eraBp00(date1, date2, rb, rp, rbp);

   return;

}

void eraPmat06(double date1, double date2, double rbp[3][3])
/*
**  - - - - - - - - - -
**   e r a P m a t 0 6
**  - - - - - - - - - -
**
**  Precession matrix (including frame bias) from GCRS to a specified
**  date, IAU 2006 model.
**
**  Given:
**     date1,date2  double          TT as a 2-part Julian Date (Note 1)
**
**  Returned:
**     rbp          double[3][3]    bias-precession matrix (Note 2)
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
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
**  2) The matrix operates in the sense V(date) = rbp * V(GCRS), where
**     the p-vector V(GCRS) is with respect to the Geocentric Celestial
**     Reference System (IAU, 2000) and the p-vector V(date) is with
**     respect to the mean equatorial triad of the given date.
**
**  Called:
**     eraPfw06     bias-precession F-W angles, IAU 2006
**     eraFw2m      F-W angles to r-matrix
**
**  References:
**
**     Capitaine, N. & Wallace, P.T., 2006, Astron.Astrophys. 450, 855
**
**     Wallace, P.T. & Capitaine, N., 2006, Astron.Astrophys. 459, 981
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double gamb, phib, psib, epsa;


/* Bias-precession Fukushima-Williams angles. */
   eraPfw06(date1, date2, &gamb, &phib, &psib, &epsa);

/* Form the matrix. */
   eraFw2m(gamb, phib, psib, epsa, rbp);

   return;

}

void eraPmat76(double date1, double date2, double rmatp[3][3])
/*
**  - - - - - - - - - -
**   e r a P m a t 7 6
**  - - - - - - - - - -
**
**  Precession matrix from J2000.0 to a specified date, IAU 1976 model.
**
**  Given:
**     date1,date2 double       ending date, TT (Note 1)
**
**  Returned:
**     rmatp       double[3][3] precession matrix, J2000.0 -> date1+date2
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
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
**  2) The matrix operates in the sense V(date) = RMATP * V(J2000),
**     where the p-vector V(J2000) is with respect to the mean
**     equatorial triad of epoch J2000.0 and the p-vector V(date)
**     is with respect to the mean equatorial triad of the given
**     date.
**
**  3) Though the matrix method itself is rigorous, the precession
**     angles are expressed through canonical polynomials which are
**     valid only for a limited time span.  In addition, the IAU 1976
**     precession rate is known to be imperfect.  The absolute accuracy
**     of the present formulation is better than 0.1 arcsec from
**     1960AD to 2040AD, better than 1 arcsec from 1640AD to 2360AD,
**     and remains below 3 arcsec for the whole of the period
**     500BC to 3000AD.  The errors exceed 10 arcsec outside the
**     range 1200BC to 3900AD, exceed 100 arcsec outside 4200BC to
**     5600AD and exceed 1000 arcsec outside 6800BC to 8200AD.
**
**  Called:
**     eraPrec76    accumulated precession angles, IAU 1976
**     eraIr        initialize r-matrix to identity
**     eraRz        rotate around Z-axis
**     eraRy        rotate around Y-axis
**     eraCr        copy r-matrix
**
**  References:
**
**     Lieske, J.H., 1979, Astron.Astrophys. 73, 282.
**      equations (6) & (7), p283.
**
**     Kaplan,G.H., 1981. USNO circular no. 163, pA2.
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double zeta, z, theta, wmat[3][3];


/* Precession Euler angles, J2000.0 to specified date. */
   eraPrec76(DJ00, 0.0, date1, date2, &zeta, &z, &theta);

/* Form the rotation matrix. */
   eraIr(  wmat);
   eraRz( -zeta, wmat);
   eraRy(  theta, wmat);
   eraRz( -z, wmat);
   eraCr( wmat, rmatp);

   return;

}

void eraPmp(double a[3], double b[3], double amb[3])
/*
**  - - - - - - -
**   e r a P m p
**  - - - - - - -
**
**  P-vector subtraction.
**
**  Given:
**     a        double[3]      first p-vector
**     b        double[3]      second p-vector
**
**  Returned:
**     amb      double[3]      a - b
**
**  Note:
**     It is permissible to re-use the same array for any of the
**     arguments.
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   amb[0] = a[0] - b[0];
   amb[1] = a[1] - b[1];
   amb[2] = a[2] - b[2];

   return;

}

void eraPn(double p[3], double *r, double u[3])
/*
**  - - - - - -
**   e r a P n
**  - - - - - -
**
**  Convert a p-vector into modulus and unit vector.
**
**  Given:
**     p        double[3]      p-vector
**
**  Returned:
**     r        double         modulus
**     u        double[3]      unit vector
**
**  Notes:
**
**  1) If p is null, the result is null.  Otherwise the result is a unit
**     vector.
**
**  2) It is permissible to re-use the same array for any of the
**     arguments.
**
**  Called:
**     eraPm        modulus of p-vector
**     eraZp        zero p-vector
**     eraSxp       multiply p-vector by scalar
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double w;


/* Obtain the modulus and test for zero. */
   w = eraPm(p);
   if (w == 0.0) {

   /* Null vector. */
      eraZp(u);

   } else {

   /* Unit vector. */
      eraSxp(1.0/w, p, u);
   }

/* Return the modulus. */
   *r = w;

   return;

}

void eraPn00(double date1, double date2, double dpsi, double deps,
             double *epsa,
             double rb[3][3], double rp[3][3], double rbp[3][3],
             double rn[3][3], double rbpn[3][3])
/*
**  - - - - - - - -
**   e r a P n 0 0
**  - - - - - - - -
**
**  Precession-nutation, IAU 2000 model:  a multi-purpose function,
**  supporting classical (equinox-based) use directly and CIO-based
**  use indirectly.
**
**  Given:
**     date1,date2  double          TT as a 2-part Julian Date (Note 1)
**     dpsi,deps    double          nutation (Note 2)
**
**  Returned:
**     epsa         double          mean obliquity (Note 3)
**     rb           double[3][3]    frame bias matrix (Note 4)
**     rp           double[3][3]    precession matrix (Note 5)
**     rbp          double[3][3]    bias-precession matrix (Note 6)
**     rn           double[3][3]    nutation matrix (Note 7)
**     rbpn         double[3][3]    GCRS-to-true matrix (Note 8)
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
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
**  2) The caller is responsible for providing the nutation components;
**     they are in longitude and obliquity, in radians and are with
**     respect to the equinox and ecliptic of date.  For high-accuracy
**     applications, free core nutation should be included as well as
**     any other relevant corrections to the position of the CIP.
**
**  3) The returned mean obliquity is consistent with the IAU 2000
**     precession-nutation models.
**
**  4) The matrix rb transforms vectors from GCRS to J2000.0 mean
**     equator and equinox by applying frame bias.
**
**  5) The matrix rp transforms vectors from J2000.0 mean equator and
**     equinox to mean equator and equinox of date by applying
**     precession.
**
**  6) The matrix rbp transforms vectors from GCRS to mean equator and
**     equinox of date by applying frame bias then precession.  It is
**     the product rp x rb.
**
**  7) The matrix rn transforms vectors from mean equator and equinox of
**     date to true equator and equinox of date by applying the nutation
**     (luni-solar + planetary).
**
**  8) The matrix rbpn transforms vectors from GCRS to true equator and
**     equinox of date.  It is the product rn x rbp, applying frame
**     bias, precession and nutation in that order.
**
**  9) It is permissible to re-use the same array in the returned
**     arguments.  The arrays are filled in the order given.
**
**  Called:
**     eraPr00      IAU 2000 precession adjustments
**     eraObl80     mean obliquity, IAU 1980
**     eraBp00      frame bias and precession matrices, IAU 2000
**     eraCr        copy r-matrix
**     eraNumat     form nutation matrix
**     eraRxr       product of two r-matrices
**
**  Reference:
**
**     Capitaine, N., Chapront, J., Lambert, S. and Wallace, P.,
**     "Expressions for the Celestial Intermediate Pole and Celestial
**     Ephemeris Origin consistent with the IAU 2000A precession-
**     nutation model", Astron.Astrophys. 400, 1145-1154 (2003)
**
**     n.b. The celestial ephemeris origin (CEO) was renamed "celestial
**          intermediate origin" (CIO) by IAU 2006 Resolution 2.
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double dpsipr, depspr, rbpw[3][3], rnw[3][3];


/* IAU 2000 precession-rate adjustments. */
   eraPr00(date1, date2, &dpsipr, &depspr);

/* Mean obliquity, consistent with IAU 2000 precession-nutation. */
   *epsa = eraObl80(date1, date2) + depspr;

/* Frame bias and precession matrices and their product. */
   eraBp00(date1, date2, rb, rp, rbpw);
   eraCr(rbpw, rbp);

/* Nutation matrix. */
   eraNumat(*epsa, dpsi, deps, rnw);
   eraCr(rnw, rn);

/* Bias-precession-nutation matrix (classical). */
   eraRxr(rnw, rbpw, rbpn);

   return;

}

void eraPn00a(double date1, double date2,
              double *dpsi, double *deps, double *epsa,
              double rb[3][3], double rp[3][3], double rbp[3][3],
              double rn[3][3], double rbpn[3][3])
/*
**  - - - - - - - - -
**   e r a P n 0 0 a
**  - - - - - - - - -
**
**  Precession-nutation, IAU 2000A model:  a multi-purpose function,
**  supporting classical (equinox-based) use directly and CIO-based
**  use indirectly.
**
**  Given:
**     date1,date2  double          TT as a 2-part Julian Date (Note 1)
**
**  Returned:
**     dpsi,deps    double          nutation (Note 2)
**     epsa         double          mean obliquity (Note 3)
**     rb           double[3][3]    frame bias matrix (Note 4)
**     rp           double[3][3]    precession matrix (Note 5)
**     rbp          double[3][3]    bias-precession matrix (Note 6)
**     rn           double[3][3]    nutation matrix (Note 7)
**     rbpn         double[3][3]    GCRS-to-true matrix (Notes 8,9)
**
**  Notes:
**
**  1)  The TT date date1+date2 is a Julian Date, apportioned in any
**      convenient way between the two arguments.  For example,
**      JD(TT)=2450123.7 could be expressed in any of these ways,
**      among others:
**
**             date1          date2
**
**          2450123.7           0.0       (JD method)
**          2451545.0       -1421.3       (J2000 method)
**          2400000.5       50123.2       (MJD method)
**          2450123.5           0.2       (date & time method)
**
**      The JD method is the most natural and convenient to use in
**      cases where the loss of several decimal digits of resolution
**      is acceptable.  The J2000 method is best matched to the way
**      the argument is handled internally and will deliver the
**      optimum resolution.  The MJD method and the date & time methods
**      are both good compromises between resolution and convenience.
**
**  2)  The nutation components (luni-solar + planetary, IAU 2000A) in
**      longitude and obliquity are in radians and with respect to the
**      equinox and ecliptic of date.  Free core nutation is omitted;
**      for the utmost accuracy, use the eraPn00  function, where the
**      nutation components are caller-specified.  For faster but
**      slightly less accurate results, use the eraPn00b function.
**
**  3)  The mean obliquity is consistent with the IAU 2000 precession.
**
**  4)  The matrix rb transforms vectors from GCRS to J2000.0 mean
**      equator and equinox by applying frame bias.
**
**  5)  The matrix rp transforms vectors from J2000.0 mean equator and
**      equinox to mean equator and equinox of date by applying
**      precession.
**
**  6)  The matrix rbp transforms vectors from GCRS to mean equator and
**      equinox of date by applying frame bias then precession.  It is
**      the product rp x rb.
**
**  7)  The matrix rn transforms vectors from mean equator and equinox
**      of date to true equator and equinox of date by applying the
**      nutation (luni-solar + planetary).
**
**  8)  The matrix rbpn transforms vectors from GCRS to true equator and
**      equinox of date.  It is the product rn x rbp, applying frame
**      bias, precession and nutation in that order.
**
**  9)  The X,Y,Z coordinates of the IAU 2000B Celestial Intermediate
**      Pole are elements (3,1-3) of the matrix rbpn.
**
**  10) It is permissible to re-use the same array in the returned
**      arguments.  The arrays are filled in the order given.
**
**  Called:
**     eraNut00a    nutation, IAU 2000A
**     eraPn00      bias/precession/nutation results, IAU 2000
**
**  Reference:
**
**     Capitaine, N., Chapront, J., Lambert, S. and Wallace, P.,
**     "Expressions for the Celestial Intermediate Pole and Celestial
**     Ephemeris Origin consistent with the IAU 2000A precession-
**     nutation model", Astron.Astrophys. 400, 1145-1154 (2003)
**
**     n.b. The celestial ephemeris origin (CEO) was renamed "celestial
**          intermediate origin" (CIO) by IAU 2006 Resolution 2.
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
/* Nutation. */
   eraNut00a(date1, date2, dpsi, deps);

/* Remaining results. */
   eraPn00(date1, date2, *dpsi, *deps, epsa, rb, rp, rbp, rn, rbpn);

   return;

}

void eraPn00b(double date1, double date2,
              double *dpsi, double *deps, double *epsa,
              double rb[3][3], double rp[3][3], double rbp[3][3],
              double rn[3][3], double rbpn[3][3])
/*
**  - - - - - - - - -
**   e r a P n 0 0 b
**  - - - - - - - - -
**
**  Precession-nutation, IAU 2000B model:  a multi-purpose function,
**  supporting classical (equinox-based) use directly and CIO-based
**  use indirectly.
**
**  Given:
**     date1,date2  double          TT as a 2-part Julian Date (Note 1)
**
**  Returned:
**     dpsi,deps    double          nutation (Note 2)
**     epsa         double          mean obliquity (Note 3)
**     rb           double[3][3]    frame bias matrix (Note 4)
**     rp           double[3][3]    precession matrix (Note 5)
**     rbp          double[3][3]    bias-precession matrix (Note 6)
**     rn           double[3][3]    nutation matrix (Note 7)
**     rbpn         double[3][3]    GCRS-to-true matrix (Notes 8,9)
**
**  Notes:
**
**  1)  The TT date date1+date2 is a Julian Date, apportioned in any
**      convenient way between the two arguments.  For example,
**      JD(TT)=2450123.7 could be expressed in any of these ways,
**      among others:
**
**             date1          date2
**
**          2450123.7           0.0       (JD method)
**          2451545.0       -1421.3       (J2000 method)
**          2400000.5       50123.2       (MJD method)
**          2450123.5           0.2       (date & time method)
**
**      The JD method is the most natural and convenient to use in
**      cases where the loss of several decimal digits of resolution
**      is acceptable.  The J2000 method is best matched to the way
**      the argument is handled internally and will deliver the
**      optimum resolution.  The MJD method and the date & time methods
**      are both good compromises between resolution and convenience.
**
**  2)  The nutation components (luni-solar + planetary, IAU 2000B) in
**      longitude and obliquity are in radians and with respect to the
**      equinox and ecliptic of date.  For more accurate results, but
**      at the cost of increased computation, use the eraPn00a function.
**      For the utmost accuracy, use the eraPn00  function, where the
**      nutation components are caller-specified.
**
**  3)  The mean obliquity is consistent with the IAU 2000 precession.
**
**  4)  The matrix rb transforms vectors from GCRS to J2000.0 mean
**      equator and equinox by applying frame bias.
**
**  5)  The matrix rp transforms vectors from J2000.0 mean equator and
**      equinox to mean equator and equinox of date by applying
**      precession.
**
**  6)  The matrix rbp transforms vectors from GCRS to mean equator and
**      equinox of date by applying frame bias then precession.  It is
**      the product rp x rb.
**
**  7)  The matrix rn transforms vectors from mean equator and equinox
**      of date to true equator and equinox of date by applying the
**      nutation (luni-solar + planetary).
**
**  8)  The matrix rbpn transforms vectors from GCRS to true equator and
**      equinox of date.  It is the product rn x rbp, applying frame
**      bias, precession and nutation in that order.
**
**  9)  The X,Y,Z coordinates of the IAU 2000B Celestial Intermediate
**      Pole are elements (3,1-3) of the matrix rbpn.
**
**  10) It is permissible to re-use the same array in the returned
**      arguments.  The arrays are filled in the stated order.
**
**  Called:
**     eraNut00b    nutation, IAU 2000B
**     eraPn00      bias/precession/nutation results, IAU 2000
**
**  Reference:
**
**     Capitaine, N., Chapront, J., Lambert, S. and Wallace, P.,
**     "Expressions for the Celestial Intermediate Pole and Celestial
**     Ephemeris Origin consistent with the IAU 2000A precession-
**     nutation model", Astron.Astrophys. 400, 1145-1154 (2003).
**
**     n.b. The celestial ephemeris origin (CEO) was renamed "celestial
**          intermediate origin" (CIO) by IAU 2006 Resolution 2.
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
/* Nutation. */
   eraNut00b(date1, date2, dpsi, deps);

/* Remaining results. */
   eraPn00(date1, date2, *dpsi, *deps, epsa, rb, rp, rbp, rn, rbpn);

   return;

}

void eraPn06(double date1, double date2, double dpsi, double deps,
             double *epsa,
             double rb[3][3], double rp[3][3], double rbp[3][3],
             double rn[3][3], double rbpn[3][3])
/*
**  - - - - - - - -
**   e r a P n 0 6
**  - - - - - - - -
**
**  Precession-nutation, IAU 2006 model:  a multi-purpose function,
**  supporting classical (equinox-based) use directly and CIO-based use
**  indirectly.
**
**  Given:
**     date1,date2  double          TT as a 2-part Julian Date (Note 1)
**     dpsi,deps    double          nutation (Note 2)
**
**  Returned:
**     epsa         double          mean obliquity (Note 3)
**     rb           double[3][3]    frame bias matrix (Note 4)
**     rp           double[3][3]    precession matrix (Note 5)
**     rbp          double[3][3]    bias-precession matrix (Note 6)
**     rn           double[3][3]    nutation matrix (Note 7)
**     rbpn         double[3][3]    GCRS-to-true matrix (Note 8)
**
**  Notes:
**
**  1)  The TT date date1+date2 is a Julian Date, apportioned in any
**      convenient way between the two arguments.  For example,
**      JD(TT)=2450123.7 could be expressed in any of these ways,
**      among others:
**
**             date1          date2
**
**          2450123.7           0.0       (JD method)
**          2451545.0       -1421.3       (J2000 method)
**          2400000.5       50123.2       (MJD method)
**          2450123.5           0.2       (date & time method)
**
**      The JD method is the most natural and convenient to use in
**      cases where the loss of several decimal digits of resolution
**      is acceptable.  The J2000 method is best matched to the way
**      the argument is handled internally and will deliver the
**      optimum resolution.  The MJD method and the date & time methods
**      are both good compromises between resolution and convenience.
**
**  2)  The caller is responsible for providing the nutation components;
**      they are in longitude and obliquity, in radians and are with
**      respect to the equinox and ecliptic of date.  For high-accuracy
**      applications, free core nutation should be included as well as
**      any other relevant corrections to the position of the CIP.
**
**  3)  The returned mean obliquity is consistent with the IAU 2006
**      precession.
**
**  4)  The matrix rb transforms vectors from GCRS to J2000.0 mean
**      equator and equinox by applying frame bias.
**
**  5)  The matrix rp transforms vectors from J2000.0 mean equator and
**      equinox to mean equator and equinox of date by applying
**      precession.
**
**  6)  The matrix rbp transforms vectors from GCRS to mean equator and
**      equinox of date by applying frame bias then precession.  It is
**      the product rp x rb.
**
**  7)  The matrix rn transforms vectors from mean equator and equinox
**      of date to true equator and equinox of date by applying the
**      nutation (luni-solar + planetary).
**
**  8)  The matrix rbpn transforms vectors from GCRS to true equator and
**      equinox of date.  It is the product rn x rbp, applying frame
**      bias, precession and nutation in that order.
**
**  9)  The X,Y,Z coordinates of the IAU 2000B Celestial Intermediate
**      Pole are elements (3,1-3) of the matrix rbpn.
**
**  10) It is permissible to re-use the same array in the returned
**      arguments.  The arrays are filled in the stated order.
**
**  Called:
**     eraPfw06     bias-precession F-W angles, IAU 2006
**     eraFw2m      F-W angles to r-matrix
**     eraCr        copy r-matrix
**     eraTr        transpose r-matrix
**     eraRxr       product of two r-matrices
**
**  References:
**
**     Capitaine, N. & Wallace, P.T., 2006, Astron.Astrophys. 450, 855
**
**     Wallace, P.T. & Capitaine, N., 2006, Astron.Astrophys. 459, 981
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double gamb, phib, psib, eps, r1[3][3], r2[3][3], rt[3][3];


/* Bias-precession Fukushima-Williams angles of J2000.0 = frame bias. */
   eraPfw06(DJM0, DJM00, &gamb, &phib, &psib, &eps);

/* B matrix. */
   eraFw2m(gamb, phib, psib, eps, r1);
   eraCr(r1, rb);

/* Bias-precession Fukushima-Williams angles of date. */
   eraPfw06(date1, date2, &gamb, &phib, &psib, &eps);

/* Bias-precession matrix. */
   eraFw2m(gamb, phib, psib, eps, r2);
   eraCr(r2, rbp);

/* Solve for precession matrix. */
   eraTr(r1, rt);
   eraRxr(r2, rt, rp);

/* Equinox-based bias-precession-nutation matrix. */
   eraFw2m(gamb, phib, psib + dpsi, eps + deps, r1);
   eraCr(r1, rbpn);

/* Solve for nutation matrix. */
   eraTr(r2, rt);
   eraRxr(r1, rt, rn);

/* Obliquity, mean of date. */
   *epsa = eps;

   return;

}

void eraPn06a(double date1, double date2,
              double *dpsi, double *deps, double *epsa,
              double rb[3][3], double rp[3][3], double rbp[3][3],
              double rn[3][3], double rbpn[3][3])
/*
**  - - - - - - - - -
**   e r a P n 0 6 a
**  - - - - - - - - -
**
**  Precession-nutation, IAU 2006/2000A models:  a multi-purpose function,
**  supporting classical (equinox-based) use directly and CIO-based use
**  indirectly.
**
**  Given:
**     date1,date2  double          TT as a 2-part Julian Date (Note 1)
**
**  Returned:
**     dpsi,deps    double          nutation (Note 2)
**     epsa         double          mean obliquity (Note 3)
**     rb           double[3][3]    frame bias matrix (Note 4)
**     rp           double[3][3]    precession matrix (Note 5)
**     rbp          double[3][3]    bias-precession matrix (Note 6)
**     rn           double[3][3]    nutation matrix (Note 7)
**     rbpn         double[3][3]    GCRS-to-true matrix (Notes 8,9)
**
**  Notes:
**
**  1)  The TT date date1+date2 is a Julian Date, apportioned in any
**      convenient way between the two arguments.  For example,
**      JD(TT)=2450123.7 could be expressed in any of these ways,
**      among others:
**
**             date1          date2
**
**          2450123.7           0.0       (JD method)
**          2451545.0       -1421.3       (J2000 method)
**          2400000.5       50123.2       (MJD method)
**          2450123.5           0.2       (date & time method)
**
**      The JD method is the most natural and convenient to use in
**      cases where the loss of several decimal digits of resolution
**      is acceptable.  The J2000 method is best matched to the way
**      the argument is handled internally and will deliver the
**      optimum resolution.  The MJD method and the date & time methods
**      are both good compromises between resolution and convenience.
**
**  2)  The nutation components (luni-solar + planetary, IAU 2000A) in
**      longitude and obliquity are in radians and with respect to the
**      equinox and ecliptic of date.  Free core nutation is omitted;
**      for the utmost accuracy, use the eraPn06 function, where the
**      nutation components are caller-specified.
**
**  3)  The mean obliquity is consistent with the IAU 2006 precession.
**
**  4)  The matrix rb transforms vectors from GCRS to mean J2000.0 by
**      applying frame bias.
**
**  5)  The matrix rp transforms vectors from mean J2000.0 to mean of
**      date by applying precession.
**
**  6)  The matrix rbp transforms vectors from GCRS to mean of date by
**      applying frame bias then precession.  It is the product rp x rb.
**
**  7)  The matrix rn transforms vectors from mean of date to true of
**      date by applying the nutation (luni-solar + planetary).
**
**  8)  The matrix rbpn transforms vectors from GCRS to true of date
**      (CIP/equinox).  It is the product rn x rbp, applying frame bias,
**      precession and nutation in that order.
**
**  9)  The X,Y,Z coordinates of the IAU 2006/2000A Celestial
**      Intermediate Pole are elements (1,1-3) of the matrix rbpn.
**
**  10) It is permissible to re-use the same array in the returned
**      arguments.  The arrays are filled in the stated order.
**
**  Called:
**     eraNut06a    nutation, IAU 2006/2000A
**     eraPn06      bias/precession/nutation results, IAU 2006
**
**  Reference:
**
**     Capitaine, N. & Wallace, P.T., 2006, Astron.Astrophys. 450, 855
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
/* Nutation. */
   eraNut06a(date1, date2, dpsi, deps);

/* Remaining results. */
   eraPn06(date1, date2, *dpsi, *deps, epsa, rb, rp, rbp, rn, rbpn);

   return;

}

void eraPnm00a(double date1, double date2, double rbpn[3][3])
/*
**  - - - - - - - - - -
**   e r a P n m 0 0 a
**  - - - - - - - - - -
**
**  Form the matrix of precession-nutation for a given date (including
**  frame bias), equinox-based, IAU 2000A model.
**
**  Given:
**     date1,date2  double     TT as a 2-part Julian Date (Note 1)
**
**  Returned:
**     rbpn         double[3][3]    classical NPB matrix (Note 2)
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
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
**  2) The matrix operates in the sense V(date) = rbpn * V(GCRS), where
**     the p-vector V(date) is with respect to the true equatorial triad
**     of date date1+date2 and the p-vector V(GCRS) is with respect to
**     the Geocentric Celestial Reference System (IAU, 2000).
**
**  3) A faster, but slightly less accurate result (about 1 mas), can be
**     obtained by using instead the eraPnm00b function.
**
**  Called:
**     eraPn00a     bias/precession/nutation, IAU 2000A
**
**  Reference:
**
**     IAU: Trans. International Astronomical Union, Vol. XXIVB;  Proc.
**     24th General Assembly, Manchester, UK.  Resolutions B1.3, B1.6.
**     (2000)
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double dpsi, deps, epsa, rb[3][3], rp[3][3], rbp[3][3], rn[3][3];


/* Obtain the required matrix (discarding other results). */
   eraPn00a(date1, date2, &dpsi, &deps, &epsa, rb, rp, rbp, rn, rbpn);

   return;

}

void eraPnm00b(double date1, double date2, double rbpn[3][3])
/*
**  - - - - - - - - - -
**   e r a P n m 0 0 b
**  - - - - - - - - - -
**
**  Form the matrix of precession-nutation for a given date (including
**  frame bias), equinox-based, IAU 2000B model.
**
**  Given:
**     date1,date2 double       TT as a 2-part Julian Date (Note 1)
**
**  Returned:
**     rbpn        double[3][3] bias-precession-nutation matrix (Note 2)
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
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
**  2) The matrix operates in the sense V(date) = rbpn * V(GCRS), where
**     the p-vector V(date) is with respect to the true equatorial triad
**     of date date1+date2 and the p-vector V(GCRS) is with respect to
**     the Geocentric Celestial Reference System (IAU, 2000).
**
**  3) The present function is faster, but slightly less accurate (about
**     1 mas), than the eraPnm00a function.
**
**  Called:
**     eraPn00b     bias/precession/nutation, IAU 2000B
**
**  Reference:
**
**     IAU: Trans. International Astronomical Union, Vol. XXIVB;  Proc.
**     24th General Assembly, Manchester, UK.  Resolutions B1.3, B1.6.
**     (2000)
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double dpsi, deps, epsa, rb[3][3], rp[3][3], rbp[3][3], rn[3][3];


/* Obtain the required matrix (discarding other results). */
   eraPn00b(date1, date2, &dpsi, &deps, &epsa, rb, rp, rbp, rn, rbpn);

   return;

}

void eraPnm06a(double date1, double date2, double rnpb[3][3])
/*
**  - - - - - - - - - -
**   e r a P n m 0 6 a
**  - - - - - - - - - -
**
**  Form the matrix of precession-nutation for a given date (including
**  frame bias), IAU 2006 precession and IAU 2000A nutation models.
**
**  Given:
**     date1,date2 double       TT as a 2-part Julian Date (Note 1)
**
**  Returned:
**     rnpb        double[3][3] bias-precession-nutation matrix (Note 2)
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
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
**  2) The matrix operates in the sense V(date) = rnpb * V(GCRS), where
**     the p-vector V(date) is with respect to the true equatorial triad
**     of date date1+date2 and the p-vector V(GCRS) is with respect to
**     the Geocentric Celestial Reference System (IAU, 2000).
**
**  Called:
**     eraPfw06     bias-precession F-W angles, IAU 2006
**     eraNut06a    nutation, IAU 2006/2000A
**     eraFw2m      F-W angles to r-matrix
**
**  Reference:
**
**     Capitaine, N. & Wallace, P.T., 2006, Astron.Astrophys. 450, 855.
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double gamb, phib, psib, epsa, dp, de;


/* Fukushima-Williams angles for frame bias and precession. */
   eraPfw06(date1, date2, &gamb, &phib, &psib, &epsa);

/* Nutation components. */
   eraNut06a(date1, date2, &dp, &de);

/* Equinox based nutation x precession x bias matrix. */
   eraFw2m(gamb, phib, psib + dp, epsa + de, rnpb);

   return;

}

void eraPnm80(double date1, double date2, double rmatpn[3][3])
/*
**  - - - - - - - - -
**   e r a P n m 8 0
**  - - - - - - - - -
**
**  Form the matrix of precession/nutation for a given date, IAU 1976
**  precession model, IAU 1980 nutation model.
**
**  Given:
**     date1,date2    double         TDB date (Note 1)
**
**  Returned:
**     rmatpn         double[3][3]   combined precession/nutation matrix
**
**  Notes:
**
**  1) The TDB date date1+date2 is a Julian Date, apportioned in any
**     convenient way between the two arguments.  For example,
**     JD(TDB)=2450123.7 could be expressed in any of these ways,
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
**  2) The matrix operates in the sense V(date) = rmatpn * V(J2000),
**     where the p-vector V(date) is with respect to the true equatorial
**     triad of date date1+date2 and the p-vector V(J2000) is with
**     respect to the mean equatorial triad of epoch J2000.0.
**
**  Called:
**     eraPmat76    precession matrix, IAU 1976
**     eraNutm80    nutation matrix, IAU 1980
**     eraRxr       product of two r-matrices
**
**  Reference:
**
**     Explanatory Supplement to the Astronomical Almanac,
**     P. Kenneth Seidelmann (ed), University Science Books (1992),
**     Section 3.3 (p145).
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double rmatp[3][3], rmatn[3][3];


/* Precession matrix, J2000.0 to date. */
   eraPmat76(date1, date2, rmatp);

/* Nutation matrix. */
   eraNutm80(date1, date2, rmatn);

/* Combine the matrices:  PN = N x P. */
   eraRxr(rmatn, rmatp, rmatpn);

   return;

}

void eraPom00(double xp, double yp, double sp, double rpom[3][3])
/*
**  - - - - - - - - - -
**   e r a P o m 0 0
**  - - - - - - - - - -
**
**  Form the matrix of polar motion for a given date, IAU 2000.
**
**  Given:
**     xp,yp    double    coordinates of the pole (radians, Note 1)
**     sp       double    the TIO locator s' (radians, Note 2)
**
**  Returned:
**     rpom     double[3][3]   polar-motion matrix (Note 3)
**
**  Notes:
**
**  1) The arguments xp and yp are the coordinates (in radians) of the
**     Celestial Intermediate Pole with respect to the International
**     Terrestrial Reference System (see IERS Conventions 2003),
**     measured along the meridians to 0 and 90 deg west respectively.
**
**  2) The argument sp is the TIO locator s', in radians, which
**     positions the Terrestrial Intermediate Origin on the equator.  It
**     is obtained from polar motion observations by numerical
**     integration, and so is in essence unpredictable.  However, it is
**     dominated by a secular drift of about 47 microarcseconds per
**     century, and so can be taken into account by using s' = -47*t,
**     where t is centuries since J2000.0.  The function eraSp00
**     implements this approximation.
**
**  3) The matrix operates in the sense V(TRS) = rpom * V(CIP), meaning
**     that it is the final rotation when computing the pointing
**     direction to a celestial source.
**
**  Called:
**     eraIr        initialize r-matrix to identity
**     eraRz        rotate around Z-axis
**     eraRy        rotate around Y-axis
**     eraRx        rotate around X-axis
**
**  Reference:
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{

/* Construct the matrix. */
   eraIr(rpom);
   eraRz(sp, rpom);
   eraRy(-xp, rpom);
   eraRx(-yp, rpom);

   return;

}

void eraPpp(double a[3], double b[3], double apb[3])
/*
**  - - - - - - -
**   e r a P p p
**  - - - - - - -
**
**  P-vector addition.
**
**  Given:
**     a        double[3]      first p-vector
**     b        double[3]      second p-vector
**
**  Returned:
**     apb      double[3]      a + b
**
**  Note:
**     It is permissible to re-use the same array for any of the
**     arguments.
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   apb[0] = a[0] + b[0];
   apb[1] = a[1] + b[1];
   apb[2] = a[2] + b[2];

   return;

}

void eraPpsp(double a[3], double s, double b[3], double apsb[3])
/*
**  - - - - - - - -
**   e r a P p s p
**  - - - - - - - -
**
**  P-vector plus scaled p-vector.
**
**  Given:
**     a      double[3]     first p-vector
**     s      double        scalar (multiplier for b)
**     b      double[3]     second p-vector
**
**  Returned:
**     apsb   double[3]     a + s*b
**
**  Note:
**     It is permissible for any of a, b and apsb to be the same array.
**
**  Called:
**     eraSxp       multiply p-vector by scalar
**     eraPpp       p-vector plus p-vector
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double sb[3];


/* s*b. */
   eraSxp(s, b, sb);

/* a + s*b. */
   eraPpp(a, sb, apsb);

   return;

}

void eraPr00(double date1, double date2, double *dpsipr, double *depspr)
/*
**  - - - - - - - -
**   e r a P r 0 0
**  - - - - - - - -
**
**  Precession-rate part of the IAU 2000 precession-nutation models
**  (part of MHB2000).
**
**  Given:
**     date1,date2    double  TT as a 2-part Julian Date (Note 1)
**
**  Returned:
**     dpsipr,depspr  double  precession corrections (Notes 2,3)
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
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
**  2) The precession adjustments are expressed as "nutation
**     components", corrections in longitude and obliquity with respect
**     to the J2000.0 equinox and ecliptic.
**
**  3) Although the precession adjustments are stated to be with respect
**     to Lieske et al. (1977), the MHB2000 model does not specify which
**     set of Euler angles are to be used and how the adjustments are to
**     be applied.  The most literal and straightforward procedure is to
**     adopt the 4-rotation epsilon_0, psi_A, omega_A, xi_A option, and
**     to add dpsipr to psi_A and depspr to both omega_A and eps_A.
**
**  4) This is an implementation of one aspect of the IAU 2000A nutation
**     model, formally adopted by the IAU General Assembly in 2000,
**     namely MHB2000 (Mathews et al. 2002).
**
**  References:
**
**     Lieske, J.H., Lederle, T., Fricke, W. & Morando, B., "Expressions
**     for the precession quantities based upon the IAU (1976) System of
**     Astronomical Constants", Astron.Astrophys., 58, 1-16 (1977)
**
**     Mathews, P.M., Herring, T.A., Buffet, B.A., "Modeling of nutation
**     and precession   New nutation series for nonrigid Earth and
**     insights into the Earth's interior", J.Geophys.Res., 107, B4,
**     2002.  The MHB2000 code itself was obtained on 9th September 2002
**     from ftp://maia.usno.navy.mil/conv2000/chapter5/IAU2000A.
**
**     Wallace, P.T., "Software for Implementing the IAU 2000
**     Resolutions", in IERS Workshop 5.1 (2002).
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double t;

/* Precession and obliquity corrections (radians per century) */
   static const double PRECOR = -0.29965 * DAS2R,
                       OBLCOR = -0.02524 * DAS2R;


/* Interval between fundamental epoch J2000.0 and given date (JC). */
   t = ((date1 - DJ00) + date2) / DJC;

/* Precession rate contributions with respect to IAU 1976/80. */
   *dpsipr = PRECOR * t;
   *depspr = OBLCOR * t;

   return;

}

void eraPrec76(double ep01, double ep02, double ep11, double ep12,
               double *zeta, double *z, double *theta)
/*
**  - - - - - - - - - -
**   e r a P r e c 7 6
**  - - - - - - - - - -
**
**  IAU 1976 precession model.
**
**  This function forms the three Euler angles which implement general
**  precession between two epochs, using the IAU 1976 model (as for
**  the FK5 catalog).
**
**  Given:
**     ep01,ep02   double    TDB starting epoch (Note 1)
**     ep11,ep12   double    TDB ending epoch (Note 1)
**
**  Returned:
**     zeta        double    1st rotation: radians cw around z
**     z           double    3rd rotation: radians cw around z
**     theta       double    2nd rotation: radians ccw around y
**
**  Notes:
**
**  1) The epochs ep01+ep02 and ep11+ep12 are Julian Dates, apportioned
**     in any convenient way between the arguments epn1 and epn2.  For
**     example, JD(TDB)=2450123.7 could be expressed in any of these
**     ways, among others:
**
**             epn1          epn2
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in cases
**     where the loss of several decimal digits of resolution is
**     acceptable.  The J2000 method is best matched to the way the
**     argument is handled internally and will deliver the optimum
**     optimum resolution.  The MJD method and the date & time methods
**     are both good compromises between resolution and convenience.
**     The two epochs may be expressed using different methods, but at
**     the risk of losing some resolution.
**
**  2) The accumulated precession angles zeta, z, theta are expressed
**     through canonical polynomials which are valid only for a limited
**     time span.  In addition, the IAU 1976 precession rate is known to
**     be imperfect.  The absolute accuracy of the present formulation
**     is better than 0.1 arcsec from 1960AD to 2040AD, better than
**     1 arcsec from 1640AD to 2360AD, and remains below 3 arcsec for
**     the whole of the period 500BC to 3000AD.  The errors exceed
**     10 arcsec outside the range 1200BC to 3900AD, exceed 100 arcsec
**     outside 4200BC to 5600AD and exceed 1000 arcsec outside 6800BC to
**     8200AD.
**
**  3) The three angles are returned in the conventional order, which
**     is not the same as the order of the corresponding Euler
**     rotations.  The precession matrix is
**     R_3(-z) x R_2(+theta) x R_3(-zeta).
**
**  Reference:
**
**     Lieske, J.H., 1979, Astron.Astrophys. 73, 282, equations
**     (6) & (7), p283.
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double t0, t, tas2r, w;


/* Interval between fundamental epoch J2000.0 and start epoch (JC). */
   t0 = ((ep01 - DJ00) + ep02) / DJC;

/* Interval over which precession required (JC). */
   t = ((ep11 - ep01) + (ep12 - ep02)) / DJC;

/* Euler angles. */
   tas2r = t * DAS2R;
   w = 2306.2181 + (1.39656 - 0.000139 * t0) * t0;

   *zeta = (w + ((0.30188 - 0.000344 * t0) + 0.017998 * t) * t) * tas2r;

   *z = (w + ((1.09468 + 0.000066 * t0) + 0.018203 * t) * t) * tas2r;

   *theta = ((2004.3109 + (-0.85330 - 0.000217 * t0) * t0)
          + ((-0.42665 - 0.000217 * t0) - 0.041833 * t) * t) * tas2r;

   return;

}

void eraPv2p(double pv[2][3], double p[3])
/*
**  - - - - - - - -
**   e r a P v 2 p
**  - - - - - - - -
**
**  Discard velocity component of a pv-vector.
**
**  Given:
**     pv      double[2][3]     pv-vector
**
**  Returned:
**     p       double[3]        p-vector
**
**  Called:
**     eraCp        copy p-vector
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   eraCp(pv[0], p);

   return;

}

void eraPv2s(double pv[2][3],
             double *theta, double *phi, double *r,
             double *td, double *pd, double *rd)
/*
**  - - - - - - - -
**   e r a P v 2 s
**  - - - - - - - -
**
**  Convert position/velocity from Cartesian to spherical coordinates.
**
**  Given:
**     pv       double[2][3]  pv-vector
**
**  Returned:
**     theta    double        longitude angle (radians)
**     phi      double        latitude angle (radians)
**     r        double        radial distance
**     td       double        rate of change of theta
**     pd       double        rate of change of phi
**     rd       double        rate of change of r
**
**  Notes:
**
**  1) If the position part of pv is null, theta, phi, td and pd
**     are indeterminate.  This is handled by extrapolating the
**     position through unit time by using the velocity part of
**     pv.  This moves the origin without changing the direction
**     of the velocity component.  If the position and velocity
**     components of pv are both null, zeroes are returned for all
**     six results.
**
**  2) If the position is a pole, theta, td and pd are indeterminate.
**     In such cases zeroes are returned for all three.
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double x, y, z, xd, yd, zd, rxy2, rxy, r2, rtrue, rw, xyp;


/* Components of position/velocity vector. */
   x  = pv[0][0];
   y  = pv[0][1];
   z  = pv[0][2];
   xd = pv[1][0];
   yd = pv[1][1];
   zd = pv[1][2];

/* Component of r in XY plane squared. */
   rxy2 = x*x + y*y;

/* Modulus squared. */
   r2 = rxy2 + z*z;

/* Modulus. */
   rtrue = sqrt(r2);

/* If null vector, move the origin along the direction of movement. */
   rw = rtrue;
   if (rtrue == 0.0) {
       x = xd;
       y = yd;
       z = zd;
       rxy2 = x*x + y*y;
       r2 = rxy2 + z*z;
       rw = sqrt(r2);
   }

/* Position and velocity in spherical coordinates. */
   rxy = sqrt(rxy2);
   xyp = x*xd + y*yd;
   if (rxy2 != 0.0) {
       *theta = atan2(y, x);
       *phi = atan2(z, rxy);
       *td = (x*yd - y*xd) / rxy2;
       *pd = (zd*rxy2 - z*xyp) / (r2*rxy);
   } else {
       *theta = 0.0;
       *phi = (z != 0.0) ? atan2(z, rxy) : 0.0;
       *td = 0.0;
       *pd = 0.0;
   }
   *r = rtrue;
   *rd = (rw != 0.0) ? (xyp + z*zd) / rw : 0.0;

   return;

}

void eraPvdpv(double a[2][3], double b[2][3], double adb[2])
/*
**  - - - - - - - - -
**   e r a P v d p v
**  - - - - - - - - -
**
**  Inner (=scalar=dot) product of two pv-vectors.
**
**  Given:
**     a        double[2][3]      first pv-vector
**     b        double[2][3]      second pv-vector
**
**  Returned:
**     adb      double[2]         a . b (see note)
**
**  Note:
**
**     If the position and velocity components of the two pv-vectors are
**     ( ap, av ) and ( bp, bv ), the result, a . b, is the pair of
**     numbers ( ap . bp , ap . bv + av . bp ).  The two numbers are the
**     dot-product of the two p-vectors and its derivative.
**
**  Called:
**     eraPdp       scalar product of two p-vectors
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double adbd, addb;


/* a . b = constant part of result. */
   adb[0] = eraPdp(a[0], b[0]);

/* a . bdot */
   adbd = eraPdp(a[0], b[1]);

/* adot . b */
   addb = eraPdp(a[1], b[0]);

/* Velocity part of result. */
   adb[1] = adbd + addb;

   return;

}

void eraPvm(double pv[2][3], double *r, double *s)
/*
**  - - - - - - -
**   e r a P v m
**  - - - - - - -
**
**  Modulus of pv-vector.
**
**  Given:
**     pv     double[2][3]   pv-vector
**
**  Returned:
**     r      double         modulus of position component
**     s      double         modulus of velocity component
**
**  Called:
**     eraPm        modulus of p-vector
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
/* Distance. */
   *r = eraPm(pv[0]);

/* Speed. */
   *s = eraPm(pv[1]);

   return;

}

void eraPvmpv(double a[2][3], double b[2][3], double amb[2][3])
/*
**  - - - - - - - - -
**   e r a P v m p v
**  - - - - - - - - -
**
**  Subtract one pv-vector from another.
**
**  Given:
**     a       double[2][3]      first pv-vector
**     b       double[2][3]      second pv-vector
**
**  Returned:
**     amb     double[2][3]      a - b
**
**  Note:
**     It is permissible to re-use the same array for any of the
**     arguments.
**
**  Called:
**     eraPmp       p-vector minus p-vector
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   eraPmp(a[0], b[0], amb[0]);
   eraPmp(a[1], b[1], amb[1]);

   return;

}

void eraPvppv(double a[2][3], double b[2][3], double apb[2][3])
/*
**  - - - - - - - - -
**   e r a P v p p v
**  - - - - - - - - -
**
**  Add one pv-vector to another.
**
**  Given:
**     a        double[2][3]      first pv-vector
**     b        double[2][3]      second pv-vector
**
**  Returned:
**     apb      double[2][3]      a + b
**
**  Note:
**     It is permissible to re-use the same array for any of the
**     arguments.
**
**  Called:
**     eraPpp       p-vector plus p-vector
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   eraPpp(a[0], b[0], apb[0]);
   eraPpp(a[1], b[1], apb[1]);

   return;

}

int eraPvstar(double pv[2][3], double *ra, double *dec,
              double *pmr, double *pmd, double *px, double *rv)
/*
**  - - - - - - - - - -
**   e r a P v s t a r
**  - - - - - - - - - -
**
**  Convert star position+velocity vector to catalog coordinates.
**
**  Given (Note 1):
**     pv     double[2][3]   pv-vector (AU, AU/day)
**
**  Returned (Note 2):
**     ra     double         right ascension (radians)
**     dec    double         declination (radians)
**     pmr    double         RA proper motion (radians/year)
**     pmd    double         Dec proper motion (radians/year)
**     px     double         parallax (arcsec)
**     rv     double         radial velocity (km/s, positive = receding)
**
**  Returned (function value):
**            int            status:
**                              0 = OK
**                             -1 = superluminal speed (Note 5)
**                             -2 = null position vector
**
**  Notes:
**
**  1) The specified pv-vector is the coordinate direction (and its rate
**     of change) for the date at which the light leaving the star
**     reached the solar-system barycenter.
**
**  2) The star data returned by this function are "observables" for an
**     imaginary observer at the solar-system barycenter.  Proper motion
**     and radial velocity are, strictly, in terms of barycentric
**     coordinate time, TCB.  For most practical applications, it is
**     permissible to neglect the distinction between TCB and ordinary
**     "proper" time on Earth (TT/TAI).  The result will, as a rule, be
**     limited by the intrinsic accuracy of the proper-motion and
**     radial-velocity data;  moreover, the supplied pv-vector is likely
**     to be merely an intermediate result (for example generated by the
**     function eraStarpv), so that a change of time unit will cancel
**     out overall.
**
**     In accordance with normal star-catalog conventions, the object's
**     right ascension and declination are freed from the effects of
**     secular aberration.  The frame, which is aligned to the catalog
**     equator and equinox, is Lorentzian and centered on the SSB.
**
**     Summarizing, the specified pv-vector is for most stars almost
**     identical to the result of applying the standard geometrical
**     "space motion" transformation to the catalog data.  The
**     differences, which are the subject of the Stumpff paper cited
**     below, are:
**
**     (i) In stars with significant radial velocity and proper motion,
**     the constantly changing light-time distorts the apparent proper
**     motion.  Note that this is a classical, not a relativistic,
**     effect.
**
**     (ii) The transformation complies with special relativity.
**
**  3) Care is needed with units.  The star coordinates are in radians
**     and the proper motions in radians per Julian year, but the
**     parallax is in arcseconds; the radial velocity is in km/s, but
**     the pv-vector result is in AU and AU/day.
**
**  4) The proper motions are the rate of change of the right ascension
**     and declination at the catalog epoch and are in radians per Julian
**     year.  The RA proper motion is in terms of coordinate angle, not
**     true angle, and will thus be numerically larger at high
**     declinations.
**
**  5) Straight-line motion at constant speed in the inertial frame is
**     assumed.  If the speed is greater than or equal to the speed of
**     light, the function aborts with an error status.
**
**  6) The inverse transformation is performed by the function eraStarpv.
**
**  Called:
**     eraPn        decompose p-vector into modulus and direction
**     eraPdp       scalar product of two p-vectors
**     eraSxp       multiply p-vector by scalar
**     eraPmp       p-vector minus p-vector
**     eraPm        modulus of p-vector
**     eraPpp       p-vector plus p-vector
**     eraPv2s      pv-vector to spherical
**     eraAnp       normalize angle into range 0 to 2pi
**
**  Reference:
**
**     Stumpff, P., 1985, Astron.Astrophys. 144, 232-240.
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double r, x[3], vr, ur[3], vt, ut[3], bett, betr, d, w, del,
          usr[3], ust[3], a, rad, decd, rd;


/* Isolate the radial component of the velocity (AU/day, inertial). */
   eraPn(pv[0], &r, x);
   vr = eraPdp(x, pv[1]);
   eraSxp(vr, x, ur);

/* Isolate the transverse component of the velocity (AU/day, inertial). */
   eraPmp(pv[1], ur, ut);
   vt = eraPm(ut);

/* Special-relativity dimensionless parameters. */
   bett = vt / DC;
   betr = vr / DC;

/* The inertial-to-observed correction terms. */
   d = 1.0 + betr;
   w = 1.0 - betr*betr - bett*bett;
   if (d == 0.0 || w < 0) return -1;
   del = sqrt(w) - 1.0;

/* Apply relativistic correction factor to radial velocity component. */
   w = (betr != 0) ? (betr - del) / (betr * d) : 1.0;
   eraSxp(w, ur, usr);

/* Apply relativistic correction factor to tangential velocity */
/* component.                                                  */
   eraSxp(1.0/d, ut, ust);

/* Combine the two to obtain the observed velocity vector (AU/day). */
   eraPpp(usr, ust, pv[1]);

/* Cartesian to spherical. */
   eraPv2s(pv, &a, dec, &r, &rad, &decd, &rd);
   if (r == 0.0) return -2;

/* Return RA in range 0 to 2pi. */
   *ra = eraAnp(a);

/* Return proper motions in radians per year. */
   *pmr = rad * DJY;
   *pmd = decd * DJY;

/* Return parallax in arcsec. */
   *px = DR2AS / r;

/* Return radial velocity in km/s. */
   *rv = 1e-3 * rd * DAU / DAYSEC;

/* OK status. */
   return 0;

}

void eraPvu(double dt, double pv[2][3], double upv[2][3])
/*
**  - - - - - - -
**   e r a P v u
**  - - - - - - -
**
**  Update a pv-vector.
**
**  Given:
**     dt       double           time interval
**     pv       double[2][3]     pv-vector
**
**  Returned:
**     upv      double[2][3]     p updated, v unchanged
**
**  Notes:
**
**  1) "Update" means "refer the position component of the vector
**     to a new date dt time units from the existing date".
**
**  2) The time units of dt must match those of the velocity.
**
**  3) It is permissible for pv and upv to be the same array.
**
**  Called:
**     eraPpsp      p-vector plus scaled p-vector
**     eraCp        copy p-vector
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   eraPpsp(pv[0], dt, pv[1], upv[0]);
   eraCp(pv[1], upv[1]);

   return;

}

void eraPvup(double dt, double pv[2][3], double p[3])
/*
**  - - - - - - - -
**   e r a P v u p
**  - - - - - - - -
**
**  Update a pv-vector, discarding the velocity component.
**
**  Given:
**     dt       double            time interval
**     pv       double[2][3]      pv-vector
**
**  Returned:
**     p        double[3]         p-vector
**
**  Notes:
**
**  1) "Update" means "refer the position component of the vector to a
**     new date dt time units from the existing date".
**
**  2) The time units of dt must match those of the velocity.
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   p[0] = pv[0][0] + dt * pv[1][0];
   p[1] = pv[0][1] + dt * pv[1][1];
   p[2] = pv[0][2] + dt * pv[1][2];

   return;

}

void eraPvxpv(double a[2][3], double b[2][3], double axb[2][3])
/*
**  - - - - - - - - -
**   e r a P v x p v
**  - - - - - - - - -
**
**  Outer (=vector=cross) product of two pv-vectors.
**
**  Given:
**     a        double[2][3]      first pv-vector
**     b        double[2][3]      second pv-vector
**
**  Returned:
**     axb      double[2][3]      a x b
**
**  Notes:
**
**  1) If the position and velocity components of the two pv-vectors are
**     ( ap, av ) and ( bp, bv ), the result, a x b, is the pair of
**     vectors ( ap x bp, ap x bv + av x bp ).  The two vectors are the
**     cross-product of the two p-vectors and its derivative.
**
**  2) It is permissible to re-use the same array for any of the
**     arguments.
**
**  Called:
**     eraCpv       copy pv-vector
**     eraPxp       vector product of two p-vectors
**     eraPpp       p-vector plus p-vector
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double wa[2][3], wb[2][3], axbd[3], adxb[3];


/* Make copies of the inputs. */
   eraCpv(a, wa);
   eraCpv(b, wb);

/* a x b = position part of result. */
   eraPxp(wa[0], wb[0], axb[0]);

/* a x bdot + adot x b = velocity part of result. */
   eraPxp(wa[0], wb[1], axbd);
   eraPxp(wa[1], wb[0], adxb);
   eraPpp(axbd, adxb, axb[1]);

   return;

}

void eraPxp(double a[3], double b[3], double axb[3])
/*
**  - - - - - - -
**   e r a P x p
**  - - - - - - -
**
**  p-vector outer (=vector=cross) product.
**
**  Given:
**     a        double[3]      first p-vector
**     b        double[3]      second p-vector
**
**  Returned:
**     axb      double[3]      a x b
**
**  Note:
**     It is permissible to re-use the same array for any of the
**     arguments.
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double xa, ya, za, xb, yb, zb;


   xa = a[0];
   ya = a[1];
   za = a[2];
   xb = b[0];
   yb = b[1];
   zb = b[2];
   axb[0] = ya*zb - za*yb;
   axb[1] = za*xb - xa*zb;
   axb[2] = xa*yb - ya*xb;

   return;

}

void eraRm2v(double r[3][3], double w[3])
/*
**  - - - - - - - -
**   e r a R m 2 v
**  - - - - - - - -
**
**  Express an r-matrix as an r-vector.
**
**  Given:
**     r        double[3][3]    rotation matrix
**
**  Returned:
**     w        double[3]       rotation vector (Note 1)
**
**  Notes:
**
**  1) A rotation matrix describes a rotation through some angle about
**     some arbitrary axis called the Euler axis.  The "rotation vector"
**     returned by this function has the same direction as the Euler axis,
**     and its magnitude is the angle in radians.  (The magnitude and
**     direction can be separated by means of the function eraPn.)
**
**  2) If r is null, so is the result.  If r is not a rotation matrix
**     the result is undefined;  r must be proper (i.e. have a positive
**     determinant) and real orthogonal (inverse = transpose).
**
**  3) The reference frame rotates clockwise as seen looking along
**     the rotation vector from the origin.
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double x, y, z, s2, c2, phi, f;


   x = r[1][2] - r[2][1];
   y = r[2][0] - r[0][2];
   z = r[0][1] - r[1][0];
   s2 = sqrt(x*x + y*y + z*z);
   if (s2 != 0) {
      c2 = r[0][0] + r[1][1] + r[2][2] - 1.0;
      phi = atan2(s2, c2);
      f =  phi / s2;
      w[0] = x * f;
      w[1] = y * f;
      w[2] = z * f;
   } else {
      w[0] = 0.0;
      w[1] = 0.0;
      w[2] = 0.0;
   }

   return;

}

void eraRv2m(double w[3], double r[3][3])
/*
**  - - - - - - - -
**   e r a R v 2 m
**  - - - - - - - -
**
**  Form the r-matrix corresponding to a given r-vector.
**
**  Given:
**     w        double[3]      rotation vector (Note 1)
**
**  Returned:
**     r        double[3][3]    rotation matrix
**
**  Notes:
**
**  1) A rotation matrix describes a rotation through some angle about
**     some arbitrary axis called the Euler axis.  The "rotation vector"
**     supplied to This function has the same direction as the Euler
**     axis, and its magnitude is the angle in radians.
**
**  2) If w is null, the unit matrix is returned.
**
**  3) The reference frame rotates clockwise as seen looking along the
**     rotation vector from the origin.
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double x, y, z, phi, s, c, f;


/* Euler angle (magnitude of rotation vector) and functions. */
   x = w[0];
   y = w[1];
   z = w[2];
   phi = sqrt(x*x + y*y + z*z);
   s = sin(phi);
   c = cos(phi);
   f = 1.0 - c;

/* Euler axis (direction of rotation vector), perhaps null. */
   if (phi != 0.0) {
       x /= phi;
       y /= phi;
       z /= phi;
   }

/* Form the rotation matrix. */
   r[0][0] = x*x*f + c;
   r[0][1] = x*y*f + z*s;
   r[0][2] = x*z*f - y*s;
   r[1][0] = y*x*f - z*s;
   r[1][1] = y*y*f + c;
   r[1][2] = y*z*f + x*s;
   r[2][0] = z*x*f + y*s;
   r[2][1] = z*y*f - x*s;
   r[2][2] = z*z*f + c;

   return;

}

void eraRx(double phi, double r[3][3])
/*
**  - - - - - -
**   e r a R x
**  - - - - - -
**
**  Rotate an r-matrix about the x-axis.
**
**  Given:
**     phi    double          angle (radians)
**
**  Given and returned:
**     r      double[3][3]    r-matrix, rotated
**
**  Notes:
**
**  1) Calling this function with positive phi incorporates in the
**     supplied r-matrix r an additional rotation, about the x-axis,
**     anticlockwise as seen looking towards the origin from positive x.
**
**  2) The additional rotation can be represented by this matrix:
**
**         (  1        0            0      )
**         (                               )
**         (  0   + cos(phi)   + sin(phi)  )
**         (                               )
**         (  0   - sin(phi)   + cos(phi)  )
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double s, c, a10, a11, a12, a20, a21, a22;


   s = sin(phi);
   c = cos(phi);

   a10 =   c*r[1][0] + s*r[2][0];
   a11 =   c*r[1][1] + s*r[2][1];
   a12 =   c*r[1][2] + s*r[2][2];
   a20 = - s*r[1][0] + c*r[2][0];
   a21 = - s*r[1][1] + c*r[2][1];
   a22 = - s*r[1][2] + c*r[2][2];

   r[1][0] = a10;
   r[1][1] = a11;
   r[1][2] = a12;
   r[2][0] = a20;
   r[2][1] = a21;
   r[2][2] = a22;

   return;

}

void eraRxp(double r[3][3], double p[3], double rp[3])
/*
**  - - - - - - -
**   e r a R x p
**  - - - - - - -
**
**  Multiply a p-vector by an r-matrix.
**
**  Given:
**     r        double[3][3]    r-matrix
**     p        double[3]       p-vector
**
**  Returned:
**     rp       double[3]       r * p
**
**  Note:
**     It is permissible for p and rp to be the same array.
**
**  Called:
**     eraCp        copy p-vector
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double w, wrp[3];
   int i, j;


/* Matrix r * vector p. */
   for (j = 0; j < 3; j++) {
       w = 0.0;
       for (i = 0; i < 3; i++) {
           w += r[j][i] * p[i];
       }
       wrp[j] = w;
   }

/* Return the result. */
   eraCp(wrp, rp);

   return;

}

void eraRxpv(double r[3][3], double pv[2][3], double rpv[2][3])
/*
**  - - - - - - - -
**   e r a R x p v
**  - - - - - - - -
**
**  Multiply a pv-vector by an r-matrix.
**
**  Given:
**     r        double[3][3]    r-matrix
**     pv       double[2][3]    pv-vector
**
**  Returned:
**     rpv      double[2][3]    r * pv
**
**  Note:
**     It is permissible for pv and rpv to be the same array.
**
**  Called:
**     eraRxp       product of r-matrix and p-vector
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   eraRxp(r, pv[0], rpv[0]);
   eraRxp(r, pv[1], rpv[1]);

   return;

}

void eraRxr(double a[3][3], double b[3][3], double atb[3][3])
/*
**  - - - - - - -
**   e r a R x r
**  - - - - - - -
**
**  Multiply two r-matrices.
**
**  Given:
**     a        double[3][3]    first r-matrix
**     b        double[3][3]    second r-matrix
**
**  Returned:
**     atb      double[3][3]    a * b
**
**  Note:
**     It is permissible to re-use the same array for any of the
**     arguments.
**
**  Called:
**     eraCr        copy r-matrix
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   int i, j, k;
   double w, wm[3][3];


   for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++) {
         w = 0.0;
         for (k = 0; k < 3; k++) {
            w +=  a[i][k] * b[k][j];
         }
         wm[i][j] = w;
      }
   }
   eraCr(wm, atb);

   return;

}

void eraRy(double theta, double r[3][3])
/*
**  - - - - - -
**   e r a R y
**  - - - - - -
**
**  Rotate an r-matrix about the y-axis.
**
**  Given:
**     theta  double          angle (radians)
**
**  Given and returned:
**     r      double[3][3]    r-matrix, rotated
**
**  Notes:
**
**  1) Calling this function with positive theta incorporates in the
**     supplied r-matrix r an additional rotation, about the y-axis,
**     anticlockwise as seen looking towards the origin from positive y.
**
**  2) The additional rotation can be represented by this matrix:
**
**         (  + cos(theta)     0      - sin(theta)  )
**         (                                        )
**         (       0           1           0        )
**         (                                        )
**         (  + sin(theta)     0      + cos(theta)  )
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double s, c, a00, a01, a02, a20, a21, a22;


   s = sin(theta);
   c = cos(theta);

   a00 = c*r[0][0] - s*r[2][0];
   a01 = c*r[0][1] - s*r[2][1];
   a02 = c*r[0][2] - s*r[2][2];
   a20 = s*r[0][0] + c*r[2][0];
   a21 = s*r[0][1] + c*r[2][1];
   a22 = s*r[0][2] + c*r[2][2];

   r[0][0] = a00;
   r[0][1] = a01;
   r[0][2] = a02;
   r[2][0] = a20;
   r[2][1] = a21;
   r[2][2] = a22;

   return;

}

void eraRz(double psi, double r[3][3])
/*
**  - - - - - -
**   e r a R z
**  - - - - - -
**
**  Rotate an r-matrix about the z-axis.
**
**  Given:
**     psi    double          angle (radians)
**
**  Given and returned:
**     r      double[3][3]    r-matrix, rotated
**
**  Notes:
**
**  1) Calling this function with positive psi incorporates in the
**     supplied r-matrix r an additional rotation, about the z-axis,
**     anticlockwise as seen looking towards the origin from positive z.
**
**  2) The additional rotation can be represented by this matrix:
**
**         (  + cos(psi)   + sin(psi)     0  )
**         (                                 )
**         (  - sin(psi)   + cos(psi)     0  )
**         (                                 )
**         (       0            0         1  )
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double s, c, a00, a01, a02, a10, a11, a12;


   s = sin(psi);
   c = cos(psi);

   a00 =   c*r[0][0] + s*r[1][0];
   a01 =   c*r[0][1] + s*r[1][1];
   a02 =   c*r[0][2] + s*r[1][2];
   a10 = - s*r[0][0] + c*r[1][0];
   a11 = - s*r[0][1] + c*r[1][1];
   a12 = - s*r[0][2] + c*r[1][2];

   r[0][0] = a00;
   r[0][1] = a01;
   r[0][2] = a02;
   r[1][0] = a10;
   r[1][1] = a11;
   r[1][2] = a12;

   return;

}

double eraS00(double date1, double date2, double x, double y)
/*
**  - - - - - - -
**   e r a S 0 0
**  - - - - - - -
**
**  The CIO locator s, positioning the Celestial Intermediate Origin on
**  the equator of the Celestial Intermediate Pole, given the CIP's X,Y
**  coordinates.  Compatible with IAU 2000A precession-nutation.
**
**  Given:
**     date1,date2   double    TT as a 2-part Julian Date (Note 1)
**     x,y           double    CIP coordinates (Note 3)
**
**  Returned (function value):
**                   double    the CIO locator s in radians (Note 2)
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
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
**  2) The CIO locator s is the difference between the right ascensions
**     of the same point in two systems:  the two systems are the GCRS
**     and the CIP,CIO, and the point is the ascending node of the
**     CIP equator.  The quantity s remains below 0.1 arcsecond
**     throughout 1900-2100.
**
**  3) The series used to compute s is in fact for s+XY/2, where X and Y
**     are the x and y components of the CIP unit vector;  this series
**     is more compact than a direct series for s would be.  This
**     function requires X,Y to be supplied by the caller, who is
**     responsible for providing values that are consistent with the
**     supplied date.
**
**  4) The model is consistent with the IAU 2000A precession-nutation.
**
**  Called:
**     eraFal03     mean anomaly of the Moon
**     eraFalp03    mean anomaly of the Sun
**     eraFaf03     mean argument of the latitude of the Moon
**     eraFad03     mean elongation of the Moon from the Sun
**     eraFaom03    mean longitude of the Moon's ascending node
**     eraFave03    mean longitude of Venus
**     eraFae03     mean longitude of Earth
**     eraFapa03    general accumulated precession in longitude
**
**  References:
**
**     Capitaine, N., Chapront, J., Lambert, S. and Wallace, P.,
**     "Expressions for the Celestial Intermediate Pole and Celestial
**     Ephemeris Origin consistent with the IAU 2000A precession-
**     nutation model", Astron.Astrophys. 400, 1145-1154 (2003)
**
**     n.b. The celestial ephemeris origin (CEO) was renamed "celestial
**          intermediate origin" (CIO) by IAU 2006 Resolution 2.
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
/* Time since J2000.0, in Julian centuries */
   double t;

/* Miscellaneous */
   int i, j;
   double a, w0, w1, w2, w3, w4, w5;

/* Fundamental arguments */
   double fa[8];

/* Returned value */
   double s;

/* --------------------- */
/* The series for s+XY/2 */
/* --------------------- */

   typedef struct {
      int nfa[8];      /* coefficients of l,l',F,D,Om,LVe,LE,pA */
      double s, c;     /* sine and cosine coefficients */
   } TERM;

/* Polynomial coefficients */
   static const double sp[] = {

   /* 1-6 */
          94.00e-6,
        3808.35e-6,
        -119.94e-6,
      -72574.09e-6,
          27.70e-6,
          15.61e-6
   };

/* Terms of order t^0 */
   static const TERM s0[] = {

   /* 1-10 */
      {{ 0,  0,  0,  0,  1,  0,  0,  0}, -2640.73e-6,   0.39e-6 },
      {{ 0,  0,  0,  0,  2,  0,  0,  0},   -63.53e-6,   0.02e-6 },
      {{ 0,  0,  2, -2,  3,  0,  0,  0},   -11.75e-6,  -0.01e-6 },
      {{ 0,  0,  2, -2,  1,  0,  0,  0},   -11.21e-6,  -0.01e-6 },
      {{ 0,  0,  2, -2,  2,  0,  0,  0},     4.57e-6,   0.00e-6 },
      {{ 0,  0,  2,  0,  3,  0,  0,  0},    -2.02e-6,   0.00e-6 },
      {{ 0,  0,  2,  0,  1,  0,  0,  0},    -1.98e-6,   0.00e-6 },
      {{ 0,  0,  0,  0,  3,  0,  0,  0},     1.72e-6,   0.00e-6 },
      {{ 0,  1,  0,  0,  1,  0,  0,  0},     1.41e-6,   0.01e-6 },
      {{ 0,  1,  0,  0, -1,  0,  0,  0},     1.26e-6,   0.01e-6 },

   /* 11-20 */
      {{ 1,  0,  0,  0, -1,  0,  0,  0},     0.63e-6,   0.00e-6 },
      {{ 1,  0,  0,  0,  1,  0,  0,  0},     0.63e-6,   0.00e-6 },
      {{ 0,  1,  2, -2,  3,  0,  0,  0},    -0.46e-6,   0.00e-6 },
      {{ 0,  1,  2, -2,  1,  0,  0,  0},    -0.45e-6,   0.00e-6 },
      {{ 0,  0,  4, -4,  4,  0,  0,  0},    -0.36e-6,   0.00e-6 },
      {{ 0,  0,  1, -1,  1, -8, 12,  0},     0.24e-6,   0.12e-6 },
      {{ 0,  0,  2,  0,  0,  0,  0,  0},    -0.32e-6,   0.00e-6 },
      {{ 0,  0,  2,  0,  2,  0,  0,  0},    -0.28e-6,   0.00e-6 },
      {{ 1,  0,  2,  0,  3,  0,  0,  0},    -0.27e-6,   0.00e-6 },
      {{ 1,  0,  2,  0,  1,  0,  0,  0},    -0.26e-6,   0.00e-6 },

   /* 21-30 */
      {{ 0,  0,  2, -2,  0,  0,  0,  0},     0.21e-6,   0.00e-6 },
      {{ 0,  1, -2,  2, -3,  0,  0,  0},    -0.19e-6,   0.00e-6 },
      {{ 0,  1, -2,  2, -1,  0,  0,  0},    -0.18e-6,   0.00e-6 },
      {{ 0,  0,  0,  0,  0,  8,-13, -1},     0.10e-6,  -0.05e-6 },
      {{ 0,  0,  0,  2,  0,  0,  0,  0},    -0.15e-6,   0.00e-6 },
      {{ 2,  0, -2,  0, -1,  0,  0,  0},     0.14e-6,   0.00e-6 },
      {{ 0,  1,  2, -2,  2,  0,  0,  0},     0.14e-6,   0.00e-6 },
      {{ 1,  0,  0, -2,  1,  0,  0,  0},    -0.14e-6,   0.00e-6 },
      {{ 1,  0,  0, -2, -1,  0,  0,  0},    -0.14e-6,   0.00e-6 },
      {{ 0,  0,  4, -2,  4,  0,  0,  0},    -0.13e-6,   0.00e-6 },

   /* 31-33 */
      {{ 0,  0,  2, -2,  4,  0,  0,  0},     0.11e-6,   0.00e-6 },
      {{ 1,  0, -2,  0, -3,  0,  0,  0},    -0.11e-6,   0.00e-6 },
      {{ 1,  0, -2,  0, -1,  0,  0,  0},    -0.11e-6,   0.00e-6 }
   };

/* Terms of order t^1 */
   static const TERM s1[] ={

   /* 1-3 */
      {{ 0,  0,  0,  0,  2,  0,  0,  0},    -0.07e-6,   3.57e-6 },
      {{ 0,  0,  0,  0,  1,  0,  0,  0},     1.71e-6,  -0.03e-6 },
      {{ 0,  0,  2, -2,  3,  0,  0,  0},     0.00e-6,   0.48e-6 }
   };

/* Terms of order t^2 */
   static const TERM s2[] ={

   /* 1-10 */
      {{ 0,  0,  0,  0,  1,  0,  0,  0},   743.53e-6,  -0.17e-6 },
      {{ 0,  0,  2, -2,  2,  0,  0,  0},    56.91e-6,   0.06e-6 },
      {{ 0,  0,  2,  0,  2,  0,  0,  0},     9.84e-6,  -0.01e-6 },
      {{ 0,  0,  0,  0,  2,  0,  0,  0},    -8.85e-6,   0.01e-6 },
      {{ 0,  1,  0,  0,  0,  0,  0,  0},    -6.38e-6,  -0.05e-6 },
      {{ 1,  0,  0,  0,  0,  0,  0,  0},    -3.07e-6,   0.00e-6 },
      {{ 0,  1,  2, -2,  2,  0,  0,  0},     2.23e-6,   0.00e-6 },
      {{ 0,  0,  2,  0,  1,  0,  0,  0},     1.67e-6,   0.00e-6 },
      {{ 1,  0,  2,  0,  2,  0,  0,  0},     1.30e-6,   0.00e-6 },
      {{ 0,  1, -2,  2, -2,  0,  0,  0},     0.93e-6,   0.00e-6 },

   /* 11-20 */
      {{ 1,  0,  0, -2,  0,  0,  0,  0},     0.68e-6,   0.00e-6 },
      {{ 0,  0,  2, -2,  1,  0,  0,  0},    -0.55e-6,   0.00e-6 },
      {{ 1,  0, -2,  0, -2,  0,  0,  0},     0.53e-6,   0.00e-6 },
      {{ 0,  0,  0,  2,  0,  0,  0,  0},    -0.27e-6,   0.00e-6 },
      {{ 1,  0,  0,  0,  1,  0,  0,  0},    -0.27e-6,   0.00e-6 },
      {{ 1,  0, -2, -2, -2,  0,  0,  0},    -0.26e-6,   0.00e-6 },
      {{ 1,  0,  0,  0, -1,  0,  0,  0},    -0.25e-6,   0.00e-6 },
      {{ 1,  0,  2,  0,  1,  0,  0,  0},     0.22e-6,   0.00e-6 },
      {{ 2,  0,  0, -2,  0,  0,  0,  0},    -0.21e-6,   0.00e-6 },
      {{ 2,  0, -2,  0, -1,  0,  0,  0},     0.20e-6,   0.00e-6 },

   /* 21-25 */
      {{ 0,  0,  2,  2,  2,  0,  0,  0},     0.17e-6,   0.00e-6 },
      {{ 2,  0,  2,  0,  2,  0,  0,  0},     0.13e-6,   0.00e-6 },
      {{ 2,  0,  0,  0,  0,  0,  0,  0},    -0.13e-6,   0.00e-6 },
      {{ 1,  0,  2, -2,  2,  0,  0,  0},    -0.12e-6,   0.00e-6 },
      {{ 0,  0,  2,  0,  0,  0,  0,  0},    -0.11e-6,   0.00e-6 }
   };

/* Terms of order t^3 */
   static const TERM s3[] ={

   /* 1-4 */
      {{ 0,  0,  0,  0,  1,  0,  0,  0},     0.30e-6, -23.51e-6 },
      {{ 0,  0,  2, -2,  2,  0,  0,  0},    -0.03e-6,  -1.39e-6 },
      {{ 0,  0,  2,  0,  2,  0,  0,  0},    -0.01e-6,  -0.24e-6 },
      {{ 0,  0,  0,  0,  2,  0,  0,  0},     0.00e-6,   0.22e-6 }
   };

/* Terms of order t^4 */
   static const TERM s4[] ={

   /* 1-1 */
      {{ 0,  0,  0,  0,  1,  0,  0,  0},    -0.26e-6,  -0.01e-6 }
   };

/* Number of terms in the series */
   const int NS0 = (int) (sizeof s0 / sizeof (TERM));
   const int NS1 = (int) (sizeof s1 / sizeof (TERM));
   const int NS2 = (int) (sizeof s2 / sizeof (TERM));
   const int NS3 = (int) (sizeof s3 / sizeof (TERM));
   const int NS4 = (int) (sizeof s4 / sizeof (TERM));

/*--------------------------------------------------------------------*/

/* Interval between fundamental epoch J2000.0 and current date (JC). */
   t = ((date1 - DJ00) + date2) / DJC;

/* Fundamental Arguments (from IERS Conventions 2003) */

/* Mean anomaly of the Moon. */
   fa[0] = eraFal03(t);

/* Mean anomaly of the Sun. */
   fa[1] = eraFalp03(t);

/* Mean longitude of the Moon minus that of the ascending node. */
   fa[2] = eraFaf03(t);

/* Mean elongation of the Moon from the Sun. */
   fa[3] = eraFad03(t);

/* Mean longitude of the ascending node of the Moon. */
   fa[4] = eraFaom03(t);

/* Mean longitude of Venus. */
   fa[5] = eraFave03(t);

/* Mean longitude of Earth. */
   fa[6] = eraFae03(t);

/* General precession in longitude. */
   fa[7] = eraFapa03(t);

/* Evaluate s. */
   w0 = sp[0];
   w1 = sp[1];
   w2 = sp[2];
   w3 = sp[3];
   w4 = sp[4];
   w5 = sp[5];

   for (i = NS0-1; i >= 0; i--) {
   a = 0.0;
   for (j = 0; j < 8; j++) {
       a += (double)s0[i].nfa[j] * fa[j];
   }
   w0 += s0[i].s * sin(a) + s0[i].c * cos(a);
   }

   for (i = NS1-1; i >= 0; i--) {
   a = 0.0;
   for (j = 0; j < 8; j++) {
       a += (double)s1[i].nfa[j] * fa[j];
   }
   w1 += s1[i].s * sin(a) + s1[i].c * cos(a);
   }

   for (i = NS2-1; i >= 0; i--) {
   a = 0.0;
   for (j = 0; j < 8; j++) {
       a += (double)s2[i].nfa[j] * fa[j];
   }
   w2 += s2[i].s * sin(a) + s2[i].c * cos(a);
   }

   for (i = NS3-1; i >= 0; i--) {
   a = 0.0;
   for (j = 0; j < 8; j++) {
       a += (double)s3[i].nfa[j] * fa[j];
   }
   w3 += s3[i].s * sin(a) + s3[i].c * cos(a);
   }

   for (i = NS4-1; i >= 0; i--) {
   a = 0.0;
   for (j = 0; j < 8; j++) {
       a += (double)s4[i].nfa[j] * fa[j];
   }
   w4 += s4[i].s * sin(a) + s4[i].c * cos(a);
   }

   s = (w0 +
       (w1 +
       (w2 +
       (w3 +
       (w4 +
        w5 * t) * t) * t) * t) * t) * DAS2R - x*y/2.0;

   return s;

}

double eraS00a(double date1, double date2)
/*
**  - - - - - - - -
**   e r a S 0 0 a
**  - - - - - - - -
**
**  The CIO locator s, positioning the Celestial Intermediate Origin on
**  the equator of the Celestial Intermediate Pole, using the IAU 2000A
**  precession-nutation model.
**
**  Given:
**     date1,date2  double    TT as a 2-part Julian Date (Note 1)
**
**  Returned (function value):
**                  double    the CIO locator s in radians (Note 2)
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
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
**  2) The CIO locator s is the difference between the right ascensions
**     of the same point in two systems.  The two systems are the GCRS
**     and the CIP,CIO, and the point is the ascending node of the
**     CIP equator.  The CIO locator s remains a small fraction of
**     1 arcsecond throughout 1900-2100.
**
**  3) The series used to compute s is in fact for s+XY/2, where X and Y
**     are the x and y components of the CIP unit vector;  this series
**     is more compact than a direct series for s would be.  The present
**     function uses the full IAU 2000A nutation model when predicting
**     the CIP position.  Faster results, with no significant loss of
**     accuracy, can be obtained via the function eraS00b, which uses
**     instead the IAU 2000B truncated model.
**
**  Called:
**     eraPnm00a     classical NPB matrix, IAU 2000A
**     eraBnp2xy     extract CIP X,Y from the BPN matrix
**     eraS00        the CIO locator s, given X,Y, IAU 2000A
**
**  References:
**
**     Capitaine, N., Chapront, J., Lambert, S. and Wallace, P.,
**     "Expressions for the Celestial Intermediate Pole and Celestial
**     Ephemeris Origin consistent with the IAU 2000A precession-
**     nutation model", Astron.Astrophys. 400, 1145-1154 (2003)
**
**     n.b. The celestial ephemeris origin (CEO) was renamed "celestial
**          intermediate origin" (CIO) by IAU 2006 Resolution 2.
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double rbpn[3][3], x, y, s;


/* Bias-precession-nutation-matrix, IAU 2000A. */
   eraPnm00a(date1, date2, rbpn);

/* Extract the CIP coordinates. */
   eraBpn2xy(rbpn, &x, &y);

/* Compute the CIO locator s, given the CIP coordinates. */
   s = eraS00(date1, date2, x, y);

   return s;

}

double eraS00b(double date1, double date2)
/*
**  - - - - - - - -
**   e r a S 0 0 b
**  - - - - - - - -
**
**  The CIO locator s, positioning the Celestial Intermediate Origin on
**  the equator of the Celestial Intermediate Pole, using the IAU 2000B
**  precession-nutation model.
**
**  Given:
**     date1,date2  double    TT as a 2-part Julian Date (Note 1)
**
**  Returned (function value):
**                  double    the CIO locator s in radians (Note 2)
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
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
**  2) The CIO locator s is the difference between the right ascensions
**     of the same point in two systems.  The two systems are the GCRS
**     and the CIP,CIO, and the point is the ascending node of the
**     CIP equator.  The CIO locator s remains a small fraction of
**     1 arcsecond throughout 1900-2100.
**
**  3) The series used to compute s is in fact for s+XY/2, where X and Y
**     are the x and y components of the CIP unit vector;  this series
**     is more compact than a direct series for s would be.  The present
**     function uses the IAU 2000B truncated nutation model when
**     predicting the CIP position.  The function eraS00a uses instead
**     the full IAU 2000A model, but with no significant increase in
**     accuracy and at some cost in speed.
**
**  Called:
**     eraPnm00b     classical NPB matrix, IAU 2000B
**     eraBnp2xy     extract CIP X,Y from the BPN matrix
**     eraS00        the CIO locator s, given X,Y, IAU 2000A
**
**  References:
**
**     Capitaine, N., Chapront, J., Lambert, S. and Wallace, P.,
**     "Expressions for the Celestial Intermediate Pole and Celestial
**     Ephemeris Origin consistent with the IAU 2000A precession-
**     nutation model", Astron.Astrophys. 400, 1145-1154 (2003)
**
**     n.b. The celestial ephemeris origin (CEO) was renamed "celestial
**          intermediate origin" (CIO) by IAU 2006 Resolution 2.
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double rbpn[3][3], x, y, s;


/* Bias-precession-nutation-matrix, IAU 2000B. */
   eraPnm00b(date1, date2, rbpn);

/* Extract the CIP coordinates. */
   eraBpn2xy(rbpn, &x, &y);

/* Compute the CIO locator s, given the CIP coordinates. */
   s = eraS00(date1, date2, x, y);

   return s;

}

double eraS06(double date1, double date2, double x, double y)
/*
**  - - - - - - -
**   e r a S 0 6
**  - - - - - - -
**
**  The CIO locator s, positioning the Celestial Intermediate Origin on
**  the equator of the Celestial Intermediate Pole, given the CIP's X,Y
**  coordinates.  Compatible with IAU 2006/2000A precession-nutation.
**
**  Given:
**     date1,date2   double    TT as a 2-part Julian Date (Note 1)
**     x,y           double    CIP coordinates (Note 3)
**
**  Returned (function value):
**                   double    the CIO locator s in radians (Note 2)
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
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
**  2) The CIO locator s is the difference between the right ascensions
**     of the same point in two systems:  the two systems are the GCRS
**     and the CIP,CIO, and the point is the ascending node of the
**     CIP equator.  The quantity s remains below 0.1 arcsecond
**     throughout 1900-2100.
**
**  3) The series used to compute s is in fact for s+XY/2, where X and Y
**     are the x and y components of the CIP unit vector;  this series
**     is more compact than a direct series for s would be.  This
**     function requires X,Y to be supplied by the caller, who is
**     responsible for providing values that are consistent with the
**     supplied date.
**
**  4) The model is consistent with the "P03" precession (Capitaine et
**     al. 2003), adopted by IAU 2006 Resolution 1, 2006, and the
**     IAU 2000A nutation (with P03 adjustments).
**
**  Called:
**     eraFal03     mean anomaly of the Moon
**     eraFalp03    mean anomaly of the Sun
**     eraFaf03     mean argument of the latitude of the Moon
**     eraFad03     mean elongation of the Moon from the Sun
**     eraFaom03    mean longitude of the Moon's ascending node
**     eraFave03    mean longitude of Venus
**     eraFae03     mean longitude of Earth
**     eraFapa03    general accumulated precession in longitude
**
**  References:
**
**     Capitaine, N., Wallace, P.T. & Chapront, J., 2003, Astron.
**     Astrophys. 432, 355
**
**     McCarthy, D.D., Petit, G. (eds.) 2004, IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
/* Time since J2000.0, in Julian centuries */
   double t;

/* Miscellaneous */
   int i, j;
   double a, w0, w1, w2, w3, w4, w5;

/* Fundamental arguments */
   double fa[8];

/* Returned value */
   double s;

/* --------------------- */
/* The series for s+XY/2 */
/* --------------------- */

   typedef struct {
      int nfa[8];      /* coefficients of l,l',F,D,Om,LVe,LE,pA */
      double s, c;     /* sine and cosine coefficients */
   } TERM;

/* Polynomial coefficients */
   static const double sp[] = {

   /* 1-6 */
          94.00e-6,
        3808.65e-6,
        -122.68e-6,
      -72574.11e-6,
          27.98e-6,
          15.62e-6
   };

/* Terms of order t^0 */
   static const TERM s0[] = {

   /* 1-10 */
      {{ 0,  0,  0,  0,  1,  0,  0,  0}, -2640.73e-6,   0.39e-6 },
      {{ 0,  0,  0,  0,  2,  0,  0,  0},   -63.53e-6,   0.02e-6 },
      {{ 0,  0,  2, -2,  3,  0,  0,  0},   -11.75e-6,  -0.01e-6 },
      {{ 0,  0,  2, -2,  1,  0,  0,  0},   -11.21e-6,  -0.01e-6 },
      {{ 0,  0,  2, -2,  2,  0,  0,  0},     4.57e-6,   0.00e-6 },
      {{ 0,  0,  2,  0,  3,  0,  0,  0},    -2.02e-6,   0.00e-6 },
      {{ 0,  0,  2,  0,  1,  0,  0,  0},    -1.98e-6,   0.00e-6 },
      {{ 0,  0,  0,  0,  3,  0,  0,  0},     1.72e-6,   0.00e-6 },
      {{ 0,  1,  0,  0,  1,  0,  0,  0},     1.41e-6,   0.01e-6 },
      {{ 0,  1,  0,  0, -1,  0,  0,  0},     1.26e-6,   0.01e-6 },

   /* 11-20 */
      {{ 1,  0,  0,  0, -1,  0,  0,  0},     0.63e-6,   0.00e-6 },
      {{ 1,  0,  0,  0,  1,  0,  0,  0},     0.63e-6,   0.00e-6 },
      {{ 0,  1,  2, -2,  3,  0,  0,  0},    -0.46e-6,   0.00e-6 },
      {{ 0,  1,  2, -2,  1,  0,  0,  0},    -0.45e-6,   0.00e-6 },
      {{ 0,  0,  4, -4,  4,  0,  0,  0},    -0.36e-6,   0.00e-6 },
      {{ 0,  0,  1, -1,  1, -8, 12,  0},     0.24e-6,   0.12e-6 },
      {{ 0,  0,  2,  0,  0,  0,  0,  0},    -0.32e-6,   0.00e-6 },
      {{ 0,  0,  2,  0,  2,  0,  0,  0},    -0.28e-6,   0.00e-6 },
      {{ 1,  0,  2,  0,  3,  0,  0,  0},    -0.27e-6,   0.00e-6 },
      {{ 1,  0,  2,  0,  1,  0,  0,  0},    -0.26e-6,   0.00e-6 },

   /* 21-30 */
      {{ 0,  0,  2, -2,  0,  0,  0,  0},     0.21e-6,   0.00e-6 },
      {{ 0,  1, -2,  2, -3,  0,  0,  0},    -0.19e-6,   0.00e-6 },
      {{ 0,  1, -2,  2, -1,  0,  0,  0},    -0.18e-6,   0.00e-6 },
      {{ 0,  0,  0,  0,  0,  8,-13, -1},     0.10e-6,  -0.05e-6 },
      {{ 0,  0,  0,  2,  0,  0,  0,  0},    -0.15e-6,   0.00e-6 },
      {{ 2,  0, -2,  0, -1,  0,  0,  0},     0.14e-6,   0.00e-6 },
      {{ 0,  1,  2, -2,  2,  0,  0,  0},     0.14e-6,   0.00e-6 },
      {{ 1,  0,  0, -2,  1,  0,  0,  0},    -0.14e-6,   0.00e-6 },
      {{ 1,  0,  0, -2, -1,  0,  0,  0},    -0.14e-6,   0.00e-6 },
      {{ 0,  0,  4, -2,  4,  0,  0,  0},    -0.13e-6,   0.00e-6 },

   /* 31-33 */
      {{ 0,  0,  2, -2,  4,  0,  0,  0},     0.11e-6,   0.00e-6 },
      {{ 1,  0, -2,  0, -3,  0,  0,  0},    -0.11e-6,   0.00e-6 },
      {{ 1,  0, -2,  0, -1,  0,  0,  0},    -0.11e-6,   0.00e-6 }
   };

/* Terms of order t^1 */
   static const TERM s1[] = {

   /* 1 - 3 */
      {{ 0,  0,  0,  0,  2,  0,  0,  0},    -0.07e-6,   3.57e-6 },
      {{ 0,  0,  0,  0,  1,  0,  0,  0},     1.73e-6,  -0.03e-6 },
      {{ 0,  0,  2, -2,  3,  0,  0,  0},     0.00e-6,   0.48e-6 }
   };

/* Terms of order t^2 */
   static const TERM s2[] = {

   /* 1-10 */
      {{ 0,  0,  0,  0,  1,  0,  0,  0},   743.52e-6,  -0.17e-6 },
      {{ 0,  0,  2, -2,  2,  0,  0,  0},    56.91e-6,   0.06e-6 },
      {{ 0,  0,  2,  0,  2,  0,  0,  0},     9.84e-6,  -0.01e-6 },
      {{ 0,  0,  0,  0,  2,  0,  0,  0},    -8.85e-6,   0.01e-6 },
      {{ 0,  1,  0,  0,  0,  0,  0,  0},    -6.38e-6,  -0.05e-6 },
      {{ 1,  0,  0,  0,  0,  0,  0,  0},    -3.07e-6,   0.00e-6 },
      {{ 0,  1,  2, -2,  2,  0,  0,  0},     2.23e-6,   0.00e-6 },
      {{ 0,  0,  2,  0,  1,  0,  0,  0},     1.67e-6,   0.00e-6 },
      {{ 1,  0,  2,  0,  2,  0,  0,  0},     1.30e-6,   0.00e-6 },
      {{ 0,  1, -2,  2, -2,  0,  0,  0},     0.93e-6,   0.00e-6 },

   /* 11-20 */
      {{ 1,  0,  0, -2,  0,  0,  0,  0},     0.68e-6,   0.00e-6 },
      {{ 0,  0,  2, -2,  1,  0,  0,  0},    -0.55e-6,   0.00e-6 },
      {{ 1,  0, -2,  0, -2,  0,  0,  0},     0.53e-6,   0.00e-6 },
      {{ 0,  0,  0,  2,  0,  0,  0,  0},    -0.27e-6,   0.00e-6 },
      {{ 1,  0,  0,  0,  1,  0,  0,  0},    -0.27e-6,   0.00e-6 },
      {{ 1,  0, -2, -2, -2,  0,  0,  0},    -0.26e-6,   0.00e-6 },
      {{ 1,  0,  0,  0, -1,  0,  0,  0},    -0.25e-6,   0.00e-6 },
      {{ 1,  0,  2,  0,  1,  0,  0,  0},     0.22e-6,   0.00e-6 },
      {{ 2,  0,  0, -2,  0,  0,  0,  0},    -0.21e-6,   0.00e-6 },
      {{ 2,  0, -2,  0, -1,  0,  0,  0},     0.20e-6,   0.00e-6 },

   /* 21-25 */
      {{ 0,  0,  2,  2,  2,  0,  0,  0},     0.17e-6,   0.00e-6 },
      {{ 2,  0,  2,  0,  2,  0,  0,  0},     0.13e-6,   0.00e-6 },
      {{ 2,  0,  0,  0,  0,  0,  0,  0},    -0.13e-6,   0.00e-6 },
      {{ 1,  0,  2, -2,  2,  0,  0,  0},    -0.12e-6,   0.00e-6 },
      {{ 0,  0,  2,  0,  0,  0,  0,  0},    -0.11e-6,   0.00e-6 }
   };

/* Terms of order t^3 */
   static const TERM s3[] = {

   /* 1-4 */
      {{ 0,  0,  0,  0,  1,  0,  0,  0},     0.30e-6, -23.42e-6 },
      {{ 0,  0,  2, -2,  2,  0,  0,  0},    -0.03e-6,  -1.46e-6 },
      {{ 0,  0,  2,  0,  2,  0,  0,  0},    -0.01e-6,  -0.25e-6 },
      {{ 0,  0,  0,  0,  2,  0,  0,  0},     0.00e-6,   0.23e-6 }
   };

/* Terms of order t^4 */
   static const TERM s4[] = {

   /* 1-1 */
      {{ 0,  0,  0,  0,  1,  0,  0,  0},    -0.26e-6,  -0.01e-6 }
   };

/* Number of terms in the series */
   static const int NS0 = (int) (sizeof s0 / sizeof (TERM));
   static const int NS1 = (int) (sizeof s1 / sizeof (TERM));
   static const int NS2 = (int) (sizeof s2 / sizeof (TERM));
   static const int NS3 = (int) (sizeof s3 / sizeof (TERM));
   static const int NS4 = (int) (sizeof s4 / sizeof (TERM));

/*--------------------------------------------------------------------*/

/* Interval between fundamental epoch J2000.0 and current date (JC). */
   t = ((date1 - DJ00) + date2) / DJC;

/* Fundamental Arguments (from IERS Conventions 2003) */

/* Mean anomaly of the Moon. */
   fa[0] = eraFal03(t);

/* Mean anomaly of the Sun. */
   fa[1] = eraFalp03(t);

/* Mean longitude of the Moon minus that of the ascending node. */
   fa[2] = eraFaf03(t);

/* Mean elongation of the Moon from the Sun. */
   fa[3] = eraFad03(t);

/* Mean longitude of the ascending node of the Moon. */
   fa[4] = eraFaom03(t);

/* Mean longitude of Venus. */
   fa[5] = eraFave03(t);

/* Mean longitude of Earth. */
   fa[6] = eraFae03(t);

/* General precession in longitude. */
   fa[7] = eraFapa03(t);

/* Evaluate s. */
   w0 = sp[0];
   w1 = sp[1];
   w2 = sp[2];
   w3 = sp[3];
   w4 = sp[4];
   w5 = sp[5];

   for (i = NS0-1; i >= 0; i--) {
   a = 0.0;
   for (j = 0; j < 8; j++) {
      a += (double)s0[i].nfa[j] * fa[j];
   }
   w0 += s0[i].s * sin(a) + s0[i].c * cos(a);
   }

   for (i = NS1-1; i >= 0; i--) {
      a = 0.0;
      for (j = 0; j < 8; j++) {
         a += (double)s1[i].nfa[j] * fa[j];
      }
      w1 += s1[i].s * sin(a) + s1[i].c * cos(a);
   }

   for (i = NS2-1; i >= 0; i--) {
      a = 0.0;
      for (j = 0; j < 8; j++) {
         a += (double)s2[i].nfa[j] * fa[j];
      }
      w2 += s2[i].s * sin(a) + s2[i].c * cos(a);
   }

   for (i = NS3-1; i >= 0; i--) {
      a = 0.0;
      for (j = 0; j < 8; j++) {
         a += (double)s3[i].nfa[j] * fa[j];
      }
      w3 += s3[i].s * sin(a) + s3[i].c * cos(a);
   }

   for (i = NS4-1; i >= 0; i--) {
      a = 0.0;
      for (j = 0; j < 8; j++) {
         a += (double)s4[i].nfa[j] * fa[j];
      }
      w4 += s4[i].s * sin(a) + s4[i].c * cos(a);
   }

   s = (w0 +
       (w1 +
       (w2 +
       (w3 +
       (w4 +
        w5 * t) * t) * t) * t) * t) * DAS2R - x*y/2.0;

   return s;

}

double eraS06a(double date1, double date2)
/*
**  - - - - - - - -
**   e r a S 0 6 a
**  - - - - - - - -
**
**  The CIO locator s, positioning the Celestial Intermediate Origin on
**  the equator of the Celestial Intermediate Pole, using the IAU 2006
**  precession and IAU 2000A nutation models.
**
**  Given:
**     date1,date2  double    TT as a 2-part Julian Date (Note 1)
**
**  Returned (function value):
**                  double    the CIO locator s in radians (Note 2)
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
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
**  2) The CIO locator s is the difference between the right ascensions
**     of the same point in two systems.  The two systems are the GCRS
**     and the CIP,CIO, and the point is the ascending node of the
**     CIP equator.  The CIO locator s remains a small fraction of
**     1 arcsecond throughout 1900-2100.
**
**  3) The series used to compute s is in fact for s+XY/2, where X and Y
**     are the x and y components of the CIP unit vector;  this series is
**     more compact than a direct series for s would be.  The present
**     function uses the full IAU 2000A nutation model when predicting
**     the CIP position.
**
**  Called:
**     eraPnm06a    classical NPB matrix, IAU 2006/2000A
**     eraBpn2xy    extract CIP X,Y coordinates from NPB matrix
**     eraS06       the CIO locator s, given X,Y, IAU 2006
**
**  References:
**
**     Capitaine, N., Chapront, J., Lambert, S. and Wallace, P.,
**     "Expressions for the Celestial Intermediate Pole and Celestial
**     Ephemeris Origin consistent with the IAU 2000A precession-
**     nutation model", Astron.Astrophys. 400, 1145-1154 (2003)
**
**     n.b. The celestial ephemeris origin (CEO) was renamed "celestial
**          intermediate origin" (CIO) by IAU 2006 Resolution 2.
**
**     Capitaine, N. & Wallace, P.T., 2006, Astron.Astrophys. 450, 855
**
**     McCarthy, D. D., Petit, G. (eds.), 2004, IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG
**
**     Wallace, P.T. & Capitaine, N., 2006, Astron.Astrophys. 459, 981
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double rnpb[3][3], x, y, s;


/* Bias-precession-nutation-matrix, IAU 20006/2000A. */
   eraPnm06a(date1, date2, rnpb);

/* Extract the CIP coordinates. */
   eraBpn2xy(rnpb, &x, &y);

/* Compute the CIO locator s, given the CIP coordinates. */
   s = eraS06(date1, date2, x, y);

   return s;

}

void eraS2c(double theta, double phi, double c[3])
/*
**  - - - - - - -
**   e r a S 2 c
**  - - - - - - -
**
**  Convert spherical coordinates to Cartesian.
**
**  Given:
**     theta    double       longitude angle (radians)
**     phi      double       latitude angle (radians)
**
**  Returned:
**     c        double[3]    direction cosines
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double cp;


   cp = cos(phi);
   c[0] = cos(theta) * cp;
   c[1] = sin(theta) * cp;
   c[2] = sin(phi);

   return;

}

void eraS2p(double theta, double phi, double r, double p[3])
/*
**  - - - - - - -
**   e r a S 2 p
**  - - - - - - -
**
**  Convert spherical polar coordinates to p-vector.
**
**  Given:
**     theta   double       longitude angle (radians)
**     phi     double       latitude angle (radians)
**     r       double       radial distance
**
**  Returned:
**     p       double[3]    Cartesian coordinates
**
**  Called:
**     eraS2c       spherical coordinates to unit vector
**     eraSxp       multiply p-vector by scalar
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double u[3];


   eraS2c(theta, phi, u);
   eraSxp(r, u, p);

   return;

}

void eraS2pv(double theta, double phi, double r,
             double td, double pd, double rd,
             double pv[2][3])
/*
**  - - - - - - - -
**   e r a S 2 p v
**  - - - - - - - -
**
**  Convert position/velocity from spherical to Cartesian coordinates.
**
**  Given:
**     theta    double          longitude angle (radians)
**     phi      double          latitude angle (radians)
**     r        double          radial distance
**     td       double          rate of change of theta
**     pd       double          rate of change of phi
**     rd       double          rate of change of r
**
**  Returned:
**     pv       double[2][3]    pv-vector
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double st, ct, sp, cp, rcp, x, y, rpd, w;


   st = sin(theta);
   ct = cos(theta);
   sp = sin(phi);
   cp = cos(phi);
   rcp = r * cp;
   x = rcp * ct;
   y = rcp * st;
   rpd = r * pd;
   w = rpd*sp - cp*rd;

   pv[0][0] = x;
   pv[0][1] = y;
   pv[0][2] = r * sp;
   pv[1][0] = -y*td - w*ct;
   pv[1][1] =  x*td - w*st;
   pv[1][2] = rpd*cp + sp*rd;

   return;

}

void eraS2xpv(double s1, double s2, double pv[2][3], double spv[2][3])
/*
**  - - - - - - - - -
**   e r a S 2 x p v
**  - - - - - - - - -
**
**  Multiply a pv-vector by two scalars.
**
**  Given:
**     s1     double         scalar to multiply position component by
**     s2     double         scalar to multiply velocity component by
**     pv     double[2][3]   pv-vector
**
**  Returned:
**     spv    double[2][3]   pv-vector: p scaled by s1, v scaled by s2
**
**  Note:
**     It is permissible for pv and spv to be the same array.
**
**  Called:
**     eraSxp       multiply p-vector by scalar
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   eraSxp(s1, pv[0], spv[0]);
   eraSxp(s2, pv[1], spv[1]);

   return;

}

double eraSepp(double a[3], double b[3])
/*
**  - - - - - - - -
**   e r a S e p p
**  - - - - - - - -
**
**  Angular separation between two p-vectors.
**
**  Given:
**     a      double[3]    first p-vector (not necessarily unit length)
**     b      double[3]    second p-vector (not necessarily unit length)
**
**  Returned (function value):
**            double       angular separation (radians, always positive)
**
**  Notes:
**
**  1) If either vector is null, a zero result is returned.
**
**  2) The angular separation is most simply formulated in terms of
**     scalar product.  However, this gives poor accuracy for angles
**     near zero and pi.  The present algorithm uses both cross product
**     and dot product, to deliver full accuracy whatever the size of
**     the angle.
**
**  Called:
**     eraPxp       vector product of two p-vectors
**     eraPm        modulus of p-vector
**     eraPdp       scalar product of two p-vectors
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double axb[3], ss, cs, s;


/* Sine of angle between the vectors, multiplied by the two moduli. */
   eraPxp(a, b, axb);
   ss = eraPm(axb);

/* Cosine of the angle, multiplied by the two moduli. */
   cs = eraPdp(a, b);

/* The angle. */
   s = ((ss != 0.0) || (cs != 0.0)) ? atan2(ss, cs) : 0.0;

   return s;

}

double eraSeps(double al, double ap, double bl, double bp)
/*
**  - - - - - - - -
**   e r a S e p s
**  - - - - - - - -
**
**  Angular separation between two sets of spherical coordinates.
**
**  Given:
**     al     double       first longitude (radians)
**     ap     double       first latitude (radians)
**     bl     double       second longitude (radians)
**     bp     double       second latitude (radians)
**
**  Returned (function value):
**            double       angular separation (radians)
**
**  Called:
**     eraS2c       spherical coordinates to unit vector
**     eraSepp      angular separation between two p-vectors
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double ac[3], bc[3], s;


/* Spherical to Cartesian. */
   eraS2c(al, ap, ac);
   eraS2c(bl, bp, bc);

/* Angle between the vectors. */
   s = eraSepp(ac, bc);

   return s;

}

double eraSp00(double date1, double date2)
/*
**  - - - - - - - -
**   e r a S p 0 0
**  - - - - - - - -
**
**  The TIO locator s', positioning the Terrestrial Intermediate Origin
**  on the equator of the Celestial Intermediate Pole.
**
**  Given:
**     date1,date2  double    TT as a 2-part Julian Date (Note 1)
**
**  Returned (function value):
**                  double    the TIO locator s' in radians (Note 2)
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
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
**  2) The TIO locator s' is obtained from polar motion observations by
**     numerical integration, and so is in essence unpredictable.
**     However, it is dominated by a secular drift of about
**     47 microarcseconds per century, which is the approximation
**     evaluated by the present function.
**
**  Reference:
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double t, sp;


/* Interval between fundamental epoch J2000.0 and current date (JC). */
   t = ((date1 - DJ00) + date2) / DJC;

/* Approximate s'. */
   sp = -47e-6 * t * DAS2R;

   return sp;

}

int eraStarpm(double ra1, double dec1,
              double pmr1, double pmd1, double px1, double rv1,
              double ep1a, double ep1b, double ep2a, double ep2b,
              double *ra2, double *dec2,
              double *pmr2, double *pmd2, double *px2, double *rv2)
/*
**  - - - - - - - - - -
**   e r a S t a r p m
**  - - - - - - - - - -
**
**  Star proper motion:  update star catalog data for space motion.
**
**  Given:
**     ra1    double     right ascension (radians), before
**     dec1   double     declination (radians), before
**     pmr1   double     RA proper motion (radians/year), before
**     pmd1   double     Dec proper motion (radians/year), before
**     px1    double     parallax (arcseconds), before
**     rv1    double     radial velocity (km/s, +ve = receding), before
**     ep1a   double     "before" epoch, part A (Note 1)
**     ep1b   double     "before" epoch, part B (Note 1)
**     ep2a   double     "after" epoch, part A (Note 1)
**     ep2b   double     "after" epoch, part B (Note 1)
**
**  Returned:
**     ra2    double     right ascension (radians), after
**     dec2   double     declination (radians), after
**     pmr2   double     RA proper motion (radians/year), after
**     pmd2   double     Dec proper motion (radians/year), after
**     px2    double     parallax (arcseconds), after
**     rv2    double     radial velocity (km/s, +ve = receding), after
**
**  Returned (function value):
**            int        status:
**                          -1 = system error (should not occur)
**                           0 = no warnings or errors
**                           1 = distance overridden (Note 6)
**                           2 = excessive velocity (Note 7)
**                           4 = solution didn't converge (Note 8)
**                        else = binary logical OR of the above warnings
**
**  Notes:
**
**  1) The starting and ending TDB dates ep1a+ep1b and ep2a+ep2b are
**     Julian Dates, apportioned in any convenient way between the two
**     parts (A and B).  For example, JD(TDB)=2450123.7 could be
**     expressed in any of these ways, among others:
**
**             epna          epnb
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
**  2) In accordance with normal star-catalog conventions, the object's
**     right ascension and declination are freed from the effects of
**     secular aberration.  The frame, which is aligned to the catalog
**     equator and equinox, is Lorentzian and centered on the SSB.
**
**     The proper motions are the rate of change of the right ascension
**     and declination at the catalog epoch and are in radians per TDB
**     Julian year.
**
**     The parallax and radial velocity are in the same frame.
**
**  3) Care is needed with units.  The star coordinates are in radians
**     and the proper motions in radians per Julian year, but the
**     parallax is in arcseconds.
**
**  4) The RA proper motion is in terms of coordinate angle, not true
**     angle.  If the catalog uses arcseconds for both RA and Dec proper
**     motions, the RA proper motion will need to be divided by cos(Dec)
**     before use.
**
**  5) Straight-line motion at constant speed, in the inertial frame,
**     is assumed.
**
**  6) An extremely small (or zero or negative) parallax is interpreted
**     to mean that the object is on the "celestial sphere", the radius
**     of which is an arbitrary (large) value (see the eraStarpv
**     function for the value used).  When the distance is overridden in
**     this way, the status, initially zero, has 1 added to it.
**
**  7) If the space velocity is a significant fraction of c (see the
**     constant VMAX in the function eraStarpv),  it is arbitrarily set
**     to zero.  When this action occurs, 2 is added to the status.
**
**  8) The relativistic adjustment carried out in the eraStarpv function
**     involves an iterative calculation.  If the process fails to
**     converge within a set number of iterations, 4 is added to the
**     status.
**
**  Called:
**     eraStarpv    star catalog data to space motion pv-vector
**     eraPvu       update a pv-vector
**     eraPdp       scalar product of two p-vectors
**     eraPvstar    space motion pv-vector to star catalog data
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double pv1[2][3], tl1, dt, pv[2][3], r2, rdv, v2, c2mv2, tl2,
          pv2[2][3];
   int j1, j2, j;


/* RA,Dec etc. at the "before" epoch to space motion pv-vector. */
   j1 = eraStarpv(ra1, dec1, pmr1, pmd1, px1, rv1, pv1);

/* Light time when observed (days). */
   tl1 = eraPm(pv1[0]) / DC;

/* Time interval, "before" to "after" (days). */
   dt = (ep2a - ep1a) + (ep2b - ep1b);

/* Move star along track from the "before" observed position to the */
/* "after" geometric position. */
   eraPvu(dt + tl1, pv1, pv);

/* From this geometric position, deduce the observed light time (days) */
/* at the "after" epoch (with theoretically unneccessary error check). */
   r2 = eraPdp(pv[0], pv[0]);
   rdv = eraPdp(pv[0], pv[1]);
   v2 = eraPdp(pv[1], pv[1]);
   c2mv2 = DC*DC - v2;
   if (c2mv2 <=  0) return -1;
   tl2 = (-rdv + sqrt(rdv*rdv + c2mv2*r2)) / c2mv2;

/* Move the position along track from the observed place at the */
/* "before" epoch to the observed place at the "after" epoch. */
   eraPvu(dt + (tl1 - tl2), pv1, pv2);

/* Space motion pv-vector to RA,Dec etc. at the "after" epoch. */
   j2 = eraPvstar(pv2, ra2, dec2, pmr2, pmd2, px2, rv2);

/* Final status. */
   j = (j2 == 0) ? j1 : -1;

   return j;

}

int eraStarpv(double ra, double dec,
              double pmr, double pmd, double px, double rv,
              double pv[2][3])
/*
**  - - - - - - - - - -
**   e r a S t a r p v
**  - - - - - - - - - -
**
**  Convert star catalog coordinates to position+velocity vector.
**
**  Given (Note 1):
**     ra     double        right ascension (radians)
**     dec    double        declination (radians)
**     pmr    double        RA proper motion (radians/year)
**     pmd    double        Dec proper motion (radians/year)
**     px     double        parallax (arcseconds)
**     rv     double        radial velocity (km/s, positive = receding)
**
**  Returned (Note 2):
**     pv     double[2][3]  pv-vector (AU, AU/day)
**
**  Returned (function value):
**            int           status:
**                              0 = no warnings
**                              1 = distance overridden (Note 6)
**                              2 = excessive speed (Note 7)
**                              4 = solution didn't converge (Note 8)
**                           else = binary logical OR of the above
**
**  Notes:
**
**  1) The star data accepted by this function are "observables" for an
**     imaginary observer at the solar-system barycenter.  Proper motion
**     and radial velocity are, strictly, in terms of barycentric
**     coordinate time, TCB.  For most practical applications, it is
**     permissible to neglect the distinction between TCB and ordinary
**     "proper" time on Earth (TT/TAI).  The result will, as a rule, be
**     limited by the intrinsic accuracy of the proper-motion and
**     radial-velocity data;  moreover, the pv-vector is likely to be
**     merely an intermediate result, so that a change of time unit
**     would cancel out overall.
**
**     In accordance with normal star-catalog conventions, the object's
**     right ascension and declination are freed from the effects of
**     secular aberration.  The frame, which is aligned to the catalog
**     equator and equinox, is Lorentzian and centered on the SSB.
**
**  2) The resulting position and velocity pv-vector is with respect to
**     the same frame and, like the catalog coordinates, is freed from
**     the effects of secular aberration.  Should the "coordinate
**     direction", where the object was located at the catalog epoch, be
**     required, it may be obtained by calculating the magnitude of the
**     position vector pv[0][0-2] dividing by the speed of light in
**     AU/day to give the light-time, and then multiplying the space
**     velocity pv[1][0-2] by this light-time and adding the result to
**     pv[0][0-2].
**
**     Summarizing, the pv-vector returned is for most stars almost
**     identical to the result of applying the standard geometrical
**     "space motion" transformation.  The differences, which are the
**     subject of the Stumpff paper referenced below, are:
**
**     (i) In stars with significant radial velocity and proper motion,
**     the constantly changing light-time distorts the apparent proper
**     motion.  Note that this is a classical, not a relativistic,
**     effect.
**
**     (ii) The transformation complies with special relativity.
**
**  3) Care is needed with units.  The star coordinates are in radians
**     and the proper motions in radians per Julian year, but the
**     parallax is in arcseconds; the radial velocity is in km/s, but
**     the pv-vector result is in AU and AU/day.
**
**  4) The RA proper motion is in terms of coordinate angle, not true
**     angle.  If the catalog uses arcseconds for both RA and Dec proper
**     motions, the RA proper motion will need to be divided by cos(Dec)
**     before use.
**
**  5) Straight-line motion at constant speed, in the inertial frame,
**     is assumed.
**
**  6) An extremely small (or zero or negative) parallax is interpreted
**     to mean that the object is on the "celestial sphere", the radius
**     of which is an arbitrary (large) value (see the constant PXMIN).
**     When the distance is overridden in this way, the status,
**     initially zero, has 1 added to it.
**
**  7) If the space velocity is a significant fraction of c (see the
**     constant VMAX), it is arbitrarily set to zero.  When this action
**     occurs, 2 is added to the status.
**
**  8) The relativistic adjustment involves an iterative calculation.
**     If the process fails to converge within a set number (IMAX) of
**     iterations, 4 is added to the status.
**
**  9) The inverse transformation is performed by the function
**     eraPvstar.
**
**  Called:
**     eraS2pv      spherical coordinates to pv-vector
**     eraPm        modulus of p-vector
**     eraZp        zero p-vector
**     eraPn        decompose p-vector into modulus and direction
**     eraPdp       scalar product of two p-vectors
**     eraSxp       multiply p-vector by scalar
**     eraPmp       p-vector minus p-vector
**     eraPpp       p-vector plus p-vector
**
**  Reference:
**
**     Stumpff, P., 1985, Astron.Astrophys. 144, 232-240.
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
/* Smallest allowed parallax */
   static const double PXMIN = 1e-7;

/* Largest allowed speed (fraction of c) */
   static const double VMAX = 0.5;

/* Maximum number of iterations for relativistic solution */
   static const int IMAX = 100;

   int i, iwarn;
   double w, r, rd, rad, decd, v, x[3], usr[3], ust[3],
          vsr, vst, betst, betsr, bett, betr,
          dd, ddel, ur[3], ut[3],
          d = 0.0, del = 0.0,       /* to prevent */
          odd = 0.0, oddel = 0.0,   /* compiler   */
          od = 0.0, odel = 0.0;     /* warnings   */


/* Distance (AU). */
   if (px >= PXMIN) {
      w = px;
      iwarn = 0;
   } else {
      w = PXMIN;
      iwarn = 1;
   }
   r = DR2AS / w;

/* Radial velocity (AU/day). */
   rd = DAYSEC * rv * 1e3 / DAU;

/* Proper motion (radian/day). */
   rad = pmr / DJY;
   decd = pmd / DJY;

/* To pv-vector (AU,AU/day). */
   eraS2pv(ra, dec, r, rad, decd, rd, pv);

/* If excessive velocity, arbitrarily set it to zero. */
   v = eraPm(pv[1]);
   if (v / DC > VMAX) {
      eraZp(pv[1]);
      iwarn += 2;
   }

/* Isolate the radial component of the velocity (AU/day). */
   eraPn(pv[0], &w, x);
   vsr = eraPdp(x, pv[1]);
   eraSxp(vsr, x, usr);

/* Isolate the transverse component of the velocity (AU/day). */
   eraPmp(pv[1], usr, ust);
   vst = eraPm(ust);

/* Special-relativity dimensionless parameters. */
   betsr = vsr / DC;
   betst = vst / DC;

/* Determine the inertial-to-observed relativistic correction terms. */
   bett = betst;
   betr = betsr;
   for (i = 0; i < IMAX; i++) {
      d = 1.0 + betr;
      del = sqrt(1.0 - betr*betr - bett*bett) - 1.0;
      betr = d * betsr + del;
      bett = d * betst;
      if (i > 0) {
         dd = fabs(d - od);
         ddel = fabs(del - odel);
         if ((i > 1) && (dd >= odd) && (ddel >= oddel)) break;
         odd = dd;
         oddel = ddel;
      }
      od = d;
      odel = del;
   }
   if (i >= IMAX) iwarn += 4;

/* Replace observed radial velocity with inertial value. */
   w = (betsr != 0.0) ? d + del / betsr : 1.0;
   eraSxp(w, usr, ur);

/* Replace observed tangential velocity with inertial value. */
   eraSxp(d, ust, ut);

/* Combine the two to obtain the inertial space velocity. */
   eraPpp(ur, ut, pv[1]);

/* Return the status. */
   return iwarn;

}

void eraSxp(double s, double p[3], double sp[3])
/*
**  - - - - - - -
**   e r a S x p
**  - - - - - - -
**
**  Multiply a p-vector by a scalar.
**
**  Given:
**     s      double        scalar
**     p      double[3]     p-vector
**
**  Returned:
**     sp     double[3]     s * p
**
**  Note:
**     It is permissible for p and sp to be the same array.
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   sp[0] = s * p[0];
   sp[1] = s * p[1];
   sp[2] = s * p[2];

   return;

}

void eraSxpv(double s, double pv[2][3], double spv[2][3])
/*
**  - - - - - - - -
**   e r a S x p v
**  - - - - - - - -
**
**  Multiply a pv-vector by a scalar.
**
**  Given:
**     s       double          scalar
**     pv      double[2][3]    pv-vector
**
**  Returned:
**     spv     double[2][3]    s * pv
**
**  Note:
**     It is permissible for pv and psv to be the same array
**
**  Called:
**     eraS2xpv     multiply pv-vector by two scalars
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   eraS2xpv(s, s, pv, spv);

   return;

}

int eraTaitt(double tai1, double tai2, double *tt1, double *tt2)
/*
**  - - - - - - - - -
**   e r a T a i t t
**  - - - - - - - - -
**
**  Time scale transformation:  International Atomic Time, TAI, to
**  Terrestrial Time, TT.
**
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
**
**  References:
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**     Explanatory Supplement to the Astronomical Almanac,
**     P. Kenneth Seidelmann (ed), University Science Books (1992)
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{

/* TT minus TAI (days). */
   static const double dtat = TTMTAI/DAYSEC;


/* Result, safeguarding precision. */
   if ( tai1 > tai2 ) {
      *tt1 = tai1;
      *tt2 = tai2 + dtat;
   } else {
      *tt1 = tai1 + dtat;
      *tt2 = tai2;
   }

/* Status (always OK). */
   return 0;

}

int eraTaiut1(double tai1, double tai2, double dta,
              double *ut11, double *ut12)
/*
**  - - - - - - - - - -
**   e r a T a i u t 1
**  - - - - - - - - - -
**
**  Time scale transformation:  International Atomic Time, TAI, to
**  Universal Time, UT1.
**
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
**
**  Reference:
**
**     Explanatory Supplement to the Astronomical Almanac,
**     P. Kenneth Seidelmann (ed), University Science Books (1992)
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double dtad;


/* Result, safeguarding precision. */
   dtad = dta / DAYSEC;
   if ( tai1 > tai2 ) {
      *ut11 = tai1;
      *ut12 = tai2 + dtad;
   } else {
      *ut11 = tai1 + dtad;
      *ut12 = tai2;
   }

/* Status (always OK). */
   return 0;

}

int eraTaiutc(double tai1, double tai2, double *utc1, double *utc2)
/*
**  - - - - - - - - - -
**   e r a T a i u t c
**  - - - - - - - - - -
**
**  Time scale transformation:  International Atomic Time, TAI, to
**  Coordinated Universal Time, UTC.
**
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
**
**  Called:
**     eraJd2cal    JD to Gregorian calendar
**     eraDat       delta(AT) = TAI-UTC
**     eraCal2jd    Gregorian calendar to JD
**
**  References:
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**     Explanatory Supplement to the Astronomical Almanac,
**     P. Kenneth Seidelmann (ed), University Science Books (1992)
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   int big1;
   int i, iy, im, id, js;
   double a1, a2, d1, dats1, d2, fd, ddats, dats2, datd, as1, as2, da;


/* Put the two parts of the TAI into big-first order. */
   big1 = ( tai1 >= tai2 );
   if ( big1 ) {
      a1 = tai1;
      a2 = tai2;
   } else {
      a1 = tai2;
      a2 = tai1;
   }

/* See if the TAI can possibly be in a leap-second day. */
   d1 = a1;
   dats1 = 0.0;
   for ( i = -1; i <= 3; i++ ) {
      d2 = a2 + (double) i;
      if ( eraJd2cal(d1, d2, &iy, &im, &id, &fd) ) return -1;
      js = eraDat(iy, im, id, 0.0, &dats2);
      if ( js < 0 ) return -1;
      if ( i == -1 ) dats1 = dats2;
      ddats = dats2 - dats1;
      datd = dats1 / DAYSEC;
      if ( fabs(ddats) >= 0.5 ) {

      /* Yes.  Get TAI for the start of the UTC day that */
      /* ends in a leap. */
         if ( eraCal2jd(iy, im, id, &d1, &d2) ) return -1;
         as1 = d1;
         as2 = d2 - 1.0 + datd;

      /* Is the TAI after this point? */
         da = a1 - as1;
         da = da + ( a2 - as2 );
         if ( da > 0 ) {

         /* Yes:  fraction of the current UTC day that has elapsed. */
            fd = da * DAYSEC / ( DAYSEC + ddats );

         /* Ramp TAI-UTC to bring about ERFA's JD(UTC) convention. */
            datd += ddats * ( fd <= 1.0 ? fd : 1.0 ) / DAYSEC;
         }

      /* Done. */
         break;
      }
      dats1 = dats2;
   }

/* Subtract the (possibly adjusted) TAI-UTC from TAI to give UTC. */
   a2 -= datd;

/* Return the UTC result, preserving the TAI order. */
   if ( big1 ) {
      *utc1 = a1;
      *utc2 = a2;
   } else {
      *utc1 = a2;
      *utc2 = a1;
   }

/* Status. */
   return js;

}

int eraTcbtdb(double tcb1, double tcb2, double *tdb1, double *tdb2)
/*
**  - - - - - - - - - -
**   e r a T c b t d b
**  - - - - - - - - - -
**
**  Time scale transformation:  Barycentric Coordinate Time, TCB, to
**  Barycentric Dynamical Time, TDB.
**
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
**
**  Reference:
**
**     IAU 2006 Resolution B3
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{

/* 1977 Jan 1 00:00:32.184 TT, as two-part JD */
   static const double t77td = DJM0 + DJM77;
   static const double t77tf = TTMTAI/DAYSEC;

/* TDB (days) at TAI 1977 Jan 1.0 */
   static const double tdb0 = TDB0/DAYSEC;

   double d;


/* Result, safeguarding precision. */
   if ( tcb1 > tcb2 ) {
      d = tcb1 - t77td;
      *tdb1 = tcb1;
      *tdb2 = tcb2 + tdb0 - ( d + ( tcb2 - t77tf ) ) * ELB;
   } else {
      d = tcb2 - t77td;
      *tdb1 = tcb1 + tdb0 - ( d + ( tcb1 - t77tf ) ) * ELB;
      *tdb2 = tcb2;
   }

/* Status (always OK). */
   return 0;

}

int eraTcgtt(double tcg1, double tcg2, double *tt1, double *tt2)
/*
**  - - - - - - - - -
**   e r a T c g t t
**  - - - - - - - - -
**
**  Time scale transformation:  Geocentric Coordinate Time, TCG, to
**  Terrestrial Time, TT.
**
**  Given:
**     tcg1,tcg2  double    TCG as a 2-part Julian Date
**
**  Returned:
**     tt1,tt2    double    TT as a 2-part Julian Date
**
**  Returned (function value):
**                int       status:  0 = OK
**
**  Note:
**
**     tcg1+tcg2 is Julian Date, apportioned in any convenient way
**     between the two arguments, for example where tcg1 is the Julian
**     Day Number and tcg22 is the fraction of a day.  The returned
**     tt1,tt2 follow suit.
**
**  References:
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),.
**     IERS Technical Note No. 32, BKG (2004)
**
**     IAU 2000 Resolution B1.9
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{

/* 1977 Jan 1 00:00:32.184 TT, as MJD */
   static const double t77t = DJM77 + TTMTAI/DAYSEC;


/* Result, safeguarding precision. */
   if ( tcg1 > tcg2 ) {
      *tt1 = tcg1;
      *tt2 = tcg2 - ( ( tcg1 - DJM0 ) + ( tcg2 - t77t ) ) * ELG;
   } else {
      *tt1 = tcg1 - ( ( tcg2 - DJM0 ) + ( tcg1 - t77t ) ) * ELG;
      *tt2 = tcg2;
   }

/* OK status. */
   return 0;

}

int eraTdbtcb(double tdb1, double tdb2, double *tcb1, double *tcb2)
/*
**  - - - - - - - - - -
**   e r a T d b t c b
**  - - - - - - - - - -
**
**  Time scale transformation:  Barycentric Dynamical Time, TDB, to
**  Barycentric Coordinate Time, TCB.
**
**  Given:
**     tdb1,tdb2  double    TDB as a 2-part Julian Date
**
**  Returned:
**     tcb1,tcb2  double    TCB as a 2-part Julian Date
**
**  Returned (function value):
**                int       status:  0 = OK
**
**  Notes:
**
**  1) tdb1+tdb2 is Julian Date, apportioned in any convenient way
**     between the two arguments, for example where tdb1 is the Julian
**     Day Number and tdb2 is the fraction of a day.  The returned
**     tcb1,tcb2 follow suit.
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
**
**  Reference:
**
**     IAU 2006 Resolution B3
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{

/* 1977 Jan 1 00:00:32.184 TT, as two-part JD */
   static const double t77td = DJM0 + DJM77;
   static const double t77tf = TTMTAI/DAYSEC;

/* TDB (days) at TAI 1977 Jan 1.0 */
   static const double tdb0 = TDB0/DAYSEC;

/* TDB to TCB rate */
   static const double elbb = ELB/(1.0-ELB);

   double d, f;


/* Result, preserving date format but safeguarding precision. */
   if ( tdb1 > tdb2 ) {
      d = t77td - tdb1;
      f  = tdb2 - tdb0;
      *tcb1 = tdb1;
      *tcb2 = f - ( d - ( f - t77tf ) ) * elbb;
   } else {
      d = t77td - tdb2;
      f  = tdb1 - tdb0;
      *tcb1 = f + ( d - ( f - t77tf ) ) * elbb;
      *tcb2 = tdb2;
   }

/* Status (always OK). */
   return 0;

}

int eraTdbtt(double tdb1, double tdb2, double dtr,
             double *tt1, double *tt2 )
/*
**  - - - - - - - - -
**   e r a T d b t t
**  - - - - - - - - -
**
**  Time scale transformation:  Barycentric Dynamical Time, TDB, to
**  Terrestrial Time, TT.
**
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
**
**  References:
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**     IAU 2006 Resolution 3
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double dtrd;


/* Result, safeguarding precision. */
   dtrd = dtr / DAYSEC;
   if ( tdb1 > tdb2 ) {
      *tt1 = tdb1;
      *tt2 = tdb2 - dtrd;
   } else {
      *tt1 = tdb1 - dtrd;
      *tt2 = tdb2;
   }

/* Status (always OK). */
   return 0;

}

int eraTf2a(char s, int ihour, int imin, double sec, double *rad)
/*
**  - - - - - - - -
**   e r a T f 2 a
**  - - - - - - - -
**
**  Convert hours, minutes, seconds to radians.
**
**  Given:
**     s         char    sign:  '-' = negative, otherwise positive
**     ihour     int     hours
**     imin      int     minutes
**     sec       double  seconds
**
**  Returned:
**     rad       double  angle in radians
**
**  Returned (function value):
**               int     status:  0 = OK
**                                1 = ihour outside range 0-23
**                                2 = imin outside range 0-59
**                                3 = sec outside range 0-59.999...
**
**  Notes:
**
**  1)  The result is computed even if any of the range checks fail.
**
**  2)  Negative ihour, imin and/or sec produce a warning status, but
**      the absolute value is used in the conversion.
**
**  3)  If there are multiple errors, the status value reflects only the
**      first, the smallest taking precedence.
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{

/* Compute the interval. */
   *rad  = ( s == '-' ? -1.0 : 1.0 ) *
           ( 60.0 * ( 60.0 * ( (double) abs(ihour) ) +
                             ( (double) abs(imin) ) ) +
                                        fabs(sec) ) * DS2R;

/* Validate arguments and return status. */
   if ( ihour < 0 || ihour > 23 ) return 1;
   if ( imin < 0 || imin > 59 ) return 2;
   if ( sec < 0.0 || sec >= 60.0 ) return 3;
   return 0;

}

int eraTf2d(char s, int ihour, int imin, double sec, double *days)
/*
**  - - - - - - - -
**   e r a T f 2 d
**  - - - - - - - -
**
**  Convert hours, minutes, seconds to days.
**
**  Given:
**     s         char    sign:  '-' = negative, otherwise positive
**     ihour     int     hours
**     imin      int     minutes
**     sec       double  seconds
**
**  Returned:
**     days      double  interval in days
**
**  Returned (function value):
**               int     status:  0 = OK
**                                1 = ihour outside range 0-23
**                                2 = imin outside range 0-59
**                                3 = sec outside range 0-59.999...
**
**  Notes:
**
**  1)  The result is computed even if any of the range checks fail.
**
**  2)  Negative ihour, imin and/or sec produce a warning status, but
**      the absolute value is used in the conversion.
**
**  3)  If there are multiple errors, the status value reflects only the
**      first, the smallest taking precedence.
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{

/* Compute the interval. */
   *days  = ( s == '-' ? -1.0 : 1.0 ) *
            ( 60.0 * ( 60.0 * ( (double) abs(ihour) ) +
                              ( (double) abs(imin) ) ) +
                                         fabs(sec) ) / DAYSEC;

/* Validate arguments and return status. */
   if ( ihour < 0 || ihour > 23 ) return 1;
   if ( imin < 0 || imin > 59 ) return 2;
   if ( sec < 0.0 || sec >= 60.0 ) return 3;
   return 0;

}

void eraTr(double r[3][3], double rt[3][3])
/*
**  - - - - - -
**   e r a T r
**  - - - - - -
**
**  Transpose an r-matrix.
**
**  Given:
**     r        double[3][3]    r-matrix
**
**  Returned:
**     rt       double[3][3]    transpose
**
**  Note:
**     It is permissible for r and rt to be the same array.
**
**  Called:
**     eraCr        copy r-matrix
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double wm[3][3];
   int i, j;


   for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++) {
         wm[i][j] = r[j][i];
      }
   }
   eraCr(wm, rt);

   return;

}

void eraTrxp(double r[3][3], double p[3], double trp[3])
/*
**  - - - - - - - -
**   e r a T r x p
**  - - - - - - - -
**
**  Multiply a p-vector by the transpose of an r-matrix.
**
**  Given:
**     r        double[3][3]   r-matrix
**     p        double[3]      p-vector
**
**  Returned:
**     trp      double[3]      r * p
**
**  Note:
**     It is permissible for p and trp to be the same array.
**
**  Called:
**     eraTr        transpose r-matrix
**     eraRxp       product of r-matrix and p-vector
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double tr[3][3];


/* Transpose of matrix r. */
   eraTr(r, tr);

/* Matrix tr * vector p -> vector trp. */
   eraRxp(tr, p, trp);

   return;

}

void eraTrxpv(double r[3][3], double pv[2][3], double trpv[2][3])
/*
**  - - - - - - - - -
**   e r a T r x p v
**  - - - - - - - - -
**
**  Multiply a pv-vector by the transpose of an r-matrix.
**
**  Given:
**     r        double[3][3]    r-matrix
**     pv       double[2][3]    pv-vector
**
**  Returned:
**     trpv     double[2][3]    r * pv
**
**  Note:
**     It is permissible for pv and trpv to be the same array.
**
**  Called:
**     eraTr        transpose r-matrix
**     eraRxpv      product of r-matrix and pv-vector
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double tr[3][3];


/* Transpose of matrix r. */
   eraTr(r, tr);

/* Matrix tr * vector pv -> vector trpv. */
   eraRxpv(tr, pv, trpv);

   return;

}

int eraTttai(double tt1, double tt2, double *tai1, double *tai2)
/*
**  - - - - - - - - -
**   e r a T t t a i
**  - - - - - - - - -
**
**  Time scale transformation:  Terrestrial Time, TT, to International
**  Atomic Time, TAI.
**
**  Given:
**     tt1,tt2    double    TT as a 2-part Julian Date
**
**  Returned:
**     tai1,tai2  double    TAI as a 2-part Julian Date
**
**  Returned (function value):
**                int       status:  0 = OK
**
**  Note:
**
**     tt1+tt2 is Julian Date, apportioned in any convenient way between
**     the two arguments, for example where tt1 is the Julian Day Number
**     and tt2 is the fraction of a day.  The returned tai1,tai2 follow
**     suit.
**
**  References:
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**     Explanatory Supplement to the Astronomical Almanac,
**     P. Kenneth Seidelmann (ed), University Science Books (1992)
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{

/* TT minus TAI (days). */
   static const double dtat = TTMTAI/DAYSEC;


/* Result, safeguarding precision. */
   if ( tt1 > tt2 ) {
      *tai1 = tt1;
      *tai2 = tt2 - dtat;
   } else {
      *tai1 = tt1 - dtat;
      *tai2 = tt2;
   }

/* Status (always OK). */
   return 0;

}

int eraTttcg(double tt1, double tt2, double *tcg1, double *tcg2)
/*
**  - - - - - - - - -
**   e r a T t t c g
**  - - - - - - - - -
**
**  Time scale transformation:  Terrestrial Time, TT, to Geocentric
**  Coordinate Time, TCG.
**
**  Given:
**     tt1,tt2    double    TT as a 2-part Julian Date
**
**  Returned:
**     tcg1,tcg2  double    TCG as a 2-part Julian Date
**
**  Returned (function value):
**                int       status:  0 = OK
**
**  Note:
**
**     tt1+tt2 is Julian Date, apportioned in any convenient way between
**     the two arguments, for example where tt1 is the Julian Day Number
**     and tt2 is the fraction of a day.  The returned tcg1,tcg2 follow
**     suit.
**
**  References:
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**     IAU 2000 Resolution B1.9
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{

/* 1977 Jan 1 00:00:32.184 TT, as MJD */
   static const double t77t = DJM77 + TTMTAI/DAYSEC;

/* TT to TCG rate */
   static const double elgg = ELG/(1.0-ELG);


/* Result, safeguarding precision. */
   if ( tt1 > tt2 ) {
      *tcg1 = tt1;
      *tcg2 = tt2 + ( ( tt1 - DJM0 ) + ( tt2 - t77t ) ) * elgg;
   } else {
      *tcg1 = tt1 + ( ( tt2 - DJM0 ) + ( tt1 - t77t ) ) * elgg;
      *tcg2 = tt2;
   }

/* Status (always OK). */
   return 0;

}

int eraTttdb(double tt1, double tt2, double dtr,
             double *tdb1, double *tdb2)
/*
**  - - - - - - - - -
**   e r a T t t d b
**  - - - - - - - - -
**
**  Time scale transformation:  Terrestrial Time, TT, to Barycentric
**  Dynamical Time, TDB.
**
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
**
**  References:
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**     IAU 2006 Resolution 3
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double dtrd;


/* Result, safeguarding precision. */
   dtrd = dtr / DAYSEC;
   if ( tt1 > tt2 ) {
      *tdb1 = tt1;
      *tdb2 = tt2 + dtrd;
   } else {
      *tdb1 = tt1 + dtrd;
      *tdb2 = tt2;
   }

/* Status (always OK). */
   return 0;

}

int eraTtut1(double tt1, double tt2, double dt,
             double *ut11, double *ut12)
/*
**  - - - - - - - - -
**   e r a T t u t 1
**  - - - - - - - - -
**
**  Time scale transformation:  Terrestrial Time, TT, to Universal Time,
**  UT1.
**
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
**
**  Reference:
**
**     Explanatory Supplement to the Astronomical Almanac,
**     P. Kenneth Seidelmann (ed), University Science Books (1992)
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double dtd;


/* Result, safeguarding precision. */
   dtd = dt / DAYSEC;
   if ( tt1 > tt2 ) {
      *ut11 = tt1;
      *ut12 = tt2 - dtd;
   } else {
      *ut11 = tt1 - dtd;
      *ut12 = tt2;
   }

/* Status (always OK). */
   return 0;

}

int eraUt1tai(double ut11, double ut12, double dta,
              double *tai1, double *tai2)
/*
**  - - - - - - - - - -
**   e r a U t 1 t a i
**  - - - - - - - - - -
**
**  Time scale transformation:  Universal Time, UT1, to International
**  Atomic Time, TAI.
**
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
**
**  Reference:
**
**     Explanatory Supplement to the Astronomical Almanac,
**     P. Kenneth Seidelmann (ed), University Science Books (1992)
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double dtad;


/* Result, safeguarding precision. */
   dtad = dta / DAYSEC;
   if ( ut11 > ut12 ) {
      *tai1 = ut11;
      *tai2 = ut12 - dtad;
   } else {
      *tai1 = ut11 - dtad;
      *tai2 = ut12;
   }

/* Status (always OK). */
   return 0;

}

int eraUt1tt(double ut11, double ut12, double dt,
             double *tt1, double *tt2)
/*
**  - - - - - - - - -
**   e r a U t 1 t t
**  - - - - - - - - -
**
**  Time scale transformation:  Universal Time, UT1, to Terrestrial
**  Time, TT.
**
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
**
**  Reference:
**
**     Explanatory Supplement to the Astronomical Almanac,
**     P. Kenneth Seidelmann (ed), University Science Books (1992)
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double dtd;


/* Result, safeguarding precision. */
   dtd = dt / DAYSEC;
   if ( ut11 > ut12 ) {
      *tt1 = ut11;
      *tt2 = ut12 + dtd;
   } else {
      *tt1 = ut11 + dtd;
      *tt2 = ut12;
   }

/* Status (always OK). */
   return 0;

}

int eraUt1utc(double ut11, double ut12, double dut1,
              double *utc1, double *utc2)
/*
**  - - - - - - - - - -
**   e r a U t 1 u t c
**  - - - - - - - - - -
**
**  Time scale transformation:  Universal Time, UT1, to Coordinated
**  Universal Time, UTC.
**
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
**
**  Called:
**     eraJd2cal    JD to Gregorian calendar
**     eraDat       delta(AT) = TAI-UTC
**     eraCal2jd    Gregorian calendar to JD
**
**  References:
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**     Explanatory Supplement to the Astronomical Almanac,
**     P. Kenneth Seidelmann (ed), University Science Books (1992)
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   int big1;
   int i, iy, im, id, js;
   double duts, u1, u2, d1, dats1, d2, fd, dats2, ddats, us1, us2, du;


/* UT1-UTC in seconds. */
   duts = dut1;

/* Put the two parts of the UT1 into big-first order. */
   big1 = ( ut11 >= ut12 );
   if ( big1 ) {
      u1 = ut11;
      u2 = ut12;
   } else {
      u1 = ut12;
      u2 = ut11;
   }

/* See if the UT1 can possibly be in a leap-second day. */
   d1 = u1;
   dats1 = 0;
   for ( i = -1; i <= 3; i++ ) {
      d2 = u2 + (double) i;
      if ( eraJd2cal(d1, d2, &iy, &im, &id, &fd) ) return -1;
      js = eraDat(iy, im, id, 0.0, &dats2);
      if ( js < 0 ) return -1;
      if ( i == - 1 ) dats1 = dats2;
      ddats = dats2 - dats1;
      if ( fabs(ddats) >= 0.5 ) {

      /* Yes, leap second nearby: ensure UT1-UTC is "before" value. */
         if ( ddats * duts >= 0 ) duts -= ddats;

      /* UT1 for the start of the UTC day that ends in a leap. */
         if ( eraCal2jd(iy, im, id, &d1, &d2) ) return -1;
         us1 = d1;
         us2 = d2 - 1.0 + duts/DAYSEC;

      /* Is the UT1 after this point? */
         du = u1 - us1;
         du += u2 - us2;
         if ( du > 0 ) {

         /* Yes:  fraction of the current UTC day that has elapsed. */
            fd = du * DAYSEC / ( DAYSEC + ddats );

         /* Ramp UT1-UTC to bring about ERFA's JD(UTC) convention. */
            duts += ddats * ( fd <= 1.0 ? fd : 1.0 );
         }

      /* Done. */
         break;
      }
      dats1 = dats2;
   }

/* Subtract the (possibly adjusted) UT1-UTC from UT1 to give UTC. */
   u2 -= duts / DAYSEC;

/* Result, safeguarding precision. */
   if ( big1 ) {
      *utc1 = u1;
      *utc2 = u2;
   } else {
      *utc1 = u2;
      *utc2 = u1;
   }

/* Status. */
   return js;

}

int eraUtctai(double utc1, double utc2, double *tai1, double *tai2)
/*
**  - - - - - - - - - -
**   e r a U t c t a i
**  - - - - - - - - - -
**
**  Time scale transformation:  Coordinated Universal Time, UTC, to
**  International Atomic Time, TAI.
**
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
**  Called:
**     eraJd2cal    JD to Gregorian calendar
**     eraDat       delta(AT) = TAI-UTC
**     eraCal2jd    Gregorian calendar to JD
**
**  References:
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**     Explanatory Supplement to the Astronomical Almanac,
**     P. Kenneth Seidelmann (ed), University Science Books (1992)
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   int big1;
   int iy, im, id, js, iyt, imt, idt;
   double u1, u2, fd, dats, fdt, datst, ddat, z1, z2, a2;


/* Put the two parts of the UTC into big-first order. */
   big1 = ( utc1 >= utc2 );
   if ( big1 ) {
      u1 = utc1;
      u2 = utc2;
   } else {
      u1 = utc2;
      u2 = utc1;
   }

/* Get TAI-UTC now. */
   if ( eraJd2cal(u1, u2, &iy, &im, &id, &fd) ) return -1;
   js = eraDat(iy, im, id, fd, &dats);
   if ( js < 0 ) return -1;

/* Get TAI-UTC tomorrow. */
   if ( eraJd2cal(u1+1.5, u2-fd, &iyt, &imt, &idt, &fdt) ) return -1;
   js = eraDat(iyt, imt, idt, fdt, &datst);
   if ( js < 0 ) return -1;

/* If today ends in a leap second, scale the fraction into SI days. */
   ddat = datst - dats;
   if ( fabs(ddat) > 0.5 ) fd += fd * ddat / DAYSEC;

/* Today's calendar date to 2-part JD. */
   if ( eraCal2jd(iy, im, id, &z1, &z2) ) return -1;

/* Assemble the TAI result, preserving the UTC split and order. */
   a2 = z1 - u1;
   a2 += z2;
   a2 += fd + dats / DAYSEC;
   if ( big1 ) {
      *tai1 = u1;
      *tai2 = a2;
   } else {
      *tai1 = a2;
      *tai2 = u1;
   }

/* Status. */
   return js;

}

int eraUtcut1(double utc1, double utc2, double dut1,
              double *ut11, double *ut12)
/*
**  - - - - - - - - - -
**   e r a U t c u t 1
**  - - - - - - - - - -
**
**  Time scale transformation:  Coordinated Universal Time, UTC, to
**  Universal Time, UT1.
**
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
**
**  References:
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**     Explanatory Supplement to the Astronomical Almanac,
**     P. Kenneth Seidelmann (ed), University Science Books (1992)
**
**  Called:
**     eraJd2cal    JD to Gregorian calendar
**     eraDat       delta(AT) = TAI-UTC
**     eraUtctai    UTC to TAI
**     eraTaiut1    TAI to UT1
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   int iy, im, id, js, jw;
   double w, dat, dta, tai1, tai2;


/* Look up TAI-UTC. */
   if ( eraJd2cal(utc1, utc2, &iy, &im, &id, &w) ) return -1;
   js = eraDat ( iy, im, id, 0.0, &dat);
   if ( js < 0 ) return -1;

/* Form UT1-TAI. */
   dta = dut1 - dat;

/* UTC to TAI to UT1. */
   jw = eraUtctai(utc1, utc2, &tai1, &tai2);
   if ( jw < 0 ) {
      return -1;
   } else if ( jw > 0 ) {
      js = jw;
   }
   if ( eraTaiut1(tai1, tai2, dta, ut11, ut12) ) return -1;

/* Status. */
   return js;

}

void eraXy06(double date1, double date2, double *x, double *y)
/*
**  - - - - - - - -
**   e r a X y 0 6
**  - - - - - - - -
**
**  X,Y coordinates of celestial intermediate pole from series based
**  on IAU 2006 precession and IAU 2000A nutation.
**
**  Given:
**     date1,date2  double     TT as a 2-part Julian Date (Note 1)
**
**  Returned:
**     x,y          double     CIP X,Y coordinates (Note 2)
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
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
**  2) The X,Y coordinates are those of the unit vector towards the
**     celestial intermediate pole.  They represent the combined effects
**     of frame bias, precession and nutation.
**
**  3) The fundamental arguments used are as adopted in IERS Conventions
**     (2003) and are from Simon et al. (1994) and Souchay et al.
**     (1999).
**
**  4) This is an alternative to the angles-based method, via the ERFA
**     function eraFw2xy and as used in eraXys06a for example.  The two
**     methods agree at the 1 microarcsecond level (at present), a
**     negligible amount compared with the intrinsic accuracy of the
**     models.  However, it would be unwise to mix the two methods
**     (angles-based and series-based) in a single application.
**
**  Called:
**     eraFal03     mean anomaly of the Moon
**     eraFalp03    mean anomaly of the Sun
**     eraFaf03     mean argument of the latitude of the Moon
**     eraFad03     mean elongation of the Moon from the Sun
**     eraFaom03    mean longitude of the Moon's ascending node
**     eraFame03    mean longitude of Mercury
**     eraFave03    mean longitude of Venus
**     eraFae03     mean longitude of Earth
**     eraFama03    mean longitude of Mars
**     eraFaju03    mean longitude of Jupiter
**     eraFasa03    mean longitude of Saturn
**     eraFaur03    mean longitude of Uranus
**     eraFane03    mean longitude of Neptune
**     eraFapa03    general accumulated precession in longitude
**
**  References:
**
**     Capitaine, N., Wallace, P.T. & Chapront, J., 2003,
**     Astron.Astrophys., 412, 567
**
**     Capitaine, N. & Wallace, P.T., 2006, Astron.Astrophys. 450, 855
**
**     McCarthy, D. D., Petit, G. (eds.), 2004, IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG
**
**     Simon, J.L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
**     Francou, G. & Laskar, J., Astron.Astrophys., 1994, 282, 663
**
**     Souchay, J., Loysel, B., Kinoshita, H., Folgueira, M., 1999,
**     Astron.Astrophys.Supp.Ser. 135, 111
**
**     Wallace, P.T. & Capitaine, N., 2006, Astron.Astrophys. 459, 981
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{

/* Maximum power of T in the polynomials for X and Y */
#define MAXPT (5)

/* Polynomial coefficients (arcsec, X then Y). */
   static const double xyp[2][MAXPT+1] = {

      {    -0.016617,
         2004.191898,
           -0.4297829,
           -0.19861834,
            0.000007578,
            0.0000059285
      },
      {    -0.006951,
           -0.025896,
          -22.4072747,
            0.00190059,
            0.001112526,
            0.0000001358
      }
   };

/* Fundamental-argument multipliers:  luni-solar terms */
   static const int mfals[][5] = {

   /* 1-10 */
      {  0,   0,   0,   0,   1 },
      {  0,   0,   2,  -2,   2 },
      {  0,   0,   2,   0,   2 },
      {  0,   0,   0,   0,   2 },
      {  0,   1,   0,   0,   0 },
      {  0,   1,   2,  -2,   2 },
      {  1,   0,   0,   0,   0 },
      {  0,   0,   2,   0,   1 },
      {  1,   0,   2,   0,   2 },
      {  0,   1,  -2,   2,  -2 },

   /* 11-20 */
      {  0,   0,   2,  -2,   1 },
      {  1,   0,  -2,   0,  -2 },
      {  1,   0,   0,  -2,   0 },
      {  1,   0,   0,   0,   1 },
      {  1,   0,   0,   0,  -1 },
      {  1,   0,  -2,  -2,  -2 },
      {  1,   0,   2,   0,   1 },
      {  2,   0,  -2,   0,  -1 },
      {  0,   0,   0,   2,   0 },
      {  0,   0,   2,   2,   2 },

   /* 21-30 */
      {  2,   0,   0,  -2,   0 },
      {  0,   2,  -2,   2,  -2 },
      {  2,   0,   2,   0,   2 },
      {  1,   0,   2,  -2,   2 },
      {  1,   0,  -2,   0,  -1 },
      {  2,   0,   0,   0,   0 },
      {  0,   0,   2,   0,   0 },
      {  0,   1,   0,   0,   1 },
      {  1,   0,   0,  -2,  -1 },
      {  0,   2,   2,  -2,   2 },

   /* 31-40 */
      {  0,   0,   2,  -2,   0 },
      {  1,   0,   0,  -2,   1 },
      {  0,   1,   0,   0,  -1 },
      {  0,   2,   0,   0,   0 },
      {  1,   0,  -2,  -2,  -1 },
      {  1,   0,   2,   2,   2 },
      {  0,   1,   2,   0,   2 },
      {  2,   0,  -2,   0,   0 },
      {  0,   0,   2,   2,   1 },
      {  0,   1,  -2,   0,  -2 },

   /* 41-50 */
      {  0,   0,   0,   2,   1 },
      {  1,   0,   2,  -2,   1 },
      {  2,   0,   0,  -2,  -1 },
      {  2,   0,   2,  -2,   2 },
      {  2,   0,   2,   0,   1 },
      {  0,   0,   0,   2,  -1 },
      {  0,   1,  -2,   2,  -1 },
      {  1,   1,   0,  -2,   0 },
      {  2,   0,   0,  -2,   1 },
      {  1,   0,   0,   2,   0 },

   /* 51-60 */
      {  0,   1,   2,  -2,   1 },
      {  1,  -1,   0,   0,   0 },
      {  0,   1,  -1,   1,  -1 },
      {  2,   0,  -2,   0,  -2 },
      {  0,   1,   0,  -2,   0 },
      {  1,   0,   0,  -1,   0 },
      {  3,   0,   2,   0,   2 },
      {  0,   0,   0,   1,   0 },
      {  1,  -1,   2,   0,   2 },
      {  1,   1,  -2,  -2,  -2 },

   /* 61-70 */
      {  1,   0,  -2,   0,   0 },
      {  2,   0,   0,   0,  -1 },
      {  0,   1,  -2,  -2,  -2 },
      {  1,   1,   2,   0,   2 },
      {  2,   0,   0,   0,   1 },
      {  1,   1,   0,   0,   0 },
      {  1,   0,  -2,   2,  -1 },
      {  1,   0,   2,   0,   0 },
      {  1,  -1,   0,  -1,   0 },
      {  1,   0,   0,   0,   2 },

   /* 71-80 */
      {  1,   0,  -1,   0,  -1 },
      {  0,   0,   2,   1,   2 },
      {  1,   0,  -2,  -4,  -2 },
      {  1,  -1,   0,  -1,  -1 },
      {  1,   0,   2,   2,   1 },
      {  0,   2,  -2,   2,  -1 },
      {  1,   0,   0,   0,  -2 },
      {  2,   0,  -2,  -2,  -2 },
      {  1,   1,   2,  -2,   2 },
      {  2,   0,  -2,  -4,  -2 },

   /* 81-90 */
      {  1,   0,  -4,   0,  -2 },
      {  2,   0,   2,  -2,   1 },
      {  1,   0,   0,  -1,  -1 },
      {  2,   0,   2,   2,   2 },
      {  3,   0,   0,   0,   0 },
      {  1,   0,   0,   2,   1 },
      {  0,   0,   2,  -2,  -1 },
      {  3,   0,   2,  -2,   2 },
      {  0,   0,   4,  -2,   2 },
      {  1,   0,   0,  -4,   0 },

   /* 91-100 */
      {  0,   1,   2,   0,   1 },
      {  2,   0,   0,  -4,   0 },
      {  1,   1,   0,  -2,  -1 },
      {  2,   0,  -2,   0,   1 },
      {  0,   0,   2,   0,  -1 },
      {  0,   1,  -2,   0,  -1 },
      {  0,   1,   0,   0,   2 },
      {  0,   0,   2,  -1,   2 },
      {  0,   0,   2,   4,   2 },
      {  2,   1,   0,  -2,   0 },

   /* 101-110 */
      {  1,   1,   0,  -2,   1 },
      {  1,  -1,   0,  -2,   0 },
      {  1,  -1,   0,  -1,  -2 },
      {  1,  -1,   0,   0,   1 },
      {  0,   1,  -2,   2,   0 },
      {  0,   1,   0,   0,  -2 },
      {  1,  -1,   2,   2,   2 },
      {  1,   0,   0,   2,  -1 },
      {  1,  -1,  -2,  -2,  -2 },
      {  3,   0,   2,   0,   1 },

   /* 111-120 */
      {  0,   1,   2,   2,   2 },
      {  1,   0,   2,  -2,   0 },
      {  1,   1,  -2,  -2,  -1 },
      {  1,   0,   2,  -4,   1 },
      {  0,   1,  -2,  -2,  -1 },
      {  2,  -1,   2,   0,   2 },
      {  0,   0,   0,   2,   2 },
      {  1,  -1,   2,   0,   1 },
      {  1,  -1,  -2,   0,  -2 },
      {  0,   1,   0,   2,   0 },

   /* 121-130 */
      {  0,   1,   2,  -2,   0 },
      {  0,   0,   0,   1,   1 },
      {  1,   0,  -2,  -2,   0 },
      {  0,   3,   2,  -2,   2 },
      {  2,   1,   2,   0,   2 },
      {  1,   1,   0,   0,   1 },
      {  2,   0,   0,   2,   0 },
      {  1,   1,   2,   0,   1 },
      {  1,   0,   0,  -2,  -2 },
      {  1,   0,  -2,   2,   0 },

   /* 131-140 */
      {  1,   0,  -1,   0,  -2 },
      {  0,   1,   0,  -2,   1 },
      {  0,   1,   0,   1,   0 },
      {  0,   0,   0,   1,  -1 },
      {  1,   0,  -2,   2,  -2 },
      {  1,  -1,   0,   0,  -1 },
      {  0,   0,   0,   4,   0 },
      {  1,  -1,   0,   2,   0 },
      {  1,   0,   2,   1,   2 },
      {  1,   0,   2,  -1,   2 },

   /* 141-150 */
      {  0,   0,   2,   1,   1 },
      {  1,   0,   0,  -2,   2 },
      {  1,   0,  -2,   0,   1 },
      {  1,   0,  -2,  -4,  -1 },
      {  0,   0,   2,   2,   0 },
      {  1,   1,   2,  -2,   1 },
      {  1,   0,  -2,   1,  -1 },
      {  0,   0,   1,   0,   1 },
      {  2,   0,  -2,  -2,  -1 },
      {  4,   0,   2,   0,   2 },

   /* 151-160 */
      {  2,  -1,   0,   0,   0 },
      {  2,   1,   2,  -2,   2 },
      {  0,   1,   2,   1,   2 },
      {  1,   0,   4,  -2,   2 },
      {  1,   1,   0,   0,  -1 },
      {  2,   0,   2,   0,   0 },
      {  2,   0,  -2,  -4,  -1 },
      {  1,   0,  -1,   0,   0 },
      {  1,   0,   0,   1,   0 },
      {  0,   1,   0,   2,   1 },

   /* 161-170 */
      {  1,   0,  -4,   0,  -1 },
      {  1,   0,   0,  -4,  -1 },
      {  2,   0,   2,   2,   1 },
      {  2,   1,   0,   0,   0 },
      {  0,   0,   2,  -3,   2 },
      {  1,   2,   0,  -2,   0 },
      {  0,   3,   0,   0,   0 },
      {  0,   0,   4,   0,   2 },
      {  0,   0,   2,  -4,   1 },
      {  2,   0,   0,  -2,  -2 },

   /* 171-180 */
      {  1,   1,  -2,  -4,  -2 },
      {  0,   1,   0,  -2,  -1 },
      {  0,   0,   0,   4,   1 },
      {  3,   0,   2,  -2,   1 },
      {  1,   0,   2,   4,   2 },
      {  1,   1,  -2,   0,  -2 },
      {  0,   0,   4,  -2,   1 },
      {  2,  -2,   0,  -2,   0 },
      {  2,   1,   0,  -2,  -1 },
      {  0,   2,   0,  -2,   0 },

   /* 181-190 */
      {  1,   0,   0,  -1,   1 },
      {  1,   1,   2,   2,   2 },
      {  3,   0,   0,   0,  -1 },
      {  2,   0,   0,  -4,  -1 },
      {  3,   0,   2,   2,   2 },
      {  0,   0,   2,   4,   1 },
      {  0,   2,  -2,  -2,  -2 },
      {  1,  -1,   0,  -2,  -1 },
      {  0,   0,   2,  -1,   1 },
      {  2,   0,   0,   2,   1 },

   /* 191-200 */
      {  1,  -1,  -2,   2,  -1 },
      {  0,   0,   0,   2,  -2 },
      {  2,   0,   0,  -4,   1 },
      {  1,   0,   0,  -4,   1 },
      {  2,   0,   2,  -4,   1 },
      {  4,   0,   2,  -2,   2 },
      {  2,   1,  -2,   0,  -1 },
      {  2,   1,  -2,  -4,  -2 },
      {  3,   0,   0,  -4,   0 },
      {  1,  -1,   2,   2,   1 },

   /* 201-210 */
      {  1,  -1,  -2,   0,  -1 },
      {  0,   2,   0,   0,   1 },
      {  1,   2,  -2,  -2,  -2 },
      {  1,   1,   0,  -4,   0 },
      {  2,   0,   0,  -2,   2 },
      {  0,   2,   2,  -2,   1 },
      {  1,   0,   2,   0,  -1 },
      {  2,   1,   0,  -2,   1 },
      {  2,  -1,  -2,   0,  -1 },
      {  1,  -1,  -2,  -2,  -1 },

   /* 211-220 */
      {  0,   1,  -2,   1,  -2 },
      {  1,   0,  -4,   2,  -2 },
      {  0,   1,   2,   2,   1 },
      {  3,   0,   0,   0,   1 },
      {  2,  -1,   2,   2,   2 },
      {  0,   1,  -2,  -4,  -2 },
      {  1,   0,  -2,  -3,  -2 },
      {  2,   0,   0,   0,   2 },
      {  1,  -1,   0,  -2,  -2 },
      {  2,   0,  -2,   2,  -1 },

   /* 221-230 */
      {  0,   2,  -2,   0,  -2 },
      {  3,   0,  -2,   0,  -1 },
      {  2,  -1,   2,   0,   1 },
      {  1,   0,  -2,  -1,  -2 },
      {  0,   0,   2,   0,   3 },
      {  2,   0,  -4,   0,  -2 },
      {  2,   1,   0,  -4,   0 },
      {  1,   1,  -2,   1,  -1 },
      {  0,   2,   2,   0,   2 },
      {  1,  -1,   2,  -2,   2 },

   /* 231-240 */
      {  1,  -1,   0,  -2,   1 },
      {  2,   1,   2,   0,   1 },
      {  1,   0,   2,  -4,   2 },
      {  1,   1,  -2,   0,  -1 },
      {  1,   1,   0,   2,   0 },
      {  1,   0,   0,  -3,   0 },
      {  2,   0,   2,  -1,   2 },
      {  0,   2,   0,   0,  -1 },
      {  2,  -1,   0,  -2,   0 },
      {  4,   0,   0,   0,   0 },

   /* 241-250 */
      {  2,   1,  -2,  -2,  -2 },
      {  0,   2,  -2,   2,   0 },
      {  1,   0,   2,   1,   1 },
      {  1,   0,  -1,   0,  -3 },
      {  3,  -1,   2,   0,   2 },
      {  2,   0,   2,  -2,   0 },
      {  1,  -2,   0,   0,   0 },
      {  2,   0,   0,   0,  -2 },
      {  1,   0,   0,   4,   0 },
      {  0,   1,   0,   1,   1 },

   /* 251-260 */
      {  1,   0,   2,   2,   0 },
      {  0,   1,   0,   2,  -1 },
      {  0,   1,   0,   1,  -1 },
      {  0,   0,   2,  -2,   3 },
      {  3,   1,   2,   0,   2 },
      {  1,   1,   2,   1,   2 },
      {  1,   1,  -2,   2,  -1 },
      {  2,  -1,   2,  -2,   2 },
      {  1,  -2,   2,   0,   2 },
      {  1,   0,   2,  -4,   0 },

   /* 261-270 */
      {  0,   0,   1,   0,   0 },
      {  1,   0,   2,  -3,   1 },
      {  1,  -2,   0,  -2,   0 },
      {  2,   0,   0,   2,  -1 },
      {  1,   1,   2,  -4,   1 },
      {  4,   0,   2,   0,   1 },
      {  0,   1,   2,   1,   1 },
      {  1,   2,   2,  -2,   2 },
      {  2,   0,   2,   1,   2 },
      {  2,   1,   2,  -2,   1 },

   /* 271-280 */
      {  1,   0,   2,  -1,   1 },
      {  1,   0,   4,  -2,   1 },
      {  1,  -1,   2,  -2,   1 },
      {  0,   1,   0,  -4,   0 },
      {  3,   0,  -2,  -2,  -2 },
      {  0,   0,   4,  -4,   2 },
      {  2,   0,  -4,  -2,  -2 },
      {  2,  -2,   0,  -2,  -1 },
      {  1,   0,   2,  -2,  -1 },
      {  2,   0,  -2,  -6,  -2 },

   /* 281-290 */
      {  1,   0,  -2,   1,  -2 },
      {  1,   0,  -2,   2,   1 },
      {  1,  -1,   0,   2,  -1 },
      {  1,   0,  -2,   1,   0 },
      {  2,  -1,   0,  -2,   1 },
      {  1,  -1,   0,   2,   1 },
      {  2,   0,  -2,  -2,   0 },
      {  1,   0,   2,  -3,   2 },
      {  0,   0,   0,   4,  -1 },
      {  2,  -1,   0,   0,   1 },

   /* 291-300 */
      {  2,   0,   4,  -2,   2 },
      {  0,   0,   2,   3,   2 },
      {  0,   1,   4,  -2,   2 },
      {  0,   1,  -2,   2,   1 },
      {  1,   1,   0,   2,   1 },
      {  1,   0,   0,   4,   1 },
      {  0,   0,   4,   0,   1 },
      {  2,   0,   0,  -3,   0 },
      {  1,   0,   0,  -1,  -2 },
      {  1,  -2,  -2,  -2,  -2 },

   /* 301-310 */
      {  3,   0,   0,   2,   0 },
      {  2,   0,   2,  -4,   2 },
      {  1,   1,  -2,  -4,  -1 },
      {  1,   0,  -2,  -6,  -2 },
      {  2,  -1,   0,   0,  -1 },
      {  2,  -1,   0,   2,   0 },
      {  0,   1,   2,  -2,  -1 },
      {  1,   1,   0,   1,   0 },
      {  1,   2,   0,  -2,  -1 },
      {  1,   0,   0,   1,  -1 },

   /* 311-320 */
      {  0,   0,   1,   0,   2 },
      {  3,   1,   2,  -2,   2 },
      {  1,   0,  -4,  -2,  -2 },
      {  1,   0,   2,   4,   1 },
      {  1,  -2,   2,   2,   2 },
      {  1,  -1,  -2,  -4,  -2 },
      {  0,   0,   2,  -4,   2 },
      {  0,   0,   2,  -3,   1 },
      {  2,   1,  -2,   0,   0 },
      {  3,   0,  -2,  -2,  -1 },

   /* 321-330 */
      {  2,   0,   2,   4,   2 },
      {  0,   0,   0,   0,   3 },
      {  2,  -1,  -2,  -2,  -2 },
      {  2,   0,   0,  -1,   0 },
      {  3,   0,   2,  -4,   2 },
      {  2,   1,   2,   2,   2 },
      {  0,   0,   3,   0,   3 },
      {  1,   1,   2,   2,   1 },
      {  2,   1,   0,   0,  -1 },
      {  1,   2,   0,  -2,   1 },

   /* 331-340 */
      {  3,   0,   2,   2,   1 },
      {  1,  -1,  -2,   2,  -2 },
      {  1,   1,   0,  -1,   0 },
      {  1,   2,   0,   0,   0 },
      {  1,   0,   4,   0,   2 },
      {  1,  -1,   2,   4,   2 },
      {  2,   1,   0,   0,   1 },
      {  1,   0,   0,   2,   2 },
      {  1,  -1,  -2,   2,   0 },
      {  0,   2,  -2,  -2,  -1 },

   /* 341-350 */
      {  2,   0,  -2,   0,   2 },
      {  5,   0,   2,   0,   2 },
      {  3,   0,  -2,  -6,  -2 },
      {  1,  -1,   2,  -1,   2 },
      {  3,   0,   0,  -4,  -1 },
      {  1,   0,   0,   1,   1 },
      {  1,   0,  -4,   2,  -1 },
      {  0,   1,   2,  -4,   1 },
      {  1,   2,   2,   0,   2 },
      {  0,   1,   0,  -2,  -2 },

   /* 351-360 */
      {  0,   0,   2,  -1,   0 },
      {  1,   0,   1,   0,   1 },
      {  0,   2,   0,  -2,   1 },
      {  3,   0,   2,   0,   0 },
      {  1,   1,  -2,   1,   0 },
      {  2,   1,  -2,  -4,  -1 },
      {  3,  -1,   0,   0,   0 },
      {  2,  -1,  -2,   0,   0 },
      {  4,   0,   2,  -2,   1 },
      {  2,   0,  -2,   2,   0 },

   /* 361-370 */
      {  1,   1,   2,  -2,   0 },
      {  1,   0,  -2,   4,  -1 },
      {  1,   0,  -2,  -2,   1 },
      {  2,   0,   2,  -4,   0 },
      {  1,   1,   0,  -2,  -2 },
      {  1,   1,  -2,  -2,   0 },
      {  1,   0,   1,  -2,   1 },
      {  2,  -1,  -2,  -4,  -2 },
      {  3,   0,  -2,   0,  -2 },
      {  0,   1,  -2,  -2,   0 },

   /* 371-380 */
      {  3,   0,   0,  -2,  -1 },
      {  1,   0,  -2,  -3,  -1 },
      {  0,   1,   0,  -4,  -1 },
      {  1,  -2,   2,  -2,   1 },
      {  0,   1,  -2,   1,  -1 },
      {  1,  -1,   0,   0,   2 },
      {  2,   0,   0,   1,   0 },
      {  1,  -2,   0,   2,   0 },
      {  1,   2,  -2,  -2,  -1 },
      {  0,   0,   4,  -4,   1 },

   /* 381-390 */
      {  0,   1,   2,   4,   2 },
      {  0,   1,  -4,   2,  -2 },
      {  3,   0,  -2,   0,   0 },
      {  2,  -1,   2,   2,   1 },
      {  0,   1,  -2,  -4,  -1 },
      {  4,   0,   2,   2,   2 },
      {  2,   0,  -2,  -3,  -2 },
      {  2,   0,   0,  -6,   0 },
      {  1,   0,   2,   0,   3 },
      {  3,   1,   0,   0,   0 },

   /* 391-400 */
      {  3,   0,   0,  -4,   1 },
      {  1,  -1,   2,   0,   0 },
      {  1,  -1,   0,  -4,   0 },
      {  2,   0,  -2,   2,  -2 },
      {  1,   1,   0,  -2,   2 },
      {  4,   0,   0,  -2,   0 },
      {  2,   2,   0,  -2,   0 },
      {  0,   1,   2,   0,   0 },
      {  1,   1,   0,  -4,   1 },
      {  1,   0,   0,  -4,  -2 },

   /* 401-410 */
      {  0,   0,   0,   1,   2 },
      {  3,   0,   0,   2,   1 },
      {  1,   1,   0,  -4,  -1 },
      {  0,   0,   2,   2,  -1 },
      {  1,   1,   2,   0,   0 },
      {  1,  -1,   2,  -4,   1 },
      {  1,   1,   0,   0,   2 },
      {  0,   0,   2,   6,   2 },
      {  4,   0,  -2,  -2,  -1 },
      {  2,   1,   0,  -4,  -1 },

   /* 411-420 */
      {  0,   0,   0,   3,   1 },
      {  1,  -1,  -2,   0,   0 },
      {  0,   0,   2,   1,   0 },
      {  1,   0,   0,   2,  -2 },
      {  3,  -1,   2,   2,   2 },
      {  3,  -1,   2,  -2,   2 },
      {  1,   0,   0,  -1,   2 },
      {  1,  -2,   2,  -2,   2 },
      {  0,   1,   0,   2,   2 },
      {  0,   1,  -2,  -1,  -2 },

   /* 421-430 */
      {  1,   1,  -2,   0,   0 },
      {  0,   2,   2,  -2,   0 },
      {  3,  -1,  -2,  -1,  -2 },
      {  1,   0,   0,  -6,   0 },
      {  1,   0,  -2,  -4,   0 },
      {  2,   1,   0,  -4,   1 },
      {  2,   0,   2,   0,  -1 },
      {  2,   0,  -4,   0,  -1 },
      {  0,   0,   3,   0,   2 },
      {  2,   1,  -2,  -2,  -1 },

   /* 431-440 */
      {  1,  -2,   0,   0,   1 },
      {  2,  -1,   0,  -4,   0 },
      {  0,   0,   0,   3,   0 },
      {  5,   0,   2,  -2,   2 },
      {  1,   2,  -2,  -4,  -2 },
      {  1,   0,   4,  -4,   2 },
      {  0,   0,   4,  -1,   2 },
      {  3,   1,   0,  -4,   0 },
      {  3,   0,   0,  -6,   0 },
      {  2,   0,   0,   2,   2 },

   /* 441-450 */
      {  2,  -2,   2,   0,   2 },
      {  1,   0,   0,  -3,   1 },
      {  1,  -2,  -2,   0,  -2 },
      {  1,  -1,  -2,  -3,  -2 },
      {  0,   0,   2,  -2,  -2 },
      {  2,   0,  -2,  -4,   0 },
      {  1,   0,  -4,   0,   0 },
      {  0,   1,   0,  -1,   0 },
      {  4,   0,   0,   0,  -1 },
      {  3,   0,   2,  -1,   2 },

   /* 451-460 */
      {  3,  -1,   2,   0,   1 },
      {  2,   0,   2,  -1,   1 },
      {  1,   2,   2,  -2,   1 },
      {  1,   1,   0,   2,  -1 },
      {  0,   2,   2,   0,   1 },
      {  3,   1,   2,   0,   1 },
      {  1,   1,   2,   1,   1 },
      {  1,   1,   0,  -1,   1 },
      {  1,  -2,   0,  -2,  -1 },
      {  4,   0,   0,  -4,   0 },

   /* 461-470 */
      {  2,   1,   0,   2,   0 },
      {  1,  -1,   0,   4,   0 },
      {  0,   1,   0,  -2,   2 },
      {  0,   0,   2,   0,  -2 },
      {  1,   0,  -1,   0,   1 },
      {  3,   0,   2,  -2,   0 },
      {  2,   0,   2,   2,   0 },
      {  1,   2,   0,  -4,   0 },
      {  1,  -1,   0,  -3,   0 },
      {  0,   1,   0,   4,   0 },

   /* 471 - 480 */
      {  0,   1,  -2,   0,   0 },
      {  2,   2,   2,  -2,   2 },
      {  0,   0,   0,   1,  -2 },
      {  0,   2,  -2,   0,  -1 },
      {  4,   0,   2,  -4,   2 },
      {  2,   0,  -4,   2,  -2 },
      {  2,  -1,  -2,   0,  -2 },
      {  1,   1,   4,  -2,   2 },
      {  1,   1,   2,  -4,   2 },
      {  1,   0,   2,   3,   2 },

   /* 481-490 */
      {  1,   0,   0,   4,  -1 },
      {  0,   0,   0,   4,   2 },
      {  2,   0,   0,   4,   0 },
      {  1,   1,  -2,   2,   0 },
      {  2,   1,   2,   1,   2 },
      {  2,   1,   2,  -4,   1 },
      {  2,   0,   2,   1,   1 },
      {  2,   0,  -4,  -2,  -1 },
      {  2,   0,  -2,  -6,  -1 },
      {  2,  -1,   2,  -1,   2 },

   /* 491-500 */
      {  1,  -2,   2,   0,   1 },
      {  1,  -2,   0,  -2,   1 },
      {  1,  -1,   0,  -4,  -1 },
      {  0,   2,   2,   2,   2 },
      {  0,   2,  -2,  -4,  -2 },
      {  0,   1,   2,   3,   2 },
      {  0,   1,   0,  -4,   1 },
      {  3,   0,   0,  -2,   1 },
      {  2,   1,  -2,   0,   1 },
      {  2,   0,   4,  -2,   1 },

   /* 501-510 */
      {  2,   0,   0,  -3,  -1 },
      {  2,  -2,   0,  -2,   1 },
      {  2,  -1,   2,  -2,   1 },
      {  1,   0,   0,  -6,  -1 },
      {  1,  -2,   0,   0,  -1 },
      {  1,  -2,  -2,  -2,  -1 },
      {  0,   1,   4,  -2,   1 },
      {  0,   0,   2,   3,   1 },
      {  2,  -1,   0,  -1,   0 },
      {  1,   3,   0,  -2,   0 },

   /* 511-520 */
      {  0,   3,   0,  -2,   0 },
      {  2,  -2,   2,  -2,   2 },
      {  0,   0,   4,  -2,   0 },
      {  4,  -1,   2,   0,   2 },
      {  2,   2,  -2,  -4,  -2 },
      {  4,   1,   2,   0,   2 },
      {  4,  -1,  -2,  -2,  -2 },
      {  2,   1,   0,  -2,  -2 },
      {  2,   1,  -2,  -6,  -2 },
      {  2,   0,   0,  -1,   1 },

   /* 521-530 */
      {  2,  -1,  -2,   2,  -1 },
      {  1,   1,  -2,   2,  -2 },
      {  1,   1,  -2,  -3,  -2 },
      {  1,   0,   3,   0,   3 },
      {  1,   0,  -2,   1,   1 },
      {  1,   0,  -2,   0,   2 },
      {  1,  -1,   2,   1,   2 },
      {  1,  -1,   0,   0,  -2 },
      {  1,  -1,  -4,   2,  -2 },
      {  0,   3,  -2,  -2,  -2 },

   /* 531-540 */
      {  0,   1,   0,   4,   1 },
      {  0,   0,   4,   2,   2 },
      {  3,   0,  -2,  -2,   0 },
      {  2,  -2,   0,   0,   0 },
      {  1,   1,   2,  -4,   0 },
      {  1,   1,   0,  -3,   0 },
      {  1,   0,   2,  -3,   0 },
      {  1,  -1,   2,  -2,   0 },
      {  0,   2,   0,   2,   0 },
      {  0,   0,   2,   4,   0 },

   /* 541-550 */
      {  1,   0,   1,   0,   0 },
      {  3,   1,   2,  -2,   1 },
      {  3,   0,   4,  -2,   2 },
      {  3,   0,   2,   1,   2 },
      {  3,   0,   0,   2,  -1 },
      {  3,   0,   0,   0,   2 },
      {  3,   0,  -2,   2,  -1 },
      {  2,   0,   4,  -4,   2 },
      {  2,   0,   2,  -3,   2 },
      {  2,   0,   0,   4,   1 },

   /* 551-560 */
      {  2,   0,   0,  -3,   1 },
      {  2,   0,  -4,   2,  -1 },
      {  2,   0,  -2,  -2,   1 },
      {  2,  -2,   2,   2,   2 },
      {  2,  -2,   0,  -2,  -2 },
      {  2,  -1,   0,   2,   1 },
      {  2,  -1,   0,   2,  -1 },
      {  1,   1,   2,   4,   2 },
      {  1,   1,   0,   1,   1 },
      {  1,   1,   0,   1,  -1 },

   /* 561-570 */
      {  1,   1,  -2,  -6,  -2 },
      {  1,   0,   0,  -3,  -1 },
      {  1,   0,  -4,  -2,  -1 },
      {  1,   0,  -2,  -6,  -1 },
      {  1,  -2,   2,   2,   1 },
      {  1,  -2,  -2,   2,  -1 },
      {  1,  -1,  -2,  -4,  -1 },
      {  0,   2,   0,   0,   2 },
      {  0,   1,   2,  -4,   2 },
      {  0,   1,  -2,   4,  -1 },

   /* 571-580 */
      {  5,   0,   0,   0,   0 },
      {  3,   0,   0,  -3,   0 },
      {  2,   2,   0,  -4,   0 },
      {  1,  -1,   2,   2,   0 },
      {  0,   1,   0,   3,   0 },
      {  4,   0,  -2,   0,  -1 },
      {  3,   0,  -2,  -6,  -1 },
      {  3,   0,  -2,  -1,  -1 },
      {  2,   1,   2,   2,   1 },
      {  2,   1,   0,   2,   1 },

   /* 581-590 */
      {  2,   0,   2,   4,   1 },
      {  2,   0,   2,  -6,   1 },
      {  2,   0,   2,  -2,  -1 },
      {  2,   0,   0,  -6,  -1 },
      {  2,  -1,  -2,  -2,  -1 },
      {  1,   2,   2,   0,   1 },
      {  1,   2,   0,   0,   1 },
      {  1,   0,   4,   0,   1 },
      {  1,   0,   2,  -6,   1 },
      {  1,   0,   2,  -4,  -1 },

   /* 591-600 */
      {  1,   0,  -1,  -2,  -1 },
      {  1,  -1,   2,   4,   1 },
      {  1,  -1,   2,  -3,   1 },
      {  1,  -1,   0,   4,   1 },
      {  1,  -1,  -2,   1,  -1 },
      {  0,   1,   2,  -2,   3 },
      {  3,   0,   0,  -2,   0 },
      {  1,   0,   1,  -2,   0 },
      {  0,   2,   0,  -4,   0 },
      {  0,   0,   2,  -4,   0 },

   /* 601-610 */
      {  0,   0,   1,  -1,   0 },
      {  0,   0,   0,   6,   0 },
      {  0,   2,   0,   0,  -2 },
      {  0,   1,  -2,   2,  -3 },
      {  4,   0,   0,   2,   0 },
      {  3,   0,   0,  -1,   0 },
      {  3,  -1,   0,   2,   0 },
      {  2,   1,   0,   1,   0 },
      {  2,   1,   0,  -6,   0 },
      {  2,  -1,   2,   0,   0 },

   /* 611-620 */
      {  1,   0,   2,  -1,   0 },
      {  1,  -1,   0,   1,   0 },
      {  1,  -1,  -2,  -2,   0 },
      {  0,   1,   2,   2,   0 },
      {  0,   0,   2,  -3,   0 },
      {  2,   2,   0,  -2,  -1 },
      {  2,  -1,  -2,   0,   1 },
      {  1,   2,   2,  -4,   1 },
      {  0,   1,   4,  -4,   2 },
      {  0,   0,   0,   3,   2 },

   /* 621-630 */
      {  5,   0,   2,   0,   1 },
      {  4,   1,   2,  -2,   2 },
      {  4,   0,  -2,  -2,   0 },
      {  3,   1,   2,   2,   2 },
      {  3,   1,   0,  -2,   0 },
      {  3,   1,  -2,  -6,  -2 },
      {  3,   0,   0,   0,  -2 },
      {  3,   0,  -2,  -4,  -2 },
      {  3,  -1,   0,  -3,   0 },
      {  3,  -1,   0,  -2,   0 },

   /* 631-640 */
      {  2,   1,   2,   0,   0 },
      {  2,   1,   2,  -4,   2 },
      {  2,   1,   2,  -2,   0 },
      {  2,   1,   0,  -3,   0 },
      {  2,   1,  -2,   0,  -2 },
      {  2,   0,   0,  -4,   2 },
      {  2,   0,   0,  -4,  -2 },
      {  2,   0,  -2,  -5,  -2 },
      {  2,  -1,   2,   4,   2 },
      {  2,  -1,   0,  -2,   2 },

   /* 641-650 */
      {  1,   3,  -2,  -2,  -2 },
      {  1,   1,   0,   0,  -2 },
      {  1,   1,   0,  -6,   0 },
      {  1,   1,  -2,   1,  -2 },
      {  1,   1,  -2,  -1,  -2 },
      {  1,   0,   2,   1,   0 },
      {  1,   0,   0,   3,   0 },
      {  1,   0,   0,  -4,   2 },
      {  1,   0,  -2,   4,  -2 },
      {  1,  -2,   0,  -1,   0 },

   /* 651-NFLS */
      {  0,   1,  -4,   2,  -1 },
      {  1,   0,  -2,   0,  -3 },
      {  0,   0,   4,  -4,   4 }
   };

/* Number of frequencies:  luni-solar */
   static const int NFLS = (int) (sizeof mfals / sizeof (int) / 5);

/* Fundamental-argument multipliers:  planetary terms */
   static const int mfapl[][14] = {

   /* 1-10 */
      {  0,  0,  1, -1,  1,  0,  0, -1,  0, -2,  5,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  0,  0,  2, -5,  0,  0, -1 },
      {  0,  0,  0,  0,  0,  0,  3, -5,  0,  0,  0,  0,  0, -2 },
      {  0,  0,  1, -1,  1,  0, -8, 12,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  2 },
      {  0,  0,  0,  0,  0,  0,  0,  4, -8,  3,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  8,-16,  4,  5,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  1,  0, -1,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  1,  0,  0, -1,  2,  0,  0,  0,  0,  0 },

   /* 11-20 */
      {  0,  0,  0,  0,  0,  0,  8,-13,  0,  0,  0,  0,  0, -1 },
      {  0,  0,  1, -1,  1,  0,  0, -1,  0,  2, -5,  0,  0,  0 },
      {  0,  0,  2, -2,  1,  0, -5,  6,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  4, -6,  0,  0,  0,  0,  0, -2 },
      {  0,  0,  0,  0,  0,  0,  0,  3,  0, -1,  0,  0,  0,  2 },
      {  0,  0,  0,  0,  0,  0,  0,  2, -8,  3,  0,  0,  0, -2 },
      {  0,  0,  0,  0,  0,  0,  2, -4,  0,  0,  0,  0,  0, -2 },
      {  0,  0,  0,  0,  0,  0,  0,  6, -8,  3,  0,  0,  0,  2 },
      {  0,  0,  0,  0,  0,  0,  0,  1, -2,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  2, -3,  0,  0,  0,  0,  0,  0 },

   /* 21-30 */
      {  0,  0,  0,  0,  0,  0,  2, -2,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  0,  2 },
      {  0,  0,  0,  0,  1,  0,  0, -4,  8, -3,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  1,  0,  0,  4, -8,  3,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  0,  0,  2, -5,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  1,  1,  0,  0,  0,  0,  0,  2 },
      {  0,  0,  1, -1,  1,  0,  0,  0, -2,  0,  0,  0,  0,  0 },
      {  2,  0,  0, -2, -1,  0,  0, -2,  0,  2,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  1 },
      {  2,  0,  0, -2,  0,  0,  0, -2,  0,  2,  0,  0,  0,  0 },

   /* 31-40 */
      {  0,  0,  0,  0,  0,  0,  0,  2,  0, -2,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  8,-13,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  1,  0,  1,  0,  0,  0,  2 },
      {  0,  0,  0,  0,  0,  0,  5, -8,  0,  0,  0,  0,  0, -2 },
      {  0,  0,  0,  0,  0,  0,  0,  2, -2,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  0,  0,  2, -5,  0,  0,  1 },
      {  2,  0,  0, -2,  0,  0,  0, -2,  0,  3,  0,  0,  0,  0 },
      {  0,  0,  1, -1,  1,  0,  0, -1,  0, -1,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  3, -4,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  1, -1,  1,  0,  0, -1,  0,  0, -1,  0,  0,  0 },

   /* 41-50 */
      {  0,  0,  0,  0,  0,  0,  0,  1,  0, -2,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  5, -7,  0,  0,  0,  0,  0, -2 },
      {  0,  0,  1, -1,  0,  0,  0,  0, -2,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  4,  0, -2,  0,  0,  0,  2 },
      {  0,  0,  0,  0,  0,  0,  8,-13,  0,  0,  0,  0,  0, -2 },
      {  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  2, -1,  0,  0,  0,  0,  0,  2 },
      {  1,  0,  0,  0,  0,  0,-18, 16,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  1, -1,  1,  0,  0, -1,  0,  2,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  2,  0,  1,  0,  0,  0,  2 },

   /* 51-60 */
      {  0,  0,  1, -1,  1,  0, -5,  7,  0,  0,  0,  0,  0,  0 },
      {  1,  0,  0,  0,  0,  0,-10,  3,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  2, -2,  0,  0, -5,  6,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  2,  0, -1,  0,  0,  0,  2 },
      {  1,  0,  2,  0,  2,  0,  0,  1,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  4, -2,  0,  0,  0,  0,  2 },
      {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  1 },
      {  1,  0, -2,  0, -2,  0,  0,  4, -8,  3,  0,  0,  0,  0 },
      {  0,  0,  1, -1,  1,  0,  0, -1,  0,  0,  2,  0,  0,  0 },
      {  0,  0,  2, -2,  1,  0, -3,  3,  0,  0,  0,  0,  0,  0 },

   /* 61-70 */
      {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  2 },
      {  0,  0,  0,  0,  0,  0,  0,  8,-16,  4,  5,  0,  0, -2 },
      {  0,  0,  1, -1,  1,  0,  0,  3, -8,  3,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  8,-11,  0,  0,  0,  0,  0, -2 },
      {  0,  0,  0,  0,  0,  0,  0,  0,  0,  3,  0,  0,  0,  2 },
      {  0,  0,  0,  0,  0,  0,  0,  8,-16,  4,  5,  0,  0,  2 },
      {  0,  0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  0, -1 },
      {  0,  0,  0,  0,  0,  0,  4, -6,  0,  0,  0,  0,  0, -1 },
      {  0,  0,  0,  0,  0,  0,  0,  1,  0, -3,  0,  0,  0, -2 },
      {  0,  0,  0,  0,  0,  0,  0,  2, -4,  0,  0,  0,  0,  0 },

   /* 71-80 */
      {  0,  0,  0,  0,  0,  0,  6, -8,  0,  0,  0,  0,  0, -2 },
      {  0,  0,  0,  0,  0,  0,  3, -2,  0,  0,  0,  0,  0,  2 },
      {  0,  0,  0,  0,  0,  0,  8,-15,  0,  0,  0,  0,  0, -2 },
      {  0,  0,  0,  0,  0,  0,  2, -5,  0,  0,  0,  0,  0, -2 },
      {  0,  0,  0,  0,  0,  0,  1, -3,  0,  0,  0,  0,  0, -2 },
      {  0,  0,  0,  0,  0,  0,  0,  3,  0, -2,  0,  0,  0,  2 },
      {  0,  0,  1, -1,  1,  0,  0, -5,  8, -3,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  1,  2,  0,  0,  0,  0,  2 },
      {  0,  0,  0,  0,  0,  0,  0,  3, -2,  0,  0,  0,  0,  2 },
      {  0,  0,  0,  0,  0,  0,  3, -5,  0,  0,  0,  0,  0,  0 },

   /* 81-90 */
      {  2,  0,  0, -2,  1,  0,  0, -2,  0,  3,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  5, -8,  0,  0,  0,  0,  0, -1 },
      {  2,  0,  0, -2,  0,  0, -3,  3,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  1,  0,  8,-13,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  1,  0,  0,  0,  0, -2,  5,  0,  0,  0 },
      {  1,  0,  0, -1,  0,  0, -3,  4,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  2 },
      {  1,  0,  0,  0, -1,  0,-18, 16,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  1,  0,  0,  0,  0,  2, -5,  0,  0,  0 },
      {  0,  0,  0,  0,  1,  0,  0,  0,  0,  1,  0,  0,  0,  0 },

   /* 91-100 */
      {  1,  0,  0, -2,  0,  0, 19,-21,  3,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  1,  0, -8, 13,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  1, -1,  1,  0,  0, -1,  0,  0,  1,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  7, -9,  0,  0,  0,  0,  0, -2 },
      {  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  2 },
      {  1,  0,  0,  0,  1,  0,-18, 16,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  2, -4,  0,  0,  0,  0,  0, -1 },
      {  0,  0,  0,  0,  0,  0,  0,  6,-16,  4,  5,  0,  0, -2 },
      {  0,  0,  0,  0,  0,  0,  4, -7,  0,  0,  0,  0,  0, -2 },
      {  0,  0,  0,  0,  0,  0,  3, -7,  0,  0,  0,  0,  0, -2 },

   /* 101-110 */
      {  0,  0,  0,  0,  0,  0,  2, -2,  0,  0,  0,  0,  0, -1 },
      {  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  1 },
      {  2,  0,  0, -2,  1,  0,  0, -2,  0,  2,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0, -1 },
      {  0,  0,  0,  0,  0,  0,  0,  3, -4,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  1, -2,  0,  0,  0,  0,  0,  0 },
      {  2,  0,  0, -2, -1,  0,  0, -2,  0,  3,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  3, -3,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  1,  0,  2,  0,  0,  0,  2 },

   /* 111-120 */
      {  0,  0,  0,  0,  1,  0,  0,  1, -2,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  2 },
      {  0,  0,  2, -2,  1,  0,  0, -2,  0,  2,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  2,  0, -3,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  3, -5,  0,  0,  0,  0,  0, -1 },
      {  0,  0,  0,  0,  0,  0,  3, -3,  0,  0,  0,  0,  0,  2 },
      {  0,  0,  0,  0,  0,  0,  4, -4,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  1, -1,  0,  0,  0, -1,  0, -1,  0,  0,  0,  0 },
      {  2,  0,  0, -2,  0,  0, -6,  8,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  1, -1,  1,  0,  0, -2,  2,  0,  0,  0,  0,  0 },

   /* 121-130 */
      {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  1 },
      {  0,  0,  1, -1,  1,  0,  0, -1,  0,  1,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  1, -2,  0,  0,  0,  0, -1 },
      {  0,  0,  0,  0,  0,  0,  0,  2, -3,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  2, -4,  0,  0,  0,  0, -2 },
      {  0,  0,  0,  0,  0,  0,  0,  1,  0,  0, -1,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  8,-10,  0,  0,  0,  0,  0, -2 },
      {  0,  0,  1, -1,  1,  0, -3,  4,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  6, -9,  0,  0,  0,  0,  0, -2 },
      {  1,  0,  0, -1,  1,  0,  0, -1,  0,  2,  0,  0,  0,  0 },

   /* 131-140 */
      {  0,  0,  0,  0,  0,  0,  5, -7,  0,  0,  0,  0,  0, -1 },
      {  0,  0,  0,  0,  0,  0,  5, -5,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0, -1 },
      {  0,  0,  0,  0,  0,  0,  3, -3,  0,  0,  0,  0,  0, -1 },
      {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  0,  4,  0,  0,  0,  0,  2 },
      {  0,  0,  0,  0,  0,  0,  0,  4,  0, -3,  0,  0,  0,  2 },
      {  0,  0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  0,  1 },
      {  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  1 },
      {  0,  0,  0,  0,  1,  0,  2, -3,  0,  0,  0,  0,  0,  0 },

   /* 141-150 */
      {  1,  0,  0, -1,  0,  0,  0, -1,  0,  1,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  1, -3,  0,  0,  0,  0,  0, -1 },
      {  0,  0,  0,  0,  0,  0,  0,  5, -4,  0,  0,  0,  0,  2 },
      {  0,  0,  0,  0,  0,  0,  0,  4, -4,  0,  0,  0,  0,  2 },
      {  0,  0,  0,  0,  0,  0,  9,-11,  0,  0,  0,  0,  0, -2 },
      {  0,  0,  0,  0,  0,  0,  2, -3,  0,  0,  0,  0,  0, -1 },
      {  0,  0,  0,  0,  0,  0,  0,  8,-15,  0,  0,  0,  0,  0 },
      {  0,  0,  1, -1,  1,  0, -4,  5,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  4, -6,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  4,  0, -1,  0,  0,  0,  2 },

   /* 151-160 */
      {  1,  0,  0, -1,  1,  0, -3,  4,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  1,  1,  1,  0,  0,  1,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  1, -1,  1,  0,  0, -1,  0, -4, 10,  0,  0,  0 },
      {  0,  0,  0,  0,  1,  0,  1, -1,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  1, -1,  0,  0,  0, -1,  0,  0, -1,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  1,  0, -3,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  3, -1,  0,  0,  0,  0,  0,  2 },
      {  0,  0,  0,  0,  0,  0,  0,  1,  0, -4,  0,  0,  0, -2 },
      {  0,  0,  0,  0,  0,  0,  0,  0,  0,  2, -5,  0,  0, -2 },
      {  0,  0,  2, -2,  1,  0, -4,  4,  0,  0,  0,  0,  0,  0 },

   /* 161-170 */
      {  0,  0,  0,  0,  0,  0,  0,  3,  0,  0, -1,  0,  0,  2 },
      {  0,  0,  0,  0,  0,  0,  0,  4, -3,  0,  0,  0,  0,  2 },
      {  0,  0,  1, -1,  1,  0,  0, -1,  0,  0,  0,  0,  2,  0 },
      {  0,  0,  0,  0,  0,  0,  4, -4,  0,  0,  0,  0,  0, -1 },
      {  0,  0,  0,  0,  0,  0,  0,  2, -4,  0,  0,  0,  0, -1 },
      {  0,  0,  0,  0,  0,  0,  5, -8,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  1, -2,  0,  0,  0,  0,  1 },
      {  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  1,  0,  0,  0 },
      {  0,  0,  2, -2,  1,  0,  0, -9, 13,  0,  0,  0,  0,  0 },
      {  2,  0,  2,  0,  2,  0,  0,  2,  0, -3,  0,  0,  0,  0 },

   /* 171-180 */
      {  0,  0,  0,  0,  0,  0,  3, -6,  0,  0,  0,  0,  0, -2 },
      {  0,  0,  1, -1,  2,  0,  0, -1,  0,  0,  2,  0,  0,  0 },
      {  1,  0,  0, -1, -1,  0, -3,  4,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  3, -6,  0,  0,  0,  0, -2 },
      {  0,  0,  0,  0,  0,  0,  6, -6,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  0,  0,  3,  0,  0,  0,  1 },
      {  1,  0,  2,  0,  1,  0,  0, -2,  0,  3,  0,  0,  0,  0 },
      {  1,  0, -2,  0, -1,  0,  0, -1,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  1,  0,  0, -2,  4,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  3, -5,  0,  0,  0,  0,  0 },

   /* 181-190 */
      {  0,  0,  0,  0,  0,  0,  2,  1,  0,  0,  0,  0,  0,  2 },
      {  0,  0,  0,  0,  0,  0,  1,  1,  0,  0,  0,  0,  0,  1 },
      {  0,  0,  2,  0,  2,  0,  0,  1,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  1, -8,  3,  0,  0,  0, -2 },
      {  0,  0,  0,  0,  0,  0,  6,-10,  0,  0,  0,  0,  0, -2 },
      {  0,  0,  0,  0,  0,  0,  0,  7, -8,  3,  0,  0,  0,  2 },
      {  0,  0,  0,  0,  1,  0, -3,  5,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  1, -1,  1,  0, -1,  0,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  1, -1,  0,  0, -5,  7,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  2,  0, -2,  0,  0,  0,  1 },

   /* 191-200 */
      {  0,  0,  0,  0,  0,  0,  0,  2,  0, -1,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  7,-10,  0,  0,  0,  0,  0, -2 },
      {  1,  0,  0, -2,  0,  0,  0, -2,  0,  2,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  1,  0,  2, -5,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  6, -8,  0,  0,  0,  0,  0, -1 },
      {  0,  0,  1, -1,  1,  0,  0, -9, 15,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  1,  0, -2,  3,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  1,  0, -1,  1,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  3, -6,  0,  0,  0,  0,  0 },

   /* 201-210 */
      {  0,  0,  0,  0,  0,  0,  0,  1, -4,  0,  0,  0,  0, -2 },
      {  0,  0,  0,  0,  0,  0,  0,  0,  3,  0,  0,  0,  0,  2 },
      {  0,  0,  0,  0,  0,  0,  0,  2,  0,  0, -1,  0,  0,  2 },
      {  2,  0,  0, -2,  1,  0, -6,  8,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  5, -5,  0,  0,  0,  0,  0, -1 },
      {  0,  0,  1, -1,  1,  0,  3, -6,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  1, -1,  1,  0, -2,  2,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  1, -1,  1,  0,  8,-14,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0 },

   /* 211-220 */
      {  0,  0,  0,  0,  1,  0,  0,  8,-15,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  4, -6,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  7, -7,  0,  0,  0,  0,  0,  0 },
      {  2,  0,  0, -2,  1,  0, -3,  3,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  3, -1,  0,  0,  0,  0,  2 },
      {  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  1,  0,  0,  2 },
      {  2,  0, -1, -1,  0,  0,  0,  3, -7,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  4, -7,  0,  0,  0,  0, -2 },
      {  0,  0,  0,  0,  0,  0,  0,  3, -3,  0,  0,  0,  0,  0 },
      {  0,  0,  1, -1,  1,  0,  0, -3,  4,  0,  0,  0,  0,  0 },

   /* 221-230 */
      {  2,  0,  0, -2,  0,  0,  0, -6,  8,  0,  0,  0,  0,  0 },
      {  2,  0,  0, -2,  0,  0,  0, -5,  6,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  1,  0,  0,  0,  0, -1,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  0,  1 },
      {  0,  0,  0,  0,  0,  0,  2,  1,  0,  0,  0,  0,  0,  1 },
      {  0,  0,  0,  0,  0,  0,  1,  2,  0,  0,  0,  0,  0,  2 },
      {  0,  0,  0,  0,  1,  0,  0,  1,  0, -1,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  3, -9,  4,  0,  0,  0,  0, -2 },
      {  0,  0,  0,  0,  0,  0,  0,  3, -5,  0,  0,  0,  0, -2 },

   /* 231-240 */
      {  0,  0,  0,  0,  0,  0,  0,  2,  0, -4,  0,  0,  0, -2 },
      {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  1 },
      {  0,  0,  0,  0,  0,  0,  7,-11,  0,  0,  0,  0,  0, -2 },
      {  0,  0,  0,  0,  0,  0,  3, -5,  4,  0,  0,  0,  0,  2 },
      {  0,  0,  1, -1,  0,  0,  0, -1,  0, -1,  1,  0,  0,  0 },
      {  2,  0,  0,  0,  0,  0,  0, -2,  0,  3,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  8,-15,  0,  0,  0,  0, -2 },
      {  0,  0,  1, -1,  2,  0,  0, -2,  2,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  3,  0,  0,  2 },
      {  0,  0,  0,  0,  0,  0,  6, -6,  0,  0,  0,  0,  0, -1 },

   /* 241-250 */
      {  0,  0,  1, -1,  1,  0,  0, -1,  0, -1,  1,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  2, -2,  0,  0,  0,  0,  0,  1 },
      {  0,  0,  0,  0,  0,  0,  0,  4, -7,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  3, -8,  3,  0,  0,  0,  0 },
      {  0,  0,  1, -1,  1,  0,  2, -4,  0, -3,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  1,  0,  3, -5,  0,  2,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  3,  0, -3,  0,  0,  0,  2 },
      {  0,  0,  2, -2,  2,  0, -8, 11,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  5, -8,  3,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  1,  0,  0, -2,  0,  0,  0 },

   /* 251-260 */
      {  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  1,  0,  0,  2 },
      {  0,  0,  0,  0,  0,  0,  0,  5, -9,  0,  0,  0,  0, -2 },
      {  0,  0,  0,  0,  0,  0,  0,  5, -5,  0,  0,  0,  0,  2 },
      {  0,  0,  0,  0,  0,  0,  7, -9,  0,  0,  0,  0,  0, -1 },
      {  0,  0,  0,  0,  0,  0,  4, -7,  0,  0,  0,  0,  0, -1 },
      {  0,  0,  0,  0,  0,  0,  2, -1,  0,  0,  0,  0,  0,  0 },
      {  1,  0, -2, -2, -2,  0,  0, -2,  0,  2,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  1,  1,  0,  0,  0,  0,  2 },
      {  0,  0,  0,  0,  0,  0,  0,  2,  0, -2,  5,  0,  0,  2 },
      {  0,  0,  0,  0,  0,  0,  3, -3,  0,  0,  0,  0,  0,  1 },

   /* 261-270 */
      {  0,  0,  0,  0,  0,  0,  0,  6,  0,  0,  0,  0,  0,  2 },
      {  0,  0,  0,  0,  0,  0,  0,  2,  0,  2, -5,  0,  0,  2 },
      {  2,  0,  0, -2, -1,  0,  0, -2,  0,  0,  5,  0,  0,  0 },
      {  2,  0,  0, -2, -1,  0, -6,  8,  0,  0,  0,  0,  0,  0 },
      {  1,  0,  0, -2,  0,  0, -3,  3,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  8, -8,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  3,  0,  2, -5,  0,  0,  2 },
      {  0,  0,  0,  0,  1,  0,  3, -7,  4,  0,  0,  0,  0,  0 },
      {  0,  0,  2, -2,  1,  0, -2,  2,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  1,  0,  0, -1,  0,  1,  0,  0,  0,  0 },

   /* 271-280 */
      {  0,  0,  1, -1,  0,  0,  0, -1,  0, -2,  5,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  3,  0, -3,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  3, -1,  0,  0,  0,  0,  0,  1 },
      {  0,  0,  0,  0,  0,  0,  2, -3,  0,  0,  0,  0,  0, -2 },
      {  0,  0,  0,  0,  0,  0,  0, 11,  0,  0,  0,  0,  0,  2 },
      {  0,  0,  0,  0,  0,  0,  0,  6,-15,  0,  0,  0,  0, -2 },
      {  0,  0,  0,  0,  0,  0,  0,  3,  0,  1,  0,  0,  0,  2 },
      {  1,  0,  0, -1,  0,  0,  0, -3,  4,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  1,  0, -3,  7, -4,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  5,  0, -2,  0,  0,  0,  2 },

   /* 281-290 */
      {  0,  0,  0,  0,  0,  0,  3, -5,  0,  0,  0,  0,  0,  1 },
      {  0,  0,  2, -2,  2,  0, -5,  6,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  2, -2,  2,  0, -3,  3,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  3,  0,  0,  0,  0,  0,  0,  2 },
      {  0,  0,  0,  0,  0,  0,  0,  6,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  4, -4,  0,  0,  0,  0,  0,  2 },
      {  0,  0,  0,  0,  0,  0,  0,  4, -8,  0,  0,  0,  0, -2 },
      {  0,  0,  0,  0,  0,  0,  0,  4, -5,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  5, -7,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  6,-11,  0,  0,  0,  0, -2 },

   /* 291-300 */
      {  0,  0,  0,  0,  0,  0,  0,  1, -3,  0,  0,  0,  0, -2 },
      {  0,  0,  1, -1,  1,  0,  0, -1,  0,  3,  0,  0,  0,  0 },
      {  0,  0,  1, -1,  0,  0,  0, -1,  0,  2,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  1, -2,  0,  0,  0,  0,  0,  1 },
      {  0,  0,  0,  0,  0,  0,  9,-12,  0,  0,  0,  0,  0, -2 },
      {  0,  0,  0,  0,  0,  0,  4, -4,  0,  0,  0,  0,  0,  1 },
      {  0,  0,  1, -1,  0,  0, -8, 12,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  1, -1,  1,  0, -2,  3,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  7, -7,  0,  0,  0,  0,  0, -1 },
      {  0,  0,  0,  0,  0,  0,  0,  3, -6,  0,  0,  0,  0, -1 },

   /* 301-310 */
      {  0,  0,  0,  0,  0,  0,  0,  6, -6,  0,  0,  0,  0,  2 },
      {  0,  0,  0,  0,  0,  1,  0, -4,  0,  0,  0,  0,  0, -2 },
      {  0,  0,  1, -1,  1,  0,  0,  1,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  6, -9,  0,  0,  0,  0,  0, -1 },
      {  0,  0,  1, -1, -1,  0,  0,  0, -2,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  1, -5,  0,  0,  0,  0, -2 },
      {  2,  0,  0, -2,  0,  0,  0, -2,  0,  3, -1,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  2,  0,  0, -2,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  5, -9,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  5, -6,  0,  0,  0,  0,  0,  2 },

   /* 311-320 */
      {  0,  0,  0,  0,  0,  0,  9, -9,  0,  0,  0,  0,  0, -1 },
      {  0,  0,  1, -1,  1,  0,  0, -1,  0,  0,  3,  0,  0,  0 },
      {  0,  0,  0,  0,  1,  0,  0,  2, -4,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  5, -3,  0,  0,  0,  0,  0,  2 },
      {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  3,  0,  0,  1 },
      {  0,  0,  1, -1,  2,  0,  0, -1,  0,  2,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  5, -9,  0,  0,  0,  0,  0, -2 },
      {  0,  0,  0,  0,  0,  0,  0,  5, -3,  0,  0,  0,  0,  2 },
      {  0,  0,  0,  0,  0,  0,  0,  0,  0,  4,  0,  0,  0,  2 },
      {  0,  0,  2,  0,  2,  0,  0,  4, -8,  3,  0,  0,  0,  0 },

   /* 321-330 */
      {  0,  0,  2,  0,  2,  0,  0, -4,  8, -3,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  5,  0, -3,  0,  0,  0,  2 },
      {  0,  0,  0,  0,  0,  0,  0,  1,  0,  1,  0,  0,  0,  0 },
      {  2,  0, -1, -1, -1,  0,  0, -1,  0,  3,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  4, -3,  0,  0,  0,  0,  0,  2 },
      {  0,  0,  0,  0,  0,  0,  4, -2,  0,  0,  0,  0,  0,  2 },
      {  0,  0,  0,  0,  0,  0,  5,-10,  0,  0,  0,  0,  0, -2 },
      {  0,  0,  0,  0,  0,  0,  8,-13,  0,  0,  0,  0,  0,  1 },
      {  0,  0,  2, -2,  1, -1,  0,  2,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  1, -1,  1,  0,  0, -1,  0,  0,  0,  2,  0,  0 },

   /* 331-340 */
      {  0,  0,  0,  0,  1,  0,  3, -5,  0,  0,  0,  0,  0,  0 },
      {  1,  0,  0, -2,  0,  0,  0, -2,  0,  3,  0,  0,  0,  0 },
      {  0,  0,  2, -2,  0,  0, -3,  3,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  9, -9,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  2,  0,  2,  0,  1, -1,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  2, -2,  1,  0,  0, -8, 11,  0,  0,  0,  0,  0 },
      {  0,  0,  2, -2,  1,  0,  0, -2,  0,  0,  2,  0,  0,  0 },
      {  0,  0,  1, -1,  1,  0,  0, -1,  0, -1,  2,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  5, -5,  0,  0,  0,  0,  0,  2 },
      {  0,  0,  0,  0,  0,  0,  2, -6,  0,  0,  0,  0,  0, -2 },

   /* 341-350 */
      {  0,  0,  0,  0,  0,  0,  0,  8,-15,  0,  0,  0,  0, -1 },
      {  0,  0,  0,  0,  0,  0,  0,  5, -2,  0,  0,  0,  0,  2 },
      {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  2 },
      {  0,  0,  0,  0,  0,  0,  0,  7,-13,  0,  0,  0,  0, -2 },
      {  0,  0,  0,  0,  0,  0,  0,  3,  0, -2,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  1,  0,  3,  0,  0,  0,  2 },
      {  0,  0,  2, -2,  1,  0,  0, -2,  0,  3,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  8, -8,  0,  0,  0,  0,  0, -1 },
      {  0,  0,  0,  0,  0,  0,  8,-10,  0,  0,  0,  0,  0, -1 },
      {  0,  0,  0,  0,  0,  0,  4, -2,  0,  0,  0,  0,  0,  1 },

   /* 351-360 */
      {  0,  0,  0,  0,  0,  0,  3, -6,  0,  0,  0,  0,  0, -1 },
      {  0,  0,  0,  0,  0,  0,  3, -4,  0,  0,  0,  0,  0, -1 },
      {  0,  0,  0,  0,  0,  0,  0,  0,  0,  2, -5,  0,  0,  2 },
      {  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  2 },
      {  0,  0,  0,  0,  0,  0,  0,  2,  0, -4,  0,  0,  0,  0 },
      {  2,  0,  0, -2, -1,  0,  0, -5,  6,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  2, -5,  0,  0,  0,  0, -2 },
      {  2,  0, -1, -1, -1,  0,  0,  3, -7,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  5, -8,  0,  0,  0,  0,  0 },
      {  0,  0,  2,  0,  2,  0, -1,  1,  0,  0,  0,  0,  0,  0 },

   /* 361-370 */
      {  2,  0,  0, -2,  0,  0,  0, -2,  0,  4, -3,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  6,-11,  0,  0,  0,  0,  0 },
      {  2,  0,  0, -2,  1,  0,  0, -6,  8,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  4, -8,  1,  5,  0,  0, -2 },
      {  0,  0,  0,  0,  0,  0,  0,  6, -5,  0,  0,  0,  0,  2 },
      {  1,  0, -2, -2, -2,  0, -3,  3,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  1, -1,  2,  0,  0,  0, -2,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  2,  0,  0,  4, -8,  3,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  2,  0,  0, -4,  8, -3,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  6,  0,  0,  0,  0,  0,  1 },

   /* 371-380 */
      {  0,  0,  0,  0,  0,  0,  0,  6, -7,  0,  0,  0,  0,  2 },
      {  0,  0,  0,  0,  0,  0,  0,  4,  0,  0, -2,  0,  0,  2 },
      {  0,  0,  0,  0,  0,  0,  0,  3,  0,  0, -2,  0,  0,  2 },
      {  0,  0,  0,  0,  0,  0,  0,  1,  0, -1,  0,  0,  0,  1 },
      {  0,  0,  0,  0,  0,  0,  0,  1, -6,  0,  0,  0,  0, -2 },
      {  0,  0,  0,  0,  0,  0,  0,  0,  0,  4, -5,  0,  0,  2 },
      {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  2 },
      {  0,  0,  0,  0,  0,  0,  3, -5,  0,  2,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  7,-13,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  2,  0, -2,  0,  0,  0,  2 },

   /* 381-390 */
      {  0,  0,  1, -1,  0,  0,  0, -1,  0,  0,  2,  0,  0,  0 },
      {  0,  0,  0,  0,  1,  0,  0, -8, 15,  0,  0,  0,  0,  0 },
      {  2,  0,  0, -2, -2,  0, -3,  3,  0,  0,  0,  0,  0,  0 },
      {  2,  0, -1, -1, -1,  0,  0, -1,  0,  2,  0,  0,  0,  0 },
      {  1,  0,  2, -2,  2,  0,  0, -2,  0,  2,  0,  0,  0,  0 },
      {  1,  0, -1,  1, -1,  0,-18, 17,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  2,  0,  2,  0,  0,  1,  0, -1,  0,  0,  0,  0 },
      {  0,  0,  2,  0,  2,  0,  0, -1,  0,  1,  0,  0,  0,  0 },
      {  0,  0,  2, -2, -1,  0, -5,  6,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  1, -1,  2,  0,  0, -1,  0,  1,  0,  0,  0,  0 },

   /* 391-400 */
      {  0,  0,  0,  0,  1,  0,  2, -2,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  8,-16,  0,  0,  0,  0,  0, -2 },
      {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  5,  0,  0,  2 },
      {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  2 },
      {  0,  0,  0,  0,  2,  0,  0, -1,  2,  0,  0,  0,  0,  0 },
      {  2,  0, -1, -1, -2,  0,  0, -1,  0,  2,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  6,-10,  0,  0,  0,  0,  0, -1 },
      {  0,  0,  1, -1,  1,  0,  0, -1,  0, -2,  4,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  2,  2,  0,  0,  0,  0,  2 },
      {  2,  0,  0, -2, -1,  0,  0, -2,  0,  4, -5,  0,  0,  0 },

   /* 401-410 */
      {  2,  0,  0, -2, -1,  0, -3,  3,  0,  0,  0,  0,  0,  0 },
      {  2,  0, -1, -1, -1,  0,  0, -1,  0,  0,  0,  0,  0,  0 },
      {  1,  0,  1, -1,  1,  0,  0, -1,  0,  0,  0,  0,  0,  0 },
      {  1,  0,  0, -1, -1,  0,  0, -2,  2,  0,  0,  0,  0,  0 },
      {  1,  0, -1, -1, -1,  0, 20,-20,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  2, -2,  1,  0,  0, -1,  0,  1,  0,  0,  0,  0 },
      {  0,  0,  1, -1,  1,  0,  1, -2,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  1, -1,  1,  0, -2,  1,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  1,  0,  5, -8,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  1,  0,  0,  0,  0,  0, -1,  0,  0,  0 },

   /* 411-420 */
      {  0,  0,  0,  0,  0,  0,  9,-11,  0,  0,  0,  0,  0, -1 },
      {  0,  0,  0,  0,  0,  0,  5, -3,  0,  0,  0,  0,  0,  1 },
      {  0,  0,  0,  0,  0,  0,  0,  1,  0, -3,  0,  0,  0, -1 },
      {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  1 },
      {  0,  0,  0,  0,  0,  0,  6, -7,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  3, -2,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  1, -2,  0,  0,  0,  0,  0, -2 },
      {  0,  0,  1, -1,  1,  0,  0, -1,  0,  0, -2,  0,  0,  0 },
      {  0,  0,  1, -1,  2,  0,  0, -1,  0, -2,  5,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  5, -7,  0,  0,  0,  0,  0 },

   /* 421-430 */
      {  0,  0,  0,  0,  0,  0,  1, -3,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  5, -8,  0,  0,  0,  0, -2 },
      {  0,  0,  0,  0,  0,  0,  0,  2, -6,  0,  0,  0,  0, -2 },
      {  1,  0,  0, -2,  0,  0, 20,-21,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  8,-12,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  5, -6,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  4, -4,  0,  0,  0,  0,  0 },
      {  0,  0,  1, -1,  2,  0,  0, -1,  0, -1,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  8,-12,  0,  0,  0,  0,  0, -2 },
      {  0,  0,  0,  0,  0,  0,  0,  9,-17,  0,  0,  0,  0,  0 },

   /* 431-440 */
      {  0,  0,  0,  0,  0,  0,  0,  5, -6,  0,  0,  0,  0,  2 },
      {  0,  0,  0,  0,  0,  0,  0,  4, -8,  1,  5,  0,  0,  2 },
      {  0,  0,  0,  0,  0,  0,  0,  4, -6,  0,  0,  0,  0, -2 },
      {  0,  0,  0,  0,  0,  0,  0,  2, -7,  0,  0,  0,  0, -2 },
      {  1,  0,  0, -1,  1,  0,  0, -3,  4,  0,  0,  0,  0,  0 },
      {  1,  0, -2,  0, -2,  0,-10,  3,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  1,  0,  0, -9, 17,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  1, -4,  0,  0,  0,  0,  0, -2 },
      {  1,  0, -2, -2, -2,  0,  0, -2,  0,  3,  0,  0,  0,  0 },
      {  1,  0, -1,  1, -1,  0,  0,  1,  0,  0,  0,  0,  0,  0 },

   /* 441-450 */
      {  0,  0,  2, -2,  2,  0,  0, -2,  0,  2,  0,  0,  0,  0 },
      {  0,  0,  1, -1,  2,  0,  0, -1,  0,  0,  1,  0,  0,  0 },
      {  0,  0,  1, -1,  2,  0, -5,  7,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  1,  0,  0,  2, -2,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  4, -5,  0,  0,  0,  0,  0, -1 },
      {  0,  0,  0,  0,  0,  0,  3, -4,  0,  0,  0,  0,  0, -2 },
      {  0,  0,  0,  0,  0,  0,  2, -4,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  5,-10,  0,  0,  0,  0, -2 },
      {  0,  0,  0,  0,  0,  0,  0,  4,  0, -4,  0,  0,  0,  2 },
      {  0,  0,  0,  0,  0,  0,  0,  2,  0, -5,  0,  0,  0, -2 },

   /* 451-460 */
      {  0,  0,  0,  0,  0,  0,  0,  1,  0, -5,  0,  0,  0, -2 },
      {  0,  0,  0,  0,  0,  0,  0,  1,  0, -2,  5,  0,  0,  2 },
      {  0,  0,  0,  0,  0,  0,  0,  1,  0, -2,  0,  0,  0, -2 },
      {  0,  0,  0,  0,  0,  0,  2, -3,  0,  0,  0,  0,  0,  1 },
      {  1,  0,  0, -2,  0,  0,  0,  1,  0, -1,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  3, -7,  4,  0,  0,  0,  0,  0 },
      {  2,  0,  2,  0,  1,  0,  0,  1,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  1, -1, -1,  0,  0, -1,  0, -1,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  1,  0,  0,  1,  0, -2,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  6,-10,  0,  0,  0,  0, -2 },

   /* 461-470 */
      {  1,  0,  0, -1,  1,  0,  0, -1,  0,  1,  0,  0,  0,  0 },
      {  0,  0,  2, -2,  1,  0,  0,  4, -8,  3,  0,  0,  0,  0 },
      {  0,  0,  2, -2,  1,  0,  0,  1,  0, -1,  0,  0,  0,  0 },
      {  0,  0,  2, -2,  1,  0,  0, -4,  8, -3,  0,  0,  0,  0 },
      {  0,  0,  2, -2,  1,  0,  0, -3,  0,  3,  0,  0,  0,  0 },
      {  0,  0,  2, -2,  1,  0, -5,  5,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  1, -1,  1,  0,  1, -3,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  1, -1,  1,  0,  0, -4,  6,  0,  0,  0,  0,  0 },
      {  0,  0,  1, -1,  1,  0,  0, -1,  0,  0,  0, -1,  0,  0 },
      {  0,  0,  1, -1,  1,  0, -5,  6,  0,  0,  0,  0,  0,  0 },

   /* 471-480 */
      {  0,  0,  0,  0,  1,  0,  3, -4,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  1,  0, -2,  2,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  7,-10,  0,  0,  0,  0,  0, -1 },
      {  0,  0,  0,  0,  0,  0,  5, -5,  0,  0,  0,  0,  0,  1 },
      {  0,  0,  0,  0,  0,  0,  4, -5,  0,  0,  0,  0,  0, -2 },
      {  0,  0,  0,  0,  0,  0,  3, -8,  0,  0,  0,  0,  0, -2 },
      {  0,  0,  0,  0,  0,  0,  2, -5,  0,  0,  0,  0,  0, -1 },
      {  0,  0,  0,  0,  0,  0,  1, -2,  0,  0,  0,  0,  0, -1 },
      {  0,  0,  0,  0,  0,  0,  0,  7, -9,  0,  0,  0,  0,  2 },
      {  0,  0,  0,  0,  0,  0,  0,  7, -8,  0,  0,  0,  0,  2 },

   /* 481-490 */
      {  0,  0,  0,  0,  0,  0,  0,  3,  0,  0,  0,  0,  0,  2 },
      {  0,  0,  0,  0,  0,  0,  0,  3, -8,  3,  0,  0,  0, -2 },
      {  0,  0,  0,  0,  0,  0,  0,  2,  0,  0, -2,  0,  0,  1 },
      {  0,  0,  0,  0,  0,  0,  0,  2, -4,  0,  0,  0,  0,  1 },
      {  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0, -1 },
      {  0,  0,  0,  0,  0,  0,  0,  1,  0, -1,  0,  0,  0, -1 },
      {  2,  0,  0, -2, -1,  0,  0, -6,  8,  0,  0,  0,  0,  0 },
      {  2,  0, -1, -1,  1,  0,  0,  3, -7,  0,  0,  0,  0,  0 },
      {  0,  0,  2, -2,  1,  0,  0, -7,  9,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  3, -5,  0,  0,  0,  0, -1 },

   /* 491-500 */
      {  0,  0,  1, -1,  2,  0, -8, 12,  0,  0,  0,  0,  0,  0 },
      {  1,  0,  0,  0,  0,  0,  0, -2,  0,  2,  0,  0,  0,  0 },
      {  1,  0,  0, -2,  0,  0,  2, -2,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  7, -8,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0 },
      {  2,  0,  0, -2,  1,  0,  0, -5,  6,  0,  0,  0,  0,  0 },
      {  2,  0,  0, -2, -1,  0,  0, -2,  0,  3, -1,  0,  0,  0 },
      {  1,  0,  1,  1,  1,  0,  0,  1,  0,  0,  0,  0,  0,  0 },
      {  1,  0,  0, -2,  1,  0,  0, -2,  0,  2,  0,  0,  0,  0 },
      {  1,  0,  0, -2, -1,  0,  0, -2,  0,  2,  0,  0,  0,  0 },

   /* 501-510 */
      {  1,  0,  0, -1, -1,  0,  0, -3,  4,  0,  0,  0,  0,  0 },
      {  1,  0, -1,  0, -1,  0, -3,  5,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  2, -2,  1,  0,  0, -4,  4,  0,  0,  0,  0,  0 },
      {  0,  0,  2, -2,  1,  0,  0, -2,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  2, -2,  1,  0, -8, 11,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  2, -2,  0,  0,  0, -9, 13,  0,  0,  0,  0,  0 },
      {  0,  0,  1,  1,  2,  0,  0,  1,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  1, -1,  1,  0,  0,  1, -4,  0,  0,  0,  0,  0 },
      {  0,  0,  1, -1,  1,  0,  0, -1,  0,  1, -3,  0,  0,  0 },
      {  0,  0,  0,  0,  1,  0,  0,  7,-13,  0,  0,  0,  0,  0 },

   /* 511-520 */
      {  0,  0,  0,  0,  1,  0,  0,  2,  0, -2,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  1,  0,  0, -2,  2,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  1,  0, -3,  4,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  1,  0, -4,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  7,-11,  0,  0,  0,  0,  0, -1 },
      {  0,  0,  0,  0,  0,  0,  6, -6,  0,  0,  0,  0,  0,  1 },
      {  0,  0,  0,  0,  0,  0,  6, -4,  0,  0,  0,  0,  0,  1 },
      {  0,  0,  0,  0,  0,  0,  5, -6,  0,  0,  0,  0,  0, -1 },
      {  0,  0,  0,  0,  0,  0,  4, -2,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  3, -4,  0,  0,  0,  0,  0,  1 },

   /* 521-530 */
      {  0,  0,  0,  0,  0,  0,  1, -4,  0,  0,  0,  0,  0, -1 },
      {  0,  0,  0,  0,  0,  0,  0,  9,-17,  0,  0,  0,  0, -2 },
      {  0,  0,  0,  0,  0,  0,  0,  7, -7,  0,  0,  0,  0,  2 },
      {  0,  0,  0,  0,  0,  0,  0,  4, -8,  3,  0,  0,  0,  1 },
      {  0,  0,  0,  0,  0,  0,  0,  4, -8,  3,  0,  0,  0, -1 },
      {  0,  0,  0,  0,  0,  0,  0,  4, -8,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  4, -7,  0,  0,  0,  0, -1 },
      {  0,  0,  0,  0,  0,  0,  0,  1,  0,  1,  0,  0,  0,  1 },
      {  0,  0,  0,  0,  0,  0,  0,  1,  0, -4,  0,  0,  0,  0 },
      {  2,  0,  0, -2,  0,  0,  0, -4,  8, -3,  0,  0,  0,  0 },

   /* 531-540 */
      {  2,  0,  0, -2,  0,  0, -2,  2,  0,  0,  0,  0,  0,  0 },
      {  1,  0,  0,  0,  0,  0,  0,  4, -8,  3,  0,  0,  0,  0 },
      {  1,  0,  0,  0,  0,  0,  0, -4,  8, -3,  0,  0,  0,  0 },
      {  1,  0,  0,  0,  0,  0, -1,  1,  0,  0,  0,  0,  0,  0 },
      {  1,  0,  0, -2,  0,  0, 17,-16,  0, -2,  0,  0,  0,  0 },
      {  1,  0,  0, -1,  0,  0,  0, -2,  2,  0,  0,  0,  0,  0 },
      {  0,  0,  2, -2,  0,  0,  0, -2,  0,  2,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  6, -9,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  4,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  3,  0, -4,  0,  0,  0,  0 },

   /* 541-550 */
      {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1, -2, -2 },
      {  0,  0,  0,  0,  0,  0,  0,  2,  1,  0,  0,  0,  0,  2 },
      {  2,  0,  0, -2,  0,  0,  0, -4,  4,  0,  0,  0,  0,  0 },
      {  2,  0,  0, -2,  0,  0,  0, -2,  0,  2,  2,  0,  0,  0 },
      {  1,  0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  0,  0 },
      {  1,  0,  0,  0,  0,  0,  0, -1,  0,  1,  0,  0,  0,  0 },
      {  1,  0,  0,  0,  0,  0, -3,  3,  0,  0,  0,  0,  0,  0 },
      {  1,  0,  0, -2,  0,  0,  1, -1,  0,  0,  0,  0,  0,  0 },
      {  1,  0,  0, -2,  0,  0,  0,  4, -8,  3,  0,  0,  0,  0 },
      {  1,  0,  0, -2,  0,  0,  0, -4,  8, -3,  0,  0,  0,  0 },

   /* 551-560 */
      {  1,  0,  0, -2,  0,  0, -2,  2,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  2, -2,  0,  0, -4,  4,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  1,  1,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  1, -1,  0,  0,  3, -6,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  1, -1,  0,  0,  0, -2,  2,  0,  0,  0,  0,  0 },
      {  0,  0,  1, -1,  0,  0,  0, -1,  0,  1,  0,  0,  0,  0 },
      {  0,  0,  1, -1,  0,  0,  0, -1,  0,  0,  1,  0,  0,  0 },
      {  0,  0,  1, -1,  0,  0, -4,  5,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  1, -1,  0,  0, -3,  4,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  2,  0,  0,  0, -1,  0,  1,  0,  0,  0,  0 },

   /* 561-570 */
      {  0,  0,  0,  0,  0,  0,  8, -9,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  3, -6,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  1,  1,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  0,  0,  3, -5,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  0,  0,  2, -2,  0,  0,  0 },
      {  2,  0, -2, -2, -2,  0,  0, -2,  0,  2,  0,  0,  0,  0 },
      {  1,  0,  0,  0,  1,  0,-10,  3,  0,  0,  0,  0,  0,  0 },
      {  1,  0,  0,  0, -1,  0,-10,  3,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  2,  0,  2,  0,  2, -3,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  2,  0,  2,  0,  2, -2,  0,  0,  0,  0,  0,  0 },

   /* 571-580 */
      {  0,  0,  2,  0,  2,  0, -2,  3,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  2,  0,  2,  0, -2,  2,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  2,  0,  0,  0,  0,  1,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  1,  0,  0, -1,  0,  2,  0,  0,  0,  0 },
      {  2,  0,  2, -2,  2,  0,  0, -2,  0,  3,  0,  0,  0,  0 },
      {  2,  0,  1, -3,  1,  0, -6,  7,  0,  0,  0,  0,  0,  0 },
      {  2,  0,  0, -2,  0,  0,  2, -5,  0,  0,  0,  0,  0,  0 },
      {  2,  0,  0, -2,  0,  0,  0, -2,  0,  5, -5,  0,  0,  0 },
      {  2,  0,  0, -2,  0,  0,  0, -2,  0,  1,  5,  0,  0,  0 },
      {  2,  0,  0, -2,  0,  0,  0, -2,  0,  0,  5,  0,  0,  0 },

   /* 581-590 */
      {  2,  0,  0, -2,  0,  0,  0, -2,  0,  0,  2,  0,  0,  0 },
      {  2,  0,  0, -2,  0,  0, -4,  4,  0,  0,  0,  0,  0,  0 },
      {  2,  0, -2,  0, -2,  0,  0,  5, -9,  0,  0,  0,  0,  0 },
      {  2,  0, -1, -1,  0,  0,  0, -1,  0,  3,  0,  0,  0,  0 },
      {  1,  0,  2,  0,  2,  0,  1, -1,  0,  0,  0,  0,  0,  0 },
      {  1,  0,  2,  0,  2,  0,  0,  4, -8,  3,  0,  0,  0,  0 },
      {  1,  0,  2,  0,  2,  0,  0, -4,  8, -3,  0,  0,  0,  0 },
      {  1,  0,  2,  0,  2,  0, -1,  1,  0,  0,  0,  0,  0,  0 },
      {  1,  0,  2, -2,  2,  0, -3,  3,  0,  0,  0,  0,  0,  0 },
      {  1,  0,  0,  0,  0,  0,  0,  1,  0, -1,  0,  0,  0,  0 },

   /* 591-600 */
      {  1,  0,  0,  0,  0,  0,  0, -2,  0,  3,  0,  0,  0,  0 },
      {  1,  0,  0, -2,  0,  0,  0,  2,  0, -2,  0,  0,  0,  0 },
      {  1,  0, -2, -2, -2,  0,  0,  1,  0, -1,  0,  0,  0,  0 },
      {  1,  0, -1,  1,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0 },
      {  1,  0, -1, -1,  0,  0,  0,  8,-15,  0,  0,  0,  0,  0 },
      {  0,  0,  2,  2,  2,  0,  0,  2,  0, -2,  0,  0,  0,  0 },
      {  0,  0,  2, -2,  1,  0,  1, -1,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  2, -2,  1,  0,  0, -2,  0,  1,  0,  0,  0,  0 },
      {  0,  0,  2, -2,  1,  0,  0,-10, 15,  0,  0,  0,  0,  0 },
      {  0,  0,  2, -2,  0, -1,  0,  2,  0,  0,  0,  0,  0,  0 },

   /* 601-610 */
      {  0,  0,  1, -1,  2,  0,  0, -1,  0,  0, -1,  0,  0,  0 },
      {  0,  0,  1, -1,  2,  0, -3,  4,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  1, -1,  1,  0, -4,  6,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  1, -1,  1,  0, -1,  2,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  1, -1,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  1, -1,  0,  0,  0, -1,  0,  0, -2,  0,  0,  0 },
      {  0,  0,  1, -1,  0,  0, -2,  2,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  1, -1,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  1, -1, -1,  0, -5,  7,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  2,  0,  0,  0,  2,  0, -2,  0,  0,  0,  0 },

   /* 611-620 */
      {  0,  0,  0,  2,  0,  0, -2,  2,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  2,  0, -3,  5,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  1,  0, -1,  2,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  9,-13,  0,  0,  0,  0,  0, -2 },
      {  0,  0,  0,  0,  0,  0,  8,-14,  0,  0,  0,  0,  0, -2 },
      {  0,  0,  0,  0,  0,  0,  8,-11,  0,  0,  0,  0,  0, -1 },
      {  0,  0,  0,  0,  0,  0,  6, -9,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  6, -8,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  6, -7,  0,  0,  0,  0,  0, -1 },
      {  0,  0,  0,  0,  0,  0,  5, -6,  0,  0,  0,  0,  0, -2 },

   /* 621-630 */
      {  0,  0,  0,  0,  0,  0,  5, -6, -4,  0,  0,  0,  0, -2 },
      {  0,  0,  0,  0,  0,  0,  5, -4,  0,  0,  0,  0,  0,  2 },
      {  0,  0,  0,  0,  0,  0,  4, -8,  0,  0,  0,  0,  0, -2 },
      {  0,  0,  0,  0,  0,  0,  4, -5,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  3, -3,  0,  2,  0,  0,  0,  2 },
      {  0,  0,  0,  0,  0,  0,  3, -1,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  0, -2 },
      {  0,  0,  0,  0,  0,  0,  0,  7,-12,  0,  0,  0,  0, -2 },
      {  0,  0,  0,  0,  0,  0,  0,  6, -9,  0,  0,  0,  0, -2 },

   /* 631-640 */
      {  0,  0,  0,  0,  0,  0,  0,  6, -8,  1,  5,  0,  0,  2 },
      {  0,  0,  0,  0,  0,  0,  0,  6, -4,  0,  0,  0,  0,  2 },
      {  0,  0,  0,  0,  0,  0,  0,  6,-10,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  5,  0, -4,  0,  0,  0,  2 },
      {  0,  0,  0,  0,  0,  0,  0,  5, -9,  0,  0,  0,  0, -1 },
      {  0,  0,  0,  0,  0,  0,  0,  5, -8,  3,  0,  0,  0,  2 },
      {  0,  0,  0,  0,  0,  0,  0,  5, -7,  0,  0,  0,  0, -2 },
      {  0,  0,  0,  0,  0,  0,  0,  5, -6,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  5,-16,  4,  5,  0,  0, -2 },
      {  0,  0,  0,  0,  0,  0,  0,  5,-13,  0,  0,  0,  0, -2 },

   /* 641-650 */
      {  0,  0,  0,  0,  0,  0,  0,  3,  0, -5,  0,  0,  0, -2 },
      {  0,  0,  0,  0,  0,  0,  0,  3, -9,  0,  0,  0,  0, -2 },
      {  0,  0,  0,  0,  0,  0,  0,  3, -7,  0,  0,  0,  0, -2 },
      {  0,  0,  0,  0,  0,  0,  0,  2,  0,  2,  0,  0,  0,  2 },
      {  0,  0,  0,  0,  0,  0,  0,  2,  0,  0, -3,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  2, -8,  1,  5,  0,  0, -2 },
      {  0,  0,  0,  0,  0,  0,  0,  1,  0,  1, -5,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  2,  0,  0,  2 },
      {  0,  0,  0,  0,  0,  0,  0,  1,  0,  0, -3,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  1,  0, -3,  5,  0,  0,  0 },

   /* 651-NFPL */
      {  0,  0,  0,  0,  0,  0,  0,  1, -3,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  0,  0,  2, -6,  3,  0, -2 },
      {  0,  0,  0,  0,  0,  0,  0,  0,  0,  1, -2,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2 },
      {  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0 }
   };

/* Number of frequencies:  planetary */
   static const int NFPL = (int) (sizeof mfapl / sizeof (int) / 14);

/* Pointers into amplitudes array, one pointer per frequency */
   static const int nc[] = {

   /* 1-100 */
       1,    21,    37,    51,    65,    79,    91,   103,   115,   127,
     139,   151,   163,   172,   184,   196,   207,   219,   231,   240,
     252,   261,   273,   285,   297,   309,   318,   327,   339,   351,
     363,   372,   384,   396,   405,   415,   423,   435,   444,   452,
     460,   467,   474,   482,   490,   498,   506,   513,   521,   528,
     536,   543,   551,   559,   566,   574,   582,   590,   597,   605,
     613,   620,   628,   636,   644,   651,   658,   666,   674,   680,
     687,   695,   702,   710,   717,   725,   732,   739,   746,   753,
     760,   767,   774,   782,   790,   798,   805,   812,   819,   826,
     833,   840,   846,   853,   860,   867,   874,   881,   888,   895,

   /* 101-200 */
     901,   908,   914,   921,   928,   934,   941,   948,   955,   962,
     969,   976,   982,   989,   996,  1003,  1010,  1017,  1024,  1031,
    1037,  1043,  1050,  1057,  1064,  1071,  1078,  1084,  1091,  1098,
    1104,  1112,  1118,  1124,  1131,  1138,  1145,  1151,  1157,  1164,
    1171,  1178,  1185,  1192,  1199,  1205,  1212,  1218,  1226,  1232,
    1239,  1245,  1252,  1259,  1266,  1272,  1278,  1284,  1292,  1298,
    1304,  1310,  1316,  1323,  1329,  1335,  1341,  1347,  1353,  1359,
    1365,  1371,  1377,  1383,  1389,  1396,  1402,  1408,  1414,  1420,
    1426,  1434,  1440,  1446,  1452,  1459,  1465,  1471,  1477,  1482,
    1488,  1493,  1499,  1504,  1509,  1514,  1520,  1527,  1532,  1538,

   /* 201-300 */
    1543,  1548,  1553,  1558,  1564,  1569,  1574,  1579,  1584,  1589,
    1594,  1596,  1598,  1600,  1602,  1605,  1608,  1610,  1612,  1617,
    1619,  1623,  1625,  1627,  1629,  1632,  1634,  1640,  1642,  1644,
    1646,  1648,  1650,  1652,  1654,  1658,  1660,  1662,  1664,  1668,
    1670,  1672,  1673,  1675,  1679,  1681,  1683,  1684,  1686,  1688,
    1690,  1693,  1695,  1697,  1701,  1703,  1705,  1707,  1709,  1711,
    1712,  1715,  1717,  1721,  1723,  1725,  1727,  1729,  1731,  1733,
    1735,  1737,  1739,  1741,  1743,  1745,  1747,  1749,  1751,  1753,
    1755,  1757,  1759,  1761,  1762,  1764,  1766,  1768,  1769,  1771,
    1773,  1775,  1777,  1779,  1781,  1783,  1785,  1787,  1788,  1790,

   /* 301-400 */
    1792,  1794,  1796,  1798,  1800,  1802,  1804,  1806,  1807,  1809,
    1811,  1815,  1817,  1819,  1821,  1823,  1825,  1827,  1829,  1831,
    1833,  1835,  1837,  1839,  1840,  1842,  1844,  1848,  1850,  1852,
    1854,  1856,  1858,  1859,  1860,  1862,  1864,  1866,  1868,  1869,
    1871,  1873,  1875,  1877,  1879,  1881,  1883,  1885,  1887,  1889,
    1891,  1892,  1896,  1898,  1900,  1901,  1903,  1905,  1907,  1909,
    1910,  1911,  1913,  1915,  1919,  1921,  1923,  1927,  1929,  1931,
    1933,  1935,  1937,  1939,  1943,  1945,  1947,  1948,  1949,  1951,
    1953,  1955,  1957,  1958,  1960,  1962,  1964,  1966,  1968,  1970,
    1971,  1973,  1974,  1975,  1977,  1979,  1980,  1981,  1982,  1984,

   /* 401-500 */
    1986,  1988,  1990,  1992,  1994,  1995,  1997,  1999,  2001,  2003,
    2005,  2007,  2008,  2009,  2011,  2013,  2015,  2017,  2019,  2021,
    2023,  2024,  2025,  2027,  2029,  2031,  2033,  2035,  2037,  2041,
    2043,  2045,  2046,  2047,  2049,  2051,  2053,  2055,  2056,  2057,
    2059,  2061,  2063,  2065,  2067,  2069,  2070,  2071,  2072,  2074,
    2076,  2078,  2080,  2082,  2084,  2086,  2088,  2090,  2092,  2094,
    2095,  2096,  2097,  2099,  2101,  2105,  2106,  2107,  2108,  2109,
    2110,  2111,  2113,  2115,  2119,  2121,  2123,  2125,  2127,  2129,
    2131,  2133,  2135,  2136,  2137,  2139,  2141,  2143,  2145,  2147,
    2149,  2151,  2153,  2155,  2157,  2159,  2161,  2163,  2165,  2167,

   /* 501-600 */
    2169,  2171,  2173,  2175,  2177,  2179,  2181,  2183,  2185,  2186,
    2187,  2188,  2192,  2193,  2195,  2197,  2199,  2201,  2203,  2205,
    2207,  2209,  2211,  2213,  2217,  2219,  2221,  2223,  2225,  2227,
    2229,  2231,  2233,  2234,  2235,  2236,  2237,  2238,  2239,  2240,
    2241,  2244,  2246,  2248,  2250,  2252,  2254,  2256,  2258,  2260,
    2262,  2264,  2266,  2268,  2270,  2272,  2274,  2276,  2278,  2280,
    2282,  2284,  2286,  2288,  2290,  2292,  2294,  2296,  2298,  2300,
    2302,  2303,  2304,  2305,  2306,  2307,  2309,  2311,  2313,  2315,
    2317,  2319,  2321,  2323,  2325,  2327,  2329,  2331,  2333,  2335,
    2337,  2341,  2343,  2345,  2347,  2349,  2351,  2352,  2355,  2356,

   /* 601-700 */
    2357,  2358,  2359,  2361,  2363,  2364,  2365,  2366,  2367,  2368,
    2369,  2370,  2371,  2372,  2373,  2374,  2376,  2378,  2380,  2382,
    2384,  2385,  2386,  2387,  2388,  2389,  2390,  2391,  2392,  2393,
    2394,  2395,  2396,  2397,  2398,  2399,  2400,  2401,  2402,  2403,
    2404,  2405,  2406,  2407,  2408,  2409,  2410,  2411,  2412,  2413,
    2414,  2415,  2417,  2418,  2430,  2438,  2445,  2453,  2460,  2468,
    2474,  2480,  2488,  2496,  2504,  2512,  2520,  2527,  2535,  2543,
    2550,  2558,  2566,  2574,  2580,  2588,  2596,  2604,  2612,  2619,
    2627,  2634,  2642,  2648,  2656,  2664,  2671,  2679,  2685,  2693,
    2701,  2709,  2717,  2725,  2733,  2739,  2747,  2753,  2761,  2769,

   /* 701-800 */
    2777,  2785,  2793,  2801,  2809,  2817,  2825,  2833,  2841,  2848,
    2856,  2864,  2872,  2878,  2884,  2892,  2898,  2906,  2914,  2922,
    2930,  2938,  2944,  2952,  2958,  2966,  2974,  2982,  2988,  2996,
    3001,  3009,  3017,  3025,  3032,  3039,  3045,  3052,  3059,  3067,
    3069,  3076,  3083,  3090,  3098,  3105,  3109,  3111,  3113,  3120,
    3124,  3128,  3132,  3136,  3140,  3144,  3146,  3150,  3158,  3161,
    3165,  3166,  3168,  3172,  3176,  3180,  3182,  3185,  3189,  3193,
    3194,  3197,  3200,  3204,  3208,  3212,  3216,  3219,  3221,  3222,
    3226,  3230,  3234,  3238,  3242,  3243,  3247,  3251,  3254,  3258,
    3262,  3266,  3270,  3274,  3275,  3279,  3283,  3287,  3289,  3293,

   /* 801-900 */
    3296,  3300,  3303,  3307,  3311,  3315,  3319,  3321,  3324,  3327,
    3330,  3334,  3338,  3340,  3342,  3346,  3350,  3354,  3358,  3361,
    3365,  3369,  3373,  3377,  3381,  3385,  3389,  3393,  3394,  3398,
    3402,  3406,  3410,  3413,  3417,  3421,  3425,  3429,  3433,  3435,
    3439,  3443,  3446,  3450,  3453,  3457,  3458,  3461,  3464,  3468,
    3472,  3476,  3478,  3481,  3485,  3489,  3493,  3497,  3501,  3505,
    3507,  3511,  3514,  3517,  3521,  3524,  3525,  3527,  3529,  3533,
    3536,  3540,  3541,  3545,  3548,  3551,  3555,  3559,  3563,  3567,
    3569,  3570,  3574,  3576,  3578,  3582,  3586,  3590,  3593,  3596,
    3600,  3604,  3608,  3612,  3616,  3620,  3623,  3626,  3630,  3632,

   /* 901-1000 */
    3636,  3640,  3643,  3646,  3648,  3652,  3656,  3660,  3664,  3667,
    3669,  3671,  3675,  3679,  3683,  3687,  3689,  3693,  3694,  3695,
    3699,  3703,  3705,  3707,  3710,  3713,  3717,  3721,  3725,  3729,
    3733,  3736,  3740,  3744,  3748,  3752,  3754,  3757,  3759,  3763,
    3767,  3770,  3773,  3777,  3779,  3783,  3786,  3790,  3794,  3798,
    3801,  3805,  3809,  3813,  3817,  3821,  3825,  3827,  3831,  3835,
    3836,  3837,  3840,  3844,  3848,  3852,  3856,  3859,  3863,  3867,
    3869,  3871,  3875,  3879,  3883,  3887,  3890,  3894,  3898,  3901,
    3905,  3909,  3913,  3917,  3921,  3922,  3923,  3924,  3926,  3930,
    3932,  3936,  3938,  3940,  3944,  3948,  3952,  3956,  3959,  3963,

   /* 1001-1100 */
    3965,  3969,  3973,  3977,  3979,  3981,  3982,  3986,  3989,  3993,
    3997,  4001,  4004,  4006,  4009,  4012,  4016,  4020,  4024,  4026,
    4028,  4032,  4036,  4040,  4044,  4046,  4050,  4054,  4058,  4060,
    4062,  4063,  4064,  4068,  4071,  4075,  4077,  4081,  4083,  4087,
    4089,  4091,  4095,  4099,  4101,  4103,  4105,  4107,  4111,  4115,
    4119,  4123,  4127,  4129,  4131,  4135,  4139,  4141,  4143,  4145,
    4149,  4153,  4157,  4161,  4165,  4169,  4173,  4177,  4180,  4183,
    4187,  4191,  4195,  4198,  4201,  4205,  4209,  4212,  4213,  4216,
    4217,  4221,  4223,  4226,  4230,  4234,  4236,  4240,  4244,  4248,
    4252,  4256,  4258,  4262,  4264,  4266,  4268,  4270,  4272,  4276,

   /* 1101-1200 */
    4279,  4283,  4285,  4287,  4289,  4293,  4295,  4299,  4300,  4301,
    4305,  4309,  4313,  4317,  4319,  4323,  4325,  4329,  4331,  4333,
    4335,  4337,  4341,  4345,  4349,  4351,  4353,  4357,  4361,  4365,
    4367,  4369,  4373,  4377,  4381,  4383,  4387,  4389,  4391,  4395,
    4399,  4403,  4407,  4411,  4413,  4414,  4415,  4418,  4419,  4421,
    4423,  4427,  4429,  4431,  4433,  4435,  4437,  4439,  4443,  4446,
    4450,  4452,  4456,  4458,  4460,  4462,  4466,  4469,  4473,  4477,
    4481,  4483,  4487,  4489,  4491,  4493,  4497,  4499,  4501,  4504,
    4506,  4510,  4513,  4514,  4515,  4518,  4521,  4522,  4525,  4526,
    4527,  4530,  4533,  4534,  4537,  4541,  4542,  4543,  4544,  4545,

   /* 1201-1300 */
    4546,  4547,  4550,  4553,  4554,  4555,  4558,  4561,  4564,  4567,
    4568,  4571,  4574,  4575,  4578,  4581,  4582,  4585,  4586,  4588,
    4590,  4592,  4596,  4598,  4602,  4604,  4608,  4612,  4613,  4616,
    4619,  4622,  4623,  4624,  4625,  4626,  4629,  4632,  4633,  4636,
    4639,  4640,  4641,  4642,  4643,  4644,  4645,  4648,  4649,  4650,
    4651,  4652,  4653,  4656,  4657,  4660,  4661,  4664,  4667,  4670,
    4671,  4674,  4675,  4676,  4677,  4678,  4681,  4682,  4683,  4684,
    4687,  4688,  4689,  4692,  4693,  4696,  4697,  4700,  4701,  4702,
    4703,  4704,  4707,  4708,  4711,  4712,  4715,  4716,  4717,  4718,
    4719,  4720,  4721,  4722,  4723,  4726,  4729,  4730,  4733,  4736,

   /* 1301-(NFLS+NFPL) */
    4737,  4740,  4741,  4742,  4745,  4746,  4749,  4752,  4753
   };

/* Amplitude coefficients (microarcsec);  indexed using the nc array. */
   static const double a[] = {

   /* 1-105 */
         -6844318.44,     9205236.26,1328.67,1538.18,      205833.11,
           153041.79,       -3309.73, 853.32,2037.98,       -2301.27,
       81.46, 120.56, -20.39, -15.22,   1.73,  -1.61,  -0.10,   0.11,
       -0.02,  -0.02,     -523908.04,      573033.42,-544.75,-458.66,
            12814.01,       11714.49, 198.97,-290.91, 155.74,-143.27,
       -2.75,  -1.03,  -1.27,  -1.16,   0.00,  -0.01,      -90552.22,
            97846.69, 111.23, 137.41,2187.91,2024.68,  41.44, -51.26,
       26.92, -24.46,  -0.46,  -0.28,  -0.22,  -0.20,       82168.76,
           -89618.24, -27.64, -29.05,       -2004.36,       -1837.32,
      -36.07,  48.00, -24.43,  22.41,   0.47,   0.24,   0.20,   0.18,
            58707.02,7387.02, 470.05,-192.40, 164.33,       -1312.21,
     -179.73, -28.93, -17.36,  -1.83,  -0.50,   3.57,   0.00,   0.13,
           -20557.78,       22438.42, -20.84, -17.40, 501.82, 459.68,
       59.20, -67.30,   6.08,  -5.61,  -1.36,  -1.19,       28288.28,
     -674.99, -34.69,  35.80, -15.07,-632.54, -11.19,   0.78,  -8.41,
        0.17,   0.01,   0.07,      -15406.85,       20069.50,  15.12,

   /* 106-219 */
       31.80, 448.76, 344.50,  -5.77,   1.41,   4.59,  -5.02,   0.17,
        0.24,      -11991.74,       12902.66,  32.46,  36.70, 288.49,
      268.14,   5.70,  -7.06,   3.57,  -3.23,  -0.06,  -0.04,
            -8584.95,       -9592.72,   4.42, -13.20,-214.50, 192.06,
       23.87,  29.83,   2.54,   2.40,   0.60,  -0.48,5095.50,
            -6918.22,   7.19,   3.92,-154.91,-113.94,   2.86,  -1.04,
       -1.52,   1.73,  -0.07,  -0.10,       -4910.93,       -5331.13,
        0.76,   0.40,-119.21, 109.81,   2.16,   3.20,   1.46,   1.33,
        0.04,  -0.02,       -6245.02,-123.48,  -6.68,  -8.20,  -2.76,
      139.64,   2.71,   0.15,   1.86,2511.85,       -3323.89,   1.07,
       -0.90, -74.33, -56.17,   1.16,  -0.01,  -0.75,   0.83,  -0.02,
       -0.04,2307.58,3143.98,  -7.52,   7.50,  70.31, -51.60,   1.46,
        0.16,  -0.69,  -0.79,   0.02,  -0.05,2372.58,2554.51,   5.93,
       -6.60,  57.12, -53.05,  -0.96,  -1.24,  -0.71,  -0.64,  -0.01,
            -2053.16,2636.13,   5.13,   7.80,  58.94,  45.91,  -0.42,
       -0.12,   0.61,  -0.66,   0.02,   0.03,       -1825.49,

   /* 220-339 */
            -2423.59,   1.23,  -2.00, -54.19,  40.82,  -1.07,  -1.02,
        0.54,   0.61,  -0.04,   0.04,2521.07,-122.28,  -5.97,   2.90,
       -2.73, -56.37,  -0.82,   0.13,  -0.75,       -1534.09,1645.01,
        6.29,   6.80,  36.78,  34.30,   0.92,  -1.25,   0.46,  -0.41,
       -0.02,  -0.01,1898.27,  47.70,  -0.72,   2.50,   1.07, -42.45,
       -0.94,   0.02,  -0.56,       -1292.02,       -1387.00,   0.00,
        0.00, -31.01,  28.89,   0.68,   0.00,   0.38,   0.35,  -0.01,
       -0.01,       -1234.96,1323.81,   5.21,   5.90,  29.60,  27.61,
        0.74,  -1.22,   0.37,  -0.33,  -0.02,  -0.01,1137.48,
            -1233.89,  -0.04,  -0.30, -27.59, -25.43,  -0.61,   1.00,
       -0.34,   0.31,   0.01,   0.01,-813.13,       -1075.60,   0.40,
        0.30, -24.05,  18.18,  -0.40,  -0.01,   0.24,   0.27,  -0.01,
        0.01,1163.22, -60.90,  -2.94,   1.30,  -1.36, -26.01,  -0.58,
        0.07,  -0.35,1029.70, -55.55,  -2.63,   1.10,  -1.25, -23.02,
       -0.52,   0.06,  -0.31,-556.26, 852.85,   3.16,  -4.48,  19.06,
       12.44,  -0.81,  -0.27,   0.17,  -0.21,   0.00,   0.02,-603.52,

   /* 340-467 */
     -800.34,   0.44,   0.10, -17.90,  13.49,  -0.08,  -0.01,   0.18,
        0.20,  -0.01,   0.01,-628.24, 684.99,  -0.64,  -0.50,  15.32,
       14.05,   3.18,  -4.19,   0.19,  -0.17,  -0.09,  -0.07,-866.48,
      -16.26,   0.52,  -1.30,  -0.36,  19.37,   0.43,  -0.01,   0.26,
     -512.37, 695.54,  -1.47,  -1.40,  15.55,  11.46,  -0.16,   0.03,
        0.15,  -0.17,   0.01,   0.01, 506.65, 643.75,   2.54,  -2.62,
       14.40, -11.33,  -0.77,  -0.06,  -0.15,  -0.16,   0.00,   0.01,
      664.57,  16.81,  -0.40,   1.00,   0.38, -14.86,  -3.71,  -0.09,
       -0.20, 405.91, 522.11,   0.99,  -1.50,  11.67,  -9.08,  -0.25,
       -0.02,  -0.12,  -0.13,-305.78, 326.60,   1.75,   1.90,   7.30,
        6.84,   0.20,  -0.04, 300.99,-325.03,  -0.44,  -0.50,  -7.27,
       -6.73,  -1.01,   0.01,   0.00,   0.08,   0.00,   0.02, 438.51,
       10.47,  -0.56,  -0.20,   0.24,  -9.81,  -0.24,   0.01,  -0.13,
     -264.02, 335.24,   0.99,   1.40,   7.49,   5.90,  -0.27,  -0.02,
      284.09, 307.03,   0.32,  -0.40,   6.87,  -6.35,  -0.99,  -0.01,
     -250.54, 327.11,   0.08,   0.40,   7.31,   5.60,  -0.30, 230.72,

   /* 468-595 */
     -304.46,   0.08,  -0.10,  -6.81,  -5.16,   0.27, 229.78, 304.17,
       -0.60,   0.50,   6.80,  -5.14,   0.33,   0.01, 256.30,-276.81,
       -0.28,  -0.40,  -6.19,  -5.73,  -0.14,   0.01,-212.82, 269.45,
        0.84,   1.20,   6.02,   4.76,   0.14,  -0.02, 196.64, 272.05,
       -0.84,   0.90,   6.08,  -4.40,   0.35,   0.02, 188.95, 272.22,
       -0.12,   0.30,   6.09,  -4.22,   0.34,-292.37,  -5.10,  -0.32,
       -0.40,  -0.11,   6.54,   0.14,   0.01, 161.79,-220.67,   0.24,
        0.10,  -4.93,  -3.62,  -0.08, 261.54, -19.94,  -0.95,   0.20,
       -0.45,  -5.85,  -0.13,   0.02, 142.16,-190.79,   0.20,   0.10,
       -4.27,  -3.18,  -0.07, 187.95,  -4.11,  -0.24,   0.30,  -0.09,
       -4.20,  -0.09,   0.01,   0.00,   0.00, -79.08, 167.90,   0.04,
        0.00,   3.75,   1.77, 121.98, 131.04,  -0.08,   0.10,   2.93,
       -2.73,  -0.06,-172.95,  -8.11,  -0.40,  -0.20,  -0.18,   3.87,
        0.09,   0.01,-160.15, -55.30, -14.04,  13.90,  -1.23,   3.58,
        0.40,   0.31,-115.40, 123.20,   0.60,   0.70,   2.75,   2.58,
        0.08,  -0.01,-168.26,  -2.00,   0.20,  -0.20,  -0.04,   3.76,

   /* 596-723 */
        0.08,-114.49, 123.20,   0.32,   0.40,   2.75,   2.56,   0.07,
       -0.01, 112.14, 120.70,   0.28,  -0.30,   2.70,  -2.51,  -0.07,
       -0.01, 161.34,   4.03,   0.20,   0.20,   0.09,  -3.61,  -0.08,
       91.31, 126.64,  -0.40,   0.40,   2.83,  -2.04,  -0.04,   0.01,
      105.29, 112.90,   0.44,  -0.50,   2.52,  -2.35,  -0.07,  -0.01,
       98.69,-106.20,  -0.28,  -0.30,  -2.37,  -2.21,  -0.06,   0.01,
       86.74,-112.94,  -0.08,  -0.20,  -2.53,  -1.94,  -0.05,-134.81,
        3.51,   0.20,  -0.20,   0.08,   3.01,   0.07,  79.03, 107.31,
       -0.24,   0.20,   2.40,  -1.77,  -0.04,   0.01, 132.81, -10.77,
       -0.52,   0.10,  -0.24,  -2.97,  -0.07,   0.01,-130.31,  -0.90,
        0.04,   0.00,   0.00,   2.91, -78.56,  85.32,   0.00,   0.00,
        1.91,   1.76,   0.04,   0.00,   0.00, -41.53,  89.10,   0.02,
        0.00,   1.99,   0.93,  66.03, -71.00,  -0.20,  -0.20,  -1.59,
       -1.48,  -0.04,  60.50,  64.70,   0.36,  -0.40,   1.45,  -1.35,
       -0.04,  -0.01, -52.27, -70.01,   0.00,   0.00,  -1.57,   1.17,
        0.03, -52.95,  66.29,   0.32,   0.40,   1.48,   1.18,   0.04,

   /* 724-851 */
       -0.01,  51.02,  67.25,   0.00,   0.00,   1.50,  -1.14,  -0.03,
      -55.66, -60.92,   0.16,  -0.20,  -1.36,   1.24,   0.03, -54.81,
      -59.20,  -0.08,   0.20,  -1.32,   1.23,   0.03,  51.32, -55.60,
        0.00,   0.00,  -1.24,  -1.15,  -0.03,  48.29,  51.80,   0.20,
       -0.20,   1.16,  -1.08,  -0.03, -45.59, -49.00,  -0.12,   0.10,
       -1.10,   1.02,   0.03,  40.54, -52.69,  -0.04,  -0.10,  -1.18,
       -0.91,  -0.02, -40.58, -49.51,  -1.00,   1.00,  -1.11,   0.91,
        0.04,   0.02, -43.76,  46.50,   0.36,   0.40,   1.04,   0.98,
        0.03,  -0.01,  62.65,  -5.00,  -0.24,   0.00,  -0.11,  -1.40,
       -0.03,   0.01, -38.57,  49.59,   0.08,   0.10,   1.11,   0.86,
        0.02, -33.22, -44.04,   0.08,  -0.10,  -0.98,   0.74,   0.02,
       37.15, -39.90,  -0.12,  -0.10,  -0.89,  -0.83,  -0.02,  36.68,
      -39.50,  -0.04,  -0.10,  -0.88,  -0.82,  -0.02, -53.22,  -3.91,
       -0.20,   0.00,  -0.09,   1.19,   0.03,  32.43, -42.19,  -0.04,
       -0.10,  -0.94,  -0.73,  -0.02, -51.00,  -2.30,  -0.12,  -0.10,
        0.00,   1.14, -29.53, -39.11,   0.04,   0.00,  -0.87,   0.66,

   /* 852-979 */
        0.02,  28.50, -38.92,  -0.08,  -0.10,  -0.87,  -0.64,  -0.02,
       26.54,  36.95,  -0.12,   0.10,   0.83,  -0.59,  -0.01,  26.54,
       34.59,   0.04,  -0.10,   0.77,  -0.59,  -0.02,  28.35, -32.55,
       -0.16,   0.20,  -0.73,  -0.63,  -0.01, -28.00,  30.40,   0.00,
        0.00,   0.68,   0.63,   0.01, -27.61,  29.40,   0.20,   0.20,
        0.66,   0.62,   0.02,  40.33,   0.40,  -0.04,   0.10,   0.00,
       -0.90, -23.28,  31.61,  -0.08,  -0.10,   0.71,   0.52,   0.01,
       37.75,   0.80,   0.04,   0.10,   0.00,  -0.84,  23.66,  25.80,
        0.00,   0.00,   0.58,  -0.53,  -0.01,  21.01, -27.91,   0.00,
        0.00,  -0.62,  -0.47,  -0.01, -34.81,   2.89,   0.04,   0.00,
        0.00,   0.78, -23.49, -25.31,   0.00,   0.00,  -0.57,   0.53,
        0.01, -23.47,  25.20,   0.16,   0.20,   0.56,   0.52,   0.02,
       19.58,  27.50,  -0.12,   0.10,   0.62,  -0.44,  -0.01, -22.67,
      -24.40,  -0.08,   0.10,  -0.55,   0.51,   0.01, -19.97,  25.00,
        0.12,   0.20,   0.56,   0.45,   0.01,  21.28, -22.80,  -0.08,
       -0.10,  -0.51,  -0.48,  -0.01, -30.47,   0.91,   0.04,   0.00,

   /* 980-1107 */
        0.00,   0.68,  18.58,  24.00,   0.04,  -0.10,   0.54,  -0.42,
       -0.01, -18.02,  24.40,  -0.04,  -0.10,   0.55,   0.40,   0.01,
       17.74,  22.50,   0.08,  -0.10,   0.50,  -0.40,  -0.01, -19.41,
       20.70,   0.08,   0.10,   0.46,   0.43,   0.01, -18.64,  20.11,
        0.00,   0.00,   0.45,   0.42,   0.01, -16.75,  21.60,   0.04,
        0.10,   0.48,   0.37,   0.01, -18.42, -20.00,   0.00,   0.00,
       -0.45,   0.41,   0.01, -26.77,   1.41,   0.08,   0.00,   0.00,
        0.60, -26.17,  -0.19,   0.00,   0.00,   0.00,   0.59, -15.52,
       20.51,   0.00,   0.00,   0.46,   0.35,   0.01, -25.42,  -1.91,
       -0.08,   0.00,  -0.04,   0.57,   0.45, -17.42,  18.10,   0.00,
        0.00,   0.40,   0.39,   0.01,  16.39, -17.60,  -0.08,  -0.10,
       -0.39,  -0.37,  -0.01, -14.37,  18.91,   0.00,   0.00,   0.42,
        0.32,   0.01,  23.39,  -2.40,  -0.12,   0.00,   0.00,  -0.52,
       14.32, -18.50,  -0.04,  -0.10,  -0.41,  -0.32,  -0.01,  15.69,
       17.08,   0.00,   0.00,   0.38,  -0.35,  -0.01, -22.99,   0.50,
        0.04,   0.00,   0.00,   0.51,   0.00,   0.00,  14.47, -17.60,

   /* 1108-1235 */
       -0.01,   0.00,  -0.39,  -0.32, -13.33,  18.40,  -0.04,  -0.10,
        0.41,   0.30,  22.47,  -0.60,  -0.04,   0.00,   0.00,  -0.50,
      -12.78, -17.41,   0.04,   0.00,  -0.39,   0.29,   0.01, -14.10,
      -15.31,   0.04,   0.00,  -0.34,   0.32,   0.01,  11.98,  16.21,
       -0.04,   0.00,   0.36,  -0.27,  -0.01,  19.65,  -1.90,  -0.08,
        0.00,   0.00,  -0.44,  19.61,  -1.50,  -0.08,   0.00,   0.00,
       -0.44,  13.41, -14.30,  -0.04,  -0.10,  -0.32,  -0.30,  -0.01,
      -13.29,  14.40,   0.00,   0.00,   0.32,   0.30,   0.01,  11.14,
      -14.40,  -0.04,   0.00,  -0.32,  -0.25,  -0.01,  12.24, -13.38,
        0.04,   0.00,  -0.30,  -0.27,  -0.01,  10.07, -13.81,   0.04,
        0.00,  -0.31,  -0.23,  -0.01,  10.46,  13.10,   0.08,  -0.10,
        0.29,  -0.23,  -0.01,  16.55,  -1.71,  -0.08,   0.00,   0.00,
       -0.37,   9.75, -12.80,   0.00,   0.00,  -0.29,  -0.22,  -0.01,
        9.11,  12.80,   0.00,   0.00,   0.29,  -0.20,   0.00,   0.00,
       -6.44, -13.80,   0.00,   0.00,  -0.31,   0.14,  -9.19, -12.00,
        0.00,   0.00,  -0.27,   0.21, -10.30,  10.90,   0.08,   0.10,

   /* 1236-1363 */
        0.24,   0.23,   0.01,  14.92,  -0.80,  -0.04,   0.00,   0.00,
       -0.33,  10.02, -10.80,   0.00,   0.00,  -0.24,  -0.22,  -0.01,
       -9.75,  10.40,   0.04,   0.00,   0.23,   0.22,   0.01,   9.67,
      -10.40,  -0.04,   0.00,  -0.23,  -0.22,  -0.01,  -8.28, -11.20,
        0.04,   0.00,  -0.25,   0.19,  13.32,  -1.41,  -0.08,   0.00,
        0.00,  -0.30,   8.27,  10.50,   0.04,   0.00,   0.23,  -0.19,
        0.00,   0.00,  13.13,   0.00,   0.00,   0.00,   0.00,  -0.29,
      -12.93,   0.70,   0.04,   0.00,   0.00,   0.29,   7.91, -10.20,
        0.00,   0.00,  -0.23,  -0.18,  -7.84, -10.00,  -0.04,   0.00,
       -0.22,   0.18,   7.44,   9.60,   0.00,   0.00,   0.21,  -0.17,
       -7.64,   9.40,   0.08,   0.10,   0.21,   0.17,   0.01, -11.38,
        0.60,   0.04,   0.00,   0.00,   0.25,  -7.48,   8.30,   0.00,
        0.00,   0.19,   0.17, -10.98,  -0.20,   0.00,   0.00,   0.00,
        0.25,  10.98,   0.20,   0.00,   0.00,   0.00,  -0.25,   7.40,
       -7.90,  -0.04,   0.00,  -0.18,  -0.17,  -6.09,   8.40,  -0.04,
        0.00,   0.19,   0.14,  -6.94,  -7.49,   0.00,   0.00,  -0.17,

   /* 1364-1491 */
        0.16,   6.92,   7.50,   0.04,   0.00,   0.17,  -0.15,   6.20,
        8.09,   0.00,   0.00,   0.18,  -0.14,  -6.12,   7.80,   0.04,
        0.00,   0.17,   0.14,   5.85,  -7.50,   0.00,   0.00,  -0.17,
       -0.13,  -6.48,   6.90,   0.08,   0.10,   0.15,   0.14,   0.01,
        6.32,   6.90,   0.00,   0.00,   0.15,  -0.14,   5.61,  -7.20,
        0.00,   0.00,  -0.16,  -0.13,   9.07,   0.00,   0.00,   0.00,
        0.00,  -0.20,   5.25,   6.90,   0.00,   0.00,   0.15,  -0.12,
       -8.47,  -0.40,   0.00,   0.00,   0.00,   0.19,   6.32,  -5.39,
       -1.11,   1.10,  -0.12,  -0.14,   0.02,   0.02,   5.73,  -6.10,
       -0.04,   0.00,  -0.14,  -0.13,   4.70,   6.60,  -0.04,   0.00,
        0.15,  -0.11,  -4.90,  -6.40,   0.00,   0.00,  -0.14,   0.11,
       -5.33,   5.60,   0.04,   0.10,   0.13,   0.12,   0.01,  -4.81,
        6.00,   0.04,   0.00,   0.13,   0.11,   5.13,   5.50,   0.04,
        0.00,   0.12,  -0.11,   4.50,   5.90,   0.00,   0.00,   0.13,
       -0.10,  -4.22,   6.10,   0.00,   0.00,   0.14,  -4.53,   5.70,
        0.00,   0.00,   0.13,   0.10,   4.18,   5.70,   0.00,   0.00,

   /* 1492-1619 */
        0.13,  -4.75,  -5.19,   0.00,   0.00,  -0.12,   0.11,  -4.06,
        5.60,   0.00,   0.00,   0.13,  -3.98,   5.60,  -0.04,   0.00,
        0.13,   4.02,  -5.40,   0.00,   0.00,  -0.12,   4.49,  -4.90,
       -0.04,   0.00,  -0.11,  -0.10,  -3.62,  -5.40,  -0.16,   0.20,
       -0.12,   0.00,   0.01,   4.38,   4.80,   0.00,   0.00,   0.11,
       -6.40,  -0.10,   0.00,   0.00,   0.00,   0.14,  -3.98,   5.00,
        0.04,   0.00,   0.11,  -3.82,  -5.00,   0.00,   0.00,  -0.11,
       -3.71,   5.07,   0.00,   0.00,   0.11,   4.14,   4.40,   0.00,
        0.00,   0.10,  -6.01,  -0.50,  -0.04,   0.00,   0.00,   0.13,
       -4.04,   4.39,   0.00,   0.00,   0.10,   3.45,  -4.72,   0.00,
        0.00,  -0.11,   3.31,   4.71,   0.00,   0.00,   0.11,   3.26,
       -4.50,   0.00,   0.00,  -0.10,  -3.26,  -4.50,   0.00,   0.00,
       -0.10,  -3.34,  -4.40,   0.00,   0.00,  -0.10,  -3.74,  -4.00,
        3.70,   4.00,   3.34,  -4.30,   3.30,  -4.30,  -3.66,   3.90,
        0.04,   3.66,   3.90,   0.04,  -3.62,  -3.90,  -3.61,   3.90,
       -0.20,   5.30,   0.00,   0.00,   0.12,   3.06,   4.30,   3.30,

   /* 1620-1747 */
        4.00,   0.40,   0.20,   3.10,   4.10,  -3.06,   3.90,  -3.30,
       -3.60,  -3.30,   3.36,   0.01,   3.14,   3.40,  -4.57,  -0.20,
        0.00,   0.00,   0.00,   0.10,  -2.70,  -3.60,   2.94,  -3.20,
       -2.90,   3.20,   2.47,  -3.40,   2.55,  -3.30,   2.80,  -3.08,
        2.51,   3.30,  -4.10,   0.30,  -0.12,  -0.10,   4.10,   0.20,
       -2.74,   3.00,   2.46,   3.23,  -3.66,   1.20,  -0.20,   0.20,
        3.74,  -0.40,  -2.51,  -2.80,  -3.74,   2.27,  -2.90,   0.00,
        0.00,  -2.50,   2.70,  -2.51,   2.60,  -3.50,   0.20,   3.38,
       -2.22,  -2.50,   3.26,  -0.40,   1.95,  -2.60,   3.22,  -0.40,
       -0.04,  -1.79,  -2.60,   1.91,   2.50,   0.74,   3.05,  -0.04,
        0.08,   2.11,  -2.30,  -2.11,   2.20,  -1.87,  -2.40,   2.03,
       -2.20,  -2.03,   2.20,   2.98,   0.00,   0.00,   2.98,  -1.71,
        2.40,   2.94,  -0.10,  -0.12,   0.10,   1.67,   2.40,  -1.79,
        2.30,  -1.79,   2.20,  -1.67,   2.20,   1.79,  -2.00,   1.87,
       -1.90,   1.63,  -2.10,  -1.59,   2.10,   1.55,  -2.10,  -1.55,
        2.10,  -2.59,  -0.20,  -1.75,  -1.90,  -1.75,   1.90,  -1.83,

   /* 1748-1875 */
       -1.80,   1.51,   2.00,  -1.51,  -2.00,   1.71,   1.80,   1.31,
        2.10,  -1.43,   2.00,   1.43,   2.00,  -2.43,  -1.51,   1.90,
       -1.47,   1.90,   2.39,   0.20,  -2.39,   1.39,   1.90,   1.39,
       -1.80,   1.47,  -1.60,   1.47,  -1.60,   1.43,  -1.50,  -1.31,
        1.60,   1.27,  -1.60,  -1.27,   1.60,   1.27,  -1.60,   2.03,
        1.35,   1.50,  -1.39,  -1.40,   1.95,  -0.20,  -1.27,   1.49,
        1.19,   1.50,   1.27,   1.40,   1.15,   1.50,   1.87,  -0.10,
       -1.12,  -1.50,   1.87,  -1.11,  -1.50,  -1.11,  -1.50,   0.00,
        0.00,   1.19,   1.40,   1.27,  -1.30,  -1.27,  -1.30,  -1.15,
        1.40,  -1.23,   1.30,  -1.23,  -1.30,   1.22,  -1.29,   1.07,
       -1.40,   1.75,  -0.20,  -1.03,  -1.40,  -1.07,   1.20,  -1.03,
        1.15,   1.07,   1.10,   1.51,  -1.03,   1.10,   1.03,  -1.10,
        0.00,   0.00,  -1.03,  -1.10,   0.91,  -1.20,  -0.88,  -1.20,
       -0.88,   1.20,  -0.95,   1.10,  -0.95,  -1.10,   1.43,  -1.39,
        0.95,  -1.00,  -0.95,   1.00,  -0.80,   1.10,   0.91,  -1.00,
       -1.35,   0.88,   1.00,  -0.83,   1.00,  -0.91,   0.90,   0.91,

   /* 1876-2003 */
        0.90,   0.88,  -0.90,  -0.76,  -1.00,  -0.76,   1.00,   0.76,
        1.00,  -0.72,   1.00,   0.84,  -0.90,   0.84,   0.90,   1.23,
        0.00,   0.00,  -0.52,  -1.10,  -0.68,   1.00,   1.19,  -0.20,
        1.19,   0.76,   0.90,   1.15,  -0.10,   1.15,  -0.10,   0.72,
       -0.90,  -1.15,  -1.15,   0.68,   0.90,  -0.68,   0.90,  -1.11,
        0.00,   0.00,   0.20,   0.79,   0.80,  -1.11,  -0.10,   0.00,
        0.00,  -0.48,  -1.00,  -0.76,  -0.80,  -0.72,  -0.80,  -1.07,
       -0.10,   0.64,   0.80,  -0.64,  -0.80,   0.64,   0.80,   0.40,
        0.60,   0.52,  -0.50,  -0.60,  -0.80,  -0.71,   0.70,  -0.99,
        0.99,   0.56,   0.80,  -0.56,   0.80,   0.68,  -0.70,   0.68,
        0.70,  -0.95,  -0.64,   0.70,   0.64,   0.70,  -0.60,   0.70,
       -0.60,  -0.70,  -0.91,  -0.10,  -0.51,   0.76,  -0.91,  -0.56,
        0.70,   0.88,   0.88,  -0.63,  -0.60,   0.55,  -0.60,  -0.80,
        0.80,  -0.80,  -0.52,   0.60,   0.52,   0.60,   0.52,  -0.60,
       -0.48,   0.60,   0.48,   0.60,   0.48,   0.60,  -0.76,   0.44,
       -0.60,   0.52,  -0.50,  -0.52,   0.50,   0.40,   0.60,  -0.40,

   /* 2004-2131 */
       -0.60,   0.40,  -0.60,   0.72,  -0.72,  -0.51,  -0.50,  -0.48,
        0.50,   0.48,  -0.50,  -0.48,   0.50,  -0.48,   0.50,   0.48,
       -0.50,  -0.48,  -0.50,  -0.68,  -0.68,   0.44,   0.50,  -0.64,
       -0.10,  -0.64,  -0.10,  -0.40,   0.50,   0.40,   0.50,   0.40,
        0.50,   0.00,   0.00,  -0.40,  -0.50,  -0.36,  -0.50,   0.36,
       -0.50,   0.60,  -0.60,   0.40,  -0.40,   0.40,   0.40,  -0.40,
        0.40,  -0.40,   0.40,  -0.56,  -0.56,   0.36,  -0.40,  -0.36,
        0.40,   0.36,  -0.40,  -0.36,  -0.40,   0.36,   0.40,   0.36,
        0.40,  -0.52,   0.52,   0.52,   0.32,   0.40,  -0.32,   0.40,
       -0.32,   0.40,  -0.32,   0.40,   0.32,  -0.40,  -0.32,  -0.40,
        0.32,  -0.40,   0.28,  -0.40,  -0.28,   0.40,   0.28,  -0.40,
        0.28,   0.40,   0.48,  -0.48,   0.48,   0.36,  -0.30,  -0.36,
       -0.30,   0.00,   0.00,   0.20,   0.40,  -0.44,   0.44,  -0.44,
       -0.44,  -0.44,  -0.44,   0.32,  -0.30,   0.32,   0.30,   0.24,
        0.30,  -0.12,  -0.10,  -0.28,   0.30,   0.28,   0.30,   0.28,
        0.30,   0.28,  -0.30,   0.28,  -0.30,   0.28,  -0.30,   0.28,

   /* 2132-2259 */
        0.30,  -0.28,   0.30,   0.40,   0.40,  -0.24,   0.30,   0.24,
       -0.30,   0.24,  -0.30,  -0.24,  -0.30,   0.24,   0.30,   0.24,
       -0.30,  -0.24,   0.30,   0.24,  -0.30,  -0.24,  -0.30,   0.24,
       -0.30,   0.24,   0.30,  -0.24,   0.30,  -0.24,   0.30,   0.20,
       -0.30,   0.20,  -0.30,   0.20,  -0.30,   0.20,   0.30,   0.20,
       -0.30,   0.20,  -0.30,   0.20,   0.30,   0.20,   0.30,  -0.20,
       -0.30,   0.20,  -0.30,   0.20,  -0.30,  -0.36,  -0.36,  -0.36,
       -0.04,   0.30,   0.12,  -0.10,  -0.32,  -0.24,   0.20,   0.24,
        0.20,   0.20,  -0.20,  -0.20,  -0.20,  -0.20,  -0.20,   0.20,
        0.20,   0.20,  -0.20,   0.20,   0.20,   0.20,   0.20,  -0.20,
       -0.20,   0.00,   0.00,  -0.20,  -0.20,  -0.20,   0.20,  -0.20,
        0.20,   0.20,  -0.20,  -0.20,  -0.20,   0.20,   0.20,   0.20,
        0.20,   0.20,  -0.20,   0.20,  -0.20,   0.28,   0.28,   0.28,
        0.28,   0.28,   0.28,  -0.28,   0.28,   0.12,   0.00,   0.24,
        0.16,  -0.20,   0.16,  -0.20,   0.16,  -0.20,   0.16,   0.20,
       -0.16,   0.20,   0.16,   0.20,  -0.16,   0.20,  -0.16,   0.20,

   /* 2260-2387 */
       -0.16,   0.20,   0.16,  -0.20,   0.16,   0.20,   0.16,  -0.20,
       -0.16,   0.20,  -0.16,  -0.20,  -0.16,   0.20,   0.16,   0.20,
        0.16,  -0.20,   0.16,  -0.20,   0.16,   0.20,   0.16,   0.20,
        0.16,   0.20,  -0.16,  -0.20,   0.16,   0.20,  -0.16,   0.20,
        0.16,   0.20,  -0.16,  -0.20,   0.16,  -0.20,   0.16,  -0.20,
       -0.16,  -0.20,   0.24,  -0.24,  -0.24,   0.24,   0.24,   0.12,
        0.20,   0.12,   0.20,  -0.12,  -0.20,   0.12,  -0.20,   0.12,
       -0.20,  -0.12,   0.20,  -0.12,   0.20,  -0.12,  -0.20,   0.12,
        0.20,   0.12,   0.20,   0.12,  -0.20,  -0.12,   0.20,   0.12,
       -0.20,  -0.12,   0.20,   0.12,   0.20,   0.00,   0.00,  -0.12,
        0.20,  -0.12,   0.20,   0.12,  -0.20,  -0.12,   0.20,   0.12,
        0.20,   0.00,  -0.21,  -0.20,   0.00,   0.00,   0.20,  -0.20,
       -0.20,  -0.20,   0.20,  -0.16,  -0.10,   0.00,   0.17,   0.16,
        0.16,   0.16,   0.16,  -0.16,   0.16,   0.16,  -0.16,   0.16,
       -0.16,   0.16,   0.12,   0.10,   0.12,  -0.10,  -0.12,   0.10,
       -0.12,   0.10,   0.12,  -0.10,  -0.12,   0.12,  -0.12,   0.12,

   /* 2388-2515 */
       -0.12,   0.12,  -0.12,  -0.12,  -0.12,  -0.12,  -0.12,  -0.12,
       -0.12,   0.12,   0.12,   0.12,   0.12,  -0.12,  -0.12,   0.12,
        0.12,   0.12,  -0.12,   0.12,  -0.12,  -0.12,  -0.12,   0.12,
       -0.12,  -0.12,   0.12,   0.00,   0.11,   0.11,-122.67, 164.70,
      203.78, 273.50,   3.58,   2.74,   6.18,  -4.56,   0.00,  -0.04,
        0.00,  -0.07,  57.44, -77.10,  95.82, 128.60,  -1.77,  -1.28,
        2.85,  -2.14,  82.14,  89.50,   0.00,   0.00,   2.00,  -1.84,
       -0.04,  47.73, -64.10,  23.79,  31.90,  -1.45,  -1.07,   0.69,
       -0.53, -46.38,  50.50,   0.00,   0.00,   1.13,   1.04,   0.02,
      -18.38,   0.00,  63.80,   0.00,   0.00,   0.41,   0.00,  -1.43,
       59.07,   0.00,   0.00,   0.00,   0.00,  -1.32,  57.28,   0.00,
        0.00,   0.00,   0.00,  -1.28, -48.65,   0.00,  -1.15,   0.00,
        0.00,   1.09,   0.00,   0.03, -18.30,  24.60, -17.30, -23.20,
        0.56,   0.41,  -0.51,   0.39, -16.91,  26.90,   8.43,  13.30,
        0.60,   0.38,   0.31,  -0.19,   1.23,  -1.70, -19.13, -25.70,
       -0.03,  -0.03,  -0.58,   0.43,  -0.72,   0.90, -17.34, -23.30,

   /* 2516-2643 */
        0.03,   0.02,  -0.52,   0.39, -19.49, -21.30,   0.00,   0.00,
       -0.48,   0.44,   0.01,  20.57, -20.10,   0.64,   0.70,  -0.45,
       -0.46,   0.00,  -0.01,   4.89,   5.90, -16.55,  19.90,   0.14,
       -0.11,   0.44,   0.37,  18.22,  19.80,   0.00,   0.00,   0.44,
       -0.41,  -0.01,   4.89,  -5.30, -16.51, -18.00,  -0.11,  -0.11,
       -0.41,   0.37, -17.86,   0.00,  17.10,   0.00,   0.00,   0.40,
        0.00,  -0.38,   0.32,   0.00,  24.42,   0.00,   0.00,  -0.01,
        0.00,  -0.55, -23.79,   0.00,   0.00,   0.00,   0.00,   0.53,
       14.72, -16.00,  -0.32,   0.00,  -0.36,  -0.33,  -0.01,   0.01,
        3.34,  -4.50,  11.86,  15.90,  -0.11,  -0.07,   0.35,  -0.27,
       -3.26,   4.40,  11.62,  15.60,   0.09,   0.07,   0.35,  -0.26,
      -19.53,   0.00,   5.09,   0.00,   0.00,   0.44,   0.00,  -0.11,
      -13.48,  14.70,   0.00,   0.00,   0.33,   0.30,   0.01,  10.86,
      -14.60,   3.18,   4.30,  -0.33,  -0.24,   0.09,  -0.07, -11.30,
      -15.10,   0.00,   0.00,  -0.34,   0.25,   0.01,   2.03,  -2.70,
       10.82,  14.50,  -0.07,  -0.05,   0.32,  -0.24,  17.46,   0.00,

   /* 2644-2771 */
        0.00,   0.00,   0.00,  -0.39,  16.43,   0.00,   0.52,   0.00,
        0.00,  -0.37,   0.00,  -0.01,   9.35,   0.00,  13.29,   0.00,
        0.00,  -0.21,   0.00,  -0.30, -10.42,  11.40,   0.00,   0.00,
        0.25,   0.23,   0.01,   0.44,   0.50, -10.38,  11.30,   0.02,
       -0.01,   0.25,   0.23, -14.64,   0.00,   0.00,   0.00,   0.00,
        0.33,   0.56,   0.80,  -8.67,  11.70,   0.02,  -0.01,   0.26,
        0.19,  13.88,   0.00,  -2.47,   0.00,   0.00,  -0.31,   0.00,
        0.06,  -1.99,   2.70,   7.72,  10.30,   0.06,   0.04,   0.23,
       -0.17,  -0.20,   0.00,  13.05,   0.00,   0.00,   0.00,   0.00,
       -0.29,   6.92,  -9.30,   3.34,   4.50,  -0.21,  -0.15,   0.10,
       -0.07,  -6.60,   0.00,  10.70,   0.00,   0.00,   0.15,   0.00,
       -0.24,  -8.04,  -8.70,   0.00,   0.00,  -0.19,   0.18, -10.58,
        0.00,  -3.10,   0.00,   0.00,   0.24,   0.00,   0.07,  -7.32,
        8.00,  -0.12,  -0.10,   0.18,   0.16,   1.63,   1.70,   6.96,
       -7.60,   0.03,  -0.04,  -0.17,  -0.16,  -3.62,   0.00,   9.86,
        0.00,   0.00,   0.08,   0.00,  -0.22,   0.20,  -0.20,  -6.88,

   /* 2772-2899 */
       -7.50,   0.00,   0.00,  -0.17,   0.15,  -8.99,   0.00,   4.02,
        0.00,   0.00,   0.20,   0.00,  -0.09,  -1.07,   1.40,  -5.69,
       -7.70,   0.03,   0.02,  -0.17,   0.13,   6.48,  -7.20,  -0.48,
       -0.50,  -0.16,  -0.14,  -0.01,   0.01,   5.57,  -7.50,   1.07,
        1.40,  -0.17,  -0.12,   0.03,  -0.02,   8.71,   0.00,   3.54,
        0.00,   0.00,  -0.19,   0.00,  -0.08,   0.40,   0.00,   9.27,
        0.00,   0.00,  -0.01,   0.00,  -0.21,  -6.13,   6.70,  -1.19,
       -1.30,   0.15,   0.14,  -0.03,   0.03,   5.21,  -5.70,  -2.51,
       -2.60,  -0.13,  -0.12,  -0.06,   0.06,   5.69,  -6.20,  -0.12,
       -0.10,  -0.14,  -0.13,  -0.01,   2.03,  -2.70,   4.53,   6.10,
       -0.06,  -0.05,   0.14,  -0.10,   5.01,   5.50,  -2.51,   2.70,
        0.12,  -0.11,   0.06,   0.06,  -1.91,   2.60,  -4.38,  -5.90,
        0.06,   0.04,  -0.13,   0.10,   4.65,  -6.30,   0.00,   0.00,
       -0.14,  -0.10,  -5.29,   5.70,   0.00,   0.00,   0.13,   0.12,
       -2.23,  -4.00,  -4.65,   4.20,  -0.09,   0.05,   0.10,   0.10,
       -4.53,   6.10,   0.00,   0.00,   0.14,   0.10,   2.47,   2.70,

   /* 2900-3027 */
       -4.46,   4.90,   0.06,  -0.06,   0.11,   0.10,  -5.05,   5.50,
        0.84,   0.90,   0.12,   0.11,   0.02,  -0.02,   4.97,  -5.40,
       -1.71,   0.00,  -0.12,  -0.11,   0.00,   0.04,  -0.99,  -1.30,
        4.22,  -5.70,  -0.03,   0.02,  -0.13,  -0.09,   0.99,   1.40,
        4.22,  -5.60,   0.03,  -0.02,  -0.13,  -0.09,  -4.69,  -5.20,
        0.00,   0.00,  -0.12,   0.10,  -3.42,   0.00,   6.09,   0.00,
        0.00,   0.08,   0.00,  -0.14,  -4.65,  -5.10,   0.00,   0.00,
       -0.11,   0.10,   0.00,   0.00,  -4.53,  -5.00,   0.00,   0.00,
       -0.11,   0.10,  -2.43,  -2.70,  -3.82,   4.20,  -0.06,   0.05,
        0.10,   0.09,   0.00,   0.00,  -4.53,   4.90,   0.00,   0.00,
        0.11,   0.10,  -4.49,  -4.90,   0.00,   0.00,  -0.11,   0.10,
        2.67,  -2.90,  -3.62,  -3.90,  -0.06,  -0.06,  -0.09,   0.08,
        3.94,  -5.30,   0.00,   0.00,  -0.12,  -3.38,   3.70,  -2.78,
       -3.10,   0.08,   0.08,  -0.07,   0.06,   3.18,  -3.50,  -2.82,
       -3.10,  -0.08,  -0.07,  -0.07,   0.06,  -5.77,   0.00,   1.87,
        0.00,   0.00,   0.13,   0.00,  -0.04,   3.54,  -4.80,  -0.64,

   /* 3028-3155 */
       -0.90,  -0.11,   0.00,  -0.02,  -3.50,  -4.70,   0.68,  -0.90,
       -0.11,   0.00,  -0.02,   5.49,   0.00,   0.00,   0.00,   0.00,
       -0.12,   1.83,  -2.50,   2.63,   3.50,  -0.06,   0.00,   0.08,
        3.02,  -4.10,   0.68,   0.90,  -0.09,   0.00,   0.02,   0.00,
        0.00,   5.21,   0.00,   0.00,   0.00,   0.00,  -0.12,  -3.54,
        3.80,   2.70,   3.60,  -1.35,   1.80,   0.08,   0.00,   0.04,
       -2.90,   3.90,   0.68,   0.90,   0.09,   0.00,   0.02,   0.80,
       -1.10,  -2.78,  -3.70,  -0.02,   0.00,  -0.08,   4.10,   0.00,
       -2.39,   0.00,   0.00,  -0.09,   0.00,   0.05,  -1.59,   2.10,
        2.27,   3.00,   0.05,   0.00,   0.07,  -2.63,   3.50,  -0.48,
       -0.60,  -2.94,  -3.20,  -2.94,   3.20,   2.27,  -3.00,  -1.11,
       -1.50,  -0.07,   0.00,  -0.03,  -0.56,  -0.80,  -2.35,   3.10,
        0.00,  -0.60,  -3.42,   1.90,  -0.12,  -0.10,   2.63,  -2.90,
        2.51,   2.80,  -0.64,   0.70,  -0.48,  -0.60,   2.19,  -2.90,
        0.24,  -0.30,   2.15,   2.90,   2.15,  -2.90,   0.52,   0.70,
        2.07,  -2.80,  -3.10,   0.00,   1.79,   0.00,   0.00,   0.07,

   /* 3156-3283 */
        0.00,  -0.04,   0.88,   0.00,  -3.46,   2.11,   2.80,  -0.36,
        0.50,   3.54,  -0.20,  -3.50,  -1.39,   1.50,  -1.91,  -2.10,
       -1.47,   2.00,   1.39,   1.90,   2.07,  -2.30,   0.91,   1.00,
        1.99,  -2.70,   3.30,   0.00,   0.60,  -0.44,  -0.70,  -1.95,
        2.60,   2.15,  -2.40,  -0.60,  -0.70,   3.30,   0.84,   0.00,
       -3.10,  -3.10,   0.00,  -0.72,  -0.32,   0.40,  -1.87,  -2.50,
        1.87,  -2.50,   0.32,   0.40,  -0.24,   0.30,  -1.87,  -2.50,
       -0.24,  -0.30,   1.87,  -2.50,  -2.70,   0.00,   1.55,   2.03,
        2.20,  -2.98,  -1.99,  -2.20,   0.12,  -0.10,  -0.40,   0.50,
        1.59,   2.10,   0.00,   0.00,  -1.79,   2.00,  -1.03,   1.40,
       -1.15,  -1.60,   0.32,   0.50,   1.39,  -1.90,   2.35,  -1.27,
        1.70,   0.60,   0.80,  -0.32,  -0.40,   1.35,  -1.80,   0.44,
        0.00,   2.23,  -0.84,   0.90,  -1.27,  -1.40,  -1.47,   1.60,
       -0.28,  -0.30,  -0.28,   0.40,  -1.27,  -1.70,   0.28,  -0.40,
       -1.43,  -1.50,   0.00,   0.00,  -1.27,  -1.70,   2.11,  -0.32,
       -0.40,  -1.23,   1.60,   1.19,  -1.30,  -0.72,  -0.80,   0.72,

   /* 3284-3411 */
       -0.80,  -1.15,  -1.30,  -1.35,  -1.50,  -1.19,  -1.60,  -0.12,
        0.20,   1.79,   0.00,  -0.88,  -0.28,   0.40,   1.11,   1.50,
       -1.83,   0.00,   0.56,  -0.12,   0.10,  -1.27,  -1.40,   0.00,
        0.00,   1.15,   1.50,  -0.12,   0.20,   1.11,   1.50,   0.36,
       -0.50,  -1.07,  -1.40,  -1.11,   1.50,   1.67,   0.00,   0.80,
       -1.11,   0.00,   1.43,   1.23,  -1.30,  -0.24,  -1.19,  -1.30,
       -0.24,   0.20,  -0.44,  -0.90,  -0.95,   1.10,   1.07,  -1.40,
        1.15,  -1.30,   1.03,  -1.10,  -0.56,  -0.60,  -0.68,   0.90,
       -0.76,  -1.00,  -0.24,  -0.30,   0.95,  -1.30,   0.56,   0.70,
        0.84,  -1.10,  -0.56,   0.00,  -1.55,   0.91,  -1.30,   0.28,
        0.30,   0.16,  -0.20,   0.95,   1.30,   0.40,  -0.50,  -0.88,
       -1.20,   0.95,  -1.10,  -0.48,  -0.50,   0.00,   0.00,  -1.07,
        1.20,   0.44,  -0.50,   0.95,   1.10,   0.00,   0.00,   0.92,
       -1.30,   0.95,   1.00,  -0.52,   0.60,   1.59,   0.24,  -0.40,
        0.91,   1.20,   0.84,  -1.10,  -0.44,  -0.60,   0.84,   1.10,
       -0.44,   0.60,  -0.44,   0.60,  -0.84,  -1.10,  -0.80,   0.00,

   /* 3412-3539 */
        1.35,   0.76,   0.20,  -0.91,  -1.00,   0.20,  -0.30,  -0.91,
       -1.20,  -0.95,   1.00,  -0.48,  -0.50,   0.88,   1.00,   0.48,
       -0.50,  -0.95,  -1.10,   0.20,  -0.20,  -0.99,   1.10,  -0.84,
        1.10,  -0.24,  -0.30,   0.20,  -0.30,   0.84,   1.10,  -1.39,
        0.00,  -0.28,  -0.16,   0.20,   0.84,   1.10,   0.00,   0.00,
        1.39,   0.00,   0.00,  -0.95,   1.00,   1.35,  -0.99,   0.00,
        0.88,  -0.52,   0.00,  -1.19,   0.20,   0.20,   0.76,  -1.00,
        0.00,   0.00,   0.76,   1.00,   0.00,   0.00,   0.76,   1.00,
       -0.76,   1.00,   0.00,   0.00,   1.23,   0.76,   0.80,  -0.32,
        0.40,  -0.72,   0.80,  -0.40,  -0.40,   0.00,   0.00,  -0.80,
       -0.90,  -0.68,   0.90,  -0.16,  -0.20,  -0.16,  -0.20,   0.68,
       -0.90,  -0.36,   0.50,  -0.56,  -0.80,   0.72,  -0.90,   0.44,
       -0.60,  -0.48,  -0.70,  -0.16,   0.00,  -1.11,   0.32,   0.00,
       -1.07,   0.60,  -0.80,  -0.28,  -0.40,  -0.64,   0.00,   0.91,
        1.11,   0.64,  -0.90,   0.76,  -0.80,   0.00,   0.00,  -0.76,
       -0.80,   1.03,   0.00,  -0.36,  -0.64,  -0.70,   0.36,  -0.40,

   /* 3540-3667 */
        1.07,   0.36,  -0.50,  -0.52,  -0.70,   0.60,   0.00,   0.88,
        0.95,   0.00,   0.48,   0.16,  -0.20,   0.60,   0.80,   0.16,
       -0.20,  -0.60,  -0.80,   0.00,  -1.00,   0.12,   0.20,   0.16,
       -0.20,   0.68,   0.70,   0.59,  -0.80,  -0.99,  -0.56,  -0.60,
        0.36,  -0.40,  -0.68,  -0.70,  -0.68,  -0.70,  -0.36,  -0.50,
       -0.44,   0.60,   0.64,   0.70,  -0.12,   0.10,  -0.52,   0.60,
        0.36,   0.40,   0.00,   0.00,   0.95,  -0.84,   0.00,   0.44,
        0.56,   0.60,   0.32,  -0.30,   0.00,   0.00,   0.60,   0.70,
        0.00,   0.00,   0.60,   0.70,  -0.12,  -0.20,   0.52,  -0.70,
        0.00,   0.00,   0.56,   0.70,  -0.12,   0.10,  -0.52,  -0.70,
        0.00,   0.00,   0.88,  -0.76,   0.00,  -0.44,   0.00,   0.00,
       -0.52,  -0.70,   0.52,  -0.70,   0.36,  -0.40,  -0.44,  -0.50,
        0.00,   0.00,   0.60,   0.60,   0.84,   0.00,   0.12,  -0.24,
        0.00,   0.80,  -0.56,   0.60,  -0.32,  -0.30,   0.48,  -0.50,
        0.28,  -0.30,  -0.48,  -0.50,   0.12,   0.20,   0.48,  -0.60,
        0.48,   0.60,  -0.12,   0.20,   0.24,   0.00,   0.76,  -0.52,

   /* 3668-3795 */
       -0.60,  -0.52,   0.60,   0.48,  -0.50,  -0.24,  -0.30,   0.12,
       -0.10,   0.48,   0.60,   0.52,  -0.20,   0.36,   0.40,  -0.44,
        0.50,  -0.24,  -0.30,  -0.48,  -0.60,  -0.44,  -0.60,  -0.12,
        0.10,   0.76,   0.76,   0.20,  -0.20,   0.48,   0.50,   0.40,
       -0.50,  -0.24,  -0.30,   0.44,  -0.60,   0.44,  -0.60,   0.36,
        0.00,  -0.64,   0.72,   0.00,  -0.12,   0.00,  -0.10,  -0.40,
       -0.60,  -0.20,  -0.20,  -0.44,   0.50,  -0.44,   0.50,   0.20,
        0.20,  -0.44,  -0.50,   0.20,  -0.20,  -0.20,   0.20,  -0.44,
       -0.50,   0.64,   0.00,   0.32,  -0.36,   0.50,  -0.20,  -0.30,
        0.12,  -0.10,   0.48,   0.50,  -0.12,   0.30,  -0.36,  -0.50,
        0.00,   0.00,   0.48,   0.50,  -0.48,   0.50,   0.68,   0.00,
       -0.12,   0.56,  -0.40,   0.44,  -0.50,  -0.12,  -0.10,   0.24,
        0.30,  -0.40,   0.40,   0.64,   0.00,  -0.24,   0.64,   0.00,
       -0.20,   0.00,   0.00,   0.44,  -0.50,   0.44,   0.50,  -0.12,
        0.20,  -0.36,  -0.50,   0.12,   0.00,   0.64,  -0.40,   0.50,
        0.00,   0.10,   0.00,   0.00,  -0.40,   0.50,   0.00,   0.00,

   /* 3796-3923 */
       -0.40,  -0.50,   0.56,   0.00,   0.28,   0.00,   0.10,   0.36,
        0.50,   0.00,  -0.10,   0.36,  -0.50,   0.36,   0.50,   0.00,
       -0.10,   0.24,  -0.20,  -0.36,  -0.40,   0.16,   0.20,   0.40,
       -0.40,   0.00,   0.00,  -0.36,  -0.50,  -0.36,  -0.50,  -0.32,
       -0.50,  -0.12,   0.10,   0.20,   0.20,  -0.36,   0.40,  -0.60,
        0.60,   0.28,   0.00,   0.52,   0.12,  -0.10,   0.40,   0.40,
        0.00,  -0.50,   0.20,  -0.20,  -0.32,   0.40,   0.16,   0.20,
       -0.16,   0.20,   0.32,   0.40,   0.56,   0.00,  -0.12,   0.32,
       -0.40,  -0.16,  -0.20,   0.00,   0.00,   0.40,   0.40,  -0.40,
       -0.40,  -0.40,   0.40,  -0.36,   0.40,   0.12,   0.10,   0.00,
        0.10,   0.36,   0.40,   0.00,  -0.10,   0.36,   0.40,  -0.36,
        0.40,   0.00,   0.10,   0.32,   0.00,   0.44,   0.12,   0.20,
        0.28,  -0.40,   0.00,   0.00,   0.36,   0.40,   0.32,  -0.40,
       -0.16,   0.12,   0.10,   0.32,  -0.40,   0.20,   0.30,  -0.24,
        0.30,   0.00,   0.10,   0.32,   0.40,   0.00,  -0.10,  -0.32,
       -0.40,  -0.32,   0.40,   0.00,   0.10,  -0.52,  -0.52,   0.52,

   /* 3924-4051 */
        0.32,  -0.40,   0.00,   0.00,   0.32,   0.40,   0.32,  -0.40,
        0.00,   0.00,  -0.32,  -0.40,  -0.32,   0.40,   0.32,   0.40,
        0.00,   0.00,   0.32,   0.40,   0.00,   0.00,  -0.32,  -0.40,
        0.00,   0.00,   0.32,   0.40,   0.16,   0.20,   0.32,  -0.30,
       -0.16,   0.00,  -0.48,  -0.20,   0.20,  -0.28,  -0.30,   0.28,
       -0.40,   0.00,   0.00,   0.28,  -0.40,   0.00,   0.00,   0.28,
       -0.40,   0.00,   0.00,  -0.28,  -0.40,   0.28,   0.40,  -0.28,
       -0.40,  -0.48,  -0.20,   0.20,   0.24,   0.30,   0.44,   0.00,
        0.16,   0.24,   0.30,   0.16,  -0.20,   0.24,   0.30,  -0.12,
        0.20,   0.20,   0.30,  -0.16,   0.20,   0.00,   0.00,   0.44,
       -0.32,   0.30,   0.24,   0.00,  -0.36,   0.36,   0.00,   0.24,
        0.12,  -0.20,   0.20,   0.30,  -0.12,   0.00,  -0.28,   0.30,
       -0.24,   0.30,   0.12,   0.10,  -0.28,  -0.30,  -0.28,   0.30,
        0.00,   0.00,  -0.28,  -0.30,   0.00,   0.00,  -0.28,  -0.30,
        0.00,   0.00,   0.28,   0.30,   0.00,   0.00,  -0.28,  -0.30,
       -0.28,   0.30,   0.00,   0.00,  -0.28,  -0.30,   0.00,   0.00,

   /* 4052-4179 */
        0.28,   0.30,   0.00,   0.00,  -0.28,   0.30,   0.28,  -0.30,
       -0.28,   0.30,   0.40,   0.40,  -0.24,   0.30,   0.00,  -0.10,
        0.16,   0.00,   0.36,  -0.20,   0.30,  -0.12,  -0.10,  -0.24,
       -0.30,   0.00,   0.00,  -0.24,   0.30,  -0.24,   0.30,   0.00,
        0.00,  -0.24,   0.30,  -0.24,   0.30,   0.24,  -0.30,   0.00,
        0.00,   0.24,  -0.30,   0.00,   0.00,   0.24,   0.30,   0.24,
       -0.30,   0.24,   0.30,  -0.24,   0.30,  -0.24,   0.30,  -0.20,
        0.20,  -0.16,  -0.20,   0.00,   0.00,  -0.32,   0.20,   0.00,
        0.10,   0.20,  -0.30,   0.20,  -0.20,   0.12,   0.20,  -0.16,
        0.20,   0.16,   0.20,   0.20,   0.30,   0.20,   0.30,   0.00,
        0.00,  -0.20,   0.30,   0.00,   0.00,   0.20,   0.30,  -0.20,
       -0.30,  -0.20,  -0.30,   0.20,  -0.30,   0.00,   0.00,   0.20,
        0.30,   0.00,   0.00,   0.20,   0.30,   0.00,   0.00,   0.20,
        0.30,   0.00,   0.00,   0.20,   0.30,   0.00,   0.00,   0.20,
       -0.30,   0.00,   0.00,  -0.20,  -0.30,   0.00,   0.00,  -0.20,
        0.30,   0.00,   0.00,  -0.20,   0.30,   0.00,   0.00,   0.36,

   /* 4180-4307 */
        0.00,   0.00,   0.36,   0.12,   0.10,  -0.24,   0.20,   0.12,
       -0.20,  -0.16,  -0.20,  -0.13,   0.10,   0.22,   0.21,   0.20,
        0.00,  -0.28,   0.32,   0.00,  -0.12,  -0.20,  -0.20,   0.12,
       -0.10,   0.12,   0.10,  -0.20,   0.20,   0.00,   0.00,  -0.32,
        0.32,   0.00,   0.00,   0.32,   0.32,   0.00,   0.00,  -0.24,
       -0.20,   0.24,   0.20,   0.20,   0.00,  -0.24,   0.00,   0.00,
       -0.24,  -0.20,   0.00,   0.00,   0.24,   0.20,  -0.24,  -0.20,
        0.00,   0.00,  -0.24,   0.20,   0.16,  -0.20,   0.12,   0.10,
        0.20,   0.20,   0.00,  -0.10,  -0.12,   0.10,  -0.16,  -0.20,
       -0.12,  -0.10,  -0.16,   0.20,   0.20,   0.20,   0.00,   0.00,
       -0.20,   0.20,  -0.20,   0.20,  -0.20,   0.20,  -0.20,   0.20,
        0.20,  -0.20,  -0.20,  -0.20,   0.00,   0.00,  -0.20,   0.20,
        0.20,   0.00,  -0.20,   0.00,   0.00,  -0.20,   0.20,  -0.20,
        0.20,  -0.20,  -0.20,  -0.20,  -0.20,   0.00,   0.00,   0.20,
        0.20,   0.20,   0.20,   0.12,  -0.20,  -0.12,  -0.10,   0.28,
       -0.28,   0.16,  -0.20,   0.00,  -0.10,   0.00,   0.10,  -0.16,

   /* 4308-4435 */
        0.20,   0.00,  -0.10,  -0.16,  -0.20,   0.00,  -0.10,   0.16,
       -0.20,   0.16,  -0.20,   0.00,   0.00,   0.16,   0.20,  -0.16,
        0.20,   0.00,   0.00,   0.16,   0.20,   0.16,  -0.20,   0.16,
       -0.20,  -0.16,   0.20,   0.16,  -0.20,   0.00,   0.00,   0.16,
        0.20,   0.00,   0.00,   0.16,   0.20,   0.00,   0.00,  -0.16,
       -0.20,   0.16,  -0.20,  -0.16,  -0.20,   0.00,   0.00,  -0.16,
       -0.20,   0.00,   0.00,  -0.16,   0.20,   0.00,   0.00,   0.16,
       -0.20,   0.16,   0.20,   0.16,   0.20,   0.00,   0.00,  -0.16,
       -0.20,   0.00,   0.00,  -0.16,  -0.20,   0.00,   0.00,   0.16,
        0.20,   0.16,   0.20,   0.00,   0.00,   0.16,   0.20,   0.16,
       -0.20,   0.16,   0.20,   0.00,   0.00,  -0.16,   0.20,   0.00,
        0.10,   0.12,  -0.20,   0.12,  -0.20,   0.00,  -0.10,   0.00,
       -0.10,   0.12,   0.20,   0.00,  -0.10,  -0.12,   0.20,  -0.15,
        0.20,  -0.24,   0.24,   0.00,   0.00,   0.24,   0.24,   0.12,
       -0.20,  -0.12,  -0.20,   0.00,   0.00,   0.12,   0.20,   0.12,
       -0.20,   0.12,   0.20,   0.12,   0.20,   0.12,   0.20,   0.12,

   /* 4436-4563 */
       -0.20,  -0.12,   0.20,   0.00,   0.00,   0.12,   0.20,   0.12,
        0.00,  -0.20,   0.00,   0.00,  -0.12,  -0.20,   0.12,  -0.20,
        0.00,   0.00,   0.12,   0.20,  -0.12,   0.20,  -0.12,   0.20,
        0.12,  -0.20,   0.00,   0.00,   0.12,   0.20,   0.20,   0.00,
        0.12,   0.00,   0.00,  -0.12,   0.20,   0.00,   0.00,  -0.12,
       -0.20,   0.00,   0.00,  -0.12,  -0.20,  -0.12,  -0.20,   0.00,
        0.00,   0.12,  -0.20,   0.12,  -0.20,   0.12,   0.20,  -0.12,
       -0.20,   0.00,   0.00,   0.12,  -0.20,   0.12,  -0.20,   0.12,
        0.20,   0.12,   0.00,   0.20,  -0.12,  -0.20,   0.00,   0.00,
        0.12,   0.20,  -0.16,   0.00,   0.16,  -0.20,   0.20,   0.00,
        0.00,  -0.20,   0.00,   0.00,  -0.20,   0.20,   0.00,   0.00,
        0.20,   0.20,  -0.20,   0.00,   0.00,  -0.20,   0.12,   0.00,
       -0.16,   0.20,   0.00,   0.00,   0.20,   0.12,  -0.10,   0.00,
        0.10,   0.16,  -0.16,  -0.16,  -0.16,  -0.16,  -0.16,   0.00,
        0.00,  -0.16,   0.00,   0.00,  -0.16,  -0.16,  -0.16,   0.00,
        0.00,  -0.16,   0.00,   0.00,   0.16,   0.00,   0.00,   0.16,

   /* 4564-4691 */
        0.00,   0.00,   0.16,   0.16,   0.00,   0.00,  -0.16,   0.00,
        0.00,  -0.16,  -0.16,   0.00,   0.00,   0.16,   0.00,   0.00,
       -0.16,  -0.16,   0.00,   0.00,  -0.16,  -0.16,   0.12,   0.10,
        0.12,  -0.10,   0.12,   0.10,   0.00,   0.00,   0.12,   0.10,
       -0.12,   0.10,   0.00,   0.00,   0.12,   0.10,   0.12,  -0.10,
        0.00,   0.00,  -0.12,  -0.10,   0.00,   0.00,   0.12,   0.10,
        0.12,   0.00,   0.00,   0.12,   0.00,   0.00,  -0.12,   0.00,
        0.00,   0.12,   0.12,   0.12,   0.12,   0.12,   0.00,   0.00,
        0.12,   0.00,   0.00,   0.12,   0.12,   0.00,   0.00,   0.12,
        0.00,   0.00,   0.12,  -0.12,  -0.12,   0.12,   0.12,  -0.12,
       -0.12,   0.00,   0.00,   0.12,  -0.12,   0.12,   0.12,  -0.12,
       -0.12,   0.00,   0.00,  -0.12,  -0.12,   0.00,   0.00,  -0.12,
        0.12,   0.00,   0.00,   0.12,   0.00,   0.00,   0.12,   0.00,
        0.00,   0.12,  -0.12,   0.00,   0.00,  -0.12,   0.12,  -0.12,
       -0.12,   0.12,   0.00,   0.00,   0.12,   0.12,   0.12,  -0.12,
        0.00,   0.00,  -0.12,  -0.12,  -0.12,   0.00,   0.00,  -0.12,

   /* 4692-NA */
       -0.12,   0.00,   0.00,   0.12,   0.12,   0.00,   0.00,  -0.12,
       -0.12,  -0.12,  -0.12,   0.12,   0.00,   0.00,   0.12,  -0.12,
        0.00,   0.00,  -0.12,  -0.12,   0.00,   0.00,   0.12,  -0.12,
       -0.12,  -0.12,  -0.12,   0.12,   0.12,  -0.12,  -0.12,   0.00,
        0.00,  -0.12,   0.00,   0.00,  -0.12,   0.12,   0.00,   0.00,
        0.12,   0.00,   0.00,  -0.12,  -0.12,   0.00,   0.00,  -0.12,
       -0.12,   0.12,   0.00,   0.00,   0.12,   0.12,   0.00,   0.00,
        0.12,   0.00,   0.00,   0.12,   0.12,   0.08,   0.00,   0.04
   };

/* Number of amplitude coefficients */
   static const int NA = (int) (sizeof a / sizeof (double));

/* Amplitude usage: X or Y, sin or cos, power of T. */
   static const int jaxy[] = {0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1};
   static const int jasc[] = {0,1,1,0,1,0,0,1,0,1,1,0,1,0,0,1,0,1,1,0};
   static const int japt[] = {0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4};

/* Miscellaneous */
   double t, w, pt[MAXPT+1], fa[14], xypr[2], xypl[2], xyls[2], arg,
          sc[2];
   int jpt, i, j, jxy, ialast, ifreq, m, ia, jsc;

/*--------------------------------------------------------------------*/

/* Interval between fundamental date J2000.0 and given date (JC). */
   t = ((date1 - DJ00) + date2) / DJC;

/* Powers of T. */
   w = 1.0;
   for (jpt = 0; jpt <= MAXPT; jpt++) {
      pt[jpt] = w;
      w *= t;
   }

/* Initialize totals in X and Y:  polynomial, luni-solar, planetary. */
   for (jxy = 0; jxy < 2; jxy++) {
      xypr[jxy] = 0.0;
      xyls[jxy] = 0.0;
      xypl[jxy] = 0.0;
   }

/* --------------------------------- */
/* Fundamental arguments (IERS 2003) */
/* --------------------------------- */

/* Mean anomaly of the Moon. */
   fa[0] = eraFal03(t);

/* Mean anomaly of the Sun. */
   fa[1] = eraFalp03(t);

/* Mean argument of the latitude of the Moon. */
   fa[2] = eraFaf03(t);

/* Mean elongation of the Moon from the Sun. */
   fa[3] = eraFad03(t);

/* Mean longitude of the ascending node of the Moon. */
   fa[4] = eraFaom03(t);

/* Planetary longitudes, Mercury through Neptune. */
   fa[5] = eraFame03(t);
   fa[6] = eraFave03(t);
   fa[7] = eraFae03(t);
   fa[8] = eraFama03(t);
   fa[9] = eraFaju03(t);
   fa[10] = eraFasa03(t);
   fa[11] = eraFaur03(t);
   fa[12] = eraFane03(t);

/* General accumulated precession in longitude. */
   fa[13] = eraFapa03(t);

/* -------------------------------------- */
/* Polynomial part of precession-nutation */
/* -------------------------------------- */

   for (jxy = 0; jxy < 2; jxy++) {
      for (j = MAXPT; j >= 0; j--) {
         xypr[jxy] += xyp[jxy][j] * pt[j];
      }
   }

/* ---------------------------------- */
/* Nutation periodic terms, planetary */
/* ---------------------------------- */

/* Work backwards through the coefficients per frequency list. */
   ialast = NA;
   for (ifreq = NFPL-1; ifreq >= 0; ifreq--) {

   /* Obtain the argument functions. */
      arg = 0.0;
      for (i = 0; i < 14; i++) {
         m = mfapl[ifreq][i];
         if (m != 0) arg += (double)m * fa[i];
      }
      sc[0] = sin(arg);
      sc[1] = cos(arg);

   /* Work backwards through the amplitudes at this frequency. */
      ia = nc[ifreq+NFLS];
      for (i = ialast; i >= ia; i--) {

      /* Coefficient number (0 = 1st). */
         j = i-ia;

      /* X or Y. */
         jxy = jaxy[j];

      /* Sin or cos. */
         jsc = jasc[j];

      /* Power of T. */
         jpt = japt[j];

      /* Accumulate the component. */
         xypl[jxy] += a[i-1] * sc[jsc] * pt[jpt];
      }
      ialast = ia-1;
   }

/* ----------------------------------- */
/* Nutation periodic terms, luni-solar */
/* ----------------------------------- */

/* Continue working backwards through the number of coefficients list. */
   for (ifreq = NFLS-1; ifreq >= 0; ifreq--) {

   /* Obtain the argument functions. */
      arg = 0.0;
      for (i = 0; i < 5; i++) {
         m = mfals[ifreq][i];
         if (m != 0) arg += (double)m * fa[i];
      }
      sc[0] = sin(arg);
      sc[1] = cos(arg);

   /* Work backwards through the amplitudes at this frequency. */
      ia = nc[ifreq];
      for (i = ialast; i >= ia; i--) {

      /* Coefficient number (0 = 1st). */
         j = i-ia;

      /* X or Y. */
         jxy = jaxy[j];

      /* Sin or cos. */
         jsc = jasc[j];

      /* Power of T. */
         jpt = japt[j];

      /* Accumulate the component. */
         xyls[jxy] += a[i-1] * sc[jsc] * pt[jpt];
      }
      ialast = ia-1;
   }

/* ------------------------------------ */
/* Results:  CIP unit vector components */
/* ------------------------------------ */

   *x = DAS2R * (xypr[0] + (xyls[0] + xypl[0]) / 1e6);
   *y = DAS2R * (xypr[1] + (xyls[1] + xypl[1]) / 1e6);

   return;

}

void eraXys00a(double date1, double date2,
               double *x, double *y, double *s)
/*
**  - - - - - - - - - -
**   e r a X y s 0 0 a
**  - - - - - - - - - -
**
**  For a given TT date, compute the X,Y coordinates of the Celestial
**  Intermediate Pole and the CIO locator s, using the IAU 2000A
**  precession-nutation model.
**
**  Given:
**     date1,date2  double   TT as a 2-part Julian Date (Note 1)
**
**  Returned:
**     x,y          double   Celestial Intermediate Pole (Note 2)
**     s            double   the CIO locator s (Note 2)
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
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
**  2) The Celestial Intermediate Pole coordinates are the x,y
**     components of the unit vector in the Geocentric Celestial
**     Reference System.
**
**  3) The CIO locator s (in radians) positions the Celestial
**     Intermediate Origin on the equator of the CIP.
**
**  4) A faster, but slightly less accurate result (about 1 mas for
**     X,Y), can be obtained by using instead the eraXys00b function.
**
**  Called:
**     eraPnm00a    classical NPB matrix, IAU 2000A
**     eraBpn2xy    extract CIP X,Y coordinates from NPB matrix
**     eraS00       the CIO locator s, given X,Y, IAU 2000A
**
**  Reference:
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double rbpn[3][3];


/* Form the bias-precession-nutation matrix, IAU 2000A. */
   eraPnm00a(date1, date2, rbpn);

/* Extract X,Y. */
   eraBpn2xy(rbpn, x, y);

/* Obtain s. */
   *s = eraS00(date1, date2, *x, *y);

   return;

}

void eraXys00b(double date1, double date2,
               double *x, double *y, double *s)
/*
**  - - - - - - - - - -
**   e r a X y s 0 0 b
**  - - - - - - - - - -
**
**  For a given TT date, compute the X,Y coordinates of the Celestial
**  Intermediate Pole and the CIO locator s, using the IAU 2000B
**  precession-nutation model.
**
**  Given:
**     date1,date2  double   TT as a 2-part Julian Date (Note 1)
**
**  Returned:
**     x,y          double   Celestial Intermediate Pole (Note 2)
**     s            double   the CIO locator s (Note 2)
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
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
**  2) The Celestial Intermediate Pole coordinates are the x,y
**     components of the unit vector in the Geocentric Celestial
**     Reference System.
**
**  3) The CIO locator s (in radians) positions the Celestial
**     Intermediate Origin on the equator of the CIP.
**
**  4) The present function is faster, but slightly less accurate (about
**     1 mas in X,Y), than the eraXys00a function.
**
**  Called:
**     eraPnm00b    classical NPB matrix, IAU 2000B
**     eraBpn2xy    extract CIP X,Y coordinates from NPB matrix
**     eraS00       the CIO locator s, given X,Y, IAU 2000A
**
**  Reference:
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double rbpn[3][3];


/* Form the bias-precession-nutation matrix, IAU 2000A. */
   eraPnm00b(date1, date2, rbpn);

/* Extract X,Y. */
   eraBpn2xy(rbpn, x, y);

/* Obtain s. */
   *s = eraS00(date1, date2, *x, *y);

   return;

}

void eraXys06a(double date1, double date2,
               double *x, double *y, double *s)
/*
**  - - - - - - - - - -
**   e r a X y s 0 6 a
**  - - - - - - - - - -
**
**  For a given TT date, compute the X,Y coordinates of the Celestial
**  Intermediate Pole and the CIO locator s, using the IAU 2006
**  precession and IAU 2000A nutation models.
**
**  Given:
**     date1,date2  double  TT as a 2-part Julian Date (Note 1)
**
**  Returned:
**     x,y          double  Celestial Intermediate Pole (Note 2)
**     s            double  the CIO locator s (Note 2)
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
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
**  2) The Celestial Intermediate Pole coordinates are the x,y components
**     of the unit vector in the Geocentric Celestial Reference System.
**
**  3) The CIO locator s (in radians) positions the Celestial
**     Intermediate Origin on the equator of the CIP.
**
**  4) Series-based solutions for generating X and Y are also available:
**     see Capitaine & Wallace (2006) and eraXy06.
**
**  Called:
**     eraPnm06a    classical NPB matrix, IAU 2006/2000A
**     eraBpn2xy    extract CIP X,Y coordinates from NPB matrix
**     eraS06       the CIO locator s, given X,Y, IAU 2006
**
**  References:
**
**     Capitaine, N. & Wallace, P.T., 2006, Astron.Astrophys. 450, 855
**
**     Wallace, P.T. & Capitaine, N., 2006, Astron.Astrophys. 459, 981
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double rbpn[3][3];


/* Form the bias-precession-nutation matrix, IAU 2000A. */
   eraPnm06a(date1, date2, rbpn);

/* Extract X,Y. */
   eraBpn2xy(rbpn, x, y);

/* Obtain s. */
   *s = eraS06(date1, date2, *x, *y);

   return;

}

void eraZp(double p[3])
/*
**  - - - - - -
**   e r a Z p
**  - - - - - -
**
**  Zero a p-vector.
**
**  Returned:
**     p        double[3]      p-vector
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   p[0] = 0.0;
   p[1] = 0.0;
   p[2] = 0.0;

   return;

}

void eraZpv(double pv[2][3])
/*
**  - - - - - - -
**   e r a Z p v
**  - - - - - - -
**
**  Zero a pv-vector.
**
**  Returned:
**     pv       double[2][3]      pv-vector
**
**  Called:
**     eraZp        zero p-vector
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   eraZp(pv[0]);
   eraZp(pv[1]);

   return;

}

void eraZr(double r[3][3])
/*
**  - - - - - -
**   e r a Z r
**  - - - - - -
**
**  Initialize an r-matrix to the null matrix.
**
**  Returned:
**     r        double[3][3]    r-matrix
**
**  Copyright (C) 2013, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   r[0][0] = 0.0;
   r[0][1] = 0.0;
   r[0][2] = 0.0;
   r[1][0] = 0.0;
   r[1][1] = 0.0;
   r[1][2] = 0.0;
   r[2][0] = 0.0;
   r[2][1] = 0.0;
   r[2][2] = 0.0;

   return;

}

/*----------------------------------------------------------------------
**  
**  
**  Copyright (C) 2013, NumFOCUS Foundation.
**  All rights reserved.
**  
**  This library is derived, with permission, from the International
**  Astronomical Union's "Standards of Fundamental Astronomy" library,
**  available from http://www.iausofa.org.
**  
**  The ERFA version is intended to retain identical functionality to
**  the SOFA library, but made distinct through different function and
**  file names, as set out in the SOFA license conditions.  The SOFA
**  original has a role as a reference standard for the IAU and IERS,
**  and consequently redistribution is permitted only in its unaltered
**  state.  The ERFA version is not subject to this restriction and
**  therefore can be included in distributions which do not support the
**  concept of "read only" software.
**  
**  Although the intent is to replicate the SOFA API (other than
**  replacement of prefix names) and results (with the exception of
**  bugs;  any that are discovered will be fixed), SOFA is not
**  responsible for any errors found in this version of the library.
**  
**  If you wish to acknowledge the SOFA heritage, please acknowledge
**  that you are using a library derived from SOFA, rather than SOFA
**  itself.
**  
**  
**  TERMS AND CONDITIONS
**  
**  Redistribution and use in source and binary forms, with or without
**  modification, are permitted provided that the following conditions
**  are met:
**  
**  1 Redistributions of source code must retain the above copyright
**    notice, this list of conditions and the following disclaimer.
**  
**  2 Redistributions in binary form must reproduce the above copyright
**    notice, this list of conditions and the following disclaimer in
**    the documentation and/or other materials provided with the
**    distribution.
**  
**  3 Neither the name of the Standards Of Fundamental Astronomy Board,
**    the International Astronomical Union nor the names of its
**    contributors may be used to endorse or promote products derived
**    from this software without specific prior written permission.
**  
**  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
**  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
**  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
**  FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE
**  COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
**  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
**  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
**  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
**  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
**  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
**  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
**  POSSIBILITY OF SUCH DAMAGE.
**  
*/
