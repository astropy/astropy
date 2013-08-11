#ifndef ERFAHDEF
#define ERFAHDEF

#include <math.h>

/*
**  - - - - - - -
**   e r f a . h
**  - - - - - - -
**
**  Prototype function declarations and macros for erfa library.
**
*/

#ifdef __cplusplus
extern "C" {
#endif

/* Astronomy/Calendars */
int eraCal2jd(int iy, int im, int id, double *djm0, double *djm);
double eraEpb(double dj1, double dj2);
void eraEpb2jd(double epb, double *djm0, double *djm);
double eraEpj(double dj1, double dj2);
void eraEpj2jd(double epj, double *djm0, double *djm);
int eraJd2cal(double dj1, double dj2,
                     int *iy, int *im, int *id, double *fd);
int eraJdcalf(int ndp, double dj1, double dj2, int iymdf[4]);

/* Astronomy/Ephemerides */
int eraEpv00(double date1, double date2,
             double pvh[2][3], double pvb[2][3]);
int eraPlan94(double date1, double date2, int np, double pv[2][3]);

/* Astronomy/FundamentalArgs */
double eraFad03(double t);
double eraFae03(double t);
double eraFaf03(double t);
double eraFaju03(double t);
double eraFal03(double t);
double eraFalp03(double t);
double eraFama03(double t);
double eraFame03(double t);
double eraFane03(double t);
double eraFaom03(double t);
double eraFapa03(double t);
double eraFasa03(double t);
double eraFaur03(double t);
double eraFave03(double t);

/* Astronomy/PrecNutPolar */
void eraBi00(double *dpsibi, double *depsbi, double *dra);
void eraBp00(double date1, double date2,
             double rb[3][3], double rp[3][3], double rbp[3][3]);
void eraBp06(double date1, double date2,
             double rb[3][3], double rp[3][3], double rbp[3][3]);
void eraBpn2xy(double rbpn[3][3], double *x, double *y);
void eraC2i00a(double date1, double date2, double rc2i[3][3]);
void eraC2i00b(double date1, double date2, double rc2i[3][3]);
void eraC2i06a(double date1, double date2, double rc2i[3][3]);
void eraC2ibpn(double date1, double date2, double rbpn[3][3],
               double rc2i[3][3]);
void eraC2ixy(double date1, double date2, double x, double y,
              double rc2i[3][3]);
void eraC2ixys(double x, double y, double s, double rc2i[3][3]);
void eraC2t00a(double tta, double ttb, double uta, double utb,
               double xp, double yp, double rc2t[3][3]);
void eraC2t00b(double tta, double ttb, double uta, double utb,
               double xp, double yp, double rc2t[3][3]);
void eraC2t06a(double tta, double ttb, double uta, double utb,
               double xp, double yp, double rc2t[3][3]);
void eraC2tcio(double rc2i[3][3], double era, double rpom[3][3],
               double rc2t[3][3]);
void eraC2teqx(double rbpn[3][3], double gst, double rpom[3][3],
               double rc2t[3][3]);
void eraC2tpe(double tta, double ttb, double uta, double utb,
              double dpsi, double deps, double xp, double yp,
              double rc2t[3][3]);
void eraC2txy(double tta, double ttb, double uta, double utb,
              double x, double y, double xp, double yp,
              double rc2t[3][3]);
double eraEo06a(double date1, double date2);
double eraEors(double rnpb[3][3], double s);
void eraFw2m(double gamb, double phib, double psi, double eps,
             double r[3][3]);
void eraFw2xy(double gamb, double phib, double psi, double eps,
              double *x, double *y);
void eraNum00a(double date1, double date2, double rmatn[3][3]);
void eraNum00b(double date1, double date2, double rmatn[3][3]);
void eraNum06a(double date1, double date2, double rmatn[3][3]);
void eraNumat(double epsa, double dpsi, double deps, double rmatn[3][3]);
void eraNut00a(double date1, double date2, double *dpsi, double *deps);
void eraNut00b(double date1, double date2, double *dpsi, double *deps);
void eraNut06a(double date1, double date2, double *dpsi, double *deps);
void eraNut80(double date1, double date2, double *dpsi, double *deps);
void eraNutm80(double date1, double date2, double rmatn[3][3]);
double eraObl06(double date1, double date2);
double eraObl80(double date1, double date2);
void eraP06e(double date1, double date2,
             double *eps0, double *psia, double *oma, double *bpa,
             double *bqa, double *pia, double *bpia,
             double *epsa, double *chia, double *za, double *zetaa,
             double *thetaa, double *pa,
             double *gam, double *phi, double *psi);
void eraPb06(double date1, double date2,
             double *bzeta, double *bz, double *btheta);
void eraPfw06(double date1, double date2,
              double *gamb, double *phib, double *psib, double *epsa);
void eraPmat00(double date1, double date2, double rbp[3][3]);
void eraPmat06(double date1, double date2, double rbp[3][3]);
void eraPmat76(double date1, double date2, double rmatp[3][3]);
void eraPn00(double date1, double date2, double dpsi, double deps,
             double *epsa,
             double rb[3][3], double rp[3][3], double rbp[3][3],
             double rn[3][3], double rbpn[3][3]);
void eraPn00a(double date1, double date2,
              double *dpsi, double *deps, double *epsa,
              double rb[3][3], double rp[3][3], double rbp[3][3],
              double rn[3][3], double rbpn[3][3]);
void eraPn00b(double date1, double date2,
              double *dpsi, double *deps, double *epsa,
              double rb[3][3], double rp[3][3], double rbp[3][3],
              double rn[3][3], double rbpn[3][3]);
void eraPn06(double date1, double date2, double dpsi, double deps,
             double *epsa,
             double rb[3][3], double rp[3][3], double rbp[3][3],
             double rn[3][3], double rbpn[3][3]);
void eraPn06a(double date1, double date2,
              double *dpsi, double *deps, double *epsa,
              double rb[3][3], double rp[3][3], double rbp[3][3],
              double rn[3][3], double rbpn[3][3]);
void eraPnm00a(double date1, double date2, double rbpn[3][3]);
void eraPnm00b(double date1, double date2, double rbpn[3][3]);
void eraPnm06a(double date1, double date2, double rnpb[3][3]);
void eraPnm80(double date1, double date2, double rmatpn[3][3]);
void eraPom00(double xp, double yp, double sp, double rpom[3][3]);
void eraPr00(double date1, double date2, double *dpsipr, double *depspr);
void eraPrec76(double ep01, double ep02, double ep11, double ep12,
               double *zeta, double *z, double *theta);
double eraS00(double date1, double date2, double x, double y);
double eraS00a(double date1, double date2);
double eraS00b(double date1, double date2);
double eraS06(double date1, double date2, double x, double y);
double eraS06a(double date1, double date2);
double eraSp00(double date1, double date2);
void eraXy06(double date1, double date2, double *x, double *y);
void eraXys00a(double date1, double date2,
               double *x, double *y, double *s);
void eraXys00b(double date1, double date2,
               double *x, double *y, double *s);
void eraXys06a(double date1, double date2,
               double *x, double *y, double *s);

/* Astronomy/RotationAndTime */
double eraEe00(double date1, double date2, double epsa, double dpsi);
double eraEe00a(double date1, double date2);
double eraEe00b(double date1, double date2);
double eraEe06a(double date1, double date2);
double eraEect00(double date1, double date2);
double eraEqeq94(double date1, double date2);
double eraEra00(double dj1, double dj2);
double eraGmst00(double uta, double utb, double tta, double ttb);
double eraGmst06(double uta, double utb, double tta, double ttb);
double eraGmst82(double dj1, double dj2);
double eraGst00a(double uta, double utb, double tta, double ttb);
double eraGst00b(double uta, double utb);
double eraGst06(double uta, double utb, double tta, double ttb,
                double rnpb[3][3]);
double eraGst06a(double uta, double utb, double tta, double ttb);
double eraGst94(double uta, double utb);

/* Astronomy/SpaceMotion */
int eraPvstar(double pv[2][3], double *ra, double *dec,
              double *pmr, double *pmd, double *px, double *rv);
int eraStarpv(double ra, double dec,
              double pmr, double pmd, double px, double rv,
              double pv[2][3]);

/* Astronomy/StarCatalogs */
void eraFk52h(double r5, double d5,
              double dr5, double dd5, double px5, double rv5,
              double *rh, double *dh,
              double *drh, double *ddh, double *pxh, double *rvh);
void eraFk5hip(double r5h[3][3], double s5h[3]);
void eraFk5hz(double r5, double d5, double date1, double date2,
              double *rh, double *dh);
void eraH2fk5(double rh, double dh,
              double drh, double ddh, double pxh, double rvh,
              double *r5, double *d5,
              double *dr5, double *dd5, double *px5, double *rv5);
void eraHfk5z(double rh, double dh, double date1, double date2,
              double *r5, double *d5, double *dr5, double *dd5);
int eraStarpm(double ra1, double dec1,
              double pmr1, double pmd1, double px1, double rv1,
              double ep1a, double ep1b, double ep2a, double ep2b,
              double *ra2, double *dec2,
              double *pmr2, double *pmd2, double *px2, double *rv2);

/* Astronomy/Geodetic/Geocentric */
int eraEform(int n, double *a, double *f);
int eraGc2gd(int n, double xyz[3],
             double *elong, double *phi, double *height);
int eraGc2gde(double a, double f, double xyz[3],
              double *elong, double *phi, double *height);
int eraGd2gc(int n, double elong, double phi, double height,
             double xyz[3]);
int eraGd2gce(double a, double f,
              double elong, double phi, double height, double xyz[3]);

/* Astronomy/Timescales */
int eraD2dtf(const char *scale, int ndp, double d1, double d2,
             int *iy, int *im, int *id, int ihmsf[4]);
int eraDat(int iy, int im, int id, double fd, double *deltat);
double eraDtdb(double date1, double date2,
               double ut, double elong, double u, double v);
int eraDtf2d(const char *scale, int iy, int im, int id,
             int ihr, int imn, double sec, double *d1, double *d2);
int eraTaitt(double tai1, double tai2, double *tt1, double *tt2);
int eraTaiut1(double tai1, double tai2, double dta,
              double *ut11, double *ut12);
int eraTaiutc(double tai1, double tai2, double *utc1, double *utc2);
int eraTcbtdb(double tcb1, double tcb2, double *tdb1, double *tdb2);
int eraTcgtt(double tcg1, double tcg2, double *tt1, double *tt2);
int eraTdbtcb(double tdb1, double tdb2, double *tcb1, double *tcb2);
int eraTdbtt(double tdb1, double tdb2, double dtr,
             double *tt1, double *tt2);
int eraTttai(double tt1, double tt2, double *tai1, double *tai2);
int eraTttcg(double tt1, double tt2, double *tcg1, double *tcg2);
int eraTttdb(double tt1, double tt2, double dtr,
             double *tdb1, double *tdb2);
int eraTtut1(double tt1, double tt2, double dt,
             double *ut11, double *ut12);
int eraUt1tai(double ut11, double ut12, double dta,
              double *tai1, double *tai2);
int eraUt1tt(double ut11, double ut12, double dt,
             double *tt1, double *tt2);
int eraUt1utc(double ut11, double ut12, double dut1,
              double *utc1, double *utc2);
int eraUtctai(double utc1, double utc2, double *tai1, double *tai2);
int eraUtcut1(double utc1, double utc2, double dut1,
              double *ut11, double *ut12);

/* VectorMatrix/AngleOps */
void eraA2af(int ndp, double angle, char *sign, int idmsf[4]);
void eraA2tf(int ndp, double angle, char *sign, int ihmsf[4]);
int eraAf2a(char s, int ideg, int iamin, double asec, double *rad);
double eraAnp(double a);
double eraAnpm(double a);
void eraD2tf(int ndp, double days, char *sign, int ihmsf[4]);
int eraTf2a(char s, int ihour, int imin, double sec, double *rad);
int eraTf2d(char s, int ihour, int imin, double sec, double *days);

/* VectorMatrix/BuildRotations */
void eraRx(double phi, double r[3][3]);
void eraRy(double theta, double r[3][3]);
void eraRz(double psi, double r[3][3]);

/* VectorMatrix/CopyExtendExtract */
void eraCp(double p[3], double c[3]);
void eraCpv(double pv[2][3], double c[2][3]);
void eraCr(double r[3][3], double c[3][3]);
void eraP2pv(double p[3], double pv[2][3]);
void eraPv2p(double pv[2][3], double p[3]);

/* VectorMatrix/Initialization */
void eraIr(double r[3][3]);
void eraZp(double p[3]);
void eraZpv(double pv[2][3]);
void eraZr(double r[3][3]);

/* VectorMatrix/MatrixOps */
void eraRxr(double a[3][3], double b[3][3], double atb[3][3]);
void eraTr(double r[3][3], double rt[3][3]);

/* VectorMatrix/MatrixVectorProducts */
void eraRxp(double r[3][3], double p[3], double rp[3]);
void eraRxpv(double r[3][3], double pv[2][3], double rpv[2][3]);
void eraTrxp(double r[3][3], double p[3], double trp[3]);
void eraTrxpv(double r[3][3], double pv[2][3], double trpv[2][3]);

/* VectorMatrix/RotationVectors */
void eraRm2v(double r[3][3], double w[3]);
void eraRv2m(double w[3], double r[3][3]);

/* VectorMatrix/SeparationAndAngle */
double eraPap(double a[3], double b[3]);
double eraPas(double al, double ap, double bl, double bp);
double eraSepp(double a[3], double b[3]);
double eraSeps(double al, double ap, double bl, double bp);

/* VectorMatrix/SphericalCartesian */
void eraC2s(double p[3], double *theta, double *phi);
void eraP2s(double p[3], double *theta, double *phi, double *r);
void eraPv2s(double pv[2][3],
             double *theta, double *phi, double *r,
             double *td, double *pd, double *rd);
void eraS2c(double theta, double phi, double c[3]);
void eraS2p(double theta, double phi, double r, double p[3]);
void eraS2pv(double theta, double phi, double r,
             double td, double pd, double rd,
             double pv[2][3]);

/* VectorMatrix/VectorOps */
double eraPdp(double a[3], double b[3]);
double eraPm(double p[3]);
void eraPmp(double a[3], double b[3], double amb[3]);
void eraPn(double p[3], double *r, double u[3]);
void eraPpp(double a[3], double b[3], double apb[3]);
void eraPpsp(double a[3], double s, double b[3], double apsb[3]);
void eraPvdpv(double a[2][3], double b[2][3], double adb[2]);
void eraPvm(double pv[2][3], double *r, double *s);
void eraPvmpv(double a[2][3], double b[2][3], double amb[2][3]);
void eraPvppv(double a[2][3], double b[2][3], double apb[2][3]);
void eraPvu(double dt, double pv[2][3], double upv[2][3]);
void eraPvup(double dt, double pv[2][3], double p[3]);
void eraPvxpv(double a[2][3], double b[2][3], double axb[2][3]);
void eraPxp(double a[3], double b[3], double axb[3]);
void eraS2xpv(double s1, double s2, double pv[2][3], double spv[2][3]);
void eraSxp(double s, double p[3], double sp[3]);
void eraSxpv(double s, double pv[2][3], double spv[2][3]);



/* Pi */
#define DPI (3.141592653589793238462643)

/* 2Pi */
#define D2PI (6.283185307179586476925287)

/* Degrees to radians */
#define DD2R (1.745329251994329576923691e-2)

/* Radians to arcseconds */
#define DR2AS (206264.8062470963551564734)

/* Arcseconds to radians */
#define DAS2R (4.848136811095359935899141e-6)

/* Seconds of time to radians */
#define DS2R (7.272205216643039903848712e-5)

/* Arcseconds in a full circle */
#define TURNAS (1296000.0)

/* Milliarcseconds to radians */
#define DMAS2R (DAS2R / 1e3)

/* Length of tropical year B1900 (days) */
#define DTY (365.242198781)

/* Seconds per day. */
#define DAYSEC (86400.0)

/* Days per Julian year */
#define DJY (365.25)

/* Days per Julian century */
#define DJC (36525.0)

/* Days per Julian millennium */
#define DJM (365250.0)

/* Reference epoch (J2000.0), Julian Date */
#define DJ00 (2451545.0)

/* Julian Date of Modified Julian Date zero */
#define DJM0 (2400000.5)

/* Reference epoch (J2000.0), Modified Julian Date */
#define DJM00 (51544.5)

/* 1977 Jan 1.0 as MJD */
#define DJM77 (43144.0)

/* TT minus TAI (s) */
#define TTMTAI (32.184)

/* AU (m) */
#define DAU (149597870e3)

/* Speed of light (AU per day) */
#define DC (DAYSEC / 499.004782)

/* L_G = 1 - d(TT)/d(TCG) */
#define ELG (6.969290134e-10)

/* L_B = 1 - d(TDB)/d(TCB), and TDB (s) at TAI 1977/1/1.0 */
#define ELB (1.550519768e-8)
#define TDB0 (-6.55e-5)

/* dint(A) - truncate to nearest whole number towards zero (double) */
#define dint(A) ((A)<0.0?ceil(A):floor(A))

/* dnint(A) - round to nearest whole number (double) */
#define dnint(A) ((A)<0.0?ceil((A)-0.5):floor((A)+0.5))

/* dsign(A,B) - magnitude of A with sign of B (double) */
#define dsign(A,B) ((B)<0.0?-fabs(A):fabs(A))

/* Reference ellipsoids */
#define WGS84 1
#define GRS80 2
#define WGS72 3

#endif

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
