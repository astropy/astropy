#include "erfa.h"

void eraAticq(double ri, double di, eraASTROM *astrom,
              double *rc, double *dc)
/*
**  - - - - - - - - -
**   e r a A t i c q
**  - - - - - - - - -
**
**  Quick CIRS RA,Dec to ICRS astrometric place, given the star-
**  independent astrometry parameters.
**
**  Use of this function is appropriate when efficiency is important and
**  where many star positions are all to be transformed for one date.
**  The star-independent astrometry parameters can be obtained by
**  calling one of the functions eraApci[13], eraApcg[13], eraApco[13]
**  or eraApcs[13].
**
**  Given:
**     ri,di  double     CIRS RA,Dec (radians)
**     astrom eraASTROM* star-independent astrometry parameters:
**      pmt    double       PM time interval (SSB, Julian years)
**      eb     double[3]    SSB to observer (vector, au)
**      eh     double[3]    Sun to observer (unit vector)
**      em     double       distance from Sun to observer (au)
**      v      double[3]    barycentric observer velocity (vector, c)
**      bm1    double       sqrt(1-|v|^2): reciprocal of Lorenz factor
**      bpn    double[3][3] bias-precession-nutation matrix
**      along  double       longitude + s' (radians)
**      xpl    double       polar motion xp wrt local meridian (radians)
**      ypl    double       polar motion yp wrt local meridian (radians)
**      sphi   double       sine of geodetic latitude
**      cphi   double       cosine of geodetic latitude
**      diurab double       magnitude of diurnal aberration vector
**      eral   double       "local" Earth rotation angle (radians)
**      refa   double       refraction constant A (radians)
**      refb   double       refraction constant B (radians)
**
**  Returned:
**     rc,dc  double     ICRS astrometric RA,Dec (radians)
**
**  Notes:
**
**  1) Only the Sun is taken into account in the light deflection
**     correction.
**
**  2) Iterative techniques are used for the aberration and light
**     deflection corrections so that the functions eraAtic13 (or
**     eraAticq) and eraAtci13 (or eraAtciq) are accurate inverses;
**     even at the edge of the Sun's disk the discrepancy is only about
**     1 nanoarcsecond.
**
**  Called:
**     eraS2c       spherical coordinates to unit vector
**     eraTrxp      product of transpose of r-matrix and p-vector
**     eraZp        zero p-vector
**     eraAb        stellar aberration
**     eraLdsun     light deflection by the Sun
**     eraC2s       p-vector to spherical
**     eraAnp       normalize angle into range +/- pi
**
**  Copyright (C) 2013-2015, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   int j, i;
   double pi[3], ppr[3], pnat[3], pco[3], w, d[3], before[3], r2, r,
          after[3];

/* CIRS RA,Dec to Cartesian. */
   eraS2c(ri, di, pi);

/* Bias-precession-nutation, giving GCRS proper direction. */
   eraTrxp(astrom->bpn, pi, ppr);

/* Aberration, giving GCRS natural direction. */
   eraZp(d);
   for (j = 0; j < 2; j++) {
      r2 = 0.0;
      for (i = 0; i < 3; i++) {
         w = ppr[i] - d[i];
         before[i] = w;
         r2 += w*w;
      }
      r = sqrt(r2);
      for (i = 0; i < 3; i++) {
         before[i] /= r;
      }
      eraAb(before, astrom->v, astrom->em, astrom->bm1, after);
      r2 = 0.0;
      for (i = 0; i < 3; i++) {
         d[i] = after[i] - before[i];
         w = ppr[i] - d[i];
         pnat[i] = w;
         r2 += w*w;
      }
      r = sqrt(r2);
      for (i = 0; i < 3; i++) {
         pnat[i] /= r;
      }
   }

/* Light deflection by the Sun, giving BCRS coordinate direction. */
   eraZp(d);
   for (j = 0; j < 5; j++) {
      r2 = 0.0;
      for (i = 0; i < 3; i++) {
         w = pnat[i] - d[i];
         before[i] = w;
         r2 += w*w;
      }
      r = sqrt(r2);
      for (i = 0; i < 3; i++) {
         before[i] /= r;
      }
      eraLdsun(before, astrom->eh, astrom->em, after);
      r2 = 0.0;
      for (i = 0; i < 3; i++) {
         d[i] = after[i] - before[i];
         w = pnat[i] - d[i];
         pco[i] = w;
         r2 += w*w;
      }
      r = sqrt(r2);
      for (i = 0; i < 3; i++) {
         pco[i] /= r;
      }
   }

/* ICRS astrometric RA,Dec. */
   eraC2s(pco, &w, dc);
   *rc = eraAnp(w);

/* Finished. */

}
/*----------------------------------------------------------------------
**  
**  
**  Copyright (C) 2013-2015, NumFOCUS Foundation.
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
