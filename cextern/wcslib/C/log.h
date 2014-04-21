/*============================================================================

  WCSLIB 4.23 - an implementation of the FITS WCS standard.
  Copyright (C) 1995-2014, Mark Calabretta

  This file is part of WCSLIB.

  WCSLIB is free software: you can redistribute it and/or modify it under the
  terms of the GNU Lesser General Public License as published by the Free
  Software Foundation, either version 3 of the License, or (at your option)
  any later version.

  WCSLIB is distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for
  more details.

  You should have received a copy of the GNU Lesser General Public License
  along with WCSLIB.  If not, see http://www.gnu.org/licenses.

  Direct correspondence concerning WCSLIB to mark@calabretta.id.au

  Author: Mark Calabretta, Australia Telescope National Facility, CSIRO.
  http://www.atnf.csiro.au/people/Mark.Calabretta
  $Id: log.h,v 4.23 2014/05/11 04:09:38 mcalabre Exp $
*=============================================================================
*
* WCSLIB 4.23 - C routines that implement logarithmic coordinate systems as
* defined by the FITS World Coordinate System (WCS) standard.  Refer to
*
*   "Representations of world coordinates in FITS",
*   Greisen, E.W., & Calabretta, M.R. 2002, A&A, 395, 1061 (Paper I)
*
*   "Representations of spectral coordinates in FITS",
*   Greisen, E.W., Calabretta, M.R., Valdes, F.G., & Allen, S.L.
*   2006, A&A, 446, 747 (Paper III)
*
* Refer to the README file provided with WCSLIB for an overview of the
* library.
*
*
* Summary of the log routines
* ---------------------------
* These routines implement the part of the FITS WCS standard that deals with
* logarithmic coordinates.  They define methods to be used for computing
* logarithmic world coordinates from intermediate world coordinates (a linear
* transformation of image pixel coordinates), and vice versa.
*
* logx2s() and logs2x() implement the WCS logarithmic coordinate
* transformations.
*
* Argument checking:
* ------------------
* The input log-coordinate values are only checked for values that would
* result in floating point exceptions and the same is true for the
* log-coordinate reference value.
*
* Accuracy:
* ---------
* No warranty is given for the accuracy of these routines (refer to the
* copyright notice); intending users must satisfy for themselves their
* adequacy for the intended purpose.  However, closure effectively to within
* double precision rounding error was demonstrated by test routine tlog.c
* which accompanies this software.
*
*
* logx2s() - Transform to logarithmic coordinates
* -----------------------------------------------
* logx2s() transforms intermediate world coordinates to logarithmic
* coordinates.
*
* Given and returned:
*   crval     double    Log-coordinate reference value (CRVALia).
*
* Given:
*   nx        int       Vector length.
*
*   sx        int       Vector stride.
*
*   slogc     int       Vector stride.
*
*   x         const double[]
*                       Intermediate world coordinates, in SI units.
*
* Returned:
*   logc      double[]  Logarithmic coordinates, in SI units.
*
*   stat      int[]     Status return value status for each vector element:
*                         0: Success.
*
* Function return value:
*             int       Status return value:
*                         0: Success.
*                         2: Invalid log-coordinate reference value.
*
*
* logs2x() - Transform logarithmic coordinates
* --------------------------------------------
* logs2x() transforms logarithmic world coordinates to intermediate world
* coordinates.
*
* Given and returned:
*   crval     double    Log-coordinate reference value (CRVALia).
*
* Given:
*   nlogc     int       Vector length.
*
*   slogc     int       Vector stride.
*
*   sx        int       Vector stride.
*
*   logc      const double[]
*                       Logarithmic coordinates, in SI units.
*
* Returned:
*   x         double[]  Intermediate world coordinates, in SI units.
*
*   stat      int[]     Status return value status for each vector element:
*                         0: Success.
*                         1: Invalid value of logc.
*
* Function return value:
*             int       Status return value:
*                         0: Success.
*                         2: Invalid log-coordinate reference value.
*                         4: One or more of the world-coordinate values
*                            are incorrect, as indicated by the stat vector.
*
*
* Global variable: const char *log_errmsg[] - Status return messages
* ------------------------------------------------------------------
* Error messages to match the status value returned from each function.
*
*===========================================================================*/

#ifndef WCSLIB_LOG
#define WCSLIB_LOG

#ifdef __cplusplus
extern "C" {
#endif

extern const char *log_errmsg[];

enum log_errmsg_enum {
  LOGERR_SUCCESS         = 0,	/* Success. */
  LOGERR_NULL_POINTER    = 1,	/* Null pointer passed. */
  LOGERR_BAD_LOG_REF_VAL = 2,	/* Invalid log-coordinate reference value. */
  LOGERR_BAD_X           = 3,	/* One or more of the x coordinates were
				   invalid. */
  LOGERR_BAD_WORLD       = 4 	/* One or more of the world coordinates were
				   invalid. */
};

int logx2s(double crval, int nx, int sx, int slogc, const double x[],
           double logc[], int stat[]);

int logs2x(double crval, int nlogc, int slogc, int sx, const double logc[],
           double x[], int stat[]);


#ifdef __cplusplus
}
#endif

#endif /* WCSLIB_LOG */
