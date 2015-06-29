/*============================================================================

  WCSLIB 5.7 - an implementation of the FITS WCS standard.
  Copyright (C) 1995-2015, Mark Calabretta

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
  $Id: wcsfix.h,v 5.7 2015/06/29 02:44:16 mcalabre Exp $
*=============================================================================
*
* WCSLIB 5.7 - C routines that implement the FITS World Coordinate System
* (WCS) standard.  Refer to the README file provided with WCSLIB for an
* overview of the library.
*
*
* Summary of the wcsfix routines
* ------------------------------
* Routines in this suite identify and translate various forms of construct
* known to occur in FITS headers that violate the FITS World Coordinate System
* (WCS) standard described in
*
=   "Representations of world coordinates in FITS",
=   Greisen, E.W., & Calabretta, M.R. 2002, A&A, 395, 1061 (WCS Paper I)
=
=   "Representations of celestial coordinates in FITS",
=   Calabretta, M.R., & Greisen, E.W. 2002, A&A, 395, 1077 (WCS Paper II)
=
=   "Representations of spectral coordinates in FITS",
=   Greisen, E.W., Calabretta, M.R., Valdes, F.G., & Allen, S.L.
=   2006, A&A, 446, 747 (WCS Paper III)
*
* Repairs effected by these routines range from the translation of
* non-standard values for standard WCS keywords, to the repair of malformed
* coordinate representations.
*
* Non-standard keyvalues:
* -----------------------
*   AIPS-convention celestial projection types, NCP and GLS, and spectral
*   types, 'FREQ-LSR', 'FELO-HEL', etc., set in CTYPEia are translated
*   on-the-fly by wcsset() but without modifying the relevant ctype[], pv[] or
*   specsys members of the wcsprm struct.  That is, only the information
*   extracted from ctype[] is translated when wcsset() fills in wcsprm::cel
*   (celprm struct) or wcsprm::spc (spcprm struct).
*
*   On the other hand, these routines do change the values of wcsprm::ctype[],
*   wcsprm::pv[], wcsprm::specsys and other wcsprm struct members as
*   appropriate to produce the same result as if the FITS header itself had
*   been translated.
*
*   Auxiliary WCS header information not used directly by WCSLIB may also be
*   translated.  For example, the older DATE-OBS date format (wcsprm::dateobs)
*   is recast to year-2000 standard form, and MJD-OBS (wcsprm::mjdobs) will be
*   deduced from it if not already set.
*
*   Certain combinations of keyvalues that result in malformed coordinate
*   systems, as described in Sect. 7.3.4 of Paper I, may also be repaired.
*   These are handled by cylfix().
*
* Non-standard keywords:
* ----------------------
*   The AIPS-convention CROTAn keywords are recognized as quasi-standard and
*   as such are accomodated by the wcsprm::crota[] and translated to
*   wcsprm::pc[][] by wcsset().  These are not dealt with here, nor are any
*   other non-standard keywords since these routines work only on the contents
*   of a wcsprm struct and do not deal with FITS headers per se.  In
*   particular, they do not identify or translate CD00i00j, PC00i00j, PROJPn,
*   EPOCH, VELREF or VSOURCEa keywords; this may be done by the FITS WCS
*   header parser supplied with WCSLIB, refer to wcshdr.h.
*
* wcsfix() and wcsfixi() apply all of the corrections handled by the following
* specific functions which may also be invoked separately:
*
*   - cdfix(): Sets the diagonal element of the CDi_ja matrix to 1.0 if all
*     CDi_ja keywords associated with a particular axis are omitted.
*
*   - datfix(): recast an older DATE-OBS date format in dateobs to year-2000
*     standard form and derive mjdobs from it if not already set.
*     Alternatively, if mjdobs is set and dateobs isn't, then derive dateobs
*     from it.
*
*   - unitfix(): translate some commonly used but non-standard unit strings in
*     the CUNITia keyvalues, e.g. 'DEG' -> 'deg'.
*
*   - spcfix(): translate AIPS-convention spectral types, 'FREQ-LSR',
*     'FELO-HEL', etc., in ctype[] as set from CTYPEia.
*
*   - celfix(): translate AIPS-convention celestial projection types, NCP and
*     GLS, in ctype[] as set from CTYPEia.
*
*   - cylfix(): fixes WCS keyvalues for malformed cylindrical projections that
*     suffer from the problem described in Sect. 7.3.4 of Paper I.
*
*
* wcsfix() - Translate a non-standard WCS struct
* ----------------------------------------------
* wcsfix() is identical to wcsfixi(), but lacks the info argument.
*
*
* wcsfixi() - Translate a non-standard WCS struct
* -----------------------------------------------
* wcsfix() applies all of the corrections handled separately by cdfix(),
* datfix(), unitfix(), spcfix(), celfix(), and cylfix().
*
* Given:
*   ctrl      int       Do potentially unsafe translations of non-standard
*                       unit strings as described in the usage notes to
*                       wcsutrn().
*
*   naxis     const int []
*                       Image axis lengths.  If this array pointer is set to
*                       zero then cylfix() will not be invoked.
*
* Given and returned:
*   wcs       struct wcsprm*
*                       Coordinate transformation parameters.
*
* Returned:
*   stat      int [NWCSFIX]
*                       Status returns from each of the functions.  Use the
*                       preprocessor macros NWCSFIX to dimension this vector
*                       and CDFIX, DATFIX, UNITFIX, SPCFIX, CELFIX, and CYLFIX
*                       to access its elements.  A status value of -2 is set
*                       for functions that were not invoked.
*
*   info      struct wcserr [NWCSFIX]
*                       Status messages from each of the functions.  Use the
*                       preprocessor macros NWCSFIX to dimension this vector
*                       and CDFIX, DATFIX, UNITFIX, SPCFIX, CELFIX, and CYLFIX
*                       to access its elements.
*
* Function return value:
*             int       Status return value:
*                         0: Success.
*                         1: One or more of the translation functions
*                            returned an error.
*
*
* cdfix() - Fix erroneously omitted CDi_ja keywords
* -------------------------------------------------
* cdfix() sets the diagonal element of the CDi_ja matrix to unity if all
* CDi_ja keywords associated with a given axis were omitted.  According to
* Paper I, if any CDi_ja keywords at all are given in a FITS header then those
* not given default to zero.  This results in a singular matrix with an
* intersecting row and column of zeros.
*
* Given and returned:
*   wcs       struct wcsprm*
*                       Coordinate transformation parameters.
*
* Function return value:
*             int       Status return value:
*                        -1: No change required (not an error).
*                         0: Success.
*                         1: Null wcsprm pointer passed.
*
*
* datfix() - Translate DATE-OBS and derive MJD-OBS or vice versa
* --------------------------------------------------------------
* datfix() translates the old DATE-OBS date format set in wcsprm::dateobs to
* year-2000 standard form (yyyy-mm-ddThh:mm:ss) and derives MJD-OBS from it if
* not already set.  Alternatively, if wcsprm::mjdobs is set and
* wcsprm::dateobs isn't, then datfix() derives wcsprm::dateobs from it.  If
* both are set but disagree by more than half a day then status 5 is returned.
*
* Given and returned:
*   wcs       struct wcsprm*
*                       Coordinate transformation parameters.  wcsprm::dateobs
*                       and/or wcsprm::mjdobs may be changed.
*
* Function return value:
*             int       Status return value:
*                        -1: No change required (not an error).
*                         0: Success.
*                         1: Null wcsprm pointer passed.
*                         5: Invalid parameter value.
*
*                       For returns > 1, a detailed error message is set in
*                       wcsprm::err if enabled, see wcserr_enable().
*
* Notes:
*   The MJD algorithms used by datfix() are from D.A. Hatcher, 1984, QJRAS,
*   25, 53-55, as modified by P.T. Wallace for use in SLALIB subroutines CLDJ
*   and DJCL.
*
*
* unitfix() - Correct aberrant CUNITia keyvalues
* ----------------------------------------------
* unitfix() applies wcsutrn() to translate non-standard CUNITia keyvalues,
* e.g. 'DEG' -> 'deg', also stripping off unnecessary whitespace.
*
* Given:
*   ctrl      int       Do potentially unsafe translations described in the
*                       usage notes to wcsutrn().
*
* Given and returned:
*   wcs       struct wcsprm*
*                       Coordinate transformation parameters.
*
* Function return value:
*             int       Status return value:
*                        -1: No change required (not an error).
*                         0: Success (an alias was applied).
*                         1: Null wcsprm pointer passed.
*
*                       When units are translated (i.e. status 0), status -2
*                       is set in the wcserr struct to allow an informative
*                       message to be returned.
*
*
* spcfix() - Translate AIPS-convention spectral types
* ---------------------------------------------------
* spcfix() translates AIPS-convention spectral coordinate types,
* '{FREQ,FELO,VELO}-{LSR,HEL,OBS}' (e.g. 'FREQ-OBS', 'FELO-HEL', 'VELO-LSR')
* set in wcsprm::ctype[], subject to VELREF set in wcsprm::velref.
*
* Note that if wcs::specsys is already set then it will not be overridden.
*
* Given and returned:
*   wcs       struct wcsprm*
*                       Coordinate transformation parameters.  wcsprm::ctype[]
*                       and/or wcsprm::specsys may be changed.
*
* Function return value:
*             int       Status return value:
*                        -1: No change required (not an error).
*                         0: Success.
*                         1: Null wcsprm pointer passed.
*                         2: Memory allocation failed.
*                         3: Linear transformation matrix is singular.
*                         4: Inconsistent or unrecognized coordinate axis
*                            types.
*                         5: Invalid parameter value.
*                         6: Invalid coordinate transformation parameters.
*                         7: Ill-conditioned coordinate transformation
*                            parameters.
*
*                       For returns > 1, a detailed error message is set in
*                       wcsprm::err if enabled, see wcserr_enable().
*
*
* celfix() - Translate AIPS-convention celestial projection types
* ---------------------------------------------------------------
* celfix() translates AIPS-convention celestial projection types, NCP and
* GLS, set in the ctype[] member of the wcsprm struct.
*
* Two additional pv[] keyvalues are created when translating NCP, and three
* are created when translating GLS with non-zero reference point.  If the pv[]
* array was initially allocated by wcsini() then the array will be expanded if
* necessary.  Otherwise, error 2 will be returned if sufficient empty slots
* are not already available for use.
*
* Given and returned:
*   wcs       struct wcsprm*
*                       Coordinate transformation parameters.  wcsprm::ctype[]
*                       and/or wcsprm::pv[] may be changed.
*
* Function return value:
*             int       Status return value:
*                        -1: No change required (not an error).
*                         0: Success.
*                         1: Null wcsprm pointer passed.
*                         2: Memory allocation failed.
*                         3: Linear transformation matrix is singular.
*                         4: Inconsistent or unrecognized coordinate axis
*                            types.
*                         5: Invalid parameter value.
*                         6: Invalid coordinate transformation parameters.
*                         7: Ill-conditioned coordinate transformation
*                            parameters.
*
*                       For returns > 1, a detailed error message is set in
*                       wcsprm::err if enabled, see wcserr_enable().
*
*
* cylfix() - Fix malformed cylindrical projections
* ------------------------------------------------
* cylfix() fixes WCS keyvalues for malformed cylindrical projections that
* suffer from the problem described in Sect. 7.3.4 of Paper I.
*
* Given:
*   naxis     const int []
*                       Image axis lengths.
*
* Given and returned:
*   wcs       struct wcsprm*
*                       Coordinate transformation parameters.
*
* Function return value:
*             int       Status return value:
*                        -1: No change required (not an error).
*                         0: Success.
*                         1: Null wcsprm pointer passed.
*                         2: Memory allocation failed.
*                         3: Linear transformation matrix is singular.
*                         4: Inconsistent or unrecognized coordinate axis
*                            types.
*                         5: Invalid parameter value.
*                         6: Invalid coordinate transformation parameters.
*                         7: Ill-conditioned coordinate transformation
*                            parameters.
*                         8: All of the corner pixel coordinates are invalid.
*                         9: Could not determine reference pixel coordinate.
*                        10: Could not determine reference pixel value.
*
*                       For returns > 1, a detailed error message is set in
*                       wcsprm::err if enabled, see wcserr_enable().
*
*
* Global variable: const char *wcsfix_errmsg[] - Status return messages
* ---------------------------------------------------------------------
* Error messages to match the status value returned from each function.
*
*===========================================================================*/

#ifndef WCSLIB_WCSFIX
#define WCSLIB_WCSFIX

#include "wcs.h"
#include "wcserr.h"

#ifdef __cplusplus
extern "C" {
#endif

#define CDFIX    0
#define DATFIX   1
#define UNITFIX  2
#define SPCFIX   3
#define CELFIX   4
#define CYLFIX   5
#define NWCSFIX  6

extern const char *wcsfix_errmsg[];
#define cylfix_errmsg wcsfix_errmsg

enum wcsfix_errmsg_enum {
  FIXERR_DATE_FIX         = -4, /* The date formatting has been fixed up. */
  FIXERR_SPC_UPDATE       = -3, /* Spectral axis type modified. */
  FIXERR_UNITS_ALIAS      = -2,	/* Units alias translation. */
  FIXERR_NO_CHANGE        = -1,	/* No change. */
  FIXERR_SUCCESS          =  0,	/* Success. */
  FIXERR_NULL_POINTER     =  1,	/* Null wcsprm pointer passed. */
  FIXERR_MEMORY           =  2,	/* Memory allocation failed. */
  FIXERR_SINGULAR_MTX     =  3,	/* Linear transformation matrix is
				   singular. */
  FIXERR_BAD_CTYPE        =  4,	/* Inconsistent or unrecognized coordinate
				   axis types. */
  FIXERR_BAD_PARAM        =  5,	/* Invalid parameter value. */
  FIXERR_BAD_COORD_TRANS  =  6,	/* Invalid coordinate transformation
				   parameters. */
  FIXERR_ILL_COORD_TRANS  =  7,	/* Ill-conditioned coordinate transformation
				   parameters. */
  FIXERR_BAD_CORNER_PIX   =  8,	/* All of the corner pixel coordinates are
				   invalid. */
  FIXERR_NO_REF_PIX_COORD =  9,	/* Could not determine reference pixel
				   coordinate. */
  FIXERR_NO_REF_PIX_VAL   = 10	/* Could not determine reference pixel
				   value. */
};

int wcsfix(int ctrl, const int naxis[], struct wcsprm *wcs, int stat[]);

int wcsfixi(int ctrl, const int naxis[], struct wcsprm *wcs, int stat[],
            struct wcserr info[]);

int cdfix(struct wcsprm *wcs);

int datfix(struct wcsprm *wcs);

int unitfix(int ctrl, struct wcsprm *wcs);

int spcfix(struct wcsprm *wcs);

int celfix(struct wcsprm *wcs);

int cylfix(const int naxis[], struct wcsprm *wcs);


#ifdef __cplusplus
}
#endif

#endif /* WCSLIB_WCSFIX */
