/*============================================================================
  WCSLIB 8.3 - an implementation of the FITS WCS standard.
  Copyright (C) 1995-2024, Mark Calabretta

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

  Author: Mark Calabretta, Australia Telescope National Facility, CSIRO.
  http://www.atnf.csiro.au/people/Mark.Calabretta
  $Id: wcsfix.h,v 8.3 2024/05/13 16:33:00 mcalabre Exp $
*=============================================================================
*
* WCSLIB 8.3 - C routines that implement the FITS World Coordinate System
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
=
=   "Representations of time coordinates in FITS -
=    Time and relative dimension in space",
=   Rots, A.H., Bunclark, P.S., Calabretta, M.R., Allen, S.L.,
=   Manchester, R.N., & Thompson, W.T. 2015, A&A, 574, A36 (WCS Paper VII)
*
* Repairs effected by these routines range from the translation of
* non-standard values for standard WCS keywords, to the repair of malformed
* coordinate representations.  Some routines are also provided to check the
* consistency of pairs of keyvalues that define the same measure in two
* different ways, for example, as a date and an MJD.
*
* A separate routine, wcspcx(), "regularizes" the linear transformation matrix
* component (PCi_j) of the coordinate transformation to make it more human-
* readable.  Where a coordinate description was constructed from CDi_j, it
* decomposes it into PCi_j + CDELTi in a meaningful way.  Optionally, it can
* also diagonalize the PCi_j matrix (as far as possible), i.e. undo a
* transposition of axes in the intermediate pixel coordinate system.
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
*   The AIPS-convention CROTAn keywords are recognized as quasi-standard
*   and as such are accomodated by wcsprm::crota[] and translated to
*   wcsprm::pc[][] by wcsset().  These are not dealt with here, nor are any
*   other non-standard keywords since these routines work only on the contents
*   of a wcsprm struct and do not deal with FITS headers per se.  In
*   particular, they do not identify or translate CD00i00j, PC00i00j, PROJPn,
*   EPOCH, VELREF or VSOURCEa keywords; this may be done by the FITS WCS
*   header parser supplied with WCSLIB, refer to wcshdr.h.
*
* wcsfix() and wcsfixi() apply all of the corrections handled by the following
* specific functions, which may also be invoked separately:
*
*   - cdfix(): Sets the diagonal element of the CDi_ja matrix to 1.0 if all
*     CDi_ja keywords associated with a particular axis are omitted.
*
*   - datfix(): recast an older DATE-OBS date format in dateobs to year-2000
*     standard form.  Derive dateref from mjdref if not already set.
*     Alternatively, if dateref is set and mjdref isn't, then derive mjdref
*     from it.  If both are set, then check consistency.  Likewise for dateobs
*     and mjdobs; datebeg and mjdbeg; dateavg and mjdavg; and dateend and
*     mjdend.
*
*   - obsfix(): if only one half of obsgeo[] is set, then derive the other
*     half from it.  If both halves are set, then check consistency.
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
* wcsfixi() applies all of the corrections handled separately by cdfix(),
* datfix(), obsfix(), unitfix(), spcfix(), celfix(), and cylfix().
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
*                       and CDFIX, DATFIX, OBSFIX, UNITFIX, SPCFIX, CELFIX,
*                       and CYLFIX to access its elements.  A status value
*                       of -2 is set for functions that were not invoked.
*
*   info      struct wcserr [NWCSFIX]
*                       Status messages from each of the functions.  Use the
*                       preprocessor macros NWCSFIX to dimension this vector
*                       and CDFIX, DATFIX, OBSFIX, UNITFIX, SPCFIX, CELFIX,
*                       and CYLFIX to access its elements.
*
*                       Note that the memory allocated by wcsfixi() for the
*                       message in each wcserr struct (wcserr::msg, if
*                       non-zero) must be freed by the user.  See
*                       wcsdealloc().
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
* CDi_ja keywords associated with a given axis were omitted.  According to WCS
* Paper I, if any CDi_ja keywords at all are given in a FITS header then those
* not given default to zero.  This results in a singular matrix with an
* intersecting row and column of zeros.
*
* cdfix() is expected to be invoked before wcsset(), which will fail if these
* errors have not been corrected.
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
* year-2000 standard form (yyyy-mm-ddThh:mm:ss).  It derives wcsprm::dateref
* from wcsprm::mjdref if not already set.  Alternatively, if dateref is set
* and mjdref isn't, then it derives mjdref from it.  If both are set but
* disagree by more than 0.001 day (86.4 seconds) then an error status is
* returned.  Likewise for wcsprm::dateobs and wcsprm::mjdobs; wcsprm::datebeg
* and wcsprm::mjdbeg; wcsprm::dateavg and wcsprm::mjdavg; and wcsprm::dateend
* and wcsprm::mjdend.
*
* If neither dateobs nor mjdobs are set, but wcsprm::jepoch (primarily) or
* wcsprm::bepoch is, then both are derived from it.  If jepoch and/or bepoch
* are set but disagree with dateobs or mjdobs by more than 0.000002 year
* (63.2 seconds), an informative message is produced.
*
* The translations done by datfix() do not affect and are not affected by
* wcsset().
*
* Given and returned:
*   wcs       struct wcsprm*
*                       Coordinate transformation parameters.
*                       wcsprm::dateref and/or wcsprm::mjdref may be changed.
*                       wcsprm::dateobs and/or wcsprm::mjdobs may be changed.
*                       wcsprm::datebeg and/or wcsprm::mjdbeg may be changed.
*                       wcsprm::dateavg and/or wcsprm::mjdavg may be changed.
*                       wcsprm::dateend and/or wcsprm::mjdend may be changed.
*
* Function return value:
*             int       Status return value:
*                        -1: No change required (not an error).
*                         0: Success.
*                         1: Null wcsprm pointer passed.
*                         5: Invalid parameter value.
*
*                       For returns >= 0, a detailed message, whether
*                       informative or an error message, may be set in
*                       wcsprm::err if enabled, see wcserr_enable(), with
*                       wcsprm::err.status set to FIXERR_DATE_FIX.
*
* Notes:
*   1: The MJD algorithms used by datfix() are from D.A. Hatcher, 1984, QJRAS,
*      25, 53-55, as modified by P.T. Wallace for use in SLALIB subroutines
*      CLDJ and DJCL.
*
*
* obsfix() - complete the OBSGEO-[XYZLBH] vector of observatory coordinates
* -------------------------------------------------------------------------
* obsfix() completes the wcsprm::obsgeo vector of observatory coordinates.
* That is, if only the (x,y,z) Cartesian coordinate triplet or the (l,b,h)
* geodetic coordinate triplet are set, then it derives the other triplet from
* it.  If both triplets are set, then it checks for consistency at the level
* of 1 metre.
*
* The operations done by obsfix() do not affect and are not affected by
* wcsset().
*
* Given:
*   ctrl      int       Flag that controls behaviour if one triplet is
*                       defined and the other is only partially defined:
*                         0: Reset only the undefined elements of an
*                            incomplete coordinate triplet.
*                         1: Reset all elements of an incomplete triplet.
*                         2: Don't make any changes, check for consistency
*                            only.  Returns an error if either of the two
*                            triplets is incomplete.
*
* Given and returned:
*   wcs       struct wcsprm*
*                       Coordinate transformation parameters.
*                       wcsprm::obsgeo may be changed.
*
* Function return value:
*             int       Status return value:
*                        -1: No change required (not an error).
*                         0: Success.
*                         1: Null wcsprm pointer passed.
*                         5: Invalid parameter value.
*
*                       For returns >= 0, a detailed message, whether
*                       informative or an error message, may be set in
*                       wcsprm::err if enabled, see wcserr_enable(), with
*                       wcsprm::err.status set to FIXERR_OBS_FIX.
*
* Notes:
*   1: While the International Terrestrial Reference System (ITRS) is based
*      solely on Cartesian coordinates, it recommends the use of the GRS80
*      ellipsoid in converting to geodetic coordinates.  However, while WCS
*      Paper III recommends ITRS Cartesian coordinates, Paper VII prescribes
*      the use of the IAU(1976) ellipsoid for geodetic coordinates, and
*      consequently that is what is used here.
*
*   2: For reference, parameters of commonly used global reference ellipsoids:
*
=          a (m)          1/f                    Standard
=        ---------  -------------  --------------------------------
=        6378140    298.2577        IAU(1976)
=        6378137    298.257222101   GRS80
=        6378137    298.257223563   WGS84
=        6378136    298.257         IERS(1989)
=        6378136.6  298.25642       IERS(2003,2010), IAU(2009/2012)
*
*      where f = (a - b) / a is the flattening, and a and b are the semi-major
*      and semi-minor radii in metres.
*
*   3: The transformation from geodetic (lng,lat,hgt) to Cartesian (x,y,z) is
*
=        x = (n + hgt)*coslng*coslat,
=        y = (n + hgt)*sinlng*coslat,
=        z = (n*(1.0 - e^2) + hgt)*sinlat,
*
*      where the "prime vertical radius", n, is a function of latitude
*
=        n = a / sqrt(1 - (e*sinlat)^2),
*
*      and a, the equatorial radius, and e^2 = (2 - f)*f, the (first)
*      eccentricity of the ellipsoid, are constants.  obsfix() inverts these
*      iteratively by writing
*
=           x = rho*coslng*coslat,
=           y = rho*sinlng*coslat,
=        zeta = rho*sinlat,
*
*      where
*
=         rho = n + hgt,
=             = sqrt(x^2 + y^2 + zeta^2),
=        zeta = z / (1 - n*e^2/rho),
*
*      and iterating over the value of zeta.  Since e is small, a good first
*      approximation is given by zeta = z.
*
*
* unitfix() - Correct aberrant CUNITia keyvalues
* ----------------------------------------------
* unitfix() applies wcsutrn() to translate non-standard CUNITia keyvalues,
* e.g. 'DEG' -> 'deg', also stripping off unnecessary whitespace.
*
* unitfix() is expected to be invoked before wcsset(), which will fail if
* non-standard CUNITia keyvalues have not been translated.
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
*                       When units are translated (i.e. 0 is returned), an
*                       informative message is set in wcsprm::err if enabled,
*                       see wcserr_enable(), with wcsprm::err.status set to
*                       FIXERR_UNITS_ALIAS.
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
* AIPS-convention spectral types set in CTYPEia are translated on-the-fly by
* wcsset() but without modifying wcsprm::ctype[] or wcsprm::specsys.  That is,
* only the information extracted from wcsprm::ctype[] is translated when
* wcsset() fills in wcsprm::spc (spcprm struct).  spcfix() modifies
* wcsprm::ctype[] so that if the header is subsequently written out, e.g. by
* wcshdo(), then it will contain translated CTYPEia keyvalues.
*
* The operations done by spcfix() do not affect and are not affected by
* wcsset().
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
*                       For returns >= 0, a detailed message, whether
*                       informative or an error message, may be set in
*                       wcsprm::err if enabled, see wcserr_enable(), with
*                       wcsprm::err.status set to FIXERR_SPC_UPDTE.
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
* AIPS-convention celestial projection types set in CTYPEia are translated
* on-the-fly by wcsset() but without modifying wcsprm::ctype[], wcsprm::pv[],
* or wcsprm::npv.  That is, only the information extracted from
* wcsprm::ctype[] is translated when wcsset() fills in wcsprm::cel (celprm
* struct).  celfix() modifies wcsprm::ctype[], wcsprm::pv[], and wcsprm::npv
* so that if the header is subsequently written out, e.g. by wcshdo(), then it
* will contain translated CTYPEia keyvalues and the relevant PVi_ma.
*
* The operations done by celfix() do not affect and are not affected by
* wcsset().  However, it uses information in the wcsprm struct provided by
* wcsset(), and will invoke it if necessary.
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
* cylfix() requires the wcsprm struct to have been set up by wcsset(), and
* will invoke it if necessary.  After modification, the struct is reset on
* return with an explicit call to wcsset().
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
* wcspcx() - regularize PCi_j
* ---------------------------
* wcspcx() "regularizes" the linear transformation matrix component of the
* coordinate transformation (PCi_ja) to make it more human-readable.
*
* Normally, upon encountering a FITS header containing a CDi_ja matrix,
* wcsset() simply treats it as PCi_ja and sets CDELTia to unity.  However,
* wcspcx() decomposes CDi_ja into PCi_ja and CDELTia in such a way that
* CDELTia form meaningful scaling parameters.  In practice, the residual
* PCi_ja matrix will often then be orthogonal, i.e. unity, or describing a
* pure rotation, axis permutation, or reflection, or a combination thereof.
*
* The decomposition is based on normalizing the length in the transformed
* system (i.e. intermediate pixel coordinates) of the orthonormal basis
* vectors of the pixel coordinate system.  This deviates slightly from the
* prescription given by Eq. (4) of WCS Paper I, namely Sum(j=1,N)(PCi_ja)² = 1,
* in replacing the sum over j with the sum over i.  Consequently, the columns
* of PCi_ja will consist of unit vectors.  In practice, especially in cubes
* and higher dimensional images, at least some pairs of these unit vectors, if
* not all, will often be orthogonal or close to orthogonal.
*
* The sign of CDELTia is chosen to make the PCi_ja matrix as close to the,
* possibly permuted, unit matrix as possible, except that where the coordinate
* description contains a pair of celestial axes, the sign of CDELTia is set
* negative for the longitude axis and positive for the latitude axis.
*
* Optionally, rows of the PCi_ja matrix may also be permuted to diagonalize
* it as far as possible, thus undoing any transposition of axes in the
* intermediate pixel coordinate system.
*
* If the coordinate description contains a celestial plane, then the angle of
* rotation of each of the basis vectors associated with the celestial axes is
* returned.  For a pure rotation the two angles should be identical.  Any
* difference between them is a measure of axis skewness.
*
* The decomposition is not performed for axes involving a sequent distortion
* function that is defined in terms of CDi_ja, such as TPV, TNX, or ZPX, which
* always are.  The independent variables of the polynomial are therefore
* intermediate world coordinates rather than intermediate pixel coordinates.
* Because sequent distortions are always applied before CDELTia, if CDi_ja was
* translated to PCi_ja plus CDELTia, then the distortion would be altered
* unless the polynomial coefficients were also adjusted to account for the
* change of scale.
*
* wcspcx() requires the wcsprm struct to have been set up by wcsset(), and
* will invoke it if necessary.  The wcsprm struct is reset on return with an
* explicit call to wcsset().
*
* Given and returned:
*   wcs       struct wcsprm*
*                       Coordinate transformation parameters.
*
* Given:
*   dopc      int       If 1, then PCi_ja and CDELTia, as given, will be
*                       recomposed according to the above prescription.  If 0,
*                       the operation is restricted to decomposing CDi_ja.
*
*   permute   int       If 1, then after decomposition (or recomposition),
*                       permute rows of PCi_ja to make the axes of the
*                       intermediate pixel coordinate system match as closely
*                       as possible those of the pixel coordinates.  That is,
*                       make it as close to a diagonal matrix as possible.
*                       However, celestial axes are special in always being
*                       paired, with the longitude axis preceding the latitude
*                       axis.
*
*                       All WCS entities indexed by i, such as CTYPEia,
*                       CRVALia, CDELTia, etc., including coordinate lookup
*                       tables, will also be permuted as necessary to account
*                       for the change to PCi_ja.  This does not apply to
*                       CRPIXja, nor prior distortion functions.  These
*                       operate on pixel coordinates, which are not affected
*                       by the permutation.
*
* Returned:
*   rotn      double[2] Rotation angle [deg] of each basis vector associated
*                       with the celestial axes.  For a pure rotation the two
*                       angles should be identical.  Any difference between
*                       them is a measure of axis skewness.
*
*                       May be set to the NULL pointer if this information is
*                       not required.
*
* Function return value:
*             int       Status return value:
*                         0: Success.
*                         1: Null wcsprm pointer passed.
*                         2: Memory allocation failed.
*                         5: CDi_j matrix not used.
*                         6: Sequent distortion function present.
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
#define OBSFIX   2
#define UNITFIX  3
#define SPCFIX   4
#define CELFIX   5
#define CYLFIX   6
#define NWCSFIX  7

extern const char *wcsfix_errmsg[];
#define cylfix_errmsg wcsfix_errmsg

enum wcsfix_errmsg_enum {
  FIXERR_OBSGEO_FIX       = -5, // Observatory coordinates amended.
  FIXERR_DATE_FIX         = -4, // Date string reformatted.
  FIXERR_SPC_UPDATE       = -3, // Spectral axis type modified.
  FIXERR_UNITS_ALIAS      = -2,	// Units alias translation.
  FIXERR_NO_CHANGE        = -1,	// No change.
  FIXERR_SUCCESS          =  0,	// Success.
  FIXERR_NULL_POINTER     =  1,	// Null wcsprm pointer passed.
  FIXERR_MEMORY           =  2,	// Memory allocation failed.
  FIXERR_SINGULAR_MTX     =  3,	// Linear transformation matrix is singular.
  FIXERR_BAD_CTYPE        =  4,	// Inconsistent or unrecognized coordinate
				// axis types.
  FIXERR_BAD_PARAM        =  5,	// Invalid parameter value.
  FIXERR_BAD_COORD_TRANS  =  6,	// Invalid coordinate transformation
				// parameters.
  FIXERR_ILL_COORD_TRANS  =  7,	// Ill-conditioned coordinate transformation
				// parameters.
  FIXERR_BAD_CORNER_PIX   =  8,	// All of the corner pixel coordinates are
				// invalid.
  FIXERR_NO_REF_PIX_COORD =  9,	// Could not determine reference pixel
				// coordinate.
  FIXERR_NO_REF_PIX_VAL   = 10	// Could not determine reference pixel value.
};

int wcsfix(int ctrl, const int naxis[], struct wcsprm *wcs, int stat[]);

int wcsfixi(int ctrl, const int naxis[], struct wcsprm *wcs, int stat[],
            struct wcserr info[]);

int cdfix(struct wcsprm *wcs);

int datfix(struct wcsprm *wcs);

int obsfix(int ctrl, struct wcsprm *wcs);

int unitfix(int ctrl, struct wcsprm *wcs);

int spcfix(struct wcsprm *wcs);

int celfix(struct wcsprm *wcs);

int cylfix(const int naxis[], struct wcsprm *wcs);

int wcspcx(struct wcsprm *wcs, int dopc, int permute, double rotn[2]);


#ifdef __cplusplus
}
#endif

#endif // WCSLIB_WCSFIX
