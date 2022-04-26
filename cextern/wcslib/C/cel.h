/*============================================================================
  WCSLIB 7.11 - an implementation of the FITS WCS standard.
  Copyright (C) 1995-2022, Mark Calabretta

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
  $Id: cel.h,v 7.11 2022/04/26 06:13:52 mcalabre Exp $
*=============================================================================
*
* WCSLIB 7.11 - C routines that implement the FITS World Coordinate System
* (WCS) standard.  Refer to the README file provided with WCSLIB for an
* overview of the library.
*
*
* Summary of the cel routines
* ---------------------------
* Routines in this suite implement the part of the FITS World Coordinate
* System (WCS) standard that deals with celestial coordinates, as described in
*
=   "Representations of world coordinates in FITS",
=   Greisen, E.W., & Calabretta, M.R. 2002, A&A, 395, 1061 (WCS Paper I)
=
=   "Representations of celestial coordinates in FITS",
=   Calabretta, M.R., & Greisen, E.W. 2002, A&A, 395, 1077 (WCS Paper II)
*
* These routines define methods to be used for computing celestial world
* coordinates from intermediate world coordinates (a linear transformation
* of image pixel coordinates), and vice versa.  They are based on the celprm
* struct which contains all information needed for the computations.  This
* struct contains some elements that must be set by the user, and others that
* are maintained by these routines, somewhat like a C++ class but with no
* encapsulation.
*
* Routine celini() is provided to initialize the celprm struct with default
* values, celfree() reclaims any memory that may have been allocated to store
* an error message, celsize() computes its total size including allocated
* memory, and celprt() prints its contents.
*
* celperr() prints the error message(s), if any, stored in a celprm struct and
* the prjprm struct that it contains.
*
* A setup routine, celset(), computes intermediate values in the celprm struct
* from parameters in it that were supplied by the user.  The struct always
* needs to be set up by celset() but it need not be called explicitly - refer
* to the explanation of celprm::flag.
*
* celx2s() and cels2x() implement the WCS celestial coordinate
* transformations.  In fact, they are high level driver routines for the lower
* level spherical coordinate rotation and projection routines described in
* sph.h and prj.h.
*
*
* celini() - Default constructor for the celprm struct
* ----------------------------------------------------
* celini() sets all members of a celprm struct to default values.  It should
* be used to initialize every celprm struct.
*
* PLEASE NOTE: If the celprm struct has already been initialized, then before
* reinitializing, it celfree() should be used to free any memory that may have
* been allocated to store an error message.  A memory leak may otherwise
* result.
*
* Returned:
*   cel       struct celprm*
*                       Celestial transformation parameters.
*
* Function return value:
*             int       Status return value:
*                         0: Success.
*                         1: Null celprm pointer passed.
*
*
* celfree() - Destructor for the celprm struct
* --------------------------------------------
* celfree() frees any memory that may have been allocated to store an error
* message in the celprm struct.
*
* Given:
*   cel       struct celprm*
*                       Celestial transformation parameters.
*
* Function return value:
*             int       Status return value:
*                         0: Success.
*                         1: Null celprm pointer passed.
*
*
* celsize() - Compute the size of a celprm struct
* -----------------------------------------------
* celsize() computes the full size of a celprm struct, including allocated
* memory.
*
* Given:
*   cel       const struct celprm*
*                       Celestial transformation parameters.
*
*                       If NULL, the base size of the struct and the allocated
*                       size are both set to zero.
*
* Returned:
*   sizes     int[2]    The first element is the base size of the struct as
*                       returned by sizeof(struct celprm).  The second element
*                       is the total allocated size, in bytes.  This figure
*                       includes memory allocated for the constituent struct,
*                       celprm::err.
*
*                       It is not an error for the struct not to have been set
*                       up via celset().
*
* Function return value:
*             int       Status return value:
*                         0: Success.
*
*
* celprt() - Print routine for the celprm struct
* ----------------------------------------------
* celprt() prints the contents of a celprm struct using wcsprintf().  Mainly
* intended for diagnostic purposes.
*
* Given:
*   cel       const struct celprm*
*                       Celestial transformation parameters.
*
* Function return value:
*             int       Status return value:
*                         0: Success.
*                         1: Null celprm pointer passed.
*
*
* celperr() - Print error messages from a celprm struct
* -----------------------------------------------------
* celperr() prints the error message(s), if any, stored in a celprm struct and
* the prjprm struct that it contains.  If there are no errors then nothing is
* printed.  It uses wcserr_prt(), q.v.
*
* Given:
*   cel       const struct celprm*
*                       Coordinate transformation parameters.
*
*   prefix    const char *
*                       If non-NULL, each output line will be prefixed with
*                       this string.
*
* Function return value:
*             int       Status return value:
*                         0: Success.
*                         1: Null celprm pointer passed.
*
*
* celset() - Setup routine for the celprm struct
* ----------------------------------------------
* celset() sets up a celprm struct according to information supplied within
* it.
*
* Note that this routine need not be called directly; it will be invoked by
* celx2s() and cels2x() if celprm::flag is anything other than a predefined
* magic value.
*
* Given and returned:
*   cel       struct celprm*
*                       Celestial transformation parameters.
*
* Function return value:
*             int       Status return value:
*                         0: Success.
*                         1: Null celprm pointer passed.
*                         2: Invalid projection parameters.
*                         3: Invalid coordinate transformation parameters.
*                         4: Ill-conditioned coordinate transformation
*                            parameters.
*
*                       For returns > 1, a detailed error message is set in
*                       celprm::err if enabled, see wcserr_enable().
*
*
* celx2s() - Pixel-to-world celestial transformation
* --------------------------------------------------
* celx2s() transforms (x,y) coordinates in the plane of projection to
* celestial coordinates (lng,lat).
*
* Given and returned:
*   cel       struct celprm*
*                       Celestial transformation parameters.
*
* Given:
*   nx,ny     int       Vector lengths.
*
*   sxy,sll   int       Vector strides.
*
*   x,y       const double[]
*                       Projected coordinates in pseudo "degrees".
*
* Returned:
*   phi,theta double[]  Longitude and latitude (phi,theta) in the native
*                       coordinate system of the projection [deg].
*
*   lng,lat   double[]  Celestial longitude and latitude (lng,lat) of the
*                       projected point [deg].
*
*   stat      int[]     Status return value for each vector element:
*                         0: Success.
*                         1: Invalid value of (x,y).
*
* Function return value:
*             int       Status return value:
*                         0: Success.
*                         1: Null celprm pointer passed.
*                         2: Invalid projection parameters.
*                         3: Invalid coordinate transformation parameters.
*                         4: Ill-conditioned coordinate transformation
*                            parameters.
*                         5: One or more of the (x,y) coordinates were
*                            invalid, as indicated by the stat vector.
*
*                       For returns > 1, a detailed error message is set in
*                       celprm::err if enabled, see wcserr_enable().
*
*
* cels2x() - World-to-pixel celestial transformation
* --------------------------------------------------
* cels2x() transforms celestial coordinates (lng,lat) to (x,y) coordinates in
* the plane of projection.
*
* Given and returned:
*   cel       struct celprm*
*                       Celestial transformation parameters.
*
* Given:
*   nlng,nlat int       Vector lengths.
*
*   sll,sxy   int       Vector strides.
*
*   lng,lat   const double[]
*                       Celestial longitude and latitude (lng,lat) of the
*                       projected point [deg].
*
* Returned:
*   phi,theta double[]  Longitude and latitude (phi,theta) in the native
*                       coordinate system of the projection [deg].
*
*   x,y       double[]  Projected coordinates in pseudo "degrees".
*
*   stat      int[]     Status return value for each vector element:
*                         0: Success.
*                         1: Invalid value of (lng,lat).
*
* Function return value:
*             int       Status return value:
*                         0: Success.
*                         1: Null celprm pointer passed.
*                         2: Invalid projection parameters.
*                         3: Invalid coordinate transformation parameters.
*                         4: Ill-conditioned coordinate transformation
*                            parameters.
*                         6: One or more of the (lng,lat) coordinates were
*                            invalid, as indicated by the stat vector.
*
*                       For returns > 1, a detailed error message is set in
*                       celprm::err if enabled, see wcserr_enable().
*
*
* celprm struct - Celestial transformation parameters
* ---------------------------------------------------
* The celprm struct contains information required to transform celestial
* coordinates.  It consists of certain members that must be set by the user
* ("given") and others that are set by the WCSLIB routines ("returned").  Some
* of the latter are supplied for informational purposes and others are for
* internal use only.
*
* Returned celprm struct members must not be modified by the user.
*
*   int flag
*     (Given and returned) This flag must be set to zero whenever any of the
*     following celprm struct members are set or changed:
*
*       - celprm::offset,
*       - celprm::phi0,
*       - celprm::theta0,
*       - celprm::ref[4],
*       - celprm::prj:
*         - prjprm::code,
*         - prjprm::r0,
*         - prjprm::pv[],
*         - prjprm::phi0,
*         - prjprm::theta0.
*
*     This signals the initialization routine, celset(), to recompute the
*     returned members of the celprm struct.  celset() will reset flag to
*     indicate that this has been done.
*
*   int offset
*     (Given) If true (non-zero), an offset will be applied to (x,y) to
*     force (x,y) = (0,0) at the fiducial point, (phi_0,theta_0).
*     Default is 0 (false).
*
*   double phi0
*     (Given) The native longitude, phi_0 [deg], and ...
*
*   double theta0
*     (Given) ... the native latitude, theta_0 [deg], of the fiducial point,
*     i.e. the point whose celestial coordinates are given in
*     celprm::ref[1:2].  If undefined (set to a magic value by prjini()) the
*     initialization routine, celset(), will set this to a projection-specific
*     default.
*
*   double ref[4]
*     (Given) The first pair of values should be set to the celestial
*     longitude and latitude of the fiducial point [deg] - typically right
*     ascension and declination.  These are given by the CRVALia keywords in
*     FITS.
*
*     (Given and returned) The second pair of values are the native longitude,
*     phi_p [deg], and latitude, theta_p [deg], of the celestial pole (the
*     latter is the same as the celestial latitude of the native pole,
*     delta_p) and these are given by the FITS keywords LONPOLEa and LATPOLEa
*     (or by PVi_2a and PVi_3a attached to the longitude axis which take
*     precedence if defined).
*
*     LONPOLEa defaults to phi_0 (see above) if the celestial latitude of the
*     fiducial point of the projection is greater than or equal to the native
*     latitude, otherwise phi_0 + 180 [deg].  (This is the condition for the
*     celestial latitude to increase in the same direction as the native
*     latitude at the fiducial point.)  ref[2] may be set to UNDEFINED (from
*     wcsmath.h) or 999.0 to indicate that the correct default should be
*     substituted.
*
*     theta_p, the native latitude of the celestial pole (or equally the
*     celestial latitude of the native pole, delta_p) is often determined
*     uniquely by CRVALia and LONPOLEa in which case LATPOLEa is ignored.
*     However, in some circumstances there are two valid solutions for theta_p
*     and LATPOLEa is used to choose between them.  LATPOLEa is set in ref[3]
*     and the solution closest to this value is used to reset ref[3].  It is
*     therefore legitimate, for example, to set ref[3] to +90.0 to choose the
*     more northerly solution - the default if the LATPOLEa keyword is omitted
*     from the FITS header.  For the special case where the fiducial point of
*     the projection is at native latitude zero, its celestial latitude is
*     zero, and LONPOLEa = +/- 90.0 then the celestial latitude of the native
*     pole is not determined by the first three reference values and LATPOLEa
*     specifies it completely.
*
*     The returned value, celprm::latpreq, specifies how LATPOLEa was actually
*     used.
*
*   struct prjprm prj
*     (Given and returned) Projection parameters described in the prologue to
*     prj.h.
*
*   double euler[5]
*     (Returned) Euler angles and associated intermediaries derived from the
*     coordinate reference values.  The first three values are the Z-, X-, and
*     Z'-Euler angles [deg], and the remaining two are the cosine and sine of
*     the X-Euler angle.
*
*   int latpreq
*     (Returned) For informational purposes, this indicates how the LATPOLEa
*     keyword was used
*       - 0: Not required, theta_p (== delta_p) was determined uniquely by the
*            CRVALia and LONPOLEa keywords.
*       - 1: Required to select between two valid solutions of theta_p.
*       - 2: theta_p was specified solely by LATPOLEa.
*
*   int isolat
*     (Returned) True if the spherical rotation preserves the magnitude of the
*     latitude, which occurs iff the axes of the native and celestial
*     coordinates are coincident.  It signals an opportunity to cache
*     intermediate calculations common to all elements in a vector
*     computation.
*
*   struct wcserr *err
*     (Returned) If enabled, when an error status is returned, this struct
*     contains detailed information about the error, see wcserr_enable().
*
*   void *padding
*     (An unused variable inserted for alignment purposes only.)
*
* Global variable: const char *cel_errmsg[] - Status return messages
* ------------------------------------------------------------------
* Status messages to match the status value returned from each function.
*
*===========================================================================*/

#ifndef WCSLIB_CEL
#define WCSLIB_CEL

#include "prj.h"

#ifdef __cplusplus
extern "C" {
#endif


extern const char *cel_errmsg[];

enum cel_errmsg_enum {
  CELERR_SUCCESS         = 0,	// Success.
  CELERR_NULL_POINTER    = 1,	// Null celprm pointer passed.
  CELERR_BAD_PARAM       = 2,	// Invalid projection parameters.
  CELERR_BAD_COORD_TRANS = 3,	// Invalid coordinate transformation
				// parameters.
  CELERR_ILL_COORD_TRANS = 4,	// Ill-conditioned coordinated transformation
				// parameters.
  CELERR_BAD_PIX         = 5,	// One or more of the (x,y) coordinates were
				// invalid.
  CELERR_BAD_WORLD       = 6 	// One or more of the (lng,lat) coordinates
				// were invalid.
};

struct celprm {
  // Initialization flag (see the prologue above).
  //--------------------------------------------------------------------------
  int    flag;			// Set to zero to force initialization.

  // Parameters to be provided (see the prologue above).
  //--------------------------------------------------------------------------
  int    offset;		// Force (x,y) = (0,0) at (phi_0,theta_0).
  double phi0, theta0;		// Native coordinates of fiducial point.
  double ref[4];		// Celestial coordinates of fiducial
                                // point and native coordinates of
                                // celestial pole.

  struct prjprm prj;		// Projection parameters (see prj.h).

  // Information derived from the parameters supplied.
  //--------------------------------------------------------------------------
  double euler[5];		// Euler angles and functions thereof.
  int    latpreq;		// LATPOLEa requirement.
  int    isolat;		// True if |latitude| is preserved.

  // Error handling
  //--------------------------------------------------------------------------
  struct wcserr *err;

  // Private
  //--------------------------------------------------------------------------
  void   *padding;		// (Dummy inserted for alignment purposes.)
};

// Size of the celprm struct in int units, used by the Fortran wrappers.
#define CELLEN (sizeof(struct celprm)/sizeof(int))


int celini(struct celprm *cel);

int celfree(struct celprm *cel);

int celsize(const struct celprm *cel, int sizes[2]);

int celprt(const struct celprm *cel);

int celperr(const struct celprm *cel, const char *prefix);

int celset(struct celprm *cel);

int celx2s(struct celprm *cel, int nx, int ny, int sxy, int sll,
           const double x[], const double y[],
           double phi[], double theta[], double lng[], double lat[],
           int stat[]);

int cels2x(struct celprm *cel, int nlng, int nlat, int sll, int sxy,
           const double lng[], const double lat[],
           double phi[], double theta[], double x[], double y[],
           int stat[]);


// Deprecated.
#define celini_errmsg cel_errmsg
#define celprt_errmsg cel_errmsg
#define celset_errmsg cel_errmsg
#define celx2s_errmsg cel_errmsg
#define cels2x_errmsg cel_errmsg

#ifdef __cplusplus
}
#endif

#endif // WCSLIB_CEL
