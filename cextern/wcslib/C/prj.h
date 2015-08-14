/*============================================================================

  WCSLIB 5.9 - an implementation of the FITS WCS standard.
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
  $Id: prj.h,v 5.9 2015/07/21 09:20:01 mcalabre Exp $
*=============================================================================
*
* WCSLIB 5.9 - C routines that implement the FITS World Coordinate System
* (WCS) standard.  Refer to the README file provided with WCSLIB for an
* overview of the library.
*
*
* Summary of the prj routines
* ---------------------------
* Routines in this suite  implement the spherical map projections defined by
* the FITS World Coordinate System (WCS) standard, as described in
*
=   "Representations of world coordinates in FITS",
=   Greisen, E.W., & Calabretta, M.R. 2002, A&A, 395, 1061 (WCS Paper I)
=
=   "Representations of celestial coordinates in FITS",
=   Calabretta, M.R., & Greisen, E.W. 2002, A&A, 395, 1077 (WCS Paper II)
=
=   "Mapping on the HEALPix grid",
=   Calabretta, M.R., & Roukema, B.F. 2007, MNRAS, 381, 865 (WCS Paper V)
=
=   "Representing the 'Butterfly' Projection in FITS -- Projection Code XPH",
=   Calabretta, M.R., & Lowe, S.R. 2013, PASA, 30, e050 (WCS Paper VI)
*
* These routines are based on the prjprm struct which contains all information
* needed for the computations.  The struct contains some members that must be
* set by the user, and others that are maintained by these routines, somewhat
* like a C++ class but with no encapsulation.
*
* Routine prjini() is provided to initialize the prjprm struct with default
* values, prjfree() reclaims any memory that may have been allocated to store
* an error message, and prjprt() prints its contents.
*
* prjperr() prints the error message(s) (if any) stored in a prjprm struct.
* prjbchk() performs bounds checking on native spherical coordinates.
*
* Setup routines for each projection with names of the form ???set(), where
* "???" is the down-cased three-letter projection code, compute intermediate
* values in the prjprm struct from parameters in it that were supplied by the
* user.  The struct always needs to be set by the projection's setup routine
* but that need not be called explicitly - refer to the explanation of
* prjprm::flag.
*
* Each map projection is implemented via separate functions for the spherical
* projection, ???s2x(), and deprojection, ???x2s().
*
* A set of driver routines, prjset(), prjx2s(), and prjs2x(), provides a
* generic interface to the specific projection routines which they invoke
* via pointers-to-functions stored in the prjprm struct.
*
* In summary, the routines are:
*   - prjini()                Initialization routine for the prjprm struct.
*   - prjfree()               Reclaim memory allocated for error messages.
*   - prjprt()                Print the prjprm struct.
*   - prjperr()               Print error message (if any).
*   - prjbchk()               Bounds checking on native coordinates.
*
*   - prjset(), prjx2s(), prjs2x():   Generic driver routines
*
*   - azpset(), azpx2s(), azps2x():   AZP (zenithal/azimuthal perspective)
*   - szpset(), szpx2s(), szps2x():   SZP (slant zenithal perspective)
*   - tanset(), tanx2s(), tans2x():   TAN (gnomonic)
*   - stgset(), stgx2s(), stgs2x():   STG (stereographic)
*   - sinset(), sinx2s(), sins2x():   SIN (orthographic/synthesis)
*   - arcset(), arcx2s(), arcs2x():   ARC (zenithal/azimuthal equidistant)
*   - zpnset(), zpnx2s(), zpns2x():   ZPN (zenithal/azimuthal polynomial)
*   - zeaset(), zeax2s(), zeas2x():   ZEA (zenithal/azimuthal equal area)
*   - airset(), airx2s(), airs2x():   AIR (Airy)
*   - cypset(), cypx2s(), cyps2x():   CYP (cylindrical perspective)
*   - ceaset(), ceax2s(), ceas2x():   CEA (cylindrical equal area)
*   - carset(), carx2s(), cars2x():   CAR (Plate carree)
*   - merset(), merx2s(), mers2x():   MER (Mercator)
*   - sflset(), sflx2s(), sfls2x():   SFL (Sanson-Flamsteed)
*   - parset(), parx2s(), pars2x():   PAR (parabolic)
*   - molset(), molx2s(), mols2x():   MOL (Mollweide)
*   - aitset(), aitx2s(), aits2x():   AIT (Hammer-Aitoff)
*   - copset(), copx2s(), cops2x():   COP (conic perspective)
*   - coeset(), coex2s(), coes2x():   COE (conic equal area)
*   - codset(), codx2s(), cods2x():   COD (conic equidistant)
*   - cooset(), coox2s(), coos2x():   COO (conic orthomorphic)
*   - bonset(), bonx2s(), bons2x():   BON (Bonne)
*   - pcoset(), pcox2s(), pcos2x():   PCO (polyconic)
*   - tscset(), tscx2s(), tscs2x():   TSC (tangential spherical cube)
*   - cscset(), cscx2s(), cscs2x():   CSC (COBE spherical cube)
*   - qscset(), qscx2s(), qscs2x():   QSC (quadrilateralized spherical cube)
*   - hpxset(), hpxx2s(), hpxs2x():   HPX (HEALPix)
*   - xphset(), xphx2s(), xphs2x():   XPH (HEALPix polar, aka "butterfly")
*
* Argument checking (projection routines):
* ----------------------------------------
* The values of phi and theta (the native longitude and latitude) normally lie
* in the range [-180,180] for phi, and [-90,90] for theta.  However, all
* projection routines will accept any value of phi and will not normalize it.
*
* The projection routines do not explicitly check that theta lies within the
* range [-90,90].  They do check for any value of theta that produces an
* invalid argument to the projection equations (e.g. leading to division by
* zero).  The projection routines for AZP, SZP, TAN, SIN, ZPN, and COP also
* return error 2 if (phi,theta) corresponds to the overlapped (far) side of
* the projection but also return the corresponding value of (x,y).  This
* strict bounds checking may be relaxed at any time by setting
* prjprm::bounds%2 to 0 (rather than 1); the projections need not be
* reinitialized.
*
* Argument checking (deprojection routines):
* ------------------------------------------
* Error checking on the projected coordinates (x,y) is limited to that
* required to ascertain whether a solution exists.  Where a solution does
* exist, an optional check is made that the value of phi and theta obtained
* lie within the ranges [-180,180] for phi, and [-90,90] for theta.  This
* check, performed by prjbchk(), is enabled by default.  It may be disabled by
* setting prjprm::bounds%4 to 0 (rather than 1); the projections need not be
* reinitialized.
*
* Accuracy:
* ---------
* No warranty is given for the accuracy of these routines (refer to the
* copyright notice); intending users must satisfy for themselves their
* adequacy for the intended purpose.  However, closure to a precision of at
* least 1E-10 degree of longitude and latitude has been verified for typical
* projection parameters on the 1 degree graticule of native longitude and
* latitude (to within 5 degrees of any latitude where the projection may
* diverge).  Refer to the tprj1.c and tprj2.c test routines that accompany
* this software.
*
*
* prjini() - Default constructor for the prjprm struct
* ----------------------------------------------------
* prjini() sets all members of a prjprm struct to default values.  It should
* be used to initialize every prjprm struct.
*
* Returned:
*   prj       struct prjprm*
*                       Projection parameters.
*
* Function return value:
*             int       Status return value:
*                         0: Success.
*                         1: Null prjprm pointer passed.
*
*
* prjfree() - Destructor for the prjprm struct
* --------------------------------------------
* prjfree() frees any memory that may have been allocated to store an error
* message in the prjprm struct.
*
* Given:
*   prj       struct prjprm*
*                       Projection parameters.
*
* Function return value:
*             int       Status return value:
*                         0: Success.
*                         1: Null prjprm pointer passed.
*
*
* prjprt() - Print routine for the prjprm struct
* ----------------------------------------------
* prjprt() prints the contents of a prjprm struct using wcsprintf().  Mainly
* intended for diagnostic purposes.
*
* Given:
*   prj       const struct prjprm*
*                       Projection parameters.
*
* Function return value:
*             int       Status return value:
*                         0: Success.
*                         1: Null prjprm pointer passed.
*
*
* prjperr() - Print error messages from a prjprm struct
* -----------------------------------------------------
* prjperr() prints the error message(s) (if any) stored in a prjprm struct.
* If there are no errors then nothing is printed.  It uses wcserr_prt(), q.v.
*
* Given:
*   prj       const struct prjprm*
*                       Projection parameters.
*
*   prefix    const char *
*                       If non-NULL, each output line will be prefixed with
*                       this string.
*
* Function return value:
*             int       Status return value:
*                         0: Success.
*                         1: Null prjprm pointer passed.
*
*
* prjbchk() - Bounds checking on native coordinates
* -------------------------------------------------
* prjbchk() performs bounds checking on native spherical coordinates.  As
* returned by the deprojection (x2s) routines, native longitude is expected
* to lie in the closed interval [-180,180], with latitude in [-90,90].
*
* A tolerance may be specified to provide a small allowance for numerical
* imprecision.  Values that lie outside the allowed range by not more than
* the specified tolerance will be adjusted back into range.
*
* If prjprm::bounds&4 is set, as it is by prjini(), then prjbchk() will be
* invoked automatically by the Cartesian-to-spherical deprojection (x2s)
* routines with an appropriate tolerance set for each projection.
*
* Given:
*   tol       double    Tolerance for the bounds check [deg].
*
*   nphi,
*   ntheta    int       Vector lengths.
*
*   spt       int       Vector stride.
*
* Given and returned:
*   phi,theta double[]  Native longitude and latitude (phi,theta) [deg].
*
* Returned:
*   stat      int[]     Status value for each vector element:
*                         0: Valid value of (phi,theta).
*                         1: Invalid value.
*
* Function return value:
*             int       Status return value:
*                         0: Success.
*                         1: One or more of the (phi,theta) coordinates
*                            were, invalid, as indicated by the stat vector.
*
*
* prjset() - Generic setup routine for the prjprm struct
* ------------------------------------------------------
* prjset() sets up a prjprm struct according to information supplied within
* it.
*
* Note that this routine need not be called directly; it will be invoked by
* prjx2s() and prjs2x() if prj.flag is anything other than a predefined magic
* value.
*
* The one important distinction between prjset() and the setup routines for
* the specific projections is that the projection code must be defined in the
* prjprm struct in order for prjset() to identify the required projection.
* Once prjset() has initialized the prjprm struct, prjx2s() and prjs2x() use
* the pointers to the specific projection and deprojection routines contained
* therein.
*
* Given and returned:
*   prj       struct prjprm*
*                       Projection parameters.
*
* Function return value:
*             int       Status return value:
*                         0: Success.
*                         1: Null prjprm pointer passed.
*                         2: Invalid projection parameters.
*
*                       For returns > 1, a detailed error message is set in
*                       prjprm::err if enabled, see wcserr_enable().
*
*
* prjx2s() - Generic Cartesian-to-spherical deprojection
* ------------------------------------------------------
* Deproject Cartesian (x,y) coordinates in the plane of projection to native
* spherical coordinates (phi,theta).
*
* The projection is that specified by prjprm::code.
*
* Given and returned:
*   prj       struct prjprm*
*                       Projection parameters.
*
* Given:
*   nx,ny     int       Vector lengths.
*
*   sxy,spt   int       Vector strides.
*
*   x,y       const double[]
*                       Projected coordinates.
*
* Returned:
*   phi,theta double[]  Longitude and latitude (phi,theta) of the projected
*                       point in native spherical coordinates [deg].
*
*   stat      int[]     Status value for each vector element:
*                         0: Success.
*                         1: Invalid value of (x,y).
*
* Function return value:
*             int       Status return value:
*                         0: Success.
*                         1: Null prjprm pointer passed.
*                         2: Invalid projection parameters.
*                         3: One or more of the (x,y) coordinates were
*                            invalid, as indicated by the stat vector.
*
*                       For returns > 1, a detailed error message is set in
*                       prjprm::err if enabled, see wcserr_enable().
*
*
* prjs2x() - Generic spherical-to-Cartesian projection
* ----------------------------------------------------
* Project native spherical coordinates (phi,theta) to Cartesian (x,y)
* coordinates in the plane of projection.
*
* The projection is that specified by prjprm::code.
*
* Given and returned:
*   prj       struct prjprm*
*                       Projection parameters.
*
* Given:
*   nphi,
*   ntheta    int       Vector lengths.
*
*   spt,sxy   int       Vector strides.
*
*   phi,theta const double[]
*                       Longitude and latitude (phi,theta) of the projected
*                       point in native spherical coordinates [deg].
*
* Returned:
*   x,y       double[]  Projected coordinates.
*
*   stat      int[]     Status value for each vector element:
*                         0: Success.
*                         1: Invalid value of (phi,theta).
*
* Function return value:
*             int       Status return value:
*                         0: Success.
*                         1: Null prjprm pointer passed.
*                         2: Invalid projection parameters.
*                         4: One or more of the (phi,theta) coordinates
*                            were, invalid, as indicated by the stat vector.
*
*                       For returns > 1, a detailed error message is set in
*                       prjprm::err if enabled, see wcserr_enable().
*
*
* ???set() - Specific setup routines for the prjprm struct
* --------------------------------------------------------
* Set up a prjprm struct for a particular projection according to information
* supplied within it.
*
* Given and returned:
*   prj       struct prjprm*
*                       Projection parameters.
*
* Function return value:
*             int       Status return value:
*                         0: Success.
*                         1: Null prjprm pointer passed.
*                         2: Invalid projection parameters.
*
*                       For returns > 1, a detailed error message is set in
*                       prjprm::err if enabled, see wcserr_enable().
*
*
* ???x2s() - Specific Cartesian-to-spherical deprojection routines
* ----------------------------------------------------------------
* Transform (x,y) coordinates in the plane of projection to native spherical
* coordinates (phi,theta).
*
* Given and returned:
*   prj       struct prjprm*
*                       Projection parameters.
*
* Given:
*   nx,ny     int       Vector lengths.
*
*   sxy,spt   int       Vector strides.
*
*   x,y       const double[]
*                       Projected coordinates.
*
* Returned:
*   phi,theta double[]  Longitude and latitude of the projected point in
*                       native spherical coordinates [deg].
*
*   stat      int[]     Status value for each vector element:
*                         0: Success.
*                         1: Invalid value of (x,y).
*
* Function return value:
*             int       Status return value:
*                         0: Success.
*                         1: Null prjprm pointer passed.
*                         2: Invalid projection parameters.
*                         3: One or more of the (x,y) coordinates were
*                            invalid, as indicated by the stat vector.
*
*                       For returns > 1, a detailed error message is set in
*                       prjprm::err if enabled, see wcserr_enable().
*
*
* ???s2x() - Specific spherical-to-Cartesian projection routines
*---------------------------------------------------------------
* Transform native spherical coordinates (phi,theta) to (x,y) coordinates in
* the plane of projection.
*
* Given and returned:
*   prj       struct prjprm*
*                       Projection parameters.
*
* Given:
*   nphi,
*   ntheta    int       Vector lengths.
*
*   spt,sxy   int       Vector strides.
*
*   phi,theta const double[]
*                       Longitude and latitude of the projected point in
*                       native spherical coordinates [deg].
*
* Returned:
*   x,y       double[]  Projected coordinates.
*
*   stat      int[]     Status value for each vector element:
*                         0: Success.
*                         1: Invalid value of (phi,theta).
*
* Function return value:
*             int       Status return value:
*                         0: Success.
*                         1: Null prjprm pointer passed.
*                         2: Invalid projection parameters.
*                         4: One or more of the (phi,theta) coordinates
*                            were, invalid, as indicated by the stat vector.
*
*                       For returns > 1, a detailed error message is set in
*                       prjprm::err if enabled, see wcserr_enable().
*
*
* prjprm struct - Projection parameters
* -------------------------------------
* The prjprm struct contains all information needed to project or deproject
* native spherical coordinates.  It consists of certain members that must be
* set by the user ("given") and others that are set by the WCSLIB routines
* ("returned").  Some of the latter are supplied for informational purposes
* while others are for internal use only.
*
*   int flag
*     (Given and returned) This flag must be set to zero whenever any of the
*     following prjprm struct members are set or changed:
*
*       - prjprm::code,
*       - prjprm::r0,
*       - prjprm::pv[],
*       - prjprm::phi0,
*       - prjprm::theta0.
*
*     This signals the initialization routine (prjset() or ???set()) to
*     recompute the returned members of the prjprm struct.  flag will then be
*     reset to indicate that this has been done.
*
*     Note that flag need not be reset when prjprm::bounds is changed.
*
*   char code[4]
*     (Given) Three-letter projection code defined by the FITS standard.
*
*   double r0
*     (Given) The radius of the generating sphere for the projection, a linear
*     scaling parameter.  If this is zero, it will be reset to its default
*     value of 180/pi (the value for FITS WCS).
*
*   double pv[30]
*     (Given) Projection parameters.  These correspond to the PVi_ma keywords
*     in FITS, so pv[0] is PVi_0a, pv[1] is PVi_1a, etc., where i denotes the
*     latitude-like axis.  Many projections use pv[1] (PVi_1a), some also use
*     pv[2] (PVi_2a) and SZP uses pv[3] (PVi_3a).  ZPN is currently the only
*     projection that uses any of the others.
*
*     Usage of the pv[] array as it applies to each projection is described in
*     the prologue to each trio of projection routines in prj.c.
*
*   double phi0
*     (Given) The native longitude, phi_0 [deg], and ...
*   double theta0
*     (Given) ... the native latitude, theta_0 [deg], of the reference point,
*     i.e. the point (x,y) = (0,0).  If undefined (set to a magic value by
*     prjini()) the initialization routine will set this to a
*     projection-specific default.
*
*   int bounds
*     (Given) Controls bounds checking.  If bounds&1 then enable strict bounds
*     checking for the spherical-to-Cartesian (s2x) transformation for the
*     AZP, SZP, TAN, SIN, ZPN, and COP projections.  If bounds&2 then enable
*     strict bounds checking for the Cartesian-to-spherical transformation
*     (x2s) for the HPX and XPH projections.  If bounds&4 then the Cartesian-
*     to-spherical transformations (x2s) will invoke prjbchk() to perform
*     bounds checking on the computed native coordinates, with a tolerance set
*     to suit each projection.  bounds is set to 7 by prjini() by default
*     which enables all checks.  Zero it to disable all checking.
*
*     It is not necessary to reset the prjprm struct (via prjset() or
*     ???set()) when prjprm::bounds is changed.
*
* The remaining members of the prjprm struct are maintained by the setup
* routines and must not be modified elsewhere:
*
*   char name[40]
*     (Returned) Long name of the projection.
*
*     Provided for information only, not used by the projection routines.
*
*   int  category
*     (Returned) Projection category matching the value of the relevant global
*     variable:
*
*     - ZENITHAL,
*     - CYLINDRICAL,
*     - PSEUDOCYLINDRICAL,
*     - CONVENTIONAL,
*     - CONIC,
*     - POLYCONIC,
*     - QUADCUBE, and
*     - HEALPIX.
*
*     The category name may be identified via the prj_categories character
*     array, e.g.
*
=       struct prjprm prj;
=         ...
=       printf("%s\n", prj_categories[prj.category]);
*
*     Provided for information only, not used by the projection routines.
*
*   int  pvrange
*     (Returned) Range of projection parameter indices: 100 times the first
*     allowed index plus the number of parameters, e.g. TAN is 0 (no
*     parameters), SZP is 103 (1 to 3), and ZPN is 30 (0 to 29).
*
*     Provided for information only, not used by the projection routines.
*
*   int  simplezen
*     (Returned) True if the projection is a radially-symmetric zenithal
*     projection.
*
*     Provided for information only, not used by the projection routines.
*
*   int  equiareal
*     (Returned) True if the projection is equal area.
*
*     Provided for information only, not used by the projection routines.
*
*   int  conformal
*     (Returned) True if the projection is conformal.
*
*     Provided for information only, not used by the projection routines.
*
*   int  global
*     (Returned) True if the projection can represent the whole sphere in a
*     finite, non-overlapped mapping.
*
*     Provided for information only, not used by the projection routines.
*
*   int  divergent
*     (Returned) True if the projection diverges in latitude.
*
*     Provided for information only, not used by the projection routines.
*
*   double x0
*     (Returned) The offset in x, and ...
*   double y0
*     (Returned) ... the offset in y used to force (x,y) = (0,0) at
*     (phi_0,theta_0).
*
*   struct wcserr *err
*     (Returned) If enabled, when an error status is returned, this struct
*     contains detailed information about the error, see wcserr_enable().
*
*   void *padding
*     (An unused variable inserted for alignment purposes only.)
*
*   double w[10]
*     (Returned) Intermediate floating-point values derived from the
*     projection parameters, cached here to save recomputation.
*
*     Usage of the w[] array as it applies to each projection is described in
*     the prologue to each trio of projection routines in prj.c.
*
*   int n
*     (Returned) Intermediate integer value (used only for the ZPN and HPX
*     projections).
*
*   int (*prjx2s)(PRJX2S_ARGS)
*     (Returned) Pointer to the spherical projection ...
*   int (*prjs2x)(PRJ_ARGS)
*     (Returned) ... and deprojection routines.
*
*
* Global variable: const char *prj_errmsg[] - Status return messages
* ------------------------------------------------------------------
* Error messages to match the status value returned from each function.
*
*===========================================================================*/

#ifndef WCSLIB_PROJ
#define WCSLIB_PROJ

#ifdef __cplusplus
extern "C" {
#endif


/* Total number of projection parameters; 0 to PVN-1. */
#define PVN 30

extern const char *prj_errmsg[];

enum prj_errmsg_enum {
  PRJERR_SUCCESS      = 0,	/* Success. */
  PRJERR_NULL_POINTER = 1,	/* Null prjprm pointer passed. */
  PRJERR_BAD_PARAM    = 2,	/* Invalid projection parameters. */
  PRJERR_BAD_PIX      = 3,	/* One or more of the (x, y) coordinates were
				   invalid. */
  PRJERR_BAD_WORLD    = 4	/* One or more of the (phi, theta) coordinates
				   were invalid. */
};

extern const int CONIC, CONVENTIONAL, CYLINDRICAL, POLYCONIC,
                 PSEUDOCYLINDRICAL, QUADCUBE, ZENITHAL, HEALPIX;
extern const char prj_categories[9][32];

extern const int  prj_ncode;
extern const char prj_codes[28][4];

#ifdef PRJX2S_ARGS
#undef PRJX2S_ARGS
#endif

#ifdef PRJS2X_ARGS
#undef PRJS2X_ARGS
#endif

/* For use in declaring deprojection function prototypes. */
#define PRJX2S_ARGS struct prjprm *prj, int nx, int ny, int sxy, int spt, \
const double x[], const double y[], double phi[], double theta[], int stat[]

/* For use in declaring projection function prototypes. */
#define PRJS2X_ARGS struct prjprm *prj, int nx, int ny, int sxy, int spt, \
const double phi[], const double theta[], double x[], double y[], int stat[]


struct prjprm {
  /* Initialization flag (see the prologue above).                          */
  /*------------------------------------------------------------------------*/
  int    flag;			/* Set to zero to force initialization.     */

  /* Parameters to be provided (see the prologue above).                    */
  /*------------------------------------------------------------------------*/
  char   code[4];		/* Three-letter projection code.            */
  double r0;			/* Radius of the generating sphere.         */
  double pv[PVN];		/* Projection parameters.                   */
  double phi0, theta0;		/* Fiducial native coordinates.             */
  int    bounds;		/* Controls bounds checking.                */

  /* Information derived from the parameters supplied.                      */
  /*------------------------------------------------------------------------*/
  char   name[40];		/* Projection name.                         */
  int    category;		/* Projection category.                     */
  int    pvrange;		/* Range of projection parameter indices.   */
  int    simplezen;		/* Is it a simple zenithal projection?      */
  int    equiareal;		/* Is it an equal area projection?          */
  int    conformal;		/* Is it a conformal projection?            */
  int    global;		/* Can it map the whole sphere?             */
  int    divergent;		/* Does the projection diverge in latitude? */
  double x0, y0;		/* Fiducial offsets.                        */

  /* Error handling                                                         */
  /*------------------------------------------------------------------------*/
  struct wcserr *err;

  /* Private                                                                */
  /*------------------------------------------------------------------------*/
  void   *padding;		/* (Dummy inserted for alignment purposes.) */
  double w[10];			/* Intermediate values.                     */
  int    m, n;			/* Intermediate values.                     */

  int (*prjx2s)(PRJX2S_ARGS);	/* Pointers to the spherical projection and */
  int (*prjs2x)(PRJS2X_ARGS);	/* deprojection functions.                  */
};

/* Size of the prjprm struct in int units, used by the Fortran wrappers. */
#define PRJLEN (sizeof(struct prjprm)/sizeof(int))


/* Use the preprocessor to help declare function prototypes (see above). */
int prjini(struct prjprm *prj);
int prjfree(struct prjprm *prj);
int prjprt(const struct prjprm *prj);
int prjperr(const struct prjprm *prj, const char *prefix);
int prjbchk(double tol, int nx, int ny, int spt, double phi[], double theta[],
           int stat[]);

int prjset(struct prjprm *prj);
int prjx2s(PRJX2S_ARGS);
int prjs2x(PRJS2X_ARGS);

int azpset(struct prjprm *prj);
int azpx2s(PRJX2S_ARGS);
int azps2x(PRJS2X_ARGS);

int szpset(struct prjprm *prj);
int szpx2s(PRJX2S_ARGS);
int szps2x(PRJS2X_ARGS);

int tanset(struct prjprm *prj);
int tanx2s(PRJX2S_ARGS);
int tans2x(PRJS2X_ARGS);

int stgset(struct prjprm *prj);
int stgx2s(PRJX2S_ARGS);
int stgs2x(PRJS2X_ARGS);

int sinset(struct prjprm *prj);
int sinx2s(PRJX2S_ARGS);
int sins2x(PRJS2X_ARGS);

int arcset(struct prjprm *prj);
int arcx2s(PRJX2S_ARGS);
int arcs2x(PRJS2X_ARGS);

int zpnset(struct prjprm *prj);
int zpnx2s(PRJX2S_ARGS);
int zpns2x(PRJS2X_ARGS);

int zeaset(struct prjprm *prj);
int zeax2s(PRJX2S_ARGS);
int zeas2x(PRJS2X_ARGS);

int airset(struct prjprm *prj);
int airx2s(PRJX2S_ARGS);
int airs2x(PRJS2X_ARGS);

int cypset(struct prjprm *prj);
int cypx2s(PRJX2S_ARGS);
int cyps2x(PRJS2X_ARGS);

int ceaset(struct prjprm *prj);
int ceax2s(PRJX2S_ARGS);
int ceas2x(PRJS2X_ARGS);

int carset(struct prjprm *prj);
int carx2s(PRJX2S_ARGS);
int cars2x(PRJS2X_ARGS);

int merset(struct prjprm *prj);
int merx2s(PRJX2S_ARGS);
int mers2x(PRJS2X_ARGS);

int sflset(struct prjprm *prj);
int sflx2s(PRJX2S_ARGS);
int sfls2x(PRJS2X_ARGS);

int parset(struct prjprm *prj);
int parx2s(PRJX2S_ARGS);
int pars2x(PRJS2X_ARGS);

int molset(struct prjprm *prj);
int molx2s(PRJX2S_ARGS);
int mols2x(PRJS2X_ARGS);

int aitset(struct prjprm *prj);
int aitx2s(PRJX2S_ARGS);
int aits2x(PRJS2X_ARGS);

int copset(struct prjprm *prj);
int copx2s(PRJX2S_ARGS);
int cops2x(PRJS2X_ARGS);

int coeset(struct prjprm *prj);
int coex2s(PRJX2S_ARGS);
int coes2x(PRJS2X_ARGS);

int codset(struct prjprm *prj);
int codx2s(PRJX2S_ARGS);
int cods2x(PRJS2X_ARGS);

int cooset(struct prjprm *prj);
int coox2s(PRJX2S_ARGS);
int coos2x(PRJS2X_ARGS);

int bonset(struct prjprm *prj);
int bonx2s(PRJX2S_ARGS);
int bons2x(PRJS2X_ARGS);

int pcoset(struct prjprm *prj);
int pcox2s(PRJX2S_ARGS);
int pcos2x(PRJS2X_ARGS);

int tscset(struct prjprm *prj);
int tscx2s(PRJX2S_ARGS);
int tscs2x(PRJS2X_ARGS);

int cscset(struct prjprm *prj);
int cscx2s(PRJX2S_ARGS);
int cscs2x(PRJS2X_ARGS);

int qscset(struct prjprm *prj);
int qscx2s(PRJX2S_ARGS);
int qscs2x(PRJS2X_ARGS);

int hpxset(struct prjprm *prj);
int hpxx2s(PRJX2S_ARGS);
int hpxs2x(PRJS2X_ARGS);

int xphset(struct prjprm *prj);
int xphx2s(PRJX2S_ARGS);
int xphs2x(PRJS2X_ARGS);


/* Deprecated. */
#define prjini_errmsg prj_errmsg
#define prjprt_errmsg prj_errmsg
#define prjset_errmsg prj_errmsg
#define prjx2s_errmsg prj_errmsg
#define prjs2x_errmsg prj_errmsg

#ifdef __cplusplus
}
#endif

#endif /* WCSLIB_PROJ */
