/*============================================================================

  WCSLIB 7.3 - an implementation of the FITS WCS standard.
  Copyright (C) 1995-2020, Mark Calabretta

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
  $Id: wcs.h,v 7.3 2020/06/03 03:37:02 mcalabre Exp $
*=============================================================================
*
* WCSLIB 7.3 - C routines that implement the FITS World Coordinate System
* (WCS) standard.  Refer to the README file provided with WCSLIB for an
* overview of the library.
*
*
* Summary of the wcs routines
* ---------------------------
* Routines in this suite implement the FITS World Coordinate System (WCS)
* standard which defines methods to be used for computing world coordinates
* from image pixel coordinates, and vice versa.  The standard, and proposed
* extensions for handling distortions, are described in
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
=   "Representations of distortions in FITS world coordinate systems",
=   Calabretta, M.R. et al. (WCS Paper IV, draft dated 2004/04/22),
=   available from http://www.atnf.csiro.au/people/Mark.Calabretta
=
=   "Mapping on the HEALPix grid",
=   Calabretta, M.R., & Roukema, B.F. 2007, MNRAS, 381, 865 (WCS Paper V)
=
=   "Representing the 'Butterfly' Projection in FITS -- Projection Code XPH",
=   Calabretta, M.R., & Lowe, S.R. 2013, PASA, 30, e050 (WCS Paper VI)
=
=   "Representations of time coordinates in FITS -
=    Time and relative dimension in space",
=   Rots, A.H., Bunclark, P.S., Calabretta, M.R., Allen, S.L.,
=   Manchester, R.N., & Thompson, W.T. 2015, A&A, 574, A36 (WCS Paper VII)
*
* These routines are based on the wcsprm struct which contains all information
* needed for the computations.  The struct contains some members that must be
* set by the user, and others that are maintained by these routines, somewhat
* like a C++ class but with no encapsulation.
*
* wcsnpv(), wcsnps(), wcsini(), wcsinit(), wcssub(), and wcsfree() are
* provided to manage the wcsprm struct and another, wcsprt(), prints its
* contents.  Refer to the description of the wcsprm struct for an explanation
* of the anticipated usage of these routines.  wcscopy(), which does a deep
* copy of one wcsprm struct to another, is defined as a preprocessor macro
* function that invokes wcssub().
*
* wcsperr() prints the error message(s) (if any) stored in a wcsprm struct,
* and the linprm, celprm, prjprm, spcprm, and tabprm structs that it contains.
*
* A setup routine, wcsset(), computes intermediate values in the wcsprm struct
* from parameters in it that were supplied by the user.  The struct always
* needs to be set up by wcsset() but this need not be called explicitly -
* refer to the explanation of wcsprm::flag.
*
* wcsp2s() and wcss2p() implement the WCS world coordinate transformations.
* In fact, they are high level driver routines for the WCS linear,
* logarithmic, celestial, spectral and tabular transformation routines
* described in lin.h, log.h, cel.h, spc.h and tab.h.
*
* Given either the celestial longitude or latitude plus an element of the
* pixel coordinate a hybrid routine, wcsmix(), iteratively solves for the
* unknown elements.
*
* wcssptr() translates the spectral axis in a wcsprm struct.  For example, a
* 'FREQ' axis may be translated into 'ZOPT-F2W' and vice versa.
*
* wcslib_version() returns the WCSLIB version number.
*
* Quadcube projections:
* ---------------------
*   The quadcube projections (TSC, CSC, QSC) may be represented in FITS in
*   either of two ways:
*
*     a: The six faces may be laid out in one plane and numbered as follows:
*
=                                 0
=
=                        4  3  2  1  4  3  2
=
=                                 5
*
*        Faces 2, 3 and 4 may appear on one side or the other (or both).  The
*        world-to-pixel routines map faces 2, 3 and 4 to the left but the
*        pixel-to-world routines accept them on either side.
*
*     b: The "COBE" convention in which the six faces are stored in a
*        three-dimensional structure using a CUBEFACE axis indexed from
*        0 to 5 as above.
*
*   These routines support both methods; wcsset() determines which is being
*   used by the presence or absence of a CUBEFACE axis in ctype[].  wcsp2s()
*   and wcss2p() translate the CUBEFACE axis representation to the single
*   plane representation understood by the lower-level WCSLIB projection
*   routines.
*
*
* wcsnpv() - Memory allocation for PVi_ma
* ---------------------------------------
* wcsnpv() sets or gets the value of NPVMAX (default 64).  This global
* variable controls the number of pvcard structs, for holding PVi_ma
* keyvalues, that wcsini() should allocate space for.  It is also used by
* wcsinit() as the default value of npvmax.
*
* PLEASE NOTE: This function is not thread-safe.
*
* Given:
*   n         int       Value of NPVMAX; ignored if < 0.  Use a value less
*                       than zero to get the current value.
*
* Function return value:
*             int       Current value of NPVMAX.
*
*
* wcsnps() - Memory allocation for PSi_ma
* ---------------------------------------
* wcsnps() sets or gets the value of NPSMAX (default 8).  This global variable
* controls the number of pscard structs, for holding PSi_ma keyvalues, that
* wcsini() should allocate space for.  It is also used by wcsinit() as the
* default value of npsmax.
*
* PLEASE NOTE: This function is not thread-safe.
*
* Given:
*   n         int       Value of NPSMAX; ignored if < 0.  Use a value less
*                       than zero to get the current value.
*
* Function return value:
*             int       Current value of NPSMAX.
*
*
* wcsini() - Default constructor for the wcsprm struct
* ----------------------------------------------------
* wcsini() is a thin wrapper on wcsinit().  It invokes it with npvmax,
* npsmax, and ndpmax set to -1 which causes it to use the values of the
* global variables NDPMAX, NPSMAX, and NDPMAX.  It is thereby potentially
* thread-unsafe if these variables are altered dynamically via wcsnpv(),
* wcsnps(), and disndp().  Use wcsinit() for a thread-safe alternative in
* this case.
*
*
* wcsinit() - Default constructor for the wcsprm struct
* -----------------------------------------------------
* wcsinit() optionally allocates memory for arrays in a wcsprm struct and sets
* all members of the struct to default values.
*
* PLEASE NOTE: every wcsprm struct should be initialized by wcsinit(),
* possibly repeatedly.  On the first invokation, and only the first
* invokation, wcsprm::flag must be set to -1 to initialize memory management,
* regardless of whether wcsinit() will actually be used to allocate memory.
*
* Given:
*   alloc     int       If true, allocate memory unconditionally for the
*                       crpix, etc. arrays.  Please note that memory is never
*                       allocated by wcsinit() for the auxprm, tabprm, nor
*                       wtbarr structs.
*
*                       If false, it is assumed that pointers to these arrays
*                       have been set by the user except if they are null
*                       pointers in which case memory will be allocated for
*                       them regardless.  (In other words, setting alloc true
*                       saves having to initalize these pointers to zero.)
*
*   naxis     int       The number of world coordinate axes.  This is used to
*                       determine the length of the various wcsprm vectors and
*                       matrices and therefore the amount of memory to
*                       allocate for them.
*
* Given and returned:
*   wcs       struct wcsprm*
*                       Coordinate transformation parameters.
*
*                       Note that, in order to initialize memory management,
*                       wcsprm::flag should be set to -1 when wcs is
*                       initialized for the first time (memory leaks may
*                       result if it had already been initialized).
*
* Given:
*   npvmax    int       The number of PVi_ma keywords to allocate space for.
*                       If set to -1, the value of the global variable NPVMAX
*                       will be used.  This is potentially thread-unsafe if
*                       wcsnpv() is being used dynamically to alter its value.
*
*   npsmax    int       The number of PSi_ma keywords to allocate space for.
*                       If set to -1, the value of the global variable NPSMAX
*                       will be used.  This is potentially thread-unsafe if
*                       wcsnps() is being used dynamically to alter its value.
*
*   ndpmax    int       The number of DPja or DQia keywords to allocate space
*                       for.  If set to -1, the value of the global variable
*                       NDPMAX will be used.  This is potentially
*                       thread-unsafe if disndp() is being used dynamically to
*                       alter its value.
*
* Function return value:
*             int       Status return value:
*                         0: Success.
*                         1: Null wcsprm pointer passed.
*                         2: Memory allocation failed.
*
*                       For returns > 1, a detailed error message is set in
*                       wcsprm::err if enabled, see wcserr_enable().
*
*
* wcsauxi() - Default constructor for the auxprm struct
* -----------------------------------------------------
* wcsauxi() optionally allocates memory for an auxprm struct, attaches it to
* wcsprm, and sets all members of the struct to default values.
*
* Given:
*   alloc     int       If true, allocate memory unconditionally for the
*                       auxprm struct.
*
*                       If false, it is assumed that wcsprm::aux has already
*                       been set to point to an auxprm struct, in which case
*                       the user is responsible for managing that memory.
*                       However, if wcsprm::aux is a null pointer, memory will
*                       be allocated regardless.  (In other words, setting
*                       alloc true saves having to initalize the pointer to
*                       zero.)
*
* Given and returned:
*   wcs       struct wcsprm*
*                       Coordinate transformation parameters.
*
* Function return value:
*             int       Status return value:
*                         0: Success.
*                         1: Null wcsprm pointer passed.
*                         2: Memory allocation failed.
*
*
* wcssub() - Subimage extraction routine for the wcsprm struct
* ------------------------------------------------------------
* wcssub() extracts the coordinate description for a subimage from a wcsprm
* struct.  It does a deep copy, using wcsinit() to allocate memory for its
* arrays if required.  Only the "information to be provided" part of the
* struct is extracted.  Consequently, wcsset() need not have been, and won't
* be invoked on the struct from which the subimage is extracted.  A call to
* wcsset() is required to set up the subimage struct.
*
* The world coordinate system of the subimage must be separable in the sense
* that the world coordinates at any point in the subimage must depend only on
* the pixel coordinates of the axes extracted.  In practice, this means that
* the linear transformation matrix of the original image must not contain
* non-zero off-diagonal terms that associate any of the subimage axes with any
* of the non-subimage axes.  Likewise, if any distortions are associated with
* the subimage axes, they must not depend on any of the axes that are not
* being extracted.
*
* Note that while the required elements of the tabprm array are extracted, the
* wtbarr array is not.  (Thus it is not appropriate to call wcssub() after
* wcstab() but before filling the tabprm structs - refer to wcshdr.h.)
*
* wcssub() can also add axes to a wcsprm struct.  The new axes will be created
* using the defaults set by wcsinit() which produce a simple, unnamed, linear
* axis with world coordinate equal to the pixel coordinate.  These default
* values can be changed afterwards, before invoking wcsset().
*
* Given:
*   alloc     int       If true, allocate memory for the crpix, etc. arrays in
*                       the destination.  Otherwise, it is assumed that
*                       pointers to these arrays have been set by the user
*                       except if they are null pointers in which case memory
*                       will be allocated for them regardless.
*
*   wcssrc    const struct wcsprm*
*                       Struct to extract from.
*
* Given and returned:
*   nsub      int*
*   axes      int[]     Vector of length *nsub containing the image axis
*                       numbers (1-relative) to extract.  Order is
*                       significant; axes[0] is the axis number of the input
*                       image that corresponds to the first axis in the
*                       subimage, etc.
*
*                       Use an axis number of 0 to create a new axis using
*                       the defaults set by wcsinit().  They can be changed
*                       later.
*
*                       nsub (the pointer) may be set to zero, and so also may
*                       *nsub, which is interpreted to mean all axes in the
*                       input image; the number of axes will be returned if
*                       nsub != 0x0.  axes itself (the pointer) may be set to
*                       zero to indicate the first *nsub axes in their
*                       original order.
*
*                       Set both nsub (or *nsub) and axes to zero to do a deep
*                       copy of one wcsprm struct to another.
*
*                       Subimage extraction by coordinate axis type may be
*                       done by setting the elements of axes[] to the
*                       following special preprocessor macro values:
*
*                         WCSSUB_LONGITUDE: Celestial longitude.
*                         WCSSUB_LATITUDE:  Celestial latitude.
*                         WCSSUB_CUBEFACE:  Quadcube CUBEFACE axis.
*                         WCSSUB_SPECTRAL:  Spectral axis.
*                         WCSSUB_STOKES:    Stokes axis.
*
*                       Refer to the notes (below) for further usage examples.
*
*                       On return, *nsub will be set to the number of axes in
*                       the subimage; this may be zero if there were no axes
*                       of the required type(s) (in which case no memory will
*                       be allocated).  axes[] will contain the axis numbers
*                       that were extracted, or 0 for newly created axes.  The
*                       vector length must be sufficient to contain all axis
*                       numbers.  No checks are performed to verify that the
*                       coordinate axes are consistent, this is done by
*                       wcsset().
*
*   wcsdst    struct wcsprm*
*                       Struct describing the subimage.  wcsprm::flag should
*                       be set to -1 if wcsdst was not previously initialized
*                       (memory leaks may result if it was previously
*                       initialized).
*
* Function return value:
*             int       Status return value:
*                         0: Success.
*                         1: Null wcsprm pointer passed.
*                         2: Memory allocation failed.
*                        12: Invalid subimage specification.
*                        13: Non-separable subimage coordinate system.
*
*                       For returns > 1, a detailed error message is set in
*                       wcsprm::err if enabled, see wcserr_enable().
*
* Notes:
*   Combinations of subimage axes of particular types may be extracted in the
*   same order as they occur in the input image by combining preprocessor
*   codes, for example
*
=     *nsub = 1;
=     axes[0] = WCSSUB_LONGITUDE | WCSSUB_LATITUDE | WCSSUB_SPECTRAL;
*
*   would extract the longitude, latitude, and spectral axes in the same order
*   as the input image.  If one of each were present, *nsub = 3 would be
*   returned.
*
*   For convenience, WCSSUB_CELESTIAL is defined as the combination
*   WCSSUB_LONGITUDE | WCSSUB_LATITUDE | WCSSUB_CUBEFACE.
*
*   The codes may also be negated to extract all but the types specified, for
*   example
*
=     *nsub = 4;
=     axes[0] = WCSSUB_LONGITUDE;
=     axes[1] = WCSSUB_LATITUDE;
=     axes[2] = WCSSUB_CUBEFACE;
=     axes[3] = -(WCSSUB_SPECTRAL | WCSSUB_STOKES);
*
*   The last of these specifies all axis types other than spectral or Stokes.
*   Extraction is done in the order specified by axes[] a longitude axis (if
*   present) would be extracted first (via axes[0]) and not subsequently (via
*   axes[3]).  Likewise for the latitude and cubeface axes in this example.
*
*   From the foregoing, it is apparent that the value of *nsub returned may be
*   less than or greater than that given.  However, it will never exceed the
*   number of axes in the input image (plus the number of newly-created axes
*   if any were specified on input).
*
*
* wcscompare() - Compare two wcsprm structs for equality
* ------------------------------------------------------
* wcscompare() compares two wcsprm structs for equality.
*
* Given:
*   cmp       int       A bit field controlling the strictness of the
*                       comparison.  When 0, all fields must be identical.
*
*                       The following constants may be or'ed together to
*                       relax the comparison:
*                         WCSCOMPARE_ANCILLARY: Ignore ancillary keywords
*                           that don't change the WCS transformation, such
*                           as DATE-OBS or EQUINOX.
*                         WCSCOMPARE_TILING: Ignore integral differences in
*                           CRPIXja.  This is the 'tiling' condition, where
*                           two WCSes cover different regions of the same
*                           map projection and align on the same map grid.
*                         WCSCOMPARE_CRPIX: Ignore any differences at all in
*                           CRPIXja.  The two WCSes cover different regions
*                           of the same map projection but may not align on
*                           the same map grid.  Overrides WCSCOMPARE_TILING.
*
*   tol       double    Tolerance for comparison of floating-point values.
*                       For example, for tol == 1e-6, all floating-point
*                       values in the structs must be equal to the first 6
*                       decimal places.  A value of 0 implies exact equality.
*
*   wcs1      const struct wcsprm*
*                       The first wcsprm struct to compare.
*
*   wcs2      const struct wcsprm*
*                       The second wcsprm struct to compare.
*
* Returned:
*   equal     int*      Non-zero when the given structs are equal.
*
* Function return value:
*             int       Status return value:
*                         0: Success.
*                         1: Null pointer passed.
*
*
* wcscopy() macro - Copy routine for the wcsprm struct
* ----------------------------------------------------
* wcscopy() does a deep copy of one wcsprm struct to another.  As of
* WCSLIB 3.6, it is implemented as a preprocessor macro that invokes
* wcssub() with the nsub and axes pointers both set to zero.
*
*
* wcsfree() - Destructor for the wcsprm struct
* --------------------------------------------
* wcsfree() frees memory allocated for the wcsprm arrays by wcsinit() and/or
* wcsset().  wcsinit() records the memory it allocates and wcsfree() will only
* attempt to free this.
*
* PLEASE NOTE: wcsfree() must not be invoked on a wcsprm struct that was not
* initialized by wcsinit().
*
* Returned:
*   wcs       struct wcsprm*
*                       Coordinate transformation parameters.
*
* Function return value:
*             int       Status return value:
*                         0: Success.
*                         1: Null wcsprm pointer passed.
*
*
* wcsprt() - Print routine for the wcsprm struct
* ----------------------------------------------
* wcsprt() prints the contents of a wcsprm struct using wcsprintf().  Mainly
* intended for diagnostic purposes.
*
* Given:
*   wcs       const struct wcsprm*
*                       Coordinate transformation parameters.
*
* Function return value:
*             int       Status return value:
*                         0: Success.
*                         1: Null wcsprm pointer passed.
*
*
* wcsperr() - Print error messages from a wcsprm struct
* -----------------------------------------------------
* wcsperr() prints the error message(s), if any, stored in a wcsprm struct,
* and the linprm, celprm, prjprm, spcprm, and tabprm structs that it contains.
* If there are no errors then nothing is printed.  It uses wcserr_prt(), q.v.
*
* Given:
*   wcs       const struct wcsprm*
*                       Coordinate transformation parameters.
*
*   prefix    const char *
*                       If non-NULL, each output line will be prefixed with
*                       this string.
*
* Function return value:
*             int       Status return value:
*                         0: Success.
*                         1: Null wcsprm pointer passed.
*
*
* wcsbchk() - Enable/disable bounds checking
* ------------------------------------------
* wcsbchk() is used to control bounds checking in the projection routines.
* Note that wcsset() always enables bounds checking.  wcsbchk() will invoke
* wcsset() on the wcsprm struct beforehand if necessary.
*
* Given and returned:
*   wcs       struct wcsprm*
*                       Coordinate transformation parameters.
*
* Given:
*   bounds    int       If bounds&1 then enable strict bounds checking for the
*                       spherical-to-Cartesian (s2x) transformation for the
*                       AZP, SZP, TAN, SIN, ZPN, and COP projections.
*
*                       If bounds&2 then enable strict bounds checking for the
*                       Cartesian-to-spherical (x2s) transformation for the
*                       HPX and XPH projections.
*
*                       If bounds&4 then enable bounds checking on the native
*                       coordinates returned by the Cartesian-to-spherical
*                       (x2s) transformations using prjchk().
*
*                       Zero it to disable all checking.
*
* Function return value:
*             int       Status return value:
*                         0: Success.
*                         1: Null wcsprm pointer passed.
*
*
* wcsset() - Setup routine for the wcsprm struct
* ----------------------------------------------
* wcsset() sets up a wcsprm struct according to information supplied within
* it (refer to the description of the wcsprm struct).
*
* wcsset() recognizes the NCP projection and converts it to the equivalent SIN
* projection and likewise translates GLS into SFL.  It also translates the
* AIPS spectral types ('FREQ-LSR', 'FELO-HEL', etc.), possibly changing the
* input header keywords wcsprm::ctype and/or wcsprm::specsys if necessary.
*
* Note that this routine need not be called directly; it will be invoked by
* wcsp2s() and wcss2p() if the wcsprm::flag is anything other than a
* predefined magic value.
*
* Given and returned:
*   wcs       struct wcsprm*
*                       Coordinate transformation parameters.
*
* Function return value:
*             int       Status return value:
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
* Notes:
*   wcsset() always enables strict bounds checking in the projection routines
*   (via a call to prjini()).  Use wcsbchk() to modify bounds-checking after
*   wcsset() is invoked.
*
*
* wcsp2s() - Pixel-to-world transformation
* ----------------------------------------
* wcsp2s() transforms pixel coordinates to world coordinates.
*
* Given and returned:
*   wcs       struct wcsprm*
*                       Coordinate transformation parameters.
*
* Given:
*   ncoord,
*   nelem     int       The number of coordinates, each of vector length
*                       nelem but containing wcs.naxis coordinate elements.
*                       Thus nelem must equal or exceed the value of the
*                       NAXIS keyword unless ncoord == 1, in which case nelem
*                       is not used.
*
*   pixcrd    const double[ncoord][nelem]
*                       Array of pixel coordinates.
*
* Returned:
*   imgcrd    double[ncoord][nelem]
*                       Array of intermediate world coordinates.  For
*                       celestial axes, imgcrd[][wcs.lng] and
*                       imgcrd[][wcs.lat] are the projected x-, and
*                       y-coordinates in pseudo "degrees".  For spectral
*                       axes, imgcrd[][wcs.spec] is the intermediate spectral
*                       coordinate, in SI units.
*
*   phi,theta double[ncoord]
*                       Longitude and latitude in the native coordinate system
*                       of the projection [deg].
*
*   world     double[ncoord][nelem]
*                       Array of world coordinates.  For celestial axes,
*                       world[][wcs.lng] and world[][wcs.lat] are the
*                       celestial longitude and latitude [deg].  For
*                       spectral axes, imgcrd[][wcs.spec] is the intermediate
*                       spectral coordinate, in SI units.
*
*   stat      int[ncoord]
*                       Status return value for each coordinate:
*                         0: Success.
*                        1+: A bit mask indicating invalid pixel coordinate
*                            element(s).
*
* Function return value:
*             int       Status return value:
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
*                         8: One or more of the pixel coordinates were
*                            invalid, as indicated by the stat vector.
*
*                       For returns > 1, a detailed error message is set in
*                       wcsprm::err if enabled, see wcserr_enable().
*
*
* wcss2p() - World-to-pixel transformation
* ----------------------------------------
* wcss2p() transforms world coordinates to pixel coordinates.
*
* Given and returned:
*   wcs       struct wcsprm*
*                       Coordinate transformation parameters.
*
* Given:
*   ncoord,
*   nelem     int       The number of coordinates, each of vector length nelem
*                       but containing wcs.naxis coordinate elements.  Thus
*                       nelem must equal or exceed the value of the NAXIS
*                       keyword unless ncoord == 1, in which case nelem is not
*                       used.
*
*   world     const double[ncoord][nelem]
*                       Array of world coordinates.  For celestial axes,
*                       world[][wcs.lng] and world[][wcs.lat] are the
*                       celestial longitude and latitude [deg]. For spectral
*                       axes, world[][wcs.spec] is the spectral coordinate, in
*                       SI units.
*
* Returned:
*   phi,theta double[ncoord]
*                       Longitude and latitude in the native coordinate
*                       system of the projection [deg].
*
*   imgcrd    double[ncoord][nelem]
*                       Array of intermediate world coordinates.  For
*                       celestial axes, imgcrd[][wcs.lng] and
*                       imgcrd[][wcs.lat] are the projected x-, and
*                       y-coordinates in pseudo "degrees".  For quadcube
*                       projections with a CUBEFACE axis the face number is
*                       also returned in imgcrd[][wcs.cubeface].  For
*                       spectral axes, imgcrd[][wcs.spec] is the intermediate
*                       spectral coordinate, in SI units.
*
*   pixcrd    double[ncoord][nelem]
*                       Array of pixel coordinates.
*
*   stat      int[ncoord]
*                       Status return value for each coordinate:
*                         0: Success.
*                        1+: A bit mask indicating invalid world coordinate
*                            element(s).
*
* Function return value:
*             int       Status return value:
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
*                         9: One or more of the world coordinates were
*                            invalid, as indicated by the stat vector.
*
*                       For returns > 1, a detailed error message is set in
*                       wcsprm::err if enabled, see wcserr_enable().
*
*
* wcsmix() - Hybrid coordinate transformation
* -------------------------------------------
* wcsmix(), given either the celestial longitude or latitude plus an element
* of the pixel coordinate, solves for the remaining elements by iterating on
* the unknown celestial coordinate element using wcss2p().  Refer also to the
* notes below.
*
* Given and returned:
*   wcs       struct wcsprm*
*                       Indices for the celestial coordinates obtained
*                       by parsing the wcsprm::ctype[].
*
* Given:
*   mixpix    int       Which element of the pixel coordinate is given.
*
*   mixcel    int       Which element of the celestial coordinate is given:
*                         1: Celestial longitude is given in
*                            world[wcs.lng], latitude returned in
*                            world[wcs.lat].
*                         2: Celestial latitude is given in
*                            world[wcs.lat], longitude returned in
*                            world[wcs.lng].
*
*   vspan     const double[2]
*                       Solution interval for the celestial coordinate [deg].
*                       The ordering of the two limits is irrelevant.
*                       Longitude ranges may be specified with any convenient
*                       normalization, for example [-120,+120] is the same as
*                       [240,480], except that the solution will be returned
*                       with the same normalization, i.e. lie within the
*                       interval specified.
*
*   vstep     const double
*                       Step size for solution search [deg].  If zero, a
*                       sensible, although perhaps non-optimal default will be
*                       used.
*
*   viter     int       If a solution is not found then the step size will be
*                       halved and the search recommenced.  viter controls how
*                       many times the step size is halved.  The allowed range
*                       is 5 - 10.
*
* Given and returned:
*   world     double[naxis]
*                       World coordinate elements.  world[wcs.lng] and
*                       world[wcs.lat] are the celestial longitude and
*                       latitude [deg].  Which is given and which returned
*                       depends on the value of mixcel.  All other elements
*                       are given.
*
* Returned:
*   phi,theta double[naxis]
*                       Longitude and latitude in the native coordinate
*                       system of the projection [deg].
*
*   imgcrd    double[naxis]
*                       Image coordinate elements.  imgcrd[wcs.lng] and
*                       imgcrd[wcs.lat] are the projected x-, and
*                       y-coordinates in pseudo "degrees".
*
* Given and returned:
*   pixcrd    double[naxis]
*                       Pixel coordinate.  The element indicated by mixpix is
*                       given and the remaining elements are returned.
*
* Function return value:
*             int       Status return value:
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
*                        10: Invalid world coordinate.
*                        11: No solution found in the specified interval.
*
*                       For returns > 1, a detailed error message is set in
*                       wcsprm::err if enabled, see wcserr_enable().
*
* Notes:
*   Initially the specified solution interval is checked to see if it's a
*   "crossing" interval.  If it isn't, a search is made for a crossing
*   solution by iterating on the unknown celestial coordinate starting at the
*   upper limit of the solution interval and decrementing by the specified
*   step size.  A crossing is indicated if the trial value of the pixel
*   coordinate steps through the value specified.  If a crossing interval is
*   found then the solution is determined by a modified form of "regula falsi"
*   division of the crossing interval.  If no crossing interval was found
*   within the specified solution interval then a search is made for a
*   "non-crossing" solution as may arise from a point of tangency.  The
*   process is complicated by having to make allowance for the discontinuities
*   that occur in all map projections.
*
*   Once one solution has been determined others may be found by subsequent
*   invokations of wcsmix() with suitably restricted solution intervals.
*
*   Note the circumstance that arises when the solution point lies at a native
*   pole of a projection in which the pole is represented as a finite curve,
*   for example the zenithals and conics.  In such cases two or more valid
*   solutions may exist but wcsmix() only ever returns one.
*
*   Because of its generality wcsmix() is very compute-intensive.  For
*   compute-limited applications more efficient special-case solvers could be
*   written for simple projections, for example non-oblique cylindrical
*   projections.
*
*
* wcssptr() - Spectral axis translation
* -------------------------------------
* wcssptr() translates the spectral axis in a wcsprm struct.  For example, a
* 'FREQ' axis may be translated into 'ZOPT-F2W' and vice versa.
*
* Given and returned:
*   wcs       struct wcsprm*
*                       Coordinate transformation parameters.
*
*   i         int*      Index of the spectral axis (0-relative).  If given < 0
*                       it will be set to the first spectral axis identified
*                       from the ctype[] keyvalues in the wcsprm struct.
*
*   ctype     char[9]   Desired spectral CTYPEia.  Wildcarding may be used as
*                       for the ctypeS2 argument to spctrn() as described in
*                       the prologue of spc.h, i.e. if the final three
*                       characters are specified as "???", or if just the
*                       eighth character is specified as '?', the correct
*                       algorithm code will be substituted and returned.
*
* Function return value:
*             int       Status return value:
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
*                        12: Invalid subimage specification (no spectral
*                            axis).
*
*                       For returns > 1, a detailed error message is set in
*                       wcsprm::err if enabled, see wcserr_enable().
*
*
* wcslib_version() - WCSLIB version number
* ----------------------------------------
* wcslib_version() returns the WCSLIB version number.
*
* The major version number changes when the ABI changes or when the license
* conditions change.  ABI changes typically result from a change to the
* contents of one of the structs.  The major version number is used to
* distinguish between incompatible versions of the sharable library.
*
* The minor version number changes with new functionality or bug fixes that do
* not involve a change in the ABI.
*
* The auxiliary version number (which is often absent) signals changes to the
* documentation, test suite, build procedures, or any other change that does
* not affect the compiled library.
*
* Returned:
*   vers[3]   int[3]    The broken-down version number:
*                         0: Major version number.
*                         1: Minor version number.
*                         2: Auxiliary version number (zero if absent).
*                       May be given as a null pointer if not required.
*
* Function return value:
*             char*     A null-terminated, statically allocated string
*                       containing the version number in the usual form, i.e.
*                       "<major>.<minor>.<auxiliary>".
*
*
* wcsprm struct - Coordinate transformation parameters
* ----------------------------------------------------
* The wcsprm struct contains information required to transform world
* coordinates.  It consists of certain members that must be set by the user
* ("given") and others that are set by the WCSLIB routines ("returned").
* While the addresses of the arrays themselves may be set by wcsinit() if it
* (optionally) allocates memory, their contents must be set by the user.
*
* Some parameters that are given are not actually required for transforming
* coordinates.  These are described as "auxiliary"; the struct simply provides
* a place to store them, though they may be used by wcshdo() in constructing a
* FITS header from a wcsprm struct.  Some of the returned values are supplied
* for informational purposes and others are for internal use only as
* indicated.
*
* In practice, it is expected that a WCS parser would scan the FITS header to
* determine the number of coordinate axes.  It would then use wcsinit() to
* allocate memory for arrays in the wcsprm struct and set default values.
* Then as it reread the header and identified each WCS keyrecord it would load
* the value into the relevant wcsprm array element.  This is essentially what
* wcspih() does - refer to the prologue of wcshdr.h.  As the final step,
* wcsset() is invoked, either directly or indirectly, to set the derived
* members of the wcsprm struct.  wcsset() strips off trailing blanks in all
* string members and null-fills the character array.
*
*   int flag
*     (Given and returned) This flag must be set to zero whenever any of the
*     following wcsprm struct members are set or changed:
*
*       - wcsprm::naxis (q.v., not normally set by the user),
*       - wcsprm::crpix,
*       - wcsprm::pc,
*       - wcsprm::cdelt,
*       - wcsprm::crval,
*       - wcsprm::cunit,
*       - wcsprm::ctype,
*       - wcsprm::lonpole,
*       - wcsprm::latpole,
*       - wcsprm::restfrq,
*       - wcsprm::restwav,
*       - wcsprm::npv,
*       - wcsprm::pv,
*       - wcsprm::nps,
*       - wcsprm::ps,
*       - wcsprm::cd,
*       - wcsprm::crota,
*       - wcsprm::altlin,
*       - wcsprm::ntab,
*       - wcsprm::nwtb,
*       - wcsprm::tab,
*       - wcsprm::wtb.
*
*     This signals the initialization routine, wcsset(), to recompute the
*     returned members of the celprm struct.  celset() will reset flag to
*     indicate that this has been done.
*
*     PLEASE NOTE: flag should be set to -1 when wcsinit() is called for the
*     first time for a particular wcsprm struct in order to initialize memory
*     management.  It must ONLY be used on the first initialization otherwise
*     memory leaks may result.
*
*   int naxis
*     (Given or returned) Number of pixel and world coordinate elements.
*
*     If wcsinit() is used to initialize the linprm struct (as would normally
*     be the case) then it will set naxis from the value passed to it as a
*     function argument.  The user should not subsequently modify it.
*
*   double *crpix
*     (Given) Address of the first element of an array of double containing
*     the coordinate reference pixel, CRPIXja.
*
*   double *pc
*     (Given) Address of the first element of the PCi_ja (pixel coordinate)
*     transformation matrix.  The expected order is
*
=       struct wcsprm wcs;
=       wcs.pc = {PC1_1, PC1_2, PC2_1, PC2_2};
*
*     This may be constructed conveniently from a 2-D array via
*
=       double m[2][2] = {{PC1_1, PC1_2},
=                         {PC2_1, PC2_2}};
*
*     which is equivalent to
*
=       double m[2][2];
=       m[0][0] = PC1_1;
=       m[0][1] = PC1_2;
=       m[1][0] = PC2_1;
=       m[1][1] = PC2_2;
*
*     The storage order for this 2-D array is the same as for the 1-D array,
*     whence
*
=       wcs.pc = *m;
*
*     would be legitimate.
*
*   double *cdelt
*     (Given) Address of the first element of an array of double containing
*     the coordinate increments, CDELTia.
*
*   double *crval
*     (Given) Address of the first element of an array of double containing
*     the coordinate reference values, CRVALia.
*
*   char (*cunit)[72]
*     (Given) Address of the first element of an array of char[72] containing
*     the CUNITia keyvalues which define the units of measurement of the
*     CRVALia, CDELTia, and CDi_ja keywords.
*
*     As CUNITia is an optional header keyword, cunit[][72] may be left blank
*     but otherwise is expected to contain a standard units specification as
*     defined by WCS Paper I.  Utility function wcsutrn(), described in
*     wcsunits.h, is available to translate commonly used non-standard units
*     specifications but this must be done as a separate step before invoking
*     wcsset().
*
*     For celestial axes, if cunit[][72] is not blank, wcsset() uses
*     wcsunits() to parse it and scale cdelt[], crval[], and cd[][*] to
*     degrees.  It then resets cunit[][72] to "deg".
*
*     For spectral axes, if cunit[][72] is not blank, wcsset() uses wcsunits()
*     to parse it and scale cdelt[], crval[], and cd[][*] to SI units.  It
*     then resets cunit[][72] accordingly.
*
*     wcsset() ignores cunit[][72] for other coordinate types; cunit[][72] may
*     be used to label coordinate values.
*
*     These variables accomodate the longest allowed string-valued FITS
*     keyword, being limited to 68 characters, plus the null-terminating
*     character.
*
*   char (*ctype)[72]
*     (Given) Address of the first element of an array of char[72] containing
*     the coordinate axis types, CTYPEia.
*
*     The ctype[][72] keyword values must be in upper case and there must be
*     zero or one pair of matched celestial axis types, and zero or one
*     spectral axis.  The ctype[][72] strings should be padded with blanks on
*     the right and null-terminated so that they are at least eight characters
*     in length.
*
*     These variables accomodate the longest allowed string-valued FITS
*     keyword, being limited to 68 characters, plus the null-terminating
*     character.
*
*   double lonpole
*     (Given and returned) The native longitude of the celestial pole, phi_p,
*     given by LONPOLEa [deg] or by PVi_2a [deg] attached to the longitude
*     axis which takes precedence if defined, and ...
*   double latpole
*     (Given and returned) ... the native latitude of the celestial pole,
*     theta_p, given by LATPOLEa [deg] or by PVi_3a [deg] attached to the
*     longitude axis which takes precedence if defined.
*
*     lonpole and latpole may be left to default to values set by wcsinit()
*     (see celprm::ref), but in any case they will be reset by wcsset() to
*     the values actually used.  Note therefore that if the wcsprm struct is
*     reused without resetting them, whether directly or via wcsinit(), they
*     will no longer have their default values.
*
*   double restfrq
*     (Given) The rest frequency [Hz], and/or ...
*   double restwav
*     (Given) ... the rest wavelength in vacuo [m], only one of which need be
*     given, the other should be set to zero.
*
*   int npv
*     (Given) The number of entries in the wcsprm::pv[] array.
*
*   int npvmax
*     (Given or returned) The length of the wcsprm::pv[] array.
*
*     npvmax will be set by wcsinit() if it allocates memory for wcsprm::pv[],
*     otherwise it must be set by the user.  See also wcsnpv().
*
*   struct pvcard *pv
*     (Given) Address of the first element of an array of length npvmax of
*     pvcard structs.
*
*     As a FITS header parser encounters each PVi_ma keyword it should load it
*     into a pvcard struct in the array and increment npv.  wcsset()
*     interprets these as required.
*
*     Note that, if they were not given, wcsset() resets the entries for
*     PVi_1a, PVi_2a, PVi_3a, and PVi_4a for longitude axis i to match
*     phi_0 and theta_0 (the native longitude and latitude of the reference
*     point), LONPOLEa and LATPOLEa respectively.
*
*   int nps
*     (Given) The number of entries in the wcsprm::ps[] array.
*
*   int npsmax
*     (Given or returned) The length of the wcsprm::ps[] array.
*
*     npsmax will be set by wcsinit() if it allocates memory for wcsprm::ps[],
*     otherwise it must be set by the user.  See also wcsnps().
*
*   struct pscard *ps
*     (Given) Address of the first element of an array of length npsmax of
*     pscard structs.
*
*     As a FITS header parser encounters each PSi_ma keyword it should load it
*     into a pscard struct in the array and increment nps.  wcsset()
*     interprets these as required (currently no PSi_ma keyvalues are
*     recognized).
*
*   double *cd
*     (Given) For historical compatibility, the wcsprm struct supports two
*     alternate specifications of the linear transformation matrix, those
*     associated with the CDi_ja keywords, and ...
*   double *crota
*     (Given) ... those associated with the CROTAi keywords.  Although these
*     may not formally co-exist with PCi_ja, the approach taken here is simply
*     to ignore them if given in conjunction with PCi_ja.
*
*   int altlin
*     (Given) altlin is a bit flag that denotes which of the PCi_ja, CDi_ja
*     and CROTAi keywords are present in the header:
*
*     - Bit 0: PCi_ja is present.
*
*     - Bit 1: CDi_ja is present.
*
*       Matrix elements in the IRAF convention are
*       equivalent to the product CDi_ja = CDELTia * PCi_ja, but the
*       defaults differ from that of the PCi_ja matrix.  If one or more
*       CDi_ja keywords are present then all unspecified CDi_ja default to
*       zero.  If no CDi_ja (or CROTAi) keywords are present, then the
*       header is assumed to be in PCi_ja form whether or not any PCi_ja
*       keywords are present since this results in an interpretation of
*       CDELTia consistent with the original FITS specification.
*
*       While CDi_ja may not formally co-exist with PCi_ja, it may co-exist
*       with CDELTia and CROTAi which are to be ignored.
*
*     - Bit 2: CROTAi is present.
*
*       In the AIPS convention, CROTAi may only be
*       associated with the latitude axis of a celestial axis pair.  It
*       specifies a rotation in the image plane that is applied AFTER the
*       CDELTia; any other CROTAi keywords are ignored.
*
*       CROTAi may not formally co-exist with PCi_ja.
*
*       CROTAi and CDELTia may formally co-exist with CDi_ja but if so are to
*       be ignored.
*
*     CDi_ja and CROTAi keywords, if found, are to be stored in the
*     wcsprm::cd and wcsprm::crota arrays which are dimensioned similarly to
*     wcsprm::pc and wcsprm::cdelt.  FITS
*     header parsers should use the following procedure:
*
*     - Whenever a PCi_ja  keyword is encountered: altlin |= 1;
*
*     - Whenever a CDi_ja  keyword is encountered: altlin |= 2;
*
*     - Whenever a CROTAi keyword is encountered: altlin |= 4;
*
*     If none of these bits are set the PCi_ja representation results, i.e.
*     wcsprm::pc and wcsprm::cdelt will be used as given.
*
*     These alternate specifications of the linear transformation matrix are
*     translated immediately to PCi_ja by wcsset() and are invisible to the
*     lower-level WCSLIB routines.  In particular, wcsset() resets
*     wcsprm::cdelt to unity if CDi_ja is present (and no PCi_ja).
*
*     If CROTAi are present but none is associated with the latitude axis
*     (and no PCi_ja or CDi_ja), then wcsset() reverts to a unity PCi_ja
*     matrix.
*
*   int velref
*     (Given) AIPS velocity code VELREF, refer to spcaips().
*
*     It is not necessary to reset the wcsprm struct (via wcsset()) when
*     wcsprm::velref is changed.
*
*   char alt[4]
*     (Given, auxiliary) Character code for alternate coordinate descriptions
*     (i.e. the 'a' in keyword names such as CTYPEia).  This is blank for the
*     primary coordinate description, or one of the 26 upper-case letters,
*     A-Z.
*
*     An array of four characters is provided for alignment purposes, only the
*     first is used.
*
*     It is not necessary to reset the wcsprm struct (via wcsset()) when
*     wcsprm::alt is changed.
*
*   int colnum
*     (Given, auxiliary) Where the coordinate representation is associated
*     with an image-array column in a FITS binary table, this variable may be
*     used to record the relevant column number.
*
*     It should be set to zero for an image header or pixel list.
*
*     It is not necessary to reset the wcsprm struct (via wcsset()) when
*     wcsprm::colnum is changed.
*
*   int *colax
*     (Given, auxiliary) Address of the first element of an array of int
*     recording the column numbers for each axis in a pixel list.
*
*     The array elements should be set to zero for an image header or image
*     array in a binary table.
*
*     It is not necessary to reset the wcsprm struct (via wcsset()) when
*     wcsprm::colax is changed.
*
*   char (*cname)[72]
*     (Given, auxiliary) The address of the first element of an array of
*     char[72] containing the coordinate axis names, CNAMEia.
*
*     These variables accomodate the longest allowed string-valued FITS
*     keyword, being limited to 68 characters, plus the null-terminating
*     character.
*
*     It is not necessary to reset the wcsprm struct (via wcsset()) when
*     wcsprm::cname is changed.
*
*   double *crder
*     (Given, auxiliary) Address of the first element of an array of double
*     recording the random error in the coordinate value, CRDERia.
*
*     It is not necessary to reset the wcsprm struct (via wcsset()) when
*     wcsprm::crder is changed.
*
*   double *csyer
*     (Given, auxiliary) Address of the first element of an array of double
*     recording the systematic error in the coordinate value, CSYERia.
*
*     It is not necessary to reset the wcsprm struct (via wcsset()) when
*     wcsprm::csyer is changed.
*
*   double *czphs
*     (Given, auxiliary) Address of the first element of an array of double
*     recording the time at the zero point of a phase axis, CZPHSia.
*
*     It is not necessary to reset the wcsprm struct (via wcsset()) when
*     wcsprm::czphs is changed.
*
*   double *cperi
*     (Given, auxiliary) Address of the first element of an array of double
*     recording the period of a phase axis, CPERIia.
*
*     It is not necessary to reset the wcsprm struct (via wcsset()) when
*     wcsprm::cperi is changed.
*
*   char wcsname[72]
*     (Given, auxiliary) The name given to the coordinate representation,
*     WCSNAMEa.  This variable accomodates the longest allowed string-valued
*     FITS keyword, being limited to 68 characters, plus the null-terminating
*     character.
*
*     It is not necessary to reset the wcsprm struct (via wcsset()) when
*     wcsprm::wcsname is changed.
*
*   char timesys[72]
*     (Given, auxiliary) TIMESYS keyvalue, being the time scale (UTC, TAI,
*     etc.) in which all other time-related auxiliary header values are
*     recorded.  Also defines the time scale for an image axis with CTYPEia
*     set to 'TIME'.
*
*     It is not necessary to reset the wcsprm struct (via wcsset()) when
*     wcsprm::timesys is changed.
*
*   char trefpos[72]
*     (Given, auxiliary) TREFPOS keyvalue, being the location in space where
*     the recorded time is valid.
*
*     It is not necessary to reset the wcsprm struct (via wcsset()) when
*     wcsprm::trefpos is changed.
*
*   char trefdir[72]
*     (Given, auxiliary) TREFDIR keyvalue, being the reference direction used
*     in calculating a pathlength delay.
*
*     It is not necessary to reset the wcsprm struct (via wcsset()) when
*     wcsprm::trefdir is changed.
*
*   char plephem[72]
*     (Given, auxiliary) PLEPHEM keyvalue, being the Solar System ephemeris
*     used for calculating a pathlength delay.
*
*     It is not necessary to reset the wcsprm struct (via wcsset()) when
*     wcsprm::plephem is changed.
*
*   char timeunit[72]
*     (Given, auxiliary) TIMEUNIT keyvalue, being the time units in which
*     the following header values are expressed: TSTART, TSTOP, TIMEOFFS,
*     TIMSYER, TIMRDER, TIMEDEL.  It also provides the default value for
*     CUNITia for time axes.
*
*     It is not necessary to reset the wcsprm struct (via wcsset()) when
*     wcsprm::timeunit is changed.
*
*   char dateref[72]
*     (Given, auxiliary) DATEREF keyvalue, being the date of a reference epoch
*     relative to which other time measurements refer.
*
*     It is not necessary to reset the wcsprm struct (via wcsset()) when
*     wcsprm::dateref is changed.
*
*   double mjdref[2]
*     (Given, auxiliary) MJDREF keyvalue, equivalent to DATEREF expressed as
*     a Modified Julian Date (MJD = JD - 2400000.5).  The value is given as
*     the sum of the two-element vector, allowing increased precision.
*
*     It is not necessary to reset the wcsprm struct (via wcsset()) when
*     wcsprm::mjdref is changed.
*
*   double timeoffs
*     (Given, auxiliary) TIMEOFFS keyvalue, being a time offset, which may be
*     used, for example, to provide a uniform clock correction for times
*     referenced to DATEREF.
*
*     It is not necessary to reset the wcsprm struct (via wcsset()) when
*     wcsprm::timeoffs is changed.
*
*   char dateobs[72]
*     (Given, auxiliary) DATE-OBS keyvalue, being the date at the start of the
*     observation unless otherwise explained in the DATE-OBS keycomment, in
*     ISO format, yyyy-mm-ddThh:mm:ss.
*
*     It is not necessary to reset the wcsprm struct (via wcsset()) when
*     wcsprm::dateobs is changed.
*
*   char datebeg[72]
*     (Given, auxiliary) DATE-BEG keyvalue, being the date at the start of the
*     observation in ISO format, yyyy-mm-ddThh:mm:ss.
*
*     It is not necessary to reset the wcsprm struct (via wcsset()) when
*     wcsprm::datebeg is changed.
*
*   char dateavg[72]
*     (Given, auxiliary) DATE-AVG keyvalue, being the date at a representative
*     mid-point of the observation in ISO format, yyyy-mm-ddThh:mm:ss.
*
*     It is not necessary to reset the wcsprm struct (via wcsset()) when
*     wcsprm::dateavg is changed.
*
*   char dateend[72]
*     (Given, auxiliary) DATE-END keyvalue, baing the date at the end of the
*     observation in ISO format, yyyy-mm-ddThh:mm:ss.
*
*     It is not necessary to reset the wcsprm struct (via wcsset()) when
*     wcsprm::dateend is changed.
*
*   double mjdobs
*     (Given, auxiliary) MJD-OBS keyvalue, equivalent to DATE-OBS expressed
*     as a Modified Julian Date (MJD = JD - 2400000.5).
*
*     It is not necessary to reset the wcsprm struct (via wcsset()) when
*     wcsprm::mjdobs is changed.
*
*   double mjdbeg
*     (Given, auxiliary) MJD-BEG keyvalue, equivalent to DATE-BEG expressed
*     as a Modified Julian Date (MJD = JD - 2400000.5).
*
*     It is not necessary to reset the wcsprm struct (via wcsset()) when
*     wcsprm::mjdbeg is changed.
*
*   double mjdavg
*     (Given, auxiliary) MJD-AVG keyvalue, equivalent to DATE-AVG expressed
*     as a Modified Julian Date (MJD = JD - 2400000.5).
*
*     It is not necessary to reset the wcsprm struct (via wcsset()) when
*     wcsprm::mjdavg is changed.
*
*   double mjdend
*     (Given, auxiliary) MJD-END keyvalue, equivalent to DATE-END expressed
*     as a Modified Julian Date (MJD = JD - 2400000.5).
*
*     It is not necessary to reset the wcsprm struct (via wcsset()) when
*     wcsprm::mjdend is changed.
*
*   double jepoch
*     (Given, auxiliary) JEPOCH keyvalue, equivalent to DATE-OBS expressed
*     as a Julian epoch.
*
*     It is not necessary to reset the wcsprm struct (via wcsset()) when
*     wcsprm::jepoch is changed.
*
*   double bepoch
*     (Given, auxiliary) BEPOCH keyvalue, equivalent to DATE-OBS expressed
*     as a Besselian epoch
*
*     It is not necessary to reset the wcsprm struct (via wcsset()) when
*     wcsprm::bepoch is changed.
*
*   double tstart
*     (Given, auxiliary) TSTART keyvalue, equivalent to DATE-BEG expressed
*     as a time in units of TIMEUNIT relative to DATEREF+TIMEOFFS.
*
*     It is not necessary to reset the wcsprm struct (via wcsset()) when
*     wcsprm::tstart is changed.
*
*   double tstop
*     (Given, auxiliary) TSTOP keyvalue, equivalent to DATE-END expressed
*     as a time in units of TIMEUNIT relative to DATEREF+TIMEOFFS.
*
*     It is not necessary to reset the wcsprm struct (via wcsset()) when
*     wcsprm::tstop is changed.
*
*   double xposure
*     (Given, auxiliary) XPOSURE keyvalue, being the effective exposure time
*     in units of TIMEUNIT.
*
*     It is not necessary to reset the wcsprm struct (via wcsset()) when
*     wcsprm::xposure is changed.
*
*   double telapse
*     (Given, auxiliary) TELAPSE keyvalue, equivalent to the elapsed time
*     between DATE-BEG and DATE-END, in units of TIMEUNIT.
*
*     It is not necessary to reset the wcsprm struct (via wcsset()) when
*     wcsprm::telapse is changed.
*
*   double timsyer
*     (Given, auxiliary) TIMSYER keyvalue, being the absolute error of the
*     time values, in units of TIMEUNIT.
*
*     It is not necessary to reset the wcsprm struct (via wcsset()) when
*     wcsprm::timsyer is changed.
*
*   double timrder
*     (Given, auxiliary) TIMRDER keyvalue, being the accuracy of time stamps
*     relative to each other, in units of TIMEUNIT.
*
*     It is not necessary to reset the wcsprm struct (via wcsset()) when
*     wcsprm::timrder is changed.
*
*   double timedel
*     (Given, auxiliary) TIMEDEL keyvalue, being the resolution of the time
*     stamps.
*
*     It is not necessary to reset the wcsprm struct (via wcsset()) when
*     wcsprm::timedel is changed.
*
*   double timepixr
*     (Given, auxiliary) TIMEPIXR keyvalue, being the relative position of the
*     time stamps in binned time intervals, a value between 0.0 and 1.0.
*
*     It is not necessary to reset the wcsprm struct (via wcsset()) when
*     wcsprm::timepixr is changed.
*
*   double obsgeo[6]
*     (Given, auxiliary) Location of the observer in a standard terrestrial
*     reference frame.  The first three give ITRS Cartesian coordinates
*     OBSGEO-X [m],   OBSGEO-Y [m],   OBSGEO-Z [m], and the second three give
*     OBSGEO-L [deg], OBSGEO-B [deg], OBSGEO-H [m], which are related through
*     a standard transformation.
*
*     It is not necessary to reset the wcsprm struct (via wcsset()) when
*     wcsprm::obsgeo is changed.
*
*   char obsorbit[72]
*     (Given, auxiliary) OBSORBIT keyvalue, being the URI, URL, or name of an
*     orbit ephemeris file giving spacecraft coordinates relating to TREFPOS.
*
*     It is not necessary to reset the wcsprm struct (via wcsset()) when
*     wcsprm::obsorbit is changed.
*
*   char radesys[72]
*     (Given, auxiliary) The equatorial or ecliptic coordinate system type,
*     RADESYSa.
*
*     It is not necessary to reset the wcsprm struct (via wcsset()) when
*     wcsprm::radesys is changed.
*
*   double equinox
*     (Given, auxiliary) The equinox associated with dynamical equatorial or
*     ecliptic coordinate systems, EQUINOXa (or EPOCH in older headers).  Not
*     applicable to ICRS equatorial or ecliptic coordinates.
*
*     It is not necessary to reset the wcsprm struct (via wcsset()) when
*     wcsprm::equinox is changed.
*
*   char specsys[72]
*     (Given, auxiliary) Spectral reference frame (standard of rest),
*     SPECSYSa.
*
*     It is not necessary to reset the wcsprm struct (via wcsset()) when
*     wcsprm::specsys is changed.
*
*   char ssysobs[72]
*     (Given, auxiliary) The spectral reference frame in which there is no
*     differential variation in the spectral coordinate across the
*     field-of-view, SSYSOBSa.
*
*     It is not necessary to reset the wcsprm struct (via wcsset()) when
*     wcsprm::ssysobs is changed.
*
*   double velosys
*     (Given, auxiliary) The relative radial velocity [m/s] between the
*     observer and the selected standard of rest in the direction of the
*     celestial reference coordinate, VELOSYSa.
*
*     It is not necessary to reset the wcsprm struct (via wcsset()) when
*     wcsprm::velosys is changed.
*
*   double zsource
*     (Given, auxiliary) The redshift, ZSOURCEa, of the source.
*
*     It is not necessary to reset the wcsprm struct (via wcsset()) when
*     wcsprm::zsource is changed.
*
*   char ssyssrc[72]
*     (Given, auxiliary) The spectral reference frame (standard of rest),
*     SSYSSRCa, in which wcsprm::zsource was measured.
*
*     It is not necessary to reset the wcsprm struct (via wcsset()) when
*     wcsprm::ssyssrc is changed.
*
*   double velangl
*     (Given, auxiliary) The angle [deg] that should be used to decompose an
*     observed velocity into radial and transverse components.
*
*     It is not necessary to reset the wcsprm struct (via wcsset()) when
*     wcsprm::velangl is changed.
*
*   struct auxprm *aux
*     (Given, auxiliary) This struct holds auxiliary coordinate system
*     information of a specialist nature.  While these parameters may be
*     widely recognized within particular fields of astronomy, they differ
*     from the above auxiliary parameters in not being defined by any of the
*     FITS WCS standards.  Collecting them together in a separate struct that
*     is allocated only when required helps to control bloat in the size of
*     the wcsprm struct.
*
*   int ntab
*     (Given) See wcsprm::tab.
*
*   int nwtb
*     (Given) See wcsprm::wtb.
*
*   struct tabprm *tab
*     (Given) Address of the first element of an array of ntab tabprm structs
*     for which memory has been allocated.  These are used to store tabular
*     transformation parameters.
*
*     Although technically wcsprm::ntab and tab are "given", they will
*     normally be set by invoking wcstab(), whether directly or indirectly.
*
*     The tabprm structs contain some members that must be supplied and others
*     that are derived.  The information to be supplied comes primarily from
*     arrays stored in one or more FITS binary table extensions.  These
*     arrays, referred to here as "wcstab arrays", are themselves located by
*     parameters stored in the FITS image header.
*
*   struct wtbarr *wtb
*     (Given) Address of the first element of an array of nwtb wtbarr structs
*     for which memory has been allocated.  These are used in extracting
*     wcstab arrays from a FITS binary table.
*
*     Although technically wcsprm::nwtb and wtb are "given", they will
*     normally be set by invoking wcstab(), whether directly or indirectly.
*
*   char lngtyp[8]
*     (Returned) Four-character WCS celestial longitude and ...
*   char lattyp[8]
*     (Returned) ... latitude axis types. e.g. "RA", "DEC", "GLON", "GLAT",
*     etc. extracted from 'RA--', 'DEC-', 'GLON', 'GLAT', etc. in the first
*     four characters of CTYPEia but with trailing dashes removed.  (Declared
*     as char[8] for alignment reasons.)
*
*   int lng
*     (Returned) Index for the longitude coordinate, and ...
*   int lat
*     (Returned) ... index for the latitude coordinate, and ...
*   int spec
*     (Returned) ... index for the spectral coordinate in the imgcrd[][] and
*     world[][] arrays in the API of wcsp2s(), wcss2p() and wcsmix().
*
*     These may also serve as indices into the pixcrd[][] array provided that
*     the PCi_ja matrix does not transpose axes.
*
*   int cubeface
*     (Returned) Index into the pixcrd[][] array for the CUBEFACE axis.  This
*     is used for quadcube projections where the cube faces are stored on a
*     separate axis (see wcs.h).
*
*   int *types
*     (Returned) Address of the first element of an array of int containing a
*     four-digit type code for each axis.
*
*     - First digit (i.e. 1000s):
*       - 0: Non-specific coordinate type.
*       - 1: Stokes coordinate.
*       - 2: Celestial coordinate (including CUBEFACE).
*       - 3: Spectral coordinate.
*
*     - Second digit (i.e. 100s):
*       - 0: Linear axis.
*       - 1: Quantized axis (STOKES, CUBEFACE).
*       - 2: Non-linear celestial axis.
*       - 3: Non-linear spectral axis.
*       - 4: Logarithmic axis.
*       - 5: Tabular axis.
*
*     - Third digit (i.e. 10s):
*       - 0: Group number, e.g. lookup table number, being an index into the
*            tabprm array (see above).
*
*     - The fourth digit is used as a qualifier depending on the axis type.
*
*       - For celestial axes:
*         - 0: Longitude coordinate.
*         - 1: Latitude coordinate.
*         - 2: CUBEFACE number.
*
*       - For lookup tables: the axis number in a multidimensional table.
*
*     CTYPEia in "4-3" form with unrecognized algorithm code will have its
*     type set to -1 and generate an error.
*
*   struct linprm lin
*     (Returned) Linear transformation parameters (usage is described in the
*     prologue to lin.h).
*
*   struct celprm cel
*     (Returned) Celestial transformation parameters (usage is described in
*     the prologue to cel.h).
*
*   struct spcprm spc
*     (Returned) Spectral transformation parameters (usage is described in the
*     prologue to spc.h).
*
*   struct wcserr *err
*     (Returned) If enabled, when an error status is returned, this struct
*     contains detailed information about the error, see wcserr_enable().
*
*   int m_flag
*     (For internal use only.)
*   int m_naxis
*     (For internal use only.)
*   double *m_crpix
*     (For internal use only.)
*   double *m_pc
*     (For internal use only.)
*   double *m_cdelt
*     (For internal use only.)
*   double *m_crval
*     (For internal use only.)
*   char (*m_cunit)[72]
*     (For internal use only.)
*   char (*m_ctype)[72]
*     (For internal use only.)
*   struct pvcard *m_pv
*     (For internal use only.)
*   struct pscard *m_ps
*     (For internal use only.)
*   double *m_cd
*     (For internal use only.)
*   double *m_crota
*     (For internal use only.)
*   int *m_colax
*     (For internal use only.)
*   char (*m_cname)[72]
*     (For internal use only.)
*   double *m_crder
*     (For internal use only.)
*   double *m_csyer
*     (For internal use only.)
*   double *m_czphs
*     (For internal use only.)
*   double *m_cperi
*     (For internal use only.)
*   struct tabprm *m_tab
*     (For internal use only.)
*   struct wtbarr *m_wtb
*     (For internal use only.)
*
*
* pvcard struct - Store for PVi_ma keyrecords
* -------------------------------------------
* The pvcard struct is used to pass the parsed contents of PVi_ma keyrecords
* to wcsset() via the wcsprm struct.
*
* All members of this struct are to be set by the user.
*
*   int i
*     (Given) Axis number (1-relative), as in the FITS PVi_ma keyword.  If
*     i == 0, wcsset() will replace it with the latitude axis number.
*
*   int m
*     (Given) Parameter number (non-negative), as in the FITS PVi_ma keyword.
*
*   double value
*     (Given) Parameter value.
*
*
* pscard struct - Store for PSi_ma keyrecords
* -------------------------------------------
* The pscard struct is used to pass the parsed contents of PSi_ma keyrecords
* to wcsset() via the wcsprm struct.
*
* All members of this struct are to be set by the user.
*
*   int i
*     (Given) Axis number (1-relative), as in the FITS PSi_ma keyword.
*
*   int m
*     (Given) Parameter number (non-negative), as in the FITS PSi_ma keyword.
*
*   char value[72]
*     (Given) Parameter value.
*
*
* auxprm struct - Additional auxiliary parameters
* -----------------------------------------------
* The auxprm struct holds auxiliary coordinate system information of a
* specialist nature.  It is anticipated that this struct will expand in future
* to accomodate additional parameters.
*
* All members of this struct are to be set by the user.
*
*   double rsun_ref
*     (Given, auxiliary) Reference radius of the Sun used in coordinate
*     calculations (m).
*
*   double dsun_obs
*     (Given, auxiliary) Distance between the centre of the Sun and the
*     observer (m).
*
*   double crln_obs
*     (Given, auxiliary) Carrington heliographic longitude of the observer
*     (deg).
*
*   double hgln_obs
*     (Given, auxiliary) Stonyhurst heliographic longitude of the observer
*     (deg).
*
*   double hglt_obs
*     (Given, auxiliary) Heliographic latitude (Carrington or Stonyhurst) of
*     the observer (deg).
*
*
* Global variable: const char *wcs_errmsg[] - Status return messages
* ------------------------------------------------------------------
* Error messages to match the status value returned from each function.
*
*===========================================================================*/

#ifndef WCSLIB_WCS
#define WCSLIB_WCS

#include "lin.h"
#include "cel.h"
#include "spc.h"

#ifdef __cplusplus
extern "C" {
#define wtbarr wtbarr_s		/* See prologue of wtbarr.h.                */
#endif

#define WCSSUB_LONGITUDE 0x1001
#define WCSSUB_LATITUDE  0x1002
#define WCSSUB_CUBEFACE  0x1004
#define WCSSUB_CELESTIAL 0x1007
#define WCSSUB_SPECTRAL  0x1008
#define WCSSUB_STOKES    0x1010


#define WCSCOMPARE_ANCILLARY 0x0001
#define WCSCOMPARE_TILING    0x0002
#define WCSCOMPARE_CRPIX     0x0004


extern const char *wcs_errmsg[];

enum wcs_errmsg_enum {
  WCSERR_SUCCESS         =  0,	/* Success. */
  WCSERR_NULL_POINTER    =  1,	/* Null wcsprm pointer passed. */
  WCSERR_MEMORY          =  2,	/* Memory allocation failed. */
  WCSERR_SINGULAR_MTX    =  3,	/* Linear transformation matrix is
				   singular. */
  WCSERR_BAD_CTYPE       =  4,	/* Inconsistent or unrecognized coordinate
				   axis type. */
  WCSERR_BAD_PARAM       =  5,	/* Invalid parameter value. */
  WCSERR_BAD_COORD_TRANS =  6,	/* Unrecognized coordinate transformation
				   parameter. */
  WCSERR_ILL_COORD_TRANS =  7,	/* Ill-conditioned coordinate transformation
				   parameter. */
  WCSERR_BAD_PIX         =  8,	/* One or more of the pixel coordinates were
				   invalid. */
  WCSERR_BAD_WORLD       =  9,	/* One or more of the world coordinates were
				   invalid. */
  WCSERR_BAD_WORLD_COORD = 10,	/* Invalid world coordinate. */
  WCSERR_NO_SOLUTION     = 11,	/* No solution found in the specified
				   interval. */
  WCSERR_BAD_SUBIMAGE    = 12,	/* Invalid subimage specification. */
  WCSERR_NON_SEPARABLE   = 13 	/* Non-separable subimage coordinate
				   system. */
};


/* Struct used for storing PVi_ma keywords. */
struct pvcard {
  int i;			/* Axis number, as in PVi_ma (1-relative).  */
  int m;			/* Parameter number, ditto  (0-relative).   */
  double value;			/* Parameter value.                         */
};

/* Size of the pvcard struct in int units, used by the Fortran wrappers. */
#define PVLEN (sizeof(struct pvcard)/sizeof(int))

/* Struct used for storing PSi_ma keywords. */
struct pscard {
  int i;			/* Axis number, as in PSi_ma (1-relative).  */
  int m;			/* Parameter number, ditto  (0-relative).   */
  char value[72];		/* Parameter value.                         */
};

/* Size of the pscard struct in int units, used by the Fortran wrappers. */
#define PSLEN (sizeof(struct pscard)/sizeof(int))

/* Struct used to hold additional auxiliary parameters.                     */
struct auxprm {
  double rsun_ref;              /* Solar radius.                            */
  double dsun_obs;              /* Distance from Sun centre to observer.    */
  double crln_obs;              /* Carrington heliographic lng of observer. */
  double hgln_obs;              /* Stonyhurst heliographic lng of observer. */
  double hglt_obs;              /* Heliographic latitude of observer.       */
};

/* Size of the auxprm struct in int units, used by the Fortran wrappers. */
#define AUXLEN (sizeof(struct auxprm)/sizeof(int))


struct wcsprm {
  /* Initialization flag (see the prologue above).                          */
  /*------------------------------------------------------------------------*/
  int    flag;			/* Set to zero to force initialization.     */

  /* FITS header keyvalues to be provided (see the prologue above).         */
  /*------------------------------------------------------------------------*/
  int    naxis;			/* Number of axes (pixel and coordinate).   */
  double *crpix;		/* CRPIXja keyvalues for each pixel axis.   */
  double *pc;			/* PCi_ja  linear transformation matrix.    */
  double *cdelt;		/* CDELTia keyvalues for each coord axis.   */
  double *crval;		/* CRVALia keyvalues for each coord axis.   */

  char   (*cunit)[72];		/* CUNITia keyvalues for each coord axis.   */
  char   (*ctype)[72];		/* CTYPEia keyvalues for each coord axis.   */

  double lonpole;		/* LONPOLEa keyvalue.                       */
  double latpole;		/* LATPOLEa keyvalue.                       */

  double restfrq;		/* RESTFRQa keyvalue.                       */
  double restwav;		/* RESTWAVa keyvalue.                       */

  int    npv;			/* Number of PVi_ma keywords, and the       */
  int    npvmax;		/* number for which space was allocated.    */
  struct pvcard *pv;		/* PVi_ma keywords for each i and m.        */

  int    nps;			/* Number of PSi_ma keywords, and the       */
  int    npsmax;		/* number for which space was allocated.    */
  struct pscard *ps;		/* PSi_ma keywords for each i and m.        */

  /* Alternative header keyvalues (see the prologue above).                 */
  /*------------------------------------------------------------------------*/
  double *cd;			/* CDi_ja linear transformation matrix.     */
  double *crota;		/* CROTAi keyvalues for each coord axis.    */
  int    altlin;		/* Alternative representations              */
				/*   Bit 0: PCi_ja  is present,             */
				/*   Bit 1: CDi_ja  is present,             */
				/*   Bit 2: CROTAi is present.              */
  int    velref;		/* AIPS velocity code, VELREF.              */

  /* Auxiliary coordinate system information of a general nature.  Not      */
  /* used by WCSLIB.  Refer to the prologue comments above for a brief      */
  /* explanation of these values.                                           */
  char   alt[4];
  int    colnum;
  int    *colax;
				/* Auxiliary coordinate axis information.   */
  char   (*cname)[72];
  double *crder;
  double *csyer;
  double *czphs;
  double *cperi;

  char   wcsname[72];
				/* Time reference system and measurement.   */
  char   timesys[72], trefpos[72], trefdir[72], plephem[72];
  char   timeunit[72];
  char   dateref[72];
  double mjdref[2];
  double timeoffs;
				/* Data timestamps and durations.           */
  char   dateobs[72], datebeg[72], dateavg[72], dateend[72];
  double mjdobs, mjdbeg, mjdavg, mjdend;
  double jepoch, bepoch;
  double tstart, tstop;
  double xposure, telapse;
				/* Timing accuracy.                         */
  double timsyer, timrder;
  double timedel, timepixr;
				/* Spatial & celestial reference frame.     */
  double obsgeo[6];
  char   obsorbit[72];
  char   radesys[72];
  double equinox;
  char   specsys[72];
  char   ssysobs[72];
  double velosys;
  double zsource;
  char   ssyssrc[72];
  double velangl;

  /* Additional auxiliary coordinate system information of a specialist     */
  /* nature.  Not used by WCSLIB.  Refer to the prologue comments above.    */
  struct auxprm *aux;

  /* Coordinate lookup tables (see the prologue above).                     */
  /*------------------------------------------------------------------------*/
  int    ntab;			/* Number of separate tables.               */
  int    nwtb;			/* Number of wtbarr structs.                */
  struct tabprm *tab;		/* Tabular transformation parameters.       */
  struct wtbarr *wtb;		/* Array of wtbarr structs.                 */

  /*------------------------------------------------------------------------*/
  /* Information derived from the FITS header keyvalues by wcsset().        */
  /*------------------------------------------------------------------------*/
  char   lngtyp[8], lattyp[8];	/* Celestial axis types, e.g. RA, DEC.      */
  int    lng, lat, spec;	/* Longitude, latitude and spectral axis    */
				/* indices (0-relative).                    */
  int    cubeface;		/* True if there is a CUBEFACE axis.        */
  int    *types;		/* Coordinate type codes for each axis.     */

  struct linprm lin;		/*    Linear transformation parameters.     */
  struct celprm cel;		/* Celestial transformation parameters.     */
  struct spcprm spc;		/*  Spectral transformation parameters.     */

  /*------------------------------------------------------------------------*/
  /*             THE REMAINDER OF THE WCSPRM STRUCT IS PRIVATE.             */
  /*------------------------------------------------------------------------*/

  /* Error handling, if enabled.                                            */
  /*------------------------------------------------------------------------*/
  struct wcserr *err;

  /* Memory management.                                                     */
  /*------------------------------------------------------------------------*/
  int    m_flag, m_naxis;
  double *m_crpix, *m_pc, *m_cdelt, *m_crval;
  char  (*m_cunit)[72], (*m_ctype)[72];
  struct pvcard *m_pv;
  struct pscard *m_ps;
  double *m_cd, *m_crota;
  int    *m_colax;
  char  (*m_cname)[72];
  double *m_crder, *m_csyer, *m_czphs, *m_cperi;
  struct auxprm *m_aux;
  struct tabprm *m_tab;
  struct wtbarr *m_wtb;
};

/* Size of the wcsprm struct in int units, used by the Fortran wrappers. */
#define WCSLEN (sizeof(struct wcsprm)/sizeof(int))


int wcsnpv(int n);

int wcsnps(int n);

int wcsini(int alloc, int naxis, struct wcsprm *wcs);

int wcsinit(int alloc, int naxis, struct wcsprm *wcs, int npvmax, int npsmax,
            int ndpmax);

int wcsauxi(int alloc, struct wcsprm *wcs);

int wcssub(int alloc, const struct wcsprm *wcssrc, int *nsub, int axes[],
           struct wcsprm *wcsdst);

int wcscompare(int cmp, double tol, const struct wcsprm *wcs1,
               const struct wcsprm *wcs2, int *equal);

int wcsfree(struct wcsprm *wcs);

int wcsprt(const struct wcsprm *wcs);

int wcsperr(const struct wcsprm *wcs, const char *prefix);

int wcsbchk(struct wcsprm *wcs, int bounds);

int wcsset(struct wcsprm *wcs);

int wcsp2s(struct wcsprm *wcs, int ncoord, int nelem, const double pixcrd[],
           double imgcrd[], double phi[], double theta[], double world[],
           int stat[]);

int wcss2p(struct wcsprm *wcs, int ncoord, int nelem, const double world[],
           double phi[], double theta[], double imgcrd[], double pixcrd[],
           int stat[]);

int wcsmix(struct wcsprm *wcs, int mixpix, int mixcel, const double vspan[],
           double vstep, int viter, double world[], double phi[],
           double theta[], double imgcrd[], double pixcrd[]);

int wcssptr(struct wcsprm *wcs, int *i, char ctype[9]);

const char* wcslib_version(int vers[3]);

/* Defined mainly for backwards compatibility, use wcssub() instead. */
#define wcscopy(alloc, wcssrc, wcsdst) wcssub(alloc, wcssrc, 0x0, 0x0, wcsdst)


/* Deprecated. */
#define wcsini_errmsg wcs_errmsg
#define wcssub_errmsg wcs_errmsg
#define wcscopy_errmsg wcs_errmsg
#define wcsfree_errmsg wcs_errmsg
#define wcsprt_errmsg wcs_errmsg
#define wcsset_errmsg wcs_errmsg
#define wcsp2s_errmsg wcs_errmsg
#define wcss2p_errmsg wcs_errmsg
#define wcsmix_errmsg wcs_errmsg

#ifdef __cplusplus
#undef wtbarr
}
#endif

#endif /* WCSLIB_WCS */
