/*============================================================================

  WCSLIB 5.3 - an implementation of the FITS WCS standard.
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
  $Id: dis.h,v 5.3 2015/04/21 02:50:51 mcalabre Exp $
*=============================================================================
*
* WCSLIB 5.3 - experimental C routines that implement proposed extensions to
* the FITS World Coordinate System (WCS) standard.  Refer to
*
*   "Representations of distortions in FITS world coordinate systems",
*   Calabretta, M.R. et al. (Paper IV), draft dated 2004/04/22 available from
*   http://www.atnf.csiro.au/people/Mark.Calabretta
*
* The distortion function component of the TPV "projection" is also supported.
* The TPV projection, originally proposed in a draft of WCS Paper II, consists
* of a TAN projection with sequent polynomial distortion, the coefficients of
* which are encoded in PVi_ma keyrecords.  Full details may be found at the
* registry of FITS conventions:
*
*   http://fits.gsfc.nasa.gov/registry/tpvwcs/tpv.html
*
* Refer to the README file provided with WCSLIB for an overview of the
* library.
*
*
* Summary of the dis routines
* ---------------------------
* These routines apply the distortion functions defined by the extension to
* the FITS WCS standard proposed in Paper IV.  They are based on the disprm
* struct which contains all information needed for the computations.  The
* struct contains some members that must be set by the user, and others that
* are maintained by these routines, somewhat like a C++ class but with no
* encapsulation.
*
* disndp(), dpfill(), disini(), discpy(), and disfree() are provided to manage
* the disprm struct, and another, disprt(), prints its contents.
*
* A setup routine, disset(), computes intermediate values in the disprm struct
* from parameters in it that were supplied by the user.  The struct always
* needs to be set up by disset(), though disset() need not be called
* explicitly - refer to the explanation of disprm::flag.
*
* disp2x() and disx2p() implement the WCS distortion functions, disp2x() using
* separate functions, such as dispoly() and tpv7(), to do the computation.
*
* An auxiliary routine, diswarp(), computes various measures of the distortion
* over a specified range of coordinates.
*
* PLEASE NOTE: Distortions are not currently handled by wcspih(), wcsbth(),
*              wcssub(), wcscompare(), or wcshdo().
*
*
* disndp() - Memory allocation for DPja and DQia
* ----------------------------------------------
* disndp() changes the value of NDPMAX (default 256).  This global variable
* controls the number of dpkey structs, for holding DPja or DQia keyvalues,
* that disini() should allocate space for.
*
* PLEASE NOTE: This function is not thread-safe.
*
* Given:
*   n         int       Value of NDPMAX; ignored if < 0.
*
* Function return value:
*             int       Current value of NDPMAX.
*
*
* dpfill() - Fill the contents of a dpkey struct
* ----------------------------------------------
* dpfill() is a utility routine to aid in filling the contents of the dpkey
* struct.  No checks are done on the validity of the inputs.
*
* WCS Paper IV specifies the syntax of a record-valued keyword as
*
=   keyword = '<field-specifier>: <float>'
* 
* However, some DPja and DQia record values, such as those of DPja.NAXES and
* DPja.AXIS.j, are intrinsically integer-valued.  While FITS header parsers
* are not expected to know in advance which of DPja and DQia are integral and
* which are floating point, if the record's value parses as an integer (i.e.
* without decimal point or exponent), then preferably enter it into the dpkey
* struct as an integer.  Either way, it doesn't matter as disset() accepts
* either data type for all record values.
*
* Given and returned:
*   dp        struct dpkey*
*                       Store for DPja and DQia keyvalues.
*
* Given:
*   keyword   const char *
*   field     const char *
*                       These arguments are concatenated with an intervening
*                       "." to construct the full record field name, i.e.
*                       including the keyword name, DPja or DQia (but
*                       excluding the colon delimiter which is NOT part of the
*                       name).  Either may be given as a NULL pointer.  Set
*                       both NULL to omit setting this component of the
*                       struct.
*
*   j         int       Axis number (1-relative), i.e. the j in DPja or
*                       i in DQia.  Can be given as 0, in which case the axis
*                       number will be obtained from the keyword component of
*                       the field name which must either have been given or
*                       preset.
*
*   type      int       Data type of the record's value
*                         0: Integer,
*                         1: Floating point.
*
*   ival      int       For type == 0, the integer value of the record.
*
*   fval      double    For type == 1, the floating point value of the record.
*
* Function return value:
*             int       Status return value:
*                         0: Success.
*
*
* disini() - Default constructor for the disprm struct
* ----------------------------------------------------
* disini() allocates memory for arrays in a disprm struct and sets all members
* of the struct to default values.  Memory is allocated for up to NDPMAX DPja
* or DQia keywords per WCS representation.  This may be changed via disndp()
* before disini() is called.
*
* PLEASE NOTE: every disprm struct must be initialized by disini(), possibly
* repeatedly.  On the first invokation, and only the first invokation,
* disprm::flag must be set to -1 to initialize memory management, regardless
* of whether disini() will actually be used to allocate memory.
*
* Given:
*   alloc     int       If true, allocate memory unconditionally for arrays in
*                       the disprm struct.
*
*                       If false, it is assumed that pointers to these arrays
*                       have been set by the user except if they are null
*                       pointers in which case memory will be allocated for
*                       them regardless.  (In other words, setting alloc true
*                       saves having to initalize these pointers to zero.)
*
*   naxis     int       The number of world coordinate axes, used to determine
*                       array sizes.
*
* Given and returned:
*   dis       struct disprm*
*                       Distortion function parameters.  Note that, in order
*                       to initialize memory management disprm::flag must be
*                       set to -1 when dis is initialized for the first time
*                       (memory leaks may result if it had already been
*                       initialized).
*
* Function return value:
*             int       Status return value:
*                         0: Success.
*                         1: Null disprm pointer passed.
*                         2: Memory allocation failed.
*
*                       For returns > 1, a detailed error message is set in
*                       disprm::err if enabled, see wcserr_enable().
*
*
* discpy() - Copy routine for the disprm struct
* ---------------------------------------------
* discpy() does a deep copy of one disprm struct to another, using disini() to
* allocate memory unconditionally for its arrays if required.  Only the
* "information to be provided" part of the struct is copied; a call to
* disset() is required to initialize the remainder.
*
* Given:
*   alloc     int       If true, allocate memory unconditionally for arrays in
*                       the destination.  Otherwise, it is assumed that
*                       pointers to these arrays have been set by the user
*                       except if they are null pointers in which case memory
*                       will be allocated for them regardless.
*
*   dissrc    const struct disprm*
*                       Struct to copy from.
*
* Given and returned:
*   disdst    struct disprm*
*                       Struct to copy to.  disprm::flag should be set to -1
*                       if disdst was not previously initialized (memory leaks
*                       may result if it was previously initialized).
*
* Function return value:
*             int       Status return value:
*                         0: Success.
*                         1: Null disprm pointer passed.
*                         2: Memory allocation failed.
*
*                       For returns > 1, a detailed error message is set in
*                       disprm::err if enabled, see wcserr_enable().
*
*
* disfree() - Destructor for the disprm struct
* --------------------------------------------
* disfree() frees memory allocated for the disprm arrays by disini().
* disini() keeps a record of the memory it allocates and disfree() will only
* attempt to free this.
*
* PLEASE NOTE: disfree() must not be invoked on a disprm struct that was not
* initialized by disini().
*
* Given:
*   dis       struct disprm*
*                       Distortion function parameters.
*
* Function return value:
*             int       Status return value:
*                         0: Success.
*                         1: Null disprm pointer passed.
*
*
* disprt() - Print routine for the disprm struct
* ----------------------------------------------
* disprt() prints the contents of a disprm struct using wcsprintf().  Mainly
* intended for diagnostic purposes.
*
* Given:
*   dis       const struct disprm*
*                       Distortion function parameters.
*
* Function return value:
*             int       Status return value:
*                         0: Success.
*                         1: Null disprm pointer passed.
*
*
* disset() - Setup routine for the disprm struct
* ----------------------------------------------
* disset(), sets up the disprm struct according to information supplied within
* it - refer to the explanation of disprm::flag.
*
* Note that this routine need not be called directly; it will be invoked by
* disp2x() and disx2p() if the disprm::flag is anything other than a
* predefined magic value.
*
* Given and returned:
*   dis       struct disprm*
*                       Distortion function parameters.
*
* Function return value:
*             int       Status return value:
*                         0: Success.
*                         1: Null disprm pointer passed.
*                         2: Memory allocation failed.
*                         3: Invalid parameter.
*
*                       For returns > 1, a detailed error message is set in
*                       disprm::err if enabled, see wcserr_enable().
*
*
* disp2x() - Apply distortion function
* ------------------------------------
* disp2x() applies the distortion functions.  By definition, the distortion
* is in the pixel-to-world direction.
*
* Depending on the point in the algorithm chain at which it is invoked,
* disp2x() may transform pixel coordinates to corrected pixel coordinates, or
* intermediate pixel coordinates to corrected intermediate pixel coordinates,
* or image coordinates to corrected image coordinates.
*
*
* Given and returned:
*   dis       struct disprm*
*                       Distortion function parameters.
*
* Given:
*   rawcrd    const double[naxis]
*                       Array of coordinates.
*
* Returned:
*   discrd    double[naxis]
*                       Array of coordinates to which the distortion functions
*                       have been applied.
*
* Function return value:
*             int       Status return value:
*                         0: Success.
*                         1: Null disprm pointer passed.
*                         2: Memory allocation failed.
*                         3: Invalid parameter.
*                         4: Distort error.
*
*                       For returns > 1, a detailed error message is set in
*                       disprm::err if enabled, see wcserr_enable().
*
*
* disx2p() - Apply de-distortion function
* ---------------------------------------
* disx2p() applies the inverse of the distortion functions.  By definition,
* the de-distortion is in the world-to-pixel direction.
*
* Depending on the point in the algorithm chain at which it is invoked,
* disx2p() may transform corrected pixel coordinates to pixel coordinates, or
* corrected intermediate pixel coordinates to intermediate pixel coordinates,
* or corrected image coordinates to image coordinates.
*
* disx2p() iteratively solves for the inverse using disp2x().  It assumes
* that the distortion is small and the functions are well-behaved, being
* continuous and with continuous derivatives.  Also that, to first order
* in the neighbourhood of the solution, discrd[j] ~= a + b*rawcrd[j], i.e.
* independent of rawcrd[i], where i != j.  This is effectively equivalent to
* assuming that the distortion functions are separable to first order.
* Furthermore, a is assumed to be small, and b close to unity.
*
* Given and returned:
*   dis       struct disprm*
*                       Distortion function parameters.
*
* Given:
*   discrd    const double[naxis]
*                       Array of coordinates.
*
* Returned:
*   rawcrd    double[naxis]
*                       Array of coordinates to which the inverse distortion
*                       functions have been applied.
*
* Function return value:
*             int       Status return value:
*                         0: Success.
*                         1: Null disprm pointer passed.
*                         2: Memory allocation failed.
*                         3: Invalid parameter.
*                         5: De-distort error.
*
*                       For returns > 1, a detailed error message is set in
*                       disprm::err if enabled, see wcserr_enable().
*
*
* diswarp() - Compute measures of distortion
* ------------------------------------------
* diswarp() computes various measures of the distortion over a specified range
* of coordinates.
*
* For prior distortions, the measures may be interpreted simply as an offset
* in pixel coordinates.  For sequent distortions, the interpretation depends
* on the nature of the linear transformation matrix (PCi_ja or CDi_ja).  If
* the latter introduces a scaling, then the measures will also be scaled.
* Note also that the image domain, which is rectangular in pixel coordinates,
* may be rotated, skewed, and/or stretched in intermediate pixel coordinates,
* and in general cannot be defined using pixblc[] and pixtrc[].
*
* PLEASE NOTE: the measures of total distortion may be essentially meaningless
* if there are multiple sequent distortions with different scaling.
*
* See also linwarp().
*
* Given and returned:
*   dis       struct disprm*
*                       Distortion function parameters.
*
* Given:
*   pixblc    const double[naxis]
*                       Start of the range of pixel coordinates (for prior
*                       distortions), or intermediate pixel coordinates (for
*                       sequent distortions).  May be specified as a NULL
*                       pointer which is interpreted as (1,1,...).
*
*   pixtrc    const double[naxis]
*                       End of the range of pixel coordinates (prior) or
*                       intermediate pixel coordinates (sequent).
*
*   pixsamp   const double[naxis]
*                       If positive or zero, the increment on the particular
*                       axis, starting at pixblc[].  Zero is interpreted as a
*                       unit increment.  pixsamp may also be specified as a
*                       NULL pointer which is interpreted as all zeroes, i.e.
*                       unit increments on all axes.
*
*                       If negative, the grid size on the particular axis (the
*                       absolute value being rounded to the nearest integer).
*                       For example, if pixsamp is (-128.0,-128.0,...) then
*                       each axis will be sampled at 128 points between
*                       pixblc[] and pixtrc[] inclusive.  Use caution when
*                       using this option on non-square images.
*
* Returned:
*   nsamp     int*      The number of pixel coordinates sampled.
*
*                       Can be specified as a NULL pointer if not required.
*
*   maxdis    double[naxis]
*                       For each individual distortion function, the
*                       maximum absolute value of the distortion.
*
*                       Can be specified as a NULL pointer if not required.
*
*   maxtot    double*   For the combination of all distortion functions, the
*                       maximum absolute value of the distortion.
*
*                       Can be specified as a NULL pointer if not required.
*
*   avgdis    double[naxis]
*                       For each individual distortion function, the
*                       mean value of the distortion.
*
*                       Can be specified as a NULL pointer if not required.
*
*   avgtot    double*   For the combination of all distortion functions, the
*                       mean value of the distortion.
*
*                       Can be specified as a NULL pointer if not required.
*
*   rmsdis    double[naxis]
*                       For each individual distortion function, the
*                       root mean square deviation of the distortion.
*
*                       Can be specified as a NULL pointer if not required.
*
*   rmstot    double*   For the combination of all distortion functions, the
*                       root mean square deviation of the distortion.
*
*                       Can be specified as a NULL pointer if not required.
*
* Function return value:
*             int       Status return value:
*                         0: Success.
*                         1: Null disprm pointer passed.
*                         2: Memory allocation failed.
*                         3: Invalid parameter.
*                         4: Distort error.
*
*
* disprm struct - Distortion parameters
* -------------------------------------
* The disprm struct contains all of the information required to apply a set of
* distortion functions.  It consists of certain members that must be set by
* the user ("given") and others that are set by the WCSLIB routines
* ("returned").  While the addresses of the arrays themselves may be set by
* disini() if it (optionally) allocates memory, their contents must be set by
* the user.
*
*   int flag
*     (Given and returned) This flag must be set to zero whenever any of the
*     following members of the disprm struct are set or modified:
*
*       - disprm::naxis,
*       - disprm::dtype,
*       - disprm::ndp,
*       - disprm::dp.
*
*     This signals the initialization routine, disset(), to recompute the
*     returned members of the disprm struct.  disset() will reset flag to
*     indicate that this has been done.
*
*     PLEASE NOTE: flag must be set to -1 when disini() is called for the
*     first time for a particular disprm struct in order to initialize memory
*     management.  It must ONLY be used on the first initialization otherwise
*     memory leaks may result.
*
*   int naxis
*     (Given or returned) Number of pixel and world coordinate elements.
*
*     If disini() is used to initialize the disprm struct (as would normally
*     be the case) then it will set naxis from the value passed to it as a
*     function argument.  The user should not subsequently modify it.
*
*   char (*dtype)[72]
*     (Given) Pointer to the first element of an array of char[72] containing
*     the name of the distortion function for each axis.
*
*   int ndp
*     (Given) The number of entries in the disprm::dp[] array.
*
*   int ndpmax
*     (Given) The length of the disprm::dp[] array.
*
*     ndpmax will be set by disini() if it allocates memory for disprm::dp[],
*     otherwise it must be set by the user.  See also disndp().
*
*   struct dpkey dp
*     (Given) Address of the first element of an array of length ndpmax of
*     dpkey structs.
*
*     As a FITS header parser encounters each DPja or DQia keyword it should
*     load it into a dpkey struct in the array and increment ndp.  However,
*     note that a single disprm struct must hold only DPja or DQia keyvalues,
*     not both.  disset() interprets them as required by the particular
*     distortion function.
*
*   double *maxdis
*     (Given) Pointer to the first element of an array of double specifying
*     the maximum absolute value of the distortion for each axis computed over
*     the whole image.
*
*     It is not necessary to reset the disprm struct (via disset()) when
*     disprm::maxdis is changed.
*
*   double totdis
*     (Given) The maximum absolute value of the combination of all distortion
*     functions specified as an offset in pixel coordinates computed over the
*     whole image.
*
*     It is not necessary to reset the disprm struct (via disset()) when
*     disprm::totdis is changed.
*
*   int **axmap
*     (Returned) Pointer to the first element of an array of int* containing
*     pointers to the first elements of the axis mapping arrays for each axis.
*
*     An axis mapping associates the independent variables of a distortion
*     function with the 1-relative image axis number.  For example, consider
*     an image with a spectrum on the first axis, followed by RA, Dec, and
*     time axes.  For a distortion in (RA,Dec) and no distortion on the
*     spectral or time axes, the axis mapping arrays, axmap[j][], would be
*
=       j=0: [0, 0, 0, 0]   ...no  distortion on spectral axis,
=         1: [2, 3, 0, 0]   ...RA  distortion depends on RA and Dec,
=         2: [3, 2, 0, 0]   ...Dec distortion depends on Dec and RA,
=         3: [0, 0, 0, 0]   ...no  distortion on time axis,
*
*     where zero indicates that there is no corresponding independent
*     variable.
*
*   int *Nhat
*     (Returned) Pointer to the first element of an array of int* containing
*     the number of coordinate axes that form the independent variables of the
*     distortion function.
*
*   double **offset
*     (Returned) Pointer to the first element of an array of double*
*     containing an offset used to renormalize the independent variables of
*     the distortion function for each axis.
*
*     The offsets are subtracted from the independent variables before
*     scaling.
*
*   double **scale
*     (Returned) Pointer to the first element of an array of double*
*     containing a scale used to renormalize the independent variables of the
*     distortion function for each axis.
*
*     The scale is applied to the independent variables after the offsets are
*     subtracted.
*
*   int **iparm
*     (Returned) Pointer to the first element of an array of int*
*     containing pointers to the first elements of the arrays of integer
*     distortion parameters for each axis.
*
*   double **dparm
*     (Returned) Pointer to the first element of an array of double*
*     containing pointers to the first elements of the arrays of floating
*     point distortion parameters for each axis.
*
*   int i_naxis
*     (Returned) Dimension of the internal arrays (normally equal to naxis).
*
*   int ndis
*     (Returned) The number of distortion functions.
*
*   struct wcserr *err
*     (Returned) If enabled, when an error status is returned, this struct
*     contains detailed information about the error, see wcserr_enable().
*
*   int (**disp2x)(DISP2X_ARGS)
*     (For internal use only.)
*   int (**disx2p)(DISX2P_ARGS)
*     (For internal use only.)
*   double *tmpmem
*     (For internal use only.)
*   int m_flag
*     (For internal use only.)
*   int m_naxis
*     (For internal use only.)
*   char (*m_dtype)[72]
*     (For internal use only.)
*   double **m_dp
*     (For internal use only.)
*   double *m_maxdis
*     (For internal use only.)
*
*
* dpkey struct - Store for DPja and DQia keyvalues
* ------------------------------------------------
* The dpkey struct is used to pass the parsed contents of DPja or DQia
* keyrecords to disset() via the disprm struct.  A disprm struct must hold
* only DPja or DQia keyvalues, not both.
*
* All members of this struct are to be set by the user.
*
*   char field[72]
*     (Given) The full field name of the record, including the keyword name.
*     Note that the colon delimiter separating the field name and the value in
*     record-valued keyvalues is not part of the field name.  For example, in
*     the following:
*
=       DP3A = 'AXIS.1: 2'
*
*     the full record field name is "DP3A.AXIS.1", and the record's value
*     is 2.
*
*   int j
*     (Given) Axis number (1-relative), i.e. the j in DPja or i in DQia.
*
*   int type
*     (Given) The data type of the record's value
*       - 0: Integer (stored as an int),
*       - 1: Floating point (stored as a double).
*
*   union value
*     (Given) A union comprised of
*       - dpkey::ival,
*       - dpkey::fval,
*
*     the record's value.
*
*
* Global variable: const char *dis_errmsg[] - Status return messages
* ------------------------------------------------------------------
* Error messages to match the status value returned from each function.
*
*===========================================================================*/

#ifndef WCSLIB_DIS
#define WCSLIB_DIS

#ifdef __cplusplus
extern "C" {
#endif


extern const char *dis_errmsg[];

enum dis_errmsg_enum {
  DISERR_SUCCESS      = 0,	/* Success. */
  DISERR_NULL_POINTER = 1,	/* Null disprm pointer passed. */
  DISERR_MEMORY       = 2,	/* Memory allocation failed. */
  DISERR_BAD_PARAM    = 3,	/* Invalid parameter value. */
  DISERR_DISTORT      = 4,	/* Distortion error. */
  DISERR_DEDISTORT    = 5	/* De-distortion error. */
};

/* For use in declaring distortion function prototypes (= DISX2P_ARGS). */
#define DISP2X_ARGS const int iparm[], const double dparm[], \
int ncrd, const double rawcrd[], double *discrd

/* For use in declaring de-distortion function prototypes (= DISP2X_ARGS). */
#define DISX2P_ARGS const int iparm[], const double dparm[], \
int ncrd, const double discrd[], double *rawcrd


/* Struct used for storing DPja and DQia keyvalues. */
struct dpkey {
  char field[72];		/* Full record field name (no colon).       */
  int j;			/* Axis number, as in DPja (1-relative).    */
  int type;			/* Data type of value.                      */
  union {
    int    i;			/* Integer record value.                    */
    double f;			/* Floating point record value.             */
  } value;			/* Record value.                            */
};

/* Size of the dpkey struct in int units, used by the Fortran wrappers. */
#define DPLEN (sizeof(struct dpkey)/sizeof(int))


struct disprm {
  /* Initialization flag (see the prologue above).                          */
  /*------------------------------------------------------------------------*/
  int flag;			/* Set to zero to force initialization.     */

  /* Parameters to be provided (see the prologue above).                    */
  /*------------------------------------------------------------------------*/
  int naxis;			/* The number of pixel coordinate elements, */
				/* given by NAXIS.                          */
  char   (*dtype)[72];		/* For each axis, the distortion type.      */
  int    ndp;			/* Number of DPja or DQia keywords, and the */
  int    ndpmax;		/* number for which space was allocated.    */
  struct dpkey *dp;		/* DPja or DQia keyvalues (not both).       */
  double *maxdis;		/* For each axis, the maximum distortion.   */
  double totdis;		/* The maximum combined distortion.         */

  /* Information derived from the parameters supplied.                      */
  /*------------------------------------------------------------------------*/
  int    **axmap;		/* For each axis, the axis mapping array.   */
  int    *Nhat;			/* For each axis, the number of coordinate  */
				/* axes that form the independent variables */
				/* of the distortion function.              */
  double **offset;		/* For each axis, renormalization offsets.  */
  double **scale;		/* For each axis, renormalization scales.   */
  int    **iparm;		/* For each axis, the array of integer      */
				/* distortion parameters.                   */
  double **dparm;		/* For each axis, the array of floating     */
				/* point distortion parameters.             */
  int    i_naxis;		/* Dimension of the internal arrays.        */
  int    ndis;			/* The number of distortion functions.      */

  /* Error handling, if enabled.                                            */
  /*------------------------------------------------------------------------*/
  struct wcserr *err;

  /* Private - the remainder are for internal use.                          */
  /*------------------------------------------------------------------------*/
  int (**disp2x)(DISP2X_ARGS);	/* For each axis, pointers to the           */
  int (**disx2p)(DISX2P_ARGS);	/* distortion function and its inverse.     */

  double *tmpmem;

  int    m_flag, m_naxis;	/* The remainder are for memory management. */
  char   (*m_dtype)[72];
  struct dpkey *m_dp;
  double *m_maxdis;
};

/* Size of the disprm struct in int units, used by the Fortran wrappers. */
#define DISLEN (sizeof(struct disprm)/sizeof(int))


int disndp(int n);

int dpfill(struct dpkey *dp, const char *keyword, const char *field, int j,
           int type, int ival, double fval);

int disini(int alloc, int naxis, struct disprm *dis);

int discpy(int alloc, const struct disprm *dissrc, struct disprm *disdst);

int disfree(struct disprm *dis);

int disprt(const struct disprm *dis);

int disset(struct disprm *dis);

int disp2x(struct disprm *dis, const double rawcrd[], double discrd[]);

int disx2p(struct disprm *dis, const double discrd[], double rawcrd[]);

int diswarp(struct disprm *dis, const double pixblc[], const double pixtrc[],
            const double pixsamp[], int *nsamp,
            double maxdis[], double *maxtot,
            double avgdis[], double *avgtot,
            double rmsdis[], double *rmstot);


/* Specialist distortion functions (internal use only). */
int polyset(int j, struct disprm *dis);
int dispoly(DISP2X_ARGS);

int tpvset(int j, struct disprm *dis);

int tpv1(DISP2X_ARGS);
int tpv2(DISP2X_ARGS);
int tpv3(DISP2X_ARGS);
int tpv4(DISP2X_ARGS);
int tpv5(DISP2X_ARGS);
int tpv6(DISP2X_ARGS);
int tpv7(DISP2X_ARGS);

#ifdef __cplusplus
}
#endif

#endif /* WCSLIB_DIS */
