/*============================================================================

  WCSLIB 4.25 - an implementation of the FITS WCS standard.
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
  $Id: tab.h,v 4.25 2014/12/14 14:29:36 mcalabre Exp $
*=============================================================================
*
* WCSLIB 4.25 - C routines that implement tabular coordinate systems as
* defined by the FITS World Coordinate System (WCS) standard.  Refer to
*
*   "Representations of world coordinates in FITS",
*   Greisen, E.W., & Calabretta, M.R. 2002, A&A, 395, 1061 (paper I)
*
*   "Representations of spectral coordinates in FITS",
*   Greisen, E.W., Calabretta, M.R., Valdes, F.G., & Allen, S.L.
*   2006, A&A, 446, 747 (Paper III)
*
* Refer to the README file provided with WCSLIB for an overview of the
* library.
*
*
* Summary of the tab routines
* ---------------------------
* These routines implement the part of the FITS WCS standard that deals with
* tabular coordinates, i.e. coordinates that are defined via a lookup table.
* They define methods to be used for computing tabular world coordinates from
* intermediate world coordinates (a linear transformation of image pixel
* coordinates), and vice versa.  They are based on the tabprm struct which
* contains all information needed for the computations.  The struct contains
* some members that must be set by the user, and others that are maintained
* by these routines, somewhat like a C++ class but with no encapsulation.
*
* tabini(), tabmem(), tabcpy(), and tabfree() are provided to manage the
* tabprm struct, and another, tabprt(), to print its contents.
*
* A setup routine, tabset(), computes intermediate values in the tabprm struct
* from parameters in it that were supplied by the user.  The struct always
* needs to be set up by tabset() but it need not be called explicitly - refer
* to the explanation of tabprm::flag.
*
* tabx2s() and tabs2x() implement the WCS tabular coordinate transformations.
*
* Accuracy:
* ---------
* No warranty is given for the accuracy of these routines (refer to the
* copyright notice); intending users must satisfy for themselves their
* adequacy for the intended purpose.  However, closure effectively to within
* double precision rounding error was demonstrated by test routine ttab.c
* which accompanies this software.
*
*
* tabini() - Default constructor for the tabprm struct
* ----------------------------------------------------
* tabini() allocates memory for arrays in a tabprm struct and sets all members
* of the struct to default values.
*
* PLEASE NOTE: every tabprm struct should be initialized by tabini(), possibly
* repeatedly.  On the first invokation, and only the first invokation, the
* flag member of the tabprm struct must be set to -1 to initialize memory
* management, regardless of whether tabini() will actually be used to allocate
* memory.
*
* Given:
*   alloc     int       If true, allocate memory unconditionally for arrays in
*                       the tabprm struct.
*
*                       If false, it is assumed that pointers to these arrays
*                       have been set by the user except if they are null
*                       pointers in which case memory will be allocated for
*                       them regardless.  (In other words, setting alloc true
*                       saves having to initalize these pointers to zero.)
*
*   M         int       The number of tabular coordinate axes.
*
*   K         const int[]
*                       Vector of length M whose elements (K_1, K_2,... K_M)
*                       record the lengths of the axes of the coordinate array
*                       and of each indexing vector.  M and K[] are used to
*                       determine the length of the various tabprm arrays and
*                       therefore the amount of memory to allocate for them.
*                       Their values are copied into the tabprm struct.
*
*                       It is permissible to set K (i.e. the address of the
*                       array) to zero which has the same effect as setting
*                       each element of K[] to zero.  In this case no memory
*                       will be allocated for the index vectors or coordinate
*                       array in the tabprm struct.  These together with the
*                       K vector must be set separately before calling
*                       tabset().
*
* Given and returned:
*   tab       struct tabprm*
*                       Tabular transformation parameters.  Note that, in
*                       order to initialize memory management tabprm::flag
*                       should be set to -1 when tab is initialized for the
*                       first time (memory leaks may result if it had already
*                       been initialized).
*
* Function return value:
*             int       Status return value:
*                         0: Success.
*                         1: Null tabprm pointer passed.
*                         2: Memory allocation failed.
*                         3: Invalid tabular parameters.
*
*                       For returns > 1, a detailed error message is set in
*                       tabprm::err if enabled, see wcserr_enable().
*
*
* tabmem() - Acquire tabular memory
* ---------------------------------
* tabmem() takes control of memory allocated by the user for arrays in the
* tabprm struct.
*
* Given and returned:
*   tab       struct tabprm*
*                       Tabular transformation parameters.
*
* Function return value:
*             int       Status return value:
*                         0: Success.
*                         1: Null tabprm pointer passed.
*                         2: Memory allocation failed.
*
*                       For returns > 1, a detailed error message is set in
*                       tabprm::err if enabled, see wcserr_enable().
*
*
* tabcpy() - Copy routine for the tabprm struct
* ---------------------------------------------
* tabcpy() does a deep copy of one tabprm struct to another, using tabini() to
* allocate memory for its arrays if required.  Only the "information to be
* provided" part of the struct is copied; a call to tabset() is required to
* set up the remainder.
*
* Given:
*   alloc     int       If true, allocate memory unconditionally for arrays in
*                       the tabprm struct.
*
*                       If false, it is assumed that pointers to these arrays
*                       have been set by the user except if they are null
*                       pointers in which case memory will be allocated for
*                       them regardless.  (In other words, setting alloc true
*                       saves having to initalize these pointers to zero.)
*
*   tabsrc    const struct tabprm*
*                       Struct to copy from.
*
* Given and returned:
*   tabdst    struct tabprm*
*                       Struct to copy to.  tabprm::flag should be set to -1
*                       if tabdst was not previously initialized (memory leaks
*                       may result if it was previously initialized).
*
* Function return value:
*             int       Status return value:
*                         0: Success.
*                         1: Null tabprm pointer passed.
*                         2: Memory allocation failed.
*
*                       For returns > 1, a detailed error message is set in
*                       tabprm::err (associated with tabdst) if enabled, see
*                       wcserr_enable().
*
*
* tabcmp() - Compare two tabprm structs for equality
* --------------------------------------------------
* tabcmp() compares two tabprm structs for equality.
*
* Given:
*   cmp       int       A bit field controlling the strictness of the
*                       comparison.  At present, this value must always be 0,
*                       indicating a strict comparison.  In the future, other
*                       options may be added.
*
*   tol       double    Tolerance for comparison of floating-point values.
*                       For example, for tol == 1e-6, all floating-point
*                       values in the structs must be equal to the first 6
*                       decimal places.  A value of 0 implies exact equality.
*
*   tab1      const struct tabprm*
*                       The first tabprm struct to compare.
*
*   tab2      const struct tabprm*
*                       The second tabprm struct to compare.
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
* tabfree() - Destructor for the tabprm struct
* --------------------------------------------
* tabfree() frees memory allocated for the tabprm arrays by tabini().
* tabini() records the memory it allocates and tabfree() will only attempt to
* free this.
*
* PLEASE NOTE: tabfree() must not be invoked on a tabprm struct that was not
* initialized by tabini().
*
* Returned:
*   tab       struct tabprm*
*                       Coordinate transformation parameters.
*
* Function return value:
*             int       Status return value:
*                         0: Success.
*                         1: Null tabprm pointer passed.
*
*
* tabprt() - Print routine for the tabprm struct
* ----------------------------------------------
* tabprt() prints the contents of a tabprm struct using wcsprintf().  Mainly
* intended for diagnostic purposes.
*
* Given:
*   tab       const struct tabprm*
*                       Tabular transformation parameters.
*
* Function return value:
*             int       Status return value:
*                         0: Success.
*                         1: Null tabprm pointer passed.
*
*
* tabset() - Setup routine for the tabprm struct
* -----------------------------------------------
* tabset() allocates memory for work arrays in the tabprm struct and sets up
* the struct according to information supplied within it.
*
* Note that this routine need not be called directly; it will be invoked by
* tabx2s() and tabs2x() if tabprm::flag is anything other than a predefined
* magic value.
*
* Given and returned:
*   tab       struct tabprm*
*                       Tabular transformation parameters.
*
* Function return value:
*             int       Status return value:
*                         0: Success.
*                         1: Null tabprm pointer passed.
*                         3: Invalid tabular parameters.
*
*                       For returns > 1, a detailed error message is set in
*                       tabprm::err if enabled, see wcserr_enable().
*
*
* tabx2s() - Pixel-to-world transformation
* ----------------------------------------
* tabx2s() transforms intermediate world coordinates to world coordinates
* using coordinate lookup.
*
* Given and returned:
*   tab       struct tabprm*
*                       Tabular transformation parameters.
*
* Given:
*   ncoord,
*   nelem     int       The number of coordinates, each of vector length
*                       nelem.
*
*   x         const double[ncoord][nelem]
*                       Array of intermediate world coordinates, SI units.
*
* Returned:
*   world     double[ncoord][nelem]
*                       Array of world coordinates, in SI units.
*
*   stat      int[ncoord]
*                       Status return value status for each coordinate:
*                         0: Success.
*                         1: Invalid intermediate world coordinate.
*
* Function return value:
*             int       Status return value:
*                         0: Success.
*                         1: Null tabprm pointer passed.
*                         3: Invalid tabular parameters.
*                         4: One or more of the x coordinates were invalid,
*                            as indicated by the stat vector.
*
*                       For returns > 1, a detailed error message is set in
*                       tabprm::err if enabled, see wcserr_enable().
*
*
* tabs2x() - World-to-pixel transformation
* ----------------------------------------
* tabs2x() transforms world coordinates to intermediate world coordinates.
*
* Given and returned:
*   tab       struct tabprm*
*                       Tabular transformation parameters.
*
* Given:
*   ncoord,
*   nelem     int       The number of coordinates, each of vector length
*                       nelem.
*   world     const double[ncoord][nelem]
*                       Array of world coordinates, in SI units.
*
* Returned:
*   x         double[ncoord][nelem]
*                       Array of intermediate world coordinates, SI units.
*   stat      int[ncoord]
*                       Status return value status for each vector element:
*                         0: Success.
*                         1: Invalid world coordinate.
*
* Function return value:
*             int       Status return value:
*                         0: Success.
*                         1: Null tabprm pointer passed.
*                         3: Invalid tabular parameters.
*                         5: One or more of the world coordinates were
*                            invalid, as indicated by the stat vector.
*
*                       For returns > 1, a detailed error message is set in
*                       tabprm::err if enabled, see wcserr_enable().
*
*
* tabprm struct - Tabular transformation parameters
* -------------------------------------------------
* The tabprm struct contains information required to transform tabular
* coordinates.  It consists of certain members that must be set by the user
* ("given") and others that are set by the WCSLIB routines ("returned").  Some
* of the latter are supplied for informational purposes while others are for
* internal use only.
*
*   int flag
*     (Given and returned) This flag must be set to zero whenever any of the
*     following tabprm structure members are set or changed:
*
*       - tabprm::M (q.v., not normally set by the user),
*       - tabprm::K (q.v., not normally set by the user),
*       - tabprm::map,
*       - tabprm::crval,
*       - tabprm::index,
*       - tabprm::coord.
*
*     This signals the initialization routine, tabset(), to recompute the
*     returned members of the tabprm struct.  tabset() will reset flag to
*     indicate that this has been done.
*
*     PLEASE NOTE: flag should be set to -1 when tabini() is called for the
*     first time for a particular tabprm struct in order to initialize memory
*     management.  It must ONLY be used on the first initialization otherwise
*     memory leaks may result.
*
*   int M
*     (Given or returned) Number of tabular coordinate axes.
*
*     If tabini() is used to initialize the tabprm struct (as would normally
*     be the case) then it will set M from the value passed to it as a
*     function argument.  The user should not subsequently modify it.
*
*   int *K
*     (Given or returned) Pointer to the first element of a vector of length
*     tabprm::M whose elements (K_1, K_2,... K_M) record the lengths of the
*     axes of the coordinate array and of each indexing vector.
*
*     If tabini() is used to initialize the tabprm struct (as would normally
*     be the case) then it will set K from the array passed to it as a
*     function argument.  The user should not subsequently modify it.
*
*   int *map
*     (Given) Pointer to the first element of a vector of length tabprm::M
*     that defines the association between axis m in the M-dimensional
*     coordinate array (1 <= m <= M) and the indices of the intermediate world
*     coordinate and world coordinate arrays, x[] and world[], in the argument
*     lists for tabx2s() and tabs2x().
*
*     When x[] and world[] contain the full complement of coordinate elements
*     in image-order, as will usually be the case, then map[m-1] == i-1 for
*     axis i in the N-dimensional image (1 <= i <= N).  In terms of the FITS
*     keywords
*
*       map[PVi_3a - 1] == i - 1.
*
*     However, a different association may result if x[], for example, only
*     contains a (relevant) subset of intermediate world coordinate elements.
*     For example, if M == 1 for an image with N > 1, it is possible to fill
*     x[] with the relevant coordinate element with nelem set to 1.  In this
*     case map[0] = 0 regardless of the value of i.
*
*   double *crval
*     (Given) Pointer to the first element of a vector of length tabprm::M
*     whose elements contain the index value for the reference pixel for each
*     of the tabular coordinate axes.
*
*   double **index
*     (Given) Pointer to the first element of a vector of length tabprm::M of
*     pointers to vectors of lengths (K_1, K_2,... K_M) of 0-relative indexes
*     (see tabprm::K).
*
*     The address of any or all of these index vectors may be set to zero,
*     i.e.
*
=       index[m] == 0;
*
*     this is interpreted as default indexing, i.e.
*
=       index[m][k] = k;
*
*   double *coord
*     (Given) Pointer to the first element of the tabular coordinate array,
*     treated as though it were defined as
*
=       double coord[K_M]...[K_2][K_1][M];
*
*     (see tabprm::K) i.e. with the M dimension varying fastest so that the
*     M elements of a coordinate vector are stored contiguously in memory.
*
*   int nc
*     (Returned) Total number of coordinate vectors in the coordinate array
*     being the product K_1 * K_2 * ... * K_M (see tabprm::K).
*
*   int padding
*     (An unused variable inserted for alignment purposes only.)
*
*   int *sense
*     (Returned) Pointer to the first element of a vector of length tabprm::M
*     whose elements indicate whether the corresponding indexing vector is
*     monotonic increasing (+1), or decreasing (-1).
*
*   int *p0
*     (Returned) Pointer to the first element of a vector of length tabprm::M
*     of interpolated indices into the coordinate array such that Upsilon_m,
*     as defined in Paper III, is equal to (p0[m] + 1) + tabprm::delta[m].
*
*   double *delta
*     (Returned) Pointer to the first element of a vector of length tabprm::M
*     of interpolated indices into the coordinate array such that Upsilon_m,
*     as defined in Paper III, is equal to (tabprm::p0[m] + 1) + delta[m].
*
*   double *extrema
*     (Returned) Pointer to the first element of an array that records the
*     minimum and maximum value of each element of the coordinate vector in
*     each row of the coordinate array, treated as though it were defined as
*
=       double extrema[K_M]...[K_2][2][M]
*
*     (see tabprm::K).  The minimum is recorded in the first element of the
*     compressed K_1 dimension, then the maximum.  This array is used by the
*     inverse table lookup function, tabs2x(), to speed up table searches.
*
*   struct wcserr *err
*     (Returned) If enabled, when an error status is returned this struct
*     contains detailed information about the error, see wcserr_enable().
*
*   int m_flag
*     (For internal use only.)
*   int m_M
*     (For internal use only.)
*   int m_N
*     (For internal use only.)
*   int set_M
*     (For internal use only.)
*   int m_K
*     (For internal use only.)
*   int m_map
*     (For internal use only.)
*   int m_crval
*     (For internal use only.)
*   int m_index
*     (For internal use only.)
*   int m_indxs
*     (For internal use only.)
*   int m_coord
*     (For internal use only.)
*
*
* Global variable: const char *tab_errmsg[] - Status return messages
* ------------------------------------------------------------------
* Error messages to match the status value returned from each function.
*
*===========================================================================*/

#ifndef WCSLIB_TAB
#define WCSLIB_TAB

#include "wcserr.h"

#ifdef __cplusplus
extern "C" {
#endif


extern const char *tab_errmsg[];

enum tab_errmsg_enum {
  TABERR_SUCCESS      = 0,	/* Success. */
  TABERR_NULL_POINTER = 1,	/* Null tabprm pointer passed. */
  TABERR_MEMORY       = 2,	/* Memory allocation failed. */
  TABERR_BAD_PARAMS   = 3,	/* Invalid tabular parameters. */
  TABERR_BAD_X        = 4,	/* One or more of the x coordinates were
				   invalid. */
  TABERR_BAD_WORLD    = 5	/* One or more of the world coordinates were
				   invalid. */
};

struct tabprm {
  /* Initialization flag (see the prologue above).                          */
  /*------------------------------------------------------------------------*/
  int    flag;			/* Set to zero to force initialization.     */

  /* Parameters to be provided (see the prologue above).                    */
  /*------------------------------------------------------------------------*/
  int    M;			/* Number of tabular coordinate axes.       */
  int    *K;			/* Vector of length M whose elements        */
				/* (K_1, K_2,... K_M) record the lengths of */
				/* the axes of the coordinate array and of  */
				/* each indexing vector.                    */
  int    *map;			/* Vector of length M usually such that     */
				/* map[m-1] == i-1 for coordinate array     */
				/* axis m and image axis i (see above).     */
  double *crval;		/* Vector of length M containing the index  */
				/* value for the reference pixel for each   */
				/* of the tabular coordinate axes.          */
  double **index;		/* Vector of pointers to M indexing vectors */
				/* of lengths (K_1, K_2,... K_M).           */
  double *coord;		/* (1+M)-dimensional tabular coordinate     */
				/* array (see above).                       */

  /* Information derived from the parameters supplied.                      */
  /*------------------------------------------------------------------------*/
  int    nc;			/* Number of coordinate vectors (of length  */
				/* M) in the coordinate array.              */
  int    padding;		/* (Dummy inserted for alignment purposes.) */
  int    *sense;		/* Vector of M flags that indicate whether  */
				/* the Mth indexing vector is monotonic     */
				/* increasing, or else decreasing.          */
  int    *p0;			/* Vector of M indices.                     */
  double *delta;		/* Vector of M increments.                  */
  double *extrema;		/* (1+M)-dimensional array of coordinate    */
				/* extrema.                                 */

  /* Error handling                                                         */
  /*------------------------------------------------------------------------*/
  struct wcserr *err;

  /* Private - the remainder are for memory management.                     */
  /*------------------------------------------------------------------------*/
  int    m_flag, m_M, m_N;
  int    set_M;
  int    *m_K, *m_map;
  double *m_crval, **m_index, **m_indxs, *m_coord;
};

/* Size of the tabprm struct in int units, used by the Fortran wrappers. */
#define TABLEN (sizeof(struct tabprm)/sizeof(int))


int tabini(int alloc, int M, const int K[], struct tabprm *tab);

int tabmem(struct tabprm *tab);

int tabcpy(int alloc, const struct tabprm *tabsrc, struct tabprm *tabdst);

int tabcmp(int cmp, double tol, const struct tabprm *tab1,
           const struct tabprm *tab2, int *equal);

int tabfree(struct tabprm *tab);

int tabprt(const struct tabprm *tab);

int tabset(struct tabprm *tab);

int tabx2s(struct tabprm *tab, int ncoord, int nelem, const double x[],
           double world[], int stat[]);

int tabs2x(struct tabprm *tab, int ncoord, int nelem, const double world[],
           double x[], int stat[]);


/* Deprecated. */
#define tabini_errmsg tab_errmsg
#define tabcpy_errmsg tab_errmsg
#define tabfree_errmsg tab_errmsg
#define tabprt_errmsg tab_errmsg
#define tabset_errmsg tab_errmsg
#define tabx2s_errmsg tab_errmsg
#define tabs2x_errmsg tab_errmsg

#ifdef __cplusplus
}
#endif

#endif /* WCSLIB_TAB */
