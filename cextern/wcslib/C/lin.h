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
  $Id: lin.h,v 4.25 2014/12/14 14:29:36 mcalabre Exp $
*=============================================================================
*
* WCSLIB 4.25 - C routines that implement the FITS World Coordinate System
* (WCS) standard.  Refer to
*
*   "Representations of world coordinates in FITS",
*   Greisen, E.W., & Calabretta, M.R. 2002, A&A, 395, 1061 (Paper I)
*
* Refer to the README file provided with WCSLIB for an overview of the
* library.
*
*
* Summary of the lin routines
* ---------------------------
* These routines apply the linear transformation defined by the FITS WCS
* standard.  They are based on the linprm struct which contains all
* information needed for the computations.  The struct contains some members
* that must be set by the user, and others that are maintained by these
* routines, somewhat like a C++ class but with no encapsulation.
*
* Three routines, linini(), lincpy(), and linfree() are provided to manage the
* linprm struct, and another, linprt(), prints its contents.
*
* A setup routine, linset(), computes intermediate values in the linprm struct
* from parameters in it that were supplied by the user.  The struct always
* needs to be set up by linset() but need not be called explicitly - refer to
* the explanation of linprm::flag.
*
* linp2x() and linx2p() implement the WCS linear transformations.
*
* An auxiliary matrix inversion routine, matinv(), is included.  It uses
* LU-triangular factorization with scaled partial pivoting.
*
*
* linini() - Default constructor for the linprm struct
* ----------------------------------------------------
* linini() allocates memory for arrays in a linprm struct and sets all members
* of the struct to default values.
*
* PLEASE NOTE: every linprm struct should be initialized by linini(), possibly
* repeatedly.  On the first invokation, and only the first invokation,
* linprm::flag must be set to -1 to initialize memory management, regardless
* of whether linini() will actually be used to allocate memory.
*
* Given:
*   alloc     int       If true, allocate memory unconditionally for arrays in
*                       the linprm struct.
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
*   lin       struct linprm*
*                       Linear transformation parameters.  Note that, in order
*                       to initialize memory management linprm::flag should be
*                       set to -1 when lin is initialized for the first time
*                       (memory leaks may result if it had already been
*                       initialized).
*
* Function return value:
*             int       Status return value:
*                         0: Success.
*                         1: Null linprm pointer passed.
*                         2: Memory allocation failed.
*
*                       For returns > 1, a detailed error message is set in
*                       linprm::err if enabled, see wcserr_enable().
*
*
* lincpy() - Copy routine for the linprm struct
* ---------------------------------------------
* lincpy() does a deep copy of one linprm struct to another, using linini() to
* allocate memory for its arrays if required.  Only the "information to be
* provided" part of the struct is copied; a call to linset() is required to
* initialize the remainder.
*
* Given:
*   alloc     int       If true, allocate memory for the crpix, pc, and cdelt
*                       arrays in the destination.  Otherwise, it is assumed
*                       that pointers to these arrays have been set by the
*                       user except if they are null pointers in which case
*                       memory will be allocated for them regardless.
*
*   linsrc    const struct linprm*
*                       Struct to copy from.
*
* Given and returned:
*   lindst    struct linprm*
*                       Struct to copy to.  linprm::flag should be set to -1
*                       if lindst was not previously initialized (memory leaks
*                       may result if it was previously initialized).
*
* Function return value:
*             int       Status return value:
*                         0: Success.
*                         1: Null linprm pointer passed.
*                         2: Memory allocation failed.
*
*                       For returns > 1, a detailed error message is set in
*                       linprm::err if enabled, see wcserr_enable().
*
*
* linfree() - Destructor for the linprm struct
* --------------------------------------------
* linfree() frees memory allocated for the linprm arrays by linini() and/or
* linset().  linini() keeps a record of the memory it allocates and linfree()
* will only attempt to free this.
*
* PLEASE NOTE: linfree() must not be invoked on a linprm struct that was not
* initialized by linini().
*
* Given:
*   lin       struct linprm*
*                       Linear transformation parameters.
*
* Function return value:
*             int       Status return value:
*                         0: Success.
*                         1: Null linprm pointer passed.
*
*
* linprt() - Print routine for the linprm struct
* ----------------------------------------------
* linprt() prints the contents of a linprm struct using wcsprintf().  Mainly
* intended for diagnostic purposes.
*
* Given:
*   lin       const struct linprm*
*                       Linear transformation parameters.
*
* Function return value:
*             int       Status return value:
*                         0: Success.
*                         1: Null linprm pointer passed.
*
*
* linset() - Setup routine for the linprm struct
* ----------------------------------------------
* linset(), if necessary, allocates memory for the linprm::piximg and
* linprm::imgpix arrays and sets up the linprm struct according to information
* supplied within it - refer to the explanation of linprm::flag.
*
* Note that this routine need not be called directly; it will be invoked by
* linp2x() and linx2p() if the linprm::flag is anything other than a
* predefined magic value.
*
* Given and returned:
*   lin       struct linprm*
*                       Linear transformation parameters.
*
* Function return value:
*             int       Status return value:
*                         0: Success.
*                         1: Null linprm pointer passed.
*                         2: Memory allocation failed.
*                         3: PCi_ja matrix is singular.
*
*                       For returns > 1, a detailed error message is set in
*                       linprm::err if enabled, see wcserr_enable().
*
*
* linp2x() - Pixel-to-world linear transformation
* -----------------------------------------------
* linp2x() transforms pixel coordinates to intermediate world coordinates.
*
* Given and returned:
*   lin       struct linprm*
*                       Linear transformation parameters.
*
* Given:
*   ncoord,
*   nelem     int       The number of coordinates, each of vector length nelem
*                       but containing lin.naxis coordinate elements.
*
*   pixcrd    const double[ncoord][nelem]
*                       Array of pixel coordinates.
*
* Returned:
*   imgcrd    double[ncoord][nelem]
*                       Array of intermediate world coordinates.
*
* Function return value:
*             int       Status return value:
*                         0: Success.
*                         1: Null linprm pointer passed.
*                         2: Memory allocation failed.
*                         3: PCi_ja matrix is singular.
*
*                       For returns > 1, a detailed error message is set in
*                       linprm::err if enabled, see wcserr_enable().
*
*
* linx2p() - World-to-pixel linear transformation
* -----------------------------------------------
* linx2p() transforms intermediate world coordinates to pixel coordinates.
*
* Given and returned:
*   lin       struct linprm*
*                       Linear transformation parameters.
*
* Given:
*   ncoord,
*   nelem     int       The number of coordinates, each of vector length nelem
*                       but containing lin.naxis coordinate elements.
*
*   imgcrd   const double[ncoord][nelem]
*                       Array of intermediate world coordinates.
*
* Returned:
*   pixcrd    double[ncoord][nelem]
*                       Array of pixel coordinates.
*
*             int       Status return value:
*                         0: Success.
*                         1: Null linprm pointer passed.
*                         2: Memory allocation failed.
*                         3: PCi_ja matrix is singular.
*
*                       For returns > 1, a detailed error message is set in
*                       linprm::err if enabled, see wcserr_enable().
*
*
* linprm struct - Linear transformation parameters
* ------------------------------------------------
* The linprm struct contains all of the information required to perform a
* linear transformation.  It consists of certain members that must be set by
* the user ("given") and others that are set by the WCSLIB routines
* ("returned").
*
*   int flag
*     (Given and returned) This flag must be set to zero whenever any of the
*     following members of the linprm struct are set or modified:
*
*       - linprm::naxis (q.v., not normally set by the user),
*       - linprm::pc,
*       - linprm::cdelt.
*
*     This signals the initialization routine, linset(), to recompute the
*     returned members of the linprm struct.  linset() will reset flag to
*     indicate that this has been done.
*
*     PLEASE NOTE: flag should be set to -1 when linini() is called for the
*     first time for a particular linprm struct in order to initialize memory
*     management.  It must ONLY be used on the first initialization otherwise
*     memory leaks may result.
*
*   int naxis
*     (Given or returned) Number of pixel and world coordinate elements.
*
*     If linini() is used to initialize the linprm struct (as would normally
*     be the case) then it will set naxis from the value passed to it as a
*     function argument.  The user should not subsequently modify it.
*
*   double *crpix
*     (Given) Pointer to the first element of an array of double containing
*     the coordinate reference pixel, CRPIXja.
*
*   double *pc
*     (Given) Pointer to the first element of the PCi_ja (pixel coordinate)
*     transformation matrix.  The expected order is
*
=       struct linprm lin;
=       lin.pc = {PC1_1, PC1_2, PC2_1, PC2_2};
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
=       lin.pc = *m;
*
*     would be legitimate.
*
*   double *cdelt
*     (Given) Pointer to the first element of an array of double containing
*     the coordinate increments, CDELTia.
*
*   int unity
*     (Returned) True if the linear transformation matrix is unity.
*
*   int padding
*     (An unused variable inserted for alignment purposes only.)
*
*   double *piximg
*     (Returned) Pointer to the first element of the matrix containing the
*     product of the CDELTia diagonal matrix and the PCi_ja matrix.
*
*   double *imgpix
*     (Returned) Pointer to the first element of the inverse of the
*     linprm::piximg matrix.
*
*   struct wcserr *err
*     (Returned) If enabled, when an error status is returned this struct
*     contains detailed information about the error, see wcserr_enable().
*
*   int i_naxis
*     (For internal use only.)
*   int m_flag
*     (For internal use only.)
*   int m_naxis
*     (For internal use only.)
*   int m_padding
*     (For internal use only.)
*   double *m_crpix
*     (For internal use only.)
*   double *m_pc
*     (For internal use only.)
*   double *m_cdelt
*     (For internal use only.)
*   void *padding2
*     (For internal use only.)
*
*
* Global variable: const char *lin_errmsg[] - Status return messages
* ------------------------------------------------------------------
* Error messages to match the status value returned from each function.
*
*===========================================================================*/

#ifndef WCSLIB_LIN
#define WCSLIB_LIN

#include "wcserr.h"

#ifdef __cplusplus
extern "C" {
#endif


extern const char *lin_errmsg[];

enum lin_errmsg_enum {
  LINERR_SUCCESS      = 0,	/* Success. */
  LINERR_NULL_POINTER = 1,	/* Null linprm pointer passed. */
  LINERR_MEMORY       = 2,	/* Memory allocation failed. */
  LINERR_SINGULAR_MTX = 3 	/* PCi_ja matrix is singular. */
};

struct linprm {
  /* Initialization flag (see the prologue above).                          */
  /*------------------------------------------------------------------------*/
  int flag;			/* Set to zero to force initialization.     */

  /* Parameters to be provided (see the prologue above).                    */
  /*------------------------------------------------------------------------*/
  int naxis;			/* The number of axes, given by NAXIS.      */
  double *crpix;		/* CRPIXja keywords for each pixel axis.    */
  double *pc;			/* PCi_ja  linear transformation matrix.    */
  double *cdelt;		/* CDELTia keywords for each coord axis.    */

  /* Information derived from the parameters supplied.                      */
  /*------------------------------------------------------------------------*/
  double *piximg;		/* Product of CDELTia and PCi_ja matrices.  */
  double *imgpix;		/* Inverse of the piximg matrix.            */
  int    unity;			/* True if the PCi_ja matrix is unity.      */

  /* Error handling                                                         */
  /*------------------------------------------------------------------------*/
  int    padding;		/* (Dummy inserted for alignment purposes.) */
  struct wcserr *err;

  /* Private - the remainder are for memory management.                     */
  /*------------------------------------------------------------------------*/
  int    i_naxis;
  int    m_flag, m_naxis, m_padding;
  double *m_crpix, *m_pc, *m_cdelt;
  void   *padding2;
};

/* Size of the linprm struct in int units, used by the Fortran wrappers. */
#define LINLEN (sizeof(struct linprm)/sizeof(int))


int linini(int alloc, int naxis, struct linprm *lin);

int lincpy(int alloc, const struct linprm *linsrc, struct linprm *lindst);

int linfree(struct linprm *lin);

int linprt(const struct linprm *lin);

int linset(struct linprm *lin);

int linp2x(struct linprm *lin, int ncoord, int nelem, const double pixcrd[],
           double imgcrd[]);

int linx2p(struct linprm *lin, int ncoord, int nelem, const double imgcrd[],
           double pixcrd[]);

int matinv(int n, const double mat[], double inv[]);


/* Deprecated. */
#define linini_errmsg lin_errmsg
#define lincpy_errmsg lin_errmsg
#define linfree_errmsg lin_errmsg
#define linprt_errmsg lin_errmsg
#define linset_errmsg lin_errmsg
#define linp2x_errmsg lin_errmsg
#define linx2p_errmsg lin_errmsg

#ifdef __cplusplus
}
#endif

#endif /* WCSLIB_LIN */
