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
  $Id: lin.h,v 7.11 2022/04/26 06:13:52 mcalabre Exp $
*=============================================================================
*
* WCSLIB 7.11 - C routines that implement the FITS World Coordinate System
* (WCS) standard.  Refer to the README file provided with WCSLIB for an
* overview of the library.
*
*
* Summary of the lin routines
* ---------------------------
* Routines in this suite apply the linear transformation defined by the FITS
* World Coordinate System (WCS) standard, as described in
*
=   "Representations of world coordinates in FITS",
=   Greisen, E.W., & Calabretta, M.R. 2002, A&A, 395, 1061 (WCS Paper I)
*
* These routines are based on the linprm struct which contains all information
* needed for the computations.  The struct contains some members that must be
* set by the user, and others that are maintained by these routines, somewhat
* like a C++ class but with no encapsulation.
*
* Six routines, linini(), lininit(), lindis(), lindist() lincpy(), and
* linfree() are provided to manage the linprm struct, linsize() computes its
* total size including allocated memory, and linprt() prints its contents.
*
* linperr() prints the error message(s) (if any) stored in a linprm struct,
* and the disprm structs that it may contain.
*
* A setup routine, linset(), computes intermediate values in the linprm struct
* from parameters in it that were supplied by the user.  The struct always
* needs to be set up by linset() but need not be called explicitly - refer to
* the explanation of linprm::flag.
*
* linp2x() and linx2p() implement the WCS linear transformations.
*
* An auxiliary routine, linwarp(), computes various measures of the distortion
* over a specified range of pixel coordinates.
*
* An auxiliary matrix inversion routine, matinv(), is included.  It uses
* LU-triangular factorization with scaled partial pivoting.
*
*
* linini() - Default constructor for the linprm struct
* ----------------------------------------------------
* linini() is a thin wrapper on lininit().  It invokes it with ndpmax set
* to -1 which causes it to use the value of the global variable NDPMAX.  It
* is thereby potentially thread-unsafe if NDPMAX is altered dynamically via
* disndp().  Use lininit() for a thread-safe alternative in this case.
*
*
* lininit() - Default constructor for the linprm struct
* -----------------------------------------------------
* lininit() allocates memory for arrays in a linprm struct and sets all
* members of the struct to default values.
*
* PLEASE NOTE: every linprm struct must be initialized by lininit(), possibly
* repeatedly.  On the first invokation, and only the first invokation,
* linprm::flag must be set to -1 to initialize memory management, regardless
* of whether lininit() will actually be used to allocate memory.
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
* Given:
*   ndpmax    int       The number of DPja or DQia keywords to allocate space
*                       for.  If set to -1, the value of the global variable
*                       NDPMAX will be used.  This is potentially
*                       thread-unsafe if disndp() is being used dynamically to
*                       alter its value.
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
* lindis() - Assign a distortion to a linprm struct
* -------------------------------------------------
* lindis() is a thin wrapper on lindist().   It invokes it with ndpmax set
* to -1 which causes the value of the global variable NDPMAX to be used (by
* disinit()).  It is thereby potentially thread-unsafe if NDPMAX is altered
* dynamically via disndp().  Use lindist() for a thread-safe alternative in
* this case.
*
*
* lindist() - Assign a distortion to a linprm struct
* --------------------------------------------------
* lindist() may be used to assign the address of a disprm struct to
* linprm::dispre or linprm::disseq.  The linprm struct must already have been
* initialized by lininit().
*
* The disprm struct must have been allocated from the heap (e.g. using
* malloc(), calloc(), etc.).  lindist() will immediately initialize it via a
* call to disini() using the value of linprm::naxis.  Subsequently, it will be
* reinitialized by calls to lininit(), and freed by linfree(), neither of
* which would happen if the disprm struct was assigned directly.
*
* If the disprm struct had previously been assigned via lindist(), it will be
* freed before reassignment.  It is also permissable for a null disprm pointer
* to be assigned to disable the distortion correction.
*
* Given:
*   sequence  int       Is it a prior or sequent distortion?
*                         1: Prior,   the assignment is to linprm::dispre.
*                         2: Sequent, the assignment is to linprm::disseq.
*
*                       Anything else is an error.
*
* Given and returned:
*   lin       struct linprm*
*                       Linear transformation parameters.
*
*   dis       struct disprm*
*                       Distortion function parameters.
*
* Given:
*   ndpmax    int       The number of DPja or DQia keywords to allocate space
*                       for.  If set to -1, the value of the global variable
*                       NDPMAX will be used.  This is potentially
*                       thread-unsafe if disndp() is being used dynamically to
*                       alter its value.
*
* Function return value:
*             int       Status return value:
*                         0: Success.
*                         1: Null linprm pointer passed.
*                         4: Invalid sequence.
*
*
* lincpy() - Copy routine for the linprm struct
* ---------------------------------------------
* lincpy() does a deep copy of one linprm struct to another, using lininit()
* to allocate memory for its arrays if required.  Only the "information to be
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
* linfree() frees memory allocated for the linprm arrays by lininit() and/or
* linset().  lininit() keeps a record of the memory it allocates and linfree()
* will only attempt to free this.
*
* PLEASE NOTE: linfree() must not be invoked on a linprm struct that was not
* initialized by lininit().
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
* linsize() - Compute the size of a linprm struct
* -----------------------------------------------
* linsize() computes the full size of a linprm struct, including allocated
* memory.
*
* Given:
*   lin       const struct linprm*
*                       Linear transformation parameters.
*
*                       If NULL, the base size of the struct and the allocated
*                       size are both set to zero.
*
* Returned:
*   sizes     int[2]    The first element is the base size of the struct as
*                       returned by sizeof(struct linprm).
*
*                       The second element is the total size of memory
*                       allocated in the struct, in bytes, assuming that the
*                       allocation was done by linini().  This figure includes
*                       memory allocated for members of constituent structs,
*                       such as linprm::dispre.
*
*                       It is not an error for the struct not to have been set
*                       up via linset(), which normally results in additional
*                       memory allocation.
*
* Function return value:
*             int       Status return value:
*                         0: Success.
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
* linperr() - Print error messages from a linprm struct
* -----------------------------------------------------
* linperr() prints the error message(s) (if any) stored in a linprm struct,
* and the disprm structs that it may contain.  If there are no errors then
* nothing is printed.  It uses wcserr_prt(), q.v.
*
* Given:
*   lin       const struct linprm*
*                       Coordinate transformation parameters.
*
*   prefix    const char *
*                       If non-NULL, each output line will be prefixed with
*                       this string.
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
* linwarp() - Compute measures of distortion
* ------------------------------------------
* linwarp() computes various measures of the distortion over a specified range
* of pixel coordinates.
*
* All distortion measures are specified as an offset in pixel coordinates,
* as given directly by prior distortions.  The offset in intermediate pixel
* coordinates given by sequent distortions is translated back to pixel
* coordinates by applying the inverse of the linear transformation matrix
* (PCi_ja or CDi_ja).  The difference may be significant if the matrix
* introduced a scaling.
*
* If all distortions are prior, then linwarp() uses diswarp(), q.v.
*
* Given and returned:
*   lin       struct linprm*
*                       Linear transformation parameters plus distortions.
*
* Given:
*   pixblc    const double[naxis]
*                       Start of the range of pixel coordinates (i.e. "bottom
*                       left-hand corner" in the conventional FITS image
*                       display orientation).  May be specified as a NULL
*                       pointer which is interpreted as (1,1,...).
*
*   pixtrc    const double[naxis]
*                       End of the range of pixel coordinates (i.e. "top
*                       right-hand corner" in the conventional FITS image
*                       display orientation).
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
*                         1: Null linprm pointer passed.
*                         2: Memory allocation failed.
*                         3: Invalid parameter.
*                         4: Distort error.
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
*       - linprm::cdelt,
*       - linprm::dispre.
*       - linprm::disseq.
*
*     This signals the initialization routine, linset(), to recompute the
*     returned members of the linprm struct.  linset() will reset flag to
*     indicate that this has been done.
*
*     PLEASE NOTE: flag should be set to -1 when lininit() is called for the
*     first time for a particular linprm struct in order to initialize memory
*     management.  It must ONLY be used on the first initialization otherwise
*     memory leaks may result.
*
*   int naxis
*     (Given or returned) Number of pixel and world coordinate elements.
*
*     If lininit() is used to initialize the linprm struct (as would normally
*     be the case) then it will set naxis from the value passed to it as a
*     function argument.  The user should not subsequently modify it.
*
*   double *crpix
*     (Given) Pointer to the first element of an array of double containing
*     the coordinate reference pixel, CRPIXja.
*
*     It is not necessary to reset the linprm struct (via linset()) when
*     linprm::crpix is changed.
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
*   struct disprm *dispre
*     (Given) Pointer to a disprm struct holding parameters for prior
*     distortion functions, or a null (0x0) pointer if there are none.
*
*     Function lindist() may be used to assign a disprm pointer to a linprm
*     struct, allowing it to take control of any memory allocated for it, as
*     in the following example:
*
=       void add_distortion(struct linprm *lin)
=       {
=         struct disprm *dispre;
=
=         dispre = malloc(sizeof(struct disprm));
=         dispre->flag = -1;
=         lindist(1, lin, dispre, ndpmax);
=           :
=          (Set up dispre.)
=           :
=
=         return;
=       }
*
*     Here, after the distortion function parameters etc. are copied into
*     dispre, dispre is assigned using lindist() which takes control of the
*     allocated memory.  It will be freed later when linfree() is invoked on
*     the linprm struct.
*
*     Consider also the following erroneous code:
*
=       void bad_code(struct linprm *lin)
=       {
=         struct disprm dispre;
=
=         dispre.flag = -1;
=         lindist(1, lin, &dispre, ndpmax);   // WRONG.
=           :
=
=         return;
=       }
*
*     Here, dispre is declared as a struct, rather than a pointer.  When the
*     function returns, dispre will go out of scope and its memory will most
*     likely be reused, thereby trashing its contents.  Later, a segfault will
*     occur when linfree() tries to free dispre's stale address.
*
*   struct disprm *disseq
*     (Given) Pointer to a disprm struct holding parameters for sequent
*     distortion functions, or a null (0x0) pointer if there are none.
*
*     Refer to the comments and examples given for disprm::dispre.
*
*   double *piximg
*     (Returned) Pointer to the first element of the matrix containing the
*     product of the CDELTia diagonal matrix and the PCi_ja matrix.
*
*   double *imgpix
*     (Returned) Pointer to the first element of the inverse of the
*     linprm::piximg matrix.
*
*   int i_naxis
*     (Returned) The dimension of linprm::piximg and linprm::imgpix (normally
*     equal to naxis).
*
*   int unity
*     (Returned) True if the linear transformation matrix is unity.
*
*   int affine
*     (Returned) True if there are no distortions.
*
*   int simple
*     (Returned) True if unity and no distortions.
*
*   struct wcserr *err
*     (Returned) If enabled, when an error status is returned, this struct
*     contains detailed information about the error, see wcserr_enable().
*
*   double *tmpcrd
*     (For internal use only.)
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
*   struct disprm *m_dispre
*     (For internal use only.)
*   struct disprm *m_disseq
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

#ifdef __cplusplus
extern "C" {
#endif


extern const char *lin_errmsg[];

enum lin_errmsg_enum {
  LINERR_SUCCESS      = 0,	// Success.
  LINERR_NULL_POINTER = 1,	// Null linprm pointer passed.
  LINERR_MEMORY       = 2,	// Memory allocation failed.
  LINERR_SINGULAR_MTX = 3,	// PCi_ja matrix is singular.
  LINERR_DISTORT_INIT = 4,	// Failed to initialise distortions.
  LINERR_DISTORT      = 5,	// Distort error.
  LINERR_DEDISTORT    = 6	// De-distort error.
};

struct linprm {
  // Initialization flag (see the prologue above).
  //--------------------------------------------------------------------------
  int flag;			// Set to zero to force initialization.

  // Parameters to be provided (see the prologue above).
  //--------------------------------------------------------------------------
  int naxis;			// The number of axes, given by NAXIS.
  double *crpix;		// CRPIXja keywords for each pixel axis.
  double *pc;			// PCi_ja  linear transformation matrix.
  double *cdelt;		// CDELTia keywords for each coord axis.
  struct disprm *dispre;	// Prior   distortion parameters, if any.
  struct disprm *disseq;	// Sequent distortion parameters, if any.

  // Information derived from the parameters supplied.
  //--------------------------------------------------------------------------
  double *piximg;		// Product of CDELTia and PCi_ja matrices.
  double *imgpix;		// Inverse of the piximg matrix.
  int    i_naxis;		// Dimension of piximg and imgpix.
  int    unity;			// True if the PCi_ja matrix is unity.
  int    affine;		// True if there are no distortions.
  int    simple;		// True if unity and no distortions.

  // Error handling, if enabled.
  //--------------------------------------------------------------------------
  struct wcserr *err;

  // Private - the remainder are for internal use.
  //--------------------------------------------------------------------------
  double *tmpcrd;

  int    m_flag, m_naxis;
  double *m_crpix, *m_pc, *m_cdelt;
  struct disprm *m_dispre, *m_disseq;
};

// Size of the linprm struct in int units, used by the Fortran wrappers.
#define LINLEN (sizeof(struct linprm)/sizeof(int))


int linini(int alloc, int naxis, struct linprm *lin);

int lininit(int alloc, int naxis, struct linprm *lin, int ndpmax);

int lindis(int sequence, struct linprm *lin, struct disprm *dis);

int lindist(int sequence, struct linprm *lin, struct disprm *dis, int ndpmax);

int lincpy(int alloc, const struct linprm *linsrc, struct linprm *lindst);

int linfree(struct linprm *lin);

int linsize(const struct linprm *lin, int sizes[2]);

int linprt(const struct linprm *lin);

int linperr(const struct linprm *lin, const char *prefix);

int linset(struct linprm *lin);

int linp2x(struct linprm *lin, int ncoord, int nelem, const double pixcrd[],
           double imgcrd[]);

int linx2p(struct linprm *lin, int ncoord, int nelem, const double imgcrd[],
           double pixcrd[]);

int linwarp(struct linprm *lin, const double pixblc[], const double pixtrc[],
            const double pixsamp[], int *nsamp,
            double maxdis[], double *maxtot,
            double avgdis[], double *avgtot,
            double rmsdis[], double *rmstot);

int matinv(int n, const double mat[], double inv[]);


// Deprecated.
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

#endif // WCSLIB_LIN
