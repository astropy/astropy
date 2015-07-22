/*============================================================================

  WCSLIB 5.8 - an implementation of the FITS WCS standard.
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
  $Id: dis.h,v 5.8 2015/07/08 11:03:59 mcalabre Exp $
*=============================================================================
*
* WCSLIB 5.8 - C routines that implement the FITS World Coordinate System
* (WCS) standard.  Refer to the README file provided with WCSLIB for an
* overview of the library.
*
*
* Summary of the dis routines
* ---------------------------
* Routines in this suite implement extensions to the FITS World Coordinate
* System (WCS) standard proposed by
*
=   "Representations of distortions in FITS world coordinate systems",
=   Calabretta, M.R. et al. (WCS Paper IV, draft dated 2004/04/22),
=   available from http://www.atnf.csiro.au/people/Mark.Calabretta
*
* In brief, a distortion function may occupy one of two positions in the WCS
* algorithm chain.  Prior distortions precede the linear transformation
* matrix, whether it be PCi_ja or CDi_ja, and sequent distortions follow it.
* WCS Paper IV defines FITS keywords used to specify parameters for predefined
* distortion functions.  The following are used for prior distortions:
*
=   CPDISja   ...(string-valued, identifies the distortion function)
=   DPja      ...(record-valued, parameters)
=   CPERRja   ...(floating-valued, maximum value)
*
* Their counterparts for sequent distortions are CQDISia, DQia, and CQERRia.
* An additional floating-valued keyword, DVERRa, records the maximum value of
* the combined distortions.
*
* DPja and DQia are "record-valued".  Syntactically, the keyvalues are
* standard FITS strings, but they are to be interpreted in a special way.
* The general form is
*
=   DPja = '<field-specifier>: <float>'
*
* where the field-specifier consists of a sequence of fields separated by
* periods, and the ': ' between the field-specifier and the floating-point
* value is part of the record syntax.  For example:
*
=   DP1 = 'AXIS.1: 1'
*
* Certain field-specifiers are defined for all distortion functions, while
* others are defined only for particular distortions.  Refer to WCS Paper IV
* for further details.  wcspih() parses all distortion keywords and loads them
* into a disprm struct for analysis by disset() which knows (or possibly does
* not know) how to interpret them.  Of the Paper IV distortion functions, only
* the general Polynomial distortion is currently implemented here.
*
* TPV: The TPV "projection"
* -------------------------
* The distortion function component of the TPV celestial "projection" is also
* supported.  The TPV projection, originally proposed in a draft of WCS Paper
* II, consists of a TAN projection with sequent polynomial distortion, the
* coefficients of which are encoded in PVi_ma keyrecords.  Full details may be
* found at the registry of FITS conventions:
*
=   http://fits.gsfc.nasa.gov/registry/tpvwcs/tpv.html
*
* Internally, wcsset() changes TPV to a TAN projection, translates the PVi_ma
* keywords to DQia and loads them into a disprm struct.  These DQia keyrecords
* have the form
*
=   DQia = 'TPV.m: <value>'
*
* where i, a, m, and the value for each DQia match each PVi_ma.  Consequently,
* WCSLIB would handle a FITS header containing these keywords, along with
* CQDISia = 'TPV' and the required DQia.NAXES and DQia.AXIS.ihat keywords.
*
* SIP: Simple Imaging Polynomial
* ------------------------------
* These routines also support the Simple Imaging Polynomial (SIP), whose
* design was influenced by early drafts of WCS Paper IV.  It is described in
* detail in
*
=   http://fits.gsfc.nasa.gov/registry/sip.html
*
* SIP, which is defined only as a prior distortion for 2D celestial images,
* has the interesting feature that it records an approximation to the inverse
* polynomial distortion function.  This is used by disx2p() to provide an
* initial estimate for its more precise iterative inversion.  The
* special-purpose keywords used by SIP are parsed and translated by wcspih()
* as follows:
*
=    A_p_q = <value>   ->   DP1 = 'SIP.FWD.p_q: <value>'
=   AP_p_q = <value>   ->   DP1 = 'SIP.REV.p_q: <value>'
=    B_p_q = <value>   ->   DP2 = 'SIP.FWD.p_q: <value>'
=   BP_p_q = <value>   ->   DP2 = 'SIP.REV.p_q: <value>'
=   A_DMAX = <value>   ->   DPERR1 = <value>
=   B_DMAX = <value>   ->   DPERR2 = <value>
*
* SIP's A_ORDER and B_ORDER keywords are not used.  WCSLIB would recognise a
* FITS header containing the above keywords, along with CPDISja = 'SIP' and
* the required DPja.NAXES keywords.
*
* DSS: Digitized Sky Survey
* -------------------------
* The Digitized Sky Survey resulted from the production of the Guide Star
* Catalogue for the Hubble Space Telescope.  Plate solutions based on a
* polynomial distortion function were encoded in FITS using non-standard
* keywords.  Sect. 5.2 of WCS Paper IV describes how DSS coordinates may be
* translated to a sequent Polynomial distortion using two auxiliary variables.
* That translation is based on optimising the non-distortion component of the
* plate solution.
*
* Following Paper IV, wcspih() translates the non-distortion component of DSS
* coordinates to standard WCS keywords (CRPIXja, PCi_ja, CRVALia, etc), and
* fills a wcsprm struct with their values.  It encodes the DSS polynomial
* coefficients as
*
=    AMDXm = <value>   ->   DQ1 = 'AMD.m: <value>'
=    AMDYm = <value>   ->   DQ2 = 'AMD.m: <value>'
*
* WCSLIB would recognise a FITS header containing the above keywords, along
* with CQDISia = 'DSS' and the required DQia.NAXES keywords.
*
* WAT: The TNX and ZPX "projections"
* ----------------------------------
* The TNX and ZPX "projections" add a polynomial distortion function to the
* standard TAN and ZPN projections respectively.  Unusually, the polynomial
* may be expressed as the sum of Chebyshev or Legendre polynomials, or as a
* simple sum of monomials, as described in
*
=   http://fits.gsfc.nasa.gov/registry/tnx/tnx-doc.html
=   http://fits.gsfc.nasa.gov/registry/zpxwcs/zpx.html
*
* The polynomial coefficients are encoded in special-purpose WATi_n keywords
* as a set of continued strings, thus providing the name for this distortion
* type.  WATi_n are parsed and translated by wcspih() into the following set:
*
=    DQi = 'WAT.POLY: <value>'
=    DQi = 'WAT.XMIN: <value>'
=    DQi = 'WAT.XMAX: <value>'
=    DQi = 'WAT.YMIN: <value>'
=    DQi = 'WAT.YMAX: <value>'
=    DQi = 'WAT.CHBY.m_n: <value>'  or
=    DQi = 'WAT.LEGR.m_n: <value>'  or
=    DQi = 'WAT.MONO.m_n: <value>'
*
* along with CQDISia = 'WAT' and the required DPja.NAXES keywords.  For ZPX,
* the ZPN projection parameters are also encoded in WATi_n, and wcspih()
* translates these to standard PVi_ma.
*
* TPD: Template Polynomial Distortion
* -----------------------------------
* The "Template Polynomial Distortion" (TPD) is a superset of the TPV, SIP,
* DSS, and WAT (TNX & ZPX) polynomial distortions that also supports 1D usage
* and inversions.  Like TPV, SIP, and DSS, the form of the polynomial is fixed
* (the "template") and only the coefficients for the required terms are set
* non-zero.  TPD generalizes TPV in going to 9th degree, SIP by accomodating
* TPV's linear and radial terms, and DSS in both respects.  While in theory
* the degree of the WAT polynomial distortion in unconstrained, in practice it
* is limited to values that can be handled by TPD.
*
* Within WCSLIB, TPV, SIP, DSS, and WAT are all implemented as special cases
* of TPD.  Indeed, TPD was developed precisely for that purpose.  WAT
* distortions expressed as the sum of Chebyshev or Legendre polynomials are
* expanded for TPD as a simple sum of monomials.  Moreover, the general
* Polynomial distortion is translated and implemented internally as TPD
* whenever possible.
*
* However, WCSLIB also recognizes 'TPD' as a distortion function in its own
* right (i.e. a recognized value of CPDISja or CQDISia), for use as both prior
* and sequent distortions.  Its DPja and DQia keyrecords have the form
*
=   DPja = 'TPD.FWD.m: <value>'
=   DPja = 'TPD.REV.m: <value>'
*
* for the forward and reverse distortion functions.  Moreover, like the
* general Polynomial distortion, TPD supports auxiliary variables, though only
* as a linear transformation of pixel coordinates (p1,p2):
*
=   x = a0 + a1*p1 + a2*p2
=   y = b0 + b1*p1 + b2*p2
*
* where the coefficients of the auxiliary variables (x,y) are recorded as
*
=   DPja = 'AUX.1.COEFF.0: a0'      ...default 0.0
=   DPja = 'AUX.1.COEFF.1: a1'      ...default 1.0
=   DPja = 'AUX.1.COEFF.2: a2'      ...default 0.0
=   DPja = 'AUX.2.COEFF.0: b0'      ...default 0.0
=   DPja = 'AUX.2.COEFF.1: b1'      ...default 0.0
=   DPja = 'AUX.2.COEFF.2: b2'      ...default 1.0
*
* Though nowhere near as powerful, in typical applications TPD is considerably
* faster than the general Polynomial distortion.  As TPD has a finite and not
* too large number of possible terms (60), the coefficients for each can be
* stored (by disset()) in a fixed location in the disprm::dparm[] array.  A
* large part of the speedup then arises from evaluating the polynomial using
* Horner's scheme.
*
* Separate implementations for polynomials of each degree, and conditionals
* for 1D polynomials and 2D polynomials with and without the radial variable,
* ensure that unused terms mostly do not impose a significant computational
* overhead.
*
* The TPD terms are as follows
*
=   0: 1     4: xx      12: xxxx      24: xxxxxx      40: xxxxxxxx
=            5: xy      13: xxxy      25: xxxxxy      41: xxxxxxxy
=   1: x     6: yy      14: xxyy      26: xxxxyy      42: xxxxxxyy
=   2: y                15: xyyy      27: xxxyyy      43: xxxxxyyy
=   3: r     7: xxx     16: yyyy      28: xxyyyy      44: xxxxyyyy
=            8: xxy                   29: xyyyyy      45: xxxyyyyy
=            9: xyy     17: xxxxx     30: yyyyyy      46: xxyyyyyy
=           10: yyy     18: xxxxy                     47: xyyyyyyy
=           11: rrr     19: xxxyy     31: xxxxxxx     48: yyyyyyyy
=                       20: xxyyy     32: xxxxxxy
=                       21: xyyyy     33: xxxxxyy     49: xxxxxxxxx
=                       22: yyyyy     34: xxxxyyy     50: xxxxxxxxy
=                       23: rrrrr     35: xxxyyyy     51: xxxxxxxyy
=                                     36: xxyyyyy     52: xxxxxxyyy
=                                     37: xyyyyyy     53: xxxxxyyyy
=                                     38: yyyyyyy     54: xxxxyyyyy
=                                     39: rrrrrrr     55: xxxyyyyyy
=                                                     56: xxyyyyyyy
=                                                     57: xyyyyyyyy
=                                                     58: yyyyyyyyy
=                                                     59: rrrrrrrrr
*
* where r = sqrt(xx + yy).  Note that even powers of r are excluded since they
* can be accomodated by powers of (xx + yy).
*
* TPV uses all terms up to 39.  The m in its PVi_ma keywords translates
* directly to the TPD coefficient number.
*
* SIP uses all terms except for 0, 3, 11, 23, 39, and 59, with terms 1 and 2
* only used for the inverse.  Its A_p_q, etc. keywords must be translated
* using a map.
*
* DSS uses terms 0, 1, 2, 4, 5, 6, 7, 8, 9, 10, 17, 19, and 21.  The presence
* of a non-zero constant term arises through the use of auxiliary variables
* with origin offset from the reference point of the TAN projection.  However,
* in the translation given by WCS Paper IV, the distortion polynomial is zero,
* or very close to zero, at the reference pixel itself.  The mapping between
* DSS's AMDXm (or AMDYm) keyvalues and TPD coefficients, while still simple,
* is not quite as straightforward as for TPV and SIP.
*
* WAT uses all but the radial terms: 3, 11, 23, 39, and 59.  While the mapping
* between WAT's monomial coefficients and TPD is fairly simple, for its
* expression in terms of a sum of Chebyshev or Legendre polynomials it is much
* less so. 
*
* Summary of the dis routines:
* ----------------------------
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
* disperr() prints the error message(s) (if any) stored in a disprm struct.
*
* wcshdo() normally writes SIP and TPV headers in their native form if at all
* possible.  However, dishdo() may be used to set a flag that tells it to
* write the header in the form of the TPD translation used internally.
*
* A setup routine, disset(), computes intermediate values in the disprm struct
* from parameters in it that were supplied by the user.  The struct always
* needs to be set up by disset(), though disset() need not be called
* explicitly - refer to the explanation of disprm::flag.
*
* disp2x() and disx2p() implement the WCS distortion functions, disp2x() using
* separate functions, such as dispoly() and tpd7(), to do the computation.
*
* An auxiliary routine, diswarp(), computes various measures of the distortion
* over a specified range of coordinates.
*
* PLEASE NOTE: Distortions are not yet handled by wcsbth(), or wcscompare().
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
*                       If j is non-zero, and keyword was given, then the
*                       value of j will be used to fill in the axis number.
*
*   type      int       Data type of the record's value
*                         0: Integer,
*                         1: Floating point.
*
*   i         int       For type == 0, the integer value of the record.
*
*   f         double    For type == 1, the floating point value of the record.
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
* disperr() - Print error messages from a disprm struct
* -----------------------------------------------------
* disperr() prints the error message(s) (if any) stored in a disprm struct.
* If there are no errors then nothing is printed.  It uses wcserr_prt(), q.v.
*
* Given:
*   dis       const struct disprm*
*                       Distortion function parameters.
*
*   prefix    const char *
*                       If non-NULL, each output line will be prefixed with
*                       this string.
*
* Function return value:
*             int       Status return value:
*                         0: Success.
*                         1: Null disprm pointer passed.
*
*
* dishdo() - write FITS headers using TPD
* ---------------------------------------
* dishdo() sets a flag that tells wcshdo() to write FITS headers in the form
* of the TPD translation used internally.  Normally SIP and TPV would be
* written in their native form if at all possible.
*
* Given and returned:
*   dis       struct disprm*
*                       Distortion function parameters.
*
* Function return value:
*             int       Status return value:
*                         0: Success.
*                         1: Null disprm pointer passed.
*                         3: No TPD translation.
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
* If disprm::disx2p() is defined, then disx2p() uses it to provide an initial
* estimate for its more precise iterative inversion.
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
*     function with the 0-relative image axis number.  For example, consider
*     an image with a spectrum on the first axis (axis 0), followed by RA
*     (axis 1), Dec (axis2), and time (axis 3) axes.  For a distortion in
*     (RA,Dec) and no distortion on the spectral or time axes, the axis
*     mapping arrays, axmap[j][], would be
*
=       j=0: [-1, -1, -1, -1]   ...no  distortion on spectral axis,
=         1: [ 1,  2, -1, -1]   ...RA  distortion depends on RA and Dec,
=         2: [ 2,  1, -1, -1]   ...Dec distortion depends on Dec and RA,
=         3: [-1, -1, -1, -1]   ...no  distortion on time axis,
*
*     where -1 indicates that there is no corresponding independent
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
*       - dpkey::i,
*       - dpkey::f,
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
#define DISP2X_ARGS int inverse, const int iparm[], const double dparm[], \
int ncrd, const double rawcrd[], double *discrd

/* For use in declaring de-distortion function prototypes (= DISP2X_ARGS). */
#define DISX2P_ARGS int inverse, const int iparm[], const double dparm[], \
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
           int type, int i, double f);

int disini(int alloc, int naxis, struct disprm *dis);

int discpy(int alloc, const struct disprm *dissrc, struct disprm *disdst);

int disfree(struct disprm *dis);

int disprt(const struct disprm *dis);

int disperr(const struct disprm *dis, const char *prefix);

int dishdo(struct disprm *dis);

int disset(struct disprm *dis);

int disp2x(struct disprm *dis, const double rawcrd[], double discrd[]);

int disx2p(struct disprm *dis, const double discrd[], double rawcrd[]);

int diswarp(struct disprm *dis, const double pixblc[], const double pixtrc[],
            const double pixsamp[], int *nsamp,
            double maxdis[], double *maxtot,
            double avgdis[], double *avgtot,
            double rmsdis[], double *rmstot);

#ifdef __cplusplus
}
#endif

#endif /* WCSLIB_DIS */
