/*============================================================================

  WCSLIB 5.2 - an implementation of the FITS WCS standard.
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
  $Id: spc.h,v 5.2 2015/04/15 12:35:07 mcalabre Exp $
*=============================================================================
*
* WCSLIB 5.2 - C routines that implement the spectral coordinate systems
* recognized by the FITS World Coordinate System (WCS) standard.  Refer to
*
*   "Representations of world coordinates in FITS",
*   Greisen, E.W., & Calabretta, M.R. 2002, A&A, 395, 1061 (Paper I)
*
*   "Representations of spectral coordinates in FITS",
*   Greisen, E.W., Calabretta, M.R., Valdes, F.G., & Allen, S.L.
*   2006, A&A, 446, 747 (Paper III)
*
* Refer to the README file provided with WCSLIB for an overview of the
* library.
*
*
* Summary of the spc routines
* ---------------------------
* These routines implement the part of the FITS WCS standard that deals with
* spectral coordinates.  They define methods to be used for computing spectral
* world coordinates from intermediate world coordinates (a linear
* transformation of image pixel coordinates), and vice versa.  They are based
* on the spcprm struct which contains all information needed for the
* computations.  The struct contains some members that must be set by the
* user, and others that are maintained by these routines, somewhat like a
* C++ class but with no encapsulation.
*
* Routine spcini() is provided to initialize the spcprm struct with default
* values, spcfree() reclaims any memory that may have been allocated to store
* an error message, and spcprt() prints its contents.
*
* A setup routine, spcset(), computes intermediate values in the spcprm struct
* from parameters in it that were supplied by the user.  The struct always
* needs to be set up by spcset() but it need not be called explicitly - refer
* to the explanation of spcprm::flag.
*
* spcx2s() and spcs2x() implement the WCS spectral coordinate transformations.
* In fact, they are high level driver routines for the lower level spectral
* coordinate transformation routines described in spx.h.
*
* A number of routines are provided to aid in analysing or synthesising sets
* of FITS spectral axis keywords:
*
*   - spctype() checks a spectral CTYPEia keyword for validity and returns
*     information derived from it.
*
*   - Spectral keyword analysis routine spcspxe() computes the values of the
*     X-type spectral variables for the S-type variables supplied.
*
*   - Spectral keyword synthesis routine, spcxpse(), computes the S-type
*     variables for the X-types supplied.
*
*   - Given a set of spectral keywords, a translation routine, spctrne(),
*     produces the corresponding set for the specified spectral CTYPEia.
*
*   - spcaips() translates AIPS-convention spectral CTYPEia and VELREF
*     keyvalues.
*
* Spectral variable types - S, P, and X:
* --------------------------------------
* A few words of explanation are necessary regarding spectral variable types
* in FITS.
*
* Every FITS spectral axis has three associated spectral variables:
*
*   S-type: the spectral variable in which coordinates are to be
*     expressed.  Each S-type is encoded as four characters and is
*     linearly related to one of four basic types as follows:
*
*     F: frequency
*       'FREQ':  frequency
*       'AFRQ':  angular frequency
*       'ENER':  photon energy
*       'WAVN':  wave number
*       'VRAD':  radio velocity
*
*     W: wavelength in vacuo
*       'WAVE':  wavelength
*       'VOPT':  optical velocity
*       'ZOPT':  redshift
*
*     A: wavelength in air
*       'AWAV':  wavelength in air
*
*     V: velocity
*       'VELO':  relativistic velocity
*       'BETA':  relativistic beta factor
*
*     The S-type forms the first four characters of the CTYPEia keyvalue,
*     and CRVALia and CDELTia are expressed as S-type quantities so that
*     they provide a first-order approximation to the S-type variable at
*     the reference point.
*
*     Note that 'AFRQ', angular frequency, is additional to the variables
*     defined in WCS Paper III.
*
*   P-type: the basic spectral variable (F, W, A, or V) with which the
*     S-type variable is associated (see list above).
*
*     For non-grism axes, the P-type is encoded as the eighth character of
*     CTYPEia.
*
*   X-type: the basic spectral variable (F, W, A, or V) for which the
*     spectral axis is linear, grisms excluded (see below).
*
*     For non-grism axes, the X-type is encoded as the sixth character of
*     CTYPEia.
*
*   Grisms: Grism axes have normal S-, and P-types but the axis is linear,
*     not in any spectral variable, but in a special "grism parameter".
*     The X-type spectral variable is either W or A for grisms in vacuo or
*     air respectively, but is encoded as 'w' or 'a' to indicate that an
*     additional transformation is required to convert to or from the
*     grism parameter.  The spectral algorithm code for grisms also has a
*     special encoding in CTYPEia, either 'GRI' (in vacuo) or 'GRA' (in air).
*
* In the algorithm chain, the non-linear transformation occurs between the
* X-type and the P-type variables; the transformation between P-type and
* S-type variables is always linear.
*
* When the P-type and X-type variables are the same, the spectral axis is
* linear in the S-type variable and the second four characters of CTYPEia
* are blank.  This can never happen for grism axes.
*
* As an example, correlating radio spectrometers always produce spectra that
* are regularly gridded in frequency; a redshift scale on such a spectrum is
* non-linear.  The required value of CTYPEia would be 'ZOPT-F2W', where the
* desired S-type is 'ZOPT' (redshift), the P-type is necessarily 'W'
* (wavelength), and the X-type is 'F' (frequency) by the nature of the
* instrument.
*
* Argument checking:
* ------------------
* The input spectral values are only checked for values that would result in
* floating point exceptions.  In particular, negative frequencies and
* wavelengths are allowed, as are velocities greater than the speed of
* light.  The same is true for the spectral parameters - rest frequency and
* wavelength.
*
* Accuracy:
* ---------
* No warranty is given for the accuracy of these routines (refer to the
* copyright notice); intending users must satisfy for themselves their
* adequacy for the intended purpose.  However, closure effectively to within
* double precision rounding error was demonstrated by test routine tspc.c
* which accompanies this software.
*
*
* spcini() - Default constructor for the spcprm struct
* ----------------------------------------------------
* spcini() sets all members of a spcprm struct to default values.  It should
* be used to initialize every spcprm struct.
*
* Given and returned:
*   spc       struct spcprm*
*                       Spectral transformation parameters.
*
* Function return value:
*             int       Status return value:
*                         0: Success.
*                         1: Null spcprm pointer passed.
*
*
* spcfree() - Destructor for the spcprm struct
* --------------------------------------------
* spcfree() frees any memory that may have been allocated to store an error
* message in the spcprm struct.
*
* Given:
*   spc       struct spcprm*
*                       Spectral transformation parameters.
*
* Function return value:
*             int       Status return value:
*                         0: Success.
*                         1: Null spcprm pointer passed.
*
*
* spcprt() - Print routine for the spcprm struct
* ----------------------------------------------
* spcprt() prints the contents of a spcprm struct using wcsprintf().  Mainly
* intended for diagnostic purposes.
*
* Given:
*   spc       const struct spcprm*
*                       Spectral transformation parameters.
*
* Function return value:
*             int       Status return value:
*                         0: Success.
*                         1: Null spcprm pointer passed.
*
*
* spcset() - Setup routine for the spcprm struct
* ----------------------------------------------
* spcset() sets up a spcprm struct according to information supplied within
* it.
*
* Note that this routine need not be called directly; it will be invoked by
* spcx2s() and spcs2x() if spcprm::flag is anything other than a predefined
* magic value.
*
* Given and returned:
*   spc       struct spcprm*
*                       Spectral transformation parameters.
*
* Function return value:
*             int       Status return value:
*                         0: Success.
*                         1: Null spcprm pointer passed.
*                         2: Invalid spectral parameters.
*
*                       For returns > 1, a detailed error message is set in
*                       spcprm::err if enabled, see wcserr_enable().
*
*
* spcx2s() - Transform to spectral coordinates
* --------------------------------------------
* spcx2s() transforms intermediate world coordinates to spectral coordinates.
*
* Given and returned:
*   spc       struct spcprm*
*                       Spectral transformation parameters.
*
* Given:
*   nx        int       Vector length.
*
*   sx        int       Vector stride.
*
*   sspec     int       Vector stride.
*
*   x         const double[]
*                       Intermediate world coordinates, in SI units.
*
* Returned:
*   spec      double[]  Spectral coordinates, in SI units.
*
*   stat      int[]     Status return value status for each vector element:
*                         0: Success.
*                         1: Invalid value of x.
*
* Function return value:
*             int       Status return value:
*                         0: Success.
*                         1: Null spcprm pointer passed.
*                         2: Invalid spectral parameters.
*                         3: One or more of the x coordinates were invalid,
*                            as indicated by the stat vector.
*
*                       For returns > 1, a detailed error message is set in
*                       spcprm::err if enabled, see wcserr_enable().
*
*
* spcs2x() - Transform spectral coordinates
* -----------------------------------------
* spcs2x() transforms spectral world coordinates to intermediate world
* coordinates.
*
* Given and returned:
*   spc       struct spcprm*
*                       Spectral transformation parameters.
*
* Given:
*   nspec     int       Vector length.
*
*   sspec     int       Vector stride.
*
*   sx        int       Vector stride.
*
*   spec      const double[]
*                       Spectral coordinates, in SI units.
*
* Returned:
*   x         double[]  Intermediate world coordinates, in SI units.
*
*   stat      int[]     Status return value status for each vector element:
*                         0: Success.
*                         1: Invalid value of spec.
*
* Function return value:
*             int       Status return value:
*                         0: Success.
*                         1: Null spcprm pointer passed.
*                         2: Invalid spectral parameters.
*                         4: One or more of the spec coordinates were
*                            invalid, as indicated by the stat vector.
*
*                       For returns > 1, a detailed error message is set in
*                       spcprm::err if enabled, see wcserr_enable().
*
*
* spctype() - Spectral CTYPEia keyword analysis
* ---------------------------------------------
* spctype() checks whether a CTYPEia keyvalue is a valid spectral axis type
* and if so returns information derived from it relating to the associated S-,
* P-, and X-type spectral variables (see explanation above).
*
* The return arguments are guaranteed not be modified if CTYPEia is not a
* valid spectral type; zero-pointers may be specified for any that are not of
* interest.
*
* A deprecated form of this function, spctyp(), lacks the wcserr** parameter.
*
* Given:
*   ctype     const char[9]
*                       The CTYPEia keyvalue, (eight characters with null
*                       termination).
*
* Returned:
*   stype     char[]    The four-letter name of the S-type spectral variable
*                       copied or translated from ctype.  If a non-zero
*                       pointer is given, the array must accomodate a null-
*                       terminated string of length 5.
*
*   scode     char[]    The three-letter spectral algorithm code copied or
*                       translated from ctype.  Logarithmic ('LOG') and
*                       tabular ('TAB') codes are also recognized.  If a
*                       non-zero pointer is given, the array must accomodate a
*                       null-terminated string of length 4.
*
*   sname     char[]    Descriptive name of the S-type spectral variable.
*                       If a non-zero pointer is given, the array must
*                       accomodate a null-terminated string of length 22.
*
*   units     char[]    SI units of the S-type spectral variable.  If a
*                       non-zero pointer is given, the array must accomodate a
*                       null-terminated string of length 8.
*
*   ptype     char*     Character code for the P-type spectral variable
*                       derived from ctype, one of 'F', 'W', 'A', or 'V'.
*
*   xtype     char*     Character code for the X-type spectral variable
*                       derived from ctype, one of 'F', 'W', 'A', or 'V'.
*                       Also, 'w' and 'a' are synonymous to 'W' and 'A' for
*                       grisms in vacuo and air respectively.  Set to 'L' or
*                       'T' for logarithmic ('LOG') and tabular ('TAB') axes.
*
*   restreq   int*      Multivalued flag that indicates whether rest
*                       frequency or wavelength is required to compute
*                       spectral variables for this CTYPEia:
*                         0: Not required.
*                         1: Required for the conversion between S- and
*                            P-types (e.g. 'ZOPT-F2W').
*                         2: Required for the conversion between P- and
*                            X-types (e.g. 'BETA-W2V').
*                         3: Required for the conversion between S- and
*                            P-types, and between P- and X-types, but not
*                            between S- and X-types (this applies only for
*                            'VRAD-V2F', 'VOPT-V2W', and 'ZOPT-V2W').
*                        Thus the rest frequency or wavelength is required for
*                        spectral coordinate computations (i.e. between S- and
*                        X-types) only if restreq%3 != 0.
*
*   err       struct wcserr **
*                       If enabled, for function return values > 1, this
*                       struct will contain a detailed error message, see
*                       wcserr_enable().  May be NULL if an error message is
*                       not desired.  Otherwise, the user is responsible for
*                       deleting the memory allocated for the wcserr struct.
*
* Function return value:
*             int       Status return value:
*                         0: Success.
*                         2: Invalid spectral parameters (not a spectral
*                            CTYPEia).
*
*
* spcspxe() - Spectral keyword analysis
* ------------------------------------
* spcspxe() analyses the CTYPEia and CRVALia FITS spectral axis keyword values
* and returns information about the associated X-type spectral variable.
*
* A deprecated form of this function, spcspx(), lacks the wcserr** parameter.
*
* Given:
*   ctypeS    const char[9]
*                       Spectral axis type, i.e. the CTYPEia keyvalue, (eight
*                       characters with null termination).  For non-grism
*                       axes, the character code for the P-type spectral
*                       variable in the algorithm code (i.e. the eighth
*                       character of CTYPEia) may be set to '?' (it will not
*                       be reset).
*
*   crvalS    double    Value of the S-type spectral variable at the reference
*                       point, i.e. the CRVALia keyvalue, SI units.
*
*   restfrq,
*   restwav   double    Rest frequency [Hz] and rest wavelength in vacuo [m],
*                       only one of which need be given, the other should be
*                       set to zero.
*
* Returned:
*   ptype     char*     Character code for the P-type spectral variable
*                       derived from ctypeS, one of 'F', 'W', 'A', or 'V'.
*
*   xtype     char*     Character code for the X-type spectral variable
*                       derived from ctypeS, one of 'F', 'W', 'A', or 'V'.
*                       Also, 'w' and 'a' are synonymous to 'W' and 'A' for
*                       grisms in vacuo and air respectively; crvalX and dXdS
*                       (see below) will conform to these.
*
*   restreq   int*      Multivalued flag that indicates whether rest frequency
*                       or wavelength is required to compute spectral
*                       variables for this CTYPEia, as for spctype().
*
*   crvalX    double*   Value of the X-type spectral variable at the reference
*                       point, SI units.
*
*   dXdS      double*   The derivative, dX/dS, evaluated at the reference
*                       point, SI units.  Multiply the CDELTia keyvalue by
*                       this to get the pixel spacing in the X-type spectral
*                       coordinate.
*
*   err       struct wcserr **
*                       If enabled, for function return values > 1, this
*                       struct will contain a detailed error message, see
*                       wcserr_enable().  May be NULL if an error message is
*                       not desired.  Otherwise, the user is responsible for
*                       deleting the memory allocated for the wcserr struct.
*
* Function return value:
*             int       Status return value:
*                         0: Success.
*                         2: Invalid spectral parameters.
*
*
* spcxpse() - Spectral keyword synthesis
* -------------------------------------
* spcxpse(), for the spectral axis type specified and the value provided for
* the X-type spectral variable at the reference point, deduces the value of
* the FITS spectral axis keyword CRVALia and also the derivative dS/dX which
* may be used to compute CDELTia.  See above for an explanation of the S-,
* P-, and X-type spectral variables.
*
* A deprecated form of this function, spcxps(), lacks the wcserr** parameter.
*
* Given:
*   ctypeS    const char[9]
*                       The required spectral axis type, i.e. the CTYPEia
*                       keyvalue, (eight characters with null termination).
*                       For non-grism axes, the character code for the P-type
*                       spectral variable in the algorithm code (i.e. the
*                       eighth character of CTYPEia) may be set to '?' (it
*                       will not be reset).
*
*   crvalX    double    Value of the X-type spectral variable at the reference
*                       point (N.B. NOT the CRVALia keyvalue), SI units.
*
*   restfrq,
*   restwav   double    Rest frequency [Hz] and rest wavelength in vacuo [m],
*                       only one of which need be given, the other should be
*                       set to zero.
*
* Returned:
*   ptype     char*     Character code for the P-type spectral variable
*                       derived from ctypeS, one of 'F', 'W', 'A', or 'V'.
*
*   xtype     char*     Character code for the X-type spectral variable
*                       derived from ctypeS, one of 'F', 'W', 'A', or 'V'.
*                       Also, 'w' and 'a' are synonymous to 'W' and 'A' for
*                       grisms; crvalX and cdeltX must conform to these.
*
*   restreq   int*      Multivalued flag that indicates whether rest frequency
*                       or wavelength is required to compute spectral
*                       variables for this CTYPEia, as for spctype().
*
*   crvalS    double*   Value of the S-type spectral variable at the reference
*                       point (i.e. the appropriate CRVALia keyvalue), SI
*                       units.
*
*   dSdX      double*   The derivative, dS/dX, evaluated at the reference
*                       point, SI units.  Multiply this by the pixel spacing
*                       in the X-type spectral coordinate to get the CDELTia
*                       keyvalue.
*
*   err       struct wcserr **
*                       If enabled, for function return values > 1, this
*                       struct will contain a detailed error message, see
*                       wcserr_enable().  May be NULL if an error message is
*                       not desired.  Otherwise, the user is responsible for
*                       deleting the memory allocated for the wcserr struct.
*
* Function return value:
*             int       Status return value:
*                         0: Success.
*                         2: Invalid spectral parameters.
*
*
* spctrne() - Spectral keyword translation
* ---------------------------------------
* spctrne() translates a set of FITS spectral axis keywords into the
* corresponding set for the specified spectral axis type.  For example, a
* 'FREQ' axis may be translated into 'ZOPT-F2W' and vice versa.
*
* A deprecated form of this function, spctrn(), lacks the wcserr** parameter.
*
* Given:
*   ctypeS1   const char[9]
*                       Spectral axis type, i.e. the CTYPEia keyvalue, (eight
*                       characters with null termination).  For non-grism
*                       axes, the character code for the P-type spectral
*                       variable in the algorithm code (i.e. the eighth
*                       character of CTYPEia) may be set to '?' (it will not
*                       be reset).
*
*   crvalS1   double    Value of the S-type spectral variable at the reference
*                       point, i.e. the CRVALia keyvalue, SI units.
*
*   cdeltS1   double    Increment of the S-type spectral variable at the
*                       reference point, SI units.
*
*   restfrq,
*   restwav   double    Rest frequency [Hz] and rest wavelength in vacuo [m],
*                       only one of which need be given, the other should be
*                       set to zero.  Neither are required if the translation
*                       is between wave-characteristic types, or between
*                       velocity-characteristic types.  E.g., required for
*                       'FREQ'     -> 'ZOPT-F2W', but not required for
*                       'VELO-F2V' -> 'ZOPT-F2W'.
*
* Given and returned:
*   ctypeS2   char[9]   Required spectral axis type (eight characters with
*                       null termination).  The first four characters are
*                       required to be given and are never modified.  The
*                       remaining four, the algorithm code, are completely
*                       determined by, and must be consistent with, ctypeS1
*                       and the first four characters of ctypeS2.  A non-zero
*                       status value will be returned if they are inconsistent
*                       (see below).  However, if the final three characters
*                       are specified as "???", or if just the eighth
*                       character is specified as '?', the correct algorithm
*                       code will be substituted (applies for grism axes as
*                       well as non-grism).
*
* Returned:
*   crvalS2   double*   Value of the new S-type spectral variable at the
*                       reference point, i.e. the new CRVALia keyvalue, SI
*                       units.
*
*   cdeltS2   double*   Increment of the new S-type spectral variable at the
*                       reference point, i.e. the new CDELTia keyvalue, SI
*                       units.
*
*   err       struct wcserr **
*                       If enabled, for function return values > 1, this
*                       struct will contain a detailed error message, see
*                       wcserr_enable().  May be NULL if an error message is
*                       not desired.  Otherwise, the user is responsible for
*                       deleting the memory allocated for the wcserr struct.
*
* Function return value:
*             int       Status return value:
*                         0: Success.
*                         2: Invalid spectral parameters.
*
*                       A status value of 2 will be returned if restfrq or
*                       restwav are not specified when required, or if ctypeS1
*                       or ctypeS2 are self-inconsistent, or have different
*                       spectral X-type variables.
*
*
* spcaips() - Translate AIPS-convention spectral keywords
* -------------------------------------------------------
* spcaips() translates AIPS-convention spectral CTYPEia and VELREF keyvalues.
*
* Given:
*   ctypeA    const char[9]
*                       CTYPEia keyvalue possibly containing an
*                       AIPS-convention spectral code (eight characters, need
*                       not be null-terminated).
*
*   velref    int       AIPS-convention VELREF code.  It has the following
*                       integer values:
*                         1: LSR kinematic, originally described simply as
*                            "LSR" without distinction between the kinematic
*                            and dynamic definitions.
*                         2: Barycentric, originally described as "HEL"
*                            meaning heliocentric.
*                         3: Topocentric, originally described as "OBS"
*                            meaning geocentric but widely interpreted as
*                            topocentric.
*                       AIPS++ extensions to VELREF are also recognized:
*                         4: LSR dynamic.
*                         5: Geocentric.
*                         6: Source rest frame.
*                         7: Galactocentric.
*
*                       For an AIPS 'VELO' axis, a radio convention velocity
*                       (VRAD) is denoted by adding 256 to VELREF, otherwise
*                       an optical velocity (VOPT) is indicated (this is not
*                       applicable to 'FREQ' or 'FELO' axes).  Setting velref
*                       to 0 or 256 chooses between optical and radio velocity
*                       without specifying a Doppler frame, provided that a
*                       frame is encoded in ctypeA.  If not, i.e. for
*                       ctypeA = 'VELO', ctype will be returned as 'VELO'.
*
*                       VELREF takes precedence over CTYPEia in defining the
*                       Doppler frame, e.g.
*
=                         ctypeA = 'VELO-HEL'
=                         velref = 1
*
*                       returns ctype = 'VOPT' with specsys set to 'LSRK'.
*
* Returned:
*   ctype     char[9]   Translated CTYPEia keyvalue, or a copy of ctypeA if no
*                       translation was performed (in which case any trailing
*                       blanks in ctypeA will be replaced with nulls).
*
*   specsys   char[9]   Doppler reference frame indicated by VELREF or else
*                       by CTYPEia with value corresponding to the SPECSYS
*                       keyvalue in the FITS WCS standard.  May be returned
*                       blank if neither specifies a Doppler frame, e.g.
*                       ctypeA = 'FELO' and velref%256 == 0.
*
* Function return value:
*             int       Status return value:
*                        -1: No translation required (not an error).
*                         0: Success.
*                         2: Invalid value of VELREF.
*
*
* spcprm struct - Spectral transformation parameters
* --------------------------------------------------
* The spcprm struct contains information required to transform spectral
* coordinates.  It consists of certain members that must be set by the user
* ("given") and others that are set by the WCSLIB routines ("returned").  Some
* of the latter are supplied for informational purposes while others are for
* internal use only.
*
*   int flag
*     (Given and returned) This flag must be set to zero whenever any of the
*     following spcprm structure members are set or changed:
*
*       - spcprm::type,
*       - spcprm::code,
*       - spcprm::crval,
*       - spcprm::restfrq,
*       - spcprm::restwav,
*       - spcprm::pv[].
*
*     This signals the initialization routine, spcset(), to recompute the
*     returned members of the spcprm struct.  spcset() will reset flag to
*     indicate that this has been done.
*
*   char type[8]
*     (Given) Four-letter spectral variable type, e.g "ZOPT" for
*     CTYPEia = 'ZOPT-F2W'.  (Declared as char[8] for alignment reasons.)
*
*   char code[4]
*     (Given) Three-letter spectral algorithm code, e.g "F2W" for
*     CTYPEia = 'ZOPT-F2W'.
*
*   double crval
*     (Given) Reference value (CRVALia), SI units.
*
*   double restfrq
*     (Given) The rest frequency [Hz], and ...
*
*   double restwav
*     (Given) ... the rest wavelength in vacuo [m], only one of which need be
*     given, the other should be set to zero.  Neither are required if the
*     X and S spectral variables are both wave-characteristic, or both
*     velocity-characteristic, types.
*
*   double pv[7]
*     (Given) Grism parameters for 'GRI' and 'GRA' algorithm codes:
*       - 0: G, grating ruling density.
*       - 1: m, interference order.
*       - 2: alpha, angle of incidence [deg].
*       - 3: n_r, refractive index at the reference wavelength, lambda_r.
*       - 4: n'_r, dn/dlambda at the reference wavelength, lambda_r (/m).
*       - 5: epsilon, grating tilt angle [deg].
*       - 6: theta, detector tilt angle [deg].
*
* The remaining members of the spcprm struct are maintained by spcset() and
* must not be modified elsewhere:
*
*   double w[6]
*     (Returned) Intermediate values:
*       - 0: Rest frequency or wavelength (SI).
*       - 1: The value of the X-type spectral variable at the reference point
*           (SI units).
*       - 2: dX/dS at the reference point (SI units).
*      The remainder are grism intermediates.
*
*   int isGrism
*     (Returned) Grism coordinates?
*       - 0: no,
*       - 1: in vacuum,
*       - 2: in air.
*
*   int padding1
*     (An unused variable inserted for alignment purposes only.)
*
*   struct wcserr *err
*     (Returned) If enabled, when an error status is returned, this struct
*     contains detailed information about the error, see wcserr_enable().
*
*   void *padding2
*     (An unused variable inserted for alignment purposes only.)
*   int (*spxX2P)(SPX_ARGS)
*     (Returned) The first and ...
*   int (*spxP2S)(SPX_ARGS)
*     (Returned) ... the second of the pointers to the transformation
*     functions in the two-step algorithm chain X -> P -> S in the
*     pixel-to-spectral direction where the non-linear transformation is from
*     X to P.  The argument list, SPX_ARGS, is defined in spx.h.
*
*   int (*spxS2P)(SPX_ARGS)
*     (Returned) The first and ...
*   int (*spxP2X)(SPX_ARGS)
*     (Returned) ... the second of the pointers to the transformation
*     functions in the two-step algorithm chain S -> P -> X in the
*     spectral-to-pixel direction where the non-linear transformation is from
*     P to X.  The argument list, SPX_ARGS, is defined in spx.h.
*
*
* Global variable: const char *spc_errmsg[] - Status return messages
* ------------------------------------------------------------------
* Error messages to match the status value returned from each function.
*
*===========================================================================*/

#ifndef WCSLIB_SPC
#define WCSLIB_SPC

#include "spx.h"

#ifdef __cplusplus
extern "C" {
#endif


extern const char *spc_errmsg[];

enum spc_errmsg_enum {
  SPCERR_NO_CHANGE       = -1,	/* No change. */
  SPCERR_SUCCESS         =  0,	/* Success. */
  SPCERR_NULL_POINTER    =  1,	/* Null spcprm pointer passed. */
  SPCERR_BAD_SPEC_PARAMS =  2,	/* Invalid spectral parameters. */
  SPCERR_BAD_X           =  3,	/* One or more of x coordinates were
				   invalid. */
  SPCERR_BAD_SPEC        =  4 	/* One or more of the spec coordinates were
				   invalid. */
};

struct spcprm {
  /* Initialization flag (see the prologue above).                          */
  /*------------------------------------------------------------------------*/
  int    flag;			/* Set to zero to force initialization.     */

  /* Parameters to be provided (see the prologue above).                    */
  /*------------------------------------------------------------------------*/
  char   type[8];		/* Four-letter spectral variable type.      */
  char   code[4];		/* Three-letter spectral algorithm code.    */

  double crval;			/* Reference value (CRVALia), SI units.     */
  double restfrq;		/* Rest frequency, Hz.                      */
  double restwav;		/* Rest wavelength, m.                      */

  double pv[7];			/* Grism parameters:                        */
				/*   0: G, grating ruling density.          */
				/*   1: m, interference order.              */
				/*   2: alpha, angle of incidence.          */
				/*   3: n_r, refractive index at lambda_r.  */
				/*   4: n'_r, dn/dlambda at lambda_r.       */
				/*   5: epsilon, grating tilt angle.        */
				/*   6: theta, detector tilt angle.         */

  /* Information derived from the parameters supplied.                      */
  /*------------------------------------------------------------------------*/
  double w[6];			/* Intermediate values.                     */
				/*   0: Rest frequency or wavelength (SI).  */
				/*   1: CRVALX (SI units).                  */
				/*   2: CDELTX/CDELTia = dX/dS (SI units).  */
				/* The remainder are grism intermediates.   */

  int    isGrism;		/* Grism coordinates?  1: vacuum, 2: air.   */
  int    padding1;		/* (Dummy inserted for alignment purposes.) */

  /* Error handling                                                         */
  /*------------------------------------------------------------------------*/
  struct wcserr *err;

  /* Private                                                                */
  /*------------------------------------------------------------------------*/
  void   *padding2;		/* (Dummy inserted for alignment purposes.) */
  int (*spxX2P)(SPX_ARGS);	/* Pointers to the transformation functions */
  int (*spxP2S)(SPX_ARGS);	/* in the two-step algorithm chain in the   */
				/* pixel-to-spectral direction.             */

  int (*spxS2P)(SPX_ARGS);	/* Pointers to the transformation functions */
  int (*spxP2X)(SPX_ARGS);	/* in the two-step algorithm chain in the   */
				/* spectral-to-pixel direction.             */
};

/* Size of the spcprm struct in int units, used by the Fortran wrappers. */
#define SPCLEN (sizeof(struct spcprm)/sizeof(int))


int spcini(struct spcprm *spc);

int spcfree(struct spcprm *spc);

int spcprt(const struct spcprm *spc);

int spcset(struct spcprm *spc);

int spcx2s(struct spcprm *spc, int nx, int sx, int sspec,
           const double x[], double spec[], int stat[]);

int spcs2x(struct spcprm *spc, int nspec, int sspec, int sx,
           const double spec[], double x[], int stat[]);

int spctype(const char ctype[9], char stype[], char scode[], char sname[],
            char units[], char *ptype, char *xtype, int *restreq,
            struct wcserr **err);

int spcspxe(const char ctypeS[9], double crvalS, double restfrq,
            double restwav, char *ptype, char *xtype, int *restreq,
            double *crvalX, double *dXdS, struct wcserr **err);

int spcxpse(const char ctypeS[9], double crvalX, double restfrq,
            double restwav, char *ptype, char *xtype, int *restreq,
            double *crvalS, double *dSdX, struct wcserr **err);

int spctrne(const char ctypeS1[9], double crvalS1, double cdeltS1,
            double restfrq, double restwav, char ctypeS2[9], double *crvalS2,
            double *cdeltS2, struct wcserr **err);

int spcaips(const char ctypeA[9], int velref, char ctype[9], char specsys[9]);


/* Deprecated. */
#define spcini_errmsg spc_errmsg
#define spcprt_errmsg spc_errmsg
#define spcset_errmsg spc_errmsg
#define spcx2s_errmsg spc_errmsg
#define spcs2x_errmsg spc_errmsg

int spctyp(const char ctype[9], char stype[], char scode[], char sname[],
           char units[], char *ptype, char *xtype, int *restreq);
int spcspx(const char ctypeS[9], double crvalS, double restfrq,
           double restwav, char *ptype, char *xtype, int *restreq,
           double *crvalX, double *dXdS);
int spcxps(const char ctypeS[9], double crvalX, double restfrq,
           double restwav, char *ptype, char *xtype, int *restreq,
           double *crvalS, double *dSdX);
int spctrn(const char ctypeS1[9], double crvalS1, double cdeltS1,
           double restfrq, double restwav, char ctypeS2[9], double *crvalS2,
           double *cdeltS2);

#ifdef __cplusplus
}
#endif

#endif /* WCSLIB_SPC */
