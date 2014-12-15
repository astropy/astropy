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
  $Id: wcsunits.h,v 4.25 2014/12/14 14:29:36 mcalabre Exp $
*=============================================================================
*
* WCSLIB 4.25 - C routines that implement the FITS World Coordinate System
* (WCS) standard.  Refer to
*
*   "Representations of world coordinates in FITS",
*   Greisen, E.W., & Calabretta, M.R. 2002, A&A, 395, 1061 (Paper I)
*
* The Flexible Image Transport System (FITS), a data format widely used in
* astronomy for data interchange and archive, is described in
*
*   "Definition of the Flexible Image Transport System (FITS), version 3.0",
*   Pence, W.D., Chiappetti, L., Page, C.G., Shaw, R.A., & Stobie, E. 2010,
*   A&A, 524, A42 - http://dx.doi.org/10.1051/0004-6361/201015362
*
* See also http://fits.gsfc.nasa.gov
*
* Refer to the README file provided with WCSLIB for an overview of the
* library.
*
*
* Summary of the wcsunits routines
* --------------------------------
* Routines in this suite deal with units specifications and conversions:
*
*   - wcsunitse(): given two unit specifications, derive the conversion from
*     one to the other.
*
*   - wcsutrne(): translates certain commonly used but non-standard unit
*     strings.  It is intended to be called before wcsulexe() which only
*     handles standard FITS units specifications.
*
*   - wcsulexe(): parses a standard FITS units specification of arbitrary
*     complexity, deriving the conversion to canonical units.
*
*
* wcsunitse() - FITS units specification conversion
* -------------------------------------------------
* wcsunitse() derives the conversion from one system of units to another.
*
* A deprecated form of this function, wcsunits(), lacks the wcserr**
* parameter.
*
* Given:
*   have      const char []
*                       FITS units specification to convert from (null-
*                       terminated), with or without surrounding square
*                       brackets (for inline specifications); text following
*                       the closing bracket is ignored.
*
*   want      const char []
*                       FITS units specification to convert to (null-
*                       terminated), with or without surrounding square
*                       brackets (for inline specifications); text following
*                       the closing bracket is ignored.
*
* Returned:
*   scale,
*   offset,
*   power     double*   Convert units using
*
=                         pow(scale*value + offset, power);
*
*                       Normally offset is zero except for log() or ln()
*                       conversions, e.g. "log(MHz)" to "ln(Hz)".  Likewise,
*                       power is normally unity except for exp() conversions,
*                       e.g. "exp(ms)" to "exp(/Hz)".  Thus conversions
*                       ordinarily consist of
*
=                         value *= scale;
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
*                          0: Success.
*                        1-9: Status return from wcsulexe().
*                         10: Non-conformant unit specifications.
*                         11: Non-conformant functions.
*
*                       scale is zeroed on return if an error occurs.
*
*
* wcsutrne() - Translation of non-standard unit specifications
* ------------------------------------------------------------
* wcsutrne() translates certain commonly used but non-standard unit strings,
* e.g. "DEG", "MHZ", "KELVIN", that are not recognized by wcsulexe(), refer to
* the notes below for a full list.  Compounds are also recognized, e.g.
* "JY/BEAM" and "KM/SEC/SEC".  Extraneous embedded blanks are removed.
*
* A deprecated form of this function, wcsutrn(), lacks the wcserr** parameter.
*
* Given:
*   ctrl      int       Although "S" is commonly used to represent seconds,
*                       its translation to "s" is potentially unsafe since the
*                       standard recognizes "S" formally as Siemens, however
*                       rarely that may be used.  The same applies to "H" for
*                       hours (Henry), and "D" for days (Debye).  This
*                       bit-flag controls what to do in such cases:
*                         1: Translate "S" to "s".
*                         2: Translate "H" to "h".
*                         4: Translate "D" to "d".
*                       Thus ctrl == 0 doesn't do any unsafe translations,
*                       whereas ctrl == 7 does all of them.
*
* Given and returned:
*   unitstr   char []   Null-terminated character array containing the units
*                       specification to be translated.
*
*                       Inline units specifications in the a FITS header
*                       keycomment are also handled.  If the first non-blank
*                       character in unitstr is '[' then the unit string is
*                       delimited by its matching ']'.  Blanks preceding '['
*                       will be stripped off, but text following the closing
*                       bracket will be preserved without modification.
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
*                        -1: No change was made, other than stripping blanks
*                            (not an error).
*                         0: Success.
*                         9: Internal parser error.
*                        12: Potentially unsafe translation, whether applied
*                            or not (see notes).
*
* Notes:
*   Translation of non-standard unit specifications: apart from leading and
*   trailing blanks, a case-sensitive match is required for the aliases listed
*   below, in particular the only recognized aliases with metric prefixes are
*   "KM", "KHZ", "MHZ", and "GHZ".  Potentially unsafe translations of "D",
*   "H", and "S", shown in parentheses, are optional.
*
=     Unit       Recognized aliases
=     ----       -------------------------------------------------------------
=     Angstrom   angstrom
=     arcmin     arcmins, ARCMIN, ARCMINS
=     arcsec     arcsecs, ARCSEC, ARCSECS
=     beam       BEAM
=     byte       Byte
=     d          day, days, (D), DAY, DAYS
=     deg        degree, degrees, DEG, DEGREE, DEGREES
=     GHz        GHZ
=     h          hr, (H), HR
=     Hz         hz, HZ
=     kHz        KHZ
=     Jy         JY
=     K          kelvin, kelvins, Kelvin, Kelvins, KELVIN, KELVINS
=     km         KM
=     m          metre, meter, metres, meters, M, METRE, METER, METRES, METERS
=     min        MIN
=     MHz        MHZ
=     Ohm        ohm
=     Pa         pascal, pascals, Pascal, Pascals, PASCAL, PASCALS
=     pixel      pixels, PIXEL, PIXELS
=     rad        radian, radians, RAD, RADIAN, RADIANS
=     s          sec, second, seconds, (S), SEC, SECOND, SECONDS
=     V          volt, volts, Volt, Volts, VOLT, VOLTS
=     yr         year, years, YR, YEAR, YEARS
*
*   The aliases "angstrom", "ohm", and "Byte" for (Angstrom, Ohm, and byte)
*   are recognized by wcsulexe() itself as an unofficial extension of the
*   standard, but they are converted to the standard form here.
*
*
* wcsulexe() - FITS units specification parser
* --------------------------------------------
* wcsulexe() parses a standard FITS units specification of arbitrary
* complexity, deriving the scale factor required to convert to canonical
* units - basically SI with degrees and "dimensionless" additions such as
* byte, pixel and count.
*
* A deprecated form of this function, wcsulex(), lacks the wcserr** parameter.
*
* Given:
*   unitstr   const char []
*                       Null-terminated character array containing the units
*                       specification, with or without surrounding square
*                       brackets (for inline specifications); text following
*                       the closing bracket is ignored.
*
* Returned:
*   func      int*      Special function type, see note 4:
*                         0: None
*                         1: log()  ...base 10
*                         2: ln()   ...base e
*                         3: exp()
*
*   scale     double*   Scale factor for the unit specification; multiply a
*                       value expressed in the given units by this factor to
*                       convert it to canonical units.
*
*   units     double[WCSUNITS_NTYPE]
*                       A units specification is decomposed into powers of 16
*                       fundamental unit types: angle, mass, length, time,
*                       count, pixel, etc.  Preprocessor macro WCSUNITS_NTYPE
*                       is defined to dimension this vector, and others such
*                       WCSUNITS_PLANE_ANGLE, WCSUNITS_LENGTH, etc. to access
*                       its elements.
*
*                       Corresponding character strings, wcsunits_types[] and
*                       wcsunits_units[], are predefined to describe each
*                       quantity and its canonical units.
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
*                         1: Invalid numeric multiplier.
*                         2: Dangling binary operator.
*                         3: Invalid symbol in INITIAL context.
*                         4: Function in invalid context.
*                         5: Invalid symbol in EXPON context.
*                         6: Unbalanced bracket.
*                         7: Unbalanced parenthesis.
*                         8: Consecutive binary operators.
*                         9: Internal parser error.
*
*                       scale and units[] are zeroed on return if an error
*                       occurs.
*
* Notes:
*   1: wcsulexe() is permissive in accepting whitespace in all contexts in a
*      units specification where it does not create ambiguity (e.g. not
*      between a metric prefix and a basic unit string), including in strings
*      like "log (m ** 2)" which is formally disallowed.
*
*   2: Supported extensions:
*      - "angstrom" (OGIP usage) is allowed in addition to "Angstrom".
*      - "ohm"      (OGIP usage) is allowed in addition to "Ohm".
*      - "Byte"   (common usage) is allowed in addition to "byte".
*
*   3: Table 6 of WCS Paper I lists eleven units for which metric prefixes are
*      allowed.  However, in this implementation only prefixes greater than
*      unity are allowed for "a" (annum), "yr" (year), "pc" (parsec), "bit",
*      and "byte", and only prefixes less than unity are allowed for "mag"
*      (stellar magnitude).
*
*      Metric prefix "P" (peta) is specifically forbidden for "a" (annum) to
*      avoid confusion with "Pa" (Pascal, not peta-annum).  Note that metric
*      prefixes are specifically disallowed for "h" (hour) and "d" (day) so
*      that "ph" (photons) cannot be interpreted as pico-hours, nor "cd"
*      (candela) as centi-days.
*
*   4: Function types log(), ln() and exp() may only occur at the start of the
*      units specification.  The scale and units[] returned for these refers
*      to the string inside the function "argument", e.g. to "MHz" in log(MHz)
*      for which a scale of 1e6 will be returned.
*
*
* Global variable: const char *wcsunits_errmsg[] - Status return messages
* -----------------------------------------------------------------------
* Error messages to match the status value returned from each function.
*
*
* Global variable: const char *wcsunits_types[] - Names of physical quantities
* ----------------------------------------------------------------------------
* Names for physical quantities to match the units vector returned by
* wcsulexe():
*   -  0: plane angle
*   -  1: solid angle
*   -  2: charge
*   -  3: mole
*   -  4: temperature
*   -  5: luminous intensity
*   -  6: mass
*   -  7: length
*   -  8: time
*   -  9: beam
*   - 10: bin
*   - 11: bit
*   - 12: count
*   - 13: stellar magnitude
*   - 14: pixel
*   - 15: solar ratio
*   - 16: voxel
*
*
* Global variable: const char *wcsunits_units[] - Names of units
* --------------------------------------------------------------
* Names for the units (SI) to match the units vector returned by wcsulexe():
*   -  0: degree
*   -  1: steradian
*   -  2: Coulomb
*   -  3: mole
*   -  4: Kelvin
*   -  5: candela
*   -  6: kilogram
*   -  7: metre
*   -  8: second
*
* The remainder are dimensionless.
*===========================================================================*/

#ifndef WCSLIB_WCSUNITS
#define WCSLIB_WCSUNITS

#include "wcserr.h"

#ifdef __cplusplus
extern "C" {
#endif


extern const char *wcsunits_errmsg[];

enum wcsunits_errmsg_enum {
  UNITSERR_SUCCESS            =  0,	/* Success. */
  UNITSERR_BAD_NUM_MULTIPLIER =  1,	/* Invalid numeric multiplier. */
  UNITSERR_DANGLING_BINOP     =  2,	/* Dangling binary operator. */
  UNITSERR_BAD_INITIAL_SYMBOL =  3,	/* Invalid symbol in INITIAL
					   context. */
  UNITSERR_FUNCTION_CONTEXT   =  4,	/* Function in invalid context. */
  UNITSERR_BAD_EXPON_SYMBOL   =  5,	/* Invalid symbol in EXPON context. */
  UNITSERR_UNBAL_BRACKET      =  6,	/* Unbalanced bracket. */
  UNITSERR_UNBAL_PAREN        =  7,	/* Unbalanced parenthesis. */
  UNITSERR_CONSEC_BINOPS      =  8,	/* Consecutive binary operators. */
  UNITSERR_PARSER_ERROR       =  9,	/* Internal parser error. */
  UNITSERR_BAD_UNIT_SPEC      = 10,	/* Non-conformant unit
					   specifications. */
  UNITSERR_BAD_FUNCS          = 11,	/* Non-conformant functions. */
  UNITSERR_UNSAFE_TRANS       = 12	/* Potentially unsafe translation. */
};

extern const char *wcsunits_types[];
extern const char *wcsunits_units[];

#define WCSUNITS_PLANE_ANGLE 0
#define WCSUNITS_SOLID_ANGLE 1
#define WCSUNITS_CHARGE      2
#define WCSUNITS_MOLE        3
#define WCSUNITS_TEMPERATURE 4
#define WCSUNITS_LUMINTEN    5
#define WCSUNITS_MASS        6
#define WCSUNITS_LENGTH      7
#define WCSUNITS_TIME        8
#define WCSUNITS_BEAM        9
#define WCSUNITS_BIN        10
#define WCSUNITS_BIT        11
#define WCSUNITS_COUNT      12
#define WCSUNITS_MAGNITUDE  13
#define WCSUNITS_PIXEL      14
#define WCSUNITS_SOLRATIO   15
#define WCSUNITS_VOXEL      16

#define WCSUNITS_NTYPE      17


int wcsunitse(const char have[], const char want[], double *scale,
              double *offset, double *power, struct wcserr **err);

int wcsutrne(int ctrl, char unitstr[], struct wcserr **err);

int wcsulexe(const char unitstr[], int *func, double *scale,
             double units[WCSUNITS_NTYPE], struct wcserr **err);

/* Deprecated. */
int wcsunits(const char have[], const char want[], double *scale,
             double *offset, double *power);
int wcsutrn(int ctrl, char unitstr[]);
int wcsulex(const char unitstr[], int *func, double *scale,
            double units[WCSUNITS_NTYPE]);

#ifdef __cplusplus
}
#endif

#endif /* WCSLIB_WCSUNITS */
