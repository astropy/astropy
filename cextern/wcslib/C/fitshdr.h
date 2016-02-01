/*============================================================================

  WCSLIB 5.10 - an implementation of the FITS WCS standard.
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
  $Id: fitshdr.h,v 5.10 2015/10/09 08:19:15 mcalabre Exp $
*=============================================================================
*
* WCSLIB 5.10 - C routines that implement the FITS World Coordinate System
* (WCS) standard.  Refer to the README file provided with WCSLIB for an
* overview of the library.
*
*
* Summary of the fitshdr routines
* -------------------------------
* The Flexible Image Transport System (FITS), is a data format widely used in
* astronomy for data interchange and archive.  It is described in
*
=   "Definition of the Flexible Image Transport System (FITS), version 3.0",
=   Pence, W.D., Chiappetti, L., Page, C.G., Shaw, R.A., & Stobie, E. 2010,
=   A&A, 524, A42 - http://dx.doi.org/10.1051/0004-6361/201015362
*
* See also http://fits.gsfc.nasa.gov
*
* fitshdr() is a generic FITS header parser provided to handle keyrecords that
* are ignored by the WCS header parsers, wcspih() and wcsbth().  Typically the
* latter may be set to remove WCS keyrecords from a header leaving fitshdr()
* to handle the remainder.
*
*
* fitshdr() - FITS header parser routine
* --------------------------------------
* fitshdr() parses a character array containing a FITS header, extracting
* all keywords and their values into an array of fitskey structs.
*
* Given:
*   header    const char []
*                       Character array containing the (entire) FITS header,
*                       for example, as might be obtained conveniently via the
*                       CFITSIO routine fits_hdr2str().
*
*                       Each header "keyrecord" (formerly "card image")
*                       consists of exactly 80 7-bit ASCII printing characters
*                       in the range 0x20 to 0x7e (which excludes NUL, BS,
*                       TAB, LF, FF and CR) especially noting that the
*                       keyrecords are NOT null-terminated.
*
*   nkeyrec   int       Number of keyrecords in header[].
*
*   nkeyids   int       Number of entries in keyids[].
*
* Given and returned:
*   keyids    struct fitskeyid []
*                       While all keywords are extracted from the header,
*                       keyids[] provides a convienient way of indexing them.
*                       The fitskeyid struct contains three members;
*                       fitskeyid::name must be set by the user while
*                       fitskeyid::count and fitskeyid::idx are returned by
*                       fitshdr().  All matched keywords will have their
*                       fitskey::keyno member negated.
*
* Returned:
*   nreject   int*      Number of header keyrecords rejected for syntax
*                       errors.
*
*   keys      struct fitskey**
*                       Pointer to an array of nkeyrec fitskey structs
*                       containing all keywords and keyvalues extracted from
*                       the header.
*
*                       Memory for the array is allocated by fitshdr() and
*                       this must be freed by the user by invoking free() on
*                       the array.
*
* Function return value:
*             int       Status return value:
*                         0: Success.
*                         1: Null fitskey pointer passed.
*                         2: Memory allocation failed.
*                         3: Fatal error returned by Flex parser.
*
* Notes:
*   1: Keyword parsing is done in accordance with the syntax defined by
*      NOST 100-2.0, noting the following points in particular:
*
*      a: Sect. 5.1.2.1 specifies that keywords be left-justified in columns
*         1-8, blank-filled with no embedded spaces, composed only of the
*         ASCII characters ABCDEFGHJKLMNOPQRSTUVWXYZ0123456789-_
*
*         fitshdr() accepts any characters in columns 1-8 but flags keywords
*         that do not conform to standard syntax.
*
*      b: Sect. 5.1.2.2 defines the "value indicator" as the characters "= "
*         occurring in columns 9 and 10.  If these are absent then the
*         keyword has no value and columns 9-80 may contain any ASCII text
*         (but see note 2 for CONTINUE keyrecords).  This is copied to the
*         comment member of the fitskey struct.
*
*      c: Sect. 5.1.2.3 states that a keyword may have a null (undefined)
*         value if the value/comment field, columns 11-80, consists entirely
*         of spaces, possibly followed by a comment.
*
*      d: Sect. 5.1.1 states that trailing blanks in a string keyvalue are
*         not significant and the parser always removes them.  A string
*         containing nothing but blanks will be replaced with a single
*         blank.
*
*         Sect. 5.2.1 also states that a quote character (') in a string
*         value is to be represented by two successive quote characters and
*         the parser removes the repeated quote.
*
*      e: The parser recognizes free-format character (NOST 100-2.0,
*         Sect. 5.2.1), integer (Sect. 5.2.3), and floating-point values
*         (Sect. 5.2.4) for all keywords.
*
*      f: Sect. 5.2.3 offers no comment on the size of an integer keyvalue
*         except indirectly in limiting it to 70 digits.  The parser will
*         translate an integer keyvalue to a 32-bit signed integer if it
*         lies in the range -2147483648 to +2147483647, otherwise it
*         interprets it as a 64-bit signed integer if possible, or else a
*         "very long" integer (see fitskey::type).
*
*      g: END not followed by 77 blanks is not considered to be a legitimate
*         end keyrecord.
*
*   2: The parser supports a generalization of the OGIP Long String Keyvalue
*      Convention (v1.0) whereby strings may be continued onto successive
*      header keyrecords.  A keyrecord contains a segment of a continued
*      string if and only if
*
*      a: it contains the pseudo-keyword CONTINUE,
*
*      b: columns 9 and 10 are both blank,
*
*      c: columns 11 to 80 contain what would be considered a valid string
*         keyvalue, including optional keycomment, if column 9 had contained
*         '=',
*
*      d: the previous keyrecord contained either a valid string keyvalue or
*         a valid CONTINUE keyrecord.
*
*      If any of these conditions is violated, the keyrecord is considered in
*      isolation.
*
*      Syntax errors in keycomments in a continued string are treated more
*      permissively than usual; the '/' delimiter may be omitted provided that
*      parsing of the string keyvalue is not compromised.  However, the
*      FITSHDR_COMMENT status bit will be set for the keyrecord (see
*      fitskey::status).
*
*      As for normal strings, trailing blanks in a continued string are not
*      significant.
*
*      In the OGIP convention "the '&' character is used as the last non-blank
*      character of the string to indicate that the string is (probably)
*      continued on the following keyword".  This additional syntax is not
*      required by fitshdr(), but if '&' does occur as the last non-blank
*      character of a continued string keyvalue then it will be removed, along
*      with any trailing blanks.  However, blanks that occur before the '&'
*      will be preserved.
*
*
* fitskeyid struct - Keyword indexing
* -----------------------------------
* fitshdr() uses the fitskeyid struct to return indexing information for
* specified keywords.  The struct contains three members, the first of which,
* fitskeyid::name, must be set by the user with the remainder returned by
* fitshdr().
*
*   char name[12]:
*     (Given) Name of the required keyword.  This is to be set by the user;
*     the '.' character may be used for wildcarding.  Trailing blanks will be
*     replaced with nulls.
*
*   int count:
*     (Returned) The number of matches found for the keyword.
*
*   int idx[2]:
*     (Returned) Indices into keys[], the array of fitskey structs returned by
*     fitshdr().  Note that these are 0-relative array indices, not keyrecord
*     numbers.
*
*     If the keyword is found in the header the first index will be set to the
*     array index of its first occurrence, otherwise it will be set to -1.
*
*     If multiples of the keyword are found, the second index will be set to
*     the array index of its last occurrence, otherwise it will be set to -1.
*
*
* fitskey struct - Keyword/value information
* ------------------------------------------
* fitshdr() returns an array of fitskey structs, each of which contains the
* result of parsing one FITS header keyrecord.  All members of the fitskey
* struct are returned by fitshdr(), none are given by the user.
*
*   int keyno
*     (Returned) Keyrecord number (1-relative) in the array passed as input to
*     fitshdr().  This will be negated if the keyword matched any specified in
*     the keyids[] index.
*
*   int keyid
*     (Returned) Index into the first entry in keyids[] with which the
*     keyrecord matches, else -1.
*
*   int status
*     (Returned) Status flag bit-vector for the header keyrecord employing the
*     following bit masks defined as preprocessor macros:
*
*       - FITSHDR_KEYWORD:    Illegal keyword syntax.
*       - FITSHDR_KEYVALUE:   Illegal keyvalue syntax.
*       - FITSHDR_COMMENT:    Illegal keycomment syntax.
*       - FITSHDR_KEYREC:     Illegal keyrecord, e.g. an END keyrecord with
*                             trailing text.
*       - FITSHDR_TRAILER:    Keyrecord following a valid END keyrecord.
*
*     The header keyrecord is syntactically correct if no bits are set.
*
*   char keyword[12]
*     (Returned) Keyword name, null-filled for keywords of less than eight
*     characters (trailing blanks replaced by nulls).
*
*     Use
*
=       sprintf(dst, "%.8s", keyword)
*
*     to copy it to a character array with null-termination, or
*
=       sprintf(dst, "%8.8s", keyword)
*
*     to blank-fill to eight characters followed by null-termination.
*
*   int type
*     (Returned) Keyvalue data type:
*       - 0: No keyvalue (both the value and type are undefined).
*       - 1: Logical, represented as int.
*       - 2: 32-bit signed integer.
*       - 3: 64-bit signed integer (see below).
*       - 4: Very long integer (see below).
*       - 5: Floating point (stored as double).
*       - 6: Integer complex (stored as double[2]).
*       - 7: Floating point complex (stored as double[2]).
*       - 8: String.
*       - 8+10*n: Continued string (described below and in fitshdr() note 2).
*
*     A negative type indicates that a syntax error was encountered when
*     attempting to parse a keyvalue of the particular type.
*
*     Comments on particular data types:
*       - 64-bit signed integers lie in the range
*
=           (-9223372036854775808 <= int64 <  -2147483648) ||
=                    (+2147483647 <  int64 <= +9223372036854775807)
*
*         A native 64-bit data type may be defined via preprocessor macro
*         WCSLIB_INT64 defined in wcsconfig.h, e.g. as 'long long int'; this
*         will be typedef'd to 'int64' here.  If WCSLIB_INT64 is not set, then
*         int64 is typedef'd to int[3] instead and fitskey::keyvalue is to be
*         computed as
*
=           ((keyvalue.k[2]) * 1000000000 +
=             keyvalue.k[1]) * 1000000000 +
=             keyvalue.k[0]
*
*         and may reported via
*
=            if (keyvalue.k[2]) {
=              printf("%d%09d%09d", keyvalue.k[2], abs(keyvalue.k[1]),
=                                   abs(keyvalue.k[0]));
=            } else {
=              printf("%d%09d", keyvalue.k[1], abs(keyvalue.k[0]));
=            }
*
*         where keyvalue.k[0] and keyvalue.k[1] range from -999999999 to
*         +999999999.
*
*       - Very long integers, up to 70 decimal digits in length, are encoded
*         in keyvalue.l as an array of int[8], each of which stores 9 decimal
*         digits.  fitskey::keyvalue is to be computed as
*
=           (((((((keyvalue.l[7]) * 1000000000 +
=                  keyvalue.l[6]) * 1000000000 +
=                  keyvalue.l[5]) * 1000000000 +
=                  keyvalue.l[4]) * 1000000000 +
=                  keyvalue.l[3]) * 1000000000 +
=                  keyvalue.l[2]) * 1000000000 +
=                  keyvalue.l[1]) * 1000000000 +
=                  keyvalue.l[0]
*
*       - Continued strings are not reconstructed, they remain split over
*         successive fitskey structs in the keys[] array returned by
*         fitshdr().  fitskey::keyvalue data type, 8 + 10n, indicates the
*         segment number, n, in the continuation.
*
*   int padding
*     (An unused variable inserted for alignment purposes only.)
*
*   union keyvalue
*     (Returned) A union comprised of
*
*       - fitskey::i,
*       - fitskey::k,
*       - fitskey::l,
*       - fitskey::f,
*       - fitskey::c,
*       - fitskey::s,
*
*     used by the fitskey struct to contain the value associated with a
*     keyword.
*
*   int i
*     (Returned) Logical (fitskey::type == 1) and 32-bit signed integer
*     (fitskey::type == 2) data types in the fitskey::keyvalue union.
*
*   int64 k
*     (Returned) 64-bit signed integer (fitskey::type == 3) data type in the
*     fitskey::keyvalue union.
*
*   int l[8]
*     (Returned) Very long integer (fitskey::type == 4) data type in the
*     fitskey::keyvalue union.
*
*   double f
*     (Returned) Floating point (fitskey::type == 5) data type in the
*     fitskey::keyvalue union.
*
*   double c[2]
*     (Returned) Integer and floating point complex (fitskey::type == 6 || 7)
*     data types in the fitskey::keyvalue union.
*
*   char s[72]
*     (Returned) Null-terminated string (fitskey::type == 8) data type in the
*     fitskey::keyvalue union.
*
*   int ulen
*     (Returned) Where a keycomment contains a units string in the standard
*     form, e.g. [m/s], the ulen member indicates its length, inclusive of
*     square brackets.  Otherwise ulen is zero.
*
*   char comment[84]
*     (Returned) Keycomment, i.e. comment associated with the keyword or, for
*     keyrecords rejected because of syntax errors, the compete keyrecord
*     itself with null-termination.
*
*     Comments are null-terminated with trailing spaces removed.  Leading
*     spaces are also removed from keycomments (i.e. those immediately
*     following the '/' character), but not from COMMENT or HISTORY keyrecords
*     or keyrecords without a value indicator ("= " in columns 9-80).
*
*
* Global variable: const char *fitshdr_errmsg[] - Status return messages
* ----------------------------------------------------------------------
* Error messages to match the status value returned from each function.
*
*===========================================================================*/

#ifndef WCSLIB_FITSHDR
#define WCSLIB_FITSHDR

#include "wcsconfig.h"

#ifdef __cplusplus
extern "C" {
#endif

#define FITSHDR_KEYWORD  0x01
#define FITSHDR_KEYVALUE 0x02
#define FITSHDR_COMMENT  0x04
#define FITSHDR_KEYREC   0x08
#define FITSHDR_CARD     0x08	/* Alias for backwards compatibility. */
#define FITSHDR_TRAILER  0x10


extern const char *fitshdr_errmsg[];

#ifdef WCSLIB_INT64
  typedef WCSLIB_INT64 int64;
#else
  typedef int int64[3];
#endif


/* Struct used for indexing the keywords. */
struct fitskeyid {
  char name[12];		/* Keyword name, null-terminated.           */
  int  count;			/* Number of occurrences of keyword.        */
  int  idx[2];			/* Indices into fitskey array.              */
};

/* Size of the fitskeyid struct in int units, used by the Fortran wrappers. */
#define KEYIDLEN (sizeof(struct fitskeyid)/sizeof(int))


/* Struct used for storing FITS keywords. */
struct fitskey {
  int  keyno;			/* Header keyrecord sequence number (1-rel).*/
  int  keyid;			/* Index into fitskeyid[].                  */
  int  status;			/* Header keyrecord status bit flags.       */
  char keyword[12];		/* Keyword name, null-filled.               */
  int  type;			/* Keyvalue type (see above).               */
  int  padding;			/* (Dummy inserted for alignment purposes.) */
  union {
    int    i;			/* 32-bit integer and logical values.       */
    int64  k;			/* 64-bit integer values.                   */
    int    l[8];		/* Very long signed integer values.         */
    double f;			/* Floating point values.                   */
    double c[2];		/* Complex values.                          */
    char   s[72];		/* String values, null-terminated.          */
  } keyvalue;			/* Keyvalue.                                */
  int  ulen;			/* Length of units string.                  */
  char comment[84];		/* Comment (or keyrecord), null-terminated. */
};

/* Size of the fitskey struct in int units, used by the Fortran wrappers. */
#define KEYLEN (sizeof(struct fitskey)/sizeof(int))


int fitshdr(const char header[], int nkeyrec, int nkeyids,
            struct fitskeyid keyids[], int *nreject, struct fitskey **keys);


#ifdef __cplusplus
}
#endif

#endif /* WCSLIB_FITSHDR */
