/*============================================================================

  WCSLIB 5.4 - an implementation of the FITS WCS standard.
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
  $Id: wcshdr.h,v 5.4.1.2 2015/04/23 11:03:11 mcalabre Exp mcalabre $
*=============================================================================
*
* WCSLIB 5.4 - C routines that implement the FITS World Coordinate System
* (WCS) standard.  Refer to the README file provided with WCSLIB for an
* overview of the library.
*
*
* Summary of the wcshdr routines
* ------------------------------
* Routines in this suite are aimed at extracting WCS information from a FITS
* file.  The information is encoded via keywords defined in
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
*
* These routines provide the high-level interface between the FITS file and
* the WCS coordinate transformation routines.
*
* Additionally, function wcshdo() is provided to write out the contents of a
* wcsprm struct as a FITS header.
*
* Briefly, the anticipated sequence of operations is as follows:
*
*   - 1: Open the FITS file and read the image or binary table header, e.g.
*        using CFITSIO routine fits_hdr2str().
*
*   - 2: Parse the header using wcspih() or wcsbth(); they will automatically
*        interpret 'TAB' header keywords using wcstab().
*
*   - 3: Allocate memory for, and read 'TAB' arrays from the binary table
*        extension, e.g. using CFITSIO routine fits_read_wcstab() - refer to
*        the prologue of getwcstab.h.  wcsset() will automatically take
*        control of this allocated memory, in particular causing it to be
*        free'd by wcsfree().
*
*   - 4: Translate non-standard WCS usage using wcsfix(), see wcsfix.h.
*
*   - 5: Initialize wcsprm struct(s) using wcsset() and calculate coordinates
*        using wcsp2s() and/or wcss2p().  Refer to the prologue of wcs.h for a
*        description of these and other high-level WCS coordinate
*        transformation routines.
*
*   - 6: Clean up by freeing memory with wcsvfree().
*
* In detail:
*
* - wcspih() is a high-level FITS WCS routine that parses an image header.  It
*   returns an array of up to 27 wcsprm structs on each of which it invokes
*   wcstab().
*
* - wcsbth() is the analogue of wcspih() for use with binary tables; it
*   handles image array and pixel list keywords.  As an extension of the FITS
*   WCS standard, it also recognizes image header keywords which may be used
*   to provide default values via an inheritance mechanism.
*
* - wcstab() assists in filling in members of the wcsprm struct associated
*   with coordinate lookup tables ('TAB').  These are based on arrays stored
*   in a FITS binary table extension (BINTABLE) that are located by PVi_ma
*   keywords in the image header.
*
* - wcsidx() and wcsbdx() are utility routines that return the index for a
*   specified alternate coordinate descriptor in the array of wcsprm structs
*   returned by wcspih() or wcsbth().
*
* - wcsvfree() deallocates memory for an array of wcsprm structs, such as
*   returned by wcspih() or wcsbth().
*
* - wcshdo() writes out a wcsprm struct as a FITS header.
*
*
* wcspih() - FITS WCS parser routine for image headers
* ----------------------------------------------------
* wcspih() is a high-level FITS WCS routine that parses an image header,
* either that of a primary HDU or of an image extension.  All WCS keywords
* defined in Papers I, II, and III are recognized, and also those used by the
* AIPS convention and certain other keywords that existed in early drafts of
* the WCS papers as explained in wcsbth() note 5.
*
* Given a character array containing a FITS image header, wcspih() identifies
* and reads all WCS keywords for the primary coordinate representation and up
* to 26 alternate representations.  It returns this information as an array of
* wcsprm structs.
*
* wcspih() invokes wcstab() on each of the wcsprm structs that it returns.
*
* Use wcsbth() in preference to wcspih() for FITS headers of unknown type;
* wcsbth() can parse image headers as well as binary table and pixel list
* headers.
*
* Given and returned:
*   header    char[]    Character array containing the (entire) FITS image
*                       header from which to identify and construct the
*                       coordinate representations, for example, as might be
*                       obtained conveniently via the CFITSIO routine
*                       fits_hdr2str().
*
*                       Each header "keyrecord" (formerly "card image")
*                       consists of exactly 80 7-bit ASCII printing characters
*                       in the range 0x20 to 0x7e (which excludes NUL, BS,
*                       TAB, LF, FF and CR) especially noting that the
*                       keyrecords are NOT null-terminated.
*
*                       For negative values of ctrl (see below), header[] is
*                       modified so that WCS keyrecords processed by wcspih()
*                       are removed from it.
*
* Given:
*   nkeyrec   int       Number of keyrecords in header[].
*
*   relax     int       Degree of permissiveness:
*                         0: Recognize only FITS keywords defined by the
*                            published WCS standard.
*                         WCSHDR_all: Admit all recognized informal
*                            extensions of the WCS standard.
*                       Fine-grained control of the degree of permissiveness
*                       is also possible as explained in wcsbth() note 5.
*
*   ctrl      int       Error reporting and other control options for invalid
*                       WCS and other header keyrecords:
*                           0: Do not report any rejected header keyrecords.
*                           1: Produce a one-line message stating the number
*                              of WCS keyrecords rejected (nreject).
*                           2: Report each rejected keyrecord and the reason
*                              why it was rejected.
*                           3: As above, but also report all non-WCS
*                              keyrecords that were discarded, and the number
*                              of coordinate representations (nwcs) found.
*                           4: As above, but also report the accepted WCS
*                              keyrecords, with a summary of the number
*                              accepted as well as rejected.
*                       The report is written to stderr by default, or the
*                       stream set by wcsprintf_set().
*
*                       For ctrl < 0, WCS keyrecords processed by wcspih()
*                       are removed from header[]:
*                          -1: Remove only valid WCS keyrecords whose values
*                              were successfully extracted, nothing is
*                              reported.
*                          -2: As above, but also remove WCS keyrecords that
*                              were rejected, reporting each one and the
*                              reason that it was rejected.
*                          -3: As above, and also report the number of
*                              coordinate representations (nwcs) found.
*                         -11: Same as -1 but preserving the basic keywords
*                              '{DATE,MJD}-{OBS,AVG}' and 'OBSGEO-{X,Y,Z}'.
*                       If any keyrecords are removed from header[] it will
*                       be null-terminated (NUL not being a legal FITS header
*                       character), otherwise it will contain its original
*                       complement of nkeyrec keyrecords and possibly not be
*                       null-terminated.
*
* Returned:
*   nreject   int*      Number of WCS keywords rejected for syntax errors,
*                       illegal values, etc.  Keywords not recognized as WCS
*                       keywords are simply ignored.  Refer also to wcsbth()
*                       note 5.
*
*   nwcs      int*      Number of coordinate representations found.
*
*   wcs       struct wcsprm**
*                       Pointer to an array of wcsprm structs containing up to
*                       27 coordinate representations.
*
*                       Memory for the array is allocated by wcspih() which
*                       also invokes wcsini() for each struct to allocate
*                       memory for internal arrays and initialize their
*                       members to default values.  Refer also to wcsbth()
*                       note 8.  Note that wcsset() is not invoked on these
*                       structs.
*
*                       This allocated memory must be freed by the user, first
*                       by invoking wcsfree() for each struct, and then by
*                       freeing the array itself.  A routine, wcsvfree(), is
*                       provided to do this (see below).
*
* Function return value:
*             int       Status return value:
*                         0: Success.
*                         1: Null wcsprm pointer passed.
*                         2: Memory allocation failed.
*                         4: Fatal error returned by Flex parser.
*
* Notes:
*   Refer to wcsbth() notes 1, 2, 3, 5, 7, and 8.
*
*
* wcsbth() - FITS WCS parser routine for binary table and image headers
* ---------------------------------------------------------------------
* wcsbth() is a high-level FITS WCS routine that parses a binary table header.
* It handles image array and pixel list WCS keywords which may be present
* together in one header.
*
* As an extension of the FITS WCS standard, wcsbth() also recognizes image
* header keywords in a binary table header.  These may be used to provide
* default values via an inheritance mechanism discussed in note 5 (c.f.
* WCSHDR_AUXIMG and WCSHDR_ALLIMG), or may instead result in wcsprm structs
* that are not associated with any particular column.  Thus wcsbth() can
* handle primary image and image extension headers in addition to binary table
* headers (it ignores NAXIS and does not rely on the presence of the TFIELDS
* keyword).
*
* All WCS keywords defined in Papers I, II, and III are recognized, and also
* those used by the AIPS convention and certain other keywords that existed in
* early drafts of the WCS papers as explained in note 5 below.
*
* wcsbth() sets the colnum or colax[] members of the wcsprm structs that it
* returns with the column number of an image array or the column numbers
* associated with each pixel coordinate element in a pixel list.  wcsprm
* structs that are not associated with any particular column, as may be
* derived from image header keywords, have colnum == 0.
*
* Note 6 below discusses the number of wcsprm structs returned by wcsbth(),
* and the circumstances in which image header keywords cause a struct to be
* created.  See also note 9 concerning the number of separate images that may
* be stored in a pixel list.
*
* The API to wcsbth() is similar to that of wcspih() except for the addition
* of extra arguments that may be used to restrict its operation.  Like
* wcspih(), wcsbth() invokes wcstab() on each of the wcsprm structs that it
* returns.
*
* Given and returned:
*   header    char[]    Character array containing the (entire) FITS binary
*                       table, primary image, or image extension header from
*                       which to identify and construct the coordinate
*                       representations, for example, as might be obtained
*                       conveniently via the CFITSIO routine fits_hdr2str().
*
*                       Each header "keyrecord" (formerly "card image")
*                       consists of exactly 80 7-bit ASCII printing
*                       characters in the range 0x20 to 0x7e (which excludes
*                       NUL, BS, TAB, LF, FF and CR) especially noting that
*                       the keyrecords are NOT null-terminated.
*
*                       For negative values of ctrl (see below), header[] is
*                       modified so that WCS keyrecords processed by wcsbth()
*                       are removed from it.
*
* Given:
*   nkeyrec   int       Number of keyrecords in header[].
*
*   relax     int       Degree of permissiveness:
*                         0: Recognize only FITS keywords defined by the
*                            published WCS standard.
*                         WCSHDR_all: Admit all recognized informal
*                            extensions of the WCS standard.
*                       Fine-grained control of the degree of permissiveness
*                       is also possible, as explained in note 5 below.
*
*   ctrl      int       Error reporting and other control options for invalid
*                       WCS and other header keyrecords:
*                           0: Do not report any rejected header keyrecords.
*                           1: Produce a one-line message stating the number
*                              of WCS keyrecords rejected (nreject).
*                           2: Report each rejected keyrecord and the reason
*                              why it was rejected.
*                           3: As above, but also report all non-WCS
*                              keyrecords that were discarded, and the number
*                              of coordinate representations (nwcs) found.
*                           4: As above, but also report the accepted WCS
*                              keyrecords, with a summary of the number
*                              accepted as well as rejected.
*                       The report is written to stderr by default, or the
*                       stream set by wcsprintf_set().
*
*                       For ctrl < 0, WCS keyrecords processed by wcsbth()
*                       are removed from header[]:
*                          -1: Remove only valid WCS keyrecords whose values
*                              were successfully extracted, nothing is
*                              reported.
*                          -2: Also remove WCS keyrecords that were rejected,
*                              reporting each one and the reason that it was
*                              rejected.
*                          -3: As above, and also report the number of
*                              coordinate representations (nwcs) found.
*                         -11: Same as -1 but preserving the basic keywords
*                              '{DATE,MJD}-{OBS,AVG}' and 'OBSGEO-{X,Y,Z}'.
*                       If any keyrecords are removed from header[] it will
*                       be null-terminated (NUL not being a legal FITS header
*                       character), otherwise it will contain its original
*                       complement of nkeyrec keyrecords and possibly not be
*                       null-terminated.
*
*   keysel    int       Vector of flag bits that may be used to restrict the
*                       keyword types considered:
*                         WCSHDR_IMGHEAD: Image header keywords.
*                         WCSHDR_BIMGARR: Binary table image array.
*                         WCSHDR_PIXLIST: Pixel list keywords.
*                       If zero, there is no restriction.
*
*                       Keywords such as EQUIna or RFRQna that are common to
*                       binary table image arrays and pixel lists (including
*                       WCSNna and TWCSna, as explained in note 4 below) are
*                       selected by both WCSHDR_BIMGARR and WCSHDR_PIXLIST.
*                       Thus if inheritance via WCSHDR_ALLIMG is enabled as
*                       discussed in note 5 and one of these shared keywords
*                       is present, then WCSHDR_IMGHEAD and WCSHDR_PIXLIST
*                       alone may be sufficient to cause the construction of
*                       coordinate descriptions for binary table image arrays.
*
*   colsel    int*      Pointer to an array of table column numbers used to
*                       restrict the keywords considered by wcsbth().
*
*                       A null pointer may be specified to indicate that there
*                       is no restriction.  Otherwise, the magnitude of
*                       cols[0] specifies the length of the array:
*                         cols[0] > 0: the columns are included,
*                         cols[0] < 0: the columns are excluded.
*
*                       For the pixel list keywords TPn_ka and TCn_ka (and
*                       TPCn_ka and TCDn_ka if WCSHDR_LONGKEY is enabled), it
*                       is an error for one column to be selected but not the
*                       other.  This is unlike the situation with invalid
*                       keyrecords, which are simply rejected, because the
*                       error is not intrinsic to the header itself but
*                       arises in the way that it is processed.
*
* Returned:
*   nreject   int*      Number of WCS keywords rejected for syntax errors,
*                       illegal values, etc.  Keywords not recognized as WCS
*                       keywords are simply ignored, refer also to note 5
*                       below.
*
*   nwcs      int*      Number of coordinate representations found.
*
*   wcs       struct wcsprm**
*                       Pointer to an array of wcsprm structs containing up
*                       to 27027 coordinate representations, refer to note 6
*                       below.
*
*                       Memory for the array is allocated by wcsbth() which
*                       also invokes wcsini() for each struct to allocate
*                       memory for internal arrays and initialize their
*                       members to default values.  Refer also to note 8
*                       below.  Note that wcsset() is not invoked on these
*                       structs.
*
*                       This allocated memory must be freed by the user, first
*                       by invoking wcsfree() for each struct, and then by
*                       freeing the array itself.  A routine, wcsvfree(), is
*                       provided to do this (see below).
*
* Function return value:
*             int       Status return value:
*                         0: Success.
*                         1: Null wcsprm pointer passed.
*                         2: Memory allocation failed.
*                         3: Invalid column selection.
*                         4: Fatal error returned by Flex parser.
*
* Notes:
*   1: wcspih() determines the number of coordinate axes independently for
*      each alternate coordinate representation (denoted by the "a" value in
*      keywords like CTYPEia) from the higher of
*
*        a: NAXIS,
*        b: WCSAXESa,
*        c: The highest axis number in any parameterized WCS keyword.  The
*           keyvalue, as well as the keyword, must be syntactically valid
*           otherwise it will not be considered.
*
*      If none of these keyword types is present, i.e. if the header only
*      contains auxiliary WCS keywords for a particular coordinate
*      representation, then no coordinate description is constructed for it.
*
*      wcsbth() is similar except that it ignores the NAXIS keyword if given
*      an image header to process.
*
*      The number of axes, which is returned as a member of the wcsprm
*      struct, may differ for different coordinate representations of the
*      same image.
*
*   2: wcspih() and wcsbth() enforce correct FITS "keyword = value" syntax
*      with regard to "= " occurring in columns 9 and 10.
*
*      However, they do recognize free-format character (NOST 100-2.0,
*      Sect. 5.2.1), integer (Sect. 5.2.3), and floating-point values
*      (Sect. 5.2.4) for all keywords.
*
*   3: Where CROTAn, CDi_ja, and PCi_ja occur together in one header wcspih()
*      and wcsbth() treat them as described in the prologue to wcs.h.
*
*   4: WCS Paper I mistakenly defined the pixel list form of WCSNAMEa as
*      TWCSna instead of WCSNna; the 'T' is meant to substitute for the axis
*      number in the binary table form of the keyword - note that keywords
*      defined in WCS Papers II and III that are not parameterized by axis
*      number have identical forms for binary tables and pixel lists.
*      Consequently wcsbth() always treats WCSNna and TWCSna as equivalent.
*
*   5: wcspih() and wcsbth() interpret the "relax" argument as a vector of
*      flag bits to provide fine-grained control over what non-standard WCS
*      keywords to accept.  The flag bits are subject to change in future and
*      should be set by using the preprocessor macros (see below) for the
*      purpose.
*
*      - WCSHDR_none: Don't accept any extensions (not even those in the
*              errata).  Treat non-conformant keywords in the same way as
*              non-WCS keywords in the header, i.e. simply ignore them.
*
*      - WCSHDR_all: Accept all extensions recognized by the parser.
*
*      - WCSHDR_reject: Reject non-standard keyrecords (that are not otherwise
*              explicitly accepted by one of the flags below).  A message will
*              optionally be printed on stderr by default, or the stream set
*              by wcsprintf_set(), as determined by the ctrl argument, and
*              nreject will be incremented.
*
*              This flag may be used to signal the presence of non-standard
*              keywords, otherwise they are simply passed over as though they
*              did not exist in the header.  It is mainly intended for testing
*              conformance of a FITS header to the WCS standard.
*
*              Keyrecords may be non-standard in several ways:
*
*                - The keyword may be syntactically valid but with keyvalue of
*                  incorrect type or invalid syntax, or the keycomment may be
*                  malformed.
*
*                - The keyword may strongly resemble a WCS keyword but not, in
*                  fact, be one because it does not conform to the standard.
*                  For example, "CRPIX01" looks like a CRPIXja keyword, but in
*                  fact the leading zero on the axis number violates the basic
*                  FITS standard.  Likewise, "LONPOLE2" is not a valid
*                  LONPOLEa keyword in the WCS standard, and indeed there is
*                  nothing the parser can sensibly do with it.
*
*                - Use of the keyword may be deprecated by the standard.  Such
*                  will be rejected if not explicitly accepted via one of the
*                  flags below.
*
*      - WCSHDR_strict: As for WCSHDR_reject, but also reject AIPS-convention
*              keywords and all other deprecated usage that is not explicitly
*              accepted.
*
*      - WCSHDR_CROTAia: Accept CROTAia (wcspih()),
*                               iCROTna (wcsbth()),
*                               TCROTna (wcsbth()).
*      - WCSHDR_EPOCHa:  Accept EPOCHa.
*      - WCSHDR_VELREFa: Accept VELREFa.
*              wcspih() always recognizes the AIPS-convention keywords,
*              CROTAn, EPOCH, and VELREF for the primary representation
*              (a = ' ') but alternates are non-standard.
*
*              wcsbth() accepts EPOCHa and VELREFa only if WCSHDR_AUXIMG is
*              also enabled.
*
*      - WCSHDR_CD00i00j: Accept CD00i00j (wcspih()).
*      - WCSHDR_PC00i00j: Accept PC00i00j (wcspih()).
*      - WCSHDR_PROJPn:   Accept PROJPn   (wcspih()).
*              These appeared in early drafts of WCS Paper I+II (before they
*              were split) and are equivalent to CDi_ja, PCi_ja, and PVi_ma
*              for the primary representation (a = ' ').  PROJPn is
*              equivalent to PVi_ma with m = n <= 9, and is associated
*              exclusively with the latitude axis.
*
*      - WCSHDR_CD0i_0ja: Accept CD0i_0ja (wcspih()).
*      - WCSHDR_PC0i_0ja: Accept PC0i_0ja (wcspih()).
*      - WCSHDR_PV0i_0ma: Accept PV0i_0ja (wcspih()).
*      - WCSHDR_PS0i_0ma: Accept PS0i_0ja (wcspih()).
*              Allow the numerical index to have a leading zero in doubly-
*              parameterized keywords, for example, PC01_01.  WCS Paper I
*              (Sects 2.1.2 & 2.1.4) explicitly disallows leading zeroes.
*              The FITS 3.0 standard document (Sect. 4.1.2.1) states that the
*              index in singly-parameterized keywords (e.g. CTYPEia) "shall
*              not have leading zeroes", and later in Sect. 8.1 that "leading
*              zeroes must not be used" on PVi_ma and PSi_ma.  However, by an
*              oversight, it is silent on PCi_ja and CDi_ja.
*
*      - WCSHDR_RADECSYS: Accept RADECSYS.  This appeared in early drafts of
*              WCS Paper I+II and was subsequently replaced by RADESYSa.
*
*              wcsbth() accepts RADECSYS only if WCSHDR_AUXIMG is also
*              enabled.
*
*      - WCSHDR_VSOURCE: Accept VSOURCEa or VSOUna (wcsbth()).  This appeared
*              in early drafts of WCS Paper III and was subsequently dropped
*              in favour of ZSOURCEa and ZSOUna.
*
*              wcsbth() accepts VSOURCEa only if WCSHDR_AUXIMG is also
*              enabled.
*
*      - WCSHDR_DOBSn (wcsbth() only): Allow DOBSn, the column-specific
*              analogue of DATE-OBS.  By an oversight this was never formally
*              defined in the standard.
*
*      - WCSHDR_LONGKEY (wcsbth() only): Accept long forms of the alternate
*              binary table and pixel list WCS keywords, i.e. with "a" non-
*              blank.  Specifically
*
#                jCRPXna  TCRPXna  :  jCRPXn  jCRPna  TCRPXn  TCRPna  CRPIXja
#                   -     TPCn_ka  :    -     ijPCna    -     TPn_ka  PCi_ja
#                   -     TCDn_ka  :    -     ijCDna    -     TCn_ka  CDi_ja
#                iCDLTna  TCDLTna  :  iCDLTn  iCDEna  TCDLTn  TCDEna  CDELTia
#                iCUNIna  TCUNIna  :  iCUNIn  iCUNna  TCUNIn  TCUNna  CUNITia
#                iCTYPna  TCTYPna  :  iCTYPn  iCTYna  TCTYPn  TCTYna  CTYPEia
#                iCRVLna  TCRVLna  :  iCRVLn  iCRVna  TCRVLn  TCRVna  CRVALia
#                iPVn_ma  TPVn_ma  :    -     iVn_ma    -     TVn_ma  PVi_ma
#                iPSn_ma  TPSn_ma  :    -     iSn_ma    -     TSn_ma  PSi_ma
*
*              where the primary and standard alternate forms together with
*              the image-header equivalent are shown rightwards of the colon.
*
*              The long form of these keywords could be described as quasi-
*              standard.  TPCn_ka, iPVn_ma, and TPVn_ma appeared by mistake
*              in the examples in WCS Paper II and subsequently these and
*              also TCDn_ka, iPSn_ma and TPSn_ma were legitimized by the
*              errata to the WCS papers.
*
*              Strictly speaking, the other long forms are non-standard and
*              in fact have never appeared in any draft of the WCS papers nor
*              in the errata.  However, as natural extensions of the primary
*              form they are unlikely to be written with any other intention.
*              Thus it should be safe to accept them provided, of course,
*              that the resulting keyword does not exceed the 8-character
*              limit.
*
*              If WCSHDR_CNAMn is enabled then also accept
*
#                iCNAMna  TCNAMna  :   ---   iCNAna    ---   TCNAna  CNAMEia
#                iCRDEna  TCRDEna  :   ---   iCRDna    ---   TCRDna  CRDERia
#                iCSYEna  TCSYEna  :   ---   iCSYna    ---   TCSYna  CSYERia
*
*              Note that CNAMEia, CRDERia, CSYERia, and their variants are
*              not used by WCSLIB but are stored in the wcsprm struct as
*              auxiliary information.
*
*      - WCSHDR_CNAMn (wcsbth() only): Accept iCNAMn, iCRDEn, iCSYEn, TCNAMn,
*              TCRDEn, and TCSYEn, i.e. with "a" blank.  While non-standard,
*              these are the obvious analogues of iCTYPn, TCTYPn, etc.
*
*      - WCSHDR_AUXIMG (wcsbth() only): Allow the image-header form of an
*              auxiliary WCS keyword with representation-wide scope to
*              provide a default value for all images.  This default may be
*              overridden by the column-specific form of the keyword.
*
*              For example, a keyword like EQUINOXa would apply to all image
*              arrays in a binary table, or all pixel list columns with
*              alternate representation "a" unless overridden by EQUIna.
*
*              Specifically the keywords are:
*
#                LATPOLEa  for LATPna
#                LONPOLEa  for LONPna
#                RESTFREQ  for RFRQna
#                RESTFRQa  for RFRQna
#                RESTWAVa  for RWAVna
*
*              whose keyvalues are actually used by WCSLIB, and also keywords
*              that provide auxiliary information that is simply stored in
*              the wcsprm struct:
*
#                EPOCH         -       ... (No column-specific form.)
#                EPOCHa        -       ... Only if WCSHDR_EPOCHa is set.
#                EQUINOXa  for EQUIna
#                RADESYSa  for RADEna
#                RADECSYS  for RADEna  ... Only if WCSHDR_RADECSYS is set.
#                SPECSYSa  for SPECna
#                SSYSOBSa  for SOBSna
#                SSYSSRCa  for SSRCna
#                VELOSYSa  for VSYSna
#                VELANGLa  for VANGna
#                VELREF        -       ... (No column-specific form.)
#                VELREFa       -       ... Only if WCSHDR_VELREFa is set.
#                VSOURCEa  for VSOUna  ... Only if WCSHDR_VSOURCE is set.
#                WCSNAMEa  for WCSNna  ... Or TWCSna (see below).
#                ZSOURCEa  for ZSOUna
*
#                DATE-AVG  for DAVGn
#                DATE-OBS  for DOBSn
#                MJD-AVG   for MJDAn
#                MJD-OBS   for MJDOBn
#                OBSGEO-X  for OBSGXn
#                OBSGEO-Y  for OBSGYn
#                OBSGEO-Z  for OBSGZn
*
*              where the image-header keywords on the left provide default
*              values for the column specific keywords on the right.
*
*              Keywords in the last group, such as MJD-OBS, apply to all
*              alternate representations, so MJD-OBS would provide a default
*              value for all images in the header.
*
*              This auxiliary inheritance mechanism applies to binary table
*              image arrays and pixel lists alike.  Most of these keywords
*              have no default value, the exceptions being LONPOLEa and
*              LATPOLEa, and also RADESYSa and EQUINOXa which provide
*              defaults for each other.  Thus the only potential difficulty
*              in using WCSHDR_AUXIMG is that of erroneously inheriting one
*              of these four keywords.
*
*              Unlike WCSHDR_ALLIMG, the existence of one (or all) of these
*              auxiliary WCS image header keywords will not by itself cause a
*              wcsprm struct to be created for alternate representation "a".
*              This is because they do not provide sufficient information to
*              create a non-trivial coordinate representation when used in
*              conjunction with the default values of those keywords, such as
*              CTYPEia, that are parameterized by axis number.
*
*      - WCSHDR_ALLIMG (wcsbth() only): Allow the image-header form of *all*
*              image header WCS keywords to provide a default value for all
*              image arrays in a binary table (n.b. not pixel list).  This
*              default may be overridden by the column-specific form of the
*              keyword.
*
*              For example, a keyword like CRPIXja would apply to all image
*              arrays in a binary table with alternate representation "a"
*              unless overridden by jCRPna.
*
*              Specifically the keywords are those listed above for
*              WCSHDR_AUXIMG plus
*
#                WCSAXESa  for WCAXna
*
*              which defines the coordinate dimensionality, and the following
*              keywords which are parameterized by axis number:
*
#                CRPIXja   for jCRPna
#                PCi_ja    for ijPCna
#                CDi_ja    for ijCDna
#                CDELTia   for iCDEna
#                CROTAi    for iCROTn
#                CROTAia        -      ... Only if WCSHDR_CROTAia is set.
#                CUNITia   for iCUNna
#                CTYPEia   for iCTYna
#                CRVALia   for iCRVna
#                PVi_ma    for iVn_ma
#                PSi_ma    for iSn_ma
*
#                CNAMEia   for iCNAna
#                CRDERia   for iCRDna
#                CSYERia   for iCSYna
*
*              where the image-header keywords on the left provide default
*              values for the column specific keywords on the right.
*
*              This full inheritance mechanism only applies to binary table
*              image arrays, not pixel lists, because in the latter case
*              there is no well-defined association between coordinate axis
*              number and column number.
*
*              Note that CNAMEia, CRDERia, CSYERia, and their variants are
*              not used by WCSLIB but are stored in the wcsprm struct as
*              auxiliary information.
*
*              Note especially that at least one wcsprm struct will be
*              returned for each "a" found in one of the image header
*              keywords listed above:
*
*              - If the image header keywords for "a" ARE NOT inherited by a
*                binary table, then the struct will not be associated with
*                any particular table column number and it is up to the user
*                to provide an association.
*
*              - If the image header keywords for "a" ARE inherited by a
*                binary table image array, then those keywords are considered
*                to be "exhausted" and do not result in a separate wcsprm
*                struct.
*
*      For example, to accept CD00i00j and PC00i00j and reject all other
*      extensions, use
*
=        relax = WCSHDR_reject | WCSHDR_CD00i00j | WCSHDR_PC00i00j;
*
*      The parser always treats EPOCH as subordinate to EQUINOXa if both are
*      present, and VSOURCEa is always subordinate to ZSOURCEa.
*
*      Likewise, VELREF is subordinate to the formalism of WCS Paper III, see
*      spcaips().
*
*      Neither wcspih() nor wcsbth() currently recognize the AIPS-convention
*      keywords ALTRPIX or ALTRVAL which effectively define an alternative
*      representation for a spectral axis.
*
*   6: Depending on what flags have been set in its "relax" argument,
*      wcsbth() could return as many as 27027 wcsprm structs:
*
*      - Up to 27 unattached representations derived from image header
*        keywords.
*
*      - Up to 27 structs for each of up to 999 columns containing an image
*        arrays.
*
*      - Up to 27 structs for a pixel list.
*
*      Note that it is considered legitimate for a column to contain an image
*      array and also form part of a pixel list, and in particular that
*      wcsbth() does not check the TFORM keyword for a pixel list column to
*      check that it is scalar.
*
*      In practice, of course, a realistic binary table header is unlikely to
*      contain more than a handful of images.
*
*      In order for wcsbth() to create a wcsprm struct for a particular
*      coordinate representation, at least one WCS keyword that defines an
*      axis number must be present, either directly or by inheritance if
*      WCSHDR_ALLIMG is set.
*
*      When the image header keywords for an alternate representation are
*      inherited by a binary table image array via WCSHDR_ALLIMG, those
*      keywords are considered to be "exhausted" and do not result in a
*      separate wcsprm struct.  Otherwise they do.
*
*   7: Neither wcspih() nor wcsbth() check for duplicated keywords, in most
*      cases they accept the last encountered.
*
*   8: wcspih() and wcsbth() use wcsnpv() and wcsnps() (refer to the prologue
*      of wcs.h) to match the size of the pv[] and ps[] arrays in the wcsprm
*      structs to the number in the header.  Consequently there are no unused
*      elements in the pv[] and ps[] arrays, indeed they will often be of
*      zero length.
*
*   9: The FITS WCS standard for pixel lists assumes that a pixel list
*      defines one and only one image, i.e. that each row of the binary table
*      refers to just one event, e.g. the detection of a single photon or
*      neutrino.
*
*      In the absence of a formal mechanism for identifying the columns
*      containing pixel coordinates (as opposed to pixel values or ancillary
*      data recorded at the time the photon or neutrino was detected),
*      Paper I discusses how the WCS keywords themselves may be used to
*      identify them.
*
*      In practice, however, pixel lists have been used to store multiple
*      images.  Besides not specifying how to identify columns, the pixel
*      list convention is also silent on the method to be used to associate
*      table columns with image axes.
*
*      wcsbth() simply collects all WCS keywords for a particular coordinate
*      representation (i.e. the "a" value in TCTYna) into one wcsprm struct.
*      However, these alternates need not be associated with the same table
*      columns and this allows a pixel list to contain up to 27 separate
*      images.  As usual, if one of these representations happened to contain
*      more than two celestial axes, for example, then an error would result
*      when wcsset() is invoked on it.  In this case the "colsel" argument
*      could be used to restrict the columns used to construct the
*      representation so that it only contained one pair of celestial axes.
*
*
* wcstab() - Tabular construction routine
* ---------------------------------------
* wcstab() assists in filling in the information in the wcsprm struct relating
* to coordinate lookup tables.
*
* Tabular coordinates ('TAB') present certain difficulties in that the main
* components of the lookup table - the multidimensional coordinate array plus
* an index vector for each dimension - are stored in a FITS binary table
* extension (BINTABLE).  Information required to locate these arrays is stored
* in PVi_ma and PSi_ma keywords in the image header.
*
* wcstab() parses the PVi_ma and PSi_ma keywords associated with each 'TAB'
* axis and allocates memory in the wcsprm struct for the required number of
* tabprm structs.  It sets as much of the tabprm struct as can be gleaned from
* the image header, and also sets up an array of wtbarr structs (described in
* the prologue of wcs.h) to assist in extracting the required arrays from the
* BINTABLE extension(s).
*
* It is then up to the user to allocate memory for, and copy arrays from the
* BINTABLE extension(s) into the tabprm structs.  A CFITSIO routine,
* fits_read_wcstab(), has been provided for this purpose, see getwcstab.h.
* wcsset() will automatically take control of this allocated memory, in
* particular causing it to be free'd by wcsfree(); the user must not attempt
* to free it after wcsset() has been called.
*
* Note that wcspih() and wcsbth() automatically invoke wcstab() on each of the
* wcsprm structs that they return.
*
* Given and returned:
*   wcs       struct wcsprm*
*                       Coordinate transformation parameters (see below).
*
*                       wcstab() sets ntab, tab, nwtb and wtb, allocating
*                       memory for the tab and wtb arrays.  This allocated
*                       memory will be free'd automatically by wcsfree().
*
* Function return value:
*             int       Status return value:
*                         0: Success.
*                         1: Null wcsprm pointer passed.
*                         2: Memory allocation failed.
*                         3: Invalid tabular parameters.
*
*                       For returns > 1, a detailed error message is set in
*                       wcsprm::err if enabled, see wcserr_enable().
*
*
* wcsidx() - Index alternate coordinate representations
* -----------------------------------------------------
* wcsidx() returns an array of 27 indices for the alternate coordinate
* representations in the array of wcsprm structs returned by wcspih().  For
* the array returned by wcsbth() it returns indices for the unattached
* (colnum == 0) representations derived from image header keywords - use
* wcsbdx() for those derived from binary table image arrays or pixel lists
* keywords.
*
* Given:
*   nwcs      int       Number of coordinate representations in the array.
*
*   wcs       const struct wcsprm**
*                       Pointer to an array of wcsprm structs returned by
*                       wcspih() or wcsbth().
*
* Returned:
*   alts      int[27]   Index of each alternate coordinate representation in
*                       the array: alts[0] for the primary, alts[1] for 'A',
*                       etc., set to -1 if not present.
*
*                       For example, if there was no 'P' representation then
*
=                         alts['P'-'A'+1] == -1;
*
*                       Otherwise, the address of its wcsprm struct would be
*
=                         wcs + alts['P'-'A'+1];
*
* Function return value:
*             int       Status return value:
*                         0: Success.
*                         1: Null wcsprm pointer passed.
*
*
* wcsbdx() - Index alternate coordinate representions
* ---------------------------------------------------
* wcsbdx() returns an array of 999 x 27 indices for the alternate coordinate
* representions for binary table image arrays xor pixel lists in the array of
* wcsprm structs returned by wcsbth().  Use wcsidx() for the unattached
* representations derived from image header keywords.
*
* Given:
*   nwcs      int       Number of coordinate representations in the array.
*
*   wcs       const struct wcsprm**
*                       Pointer to an array of wcsprm structs returned by
*                       wcsbth().
*
*   type      int       Select the type of coordinate representation:
*                         0: binary table image arrays,
*                         1: pixel lists.
*
* Returned:
*   alts      short[1000][28]
*                       Index of each alternate coordinate represention in the
*                       array: alts[col][0] for the primary, alts[col][1] for
*                       'A', to alts[col][26] for 'Z', where col is the
*                       1-relative column number, and col == 0 is used for
*                       unattached image headers.  Set to -1 if not present.
*
*                       alts[col][27] counts the number of coordinate
*                       representations of the chosen type for each column.
*
*                       For example, if there was no 'P' represention for
*                       column 13 then
*
=                         alts[13]['P'-'A'+1] == -1;
*
*                       Otherwise, the address of its wcsprm struct would be
*
=                         wcs + alts[13]['P'-'A'+1];
*
* Function return value:
*             int       Status return value:
*                         0: Success.
*                         1: Null wcsprm pointer passed.
*
*
* wcsvfree() - Free the array of wcsprm structs
* ---------------------------------------------
* wcsvfree() frees the memory allocated by wcspih() or wcsbth() for the array
* of wcsprm structs, first invoking wcsfree() on each of the array members.
*
* Given and returned:
*   nwcs      int*      Number of coordinate representations found; set to 0
*                       on return.
*
*   wcs       struct wcsprm**
*                       Pointer to the array of wcsprm structs; set to 0x0 on
*                       return.
*
* Function return value:
*             int       Status return value:
*                         0: Success.
*                         1: Null wcsprm pointer passed.
*
*
* wcshdo() - Write out a wcsprm struct as a FITS header
* -----------------------------------------------------
* wcshdo() translates a wcsprm struct into a FITS header.  If the colnum
* member of the struct is non-zero then a binary table image array header will
* be produced.  Otherwise, if the colax[] member of the struct is set non-zero
* then a pixel list header will be produced.  Otherwise, a primary image or
* image extension header will be produced.
*
* If the struct was originally constructed from a header, e.g. by wcspih(),
* the output header will almost certainly differ in a number of respects:
*
*   - The output header only contains WCS-related keywords.  In particular, it
*     does not contain syntactically-required keywords such as SIMPLE, NAXIS,
*     BITPIX, or END.
*
*   - Deprecated (e.g. CROTAn) or non-standard usage will be translated to
*     standard (this is partially dependent on whether wcsfix() was applied).
*
*   - Quantities will be converted to the units used internally, basically SI
*     with the addition of degrees.
*
*   - Floating-point quantities may be given to a different decimal precision.
*
*   - Elements of the PCi_ja matrix will be written if and only if they differ
*     from the unit matrix.  Thus, if the matrix is unity then no elements
*     will be written.
*
*   - Additional keywords such as WCSAXESa, CUNITia, LONPOLEa and LATPOLEa may
*     appear.
*
*   - The original keycomments will be lost, although wcshdo() tries hard to
*     write meaningful comments.
*
*   - Keyword order may be changed.
*
* Keywords can be translated between the image array, binary table, and pixel
* lists forms by manipulating the colnum or colax[] members of the wcsprm
* struct.
*
* Given:
*   relax     int       Degree of permissiveness:
*                          0: Recognize only FITS keywords defined by the
*                             published WCS standard.
*                         -1: Admit all informal extensions of the WCS
*                             standard.
*                       Fine-grained control of the degree of permissiveness
*                       is also possible as explained in the notes below.
*
* Given and returned:
*   wcs       struct wcsprm*
*                       Pointer to a wcsprm struct containing coordinate
*                       transformation parameters.  Will be initialized if
*                       necessary.
*
* Returned:
*   nkeyrec   int*      Number of FITS header keyrecords returned in the
*                       "header" array.
*
*   header    char**    Pointer to an array of char holding the header.
*                       Storage for the array is allocated by wcshdo() in
*                       blocks of 2880 bytes (32 x 80-character keyrecords)
*                       and must be free'd by the user to avoid memory leaks.
*
*                       Each keyrecord is 80 characters long and is *NOT*
*                       null-terminated, so the first keyrecord starts at
*                       (*header)[0], the second at (*header)[80], etc.
*
* Function return value:
*             int       Status return value (associated with wcs_errmsg[]):
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
*   wcshdo() interprets the "relax" argument as a vector of flag bits to
*   provide fine-grained control over what non-standard WCS keywords to write.
*   The flag bits are subject to change in future and should be set by using
*   the preprocessor macros (see below) for the purpose.
*
*   - WCSHDO_none: Don't use any extensions.
*
*   - WCSHDO_all: Write all recognized extensions, equivalent to setting each
*           flag bit.
*
*   - WCSHDO_safe: Write all extensions that are considered to be safe and
*           recommended.
*
*   - WCSHDO_DOBSn: Write DOBSn, the column-specific analogue of DATE-OBS for
*           use in binary tables and pixel lists.  WCS Paper III introduced
*           DATE-AVG and DAVGn but by an oversight DOBSn (the obvious analogy)
*           was never formally defined by the standard.  The alternative to
*           using DOBSn is to write DATE-OBS which applies to the whole table.
*           This usage is considered to be safe and is recommended.
*
*   - WCSHDO_TPCn_ka: WCS Paper I defined
*
*           - TPn_ka and TCn_ka for pixel lists
*
*           but WCS Paper II uses TPCn_ka in one example and subsequently the
*           errata for the WCS papers legitimized the use of
*
*           - TPCn_ka and TCDn_ka for pixel lists
*
*           provided that the keyword does not exceed eight characters.  This
*           usage is considered to be safe and is recommended because of the
*           non-mnemonic terseness of the shorter forms.
*
*   - WCSHDO_PVn_ma: WCS Paper I defined
*
*           - iVn_ma and iSn_ma for bintables and
*           - TVn_ma and TSn_ma for pixel lists
*
*           but WCS Paper II uses iPVn_ma and TPVn_ma in the examples and
*           subsequently the errata for the WCS papers legitimized the use of
*
*           - iPVn_ma and iPSn_ma for bintables and
*           - TPVn_ma and TPSn_ma for pixel lists
*
*           provided that the keyword does not exceed eight characters.  This
*           usage is considered to be safe and is recommended because of the
*           non-mnemonic terseness of the shorter forms.
*
*   - WCSHDO_CRPXna: For historical reasons WCS Paper I defined
*
*           - jCRPXn, iCDLTn, iCUNIn, iCTYPn, and iCRVLn for bintables and
*           - TCRPXn, TCDLTn, TCUNIn, TCTYPn, and TCRVLn for pixel lists
*
*           for use without an alternate version specifier.  However, because
*           of the eight-character keyword constraint, in order to accommodate
*           column numbers greater than 99 WCS Paper I also defined
*
*           - jCRPna, iCDEna, iCUNna, iCTYna and iCRVna for bintables and
*           - TCRPna, TCDEna, TCUNna, TCTYna and TCRVna for pixel lists
*
*           for use with an alternate version specifier (the "a").  Like the
*           PC, CD, PV, and PS keywords there is an obvious tendency to
*           confuse these two forms for column numbers up to 99.  It is very
*           unlikely that any parser would reject keywords in the first set
*           with a non-blank alternate version specifier so this usage is
*           considered to be safe and is recommended.
*
*   - WCSHDO_CNAMna: WCS Papers I and III defined
*
*           - iCNAna,  iCRDna,  and iCSYna  for bintables and
*           - TCNAna,  TCRDna,  and TCSYna  for pixel lists
*
*           By analogy with the above, the long forms would be
*
*           - iCNAMna, iCRDEna, and iCSYEna for bintables and
*           - TCNAMna, TCRDEna, and TCSYEna for pixel lists
*
*           Note that these keywords provide auxiliary information only, none
*           of them are needed to compute world coordinates.  This usage is
*           potentially unsafe and is not recommended at this time.
*
*   - WCSHDO_WCSNna: In light of wcsbth() note 4, write WCSNna instead of
*           TWCSna for pixel lists.  While wcsbth() treats WCSNna and TWCSna
*           as equivalent, other parsers may not.  Consequently, this usage
*           is potentially unsafe and is not recommended at this time.
*
*
* Global variable: const char *wcshdr_errmsg[] - Status return messages
* ---------------------------------------------------------------------
* Error messages to match the status value returned from each function.
* Use wcs_errmsg[] for status returns from wcshdo().
*
*===========================================================================*/

#ifndef WCSLIB_WCSHDR
#define WCSLIB_WCSHDR

#include "wcs.h"

#ifdef __cplusplus
extern "C" {
#endif

#define WCSHDR_none     0x00000000
#define WCSHDR_all      0x000FFFFF
#define WCSHDR_reject   0x10000000
#define WCSHDR_strict   0x20000000

#define WCSHDR_CROTAia  0x00000001
#define WCSHDR_EPOCHa   0x00000002
#define WCSHDR_VELREFa  0x00000004
#define WCSHDR_CD00i00j 0x00000008
#define WCSHDR_PC00i00j 0x00000010
#define WCSHDR_PROJPn   0x00000020
#define WCSHDR_CD0i_0ja 0x00000040
#define WCSHDR_PC0i_0ja 0x00000080
#define WCSHDR_PV0i_0ma 0x00000100
#define WCSHDR_PS0i_0ma 0x00000200
#define WCSHDR_RADECSYS 0x00000400
#define WCSHDR_VSOURCE  0x00000800
#define WCSHDR_DOBSn    0x00001000
#define WCSHDR_LONGKEY  0x00002000
#define WCSHDR_CNAMn    0x00004000
#define WCSHDR_AUXIMG   0x00008000
#define WCSHDR_ALLIMG   0x00010000

#define WCSHDR_IMGHEAD  0x00100000
#define WCSHDR_BIMGARR  0x00200000
#define WCSHDR_PIXLIST  0x00400000

#define WCSHDO_none     0x00
#define WCSHDO_all      0xFF
#define WCSHDO_safe     0x0F
#define WCSHDO_DOBSn    0x01
#define WCSHDO_TPCn_ka  0x02
#define WCSHDO_PVn_ma   0x04
#define WCSHDO_CRPXna   0x08
#define WCSHDO_CNAMna   0x10
#define WCSHDO_WCSNna   0x20


extern const char *wcshdr_errmsg[];

enum wcshdr_errmsg_enum {
  WCSHDRERR_SUCCESS            = 0,	/* Success. */
  WCSHDRERR_NULL_POINTER       = 1,	/* Null wcsprm pointer passed. */
  WCSHDRERR_MEMORY             = 2,	/* Memory allocation failed. */
  WCSHDRERR_BAD_COLUMN         = 3,	/* Invalid column selection. */
  WCSHDRERR_PARSER             = 4,	/* Fatal error returned by Flex
					   parser. */
  WCSHDRERR_BAD_TABULAR_PARAMS = 5 	/* Invalid tabular parameters. */
};

int wcspih(char *header, int nkeyrec, int relax, int ctrl, int *nreject,
           int *nwcs, struct wcsprm **wcs);

int wcsbth(char *header, int nkeyrec, int relax, int ctrl, int keysel,
           int *colsel, int *nreject, int *nwcs, struct wcsprm **wcs);

int wcstab(struct wcsprm *wcs);

int wcsidx(int nwcs, struct wcsprm **wcs, int alts[27]);

int wcsbdx(int nwcs, struct wcsprm **wcs, int type, short alts[1000][28]);

int wcsvfree(int *nwcs, struct wcsprm **wcs);

int wcshdo(int relax, struct wcsprm *wcs, int *nkeyrec, char **header);


#ifdef __cplusplus
}
#endif

#endif /* WCSLIB_WCSHDR */
