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
  $Id: wcsbth.l,v 4.25 2014/12/14 14:29:36 mcalabre Exp $
*=============================================================================
*
* wcsbth.l is a Flex description file containing the definition of a lexical
* scanner for parsing the WCS keyrecords for one or more image arrays and/or
* pixel lists in a FITS binary table header.  It can also handle primary image
* and image extension headers.
*
* wcsbth.l requires Flex v2.5.4 or later.  Refer to wcshdr.h for a description
* of the user interface and operating notes.
*
* Implementation notes
* --------------------
* wcsbth() may be invoked with an option that causes it to recognise the
* image-header form of WCS keywords as defaults for each alternate coordinate
* representation (up to 27).  By design, with this option enabled wcsbth() can
* also handle primary image and image extension headers, effectively treating
* them as a single-column binary table though with WCS keywords of a different
* form.
*
* NAXIS is always 2 for binary tables, it refers to the two-dimensional nature
* of the table.  Thus NAXIS does not count the number of image axes in either
* image arrays or pixels lists and for the latter there is not even a formal
* equivalent of WCSAXESa.  Hence NAXIS is always ignored and a first pass
* through the header is required to determine the number of images, the number
* of alternate coordinate representations for each image (up to 27), and the
* number of coordinate axes in each representation; this pass also counts the
* number of iPVn_ma and iPSn_ma or TVk_ma and TSk_ma keywords in each
* representation.
*
* On completion of the first pass, the association between column number and
* axis number is defined for each representation of a pixel list.  Memory is
* allocated for an array of the required number of wcsprm structs and each of
* these is initialized appropriately.  These structs are filled in the second
* pass.
*
* It is permissible for a scalar table column to contain degenerate (single-
* point) image arrays and simultaneously form one axis of a pixel list.
*
* The parser does not check for duplicated keywords, for most keywords it
* accepts the last encountered.
*
* wcsbth() does not currently handle the Green Bank convention.
*
*===========================================================================*/

/* Options. */
%option full
%option never-interactive
%option noyywrap
%option outfile="wcsbth.c"
%option prefix="wcsbth"

/* Indices for parameterized keywords. */
I0	[0-9]
I1	[1-9]
I2	[1-9][0-9]
I3	[1-9][0-9]{2}
I4	[1-9][0-9]{3}

/* Alternate coordinate system identifier. */
ALT	[ A-Z]

/* Keyvalue data types. */
INTEGER	[+-]?[0-9]+
FLOAT	[+-]?([0-9]+\.?[0-9]*|\.[0-9]+)([eE][+-]?[0-9]+)?
STRING	'([^']|'')*'

/* Exclusive start states. */
%x CCCCCia   iCCCna iCCCCn    TCCCna TCCCCn
%x CCi_ja    ijCCna           TCn_ka TCCn_ka
%x CROTAi           iCROTn           TCROTn
%x CCi_ma    iCn_ma iCCn_ma   TCn_ma TCCn_ma
%x PROJPm
%x CCCCCCCC CCCCCCCa
%x CCCCna   CCCCCna
%x CCCCn    CCCCCn
%x VALUE INTEGER_VAL FLOAT_VAL STRING_VAL
%x COMMENT DISCARD ERROR FLUSH

%{
#include <math.h>
#include <setjmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "wcs.h"
#include "wcshdr.h"
#include "wcsmath.h"
#include "wcsprintf.h"
#include "wcsutil.h"

			/* Codes used for keyvalue data types. */
#define INTEGER 0
#define FLOAT   1
#define STRING  2

			/* Bit masks used for keyword types:        */
#define IMGAUX  0x1	/* Auxiliary image header, e.g. LONPOLEa or */
			/* DATE-OBS.                                */
#define IMGAXIS 0x2	/* Image header with axis number, e.g.      */
			/* CTYPEia.                                 */
#define IMGHEAD 0x3	/* Image header of either type.             */
#define BIMGARR 0x4	/* Binary table image array with axis       */
			/* number, e.g. iCTYna.                     */
#define PIXLIST 0x8	/* Pixel list, e.g. TCTYna.                 */
#define BINTAB  0xC	/* Shared binary table image array (without */
			/* axis number) or pixel list, e.g. LONPna  */
			/* or OBSGXn.                               */

#define YY_DECL int wcsbth(char *header, int nkeyrec, int relax, int ctrl, \
                           int keysel, int *colsel, int *nreject, int *nwcs, \
			   struct wcsprm **wcs)

#define YY_INPUT(inbuff, count, bufsize) \
	{ \
	  if (wcsbth_nkeyrec) { \
	    strncpy(inbuff, wcsbth_hdr, 80); \
	    inbuff[80] = '\n'; \
	    wcsbth_hdr += 80; \
	    wcsbth_nkeyrec--; \
	    count = 81; \
	  } else { \
	    count = YY_NULL; \
	  } \
	}

/* A convenience macro to get around incompatibilities between unput() and
   yyless(): put yytext followed by a blank back onto the input stream. */
#define WCSBTH_PUTBACK \
  sprintf(stmp, "%s ", yytext); \
  itmp = strlen(stmp); \
  while (itmp) unput(stmp[--itmp]);

/* These global variables are required by YY_INPUT. */
char *wcsbth_hdr;
int  wcsbth_nkeyrec;

/* Used in preempting the call to exit() by yy_fatal_error(). */
jmp_buf wcsbth_abort_jmp_env;
#define exit(status) longjmp(wcsbth_abort_jmp_env, status)

/* Struct used internally for header bookkeeping. */
struct wcsbth_alts {
  int ncol, ialt, icol, imgherit;
  short int (*arridx)[27];
  short int pixidx[27];
  short int pad1;
  unsigned int *pixlist;

  unsigned char (*npv)[27];
  unsigned char (*nps)[27];
  unsigned char pixnpv[27];
  unsigned char pixnps[27];
  unsigned char pad2[2];
};

int wcsbth_pass1(int keytype, int i, int j, int n, int k, char a, char ptype,
        struct wcsbth_alts *alts);
int wcsbth_init1(struct wcsbth_alts *alts, int *nwcs, struct wcsprm **wcs);

struct wcsprm *wcsbth_idx(struct wcsprm *wcs, struct wcsbth_alts *alts,
        int keytype, int n, char a);
int wcsbth_colax(struct wcsprm *wcs, struct wcsbth_alts *alts, int k, char a);

int wcsbth_epoch(void *wptr);
int wcsbth_vsource(void *wptr);

int wcsbth_final(struct wcsbth_alts *alts, int *nwcs, struct wcsprm **wcs);

%}

%%
	/* Keyword indices, as used in the WCS papers, e.g. iVn_ma, TPn_ka. */
	char a;
	int  i, j, k, m, n;
	
	char *cptr, *errmsg, errtxt[80], exclude[1000], *extkey, *hptr, ptype,
	     stmp[16];
	int  altlin, ialt, icol, incl, ipass, ipx, itmp, ix, jx, keytype,
	     nsel, npass, status, valtype, voff;
	void   *vptr, *wptr;
	struct wcsbth_alts alts;
	struct wcsprm *wcsp, wcstem;
	int (*special)(void *);
	int yylex_destroy(void);
	
	/* The data structures produced. */
	*nwcs = 0;
	*wcs  = 0x0;
	
	/* Parameters used to implement YY_INPUT. */
	wcsbth_hdr = header;
	wcsbth_nkeyrec = nkeyrec;
	
	/* Our handle on the input stream. */
	hptr = header;
	*nreject = 0;
	
	/* Keyword parameters. */
	i = j = 0;
	n = k = 0;
	m = 0;
	a = ' ';
	
	/* Header bookkeeping. */
	alts.ncol = 0;
	alts.arridx  = 0x0;
	alts.pixlist = 0x0;
	alts.npv = 0x0;
	alts.nps = 0x0;
	
	for (ialt = 0; ialt < 27; ialt++) {
	  alts.pixidx[ialt] = 0;
	  alts.pixnpv[ialt] = 0;
	  alts.pixnps[ialt] = 0;
	}
	
	/* For decoding the keyvalue. */
	keytype =  0;
	valtype = -1;
	vptr    = 0x0;
	
	/* For keywords that require special handling. */
	altlin = 0;
	ptype  = ' ';
	special = 0x0;
	
	/* Selection by column number. */
	nsel = colsel ? colsel[0] : 0;
	incl = (nsel > 0);
	for (icol = 0; icol < 1000; icol++) {
	  exclude[icol] = incl;
	}
	for (icol = 1; icol <= abs(nsel); icol++) {
	  itmp = colsel[icol];
	  if (0 < itmp && itmp < 1000) {
	    exclude[itmp] = !incl;
	  }
	}
	exclude[0] = 0;
	
	/* Selection by keyword type. */
	itmp = keysel;
	keysel = 0;
	if (itmp) {
	  if (itmp & WCSHDR_IMGHEAD) keysel |= IMGHEAD;
	  if (itmp & WCSHDR_BIMGARR) keysel |= BIMGARR;
	  if (itmp & WCSHDR_PIXLIST) keysel |= PIXLIST;
	}
	if (keysel == 0) {
	  keysel = IMGHEAD | BINTAB;
	}
	
	/* Control variables. */
	ipass = 1;
	npass = 2;
	
	/* Return here via longjmp() invoked by yy_fatal_error(). */
	if (setjmp(wcsbth_abort_jmp_env)) {
	  return 4;
	}
	
	BEGIN(INITIAL);


^TFIELDS" = "" "*{INTEGER} {
	  if (ipass == 1) {
	    if (alts.ncol == 0) {
	      sscanf(yytext, "TFIELDS = %d", &(alts.ncol));
	      BEGIN(FLUSH);
	    } else {
	      errmsg = "Duplicate or out-of-sequence TFIELDS keyword";
	      BEGIN(ERROR);
	    }
	
	  } else {
	    BEGIN(FLUSH);
	  }
	}

^WCSAXES{ALT}=" "" "*{INTEGER} {
	  keytype = IMGAXIS;
	
	  if (!(keytype & keysel)) {
	    /* Ignore this key type. */
	    BEGIN(DISCARD);
	
	  } else {
	    if (relax & WCSHDR_ALLIMG) {
	      if (ipass == 1) {
	        sscanf(yytext, "WCSAXES%c= %d", &a, &i);
	        wcsbth_pass1(IMGAXIS, i, 0, 0, 0, a, ' ', &alts);
	      }
	
	      BEGIN(FLUSH);
	
	    } else if (relax & WCSHDR_reject) {
	      errmsg = "Image-header keyword WCSAXESa in binary table";
	      BEGIN(ERROR);
	
	    } else {
	      /* Pretend we don't recognize it. */
	      BEGIN(DISCARD);
	    }
	  }
	}

^WCAX{I1}{ALT}"  = "" "*{INTEGER} |
^WCAX{I2}{ALT}" = "" "*{INTEGER}  |
^WCAX{I3}{ALT}"= "" "*{INTEGER} {
	  keytype = BIMGARR;
	
	  /* Note that a blank in the sscanf() format string matches zero or
	     more of them in the input. */
	  sscanf(yytext, "WCAX%d%c = %d", &n, &a, &i);
	
	  if (!(keytype & keysel) || exclude[n]) {
	    /* Ignore this key type or column. */
	    BEGIN(DISCARD);
	  } else {
	    if (ipass == 1) {
	      wcsbth_pass1(BIMGARR, i, 0, n, 0, a, ' ', &alts);
	    }
	    BEGIN(FLUSH);
	  }
	}

^WCST{I1}{ALT}"  = "" "*{STRING} |
^WCST{I2}{ALT}" = "" "*{STRING} |
^WCST{I3}{ALT}"= "" "*{STRING} {
	  /* Cross-reference supplier. */
	  keytype = BIMGARR;
	  errmsg = "Cross-references are not currently implemented";
	  BEGIN(ERROR);
	}

^WCSX{I1}{ALT}"  = "" "*{STRING} |
^WCSX{I2}{ALT}" = "" "*{STRING} |
^WCSX{I3}{ALT}"= "" "*{STRING} {
	  /* Cross-reference consumer. */
	  keytype = BIMGARR;
	  errmsg = "Cross-references are not currently implemented";
	  BEGIN(ERROR);
	}

^CRPIX	{
	  valtype = FLOAT;
	  vptr = &(wcstem.crpix);
	
	  extkey = "CRPIXja";
	  BEGIN(CCCCCia);
	}

^{I1}CRP  |
^{I1}CRPX {
	  valtype = FLOAT;
	  vptr = &(wcstem.crpix);
	
	  sscanf(yytext, "%d", &i);
	
	  if (yyleng == 4) {
	    BEGIN(iCCCna);
	  } else {
	    extkey = "jCRPXn";
	    BEGIN(iCCCCn);
	  }
	}

^TCRP	|
^TCRPX	{
	  valtype = FLOAT;
	  vptr = &(wcstem.crpix);
	
	  if (yyleng == 4) {
	    BEGIN(TCCCna);
	  } else {
	    extkey = "TCRPXn";
	    BEGIN(TCCCCn);
	  }
	}

^PC	{
	  valtype = FLOAT;
	  vptr = &(wcstem.pc);
	  altlin = 1;
	
	  extkey = "PCi_ja";
	  BEGIN(CCi_ja);
	}

^{I2}PC	{
	  valtype = FLOAT;
	  vptr = &(wcstem.pc);
	  altlin = 1;
	
	  sscanf(yytext, "%1d%1d", &i, &j);
	
	  BEGIN(ijCCna);
	}

^TP	|
^TPC	{
	  valtype = FLOAT;
	  vptr = &(wcstem.pc);
	  altlin = 1;
	
	  if (yyleng == 2) {
	    BEGIN(TCn_ka);
	  } else {
	    extkey = "TPCn_ka";
	    BEGIN(TCCn_ka);
	  }
	}

^CD	{
	  valtype = FLOAT;
	  vptr = &(wcstem.cd);
	  altlin = 2;
	
	  extkey = "CDi_ja";
	  BEGIN(CCi_ja);
	}

^{I2}CD	{
	  valtype = FLOAT;
	  vptr = &(wcstem.cd);
	  altlin = 2;
	
	  sscanf(yytext, "%1d%1d", &i, &j);
	
	  BEGIN(ijCCna);
	}

^TC	|
^TCD	{
	  valtype = FLOAT;
	  vptr = &(wcstem.cd);
	  altlin = 2;
	
	  if (yyleng == 2) {
	    BEGIN(TCn_ka);
	  } else {
	    extkey = "TCDn_ka";
	    BEGIN(TCCn_ka);
	  }
	}

^CDELT	{
	  valtype = FLOAT;
	  vptr = &(wcstem.cdelt);
	
	  extkey = "CDELTia";
	  BEGIN(CCCCCia);
	}

^{I1}CDE  |
^{I1}CDLT {
	  valtype = FLOAT;
	  vptr = &(wcstem.cdelt);
	
	  sscanf(yytext, "%d", &i);
	
	  if (yyleng == 4) {
	    BEGIN(iCCCna);
	  } else {
	    extkey = "iCDLTn";
	    BEGIN(iCCCCn);
	  }
	}

^TCDE	|
^TCDLT	{
	  valtype = FLOAT;
	  vptr = &(wcstem.cdelt);
	
	  if (yyleng == 4) {
	    BEGIN(TCCCna);
	  } else {
	    extkey = "TCDLTn";
	    BEGIN(TCCCCn);
	  }
	}

^CROTA	{
	  valtype = FLOAT;
	  vptr = &(wcstem.crota);
	  altlin = 4;
	
	  extkey = "CROTAi";
	  BEGIN(CROTAi);
	}

^{I1}CROT {
	  valtype = FLOAT;
	  vptr = &(wcstem.crota);
	  altlin = 4;
	
	  sscanf(yytext, "%d", &i);
	
	  extkey = "iCROTn";
	  BEGIN(iCROTn);
	}

^TCROT	{
	  valtype = FLOAT;
	  vptr = &(wcstem.crota);
	  altlin = 4;
	
	  extkey = "TCROTn";
	  BEGIN(TCROTn);
	}

^CUNIT	{
	  valtype = STRING;
	  vptr = &(wcstem.cunit);
	
	  extkey = "CUNITia";
	  BEGIN(CCCCCia);
	}

^{I1}CUN  |
^{I1}CUNI {
	  valtype = STRING;
	  vptr = &(wcstem.cunit);
	
	  sscanf(yytext, "%d", &i);
	
	  if (yyleng == 4) {
	    BEGIN(iCCCna);
	  } else {
	    extkey = "iCUNIn";
	    BEGIN(iCCCCn);
	  }
	}

^TCUN	|
^TCUNI	{
	  valtype = STRING;
	  vptr = &(wcstem.cunit);
	
	  if (yyleng == 4) {
	    BEGIN(TCCCna);
	  } else {
	    extkey = "TCUNIn";
	    BEGIN(TCCCCn);
	  }
	}

^CTYPE	{
	  valtype = STRING;
	  vptr = &(wcstem.ctype);
	
	  extkey = "CTYPEia";
	  BEGIN(CCCCCia);
	}

^{I1}CTY  |
^{I1}CTYP {
	  valtype = STRING;
	  vptr = &(wcstem.ctype);
	
	  sscanf(yytext, "%d", &i);
	
	  if (yyleng == 4) {
	    BEGIN(iCCCna);
	  } else {
	    extkey = "iCTYPn";
	    BEGIN(iCCCCn);
	  }
	}

^TCTY	|
^TCTYP	{
	  valtype = STRING;
	  vptr = &(wcstem.ctype);
	
	  if (yyleng == 4) {
	    BEGIN(TCCCna);
	  } else {
	    extkey = "TCTYPn";
	    BEGIN(TCCCCn);
	  }
	}

^CRVAL	{
	  valtype = FLOAT;
	  vptr = &(wcstem.crval);
	
	  extkey = "CRVALia";
	  BEGIN(CCCCCia);
	}

^{I1}CRV  |
^{I1}CRVL {
	  valtype = FLOAT;
	  vptr = &(wcstem.crval);
	
	  sscanf(yytext, "%d", &i);
	
	  if (yyleng == 4) {
	    BEGIN(iCCCna);
	  } else {
	    extkey = "iCRVLn";
	    BEGIN(iCCCCn);
	  }
	}

^TCRV	|
^TCRVL	{
	  valtype = FLOAT;
	  vptr = &(wcstem.crval);
	
	  if (yyleng == 4) {
	    BEGIN(TCCCna);
	  } else {
	    extkey = "TCRVLn";
	    BEGIN(TCCCCn);
	  }
	}

^LONPOLE |
^LONP	{
	  valtype = FLOAT;
	  vptr = &(wcstem.lonpole);
	
	  if (yyleng == 7) {
	    extkey = "LONPOLEa";
	    BEGIN(CCCCCCCa);
	  } else {
	    BEGIN(CCCCna);
	  }
	}

^LATPOLE |
^LATP	{
	  valtype = FLOAT;
	  vptr = &(wcstem.latpole);
	
	  if (yyleng == 7) {
	    extkey = "LATPOLEa";
	    BEGIN(CCCCCCCa);
	  } else {
	    BEGIN(CCCCna);
	  }
	}

^RESTFRQ  |
^RESTFREQ |
^RFRQ	{
	  valtype = FLOAT;
	  vptr = &(wcstem.restfrq);
	
	  if (yyleng == 8) {
	    unput(' ');
	    extkey = "RESTFREQ";
	    BEGIN(CCCCCCCa);
	  } else if (yyleng == 7) {
	    extkey = "RESTFRQa";
	    BEGIN(CCCCCCCa);
	  } else {
	    BEGIN(CCCCna);
	  }
	}

^RESTWAV |
^RWAV	{
	  valtype = FLOAT;
	  vptr = &(wcstem.restwav);
	
	  if (yyleng == 7) {
	    extkey = "RESTWAVa";
	    BEGIN(CCCCCCCa);
	  } else {
	    BEGIN(CCCCna);
	  }
	}

^PV	{
	  valtype = FLOAT;
	  vptr = &(wcstem.pv);
	  ptype = 'v';
	
	  extkey = "PVi_ma";
	  BEGIN(CCi_ma);
	}

^{I1}V	|
^{I1}PV	{
	  valtype = FLOAT;
	  vptr = &(wcstem.pv);
	  ptype = 'v';
	
	  sscanf(yytext, "%d", &i);
	
	  if (yyleng == 2) {
	    BEGIN(iCn_ma);
	  } else {
	    extkey = "iPVn_ma";
	    BEGIN(iCCn_ma);
	  }
	}

^TV	|
^TPV	{
	  valtype = FLOAT;
	  vptr = &(wcstem.pv);
	  ptype = 'v';
	
	  if (yyleng == 2) {
	    BEGIN(TCn_ma);
	  } else {
	    extkey = "TPVn_ma";
	    BEGIN(TCCn_ma);
	  }
	}

^PROJP	{
	  valtype = FLOAT;
	  vptr = &(wcstem.pv);
	  ptype = 'v';
	
	  BEGIN(PROJPm);
	}

^PS	{
	  valtype = STRING;
	  vptr = &(wcstem.ps);
	  ptype = 's';
	
	  extkey = "PSi_ma";
	  BEGIN(CCi_ma);
	}

^{I1}S	|
^{I1}PS	{
	  valtype = STRING;
	  vptr = &(wcstem.ps);
	  ptype = 's';
	
	  sscanf(yytext, "%d", &i);
	
	  if (yyleng == 2) {
	    BEGIN(iCn_ma);
	  } else {
	    extkey = "iPSn_ma";
	    BEGIN(iCCn_ma);
	  }
	}

^TS	|
^TPS	{
	  valtype = STRING;
	  vptr = &(wcstem.ps);
	  ptype = 's';
	
	  if (yyleng == 2) {
	    BEGIN(TCn_ma);
	  } else {
	    extkey = "TPSn_ma";
	    BEGIN(TCCn_ma);
	  }
	}

^CNAME	{
	  valtype = STRING;
	  vptr = &(wcstem.cname);
	
	  extkey = "CNAMEia";
	  BEGIN(CCCCCia);
	}

^{I1}CNA  |
^{I1}CNAM {
	  valtype = STRING;
	  vptr = &(wcstem.cname);
	
	  sscanf(yytext, "%d", &i);
	
	  if (yyleng == 4) {
	    BEGIN(iCCCna);
	  } else {
	    if (!(relax & WCSHDR_CNAMn)) vptr = 0x0;
	    extkey = "iCNAMn";
	    BEGIN(iCCCCn);
	  }
	}

^TCNA	|
^TCNAM	{
	  valtype = STRING;
	  vptr = &(wcstem.cname);
	
	  if (yyleng == 4) {
	    BEGIN(TCCCna);
	  } else {
	    if (!(relax & WCSHDR_CNAMn)) vptr = 0x0;
	    extkey = "TCNAMn";
	    BEGIN(TCCCCn);
	  }
	}

^CRDER	{
	  valtype = FLOAT;
	  vptr = &(wcstem.crder);
	
	  extkey = "CRDERia";
	  BEGIN(CCCCCia);
	}

^{I1}CRD |
^{I1}CRDE {
	  valtype = FLOAT;
	  vptr = &(wcstem.crder);
	
	  sscanf(yytext, "%d", &i);
	
	  if (yyleng == 4) {
	    BEGIN(iCCCna);
	  } else {
	    if (!(relax & WCSHDR_CNAMn)) vptr = 0x0;
	    extkey = "iCRDEn";
	    BEGIN(iCCCCn);
	  }
	}

^TCRD	|
^TCRDE	{
	  valtype = FLOAT;
	  vptr = &(wcstem.crder);
	
	  if (yyleng == 4) {
	    BEGIN(TCCCna);
	  } else {
	    if (!(relax & WCSHDR_CNAMn)) vptr = 0x0;
	    extkey = "TCRDEn";
	    BEGIN(TCCCCn);
	  }
	}

^CSYER	{
	  valtype = FLOAT;
	  vptr = &(wcstem.csyer);
	
	  extkey = "CSYERia";
	  BEGIN(CCCCCia);
	}

^{I1}CSY  |
^{I1}CSYE {
	  valtype = FLOAT;
	  vptr = &(wcstem.csyer);
	
	  sscanf(yytext, "%d", &i);
	
	  if (yyleng == 4) {
	    BEGIN(iCCCna);
	  } else {
	    if (!(relax & WCSHDR_CNAMn)) vptr = 0x0;
	    extkey = "iCSYEn";
	    BEGIN(iCCCCn);
	  }
	}

^TCSY	|
^TCSYE	{
	  valtype = FLOAT;
	  vptr = &(wcstem.csyer);
	
	  if (yyleng == 4) {
	    BEGIN(TCCCna);
	  } else {
	    if (!(relax & WCSHDR_CNAMn)) vptr = 0x0;
	    extkey = "TCSYEn";
	    BEGIN(TCCCCn);
	  }
	}

^DATE-AVG |
^DAVG   {
	  valtype = STRING;
	  vptr = wcstem.dateavg;
	
	  if (yyleng == 8) {
	    extkey = "DATE-AVG";
	    BEGIN(CCCCCCCC);
	  } else {
	    BEGIN(CCCCn);
	  }
	}

^DATE-OBS {
	  valtype = STRING;
	  vptr = wcstem.dateobs;
	
	  extkey = "DATE-OBS";
	  BEGIN(CCCCCCCC);
	}

^DOBS{I1}"   " |
^DOBS{I2}"  "  |
^DOBS{I3}" " {
	  if (relax & WCSHDR_DOBSn) {
	    valtype = STRING;
	    vptr = wcstem.dateobs;
	
	    yyless(4);
	    BEGIN(CCCCn);
	
	  } else {
	    keytype = BINTAB;
	    if (relax & WCSHDR_reject) {
	      errmsg = "DOBSna keyword is non-standard";
	      BEGIN(ERROR);
	    } else {
	      BEGIN(DISCARD);
	    }
	  }
	}

^EPOCH{ALT}"  " {
	  sscanf(yytext, "EPOCH%c", &a);
	
	  if (a == ' ' || (relax & WCSHDR_EPOCHa)) {
	    valtype = FLOAT;
	    vptr = &(wcstem.equinox);
	    special = wcsbth_epoch;
	
	    unput(a);
	    extkey = "EPOCH";
	    BEGIN(CCCCCCCa);
	
	  } else {
	    keytype = IMGAUX;
	    if (relax & WCSHDR_reject) {
	      errmsg = "EPOCH keyword may not have an alternate version code";
	      BEGIN(ERROR);
	    } else {
	      BEGIN(DISCARD);
	    }
	  }
	}

^EQUINOX |
^EQUI	{
	  valtype = FLOAT;
	  vptr = &(wcstem.equinox);
	
	  if (yyleng == 7) {
	    extkey = "EQUINOXa";
	    BEGIN(CCCCCCCa);
	  } else {
	    BEGIN(CCCCna);
	  }
	}

^MJD-AVG" " |
^MJDA	{
	  valtype = FLOAT;
	  vptr = &(wcstem.mjdavg);
	
	  if (yyleng == 8) {
	    extkey = "MJD-AVG";
	    BEGIN(CCCCCCCC);
	  } else {
	    BEGIN(CCCCn);
	  }
	}

^MJD-OBS" " |
^MJDOB	{
	  valtype = FLOAT;
	  vptr = &(wcstem.mjdobs);
	
	  if (yyleng == 8) {
	    extkey = "MJD-OBS";
	    BEGIN(CCCCCCCC);
	  } else {
	    BEGIN(CCCCCn);
	  }
	}

^OBSGEO-X |
^OBSGX	{
	  valtype = FLOAT;
	  vptr = wcstem.obsgeo;
	
	  if (yyleng == 8) {
	    extkey = "OBSGEO-X";
	    BEGIN(CCCCCCCC);
	  } else {
	    BEGIN(CCCCCn);
	  }
	}

^OBSGEO-Y |
^OBSGY	{
	  valtype = FLOAT;
	  vptr = wcstem.obsgeo + 1;
	
	  if (yyleng == 8) {
	    extkey = "OBSGEO-Y";
	    BEGIN(CCCCCCCC);
	  } else {
	    BEGIN(CCCCCn);
	  }
	}

^OBSGEO-Z |
^OBSGZ	{
	  valtype = FLOAT;
	  vptr = wcstem.obsgeo + 2;
	
	  if (yyleng == 8) {
	    extkey = "OBSGEO-Z";
	    BEGIN(CCCCCCCC);
	  } else {
	    BEGIN(CCCCCn);
	  }
	}

^RADESYS |
^RADE	{
	  valtype = STRING;
	  vptr = wcstem.radesys;
	
	  if (yyleng == 7) {
	    extkey = "RADESYSa";
	    BEGIN(CCCCCCCa);
	  } else {
	    BEGIN(CCCCna);
	  }
	}

^RADECSYS {
	  if (relax & WCSHDR_RADECSYS) {
	    valtype = STRING;
	    vptr = wcstem.radesys;
	
	    unput(' ');
	    extkey = "RADECSYS";
	    BEGIN(CCCCCCCa);
	
	  } else {
	    keytype = IMGAUX;
	    if (relax & WCSHDR_reject) {
	      errmsg = "RADECSYS keyword is non-standard";
	      BEGIN(ERROR);
	    } else {
	      BEGIN(DISCARD);
	    }
	  }
	}

^SPECSYS |
^SPEC	{
	  valtype = STRING;
	  vptr = wcstem.specsys;
	
	  if (yyleng == 7) {
	    extkey = "SPECSYSa";
	    BEGIN(CCCCCCCa);
	  } else {
	    BEGIN(CCCCna);
	  }
	}

^SSYSOBS |
^SOBS	{
	  valtype = STRING;
	  vptr = wcstem.ssysobs;
	
	  if (yyleng == 7) {
	    extkey = "SSYSOBSa";
	    BEGIN(CCCCCCCa);
	  } else {
	    BEGIN(CCCCna);
	  }
	}

^SSYSSRC |
^SSRC	{
	  valtype = STRING;
	  vptr = wcstem.ssyssrc;
	
	  if (yyleng == 7) {
	    extkey = "SSYSSRCa";
	    BEGIN(CCCCCCCa);
	  } else {
	    BEGIN(CCCCna);
	  }
	}

^VELOSYS |
^VSYS	{
	  valtype = FLOAT;
	  vptr = &(wcstem.velosys);
	
	  if (yyleng == 7) {
	    extkey = "VELOSYSa";
	    BEGIN(CCCCCCCa);
	  } else {
	    BEGIN(CCCCna);
	  }
	}

^VELANGL |
^VANG	{
	  valtype = FLOAT;
	  vptr = &(wcstem.velangl);
	
	  if (yyleng == 7) {
	    extkey = "VELANGLa";
	    BEGIN(CCCCCCCa);
	  } else {
	    BEGIN(CCCCna);
	  }
	}

^VELREF{ALT}" " {
	  sscanf(yytext, "VELREF%c", &a);
	
	  if (a == ' ' || (relax & WCSHDR_VELREFa)) {
	    valtype = INTEGER;
	    vptr = &(wcstem.velref);
	
	    unput(a);
	    extkey = "VELREF";
	    BEGIN(CCCCCCCa);
	
	  } else {
	    keytype = IMGAUX;
	    if (relax & WCSHDR_reject) {
	      errmsg = "VELREF keyword may not have an alternate version code";
	      BEGIN(ERROR);
	    } else {
	      BEGIN(DISCARD);
	    }
	  }
	}

^VSOURCE{ALT} {
	  if (relax & WCSHDR_VSOURCE) {
	    valtype = FLOAT;
	    vptr = &(wcstem.zsource);
	    special = wcsbth_vsource;
	
	    yyless(7);
	    extkey = "VSOURCEa";
	    BEGIN(CCCCCCCa);
	
	  } else {
	    keytype = IMGAUX;
	    if (relax & WCSHDR_reject) {
	      errmsg = "VSOURCEa keyword is deprecated";
	      BEGIN(ERROR);
	    } else {
	      BEGIN(DISCARD);
	    }
	  }
	}

^VSOU{I1}{ALT}"  " |
^VSOU{I2}{ALT}" "  |
^VSOU{I3}{ALT} {
	  if (relax & WCSHDR_VSOURCE) {
	    valtype = FLOAT;
	    vptr = &(wcstem.zsource);
	    special = wcsbth_vsource;
	
	    yyless(4);
	    BEGIN(CCCCna);
	
	  } else {
	    keytype = BINTAB;
	    if (relax & WCSHDR_reject) {
	      errmsg = "VSOUna keyword is deprecated";
	      BEGIN(ERROR);
	    } else {
	      /* Pretend we don't recognize it. */
	      BEGIN(DISCARD);
	    }
	  }
	}

^WCSNAME |
^WCSN	 |
^TWCS	{
	  valtype = STRING;
	  vptr = wcstem.wcsname;
	
	  if (yyleng == 7) {
	    extkey = "WCSNAMEa";
	    BEGIN(CCCCCCCa);
	  } else {
	    BEGIN(CCCCna);
	  }
	}

^ZSOURCE |
^ZSOU	{
	  valtype = FLOAT;
	  vptr = &(wcstem.zsource);
	
	  if (yyleng == 7) {
	    extkey = "ZSOURCEa";
	    BEGIN(CCCCCCCa);
	  } else {
	    BEGIN(CCCCna);
	  }
	}

^END" "{77} {
	  yyless(0);
	  if (wcsbth_nkeyrec) {
	    wcsbth_nkeyrec = 0;
	    errmsg = "Keyrecords following the END keyrecord were ignored";
	    BEGIN(ERROR);
	  } else {
	    BEGIN(DISCARD);
	  }
	}

^.	{
	  yyless(0);
	  BEGIN(DISCARD);
	}

<CCCCCia>{I1}{ALT}" " |
<CCCCCia>{I2}{ALT} {
	  /* Image-header keyword. */
	  keytype = IMGAXIS;
	  if (relax & WCSHDR_ALLIMG) {
	    sscanf(yytext, "%d%c", &i, &a);
	    BEGIN(VALUE);
	
	  } else if (relax & WCSHDR_reject) {
	    errmsg = errtxt;
	    sprintf(errmsg,
	      "Image-header keyword %s in binary table", extkey);
	    BEGIN(ERROR);
	
	  } else {
	    /* Pretend we don't recognize it. */
	    BEGIN(DISCARD);
	  }
	}

<CCCCCia>{I3} {
	  /* Invalid axis number in image-header keyword. */
	  keytype = IMGAXIS;
	  if (relax & WCSHDR_ALLIMG) {
	    /* Will also be flagged by <VALUE> as invalid. */
	    sscanf(yytext, "%3d", &i);
	    BEGIN(VALUE);
	
	  } else {
	    /* Pretend we don't recognize it. */
	    BEGIN(DISCARD);
	  }
	}

<iCCCCn>{I1}"  " |
<iCCCCn>{I2}" "  |
<iCCCCn>{I3}     |
<TCCCCn>{I1}"  " |
<TCCCCn>{I2}" "  |
<TCCCCn>{I3} {
	  if (vptr) {
	    WCSBTH_PUTBACK;
	    BEGIN((YY_START == iCCCCn) ? iCCCna : TCCCna);
	  } else {
	    keytype = (YY_START == iCCCCn) ? BIMGARR : PIXLIST;
	    if (relax & WCSHDR_reject) {
	      errmsg = errtxt;
	      sprintf(errmsg, "%s keyword is non-standard", extkey);
	      BEGIN(ERROR);
	    } else {
	      BEGIN(DISCARD);
	    }
	  }
	}

<iCCCCn>{I1}[A-Z]" " |
<iCCCCn>{I2}[A-Z]    |
<TCCCCn>{I1}[A-Z]" " |
<TCCCCn>{I2}[A-Z] {
	  if (vptr && (relax & WCSHDR_LONGKEY)) {
	    WCSBTH_PUTBACK;
	    BEGIN((YY_START == iCCCCn) ? iCCCna : TCCCna);
	
	  } else {
	    keytype = (YY_START == iCCCna) ? BIMGARR : PIXLIST;
	    if (relax & WCSHDR_reject) {
	      errmsg = errtxt;
	      if (!vptr) {
	        sprintf(errmsg, "%s keyword is non-standard", extkey);
	      } else {
	        sprintf(errmsg,
	          "%s keyword may not have an alternate version code", extkey);
	      }
	      BEGIN(ERROR);
	
	    } else {
	      /* Pretend we don't recognize it. */
	      BEGIN(DISCARD);
	    }
	  }
	}

<iCCCCn>. |
<TCCCCn>. {
	  BEGIN(DISCARD);
	}

<iCCCna>{I1}{ALT}"  " |
<iCCCna>{I2}{ALT}" "  |
<iCCCna>{I3}{ALT}     |
<TCCCna>{I1}{ALT}"  " |
<TCCCna>{I2}{ALT}" "  |
<TCCCna>{I3}{ALT} {
	  sscanf(yytext, "%d%c", &n, &a);
	  if (YY_START == TCCCna) i = wcsbth_colax(*wcs, &alts, n, a);
	  keytype = (YY_START == iCCCna) ? BIMGARR : PIXLIST;
	  BEGIN(VALUE);
	}

<iCCCna>. |
<TCCCna>. {
	  BEGIN(DISCARD);
	}

<CCi_ja>{I1}_{I1}{ALT}"  " |
<CCi_ja>{I1}_{I2}{ALT}" " |
<CCi_ja>{I2}_{I1}{ALT}" " |
<CCi_ja>{I2}_{I2}{ALT} {
	  /* Image-header keyword. */
	  if (relax & WCSHDR_ALLIMG) {
	    sscanf(yytext, "%d_%d%c", &i, &j, &a);
	    keytype = IMGAXIS;
	    BEGIN(VALUE);
	
	  } else if (relax & WCSHDR_reject) {
	    errmsg = errtxt;
	    sprintf(errmsg,
	      "Image-header keyword %s in binary table", extkey);
	    BEGIN(ERROR);
	
	  } else {
	    /* Pretend we don't recognize it. */
	    BEGIN(DISCARD);
	  }
	}

<CCi_ja>{I1}_{I3}{ALT} |
<CCi_ja>{I3}_{I1}{ALT} |
<CCi_ja>{I1}_{I4} |
<CCi_ja>{I2}_{I3} |
<CCi_ja>{I3}_{I2} |
<CCi_ja>{I4}_{I1} {
	  /* Invalid axis number in image-header keyword. */
	  if (relax & WCSHDR_ALLIMG) {
	    /* Will be flagged by <VALUE> as invalid. */
	    sscanf(yytext, "%d_%d", &i, &j);
	    keytype = IMGAXIS;
	    BEGIN(VALUE);
	
	  } else {
	    /* Pretend we don't recognize it. */
	    BEGIN(DISCARD);
	  }
	}

<CCi_ja>{I0}{6} {
	  /* This covers the defunct forms CD00i00j and PC00i00j. */
	  if (((relax & WCSHDR_PC00i00j) && (altlin == 1)) ||
	      ((relax & WCSHDR_CD00i00j) && (altlin == 2))) {
	    sscanf(yytext, "%3d%3d", &i, &j);
	    a = ' ';
	    keytype = IMGAXIS;
	    BEGIN(VALUE);
	
	  } else if (relax & WCSHDR_reject) {
	    errmsg = errtxt;
	    sprintf(errmsg, "Defunct form of %si_ja keyword",
	                     (altlin==1) ? "PC" : "CD");
	    BEGIN(ERROR);
	
	  } else {
	    /* Pretend we don't recognize it. */
	    BEGIN(DISCARD);
	  }
	}

<CCi_ja>. {
	  BEGIN(DISCARD);
	}

<ijCCna>{I1}{ALT}"  " |
<ijCCna>{I2}{ALT}" "  |
<ijCCna>{I3}{ALT} {
	  sscanf(yytext, "%d%c", &n, &a);
	  keytype = BIMGARR;
	  BEGIN(VALUE);
	}

<TCCn_ka>{I1}_{I1}{ALT}" " |
<TCCn_ka>{I1}_{I2}{ALT} |
<TCCn_ka>{I2}_{I1}{ALT} |
<TCCn_ka>{I1}_{I3} |
<TCCn_ka>{I2}_{I2} |
<TCCn_ka>{I3}_{I1} {
	  if (relax & WCSHDR_LONGKEY) {
	    WCSBTH_PUTBACK;
	    BEGIN(TCn_ka);
	
	  } else if (relax & WCSHDR_reject) {
	    errmsg = errtxt;
	    sprintf(errmsg, "%s keyword is non-standard", extkey);
	    BEGIN(ERROR);
	
	  } else {
	    /* Pretend we don't recognize it. */
	    BEGIN(DISCARD);
	  }
	}

<TCCn_ka>. {
	  BEGIN(DISCARD);
	}

<TCn_ka>{I1}_{I1}{ALT}"  " |
<TCn_ka>{I1}_{I2}{ALT}" " |
<TCn_ka>{I2}_{I1}{ALT}" " |
<TCn_ka>{I1}_{I3}{ALT} |
<TCn_ka>{I2}_{I2}{ALT} |
<TCn_ka>{I3}_{I1}{ALT} {
	  sscanf(yytext, "%d_%d%c", &n, &k, &a);
	  i = wcsbth_colax(*wcs, &alts, n, a);
	  j = wcsbth_colax(*wcs, &alts, k, a);
	  keytype = PIXLIST;
	  BEGIN(VALUE);
	}

<TCn_ka>{I1}_{I4} |
<TCn_ka>{I2}_{I3} |
<TCn_ka>{I3}_{I2} |
<TCn_ka>{I4}_{I1} {
	  sscanf(yytext, "%d_%d", &n, &k);
	  a = ' ';
	  i = wcsbth_colax(*wcs, &alts, n, a);
	  j = wcsbth_colax(*wcs, &alts, k, a);
	  keytype = PIXLIST;
	  BEGIN(VALUE);
	}

<TCn_ka>. {
	  BEGIN(DISCARD);
	}

<CROTAi>{I1}"  " |
<CROTAi>{I2}" " {
	  yyless(0);
	  BEGIN(CCCCCia);
	}

<CROTAi>{I1}[A-Z]" " |
<CROTAi>{I2}[A-Z] {
	  if (relax & WCSHDR_CROTAia) {
	    yyless(0);
	    BEGIN(CCCCCia);
	
	  } else if (relax & WCSHDR_reject) {
	    errmsg = "CROTAn keyword may not have an alternate version code";
	    BEGIN(ERROR);
	
	  } else {
	    /* Pretend we don't recognize it. */
	    BEGIN(DISCARD);
	  }
	}

<CROTAi>. {
	  BEGIN(DISCARD);
	}

<iCROTn>{I1}"  " |
<iCROTn>{I2}" "  |
<iCROTn>{I3}     |
<TCROTn>{I1}"  " |
<TCROTn>{I2}" "  |
<TCROTn>{I3} {
	  WCSBTH_PUTBACK;
	  BEGIN((YY_START == iCROTn) ? iCCCna : TCCCna);
	}

<iCROTn>{I1}[A-Z]" " |
<iCROTn>{I2}[A-Z]    |
<TCROTn>{I1}[A-Z]" " |
<TCROTn>{I2}[A-Z] {
	  if (relax & WCSHDR_CROTAia) {
	    WCSBTH_PUTBACK;
	    BEGIN((YY_START == iCROTn) ? iCCCna : TCCCna);
	
	  } else if (relax & WCSHDR_reject) {
	    errmsg = errtxt;
	    sprintf(errmsg,
	      "%s keyword may not have an alternate version code", extkey);
	    BEGIN(ERROR);
	
	  } else {
	    /* Pretend we don't recognize it. */
	    BEGIN(DISCARD);
	  }
	}

<iCROTn>. |
<TCROTn>. {
	  BEGIN(DISCARD);
	}

<CCCCCCCa>{ALT} |
<CCCCCCCC>. {
	  /* Image-header keyword. */
	  if (relax & (WCSHDR_AUXIMG | WCSHDR_ALLIMG)) {
	    if (YY_START == CCCCCCCa) {
	      sscanf(yytext, "%c", &a);
	    } else {
	      a = 0;
	      unput(yytext[0]);
	    }
	    keytype = IMGAUX;
	    BEGIN(VALUE);
	
	  } else if (relax & WCSHDR_reject) {
	    errmsg = errtxt;
	    sprintf(errmsg,
	      "Image-header keyword %s in binary table", extkey);
	    BEGIN(ERROR);
	
	  } else {
	    /* Pretend we don't recognize it. */
	    BEGIN(DISCARD);
	  }
	}

<CCCCCCCa>. {
	  BEGIN(DISCARD);
	}

<CCCCna>{I1}{ALT}"  " |
<CCCCna>{I2}{ALT}" "  |
<CCCCna>{I3}{ALT}     |
<CCCCCna>{I1}{ALT}" " |
<CCCCCna>{I2}{ALT} {
	  sscanf(yytext, "%d%c", &n, &a);
	  keytype = BINTAB;
	  BEGIN(VALUE);
	}

<CCCCCna>{I3} {
	  sscanf(yytext, "%d", &n);
	  a = ' ';
	  keytype = BINTAB;
	  BEGIN(VALUE);
	}

<CCCCna>. |
<CCCCCna>. {
	  BEGIN(DISCARD);
	}

<CCCCn>{I1}"   " |
<CCCCn>{I2}"  "  |
<CCCCn>{I3}" "   |
<CCCCn>{I4}      |
<CCCCCn>{I1}"  " |
<CCCCCn>{I2}" "  |
<CCCCCn>{I3} {
	  sscanf(yytext, "%d", &n);
	  a = 0;
	  keytype = BINTAB;
	  BEGIN(VALUE);
	}

<CCCCn>. |
<CCCCCn>. {
	  BEGIN(DISCARD);
	}

<CCi_ma>{I1}_{I0}{ALT}"  " |
<CCi_ma>{I1}_{I2}{ALT}" " |
<CCi_ma>{I2}_{I0}{ALT}" " |
<CCi_ma>{I2}_{I2}{ALT} {
	  /* Image-header keyword. */
	  if (relax & WCSHDR_ALLIMG) {
	    sscanf(yytext, "%d_%d%c", &i, &m, &a);
	    keytype = IMGAXIS;
	    BEGIN(VALUE);
	
	  } else if (relax & WCSHDR_reject) {
	    errmsg = errtxt;
	    sprintf(errmsg,
	      "Image-header keyword %s in binary table", extkey);
	    BEGIN(ERROR);
	
	  } else {
	    /* Pretend we don't recognize it. */
	    BEGIN(DISCARD);
	  }
	}

<CCi_ma>{I1}_{I3}{ALT} |
<CCi_ma>{I3}_{I0}{ALT} |
<CCi_ma>{I1}_{I4} |
<CCi_ma>{I2}_{I3} |
<CCi_ma>{I3}_{I2} |
<CCi_ma>{I4}_{I0} {
	  /* Invalid parameter in image-header keyword. */
	  if (relax & WCSHDR_ALLIMG) {
	    /* Will be flagged by <VALUE> as invalid. */
	    sscanf(yytext, "%d_%d", &i, &m);
	    keytype = IMGAXIS;
	    BEGIN(VALUE);
	
	  } else {
	    /* Pretend we don't recognize it. */
	    BEGIN(DISCARD);
	  }
	}

<CCi_ma>. {
	  BEGIN(DISCARD);
	}

<iCCn_ma>{I1}_{I0}{ALT}" " |
<iCCn_ma>{I1}_{I2}{ALT}    |
<iCCn_ma>{I1}_{I3}         |
<iCCn_ma>{I2}_{I0}{ALT}    |
<iCCn_ma>{I2}_{I2}         |
<iCCn_ma>{I3}_{I0}         |
<TCCn_ma>{I1}_{I0}{ALT}" " |
<TCCn_ma>{I1}_{I2}{ALT}    |
<TCCn_ma>{I1}_{I3}         |
<TCCn_ma>{I2}_{I0}{ALT}    |
<TCCn_ma>{I2}_{I2}         |
<TCCn_ma>{I3}_{I0} {
	  if (relax & WCSHDR_LONGKEY) {
	    WCSBTH_PUTBACK;
	    BEGIN((YY_START == iCCn_ma) ? iCn_ma : TCn_ma);
	
	  } else if (relax & WCSHDR_reject) {
	    errmsg = errtxt;
	    sprintf(errmsg, "%s keyword is non-standard", extkey);
	    BEGIN(ERROR);
	
	  } else {
	    /* Pretend we don't recognize it. */
	    BEGIN(DISCARD);
	  }
	}

<iCCn_ma>. |
<TCCn_ma>. {
	  BEGIN(DISCARD);
	}

<iCn_ma>{I1}_{I0}{ALT}"  " |
<iCn_ma>{I1}_{I2}{ALT}" "  |
<iCn_ma>{I1}_{I3}{ALT}     |
<iCn_ma>{I2}_{I0}{ALT}" "  |
<iCn_ma>{I2}_{I2}{ALT}     |
<iCn_ma>{I3}_{I0}{ALT}     |
<TCn_ma>{I1}_{I0}{ALT}"  " |
<TCn_ma>{I1}_{I2}{ALT}" "  |
<TCn_ma>{I1}_{I3}{ALT}     |
<TCn_ma>{I2}_{I0}{ALT}" "  |
<TCn_ma>{I2}_{I2}{ALT}     |
<TCn_ma>{I3}_{I0}{ALT} {
	  sscanf(yytext, "%d_%d%c", &n, &m, &a);
	  if (YY_START == TCn_ma) i = wcsbth_colax(*wcs, &alts, n, a);
	  keytype = (YY_START == iCn_ma) ? BIMGARR : PIXLIST;
	  BEGIN(VALUE);
	}

<iCn_ma>{I1}_{I4} |
<iCn_ma>{I2}_{I3} |
<iCn_ma>{I3}_{I2} |
<iCn_ma>{I4}_{I0} |
<TCn_ma>{I1}_{I4} |
<TCn_ma>{I2}_{I3} |
<TCn_ma>{I3}_{I2} |
<TCn_ma>{I4}_{I0} {
	  /* Invalid combinations will be flagged by <VALUE>. */
	  sscanf(yytext, "%d_%d", &n, &m);
	  a = ' ';
	  if (YY_START == TCn_ma) i = wcsbth_colax(*wcs, &alts, n, a);
	  keytype = (YY_START == iCn_ma) ? BIMGARR : PIXLIST;
	  BEGIN(VALUE);
	}

<iCn_ma>. |
<TCn_ma>. {
	  BEGIN(DISCARD);
	}

<PROJPm>{I0}"  " {
	  if (relax & WCSHDR_PROJPn) {
	    sscanf(yytext, "%d", &m);
	    i = 0;
	    a = ' ';
	    keytype = IMGAXIS;
	    BEGIN(VALUE);
	
	  } else if (relax & WCSHDR_reject) {
	    errmsg = "PROJPn keyword is defunct";
	    BEGIN(ERROR);
	
	  } else {
	    /* Pretend we don't recognize it. */
	    BEGIN(DISCARD);
	  }
	}

<PROJPm>. {
	  BEGIN(DISCARD);
	}

<VALUE>=" "+ {
	  /* Do checks on i, j, m, n, k. */
	  if (!(keytype & keysel)) {
	    /* Selection by keyword type. */
	    BEGIN(DISCARD);
	
	  } else if (exclude[n] || exclude[k]) {
	    /* One or other column is not selected. */
	    if (k && (exclude[n] != exclude[k])) {
	      /* For keywords such as TCn_ka, both columns must be excluded.
	         User error, so return immediately. */
	      yylex_destroy();
	      return 3;
	
	    } else {
	      BEGIN(DISCARD);
	    }
	
	  } else if (i > 99 || j > 99 || m > 99 || n > 999 || k > 999) {
	    if (relax & WCSHDR_reject) {
	      errmsg = errtxt;
	      if (i > 99 || j > 99) {
	        sprintf(errmsg, "Axis number exceeds 99");
	      } else if (m > 99) {
	        sprintf(errmsg, "Parameter number exceeds 99");
	      } else if (n > 999 || k > 999) {
	        sprintf(errmsg, "Column number exceeds 999");
	      }
	      BEGIN(ERROR);
	
	    } else {
	      /* Pretend we don't recognize it. */
	      BEGIN(DISCARD);
	    }
	
	  } else if (ipass == 2 && npass == 3 && (keytype & BINTAB)) {
	    /* Skip keyvalues that won't be inherited. */
	    BEGIN(FLUSH);
	
	  } else if (ipass == 3 && (keytype & IMGHEAD)) {
	    /* IMGHEAD keytypes are always dealt with on the second pass. */
	    BEGIN(FLUSH);
	
	  } else if (vptr) {
	    alts.icol = 0;
	    alts.ialt = 0;
	    voff = (char *)vptr - (char *)(&wcstem);
	
	    if (valtype == INTEGER) {
	      BEGIN(INTEGER_VAL);
	    } else if (valtype == FLOAT) {
	      BEGIN(FLOAT_VAL);
	    } else if (valtype == STRING) {
	      BEGIN(STRING_VAL);
	    } else {
	      errmsg = errtxt;
	      sprintf(errmsg, "Internal parser ERROR, bad data type: %d",
	        valtype);
	      BEGIN(ERROR);
	    }
	
	  } else {
	    errmsg = "Internal parser ERROR, null pointer";
	    BEGIN(ERROR);
	  }
	}

<VALUE>. {
	  errmsg = "Invalid KEYWORD = VALUE syntax";
	  BEGIN(ERROR);
	}

<INTEGER_VAL>{INTEGER} {
	  if (ipass == 1) {
	    /* Do first-pass bookkeeping. */
	    wcsbth_pass1(keytype, i, j, n, k, a, ptype, &alts);
	    BEGIN(FLUSH);
	
	  } else {
	    /* Update each coordinate representation. */
	    while ((wcsp = wcsbth_idx(*wcs, &alts, keytype, n, a))) {
	      wptr = (void *)((char *)wcsp + voff);
	
	      /* Read the keyvalue. */
	      if (special) {
	        special(wptr);
	      } else {
	        sscanf(yytext, "%d", (int *)wptr);
	      }
	    }
	
	    BEGIN(COMMENT);
	  }
	}

<INTEGER_VAL>. {
	  errmsg = "An integer value was expected";
	  BEGIN(ERROR);
	}

<FLOAT_VAL>{FLOAT} {
	  if (ipass == 1) {
	    /* Do first-pass bookkeeping. */
	    wcsbth_pass1(keytype, i, j, n, k, a, ptype, &alts);
	    BEGIN(FLUSH);
	
	  } else {
	    /* Update each coordinate representation. */
	    while ((wcsp = wcsbth_idx(*wcs, &alts, keytype, n, a))) {
	      wptr = (void *)((char *)wcsp + voff);
	
	      /* Apply keyword parameterization. */
	      if (ptype == 'v') {
	        ipx = wcsp->npv++;
	        wcsp->pv[ipx].i = i;
	        wcsp->pv[ipx].m = m;
	        wptr = &(wcsp->pv[ipx].value);
	
	      } else if (j) {
	        /* Is the de-reference necessary? */
	        wptr = *((double **)wptr) + (i - 1)*(wcsp->naxis) + (j - 1);
	
	      } else if (i) {
	        wptr = *((double **)wptr) + (i - 1);
	      }
	
	      /* Read the keyvalue. */
	      if (special) {
	        special(wptr);
	      } else {
	        wcsutil_str2double(yytext, "%lf", (double *)wptr);
	      }
	
	      /* Flag the presence of PC, or CD and/or CROTA. */
	      if (altlin) {
	        wcsp->altlin |= altlin;
	        altlin = 0;
	      }
	    }
	
	    BEGIN(COMMENT);
	  }
	}

<FLOAT_VAL>. {
	  errmsg = "A floating-point value was expected";
	  BEGIN(ERROR);
	}

<STRING_VAL>{STRING} {
	  if (ipass == 1) {
	    /* Do first-pass bookkeeping. */
	    wcsbth_pass1(keytype, i, j, n, k, a, ptype, &alts);
	    BEGIN(FLUSH);
	
	  } else {
	    /* Update each coordinate representation. */
	    while ((wcsp = wcsbth_idx(*wcs, &alts, keytype, n, a))) {
	      wptr = (void *)((char *)wcsp + voff);
	
	      /* Apply keyword parameterization. */
	      if (ptype == 's') {
	        ipx = wcsp->nps++;
	        wcsp->ps[ipx].i = i;
	        wcsp->ps[ipx].m = m;
	        wptr = wcsp->ps[ipx].value;
	
	      } else if (j) {
	        wptr = *((char (**)[72])wptr) +
	                (i - 1)*(wcsp->naxis) + (j - 1);
	
	      } else if (i) {
	        wptr = *((char (**)[72])wptr) + (i - 1);
	      }
	
	      /* Read the keyvalue. */
	      cptr = (char *)wptr;
	      strcpy(cptr, yytext+1);
	
	      /* Squeeze out repeated quotes. */
	      ix = 0;
	      for (jx = 0; jx < 72; jx++) {
	        if (ix < jx) {
	          cptr[ix] = cptr[jx];
	        }
	
	        if (cptr[jx] == '\0') {
	          if (ix) cptr[ix-1] = '\0';
	          break;
	        } else if (cptr[jx] == '\'' && cptr[jx+1] == '\'') {
	          jx++;
	        }
	
	        ix++;
	      }
	    }
	
	    BEGIN(COMMENT);
	  }
	}

<STRING_VAL>. {
	  errmsg = "A string value was expected";
	  BEGIN(ERROR);
	}

<COMMENT>" "*\/.* |
<COMMENT>" "* {
	  BEGIN(FLUSH);
	}

<COMMENT>. {
	  errmsg = "Malformed keycomment";
	  BEGIN(ERROR);
	}

<DISCARD>.* {
	  if (ipass == npass) {
	    if (ctrl < 0) {
	      /* Preserve discards. */
	      if (hptr < wcsbth_hdr-80) {
	        strncpy(hptr, wcsbth_hdr-80, 80);
	      }
	      hptr += 80;
	
	    } else if (ctrl > 2) {
	      wcsfprintf(stderr, "%.80s\n  Discarded.\n", wcsbth_hdr-80);
	    }
	  }
	
	  BEGIN(FLUSH);
	}

<ERROR>.* {
	  (*nreject)++;
	  if (ipass == npass) {
	    if (ctrl == -1) {
	      if (hptr < wcsbth_hdr-80) {
	        /* Preserve rejects. */
	        strncpy(hptr, wcsbth_hdr-80, 80);
	      }
	      hptr += 80;
	    }
	
	    if (abs(ctrl) > 1) {
	      wcsfprintf(stderr, "%.80s\n%4d: %s.\n", wcsbth_hdr-80, *nreject,
	        errmsg);
	    }
	  }
	
	  BEGIN(FLUSH);
	}

<FLUSH>.*\n {
	  /* Throw away the rest of the line and reset for the next one. */
	  i = j = 0;
	  n = k = 0;
	  m = 0;
	  a = ' ';
	
	  keytype =  0;
	  valtype = -1;
	  vptr    = 0x0;
	
	  altlin = 0;
	  ptype  = ' ';
	  special = 0x0;
	  BEGIN(INITIAL);
	}

<<EOF>>	 {
	  /* End-of-input. */
	  if (ipass == 1) {
	    if ((status = wcsbth_init1(&alts, nwcs, wcs)) || *nwcs == 0) {
	      yylex_destroy();
	      return status;
	    }
	
	    if (alts.imgherit) npass = 3;
	
	    if (abs(ctrl) > 2) {
	      if (*nwcs == 1) {
	        wcsfprintf(stderr, "Found one coordinate representation.\n");
	      } else {
	        wcsfprintf(stderr, "Found %d coordinate representations.\n",
	          *nwcs);
	      }
	    }
	  }
	
	  if (ipass++ < npass) {
	    wcsbth_hdr = header;
	    wcsbth_nkeyrec = nkeyrec;
	    *nreject = 0;
	
	    i = j = 0;
	    k = n = 0;
	    m = 0;
	    a = ' ';
	
	    keytype =  0;
	    valtype = -1;
	    vptr    = 0x0;
	
	    altlin = 0;
	    ptype  = ' ';
	    special = 0x0;
	
	    yyrestart(yyin);
	
	  } else {
	    yylex_destroy();
	
	    if (ctrl < 0) {
	      *hptr = '\0';
	    } else if (ctrl == 1) {
	      wcsfprintf(stderr, "%d WCS keyrecords were rejected.\n",
	        *nreject);
	    }
	
	    return wcsbth_final(&alts, nwcs, wcs);
	  }
	}

%%

/*----------------------------------------------------------------------------
* Perform first-pass tasks:
*
* 1) Count the number of coordinate axes in each of the 27 possible alternate
*    image-header coordinate representations.  Also count the number of PVi_ma
*    and PSi_ma keywords in each representation.
*
* 2) Determine the number of binary table columns that have an image array
*    with a coordinate representation (up to 999), and count the number of
*    coordinate axes in each of the 27 possible alternates.  Also count the
*    number of iVn_ma and iSn_ma keywords in each representation.
*
* 3) Determine the number of alternate pixel list coordinate representations
*    (up to 27) and the table columns associated with each.  Also count the
*    number of TVn_ma and TSn_ma keywords in each representation.
*
* In the first pass alts->arridx[icol][27] is used to determine the number of
* axes in each of 27 possible image-header coordinate descriptions (icol == 0)
* and each of the 27 possible coordinate representations for an image array in
* each column.
*
* The elements of alts->pixlist[icol] are used as bit arrays to flag which of
* the 27 possible pixel list coordinate representations are associated with
* each table column.
*---------------------------------------------------------------------------*/

int wcsbth_pass1(
  int keytype,
  int i,
  int j,
  int n,
  int k,
  char a,
  char ptype,
  struct wcsbth_alts *alts)

{
  int ialt, icol, mask, ncol;

  if (a == 0) {
    /* Keywords such as DATE-OBS go along for the ride. */
    return 0;
  }

  ncol = alts->ncol;

  /* Do we need to allocate memory for alts? */
  if (alts->arridx == 0x0) {
    if (ncol == 0) {
      /* Can only happen if TFIELDS is missing or out-of-sequence.  If n and
         k are both zero then we may be processing an image header so leave
         ncol alone - the array will be realloc'd later if required. */
      if (n || k) {
        /* The header is mangled, assume the worst. */
        ncol = 999;
      }
    }

    if (!(alts->arridx  =  calloc((1 + ncol)*27, sizeof(short int))) ||
        !(alts->npv     =  calloc((1 + ncol)*27, sizeof(unsigned char)))  ||
        !(alts->nps     =  calloc((1 + ncol)*27, sizeof(unsigned char)))  ||
        !(alts->pixlist =  calloc((1 + ncol),    sizeof(unsigned int)))) {
      if (alts->arridx)  free(alts->arridx);
      if (alts->npv)     free(alts->npv);
      if (alts->nps)     free(alts->nps);
      if (alts->pixlist) free(alts->pixlist);
      return 2;
    }

    alts->ncol = ncol;

  } else if (n > ncol || k > ncol) {
    /* Can only happen if TFIELDS or the WCS keyword is wrong; carry on. */
    ncol = 999;
    if (!(alts->arridx  = realloc(alts->arridx,
                                    27*(1 + ncol)*sizeof(short int))) ||
        !(alts->npv     = realloc(alts->npv,
                                    27*(1 + ncol)*sizeof(unsigned char)))  ||
        !(alts->nps     = realloc(alts->nps,
                                    27*(1 + ncol)*sizeof(unsigned char)))  ||
        !(alts->pixlist = realloc(alts->pixlist,
                                       (1 + ncol)*sizeof(unsigned int)))) {
      if (alts->arridx)  free(alts->arridx);
      if (alts->npv)     free(alts->npv);
      if (alts->nps)     free(alts->nps);
      if (alts->pixlist) free(alts->pixlist);
      return 2;
    }

    /* Since realloc() doesn't initialize the extra memory. */
    for (icol = (1 + alts->ncol); icol < (1 + ncol); icol++) {
      for (ialt = 0; ialt < 27; ialt++) {
        alts->arridx[icol][ialt] = 0;
        alts->npv[icol][ialt] = 0;
        alts->nps[icol][ialt] = 0;
        alts->pixlist[icol]   = 0;
      }
    }

    alts->ncol = ncol;
  }

  ialt = 0;
  if (a != ' ') {
    ialt = a - 'A' + 1;
  }

  /* A BINTAB keytype such as LONPna, in conjunction with an IMGAXIS keytype
     causes a table column to be recognized as an image array. */
  if (keytype & IMGHEAD || keytype & BIMGARR) {
    /* n == 0 is expected for IMGHEAD keywords. */
    if (i == 0 && j == 0) {
      if (alts->arridx[n][ialt] == 0) {
        /* Flag that an auxiliary keyword was seen. */
        alts->arridx[n][ialt] = -1;
      }

    } else {
      /* Record the maximum axis number found. */
      if (alts->arridx[n][ialt] < i) {
        alts->arridx[n][ialt] = i;
      }

      if (alts->arridx[n][ialt] < j) {
        alts->arridx[n][ialt] = j;
      }
    }

    if (ptype == 'v') {
      alts->npv[n][ialt]++;
    } else if (ptype == 's') {
      alts->nps[n][ialt]++;
    }
  }

  /* BINTAB keytypes, which apply both to pixel lists as well as binary table
     image arrays, never contribute to recognizing a table column as a pixel
     list axis.  A PIXLIST keytype is required for that. */
  if (keytype == PIXLIST) {
    mask = 1 << ialt;

    /* n > 0 for PIXLIST keytypes. */
    alts->pixlist[n] |= mask;
    if (k) alts->pixlist[k] |= mask;

    /* Used as a flag over all columns. */
    alts->pixlist[0] |= mask;

    if (ptype == 'v') {
      alts->pixnpv[ialt]++;
    } else if (ptype == 's') {
      alts->pixnps[ialt]++;
    }
  }

  return 0;
}


/*----------------------------------------------------------------------------
* Perform initializations at the end of the first pass:
*
* 1) Determine the required number of wcsprm structs, allocate memory for
*    an array of them and initialize each one.
*---------------------------------------------------------------------------*/

int wcsbth_init1(
  struct wcsbth_alts *alts,
  int *nwcs,
  struct wcsprm **wcs)

{
  int ialt, icol, inherit[27], ix, mask, ncol, npsmax, npvmax, status = 0;
  struct wcsprm *wcsp;


  if (alts->arridx == 0x0) {
    *nwcs = 0;
    return 0;
  }

  /* Determine the number of axes in each pixel list representation. */
  ncol = alts->ncol;
  for (ialt = 0, mask = 1; ialt < 27; ialt++, mask <<= 1) {
    alts->pixidx[ialt] = 0;

    if (alts->pixlist[0] | mask) {
      for (icol = 1; icol <= ncol; icol++) {
        if (alts->pixlist[icol] & mask) {
          alts->pixidx[ialt]++;
        }
      }
    }
  }

  /* Find the total number of coordinate representations. */
  *nwcs = 0;
  alts->imgherit = 0;
  for (ialt = 0; ialt < 27; ialt++) {
    inherit[ialt] = 0;

    for (icol = 1; icol <= ncol; icol++) {
      if (alts->arridx[icol][ialt] < 0) {
        /* No BIMGARR keytype but there's at least one BINTAB. */
        if (alts->arridx[0][ialt] > 0) {
          /* There is an IMGAXIS keytype that we will inherit, so count this
             representation. */
          alts->arridx[icol][ialt] = alts->arridx[0][ialt];
        } else {
          alts->arridx[icol][ialt] = 0;
        }
      }

      if (alts->arridx[icol][ialt]) {
        if (alts->arridx[0][ialt]) {
          /* All IMGHEAD keywords are inherited for this ialt. */
          inherit[ialt] = 1;

          if (alts->arridx[icol][ialt] < alts->arridx[0][ialt]) {
            /* The extra axes are also inherited. */
            alts->arridx[icol][ialt] = alts->arridx[0][ialt];
          }
        }

        (*nwcs)++;
      }
    }

    /* Count every "a" found in any IMGHEAD keyword... */
    if (alts->arridx[0][ialt]) {
      if (inherit[ialt]) {
        /* ...but not if the IMGHEAD keywords will be inherited. */
        alts->arridx[0][ialt] = 0;
        alts->imgherit = 1;
      } else {
        (*nwcs)++;
      }
    }

    /* We need a struct for every "a" found in a PIXLIST keyword. */
    if (alts->pixidx[ialt]) {
      (*nwcs)++;
    }
  }


  if (*nwcs) {
    /* Allocate memory for the required number of wcsprm structs. */
    if (!(*wcs = calloc(*nwcs, sizeof(struct wcsprm)))) {
      return 2;
    }

    /* Record the current values of NPVMAX and NPSMAX. */
    npvmax = wcsnpv(-1);
    npsmax = wcsnps(-1);

    /* Initialize each wcsprm struct. */
    wcsp = *wcs;
    *nwcs = 0;
    for (icol = 0; icol <= ncol; icol++) {
      for (ialt = 0; ialt < 27; ialt++) {
        if (alts->arridx[icol][ialt]) {
          /* Image-header representations that are not for inheritance
             (icol == 0) or binary table image array representations. */
          wcsp->flag = -1;
          wcsnpv(alts->npv[icol][ialt]);
          wcsnps(alts->nps[icol][ialt]);
          if ((status = wcsini(1, (int)(alts->arridx[icol][ialt]), wcsp))) {
            wcsvfree(nwcs, wcs);
            break;
          }

          /* Record the alternate version code. */
          if (ialt) {
            wcsp->alt[0] = 'A' + ialt - 1;
          }

          /* Record the table column number. */
          wcsp->colnum = icol;

          /* On the second pass alts->arridx[icol][27] indexes the array of
             wcsprm structs. */
          alts->arridx[icol][ialt] = (*nwcs)++;

          wcsp++;

        } else {
          /* Signal that this column has no WCS for this "a". */
          alts->arridx[icol][ialt] = -1;
        }
      }
    }

    for (ialt = 0; ialt < 27; ialt++) {
      if (alts->pixidx[ialt]) {
        /* Pixel lists representations. */
        wcsp->flag = -1;
        wcsnpv(alts->pixnpv[ialt]);
        wcsnps(alts->pixnps[ialt]);
        if ((status = wcsini(1, (int)(alts->pixidx[ialt]), wcsp))) {
          wcsvfree(nwcs, wcs);
          break;
        }

        /* Record the alternate version code. */
        if (ialt) {
          wcsp->alt[0] = 'A' + ialt - 1;
        }

        /* Record the pixel list column numbers. */
        mask = (1 << ialt);
        for (icol = 1, ix = 0; icol <= ncol; icol++) {
          if (alts->pixlist[icol] & mask) {
            wcsp->colax[ix++] = icol;
          }
        }

        /* alts->pixidx[] indexes the array of wcsprm structs. */
        alts->pixidx[ialt] = (*nwcs)++;

        wcsp++;

      } else {
        /* Signal that this column is not a pixel list axis for this "a". */
        alts->pixidx[ialt] = -1;
      }
    }

    /* Restore the original values of NPVMAX and NPSMAX. */
    wcsnpv(npvmax);
    wcsnps(npsmax);
  }

  return status;
}


/*----------------------------------------------------------------------------
* Return a pointer to the next wcsprm struct for a particular column number
* and alternate.
*---------------------------------------------------------------------------*/

struct wcsprm *wcsbth_idx(
  struct wcsprm *wcs,
  struct wcsbth_alts *alts,
  int  keytype,
  int  n,
  char a)

{
  const char as[] = " ABCDEFGHIJKLMNOPQRSTUVWXYZ";
  int iwcs;

  if (!wcs) return 0x0;

  iwcs = -1;
  for (; iwcs < 0 && alts->ialt < 27; alts->ialt++) {
    /* Note that a == 0 applies to every alternate, otherwise this
       loop simply determines the appropriate value of alts->ialt. */
    if (a && a != as[alts->ialt]) continue;

    if (keytype & (IMGHEAD | BIMGARR)) {
      for (; iwcs < 0 && alts->icol <= alts->ncol; alts->icol++) {
        /* Image header keywords, n == 0, apply to all columns, otherwise this
           loop simply determines the appropriate value of alts->icol. */
        if (n && n != alts->icol) continue;
        iwcs = alts->arridx[alts->icol][alts->ialt];
      }

      /* Break out of the loop to stop alts->ialt from being incremented. */
      if (iwcs >= 0) break;

      /* Start from scratch for the next alts->ialt. */
      alts->icol = 0;
    }

    if (keytype & (IMGAUX | PIXLIST)) {
      iwcs = alts->pixidx[alts->ialt];
    }
  }

  return (iwcs >= 0) ? (wcs + iwcs) : 0x0;
}


/*----------------------------------------------------------------------------
* Return the axis number associated with the specified column number in a
* particular pixel list coordinate representation.
*---------------------------------------------------------------------------*/

int wcsbth_colax(
  struct wcsprm *wcs,
  struct wcsbth_alts *alts,
  int n,
  char a)

{
  int ix;
  struct wcsprm *wcsp;

  if (!wcs) return 0;

  wcsp = wcs;
  if (a != ' ') {
    wcsp += alts->pixidx[a-'A'+1];
  }

  for (ix = 0; ix < wcsp->naxis; ix++) {
    if (wcsp->colax[ix] == n) {
      return ++ix;
    }
  }

  return 0;
}


/*----------------------------------------------------------------------------
* Interpret EPOCH keywords.
*---------------------------------------------------------------------------*/

int wcsbth_epoch(void *wptr)

{
  double *equinox;

  /* If EQUINOXa is currently undefined then set it from EPOCHa. */
  equinox = (double *)wptr;
  if (undefined(*equinox)) {
    wcsutil_str2double(yytext, "%lf", equinox);
  }

  return 0;
}


/*----------------------------------------------------------------------------
* Interpret VSOURCE keywords.
*---------------------------------------------------------------------------*/

int wcsbth_vsource(void *wptr)

{
  double beta, c = 299792458.0, vsource, *zsource;

  /* If ZSOURCEa is currently undefined then set it from VSOURCEa. */
  zsource = (double *)wptr;
  if (undefined(*zsource)) {
    wcsutil_str2double(yytext, "%lf", &vsource);

    /* Convert relativistic Doppler velocity to redshift. */
    beta = vsource/c;
    *zsource = (1.0 + beta)/sqrt(1.0 - beta*beta) - 1.0;
  }

  return 0;
}


/*----------------------------------------------------------------------------
* Tie up loose ends.
*---------------------------------------------------------------------------*/

int wcsbth_final(
  struct wcsbth_alts *alts,
  int *nwcs,
  struct wcsprm **wcs)

{
  int ialt, status;

  if (alts->arridx)  free(alts->arridx);
  if (alts->npv)     free(alts->npv);
  if (alts->nps)     free(alts->nps);
  if (alts->pixlist) free(alts->pixlist);

  for (ialt = 0; ialt < *nwcs; ialt++) {
    /* Interpret -TAB header keywords. */
    if ((status = wcstab(*wcs+ialt))) {
       wcsvfree(nwcs, wcs);
       return status;
    }
  }

  return 0;
}
