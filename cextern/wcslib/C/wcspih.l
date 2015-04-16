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
  $Id: wcspih.l,v 5.2 2015/04/15 12:35:07 mcalabre Exp $
*=============================================================================
*
* wcspih.l is a Flex description file containing the definition of a lexical
* scanner for parsing the WCS keyrecords from a FITS primary image or image
* extension header.
*
* wcspih.l requires Flex v2.5.4 or later.  Refer to wcshdr.h for a description
* of the user interface and operating notes.
*
* Implementation notes
* --------------------
* Use of the WCSAXESa keyword is not mandatory.  Its default value is "the
* larger of NAXIS and the largest index of these keywords [i.e. CRPIXj, PCi_j
* or CDi_j, CDELTi, CTYPEi, CRVALi, and CUNITi] found in the FITS header".
* Consequently the definition of WCSAXESa effectively invalidates the use of
* NAXIS for determining the number of coordinate axes and forces a preliminary
* pass through the header to determine the "largest index" in headers where
* WCSAXESa was omitted.
*
* Furthermore, since the use of WCSAXESa is optional, there is no way to
* determine the number of coordinate representations (the "a" value) other
* than by parsing all of the WCS keywords in the header; even if WCSAXESa was
* specified for some representations it cannot be known in advance whether it
* was specified for all of those present in the header.
*
* Hence the definition of WCSAXESa forces the scanner to be implemented in two
* passes.  The first pass is used to determine the number of coordinate
* representations (up to 27) and the number of coordinate axes in each.
* Effectively WCSAXESa is ignored unless it exceeds the "largest index" in
* which case the keywords for the extra axes assume their default values.  The
* number of PVi_ma and PSi_ma keywords in each representation is also counted
* in the first pass.
*
* On completion of the first pass, memory is allocated for an array of the
* required number of wcsprm structs and each of these is initialized
* appropriately.  These structs are filled in the second pass.
*
* The parser does not check for duplicated keywords, it accepts the last
* encountered.
*
*===========================================================================*/

/* Options. */
%option full
%option never-interactive
%option noyywrap
%option outfile="wcspih.c"
%option prefix="wcspih"

/* Indices for parameterized keywords. */
Z1	[0-9]
Z2	[0-9]{2}
Z3	[0-9]{3}
Z4	[0-9]{4}

I1	[1-9]
I2	[1-9][0-9]
I3	[1-9][0-9]{2}
I4	[1-9][0-9]{3}

/* Alternate coordinate system identifier. */
ALT	[ A-Z]

/* Keyvalue data types. */
INTEGER	[+-]?[0-9]+
FLOAT	[+-]?([0-9]+\.?[0-9]*|\.[0-9]+)([eEdD][+-]?[0-9]+)?
STRING	'([^']|'')*'

/* Inline comment syntax. */
INLINE " "*(\/.*)?

/* Exclusive start states. */
%x CROTAi PROJPn
%x CCCCCia CCi_ja CCi_ma CCCCCCCa CCCCCCCC
%x VALUE
%x INTEGER_VAL FLOAT_VAL STRING_VAL
%x COMMENT
%x DISCARD ERROR FLUSH

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

#define INTEGER 0
#define FLOAT   1
#define STRING  2

#define YY_DECL int wcspih(char *header, int nkeyrec, int relax, int ctrl, \
                           int *nreject, int *nwcs, struct wcsprm **wcs)

#define YY_INPUT(inbuff, count, bufsize) \
	{ \
	  if (wcspih_nkeyrec) { \
	    strncpy(inbuff, wcspih_hdr, 80); \
	    inbuff[80] = '\n'; \
	    wcspih_hdr += 80; \
	    wcspih_nkeyrec--; \
	    count = 81; \
	  } else { \
	    count = YY_NULL; \
	  } \
	}

/* These global variables are required by YY_INPUT. */
char *wcspih_hdr;
int  wcspih_nkeyrec;

int wcspih_final(int alts[], int *nwcs, struct wcsprm **wcs);
int wcspih_inits(int naxis, int alts[], int npv[], int nps[], int *nwcs,
        struct wcsprm **wcs);
void wcspih_naxes(int naxis, int i, int j, char a, int alts[], int *npptr);

int wcspih_epoch(double *wptr, double epoch);
int wcspih_vsource(double *wptr, double vsource);

/* Used in preempting the call to exit() by yy_fatal_error(). */
jmp_buf wcspih_abort_jmp_env;
#define exit(status) longjmp(wcspih_abort_jmp_env, status)

%}

%%
	/* Keyword indices, as used in the WCS papers, e.g. PCi_ja, PVi_ma. */
	char a;
	int  i, j, m;
	
	char *cptr, *errmsg, errtxt[80], *hptr, *keep, *keyname, *keyrec,
	     strtmp[80];
	int  altlin, alts[27], gotone, ialt, inttmp, ipass, ipx, ix, jx,
	     naxis, nother, *npptr, nps[27], npass, npv[27], nvalid, status,
	     valtype, voff;
	double dbltmp;
	void *vptr, *wptr;
	struct wcsprm *wcsp, wcstem;
	int (*special)(double *, double);
	int yylex_destroy(void);
	
	naxis = 0;
	for (ialt = 0; ialt < 27; ialt++) {
	  alts[ialt] = 0;
	  npv[ialt] = 0;
	  nps[ialt] = 0;
	}
	
	/* Parameters used to implement YY_INPUT. */
	wcspih_hdr = header;
	wcspih_nkeyrec = nkeyrec;
	
	/* Our handle on the input stream. */
	keyrec = header;
	hptr = header;
	keep = 0x0;
	
	/* For keeping tallies of keywords found. */
	*nreject = 0;
	nvalid = 0;
	nother = 0;
	
	/* If strict, then also reject. */
	if (relax & WCSHDR_strict) relax |= WCSHDR_reject;
	
	/* Keyword parameters. */
	i = j = 0;
	m = 0;
	a = ' ';
	
	/* For decoding the keyvalue. */
	valtype = -1;
	vptr    = 0x0;
	
	/* For keywords that require special handling. */
	altlin  = 0;
	npptr   = 0x0;
	special = 0x0;
	
	/* The data structures produced. */
	*nwcs = 0;
	*wcs  = 0x0;
	
	/* Control variables. */
	ipass = 1;
	npass = 2;
	
	/* Return here via longjmp() invoked by yy_fatal_error(). */
	if (setjmp(wcspih_abort_jmp_env)) {
	  return 3;
	}
	
	BEGIN(INITIAL);


^NAXIS"   = "" "*{INTEGER}{INLINE} {
	  keyname = "NAXISn";
	
	  if (ipass == 1) {
	    sscanf(yytext, "NAXIS   = %d", &naxis);
	    if (naxis < 0) naxis = 0;
	    BEGIN(FLUSH);
	
	  } else {
	    sscanf(yytext, "NAXIS   = %d", &i);
	
	    if (i < 0) {
	      errmsg = "negative value of NAXIS ignored";
	      BEGIN(ERROR);
	    } else {
	      BEGIN(DISCARD);
	    }
	  }
	}

^WCSAXES{ALT}=" "" "*{INTEGER} {
	  sscanf(yytext, "WCSAXES%c= %d", &a, &i);
	
	  if (i < 0) {
	    errmsg = "negative value of WCSAXESa ignored";
	    BEGIN(ERROR);
	
	  } else {
	    valtype = INTEGER;
	    vptr    = 0x0;
	
	    keyname = "WCSAXESa";
	    BEGIN(COMMENT);
	  }
	}

^CRPIX	{
	  valtype = FLOAT;
	  vptr    = &(wcstem.crpix);
	
	  keyname = "CRPIXja";
	  BEGIN(CCCCCia);
	}

^PC	{
	  valtype = FLOAT;
	  vptr    = &(wcstem.pc);
	  altlin = 1;
	
	  keyname = "PCi_ja";
	  BEGIN(CCi_ja);
	}

^CD	{
	  valtype = FLOAT;
	  vptr    = &(wcstem.cd);
	  altlin = 2;
	
	  keyname = "CDi_ja";
	  BEGIN(CCi_ja);
	}

^CDELT	{
	  valtype = FLOAT;
	  vptr    = &(wcstem.cdelt);
	
	  keyname = "CDELTia";
	  BEGIN(CCCCCia);
	}

^CROTA	{
	  valtype = FLOAT;
	  vptr    = &(wcstem.crota);
	  altlin = 4;
	
	  keyname = "CROTAn";
	  BEGIN(CROTAi);
	}

^CUNIT	{
	  valtype = STRING;
	  vptr    = &(wcstem.cunit);
	
	  keyname = "CUNITia";
	  BEGIN(CCCCCia);
	}

^CTYPE	{
	  valtype = STRING;
	  vptr    = &(wcstem.ctype);
	
	  keyname = "CTYPEia";
	  BEGIN(CCCCCia);
	}

^CRVAL	{
	  valtype = FLOAT;
	  vptr    = &(wcstem.crval);
	
	  keyname = "CRVALia";
	  BEGIN(CCCCCia);
	}

^LONPOLE {
	  valtype = FLOAT;
	  vptr    = &(wcstem.lonpole);
	
	  keyname = "LONPOLEa";
	  BEGIN(CCCCCCCa);
	}

^LATPOLE {
	  valtype = FLOAT;
	  vptr    = &(wcstem.latpole);
	
	  keyname = "LATPOLEa";
	  BEGIN(CCCCCCCa);
	}

^RESTFRQ {
	  valtype = FLOAT;
	  vptr    = &(wcstem.restfrq);
	
	  keyname = "RESTFRQa";
	  BEGIN(CCCCCCCa);
	}

^RESTFREQ {
	  if (relax & WCSHDR_strict) {
	    errmsg = "the RESTFREQ keyword is deprecated, use RESTFRQa";
	    BEGIN(ERROR);
	
	  } else {
	    valtype = FLOAT;
	    vptr    = &(wcstem.restfrq);
	
	    unput(' ');
	
	    keyname = "RESTFREQ";
	    BEGIN(CCCCCCCa);
	  }
	}

^RESTWAV {
	  valtype = FLOAT;
	  vptr    = &(wcstem.restwav);
	
	  keyname = "RESTWAVa";
	  BEGIN(CCCCCCCa);
	}

^PV	{
	  valtype = FLOAT;
	  vptr    = &(wcstem.pv);
	  npptr = npv;
	
	  keyname = "PVi_ma";
	  BEGIN(CCi_ma);
	}

^PROJP	{
	  valtype = FLOAT;
	  vptr    = &(wcstem.pv);
	  npptr = npv;
	
	  keyname = "PROJPn";
	  BEGIN(PROJPn);
	}

^PS	{
	  valtype = STRING;
	  vptr    = &(wcstem.ps);
	  npptr = nps;
	
	  keyname = "PSi_ma";
	  BEGIN(CCi_ma);
	}

^CNAME	{
	  valtype = STRING;
	  vptr    = &(wcstem.cname);
	
	  keyname = "CNAMEia";
	  BEGIN(CCCCCia);
	}

^CRDER	{
	  valtype = FLOAT;
	  vptr    = &(wcstem.crder);
	
	  keyname = "CRDERia";
	  BEGIN(CCCCCia);
	}

^CSYER	{
	  valtype = FLOAT;
	  vptr    = &(wcstem.csyer);
	
	  keyname = "CSYERia";
	  BEGIN(CCCCCia);
	}

^DATE-AVG {
	  valtype = STRING;
	  vptr    = wcstem.dateavg;
	  if (ctrl < -10) keep = keyrec;
	
	  keyname = "DATE-AVG";
	  BEGIN(CCCCCCCC);
	}

^DATE-OBS {
	  valtype = STRING;
	  vptr    = wcstem.dateobs;
	  if (ctrl < -10) keep = keyrec;
	
	  keyname = "DATE-OBS";
	  BEGIN(CCCCCCCC);
	}

^EPOCH{ALT}"  " {
	  sscanf(yytext, "EPOCH%c", &a);
	
	  if (relax & WCSHDR_strict) {
	    errmsg = "the EPOCH keyword is deprecated, use EQUINOXa";
	    BEGIN(ERROR);
	
	  } else if (a == ' ' || relax & WCSHDR_EPOCHa) {
	    valtype = FLOAT;
	    vptr    = &(wcstem.equinox);
	    special = wcspih_epoch;
	
	    unput(a);
	
	    keyname = "EPOCH";
	    BEGIN(CCCCCCCa);
	
	  } else if (relax & WCSHDR_reject) {
	    errmsg = "EPOCH keyword may not have an alternate version code";
	    BEGIN(ERROR);
	
	  } else {
	    BEGIN(DISCARD);
	  }
	}

^EQUINOX {
	  valtype = FLOAT;
	  vptr    = &(wcstem.equinox);
	
	  keyname = "EQUINOXa";
	  BEGIN(CCCCCCCa);
	}

^MJD-AVG" " {
	  valtype = FLOAT;
	  vptr    = &(wcstem.mjdavg);
	
	  if (ctrl < -10) keep = keyrec;
	
	  keyname = "MJD-AVG";
	  BEGIN(CCCCCCCC);
	}

^MJD-OBS" " {
	  valtype = FLOAT;
	  vptr    = &(wcstem.mjdobs);
	
	  if (ctrl < -10) keep = keyrec;
	
	  keyname = "MJD-OBS";
	  BEGIN(CCCCCCCC);
	}

^OBSGEO-X {
	  valtype = FLOAT;
	  vptr    = wcstem.obsgeo;
	
	  if (ctrl < -10) keep = keyrec;
	
	  keyname = "OBSGEO-X";
	  BEGIN(CCCCCCCC);
	}

^OBSGEO-Y {
	  valtype = FLOAT;
	  vptr    = wcstem.obsgeo + 1;
	
	  if (ctrl < -10) keep = keyrec;
	
	  keyname = "OBSGEO-Y";
	  BEGIN(CCCCCCCC);
	}

^OBSGEO-Z {
	  valtype = FLOAT;
	  vptr    = wcstem.obsgeo + 2;
	
	  if (ctrl < -10) keep = keyrec;
	
	  keyname = "OBSGEO-Z";
	  BEGIN(CCCCCCCC);
	}

^RADESYS {
	  valtype = STRING;
	  vptr    = wcstem.radesys;
	
	  keyname = "RADESYSa";
	  BEGIN(CCCCCCCa);
	}

^RADECSYS {
	  if (relax & WCSHDR_RADECSYS) {
	    valtype = STRING;
	    vptr    = wcstem.radesys;
	
	    unput(' ');
	
	    keyname = "RADECSYS";
	    BEGIN(CCCCCCCa);
	
	  } else if (relax & WCSHDR_reject) {
	    errmsg = "the RADECSYS keyword is deprecated, use RADESYSa";
	    BEGIN(ERROR);
	
	  } else {
	    BEGIN(DISCARD);
	  }
	}

^SPECSYS {
	  valtype = STRING;
	  vptr    = wcstem.specsys;
	
	  keyname = "SPECSYSa";
	  BEGIN(CCCCCCCa);
	}

^SSYSOBS {
	  valtype = STRING;
	  vptr    = wcstem.ssysobs;
	
	  keyname = "SSYSOBSa";
	  BEGIN(CCCCCCCa);
	}

^SSYSSRC {
	  valtype = STRING;
	  vptr    = wcstem.ssyssrc;
	
	  keyname = "SSYSSRCa";
	  BEGIN(CCCCCCCa);
	}

^VELANGL {
	  valtype = FLOAT;
	  vptr    = &(wcstem.velangl);
	
	  keyname = "VELANGLa";
	  BEGIN(CCCCCCCa);
	}

^VELOSYS {
	  valtype = FLOAT;
	  vptr    = &(wcstem.velosys);
	
	  keyname = "VELOSYSa";
	  BEGIN(CCCCCCCa);
	}

^VELREF{ALT}" " {
	  sscanf(yytext, "VELREF%c", &a);
	
	  if (relax & WCSHDR_strict) {
	    errmsg = "the VELREF keyword is deprecated, use SPECSYSa";
	    BEGIN(ERROR);
	
	  } else if (a == ' ' || relax & WCSHDR_VELREFa) {
	    valtype = INTEGER;
	    vptr    = &(wcstem.velref);
	
	    unput(a);
	
	    keyname = "VELREF";
	    BEGIN(CCCCCCCa);
	
	  } else if (relax & WCSHDR_reject) {
	    errmsg = "VELREF keyword may not have an alternate version code";
	    BEGIN(ERROR);
	
	  } else {
	    BEGIN(DISCARD);
	  }
	}

^VSOURCE{ALT} {
	  if (relax & WCSHDR_VSOURCE) {
	    valtype = FLOAT;
	    vptr    = &(wcstem.zsource);
	    special = wcspih_vsource;
	
	    yyless(7);
	
	    keyname = "VSOURCEa";
	    BEGIN(CCCCCCCa);
	
	  } else if (relax & WCSHDR_reject) {
	    errmsg = "the VSOURCEa keyword is deprecated, use ZSOURCEa";
	    BEGIN(ERROR);
	
	  } else {
	    BEGIN(DISCARD);
	  }
	}

^WCSNAME {
	  valtype = STRING;
	  vptr    = wcstem.wcsname;
	
	  keyname = "WCSNAMEa";
	  BEGIN(CCCCCCCa);
	}

^ZSOURCE {
	  valtype = FLOAT;
	  vptr    = &(wcstem.zsource);
	
	  keyname = "ZSOURCEa";
	  BEGIN(CCCCCCCa);
	}

^END" "{77} {
	  if (wcspih_nkeyrec) {
	    wcspih_nkeyrec = 0;
	    errmsg = "keyrecords following the END keyrecord were ignored";
	    BEGIN(ERROR);
	  } else {
	    BEGIN(DISCARD);
	  }
	}

^.	{
	  BEGIN(DISCARD);
	}

<CCCCCia>{I1}{ALT}" " |
<CCCCCia>{I2}{ALT} {
	  sscanf(yytext, "%d%c", &i, &a);
	  BEGIN(VALUE);
	}

<CCCCCia>0{I1}{ALT} |
<CCCCCia>00{I1} {
	  if (relax & WCSHDR_reject) {
	    /* Violates the basic FITS standard. */
	    errmsg = "indices in parameterized keywords must not have "
	             "leading zeroes";
	    BEGIN(ERROR);
	
	  } else {
	    /* Pretend we don't recognize it. */
	    BEGIN(DISCARD);
	  }
	}

<CCCCCia>0{ALT}" " |
<CCCCCia>00{ALT} |
<CCCCCia>{Z3} {
	  /* Anything that has fallen through to this point must contain */
	  /* an invalid axis number. */
	  errmsg = "axis number must exceed 0";
	  BEGIN(ERROR);
	}

<CCCCCia>. {
	  if (relax & WCSHDR_reject) {
	    /* Looks too much like a FITS WCS keyword not to flag it. */
	    errmsg = errtxt;
	    sprintf(errmsg, "keyword looks very much like %s but isn't",
	      keyname);
	    BEGIN(ERROR);
	
	  } else {
	    /* Let it go. */
	    BEGIN(DISCARD);
	  }
	}

<CCi_ja>{I1}_{I1}{ALT}"  " |
<CCi_ja>{I1}_{I2}{ALT}" " |
<CCi_ja>{I2}_{I1}{ALT}" " |
<CCi_ja>{I2}_{I2}{ALT} {
	  sscanf(yytext, "%d_%d%c", &i, &j, &a);
	  BEGIN(VALUE);
	}


<CCi_ja>0{I1}_{I1}{ALT}" " |
<CCi_ja>{I1}_0{I1}{ALT}" " |
<CCi_ja>00{I1}_{I1}{ALT} |
<CCi_ja>0{I1}_0{I1}{ALT} |
<CCi_ja>{I1}_00{I1}{ALT} |
<CCi_ja>000{I1}_{I1} |
<CCi_ja>00{I1}_0{I1} |
<CCi_ja>0{I1}_00{I1} |
<CCi_ja>{I1}_000{I1} |
<CCi_ja>0{I1}_{I2}{ALT} |
<CCi_ja>{I1}_0{I2}{ALT} |
<CCi_ja>00{I1}_{I2} |
<CCi_ja>0{I1}_0{I2} |
<CCi_ja>{I1}_00{I2} |
<CCi_ja>0{I2}_{I1}{ALT} |
<CCi_ja>{I2}_0{I1}{ALT} |
<CCi_ja>00{I2}_{I1} |
<CCi_ja>0{I2}_0{I1} |
<CCi_ja>{I2}_00{I1} |
<CCi_ja>0{I2}_{I2} |
<CCi_ja>{I2}_0{I2} {
	  if (((altlin == 1) && (relax & WCSHDR_PC0i_0ja)) ||
	      ((altlin == 2) && (relax & WCSHDR_CD0i_0ja))) {
	    sscanf(yytext, "%d_%d%c", &i, &j, &a);
	    BEGIN(VALUE);
	
	  } else if (relax & WCSHDR_reject) {
	    errmsg = "indices in parameterized keywords must not have "
	             "leading zeroes";
	    BEGIN(ERROR);
	
	  } else {
	    /* Pretend we don't recognize it. */
	    BEGIN(DISCARD);
	  }
	}

<CCi_ja>{Z1}_{Z1}{ALT}"  " |
<CCi_ja>{Z2}_{Z1}{ALT}" " |
<CCi_ja>{Z1}_{Z2}{ALT}" " |
<CCi_ja>{Z3}_{Z1}{ALT} |
<CCi_ja>{Z2}_{Z2}{ALT} |
<CCi_ja>{Z1}_{Z3}{ALT} |
<CCi_ja>{Z4}_{Z1} |
<CCi_ja>{Z3}_{Z2} |
<CCi_ja>{Z2}_{Z3} |
<CCi_ja>{Z1}_{Z4} {
	  /* Anything that has fallen through to this point must contain */
	  /* an invalid axis number. */
	  errmsg = "axis number must exceed 0";
	  BEGIN(ERROR);
	}

<CCi_ja>{Z1}-{Z1}{ALT}"  " |
<CCi_ja>{Z2}-{Z1}{ALT}" " |
<CCi_ja>{Z1}-{Z2}{ALT}" " |
<CCi_ja>{Z3}-{Z1}{ALT} |
<CCi_ja>{Z2}-{Z2}{ALT} |
<CCi_ja>{Z1}-{Z3}{ALT} |
<CCi_ja>{Z4}-{Z1} |
<CCi_ja>{Z3}-{Z2} |
<CCi_ja>{Z2}-{Z3} |
<CCi_ja>{Z1}-{Z4} {
	  errmsg = errtxt;
	  sprintf(errmsg, "%s keyword must use an underscore, not a dash",
	    keyname);
	  BEGIN(ERROR);
	}

<CCi_ja>{Z1}{6} {
	  /* This covers the defunct forms CD00i00j and PC00i00j. */
	  if (((altlin == 1) && (relax & WCSHDR_PC00i00j)) ||
	      ((altlin == 2) && (relax & WCSHDR_CD00i00j))) {
	    sscanf(yytext, "%3d%3d", &i, &j);
	    a = ' ';
	    BEGIN(VALUE);
	
	  } else if (relax & WCSHDR_reject) {
	    errmsg = errtxt;
	    sprintf(errmsg,
	      "this form of the %s keyword is deprecated, use %s",
	      keyname, keyname);
	    BEGIN(ERROR);
	
	  } else {
	    /* Pretend we don't recognize it. */
	    BEGIN(DISCARD);
	  }
	}

<CCi_ja>. {
	  BEGIN(DISCARD);
	}

<CROTAi>{Z1}{ALT}" " |
<CROTAi>{Z2}{ALT} |
<CROTAi>{Z3} {
	  a = ' ';
	  sscanf(yytext, "%d%c", &i, &a);
	
	  if (relax & WCSHDR_strict) {
	    errmsg = "the CROTAn keyword is deprecated, use PCi_ja";
	    BEGIN(ERROR);
	
	  } else if (a == ' ' || relax & WCSHDR_CROTAia) {
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
	  yyless(0);
	  BEGIN(CCCCCia);
	}

<CCCCCCCa>{ALT} |
<CCCCCCCC>. {
	  if (YY_START == CCCCCCCa) {
	    sscanf(yytext, "%c", &a);
	  } else {
	    unput(yytext[0]);
	    a = 0;
	  }
	
	  BEGIN(VALUE);
	}

<CCCCCCCa>. {
	  if (relax & WCSHDR_reject) {
	    /* Looks too much like a FITS WCS keyword not to flag it. */
	    errmsg = errtxt;
	    sprintf(errmsg, "invalid alternate code, keyword resembles %s "
	      "but isn't", keyname);
	    BEGIN(ERROR);
	
	  } else {
	    /* Pretend we don't recognize it. */
	    BEGIN(DISCARD);
	  }
	}

<CCi_ma>{I1}_{Z1}{ALT}"  " |
<CCi_ma>{I1}_{I2}{ALT}" " |
<CCi_ma>{I2}_{Z1}{ALT}" " |
<CCi_ma>{I2}_{I2}{ALT} {
	  sscanf(yytext, "%d_%d%c", &i, &m, &a);
	  BEGIN(VALUE);
	}

<CCi_ma>0{I1}_{Z1}{ALT}" " |
<CCi_ma>{I1}_0{Z1}{ALT}" " |
<CCi_ma>00{I1}_{Z1}{ALT} |
<CCi_ma>0{I1}_0{Z1}{ALT} |
<CCi_ma>{I1}_00{Z1}{ALT} |
<CCi_ma>000{I1}_{Z1} |
<CCi_ma>00{I1}_0{Z1} |
<CCi_ma>0{I1}_00{Z1} |
<CCi_ma>{I1}_000{Z1} |
<CCi_ma>0{I1}_{I2}{ALT} |
<CCi_ma>{I1}_0{I2}{ALT} |
<CCi_ma>00{I1}_{I2} |
<CCi_ma>0{I1}_0{I2} |
<CCi_ma>{I1}_00{I2} |
<CCi_ma>0{I2}_{Z1}{ALT} |
<CCi_ma>{I2}_0{Z1}{ALT} |
<CCi_ma>00{I2}_{Z1} |
<CCi_ma>0{I2}_0{Z1} |
<CCi_ma>{I2}_00{Z1} |
<CCi_ma>0{I2}_{I2} |
<CCi_ma>{I2}_0{I2} {
	  if (((valtype == FLOAT)  && (relax & WCSHDR_PV0i_0ma)) ||
	      ((valtype == STRING) && (relax & WCSHDR_PS0i_0ma))) {
	    sscanf(yytext, "%d_%d%c", &i, &m, &a);
	    BEGIN(VALUE);
	
	  } else if (relax & WCSHDR_reject) {
	    errmsg = "indices in parameterized keywords must not have "
	             "leading zeroes";
	    BEGIN(ERROR);
	
	  } else {
	    /* Pretend we don't recognize it. */
	    BEGIN(DISCARD);
	  }
	}

<CCi_ma>{Z1}_{Z1}{ALT}"  " |
<CCi_ma>{Z2}_{Z1}{ALT}" " |
<CCi_ma>{Z1}_{Z2}{ALT}" " |
<CCi_ma>{Z3}_{Z1}{ALT} |
<CCi_ma>{Z2}_{Z2}{ALT} |
<CCi_ma>{Z1}_{Z3}{ALT} |
<CCi_ma>{Z4}_{Z1} |
<CCi_ma>{Z3}_{Z2} |
<CCi_ma>{Z2}_{Z3} |
<CCi_ma>{Z1}_{Z4} {
	  /* Anything that has fallen through to this point must contain */
	  /* an invalid axis number. */
	  errmsg = "axis number must exceed 0";
	  BEGIN(ERROR);
	}

<CCi_ma>{Z1}-{Z1}{ALT}"  " |
<CCi_ma>{Z2}-{Z1}{ALT}" " |
<CCi_ma>{Z1}-{Z2}{ALT}" " |
<CCi_ma>{Z3}-{Z1}{ALT} |
<CCi_ma>{Z2}-{Z2}{ALT} |
<CCi_ma>{Z1}-{Z3}{ALT} |
<CCi_ma>{Z4}-{Z1} |
<CCi_ma>{Z3}-{Z2} |
<CCi_ma>{Z2}-{Z3} |
<CCi_ma>{Z1}-{Z4} {
	  errmsg = errtxt;
	  sprintf(errmsg, "%s keyword must use an underscore, not a dash",
	    keyname);
	  BEGIN(ERROR);
	}

<CCi_ma>. {
	  BEGIN(DISCARD);
	}

<PROJPn>{Z1}"  " {
	  if (relax & WCSHDR_PROJPn) {
	    sscanf(yytext, "%d", &m);
	    i = 0;
	    a = ' ';
	    BEGIN(VALUE);
	
	  } else if (relax & WCSHDR_reject) {
	    errmsg = "the PROJPn keyword is deprecated, use PVi_ma";
	    BEGIN(ERROR);
	
	  } else {
	    BEGIN(DISCARD);
	  }
	}

<PROJPn>{Z2}" " |
<PROJPn>{Z3} {
	  if (relax & (WCSHDR_PROJPn | WCSHDR_reject)) {
	    errmsg = "invalid PROJPn keyword";
	    BEGIN(ERROR);
	
	  } else {
	    BEGIN(DISCARD);
	  }
	}

<PROJPn>. {
	  BEGIN(DISCARD);
	}

<VALUE>=" "+ {
	  /* Do checks on i, j & m. */
	  if (99 < i || 99 < j || 99 < m) {
	    if (relax & WCSHDR_reject) {
	      if (99 < i || 99 < j) {
	        errmsg = "axis number exceeds 99";
	      } else if (m > 99) {
	        errmsg = "parameter number exceeds 99";
	      }
	      BEGIN(ERROR);
	
	    } else {
	      /* Pretend we don't recognize it. */
	      BEGIN(DISCARD);
	    }
	
	  } else {
	    if (valtype == INTEGER) {
	      BEGIN(INTEGER_VAL);
	    } else if (valtype == FLOAT) {
	      BEGIN(FLOAT_VAL);
	    } else if (valtype == STRING) {
	      BEGIN(STRING_VAL);
	    } else {
	      errmsg = errtxt;
	      sprintf(errmsg, "internal parser ERROR, bad data type: %d",
	        valtype);
	      BEGIN(ERROR);
	    }
	  }
	}

<VALUE>. {
	  errmsg = "invalid KEYWORD = VALUE syntax";
	  BEGIN(ERROR);
	}

<INTEGER_VAL>{INTEGER} {
	  if (ipass == 1) {
	    BEGIN(COMMENT);
	
	  } else {
	    /* Read the keyvalue. */
	    sscanf(yytext, "%d", &inttmp);
	
	    BEGIN(COMMENT);
	  }
	}

<INTEGER_VAL>. {
	  errmsg = "an integer value was expected";
	  BEGIN(ERROR);
	}

<FLOAT_VAL>{FLOAT} {
	  if (ipass == 1) {
	    BEGIN(COMMENT);
	
	  } else {
	    /* Read the keyvalue. */
	    wcsutil_str2double(yytext, "%lf", &dbltmp);
	
	    BEGIN(COMMENT);
	  }
	}

<FLOAT_VAL>. {
	  errmsg = "a floating-point value was expected";
	  BEGIN(ERROR);
	}

<STRING_VAL>{STRING} {
	  if (ipass == 1) {
	    BEGIN(COMMENT);
	
	  } else {
	    /* Read the keyvalue. */
	    strcpy(strtmp, yytext+1);
	
	    /* Squeeze out repeated quotes. */
	    ix = 0;
	    for (jx = 0; jx < 72; jx++) {
	      if (ix < jx) {
	        strtmp[ix] = strtmp[jx];
	      }
	
	      if (strtmp[jx] == '\0') {
	        if (ix) strtmp[ix-1] = '\0';
	        break;
	      } else if (strtmp[jx] == '\'' && strtmp[jx+1] == '\'') {
	        jx++;
	      }
	
	      ix++;
	    }
	
	    BEGIN(COMMENT);
	  }
	}

<STRING_VAL>. {
	  errmsg = "a string value was expected";
	  BEGIN(ERROR);
	}

<COMMENT>{INLINE}$ {
	  if (ipass == 1) {
	    /* Do first-pass bookkeeping. */
	    wcspih_naxes(naxis, i, j, a, alts, npptr);
	    BEGIN(FLUSH);
	
	  } else if (*wcs) {
	    /* Store the value now that the keyrecord has been validated. */
	    gotone = 0;
	    for (ialt = 0; ialt < *nwcs; ialt++) {
	      /* The loop here is for keywords that apply */
	      /* to every alternate; these have a == 0. */
	      if (a >= 'A') {
	        ialt = alts[a-'A'+1];
	        if (ialt < 0) break;
	      }
	      gotone = 1;
	
	      if (vptr) {
	        wcsp = *wcs + ialt;
	        voff = (char *)vptr - (char *)(&wcstem);
	        wptr = (void *)((char *)wcsp + voff);
	
	        if (valtype == INTEGER) {
	          *((int *)wptr) = inttmp;
	
	        } else if (valtype == FLOAT) {
	          /* Apply keyword parameterization. */
	          if (npptr == npv) {
	            ipx = (wcsp->npv)++;
	            wcsp->pv[ipx].i = i;
	            wcsp->pv[ipx].m = m;
	            wptr = &(wcsp->pv[ipx].value);
	
	          } else if (j) {
	            wptr = *((double **)wptr) + (i - 1)*(wcsp->naxis)
	                                      + (j - 1);
	
	          } else if (i) {
	            wptr = *((double **)wptr) + (i - 1);
	          }
	
	          if (special) {
	            special(wptr, dbltmp);
	          } else {
	            *((double *)wptr) = dbltmp;
	          }
	
	          /* Flag the presence of PCi_ja, or CDi_ja and/or CROTAia. */
	          if (altlin) {
	            wcsp->altlin |= altlin;
	            altlin = 0;
	          }
	
	        } else if (valtype == STRING) {
	          /* Apply keyword parameterization. */
	          if (npptr == nps) {
	            ipx = (wcsp->nps)++;
	            wcsp->ps[ipx].i = i;
	            wcsp->ps[ipx].m = m;
	            wptr = wcsp->ps[ipx].value;
	
	          } else if (j) {
	            wptr = *((char (**)[72])wptr) +
	                    (i - 1)*(wcsp->naxis) + (j - 1);
	
	          } else if (i) {
	            wptr = *((char (**)[72])wptr) + (i - 1);
	          }
	
	          cptr = (char *)wptr;
	          strcpy(cptr, strtmp);
	        }
	      }
	
	      if (a) break;
	    }
	
	    if (gotone) {
	      nvalid++;
	      if (ctrl == 4) {
	        wcsfprintf(stderr,
	          "%.80s\n  Accepted (%d) as a valid WCS keyrecord.\n",
	          keyrec, nvalid);
	      }
	
	      BEGIN(FLUSH);
	
	    } else {
	      errmsg = "syntactically valid WCS keyrecord has no effect";
	      BEGIN(ERROR);
	    }
	
	  } else {
	    BEGIN(FLUSH);
	  }
	}

<COMMENT>.*" "*\/.*$ {
	  errmsg = "invalid keyvalue";
	  BEGIN(ERROR);
	}

<COMMENT>[^ \/\n]*{INLINE}$ {
	  errmsg = "invalid keyvalue";
	  BEGIN(ERROR);
	}

<COMMENT>" "+[^\/\n].*{INLINE}$ {
	  errmsg = "invalid keyvalue or malformed keycomment";
	  BEGIN(ERROR);
	}

<COMMENT>.*$ {
	  errmsg = "malformed keycomment";
	  BEGIN(ERROR);
	}

<DISCARD>.*$ {
	  if (ipass == npass) {
	    if (ctrl < 0) {
	      /* Preserve discards. */
	      keep = keyrec;
	
	    } else if (2 < ctrl) {
	      nother++;
	      wcsfprintf(stderr, "%.80s\n  Not a recognized WCS keyword.\n",
	        keyrec);
	    }
	  }
	  BEGIN(FLUSH);
	}

<ERROR>.*$ {
	  if (ipass == npass) {
	    (*nreject)++;
	
	    if (ctrl%10 == -1) {
	      /* Preserve rejects. */
	      keep = keyrec;
	    }
	
	    if (1 < abs(ctrl%10)) {
	      wcsfprintf(stderr, "%.80s\n  Rejected (%d), %s.\n",
	        keyrec, *nreject, errmsg);
	    }
	  }
	  BEGIN(FLUSH);
	}

<FLUSH>.*\n {
	  if (ipass == npass && keep) {
	    if (hptr < keep) {
	      strncpy(hptr, keep, 80);
	    }
	    hptr += 80;
	  }
	
	  /* Throw away the rest of the line and reset for the next one. */
	  i = j = 0;
	  m = 0;
	  a = ' ';
	
	  keyrec += 80;
	
	  valtype = -1;
	  vptr    = 0x0;
	  keep    = 0x0;
	
	  altlin  = 0;
	  npptr   = 0x0;
	  special = 0x0;
	
	  BEGIN(INITIAL);
	}

<<EOF>>	 {
	  /* End-of-input. */
	  if (ipass == 1) {
	    if ((status = wcspih_inits(naxis, alts, npv, nps, nwcs, wcs)) ||
	        (*nwcs == 0 && ctrl == 0)) {
	      yylex_destroy();
	      return status;
	    }
	
	    if (2 < abs(ctrl%10)) {
	      if (*nwcs == 1) {
	        if (strcmp(wcs[0]->wcsname, "DEFAULTS") != 0) {
	          wcsfprintf(stderr, "Found one coordinate representation.\n");
	        }
	      } else {
	        wcsfprintf(stderr, "Found %d coordinate representations.\n",
	          *nwcs);
	      }
	    }
	  }
	
	  if (ipass++ < npass) {
	    wcspih_hdr = header;
	    wcspih_nkeyrec = nkeyrec;
	    keyrec = header;
	    *nreject = 0;
	
	    i = j = 0;
	    m = 0;
	    a = ' ';
	
	    valtype = -1;
	    vptr    = 0x0;
	
	    altlin  = 0;
	    npptr   = 0x0;
	    special = 0x0;
	
	    yyrestart(yyin);
	
	  } else {
	    yylex_destroy();
	
	    if (ctrl < 0) {
	      *hptr = '\0';
	    } else if (ctrl == 1) {
	      wcsfprintf(stderr, "%d WCS keyrecord%s rejected.\n",
	        *nreject, (*nreject==1)?" was":"s were");
	    } else if (ctrl == 4) {
	      wcsfprintf(stderr, "\n");
	      wcsfprintf(stderr, "%5d keyrecord%s rejected for syntax or "
	        "other errors,\n", *nreject, (*nreject==1)?" was":"s were");
	      wcsfprintf(stderr, "%5d %s recognized as syntactically valid, "
	        "and\n", nvalid, (nvalid==1)?"was":"were");
	      wcsfprintf(stderr, "%5d other%s were not recognized as WCS "
	        "keyrecords.\n", nother, (nother==1)?"":"s");
	    }
	
	    return wcspih_final(alts, nwcs, wcs);
	  }
	}

%%

/*----------------------------------------------------------------------------
* Determine the number of coordinate representations (up to 27) and the
* number of coordinate axes in each, and count the number of PVi_ma and
* PSi_ma keywords in each representation.
*---------------------------------------------------------------------------*/

void wcspih_naxes(int naxis, int i, int j, char a, int alts[], int *npptr)

{
  /* On the first pass alts[] is used to determine the number of axes */
  /* for each of the 27 possible alternate coordinate descriptions.   */
  int ialt, *ip;

  if (a == 0) {
    return;
  }

  ialt = 0;
  if (a != ' ') {
    ialt = a - 'A' + 1;
  }

  ip = alts + ialt;

  if (*ip < naxis) {
    *ip = naxis;
  }

  /* i or j can be greater than naxis. */
  if (*ip < i) {
    *ip = i;
  }

  if (*ip < j) {
    *ip = j;
  }

  if (npptr) {
    npptr[ialt]++;
  }
}


/*----------------------------------------------------------------------------
* Allocate memory for an array of the required number of wcsprm structs and
* initialize each of them.
*---------------------------------------------------------------------------*/

int wcspih_inits(
  int naxis,
  int alts[],
  int npv[],
  int nps[],
  int *nwcs,
  struct wcsprm **wcs)

{
  int ialt, defaults, npsmax, npvmax, status = 0;
  struct wcsprm *wcsp;

  /* Find the number of coordinate descriptions. */
  *nwcs = 0;
  for (ialt = 0; ialt < 27; ialt++) {
    if (alts[ialt]) (*nwcs)++;
  }

  if ((defaults = !(*nwcs) && naxis)) {
    /* NAXIS is non-zero but there were no WCS keywords with an alternate
       version code; create a default WCS with blank alternate version. */
    wcspih_naxes(naxis, 0, 0, ' ', alts, 0x0);
    *nwcs = 1;
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
    for (ialt = 0; ialt < 27; ialt++) {
      if (alts[ialt]) {
        wcsp->flag = -1;
        wcsnpv(npv[ialt]);
        wcsnps(nps[ialt]);
        if ((status = wcsini(1, alts[ialt], wcsp))) {
          wcsvfree(nwcs, wcs);
          break;
        }

        /* Record the alternate version code. */
        if (ialt) {
          wcsp->alt[0] = 'A' + ialt - 1;
        }

        /* Record in wcsname whether this is a default description. */
        if (defaults) {
          strcpy(wcsp->wcsname, "DEFAULTS");
        }

        /* On the second pass alts[] indexes the array of wcsprm structs. */
        alts[ialt] = (*nwcs)++;

        wcsp++;

      } else {
        /* Signal that there is no wcsprm for this alt. */
        alts[ialt] = -1;
      }
    }

    /* Restore the original values of NPVMAX and NPSMAX. */
    wcsnpv(npvmax);
    wcsnps(npsmax);
  }

  return status;
}


/*----------------------------------------------------------------------------
* Interpret EPOCHa keywords.
*---------------------------------------------------------------------------*/

int wcspih_epoch(double *equinox, double epoch)

{
  /* If EQUINOXa is currently undefined then set it from EPOCHa. */
  if (undefined(*equinox)) {
    *equinox = epoch;
  }

  return 0;
}


/*----------------------------------------------------------------------------
* Interpret VSOURCEa keywords.
*---------------------------------------------------------------------------*/

int wcspih_vsource(double *zsource, double vsource)

{
  double beta, c = 299792458.0;

  /* If ZSOURCEa is currently undefined then set it from VSOURCEa. */
  if (undefined(*zsource)) {
    /* Convert relativistic Doppler velocity to redshift. */
    beta = vsource/c;
    *zsource = (1.0 + beta)/sqrt(1.0 - beta*beta) - 1.0;
  }

  return 0;
}


/*----------------------------------------------------------------------------
* Interpret special keywords encountered for each coordinate representation.
*---------------------------------------------------------------------------*/

int wcspih_final(
  int alts[],
  int *nwcs,
  struct wcsprm **wcs)

{
  int ialt, status;

  for (ialt = 0; ialt < *nwcs; ialt++) {
    /* Interpret -TAB header keywords. */
    if ((status = wcstab(*wcs+ialt))) {
       wcsvfree(nwcs, wcs);
       return status;
    }
  }

  return 0;
}
