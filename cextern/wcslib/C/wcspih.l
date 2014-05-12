/*============================================================================

  WCSLIB 4.23 - an implementation of the FITS WCS standard.
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
  $Id: wcspih.l,v 4.23 2014/05/11 04:09:38 mcalabre Exp $
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

int wcspih_final(int alts[], double epoch[], double vsource[], int *nwcs,
        struct wcsprm **wcs);
int wcspih_inits(int naxis, int alts[], int npv[], int nps[], int *nwcs,
        struct wcsprm **wcs);
void wcspih_naxes(int naxis, int i, int j, char a, int alts[], int *npptr);

/* Used in preempting the call to exit() by yy_fatal_error(). */
jmp_buf wcspih_abort_jmp_env;
#define exit(status) longjmp(wcspih_abort_jmp_env, status)

%}

%%
	/* Keyword indices, as used in the WCS papers, e.g. PCi_ja, PVi_ma. */
	char a;
	int  i, j, m;
	
	char *cptr, *errmsg, errtxt[80], *hptr, *keep;
	int  altlin, alts[27], ialt, idx, ipx, ix, jx, naxis, *npptr,
	     nps[27], npv[27], pass, status, valtype, voff;
	double epoch[27], vsource[27];
	void *vptr, *wptr;
	struct wcsprm *wcsp;
	int yylex_destroy(void);
	
	naxis = 0;
	for (ialt = 0; ialt < 27; ialt++) {
	  alts[ialt] = 0;
	  npv[ialt] = 0;
	  nps[ialt] = 0;
	  epoch[ialt]   = UNDEFINED;
	  vsource[ialt] = UNDEFINED;
	}
	
	/* Parameters used to implement YY_INPUT. */
	wcspih_hdr = header;
	wcspih_nkeyrec = nkeyrec;
	
	/* Our handle on the input stream. */
	hptr = header;
	keep = 0x0;
	*nreject = 0;
	
	/* Keyword parameters. */
	i = j = m = 0;
	a = ' ';
	
	/* For decoding the keyvalue. */
	valtype = -1;
	idx     = -1;
	vptr    = 0x0;
	
	/* For keywords that require special handling. */
	altlin = 0;
	npptr  = 0x0;
	
	/* The data structures produced. */
	*nwcs = 0;
	*wcs  = 0x0;
	
	pass = 1;
	
	/* Return here via longjmp() invoked by yy_fatal_error(). */
	if (setjmp(wcspih_abort_jmp_env)) {
	  return 3;
	}
	
	BEGIN(INITIAL);


^NAXIS"   = "" "*{INTEGER} {
	  if (pass == 1) {
	    sscanf(yytext, "NAXIS   = %d", &naxis);
	  }
	
	  if (naxis < 0) {
	    errmsg = errtxt;
	    sprintf(errmsg, "Negative value of NAXIS ignored: %d", naxis);
	    naxis = 0;
	    BEGIN(ERROR);
	  } else {
	    BEGIN(DISCARD);
	  }
	}

^WCSAXES{ALT}=" "" "*{INTEGER} {
	  if (pass == 1) {
	    sscanf(yytext, "WCSAXES%c= %d", &a, &i);
	    wcspih_naxes(naxis, i, 0, a, alts, 0);
	  }
	  BEGIN(FLUSH);
	}

^CRPIX	{
	  valtype = FLOAT;
	  if (pass == 2) vptr = &((*wcs)->crpix);
	  BEGIN(CCCCCia);
	}

^PC	{
	  valtype = FLOAT;
	  if (pass == 2) vptr = &((*wcs)->pc);
	  altlin = 1;
	  BEGIN(CCi_ja);
	}

^CD	{
	  valtype = FLOAT;
	  if (pass == 2) vptr = &((*wcs)->cd);
	  altlin = 2;
	  BEGIN(CCi_ja);
	}

^CDELT	{
	  valtype = FLOAT;
	  if (pass == 2) vptr = &((*wcs)->cdelt);
	  BEGIN(CCCCCia);
	}

^CROTA	{
	  valtype = FLOAT;
	  if (pass == 2) vptr = &((*wcs)->crota);
	  altlin = 4;
	  BEGIN(CROTAi);
	}

^CUNIT	{
	  valtype = STRING;
	  if (pass == 2) vptr = &((*wcs)->cunit);
	  BEGIN(CCCCCia);
	}

^CTYPE	{
	  valtype = STRING;
	  if (pass == 2) vptr = &((*wcs)->ctype);
	  BEGIN(CCCCCia);
	}

^CRVAL	{
	  valtype = FLOAT;
	  if (pass == 2) vptr = &((*wcs)->crval);
	  BEGIN(CCCCCia);
	}

^LONPOLE {
	  valtype = FLOAT;
	  if (pass == 2) vptr = &((*wcs)->lonpole);
	  BEGIN(CCCCCCCa);
	}

^LATPOLE {
	  valtype = FLOAT;
	  if (pass == 2) vptr = &((*wcs)->latpole);
	  BEGIN(CCCCCCCa);
	}

^RESTFRQ {
	  valtype = FLOAT;
	  if (pass == 2) vptr = &((*wcs)->restfrq);
	  BEGIN(CCCCCCCa);
	}

^RESTFREQ {
	  valtype = FLOAT;
	  if (pass == 2) vptr = &((*wcs)->restfrq);
	  unput(' ');
	  BEGIN(CCCCCCCa);
	}

^RESTWAV {
	  valtype = FLOAT;
	  if (pass == 2) vptr = &((*wcs)->restwav);
	  BEGIN(CCCCCCCa);
	}

^PV	{
	  valtype = FLOAT;
	  if (pass == 2) vptr = &((*wcs)->pv);
	  npptr = npv;
	  BEGIN(CCi_ma);
	}

^PROJP	{
	  valtype = FLOAT;
	  if (pass == 2) vptr = &((*wcs)->pv);
	  npptr = npv;
	  BEGIN(PROJPn);
	}

^PS	{
	  valtype = STRING;
	  if (pass == 2) vptr = &((*wcs)->ps);
	  npptr = nps;
	  BEGIN(CCi_ma);
	}

^CNAME	{
	  valtype = STRING;
	  if (pass == 2) vptr = &((*wcs)->cname);
	  BEGIN(CCCCCia);
	}

^CRDER	{
	  valtype = FLOAT;
	  if (pass == 2) vptr = &((*wcs)->crder);
	  BEGIN(CCCCCia);
	}

^CSYER	{
	  valtype = FLOAT;
	  if (pass == 2) vptr = &((*wcs)->csyer);
	  BEGIN(CCCCCia);
	}

^DATE-AVG {
	  valtype = STRING;
	  if (pass == 2) vptr = (*wcs)->dateavg;
	  if (ctrl < -10) keep = wcspih_hdr - 80;
	  BEGIN(CCCCCCCC);
	}

^DATE-OBS {
	  valtype = STRING;
	  if (pass == 2) vptr = (*wcs)->dateobs;
	  if (ctrl < -10) keep = wcspih_hdr - 80;
	  BEGIN(CCCCCCCC);
	}

^EPOCH{ALT}"  " {
	  sscanf(yytext, "EPOCH%c", &a);
	
	  if (a == ' ' || relax & WCSHDR_EPOCHa) {
	    valtype = FLOAT;
	    if (pass == 2) {
	      vptr = epoch;
	      if (a >= 'A') {
	        vptr = (void *)((double *)vptr + alts[a-'A'+1]);
	      }
	    }
	
	    unput(' ');
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
	  if (pass == 2) vptr = &((*wcs)->equinox);
	  BEGIN(CCCCCCCa);
	}

^MJD-AVG" " {
	  valtype = FLOAT;
	  if (pass == 2) vptr = &((*wcs)->mjdavg);
	  if (ctrl < -10) keep = wcspih_hdr - 80;
	  BEGIN(CCCCCCCC);
	}

^MJD-OBS" " {
	  valtype = FLOAT;
	  if (pass == 2) vptr = &((*wcs)->mjdobs);
	  if (ctrl < -10) keep = wcspih_hdr - 80;
	  BEGIN(CCCCCCCC);
	}

^OBSGEO-X {
	  valtype = FLOAT;
	  if (pass == 2) vptr = (*wcs)->obsgeo;
	  if (ctrl < -10) keep = wcspih_hdr - 80;
	  BEGIN(CCCCCCCC);
	}

^OBSGEO-Y {
	  valtype = FLOAT;
	  if (pass == 2) vptr = (*wcs)->obsgeo + 1;
	  if (ctrl < -10) keep = wcspih_hdr - 80;
	  BEGIN(CCCCCCCC);
	}

^OBSGEO-Z {
	  valtype = FLOAT;
	  if (pass == 2) vptr = (*wcs)->obsgeo + 2;
	  if (ctrl < -10) keep = wcspih_hdr - 80;
	  BEGIN(CCCCCCCC);
	}

^RADESYS {
	  valtype = STRING;
	  if (pass == 2) vptr = (*wcs)->radesys;
	  BEGIN(CCCCCCCa);
	}

^RADECSYS {
	  if (relax & WCSHDR_RADECSYS) {
	    valtype = STRING;
	    if (pass == 2) vptr = (*wcs)->radesys;
	    unput(' ');
	    BEGIN(CCCCCCCa);
	
	  } else if (relax & WCSHDR_reject) {
	    errmsg = "RADECSYS is non-standard, use RADESYSa";
	    BEGIN(ERROR);
	
	  } else {
	    BEGIN(DISCARD);
	  }
	}

^SPECSYS {
	  valtype = STRING;
	  if (pass == 2) vptr = (*wcs)->specsys;
	  BEGIN(CCCCCCCa);
	}

^SSYSOBS {
	  valtype = STRING;
	  if (pass == 2) vptr = (*wcs)->ssysobs;
	  BEGIN(CCCCCCCa);
	}

^SSYSSRC {
	  valtype = STRING;
	  if (pass == 2) vptr = (*wcs)->ssyssrc;
	  BEGIN(CCCCCCCa);
	}

^VELANGL {
	  valtype = FLOAT;
	  if (pass == 2) vptr = &((*wcs)->velangl);
	  BEGIN(CCCCCCCa);
	}

^VELOSYS {
	  valtype = FLOAT;
	  if (pass == 2) vptr = &((*wcs)->velosys);
	  BEGIN(CCCCCCCa);
	}

^VELREF{ALT}" " {
	  sscanf(yytext, "VELREF%c", &a);
	
	  if (a == ' ' || relax & WCSHDR_VELREFa) {
	    valtype = INTEGER;
	    if (pass == 2) vptr = &((*wcs)->velref);
	
	    unput(a);
	    BEGIN(CCCCCCCa);
	
	  } else if (relax & WCSHDR_reject) {
	    errmsg = "VELREF keyword may not have an alternate version code";
	    BEGIN(ERROR);
	
	  } else {
	    BEGIN(DISCARD);
	  }
	}

^VSOURCE{ALT} {
	  sscanf(yytext, "VSOURCE%c", &a);
	
	  if (relax & WCSHDR_VSOURCE) {
	    valtype = FLOAT;
	    if (pass == 2) {
	      vptr = vsource;
	      if (a >= 'A') {
	        vptr = (void *)((double *)vptr + alts[a-'A'+1]);
	      }
	    }
	
	    unput(' ');
	    BEGIN(CCCCCCCa);
	
	  } else if (relax & WCSHDR_reject) {
	    errmsg = "Deprecated VSOURCEa keyword rejected";
	    BEGIN(ERROR);
	
	  } else {
	    BEGIN(DISCARD);
	  }
	}

^WCSNAME {
	  valtype = STRING;
	  if (pass == 2) vptr = (*wcs)->wcsname;
	  BEGIN(CCCCCCCa);
	}

^ZSOURCE {
	  valtype = FLOAT;
	  if (pass == 2) vptr = &((*wcs)->zsource);
	  BEGIN(CCCCCCCa);
	}

^END" "{77} {
	  yyless(0);
	  if (wcspih_nkeyrec) {
	    wcspih_nkeyrec = 0;
	    errmsg = "Keyrecords following the END keyrecord were ignored";
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
	  idx = i - 1;
	  BEGIN(VALUE);
	}

<CCCCCia>{I3} {
	  /* Invalid axis number will be caught by <VALUE>. */
	  sscanf(yytext, "%3d", &i);
	  BEGIN(VALUE);
	}

<CCCCCia>. {
	  BEGIN(DISCARD);
	}

<CCi_ja>{I1}_{I1}{ALT}"  " |
<CCi_ja>{I1}_{I2}{ALT}" " |
<CCi_ja>{I2}_{I1}{ALT}" " |
<CCi_ja>{I2}_{I2}{ALT} {
	  sscanf(yytext, "%d_%d%c", &i, &j, &a);
	  if (pass == 2) {
	    wcsp = *wcs;
	    if (a != ' ') {
	      wcsp += alts[a-'A'+1];
	    }
	
	    idx = (i-1)*(wcsp->naxis) + j - 1;
	  }
	  BEGIN(VALUE);
	}

<CCi_ja>{I1}_{I3}{ALT} |
<CCi_ja>{I3}_{I1}{ALT} |
<CCi_ja>{I1}_{I4} |
<CCi_ja>{I2}_{I3} |
<CCi_ja>{I3}_{I2} |
<CCi_ja>{I4}_{I1} {
	  /* Invalid axis numbers will be caught by <VALUE>. */
	  sscanf(yytext, "%d_%d", &i, &j);
	  BEGIN(VALUE);
	}

<CCi_ja>{I0}{6} {
	  /* This covers the defunct forms CD00i00j and PC00i00j. */
	  if (((relax & WCSHDR_PC00i00j) && (altlin == 1)) ||
	      ((relax & WCSHDR_CD00i00j) && (altlin == 2))) {
	    sscanf(yytext, "%3d%3d", &i, &j);
	    a = ' ';
	    if (pass == 2) {
	      idx = (i-1)*((*wcs)->naxis) + j - 1;
	    }
	    BEGIN(VALUE);
	
	  } else if (relax & WCSHDR_reject) {
	    errmsg = errtxt;
	    sprintf(errmsg, "Defunct form of %si_ja keyword",
	                     (altlin==1) ? "PC" : "CD");
	    BEGIN(ERROR);
	
	  } else {
	    BEGIN(DISCARD);
	  }
	}

<CCi_ja>. {
	  BEGIN(DISCARD);
	}

<CROTAi>{I1}{ALT}" " |
<CROTAi>{I2}{ALT} {
	  sscanf(yytext, "%d%c", &i, &a);
	  if (a == ' ' || relax & WCSHDR_CROTAia) {
	    idx = i - 1;
	    BEGIN(VALUE);
	
	  } else if (relax & WCSHDR_reject) {
	    errmsg = "CROTAn keyword may not have an alternate version code";
	    BEGIN(ERROR);
	
	  } else {
	    BEGIN(DISCARD);
	  }
	}

<CROTAi>{I3} {
	  sscanf(yytext, "%d", &i);
	  a = ' ';
	  idx = i - 1;
	  BEGIN(VALUE);
	}

<CROTAi>. {
	  BEGIN(DISCARD);
	}

<CCCCCCCa>{ALT} |
<CCCCCCCC>. {
	  idx = -1;
	
	  if (YY_START == CCCCCCCa) {
	    sscanf(yytext, "%c", &a);
	  } else {
	    unput(yytext[0]);
	    a = 0;
	  }
	  BEGIN(VALUE);
	}

<CCCCCCCa>. {
	  BEGIN(DISCARD);
	}

<CCi_ma>{I1}_{I0}{ALT}"  " |
<CCi_ma>{I1}_{I2}{ALT}" " |
<CCi_ma>{I2}_{I0}{ALT}" " |
<CCi_ma>{I2}_{I2}{ALT} {
	  sscanf(yytext, "%d_%d%c", &i, &m, &a);
	  idx = -1;
	  BEGIN(VALUE);
	}

<CCi_ma>{I1}_{I3}{ALT} |
<CCi_ma>{I3}_{I0}{ALT} |
<CCi_ma>{I1}_{I4} |
<CCi_ma>{I2}_{I3} |
<CCi_ma>{I3}_{I2} |
<CCi_ma>{I4}_{I0} {
	  /* Invalid parameters will be caught by <VALUE>. */
	  sscanf(yytext, "%d_%d", &i, &m);
	  BEGIN(VALUE);
	}

<CCi_ma>. {
	  BEGIN(DISCARD);
	}

<PROJPn>{I0}"  " {
	  if (relax & WCSHDR_PROJPn) {
	    sscanf(yytext, "%d", &m);
	    i = 0;
	    a = ' ';
	    idx = -1;
	    BEGIN(VALUE);
	
	  } else if (relax & WCSHDR_reject) {
	    errmsg = "Defunct PROJPn keyword rejected";
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
	  if (i > 99 || j > 99 || m > 99) {
	    if (relax & WCSHDR_reject) {
	      errmsg = errtxt;
	      if (i > 99 || j > 99) {
	        sprintf(errmsg, "Axis number exceeds 99");
	      } else if (m > 99) {
	        sprintf(errmsg, "Parameter number exceeds 99");
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
	      sprintf(errmsg, "Internal parser ERROR, bad data type: %d",
	        valtype);
	      BEGIN(ERROR);
	    }
	  }
	}

<VALUE>. {
	  errmsg = "Invalid KEYWORD = VALUE syntax";
	  BEGIN(ERROR);
	}

<INTEGER_VAL>{INTEGER} {
	  if (pass == 1) {
	    wcspih_naxes(naxis, i, j, a, alts, npptr);
	    BEGIN(FLUSH);
	
	  } else {
	    if (vptr) {
	      /* Determine the coordinate representation. */
	      for (ialt = 0; ialt < *nwcs; ialt++) {
	        /* The loop here is for keywords that apply */
	        /* to every alternate; these have a == 0. */
	        if (a >= 'A') {
	          ialt = alts[a-'A'+1];
	        }
	
	        wptr = vptr;
	        if (ialt) {
	          voff = (char *)(*wcs+ialt) - (char *)(*wcs);
	          wptr = (void *)((char *)vptr + voff);
	        }
	
	        /* Apply keyword parameterization. */
	        if (idx >= 0) {
	          wptr = *((int **)wptr) + idx;
	        }
	
	        /* Read the keyvalue. */
	        sscanf(yytext, "%d", (int *)wptr);
	
	        if (a) break;
	      }
	
	      BEGIN(COMMENT);
	
	    } else {
	      errmsg = "Internal parser ERROR, null int pointer";
	      BEGIN(ERROR);
	    }
	  }
	}

<INTEGER_VAL>. {
	  errmsg = "An integer value was expected";
	  BEGIN(ERROR);
	}

<FLOAT_VAL>{FLOAT} {
	  if (pass == 1) {
	    wcspih_naxes(naxis, i, j, a, alts, npptr);
	    BEGIN(FLUSH);
	
	  } else {
	    if (vptr) {
	      /* Determine the coordinate representation. */
	      for (ialt = 0; ialt < *nwcs; ialt++) {
	        /* The loop here is for keywords like MJD-OBS that */
	        /* apply to every alternate; these have a == 0.  */
	        if (a >= 'A') {
	          ialt = alts[a-'A'+1];
	        }
	
	        wptr = vptr;
	        if (ialt) {
	          voff = (char *)(*wcs+ialt) - (char *)(*wcs);
	          wptr = (void *)((char *)vptr + voff);
	        }
	
	        /* Apply keyword parameterization. */
	        if (idx >= 0) {
	          wptr = *((double **)wptr) + idx;
	
	        } else if (npptr == npv) {
	          ipx = (*wcs+ialt)->npv++;
	          (*wcs+ialt)->pv[ipx].i = i;
	          (*wcs+ialt)->pv[ipx].m = m;
	          wptr = &((*wcs+ialt)->pv[ipx].value);
	        }
	
	        /* Read the keyvalue. */
	        wcsutil_str2double(yytext, "%lf", (double *)wptr);
	
	        /* Flag the presence of PCi_ja, or CDi_ja and/or CROTAia. */
	        if (altlin) {
	          (*wcs+ialt)->altlin |= altlin;
	          altlin = 0;
	        }
	
	        if (a) break;
	      }
	
	      BEGIN(COMMENT);
	
	    } else {
	      errmsg = "Internal parser ERROR, null float pointer";
	      BEGIN(ERROR);
	    }
	  }
	}

<FLOAT_VAL>. {
	  errmsg = "A floating-point value was expected";
	  BEGIN(ERROR);
	}

<STRING_VAL>{STRING} {
	  if (pass == 1) {
	    wcspih_naxes(naxis, i, j, a, alts, npptr);
	    BEGIN(FLUSH);
	
	  } else {
	    if (vptr) {
	      /* Determine the coordinate representation. */
	      for (ialt = 0; ialt < *nwcs; ialt++) {
	        /* The loop here is for keywords like DATE-OBS that */
	        /* apply to every alternate; these have a == 0.   */
	        if (a >= 'A') {
	          ialt = alts[a-'A'+1];
	        }
	
	        wptr = vptr;
	        if (ialt) {
	          voff = (char *)(*wcs+ialt) - (char *)(*wcs);
	          wptr = (void *)((char *)vptr + voff);
	        }
	
	        /* Apply keyword parameterization. */
	        if (idx >= 0) {
	          wptr = *((char (**)[72])wptr) + idx;
	
	        } else if (npptr == nps) {
	          ipx = (*wcs+ialt)->nps++;
	          (*wcs+ialt)->ps[ipx].i = i;
	          (*wcs+ialt)->ps[ipx].m = m;
	          wptr = (*wcs+ialt)->ps[ipx].value;
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
	
	        if (a) break;
	      }
	
	      BEGIN(COMMENT);
	
	    } else {
	      errmsg = "Internal parser ERROR, null string pointer";
	      BEGIN(ERROR);
	    }
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
	  if (pass == 2) {
	    if (ctrl < 0) {
	      /* Preserve discards. */
	      keep = wcspih_hdr - 80;
	
	    } else if (ctrl > 2) {
	      wcsfprintf(stderr, "%.80s\n  Discarded.\n", wcspih_hdr-80);
	    }
	  }
	  BEGIN(FLUSH);
	}

<ERROR>.* {
	  (*nreject)++;
	  if (pass == 2) {
	    if (ctrl%10 == -1) {
	      /* Preserve rejects. */
	      keep = wcspih_hdr - 80;
	    }
	
	    if (abs(ctrl%10) > 1) {
	      wcsfprintf(stderr, "%.80s\n%4d: %s.\n", wcspih_hdr-80, *nreject,
	        errmsg);
	    }
	  }
	  BEGIN(FLUSH);
	}

<FLUSH>.*\n {
	  if (pass == 2 && keep) {
	    if (hptr < keep) {
	      strncpy(hptr, keep, 80);
	    }
	    hptr += 80;
	  }
	
	  i = j = m = 0;
	  a = ' ';
	  valtype = -1;
	  keep = 0x0;
	  altlin = 0;
	  npptr = 0x0;
	  BEGIN(INITIAL);
	}

<<EOF>>	 {
	  /* End-of-input. */
	  if (pass == 1) {
	    if ((status = wcspih_inits(naxis, alts, npv, nps, nwcs, wcs)) ||
	         *nwcs == 0) {
	      yylex_destroy();
	      return status;
	    }
	
	    if (abs(ctrl%10) > 2) {
	      if (*nwcs == 1) {
	        if (strcmp(wcs[0]->wcsname, "DEFAULTS") != 0) {
	          wcsfprintf(stderr, "Found one coordinate representation.\n");
	        }
	      } else {
	        wcsfprintf(stderr, "Found %d coordinate representations.\n",
	          *nwcs);
	      }
	    }
	
	    wcspih_hdr = header;
	    wcspih_nkeyrec = nkeyrec;
	    *nreject = 0;
	
	    pass = 2;
	    i = j = m = 0;
	    a = ' ';
	    valtype = -1;
	
	    yyrestart(yyin);
	
	  } else {
	    yylex_destroy();
	
	    if (ctrl < 0) {
	      *hptr = '\0';
	    } else if (ctrl == 1) {
	      wcsfprintf(stderr, "%d WCS keyrecords were rejected.\n",
	        *nreject);
	    }
	
	    return wcspih_final(alts, epoch, vsource, nwcs, wcs);
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
      }
    }

    /* Restore the original values of NPVMAX and NPSMAX. */
    wcsnpv(npvmax);
    wcsnps(npsmax);
  }

  return status;
}


/*----------------------------------------------------------------------------
* Interpret special keywords encountered for each coordinate representation.
*---------------------------------------------------------------------------*/

int wcspih_final(
  int alts[],
  double epoch[],
  double vsource[],
  int *nwcs,
  struct wcsprm **wcs)

{
  int ialt, status;
  double beta, c = 299792458.0;

  for (ialt = 0; ialt < *nwcs; ialt++) {
    /* Check for EPOCH overriding EQUINOXa. */
    if (undefined((*wcs+ialt)->equinox) && !undefined(epoch[ialt])) {
      /* Set EQUINOXa. */
      (*wcs+ialt)->equinox = epoch[ialt];
    }

    /* Check for VSOURCEa overriding ZSOURCEa. */
    if (undefined((*wcs+ialt)->zsource) && !undefined(vsource[ialt])) {
      /* Convert relativistic Doppler velocity to redshift. */
      beta = vsource[ialt]/c;
      (*wcs+ialt)->zsource = (1.0+beta)/sqrt(1.0 - beta*beta) - 1.0;
    }

    /* Interpret -TAB header keywords. */
    if ((status = wcstab(*wcs+ialt))) {
       wcsvfree(nwcs, wcs);
       return status;
    }
  }

  return 0;
}
