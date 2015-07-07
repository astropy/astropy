/*============================================================================

  WCSLIB 5.7 - an implementation of the FITS WCS standard.
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
  $Id: fitshdr.l,v 5.7 2015/06/29 02:44:16 mcalabre Exp $
*=============================================================================
*
* fitshdr.l is a Flex description file containing a lexical scanner
* definition for extracting keywords and keyvalues from a FITS header.
*
* It requires Flex v2.5.4 or later.
*
* Refer to fitshdr.h for a description of the user interface and operating
* notes.
*
*===========================================================================*/

/* Options. */
%option full
%option never-interactive
%option nounput
%option noyywrap
%option outfile="fitshdr.c"
%option prefix="fitshdr"

/* Keywords. */
KEYCHR	[-_A-Z0-9]
KW1	{KEYCHR}{1}" "{7}
KW2	{KEYCHR}{2}" "{6}
KW3	{KEYCHR}{3}" "{5}
KW4	{KEYCHR}{4}" "{4}
KW5	{KEYCHR}{5}" "{3}
KW6	{KEYCHR}{6}" "{2}
KW7	{KEYCHR}{7}" "{1}
KW8	{KEYCHR}{8}
KEYWORD	({KW1}|{KW2}|{KW3}|{KW4}|{KW5}|{KW6}|{KW7}|{KW8})

/* Keyvalue data types. */
LOGICAL	[TF]
INT32	[+-]?0*[0-9]{1,9}
INT64	[+-]?0*[0-9]{10,18}
INTVL	[+-]?0*[0-9]{19,}
INTEGER	[+-]?[0-9]+
FLOAT	[+-]?([0-9]+\.?[0-9]*|\.[0-9]+)([eEdD][+-]?[0-9]+)?
ICOMPLX	\(" "*{INTEGER}" "*," "*{INTEGER}" "*\)
FCOMPLX	\(" "*{FLOAT}" "*," "*{FLOAT}" "*\)
STRING	'([^']|'')*'

/* Characters forming standard unit strings (jwBIQX are not used). */
UNITSTR \[[-+*/^(). 0-9a-zA-Z]+\]

/* Exclusive start states. */
%x VALUE INLINE UNITS COMMENT ERROR FLUSH

%{
#include <math.h>
#include <limits.h>
#include <setjmp.h>
#include <stdlib.h>
#include <string.h>

#include "fitshdr.h"
#include "wcsutil.h"

#define YY_DECL int fitshdr(const char header[], int nkeyrec, int nkeyids, \
                            struct fitskeyid keyids[], int *nreject, \
                            struct fitskey **keys)

#define YY_INPUT(inbuff, count, bufsize) \
	{ \
	  if (fitshdr_nkeyrec) { \
	    strncpy(inbuff, fitshdr_hdr, 80); \
	    inbuff[80] = '\n'; \
	    fitshdr_hdr += 80; \
	    fitshdr_nkeyrec--; \
	    count = 81; \
	  } else { \
	    count = YY_NULL; \
	  } \
	}

/* These global variables are required by YY_INPUT. */
const char *fitshdr_hdr;
int  fitshdr_nkeyrec;

/* Used in preempting the call to exit() by yy_fatal_error(). */
jmp_buf fitshdr_abort_jmp_env;
#define exit(status) longjmp(fitshdr_abort_jmp_env, status)

/* Map status return value to message. */
const char *fitshdr_errmsg[] = {
   "Success",
   "Null fitskey pointer-pointer passed",
   "Memory allocation failed",
   "Fatal error returned by Flex parser"};

%}

%%
	char *cptr, ctmp[72];
	int  blank, continuation, end, j, k, keyno;
	double dtmp;
	struct fitskey *kptr;
	struct fitskeyid *iptr;
	void nullfill(char cptr[], int len);
	int  yylex_destroy(void);
	
	fitshdr_hdr = header;
	fitshdr_nkeyrec = nkeyrec;
	
	*nreject = 0;
	keyno = 0;
	
	if (keys == 0x0) {
	  return 1;
	}
	
	/* Allocate memory for the required number of fitskey structs. */
	/* Recall that calloc() initializes allocated memory to zero.  */
	if (!(kptr = *keys = calloc(nkeyrec, sizeof(struct fitskey)))) {
	  return 2;
	}
	
	/* Initialize keyids[]. */
	iptr = keyids;
	for (j = 0; j < nkeyids; j++, iptr++) {
	  iptr->count  = 0;
	  iptr->idx[0] = -1;
	  iptr->idx[1] = -1;
	}
	
	blank = 0;
	continuation = 0;
	end = 0;
	
	/* Return here via longjmp() invoked by yy_fatal_error(). */
	if (setjmp(fitshdr_abort_jmp_env)) {
	  return 3;
	}
	
	BEGIN(INITIAL);

^" "{80} {
	  /* A completely blank keyrecord. */
	  strncpy(kptr->keyword, yytext, 8);
	  yyless(0);
	  blank = 1;
	  BEGIN(COMMENT);
	}

^(COMMENT|HISTORY|" "{8}) {
	  strncpy(kptr->keyword, yytext, 8);
	  BEGIN(COMMENT);
	}

^END" "{77} {
	  strncpy(kptr->keyword, yytext, 8);
	  end = 1;
	  BEGIN(FLUSH);
	}

^END" "{5}=" "+ {
	  /* Illegal END keyrecord. */
	  strncpy(kptr->keyword, yytext, 8);
	  kptr->status |= FITSHDR_KEYREC;
	  BEGIN(VALUE);
	}

^END" "{5} {
	  /* Illegal END keyrecord. */
	  strncpy(kptr->keyword, yytext, 8);
	  kptr->status |= FITSHDR_KEYREC;
	  BEGIN(COMMENT);
	}

^{KEYWORD}=" "+ {
	  strncpy(kptr->keyword, yytext, 8);
	  BEGIN(VALUE);
	}

^CONTINUE"  "+{STRING} {
	  /* Continued string keyvalue. */
	  strncpy(kptr->keyword, yytext, 8);
	
	  if (keyno > 0 && (kptr-1)->type%10 == 8) {
	    /* Put back the string keyvalue. */
	    for (k = 10; yytext[k] != '\''; k++);
	    yyless(k);
	    continuation = 1;
	    BEGIN(VALUE);
	
	  } else {
	    /* Not a valid continuation. */
	    yyless(8);
	    BEGIN(COMMENT);
	  }
	}

^{KEYWORD} {
	  /* Keyword without value. */
	  strncpy(kptr->keyword, yytext, 8);
	  BEGIN(COMMENT);
	}

^.{8}=" "+ {
	  /* Illegal keyword, carry on regardless. */
	  strncpy(kptr->keyword, yytext, 8);
	  kptr->status |= FITSHDR_KEYWORD;
	  BEGIN(VALUE);
	}

^.{8}	{
	  /* Illegal keyword, carry on regardless. */
	  strncpy(kptr->keyword, yytext, 8);
	  kptr->status |= FITSHDR_KEYWORD;
	  BEGIN(COMMENT);
	}

<VALUE>" "*/\/ {
	  /* Null keyvalue. */
	  BEGIN(INLINE);
	}

<VALUE>{LOGICAL} {
	  /* Logical keyvalue. */
	  kptr->type = 1;
	  kptr->keyvalue.i = (*yytext == 'T');
	  BEGIN(INLINE);
	}

<VALUE>{INT32} {
	  /* 32-bit signed integer keyvalue. */
	  kptr->type = 2;
	  if (sscanf(yytext, "%d", &(kptr->keyvalue.i)) < 1) {
	    kptr->status |= FITSHDR_KEYVALUE;
	    BEGIN(ERROR);
	  }
	
	  BEGIN(INLINE);
	}

<VALUE>{INT64} {
	  /* 64-bit signed integer keyvalue (up to 18 digits). */
	  if (wcsutil_str2double(yytext, "%lf", &dtmp)) {
	    kptr->status |= FITSHDR_KEYVALUE;
	    BEGIN(ERROR);
	
	  } else if (INT_MIN <= dtmp && dtmp <= INT_MAX) {
	    /* Can be accomodated as a 32-bit signed integer. */
	    kptr->type = 2;
	    if (sscanf(yytext, "%d", &(kptr->keyvalue.i)) < 1) {
	      kptr->status |= FITSHDR_KEYVALUE;
	      BEGIN(ERROR);
	    }
	
	  } else {
	    /* 64-bit signed integer. */
	    kptr->type = 3;
	    #ifdef WCSLIB_INT64
	      /* Native 64-bit integer is available. */
	      if (sscanf(yytext, "%lld", &(kptr->keyvalue.k)) < 1) {
	        kptr->status |= FITSHDR_KEYVALUE;
	        BEGIN(ERROR);
	      }
	    #else
	      /* 64-bit integer (up to 18 digits) implemented as int[3]. */
	      kptr->keyvalue.k[2] = 0;
	
	      sprintf(ctmp, "%%%dd%%9d", yyleng-9);
	      if (sscanf(yytext, ctmp, kptr->keyvalue.k+1,
	                 kptr->keyvalue.k) < 1) {
	        kptr->status |= FITSHDR_KEYVALUE;
	        BEGIN(ERROR);
	      } else if (*yytext == '-') {
	        kptr->keyvalue.k[0] *= -1;
	      }
	    #endif
	  }
	
	  BEGIN(INLINE);
	}

<VALUE>{INTVL} {
	  /* Very long integer keyvalue (and 19-digit int64). */
	  kptr->type = 4;
	  strcpy(ctmp, yytext);
	  k = yyleng;
	  for (j = 0; j < 8; j++) {
	    /* Read it backwards. */
	    k -= 9;
	    if (k < 0) k = 0;
	    if (sscanf(ctmp+k, "%d", kptr->keyvalue.l+j) < 1) {
	      kptr->status |= FITSHDR_KEYVALUE;
	      BEGIN(ERROR);
	    }
	    if (*yytext == '-') {
	      kptr->keyvalue.l[j] = -abs(kptr->keyvalue.l[j]);
	    }
	
	    if (k == 0) break;
	    ctmp[k] = '\0';
	  }
	
	  /* Can it be accomodated as a 64-bit signed integer? */
	  if (j == 2 && abs(kptr->keyvalue.l[2]) <=  9 &&
	                abs(kptr->keyvalue.l[1]) <=  223372036 &&
	                    kptr->keyvalue.l[0]  <=  854775807 &&
	                    kptr->keyvalue.l[0]  >= -854775808) {
	    kptr->type = 3;
	
	    #ifdef WCSLIB_INT64
	      /* Native 64-bit integer is available. */
	      kptr->keyvalue.l[2] = 0;
	      if (sscanf(yytext, "%lld", &(kptr->keyvalue.k)) < 1) {
	        kptr->status |= FITSHDR_KEYVALUE;
	        BEGIN(ERROR);
	      }
	    #endif
	  }
	
	  BEGIN(INLINE);
	}

<VALUE>{FLOAT} {
	  /* Float keyvalue. */
	  kptr->type = 5;
	  if (wcsutil_str2double(yytext, "%lf", &(kptr->keyvalue.f))) {
	    kptr->status |= FITSHDR_KEYVALUE;
	    BEGIN(ERROR);
	  }
	
	  BEGIN(INLINE);
	}

<VALUE>{ICOMPLX} {
	  /* Integer complex keyvalue. */
	  kptr->type = 6;
	  if (sscanf(yytext, "(%lf,%lf)", kptr->keyvalue.c,
	      kptr->keyvalue.c+1) < 2) {
	    kptr->status |= FITSHDR_KEYVALUE;
	    BEGIN(ERROR);
	  }
	
	  BEGIN(INLINE);
	}

<VALUE>{FCOMPLX} {
	  /* Floating point complex keyvalue. */
	  kptr->type = 7;
	
	  for (cptr = ctmp, k = 1; yytext[k] != ','; cptr++, k++) {
	    *cptr = yytext[k];
	  }
	  *cptr = '\0';
	
	  if (wcsutil_str2double(ctmp, "%lf", kptr->keyvalue.c)) {
	    kptr->status |= FITSHDR_KEYVALUE;
	    BEGIN(ERROR);
	  }
	
	  for (cptr = ctmp, k++; yytext[k] != ')'; cptr++, k++) {
	    *cptr = yytext[k];
	  }
	  *cptr = '\0';
	
	  if (wcsutil_str2double(ctmp, "%lf", kptr->keyvalue.c+1)) {
	    kptr->status |= FITSHDR_KEYVALUE;
	    BEGIN(ERROR);
	  }
	
	  BEGIN(INLINE);
	}

<VALUE>{STRING} {
	  /* String keyvalue. */
	  kptr->type = 8;
	  cptr = kptr->keyvalue.s;
	  strcpy(cptr, yytext+1);
	
	  /* Squeeze out repeated quotes. */
	  k = 0;
	  for (j = 0; j < 72; j++) {
	    if (k < j) {
	      cptr[k] = cptr[j];
	    }
	
	    if (cptr[j] == '\0') {
	      if (k) cptr[k-1] = '\0';
	      break;
	    } else if (cptr[j] == '\'' && cptr[j+1] == '\'') {
	      j++;
	    }
	
	    k++;
	  }
	
	  if (*cptr) {
	    /* Retain the initial blank in all-blank strings. */
	    nullfill(cptr+1, 71);
	  } else {
	    nullfill(cptr, 72);
	  }
	
	  BEGIN(INLINE);
	}

<VALUE>. {
	  kptr->status |= FITSHDR_KEYVALUE;
	  BEGIN(ERROR);
	}

<INLINE>" "*$ {
	  BEGIN(FLUSH);
	}

<INLINE>" "*\/" "*$ {
	  BEGIN(FLUSH);
	}

<INLINE>" "*\/" "* {
	  BEGIN(UNITS);
	}

<INLINE>" " {
	  kptr->status |= FITSHDR_COMMENT;
	  BEGIN(ERROR);
	}

<INLINE>. {
	  /* Keyvalue parsing must now also be suspect. */
	  kptr->status |= FITSHDR_COMMENT;
	  kptr->type = 0;
	  BEGIN(ERROR);
	}

<UNITS>{UNITSTR} {
	  kptr->ulen = yyleng;
	  yymore();
	  BEGIN(COMMENT);
	}

<UNITS>. {
	  yymore();
	  BEGIN(COMMENT);
	}

<COMMENT>.* {
	  strcpy(kptr->comment, yytext);
	  nullfill(kptr->comment, 84);
	  BEGIN(FLUSH);
	}

<ERROR>.* {
	  if (!continuation) kptr->type = -abs(kptr->type);
	
	  sprintf(kptr->comment, "%.80s", fitshdr_hdr-80);
	  kptr->comment[80] = '\0';
	  nullfill(kptr->comment+80, 4);
	
	  BEGIN(FLUSH);
	}

<FLUSH>.*\n {
	  /* Discard the rest of the input line. */
	  kptr->keyno = ++keyno;
	
	  /* Null-fill the keyword. */
	  kptr->keyword[8] = '\0';
	  nullfill(kptr->keyword, 12);
	
	  /* Do indexing. */
	  iptr = keyids;
	  kptr->keyid = -1;
	  for (j = 0; j < nkeyids; j++, iptr++) {
	    cptr = iptr->name;
	    cptr[8] = '\0';
	    nullfill(cptr, 12);
	    for (k = 0; k < 8; k++, cptr++) {
	      if (*cptr != '.' && *cptr != kptr->keyword[k]) break;
	    }
	
	    if (k == 8) {
	      /* Found a match. */
	      iptr->count++;
	      if (iptr->idx[0] == -1) {
	        iptr->idx[0] = keyno-1;
	      } else {
	        iptr->idx[1] = keyno-1;
	      }
	
	      kptr->keyno = -abs(kptr->keyno);
	      if (kptr->keyid < 0) kptr->keyid = j;
	    }
	  }
	
	  /* Deal with continued strings. */
	  if (continuation) {
	    /* Tidy up the previous string keyvalue. */
	    if ((kptr-1)->type == 8) (kptr-1)->type += 10;
	    cptr = (kptr-1)->keyvalue.s;
	    if (cptr[strlen(cptr)-1] == '&') cptr[strlen(cptr)-1] = '\0';
	
	    kptr->type = (kptr-1)->type + 10;
	  }
	
	  /* Check for keyrecords following the END keyrecord. */
	  if (end && (end++ > 1) && !blank) {
	    kptr->status |= FITSHDR_TRAILER;
	  }
	  if (kptr->status) (*nreject)++;
	
	  kptr++;
	  blank = 0;
	  continuation = 0;
	
	  BEGIN(INITIAL);
	}

<<EOF>>	{
	  /* End-of-input. */
	  yylex_destroy();
	  return 0;
	}

%%

/*--------------------------------------------------------------------------*/

void nullfill(char cptr[], int len)

{
  int j, k;

  /* Null-fill the string. */
  for (j = 0; j < len; j++) {
    if (cptr[j] == '\0') {
      for (k = j+1; k < len; k++) {
        cptr[k] = '\0';
      }
      break;
    }
  }

  for (k = j-1; k >= 0; k--) {
    if (cptr[k] != ' ') break;
    cptr[k] = '\0';
  }

  return;
}
