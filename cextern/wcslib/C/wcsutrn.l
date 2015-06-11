/*============================================================================

  WCSLIB 5.3 - an implementation of the FITS WCS standard.
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
  $Id: wcsutrn.l,v 5.3 2015/04/21 02:50:51 mcalabre Exp $
*=============================================================================
*
* wcsutrn.l is a Flex description file containing the definition of a lexical
* scanner that translates non-standard FITS units specifications.
*
* It requires Flex v2.5.4 or later.
*
* Refer to wcsunits.h for a description of the user interface and operating
* notes.
*
*===========================================================================*/

/* Options. */
%option full
%option never-interactive
%option noyywrap
%option outfile="wcsutrn.c"
%option prefix="wcsutrn"

/* Exclusive start states. */
%x NEXT FLUSH

%{
/* To get the prototype for fileno() from stdio.h when gcc is invoked with
 * -std=c89 (same as -ansi) or -std=c99 since we do not define YY_INPUT. */
#define _POSIX_SOURCE 1

#include <setjmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "wcserr.h"
#include "wcsunits.h"

#define YY_DECL int wcsutrne(int ctrl, char unitstr[], struct wcserr **err)

/* Used in preempting the call to exit() by yy_fatal_error(). */
jmp_buf wcsutrn_abort_jmp_env;
#define exit(status) longjmp(wcsutrn_abort_jmp_env, status)

%}

%%
	static const char *function = "wcsutrne";
	
	char orig[80], subs[80];
	int bracket = 0;
	int unsafe  = 0;
	int status  = -1;
	YY_BUFFER_STATE inbuff;
	int yylex_destroy(void);
	
	if (err) *err = 0x0;
	
	*orig = '\0';
	*subs = '\0';
	
	inbuff = yy_scan_string(unitstr);
	*unitstr = '\0';
	
	/* Return here via longjmp() invoked by yy_fatal_error(). */
	if (setjmp(wcsutrn_abort_jmp_env)) {
	  return wcserr_set(WCSERR_SET(UNITSERR_PARSER_ERROR),
	    "Internal units translator error parsing '%s'", unitstr);
	}
	
	BEGIN(INITIAL);
	
	#ifdef DEBUG
	fprintf(stderr, "\n%s ->\n", unitstr);
	#endif

^" "*"[" {
	  /* Looks like a keycomment. */
	  strcat(unitstr, "[");
	  bracket = 1;
	}

" "+	  /* Discard leading whitespace. */

[^A-Za-z] {
	  /* Non-alphabetic character. */
	  strcat(unitstr, yytext);
	  if (bracket && *yytext == ']') {
	    BEGIN(FLUSH);
	  }
	}

angstrom {
	  strcpy(orig, yytext);
	  strcpy(subs, "Angstrom");
	  BEGIN(NEXT);
	}

arcmins|ARCMINS? {
	  strcpy(orig, yytext);
	  strcpy(subs, "arcmin");
	  BEGIN(NEXT);
	}

arcsecs|ARCSECS? {
	  strcpy(orig, yytext);
	  strcpy(subs, "arcsec");
	  BEGIN(NEXT);
	}

BEAM	{
	  strcpy(orig, yytext);
	  strcpy(subs, "beam");
	  BEGIN(NEXT);
	}

Byte	{
	  strcpy(orig, yytext);
	  strcpy(subs, "byte");
	  BEGIN(NEXT);
	}

days?|DAYS? {
	  strcpy(orig, yytext);
	  strcpy(subs, "d");
	  BEGIN(NEXT);
	}

D	{
	  unsafe = 1;
	  strcpy(orig, yytext);
	  strcpy(subs, (ctrl & 4) ? "d" : "D");
	  BEGIN(NEXT);
	}

degrees?|DEG|DEGREES? {
	  strcpy(orig, yytext);
	  strcpy(subs, "deg");
	  BEGIN(NEXT);
	}

GHZ	{
	  strcpy(orig, yytext);
	  strcpy(subs, "GHz");
	  BEGIN(NEXT);
	}

hr|HR	{
	  strcpy(orig, yytext);
	  strcpy(subs, "h");
	  BEGIN(NEXT);
	}

H	{
	  unsafe = 1;
	  strcpy(orig, yytext);
	  strcpy(subs, (ctrl & 2) ? "h" : "H");
	  BEGIN(NEXT);
	}

hz|HZ	{
	  strcpy(orig, yytext);
	  strcpy(subs, "Hz");
	  BEGIN(NEXT);
	}

KHZ	{
	  strcpy(orig, yytext);
	  strcpy(subs, "kHz");
	  BEGIN(NEXT);
	}

JY	{
	  strcpy(orig, yytext);
	  strcpy(subs, "Jy");
	  BEGIN(NEXT);
	}

[kK]elvins?|KELVINS? {
	  strcpy(orig, yytext);
	  strcpy(subs, "K");
	  BEGIN(NEXT);
	}

KM	{
	  strcpy(orig, yytext);
	  strcpy(subs, "km");
	  BEGIN(NEXT);
	}

metres?|meters?|M|METRES?|METERS? {
	  strcpy(orig, yytext);
	  strcpy(subs, "m");
	  BEGIN(NEXT);
	}

MIN	{
	  strcpy(orig, yytext);
	  strcpy(subs, "min");
	  BEGIN(NEXT);
	}

MHZ	{
	  strcpy(orig, yytext);
	  strcpy(subs, "MHz");
	  BEGIN(NEXT);
	}

Ohm	{
	  strcpy(orig, yytext);
	  strcpy(subs, "ohm");
	  BEGIN(NEXT);
	}

[pP]ascals?|PASCALS? {
	  strcpy(orig, yytext);
	  strcpy(subs, "Pa");
	  BEGIN(NEXT);
	}

pixels|PIXELS? {
	  strcpy(orig, yytext);
	  strcpy(subs, "pixel");
	  BEGIN(NEXT);
	}

radians?|RAD|RADIANS? {
	  strcpy(orig, yytext);
	  strcpy(subs, "rad");
	  BEGIN(NEXT);
	}

sec|seconds?|SEC|SECONDS? {
	  strcpy(orig, yytext);
	  strcpy(subs, "s");
	  BEGIN(NEXT);
	}

S	{
	  unsafe = 1;
	  strcpy(orig, yytext);
	  strcpy(subs, (ctrl & 1) ? "s" : "S");
	  BEGIN(NEXT);
	}

[vV]olts?|VOLTS? {
	  strcpy(orig, yytext);
	  strcpy(subs, "V");
	  BEGIN(NEXT);
	}

years?|YR|YEARS? {
	  strcpy(orig, yytext);
	  strcpy(subs, "yr");
	  BEGIN(NEXT);
	}

[A-Za-z]+ {
	  /* Not a recognized alias. */
	  strcpy(orig, yytext);
	  strcpy(subs, orig);
	  BEGIN(NEXT);
	}

<NEXT>[A-Za-z]+ {
	  /* Reject the alias match. */
	  strcat(orig, yytext);
	  strcpy(subs, orig);
	}

<NEXT>" "+[^A-Za-z] {
	  /* Discard separating whitespace. */
	  unput(yytext[yyleng-1]);
	}

<NEXT>" "+[A-Za-z] {
	  /* Compress separating whitespace. */
	  strcat(unitstr, subs);
	  strcat(unitstr, " ");
	  if (strcmp(orig, subs)) status = 0;
	  unput(yytext[yyleng-1]);
	  *subs = '\0';
	  BEGIN(INITIAL);
	}

<NEXT>.	{
	  /* Copy anything else unchanged. */
	  strcat(unitstr, subs);
	  if (strcmp(orig, subs)) status = 0;
	  unput(*yytext);
	  *subs = '\0';
	  BEGIN(INITIAL);
	}

<FLUSH>.* {
	  /* Copy out remaining input. */
	  strcat(unitstr, yytext);
	}

<<EOF>>	{
	  /* End-of-string. */
	  if (*subs) {
	    strcat(unitstr, subs);
	    if (strcmp(orig, subs)) status = 0;
	  }
	
	  yylex_destroy();
	  if (unsafe) {
	    return wcserr_set(WCSERR_SET(UNITSERR_UNSAFE_TRANS),
	      "Unsafe unit translation in '%s'", unitstr);
	  }
	  return status;
	}

%%
