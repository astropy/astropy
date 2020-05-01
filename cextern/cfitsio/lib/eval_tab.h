/* A Bison parser, made by GNU Bison 3.0.5.  */

/* Bison interface for Yacc-like parsers in C

   Copyright (C) 1984, 1989-1990, 2000-2015, 2018 Free Software Foundation, Inc.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.

   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

#ifndef FF_FF_Y_TAB_H_INCLUDED
# define FF_FF_Y_TAB_H_INCLUDED
/* Debug traces.  */
#ifndef FFDEBUG
# define FFDEBUG 0
#endif
#if FFDEBUG
extern int ffdebug;
#endif

/* Token type.  */
#ifndef FFTOKENTYPE
# define FFTOKENTYPE
  enum fftokentype
  {
    BOOLEAN = 258,
    LONG = 259,
    DOUBLE = 260,
    STRING = 261,
    BITSTR = 262,
    FUNCTION = 263,
    BFUNCTION = 264,
    IFUNCTION = 265,
    GTIFILTER = 266,
    REGFILTER = 267,
    COLUMN = 268,
    BCOLUMN = 269,
    SCOLUMN = 270,
    BITCOL = 271,
    ROWREF = 272,
    NULLREF = 273,
    SNULLREF = 274,
    OR = 275,
    AND = 276,
    EQ = 277,
    NE = 278,
    GT = 279,
    LT = 280,
    LTE = 281,
    GTE = 282,
    XOR = 283,
    POWER = 284,
    NOT = 285,
    INTCAST = 286,
    FLTCAST = 287,
    UMINUS = 288,
    ACCUM = 289,
    DIFF = 290
  };
#endif
/* Tokens.  */
#define BOOLEAN 258
#define LONG 259
#define DOUBLE 260
#define STRING 261
#define BITSTR 262
#define FUNCTION 263
#define BFUNCTION 264
#define IFUNCTION 265
#define GTIFILTER 266
#define REGFILTER 267
#define COLUMN 268
#define BCOLUMN 269
#define SCOLUMN 270
#define BITCOL 271
#define ROWREF 272
#define NULLREF 273
#define SNULLREF 274
#define OR 275
#define AND 276
#define EQ 277
#define NE 278
#define GT 279
#define LT 280
#define LTE 281
#define GTE 282
#define XOR 283
#define POWER 284
#define NOT 285
#define INTCAST 286
#define FLTCAST 287
#define UMINUS 288
#define ACCUM 289
#define DIFF 290

/* Value type.  */
#if ! defined FFSTYPE && ! defined FFSTYPE_IS_DECLARED

union FFSTYPE
{
#line 192 "eval.y" /* yacc.c:1910  */

    int    Node;        /* Index of Node */
    double dbl;         /* real value    */
    long   lng;         /* integer value */
    char   log;         /* logical value */
    char   str[MAX_STRLEN];    /* string value  */

#line 132 "y.tab.h" /* yacc.c:1910  */
};

typedef union FFSTYPE FFSTYPE;
# define FFSTYPE_IS_TRIVIAL 1
# define FFSTYPE_IS_DECLARED 1
#endif


extern FFSTYPE fflval;

int ffparse (void);

#endif /* !FF_FF_Y_TAB_H_INCLUDED  */
