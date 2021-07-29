/* A Bison parser, made by GNU Bison 3.7.4.  */

/* Bison interface for Yacc-like parsers in C

   Copyright (C) 1984, 1989-1990, 2000-2015, 2018-2020 Free Software Foundation,
   Inc.

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

/* DO NOT RELY ON FEATURES THAT ARE NOT DOCUMENTED in the manual,
   especially those whose name start with FF_ or ff_.  They are
   private implementation details that can be changed or removed.  */

#ifndef FF_FF_Y_TAB_H_INCLUDED
# define FF_FF_Y_TAB_H_INCLUDED
/* Debug traces.  */
#ifndef FFDEBUG
# define FFDEBUG 0
#endif
#if FFDEBUG
extern int ffdebug;
#endif

/* Token kinds.  */
#ifndef FFTOKENTYPE
# define FFTOKENTYPE
  enum fftokentype
  {
    FFEMPTY = -2,
    FFEOF = 0,                     /* "end of file"  */
    FFerror = 256,                 /* error  */
    FFUNDEF = 257,                 /* "invalid token"  */
    BOOLEAN = 258,                 /* BOOLEAN  */
    LONG = 259,                    /* LONG  */
    DOUBLE = 260,                  /* DOUBLE  */
    STRING = 261,                  /* STRING  */
    BITSTR = 262,                  /* BITSTR  */
    FUNCTION = 263,                /* FUNCTION  */
    BFUNCTION = 264,               /* BFUNCTION  */
    IFUNCTION = 265,               /* IFUNCTION  */
    GTIFILTER = 266,               /* GTIFILTER  */
    GTIOVERLAP = 267,              /* GTIOVERLAP  */
    REGFILTER = 268,               /* REGFILTER  */
    COLUMN = 269,                  /* COLUMN  */
    BCOLUMN = 270,                 /* BCOLUMN  */
    SCOLUMN = 271,                 /* SCOLUMN  */
    BITCOL = 272,                  /* BITCOL  */
    ROWREF = 273,                  /* ROWREF  */
    NULLREF = 274,                 /* NULLREF  */
    SNULLREF = 275,                /* SNULLREF  */
    OR = 276,                      /* OR  */
    AND = 277,                     /* AND  */
    EQ = 278,                      /* EQ  */
    NE = 279,                      /* NE  */
    GT = 280,                      /* GT  */
    LT = 281,                      /* LT  */
    LTE = 282,                     /* LTE  */
    GTE = 283,                     /* GTE  */
    XOR = 284,                     /* XOR  */
    POWER = 285,                   /* POWER  */
    NOT = 286,                     /* NOT  */
    INTCAST = 287,                 /* INTCAST  */
    FLTCAST = 288,                 /* FLTCAST  */
    UMINUS = 289,                  /* UMINUS  */
    ACCUM = 290,                   /* ACCUM  */
    DIFF = 291                     /* DIFF  */
  };
  typedef enum fftokentype fftoken_kind_t;
#endif
/* Token kinds.  */
#define FFEMPTY -2
#define FFEOF 0
#define FFerror 256
#define FFUNDEF 257
#define BOOLEAN 258
#define LONG 259
#define DOUBLE 260
#define STRING 261
#define BITSTR 262
#define FUNCTION 263
#define BFUNCTION 264
#define IFUNCTION 265
#define GTIFILTER 266
#define GTIOVERLAP 267
#define REGFILTER 268
#define COLUMN 269
#define BCOLUMN 270
#define SCOLUMN 271
#define BITCOL 272
#define ROWREF 273
#define NULLREF 274
#define SNULLREF 275
#define OR 276
#define AND 277
#define EQ 278
#define NE 279
#define GT 280
#define LT 281
#define LTE 282
#define GTE 283
#define XOR 284
#define POWER 285
#define NOT 286
#define INTCAST 287
#define FLTCAST 288
#define UMINUS 289
#define ACCUM 290
#define DIFF 291

/* Value type.  */
#if ! defined FFSTYPE && ! defined FFSTYPE_IS_DECLARED
union FFSTYPE
{
#line 199 "eval.y"

    int    Node;        /* Index of Node */
    double dbl;         /* real value    */
    long   lng;         /* integer value */
    char   log;         /* logical value */
    char   str[MAX_STRLEN];    /* string value  */

#line 147 "y.tab.h"

};
typedef union FFSTYPE FFSTYPE;
# define FFSTYPE_IS_TRIVIAL 1
# define FFSTYPE_IS_DECLARED 1
#endif


extern FFSTYPE fflval;

int ffparse (void);

#endif /* !FF_FF_Y_TAB_H_INCLUDED  */
