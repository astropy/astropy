/* A Bison parser, made by GNU Bison 3.8.  */

/* Bison interface for Yacc-like parsers in C

   Copyright (C) 1984, 1989-1990, 2000-2015, 2018-2021 Free Software Foundation,
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
   along with this program.  If not, see <https://www.gnu.org/licenses/>.  */

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
   especially those whose name start with YY_ or yy_.  They are
   private implementation details that can be changed or removed.  */

#ifndef YY_FITS_PARSER_YY_EVAL_TAB_H_INCLUDED
# define YY_FITS_PARSER_YY_EVAL_TAB_H_INCLUDED
/* Debug traces.  */
#ifndef FITS_PARSER_YYDEBUG
# if defined YYDEBUG
#if YYDEBUG
#   define FITS_PARSER_YYDEBUG 1
#  else
#   define FITS_PARSER_YYDEBUG 0
#  endif
# else /* ! defined YYDEBUG */
#  define FITS_PARSER_YYDEBUG 0
# endif /* ! defined YYDEBUG */
#endif  /* ! defined FITS_PARSER_YYDEBUG */
#if FITS_PARSER_YYDEBUG
extern int fits_parser_yydebug;
#endif

/* Token kinds.  */
#ifndef FITS_PARSER_YYTOKENTYPE
# define FITS_PARSER_YYTOKENTYPE
  enum fits_parser_yytokentype
  {
    FITS_PARSER_YYEMPTY = -2,
    FITS_PARSER_YYEOF = 0,         /* "end of file"  */
    FITS_PARSER_YYerror = 256,     /* error  */
    FITS_PARSER_YYUNDEF = 257,     /* "invalid token"  */
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
    GTIFIND = 268,                 /* GTIFIND  */
    REGFILTER = 269,               /* REGFILTER  */
    COLUMN = 270,                  /* COLUMN  */
    BCOLUMN = 271,                 /* BCOLUMN  */
    SCOLUMN = 272,                 /* SCOLUMN  */
    BITCOL = 273,                  /* BITCOL  */
    ROWREF = 274,                  /* ROWREF  */
    NULLREF = 275,                 /* NULLREF  */
    SNULLREF = 276,                /* SNULLREF  */
    OR = 277,                      /* OR  */
    AND = 278,                     /* AND  */
    EQ = 279,                      /* EQ  */
    NE = 280,                      /* NE  */
    GT = 281,                      /* GT  */
    LT = 282,                      /* LT  */
    LTE = 283,                     /* LTE  */
    GTE = 284,                     /* GTE  */
    XOR = 285,                     /* XOR  */
    POWER = 286,                   /* POWER  */
    NOT = 287,                     /* NOT  */
    INTCAST = 288,                 /* INTCAST  */
    FLTCAST = 289,                 /* FLTCAST  */
    UMINUS = 290,                  /* UMINUS  */
    ACCUM = 291,                   /* ACCUM  */
    DIFF = 292                     /* DIFF  */
  };
  typedef enum fits_parser_yytokentype fits_parser_yytoken_kind_t;
#endif

/* Value type.  */
#if ! defined FITS_PARSER_YYSTYPE && ! defined FITS_PARSER_YYSTYPE_IS_DECLARED
union FITS_PARSER_YYSTYPE
{
#line 212 "eval.y"

    int    Node;        /* Index of Node */
    double dbl;         /* real value    */
    long   lng;         /* integer value */
    char   log;         /* logical value */
    char   str[MAX_STRLEN];    /* string value  */

#line 117 "eval_tab.h"

};
typedef union FITS_PARSER_YYSTYPE FITS_PARSER_YYSTYPE;
# define FITS_PARSER_YYSTYPE_IS_TRIVIAL 1
# define FITS_PARSER_YYSTYPE_IS_DECLARED 1
#endif




int fits_parser_yyparse (yyscan_t scanner, ParseData *lParse);


#endif /* !YY_FITS_PARSER_YY_EVAL_TAB_H_INCLUDED  */
