/* A Bison parser, made by GNU Bison 3.8.  */

/* Bison implementation for Yacc-like parsers in C

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

/* C LALR(1) parser skeleton written by Richard Stallman, by
   simplifying the original so-called "semantic" parser.  */

/* DO NOT RELY ON FEATURES THAT ARE NOT DOCUMENTED in the manual,
   especially those whose name start with YY_ or yy_.  They are
   private implementation details that can be changed or removed.  */

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Identify Bison output, and Bison version.  */
#define YYBISON 30800

/* Bison version string.  */
#define YYBISON_VERSION "3.8"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 2

/* Push parsers.  */
#define YYPUSH 0

/* Pull parsers.  */
#define YYPULL 1

/* Substitute the type names.  */
#define YYSTYPE         FITS_PARSER_YYSTYPE
/* Substitute the variable and function names.  */
#define yyparse         fits_parser_yyparse
#define yylex           fits_parser_yylex
#define yyerror         fits_parser_yyerror
#define yydebug         fits_parser_yydebug
#define yynerrs         fits_parser_yynerrs

/* First part of user prologue.  */
#line 16 "eval.y"

/* This file is one of 3 files containing code which parses an          */
/* arithmetic expression and evaluates it in the context of an input    */
/* FITS file table extension.  The CFITSIO lexical parser is divided    */
/* into the following 3 parts/files: the CFITSIO "front-end",           */
/* eval_f.c, contains the interface between the user/CFITSIO and the    */
/* real core of the parser; the FLEX interpreter, eval_l.c, takes the   */
/* input string and parses it into tokens and identifies the FITS       */
/* information required to evaluate the expression (ie, keywords and    */
/* columns); and, the BISON grammar and evaluation routines, eval_y.c,  */
/* receives the FLEX output and determines and performs the actual      */
/* operations.  The files eval_l.c and eval_y.c are produced from       */
/* running flex and bison on the files eval.l and eval.y, respectively. */
/* (flex and bison are available from any GNU archive: see www.gnu.org) */
/*                                                                      */
/* The grammar rules, rather than evaluating the expression in situ,    */
/* builds a tree, or Nodal, structure mapping out the order of          */
/* operations and expression dependencies.  This "compilation" process  */
/* allows for much faster processing of multiple rows.  This technique  */
/* was developed by Uwe Lammers of the XMM Science Analysis System,     */
/* although the CFITSIO implementation is entirely code original.       */
/*                                                                      */
/*                                                                      */
/* Modification History:                                                */
/*                                                                      */
/*   Kent Blackburn      c1992  Original parser code developed for the  */
/*                              FTOOLS software package, in particular, */
/*                              the fselect task.                       */
/*   Kent Blackburn      c1995  BIT column support added                */
/*   Peter D Wilson   Feb 1998  Vector column support added             */
/*   Peter D Wilson   May 1998  Ported to CFITSIO library.  User        */
/*                              interface routines written, in essence  */
/*                              making fselect, fcalc, and maketime     */
/*                              capabilities available to all tools     */
/*                              via single function calls.              */
/*   Peter D Wilson   Jun 1998  Major rewrite of parser core, so as to  */
/*                              create a run-time evaluation tree,      */
/*                              inspired by the work of Uwe Lammers,    */
/*                              resulting in a speed increase of        */
/*                              10-100 times.                           */
/*   Peter D Wilson   Jul 1998  gtifilter(a,b,c,d) function added       */
/*   Peter D Wilson   Aug 1998  regfilter(a,b,c,d) function added       */
/*   Peter D Wilson   Jul 1999  Make parser fitsfile-independent,       */
/*                              allowing a purely vector-based usage    */
/*  Craig B Markwardt Jun 2004  Add MEDIAN() function                   */
/*  Craig B Markwardt Jun 2004  Add SUM(), and MIN/MAX() for bit arrays */
/*  Craig B Markwardt Jun 2004  Allow subscripting of nX bit arrays     */
/*  Craig B Markwardt Jun 2004  Implement statistical functions         */
/*                              NVALID(), AVERAGE(), and STDDEV()       */
/*                              for integer and floating point vectors  */
/*  Craig B Markwardt Jun 2004  Use NULL values for range errors instead*/
/*                              of throwing a parse error               */
/*  Craig B Markwardt Oct 2004  Add ACCUM() and SEQDIFF() functions     */
/*  Craig B Markwardt Feb 2005  Add ANGSEP() function                   */
/*  Craig B Markwardt Aug 2005  CIRCLE, BOX, ELLIPSE, NEAR and REGFILTER*/
/*                              functions now accept vector arguments   */
/*  Craig B Markwardt Sum 2006  Add RANDOMN() and RANDOMP() functions   */
/*  Craig B Markwardt Mar 2007  Allow arguments to RANDOM and RANDOMN to*/
/*                              determine the output dimensions         */
/*  Craig B Markwardt Aug 2009  Add substring STRMID() and string search*/
/*                              STRSTR() functions; more overflow checks*/
/*  Craig B Markwardt Dec 2019  Add bit/hex/oct literal strings and     */
/*                              bitwise operatiosn between integers     */
/*  Craig B Markwardt Mar 2021  Add SETNULL() function                  */
/*                                                                      */
/************************************************************************/

#define  APPROX 1.0e-7
#include "eval_defs.h"
#include "region.h"
#include <time.h>

#include <stdlib.h>

#ifndef alloca
#define alloca malloc
#endif

/* Random number generators for various distributions */
#include "simplerng.h"

   /*  Shrink the initial stack depth to keep local data <32K (mac limit)  */
   /*  yacc will allocate more space if needed, though.                    */
#define  YYINITDEPTH   100

/***************************************************************/
/*  Replace Bison's BACKUP macro with one that fixes a bug --  */
/*  must update state after popping the stack -- and allows    */
/*  popping multiple terms at one time.                        */
/***************************************************************/

#define YYNEWBACKUP(token, value)                               \
   do								\
     if (yychar == YYEMPTY )   					\
       { yychar = (token);                                      \
         memcpy( &yylval, &(value), sizeof(value) );            \
         yychar1 = YYTRANSLATE (yychar);			\
         while (yylen--) YYPOPSTACK;				\
         yystate = *yyssp;					\
         goto yybackup;						\
       }							\
     else							\
       { yyerror ("syntax error: cannot back up"); YYERROR; }	\
   while (0)

/***************************************************************/
/*  Useful macros for accessing/testing Nodes                  */
/***************************************************************/

#define TEST(a)        if( (a)<0 ) YYERROR
#define SIZE(a)        lParse->Nodes[ a ].value.nelem
#define TYPE(a)        lParse->Nodes[ a ].type
#define OPER(a)        lParse->Nodes[ a ].operation
#define PROMOTE(a,b)   if( TYPE(a) > TYPE(b) )                  \
                          b = New_Unary( lParse, TYPE(a), 0, b );       \
                       else if( TYPE(a) < TYPE(b) )             \
	                  a = New_Unary( lParse, TYPE(b), 0, a );

/*****  Internal functions  *****/

#ifdef __cplusplus
extern "C" {
#endif

static int  Alloc_Node    ( ParseData * );
static void Free_Last_Node( ParseData * );
static void Evaluate_Node ( ParseData *, int thisNode );

static int  New_Const ( ParseData *, int returnType, void *value, long len );
static int  New_Column( ParseData *, int ColNum );
static int  New_Offset( ParseData *, int ColNum, int offset );
static int  New_Unary ( ParseData *, int returnType, int Op, int Node1 );
static int  New_BinOp ( ParseData *, int returnType, int Node1, int Op, int Node2 );
static int  New_Func  ( ParseData *, int returnType, funcOp Op, int nNodes,
			int Node1, int Node2, int Node3, int Node4, 
			int Node5, int Node6, int Node7 );
static int  New_FuncSize( ParseData *, int returnType, funcOp Op, int nNodes,
			int Node1, int Node2, int Node3, int Node4, 
			  int Node5, int Node6, int Node7, int Size);
static int  New_Deref ( ParseData *, int Var,  int nDim,
			int Dim1, int Dim2, int Dim3, int Dim4, int Dim5 );
static int  New_GTI   ( ParseData *, funcOp Op, char *fname, int Node1, int Node2, char *start, char *stop );
static int  New_REG   ( ParseData *, char *fname, int NodeX, int NodeY, char *colNames );
static int  New_Vector( ParseData *, int subNode );
static int  Close_Vec ( ParseData *, int vecNode );
static int  New_Array(  ParseData *, int valueNode, int dimNode );
static int  Locate_Col( ParseData *, Node *this );
static int  Test_Dims ( ParseData *, int Node1, int Node2 );
static void Copy_Dims ( ParseData *, int Node1, int Node2 );

static void Allocate_Ptrs( ParseData *, Node *this );
static void Do_Unary     ( ParseData *, Node *this );
static void Do_Offset    ( ParseData *, Node *this );
static void Do_BinOp_bit ( ParseData *, Node *this );
static void Do_BinOp_str ( ParseData *, Node *this );
static void Do_BinOp_log ( ParseData *, Node *this );
static void Do_BinOp_lng ( ParseData *, Node *this );
static void Do_BinOp_dbl ( ParseData *, Node *this );
static void Do_Func      ( ParseData *, Node *this );
static void Do_Deref     ( ParseData *, Node *this );
static void Do_GTI       ( ParseData *, Node *this );
static void Do_GTI_Over  ( ParseData *, Node *this );
static void Do_REG       ( ParseData *, Node *this );
static void Do_Vector    ( ParseData *, Node *this );
static void Do_Array     ( ParseData *, Node *this );

static long Search_GTI   ( double evtTime, long nGTI, double *start,
			   double *stop, int ordered, long *nextGTI );
static double GTI_Over(double evtStart, double evtStop,
		       long nGTI, double *start, double *stop,
		       long *gtiout);

static char  saobox (double xcen, double ycen, double xwid, double ywid,
		     double rot,  double xcol, double ycol);
static char  ellipse(double xcen, double ycen, double xrad, double yrad,
		     double rot, double xcol, double ycol);
static char  circle (double xcen, double ycen, double rad,
		     double xcol, double ycol);
static char  bnear  (double x, double y, double tolerance);
static char  bitcmp (char *bitstrm1, char *bitstrm2);
static char  bitlgte(char *bits1, int oper, char *bits2);

static void  bitand(char *result, char *bitstrm1, char *bitstrm2);
static void  bitor (char *result, char *bitstrm1, char *bitstrm2);
static void  bitnot(char *result, char *bits);
static int cstrmid(ParseData *lParse, char *dest_str, int dest_len,
		   char *src_str,  int src_len, int pos);

static void yyerror(yyscan_t scanner, ParseData *lParse, char *s);

#ifdef __cplusplus
    }
#endif


#line 273 "eval_y.c"

# ifndef YY_CAST
#  ifdef __cplusplus
#   define YY_CAST(Type, Val) static_cast<Type> (Val)
#   define YY_REINTERPRET_CAST(Type, Val) reinterpret_cast<Type> (Val)
#  else
#   define YY_CAST(Type, Val) ((Type) (Val))
#   define YY_REINTERPRET_CAST(Type, Val) ((Type) (Val))
#  endif
# endif
# ifndef YY_NULLPTR
#  if defined __cplusplus
#   if 201103L <= __cplusplus
#    define YY_NULLPTR nullptr
#   else
#    define YY_NULLPTR 0
#   endif
#  else
#   define YY_NULLPTR ((void*)0)
#  endif
# endif

#include "eval_tab.h"
/* Symbol kind.  */
enum yysymbol_kind_t
{
  YYSYMBOL_YYEMPTY = -2,
  YYSYMBOL_YYEOF = 0,                      /* "end of file"  */
  YYSYMBOL_YYerror = 1,                    /* error  */
  YYSYMBOL_YYUNDEF = 2,                    /* "invalid token"  */
  YYSYMBOL_BOOLEAN = 3,                    /* BOOLEAN  */
  YYSYMBOL_LONG = 4,                       /* LONG  */
  YYSYMBOL_DOUBLE = 5,                     /* DOUBLE  */
  YYSYMBOL_STRING = 6,                     /* STRING  */
  YYSYMBOL_BITSTR = 7,                     /* BITSTR  */
  YYSYMBOL_FUNCTION = 8,                   /* FUNCTION  */
  YYSYMBOL_BFUNCTION = 9,                  /* BFUNCTION  */
  YYSYMBOL_IFUNCTION = 10,                 /* IFUNCTION  */
  YYSYMBOL_GTIFILTER = 11,                 /* GTIFILTER  */
  YYSYMBOL_GTIOVERLAP = 12,                /* GTIOVERLAP  */
  YYSYMBOL_GTIFIND = 13,                   /* GTIFIND  */
  YYSYMBOL_REGFILTER = 14,                 /* REGFILTER  */
  YYSYMBOL_COLUMN = 15,                    /* COLUMN  */
  YYSYMBOL_BCOLUMN = 16,                   /* BCOLUMN  */
  YYSYMBOL_SCOLUMN = 17,                   /* SCOLUMN  */
  YYSYMBOL_BITCOL = 18,                    /* BITCOL  */
  YYSYMBOL_ROWREF = 19,                    /* ROWREF  */
  YYSYMBOL_NULLREF = 20,                   /* NULLREF  */
  YYSYMBOL_SNULLREF = 21,                  /* SNULLREF  */
  YYSYMBOL_22_ = 22,                       /* ','  */
  YYSYMBOL_23_ = 23,                       /* '='  */
  YYSYMBOL_24_ = 24,                       /* ':'  */
  YYSYMBOL_25_ = 25,                       /* '{'  */
  YYSYMBOL_26_ = 26,                       /* '}'  */
  YYSYMBOL_27_ = 27,                       /* '?'  */
  YYSYMBOL_OR = 28,                        /* OR  */
  YYSYMBOL_AND = 29,                       /* AND  */
  YYSYMBOL_EQ = 30,                        /* EQ  */
  YYSYMBOL_NE = 31,                        /* NE  */
  YYSYMBOL_32_ = 32,                       /* '~'  */
  YYSYMBOL_GT = 33,                        /* GT  */
  YYSYMBOL_LT = 34,                        /* LT  */
  YYSYMBOL_LTE = 35,                       /* LTE  */
  YYSYMBOL_GTE = 36,                       /* GTE  */
  YYSYMBOL_37_ = 37,                       /* '+'  */
  YYSYMBOL_38_ = 38,                       /* '-'  */
  YYSYMBOL_39_ = 39,                       /* '%'  */
  YYSYMBOL_40_ = 40,                       /* '*'  */
  YYSYMBOL_41_ = 41,                       /* '/'  */
  YYSYMBOL_42_ = 42,                       /* '|'  */
  YYSYMBOL_43_ = 43,                       /* '&'  */
  YYSYMBOL_XOR = 44,                       /* XOR  */
  YYSYMBOL_POWER = 45,                     /* POWER  */
  YYSYMBOL_NOT = 46,                       /* NOT  */
  YYSYMBOL_INTCAST = 47,                   /* INTCAST  */
  YYSYMBOL_FLTCAST = 48,                   /* FLTCAST  */
  YYSYMBOL_UMINUS = 49,                    /* UMINUS  */
  YYSYMBOL_50_ = 50,                       /* '['  */
  YYSYMBOL_ACCUM = 51,                     /* ACCUM  */
  YYSYMBOL_DIFF = 52,                      /* DIFF  */
  YYSYMBOL_53_n_ = 53,                     /* '\n'  */
  YYSYMBOL_54_ = 54,                       /* ']'  */
  YYSYMBOL_55_ = 55,                       /* '('  */
  YYSYMBOL_56_ = 56,                       /* ')'  */
  YYSYMBOL_YYACCEPT = 57,                  /* $accept  */
  YYSYMBOL_lines = 58,                     /* lines  */
  YYSYMBOL_line = 59,                      /* line  */
  YYSYMBOL_bvector = 60,                   /* bvector  */
  YYSYMBOL_vector = 61,                    /* vector  */
  YYSYMBOL_expr = 62,                      /* expr  */
  YYSYMBOL_bexpr = 63,                     /* bexpr  */
  YYSYMBOL_bits = 64,                      /* bits  */
  YYSYMBOL_sexpr = 65                      /* sexpr  */
};
typedef enum yysymbol_kind_t yysymbol_kind_t;




#ifdef short
# undef short
#endif

/* On compilers that do not define __PTRDIFF_MAX__ etc., make sure
   <limits.h> and (if available) <stdint.h> are included
   so that the code can choose integer types of a good width.  */

#ifndef __PTRDIFF_MAX__
# include <limits.h> /* INFRINGES ON USER NAME SPACE */
# if defined __STDC_VERSION__ && 199901 <= __STDC_VERSION__
#  include <stdint.h> /* INFRINGES ON USER NAME SPACE */
#  define YY_STDINT_H
# endif
#endif

/* Narrow types that promote to a signed type and that can represent a
   signed or unsigned integer of at least N bits.  In tables they can
   save space and decrease cache pressure.  Promoting to a signed type
   helps avoid bugs in integer arithmetic.  */

#ifdef __INT_LEAST8_MAX__
typedef __INT_LEAST8_TYPE__ yytype_int8;
#elif defined YY_STDINT_H
typedef int_least8_t yytype_int8;
#else
typedef signed char yytype_int8;
#endif

#ifdef __INT_LEAST16_MAX__
typedef __INT_LEAST16_TYPE__ yytype_int16;
#elif defined YY_STDINT_H
typedef int_least16_t yytype_int16;
#else
typedef short yytype_int16;
#endif

/* Work around bug in HP-UX 11.23, which defines these macros
   incorrectly for preprocessor constants.  This workaround can likely
   be removed in 2023, as HPE has promised support for HP-UX 11.23
   (aka HP-UX 11i v2) only through the end of 2022; see Table 2 of
   <https://h20195.www2.hpe.com/V2/getpdf.aspx/4AA4-7673ENW.pdf>.  */
#ifdef __hpux
# undef UINT_LEAST8_MAX
# undef UINT_LEAST16_MAX
# define UINT_LEAST8_MAX 255
# define UINT_LEAST16_MAX 65535
#endif

#if defined __UINT_LEAST8_MAX__ && __UINT_LEAST8_MAX__ <= __INT_MAX__
typedef __UINT_LEAST8_TYPE__ yytype_uint8;
#elif (!defined __UINT_LEAST8_MAX__ && defined YY_STDINT_H \
       && UINT_LEAST8_MAX <= INT_MAX)
typedef uint_least8_t yytype_uint8;
#elif !defined __UINT_LEAST8_MAX__ && UCHAR_MAX <= INT_MAX
typedef unsigned char yytype_uint8;
#else
typedef short yytype_uint8;
#endif

#if defined __UINT_LEAST16_MAX__ && __UINT_LEAST16_MAX__ <= __INT_MAX__
typedef __UINT_LEAST16_TYPE__ yytype_uint16;
#elif (!defined __UINT_LEAST16_MAX__ && defined YY_STDINT_H \
       && UINT_LEAST16_MAX <= INT_MAX)
typedef uint_least16_t yytype_uint16;
#elif !defined __UINT_LEAST16_MAX__ && USHRT_MAX <= INT_MAX
typedef unsigned short yytype_uint16;
#else
typedef int yytype_uint16;
#endif

#ifndef YYPTRDIFF_T
# if defined __PTRDIFF_TYPE__ && defined __PTRDIFF_MAX__
#  define YYPTRDIFF_T __PTRDIFF_TYPE__
#  define YYPTRDIFF_MAXIMUM __PTRDIFF_MAX__
# elif defined PTRDIFF_MAX
#  ifndef ptrdiff_t
#   include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  endif
#  define YYPTRDIFF_T ptrdiff_t
#  define YYPTRDIFF_MAXIMUM PTRDIFF_MAX
# else
#  define YYPTRDIFF_T long
#  define YYPTRDIFF_MAXIMUM LONG_MAX
# endif
#endif

#ifndef YYSIZE_T
# ifdef __SIZE_TYPE__
#  define YYSIZE_T __SIZE_TYPE__
# elif defined size_t
#  define YYSIZE_T size_t
# elif defined __STDC_VERSION__ && 199901 <= __STDC_VERSION__
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# else
#  define YYSIZE_T unsigned
# endif
#endif

#define YYSIZE_MAXIMUM                                  \
  YY_CAST (YYPTRDIFF_T,                                 \
           (YYPTRDIFF_MAXIMUM < YY_CAST (YYSIZE_T, -1)  \
            ? YYPTRDIFF_MAXIMUM                         \
            : YY_CAST (YYSIZE_T, -1)))

#define YYSIZEOF(X) YY_CAST (YYPTRDIFF_T, sizeof (X))


/* Stored state numbers (used for stacks). */
typedef yytype_int16 yy_state_t;

/* State numbers in computations.  */
typedef int yy_state_fast_t;

#ifndef YY_
# if defined YYENABLE_NLS && YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> /* INFRINGES ON USER NAME SPACE */
#   define YY_(Msgid) dgettext ("bison-runtime", Msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(Msgid) Msgid
# endif
#endif


#ifndef YY_ATTRIBUTE_PURE
# if defined __GNUC__ && 2 < __GNUC__ + (96 <= __GNUC_MINOR__)
#  define YY_ATTRIBUTE_PURE __attribute__ ((__pure__))
# else
#  define YY_ATTRIBUTE_PURE
# endif
#endif

#ifndef YY_ATTRIBUTE_UNUSED
# if defined __GNUC__ && 2 < __GNUC__ + (7 <= __GNUC_MINOR__)
#  define YY_ATTRIBUTE_UNUSED __attribute__ ((__unused__))
# else
#  define YY_ATTRIBUTE_UNUSED
# endif
#endif

/* Suppress unused-variable warnings by "using" E.  */
#if ! defined lint || defined __GNUC__
# define YY_USE(E) ((void) (E))
#else
# define YY_USE(E) /* empty */
#endif

/* Suppress an incorrect diagnostic about yylval being uninitialized.  */
#if defined __GNUC__ && ! defined __ICC && 406 <= __GNUC__ * 100 + __GNUC_MINOR__
# if __GNUC__ * 100 + __GNUC_MINOR__ < 407
#  define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN                           \
    _Pragma ("GCC diagnostic push")                                     \
    _Pragma ("GCC diagnostic ignored \"-Wuninitialized\"")
# else
#  define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN                           \
    _Pragma ("GCC diagnostic push")                                     \
    _Pragma ("GCC diagnostic ignored \"-Wuninitialized\"")              \
    _Pragma ("GCC diagnostic ignored \"-Wmaybe-uninitialized\"")
# endif
# define YY_IGNORE_MAYBE_UNINITIALIZED_END      \
    _Pragma ("GCC diagnostic pop")
#else
# define YY_INITIAL_VALUE(Value) Value
#endif
#ifndef YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_END
#endif
#ifndef YY_INITIAL_VALUE
# define YY_INITIAL_VALUE(Value) /* Nothing. */
#endif

#if defined __cplusplus && defined __GNUC__ && ! defined __ICC && 6 <= __GNUC__
# define YY_IGNORE_USELESS_CAST_BEGIN                          \
    _Pragma ("GCC diagnostic push")                            \
    _Pragma ("GCC diagnostic ignored \"-Wuseless-cast\"")
# define YY_IGNORE_USELESS_CAST_END            \
    _Pragma ("GCC diagnostic pop")
#endif
#ifndef YY_IGNORE_USELESS_CAST_BEGIN
# define YY_IGNORE_USELESS_CAST_BEGIN
# define YY_IGNORE_USELESS_CAST_END
#endif


#define YY_ASSERT(E) ((void) (0 && (E)))

#if !defined yyoverflow

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# ifdef YYSTACK_USE_ALLOCA
#  if YYSTACK_USE_ALLOCA
#   ifdef __GNUC__
#    define YYSTACK_ALLOC __builtin_alloca
#   elif defined __BUILTIN_VA_ARG_INCR
#    include <alloca.h> /* INFRINGES ON USER NAME SPACE */
#   elif defined _AIX
#    define YYSTACK_ALLOC __alloca
#   elif defined _MSC_VER
#    include <malloc.h> /* INFRINGES ON USER NAME SPACE */
#    define alloca _alloca
#   else
#    define YYSTACK_ALLOC alloca
#    if ! defined _ALLOCA_H && ! defined EXIT_SUCCESS
#     include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
      /* Use EXIT_SUCCESS as a witness for stdlib.h.  */
#     ifndef EXIT_SUCCESS
#      define EXIT_SUCCESS 0
#     endif
#    endif
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's 'empty if-body' warning.  */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (0)
#  ifndef YYSTACK_ALLOC_MAXIMUM
    /* The OS might guarantee only one guard page at the bottom of the stack,
       and a page size can be as small as 4096 bytes.  So we cannot safely
       invoke alloca (N) if N exceeds 4096.  Use a slightly smaller number
       to allow for a few compiler-allocated temporary stack slots.  */
#   define YYSTACK_ALLOC_MAXIMUM 4032 /* reasonable circa 2006 */
#  endif
# else
#  define YYSTACK_ALLOC YYMALLOC
#  define YYSTACK_FREE YYFREE
#  ifndef YYSTACK_ALLOC_MAXIMUM
#   define YYSTACK_ALLOC_MAXIMUM YYSIZE_MAXIMUM
#  endif
#  if (defined __cplusplus && ! defined EXIT_SUCCESS \
       && ! ((defined YYMALLOC || defined malloc) \
             && (defined YYFREE || defined free)))
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   ifndef EXIT_SUCCESS
#    define EXIT_SUCCESS 0
#   endif
#  endif
#  ifndef YYMALLOC
#   define YYMALLOC malloc
#   if ! defined malloc && ! defined EXIT_SUCCESS
void *malloc (YYSIZE_T); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
#  ifndef YYFREE
#   define YYFREE free
#   if ! defined free && ! defined EXIT_SUCCESS
void free (void *); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
# endif
#endif /* !defined yyoverflow */

#if (! defined yyoverflow \
     && (! defined __cplusplus \
         || (defined FITS_PARSER_YYSTYPE_IS_TRIVIAL && FITS_PARSER_YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  yy_state_t yyss_alloc;
  YYSTYPE yyvs_alloc;
};

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (YYSIZEOF (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (YYSIZEOF (yy_state_t) + YYSIZEOF (YYSTYPE)) \
      + YYSTACK_GAP_MAXIMUM)

# define YYCOPY_NEEDED 1

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack_alloc, Stack)                           \
    do                                                                  \
      {                                                                 \
        YYPTRDIFF_T yynewbytes;                                         \
        YYCOPY (&yyptr->Stack_alloc, Stack, yysize);                    \
        Stack = &yyptr->Stack_alloc;                                    \
        yynewbytes = yystacksize * YYSIZEOF (*Stack) + YYSTACK_GAP_MAXIMUM; \
        yyptr += yynewbytes / YYSIZEOF (*yyptr);                        \
      }                                                                 \
    while (0)

#endif

#if defined YYCOPY_NEEDED && YYCOPY_NEEDED
/* Copy COUNT objects from SRC to DST.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined __GNUC__ && 1 < __GNUC__
#   define YYCOPY(Dst, Src, Count) \
      __builtin_memcpy (Dst, Src, YY_CAST (YYSIZE_T, (Count)) * sizeof (*(Src)))
#  else
#   define YYCOPY(Dst, Src, Count)              \
      do                                        \
        {                                       \
          YYPTRDIFF_T yyi;                      \
          for (yyi = 0; yyi < (Count); yyi++)   \
            (Dst)[yyi] = (Src)[yyi];            \
        }                                       \
      while (0)
#  endif
# endif
#endif /* !YYCOPY_NEEDED */

/* YYFINAL -- State number of the termination state.  */
#define YYFINAL  2
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   1776

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  57
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  9
/* YYNRULES -- Number of rules.  */
#define YYNRULES  135
/* YYNSTATES -- Number of states.  */
#define YYNSTATES  322

/* YYMAXUTOK -- Last valid token kind.  */
#define YYMAXUTOK   292


/* YYTRANSLATE(TOKEN-NUM) -- Symbol number corresponding to TOKEN-NUM
   as returned by yylex, with out-of-bounds checking.  */
#define YYTRANSLATE(YYX)                                \
  (0 <= (YYX) && (YYX) <= YYMAXUTOK                     \
   ? YY_CAST (yysymbol_kind_t, yytranslate[YYX])        \
   : YYSYMBOL_YYUNDEF)

/* YYTRANSLATE[TOKEN-NUM] -- Symbol number corresponding to TOKEN-NUM
   as returned by yylex.  */
static const yytype_int8 yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
      53,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,    39,    43,     2,
      55,    56,    40,    37,    22,    38,     2,    41,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,    24,     2,
       2,    23,     2,    27,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,    50,     2,    54,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,    25,    42,    26,    32,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     1,     2,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    21,    28,    29,    30,
      31,    33,    34,    35,    36,    44,    45,    46,    47,    48,
      49,    51,    52
};

#if FITS_PARSER_YYDEBUG
/* YYRLINE[YYN] -- Source line where rule number YYN was defined.  */
static const yytype_int16 yyrline[] =
{
       0,   266,   266,   267,   270,   271,   277,   283,   289,   295,
     298,   300,   313,   315,   328,   339,   353,   357,   361,   365,
     367,   376,   379,   382,   391,   393,   395,   397,   399,   401,
     404,   408,   410,   412,   414,   423,   425,   427,   430,   433,
     436,   439,   442,   451,   460,   469,   472,   474,   476,   478,
     482,   486,   505,   524,   543,   554,   568,   617,   629,   660,
     774,   782,   885,   909,   911,   913,   915,   917,   919,   921,
     923,   925,   929,   931,   933,   942,   945,   948,   951,   954,
     957,   960,   963,   966,   969,   972,   975,   978,   981,   984,
     987,   990,   993,   996,   999,  1001,  1003,  1005,  1008,  1015,
    1032,  1045,  1058,  1069,  1085,  1109,  1137,  1174,  1178,  1182,
    1185,  1190,  1193,  1198,  1202,  1206,  1209,  1214,  1218,  1221,
    1225,  1227,  1229,  1231,  1233,  1235,  1237,  1241,  1244,  1246,
    1255,  1257,  1259,  1268,  1287,  1306
};
#endif

/** Accessing symbol of state STATE.  */
#define YY_ACCESSING_SYMBOL(State) YY_CAST (yysymbol_kind_t, yystos[State])

#if FITS_PARSER_YYDEBUG || 0
/* The user-facing name of the symbol whose (internal) number is
   YYSYMBOL.  No bounds checking.  */
static const char *yysymbol_name (yysymbol_kind_t yysymbol) YY_ATTRIBUTE_UNUSED;

/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "\"end of file\"", "error", "\"invalid token\"", "BOOLEAN", "LONG",
  "DOUBLE", "STRING", "BITSTR", "FUNCTION", "BFUNCTION", "IFUNCTION",
  "GTIFILTER", "GTIOVERLAP", "GTIFIND", "REGFILTER", "COLUMN", "BCOLUMN",
  "SCOLUMN", "BITCOL", "ROWREF", "NULLREF", "SNULLREF", "','", "'='",
  "':'", "'{'", "'}'", "'?'", "OR", "AND", "EQ", "NE", "'~'", "GT", "LT",
  "LTE", "GTE", "'+'", "'-'", "'%'", "'*'", "'/'", "'|'", "'&'", "XOR",
  "POWER", "NOT", "INTCAST", "FLTCAST", "UMINUS", "'['", "ACCUM", "DIFF",
  "'\\n'", "']'", "'('", "')'", "$accept", "lines", "line", "bvector",
  "vector", "expr", "bexpr", "bits", "sexpr", YY_NULLPTR
};

static const char *
yysymbol_name (yysymbol_kind_t yysymbol)
{
  return yytname[yysymbol];
}
#endif

#define YYPACT_NINF (-41)

#define yypact_value_is_default(Yyn) \
  ((Yyn) == YYPACT_NINF)

#define YYTABLE_NINF (-1)

#define yytable_value_is_error(Yyn) \
  0

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
static const yytype_int16 yypact[] =
{
     -41,   316,   -41,   -40,   -41,   -41,   -41,   -41,   -41,   369,
     423,   423,    -5,    15,    -4,    27,    36,    38,    40,    41,
     -41,   -41,   -41,   423,   423,   423,   423,   423,   423,   -41,
     423,   -41,    -7,    10,  1226,    81,  1646,    83,   -41,   -41,
     450,   116,   309,    12,   479,   185,   152,   222,  1593,  1673,
    1675,   -19,   -41,    13,   -18,   -41,     6,   423,   423,   423,
     423,  1593,  1673,  1684,    17,    17,    19,    24,    17,    19,
      17,    19,   710,  1253,  1611,   365,   423,   -41,   423,   -41,
     423,   423,   423,   423,   423,   423,   423,   423,   423,   423,
     423,   423,   423,   423,   423,   423,   423,   423,   -41,   423,
     423,   423,   423,   423,   423,   423,   -41,    -2,    -2,    -2,
      -2,    -2,    -2,    -2,    -2,    -2,   423,   -41,   423,   423,
     423,   423,   423,   423,   423,   -41,   423,   -41,   423,   -41,
     -41,   423,   -41,   423,   -41,   -41,   -41,   423,   423,   -41,
     423,   423,   -41,   423,   -41,  1455,  1478,  1501,  1524,   -41,
     -41,   -41,   -41,  1593,  1673,  1593,  1673,  1547,  1712,  1712,
    1712,  1726,  1726,  1726,  1726,   368,   368,   368,    28,    19,
      28,     5,     5,     5,     5,   851,  1570,   425,   260,   128,
     -20,    14,    14,    28,   876,    -2,    -2,   -25,   -25,   -25,
     -25,   -25,   -25,   -36,    24,    24,   901,   140,   140,    39,
      39,    39,    39,   -41,   508,   738,  1258,  1288,  1629,  1312,
    1638,   537,  1336,   566,  1360,   -41,   -41,   -41,   -41,   423,
     423,   -41,   423,   423,   423,   423,   -41,    24,   189,   423,
     -41,   423,   -41,   -41,   -41,   423,   -41,   423,   -41,    93,
     -41,   423,    94,   -41,   423,  1694,   926,  1694,  1673,  1694,
    1673,  1684,   951,   976,  1384,   766,   595,    79,   624,    80,
     653,   423,   -41,   423,   -41,   423,   -41,   423,   -41,   423,
     -41,   100,   101,   -41,   117,   118,   -41,  1001,  1026,  1051,
     794,  1408,    72,   111,    85,    99,   423,   -41,   423,   -41,
     423,   -41,   -41,   423,   -41,   129,   -41,   -41,  1076,  1101,
    1126,   682,   104,   423,   -41,   423,   -41,   423,   -41,   423,
     -41,   -41,  1151,  1176,  1201,  1432,   -41,   -41,   -41,   423,
     822,   -41
};

/* YYDEFACT[STATE-NUM] -- Default reduction number in state STATE-NUM.
   Performed when YYTABLE does not specify something else to do.  Zero
   means the default is an error.  */
static const yytype_uint8 yydefact[] =
{
       2,     0,     1,     0,    72,    31,    32,   127,    18,     0,
       0,     0,     0,     0,     0,     0,    33,    73,   128,    19,
      35,    36,   130,     0,     0,     0,     0,     0,     0,     4,
       0,     3,     0,     0,     0,     0,     0,     0,     9,    54,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   107,     0,     0,   113,     0,     0,     0,     0,
       0,    12,    10,     0,    46,    47,   125,    29,    68,    69,
      70,    71,     0,     0,     0,     0,     0,    17,     0,    16,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     5,     0,
       0,     0,     0,     0,     0,     0,     6,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     8,     0,     0,
       0,     0,     0,     0,     0,     7,     0,    59,     0,    55,
      58,     0,    57,     0,   100,   101,   102,     0,     0,   108,
       0,     0,   114,     0,   117,     0,     0,     0,     0,    48,
     126,    30,   131,    15,    11,    13,    14,     0,    86,    87,
      85,    81,    82,    84,    83,    38,    39,    37,    40,    49,
      41,    43,    42,    44,    45,     0,     0,     0,     0,    95,
      94,    96,    97,    50,     0,     0,     0,    75,    76,    79,
      77,    78,    80,    23,    22,    21,     0,    88,    89,    90,
      92,    93,    91,   132,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    34,    74,   129,    20,     0,
       0,    63,     0,     0,     0,     0,   120,    29,     0,     0,
      24,     0,    61,    56,   103,     0,   134,     0,    60,     0,
     109,     0,     0,   115,     0,    98,     0,    51,    53,    52,
      99,   133,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,    64,     0,   121,     0,    25,     0,   135,     0,
     104,     0,     0,   111,     0,     0,   118,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,    65,     0,   122,
       0,    26,    62,     0,   110,     0,   116,   119,     0,     0,
       0,     0,     0,     0,    66,     0,   123,     0,    27,     0,
     105,   112,     0,     0,     0,     0,    67,   124,    28,     0,
       0,   106
};

/* YYPGOTO[NTERM-NUM].  */
static const yytype_int16 yypgoto[] =
{
     -41,   -41,   -41,   -41,   -41,    -1,   170,    96,    30
};

/* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int8 yydefgoto[] =
{
       0,     1,    31,    32,    33,    48,    49,    46,    63
};

/* YYTABLE[YYPACT[STATE-NUM]] -- What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule whose
   number is the opposite.  If YYTABLE_NINF, syntax error.  */
static const yytype_int16 yytable[] =
{
      34,    51,    54,   138,   141,     8,   114,   115,    40,    44,
     102,   103,   113,    38,   116,    76,    19,   114,   115,    77,
     104,    53,    61,    64,    65,   116,    68,    70,   143,    72,
     105,    37,    78,    56,   131,   140,    79,   139,   142,    43,
      47,    50,   118,   119,   185,   120,   121,   122,   123,   124,
      96,    52,    55,   186,   104,    97,   145,   146,   147,   148,
      75,    57,   144,    58,   105,    59,    60,    97,   132,   105,
      93,    94,    95,    96,   116,   153,   124,   155,    97,   157,
     158,   159,   160,   161,   162,   163,   164,   165,   166,   167,
     168,   170,   171,   172,   173,   174,   175,    36,   176,   257,
     259,   271,   274,   183,   184,    42,   282,   283,    99,   100,
     101,   102,   103,   118,   119,   196,   120,   121,   122,   123,
     124,   104,    67,   284,   285,   204,    74,   205,   294,   178,
     207,   105,   209,   295,   106,   302,   125,   211,   128,   212,
     213,   296,   214,    99,   100,   101,   102,   103,   197,   198,
     199,   200,   201,   202,   203,   297,   104,   101,   102,   103,
     311,   208,     0,     0,     0,     0,   105,   210,   104,     0,
       0,    35,   129,   120,   121,   122,   123,   124,   105,    41,
      45,     0,   107,   108,     0,   109,   110,   111,   112,   113,
       0,     0,     0,    62,   114,   115,    66,    69,    71,     0,
      73,     0,   116,   187,   188,   189,   190,   191,   192,   193,
     194,   195,    99,   100,   101,   102,   103,     0,   245,   246,
       0,   247,   249,     0,   252,   104,   113,     0,   253,     0,
     254,   114,   115,     0,   255,   105,   256,     0,     0,   116,
     258,   135,     0,   260,     0,   151,   154,     0,   156,     0,
       0,     0,   118,   119,   251,   120,   121,   122,   123,   124,
     277,   169,   278,     0,   279,     0,   280,     0,   281,   177,
     179,   180,   181,   182,     0,     0,     0,     0,   136,     0,
       0,   227,   228,     0,   224,   298,     0,   299,     0,   300,
     118,   119,   301,   120,   121,   122,   123,   124,   206,     0,
       0,     0,   312,     0,   313,     0,   314,     0,   315,     0,
       0,     0,     0,     0,     0,     0,     2,     3,   320,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    21,    22,     0,   107,
     108,    23,   109,   110,   111,   112,   113,     0,     0,     0,
       0,   114,   115,    24,    25,     0,     0,     0,     0,   116,
       0,     0,    26,    27,    28,   130,     0,     0,     0,    29,
       0,    30,     4,     5,     6,     7,     8,     9,    10,    11,
      12,    13,    14,    15,    16,    17,    18,    19,    20,    21,
      22,     0,   248,   250,    23,   118,   119,     0,   120,   121,
     122,   123,   124,     0,     0,     0,    24,    25,    91,    92,
      93,    94,    95,    96,     0,    26,    27,    28,    97,     0,
       0,   152,     0,     0,    30,    39,     4,     5,     6,     7,
       8,     9,    10,    11,    12,    13,    14,    15,    16,    17,
      18,    19,    20,    21,    22,     0,     0,     0,    23,   223,
       0,     0,    99,   100,   101,   102,   103,     0,     0,     0,
      24,    25,     0,     0,     0,   104,     0,     0,     0,    26,
      27,    28,   126,    80,     0,   105,     0,     0,    30,     0,
      81,    82,    83,    84,    85,    86,    87,    88,    89,    90,
      91,    92,    93,    94,    95,    96,     0,     0,     0,     0,
      97,   133,    80,     0,     0,     0,   127,     0,     0,    81,
      82,    83,    84,    85,    86,    87,    88,    89,    90,    91,
      92,    93,    94,    95,    96,     0,     0,     0,     0,    97,
     231,    80,     0,     0,     0,   134,     0,     0,    81,    82,
      83,    84,    85,    86,    87,    88,    89,    90,    91,    92,
      93,    94,    95,    96,     0,     0,     0,     0,    97,   239,
      80,     0,     0,     0,   232,     0,     0,    81,    82,    83,
      84,    85,    86,    87,    88,    89,    90,    91,    92,    93,
      94,    95,    96,     0,     0,     0,     0,    97,   242,    80,
       0,     0,     0,   240,     0,     0,    81,    82,    83,    84,
      85,    86,    87,    88,    89,    90,    91,    92,    93,    94,
      95,    96,     0,     0,     0,     0,    97,   269,    80,     0,
       0,     0,   243,     0,     0,    81,    82,    83,    84,    85,
      86,    87,    88,    89,    90,    91,    92,    93,    94,    95,
      96,     0,     0,     0,     0,    97,   272,    80,     0,     0,
       0,   270,     0,     0,    81,    82,    83,    84,    85,    86,
      87,    88,    89,    90,    91,    92,    93,    94,    95,    96,
       0,     0,     0,     0,    97,   275,    80,     0,     0,     0,
     273,     0,     0,    81,    82,    83,    84,    85,    86,    87,
      88,    89,    90,    91,    92,    93,    94,    95,    96,     0,
       0,     0,     0,    97,   309,    80,     0,     0,     0,   276,
       0,     0,    81,    82,    83,    84,    85,    86,    87,    88,
      89,    90,    91,    92,    93,    94,    95,    96,     0,     0,
       0,     0,    97,    80,     0,     0,     0,     0,   310,     0,
      81,    82,    83,    84,    85,    86,    87,    88,    89,    90,
      91,    92,    93,    94,    95,    96,     0,     0,     0,     0,
      97,    80,     0,     0,     0,     0,   149,     0,    81,    82,
      83,    84,    85,    86,    87,    88,    89,    90,    91,    92,
      93,    94,    95,    96,     0,     0,     0,     0,    97,    80,
       0,     0,     0,     0,   233,     0,    81,    82,    83,    84,
      85,    86,    87,    88,    89,    90,    91,    92,    93,    94,
      95,    96,     0,     0,     0,     0,    97,    80,     0,     0,
       0,     0,   268,     0,    81,    82,    83,    84,    85,    86,
      87,    88,    89,    90,    91,    92,    93,    94,    95,    96,
       0,     0,     0,     0,    97,    80,     0,     0,     0,     0,
     292,     0,    81,    82,    83,    84,    85,    86,    87,    88,
      89,    90,    91,    92,    93,    94,    95,    96,     0,     0,
       0,     0,    97,   220,    80,     0,     0,     0,   321,     0,
       0,    81,    82,    83,    84,    85,    86,    87,    88,    89,
      90,    91,    92,    93,    94,    95,    96,     0,   225,    80,
       0,    97,     0,     0,     0,   221,    81,    82,    83,    84,
      85,    86,    87,    88,    89,    90,    91,    92,    93,    94,
      95,    96,     0,   229,    80,     0,    97,     0,     0,     0,
     226,    81,    82,    83,    84,    85,    86,    87,    88,    89,
      90,    91,    92,    93,    94,    95,    96,     0,   261,    80,
       0,    97,     0,     0,     0,   230,    81,    82,    83,    84,
      85,    86,    87,    88,    89,    90,    91,    92,    93,    94,
      95,    96,     0,   263,    80,     0,    97,     0,     0,     0,
     262,    81,    82,    83,    84,    85,    86,    87,    88,    89,
      90,    91,    92,    93,    94,    95,    96,     0,   265,    80,
       0,    97,     0,     0,     0,   264,    81,    82,    83,    84,
      85,    86,    87,    88,    89,    90,    91,    92,    93,    94,
      95,    96,     0,   286,    80,     0,    97,     0,     0,     0,
     266,    81,    82,    83,    84,    85,    86,    87,    88,    89,
      90,    91,    92,    93,    94,    95,    96,     0,   288,    80,
       0,    97,     0,     0,     0,   287,    81,    82,    83,    84,
      85,    86,    87,    88,    89,    90,    91,    92,    93,    94,
      95,    96,     0,   290,    80,     0,    97,     0,     0,     0,
     289,    81,    82,    83,    84,    85,    86,    87,    88,    89,
      90,    91,    92,    93,    94,    95,    96,     0,   303,    80,
       0,    97,     0,     0,     0,   291,    81,    82,    83,    84,
      85,    86,    87,    88,    89,    90,    91,    92,    93,    94,
      95,    96,     0,   305,    80,     0,    97,     0,     0,     0,
     304,    81,    82,    83,    84,    85,    86,    87,    88,    89,
      90,    91,    92,    93,    94,    95,    96,     0,   307,    80,
       0,    97,     0,     0,     0,   306,    81,    82,    83,    84,
      85,    86,    87,    88,    89,    90,    91,    92,    93,    94,
      95,    96,     0,     0,    80,     0,    97,     0,     0,     0,
     308,    81,    82,    83,    84,    85,    86,    87,    88,    89,
      90,    91,    92,    93,    94,    95,    96,     0,     0,    80,
       0,    97,     0,     0,     0,   316,    81,    82,    83,    84,
      85,    86,    87,    88,    89,    90,    91,    92,    93,    94,
      95,    96,     0,     0,    80,     0,    97,     0,     0,     0,
     317,    81,    82,    83,    84,    85,    86,    87,    88,    89,
      90,    91,    92,    93,    94,    95,    96,     0,     0,    80,
       0,    97,     0,     0,     0,   318,    81,    82,    83,    84,
      85,    86,    87,    88,    89,    90,    91,    92,    93,    94,
      95,    96,     0,     0,     0,     0,    97,     0,     0,    98,
      99,   100,   101,   102,   103,    99,   100,   101,   102,   103,
       0,     0,     0,   104,     0,     0,     0,     0,   104,     0,
       0,     0,     0,   105,     0,     0,     0,     0,   105,   150,
     235,    80,     0,     0,   234,     0,     0,     0,    81,    82,
      83,    84,    85,    86,    87,    88,    89,    90,    91,    92,
      93,    94,    95,    96,   237,    80,     0,     0,    97,     0,
       0,     0,    81,    82,    83,    84,    85,    86,    87,    88,
      89,    90,    91,    92,    93,    94,    95,    96,   241,    80,
       0,     0,    97,     0,     0,     0,    81,    82,    83,    84,
      85,    86,    87,    88,    89,    90,    91,    92,    93,    94,
      95,    96,   244,    80,     0,     0,    97,     0,     0,     0,
      81,    82,    83,    84,    85,    86,    87,    88,    89,    90,
      91,    92,    93,    94,    95,    96,   267,    80,     0,     0,
      97,     0,     0,     0,    81,    82,    83,    84,    85,    86,
      87,    88,    89,    90,    91,    92,    93,    94,    95,    96,
     293,    80,     0,     0,    97,     0,     0,     0,    81,    82,
      83,    84,    85,    86,    87,    88,    89,    90,    91,    92,
      93,    94,    95,    96,   319,    80,     0,     0,    97,     0,
       0,     0,    81,    82,    83,    84,    85,    86,    87,    88,
      89,    90,    91,    92,    93,    94,    95,    96,    80,     0,
       0,   215,    97,     0,     0,    81,    82,    83,    84,    85,
      86,    87,    88,    89,    90,    91,    92,    93,    94,    95,
      96,    80,     0,     0,   216,    97,     0,     0,    81,    82,
      83,    84,    85,    86,    87,    88,    89,    90,    91,    92,
      93,    94,    95,    96,    80,     0,     0,   217,    97,     0,
       0,    81,    82,    83,    84,    85,    86,    87,    88,    89,
      90,    91,    92,    93,    94,    95,    96,    80,     0,     0,
     218,    97,     0,     0,    81,    82,    83,    84,    85,    86,
      87,    88,    89,    90,    91,    92,    93,    94,    95,    96,
      80,   219,     0,     0,    97,     0,     0,    81,    82,    83,
      84,    85,    86,    87,    88,    89,    90,    91,    92,    93,
      94,    95,    96,    80,   222,     0,     0,    97,     0,     0,
      81,    82,    83,    84,    85,    86,    87,    88,    89,    90,
      91,    92,    93,    94,    95,    96,    80,     0,     0,     0,
      97,     0,     0,    81,    82,    83,    84,    85,    86,    87,
      88,    89,    90,    91,    92,    93,    94,    95,    96,     0,
       0,   107,   108,    97,   109,   110,   111,   112,   113,     0,
       0,     0,     0,   114,   115,     0,     0,     0,     0,   118,
     119,   116,   120,   121,   122,   123,   124,   151,   118,   119,
       0,   120,   121,   122,   123,   124,   107,   108,     0,   109,
     110,   111,   112,   113,     0,   236,     0,     0,   114,   115,
       0,     0,     0,     0,   238,     0,   116,   137,     0,   117,
      99,   100,   101,   102,   103,   118,   119,     0,   120,   121,
     122,   123,   124,   104,   118,   119,     0,   120,   121,   122,
     123,   124,     0,   105,    81,    82,    83,    84,    85,    86,
      87,    88,    89,    90,    91,    92,    93,    94,    95,    96,
       0,     0,     0,     0,    97,    84,    85,    86,    87,    88,
      89,    90,    91,    92,    93,    94,    95,    96,     0,     0,
       0,     0,    97,    88,    89,    90,    91,    92,    93,    94,
      95,    96,     0,     0,     0,     0,    97
};

static const yytype_int16 yycheck[] =
{
       1,     6,     6,    22,    22,     7,    42,    43,     9,    10,
      30,    31,    37,    53,    50,    22,    18,    42,    43,    26,
      40,     6,    23,    24,    25,    50,    27,    28,    22,    30,
      50,     1,    22,     6,    22,    22,    26,    56,    56,     9,
      10,    11,    30,    31,    46,    33,    34,    35,    36,    37,
      45,    56,    56,    55,    40,    50,    57,    58,    59,    60,
      30,    25,    56,    25,    50,    25,    25,    50,    56,    50,
      42,    43,    44,    45,    50,    76,    37,    78,    50,    80,
      81,    82,    83,    84,    85,    86,    87,    88,    89,    90,
      91,    92,    93,    94,    95,    96,    97,     1,    99,     6,
       6,    22,    22,   104,   105,     9,     6,     6,    27,    28,
      29,    30,    31,    30,    31,   116,    33,    34,    35,    36,
      37,    40,    26,     6,     6,   126,    30,   128,    56,    99,
     131,    50,   133,    22,    53,     6,    53,   138,    22,   140,
     141,    56,   143,    27,    28,    29,    30,    31,   118,   119,
     120,   121,   122,   123,   124,    56,    40,    29,    30,    31,
      56,   131,    -1,    -1,    -1,    -1,    50,   137,    40,    -1,
      -1,     1,    56,    33,    34,    35,    36,    37,    50,     9,
      10,    -1,    30,    31,    -1,    33,    34,    35,    36,    37,
      -1,    -1,    -1,    23,    42,    43,    26,    27,    28,    -1,
      30,    -1,    50,   107,   108,   109,   110,   111,   112,   113,
     114,   115,    27,    28,    29,    30,    31,    -1,   219,   220,
      -1,   222,   223,    -1,   225,    40,    37,    -1,   229,    -1,
     231,    42,    43,    -1,   235,    50,   237,    -1,    -1,    50,
     241,    56,    -1,   244,    -1,    56,    76,    -1,    78,    -1,
      -1,    -1,    30,    31,   224,    33,    34,    35,    36,    37,
     261,    91,   263,    -1,   265,    -1,   267,    -1,   269,    99,
     100,   101,   102,   103,    -1,    -1,    -1,    -1,    56,    -1,
      -1,   185,   186,    -1,    24,   286,    -1,   288,    -1,   290,
      30,    31,   293,    33,    34,    35,    36,    37,   128,    -1,
      -1,    -1,   303,    -1,   305,    -1,   307,    -1,   309,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,     0,     1,   319,     3,
       4,     5,     6,     7,     8,     9,    10,    11,    12,    13,
      14,    15,    16,    17,    18,    19,    20,    21,    -1,    30,
      31,    25,    33,    34,    35,    36,    37,    -1,    -1,    -1,
      -1,    42,    43,    37,    38,    -1,    -1,    -1,    -1,    50,
      -1,    -1,    46,    47,    48,    56,    -1,    -1,    -1,    53,
      -1,    55,     3,     4,     5,     6,     7,     8,     9,    10,
      11,    12,    13,    14,    15,    16,    17,    18,    19,    20,
      21,    -1,   222,   223,    25,    30,    31,    -1,    33,    34,
      35,    36,    37,    -1,    -1,    -1,    37,    38,    40,    41,
      42,    43,    44,    45,    -1,    46,    47,    48,    50,    -1,
      -1,    56,    -1,    -1,    55,    56,     3,     4,     5,     6,
       7,     8,     9,    10,    11,    12,    13,    14,    15,    16,
      17,    18,    19,    20,    21,    -1,    -1,    -1,    25,    24,
      -1,    -1,    27,    28,    29,    30,    31,    -1,    -1,    -1,
      37,    38,    -1,    -1,    -1,    40,    -1,    -1,    -1,    46,
      47,    48,    22,    23,    -1,    50,    -1,    -1,    55,    -1,
      30,    31,    32,    33,    34,    35,    36,    37,    38,    39,
      40,    41,    42,    43,    44,    45,    -1,    -1,    -1,    -1,
      50,    22,    23,    -1,    -1,    -1,    56,    -1,    -1,    30,
      31,    32,    33,    34,    35,    36,    37,    38,    39,    40,
      41,    42,    43,    44,    45,    -1,    -1,    -1,    -1,    50,
      22,    23,    -1,    -1,    -1,    56,    -1,    -1,    30,    31,
      32,    33,    34,    35,    36,    37,    38,    39,    40,    41,
      42,    43,    44,    45,    -1,    -1,    -1,    -1,    50,    22,
      23,    -1,    -1,    -1,    56,    -1,    -1,    30,    31,    32,
      33,    34,    35,    36,    37,    38,    39,    40,    41,    42,
      43,    44,    45,    -1,    -1,    -1,    -1,    50,    22,    23,
      -1,    -1,    -1,    56,    -1,    -1,    30,    31,    32,    33,
      34,    35,    36,    37,    38,    39,    40,    41,    42,    43,
      44,    45,    -1,    -1,    -1,    -1,    50,    22,    23,    -1,
      -1,    -1,    56,    -1,    -1,    30,    31,    32,    33,    34,
      35,    36,    37,    38,    39,    40,    41,    42,    43,    44,
      45,    -1,    -1,    -1,    -1,    50,    22,    23,    -1,    -1,
      -1,    56,    -1,    -1,    30,    31,    32,    33,    34,    35,
      36,    37,    38,    39,    40,    41,    42,    43,    44,    45,
      -1,    -1,    -1,    -1,    50,    22,    23,    -1,    -1,    -1,
      56,    -1,    -1,    30,    31,    32,    33,    34,    35,    36,
      37,    38,    39,    40,    41,    42,    43,    44,    45,    -1,
      -1,    -1,    -1,    50,    22,    23,    -1,    -1,    -1,    56,
      -1,    -1,    30,    31,    32,    33,    34,    35,    36,    37,
      38,    39,    40,    41,    42,    43,    44,    45,    -1,    -1,
      -1,    -1,    50,    23,    -1,    -1,    -1,    -1,    56,    -1,
      30,    31,    32,    33,    34,    35,    36,    37,    38,    39,
      40,    41,    42,    43,    44,    45,    -1,    -1,    -1,    -1,
      50,    23,    -1,    -1,    -1,    -1,    56,    -1,    30,    31,
      32,    33,    34,    35,    36,    37,    38,    39,    40,    41,
      42,    43,    44,    45,    -1,    -1,    -1,    -1,    50,    23,
      -1,    -1,    -1,    -1,    56,    -1,    30,    31,    32,    33,
      34,    35,    36,    37,    38,    39,    40,    41,    42,    43,
      44,    45,    -1,    -1,    -1,    -1,    50,    23,    -1,    -1,
      -1,    -1,    56,    -1,    30,    31,    32,    33,    34,    35,
      36,    37,    38,    39,    40,    41,    42,    43,    44,    45,
      -1,    -1,    -1,    -1,    50,    23,    -1,    -1,    -1,    -1,
      56,    -1,    30,    31,    32,    33,    34,    35,    36,    37,
      38,    39,    40,    41,    42,    43,    44,    45,    -1,    -1,
      -1,    -1,    50,    22,    23,    -1,    -1,    -1,    56,    -1,
      -1,    30,    31,    32,    33,    34,    35,    36,    37,    38,
      39,    40,    41,    42,    43,    44,    45,    -1,    22,    23,
      -1,    50,    -1,    -1,    -1,    54,    30,    31,    32,    33,
      34,    35,    36,    37,    38,    39,    40,    41,    42,    43,
      44,    45,    -1,    22,    23,    -1,    50,    -1,    -1,    -1,
      54,    30,    31,    32,    33,    34,    35,    36,    37,    38,
      39,    40,    41,    42,    43,    44,    45,    -1,    22,    23,
      -1,    50,    -1,    -1,    -1,    54,    30,    31,    32,    33,
      34,    35,    36,    37,    38,    39,    40,    41,    42,    43,
      44,    45,    -1,    22,    23,    -1,    50,    -1,    -1,    -1,
      54,    30,    31,    32,    33,    34,    35,    36,    37,    38,
      39,    40,    41,    42,    43,    44,    45,    -1,    22,    23,
      -1,    50,    -1,    -1,    -1,    54,    30,    31,    32,    33,
      34,    35,    36,    37,    38,    39,    40,    41,    42,    43,
      44,    45,    -1,    22,    23,    -1,    50,    -1,    -1,    -1,
      54,    30,    31,    32,    33,    34,    35,    36,    37,    38,
      39,    40,    41,    42,    43,    44,    45,    -1,    22,    23,
      -1,    50,    -1,    -1,    -1,    54,    30,    31,    32,    33,
      34,    35,    36,    37,    38,    39,    40,    41,    42,    43,
      44,    45,    -1,    22,    23,    -1,    50,    -1,    -1,    -1,
      54,    30,    31,    32,    33,    34,    35,    36,    37,    38,
      39,    40,    41,    42,    43,    44,    45,    -1,    22,    23,
      -1,    50,    -1,    -1,    -1,    54,    30,    31,    32,    33,
      34,    35,    36,    37,    38,    39,    40,    41,    42,    43,
      44,    45,    -1,    22,    23,    -1,    50,    -1,    -1,    -1,
      54,    30,    31,    32,    33,    34,    35,    36,    37,    38,
      39,    40,    41,    42,    43,    44,    45,    -1,    22,    23,
      -1,    50,    -1,    -1,    -1,    54,    30,    31,    32,    33,
      34,    35,    36,    37,    38,    39,    40,    41,    42,    43,
      44,    45,    -1,    -1,    23,    -1,    50,    -1,    -1,    -1,
      54,    30,    31,    32,    33,    34,    35,    36,    37,    38,
      39,    40,    41,    42,    43,    44,    45,    -1,    -1,    23,
      -1,    50,    -1,    -1,    -1,    54,    30,    31,    32,    33,
      34,    35,    36,    37,    38,    39,    40,    41,    42,    43,
      44,    45,    -1,    -1,    23,    -1,    50,    -1,    -1,    -1,
      54,    30,    31,    32,    33,    34,    35,    36,    37,    38,
      39,    40,    41,    42,    43,    44,    45,    -1,    -1,    23,
      -1,    50,    -1,    -1,    -1,    54,    30,    31,    32,    33,
      34,    35,    36,    37,    38,    39,    40,    41,    42,    43,
      44,    45,    -1,    -1,    -1,    -1,    50,    -1,    -1,    53,
      27,    28,    29,    30,    31,    27,    28,    29,    30,    31,
      -1,    -1,    -1,    40,    -1,    -1,    -1,    -1,    40,    -1,
      -1,    -1,    -1,    50,    -1,    -1,    -1,    -1,    50,    56,
      22,    23,    -1,    -1,    56,    -1,    -1,    -1,    30,    31,
      32,    33,    34,    35,    36,    37,    38,    39,    40,    41,
      42,    43,    44,    45,    22,    23,    -1,    -1,    50,    -1,
      -1,    -1,    30,    31,    32,    33,    34,    35,    36,    37,
      38,    39,    40,    41,    42,    43,    44,    45,    22,    23,
      -1,    -1,    50,    -1,    -1,    -1,    30,    31,    32,    33,
      34,    35,    36,    37,    38,    39,    40,    41,    42,    43,
      44,    45,    22,    23,    -1,    -1,    50,    -1,    -1,    -1,
      30,    31,    32,    33,    34,    35,    36,    37,    38,    39,
      40,    41,    42,    43,    44,    45,    22,    23,    -1,    -1,
      50,    -1,    -1,    -1,    30,    31,    32,    33,    34,    35,
      36,    37,    38,    39,    40,    41,    42,    43,    44,    45,
      22,    23,    -1,    -1,    50,    -1,    -1,    -1,    30,    31,
      32,    33,    34,    35,    36,    37,    38,    39,    40,    41,
      42,    43,    44,    45,    22,    23,    -1,    -1,    50,    -1,
      -1,    -1,    30,    31,    32,    33,    34,    35,    36,    37,
      38,    39,    40,    41,    42,    43,    44,    45,    23,    -1,
      -1,    26,    50,    -1,    -1,    30,    31,    32,    33,    34,
      35,    36,    37,    38,    39,    40,    41,    42,    43,    44,
      45,    23,    -1,    -1,    26,    50,    -1,    -1,    30,    31,
      32,    33,    34,    35,    36,    37,    38,    39,    40,    41,
      42,    43,    44,    45,    23,    -1,    -1,    26,    50,    -1,
      -1,    30,    31,    32,    33,    34,    35,    36,    37,    38,
      39,    40,    41,    42,    43,    44,    45,    23,    -1,    -1,
      26,    50,    -1,    -1,    30,    31,    32,    33,    34,    35,
      36,    37,    38,    39,    40,    41,    42,    43,    44,    45,
      23,    24,    -1,    -1,    50,    -1,    -1,    30,    31,    32,
      33,    34,    35,    36,    37,    38,    39,    40,    41,    42,
      43,    44,    45,    23,    24,    -1,    -1,    50,    -1,    -1,
      30,    31,    32,    33,    34,    35,    36,    37,    38,    39,
      40,    41,    42,    43,    44,    45,    23,    -1,    -1,    -1,
      50,    -1,    -1,    30,    31,    32,    33,    34,    35,    36,
      37,    38,    39,    40,    41,    42,    43,    44,    45,    -1,
      -1,    30,    31,    50,    33,    34,    35,    36,    37,    -1,
      -1,    -1,    -1,    42,    43,    -1,    -1,    -1,    -1,    30,
      31,    50,    33,    34,    35,    36,    37,    56,    30,    31,
      -1,    33,    34,    35,    36,    37,    30,    31,    -1,    33,
      34,    35,    36,    37,    -1,    56,    -1,    -1,    42,    43,
      -1,    -1,    -1,    -1,    56,    -1,    50,    22,    -1,    53,
      27,    28,    29,    30,    31,    30,    31,    -1,    33,    34,
      35,    36,    37,    40,    30,    31,    -1,    33,    34,    35,
      36,    37,    -1,    50,    30,    31,    32,    33,    34,    35,
      36,    37,    38,    39,    40,    41,    42,    43,    44,    45,
      -1,    -1,    -1,    -1,    50,    33,    34,    35,    36,    37,
      38,    39,    40,    41,    42,    43,    44,    45,    -1,    -1,
      -1,    -1,    50,    37,    38,    39,    40,    41,    42,    43,
      44,    45,    -1,    -1,    -1,    -1,    50
};

/* YYSTOS[STATE-NUM] -- The symbol kind of the accessing symbol of
   state STATE-NUM.  */
static const yytype_int8 yystos[] =
{
       0,    58,     0,     1,     3,     4,     5,     6,     7,     8,
       9,    10,    11,    12,    13,    14,    15,    16,    17,    18,
      19,    20,    21,    25,    37,    38,    46,    47,    48,    53,
      55,    59,    60,    61,    62,    63,    64,    65,    53,    56,
      62,    63,    64,    65,    62,    63,    64,    65,    62,    63,
      65,     6,    56,     6,     6,    56,     6,    25,    25,    25,
      25,    62,    63,    65,    62,    62,    63,    64,    62,    63,
      62,    63,    62,    63,    64,    65,    22,    26,    22,    26,
      23,    30,    31,    32,    33,    34,    35,    36,    37,    38,
      39,    40,    41,    42,    43,    44,    45,    50,    53,    27,
      28,    29,    30,    31,    40,    50,    53,    30,    31,    33,
      34,    35,    36,    37,    42,    43,    50,    53,    30,    31,
      33,    34,    35,    36,    37,    53,    22,    56,    22,    56,
      56,    22,    56,    22,    56,    56,    56,    22,    22,    56,
      22,    22,    56,    22,    56,    62,    62,    62,    62,    56,
      56,    56,    56,    62,    63,    62,    63,    62,    62,    62,
      62,    62,    62,    62,    62,    62,    62,    62,    62,    63,
      62,    62,    62,    62,    62,    62,    62,    63,    65,    63,
      63,    63,    63,    62,    62,    46,    55,    64,    64,    64,
      64,    64,    64,    64,    64,    64,    62,    65,    65,    65,
      65,    65,    65,    65,    62,    62,    63,    62,    65,    62,
      65,    62,    62,    62,    62,    26,    26,    26,    26,    24,
      22,    54,    24,    24,    24,    22,    54,    64,    64,    22,
      54,    22,    56,    56,    56,    22,    56,    22,    56,    22,
      56,    22,    22,    56,    22,    62,    62,    62,    63,    62,
      63,    65,    62,    62,    62,    62,    62,     6,    62,     6,
      62,    22,    54,    22,    54,    22,    54,    22,    56,    22,
      56,    22,    22,    56,    22,    22,    56,    62,    62,    62,
      62,    62,     6,     6,     6,     6,    22,    54,    22,    54,
      22,    54,    56,    22,    56,    22,    56,    56,    62,    62,
      62,    62,     6,    22,    54,    22,    54,    22,    54,    22,
      56,    56,    62,    62,    62,    62,    54,    54,    54,    22,
      62,    56
};

/* YYR1[RULE-NUM] -- Symbol kind of the left-hand side of rule RULE-NUM.  */
static const yytype_int8 yyr1[] =
{
       0,    57,    58,    58,    59,    59,    59,    59,    59,    59,
      60,    60,    61,    61,    61,    61,    62,    63,    64,    64,
      64,    64,    64,    64,    64,    64,    64,    64,    64,    64,
      64,    62,    62,    62,    62,    62,    62,    62,    62,    62,
      62,    62,    62,    62,    62,    62,    62,    62,    62,    62,
      62,    62,    62,    62,    62,    62,    62,    62,    62,    62,
      62,    62,    62,    62,    62,    62,    62,    62,    62,    62,
      62,    62,    63,    63,    63,    63,    63,    63,    63,    63,
      63,    63,    63,    63,    63,    63,    63,    63,    63,    63,
      63,    63,    63,    63,    63,    63,    63,    63,    63,    63,
      63,    63,    63,    63,    63,    63,    63,    63,    63,    63,
      63,    63,    63,    63,    63,    63,    63,    63,    63,    63,
      63,    63,    63,    63,    63,    63,    63,    65,    65,    65,
      65,    65,    65,    65,    65,    65
};

/* YYR2[RULE-NUM] -- Number of symbols on the right-hand side of rule RULE-NUM.  */
static const yytype_int8 yyr2[] =
{
       0,     2,     0,     2,     1,     2,     2,     2,     2,     2,
       2,     3,     2,     3,     3,     3,     2,     2,     1,     1,
       4,     3,     3,     3,     4,     6,     8,    10,    12,     2,
       3,     1,     1,     1,     4,     1,     1,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     2,     2,     3,     3,
       3,     5,     5,     5,     2,     3,     5,     3,     3,     3,
       5,     5,     9,     4,     6,     8,    10,    12,     2,     2,
       2,     2,     1,     1,     4,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     5,     5,
       3,     3,     3,     5,     7,    11,    15,     2,     3,     5,
       9,     7,    11,     2,     3,     5,     9,     3,     7,     9,
       4,     6,     8,    10,    12,     2,     3,     1,     1,     4,
       1,     3,     3,     5,     5,     7
};


enum { YYENOMEM = -2 };

#define yyerrok         (yyerrstatus = 0)
#define yyclearin       (yychar = FITS_PARSER_YYEMPTY)

#define YYACCEPT        goto yyacceptlab
#define YYABORT         goto yyabortlab
#define YYERROR         goto yyerrorlab
#define YYNOMEM         goto yyexhaustedlab


#define YYRECOVERING()  (!!yyerrstatus)

#define YYBACKUP(Token, Value)                                    \
  do                                                              \
    if (yychar == FITS_PARSER_YYEMPTY)                                        \
      {                                                           \
        yychar = (Token);                                         \
        yylval = (Value);                                         \
        YYPOPSTACK (yylen);                                       \
        yystate = *yyssp;                                         \
        goto yybackup;                                            \
      }                                                           \
    else                                                          \
      {                                                           \
        yyerror (scanner, lParse, YY_("syntax error: cannot back up")); \
        YYERROR;                                                  \
      }                                                           \
  while (0)

/* Backward compatibility with an undocumented macro.
   Use FITS_PARSER_YYerror or FITS_PARSER_YYUNDEF. */
#define YYERRCODE FITS_PARSER_YYUNDEF


/* Enable debugging if requested.  */
#if FITS_PARSER_YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)                        \
do {                                            \
  if (yydebug)                                  \
    YYFPRINTF Args;                             \
} while (0)




# define YY_SYMBOL_PRINT(Title, Kind, Value, Location)                    \
do {                                                                      \
  if (yydebug)                                                            \
    {                                                                     \
      YYFPRINTF (stderr, "%s ", Title);                                   \
      yy_symbol_print (stderr,                                            \
                  Kind, Value, scanner, lParse); \
      YYFPRINTF (stderr, "\n");                                           \
    }                                                                     \
} while (0)


/*-----------------------------------.
| Print this symbol's value on YYO.  |
`-----------------------------------*/

static void
yy_symbol_value_print (FILE *yyo,
                       yysymbol_kind_t yykind, YYSTYPE const * const yyvaluep, yyscan_t scanner, ParseData *lParse)
{
  FILE *yyoutput = yyo;
  YY_USE (yyoutput);
  YY_USE (scanner);
  YY_USE (lParse);
  if (!yyvaluep)
    return;
  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  YY_USE (yykind);
  YY_IGNORE_MAYBE_UNINITIALIZED_END
}


/*---------------------------.
| Print this symbol on YYO.  |
`---------------------------*/

static void
yy_symbol_print (FILE *yyo,
                 yysymbol_kind_t yykind, YYSTYPE const * const yyvaluep, yyscan_t scanner, ParseData *lParse)
{
  YYFPRINTF (yyo, "%s %s (",
             yykind < YYNTOKENS ? "token" : "nterm", yysymbol_name (yykind));

  yy_symbol_value_print (yyo, yykind, yyvaluep, scanner, lParse);
  YYFPRINTF (yyo, ")");
}

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

static void
yy_stack_print (yy_state_t *yybottom, yy_state_t *yytop)
{
  YYFPRINTF (stderr, "Stack now");
  for (; yybottom <= yytop; yybottom++)
    {
      int yybot = *yybottom;
      YYFPRINTF (stderr, " %d", yybot);
    }
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)                            \
do {                                                            \
  if (yydebug)                                                  \
    yy_stack_print ((Bottom), (Top));                           \
} while (0)


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

static void
yy_reduce_print (yy_state_t *yyssp, YYSTYPE *yyvsp,
                 int yyrule, yyscan_t scanner, ParseData *lParse)
{
  int yylno = yyrline[yyrule];
  int yynrhs = yyr2[yyrule];
  int yyi;
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %d):\n",
             yyrule - 1, yylno);
  /* The symbols being reduced.  */
  for (yyi = 0; yyi < yynrhs; yyi++)
    {
      YYFPRINTF (stderr, "   $%d = ", yyi + 1);
      yy_symbol_print (stderr,
                       YY_ACCESSING_SYMBOL (+yyssp[yyi + 1 - yynrhs]),
                       &yyvsp[(yyi + 1) - (yynrhs)], scanner, lParse);
      YYFPRINTF (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)          \
do {                                    \
  if (yydebug)                          \
    yy_reduce_print (yyssp, yyvsp, Rule, scanner, lParse); \
} while (0)

/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int yydebug;
#else /* !FITS_PARSER_YYDEBUG */
# define YYDPRINTF(Args) ((void) 0)
# define YY_SYMBOL_PRINT(Title, Kind, Value, Location)
# define YY_STACK_PRINT(Bottom, Top)
# define YY_REDUCE_PRINT(Rule)
#endif /* !FITS_PARSER_YYDEBUG */


/* YYINITDEPTH -- initial size of the parser's stacks.  */
#ifndef YYINITDEPTH
# define YYINITDEPTH 200
#endif

/* YYMAXDEPTH -- maximum size the stacks can grow to (effective only
   if the built-in stack extension method is used).

   Do not make this value too large; the results are undefined if
   YYSTACK_ALLOC_MAXIMUM < YYSTACK_BYTES (YYMAXDEPTH)
   evaluated with infinite-precision integer arithmetic.  */

#ifndef YYMAXDEPTH
# define YYMAXDEPTH 10000
#endif






/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

static void
yydestruct (const char *yymsg,
            yysymbol_kind_t yykind, YYSTYPE *yyvaluep, yyscan_t scanner, ParseData *lParse)
{
  YY_USE (yyvaluep);
  YY_USE (scanner);
  YY_USE (lParse);
  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yykind, yyvaluep, yylocationp);

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  YY_USE (yykind);
  YY_IGNORE_MAYBE_UNINITIALIZED_END
}






/*----------.
| yyparse.  |
`----------*/

int
yyparse (yyscan_t scanner, ParseData *lParse)
{
/* Lookahead token kind.  */
int yychar;


/* The semantic value of the lookahead symbol.  */
/* Default value used for initialization, for pacifying older GCCs
   or non-GCC compilers.  */
YY_INITIAL_VALUE (static YYSTYPE yyval_default;)
YYSTYPE yylval YY_INITIAL_VALUE (= yyval_default);

    /* Number of syntax errors so far.  */
    int yynerrs = 0;

    yy_state_fast_t yystate = 0;
    /* Number of tokens to shift before error messages enabled.  */
    int yyerrstatus = 0;

    /* Refer to the stacks through separate pointers, to allow yyoverflow
       to reallocate them elsewhere.  */

    /* Their size.  */
    YYPTRDIFF_T yystacksize = YYINITDEPTH;

    /* The state stack: array, bottom, top.  */
    yy_state_t yyssa[YYINITDEPTH];
    yy_state_t *yyss = yyssa;
    yy_state_t *yyssp = yyss;

    /* The semantic value stack: array, bottom, top.  */
    YYSTYPE yyvsa[YYINITDEPTH];
    YYSTYPE *yyvs = yyvsa;
    YYSTYPE *yyvsp = yyvs;

  int yyn;
  /* The return value of yyparse.  */
  int yyresult;
  /* Lookahead symbol kind.  */
  yysymbol_kind_t yytoken = YYSYMBOL_YYEMPTY;
  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;



#define YYPOPSTACK(N)   (yyvsp -= (N), yyssp -= (N))

  /* The number of symbols on the RHS of the reduced rule.
     Keep to zero when no symbol should be popped.  */
  int yylen = 0;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yychar = FITS_PARSER_YYEMPTY; /* Cause a token to be read.  */

  goto yysetstate;


/*------------------------------------------------------------.
| yynewstate -- push a new state, which is found in yystate.  |
`------------------------------------------------------------*/
yynewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed.  So pushing a state here evens the stacks.  */
  yyssp++;


/*--------------------------------------------------------------------.
| yysetstate -- set current state (the top of the stack) to yystate.  |
`--------------------------------------------------------------------*/
yysetstate:
  YYDPRINTF ((stderr, "Entering state %d\n", yystate));
  YY_ASSERT (0 <= yystate && yystate < YYNSTATES);
  YY_IGNORE_USELESS_CAST_BEGIN
  *yyssp = YY_CAST (yy_state_t, yystate);
  YY_IGNORE_USELESS_CAST_END
  YY_STACK_PRINT (yyss, yyssp);

  if (yyss + yystacksize - 1 <= yyssp)
#if !defined yyoverflow && !defined YYSTACK_RELOCATE
    YYNOMEM;
#else
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYPTRDIFF_T yysize = yyssp - yyss + 1;

# if defined yyoverflow
      {
        /* Give user a chance to reallocate the stack.  Use copies of
           these so that the &'s don't force the real ones into
           memory.  */
        yy_state_t *yyss1 = yyss;
        YYSTYPE *yyvs1 = yyvs;

        /* Each stack pointer address is followed by the size of the
           data in use in that stack, in bytes.  This used to be a
           conditional around just the two extra args, but that might
           be undefined if yyoverflow is a macro.  */
        yyoverflow (YY_("memory exhausted"),
                    &yyss1, yysize * YYSIZEOF (*yyssp),
                    &yyvs1, yysize * YYSIZEOF (*yyvsp),
                    &yystacksize);
        yyss = yyss1;
        yyvs = yyvs1;
      }
# else /* defined YYSTACK_RELOCATE */
      /* Extend the stack our own way.  */
      if (YYMAXDEPTH <= yystacksize)
        YYNOMEM;
      yystacksize *= 2;
      if (YYMAXDEPTH < yystacksize)
        yystacksize = YYMAXDEPTH;

      {
        yy_state_t *yyss1 = yyss;
        union yyalloc *yyptr =
          YY_CAST (union yyalloc *,
                   YYSTACK_ALLOC (YY_CAST (YYSIZE_T, YYSTACK_BYTES (yystacksize))));
        if (! yyptr)
          YYNOMEM;
        YYSTACK_RELOCATE (yyss_alloc, yyss);
        YYSTACK_RELOCATE (yyvs_alloc, yyvs);
#  undef YYSTACK_RELOCATE
        if (yyss1 != yyssa)
          YYSTACK_FREE (yyss1);
      }
# endif

      yyssp = yyss + yysize - 1;
      yyvsp = yyvs + yysize - 1;

      YY_IGNORE_USELESS_CAST_BEGIN
      YYDPRINTF ((stderr, "Stack size increased to %ld\n",
                  YY_CAST (long, yystacksize)));
      YY_IGNORE_USELESS_CAST_END

      if (yyss + yystacksize - 1 <= yyssp)
        YYABORT;
    }
#endif /* !defined yyoverflow && !defined YYSTACK_RELOCATE */


  if (yystate == YYFINAL)
    YYACCEPT;

  goto yybackup;


/*-----------.
| yybackup.  |
`-----------*/
yybackup:
  /* Do appropriate processing given the current state.  Read a
     lookahead token if we need one and don't already have one.  */

  /* First try to decide what to do without reference to lookahead token.  */
  yyn = yypact[yystate];
  if (yypact_value_is_default (yyn))
    goto yydefault;

  /* Not known => get a lookahead token if don't already have one.  */

  /* YYCHAR is either empty, or end-of-input, or a valid lookahead.  */
  if (yychar == FITS_PARSER_YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token\n"));
      yychar = yylex (&yylval, scanner);
    }

  if (yychar <= FITS_PARSER_YYEOF)
    {
      yychar = FITS_PARSER_YYEOF;
      yytoken = YYSYMBOL_YYEOF;
      YYDPRINTF ((stderr, "Now at end of input.\n"));
    }
  else if (yychar == FITS_PARSER_YYerror)
    {
      /* The scanner already issued an error message, process directly
         to error recovery.  But do not keep the error token as
         lookahead, it is too special and may lead us to an endless
         loop in error recovery. */
      yychar = FITS_PARSER_YYUNDEF;
      yytoken = YYSYMBOL_YYerror;
      goto yyerrlab1;
    }
  else
    {
      yytoken = YYTRANSLATE (yychar);
      YY_SYMBOL_PRINT ("Next token is", yytoken, &yylval, &yylloc);
    }

  /* If the proper action on seeing token YYTOKEN is to reduce or to
     detect an error, take that action.  */
  yyn += yytoken;
  if (yyn < 0 || YYLAST < yyn || yycheck[yyn] != yytoken)
    goto yydefault;
  yyn = yytable[yyn];
  if (yyn <= 0)
    {
      if (yytable_value_is_error (yyn))
        goto yyerrlab;
      yyn = -yyn;
      goto yyreduce;
    }

  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  /* Shift the lookahead token.  */
  YY_SYMBOL_PRINT ("Shifting", yytoken, &yylval, &yylloc);
  yystate = yyn;
  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END

  /* Discard the shifted token.  */
  yychar = FITS_PARSER_YYEMPTY;
  goto yynewstate;


/*-----------------------------------------------------------.
| yydefault -- do the default action for the current state.  |
`-----------------------------------------------------------*/
yydefault:
  yyn = yydefact[yystate];
  if (yyn == 0)
    goto yyerrlab;
  goto yyreduce;


/*-----------------------------.
| yyreduce -- do a reduction.  |
`-----------------------------*/
yyreduce:
  /* yyn is the number of a rule to reduce with.  */
  yylen = yyr2[yyn];

  /* If YYLEN is nonzero, implement the default value of the action:
     '$$ = $1'.

     Otherwise, the following line sets YYVAL to garbage.
     This behavior is undocumented and Bison
     users should not rely upon it.  Assigning to YYVAL
     unconditionally makes the parser a bit smaller, and it avoids a
     GCC warning that YYVAL may be used uninitialized.  */
  yyval = yyvsp[1-yylen];


  YY_REDUCE_PRINT (yyn);
  switch (yyn)
    {
  case 4: /* line: '\n'  */
#line 270 "eval.y"
                     {}
#line 1821 "eval_y.c"
    break;

  case 5: /* line: expr '\n'  */
#line 272 "eval.y"
                { if( (yyvsp[-1].Node)<0 ) {
		     yyerror(scanner, lParse, "Couldn't build node structure: out of memory?");
		     YYERROR;  }
                  lParse->resultNode = (yyvsp[-1].Node);
		}
#line 1831 "eval_y.c"
    break;

  case 6: /* line: bexpr '\n'  */
#line 278 "eval.y"
                { if( (yyvsp[-1].Node)<0 ) {
		     yyerror(scanner, lParse, "Couldn't build node structure: out of memory?");
		     YYERROR;  }
                  lParse->resultNode = (yyvsp[-1].Node);
		}
#line 1841 "eval_y.c"
    break;

  case 7: /* line: sexpr '\n'  */
#line 284 "eval.y"
                { if( (yyvsp[-1].Node)<0 ) {
		     yyerror(scanner, lParse, "Couldn't build node structure: out of memory?");
		     YYERROR;  } 
                  lParse->resultNode = (yyvsp[-1].Node);
		}
#line 1851 "eval_y.c"
    break;

  case 8: /* line: bits '\n'  */
#line 290 "eval.y"
                { if( (yyvsp[-1].Node)<0 ) {
		     yyerror(scanner, lParse, "Couldn't build node structure: out of memory?");
		     YYERROR;  }
                  lParse->resultNode = (yyvsp[-1].Node);
		}
#line 1861 "eval_y.c"
    break;

  case 9: /* line: error '\n'  */
#line 295 "eval.y"
                     {  yyerrok;  }
#line 1867 "eval_y.c"
    break;

  case 10: /* bvector: '{' bexpr  */
#line 299 "eval.y"
                { (yyval.Node) = New_Vector(lParse,  (yyvsp[0].Node) ); TEST((yyval.Node)); }
#line 1873 "eval_y.c"
    break;

  case 11: /* bvector: bvector ',' bexpr  */
#line 301 "eval.y"
                {
                  if( lParse->Nodes[(yyvsp[-2].Node)].nSubNodes >= MAXSUBS ) {
		     (yyvsp[-2].Node) = Close_Vec(lParse,  (yyvsp[-2].Node) ); TEST((yyvsp[-2].Node));
		     (yyval.Node) = New_Vector(lParse,  (yyvsp[-2].Node) ); TEST((yyval.Node));
                  } else {
                     (yyval.Node) = (yyvsp[-2].Node);
                  }
		  lParse->Nodes[(yyval.Node)].SubNodes[ lParse->Nodes[(yyval.Node)].nSubNodes++ ]
		     = (yyvsp[0].Node);
                }
#line 1888 "eval_y.c"
    break;

  case 12: /* vector: '{' expr  */
#line 314 "eval.y"
                { (yyval.Node) = New_Vector(lParse,  (yyvsp[0].Node) ); TEST((yyval.Node)); }
#line 1894 "eval_y.c"
    break;

  case 13: /* vector: vector ',' expr  */
#line 316 "eval.y"
                {
                  if( TYPE((yyvsp[-2].Node)) < TYPE((yyvsp[0].Node)) )
                     TYPE((yyvsp[-2].Node)) = TYPE((yyvsp[0].Node));
                  if( lParse->Nodes[(yyvsp[-2].Node)].nSubNodes >= MAXSUBS ) {
		     (yyvsp[-2].Node) = Close_Vec(lParse,  (yyvsp[-2].Node) ); TEST((yyvsp[-2].Node));
		     (yyval.Node) = New_Vector(lParse,  (yyvsp[-2].Node) ); TEST((yyval.Node));
                  } else {
                     (yyval.Node) = (yyvsp[-2].Node);
                  }
		  lParse->Nodes[(yyval.Node)].SubNodes[ lParse->Nodes[(yyval.Node)].nSubNodes++ ]
		     = (yyvsp[0].Node);
                }
#line 1911 "eval_y.c"
    break;

  case 14: /* vector: vector ',' bexpr  */
#line 329 "eval.y"
                {
                  if( lParse->Nodes[(yyvsp[-2].Node)].nSubNodes >= MAXSUBS ) {
		     (yyvsp[-2].Node) = Close_Vec(lParse,  (yyvsp[-2].Node) ); TEST((yyvsp[-2].Node));
		     (yyval.Node) = New_Vector(lParse,  (yyvsp[-2].Node) ); TEST((yyval.Node));
                  } else {
                     (yyval.Node) = (yyvsp[-2].Node);
                  }
		  lParse->Nodes[(yyval.Node)].SubNodes[ lParse->Nodes[(yyval.Node)].nSubNodes++ ]
		     = (yyvsp[0].Node);
                }
#line 1926 "eval_y.c"
    break;

  case 15: /* vector: bvector ',' expr  */
#line 340 "eval.y"
                {
                  TYPE((yyvsp[-2].Node)) = TYPE((yyvsp[0].Node));
                  if( lParse->Nodes[(yyvsp[-2].Node)].nSubNodes >= MAXSUBS ) {
		     (yyvsp[-2].Node) = Close_Vec(lParse,  (yyvsp[-2].Node) ); TEST((yyvsp[-2].Node));
		     (yyval.Node) = New_Vector(lParse,  (yyvsp[-2].Node) ); TEST((yyval.Node));
                  } else {
                     (yyval.Node) = (yyvsp[-2].Node);
                  }
		  lParse->Nodes[(yyval.Node)].SubNodes[ lParse->Nodes[(yyval.Node)].nSubNodes++ ]
		     = (yyvsp[0].Node);
                }
#line 1942 "eval_y.c"
    break;

  case 16: /* expr: vector '}'  */
#line 354 "eval.y"
                { (yyval.Node) = Close_Vec(lParse,  (yyvsp[-1].Node) ); TEST((yyval.Node)); }
#line 1948 "eval_y.c"
    break;

  case 17: /* bexpr: bvector '}'  */
#line 358 "eval.y"
                { (yyval.Node) = Close_Vec(lParse,  (yyvsp[-1].Node) ); TEST((yyval.Node)); }
#line 1954 "eval_y.c"
    break;

  case 18: /* bits: BITSTR  */
#line 362 "eval.y"
                {
                  (yyval.Node) = New_Const(lParse,  BITSTR, (yyvsp[0].str), strlen((yyvsp[0].str))+1 ); TEST((yyval.Node));
		  SIZE((yyval.Node)) = strlen((yyvsp[0].str)); }
#line 1962 "eval_y.c"
    break;

  case 19: /* bits: BITCOL  */
#line 366 "eval.y"
                { (yyval.Node) = New_Column(lParse,  (yyvsp[0].lng) ); TEST((yyval.Node)); }
#line 1968 "eval_y.c"
    break;

  case 20: /* bits: BITCOL '{' expr '}'  */
#line 368 "eval.y"
                {
                  if( TYPE((yyvsp[-1].Node)) != LONG
		      || OPER((yyvsp[-1].Node)) != CONST_OP ) {
		     yyerror(scanner, lParse, "Offset argument must be a constant integer");
		     YYERROR;
		  }
                  (yyval.Node) = New_Offset(lParse,  (yyvsp[-3].lng), (yyvsp[-1].Node) ); TEST((yyval.Node));
                }
#line 1981 "eval_y.c"
    break;

  case 21: /* bits: bits '&' bits  */
#line 377 "eval.y"
                { (yyval.Node) = New_BinOp(lParse,  BITSTR, (yyvsp[-2].Node), '&', (yyvsp[0].Node) ); TEST((yyval.Node));
                  SIZE((yyval.Node)) = ( SIZE((yyvsp[-2].Node))>SIZE((yyvsp[0].Node)) ? SIZE((yyvsp[-2].Node)) : SIZE((yyvsp[0].Node)) );  }
#line 1988 "eval_y.c"
    break;

  case 22: /* bits: bits '|' bits  */
#line 380 "eval.y"
                { (yyval.Node) = New_BinOp(lParse,  BITSTR, (yyvsp[-2].Node), '|', (yyvsp[0].Node) ); TEST((yyval.Node));
                  SIZE((yyval.Node)) = ( SIZE((yyvsp[-2].Node))>SIZE((yyvsp[0].Node)) ? SIZE((yyvsp[-2].Node)) : SIZE((yyvsp[0].Node)) );  }
#line 1995 "eval_y.c"
    break;

  case 23: /* bits: bits '+' bits  */
#line 383 "eval.y"
                { 
		  if (SIZE((yyvsp[-2].Node))+SIZE((yyvsp[0].Node)) >= MAX_STRLEN) {
		    yyerror(scanner, lParse, "Combined bit string size exceeds " MAX_STRLEN_S " bits");
		    YYERROR;
		  }
		  (yyval.Node) = New_BinOp(lParse,  BITSTR, (yyvsp[-2].Node), '+', (yyvsp[0].Node) ); TEST((yyval.Node));
                  SIZE((yyval.Node)) = SIZE((yyvsp[-2].Node)) + SIZE((yyvsp[0].Node)); 
		}
#line 2008 "eval_y.c"
    break;

  case 24: /* bits: bits '[' expr ']'  */
#line 392 "eval.y"
                { (yyval.Node) = New_Deref(lParse,  (yyvsp[-3].Node), 1, (yyvsp[-1].Node),  0,  0,  0,   0 ); TEST((yyval.Node)); }
#line 2014 "eval_y.c"
    break;

  case 25: /* bits: bits '[' expr ',' expr ']'  */
#line 394 "eval.y"
                { (yyval.Node) = New_Deref(lParse,  (yyvsp[-5].Node), 2, (yyvsp[-3].Node), (yyvsp[-1].Node),  0,  0,   0 ); TEST((yyval.Node)); }
#line 2020 "eval_y.c"
    break;

  case 26: /* bits: bits '[' expr ',' expr ',' expr ']'  */
#line 396 "eval.y"
                { (yyval.Node) = New_Deref(lParse,  (yyvsp[-7].Node), 3, (yyvsp[-5].Node), (yyvsp[-3].Node), (yyvsp[-1].Node),  0,   0 ); TEST((yyval.Node)); }
#line 2026 "eval_y.c"
    break;

  case 27: /* bits: bits '[' expr ',' expr ',' expr ',' expr ']'  */
#line 398 "eval.y"
                { (yyval.Node) = New_Deref(lParse,  (yyvsp[-9].Node), 4, (yyvsp[-7].Node), (yyvsp[-5].Node), (yyvsp[-3].Node), (yyvsp[-1].Node),   0 ); TEST((yyval.Node)); }
#line 2032 "eval_y.c"
    break;

  case 28: /* bits: bits '[' expr ',' expr ',' expr ',' expr ',' expr ']'  */
#line 400 "eval.y"
                { (yyval.Node) = New_Deref(lParse,  (yyvsp[-11].Node), 5, (yyvsp[-9].Node), (yyvsp[-7].Node), (yyvsp[-5].Node), (yyvsp[-3].Node), (yyvsp[-1].Node) ); TEST((yyval.Node)); }
#line 2038 "eval_y.c"
    break;

  case 29: /* bits: NOT bits  */
#line 402 "eval.y"
                { (yyval.Node) = New_Unary(lParse,  BITSTR, NOT, (yyvsp[0].Node) ); TEST((yyval.Node));     }
#line 2044 "eval_y.c"
    break;

  case 30: /* bits: '(' bits ')'  */
#line 405 "eval.y"
                { (yyval.Node) = (yyvsp[-1].Node); }
#line 2050 "eval_y.c"
    break;

  case 31: /* expr: LONG  */
#line 409 "eval.y"
                { (yyval.Node) = New_Const(lParse,  LONG,   &((yyvsp[0].lng)), sizeof(long)   ); TEST((yyval.Node)); }
#line 2056 "eval_y.c"
    break;

  case 32: /* expr: DOUBLE  */
#line 411 "eval.y"
                { (yyval.Node) = New_Const(lParse,  DOUBLE, &((yyvsp[0].dbl)), sizeof(double) ); TEST((yyval.Node)); }
#line 2062 "eval_y.c"
    break;

  case 33: /* expr: COLUMN  */
#line 413 "eval.y"
                { (yyval.Node) = New_Column(lParse,  (yyvsp[0].lng) ); TEST((yyval.Node)); }
#line 2068 "eval_y.c"
    break;

  case 34: /* expr: COLUMN '{' expr '}'  */
#line 415 "eval.y"
                {
                  if( TYPE((yyvsp[-1].Node)) != LONG
		      || OPER((yyvsp[-1].Node)) != CONST_OP ) {
		     yyerror(scanner, lParse, "Offset argument must be a constant integer");
		     YYERROR;
		  }
                  (yyval.Node) = New_Offset(lParse,  (yyvsp[-3].lng), (yyvsp[-1].Node) ); TEST((yyval.Node));
                }
#line 2081 "eval_y.c"
    break;

  case 35: /* expr: ROWREF  */
#line 424 "eval.y"
                { (yyval.Node) = New_Func(lParse,  LONG, row_fct,  0, 0, 0, 0, 0, 0, 0, 0 ); }
#line 2087 "eval_y.c"
    break;

  case 36: /* expr: NULLREF  */
#line 426 "eval.y"
                { (yyval.Node) = New_Func(lParse,  LONG, null_fct, 0, 0, 0, 0, 0, 0, 0, 0 ); }
#line 2093 "eval_y.c"
    break;

  case 37: /* expr: expr '%' expr  */
#line 428 "eval.y"
                { PROMOTE((yyvsp[-2].Node),(yyvsp[0].Node)); (yyval.Node) = New_BinOp(lParse,  TYPE((yyvsp[-2].Node)), (yyvsp[-2].Node), '%', (yyvsp[0].Node) );
		  TEST((yyval.Node));                                                }
#line 2100 "eval_y.c"
    break;

  case 38: /* expr: expr '+' expr  */
#line 431 "eval.y"
                { PROMOTE((yyvsp[-2].Node),(yyvsp[0].Node)); (yyval.Node) = New_BinOp(lParse,  TYPE((yyvsp[-2].Node)), (yyvsp[-2].Node), '+', (yyvsp[0].Node) );
		  TEST((yyval.Node));                                                }
#line 2107 "eval_y.c"
    break;

  case 39: /* expr: expr '-' expr  */
#line 434 "eval.y"
                { PROMOTE((yyvsp[-2].Node),(yyvsp[0].Node)); (yyval.Node) = New_BinOp(lParse,  TYPE((yyvsp[-2].Node)), (yyvsp[-2].Node), '-', (yyvsp[0].Node) ); 
		  TEST((yyval.Node));                                                }
#line 2114 "eval_y.c"
    break;

  case 40: /* expr: expr '*' expr  */
#line 437 "eval.y"
                { PROMOTE((yyvsp[-2].Node),(yyvsp[0].Node)); (yyval.Node) = New_BinOp(lParse,  TYPE((yyvsp[-2].Node)), (yyvsp[-2].Node), '*', (yyvsp[0].Node) ); 
		  TEST((yyval.Node));                                                }
#line 2121 "eval_y.c"
    break;

  case 41: /* expr: expr '/' expr  */
#line 440 "eval.y"
                { PROMOTE((yyvsp[-2].Node),(yyvsp[0].Node)); (yyval.Node) = New_BinOp(lParse,  TYPE((yyvsp[-2].Node)), (yyvsp[-2].Node), '/', (yyvsp[0].Node) ); 
		  TEST((yyval.Node));                                                }
#line 2128 "eval_y.c"
    break;

  case 42: /* expr: expr '&' expr  */
#line 443 "eval.y"
                { 
                   if (TYPE((yyvsp[-2].Node)) != LONG ||
		       TYPE((yyvsp[0].Node)) != LONG) {
                     yyerror(scanner, lParse, "Bitwise operations with incompatible types; only (bit OP bit) and (int OP int) are allowed");
                      YYERROR;
                   }
                   (yyval.Node) = New_BinOp(lParse,  TYPE((yyvsp[-2].Node)), (yyvsp[-2].Node), '&', (yyvsp[0].Node) );
                }
#line 2141 "eval_y.c"
    break;

  case 43: /* expr: expr '|' expr  */
#line 452 "eval.y"
                { 
                   if (TYPE((yyvsp[-2].Node)) != LONG ||
		       TYPE((yyvsp[0].Node)) != LONG) {
                     yyerror(scanner, lParse, "Bitwise operations with incompatible types; only (bit OP bit) and (int OP int) are allowed");
                      YYERROR;
                   }
                   (yyval.Node) = New_BinOp(lParse,  TYPE((yyvsp[-2].Node)), (yyvsp[-2].Node), '|', (yyvsp[0].Node) );
                }
#line 2154 "eval_y.c"
    break;

  case 44: /* expr: expr XOR expr  */
#line 461 "eval.y"
                { 
                   if (TYPE((yyvsp[-2].Node)) != LONG ||
		       TYPE((yyvsp[0].Node)) != LONG) {
                     yyerror(scanner, lParse, "Bitwise operations with incompatible types; only (bit OP bit) and (int OP int) are allowed");
                      YYERROR;
                   }
                   (yyval.Node) = New_BinOp(lParse,  TYPE((yyvsp[-2].Node)), (yyvsp[-2].Node), '^', (yyvsp[0].Node) );
                }
#line 2167 "eval_y.c"
    break;

  case 45: /* expr: expr POWER expr  */
#line 470 "eval.y"
                { PROMOTE((yyvsp[-2].Node),(yyvsp[0].Node)); (yyval.Node) = New_BinOp(lParse,  TYPE((yyvsp[-2].Node)), (yyvsp[-2].Node), POWER, (yyvsp[0].Node) );
		  TEST((yyval.Node));                                                }
#line 2174 "eval_y.c"
    break;

  case 46: /* expr: '+' expr  */
#line 473 "eval.y"
                { (yyval.Node) = (yyvsp[0].Node); }
#line 2180 "eval_y.c"
    break;

  case 47: /* expr: '-' expr  */
#line 475 "eval.y"
                { (yyval.Node) = New_Unary(lParse,  TYPE((yyvsp[0].Node)), UMINUS, (yyvsp[0].Node) ); TEST((yyval.Node)); }
#line 2186 "eval_y.c"
    break;

  case 48: /* expr: '(' expr ')'  */
#line 477 "eval.y"
                { (yyval.Node) = (yyvsp[-1].Node); }
#line 2192 "eval_y.c"
    break;

  case 49: /* expr: expr '*' bexpr  */
#line 479 "eval.y"
                { (yyvsp[0].Node) = New_Unary(lParse,  TYPE((yyvsp[-2].Node)), 0, (yyvsp[0].Node) );
                  (yyval.Node) = New_BinOp(lParse,  TYPE((yyvsp[-2].Node)), (yyvsp[-2].Node), '*', (yyvsp[0].Node) ); 
		  TEST((yyval.Node));                                }
#line 2200 "eval_y.c"
    break;

  case 50: /* expr: bexpr '*' expr  */
#line 483 "eval.y"
                { (yyvsp[-2].Node) = New_Unary(lParse,  TYPE((yyvsp[0].Node)), 0, (yyvsp[-2].Node) );
                  (yyval.Node) = New_BinOp(lParse,  TYPE((yyvsp[0].Node)), (yyvsp[-2].Node), '*', (yyvsp[0].Node) );
                  TEST((yyval.Node));                                }
#line 2208 "eval_y.c"
    break;

  case 51: /* expr: bexpr '?' expr ':' expr  */
#line 487 "eval.y"
                {
                  PROMOTE((yyvsp[-2].Node),(yyvsp[0].Node));
                  if( ! Test_Dims( lParse, (yyvsp[-2].Node),(yyvsp[0].Node)) ) {
                     yyerror(scanner, lParse, "Incompatible dimensions in '?:' arguments");
		     YYERROR;
                  }
                  (yyval.Node) = New_Func(lParse,  0, ifthenelse_fct, 3, (yyvsp[-2].Node), (yyvsp[0].Node), (yyvsp[-4].Node),
                                 0, 0, 0, 0 );
                  TEST((yyval.Node));
                  if( SIZE((yyvsp[-2].Node))<SIZE((yyvsp[0].Node)) )  Copy_Dims( lParse,(yyval.Node), (yyvsp[0].Node));
                  TYPE((yyvsp[-4].Node)) = TYPE((yyvsp[-2].Node));
                  if( ! Test_Dims( lParse, (yyvsp[-4].Node),(yyval.Node)) ) {
                     yyerror(scanner, lParse, "Incompatible dimensions in '?:' condition");
		     YYERROR;
                  }
                  TYPE((yyvsp[-4].Node)) = BOOLEAN;
                  if( SIZE((yyval.Node))<SIZE((yyvsp[-4].Node)) )  Copy_Dims( lParse,(yyval.Node), (yyvsp[-4].Node));
                }
#line 2231 "eval_y.c"
    break;

  case 52: /* expr: bexpr '?' bexpr ':' expr  */
#line 506 "eval.y"
                {
                  PROMOTE((yyvsp[-2].Node),(yyvsp[0].Node));
                  if( ! Test_Dims( lParse, (yyvsp[-2].Node),(yyvsp[0].Node)) ) {
                     yyerror(scanner, lParse, "Incompatible dimensions in '?:' arguments");
		     YYERROR;
                  }
                  (yyval.Node) = New_Func(lParse,  0, ifthenelse_fct, 3, (yyvsp[-2].Node), (yyvsp[0].Node), (yyvsp[-4].Node),
                                 0, 0, 0, 0 );
                  TEST((yyval.Node));
                  if( SIZE((yyvsp[-2].Node))<SIZE((yyvsp[0].Node)) )  Copy_Dims( lParse,(yyval.Node), (yyvsp[0].Node));
                  TYPE((yyvsp[-4].Node)) = TYPE((yyvsp[-2].Node));
                  if( ! Test_Dims( lParse, (yyvsp[-4].Node),(yyval.Node)) ) {
                     yyerror(scanner, lParse, "Incompatible dimensions in '?:' condition");
		     YYERROR;
                  }
                  TYPE((yyvsp[-4].Node)) = BOOLEAN;
                  if( SIZE((yyval.Node))<SIZE((yyvsp[-4].Node)) )  Copy_Dims( lParse,(yyval.Node), (yyvsp[-4].Node));
                }
#line 2254 "eval_y.c"
    break;

  case 53: /* expr: bexpr '?' expr ':' bexpr  */
#line 525 "eval.y"
                {
                  PROMOTE((yyvsp[-2].Node),(yyvsp[0].Node));
                  if( ! Test_Dims( lParse, (yyvsp[-2].Node),(yyvsp[0].Node)) ) {
                     yyerror(scanner, lParse, "Incompatible dimensions in '?:' arguments");
		     YYERROR;
                  }
                  (yyval.Node) = New_Func(lParse,  0, ifthenelse_fct, 3, (yyvsp[-2].Node), (yyvsp[0].Node), (yyvsp[-4].Node),
                                 0, 0, 0, 0 );
                  TEST((yyval.Node));
                  if( SIZE((yyvsp[-2].Node))<SIZE((yyvsp[0].Node)) )  Copy_Dims( lParse,(yyval.Node), (yyvsp[0].Node));
                  TYPE((yyvsp[-4].Node)) = TYPE((yyvsp[-2].Node));
                  if( ! Test_Dims( lParse, (yyvsp[-4].Node),(yyval.Node)) ) {
                     yyerror(scanner, lParse, "Incompatible dimensions in '?:' condition");
		     YYERROR;
                  }
                  TYPE((yyvsp[-4].Node)) = BOOLEAN;
                  if( SIZE((yyval.Node))<SIZE((yyvsp[-4].Node)) )  Copy_Dims( lParse,(yyval.Node), (yyvsp[-4].Node));
                }
#line 2277 "eval_y.c"
    break;

  case 54: /* expr: FUNCTION ')'  */
#line 544 "eval.y"
                { if (FSTRCMP((yyvsp[-1].str),"RANDOM(") == 0) {  /* Scalar RANDOM() */
                     (yyval.Node) = New_Func(lParse,  DOUBLE, rnd_fct, 0, 0, 0, 0, 0, 0, 0, 0 );
		  } else if (FSTRCMP((yyvsp[-1].str),"RANDOMN(") == 0) {/*Scalar RANDOMN()*/
		     (yyval.Node) = New_Func(lParse,  DOUBLE, gasrnd_fct, 0, 0, 0, 0, 0, 0, 0, 0 );
                  } else {
                     yyerror(scanner, lParse, "Function() not supported");
		     YYERROR;
		  }
                  TEST((yyval.Node)); 
                }
#line 2292 "eval_y.c"
    break;

  case 55: /* expr: FUNCTION bexpr ')'  */
#line 555 "eval.y"
                { if (FSTRCMP((yyvsp[-2].str),"SUM(") == 0) {
		     (yyval.Node) = New_Func(lParse,  LONG, sum_fct, 1, (yyvsp[-1].Node), 0, 0, 0, 0, 0, 0 );
                  } else if (FSTRCMP((yyvsp[-2].str),"NELEM(") == 0) {
                     (yyval.Node) = New_Const(lParse,  LONG, &( SIZE((yyvsp[-1].Node)) ), sizeof(long) );
                  } else if (FSTRCMP((yyvsp[-2].str),"ACCUM(") == 0) {
		    long zero = 0;
		    (yyval.Node) = New_BinOp(lParse,  LONG , (yyvsp[-1].Node), ACCUM, New_Const(lParse,  LONG, &zero, sizeof(zero) ));
		  } else {
                     yyerror(scanner, lParse, "Function(bool) not supported");
		     YYERROR;
		  }
                  TEST((yyval.Node)); 
		}
#line 2310 "eval_y.c"
    break;

  case 56: /* expr: FUNCTION bexpr ',' expr ')'  */
#line 569 "eval.y"
                { if (FSTRCMP((yyvsp[-4].str),"AXISELEM(") == 0) {  /* AXISELEM(V,n) */
		     if (OPER((yyvsp[-1].Node)) != CONST_OP
			 || SIZE((yyvsp[-1].Node)) != 1) {
		       yyerror(scanner, lParse, "AXISELEM second argument must be a scalar constant");
		       YYERROR;
		     }
		     if (OPER((yyvsp[-3].Node)) == CONST_OP) {
		       long one = 1;
		       (yyval.Node) = New_Const(lParse,  LONG, &one, sizeof(one) );
		     } else {
		       if ( TYPE((yyvsp[-1].Node)) != LONG ) (yyvsp[-1].Node) = New_Unary(lParse, LONG, 0, (yyvsp[-1].Node));
		       (yyval.Node) = New_Func(lParse, 0, axiselem_fct, 2, (yyvsp[-3].Node), (yyvsp[-1].Node), 0, 0, 0, 0, 0 );
		       TEST((yyval.Node));
		       TYPE((yyval.Node)) = LONG;
		     }
		   } else if (FSTRCMP((yyvsp[-4].str),"NAXES(") == 0) {  /* NAXES(V,n) */
		     if (OPER((yyvsp[-1].Node)) != CONST_OP
			 || SIZE((yyvsp[-1].Node)) != 1) {
		       yyerror(scanner, lParse, "NAXES second argument must be a scalar constant");
		       YYERROR;
		     }
		     if (OPER((yyvsp[-3].Node)) == CONST_OP) { /* if V is constant, return 1 in every case */
		       long one = 1;
		       (yyval.Node) = New_Const(lParse,  LONG, &one, sizeof(one) );
		     } else {                    /* determine now the dimension of the expression */
		       long iaxis;
		       int naxis;
		       if ( TYPE((yyvsp[-1].Node)) != LONG ) (yyvsp[-1].Node) = New_Unary(lParse, LONG, 0, (yyvsp[-1].Node));
		       /* Since it is already constant, we can extract long value directly */
		       iaxis = (lParse->Nodes[(yyvsp[-1].Node)].value.data.lng);
		       naxis = lParse->Nodes[(yyvsp[-3].Node)].value.naxis;

		       if (iaxis == 0)          iaxis = naxis;   /* NAXIS(V,0) = NAXIS */
		       else if (iaxis <= naxis) iaxis = lParse->Nodes[(yyvsp[-3].Node)].value.naxes[iaxis-1]; /* NAXIS(V,n) = NAXISn */
		       else                     iaxis = 1;       /* Out of bounds use 1 */

		       (yyval.Node) = New_Const(lParse,  LONG, &iaxis, sizeof(iaxis) );
		       TEST((yyval.Node));
		     }
		   } else if (FSTRCMP((yyvsp[-4].str),"ARRAY(") == 0) {  /* NAXES(bexpr,n) */
		     (yyval.Node) = New_Array(lParse, (yyvsp[-3].Node), (yyvsp[-1].Node));
		     TEST((yyval.Node));
		  } else {
                     yyerror(scanner, lParse, "Function(bool,expr) not supported");
		     YYERROR;
		  }
                  TEST((yyval.Node)); 
		}
#line 2363 "eval_y.c"
    break;

  case 57: /* expr: FUNCTION sexpr ')'  */
#line 618 "eval.y"
                { if (FSTRCMP((yyvsp[-2].str),"NELEM(") == 0) {
                     (yyval.Node) = New_Const(lParse,  LONG, &( SIZE((yyvsp[-1].Node)) ), sizeof(long) );
		  } else if (FSTRCMP((yyvsp[-2].str),"NVALID(") == 0) {
		     (yyval.Node) = New_Func(lParse,  LONG, nonnull_fct, 1, (yyvsp[-1].Node),
				    0, 0, 0, 0, 0, 0 );
		  } else {
                     yyerror(scanner, lParse, "Function(str) not supported");
		     YYERROR;
		  }
                  TEST((yyval.Node)); 
		}
#line 2379 "eval_y.c"
    break;

  case 58: /* expr: FUNCTION bits ')'  */
#line 630 "eval.y"
                { if (FSTRCMP((yyvsp[-2].str),"NELEM(") == 0) {
                     (yyval.Node) = New_Const(lParse,  LONG, &( SIZE((yyvsp[-1].Node)) ), sizeof(long) );
		} else if (FSTRCMP((yyvsp[-2].str),"NVALID(") == 0) { /* Bit arrays do not have NULL */
                     (yyval.Node) = New_Const(lParse,  LONG, &( SIZE((yyvsp[-1].Node)) ), sizeof(long) );
		} else if (FSTRCMP((yyvsp[-2].str),"SUM(") == 0) {
		     (yyval.Node) = New_Func(lParse,  LONG, sum_fct, 1, (yyvsp[-1].Node),
				    0, 0, 0, 0, 0, 0 );
		} else if (FSTRCMP((yyvsp[-2].str),"MIN(") == 0) {
		     (yyval.Node) = New_Func(lParse,  TYPE((yyvsp[-1].Node)),  /* Force 1D result */
				    min1_fct, 1, (yyvsp[-1].Node), 0, 0, 0, 0, 0, 0 );
		     /* Note: $2 is a vector so the result can never
		        be a constant.  Therefore it will never be set
		        inside New_Func(), and it is safe to set SIZE() */
		     SIZE((yyval.Node)) = 1;
		} else if (FSTRCMP((yyvsp[-2].str),"ACCUM(") == 0) {
		    long zero = 0;
		    (yyval.Node) = New_BinOp(lParse,  LONG , (yyvsp[-1].Node), ACCUM, New_Const(lParse,  LONG, &zero, sizeof(zero) ));
		} else if (FSTRCMP((yyvsp[-2].str),"MAX(") == 0) {
		     (yyval.Node) = New_Func(lParse,  TYPE((yyvsp[-1].Node)),  /* Force 1D result */
				    max1_fct, 1, (yyvsp[-1].Node), 0, 0, 0, 0, 0, 0 );
		     /* Note: $2 is a vector so the result can never
		        be a constant.  Therefore it will never be set
		        inside New_Func(), and it is safe to set SIZE() */
		     SIZE((yyval.Node)) = 1;
		} else {
                     yyerror(scanner, lParse, "Function(bits) not supported");
		     YYERROR;
		  }
                  TEST((yyval.Node)); 
		}
#line 2414 "eval_y.c"
    break;

  case 59: /* expr: FUNCTION expr ')'  */
#line 661 "eval.y"
                { if (FSTRCMP((yyvsp[-2].str),"SUM(") == 0)
		     (yyval.Node) = New_Func(lParse,  TYPE((yyvsp[-1].Node)), sum_fct, 1, (yyvsp[-1].Node),
				    0, 0, 0, 0, 0, 0 );
		  else if (FSTRCMP((yyvsp[-2].str),"AVERAGE(") == 0)
		     (yyval.Node) = New_Func(lParse,  DOUBLE, average_fct, 1, (yyvsp[-1].Node),
				    0, 0, 0, 0, 0, 0 );
		  else if (FSTRCMP((yyvsp[-2].str),"STDDEV(") == 0)
		     (yyval.Node) = New_Func(lParse,  DOUBLE, stddev_fct, 1, (yyvsp[-1].Node),
				    0, 0, 0, 0, 0, 0 );
		  else if (FSTRCMP((yyvsp[-2].str),"MEDIAN(") == 0)
		     (yyval.Node) = New_Func(lParse,  TYPE((yyvsp[-1].Node)), median_fct, 1, (yyvsp[-1].Node),
				    0, 0, 0, 0, 0, 0 );
		  else if (FSTRCMP((yyvsp[-2].str),"NELEM(") == 0)
                     (yyval.Node) = New_Const(lParse,  LONG, &( SIZE((yyvsp[-1].Node)) ), sizeof(long) );
		  else if (FSTRCMP((yyvsp[-2].str),"NVALID(") == 0)
		     (yyval.Node) = New_Func(lParse,  LONG, nonnull_fct, 1, (yyvsp[-1].Node),
				    0, 0, 0, 0, 0, 0 );
		  else if   ((FSTRCMP((yyvsp[-2].str),"ACCUM(") == 0) && (TYPE((yyvsp[-1].Node)) == LONG)) {
		    long zero = 0;
		    (yyval.Node) = New_BinOp(lParse,  LONG ,   (yyvsp[-1].Node), ACCUM, New_Const(lParse,  LONG,   &zero, sizeof(zero) ));
		  } else if ((FSTRCMP((yyvsp[-2].str),"ACCUM(") == 0) && (TYPE((yyvsp[-1].Node)) == DOUBLE)) {
		    double zero = 0;
		    (yyval.Node) = New_BinOp(lParse,  DOUBLE , (yyvsp[-1].Node), ACCUM, New_Const(lParse,  DOUBLE, &zero, sizeof(zero) ));
		  } else if ((FSTRCMP((yyvsp[-2].str),"SEQDIFF(") == 0) && (TYPE((yyvsp[-1].Node)) == LONG)) {
		    long zero = 0;
		    (yyval.Node) = New_BinOp(lParse,  LONG ,   (yyvsp[-1].Node), DIFF, New_Const(lParse,  LONG,   &zero, sizeof(zero) ));
		  } else if ((FSTRCMP((yyvsp[-2].str),"SEQDIFF(") == 0) && (TYPE((yyvsp[-1].Node)) == DOUBLE)) {
		    double zero = 0;
		    (yyval.Node) = New_BinOp(lParse,  DOUBLE , (yyvsp[-1].Node), DIFF, New_Const(lParse,  DOUBLE, &zero, sizeof(zero) ));
		  } else if (FSTRCMP((yyvsp[-2].str),"ABS(") == 0)
		     (yyval.Node) = New_Func(lParse,  0, abs_fct, 1, (yyvsp[-1].Node), 0, 0, 0, 0, 0, 0 );
 		  else if (FSTRCMP((yyvsp[-2].str),"MIN(") == 0)
		     (yyval.Node) = New_Func(lParse,  TYPE((yyvsp[-1].Node)),  /* Force 1D result */
				    min1_fct, 1, (yyvsp[-1].Node), 0, 0, 0, 0, 0, 0 );
		  else if (FSTRCMP((yyvsp[-2].str),"MAX(") == 0)
		     (yyval.Node) = New_Func(lParse,  TYPE((yyvsp[-1].Node)),  /* Force 1D result */
				    max1_fct, 1, (yyvsp[-1].Node), 0, 0, 0, 0, 0, 0 );
		  else if (FSTRCMP((yyvsp[-2].str),"RANDOM(") == 0) { /* Vector RANDOM() */
                     (yyval.Node) = New_Func(lParse,  0, rnd_fct, 1, (yyvsp[-1].Node), 0, 0, 0, 0, 0, 0 );
		     TEST((yyval.Node));
		     TYPE((yyval.Node)) = DOUBLE;
		  } else if (FSTRCMP((yyvsp[-2].str),"RANDOMN(") == 0) {
		     (yyval.Node) = New_Func(lParse,  0, gasrnd_fct, 1, (yyvsp[-1].Node), 0, 0, 0, 0, 0, 0 );
		     TEST((yyval.Node));
		     TYPE((yyval.Node)) = DOUBLE;
		  } else if (FSTRCMP((yyvsp[-2].str),"ELEMENTNUM(") == 0) {
		     if (OPER((yyvsp[-1].Node)) == CONST_OP) {
		       long one = 1;
		       (yyval.Node) = New_Const(lParse,  LONG, &one, sizeof(one) );
		     } else {
		       (yyval.Node) = New_Func(lParse,  0, elemnum_fct, 1, (yyvsp[-1].Node), 0, 0, 0, 0, 0, 0 );
		       TEST((yyval.Node));
		       TYPE((yyval.Node)) = LONG;
		     }
		  } else if (FSTRCMP((yyvsp[-2].str),"NAXIS(") == 0) {  /* NAXIS(V) */
		     if (OPER((yyvsp[-1].Node)) == CONST_OP) { /* if V is constant, return 1 in every case */
		       long one = 1;
		       (yyval.Node) = New_Const(lParse,  LONG, &one, sizeof(one) );
		     } else {                    /* determine now the dimension of the expression */
		       long naxis = lParse->Nodes[(yyvsp[-1].Node)].value.naxis;

		       (yyval.Node) = New_Const(lParse,  LONG, &naxis, sizeof(naxis) );
		       TEST((yyval.Node));
		     }
                  } 
  		  else {  /*  These all take DOUBLE arguments  */
		     if( TYPE((yyvsp[-1].Node)) != DOUBLE ) (yyvsp[-1].Node) = New_Unary(lParse,  DOUBLE, 0, (yyvsp[-1].Node) );
                     if (FSTRCMP((yyvsp[-2].str),"SIN(") == 0)
			(yyval.Node) = New_Func(lParse,  0, sin_fct,  1, (yyvsp[-1].Node), 0, 0, 0, 0, 0, 0 );
		     else if (FSTRCMP((yyvsp[-2].str),"COS(") == 0)
			(yyval.Node) = New_Func(lParse,  0, cos_fct,  1, (yyvsp[-1].Node), 0, 0, 0, 0, 0, 0 );
		     else if (FSTRCMP((yyvsp[-2].str),"TAN(") == 0)
			(yyval.Node) = New_Func(lParse,  0, tan_fct,  1, (yyvsp[-1].Node), 0, 0, 0, 0, 0, 0 );
		     else if (FSTRCMP((yyvsp[-2].str),"ARCSIN(") == 0
			      || FSTRCMP((yyvsp[-2].str),"ASIN(") == 0)
			(yyval.Node) = New_Func(lParse,  0, asin_fct, 1, (yyvsp[-1].Node), 0, 0, 0, 0, 0, 0 );
		     else if (FSTRCMP((yyvsp[-2].str),"ARCCOS(") == 0
			      || FSTRCMP((yyvsp[-2].str),"ACOS(") == 0)
			(yyval.Node) = New_Func(lParse,  0, acos_fct, 1, (yyvsp[-1].Node), 0, 0, 0, 0, 0, 0 );
		     else if (FSTRCMP((yyvsp[-2].str),"ARCTAN(") == 0
			      || FSTRCMP((yyvsp[-2].str),"ATAN(") == 0)
			(yyval.Node) = New_Func(lParse,  0, atan_fct, 1, (yyvsp[-1].Node), 0, 0, 0, 0, 0, 0 );
		     else if (FSTRCMP((yyvsp[-2].str),"SINH(") == 0)
			(yyval.Node) = New_Func(lParse,  0, sinh_fct,  1, (yyvsp[-1].Node), 0, 0, 0, 0, 0, 0 );
		     else if (FSTRCMP((yyvsp[-2].str),"COSH(") == 0)
			(yyval.Node) = New_Func(lParse,  0, cosh_fct,  1, (yyvsp[-1].Node), 0, 0, 0, 0, 0, 0 );
		     else if (FSTRCMP((yyvsp[-2].str),"TANH(") == 0)
			(yyval.Node) = New_Func(lParse,  0, tanh_fct,  1, (yyvsp[-1].Node), 0, 0, 0, 0, 0, 0 );
		     else if (FSTRCMP((yyvsp[-2].str),"EXP(") == 0)
			(yyval.Node) = New_Func(lParse,  0, exp_fct,  1, (yyvsp[-1].Node), 0, 0, 0, 0, 0, 0 );
		     else if (FSTRCMP((yyvsp[-2].str),"LOG(") == 0)
			(yyval.Node) = New_Func(lParse,  0, log_fct,  1, (yyvsp[-1].Node), 0, 0, 0, 0, 0, 0 );
		     else if (FSTRCMP((yyvsp[-2].str),"LOG10(") == 0)
			(yyval.Node) = New_Func(lParse,  0, log10_fct, 1, (yyvsp[-1].Node), 0, 0, 0, 0, 0, 0 );
		     else if (FSTRCMP((yyvsp[-2].str),"SQRT(") == 0)
			(yyval.Node) = New_Func(lParse,  0, sqrt_fct, 1, (yyvsp[-1].Node), 0, 0, 0, 0, 0, 0 );
		     else if (FSTRCMP((yyvsp[-2].str),"ROUND(") == 0)
			(yyval.Node) = New_Func(lParse,  0, round_fct, 1, (yyvsp[-1].Node), 0, 0, 0, 0, 0, 0 );
		     else if (FSTRCMP((yyvsp[-2].str),"FLOOR(") == 0)
			(yyval.Node) = New_Func(lParse,  0, floor_fct, 1, (yyvsp[-1].Node), 0, 0, 0, 0, 0, 0 );
		     else if (FSTRCMP((yyvsp[-2].str),"CEIL(") == 0)
			(yyval.Node) = New_Func(lParse,  0, ceil_fct, 1, (yyvsp[-1].Node), 0, 0, 0, 0, 0, 0 );
		     else if (FSTRCMP((yyvsp[-2].str),"RANDOMP(") == 0) {
		       (yyval.Node) = New_Func(lParse,  0, poirnd_fct, 1, (yyvsp[-1].Node), 
				      0, 0, 0, 0, 0, 0 );
		       TYPE((yyval.Node)) = LONG;
		     } else {
			yyerror(scanner, lParse, "Function(expr) not supported");
			YYERROR;
		     }
		  }
                  TEST((yyval.Node)); 
                }
#line 2532 "eval_y.c"
    break;

  case 60: /* expr: IFUNCTION sexpr ',' sexpr ')'  */
#line 775 "eval.y"
                { 
		  if (FSTRCMP((yyvsp[-4].str),"STRSTR(") == 0) {
		    (yyval.Node) = New_Func(lParse,  LONG, strpos_fct, 2, (yyvsp[-3].Node), (yyvsp[-1].Node), 0, 
				   0, 0, 0, 0 );
		    TEST((yyval.Node));
		  }
                }
#line 2544 "eval_y.c"
    break;

  case 61: /* expr: FUNCTION expr ',' expr ')'  */
#line 783 "eval.y"
                { 
		   if (FSTRCMP((yyvsp[-4].str),"DEFNULL(") == 0) {
		      if( SIZE((yyvsp[-3].Node))>=SIZE((yyvsp[-1].Node)) && Test_Dims( lParse,  (yyvsp[-3].Node), (yyvsp[-1].Node) ) ) {
			 PROMOTE((yyvsp[-3].Node),(yyvsp[-1].Node));
			 (yyval.Node) = New_Func(lParse,  0, defnull_fct, 2, (yyvsp[-3].Node), (yyvsp[-1].Node), 0,
					0, 0, 0, 0 );
			 TEST((yyval.Node)); 
		      } else {
			 yyerror(scanner, lParse, "Dimensions of DEFNULL arguments "
				 "are not compatible");
			 YYERROR;
		      }
		   } else if (FSTRCMP((yyvsp[-4].str),"ARCTAN2(") == 0) {
		     if( TYPE((yyvsp[-3].Node)) != DOUBLE ) (yyvsp[-3].Node) = New_Unary(lParse,  DOUBLE, 0, (yyvsp[-3].Node) );
		     if( TYPE((yyvsp[-1].Node)) != DOUBLE ) (yyvsp[-1].Node) = New_Unary(lParse,  DOUBLE, 0, (yyvsp[-1].Node) );
		     if( Test_Dims( lParse,  (yyvsp[-3].Node), (yyvsp[-1].Node) ) ) {
			(yyval.Node) = New_Func(lParse,  0, atan2_fct, 2, (yyvsp[-3].Node), (yyvsp[-1].Node), 0, 0, 0, 0, 0 );
			TEST((yyval.Node)); 
			if( SIZE((yyvsp[-3].Node))<SIZE((yyvsp[-1].Node)) ) Copy_Dims( lParse,(yyval.Node), (yyvsp[-1].Node));
		     } else {
			yyerror(scanner, lParse, "Dimensions of arctan2 arguments "
				"are not compatible");
			YYERROR;
		     }
		   } else if (FSTRCMP((yyvsp[-4].str),"MIN(") == 0) {
		      PROMOTE( (yyvsp[-3].Node), (yyvsp[-1].Node) );
		      if( Test_Dims( lParse,  (yyvsp[-3].Node), (yyvsp[-1].Node) ) ) {
			(yyval.Node) = New_Func(lParse,  0, min2_fct, 2, (yyvsp[-3].Node), (yyvsp[-1].Node), 0, 0, 0, 0, 0 );
			TEST((yyval.Node));
			if( SIZE((yyvsp[-3].Node))<SIZE((yyvsp[-1].Node)) ) Copy_Dims( lParse,(yyval.Node), (yyvsp[-1].Node));
		      } else {
			yyerror(scanner, lParse, "Dimensions of min(a,b) arguments "
				"are not compatible");
			YYERROR;
		      }
		   } else if (FSTRCMP((yyvsp[-4].str),"MAX(") == 0) {
		      PROMOTE( (yyvsp[-3].Node), (yyvsp[-1].Node) );
		      if( Test_Dims( lParse,  (yyvsp[-3].Node), (yyvsp[-1].Node) ) ) {
			(yyval.Node) = New_Func(lParse,  0, max2_fct, 2, (yyvsp[-3].Node), (yyvsp[-1].Node), 0, 0, 0, 0, 0 );
			TEST((yyval.Node));
			if( SIZE((yyvsp[-3].Node))<SIZE((yyvsp[-1].Node)) ) Copy_Dims( lParse,(yyval.Node), (yyvsp[-1].Node));
		      } else {
			yyerror(scanner, lParse, "Dimensions of max(a,b) arguments "
				"are not compatible");
			YYERROR;
		      }
		   } else if (FSTRCMP((yyvsp[-4].str),"SETNULL(") == 0) {
		     if (OPER((yyvsp[-3].Node)) != CONST_OP
			 || SIZE((yyvsp[-3].Node)) != 1) {
		       yyerror(scanner, lParse, "SETNULL first argument must be a scalar constant");
		       YYERROR;
		     }
		     /* Make sure first arg is same type as second arg */
		     if ( TYPE((yyvsp[-3].Node)) != TYPE((yyvsp[-1].Node)) ) (yyvsp[-3].Node) = New_Unary(lParse,  TYPE((yyvsp[-1].Node)), 0, (yyvsp[-3].Node) );
		     (yyval.Node) = New_Func(lParse,  0, setnull_fct, 2, (yyvsp[-1].Node), (yyvsp[-3].Node), 0, 0, 0, 0, 0 );
		   } else if (FSTRCMP((yyvsp[-4].str),"AXISELEM(") == 0) {  /* AXISELEM(V,n) */
		     if (OPER((yyvsp[-1].Node)) != CONST_OP
			 || SIZE((yyvsp[-1].Node)) != 1) {
		       yyerror(scanner, lParse, "AXISELEM second argument must be a scalar constant");
		       YYERROR;
		     }
		     if (OPER((yyvsp[-3].Node)) == CONST_OP) {
		       long one = 1;
		       (yyval.Node) = New_Const(lParse,  LONG, &one, sizeof(one) );
		     } else {
		       if ( TYPE((yyvsp[-1].Node)) != LONG ) (yyvsp[-1].Node) = New_Unary(lParse, LONG, 0, (yyvsp[-1].Node));
		       (yyval.Node) = New_Func(lParse, 0, axiselem_fct, 2, (yyvsp[-3].Node), (yyvsp[-1].Node), 0, 0, 0, 0, 0 );
		       TEST((yyval.Node));
		       TYPE((yyval.Node)) = LONG;
		     }
		   } else if (FSTRCMP((yyvsp[-4].str),"NAXES(") == 0) {  /* NAXES(V,n) */
		     if (OPER((yyvsp[-1].Node)) != CONST_OP
			 || SIZE((yyvsp[-1].Node)) != 1) {
		       yyerror(scanner, lParse, "NAXES second argument must be a scalar constant");
		       YYERROR;
		     }
		     if (OPER((yyvsp[-3].Node)) == CONST_OP) { /* if V is constant, return 1 in every case */
		       long one = 1;
		       (yyval.Node) = New_Const(lParse,  LONG, &one, sizeof(one) );
		     } else {                    /* determine now the dimension of the expression */
		       long iaxis;
		       int naxis;
		       if ( TYPE((yyvsp[-1].Node)) != LONG ) (yyvsp[-1].Node) = New_Unary(lParse, LONG, 0, (yyvsp[-1].Node));
		       /* Since it is already constant, we can extract long value directly */
		       iaxis = (lParse->Nodes[(yyvsp[-1].Node)].value.data.lng);
		       naxis = lParse->Nodes[(yyvsp[-3].Node)].value.naxis;

		       if (iaxis == 0)          iaxis = naxis;   /* NAXIS(V,0) = NAXIS */
		       else if (iaxis <= naxis) iaxis = lParse->Nodes[(yyvsp[-3].Node)].value.naxes[iaxis-1]; /* NAXIS(V,n) = NAXISn */
		       else                     iaxis = 1;       /* Out of bounds use 1 */

		       (yyval.Node) = New_Const(lParse,  LONG, &iaxis, sizeof(iaxis) );
		       TEST((yyval.Node));
		     }
		   } else if (FSTRCMP((yyvsp[-4].str),"ARRAY(") == 0) {  /* NAXES(expr,n) */
		     (yyval.Node) = New_Array(lParse, (yyvsp[-3].Node), (yyvsp[-1].Node));
		     TEST((yyval.Node));
		   } else {
		      yyerror(scanner, lParse, "Function(expr,expr) not supported");
		      YYERROR;
		   }
                }
#line 2651 "eval_y.c"
    break;

  case 62: /* expr: FUNCTION expr ',' expr ',' expr ',' expr ')'  */
#line 886 "eval.y"
                { 
		  if (FSTRCMP((yyvsp[-8].str),"ANGSEP(") == 0) {
		    if( TYPE((yyvsp[-7].Node)) != DOUBLE ) (yyvsp[-7].Node) = New_Unary(lParse,  DOUBLE, 0, (yyvsp[-7].Node) );
		    if( TYPE((yyvsp[-5].Node)) != DOUBLE ) (yyvsp[-5].Node) = New_Unary(lParse,  DOUBLE, 0, (yyvsp[-5].Node) );
		    if( TYPE((yyvsp[-3].Node)) != DOUBLE ) (yyvsp[-3].Node) = New_Unary(lParse,  DOUBLE, 0, (yyvsp[-3].Node) );
		    if( TYPE((yyvsp[-1].Node)) != DOUBLE ) (yyvsp[-1].Node) = New_Unary(lParse,  DOUBLE, 0, (yyvsp[-1].Node) );
		    if( Test_Dims( lParse,  (yyvsp[-7].Node), (yyvsp[-5].Node) ) && Test_Dims( lParse,  (yyvsp[-5].Node), (yyvsp[-3].Node) ) && 
			Test_Dims( lParse,  (yyvsp[-3].Node), (yyvsp[-1].Node) ) ) {
		      (yyval.Node) = New_Func(lParse,  0, angsep_fct, 4, (yyvsp[-7].Node), (yyvsp[-5].Node), (yyvsp[-3].Node), (yyvsp[-1].Node),0,0,0 );
		      TEST((yyval.Node)); 
		      if( SIZE((yyvsp[-7].Node))<SIZE((yyvsp[-5].Node)) ) Copy_Dims( lParse,(yyval.Node), (yyvsp[-5].Node));
		      if( SIZE((yyvsp[-5].Node))<SIZE((yyvsp[-3].Node)) ) Copy_Dims( lParse,(yyval.Node), (yyvsp[-3].Node));
		      if( SIZE((yyvsp[-3].Node))<SIZE((yyvsp[-1].Node)) ) Copy_Dims( lParse,(yyval.Node), (yyvsp[-1].Node));
		    } else {
		      yyerror(scanner, lParse, "Dimensions of ANGSEP arguments "
			      "are not compatible");
		      YYERROR;
		    }
		   } else {
		      yyerror(scanner, lParse, "Function(expr,expr,expr,expr) not supported");
		      YYERROR;
		   }
                }
#line 2679 "eval_y.c"
    break;

  case 63: /* expr: expr '[' expr ']'  */
#line 910 "eval.y"
                { (yyval.Node) = New_Deref(lParse,  (yyvsp[-3].Node), 1, (yyvsp[-1].Node),  0,  0,  0,   0 ); TEST((yyval.Node)); }
#line 2685 "eval_y.c"
    break;

  case 64: /* expr: expr '[' expr ',' expr ']'  */
#line 912 "eval.y"
                { (yyval.Node) = New_Deref(lParse,  (yyvsp[-5].Node), 2, (yyvsp[-3].Node), (yyvsp[-1].Node),  0,  0,   0 ); TEST((yyval.Node)); }
#line 2691 "eval_y.c"
    break;

  case 65: /* expr: expr '[' expr ',' expr ',' expr ']'  */
#line 914 "eval.y"
                { (yyval.Node) = New_Deref(lParse,  (yyvsp[-7].Node), 3, (yyvsp[-5].Node), (yyvsp[-3].Node), (yyvsp[-1].Node),  0,   0 ); TEST((yyval.Node)); }
#line 2697 "eval_y.c"
    break;

  case 66: /* expr: expr '[' expr ',' expr ',' expr ',' expr ']'  */
#line 916 "eval.y"
                { (yyval.Node) = New_Deref(lParse,  (yyvsp[-9].Node), 4, (yyvsp[-7].Node), (yyvsp[-5].Node), (yyvsp[-3].Node), (yyvsp[-1].Node),   0 ); TEST((yyval.Node)); }
#line 2703 "eval_y.c"
    break;

  case 67: /* expr: expr '[' expr ',' expr ',' expr ',' expr ',' expr ']'  */
#line 918 "eval.y"
                { (yyval.Node) = New_Deref(lParse,  (yyvsp[-11].Node), 5, (yyvsp[-9].Node), (yyvsp[-7].Node), (yyvsp[-5].Node), (yyvsp[-3].Node), (yyvsp[-1].Node) ); TEST((yyval.Node)); }
#line 2709 "eval_y.c"
    break;

  case 68: /* expr: INTCAST expr  */
#line 920 "eval.y"
                { (yyval.Node) = New_Unary(lParse,  LONG,   INTCAST, (yyvsp[0].Node) );  TEST((yyval.Node));  }
#line 2715 "eval_y.c"
    break;

  case 69: /* expr: INTCAST bexpr  */
#line 922 "eval.y"
                { (yyval.Node) = New_Unary(lParse,  LONG,   INTCAST, (yyvsp[0].Node) );  TEST((yyval.Node));  }
#line 2721 "eval_y.c"
    break;

  case 70: /* expr: FLTCAST expr  */
#line 924 "eval.y"
                { (yyval.Node) = New_Unary(lParse,  DOUBLE, FLTCAST, (yyvsp[0].Node) );  TEST((yyval.Node));  }
#line 2727 "eval_y.c"
    break;

  case 71: /* expr: FLTCAST bexpr  */
#line 926 "eval.y"
                { (yyval.Node) = New_Unary(lParse,  DOUBLE, FLTCAST, (yyvsp[0].Node) );  TEST((yyval.Node));  }
#line 2733 "eval_y.c"
    break;

  case 72: /* bexpr: BOOLEAN  */
#line 930 "eval.y"
                { (yyval.Node) = New_Const(lParse,  BOOLEAN, &((yyvsp[0].log)), sizeof(char) ); TEST((yyval.Node)); }
#line 2739 "eval_y.c"
    break;

  case 73: /* bexpr: BCOLUMN  */
#line 932 "eval.y"
                { (yyval.Node) = New_Column(lParse,  (yyvsp[0].lng) ); TEST((yyval.Node)); }
#line 2745 "eval_y.c"
    break;

  case 74: /* bexpr: BCOLUMN '{' expr '}'  */
#line 934 "eval.y"
                {
                  if( TYPE((yyvsp[-1].Node)) != LONG
		      || OPER((yyvsp[-1].Node)) != CONST_OP ) {
		     yyerror(scanner, lParse, "Offset argument must be a constant integer");
		     YYERROR;
		  }
                  (yyval.Node) = New_Offset(lParse,  (yyvsp[-3].lng), (yyvsp[-1].Node) ); TEST((yyval.Node));
                }
#line 2758 "eval_y.c"
    break;

  case 75: /* bexpr: bits EQ bits  */
#line 943 "eval.y"
                { (yyval.Node) = New_BinOp(lParse,  BOOLEAN, (yyvsp[-2].Node), EQ,  (yyvsp[0].Node) ); TEST((yyval.Node));
		  SIZE((yyval.Node)) = 1;                                     }
#line 2765 "eval_y.c"
    break;

  case 76: /* bexpr: bits NE bits  */
#line 946 "eval.y"
                { (yyval.Node) = New_BinOp(lParse,  BOOLEAN, (yyvsp[-2].Node), NE,  (yyvsp[0].Node) ); TEST((yyval.Node)); 
		  SIZE((yyval.Node)) = 1;                                     }
#line 2772 "eval_y.c"
    break;

  case 77: /* bexpr: bits LT bits  */
#line 949 "eval.y"
                { (yyval.Node) = New_BinOp(lParse,  BOOLEAN, (yyvsp[-2].Node), LT,  (yyvsp[0].Node) ); TEST((yyval.Node)); 
		  SIZE((yyval.Node)) = 1;                                     }
#line 2779 "eval_y.c"
    break;

  case 78: /* bexpr: bits LTE bits  */
#line 952 "eval.y"
                { (yyval.Node) = New_BinOp(lParse,  BOOLEAN, (yyvsp[-2].Node), LTE, (yyvsp[0].Node) ); TEST((yyval.Node)); 
		  SIZE((yyval.Node)) = 1;                                     }
#line 2786 "eval_y.c"
    break;

  case 79: /* bexpr: bits GT bits  */
#line 955 "eval.y"
                { (yyval.Node) = New_BinOp(lParse,  BOOLEAN, (yyvsp[-2].Node), GT,  (yyvsp[0].Node) ); TEST((yyval.Node)); 
		  SIZE((yyval.Node)) = 1;                                     }
#line 2793 "eval_y.c"
    break;

  case 80: /* bexpr: bits GTE bits  */
#line 958 "eval.y"
                { (yyval.Node) = New_BinOp(lParse,  BOOLEAN, (yyvsp[-2].Node), GTE, (yyvsp[0].Node) ); TEST((yyval.Node)); 
		  SIZE((yyval.Node)) = 1;                                     }
#line 2800 "eval_y.c"
    break;

  case 81: /* bexpr: expr GT expr  */
#line 961 "eval.y"
                { PROMOTE((yyvsp[-2].Node),(yyvsp[0].Node)); (yyval.Node) = New_BinOp(lParse,  BOOLEAN, (yyvsp[-2].Node), GT,  (yyvsp[0].Node) );
                  TEST((yyval.Node));                                               }
#line 2807 "eval_y.c"
    break;

  case 82: /* bexpr: expr LT expr  */
#line 964 "eval.y"
                { PROMOTE((yyvsp[-2].Node),(yyvsp[0].Node)); (yyval.Node) = New_BinOp(lParse,  BOOLEAN, (yyvsp[-2].Node), LT,  (yyvsp[0].Node) );
                  TEST((yyval.Node));                                               }
#line 2814 "eval_y.c"
    break;

  case 83: /* bexpr: expr GTE expr  */
#line 967 "eval.y"
                { PROMOTE((yyvsp[-2].Node),(yyvsp[0].Node)); (yyval.Node) = New_BinOp(lParse,  BOOLEAN, (yyvsp[-2].Node), GTE, (yyvsp[0].Node) );
                  TEST((yyval.Node));                                               }
#line 2821 "eval_y.c"
    break;

  case 84: /* bexpr: expr LTE expr  */
#line 970 "eval.y"
                { PROMOTE((yyvsp[-2].Node),(yyvsp[0].Node)); (yyval.Node) = New_BinOp(lParse,  BOOLEAN, (yyvsp[-2].Node), LTE, (yyvsp[0].Node) );
                  TEST((yyval.Node));                                               }
#line 2828 "eval_y.c"
    break;

  case 85: /* bexpr: expr '~' expr  */
#line 973 "eval.y"
                { PROMOTE((yyvsp[-2].Node),(yyvsp[0].Node)); (yyval.Node) = New_BinOp(lParse,  BOOLEAN, (yyvsp[-2].Node), '~', (yyvsp[0].Node) );
                  TEST((yyval.Node));                                               }
#line 2835 "eval_y.c"
    break;

  case 86: /* bexpr: expr EQ expr  */
#line 976 "eval.y"
                { PROMOTE((yyvsp[-2].Node),(yyvsp[0].Node)); (yyval.Node) = New_BinOp(lParse,  BOOLEAN, (yyvsp[-2].Node), EQ,  (yyvsp[0].Node) );
                  TEST((yyval.Node));                                               }
#line 2842 "eval_y.c"
    break;

  case 87: /* bexpr: expr NE expr  */
#line 979 "eval.y"
                { PROMOTE((yyvsp[-2].Node),(yyvsp[0].Node)); (yyval.Node) = New_BinOp(lParse,  BOOLEAN, (yyvsp[-2].Node), NE,  (yyvsp[0].Node) );
                  TEST((yyval.Node));                                               }
#line 2849 "eval_y.c"
    break;

  case 88: /* bexpr: sexpr EQ sexpr  */
#line 982 "eval.y"
                { (yyval.Node) = New_BinOp(lParse,  BOOLEAN, (yyvsp[-2].Node), EQ,  (yyvsp[0].Node) ); TEST((yyval.Node));
                  SIZE((yyval.Node)) = 1; }
#line 2856 "eval_y.c"
    break;

  case 89: /* bexpr: sexpr NE sexpr  */
#line 985 "eval.y"
                { (yyval.Node) = New_BinOp(lParse,  BOOLEAN, (yyvsp[-2].Node), NE,  (yyvsp[0].Node) ); TEST((yyval.Node));
                  SIZE((yyval.Node)) = 1; }
#line 2863 "eval_y.c"
    break;

  case 90: /* bexpr: sexpr GT sexpr  */
#line 988 "eval.y"
                { (yyval.Node) = New_BinOp(lParse,  BOOLEAN, (yyvsp[-2].Node), GT,  (yyvsp[0].Node) ); TEST((yyval.Node));
                  SIZE((yyval.Node)) = 1; }
#line 2870 "eval_y.c"
    break;

  case 91: /* bexpr: sexpr GTE sexpr  */
#line 991 "eval.y"
                { (yyval.Node) = New_BinOp(lParse,  BOOLEAN, (yyvsp[-2].Node), GTE, (yyvsp[0].Node) ); TEST((yyval.Node));
                  SIZE((yyval.Node)) = 1; }
#line 2877 "eval_y.c"
    break;

  case 92: /* bexpr: sexpr LT sexpr  */
#line 994 "eval.y"
                { (yyval.Node) = New_BinOp(lParse,  BOOLEAN, (yyvsp[-2].Node), LT,  (yyvsp[0].Node) ); TEST((yyval.Node));
                  SIZE((yyval.Node)) = 1; }
#line 2884 "eval_y.c"
    break;

  case 93: /* bexpr: sexpr LTE sexpr  */
#line 997 "eval.y"
                { (yyval.Node) = New_BinOp(lParse,  BOOLEAN, (yyvsp[-2].Node), LTE, (yyvsp[0].Node) ); TEST((yyval.Node));
                  SIZE((yyval.Node)) = 1; }
#line 2891 "eval_y.c"
    break;

  case 94: /* bexpr: bexpr AND bexpr  */
#line 1000 "eval.y"
                { (yyval.Node) = New_BinOp(lParse,  BOOLEAN, (yyvsp[-2].Node), AND, (yyvsp[0].Node) ); TEST((yyval.Node)); }
#line 2897 "eval_y.c"
    break;

  case 95: /* bexpr: bexpr OR bexpr  */
#line 1002 "eval.y"
                { (yyval.Node) = New_BinOp(lParse,  BOOLEAN, (yyvsp[-2].Node), OR,  (yyvsp[0].Node) ); TEST((yyval.Node)); }
#line 2903 "eval_y.c"
    break;

  case 96: /* bexpr: bexpr EQ bexpr  */
#line 1004 "eval.y"
                { (yyval.Node) = New_BinOp(lParse,  BOOLEAN, (yyvsp[-2].Node), EQ,  (yyvsp[0].Node) ); TEST((yyval.Node)); }
#line 2909 "eval_y.c"
    break;

  case 97: /* bexpr: bexpr NE bexpr  */
#line 1006 "eval.y"
                { (yyval.Node) = New_BinOp(lParse,  BOOLEAN, (yyvsp[-2].Node), NE,  (yyvsp[0].Node) ); TEST((yyval.Node)); }
#line 2915 "eval_y.c"
    break;

  case 98: /* bexpr: expr '=' expr ':' expr  */
#line 1009 "eval.y"
                { PROMOTE((yyvsp[-4].Node),(yyvsp[-2].Node)); PROMOTE((yyvsp[-4].Node),(yyvsp[0].Node)); PROMOTE((yyvsp[-2].Node),(yyvsp[0].Node));
		  (yyvsp[-2].Node) = New_BinOp(lParse,  BOOLEAN, (yyvsp[-2].Node), LTE, (yyvsp[-4].Node) );
                  (yyvsp[0].Node) = New_BinOp(lParse,  BOOLEAN, (yyvsp[-4].Node), LTE, (yyvsp[0].Node) );
                  (yyval.Node) = New_BinOp(lParse,  BOOLEAN, (yyvsp[-2].Node), AND, (yyvsp[0].Node) );
                  TEST((yyval.Node));                                         }
#line 2925 "eval_y.c"
    break;

  case 99: /* bexpr: bexpr '?' bexpr ':' bexpr  */
#line 1016 "eval.y"
                {
                  if( ! Test_Dims( lParse, (yyvsp[-2].Node),(yyvsp[0].Node)) ) {
                     yyerror(scanner, lParse, "Incompatible dimensions in '?:' arguments");
		     YYERROR;
                  }
                  (yyval.Node) = New_Func(lParse,  0, ifthenelse_fct, 3, (yyvsp[-2].Node), (yyvsp[0].Node), (yyvsp[-4].Node),
                                 0, 0, 0, 0 );
                  TEST((yyval.Node));
                  if( SIZE((yyvsp[-2].Node))<SIZE((yyvsp[0].Node)) )  Copy_Dims( lParse,(yyval.Node), (yyvsp[0].Node));
                  if( ! Test_Dims( lParse, (yyvsp[-4].Node),(yyval.Node)) ) {
                     yyerror(scanner, lParse, "Incompatible dimensions in '?:' condition");
		     YYERROR;
                  }
                  if( SIZE((yyval.Node))<SIZE((yyvsp[-4].Node)) )  Copy_Dims( lParse,(yyval.Node), (yyvsp[-4].Node));
                }
#line 2945 "eval_y.c"
    break;

  case 100: /* bexpr: BFUNCTION expr ')'  */
#line 1033 "eval.y"
                {
		   if (FSTRCMP((yyvsp[-2].str),"ISNULL(") == 0) {
		      (yyval.Node) = New_Func(lParse,  0, isnull_fct, 1, (yyvsp[-1].Node), 0, 0,
				     0, 0, 0, 0 );
		      TEST((yyval.Node)); 
                      /* Use expression's size, but return BOOLEAN */
		      TYPE((yyval.Node)) = BOOLEAN;
		   } else {
		      yyerror(scanner, lParse, "Boolean Function(expr) not supported");
		      YYERROR;
		   }
		}
#line 2962 "eval_y.c"
    break;

  case 101: /* bexpr: BFUNCTION bexpr ')'  */
#line 1046 "eval.y"
                {
		   if (FSTRCMP((yyvsp[-2].str),"ISNULL(") == 0) {
		      (yyval.Node) = New_Func(lParse,  0, isnull_fct, 1, (yyvsp[-1].Node), 0, 0,
				     0, 0, 0, 0 );
		      TEST((yyval.Node)); 
                      /* Use expression's size, but return BOOLEAN */
		      TYPE((yyval.Node)) = BOOLEAN;
		   } else {
		      yyerror(scanner, lParse, "Boolean Function(expr) not supported");
		      YYERROR;
		   }
		}
#line 2979 "eval_y.c"
    break;

  case 102: /* bexpr: BFUNCTION sexpr ')'  */
#line 1059 "eval.y"
                {
		   if (FSTRCMP((yyvsp[-2].str),"ISNULL(") == 0) {
		      (yyval.Node) = New_Func(lParse,  BOOLEAN, isnull_fct, 1, (yyvsp[-1].Node), 0, 0,
				     0, 0, 0, 0 );
		      TEST((yyval.Node)); 
		   } else {
		      yyerror(scanner, lParse, "Boolean Function(expr) not supported");
		      YYERROR;
		   }
		}
#line 2994 "eval_y.c"
    break;

  case 103: /* bexpr: FUNCTION bexpr ',' bexpr ')'  */
#line 1070 "eval.y"
                {
		   if (FSTRCMP((yyvsp[-4].str),"DEFNULL(") == 0) {
		      if( SIZE((yyvsp[-3].Node))>=SIZE((yyvsp[-1].Node)) && Test_Dims( lParse,  (yyvsp[-3].Node), (yyvsp[-1].Node) ) ) {
			 (yyval.Node) = New_Func(lParse,  0, defnull_fct, 2, (yyvsp[-3].Node), (yyvsp[-1].Node), 0,
					0, 0, 0, 0 );
			 TEST((yyval.Node)); 
		      } else {
			 yyerror(scanner, lParse, "Dimensions of DEFNULL arguments are not compatible");
			 YYERROR;
		      }
		   } else {
		      yyerror(scanner, lParse, "Boolean Function(expr,expr) not supported");
		      YYERROR;
		   }
		}
#line 3014 "eval_y.c"
    break;

  case 104: /* bexpr: BFUNCTION expr ',' expr ',' expr ')'  */
#line 1086 "eval.y"
                {
		   if( TYPE((yyvsp[-5].Node)) != DOUBLE ) (yyvsp[-5].Node) = New_Unary(lParse,  DOUBLE, 0, (yyvsp[-5].Node) );
		   if( TYPE((yyvsp[-3].Node)) != DOUBLE ) (yyvsp[-3].Node) = New_Unary(lParse,  DOUBLE, 0, (yyvsp[-3].Node) );
		   if( TYPE((yyvsp[-1].Node)) != DOUBLE ) (yyvsp[-1].Node) = New_Unary(lParse,  DOUBLE, 0, (yyvsp[-1].Node) );
		   if( ! (Test_Dims( lParse,  (yyvsp[-5].Node), (yyvsp[-3].Node) ) && Test_Dims( lParse,  (yyvsp[-3].Node), (yyvsp[-1].Node) ) ) ) {
		       yyerror(scanner, lParse, "Dimensions of NEAR arguments "
			       "are not compatible");
		       YYERROR;
		   } else {
		     if (FSTRCMP((yyvsp[-6].str),"NEAR(") == 0) {
		       (yyval.Node) = New_Func(lParse,  BOOLEAN, near_fct, 3, (yyvsp[-5].Node), (yyvsp[-3].Node), (yyvsp[-1].Node),
				      0, 0, 0, 0 );
		     } else {
		       yyerror(scanner, lParse, "Boolean Function not supported");
		       YYERROR;
		     }
		     TEST((yyval.Node)); 

		     if( SIZE((yyval.Node))<SIZE((yyvsp[-5].Node)) )  Copy_Dims( lParse,(yyval.Node), (yyvsp[-5].Node));
		     if( SIZE((yyvsp[-5].Node))<SIZE((yyvsp[-3].Node)) )  Copy_Dims( lParse,(yyval.Node), (yyvsp[-3].Node));
		     if( SIZE((yyvsp[-3].Node))<SIZE((yyvsp[-1].Node)) )  Copy_Dims( lParse,(yyval.Node), (yyvsp[-1].Node));
		   }
		}
#line 3042 "eval_y.c"
    break;

  case 105: /* bexpr: BFUNCTION expr ',' expr ',' expr ',' expr ',' expr ')'  */
#line 1110 "eval.y"
                {
		   if( TYPE((yyvsp[-9].Node)) != DOUBLE ) (yyvsp[-9].Node) = New_Unary(lParse,  DOUBLE, 0, (yyvsp[-9].Node) );
		   if( TYPE((yyvsp[-7].Node)) != DOUBLE ) (yyvsp[-7].Node) = New_Unary(lParse,  DOUBLE, 0, (yyvsp[-7].Node) );
		   if( TYPE((yyvsp[-5].Node)) != DOUBLE ) (yyvsp[-5].Node) = New_Unary(lParse,  DOUBLE, 0, (yyvsp[-5].Node) );
		   if( TYPE((yyvsp[-3].Node)) != DOUBLE ) (yyvsp[-3].Node) = New_Unary(lParse,  DOUBLE, 0, (yyvsp[-3].Node) );
		   if( TYPE((yyvsp[-1].Node))!= DOUBLE ) (yyvsp[-1].Node)= New_Unary(lParse,  DOUBLE, 0, (yyvsp[-1].Node));
		   if( ! (Test_Dims( lParse,  (yyvsp[-9].Node), (yyvsp[-7].Node) ) && Test_Dims( lParse,  (yyvsp[-7].Node), (yyvsp[-5].Node) ) && 
			  Test_Dims( lParse,  (yyvsp[-5].Node), (yyvsp[-3].Node) ) && Test_Dims( lParse,  (yyvsp[-3].Node), (yyvsp[-1].Node) )) ) {
		     yyerror(scanner, lParse, "Dimensions of CIRCLE arguments "
			     "are not compatible");
		     YYERROR;
		   } else {
		     if (FSTRCMP((yyvsp[-10].str),"CIRCLE(") == 0) {
		       (yyval.Node) = New_Func(lParse,  BOOLEAN, circle_fct, 5, (yyvsp[-9].Node), (yyvsp[-7].Node), (yyvsp[-5].Node), (yyvsp[-3].Node),
				      (yyvsp[-1].Node), 0, 0 );
		     } else {
		       yyerror(scanner, lParse, "Boolean Function not supported");
		       YYERROR;
		     }
		     TEST((yyval.Node)); 
		     if( SIZE((yyval.Node))<SIZE((yyvsp[-9].Node)) )  Copy_Dims( lParse,(yyval.Node), (yyvsp[-9].Node));
		     if( SIZE((yyvsp[-9].Node))<SIZE((yyvsp[-7].Node)) )  Copy_Dims( lParse,(yyval.Node), (yyvsp[-7].Node));
		     if( SIZE((yyvsp[-7].Node))<SIZE((yyvsp[-5].Node)) )  Copy_Dims( lParse,(yyval.Node), (yyvsp[-5].Node));
		     if( SIZE((yyvsp[-5].Node))<SIZE((yyvsp[-3].Node)) )  Copy_Dims( lParse,(yyval.Node), (yyvsp[-3].Node));
		     if( SIZE((yyvsp[-3].Node))<SIZE((yyvsp[-1].Node)) ) Copy_Dims( lParse,(yyval.Node), (yyvsp[-1].Node));
		   }
		}
#line 3074 "eval_y.c"
    break;

  case 106: /* bexpr: BFUNCTION expr ',' expr ',' expr ',' expr ',' expr ',' expr ',' expr ')'  */
#line 1138 "eval.y"
                {
		   if( TYPE((yyvsp[-13].Node)) != DOUBLE ) (yyvsp[-13].Node) = New_Unary(lParse,  DOUBLE, 0, (yyvsp[-13].Node) );
		   if( TYPE((yyvsp[-11].Node)) != DOUBLE ) (yyvsp[-11].Node) = New_Unary(lParse,  DOUBLE, 0, (yyvsp[-11].Node) );
		   if( TYPE((yyvsp[-9].Node)) != DOUBLE ) (yyvsp[-9].Node) = New_Unary(lParse,  DOUBLE, 0, (yyvsp[-9].Node) );
		   if( TYPE((yyvsp[-7].Node)) != DOUBLE ) (yyvsp[-7].Node) = New_Unary(lParse,  DOUBLE, 0, (yyvsp[-7].Node) );
		   if( TYPE((yyvsp[-5].Node))!= DOUBLE ) (yyvsp[-5].Node)= New_Unary(lParse,  DOUBLE, 0, (yyvsp[-5].Node));
		   if( TYPE((yyvsp[-3].Node))!= DOUBLE ) (yyvsp[-3].Node)= New_Unary(lParse,  DOUBLE, 0, (yyvsp[-3].Node));
		   if( TYPE((yyvsp[-1].Node))!= DOUBLE ) (yyvsp[-1].Node)= New_Unary(lParse,  DOUBLE, 0, (yyvsp[-1].Node));
		   if( ! (Test_Dims( lParse,  (yyvsp[-13].Node), (yyvsp[-11].Node) ) && Test_Dims( lParse,  (yyvsp[-11].Node), (yyvsp[-9].Node) ) && 
			  Test_Dims( lParse,  (yyvsp[-9].Node), (yyvsp[-7].Node) ) && Test_Dims( lParse,  (yyvsp[-7].Node), (yyvsp[-5].Node) ) &&
			  Test_Dims( lParse, (yyvsp[-5].Node),(yyvsp[-3].Node) ) && Test_Dims( lParse, (yyvsp[-3].Node), (yyvsp[-1].Node) ) ) ) {
		     yyerror(scanner, lParse, "Dimensions of BOX or ELLIPSE arguments "
			     "are not compatible");
		     YYERROR;
		   } else {
		     if (FSTRCMP((yyvsp[-14].str),"BOX(") == 0) {
		       (yyval.Node) = New_Func(lParse,  BOOLEAN, box_fct, 7, (yyvsp[-13].Node), (yyvsp[-11].Node), (yyvsp[-9].Node), (yyvsp[-7].Node),
				      (yyvsp[-5].Node), (yyvsp[-3].Node), (yyvsp[-1].Node) );
		     } else if (FSTRCMP((yyvsp[-14].str),"ELLIPSE(") == 0) {
		       (yyval.Node) = New_Func(lParse,  BOOLEAN, elps_fct, 7, (yyvsp[-13].Node), (yyvsp[-11].Node), (yyvsp[-9].Node), (yyvsp[-7].Node),
				      (yyvsp[-5].Node), (yyvsp[-3].Node), (yyvsp[-1].Node) );
		     } else {
		       yyerror(scanner, lParse, "SAO Image Function not supported");
		       YYERROR;
		     }
		     TEST((yyval.Node)); 
		     if( SIZE((yyval.Node))<SIZE((yyvsp[-13].Node)) )  Copy_Dims( lParse,(yyval.Node), (yyvsp[-13].Node));
		     if( SIZE((yyvsp[-13].Node))<SIZE((yyvsp[-11].Node)) )  Copy_Dims( lParse,(yyval.Node), (yyvsp[-11].Node));
		     if( SIZE((yyvsp[-11].Node))<SIZE((yyvsp[-9].Node)) )  Copy_Dims( lParse,(yyval.Node), (yyvsp[-9].Node));
		     if( SIZE((yyvsp[-9].Node))<SIZE((yyvsp[-7].Node)) )  Copy_Dims( lParse,(yyval.Node), (yyvsp[-7].Node));
		     if( SIZE((yyvsp[-7].Node))<SIZE((yyvsp[-5].Node)) ) Copy_Dims( lParse,(yyval.Node), (yyvsp[-5].Node));
		     if( SIZE((yyvsp[-5].Node))<SIZE((yyvsp[-3].Node)) ) Copy_Dims( lParse,(yyval.Node), (yyvsp[-3].Node));
		     if( SIZE((yyvsp[-3].Node))<SIZE((yyvsp[-1].Node)) ) Copy_Dims( lParse,(yyval.Node), (yyvsp[-1].Node));
		   }
		}
#line 3114 "eval_y.c"
    break;

  case 107: /* bexpr: GTIFILTER ')'  */
#line 1175 "eval.y"
                { /* Use defaults for all elements */
		   (yyval.Node) = New_GTI(lParse, gtifilt_fct,  "", -99, -99, "*START*", "*STOP*" );
                   TEST((yyval.Node));                                        }
#line 3122 "eval_y.c"
    break;

  case 108: /* bexpr: GTIFILTER STRING ')'  */
#line 1179 "eval.y"
                { /* Use defaults for all except filename */
		  (yyval.Node) = New_GTI(lParse, gtifilt_fct,  (yyvsp[-1].str), -99, -99, "*START*", "*STOP*" );
                   TEST((yyval.Node));                                        }
#line 3130 "eval_y.c"
    break;

  case 109: /* bexpr: GTIFILTER STRING ',' expr ')'  */
#line 1183 "eval.y"
                {  (yyval.Node) = New_GTI(lParse, gtifilt_fct,  (yyvsp[-3].str), (yyvsp[-1].Node), -99, "*START*", "*STOP*" );
                   TEST((yyval.Node));                                        }
#line 3137 "eval_y.c"
    break;

  case 110: /* bexpr: GTIFILTER STRING ',' expr ',' STRING ',' STRING ')'  */
#line 1186 "eval.y"
                {  (yyval.Node) = New_GTI(lParse, gtifilt_fct,  (yyvsp[-7].str), (yyvsp[-5].Node), -99, (yyvsp[-3].str), (yyvsp[-1].str) );
                   TEST((yyval.Node));                                        }
#line 3144 "eval_y.c"
    break;

  case 111: /* bexpr: GTIOVERLAP STRING ',' expr ',' expr ')'  */
#line 1191 "eval.y"
                {  (yyval.Node) = New_GTI(lParse, gtiover_fct,  (yyvsp[-5].str), (yyvsp[-3].Node), (yyvsp[-1].Node), "*START*", "*STOP*");
                   TEST((yyval.Node));                                        }
#line 3151 "eval_y.c"
    break;

  case 112: /* bexpr: GTIOVERLAP STRING ',' expr ',' expr ',' STRING ',' STRING ')'  */
#line 1194 "eval.y"
                {  (yyval.Node) = New_GTI(lParse, gtiover_fct,  (yyvsp[-9].str), (yyvsp[-7].Node), (yyvsp[-5].Node), (yyvsp[-3].str), (yyvsp[-1].str) );
                   TEST((yyval.Node));                                        }
#line 3158 "eval_y.c"
    break;

  case 113: /* bexpr: GTIFIND ')'  */
#line 1199 "eval.y"
                { /* Use defaults for all elements */
		   (yyval.Node) = New_GTI(lParse, gtifind_fct,  "", -99, -99, "*START*", "*STOP*" );
                   TEST((yyval.Node));                                        }
#line 3166 "eval_y.c"
    break;

  case 114: /* bexpr: GTIFIND STRING ')'  */
#line 1203 "eval.y"
                { /* Use defaults for all except filename */
		  (yyval.Node) = New_GTI(lParse, gtifind_fct,  (yyvsp[-1].str), -99, -99, "*START*", "*STOP*" );
                   TEST((yyval.Node));                                        }
#line 3174 "eval_y.c"
    break;

  case 115: /* bexpr: GTIFIND STRING ',' expr ')'  */
#line 1207 "eval.y"
                {  (yyval.Node) = New_GTI(lParse, gtifind_fct,  (yyvsp[-3].str), (yyvsp[-1].Node), -99, "*START*", "*STOP*" );
                   TEST((yyval.Node));                                        }
#line 3181 "eval_y.c"
    break;

  case 116: /* bexpr: GTIFIND STRING ',' expr ',' STRING ',' STRING ')'  */
#line 1210 "eval.y"
                {  (yyval.Node) = New_GTI(lParse, gtifind_fct,  (yyvsp[-7].str), (yyvsp[-5].Node), -99, (yyvsp[-3].str), (yyvsp[-1].str) );
                   TEST((yyval.Node));                                        }
#line 3188 "eval_y.c"
    break;

  case 117: /* bexpr: REGFILTER STRING ')'  */
#line 1215 "eval.y"
                { /* Use defaults for all except filename */
                   (yyval.Node) = New_REG(lParse,  (yyvsp[-1].str), -99, -99, "" );
                   TEST((yyval.Node));                                        }
#line 3196 "eval_y.c"
    break;

  case 118: /* bexpr: REGFILTER STRING ',' expr ',' expr ')'  */
#line 1219 "eval.y"
                {  (yyval.Node) = New_REG(lParse,  (yyvsp[-5].str), (yyvsp[-3].Node), (yyvsp[-1].Node), "" );
                   TEST((yyval.Node));                                        }
#line 3203 "eval_y.c"
    break;

  case 119: /* bexpr: REGFILTER STRING ',' expr ',' expr ',' STRING ')'  */
#line 1222 "eval.y"
                {  (yyval.Node) = New_REG(lParse,  (yyvsp[-7].str), (yyvsp[-5].Node), (yyvsp[-3].Node), (yyvsp[-1].str) );
                   TEST((yyval.Node));                                        }
#line 3210 "eval_y.c"
    break;

  case 120: /* bexpr: bexpr '[' expr ']'  */
#line 1226 "eval.y"
                { (yyval.Node) = New_Deref(lParse,  (yyvsp[-3].Node), 1, (yyvsp[-1].Node),  0,  0,  0,   0 ); TEST((yyval.Node)); }
#line 3216 "eval_y.c"
    break;

  case 121: /* bexpr: bexpr '[' expr ',' expr ']'  */
#line 1228 "eval.y"
                { (yyval.Node) = New_Deref(lParse,  (yyvsp[-5].Node), 2, (yyvsp[-3].Node), (yyvsp[-1].Node),  0,  0,   0 ); TEST((yyval.Node)); }
#line 3222 "eval_y.c"
    break;

  case 122: /* bexpr: bexpr '[' expr ',' expr ',' expr ']'  */
#line 1230 "eval.y"
                { (yyval.Node) = New_Deref(lParse,  (yyvsp[-7].Node), 3, (yyvsp[-5].Node), (yyvsp[-3].Node), (yyvsp[-1].Node),  0,   0 ); TEST((yyval.Node)); }
#line 3228 "eval_y.c"
    break;

  case 123: /* bexpr: bexpr '[' expr ',' expr ',' expr ',' expr ']'  */
#line 1232 "eval.y"
                { (yyval.Node) = New_Deref(lParse,  (yyvsp[-9].Node), 4, (yyvsp[-7].Node), (yyvsp[-5].Node), (yyvsp[-3].Node), (yyvsp[-1].Node),   0 ); TEST((yyval.Node)); }
#line 3234 "eval_y.c"
    break;

  case 124: /* bexpr: bexpr '[' expr ',' expr ',' expr ',' expr ',' expr ']'  */
#line 1234 "eval.y"
                { (yyval.Node) = New_Deref(lParse,  (yyvsp[-11].Node), 5, (yyvsp[-9].Node), (yyvsp[-7].Node), (yyvsp[-5].Node), (yyvsp[-3].Node), (yyvsp[-1].Node) ); TEST((yyval.Node)); }
#line 3240 "eval_y.c"
    break;

  case 125: /* bexpr: NOT bexpr  */
#line 1236 "eval.y"
                { (yyval.Node) = New_Unary(lParse,  BOOLEAN, NOT, (yyvsp[0].Node) ); TEST((yyval.Node)); }
#line 3246 "eval_y.c"
    break;

  case 126: /* bexpr: '(' bexpr ')'  */
#line 1238 "eval.y"
                { (yyval.Node) = (yyvsp[-1].Node); }
#line 3252 "eval_y.c"
    break;

  case 127: /* sexpr: STRING  */
#line 1242 "eval.y"
                { (yyval.Node) = New_Const(lParse,  STRING, (yyvsp[0].str), strlen((yyvsp[0].str))+1 ); TEST((yyval.Node));
                  SIZE((yyval.Node)) = strlen((yyvsp[0].str)); }
#line 3259 "eval_y.c"
    break;

  case 128: /* sexpr: SCOLUMN  */
#line 1245 "eval.y"
                { (yyval.Node) = New_Column(lParse,  (yyvsp[0].lng) ); TEST((yyval.Node)); }
#line 3265 "eval_y.c"
    break;

  case 129: /* sexpr: SCOLUMN '{' expr '}'  */
#line 1247 "eval.y"
                {
                  if( TYPE((yyvsp[-1].Node)) != LONG
		      || OPER((yyvsp[-1].Node)) != CONST_OP ) {
		     yyerror(scanner, lParse, "Offset argument must be a constant integer");
		     YYERROR;
		  }
                  (yyval.Node) = New_Offset(lParse,  (yyvsp[-3].lng), (yyvsp[-1].Node) ); TEST((yyval.Node));
                }
#line 3278 "eval_y.c"
    break;

  case 130: /* sexpr: SNULLREF  */
#line 1256 "eval.y"
                { (yyval.Node) = New_Func(lParse,  STRING, null_fct, 0, 0, 0, 0, 0, 0, 0, 0 ); }
#line 3284 "eval_y.c"
    break;

  case 131: /* sexpr: '(' sexpr ')'  */
#line 1258 "eval.y"
                { (yyval.Node) = (yyvsp[-1].Node); }
#line 3290 "eval_y.c"
    break;

  case 132: /* sexpr: sexpr '+' sexpr  */
#line 1260 "eval.y"
                { 
		  if (SIZE((yyvsp[-2].Node))+SIZE((yyvsp[0].Node)) >= MAX_STRLEN) {
		    yyerror(scanner, lParse, "Combined string size exceeds " MAX_STRLEN_S " characters");
		    YYERROR;
		  }
		  (yyval.Node) = New_BinOp(lParse,  STRING, (yyvsp[-2].Node), '+', (yyvsp[0].Node) );  TEST((yyval.Node));
		  SIZE((yyval.Node)) = SIZE((yyvsp[-2].Node)) + SIZE((yyvsp[0].Node));
		}
#line 3303 "eval_y.c"
    break;

  case 133: /* sexpr: bexpr '?' sexpr ':' sexpr  */
#line 1269 "eval.y"
                {
		  int outSize;
                  if( SIZE((yyvsp[-4].Node))!=1 ) {
                     yyerror(scanner, lParse, "Cannot have a vector string column");
		     YYERROR;
                  }
		  /* Since the output can be calculated now, as a constant
		     scalar, we must precalculate the output size, in
		     order to avoid an overflow. */
		  outSize = SIZE((yyvsp[-2].Node));
		  if (SIZE((yyvsp[0].Node)) > outSize) outSize = SIZE((yyvsp[0].Node));
                  (yyval.Node) = New_FuncSize(lParse,  0, ifthenelse_fct, 3, (yyvsp[-2].Node), (yyvsp[0].Node), (yyvsp[-4].Node),
				     0, 0, 0, 0, outSize);
		  
                  TEST((yyval.Node));
                  if( SIZE((yyvsp[-2].Node))<SIZE((yyvsp[0].Node)) )  Copy_Dims( lParse,(yyval.Node), (yyvsp[0].Node));
                }
#line 3325 "eval_y.c"
    break;

  case 134: /* sexpr: FUNCTION sexpr ',' sexpr ')'  */
#line 1288 "eval.y"
                { 
		  if (FSTRCMP((yyvsp[-4].str),"DEFNULL(") == 0) {
		     int outSize;
		     /* Since the output can be calculated now, as a constant
			scalar, we must precalculate the output size, in
			order to avoid an overflow. */
		     outSize = SIZE((yyvsp[-3].Node));
		     if (SIZE((yyvsp[-1].Node)) > outSize) outSize = SIZE((yyvsp[-1].Node));
		     
		     (yyval.Node) = New_FuncSize(lParse,  0, defnull_fct, 2, (yyvsp[-3].Node), (yyvsp[-1].Node), 0,
					0, 0, 0, 0, outSize );
		     TEST((yyval.Node)); 
		     if( SIZE((yyvsp[-1].Node))>SIZE((yyvsp[-3].Node)) ) SIZE((yyval.Node)) = SIZE((yyvsp[-1].Node));
		  } else {
		     yyerror(scanner, lParse, "Function(string,string) not supported");
		     YYERROR;
		  }
		}
#line 3348 "eval_y.c"
    break;

  case 135: /* sexpr: FUNCTION sexpr ',' expr ',' expr ')'  */
#line 1307 "eval.y"
                { 
		  if (FSTRCMP((yyvsp[-6].str),"STRMID(") == 0) {
		    int len;
		    if( TYPE((yyvsp[-3].Node)) != LONG || SIZE((yyvsp[-3].Node)) != 1 ||
			TYPE((yyvsp[-1].Node)) != LONG || SIZE((yyvsp[-1].Node)) != 1) {
		      yyerror(scanner, lParse, "When using STRMID(S,P,N), P and N must be integers (and not vector columns)");
		      YYERROR;
		    }
		    if (OPER((yyvsp[-1].Node)) == CONST_OP) {
		      /* Constant value: use that directly */
		      len = (lParse->Nodes[(yyvsp[-1].Node)].value.data.lng);
		    } else {
		      /* Variable value: use the maximum possible (from $2) */
		      len = SIZE((yyvsp[-5].Node));
		    }
		    if (len <= 0 || len >= MAX_STRLEN) {
		      yyerror(scanner, lParse, "STRMID(S,P,N), N must be 1-" MAX_STRLEN_S);
		      YYERROR;
		    }
		    (yyval.Node) = New_FuncSize(lParse,  0, strmid_fct, 3, (yyvsp[-5].Node), (yyvsp[-3].Node),(yyvsp[-1].Node),0,0,0,0,len);
		    TEST((yyval.Node));
		  } else {
		     yyerror(scanner, lParse, "Function(string,expr,expr) not supported");
		     YYERROR;
		  }
		}
#line 3379 "eval_y.c"
    break;


#line 3383 "eval_y.c"

      default: break;
    }
  /* User semantic actions sometimes alter yychar, and that requires
     that yytoken be updated with the new translation.  We take the
     approach of translating immediately before every use of yytoken.
     One alternative is translating here after every semantic action,
     but that translation would be missed if the semantic action invokes
     YYABORT, YYACCEPT, or YYERROR immediately after altering yychar or
     if it invokes YYBACKUP.  In the case of YYABORT or YYACCEPT, an
     incorrect destructor might then be invoked immediately.  In the
     case of YYERROR or YYBACKUP, subsequent parser actions might lead
     to an incorrect destructor call or verbose syntax error message
     before the lookahead is translated.  */
  YY_SYMBOL_PRINT ("-> $$ =", YY_CAST (yysymbol_kind_t, yyr1[yyn]), &yyval, &yyloc);

  YYPOPSTACK (yylen);
  yylen = 0;

  *++yyvsp = yyval;

  /* Now 'shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */
  {
    const int yylhs = yyr1[yyn] - YYNTOKENS;
    const int yyi = yypgoto[yylhs] + *yyssp;
    yystate = (0 <= yyi && yyi <= YYLAST && yycheck[yyi] == *yyssp
               ? yytable[yyi]
               : yydefgoto[yylhs]);
  }

  goto yynewstate;


/*--------------------------------------.
| yyerrlab -- here on detecting error.  |
`--------------------------------------*/
yyerrlab:
  /* Make sure we have latest lookahead translation.  See comments at
     user semantic actions for why this is necessary.  */
  yytoken = yychar == FITS_PARSER_YYEMPTY ? YYSYMBOL_YYEMPTY : YYTRANSLATE (yychar);
  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;
      yyerror (scanner, lParse, YY_("syntax error"));
    }

  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse lookahead token after an
         error, discard it.  */

      if (yychar <= FITS_PARSER_YYEOF)
        {
          /* Return failure if at end of input.  */
          if (yychar == FITS_PARSER_YYEOF)
            YYABORT;
        }
      else
        {
          yydestruct ("Error: discarding",
                      yytoken, &yylval, scanner, lParse);
          yychar = FITS_PARSER_YYEMPTY;
        }
    }

  /* Else will try to reuse lookahead token after shifting the error
     token.  */
  goto yyerrlab1;


/*---------------------------------------------------.
| yyerrorlab -- error raised explicitly by YYERROR.  |
`---------------------------------------------------*/
yyerrorlab:
  /* Pacify compilers when the user code never invokes YYERROR and the
     label yyerrorlab therefore never appears in user code.  */
  if (0)
    YYERROR;
  ++yynerrs;

  /* Do not reclaim the symbols of the rule whose action triggered
     this YYERROR.  */
  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);
  yystate = *yyssp;
  goto yyerrlab1;


/*-------------------------------------------------------------.
| yyerrlab1 -- common code for both syntax error and YYERROR.  |
`-------------------------------------------------------------*/
yyerrlab1:
  yyerrstatus = 3;      /* Each real token shifted decrements this.  */

  /* Pop stack until we find a state that shifts the error token.  */
  for (;;)
    {
      yyn = yypact[yystate];
      if (!yypact_value_is_default (yyn))
        {
          yyn += YYSYMBOL_YYerror;
          if (0 <= yyn && yyn <= YYLAST && yycheck[yyn] == YYSYMBOL_YYerror)
            {
              yyn = yytable[yyn];
              if (0 < yyn)
                break;
            }
        }

      /* Pop the current state because it cannot handle the error token.  */
      if (yyssp == yyss)
        YYABORT;


      yydestruct ("Error: popping",
                  YY_ACCESSING_SYMBOL (yystate), yyvsp, scanner, lParse);
      YYPOPSTACK (1);
      yystate = *yyssp;
      YY_STACK_PRINT (yyss, yyssp);
    }

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END


  /* Shift the error token.  */
  YY_SYMBOL_PRINT ("Shifting", YY_ACCESSING_SYMBOL (yyn), yyvsp, yylsp);

  yystate = yyn;
  goto yynewstate;


/*-------------------------------------.
| yyacceptlab -- YYACCEPT comes here.  |
`-------------------------------------*/
yyacceptlab:
  yyresult = 0;
  goto yyreturnlab;


/*-----------------------------------.
| yyabortlab -- YYABORT comes here.  |
`-----------------------------------*/
yyabortlab:
  yyresult = 1;
  goto yyreturnlab;


/*-----------------------------------------------------------.
| yyexhaustedlab -- YYNOMEM (memory exhaustion) comes here.  |
`-----------------------------------------------------------*/
yyexhaustedlab:
  yyerror (scanner, lParse, YY_("memory exhausted"));
  yyresult = 2;
  goto yyreturnlab;


/*----------------------------------------------------------.
| yyreturnlab -- parsing is finished, clean up and return.  |
`----------------------------------------------------------*/
yyreturnlab:
  if (yychar != FITS_PARSER_YYEMPTY)
    {
      /* Make sure we have latest lookahead translation.  See comments at
         user semantic actions for why this is necessary.  */
      yytoken = YYTRANSLATE (yychar);
      yydestruct ("Cleanup: discarding lookahead",
                  yytoken, &yylval, scanner, lParse);
    }
  /* Do not reclaim the symbols of the rule whose action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK (yylen);
  YY_STACK_PRINT (yyss, yyssp);
  while (yyssp != yyss)
    {
      yydestruct ("Cleanup: popping",
                  YY_ACCESSING_SYMBOL (+*yyssp), yyvsp, scanner, lParse);
      YYPOPSTACK (1);
    }
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif

  return yyresult;
}

#line 1336 "eval.y"


/*************************************************************************/
/*  Start of "New" routines which build the expression Nodal structure   */
/*************************************************************************/

static int Alloc_Node( ParseData *lParse )
{
                      /* Use this for allocation to guarantee *Nodes */
   Node *newNodePtr;  /* survives on failure, making it still valid  */
                      /* while working our way out of this error     */

   if( lParse->nNodes == lParse->nNodesAlloc ) {
      if( lParse->Nodes ) {
	 lParse->nNodesAlloc += lParse->nNodesAlloc;
	 newNodePtr = (Node *)realloc( lParse->Nodes,
				       sizeof(Node)*lParse->nNodesAlloc );
      } else {
	 lParse->nNodesAlloc = 100;
	 newNodePtr = (Node *)malloc ( sizeof(Node)*lParse->nNodesAlloc );
      }	 

      if( newNodePtr ) {
	 lParse->Nodes = newNodePtr;
      } else {
	 lParse->status = MEMORY_ALLOCATION;
	 return( -1 );
      }
   }

   return ( lParse->nNodes++ );
}

static void Free_Last_Node( ParseData *lParse )
{
   if( lParse->nNodes ) lParse->nNodes--;
}

static int New_Const( ParseData *lParse, int returnType, void *value, long len )
{
   Node *this;
   int n;

   n = Alloc_Node(lParse);
   if( n>=0 ) {
      this             = lParse->Nodes + n;
      this->operation  = CONST_OP;             /* Flag a constant */
      this->DoOp       = NULL;
      this->nSubNodes  = 0;
      this->type       = returnType;
      memcpy( &(this->value.data), value, len );
      this->value.undef = NULL;
      this->value.nelem = 1;
      this->value.naxis = 1;
      this->value.naxes[0] = 1;
   }
   return(n);
}

static int New_Column( ParseData *lParse, int ColNum )
{
   Node *this;
   int  n, i;

   n = Alloc_Node(lParse);
   if( n>=0 ) {
      this              = lParse->Nodes + n;
      this->operation   = -ColNum;
      this->DoOp        = NULL;
      this->nSubNodes   = 0;
      this->type        = lParse->varData[ColNum].type;
      this->value.nelem = lParse->varData[ColNum].nelem;
      this->value.naxis = lParse->varData[ColNum].naxis;
      for( i=0; i<lParse->varData[ColNum].naxis; i++ )
	 this->value.naxes[i] = lParse->varData[ColNum].naxes[i];
   }
   return(n);
}

static int New_Offset( ParseData *lParse, int ColNum, int offsetNode )
{
   Node *this;
   int  n, i, colNode;

   colNode = New_Column( lParse, ColNum );
   if( colNode<0 ) return(-1);

   n = Alloc_Node(lParse);
   if( n>=0 ) {
      this              = lParse->Nodes + n;
      this->operation   = '{';
      this->DoOp        = Do_Offset;
      this->nSubNodes   = 2;
      this->SubNodes[0] = colNode;
      this->SubNodes[1] = offsetNode;
      this->type        = lParse->varData[ColNum].type;
      this->value.nelem = lParse->varData[ColNum].nelem;
      this->value.naxis = lParse->varData[ColNum].naxis;
      for( i=0; i<lParse->varData[ColNum].naxis; i++ )
	 this->value.naxes[i] = lParse->varData[ColNum].naxes[i];
   }
   return(n);
}

static int New_Unary( ParseData *lParse, int returnType, int Op, int Node1 )
{
   Node *this, *that;
   int  i,n;

   if( Node1<0 ) return(-1);
   that = lParse->Nodes + Node1;

   if( !Op ) Op = returnType;

   if( (Op==DOUBLE || Op==FLTCAST) && that->type==DOUBLE  ) return( Node1 );
   if( (Op==LONG   || Op==INTCAST) && that->type==LONG    ) return( Node1 );
   if( (Op==BOOLEAN              ) && that->type==BOOLEAN ) return( Node1 );
   
   n = Alloc_Node(lParse);
   if( n>=0 ) {
      this              = lParse->Nodes + n;
      this->operation   = Op;
      this->DoOp        = Do_Unary;
      this->nSubNodes   = 1;
      this->SubNodes[0] = Node1;
      this->type        = returnType;

      that              = lParse->Nodes + Node1; /* Reset in case .Nodes mv'd */
      this->value.nelem = that->value.nelem;
      this->value.naxis = that->value.naxis;
      for( i=0; i<that->value.naxis; i++ )
	 this->value.naxes[i] = that->value.naxes[i];

      if( that->operation==CONST_OP ) this->DoOp( lParse, this );
   }
   return( n );
}

static int New_BinOp( ParseData *lParse, int returnType, int Node1, int Op, int Node2 )
{
   Node *this,*that1,*that2;
   int  n,i,constant;

   if( Node1<0 || Node2<0 ) return(-1);

   n = Alloc_Node(lParse);
   if( n>=0 ) {
      this             = lParse->Nodes + n;
      this->operation  = Op;
      this->nSubNodes  = 2;
      this->SubNodes[0]= Node1;
      this->SubNodes[1]= Node2;
      this->type       = returnType;

      that1            = lParse->Nodes + Node1;
      that2            = lParse->Nodes + Node2;
      constant         = (that1->operation==CONST_OP
                          && that2->operation==CONST_OP);
      if( that1->type!=STRING && that1->type!=BITSTR )
	if( !Test_Dims( lParse, Node1, Node2 ) ) {
	    Free_Last_Node(lParse);
	    yyerror(0, lParse, "Array sizes/dims do not match for binary operator");
	    return(-1);
	 }
      if( that1->value.nelem == 1 ) that1 = that2;

      this->value.nelem = that1->value.nelem;
      this->value.naxis = that1->value.naxis;
      for( i=0; i<that1->value.naxis; i++ )
	 this->value.naxes[i] = that1->value.naxes[i];

      if ( Op == ACCUM && that1->type == BITSTR ) {
	/* ACCUM is rank-reducing on bit strings */
	this->value.nelem = 1;
	this->value.naxis = 1;
	this->value.naxes[0] = 1;
      }

      /*  Both subnodes should be of same time  */
      switch( that1->type ) {
      case BITSTR:  this->DoOp = Do_BinOp_bit;  break;
      case STRING:  this->DoOp = Do_BinOp_str;  break;
      case BOOLEAN: this->DoOp = Do_BinOp_log;  break;
      case LONG:    this->DoOp = Do_BinOp_lng;  break;
      case DOUBLE:  this->DoOp = Do_BinOp_dbl;  break;
      }
      if( constant ) this->DoOp( lParse, this );
   }
   return( n );
}

static int New_Func( ParseData *lParse,
		     int returnType, funcOp Op, int nNodes,
		     int Node1, int Node2, int Node3, int Node4, 
		     int Node5, int Node6, int Node7 )
{
  return New_FuncSize(lParse,
		      returnType, Op, nNodes,
		      Node1, Node2, Node3, Node4, 
		      Node5, Node6, Node7, 0);
}

static int New_FuncSize( ParseData *lParse,
			 int returnType, funcOp Op, int nNodes,
			 int Node1, int Node2, int Node3, int Node4, 
			 int Node5, int Node6, int Node7, int Size )
/* If returnType==0 , use Node1's type and vector sizes as returnType, */
/* else return a single value of type returnType                       */
{
   Node *this, *that;
   int  i,n,constant;

   if( Node1<0 || Node2<0 || Node3<0 || Node4<0 || 
       Node5<0 || Node6<0 || Node7<0 ) return(-1);

   n = Alloc_Node(lParse);
   if( n>=0 ) {
      this              = lParse->Nodes + n;
      this->operation   = (int)Op;
      this->DoOp        = Do_Func;
      this->nSubNodes   = nNodes;
      this->SubNodes[0] = Node1;
      this->SubNodes[1] = Node2;
      this->SubNodes[2] = Node3;
      this->SubNodes[3] = Node4;
      this->SubNodes[4] = Node5;
      this->SubNodes[5] = Node6;
      this->SubNodes[6] = Node7;
      i = constant = nNodes;    /* Functions with zero params are not const */
      if (Op == poirnd_fct) constant = 0; /* Nor is Poisson deviate */

      while( i-- )
	constant = ( constant && OPER(this->SubNodes[i]) == CONST_OP );
      
      if( returnType ) {
	 this->type           = returnType;
	 this->value.nelem    = 1;
	 this->value.naxis    = 1;
	 this->value.naxes[0] = 1;
      } else {
	 that              = lParse->Nodes + Node1;
	 this->type        = that->type;
	 this->value.nelem = that->value.nelem;
	 this->value.naxis = that->value.naxis;
	 for( i=0; i<that->value.naxis; i++ )
	    this->value.naxes[i] = that->value.naxes[i];
      }
      /* Force explicit size before evaluating */
      if (Size > 0) this->value.nelem = Size;

      if( constant ) this->DoOp( lParse, this );
   }
   return( n );
}

static int New_Deref( ParseData *lParse, int Var,  int nDim,
		      int Dim1, int Dim2, int Dim3, int Dim4, int Dim5 )
{
   int n, idx, constant;
   long elem=0;
   Node *this, *theVar, *theDim[MAXDIMS];

   if( Var<0 || Dim1<0 || Dim2<0 || Dim3<0 || Dim4<0 || Dim5<0 ) return(-1);

   theVar = lParse->Nodes + Var;
   if( theVar->operation==CONST_OP || theVar->value.nelem==1 ) {
      yyerror(0, lParse, "Cannot index a scalar value");
      return(-1);
   }

   n = Alloc_Node(lParse);
   if( n>=0 ) {
      this              = lParse->Nodes + n;
      this->nSubNodes   = nDim+1;
      theVar            = lParse->Nodes + (this->SubNodes[0]=Var);
      theDim[0]         = lParse->Nodes + (this->SubNodes[1]=Dim1);
      theDim[1]         = lParse->Nodes + (this->SubNodes[2]=Dim2);
      theDim[2]         = lParse->Nodes + (this->SubNodes[3]=Dim3);
      theDim[3]         = lParse->Nodes + (this->SubNodes[4]=Dim4);
      theDim[4]         = lParse->Nodes + (this->SubNodes[5]=Dim5);
      constant          = theVar->operation==CONST_OP;
      for( idx=0; idx<nDim; idx++ )
	 constant = (constant && theDim[idx]->operation==CONST_OP);

      for( idx=0; idx<nDim; idx++ )
	 if( theDim[idx]->value.nelem>1 ) {
	    Free_Last_Node(lParse);
	    yyerror(0, lParse, "Cannot use an array as an index value");
	    return(-1);
	 } else if( theDim[idx]->type!=LONG ) {
	    Free_Last_Node(lParse);
	    yyerror(0, lParse, "Index value must be an integer type");
	    return(-1);
	 }

      this->operation   = '[';
      this->DoOp        = Do_Deref;
      this->type        = theVar->type;

      if( theVar->value.naxis == nDim ) { /* All dimensions specified */
	 this->value.nelem    = 1;
	 this->value.naxis    = 1;
	 this->value.naxes[0] = 1;
      } else if( nDim==1 ) { /* Dereference only one dimension */
	 elem=1;
	 this->value.naxis = theVar->value.naxis-1;
	 for( idx=0; idx<this->value.naxis; idx++ ) {
	    elem *= ( this->value.naxes[idx] = theVar->value.naxes[idx] );
	 }
	 this->value.nelem = elem;
      } else {
	 Free_Last_Node(lParse);
	 yyerror(0, lParse, "Must specify just one or all indices for vector");
	 return(-1);
      }
      if( constant ) this->DoOp( lParse, this );
   }
   return(n);
}

extern int fits_parser_yyGetVariable( ParseData *lParse, char *varName, YYSTYPE *varVal );

static int New_GTI( ParseData *lParse, funcOp Op, char *fname, int Node1, int Node2, char *start, char *stop )
{
   fitsfile *fptr;
   Node *this, *that0, *that1, *that2;
   int  type,i,n, startCol, stopCol, Node0;
   int  hdutype, hdunum, evthdu, samefile, extvers, movetotype, tstat;
   char extname[100];
   long nrows;
   double timeZeroI[2], timeZeroF[2], dt, timeSpan;
   char xcol[20], xexpr[20];
   YYSTYPE colVal;

   if( (Op == gtifilt_fct || Op == gtifind_fct) && Node1==-99 ) {
      type = fits_parser_yyGetVariable( lParse,  "TIME", &colVal );
      if( type==COLUMN ) {
	 Node1 = New_Column( lParse, (int)colVal.lng );
      } else {
	 yyerror(0, lParse, "Could not build TIME column for GTIFILTER/GTIFIND");
	 return(-1);
      }
   }

   if (Op == gtiover_fct) {
     if (Node1 == -99 || Node2 == -99) {
       yyerror(0, lParse, "startExpr and stopExpr values must be defined for GTIOVERLAP");
       return(-1);
     }
     /* Also case TIME_STOP to double precision */
     Node2 = New_Unary( lParse, DOUBLE, 0, Node2 );
     if (Node2 < 0) return(-1);

   }

   /* Type cast TIME to double precision */
   Node1 = New_Unary( lParse, DOUBLE, 0, Node1 );
   Node0 = Alloc_Node(lParse); /* This will hold the START/STOP times */
   if( Node1<0 || Node0<0 ) return(-1);

   /*  Record current HDU number in case we need to move within this file  */

   fptr = lParse->def_fptr;
   ffghdn( fptr, &evthdu );

   /*  Look for TIMEZERO keywords in current extension  */

   tstat = 0;
   if( ffgkyd( fptr, "TIMEZERO", timeZeroI, NULL, &tstat ) ) {
      tstat = 0;
      if( ffgkyd( fptr, "TIMEZERI", timeZeroI, NULL, &tstat ) ) {
	 timeZeroI[0] = timeZeroF[0] = 0.0;
      } else if( ffgkyd( fptr, "TIMEZERF", timeZeroF, NULL, &tstat ) ) {
	 timeZeroF[0] = 0.0;
      }
   } else {
      timeZeroF[0] = 0.0;
   }

   /*  Resolve filename parameter  */

   switch( fname[0] ) {
   case '\0':
      samefile = 1;
      hdunum = 1;
      break;
   case '[':
      samefile = 1;
      i = 1;
      while( fname[i] != '\0' && fname[i] != ']' ) i++;
      if( fname[i] ) {
	 fname[i] = '\0';
	 fname++;
	 ffexts( fname, &hdunum, extname, &extvers, &movetotype,
		 xcol, xexpr, &lParse->status );
         if( *extname ) {
	    ffmnhd( fptr, movetotype, extname, extvers, &lParse->status );
	    ffghdn( fptr, &hdunum );
	 } else if( hdunum ) {
	    ffmahd( fptr, ++hdunum, &hdutype, &lParse->status );
	 } else if( !lParse->status ) {
	    yyerror(0, lParse, "Cannot use primary array for GTI filter");
	    return( -1 );
	 }
      } else {
	 yyerror(0, lParse, "File extension specifier lacks closing ']'");
	 return( -1 );
      }
      break;
   case '+':
      samefile = 1;
      hdunum = atoi( fname ) + 1;
      if( hdunum>1 )
	 ffmahd( fptr, hdunum, &hdutype, &lParse->status );
      else {
	 yyerror(0, lParse, "Cannot use primary array for GTI filter / GTIFIND");
	 return( -1 );
      }
      break;
   default:
      samefile = 0;
      if( ! ffopen( &fptr, fname, READONLY, &lParse->status ) )
	 ffghdn( fptr, &hdunum );
      break;
   }
   if( lParse->status ) return(-1);

   /*  If at primary, search for GTI extension  */

   if( hdunum==1 ) {
      while( 1 ) {
	 hdunum++;
	 if( ffmahd( fptr, hdunum, &hdutype, &lParse->status ) ) break;
	 if( hdutype==IMAGE_HDU ) continue;
	 tstat = 0;
	 if( ffgkys( fptr, "EXTNAME", extname, NULL, &tstat ) ) continue;
	 ffupch( extname );
	 if( strstr( extname, "GTI" ) ) break;
      }
      if( lParse->status ) {
	 if( lParse->status==END_OF_FILE )
	    yyerror(0, lParse, "GTI extension not found in this file");
	 return(-1);
      }
   }

   /*  Locate START/STOP Columns  */

   ffgcno( fptr, CASEINSEN, start, &startCol, &lParse->status );
   ffgcno( fptr, CASEINSEN, stop,  &stopCol,  &lParse->status );
   if( lParse->status ) return(-1);

   /*  Look for TIMEZERO keywords in GTI extension  */

   tstat = 0;
   if( ffgkyd( fptr, "TIMEZERO", timeZeroI+1, NULL, &tstat ) ) {
      tstat = 0;
      if( ffgkyd( fptr, "TIMEZERI", timeZeroI+1, NULL, &tstat ) ) {
	 timeZeroI[1] = timeZeroF[1] = 0.0;
      } else if( ffgkyd( fptr, "TIMEZERF", timeZeroF+1, NULL, &tstat ) ) {
	 timeZeroF[1] = 0.0;
      }
   } else {
      timeZeroF[1] = 0.0;
   }

   n = Alloc_Node(lParse);
   if( n >= 0 ) {
      this                 = lParse->Nodes + n;
      this->SubNodes[1]    = Node1;
      this->operation      = (int) Op;
      if (Op == gtifilt_fct) {
	this->nSubNodes      = 2;
	this->DoOp           = Do_GTI;
	this->type           = BOOLEAN;
      } else if (Op == gtifind_fct) {
	this->nSubNodes      = 2;
	this->DoOp           = Do_GTI;
	this->type           = LONG;
      } else {
	this->nSubNodes      = 3;
	this->DoOp           = Do_GTI_Over;
	this->type           = DOUBLE;
      }
      that1                = lParse->Nodes + Node1;
      this->value.nelem    = that1->value.nelem;
      this->value.naxis    = that1->value.naxis;
      for( i=0; i < that1->value.naxis; i++ )
	 this->value.naxes[i] = that1->value.naxes[i];
      if (Op == gtiover_fct) {
	this->SubNodes[2]  = Node2;
	that2 = lParse->Nodes + Node2;
	if (that1->value.nelem != that2->value.nelem) {
	  yyerror(0, lParse, "Dimensions of TIME and TIME_STOP must match for GTIOVERLAP");
	  return(-1);
	}
      }

      /* Init START/STOP node to be treated as a "constant" */

      this->SubNodes[0]    = Node0;
      that0                = lParse->Nodes + Node0;
      that0->operation     = CONST_OP;
      that0->DoOp          = NULL;
      that0->value.data.ptr= NULL;

      /*  Read in START/STOP times  */

      if( ffgkyj( fptr, "NAXIS2", &nrows, NULL, &lParse->status ) )
	 return(-1);
      that0->value.nelem = nrows;
      if( nrows ) {

	 that0->value.data.dblptr = (double*)malloc( 2*nrows*sizeof(double) );
	 if( !that0->value.data.dblptr ) {
	    lParse->status = MEMORY_ALLOCATION;
	    return(-1);
	 }
	 
	 ffgcvd( fptr, startCol, 1L, 1L, nrows, 0.0,
		 that0->value.data.dblptr, &i, &lParse->status );
	 ffgcvd( fptr, stopCol, 1L, 1L, nrows, 0.0,
		 that0->value.data.dblptr+nrows, &i, &lParse->status );
	 if( lParse->status ) {
	    free( that0->value.data.dblptr );
	    return(-1);
	 }

	 /*  Test for fully time-ordered GTI... both START && STOP  */

	 that0->type = 1; /*  Assume yes  */
	 i = nrows;
	 while( --i )
	    if(    that0->value.data.dblptr[i-1]
                   >= that0->value.data.dblptr[i]
		|| that0->value.data.dblptr[i-1+nrows]
		   >= that0->value.data.dblptr[i+nrows] ) {
	       that0->type = 0;
	       break;
	    }

	 /* GTIOVERLAP() requires ordered GTI */
	 if (that0->type != 1 && Op == gtiover_fct) {
	   yyerror(0, lParse, "Input GTI must be time-ordered for GTIOVERLAP");
	   return(-1);
	 }
	 
	 /*  Handle TIMEZERO offset, if any  */
	 
	 dt = (timeZeroI[1] - timeZeroI[0]) + (timeZeroF[1] - timeZeroF[0]);
	 timeSpan = that0->value.data.dblptr[nrows+nrows-1]
	    - that0->value.data.dblptr[0];
	 if (timeSpan == 0) timeSpan = 1.0;
	 
	 if( fabs( dt / timeSpan ) > 1e-12 ) {
	    for( i=0; i<(nrows+nrows); i++ )
	       that0->value.data.dblptr[i] += dt;
	 }
      }
      /* If Node1 is constant (gtifilt_fct) or
	 Node1 and Node2 are constant (gtiover_fct), then evaluate now */
      if( OPER(Node1)==CONST_OP && (Op == gtifilt_fct || OPER(Node2)==CONST_OP)) {
	this->DoOp( lParse, this );
      }
   }

   if( samefile )
      ffmahd( fptr, evthdu, &hdutype, &lParse->status );
   else
      ffclos( fptr, &lParse->status );

   return( n );
}

static int New_REG( ParseData *lParse, char *fname, int NodeX, int NodeY, char *colNames )
{
   Node *this, *that0;
   int  type, n, Node0;
   int  Xcol, Ycol, tstat;
   WCSdata wcs;
   SAORegion *Rgn;
   char *cX, *cY;
   YYSTYPE colVal;

   if( NodeX==-99 ) {
      type = fits_parser_yyGetVariable( lParse,  "X", &colVal );
      if( type==COLUMN ) {
	 NodeX = New_Column( lParse, (int)colVal.lng );
      } else {
	 yyerror(0, lParse, "Could not build X column for REGFILTER");
	 return(-1);
      }
   }
   if( NodeY==-99 ) {
      type = fits_parser_yyGetVariable( lParse, "Y", &colVal );
      if( type==COLUMN ) {
 	 NodeY = New_Column( lParse, (int)colVal.lng );
      } else {
	 yyerror(0, lParse, "Could not build Y column for REGFILTER");
	 return(-1);
      }
   }
   NodeX = New_Unary( lParse, DOUBLE, 0, NodeX );
   NodeY = New_Unary( lParse, DOUBLE, 0, NodeY );
   Node0 = Alloc_Node(lParse); /* This will hold the Region Data */
   if( NodeX<0 || NodeY<0 || Node0<0 ) return(-1);

   if( ! (Test_Dims( lParse, NodeX, NodeY ) ) ) {
     yyerror(0, lParse, "Dimensions of REGFILTER arguments are not compatible");
     return (-1);
   }

   n = Alloc_Node(lParse);
   if( n >= 0 ) {
      this                 = lParse->Nodes + n;
      this->nSubNodes      = 3;
      this->SubNodes[0]    = Node0;
      this->SubNodes[1]    = NodeX;
      this->SubNodes[2]    = NodeY;
      this->operation      = (int)regfilt_fct;
      this->DoOp           = Do_REG;
      this->type           = BOOLEAN;
      this->value.nelem    = 1;
      this->value.naxis    = 1;
      this->value.naxes[0] = 1;
      
      Copy_Dims(lParse, n, NodeX);
      if( SIZE(NodeX)<SIZE(NodeY) )  Copy_Dims(lParse, n, NodeY);

      /* Init Region node to be treated as a "constant" */

      that0                = lParse->Nodes + Node0;
      that0->operation     = CONST_OP;
      that0->DoOp          = NULL;

      /*  Identify what columns to use for WCS information  */

      Xcol = Ycol = 0;
      if( *colNames ) {
	 /*  Use the column names in this string for WCS info  */
	 while( *colNames==' ' ) colNames++;
	 cX = cY = colNames;
	 while( *cY && *cY!=' ' && *cY!=',' ) cY++;
	 if( *cY )
	    *(cY++) = '\0';
	 while( *cY==' ' ) cY++;
	 if( !*cY ) {
	    yyerror(0, lParse, "Could not extract valid pair of column names from REGFILTER");
	    Free_Last_Node(lParse);
	    return( -1 );
	 }
	 fits_get_colnum( lParse->def_fptr, CASEINSEN, cX, &Xcol,
			  &lParse->status );
	 fits_get_colnum( lParse->def_fptr, CASEINSEN, cY, &Ycol,
			  &lParse->status );
	 if( lParse->status ) {
	    yyerror(0, lParse, "Could not locate columns indicated for WCS info");
	    Free_Last_Node(lParse);
	    return( -1 );
	 }

      } else {
	 /*  Try to find columns used in X/Y expressions  */
	 Xcol = Locate_Col( lParse, lParse->Nodes + NodeX );
	 Ycol = Locate_Col( lParse, lParse->Nodes + NodeY );
	 if( Xcol<0 || Ycol<0 ) {
	    yyerror(0, lParse, "Found multiple X/Y column references in REGFILTER");
	    Free_Last_Node(lParse);
	    return( -1 );
	 }
      }

      /*  Now, get the WCS info, if it exists, from the indicated columns  */
      wcs.exists = 0;
      if( Xcol>0 && Ycol>0 ) {
	 tstat = 0;
	 ffgtcs( lParse->def_fptr, Xcol, Ycol,
		 &wcs.xrefval, &wcs.yrefval,
		 &wcs.xrefpix, &wcs.yrefpix,
		 &wcs.xinc,    &wcs.yinc,
		 &wcs.rot,      wcs.type,
		 &tstat );
	 if( tstat==NO_WCS_KEY ) {
	    wcs.exists = 0;
	 } else if( tstat ) {
	    lParse->status = tstat;
	    Free_Last_Node(lParse);
	    return( -1 );
	 } else {
	    wcs.exists = 1;
	 }
      }

      /*  Read in Region file  */

      fits_read_rgnfile( fname, &wcs, &Rgn, &lParse->status );
      if( lParse->status ) {
	 Free_Last_Node(lParse);
	 return( -1 );
      }

      that0->value.data.ptr = Rgn;

      if( OPER(NodeX)==CONST_OP && OPER(NodeY)==CONST_OP )
	 this->DoOp( lParse, this );
   }

   return( n );
}

static int New_Vector( ParseData *lParse, int subNode )
{
   Node *this, *that;
   int n;

   n = Alloc_Node(lParse);
   if( n >= 0 ) {
      this              = lParse->Nodes + n;
      that              = lParse->Nodes + subNode;
      this->type        = that->type;
      this->nSubNodes   = 1;
      this->SubNodes[0] = subNode;
      this->operation   = '{';
      this->DoOp        = Do_Vector;
   }

   return( n );
}

static int Close_Vec( ParseData *lParse, int vecNode )
{
   Node *this;
   int n, nelem=0;

   this = lParse->Nodes + vecNode;
   for( n=0; n < this->nSubNodes; n++ ) {
      if( TYPE( this->SubNodes[n] ) != this->type ) {
	 this->SubNodes[n] = New_Unary( lParse, this->type, 0, this->SubNodes[n] );
	 if( this->SubNodes[n]<0 ) return(-1);
      }
      nelem += SIZE(this->SubNodes[n]);
   }
   this->value.naxis    = 1;
   this->value.nelem    = nelem;
   this->value.naxes[0] = nelem;

   return( vecNode );
}

static int New_Array( ParseData *lParse, int valueNode, int dimNode )
{
  Node *dims;
  long naxis, nelem;
  long naxes[MAXDIMS];
  Node *this;
  int  n,i;

   if( valueNode<0 || dimNode<0 ) return(-1);

   /* Check that dimensions are {a,b,c,d}
        - vector
	- every element is constant integer
	- 5 or fewer dimensions 
   */

   if (SIZE(valueNode) > 1) {
     yyerror(0, lParse, "ARRAY(V,n) value V must have vector dimension of 1");
     return (-1);
   }

   dims = &(lParse->Nodes[dimNode]);
   for (i=0; i<MAXDIMS; i++) naxes[i] = 1;

   if (OPER(dimNode) == CONST_OP) { /* ARRAY(V,n) is a constant integer */
     if ( TYPE(dimNode) != LONG ) dimNode = New_Unary(lParse, LONG, 0, dimNode);
     if (dimNode < 0) return (-1);
     naxis = 1;
     naxes[0] = lParse->Nodes[dimNode].value.data.lng;

   } else if (OPER(dimNode) == '{') { /* ARRAY(V,{a,b,c,d,e}) up to 5 dimensions */
     if (dims->nSubNodes > MAXDIMS) {
       yyerror(0, lParse, "ARRAY(V,{...}) number of dimensions must not exceed 5");
       return (-1);
     }
     naxis = dims->nSubNodes;
     for (i=0; i<dims->nSubNodes; i++) {
       if ( TYPE(dims->SubNodes[i]) != LONG ) {
	 dims->SubNodes[i] = New_Unary(lParse, LONG, 0, dims->SubNodes[i]);
	 if (dims->SubNodes[i] < 0) return (-1);
       }
       naxes[i] = lParse->Nodes[ dims->SubNodes[i] ].value.data.lng;
     }
   } else {
     yyerror(0, lParse, "ARRAY(V,dims) dims must be either integer or const vector");
     return (-1);
   }

   nelem = 1;
   for (i=0; i<naxis; i++) {
     if (naxes[i] <= 0) {
       yyerror(0, lParse, "ARRAY(V,dims) must have positive dimensions");
       return (-1);
     }
     nelem *= naxes[i];
   }

   n = Alloc_Node(lParse);
   if( n>=0 ) {
      this             = lParse->Nodes + n;
      this->operation  = array_fct;
      this->nSubNodes  = 1;
      this->SubNodes[0]= valueNode;
      this->type       = TYPE(valueNode);

      this->value.nelem = nelem;
      this->value.naxis = naxis;
      for( i=0; i<naxis; i++ )
	this->value.naxes[i] = naxes[i];

      this->DoOp = Do_Array;
   }
   return( n );
}

static int Locate_Col( ParseData *lParse, Node *this )
/*  Locate the TABLE column number of any columns in "this" calculation.  */
/*  Return ZERO if none found, or negative if more than 1 found.          */
{
   Node *that;
   int  i, col=0, newCol, nfound=0;
   
   if( this->nSubNodes==0
       && this->operation<=0 && this->operation!=CONST_OP )
      return lParse->colData[ - this->operation].colnum;

   for( i=0; i<this->nSubNodes; i++ ) {
      that = lParse->Nodes + this->SubNodes[i];
      if( that->operation>0 ) {
	 newCol = Locate_Col( lParse, that );
	 if( newCol<=0 ) {
	    nfound += -newCol;
	 } else {
	    if( !nfound ) {
	       col = newCol;
	       nfound++;
	    } else if( col != newCol ) {
	       nfound++;
	    }
	 }
      } else if( that->operation!=CONST_OP ) {
	 /*  Found a Column  */
	 newCol = lParse->colData[- that->operation].colnum;
	 if( !nfound ) {
	    col = newCol;
	    nfound++;
	 } else if( col != newCol ) {
	    nfound++;
	 }
      }
   }
   if( nfound!=1 )
      return( - nfound );
   else
      return( col );
}

static int Test_Dims( ParseData *lParse, int Node1, int Node2 )
{
   Node *that1, *that2;
   int valid, i;

   if( Node1<0 || Node2<0 ) return(0);

   that1 = lParse->Nodes + Node1;
   that2 = lParse->Nodes + Node2;

   if( that1->value.nelem==1 || that2->value.nelem==1 )
      valid = 1;
   else if( that1->type==that2->type
	    && that1->value.nelem==that2->value.nelem
	    && that1->value.naxis==that2->value.naxis ) {
      valid = 1;
      for( i=0; i<that1->value.naxis; i++ ) {
	 if( that1->value.naxes[i]!=that2->value.naxes[i] )
	    valid = 0;
      }
   } else
      valid = 0;
   return( valid );
}   

static void Copy_Dims( ParseData *lParse, int Node1, int Node2 )
{
   Node *that1, *that2;
   int i;

   if( Node1<0 || Node2<0 ) return;

   that1 = lParse->Nodes + Node1;
   that2 = lParse->Nodes + Node2;

   that1->value.nelem = that2->value.nelem;
   that1->value.naxis = that2->value.naxis;
   for( i=0; i<that2->value.naxis; i++ )
      that1->value.naxes[i] = that2->value.naxes[i];
}

/********************************************************************/
/*    Routines for actually evaluating the expression start here    */
/********************************************************************/

void Evaluate_Parser( ParseData *lParse, long firstRow, long nRows )
    /***********************************************************************/
    /*  Reset the parser for processing another batch of data...           */
    /*    firstRow:  Row number of the first element to evaluate           */
    /*    nRows:     Number of rows to be processed                        */
    /*  Initialize each COLUMN node so that its UNDEF and DATA pointers    */
    /*  point to the appropriate column arrays.                            */
    /*  Finally, call Evaluate_Node for final node.                        */
    /***********************************************************************/
{
   int     i, column;
   long    offset, rowOffset;
   static int rand_initialized = 0;

   /* Initialize the random number generator once and only once */
   if (rand_initialized == 0) {
     simplerng_srand( (unsigned int) time(NULL) );
     rand_initialized = 1;
   }

   lParse->firstRow = firstRow;
   lParse->nRows    = nRows;

   /*  Reset Column Nodes' pointers to point to right data and UNDEF arrays  */

   rowOffset = firstRow - lParse->firstDataRow;
   for( i=0; i<lParse->nNodes; i++ ) {
     if(    OPER(i) >  0 || OPER(i) == CONST_OP ) continue;

      column = -OPER(i);
      offset = lParse->varData[column].nelem * rowOffset;

      lParse->Nodes[i].value.undef = lParse->varData[column].undef + offset;

      switch( lParse->Nodes[i].type ) {
      case BITSTR:
	 lParse->Nodes[i].value.data.strptr =
	    (char**)lParse->varData[column].data + rowOffset;
	 lParse->Nodes[i].value.undef       = NULL;
	 break;
      case STRING:
	 lParse->Nodes[i].value.data.strptr = 
	    (char**)lParse->varData[column].data + rowOffset;
	 lParse->Nodes[i].value.undef = lParse->varData[column].undef + rowOffset;
	 break;
      case BOOLEAN:
	 lParse->Nodes[i].value.data.logptr = 
	    (char*)lParse->varData[column].data + offset;
	 break;
      case LONG:
	 lParse->Nodes[i].value.data.lngptr = 
	    (long*)lParse->varData[column].data + offset;
	 break;
      case DOUBLE:
	 lParse->Nodes[i].value.data.dblptr = 
	    (double*)lParse->varData[column].data + offset;
	 break;
      }
   }

   Evaluate_Node( lParse, lParse->resultNode );
}

static void Evaluate_Node( ParseData *lParse, int thisNode )
    /**********************************************************************/
    /*  Recursively evaluate thisNode's subNodes, then call one of the    */
    /*  Do_<Action> functions pointed to by thisNode's DoOp element.      */
    /**********************************************************************/
{
   Node *this;
   int i;
   
   if( lParse->status ) return;

   this = lParse->Nodes + thisNode;
   if( this->operation>0 ) {  /* <=0 indicate constants and columns */
      i = this->nSubNodes;
      while( i-- ) {
	 Evaluate_Node( lParse, this->SubNodes[i] );
	 if( lParse->status ) return;
      }
      this->DoOp( lParse, this );
   }
}

static void Allocate_Ptrs( ParseData *lParse, Node *this )
{
   long elem, row, size;

   if( this->type==BITSTR || this->type==STRING ) {

      this->value.data.strptr = (char**)malloc( lParse->nRows
						* sizeof(char*) );
      if( this->value.data.strptr ) {
	 this->value.data.strptr[0] = (char*)malloc( lParse->nRows
						     * (this->value.nelem+2)
						     * sizeof(char) );
	 if( this->value.data.strptr[0] ) {
	    row = 0;
	    while( (++row)<lParse->nRows ) {
	       this->value.data.strptr[row] =
		  this->value.data.strptr[row-1] + this->value.nelem+1;
	    }
	    if( this->type==STRING ) {
	       this->value.undef = this->value.data.strptr[row-1]
                                   + this->value.nelem+1;
	    } else {
	       this->value.undef = NULL;  /* BITSTRs don't use undef array */
	    }
	 } else {
	    lParse->status = MEMORY_ALLOCATION;
	    free( this->value.data.strptr );
	 }
      } else {
	 lParse->status = MEMORY_ALLOCATION;
      }

   } else {

      elem = this->value.nelem * lParse->nRows;
      switch( this->type ) {
      case DOUBLE:  size = sizeof( double ); break;
      case LONG:    size = sizeof( long   ); break;
      case BOOLEAN: size = sizeof( char   ); break;
      default:      size = 1;                break;
      }

      this->value.data.ptr = calloc(size+1, elem);

      if( this->value.data.ptr==NULL ) {
	 lParse->status = MEMORY_ALLOCATION;
      } else {
	 this->value.undef = (char *)this->value.data.ptr + elem*size;
      }
   }
}

static void Do_Unary( ParseData *lParse, Node *this )
{
   Node *that;
   long elem;

   that = lParse->Nodes + this->SubNodes[0];

   if( that->operation==CONST_OP ) {  /* Operating on a constant! */
      switch( this->operation ) {
      case DOUBLE:
      case FLTCAST:
	 if( that->type==LONG )
	    this->value.data.dbl = (double)that->value.data.lng;
	 else if( that->type==BOOLEAN )
	    this->value.data.dbl = ( that->value.data.log ? 1.0 : 0.0 );
	 break;
      case LONG:
      case INTCAST:
	 if( that->type==DOUBLE )
	    this->value.data.lng = (long)that->value.data.dbl;
	 else if( that->type==BOOLEAN )
	    this->value.data.lng = ( that->value.data.log ? 1L : 0L );
	 break;
      case BOOLEAN:
	 if( that->type==DOUBLE )
	    this->value.data.log = ( that->value.data.dbl != 0.0 );
	 else if( that->type==LONG )
	    this->value.data.log = ( that->value.data.lng != 0L );
	 break;
      case UMINUS:
	 if( that->type==DOUBLE )
	    this->value.data.dbl = - that->value.data.dbl;
	 else if( that->type==LONG )
	    this->value.data.lng = - that->value.data.lng;
	 break;
      case NOT:
	 if( that->type==BOOLEAN )
	    this->value.data.log = ( ! that->value.data.log );
	 else if( that->type==BITSTR )
	    bitnot( this->value.data.str, that->value.data.str );
	 break;
      }
      this->operation = CONST_OP;

   } else {

      Allocate_Ptrs( lParse, this );

      if( !lParse->status ) {

	 if( this->type!=BITSTR ) {
	    elem = lParse->nRows;
	    if( this->type!=STRING )
	       elem *= this->value.nelem;
	    while( elem-- )
	       this->value.undef[elem] = that->value.undef[elem];
	 }

	 elem = lParse->nRows * this->value.nelem;

	 switch( this->operation ) {

	 case BOOLEAN:
	    if( that->type==DOUBLE )
	       while( elem-- )
		  this->value.data.logptr[elem] =
		     ( that->value.data.dblptr[elem] != 0.0 );
	    else if( that->type==LONG )
	       while( elem-- )
		  this->value.data.logptr[elem] =
		     ( that->value.data.lngptr[elem] != 0L );
	    break;

	 case DOUBLE:
	 case FLTCAST:
	    if( that->type==LONG )
	       while( elem-- )
		  this->value.data.dblptr[elem] =
		     (double)that->value.data.lngptr[elem];
	    else if( that->type==BOOLEAN )
	       while( elem-- )
		  this->value.data.dblptr[elem] =
		     ( that->value.data.logptr[elem] ? 1.0 : 0.0 );
	    break;

	 case LONG:
	 case INTCAST:
	    if( that->type==DOUBLE )
	       while( elem-- )
		  this->value.data.lngptr[elem] =
		     (long)that->value.data.dblptr[elem];
	    else if( that->type==BOOLEAN )
	       while( elem-- )
		  this->value.data.lngptr[elem] =
		     ( that->value.data.logptr[elem] ? 1L : 0L );
	    break;

	 case UMINUS:
	    if( that->type==DOUBLE ) {
	       while( elem-- )
		  this->value.data.dblptr[elem] =
		     - that->value.data.dblptr[elem];
	    } else if( that->type==LONG ) {
	       while( elem-- )
		  this->value.data.lngptr[elem] =
		     - that->value.data.lngptr[elem];
	    }
	    break;

	 case NOT:
	    if( that->type==BOOLEAN ) {
	       while( elem-- )
		  this->value.data.logptr[elem] =
		     ( ! that->value.data.logptr[elem] );
	    } else if( that->type==BITSTR ) {
	       elem = lParse->nRows;
	       while( elem-- )
		  bitnot( this->value.data.strptr[elem],
			  that->value.data.strptr[elem] );
	    }
	    break;
	 }
      }
   }

   if( that->operation>0 ) {
      free( that->value.data.ptr );
   }
}

static void Do_Offset( ParseData *lParse, Node *this )
{
   Node *col;
   long fRow, nRowOverlap, nRowReload, rowOffset;
   long nelem, elem, offset, nRealElem;
   int status;

   col       = lParse->Nodes + this->SubNodes[0];
   rowOffset = lParse->Nodes[  this->SubNodes[1] ].value.data.lng;

   Allocate_Ptrs( lParse, this );

   fRow   = lParse->firstRow + rowOffset;
   if( this->type==STRING || this->type==BITSTR )
      nRealElem = 1;
   else
      nRealElem = this->value.nelem;

   nelem = nRealElem;

   if( fRow < lParse->firstDataRow ) {

      /* Must fill in data at start of array */

      nRowReload = lParse->firstDataRow - fRow;
      if( nRowReload > lParse->nRows ) nRowReload = lParse->nRows;
      nRowOverlap = lParse->nRows - nRowReload;

      offset = 0;

      /*  NULLify any values falling out of bounds  */

      while( fRow<1 && nRowReload>0 ) {
	 if( this->type == BITSTR ) {
	    nelem = this->value.nelem;
	    this->value.data.strptr[offset][ nelem ] = '\0';
	    while( nelem-- ) this->value.data.strptr[offset][nelem] = '0';
	    offset++;
	 } else {
	    while( nelem-- )
	       this->value.undef[offset++] = 1;
	 }
	 nelem = nRealElem;
	 fRow++;
	 nRowReload--;
      }

   } else if( fRow + lParse->nRows > lParse->firstDataRow + lParse->nDataRows ) {

      /* Must fill in data at end of array */

      nRowReload = (fRow+lParse->nRows) - (lParse->firstDataRow+lParse->nDataRows);
      if( nRowReload>lParse->nRows ) {
	 nRowReload = lParse->nRows;
      } else {
	 fRow = lParse->firstDataRow + lParse->nDataRows;
      }
      nRowOverlap = lParse->nRows - nRowReload;

      offset = nRowOverlap * nelem;

      /*  NULLify any values falling out of bounds  */

      elem = lParse->nRows * nelem;
      while( fRow+nRowReload>lParse->totalRows && nRowReload>0 ) {
	 if( this->type == BITSTR ) {
	    nelem = this->value.nelem;
	    elem--;
	    this->value.data.strptr[elem][ nelem ] = '\0';
	    while( nelem-- ) this->value.data.strptr[elem][nelem] = '0';
	 } else {
	    while( nelem-- )
	       this->value.undef[--elem] = 1;
	 }
	 nelem = nRealElem;
	 nRowReload--;
      }

   } else {

      nRowReload  = 0;
      nRowOverlap = lParse->nRows;
      offset      = 0;

   }

   if( nRowReload>0 ) {
      switch( this->type ) {
      case BITSTR:
      case STRING:
	 status = (*lParse->loadData)( lParse, -col->operation, fRow, nRowReload,
				      this->value.data.strptr+offset,
				      this->value.undef+offset );
	 break;
      case BOOLEAN:
	 status = (*lParse->loadData)( lParse, -col->operation, fRow, nRowReload,
				      this->value.data.logptr+offset,
				      this->value.undef+offset );
	 break;
      case LONG:
	 status = (*lParse->loadData)( lParse, -col->operation, fRow, nRowReload,
				      this->value.data.lngptr+offset,
				      this->value.undef+offset );
	 break;
      case DOUBLE:
	 status = (*lParse->loadData)( lParse, -col->operation, fRow, nRowReload,
				      this->value.data.dblptr+offset,
				      this->value.undef+offset );
	 break;
      }
   }

   /*  Now copy over the overlapping region, if any  */

   if( nRowOverlap <= 0 ) return;

   if( rowOffset>0 )
      elem = nRowOverlap * nelem;
   else
      elem = lParse->nRows * nelem;

   offset = nelem * rowOffset;
   while( nRowOverlap-- && !lParse->status ) {
      while( nelem-- && !lParse->status ) {
	 elem--;
	 if( this->type != BITSTR )
	    this->value.undef[elem] = col->value.undef[elem+offset];
	 switch( this->type ) {
	 case BITSTR:
	    strcpy( this->value.data.strptr[elem       ],
                     col->value.data.strptr[elem+offset] );
	    break;
	 case STRING:
	    strcpy( this->value.data.strptr[elem       ],
                     col->value.data.strptr[elem+offset] );
	    break;
	 case BOOLEAN:
	    this->value.data.logptr[elem] = col->value.data.logptr[elem+offset];
	    break;
	 case LONG:
	    this->value.data.lngptr[elem] = col->value.data.lngptr[elem+offset];
	    break;
	 case DOUBLE:
	    this->value.data.dblptr[elem] = col->value.data.dblptr[elem+offset];
	    break;
	 }
      }
      nelem = nRealElem;
   }
}

static void Do_BinOp_bit( ParseData *lParse, Node *this )
{
   Node *that1, *that2;
   char *sptr1=NULL, *sptr2=NULL;
   int  const1, const2;
   long rows;

   that1 = lParse->Nodes + this->SubNodes[0];
   that2 = lParse->Nodes + this->SubNodes[1];

   const1 = ( that1->operation==CONST_OP );
   const2 = ( that2->operation==CONST_OP );
   sptr1  = ( const1 ? that1->value.data.str : NULL );
   sptr2  = ( const2 ? that2->value.data.str : NULL );

   if( const1 && const2 ) {
      switch( this->operation ) {
      case NE:
	 this->value.data.log = !bitcmp( sptr1, sptr2 );
	 break;
      case EQ:
	 this->value.data.log =  bitcmp( sptr1, sptr2 );
	 break;
      case GT:
      case LT:
      case LTE:
      case GTE:
	 this->value.data.log = bitlgte( sptr1, this->operation, sptr2 );
	 break;
      case '|': 
	 bitor( this->value.data.str, sptr1, sptr2 );
	 break;
      case '&': 
	 bitand( this->value.data.str, sptr1, sptr2 );
	 break;
      case '+':
	 strcpy( this->value.data.str, sptr1 );
	 strcat( this->value.data.str, sptr2 );
	 break;
      case ACCUM:
	this->value.data.lng = 0;
	while( *sptr1 ) {
	  if ( *sptr1 == '1' ) this->value.data.lng ++;
	  sptr1 ++;
	}
	break;
	
      }
      this->operation = CONST_OP;

   } else {

      Allocate_Ptrs( lParse, this );
     
      if( !lParse->status ) {
	 rows  = lParse->nRows;
	 switch( this->operation ) {

	    /*  BITSTR comparisons  */

	 case NE:
	 case EQ:
	 case GT:
	 case LT:
	 case LTE:
	 case GTE:
	    while( rows-- ) {
	       if( !const1 )
		  sptr1 = that1->value.data.strptr[rows];
	       if( !const2 )
		  sptr2 = that2->value.data.strptr[rows];
	       switch( this->operation ) {
	       case NE:  this->value.data.logptr[rows] = 
                                                      !bitcmp( sptr1, sptr2 );
                         break;
	       case EQ:  this->value.data.logptr[rows] = 
                                                       bitcmp( sptr1, sptr2 );
                         break;
	       case GT:
	       case LT:
	       case LTE:
	       case GTE: this->value.data.logptr[rows] = 
                                     bitlgte( sptr1, this->operation, sptr2 );
	                 break;
	       }
	       this->value.undef[rows] = 0;
	    }
	    break;
	 
	    /*  BITSTR AND/ORs ...  no UNDEFS in or out */
      
	 case '|': 
	 case '&': 
	 case '+':
	    while( rows-- ) {
	       if( !const1 )
		  sptr1 = that1->value.data.strptr[rows];
	       if( !const2 )
		  sptr2 = that2->value.data.strptr[rows];
	       if( this->operation=='|' )
		  bitor(  this->value.data.strptr[rows], sptr1, sptr2 );
	       else if( this->operation=='&' )
		  bitand( this->value.data.strptr[rows], sptr1, sptr2 );
	       else {
		  strcpy( this->value.data.strptr[rows], sptr1 );
		  strcat( this->value.data.strptr[rows], sptr2 );
	       }
	    }
	    break;

	    /* Accumulate 1 bits */
	 case ACCUM:
	   { 
	     long i, previous, curr;

	     previous = that2->value.data.lng;
	     
	      /* Cumulative sum of this chunk */
	     for (i=0; i<rows; i++) {
	       sptr1 = that1->value.data.strptr[i];
	       for (curr = 0; *sptr1; sptr1 ++) {
		 if ( *sptr1 == '1' ) curr ++;
	       }
	       previous += curr;
	       this->value.data.lngptr[i] = previous;
	       this->value.undef[i] = 0;
	     }
	     
	      /* Store final cumulant for next pass */
	     that2->value.data.lng = previous;
	   }
	 }
      }
   }

   if( that1->operation>0 ) {
      free( that1->value.data.strptr[0] );
      free( that1->value.data.strptr    );
   }
   if( that2->operation>0 ) {
      free( that2->value.data.strptr[0] );
      free( that2->value.data.strptr    );
   }
}

static void Do_BinOp_str( ParseData *lParse, Node *this )
{
   Node *that1, *that2;
   char *sptr1, *sptr2, null1=0, null2=0;
   int const1, const2, val;
   long rows;

   that1 = lParse->Nodes + this->SubNodes[0];
   that2 = lParse->Nodes + this->SubNodes[1];

   const1 = ( that1->operation==CONST_OP );
   const2 = ( that2->operation==CONST_OP );
   sptr1  = ( const1 ? that1->value.data.str : NULL );
   sptr2  = ( const2 ? that2->value.data.str : NULL );

   if( const1 && const2 ) {  /*  Result is a constant  */
      switch( this->operation ) {

	 /*  Compare Strings  */

      case NE:
      case EQ:
	 val = ( FSTRCMP( sptr1, sptr2 ) == 0 );
	 this->value.data.log = ( this->operation==EQ ? val : !val );
	 break;
      case GT:
	 this->value.data.log = ( FSTRCMP( sptr1, sptr2 ) > 0 );
	 break;
      case LT:
	 this->value.data.log = ( FSTRCMP( sptr1, sptr2 ) < 0 );
	 break;
      case GTE:
	 this->value.data.log = ( FSTRCMP( sptr1, sptr2 ) >= 0 );
	 break;
      case LTE:
	 this->value.data.log = ( FSTRCMP( sptr1, sptr2 ) <= 0 );
	 break;

	 /*  Concat Strings  */

      case '+':
	 strcpy( this->value.data.str, sptr1 );
	 strcat( this->value.data.str, sptr2 );
	 break;
      }
      this->operation = CONST_OP;

   } else {  /*  Not a constant  */

     Allocate_Ptrs( lParse, this );

      if( !lParse->status ) {

	 rows = lParse->nRows;
	 switch( this->operation ) {

	    /*  Compare Strings  */

	 case NE:
	 case EQ:
	    while( rows-- ) {
	       if( !const1 ) null1 = that1->value.undef[rows];
	       if( !const2 ) null2 = that2->value.undef[rows];
	       this->value.undef[rows] = (null1 || null2);
	       if( ! this->value.undef[rows] ) {
		  if( !const1 ) sptr1  = that1->value.data.strptr[rows];
		  if( !const2 ) sptr2  = that2->value.data.strptr[rows];
		  val = ( FSTRCMP( sptr1, sptr2 ) == 0 );
		  this->value.data.logptr[rows] =
		     ( this->operation==EQ ? val : !val );
	       }
	    }
	    break;
	    
	 case GT:
	 case LT:
	    while( rows-- ) {
	       if( !const1 ) null1 = that1->value.undef[rows];
	       if( !const2 ) null2 = that2->value.undef[rows];
	       this->value.undef[rows] = (null1 || null2);
	       if( ! this->value.undef[rows] ) {
		  if( !const1 ) sptr1  = that1->value.data.strptr[rows];
		  if( !const2 ) sptr2  = that2->value.data.strptr[rows];
		  val = ( FSTRCMP( sptr1, sptr2 ) );
		  this->value.data.logptr[rows] =
		     ( this->operation==GT ? val>0 : val<0 );
	       }
	    }
	    break;

	 case GTE:
	 case LTE:
	    while( rows-- ) {
	       if( !const1 ) null1 = that1->value.undef[rows];
	       if( !const2 ) null2 = that2->value.undef[rows];
	       this->value.undef[rows] = (null1 || null2);
	       if( ! this->value.undef[rows] ) {
		  if( !const1 ) sptr1  = that1->value.data.strptr[rows];
		  if( !const2 ) sptr2  = that2->value.data.strptr[rows];
		  val = ( FSTRCMP( sptr1, sptr2 ) );
		  this->value.data.logptr[rows] =
		     ( this->operation==GTE ? val>=0 : val<=0 );
	       }
	    }
	    break;

	    /*  Concat Strings  */
	    
	 case '+':
	    while( rows-- ) {
	       if( !const1 ) null1 = that1->value.undef[rows];
	       if( !const2 ) null2 = that2->value.undef[rows];
	       this->value.undef[rows] = (null1 || null2);
	       if( ! this->value.undef[rows] ) {
		  if( !const1 ) sptr1  = that1->value.data.strptr[rows];
		  if( !const2 ) sptr2  = that2->value.data.strptr[rows];
		  strcpy( this->value.data.strptr[rows], sptr1 );
		  strcat( this->value.data.strptr[rows], sptr2 );
	       }
	    }
	    break;
	 }
      }
   }

   if( that1->operation>0 ) {
      free( that1->value.data.strptr[0] );
      free( that1->value.data.strptr );
   }
   if( that2->operation>0 ) {
      free( that2->value.data.strptr[0] );
      free( that2->value.data.strptr );
   }
}

static void Do_BinOp_log( ParseData *lParse, Node *this )
{
   Node *that1, *that2;
   int vector1, vector2;
   char val1=0, val2=0, null1=0, null2=0;
   long rows, nelem, elem;

   that1 = lParse->Nodes + this->SubNodes[0];
   that2 = lParse->Nodes + this->SubNodes[1];

   vector1 = ( that1->operation!=CONST_OP );
   if( vector1 )
      vector1 = that1->value.nelem;
   else {
      val1  = that1->value.data.log;
   }

   vector2 = ( that2->operation!=CONST_OP );
   if( vector2 )
      vector2 = that2->value.nelem;
   else {
      val2  = that2->value.data.log;
   }

   if( !vector1 && !vector2 ) {  /*  Result is a constant  */
      switch( this->operation ) {
      case OR:
	 this->value.data.log = (val1 || val2);
	 break;
      case AND:
	 this->value.data.log = (val1 && val2);
	 break;
      case EQ:
	 this->value.data.log = ( (val1 && val2) || (!val1 && !val2) );
	 break;
      case NE:
	 this->value.data.log = ( (val1 && !val2) || (!val1 && val2) );
	 break;
      case ACCUM:
	 this->value.data.lng = val1;
	 break;
      }
      this->operation=CONST_OP;
   } else if (this->operation == ACCUM) {
      long i, previous, curr;
      rows  = lParse->nRows;
      nelem = this->value.nelem;
      elem  = this->value.nelem * rows;
      
      Allocate_Ptrs( lParse, this );
      
      if( !lParse->status ) {
	previous = that2->value.data.lng;
	
	/* Cumulative sum of this chunk */
	for (i=0; i<elem; i++) {
	  if (!that1->value.undef[i]) {
	    curr = that1->value.data.logptr[i];
	    previous += curr;
	  }
	  this->value.data.lngptr[i] = previous;
	  this->value.undef[i] = 0;
	}
	
	/* Store final cumulant for next pass */
	that2->value.data.lng = previous;
      }
      
   } else {
      rows  = lParse->nRows;
      nelem = this->value.nelem;
      elem  = this->value.nelem * rows;

      Allocate_Ptrs( lParse, this );

      if( !lParse->status ) {
	
	 if (this->operation == ACCUM) {
	   long i, previous, curr;
	   
	   previous = that2->value.data.lng;
	   
	   /* Cumulative sum of this chunk */
	   for (i=0; i<elem; i++) {
	     if (!that1->value.undef[i]) {
	       curr = that1->value.data.logptr[i];
	       previous += curr;
	     }
	     this->value.data.lngptr[i] = previous;
	     this->value.undef[i] = 0;
	   }
	   
	   /* Store final cumulant for next pass */
	   that2->value.data.lng = previous;
	 }
	
	 while( rows-- ) {
	    while( nelem-- ) {
	       elem--;

	       if( vector1>1 ) {
		  val1  = that1->value.data.logptr[elem];
		  null1 = that1->value.undef[elem];
	       } else if( vector1 ) {
		  val1  = that1->value.data.logptr[rows];
		  null1 = that1->value.undef[rows];
	       }

	       if( vector2>1 ) {
		  val2  = that2->value.data.logptr[elem];
		  null2 = that2->value.undef[elem];
	       } else if( vector2 ) {
		  val2  = that2->value.data.logptr[rows];
		  null2 = that2->value.undef[rows];
	       }

	       this->value.undef[elem] = (null1 || null2);
	       switch( this->operation ) {

	       case OR:
		  /*  This is more complicated than others to suppress UNDEFs */
		  /*  in those cases where the other argument is DEF && TRUE  */

		  if( !null1 && !null2 ) {
		     this->value.data.logptr[elem] = (val1 || val2);
		  } else if( (null1 && !null2 && val2)
			     || ( !null1 && null2 && val1 ) ) {
		     this->value.data.logptr[elem] = 1;
		     this->value.undef[elem] = 0;
		  }
		  break;

	       case AND:
		  /*  This is more complicated than others to suppress UNDEFs */
		  /*  in those cases where the other argument is DEF && FALSE */

		  if( !null1 && !null2 ) {
		     this->value.data.logptr[elem] = (val1 && val2);
		  } else if( (null1 && !null2 && !val2)
			     || ( !null1 && null2 && !val1 ) ) {
		     this->value.data.logptr[elem] = 0;
		     this->value.undef[elem] = 0;
		  }
		  break;

	       case EQ:
		  this->value.data.logptr[elem] = 
		     ( (val1 && val2) || (!val1 && !val2) );
		  break;

	       case NE:
		  this->value.data.logptr[elem] =
		     ( (val1 && !val2) || (!val1 && val2) );
		  break;
	       }
	    }
	    nelem = this->value.nelem;
	 }
      }
   }

   if( that1->operation>0 ) {
      free( that1->value.data.ptr );
   }
   if( that2->operation>0 ) {
      free( that2->value.data.ptr );
   }
}

static void Do_BinOp_lng( ParseData *lParse, Node *this )
{
   Node *that1, *that2;
   int  vector1, vector2;
   long val1=0, val2=0;
   char null1=0, null2=0;
   long rows, nelem, elem;

   that1 = lParse->Nodes + this->SubNodes[0];
   that2 = lParse->Nodes + this->SubNodes[1];

   vector1 = ( that1->operation!=CONST_OP );
   if( vector1 )
      vector1 = that1->value.nelem;
   else {
      val1  = that1->value.data.lng;
   }

   vector2 = ( that2->operation!=CONST_OP );
   if( vector2 )
      vector2 = that2->value.nelem;
   else {
      val2  = that2->value.data.lng;
   }

   if( !vector1 && !vector2 ) {  /*  Result is a constant  */

      switch( this->operation ) {
      case '~':   /* Treat as == for LONGS */
      case EQ:    this->value.data.log = (val1 == val2);   break;
      case NE:    this->value.data.log = (val1 != val2);   break;
      case GT:    this->value.data.log = (val1 >  val2);   break;
      case LT:    this->value.data.log = (val1 <  val2);   break;
      case LTE:   this->value.data.log = (val1 <= val2);   break;
      case GTE:   this->value.data.log = (val1 >= val2);   break;

      case '+':   this->value.data.lng = (val1  + val2);   break;
      case '-':   this->value.data.lng = (val1  - val2);   break;
      case '*':   this->value.data.lng = (val1  * val2);   break;

      case '&':   this->value.data.lng = (val1  & val2);   break;
      case '|':   this->value.data.lng = (val1  | val2);   break;
      case '^':   this->value.data.lng = (val1  ^ val2);   break;

      case '%':
	 if( val2 ) this->value.data.lng = (val1 % val2);
	 else       yyerror(0, lParse, "Divide by Zero");
	 break;
      case '/': 
	 if( val2 ) this->value.data.lng = (val1 / val2); 
	 else       yyerror(0, lParse, "Divide by Zero");
	 break;
      case POWER:
	 this->value.data.lng = (long)pow((double)val1,(double)val2);
	 break;
      case ACCUM:
	 this->value.data.lng = val1;
	 break;
      case DIFF:
	 this->value.data.lng = 0;
	 break;
      }
      this->operation=CONST_OP;

   } else if ((this->operation == ACCUM) || (this->operation == DIFF)) {
      long i, previous, curr;
      long undef;
      rows  = lParse->nRows;
      nelem = this->value.nelem;
      elem  = this->value.nelem * rows;
      
      Allocate_Ptrs( lParse, this );
      
      if( !lParse->status ) {
	previous = that2->value.data.lng;
	undef    = (long) that2->value.undef;
	
	if (this->operation == ACCUM) {
	  /* Cumulative sum of this chunk */
	  for (i=0; i<elem; i++) {
	    if (!that1->value.undef[i]) {
	      curr = that1->value.data.lngptr[i];
	      previous += curr;
	    }
	    this->value.data.lngptr[i] = previous;
	    this->value.undef[i] = 0;
	  }
	} else {
	  /* Sequential difference for this chunk */
	  for (i=0; i<elem; i++) {
	    curr = that1->value.data.lngptr[i];
	    if (that1->value.undef[i] || undef) {
	      /* Either this, or previous, value was undefined */
	      this->value.data.lngptr[i] = 0;
	      this->value.undef[i] = 1;
	    } else {
	      /* Both defined, we are okay! */
	      this->value.data.lngptr[i] = curr - previous;
	      this->value.undef[i] = 0;
	    }

	    previous = curr;
	    undef = that1->value.undef[i];
	  }
	}	  
	
	/* Store final cumulant for next pass */
	that2->value.data.lng = previous;
	that2->value.undef    = (char *) undef; /* XXX evil, but no harm here */
      }
      
   } else {

      rows  = lParse->nRows;
      nelem = this->value.nelem;
      elem  = this->value.nelem * rows;

      Allocate_Ptrs( lParse, this );

      while( rows-- && !lParse->status ) {
	 while( nelem-- && !lParse->status ) {
	    elem--;

	    if( vector1>1 ) {
	       val1  = that1->value.data.lngptr[elem];
	       null1 = that1->value.undef[elem];
	    } else if( vector1 ) {
	       val1  = that1->value.data.lngptr[rows];
	       null1 = that1->value.undef[rows];
	    }

	    if( vector2>1 ) {
	       val2  = that2->value.data.lngptr[elem];
	       null2 = that2->value.undef[elem];
	    } else if( vector2 ) {
	       val2  = that2->value.data.lngptr[rows];
	       null2 = that2->value.undef[rows];
	    }

	    this->value.undef[elem] = (null1 || null2);
	    switch( this->operation ) {
	    case '~':   /* Treat as == for LONGS */
	    case EQ:   this->value.data.logptr[elem] = (val1 == val2);   break;
	    case NE:   this->value.data.logptr[elem] = (val1 != val2);   break;
	    case GT:   this->value.data.logptr[elem] = (val1 >  val2);   break;
	    case LT:   this->value.data.logptr[elem] = (val1 <  val2);   break;
	    case LTE:  this->value.data.logptr[elem] = (val1 <= val2);   break;
	    case GTE:  this->value.data.logptr[elem] = (val1 >= val2);   break;
	       
	    case '+':  this->value.data.lngptr[elem] = (val1  + val2);   break;
	    case '-':  this->value.data.lngptr[elem] = (val1  - val2);   break;
	    case '*':  this->value.data.lngptr[elem] = (val1  * val2);   break;

	    case '&':  this->value.data.lngptr[elem] = (val1  & val2);   break;
	    case '|':  this->value.data.lngptr[elem] = (val1  | val2);   break;
	    case '^':  this->value.data.lngptr[elem] = (val1  ^ val2);   break;

	    case '%':   
	       if( val2 ) this->value.data.lngptr[elem] = (val1 % val2);
	       else {
		 this->value.data.lngptr[elem] = 0;
		 this->value.undef[elem] = 1;
	       }
	       break;
	    case '/': 
	       if( val2 ) this->value.data.lngptr[elem] = (val1 / val2); 
	       else {
		 this->value.data.lngptr[elem] = 0;
		 this->value.undef[elem] = 1;
	       }
	       break;
	    case POWER:
	       this->value.data.lngptr[elem] = (long)pow((double)val1,(double)val2);
	       break;
	    }
	 }
	 nelem = this->value.nelem;
      }
   }

   if( that1->operation>0 ) {
      free( that1->value.data.ptr );
   }
   if( that2->operation>0 ) {
      free( that2->value.data.ptr );
   }
}

static void Do_BinOp_dbl( ParseData *lParse, Node *this )
{
   Node   *that1, *that2;
   int    vector1, vector2;
   double val1=0.0, val2=0.0;
   char   null1=0, null2=0;
   long   rows, nelem, elem;

   that1 = lParse->Nodes + this->SubNodes[0];
   that2 = lParse->Nodes + this->SubNodes[1];

   vector1 = ( that1->operation!=CONST_OP );
   if( vector1 )
      vector1 = that1->value.nelem;
   else {
      val1  = that1->value.data.dbl;
   }

   vector2 = ( that2->operation!=CONST_OP );
   if( vector2 )
      vector2 = that2->value.nelem;
   else {
      val2  = that2->value.data.dbl;
   } 

   if( !vector1 && !vector2 ) {  /*  Result is a constant  */

      switch( this->operation ) {
      case '~':   this->value.data.log = ( fabs(val1-val2) < APPROX );   break;
      case EQ:    this->value.data.log = (val1 == val2);   break;
      case NE:    this->value.data.log = (val1 != val2);   break;
      case GT:    this->value.data.log = (val1 >  val2);   break;
      case LT:    this->value.data.log = (val1 <  val2);   break;
      case LTE:   this->value.data.log = (val1 <= val2);   break;
      case GTE:   this->value.data.log = (val1 >= val2);   break;

      case '+':   this->value.data.dbl = (val1  + val2);   break;
      case '-':   this->value.data.dbl = (val1  - val2);   break;
      case '*':   this->value.data.dbl = (val1  * val2);   break;

      case '%':
	 if( val2 ) this->value.data.dbl = val1 - val2*((int)(val1/val2));
	 else       yyerror(0, lParse, "Divide by Zero");
	 break;
      case '/': 
	 if( val2 ) this->value.data.dbl = (val1 / val2); 
	 else       yyerror(0, lParse, "Divide by Zero");
	 break;
      case POWER:
	 this->value.data.dbl = (double)pow(val1,val2);
	 break;
      case ACCUM:
	 this->value.data.dbl = val1;
	 break;
      case DIFF:
	this->value.data.dbl = 0;
	 break;
      }
      this->operation=CONST_OP;

   } else if ((this->operation == ACCUM) || (this->operation == DIFF)) {
      long i;
      long undef;
      double previous, curr;
      rows  = lParse->nRows;
      nelem = this->value.nelem;
      elem  = this->value.nelem * rows;
      
      Allocate_Ptrs( lParse, this );
      
      if( !lParse->status ) {
	previous = that2->value.data.dbl;
	undef    = (long) that2->value.undef;
	
	if (this->operation == ACCUM) {
	  /* Cumulative sum of this chunk */
	  for (i=0; i<elem; i++) {
	    if (!that1->value.undef[i]) {
	      curr = that1->value.data.dblptr[i];
	      previous += curr;
	    }
	    this->value.data.dblptr[i] = previous;
	    this->value.undef[i] = 0;
	  }
	} else {
	  /* Sequential difference for this chunk */
	  for (i=0; i<elem; i++) {
	    curr = that1->value.data.dblptr[i];
	    if (that1->value.undef[i] || undef) {
	      /* Either this, or previous, value was undefined */
	      this->value.data.dblptr[i] = 0;
	      this->value.undef[i] = 1;
	    } else {
	      /* Both defined, we are okay! */
	      this->value.data.dblptr[i] = curr - previous;
	      this->value.undef[i] = 0;
	    }

	    previous = curr;
	    undef = that1->value.undef[i];
	  }
	}	  
	
	/* Store final cumulant for next pass */
	that2->value.data.dbl = previous;
	that2->value.undef    = (char *) undef; /* XXX evil, but no harm here */
      }
      
   } else {

      rows  = lParse->nRows;
      nelem = this->value.nelem;
      elem  = this->value.nelem * rows;

      Allocate_Ptrs( lParse, this );

      while( rows-- && !lParse->status ) {
	 while( nelem-- && !lParse->status ) {
	    elem--;

	    if( vector1>1 ) {
	       val1  = that1->value.data.dblptr[elem];
	       null1 = that1->value.undef[elem];
	    } else if( vector1 ) {
	       val1  = that1->value.data.dblptr[rows];
	       null1 = that1->value.undef[rows];
	    }

	    if( vector2>1 ) {
	       val2  = that2->value.data.dblptr[elem];
	       null2 = that2->value.undef[elem];
	    } else if( vector2 ) {
	       val2  = that2->value.data.dblptr[rows];
	       null2 = that2->value.undef[rows];
	    }

	    this->value.undef[elem] = (null1 || null2);
	    switch( this->operation ) {
	    case '~':   this->value.data.logptr[elem] =
                                          ( fabs(val1-val2) < APPROX );   break;
	    case EQ:    this->value.data.logptr[elem] = (val1 == val2);   break;
	    case NE:    this->value.data.logptr[elem] = (val1 != val2);   break;
	    case GT:    this->value.data.logptr[elem] = (val1 >  val2);   break;
	    case LT:    this->value.data.logptr[elem] = (val1 <  val2);   break;
	    case LTE:   this->value.data.logptr[elem] = (val1 <= val2);   break;
	    case GTE:   this->value.data.logptr[elem] = (val1 >= val2);   break;
	       
	    case '+':   this->value.data.dblptr[elem] = (val1  + val2);   break;
	    case '-':   this->value.data.dblptr[elem] = (val1  - val2);   break;
	    case '*':   this->value.data.dblptr[elem] = (val1  * val2);   break;

	    case '%':
	       if( val2 ) this->value.data.dblptr[elem] =
                                val1 - val2*((int)(val1/val2));
	       else {
		 this->value.data.dblptr[elem] = 0.0;
		 this->value.undef[elem] = 1;
	       }
	       break;
	    case '/': 
	       if( val2 ) this->value.data.dblptr[elem] = (val1 / val2); 
	       else {
		 this->value.data.dblptr[elem] = 0.0;
		 this->value.undef[elem] = 1;
	       }
	       break;
	    case POWER:
	       this->value.data.dblptr[elem] = (double)pow(val1,val2);
	       break;
	    }
	 }
	 nelem = this->value.nelem;
      }
   }

   if( that1->operation>0 ) {
      free( that1->value.data.ptr );
   }
   if( that2->operation>0 ) {
      free( that2->value.data.ptr );
   }
}

/*
 *  This Quickselect routine is based on the algorithm described in
 *  "Numerical recipes in C", Second Edition,
 *  Cambridge University Press, 1992, Section 8.5, ISBN 0-521-43108-5
 *  This code by Nicolas Devillard - 1998. Public domain.
 * http://ndevilla.free.fr/median/median/src/quickselect.c
 */

#define ELEM_SWAP(a,b) { register long t=(a);(a)=(b);(b)=t; }

/* 
 * qselect_median_lng - select the median value of a long array
 *
 * This routine selects the median value of the long integer array
 * arr[].  If there are an even number of elements, the "lower median"
 * is selected.
 *
 * The array arr[] is scrambled, so users must operate on a scratch
 * array if they wish the values to be preserved.
 *
 * long arr[] - array of values
 * int n - number of elements in arr
 *
 * RETURNS: the lower median value of arr[]
 *
 */
long qselect_median_lng(long arr[], int n)
{
    int low, high ;
    int median;
    int middle, ll, hh;

    low = 0 ; high = n-1 ; median = (low + high) / 2;
    for (;;) {

        if (high <= low) { /* One element only */
	  return arr[median];	  
	}

        if (high == low + 1) {  /* Two elements only */
            if (arr[low] > arr[high])
                ELEM_SWAP(arr[low], arr[high]) ;
	    return arr[median];
        }

    /* Find median of low, middle and high items; swap into position low */
    middle = (low + high) / 2;
    if (arr[middle] > arr[high])    ELEM_SWAP(arr[middle], arr[high]) ;
    if (arr[low] > arr[high])       ELEM_SWAP(arr[low], arr[high]) ;
    if (arr[middle] > arr[low])     ELEM_SWAP(arr[middle], arr[low]) ;

    /* Swap low item (now in position middle) into position (low+1) */
    ELEM_SWAP(arr[middle], arr[low+1]) ;

    /* Nibble from each end towards middle, swapping items when stuck */
    ll = low + 1;
    hh = high;
    for (;;) {
        do ll++; while (arr[low] > arr[ll]) ;
        do hh--; while (arr[hh]  > arr[low]) ;

        if (hh < ll)
        break;

        ELEM_SWAP(arr[ll], arr[hh]) ;
    }

    /* Swap middle item (in position low) back into correct position */
    ELEM_SWAP(arr[low], arr[hh]) ;

    /* Re-set active partition */
    if (hh <= median)
        low = ll;
        if (hh >= median)
        high = hh - 1;
    }
}

#undef ELEM_SWAP

#define ELEM_SWAP(a,b) { register double t=(a);(a)=(b);(b)=t; }

/* 
 * qselect_median_dbl - select the median value of a double array
 *
 * This routine selects the median value of the double array
 * arr[].  If there are an even number of elements, the "lower median"
 * is selected.
 *
 * The array arr[] is scrambled, so users must operate on a scratch
 * array if they wish the values to be preserved.
 *
 * double arr[] - array of values
 * int n - number of elements in arr
 *
 * RETURNS: the lower median value of arr[]
 *
 */
double qselect_median_dbl(double arr[], int n)
{
    int low, high ;
    int median;
    int middle, ll, hh;

    low = 0 ; high = n-1 ; median = (low + high) / 2;
    for (;;) {
        if (high <= low) { /* One element only */
            return arr[median] ;
	}

        if (high == low + 1) {  /* Two elements only */
            if (arr[low] > arr[high])
                ELEM_SWAP(arr[low], arr[high]) ;
            return arr[median] ;
        }

    /* Find median of low, middle and high items; swap into position low */
    middle = (low + high) / 2;
    if (arr[middle] > arr[high])    ELEM_SWAP(arr[middle], arr[high]) ;
    if (arr[low] > arr[high])       ELEM_SWAP(arr[low], arr[high]) ;
    if (arr[middle] > arr[low])     ELEM_SWAP(arr[middle], arr[low]) ;

    /* Swap low item (now in position middle) into position (low+1) */
    ELEM_SWAP(arr[middle], arr[low+1]) ;

    /* Nibble from each end towards middle, swapping items when stuck */
    ll = low + 1;
    hh = high;
    for (;;) {
        do ll++; while (arr[low] > arr[ll]) ;
        do hh--; while (arr[hh]  > arr[low]) ;

        if (hh < ll)
        break;

        ELEM_SWAP(arr[ll], arr[hh]) ;
    }

    /* Swap middle item (in position low) back into correct position */
    ELEM_SWAP(arr[low], arr[hh]) ;

    /* Re-set active partition */
    if (hh <= median)
        low = ll;
        if (hh >= median)
        high = hh - 1;
    }
}

#undef ELEM_SWAP

/*
 * angsep_calc - compute angular separation between celestial coordinates
 *   
 * This routine computes the angular separation between to coordinates
 * on the celestial sphere (i.e. RA and Dec).  Note that all units are
 * in DEGREES, unlike the other trig functions in the calculator.
 *
 * double ra1, dec1 - RA and Dec of the first position in degrees
 * double ra2, dec2 - RA and Dec of the second position in degrees
 * 
 * RETURNS: (double) angular separation in degrees
 *
 */
double angsep_calc(double ra1, double dec1, double ra2, double dec2)
{
/*  double cd;  */
  static double deg = 0;
  double a, sdec, sra;
  
  if (deg == 0) deg = ((double)4)*atan((double)1)/((double)180);
  /* deg = 1.0; **** UNCOMMENT IF YOU WANT RADIANS */

  /* The algorithm is the law of Haversines.  This algorithm is
     stable even when the points are close together.  The normal
     Law of Cosines fails for angles around 0.1 arcsec. */

  sra  = sin( (ra2 - ra1)*deg / 2 );
  sdec = sin( (dec2 - dec1)*deg / 2);
  a = sdec*sdec + cos(dec1*deg)*cos(dec2*deg)*sra*sra;

  /* Sanity checking to avoid a range error in the sqrt()'s below */
  if (a < 0) { a = 0; }
  if (a > 1) { a = 1; }

  return 2.0*atan2(sqrt(a), sqrt(1.0 - a)) / deg;
}

static void Do_Func( ParseData *lParse, Node *this )
{
   Node *theParams[MAXSUBS];
   int  vector[MAXSUBS], allConst;
   lval pVals[MAXSUBS];
   char pNull[MAXSUBS];
   long   ival;
   double dval;
   int  i, valInit;
   long row, elem, nelem;

   i = this->nSubNodes;
   allConst = 1;
   while( i-- ) {
      theParams[i] = lParse->Nodes + this->SubNodes[i];
      vector[i]   = ( theParams[i]->operation!=CONST_OP );
      if( vector[i] ) {
	 allConst = 0;
	 vector[i] = theParams[i]->value.nelem;
      } else {
	 if( theParams[i]->type==DOUBLE ) {
	    pVals[i].data.dbl = theParams[i]->value.data.dbl;
	 } else if( theParams[i]->type==LONG ) {
	    pVals[i].data.lng = theParams[i]->value.data.lng;
	 } else if( theParams[i]->type==BOOLEAN ) {
	    pVals[i].data.log = theParams[i]->value.data.log;
	 } else
	    strcpy(pVals[i].data.str, theParams[i]->value.data.str);
	 pNull[i] = 0;
      }
   }

   if( this->nSubNodes==0 ) allConst = 0; /* These do produce scalars */
   /* Random numbers are *never* constant !! */
   if( this->operation == poirnd_fct ) allConst = 0;
   if( this->operation == gasrnd_fct ) allConst = 0;
   if( this->operation == rnd_fct ) allConst = 0;

   if( allConst ) {

      switch( this->operation ) {

	    /* Non-Trig single-argument functions */

	 case sum_fct:
	    if( theParams[0]->type==BOOLEAN )
	       this->value.data.lng = ( pVals[0].data.log ? 1 : 0 );
	    else if( theParams[0]->type==LONG )
	       this->value.data.lng = pVals[0].data.lng;
	    else if( theParams[0]->type==DOUBLE )
	       this->value.data.dbl = pVals[0].data.dbl;
	    else if( theParams[0]->type==BITSTR )
	      strcpy(this->value.data.str, pVals[0].data.str);
	    break;
         case average_fct:
	    if( theParams[0]->type==LONG )
	       this->value.data.dbl = pVals[0].data.lng;
	    else if( theParams[0]->type==DOUBLE )
	       this->value.data.dbl = pVals[0].data.dbl;
	    break;
         case stddev_fct:
	    this->value.data.dbl = 0;  /* Standard deviation of a constant = 0 */
	    break;
	 case median_fct:
	    if( theParams[0]->type==BOOLEAN )
	       this->value.data.lng = ( pVals[0].data.log ? 1 : 0 );
	    else if( theParams[0]->type==LONG )
	       this->value.data.lng = pVals[0].data.lng;
	    else
	       this->value.data.dbl = pVals[0].data.dbl;
	    break;

	 case poirnd_fct:
	    if( theParams[0]->type==DOUBLE )
	      this->value.data.lng = simplerng_getpoisson(pVals[0].data.dbl);
	    else
	      this->value.data.lng = simplerng_getpoisson(pVals[0].data.lng);
	    break;

	 case abs_fct:
	    if( theParams[0]->type==DOUBLE ) {
	       dval = pVals[0].data.dbl;
	       this->value.data.dbl = (dval>0.0 ? dval : -dval);
	    } else {
	       ival = pVals[0].data.lng;
	       this->value.data.lng = (ival> 0  ? ival : -ival);
	    }
	    break;

            /* Special Null-Handling Functions */

         case nonnull_fct:
	    this->value.data.lng = 1; /* Constants are always 1-element and defined */
	    break;
         case isnull_fct:  /* Constants are always defined */
	    this->value.data.log = 0;
	    break;
         case defnull_fct:
	    if( this->type==BOOLEAN )
	       this->value.data.log = pVals[0].data.log;
            else if( this->type==LONG )
	       this->value.data.lng = pVals[0].data.lng;
            else if( this->type==DOUBLE )
	       this->value.data.dbl = pVals[0].data.dbl;
            else if( this->type==STRING )
	       strcpy(this->value.data.str,pVals[0].data.str);
	    break;
        case setnull_fct: /* Only defined for numeric expressions */
            if( this->type==LONG )
 	      this->value.data.lng = pVals[0].data.lng;
            else if( this->type==DOUBLE )
	       this->value.data.dbl = pVals[0].data.dbl;
	    break;

	    /* Math functions with 1 double argument */

	 case sin_fct:
	    this->value.data.dbl = sin( pVals[0].data.dbl );
	    break;
	 case cos_fct:
	    this->value.data.dbl = cos( pVals[0].data.dbl );
	    break;
	 case tan_fct:
	    this->value.data.dbl = tan( pVals[0].data.dbl );
	    break;
	 case asin_fct:
	    dval = pVals[0].data.dbl;
	    if( dval<-1.0 || dval>1.0 )
	       yyerror(0, lParse, "Out of range argument to arcsin");
	    else
	       this->value.data.dbl = asin( dval );
	    break;
	 case acos_fct:
	    dval = pVals[0].data.dbl;
	    if( dval<-1.0 || dval>1.0 )
	       yyerror(0, lParse, "Out of range argument to arccos");
	    else
	       this->value.data.dbl = acos( dval );
	    break;
	 case atan_fct:
	    this->value.data.dbl = atan( pVals[0].data.dbl );
	    break;
	 case sinh_fct:
	    this->value.data.dbl = sinh( pVals[0].data.dbl );
	    break;
	 case cosh_fct:
	    this->value.data.dbl = cosh( pVals[0].data.dbl );
	    break;
	 case tanh_fct:
	    this->value.data.dbl = tanh( pVals[0].data.dbl );
	    break;
	 case exp_fct:
	    this->value.data.dbl = exp( pVals[0].data.dbl );
	    break;
	 case log_fct:
	    dval = pVals[0].data.dbl;
	    if( dval<=0.0 )
	       yyerror(0, lParse, "Out of range argument to log");
	    else
	       this->value.data.dbl = log( dval );
	    break;
	 case log10_fct:
	    dval = pVals[0].data.dbl;
	    if( dval<=0.0 )
	       yyerror(0, lParse, "Out of range argument to log10");
	    else
	       this->value.data.dbl = log10( dval );
	    break;
	 case sqrt_fct:
	    dval = pVals[0].data.dbl;
	    if( dval<0.0 )
	       yyerror(0, lParse, "Out of range argument to sqrt");
	    else
	       this->value.data.dbl = sqrt( dval );
	    break;
	 case ceil_fct:
	    this->value.data.dbl = ceil( pVals[0].data.dbl );
	    break;
	 case floor_fct:
	    this->value.data.dbl = floor( pVals[0].data.dbl );
	    break;
	 case round_fct:
	    this->value.data.dbl = floor( pVals[0].data.dbl + 0.5 );
	    break;

	    /* Two-argument Trig Functions */

	 case atan2_fct:
	    this->value.data.dbl =
	       atan2( pVals[0].data.dbl, pVals[1].data.dbl );
	    break;

	    /* Four-argument ANGSEP function */
         case angsep_fct:
	    this->value.data.dbl = 
	      angsep_calc(pVals[0].data.dbl, pVals[1].data.dbl,
			  pVals[2].data.dbl, pVals[3].data.dbl);

	    /*  Min/Max functions taking 1 or 2 arguments  */

         case min1_fct:
	    /* No constant vectors! */
	    if( this->type == DOUBLE )
	       this->value.data.dbl = pVals[0].data.dbl;
	    else if( this->type == LONG )
	       this->value.data.lng = pVals[0].data.lng;
	    else if( this->type == BITSTR )
	      strcpy(this->value.data.str, pVals[0].data.str);
	    break;
         case min2_fct:
	    if( this->type == DOUBLE )
	       this->value.data.dbl =
		  minvalue( pVals[0].data.dbl, pVals[1].data.dbl );
	    else if( this->type == LONG )
	       this->value.data.lng =
		  minvalue( pVals[0].data.lng, pVals[1].data.lng );
	    break;
         case max1_fct:
	    /* No constant vectors! */
	    if( this->type == DOUBLE )
	       this->value.data.dbl = pVals[0].data.dbl;
	    else if( this->type == LONG )
	       this->value.data.lng = pVals[0].data.lng;
	    else if( this->type == BITSTR )
	      strcpy(this->value.data.str, pVals[0].data.str);
	    break;
         case max2_fct:
	    if( this->type == DOUBLE )
	       this->value.data.dbl =
		  maxvalue( pVals[0].data.dbl, pVals[1].data.dbl );
	    else if( this->type == LONG )
	       this->value.data.lng =
		  maxvalue( pVals[0].data.lng, pVals[1].data.lng );
	    break;

	    /* Boolean SAO region Functions... scalar or vector dbls */

	 case near_fct:
	    this->value.data.log = bnear( pVals[0].data.dbl, pVals[1].data.dbl,
					  pVals[2].data.dbl );
	    break;
	 case circle_fct:
	    this->value.data.log = circle( pVals[0].data.dbl, pVals[1].data.dbl,
					   pVals[2].data.dbl, pVals[3].data.dbl,
					   pVals[4].data.dbl );
	    break;
	 case box_fct:
	    this->value.data.log = saobox( pVals[0].data.dbl, pVals[1].data.dbl,
					   pVals[2].data.dbl, pVals[3].data.dbl,
					   pVals[4].data.dbl, pVals[5].data.dbl,
					   pVals[6].data.dbl );
	    break;
	 case elps_fct:
	    this->value.data.log =
                               ellipse( pVals[0].data.dbl, pVals[1].data.dbl,
					pVals[2].data.dbl, pVals[3].data.dbl,
					pVals[4].data.dbl, pVals[5].data.dbl,
					pVals[6].data.dbl );
	    break;

            /* C Conditional expression:  bool ? expr : expr */

         case ifthenelse_fct:
            switch( this->type ) {
            case BOOLEAN:
               this->value.data.log = ( pVals[2].data.log ?
                                        pVals[0].data.log : pVals[1].data.log );
               break;
            case LONG:
               this->value.data.lng = ( pVals[2].data.log ?
                                        pVals[0].data.lng : pVals[1].data.lng );
               break;
            case DOUBLE:
               this->value.data.dbl = ( pVals[2].data.log ?
                                        pVals[0].data.dbl : pVals[1].data.dbl );
               break;
            case STRING:
	       strcpy(this->value.data.str, ( pVals[2].data.log ?
                                              pVals[0].data.str :
                                              pVals[1].data.str ) );
               break;
            }
            break;

	    /* String functions */
         case strmid_fct:
	   cstrmid(lParse, 
		   this->value.data.str, this->value.nelem, 
		   pVals[0].data.str,    pVals[0].nelem,
		   pVals[1].data.lng);
	   break;
         case strpos_fct:
	   {
	     char *res = strstr(pVals[0].data.str, pVals[1].data.str);
	     if (res == NULL) {
	       this->value.data.lng = 0; 
	     } else {
	       this->value.data.lng = (res - pVals[0].data.str) + 1;
	     }
	     break;
	   }

      }
      this->operation = CONST_OP;

   } else {

     Allocate_Ptrs( lParse, this );

      row  = lParse->nRows;
      elem = row * this->value.nelem;

      if( !lParse->status ) {
	 switch( this->operation ) {

	    /* Special functions with no arguments */

	 case row_fct:
	    while( row-- ) {
	       this->value.data.lngptr[row] = lParse->firstRow + row;
	       this->value.undef[row] = 0;
	    }
	    break;
	 case null_fct:
            if( this->type==LONG ) {
               while( row-- ) {
                  this->value.data.lngptr[row] = 0;
                  this->value.undef[row] = 1;
               }
            } else if( this->type==STRING ) {
               while( row-- ) {
                  this->value.data.strptr[row][0] = '\0';
                  this->value.undef[row] = 1;
               }
            }
	    break;
	 case axiselem_fct:
	   {
	     long ielem;
	     long iaxis[MAXDIMS] = {1, 1, 1, 1, 1};
	     long ipos = pVals[1].data.lng - 1; /* This should be a constant long value */
	     int naxis = this->value.naxis;
	     int j;
	     if (ipos < 0 || ipos >= MAXDIMS) {
	         yyerror(0, lParse, "AXISELEM(V,n) n value exceeded maximum dimension");
		 free( this->value.data.ptr );
		 break;
	     }

	     for (ielem = 0; ielem<elem; ielem++) {
	       this->value.data.lngptr[ielem] = iaxis[ipos];
	       this->value.undef[ielem] = 0;
	       iaxis[0]++;
	       for (j = 0; j < naxis; j++) {
		 if (iaxis[j] > this->value.naxes[j]) { 
		   iaxis[j] = 1; 
		   if (j < (naxis-1)) iaxis[j+1]++;
		 } else {
		   break;
		 }
	       }

	     }
	   }
	   break;
	 case elemnum_fct:
	   {
	     long ielem;
	     long elemnum = 1;
	     int j;

	     for (ielem = 0; ielem<elem; ielem++) {
	       this->value.data.lngptr[ielem] = elemnum;
	       this->value.undef[ielem] = 0;
	       elemnum ++;
	       if (elemnum > this->value.nelem) elemnum = 1;
	     }
	   }
	   break;
	 case rnd_fct:
	   while( elem-- ) {
	     this->value.data.dblptr[elem] = simplerng_getuniform();
	     this->value.undef[elem] = 0;
	    }
	    break;

	 case gasrnd_fct:
	    while( elem-- ) {
	       this->value.data.dblptr[elem] = simplerng_getnorm();
	       this->value.undef[elem] = 0;
	    }
	    break;

	 case poirnd_fct:
	   if( theParams[0]->type==DOUBLE ) {
	      if (theParams[0]->operation == CONST_OP) {
		while( elem-- ) {
		  this->value.undef[elem] = (pVals[0].data.dbl < 0);
		  if (! this->value.undef[elem]) {
		    this->value.data.lngptr[elem] = simplerng_getpoisson(pVals[0].data.dbl);
		  }
		} 
	      } else {
		while( elem-- ) {
		  this->value.undef[elem] = theParams[0]->value.undef[elem];
		  if (theParams[0]->value.data.dblptr[elem] < 0) 
		    this->value.undef[elem] = 1;
		  if (! this->value.undef[elem]) {
		    this->value.data.lngptr[elem] = 
		      simplerng_getpoisson(theParams[0]->value.data.dblptr[elem]);
		  }
		} /* while */
	      } /* ! CONST_OP */
	   } else {
	     /* LONG */
	      if (theParams[0]->operation == CONST_OP) {
		while( elem-- ) {
		  this->value.undef[elem] = (pVals[0].data.lng < 0);
		  if (! this->value.undef[elem]) {
		    this->value.data.lngptr[elem] = simplerng_getpoisson(pVals[0].data.lng);
		  }
		} 
	      } else {
		while( elem-- ) {
		  this->value.undef[elem] = theParams[0]->value.undef[elem];
		  if (theParams[0]->value.data.lngptr[elem] < 0) 
		    this->value.undef[elem] = 1;
		  if (! this->value.undef[elem]) {
		    this->value.data.lngptr[elem] = 
		      simplerng_getpoisson(theParams[0]->value.data.lngptr[elem]);
		  }
		} /* while */
	      } /* ! CONST_OP */
	   } /* END LONG */
	   break;


	    /* Non-Trig single-argument functions */
	    
	 case sum_fct:
	    elem = row * theParams[0]->value.nelem;
	    if( theParams[0]->type==BOOLEAN ) {
	       while( row-- ) {
		  this->value.data.lngptr[row] = 0;
		  /* Default is UNDEF until a defined value is found */
		  this->value.undef[row] = 1;
		  nelem = theParams[0]->value.nelem;
		  while( nelem-- ) {
		     elem--;
		     if ( ! theParams[0]->value.undef[elem] ) {
		       this->value.data.lngptr[row] +=
			 ( theParams[0]->value.data.logptr[elem] ? 1 : 0 );
		       this->value.undef[row] = 0;
		     }
		  }
	       }
	    } else if( theParams[0]->type==LONG ) {
	       while( row-- ) {
		  this->value.data.lngptr[row] = 0;
		  /* Default is UNDEF until a defined value is found */
		  this->value.undef[row] = 1;
		  nelem = theParams[0]->value.nelem;
		  while( nelem-- ) {
		     elem--;
		     if ( ! theParams[0]->value.undef[elem] ) {
		       this->value.data.lngptr[row] +=
			 theParams[0]->value.data.lngptr[elem];
		       this->value.undef[row] = 0;
		     }
		  }
	       }		  
	    } else if( theParams[0]->type==DOUBLE ){
	       while( row-- ) {
		  this->value.data.dblptr[row] = 0.0;
		  /* Default is UNDEF until a defined value is found */
		  this->value.undef[row] = 1;
		  nelem = theParams[0]->value.nelem;
		  while( nelem-- ) {
		     elem--;
		     if ( ! theParams[0]->value.undef[elem] ) {
		       this->value.data.dblptr[row] +=
			 theParams[0]->value.data.dblptr[elem];
		       this->value.undef[row] = 0;
		     }
		  }
	       }		  
	    } else { /* BITSTR */
	       nelem = theParams[0]->value.nelem;
	       while( row-- ) {
		  char *sptr1 = theParams[0]->value.data.strptr[row];
		  this->value.data.lngptr[row] = 0;
		  this->value.undef[row] = 0;
		  while (*sptr1) {
		    if (*sptr1 == '1') this->value.data.lngptr[row] ++;
		    sptr1++;
		  }
	       }		  
	    }
	    break;

	 case average_fct:
	    elem = row * theParams[0]->value.nelem;
	    if( theParams[0]->type==LONG ) {
	       while( row-- ) {
		  int count = 0;
		  this->value.data.dblptr[row] = 0;
		  nelem = theParams[0]->value.nelem;
		  while( nelem-- ) {
		     elem--;
		     if (theParams[0]->value.undef[elem] == 0) {
		       this->value.data.dblptr[row] +=
			 theParams[0]->value.data.lngptr[elem];
		       count ++;
		     }
		  }
		  if (count == 0) {
		    this->value.undef[row] = 1;
		  } else {
		    this->value.undef[row] = 0;
		    this->value.data.dblptr[row] /= count;
		  }
	       }		  
	    } else if( theParams[0]->type==DOUBLE ){
	       while( row-- ) {
		  int count = 0;
		  this->value.data.dblptr[row] = 0;
		  nelem = theParams[0]->value.nelem;
		  while( nelem-- ) {
		     elem--;
		     if (theParams[0]->value.undef[elem] == 0) {
		       this->value.data.dblptr[row] +=
			 theParams[0]->value.data.dblptr[elem];
		       count ++;
		     }
		  }
		  if (count == 0) {
		    this->value.undef[row] = 1;
		  } else {
		    this->value.undef[row] = 0;
		    this->value.data.dblptr[row] /= count;
		  }
	       }		  
	    }
	    break;
	 case stddev_fct:
	    elem = row * theParams[0]->value.nelem;
	    if( theParams[0]->type==LONG ) {

	       /* Compute the mean value */
	       while( row-- ) {
		  int count = 0;
		  double sum = 0, sum2 = 0;

		  nelem = theParams[0]->value.nelem;
		  while( nelem-- ) {
		     elem--;
		     if (theParams[0]->value.undef[elem] == 0) {
		       sum += theParams[0]->value.data.lngptr[elem];
		       count ++;
		     }
		  }
		  if (count > 1) {
		    sum /= count;

		    /* Compute the sum of squared deviations */
		    nelem = theParams[0]->value.nelem;
		    elem += nelem;  /* Reset elem for second pass */
		    while( nelem-- ) {
		      elem--;
		      if (theParams[0]->value.undef[elem] == 0) {
			double dx = (theParams[0]->value.data.lngptr[elem] - sum);
			sum2 += (dx*dx);
		      }
		    }

		    sum2 /= (double)count-1;

		    this->value.undef[row] = 0;
		    this->value.data.dblptr[row] = sqrt(sum2);
		  } else {
		    this->value.undef[row] = 0;       /* STDDEV => 0 */
		    this->value.data.dblptr[row] = 0;
		  }
	       }
	    } else if( theParams[0]->type==DOUBLE ){

	       /* Compute the mean value */
	       while( row-- ) {
		  int count = 0;
		  double sum = 0, sum2 = 0;

		  nelem = theParams[0]->value.nelem;
		  while( nelem-- ) {
		     elem--;
		     if (theParams[0]->value.undef[elem] == 0) {
		       sum += theParams[0]->value.data.dblptr[elem];
		       count ++;
		     }
		  }
		  if (count > 1) {
		    sum /= count;

		    /* Compute the sum of squared deviations */
		    nelem = theParams[0]->value.nelem;
		    elem += nelem;  /* Reset elem for second pass */
		    while( nelem-- ) {
		      elem--;
		      if (theParams[0]->value.undef[elem] == 0) {
			double dx = (theParams[0]->value.data.dblptr[elem] - sum);
			sum2 += (dx*dx);
		      }
		    }

		    sum2 /= (double)count-1;

		    this->value.undef[row] = 0;
		    this->value.data.dblptr[row] = sqrt(sum2);
		  } else {
		    this->value.undef[row] = 0;       /* STDDEV => 0 */
		    this->value.data.dblptr[row] = 0;
		  }
	       }
	    }
	    break;

	 case median_fct:
	   elem = row * theParams[0]->value.nelem;
	   nelem = theParams[0]->value.nelem;
	   if( theParams[0]->type==LONG ) {
	       long *dptr = theParams[0]->value.data.lngptr;
	       char *uptr = theParams[0]->value.undef;
	       long *mptr = (long *) malloc(sizeof(long)*nelem);
	       int irow;

	       /* Allocate temporary storage for this row, since the
                  quickselect function will scramble the contents */
	       if (mptr == 0) {
		 yyerror(0, lParse, "Could not allocate temporary memory in median function");
		 free( this->value.data.ptr );
		 break;
	       }

	       for (irow=0; irow<row; irow++) {
		  long *p = mptr;
		  int nelem1 = nelem;


		  while ( nelem1-- ) { 
		    if (*uptr == 0) {
		      *p++ = *dptr;   /* Only advance the dest pointer if we copied */
		    }
		    dptr ++;  /* Advance the source pointer ... */
		    uptr ++;  /* ... and source "undef" pointer */
		  }
		  
		  nelem1 = (p - mptr);  /* Number of accepted data points */
		  if (nelem1 > 0) {
		    this->value.undef[irow] = 0;
		    this->value.data.lngptr[irow] = qselect_median_lng(mptr, nelem1);
		  } else {
		    this->value.undef[irow] = 1;
		    this->value.data.lngptr[irow] = 0;
		  }
		    
	       }		  

	       free(mptr);
	    } else {
	       double *dptr = theParams[0]->value.data.dblptr;
	       char   *uptr = theParams[0]->value.undef;
	       double *mptr = (double *) malloc(sizeof(double)*nelem);
	       int irow;

	       /* Allocate temporary storage for this row, since the
                  quickselect function will scramble the contents */
	       if (mptr == 0) {
		 yyerror(0, lParse, "Could not allocate temporary memory in median function");
		 free( this->value.data.ptr );
		 break;
	       }

	       for (irow=0; irow<row; irow++) {
		  double *p = mptr;
		  int nelem1 = nelem;

		  while ( nelem1-- ) { 
		    if (*uptr == 0) {
		      *p++ = *dptr;   /* Only advance the dest pointer if we copied */
		    }
		    dptr ++;  /* Advance the source pointer ... */
		    uptr ++;  /* ... and source "undef" pointer */
		  }

		  nelem1 = (p - mptr);  /* Number of accepted data points */
		  if (nelem1 > 0) {
		    this->value.undef[irow] = 0;
		    this->value.data.dblptr[irow] = qselect_median_dbl(mptr, nelem1);
		  } else {
		    this->value.undef[irow] = 1;
		    this->value.data.dblptr[irow] = 0;
		  }

	       }
	       free(mptr);
	    }
	    break;
	 case abs_fct:
	    if( theParams[0]->type==DOUBLE )
	       while( elem-- ) {
		  dval = theParams[0]->value.data.dblptr[elem];
		  this->value.data.dblptr[elem] = (dval>0.0 ? dval : -dval);
		  this->value.undef[elem] = theParams[0]->value.undef[elem];
	       }
	    else
	       while( elem-- ) {
		  ival = theParams[0]->value.data.lngptr[elem];
		  this->value.data.lngptr[elem] = (ival> 0  ? ival : -ival);
		  this->value.undef[elem] = theParams[0]->value.undef[elem];
	       }
	    break;

            /* Special Null-Handling Functions */

	 case nonnull_fct:
	   nelem = theParams[0]->value.nelem;
	   if ( theParams[0]->type==STRING ) nelem = 1;
	   elem = row * nelem;
	   while( row-- ) {
	     int nelem1 = nelem;

	     this->value.undef[row] = 0;        /* Initialize to 0 (defined) */
	     this->value.data.lngptr[row] = 0;
	     while( nelem1-- ) {	
	       elem --;
	       if ( theParams[0]->value.undef[elem] == 0 ) this->value.data.lngptr[row] ++;
	     }
	   }
	   break;
	 case isnull_fct:
	    if( theParams[0]->type==STRING ) elem = row;
	    while( elem-- ) {
	       this->value.data.logptr[elem] = theParams[0]->value.undef[elem];
	       this->value.undef[elem] = 0;
	    }
	    break;
         case defnull_fct:
	    switch( this->type ) {
	    case BOOLEAN:
	       while( row-- ) {
		  nelem = this->value.nelem;
		  while( nelem-- ) {
		     elem--;
		     i=2; while( i-- )
			if( vector[i]>1 ) {
			   pNull[i] = theParams[i]->value.undef[elem];
			   pVals[i].data.log =
			      theParams[i]->value.data.logptr[elem];
			} else if( vector[i] ) {
			   pNull[i] = theParams[i]->value.undef[row];
			   pVals[i].data.log =
			      theParams[i]->value.data.logptr[row];
			}
		     if( pNull[0] ) {
			this->value.undef[elem] = pNull[1];
			this->value.data.logptr[elem] = pVals[1].data.log;
		     } else {
			this->value.undef[elem] = 0;
			this->value.data.logptr[elem] = pVals[0].data.log;
		     }
		  }
	       }
	       break;
	    case LONG:
	       while( row-- ) {
		  nelem = this->value.nelem;
		  while( nelem-- ) {
		     elem--;
		     i=2; while( i-- )
			if( vector[i]>1 ) {
			   pNull[i] = theParams[i]->value.undef[elem];
			   pVals[i].data.lng =
			      theParams[i]->value.data.lngptr[elem];
			} else if( vector[i] ) {
			   pNull[i] = theParams[i]->value.undef[row];
			   pVals[i].data.lng =
			      theParams[i]->value.data.lngptr[row];
			}
		     if( pNull[0] ) {
			this->value.undef[elem] = pNull[1];
			this->value.data.lngptr[elem] = pVals[1].data.lng;
		     } else {
			this->value.undef[elem] = 0;
			this->value.data.lngptr[elem] = pVals[0].data.lng;
		     }
		  }
	       }
	       break;
	    case DOUBLE:
	       while( row-- ) {
		  nelem = this->value.nelem;
		  while( nelem-- ) {
		     elem--;
		     i=2; while( i-- )
			if( vector[i]>1 ) {
			   pNull[i] = theParams[i]->value.undef[elem];
			   pVals[i].data.dbl =
			      theParams[i]->value.data.dblptr[elem];
			} else if( vector[i] ) {
			   pNull[i] = theParams[i]->value.undef[row];
			   pVals[i].data.dbl =
			      theParams[i]->value.data.dblptr[row];
			}
		     if( pNull[0] ) {
			this->value.undef[elem] = pNull[1];
			this->value.data.dblptr[elem] = pVals[1].data.dbl;
		     } else {
			this->value.undef[elem] = 0;
			this->value.data.dblptr[elem] = pVals[0].data.dbl;
		     }
		  }
	       }
	       break;
	    case STRING:
	       while( row-- ) {
		  i=2; while( i-- )
		     if( vector[i] ) {
			pNull[i] = theParams[i]->value.undef[row];
			strcpy(pVals[i].data.str,
			       theParams[i]->value.data.strptr[row]);
		     }
		  if( pNull[0] ) {
		     this->value.undef[row] = pNull[1];
		     strcpy(this->value.data.strptr[row],pVals[1].data.str);
		  } else {
		     this->value.undef[elem] = 0;
		     strcpy(this->value.data.strptr[row],pVals[0].data.str);
		  }
	       }
	    }
	    break;
         case setnull_fct:
	    switch( this->type ) {
	    case LONG:
	      while( elem-- ) {
		if ( theParams[1]->value.data.lng == 
		     theParams[0]->value.data.lngptr[elem] ) {
		  this->value.data.lngptr[elem] = 0;
		  this->value.undef[elem] = 1;
		} else {
		  this->value.data.lngptr[elem] = theParams[0]->value.data.lngptr[elem];
		  this->value.undef[elem] = theParams[0]->value.undef[elem];
		}
	      }
	      break;
	    case DOUBLE:
	      while( elem-- ) {
		if ( theParams[1]->value.data.dbl == 
		     theParams[0]->value.data.dblptr[elem] ) {
		  this->value.data.dblptr[elem] = 0;
		  this->value.undef[elem] = 1;
		} else {
		  this->value.data.dblptr[elem] = theParams[0]->value.data.dblptr[elem];
		  this->value.undef[elem] = theParams[0]->value.undef[elem];
		}
	      }
	      break;
	    }
	    break;

	    /* Math functions with 1 double argument */

	 case sin_fct:
	    while( elem-- )
	       if( !(this->value.undef[elem] = theParams[0]->value.undef[elem]) ) {
		  this->value.data.dblptr[elem] = 
		     sin( theParams[0]->value.data.dblptr[elem] );
	       }
	    break;
	 case cos_fct:
	    while( elem-- )
	       if( !(this->value.undef[elem] = theParams[0]->value.undef[elem]) ) {
		  this->value.data.dblptr[elem] = 
		     cos( theParams[0]->value.data.dblptr[elem] );
	       }
	    break;
	 case tan_fct:
	    while( elem-- )
	       if( !(this->value.undef[elem] = theParams[0]->value.undef[elem]) ) {
		  this->value.data.dblptr[elem] = 
		     tan( theParams[0]->value.data.dblptr[elem] );
	       }
	    break;
	 case asin_fct:
	    while( elem-- )
	       if( !(this->value.undef[elem] = theParams[0]->value.undef[elem]) ) {
		  dval = theParams[0]->value.data.dblptr[elem];
		  if( dval<-1.0 || dval>1.0 ) {
		     this->value.data.dblptr[elem] = 0.0;
		     this->value.undef[elem] = 1;
		  } else
		     this->value.data.dblptr[elem] = asin( dval );
	       }
	    break;
	 case acos_fct:
	    while( elem-- )
	       if( !(this->value.undef[elem] = theParams[0]->value.undef[elem]) ) {
		  dval = theParams[0]->value.data.dblptr[elem];
		  if( dval<-1.0 || dval>1.0 ) {
		     this->value.data.dblptr[elem] = 0.0;
		     this->value.undef[elem] = 1;
		  } else
		     this->value.data.dblptr[elem] = acos( dval );
	       }
	    break;
	 case atan_fct:
	    while( elem-- )
	       if( !(this->value.undef[elem] = theParams[0]->value.undef[elem]) ) {
		  dval = theParams[0]->value.data.dblptr[elem];
		  this->value.data.dblptr[elem] = atan( dval );
	       }
	    break;
	 case sinh_fct:
	    while( elem-- )
	       if( !(this->value.undef[elem] = theParams[0]->value.undef[elem]) ) {
		  this->value.data.dblptr[elem] = 
		     sinh( theParams[0]->value.data.dblptr[elem] );
	       }
	    break;
	 case cosh_fct:
	    while( elem-- )
	       if( !(this->value.undef[elem] = theParams[0]->value.undef[elem]) ) {
		  this->value.data.dblptr[elem] = 
		     cosh( theParams[0]->value.data.dblptr[elem] );
	       }
	    break;
	 case tanh_fct:
	    while( elem-- )
	       if( !(this->value.undef[elem] = theParams[0]->value.undef[elem]) ) {
		  this->value.data.dblptr[elem] = 
		     tanh( theParams[0]->value.data.dblptr[elem] );
	       }
	    break;
	 case exp_fct:
	    while( elem-- )
	       if( !(this->value.undef[elem] = theParams[0]->value.undef[elem]) ) {
		  dval = theParams[0]->value.data.dblptr[elem];
		  this->value.data.dblptr[elem] = exp( dval );
	       }
	    break;
	 case log_fct:
	    while( elem-- )
	       if( !(this->value.undef[elem] = theParams[0]->value.undef[elem]) ) {
		  dval = theParams[0]->value.data.dblptr[elem];
		  if( dval<=0.0 ) {
		     this->value.data.dblptr[elem] = 0.0;
		     this->value.undef[elem] = 1;
		  } else
		     this->value.data.dblptr[elem] = log( dval );
	       }
	    break;
	 case log10_fct:
	    while( elem-- )
	       if( !(this->value.undef[elem] = theParams[0]->value.undef[elem]) ) {
		  dval = theParams[0]->value.data.dblptr[elem];
		  if( dval<=0.0 ) {
		     this->value.data.dblptr[elem] = 0.0;
		     this->value.undef[elem] = 1;
		  } else
		     this->value.data.dblptr[elem] = log10( dval );
	       }
	    break;
	 case sqrt_fct:
	    while( elem-- )
	       if( !(this->value.undef[elem] = theParams[0]->value.undef[elem]) ) {
		  dval = theParams[0]->value.data.dblptr[elem];
		  if( dval<0.0 ) {
		     this->value.data.dblptr[elem] = 0.0;
		     this->value.undef[elem] = 1;
		  } else
		     this->value.data.dblptr[elem] = sqrt( dval );
	       }
	    break;
	 case ceil_fct:
	    while( elem-- )
	       if( !(this->value.undef[elem] = theParams[0]->value.undef[elem]) ) {
		  this->value.data.dblptr[elem] = 
		     ceil( theParams[0]->value.data.dblptr[elem] );
	       }
	    break;
	 case floor_fct:
	    while( elem-- )
	       if( !(this->value.undef[elem] = theParams[0]->value.undef[elem]) ) {
		  this->value.data.dblptr[elem] = 
		     floor( theParams[0]->value.data.dblptr[elem] );
	       }
	    break;
	 case round_fct:
	    while( elem-- )
	       if( !(this->value.undef[elem] = theParams[0]->value.undef[elem]) ) {
		  this->value.data.dblptr[elem] = 
		     floor( theParams[0]->value.data.dblptr[elem] + 0.5);
	       }
	    break;

	    /* Two-argument Trig Functions */
	    
	 case atan2_fct:
	    while( row-- ) {
	       nelem = this->value.nelem;
	       while( nelem-- ) {
		  elem--;
		  i=2; while( i-- )
		     if( vector[i]>1 ) {
			pVals[i].data.dbl =
			   theParams[i]->value.data.dblptr[elem];
			pNull[i] = theParams[i]->value.undef[elem];
		     } else if( vector[i] ) {
			pVals[i].data.dbl =
			   theParams[i]->value.data.dblptr[row];
			pNull[i] = theParams[i]->value.undef[row];
		     }
		  if( !(this->value.undef[elem] = (pNull[0] || pNull[1]) ) )
		     this->value.data.dblptr[elem] =
			atan2( pVals[0].data.dbl, pVals[1].data.dbl );
	       }
	    }
	    break;

	    /* Four-argument ANGSEP Function */
	    
	 case angsep_fct:
	    while( row-- ) {
	       nelem = this->value.nelem;
	       while( nelem-- ) {
		  elem--;
		  i=4; while( i-- )
		     if( vector[i]>1 ) {
			pVals[i].data.dbl =
			   theParams[i]->value.data.dblptr[elem];
			pNull[i] = theParams[i]->value.undef[elem];
		     } else if( vector[i] ) {
			pVals[i].data.dbl =
			   theParams[i]->value.data.dblptr[row];
			pNull[i] = theParams[i]->value.undef[row];
		     }
		  if( !(this->value.undef[elem] = (pNull[0] || pNull[1] ||
						   pNull[2] || pNull[3]) ) )
		     this->value.data.dblptr[elem] =
		       angsep_calc(pVals[0].data.dbl, pVals[1].data.dbl,
				   pVals[2].data.dbl, pVals[3].data.dbl);
	       }
	    }
	    break;



	    /*  Min/Max functions taking 1 or 2 arguments  */

         case min1_fct:
	    elem = row * theParams[0]->value.nelem;
	    if( this->type==LONG ) {
	       long minVal=0;
	       while( row-- ) {
		  valInit = 1;
		  this->value.undef[row] = 1;
		  nelem = theParams[0]->value.nelem;
		  while( nelem-- ) {
		     elem--;
		     if ( !theParams[0]->value.undef[elem] ) {
		       if ( valInit ) {
			 valInit = 0;
			 minVal  = theParams[0]->value.data.lngptr[elem];
		       } else {
			 minVal  = minvalue( minVal,
					     theParams[0]->value.data.lngptr[elem] );
		       }
		       this->value.undef[row] = 0;
		     }
		  }  
		  this->value.data.lngptr[row] = minVal;
	       }		  
	    } else if( this->type==DOUBLE ) {
	       double minVal=0.0;
	       while( row-- ) {
		  valInit = 1;
		  this->value.undef[row] = 1;
		  nelem = theParams[0]->value.nelem;
		  while( nelem-- ) {
		     elem--;
		     if ( !theParams[0]->value.undef[elem] ) {
		       if ( valInit ) {
			 valInit = 0;
			 minVal  = theParams[0]->value.data.dblptr[elem];
		       } else {
			 minVal  = minvalue( minVal,
					     theParams[0]->value.data.dblptr[elem] );
		       }
		       this->value.undef[row] = 0;
		     }
		  }  
		  this->value.data.dblptr[row] = minVal;
	       }		  
	    } else if( this->type==BITSTR ) {
	       char minVal;
	       while( row-- ) {
		  char *sptr1 = theParams[0]->value.data.strptr[row];
		  minVal = '1';
		  while (*sptr1) {
		    if (*sptr1 == '0') minVal = '0';
		    sptr1++;
		  }
		  this->value.data.strptr[row][0] = minVal;
		  this->value.data.strptr[row][1] = 0;     /* Null terminate */
	       }		  
	    }
	    break;
         case min2_fct:
	    if( this->type==LONG ) {
	       while( row-- ) {
		  nelem = this->value.nelem;
		  while( nelem-- ) {
		     elem--;
		     i=2; while( i-- )
			if( vector[i]>1 ) {
			   pVals[i].data.lng =
			      theParams[i]->value.data.lngptr[elem];
			   pNull[i] = theParams[i]->value.undef[elem];
			} else if( vector[i] ) {
			   pVals[i].data.lng =
			      theParams[i]->value.data.lngptr[row];
			   pNull[i] = theParams[i]->value.undef[row];
			}
		     if( pNull[0] && pNull[1] ) {
		       this->value.undef[elem] = 1;
		       this->value.data.lngptr[elem] = 0;
		     } else if (pNull[0]) {
		       this->value.undef[elem] = 0;
		       this->value.data.lngptr[elem] = pVals[1].data.lng;
		     } else if (pNull[1]) {
		       this->value.undef[elem] = 0;
		       this->value.data.lngptr[elem] = pVals[0].data.lng;
		     } else {
		       this->value.undef[elem] = 0;
		       this->value.data.lngptr[elem] =
			 minvalue( pVals[0].data.lng, pVals[1].data.lng );
		     }
		  }
	       }
	    } else if( this->type==DOUBLE ) {
	       while( row-- ) {
		  nelem = this->value.nelem;
		  while( nelem-- ) {
		     elem--;
		     i=2; while( i-- )
			if( vector[i]>1 ) {
			   pVals[i].data.dbl =
			      theParams[i]->value.data.dblptr[elem];
			   pNull[i] = theParams[i]->value.undef[elem];
			} else if( vector[i] ) {
			   pVals[i].data.dbl =
			      theParams[i]->value.data.dblptr[row];
			   pNull[i] = theParams[i]->value.undef[row];
			}
		     if( pNull[0] && pNull[1] ) {
		       this->value.undef[elem] = 1;
		       this->value.data.dblptr[elem] = 0;
		     } else if (pNull[0]) {
		       this->value.undef[elem] = 0;
		       this->value.data.dblptr[elem] = pVals[1].data.dbl;
		     } else if (pNull[1]) {
		       this->value.undef[elem] = 0;
		       this->value.data.dblptr[elem] = pVals[0].data.dbl;
		     } else {
		       this->value.undef[elem] = 0;
		       this->value.data.dblptr[elem] =
			 minvalue( pVals[0].data.dbl, pVals[1].data.dbl );
		     }
		  }
 	       }
	    }
	    break;

         case max1_fct:
	    elem = row * theParams[0]->value.nelem;
	    if( this->type==LONG ) {
	       long maxVal=0;
	       while( row-- ) {
		  valInit = 1;
		  this->value.undef[row] = 1;
		  nelem = theParams[0]->value.nelem;
		  while( nelem-- ) {
		     elem--;
		     if ( !theParams[0]->value.undef[elem] ) {
		       if ( valInit ) {
			 valInit = 0;
			 maxVal  = theParams[0]->value.data.lngptr[elem];
		       } else {
			 maxVal  = maxvalue( maxVal,
					     theParams[0]->value.data.lngptr[elem] );
		       }
		       this->value.undef[row] = 0;
		     }
		  }
		  this->value.data.lngptr[row] = maxVal;
	       }		  
	    } else if( this->type==DOUBLE ) {
	       double maxVal=0.0;
	       while( row-- ) {
		  valInit = 1;
		  this->value.undef[row] = 1;
		  nelem = theParams[0]->value.nelem;
		  while( nelem-- ) {
		     elem--;
		     if ( !theParams[0]->value.undef[elem] ) {
		       if ( valInit ) {
			 valInit = 0;
			 maxVal  = theParams[0]->value.data.dblptr[elem];
		       } else {
			 maxVal  = maxvalue( maxVal,
					     theParams[0]->value.data.dblptr[elem] );
		       }
		       this->value.undef[row] = 0;
		     }
		  }
		  this->value.data.dblptr[row] = maxVal;
	       }		  
	    } else if( this->type==BITSTR ) {
	       char maxVal;
	       while( row-- ) {
		  char *sptr1 = theParams[0]->value.data.strptr[row];
		  maxVal = '0';
		  while (*sptr1) {
		    if (*sptr1 == '1') maxVal = '1';
		    sptr1++;
		  }
		  this->value.data.strptr[row][0] = maxVal;
		  this->value.data.strptr[row][1] = 0;     /* Null terminate */
	       }		  
	    }
	    break;
         case max2_fct:
	    if( this->type==LONG ) {
	       while( row-- ) {
		  nelem = this->value.nelem;
		  while( nelem-- ) {
		     elem--;
		     i=2; while( i-- )
			if( vector[i]>1 ) {
			   pVals[i].data.lng =
			      theParams[i]->value.data.lngptr[elem];
			   pNull[i] = theParams[i]->value.undef[elem];
			} else if( vector[i] ) {
			   pVals[i].data.lng =
			      theParams[i]->value.data.lngptr[row];
			   pNull[i] = theParams[i]->value.undef[row];
			}
		     if( pNull[0] && pNull[1] ) {
		       this->value.undef[elem] = 1;
		       this->value.data.lngptr[elem] = 0;
		     } else if (pNull[0]) {
		       this->value.undef[elem] = 0;
		       this->value.data.lngptr[elem] = pVals[1].data.lng;
		     } else if (pNull[1]) {
		       this->value.undef[elem] = 0;
		       this->value.data.lngptr[elem] = pVals[0].data.lng;
		     } else {
		       this->value.undef[elem] = 0;
		       this->value.data.lngptr[elem] =
			 maxvalue( pVals[0].data.lng, pVals[1].data.lng );
		     }
		  }
	       }
	    } else if( this->type==DOUBLE ) {
	       while( row-- ) {
		  nelem = this->value.nelem;
		  while( nelem-- ) {
		     elem--;
		     i=2; while( i-- )
			if( vector[i]>1 ) {
			   pVals[i].data.dbl =
			      theParams[i]->value.data.dblptr[elem];
			   pNull[i] = theParams[i]->value.undef[elem];
			} else if( vector[i] ) {
			   pVals[i].data.dbl =
			      theParams[i]->value.data.dblptr[row];
			   pNull[i] = theParams[i]->value.undef[row];
			}
		     if( pNull[0] && pNull[1] ) {
		       this->value.undef[elem] = 1;
		       this->value.data.dblptr[elem] = 0;
		     } else if (pNull[0]) {
		       this->value.undef[elem] = 0;
		       this->value.data.dblptr[elem] = pVals[1].data.dbl;
		     } else if (pNull[1]) {
		       this->value.undef[elem] = 0;
		       this->value.data.dblptr[elem] = pVals[0].data.dbl;
		     } else {
		       this->value.undef[elem] = 0;
		       this->value.data.dblptr[elem] =
			 maxvalue( pVals[0].data.dbl, pVals[1].data.dbl );
		     }
		  }
	       }
	    }
	    break;

	    /* Boolean SAO region Functions... scalar or vector dbls */

	 case near_fct:
	    while( row-- ) {
	       nelem = this->value.nelem;
	       while( nelem-- ) {
		  elem--;
		  i=3; while( i-- )
		     if( vector[i]>1 ) {
			pVals[i].data.dbl =
			   theParams[i]->value.data.dblptr[elem];
			pNull[i] = theParams[i]->value.undef[elem];
		     } else if( vector[i] ) {
			pVals[i].data.dbl =
			   theParams[i]->value.data.dblptr[row];
			pNull[i] = theParams[i]->value.undef[row];
		     }
		  if( !(this->value.undef[elem] = (pNull[0] || pNull[1] ||
						   pNull[2]) ) )
		    this->value.data.logptr[elem] =
		      bnear( pVals[0].data.dbl, pVals[1].data.dbl,
			     pVals[2].data.dbl );
	       }
	    }
	    break;

	 case circle_fct:
	    while( row-- ) {
	       nelem = this->value.nelem;
	       while( nelem-- ) {
		  elem--;
		  i=5; while( i-- )
		     if( vector[i]>1 ) {
			pVals[i].data.dbl =
			   theParams[i]->value.data.dblptr[elem];
			pNull[i] = theParams[i]->value.undef[elem];
		     } else if( vector[i] ) {
			pVals[i].data.dbl =
			   theParams[i]->value.data.dblptr[row];
			pNull[i] = theParams[i]->value.undef[row];
		     }
		  if( !(this->value.undef[elem] = (pNull[0] || pNull[1] ||
						   pNull[2] || pNull[3] ||
						   pNull[4]) ) )
		    this->value.data.logptr[elem] =
		     circle( pVals[0].data.dbl, pVals[1].data.dbl,
			     pVals[2].data.dbl, pVals[3].data.dbl,
			     pVals[4].data.dbl );
	       }
	    }
	    break;

	 case box_fct:
	    while( row-- ) {
	       nelem = this->value.nelem;
	       while( nelem-- ) {
		  elem--;
		  i=7; while( i-- )
		     if( vector[i]>1 ) {
			pVals[i].data.dbl =
			   theParams[i]->value.data.dblptr[elem];
			pNull[i] = theParams[i]->value.undef[elem];
		     } else if( vector[i] ) {
			pVals[i].data.dbl =
			   theParams[i]->value.data.dblptr[row];
			pNull[i] = theParams[i]->value.undef[row];
		     }
		  if( !(this->value.undef[elem] = (pNull[0] || pNull[1] ||
						   pNull[2] || pNull[3] ||
						   pNull[4] || pNull[5] ||
						   pNull[6] ) ) )
		    this->value.data.logptr[elem] =
		     saobox( pVals[0].data.dbl, pVals[1].data.dbl,
			     pVals[2].data.dbl, pVals[3].data.dbl,
			     pVals[4].data.dbl, pVals[5].data.dbl,
			     pVals[6].data.dbl );	
	       }
	    }
	    break;

	 case elps_fct:
	    while( row-- ) {
	       nelem = this->value.nelem;
	       while( nelem-- ) {
		  elem--;
		  i=7; while( i-- )
		     if( vector[i]>1 ) {
			pVals[i].data.dbl =
			   theParams[i]->value.data.dblptr[elem];
			pNull[i] = theParams[i]->value.undef[elem];
		     } else if( vector[i] ) {
			pVals[i].data.dbl =
			   theParams[i]->value.data.dblptr[row];
			pNull[i] = theParams[i]->value.undef[row];
		     }
		  if( !(this->value.undef[elem] = (pNull[0] || pNull[1] ||
						   pNull[2] || pNull[3] ||
						   pNull[4] || pNull[5] ||
						   pNull[6] ) ) )
		    this->value.data.logptr[elem] =
		     ellipse( pVals[0].data.dbl, pVals[1].data.dbl,
			      pVals[2].data.dbl, pVals[3].data.dbl,
			      pVals[4].data.dbl, pVals[5].data.dbl,
			      pVals[6].data.dbl );
	       }
	    }
	    break;

            /* C Conditional expression:  bool ? expr : expr */

         case ifthenelse_fct:
            switch( this->type ) {
            case BOOLEAN:
	       while( row-- ) {
		  nelem = this->value.nelem;
		  while( nelem-- ) {
		     elem--;
                     if( vector[2]>1 ) {
                        pVals[2].data.log =
                           theParams[2]->value.data.logptr[elem];
                        pNull[2] = theParams[2]->value.undef[elem];
                     } else if( vector[2] ) {
                        pVals[2].data.log =
                           theParams[2]->value.data.logptr[row];
                        pNull[2] = theParams[2]->value.undef[row];
                     }
		     i=2; while( i-- )
			if( vector[i]>1 ) {
			   pVals[i].data.log =
			      theParams[i]->value.data.logptr[elem];
			   pNull[i] = theParams[i]->value.undef[elem];
			} else if( vector[i] ) {
			   pVals[i].data.log =
			      theParams[i]->value.data.logptr[row];
			   pNull[i] = theParams[i]->value.undef[row];
			}
		     if( !(this->value.undef[elem] = pNull[2]) ) {
                        if( pVals[2].data.log ) {
                           this->value.data.logptr[elem] = pVals[0].data.log;
                           this->value.undef[elem]       = pNull[0];
                        } else {
                           this->value.data.logptr[elem] = pVals[1].data.log;
                           this->value.undef[elem]       = pNull[1];
                        }
                     }
		  }
	       }
               break;
            case LONG:
	       while( row-- ) {
		  nelem = this->value.nelem;
		  while( nelem-- ) {
		     elem--;
                     if( vector[2]>1 ) {
                        pVals[2].data.log =
                           theParams[2]->value.data.logptr[elem];
                        pNull[2] = theParams[2]->value.undef[elem];
                     } else if( vector[2] ) {
                        pVals[2].data.log =
                           theParams[2]->value.data.logptr[row];
                        pNull[2] = theParams[2]->value.undef[row];
                     }
		     i=2; while( i-- )
			if( vector[i]>1 ) {
			   pVals[i].data.lng =
			      theParams[i]->value.data.lngptr[elem];
			   pNull[i] = theParams[i]->value.undef[elem];
			} else if( vector[i] ) {
			   pVals[i].data.lng =
			      theParams[i]->value.data.lngptr[row];
			   pNull[i] = theParams[i]->value.undef[row];
			}
		     if( !(this->value.undef[elem] = pNull[2]) ) {
                        if( pVals[2].data.log ) {
                           this->value.data.lngptr[elem] = pVals[0].data.lng;
                           this->value.undef[elem]       = pNull[0];
                        } else {
                           this->value.data.lngptr[elem] = pVals[1].data.lng;
                           this->value.undef[elem]       = pNull[1];
                        }
                     }
		  }
	       }
               break;
            case DOUBLE:
	       while( row-- ) {
		  nelem = this->value.nelem;
		  while( nelem-- ) {
		     elem--;
                     if( vector[2]>1 ) {
                        pVals[2].data.log =
                           theParams[2]->value.data.logptr[elem];
                        pNull[2] = theParams[2]->value.undef[elem];
                     } else if( vector[2] ) {
                        pVals[2].data.log =
                           theParams[2]->value.data.logptr[row];
                        pNull[2] = theParams[2]->value.undef[row];
                     }
		     i=2; while( i-- )
			if( vector[i]>1 ) {
			   pVals[i].data.dbl =
			      theParams[i]->value.data.dblptr[elem];
			   pNull[i] = theParams[i]->value.undef[elem];
			} else if( vector[i] ) {
			   pVals[i].data.dbl =
			      theParams[i]->value.data.dblptr[row];
			   pNull[i] = theParams[i]->value.undef[row];
			}
		     if( !(this->value.undef[elem] = pNull[2]) ) {
                        if( pVals[2].data.log ) {
                           this->value.data.dblptr[elem] = pVals[0].data.dbl;
                           this->value.undef[elem]       = pNull[0];
                        } else {
                           this->value.data.dblptr[elem] = pVals[1].data.dbl;
                           this->value.undef[elem]       = pNull[1];
                        }
                     }
		  }
	       }
               break;
            case STRING:
	       while( row-- ) {
                  if( vector[2] ) {
                     pVals[2].data.log = theParams[2]->value.data.logptr[row];
                     pNull[2] = theParams[2]->value.undef[row];
                  }
                  i=2; while( i-- )
                     if( vector[i] ) {
                        strcpy( pVals[i].data.str,
                                theParams[i]->value.data.strptr[row] );
                        pNull[i] = theParams[i]->value.undef[row];
                     }
                  if( !(this->value.undef[row] = pNull[2]) ) {
                     if( pVals[2].data.log ) {
                        strcpy( this->value.data.strptr[row],
                                pVals[0].data.str );
                        this->value.undef[row]       = pNull[0];
                     } else {
                        strcpy( this->value.data.strptr[row],
                                pVals[1].data.str );
                        this->value.undef[row]       = pNull[1];
                     }
                  } else {
                     this->value.data.strptr[row][0] = '\0';
                  }
	       }
               break;

            }
            break;

	    /* String functions */
            case strmid_fct:
	      {
		int strconst = theParams[0]->operation == CONST_OP;
		int posconst = theParams[1]->operation == CONST_OP;
		int lenconst = theParams[2]->operation == CONST_OP;
		int dest_len = this->value.nelem;
		int src_len  = theParams[0]->value.nelem;

		while (row--) {
		  int pos;
		  int len;
		  char *str;
		  int undef = 0;

		  if (posconst) {
		    pos = theParams[1]->value.data.lng;
		  } else {
		    pos = theParams[1]->value.data.lngptr[row];
		    if (theParams[1]->value.undef[row]) undef = 1;
		  }
		  if (strconst) {
		    str = theParams[0]->value.data.str;
		    if (src_len == 0) src_len = strlen(str);
		  } else {
		    str = theParams[0]->value.data.strptr[row];
		    if (theParams[0]->value.undef[row]) undef = 1;
		  }
		  if (lenconst) {
		    len = dest_len;
		  } else {
		    len = theParams[2]->value.data.lngptr[row];
		    if (theParams[2]->value.undef[row]) undef = 1;
		  }
		  this->value.data.strptr[row][0] = '\0';
		  if (pos == 0) undef = 1;
		  if (! undef ) {
		    if (cstrmid(lParse,
				this->value.data.strptr[row], len,
				str, src_len, pos) < 0) break;
		  }
		  this->value.undef[row] = undef;
		}
	      }		      
	      break;

	    /* String functions */
            case strpos_fct:
	      {
		int const1 = theParams[0]->operation == CONST_OP;
		int const2 = theParams[1]->operation == CONST_OP;

		while (row--) {
		  char *str1, *str2;
		  int undef = 0;

		  if (const1) {
		    str1 = theParams[0]->value.data.str;
		  } else {
		    str1 = theParams[0]->value.data.strptr[row];
		    if (theParams[0]->value.undef[row]) undef = 1;
		  }
		  if (const2) {
		    str2 = theParams[1]->value.data.str;
		  } else {
		    str2 = theParams[1]->value.data.strptr[row];
		    if (theParams[1]->value.undef[row]) undef = 1;
		  }
		  this->value.data.lngptr[row] = 0;
		  if (! undef ) {
		    char *res = strstr(str1, str2);
		    if (res == NULL) {
		      undef = 1;
		      this->value.data.lngptr[row] = 0; 
		    } else {
		      this->value.data.lngptr[row] = (res - str1) + 1;
		    }
		  }
		  this->value.undef[row] = undef;
		}
	      }
	      break;

		    
	 } /* End switch(this->operation) */
      } /* End if (!lParse->status) */
   } /* End non-constant operations */

   i = this->nSubNodes;
   while( i-- ) {
      if( theParams[i]->operation>0 ) {
	 /*  Currently only numeric params allowed  */
	 free( theParams[i]->value.data.ptr );
      }
   }
}

static void Do_Deref( ParseData *lParse, Node *this )
{
   Node *theVar, *theDims[MAXDIMS];
   int  isConst[MAXDIMS], allConst;
   long dimVals[MAXDIMS];
   int  i, nDims;
   long row, elem, dsize;

   theVar = lParse->Nodes + this->SubNodes[0];

   i = nDims = this->nSubNodes-1;
   allConst = 1;
   while( i-- ) {
      theDims[i] = lParse->Nodes + this->SubNodes[i+1];
      isConst[i] = ( theDims[i]->operation==CONST_OP );
      if( isConst[i] )
	 dimVals[i] = theDims[i]->value.data.lng;
      else
	 allConst = 0;
   }

   if( this->type==DOUBLE ) {
      dsize = sizeof( double );
   } else if( this->type==LONG ) {
      dsize = sizeof( long );
   } else if( this->type==BOOLEAN ) {
      dsize = sizeof( char );
   } else
      dsize = 0;

   Allocate_Ptrs( lParse, this );

   if( !lParse->status ) {

      if( allConst && theVar->value.naxis==nDims ) {

	 /* Dereference completely using constant indices */

	 elem = 0;
	 i    = nDims;
	 while( i-- ) {
	    if( dimVals[i]<1 || dimVals[i]>theVar->value.naxes[i] ) break;
	    elem = theVar->value.naxes[i]*elem + dimVals[i]-1;
	 }
	 if( i<0 ) {
	    for( row=0; row<lParse->nRows; row++ ) {
	       if( this->type==STRING )
		 this->value.undef[row] = theVar->value.undef[row];
	       else if( this->type==BITSTR ) 
		 this->value.undef;  /* Dummy - BITSTRs do not have undefs */
	       else 
		 this->value.undef[row] = theVar->value.undef[elem];

	       if( this->type==DOUBLE )
		  this->value.data.dblptr[row] = 
		     theVar->value.data.dblptr[elem];
	       else if( this->type==LONG )
		  this->value.data.lngptr[row] = 
		     theVar->value.data.lngptr[elem];
	       else if( this->type==BOOLEAN )
		  this->value.data.logptr[row] = 
		     theVar->value.data.logptr[elem];
	       else {
		 /* XXX Note, the below expression uses knowledge of
                    the layout of the string format, namely (nelem+1)
                    characters per string, followed by (nelem+1)
                    "undef" values. */
		  this->value.data.strptr[row][0] = 
		     theVar->value.data.strptr[0][elem+row];
		  this->value.data.strptr[row][1] = 0;  /* Null terminate */
	       }
	       elem += theVar->value.nelem;
	    }
	 } else {
	    yyerror(0, lParse, "Index out of range");
	    free( this->value.data.ptr );
	 }
	 
      } else if( allConst && nDims==1 ) {
	 
	 /* Reduce dimensions by 1, using a constant index */
	 
	 if( dimVals[0] < 1 ||
	     dimVals[0] > theVar->value.naxes[ theVar->value.naxis-1 ] ) {
	    yyerror(0, lParse, "Index out of range");
	    free( this->value.data.ptr );
	 } else if ( this->type == BITSTR || this->type == STRING ) {
	    elem = this->value.nelem * (dimVals[0]-1);
	    for( row=0; row<lParse->nRows; row++ ) {
	      if (this->value.undef) 
		this->value.undef[row] = theVar->value.undef[row];
	      memcpy( (char*)this->value.data.strptr[0]
		      + row*sizeof(char)*(this->value.nelem+1),
		      (char*)theVar->value.data.strptr[0] + elem*sizeof(char),
		      this->value.nelem * sizeof(char) );
	      /* Null terminate */
	      this->value.data.strptr[row][this->value.nelem] = 0;
	      elem += theVar->value.nelem+1;
	    }	       
	 } else {
	    elem = this->value.nelem * (dimVals[0]-1);
	    for( row=0; row<lParse->nRows; row++ ) {
	       memcpy( this->value.undef + row*this->value.nelem,
		       theVar->value.undef + elem,
		       this->value.nelem * sizeof(char) );
	       memcpy( (char*)this->value.data.ptr
		       + row*dsize*this->value.nelem,
		       (char*)theVar->value.data.ptr + elem*dsize,
		       this->value.nelem * dsize );
	       elem += theVar->value.nelem;
	    }	       
	 }
      
      } else if( theVar->value.naxis==nDims ) {

	 /* Dereference completely using an expression for the indices */

	 for( row=0; row<lParse->nRows; row++ ) {

	    for( i=0; i<nDims; i++ ) {
	       if( !isConst[i] ) {
		  if( theDims[i]->value.undef[row] ) {
		     yyerror(0, lParse, "Null encountered as vector index");
		     free( this->value.data.ptr );
		     break;
		  } else
		     dimVals[i] = theDims[i]->value.data.lngptr[row];
	       }
	    }
	    if( lParse->status ) break;

	    elem = 0;
	    i    = nDims;
	    while( i-- ) {
	       if( dimVals[i]<1 || dimVals[i]>theVar->value.naxes[i] ) break;
	       elem = theVar->value.naxes[i]*elem + dimVals[i]-1;
	    }
	    if( i<0 ) {
	       elem += row*theVar->value.nelem;

	       if( this->type==STRING )
		 this->value.undef[row] = theVar->value.undef[row];
	       else if( this->type==BITSTR ) 
		 this->value.undef;  /* Dummy - BITSTRs do not have undefs */
	       else 
		 this->value.undef[row] = theVar->value.undef[elem];

	       if( this->type==DOUBLE )
		  this->value.data.dblptr[row] = 
		     theVar->value.data.dblptr[elem];
	       else if( this->type==LONG )
		  this->value.data.lngptr[row] = 
		     theVar->value.data.lngptr[elem];
	       else if( this->type==BOOLEAN )
		  this->value.data.logptr[row] = 
		     theVar->value.data.logptr[elem];
	       else {
		 /* XXX Note, the below expression uses knowledge of
                    the layout of the string format, namely (nelem+1)
                    characters per string, followed by (nelem+1)
                    "undef" values. */
		  this->value.data.strptr[row][0] = 
		     theVar->value.data.strptr[0][elem+row];
		  this->value.data.strptr[row][1] = 0;  /* Null terminate */
	       }
	    } else {
	       yyerror(0, lParse, "Index out of range");
	       free( this->value.data.ptr );
	    }
	 }

      } else {

	 /* Reduce dimensions by 1, using a nonconstant expression */

	 for( row=0; row<lParse->nRows; row++ ) {

	    /* Index cannot be a constant */

	    if( theDims[0]->value.undef[row] ) {
	       yyerror(0, lParse, "Null encountered as vector index");
	       free( this->value.data.ptr );
	       break;
	    } else
	       dimVals[0] = theDims[0]->value.data.lngptr[row];

	    if( dimVals[0] < 1 ||
		dimVals[0] > theVar->value.naxes[ theVar->value.naxis-1 ] ) {
	       yyerror(0, lParse, "Index out of range");
	       free( this->value.data.ptr );
	    } else if ( this->type == BITSTR || this->type == STRING ) {
	      elem = this->value.nelem * (dimVals[0]-1);
	      elem += row*(theVar->value.nelem+1);
	      if (this->value.undef) 
		this->value.undef[row] = theVar->value.undef[row];
	      memcpy( (char*)this->value.data.strptr[0]
		      + row*sizeof(char)*(this->value.nelem+1),
		      (char*)theVar->value.data.strptr[0] + elem*sizeof(char),
		      this->value.nelem * sizeof(char) );
	      /* Null terminate */
	      this->value.data.strptr[row][this->value.nelem] = 0;
	    } else {
	       elem  = this->value.nelem * (dimVals[0]-1);
	       elem += row*theVar->value.nelem;
	       memcpy( this->value.undef + row*this->value.nelem,
		       theVar->value.undef + elem,
		       this->value.nelem * sizeof(char) );
	       memcpy( (char*)this->value.data.ptr
		       + row*dsize*this->value.nelem,
		       (char*)theVar->value.data.ptr + elem*dsize,
		       this->value.nelem * dsize );
	    }
	 }
      }
   }

   if( theVar->operation>0 ) {
     if (theVar->type == STRING || theVar->type == BITSTR) 
       free(theVar->value.data.strptr[0] );
     else 
       free( theVar->value.data.ptr );
   }
   for( i=0; i<nDims; i++ )
      if( theDims[i]->operation>0 ) {
	 free( theDims[i]->value.data.ptr );
      }
}

static void Do_GTI( ParseData *lParse, Node *this )
{
   Node *theExpr, *theTimes;
   double *start, *stop, *times;
   long elem, nGTI, gti;
   int ordered;
   int dorow = (this->operation == gtifind_fct);

   theTimes = lParse->Nodes + this->SubNodes[0];
   theExpr  = lParse->Nodes + this->SubNodes[1];

   nGTI    = theTimes->value.nelem;
   start   = theTimes->value.data.dblptr;
   stop    = theTimes->value.data.dblptr + nGTI;
   ordered = theTimes->type;

   if( theExpr->operation==CONST_OP ) {
      gti = Search_GTI( theExpr->value.data.dbl, nGTI, start, stop, ordered, 0 );
      if (dorow) {
	this->value.data.lng = (gti >= 0) ? (gti+1) : -1;
      } else {
	this->value.data.log = (gti>=0);
      }
      this->operation      = CONST_OP;

   } else {

      Allocate_Ptrs( lParse, this );

      times = theExpr->value.data.dblptr;
      if( !lParse->status ) {

	 elem = lParse->nRows * this->value.nelem;
	 if( nGTI ) {
	    gti = -1;
	    while( elem-- ) {
	       if( (this->value.undef[elem] = theExpr->value.undef[elem]) )
		  continue;

            /*  Before searching entire GTI, check the GTI found last time  */
	       if( gti<0 || times[elem]<start[gti] || times[elem]>stop[gti] ) {
		 gti = Search_GTI( times[elem], nGTI, start, stop, ordered, 0 );
	       }
	       if (dorow) {
		 this->value.data.lngptr[elem] = ( gti >= 0 ) ? (gti + 1) : (-1);
		 this->value.undef[elem]  = ( gti >= 0 ) ? 0 : 1;
	       } else {
		 this->value.data.logptr[elem] = ( gti>=0 );
	       }
	    }
	 } else { /* nGTI == 0 */

	   if (dorow) { /* no good times so all values are undef */
	     while( elem-- ) {
	       this->value.undef[elem]       = 1;
	     }
	   } else {    /* no good times so all logicals are 0 */
	     while( elem-- ) {
	       this->value.data.logptr[elem] = 0;
	       this->value.undef[elem]       = 0;
	     }
	   }
	   
	 }
      }
   }

   if( theExpr->operation>0 )
      free( theExpr->value.data.ptr );
}

static void Do_GTI_Over( ParseData *lParse, Node *this )
{
   Node *theTimes, *theStart, *theStop;
   double *gtiStart, *gtiStop;
   double *evtStart, *evtStop;
   long elem, nGTI, gti, nextGTI;
   int ordered;

   theTimes = lParse->Nodes + this->SubNodes[0]; /* GTI times */
   theStop  = lParse->Nodes + this->SubNodes[2]; /* User start time */
   theStart = lParse->Nodes + this->SubNodes[1]; /* User stop time */

   nGTI     = theTimes->value.nelem;
   gtiStart = theTimes->value.data.dblptr;        /* GTI start */
   gtiStop  = theTimes->value.data.dblptr + nGTI; /* GTI stop */

   if( theStart->operation==CONST_OP && theStop->operation==CONST_OP) {

      this->value.data.dbl = 
	(GTI_Over( theStart->value.data.dbl, theStop->value.data.dbl,
		   nGTI, gtiStart, gtiStop, &gti));
      this->operation      = CONST_OP;

   } else {
      char undefStart = 0, undefStop = 0; /* Input values are undef? */
      double uStart, uStop;       /* User start/stop values */
      if (theStart->operation==CONST_OP) uStart = theStart->value.data.dbl;
      if (theStop ->operation==CONST_OP) uStop  = theStop ->value.data.dbl;

      Allocate_Ptrs( lParse, this );

      evtStart = theStart->value.data.dblptr;
      evtStop  = theStop ->value.data.dblptr;
      if( !lParse->status ) {

	 elem = lParse->nRows * this->value.nelem;
	 if( nGTI ) {
	    double toverlap = 0.0;
	    gti = -1;
	    while( elem-- ) {
	      if (theStart->operation!=CONST_OP) {
		undefStart = theStart->value.undef[elem];
		uStart     = evtStart[elem];
	      }
	      if (theStop->operation!=CONST_OP) {
		undefStop  = theStop ->value.undef[elem];
		uStop      = evtStop[elem];
	      }
	      /* This works because at least one of the values is not const */
	      if( (this->value.undef[elem] = (undefStart||undefStop)) )
		  continue;

            /*  Before searching entire GTI, check the GTI found last time  */
	       if( gti<0 || 
		   uStart<gtiStart[gti] || uStart>gtiStop[gti] ||
		   uStop <gtiStart[gti] || uStop >gtiStop[gti]) {
		 /* Nope, need to recalculate */
		 toverlap = GTI_Over(uStart, uStop, 
				     nGTI, gtiStart, gtiStop, 
				     &gti);
	       } else {
		 /* We are in same GTI, the overlap is just stop-start of user range */
		 toverlap = (uStop-uStart);
	       }

	       /* This works because at least one of the values is not const */
	       this->value.data.dblptr[elem] = toverlap;
	    }
	 } else
	    /* nGTI == 0; there is no overlap so set all values to 0.0 */
	    while( elem-- ) {
	       this->value.data.dblptr[elem] = 0.0;
	       this->value.undef[elem]       = 0;
	    }
      }
   }

   if( theStart->operation>0 ) {
     free( theStart->value.data.ptr );
   }
   if( theStop->operation>0 ) {
     free( theStop->value.data.ptr );
   }
}

static double GTI_Over(double evtStart, double evtStop,
		       long nGTI, double *start, double *stop,
		       long *gtiout)
{
  long gti1, gti2, nextGTI1, nextGTI2;
  long gti, nMax;
  double overlap = 0.0;

  *gtiout = -1L;
  /* Zero or negative bin size */
  if (evtStop <= evtStart) return 0.0;

  /* Locate adjacent GTIs for evtStart and evtStop */
  gti1 = Search_GTI(evtStart, nGTI, start, stop, 1, &nextGTI1);
  gti2 = Search_GTI(evtStop,  nGTI, start, stop, 1, &nextGTI2);

  /* evtStart is in gti1, we return that for future processing */
  if (gti1 >= 0) *gtiout = gti1;

  /* Both evtStart/evtStop are beyond the last GTI */
  if (nextGTI1 < 0 && nextGTI2 < 0) return 0.0;

  /* Both evtStart/evtStop are in the same gap between GTIs */
  if (gti1 < 0 && gti2 < 0 && nextGTI1 == nextGTI2) return 0.0;

  /* Both evtStart/evtStop are in the same GTI */
  if (gti1 >= 0 && gti1 == gti2) return (evtStop-evtStart);

  /* Count through the remaining GTIs; there will be at least one */
  /* The largest GTI to consider is either nextGTI2-1, if it exists,
     or nGTI-1 */
  if (nextGTI2 < 0) nMax = nGTI-1;
  else if (gti2 >= 0) nMax = nextGTI2;
  else nMax = nextGTI2-1;
  for (gti = nextGTI1; gti <= nMax; gti++) {
    double starti = start[gti], stopi = stop[gti];
    /* Trim the GTI by actual evtStart/Stop times */
    if (evtStart > starti) starti = evtStart;
    if (evtStop  < stopi ) stopi  = evtStop;
    overlap += (stopi - starti);
  }
    
  return overlap;
}

/*
 * Search_GTI - search GTI for requested evtTime
 * 
 * double evtTime - requested event time
 * long nGTI - number of entries in start[] and stop[]
 * double start[], stop[] - start and stop of each GTI
 * int ordered - set to 1 if time-ordered
 * long *nextGTI0 - upon return, *nextGTI0 is either
 *                   the GTI evtTime is inside
 *                   the next GTI if evtTime is not inside
 *                   -1L if there is no next GTI
 *                   not set if nextGTI0 is a null pointer
 *
 * NOTE: for *nextGTI to be well-defined, the GTI must
 *   be ordered.  This is true when called by Do_GTI.
 *
 * RETURNS: gti index that evtTime is located inside, or -1L
 */
static long Search_GTI( double evtTime, long nGTI, double *start,
			double *stop, int ordered, long *nextGTI0 )
{
   long gti, nextGTI = -1L, step;
                             
   if( ordered && nGTI>15 ) { /*  If time-ordered and lots of GTIs,   */
                              /*  use "FAST" Binary search algorithm  */
      if( evtTime>=start[0] && evtTime<=stop[nGTI-1] ) {
	 gti = step = (nGTI >> 1);
	 while(1) {
	    if( step>1L ) step >>= 1;
	    
	    if( evtTime>stop[gti] ) {
	       if( evtTime>=start[gti+1] )
		  gti += step;
	       else {
		  nextGTI = gti+1;
		  gti = -1L;
		  break;
	       }
	    } else if( evtTime<start[gti] ) {
	       if( evtTime<=stop[gti-1] )
		  gti -= step;
	       else {
		  nextGTI = gti;
		  gti = -1L;
		  break;
	       }
	    } else {
	       nextGTI = gti;
	       break;
	    }
	 }
      } else {
	 if (start[0] > evtTime) nextGTI = 0;
	 gti = -1L;
      }
      
   } else { /*  Use "SLOW" linear search.  Not required to be 
	        ordered, so we have to search the whole table
		no matter what.
	    */
      gti = nGTI;
      while( gti-- ) {
	if( stop[gti] >= evtTime ) nextGTI = gti;
	if( evtTime>=start[gti] && evtTime<=stop[gti] )
	    break;
      }
   }

   if (nextGTI >= nGTI) nextGTI = -1;
   if (nextGTI0) *nextGTI0 = nextGTI;

   return( gti );
}

static void Do_REG( ParseData *lParse, Node *this )
{
   Node *theRegion, *theX, *theY;
   double Xval=0.0, Yval=0.0;
   char   Xnull=0, Ynull=0;
   int    Xvector, Yvector;
   long   nelem, elem, rows;

   theRegion = lParse->Nodes + this->SubNodes[0];
   theX      = lParse->Nodes + this->SubNodes[1];
   theY      = lParse->Nodes + this->SubNodes[2];

   Xvector = ( theX->operation!=CONST_OP );
   if( Xvector )
      Xvector = theX->value.nelem;
   else {
      Xval  = theX->value.data.dbl;
   }

   Yvector = ( theY->operation!=CONST_OP );
   if( Yvector )
      Yvector = theY->value.nelem;
   else {
      Yval  = theY->value.data.dbl;
   } 

   if( !Xvector && !Yvector ) {

      this->value.data.log =
	 ( fits_in_region( Xval, Yval, (SAORegion *)theRegion->value.data.ptr )
	   != 0 );
      this->operation      = CONST_OP;

   } else {

      Allocate_Ptrs( lParse, this );

      if( !lParse->status ) {

	 rows  = lParse->nRows;
	 nelem = this->value.nelem;
	 elem  = rows*nelem;

	 while( rows-- ) {
	    while( nelem-- ) {
	       elem--;

	       if( Xvector>1 ) {
		  Xval  = theX->value.data.dblptr[elem];
		  Xnull = theX->value.undef[elem];
	       } else if( Xvector ) {
		  Xval  = theX->value.data.dblptr[rows];
		  Xnull = theX->value.undef[rows];
	       }

	       if( Yvector>1 ) {
		  Yval  = theY->value.data.dblptr[elem];
		  Ynull = theY->value.undef[elem];
	       } else if( Yvector ) {
		  Yval  = theY->value.data.dblptr[rows];
		  Ynull = theY->value.undef[rows];
	       }

	       this->value.undef[elem] = ( Xnull || Ynull );
	       if( this->value.undef[elem] )
		  continue;

	       this->value.data.logptr[elem] = 
		  ( fits_in_region( Xval, Yval,
				    (SAORegion *)theRegion->value.data.ptr )
		    != 0 );
	    }
	    nelem = this->value.nelem;
	 }
      }
   }

   if( theX->operation>0 )
      free( theX->value.data.ptr );
   if( theY->operation>0 )
      free( theY->value.data.ptr );
}

static void Do_Vector( ParseData *lParse, Node *this )
{
   Node *that;
   long row, elem, idx, jdx, offset=0;
   int node;

   Allocate_Ptrs( lParse, this );

   if( !lParse->status ) {

      for( node=0; node<this->nSubNodes; node++ ) {

	 that = lParse->Nodes + this->SubNodes[node];

	 if( that->operation == CONST_OP ) {

	    idx = lParse->nRows*this->value.nelem + offset;
	    while( (idx-=this->value.nelem)>=0 ) {
	       
	       this->value.undef[idx] = 0;

	       switch( this->type ) {
	       case BOOLEAN:
		  this->value.data.logptr[idx] = that->value.data.log;
		  break;
	       case LONG:
		  this->value.data.lngptr[idx] = that->value.data.lng;
		  break;
	       case DOUBLE:
		  this->value.data.dblptr[idx] = that->value.data.dbl;
		  break;
	       }
	    }
	    
	 } else {
	       
	    row  = lParse->nRows;
	    idx  = row * that->value.nelem;
	    while( row-- ) {
	       elem = that->value.nelem;
	       jdx = row*this->value.nelem + offset;
	       while( elem-- ) {
		  this->value.undef[jdx+elem] =
		     that->value.undef[--idx];

		  switch( this->type ) {
		  case BOOLEAN:
		     this->value.data.logptr[jdx+elem] =
			that->value.data.logptr[idx];
		     break;
		  case LONG:
		     this->value.data.lngptr[jdx+elem] =
			that->value.data.lngptr[idx];
		     break;
		  case DOUBLE:
		     this->value.data.dblptr[jdx+elem] =
			that->value.data.dblptr[idx];
		     break;
		  }
	       }
	    }
	 }
	 offset += that->value.nelem;
      }

   }

   for( node=0; node < this->nSubNodes; node++ )
     if( OPER(this->SubNodes[node])>0 )
       free( lParse->Nodes[this->SubNodes[node]].value.data.ptr );
}

static void Do_Array( ParseData *lParse, Node *this )
{
   Node *that;
   long row, elem, idx, jdx, offset=0;
   int node;

   Allocate_Ptrs( lParse, this );

   if( !lParse->status ) {

     /* This is the item to be replicated */
     that = lParse->Nodes + this->SubNodes[0];

     if( that->operation == CONST_OP ) {

       idx = lParse->nRows*this->value.nelem + offset;
       while( (idx--)>=0 ) {
	       
	 this->value.undef[idx] = 0;

	 switch( this->type ) {
	 case BOOLEAN:
	   this->value.data.logptr[idx] = that->value.data.log;
	   break;
	 case LONG:
	   this->value.data.lngptr[idx] = that->value.data.lng;
	   break;
	 case DOUBLE:
	   this->value.data.dblptr[idx] = that->value.data.dbl;
	   break;
	 }
       }
       
     } else {
       
       row  = lParse->nRows;
       idx  = row * this->value.nelem - 1;
       while( row-- ) {
	 elem = this->value.nelem;
	 while( elem-- ) {
	   this->value.undef[idx] = that->value.undef[row];

	   switch( this->type ) {
	   case BOOLEAN:
	     this->value.data.logptr[idx] = that->value.data.logptr[row];
	     break;
	   case LONG:
	     this->value.data.lngptr[idx] = that->value.data.lngptr[row];
	     break;
	   case DOUBLE:
	     this->value.data.dblptr[idx] = that->value.data.dblptr[row];
	     break;
	   }
	   idx--;
	 }
       }

     } /* not constant */

     if( OPER(this->SubNodes[0])>0 )
       free( lParse->Nodes[this->SubNodes[0]].value.data.ptr );

   }

}

/*****************************************************************************/
/*  Utility routines which perform the calculations on bits and SAO regions  */
/*****************************************************************************/

static char bitlgte(char *bits1, int oper, char *bits2)
{
 int val1, val2, nextbit;
 char result;
 int i, l1, l2, length, ldiff;
 char *stream=0;
 char chr1, chr2;

 l1 = strlen(bits1);
 l2 = strlen(bits2);
 length = (l1 > l2) ? l1 : l2;
 stream = (char *)malloc(sizeof(char)*(length+1));
 if (l1 < l2)
   {
    ldiff = l2 - l1;
    i=0;
    while( ldiff-- ) stream[i++] = '0';
    while( l1--    ) stream[i++] = *(bits1++);
    stream[i] = '\0';
    bits1 = stream;
   }
 else if (l2 < l1)
   {
    ldiff = l1 - l2;
    i=0;
    while( ldiff-- ) stream[i++] = '0';
    while( l2--    ) stream[i++] = *(bits2++);
    stream[i] = '\0';
    bits2 = stream;
   }

 val1 = val2 = 0;
 nextbit = 1;

 while( length-- )
    {
     chr1 = bits1[length];
     chr2 = bits2[length];
     if ((chr1 != 'x')&&(chr1 != 'X')&&(chr2 != 'x')&&(chr2 != 'X'))
       {
        if (chr1 == '1') val1 += nextbit;
        if (chr2 == '1') val2 += nextbit;
        nextbit *= 2;
       }
    }
 result = 0;
 switch (oper)
       {
        case LT:
             if (val1 < val2) result = 1;
             break;
        case LTE:
             if (val1 <= val2) result = 1;
             break;
        case GT:
             if (val1 > val2) result = 1;
             break;
        case GTE:
             if (val1 >= val2) result = 1;
             break;
       }
 free(stream);
 return (result);
}

static void bitand(char *result,char *bitstrm1,char *bitstrm2)
{
 int i, l1, l2, ldiff, largestStream;
 char *stream=0;
 char chr1, chr2;

 l1 = strlen(bitstrm1);
 l2 = strlen(bitstrm2);
 largestStream = (l1 > l2) ? l1 : l2;
 stream = (char *)malloc(sizeof(char)*(largestStream+1));
 if (l1 < l2)
   {
    ldiff = l2 - l1;
    i=0;
    while( ldiff-- ) stream[i++] = '0';
    while( l1--    ) stream[i++] = *(bitstrm1++);
    stream[i] = '\0';
    bitstrm1 = stream;
   }
 else if (l2 < l1)
   {
    ldiff = l1 - l2;
    i=0;
    while( ldiff-- ) stream[i++] = '0';
    while( l2--    ) stream[i++] = *(bitstrm2++);
    stream[i] = '\0';
    bitstrm2 = stream;
   }
 while ( (chr1 = *(bitstrm1++)) ) 
    {
       chr2 = *(bitstrm2++);
       if ((chr1 == 'x') || (chr2 == 'x'))
          *result = 'x';
       else if ((chr1 == '1') && (chr2 == '1'))
          *result = '1';
       else
          *result = '0';
       result++;
    }
 free(stream);
 *result = '\0';
}

static void bitor(char *result,char *bitstrm1,char *bitstrm2)
{
 int i, l1, l2, ldiff, largestStream;
 char *stream=0;
 char chr1, chr2;

 l1 = strlen(bitstrm1);
 l2 = strlen(bitstrm2);
 largestStream = (l1 > l2) ? l1 : l2;
 stream = (char *)malloc(sizeof(char)*(largestStream+1));
 if (l1 < l2)
   {
    ldiff = l2 - l1;
    i=0;
    while( ldiff-- ) stream[i++] = '0';
    while( l1--    ) stream[i++] = *(bitstrm1++);
    stream[i] = '\0';
    bitstrm1 = stream;
   }
 else if (l2 < l1)
   {
    ldiff = l1 - l2;
    i=0;
    while( ldiff-- ) stream[i++] = '0';
    while( l2--    ) stream[i++] = *(bitstrm2++);
    stream[i] = '\0';
    bitstrm2 = stream;
   }
 while ( (chr1 = *(bitstrm1++)) ) 
    {
       chr2 = *(bitstrm2++);
       if ((chr1 == '1') || (chr2 == '1'))
          *result = '1';
       else if ((chr1 == '0') || (chr2 == '0'))
          *result = '0';
       else
          *result = 'x';
       result++;
    }
 free(stream);
 *result = '\0';
}

static void bitnot(char *result,char *bits)
{
   int length;
   char chr;

   length = strlen(bits);
   while( length-- ) {
      chr = *(bits++);
      *(result++) = ( chr=='1' ? '0' : ( chr=='0' ? '1' : chr ) );
   }
   *result = '\0';
}

static char bitcmp(char *bitstrm1, char *bitstrm2)
{
 int i, l1, l2, ldiff, largestStream;
 char *stream=0;
 char chr1, chr2;

 l1 = strlen(bitstrm1);
 l2 = strlen(bitstrm2);
 largestStream = (l1 > l2) ? l1 : l2;
 stream = (char *)malloc(sizeof(char)*(largestStream+1));
 if (l1 < l2)
   {
    ldiff = l2 - l1;
    i=0;
    while( ldiff-- ) stream[i++] = '0';
    while( l1--    ) stream[i++] = *(bitstrm1++);
    stream[i] = '\0';
    bitstrm1 = stream;
   }
 else if (l2 < l1)
   {
    ldiff = l1 - l2;
    i=0;
    while( ldiff-- ) stream[i++] = '0';
    while( l2--    ) stream[i++] = *(bitstrm2++);
    stream[i] = '\0';
    bitstrm2 = stream;
   }
 while( (chr1 = *(bitstrm1++)) )
    {
       chr2 = *(bitstrm2++);
       if ( ((chr1 == '0') && (chr2 == '1'))
	    || ((chr1 == '1') && (chr2 == '0')) )
       {
          free(stream);
	  return( 0 );
       }
    }
 free(stream);
 return( 1 );
}

static char bnear(double x, double y, double tolerance)
{
 if (fabs(x - y) < tolerance)
   return ( 1 );
 else
   return ( 0 );
}

static char saobox(double xcen, double ycen, double xwid, double ywid,
		   double rot,  double xcol, double ycol)
{
 double x,y,xprime,yprime,xmin,xmax,ymin,ymax,theta;

 theta = (rot / 180.0) * myPI;
 xprime = xcol - xcen;
 yprime = ycol - ycen;
 x =  xprime * cos(theta) + yprime * sin(theta);
 y = -xprime * sin(theta) + yprime * cos(theta);
 xmin = - 0.5 * xwid; xmax = 0.5 * xwid;
 ymin = - 0.5 * ywid; ymax = 0.5 * ywid;
 if ((x >= xmin) && (x <= xmax) && (y >= ymin) && (y <= ymax))
   return ( 1 );
 else
   return ( 0 );
}

static char circle(double xcen, double ycen, double rad,
		   double xcol, double ycol)
{
 double r2,dx,dy,dlen;

 dx = xcol - xcen;
 dy = ycol - ycen;
 dx *= dx; dy *= dy;
 dlen = dx + dy;
 r2 = rad * rad;
 if (dlen <= r2)
   return ( 1 );
 else
   return ( 0 );
}

static char ellipse(double xcen, double ycen, double xrad, double yrad,
		    double rot, double xcol, double ycol)
{
 double x,y,xprime,yprime,dx,dy,dlen,theta;

 theta = (rot / 180.0) * myPI;
 xprime = xcol - xcen;
 yprime = ycol - ycen;
 x =  xprime * cos(theta) + yprime * sin(theta);
 y = -xprime * sin(theta) + yprime * cos(theta);
 dx = x / xrad; dy = y / yrad;
 dx *= dx; dy *= dy;
 dlen = dx + dy;
 if (dlen <= 1.0)
   return ( 1 );
 else
   return ( 0 );
}

/*
 * Extract substring
 */
 int cstrmid(ParseData *lParse, char *dest_str, int dest_len,
	    char *src_str,  int src_len,
	    int pos)
{
  /* char fill_char = ' '; */
  char fill_char = '\0';
  if (src_len == 0) { src_len = strlen(src_str); } /* .. if constant */

  /* Fill destination with blanks */
  if (pos < 0) { 
    yyerror(0, lParse, "STRMID(S,P,N) P must be 0 or greater");
    return -1;
  }
  if (pos > src_len || pos == 0) {
    /* pos==0: blank string requested */
    memset(dest_str, fill_char, dest_len);
  } else if (pos+dest_len > src_len) {
    /* Copy a subset */
    int nsub = src_len-pos+1;
    int npad = dest_len - nsub;
    memcpy(dest_str, src_str+pos-1, nsub);
    /* Fill remaining string with blanks */
    memset(dest_str+nsub, fill_char, npad);
  } else {
    /* Full string copy */
    memcpy(dest_str, src_str+pos-1, dest_len);
  }
  dest_str[dest_len] = '\0'; /* Null-terminate */

  return 0;
}


static void yyerror(yyscan_t scanner, ParseData *lParse, char *s)
{
    char msg[80];

    if( !lParse->status ) lParse->status = PARSE_SYNTAX_ERR;

    strncpy(msg, s, 80);
    msg[79] = '\0';
    ffpmsg(msg);
}
