
/*  A Bison parser, made from eval.y
 by  GNU Bison version 1.25
  */

#define FFBISON 1  /* Identify Bison output.  */

#define	BOOLEAN	258
#define	LONG	259
#define	DOUBLE	260
#define	STRING	261
#define	BITSTR	262
#define	FUNCTION	263
#define	BFUNCTION	264
#define	IFUNCTION	265
#define	GTIFILTER	266
#define	REGFILTER	267
#define	COLUMN	268
#define	BCOLUMN	269
#define	SCOLUMN	270
#define	BITCOL	271
#define	ROWREF	272
#define	NULLREF	273
#define	SNULLREF	274
#define	OR	275
#define	AND	276
#define	EQ	277
#define	NE	278
#define	GT	279
#define	LT	280
#define	LTE	281
#define	GTE	282
#define	POWER	283
#define	NOT	284
#define	INTCAST	285
#define	FLTCAST	286
#define	UMINUS	287
#define	ACCUM	288
#define	DIFF	289

#line 1 "eval.y"

/************************************************************************/
/*                                                                      */
/*                       CFITSIO Lexical Parser                         */
/*                                                                      */
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

   /*  Shrink the initial stack depth to keep local data <32K (mac limit)  */
   /*  yacc will allocate more space if needed, though.                    */
#define  FFINITDEPTH   100

/***************************************************************/
/*  Replace Bison's BACKUP macro with one that fixes a bug --  */
/*  must update state after popping the stack -- and allows    */
/*  popping multiple terms at one time.                        */
/***************************************************************/

#define FFNEWBACKUP(token, value)                               \
   do								\
     if (ffchar == FFEMPTY )   					\
       { ffchar = (token);                                      \
         memcpy( &fflval, &(value), sizeof(value) );            \
         ffchar1 = FFTRANSLATE (ffchar);			\
         while (fflen--) FFPOPSTACK;				\
         ffstate = *ffssp;					\
         goto ffbackup;						\
       }							\
     else							\
       { fferror ("syntax error: cannot back up"); FFERROR; }	\
   while (0)

/***************************************************************/
/*  Useful macros for accessing/testing Nodes                  */
/***************************************************************/

#define TEST(a)        if( (a)<0 ) FFERROR
#define SIZE(a)        gParse.Nodes[ a ].value.nelem
#define TYPE(a)        gParse.Nodes[ a ].type
#define OPER(a)        gParse.Nodes[ a ].operation
#define PROMOTE(a,b)   if( TYPE(a) > TYPE(b) )                  \
                          b = New_Unary( TYPE(a), 0, b );       \
                       else if( TYPE(a) < TYPE(b) )             \
	                  a = New_Unary( TYPE(b), 0, a );

/*****  Internal functions  *****/

#ifdef __cplusplus
extern "C" {
#endif

static int  Alloc_Node    ( void );
static void Free_Last_Node( void );
static void Evaluate_Node ( int thisNode );

static int  New_Const ( int returnType, void *value, long len );
static int  New_Column( int ColNum );
static int  New_Offset( int ColNum, int offset );
static int  New_Unary ( int returnType, int Op, int Node1 );
static int  New_BinOp ( int returnType, int Node1, int Op, int Node2 );
static int  New_Func  ( int returnType, funcOp Op, int nNodes,
			int Node1, int Node2, int Node3, int Node4, 
			int Node5, int Node6, int Node7 );
static int  New_FuncSize( int returnType, funcOp Op, int nNodes,
			int Node1, int Node2, int Node3, int Node4, 
			  int Node5, int Node6, int Node7, int Size);
static int  New_Deref ( int Var,  int nDim,
			int Dim1, int Dim2, int Dim3, int Dim4, int Dim5 );
static int  New_GTI   ( char *fname, int Node1, char *start, char *stop );
static int  New_REG   ( char *fname, int NodeX, int NodeY, char *colNames );
static int  New_Vector( int subNode );
static int  Close_Vec ( int vecNode );
static int  Locate_Col( Node *this );
static int  Test_Dims ( int Node1, int Node2 );
static void Copy_Dims ( int Node1, int Node2 );

static void Allocate_Ptrs( Node *this );
static void Do_Unary     ( Node *this );
static void Do_Offset    ( Node *this );
static void Do_BinOp_bit ( Node *this );
static void Do_BinOp_str ( Node *this );
static void Do_BinOp_log ( Node *this );
static void Do_BinOp_lng ( Node *this );
static void Do_BinOp_dbl ( Node *this );
static void Do_Func      ( Node *this );
static void Do_Deref     ( Node *this );
static void Do_GTI       ( Node *this );
static void Do_REG       ( Node *this );
static void Do_Vector    ( Node *this );

static long Search_GTI   ( double evtTime, long nGTI, double *start,
			   double *stop, int ordered );

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
static int cstrmid(char *dest_str, int dest_len,
		   char *src_str,  int src_len, int pos);

static void  fferror(char *msg);

#ifdef __cplusplus
    }
#endif


#line 189 "eval.y"
typedef union {
    int    Node;        /* Index of Node */
    double dbl;         /* real value    */
    long   lng;         /* integer value */
    char   log;         /* logical value */
    char   str[MAX_STRLEN];    /* string value  */
} FFSTYPE;
#include <stdio.h>

#ifndef __cplusplus
#ifndef __STDC__
#define const
#endif
#endif



#define	FFFINAL		290
#define	FFFLAG		-32768
#define	FFNTBASE	54

#define FFTRANSLATE(x) ((unsigned)(x) <= 289 ? fftranslate[x] : 62)

static const char fftranslate[] = {     0,
     2,     2,     2,     2,     2,     2,     2,     2,     2,    50,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,    37,    41,     2,    52,
    53,    38,    35,    20,    36,     2,    39,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,    22,     2,     2,
    21,     2,    25,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
    47,     2,    51,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,    23,    40,    24,    30,     2,     2,     2,     2,
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
     2,     2,     2,     2,     2,     1,     2,     3,     4,     5,
     6,     7,     8,     9,    10,    11,    12,    13,    14,    15,
    16,    17,    18,    19,    26,    27,    28,    29,    31,    32,
    33,    34,    42,    43,    44,    45,    46,    48,    49
};

#if FFDEBUG != 0
static const short ffprhs[] = {     0,
     0,     1,     4,     6,     9,    12,    15,    18,    21,    24,
    28,    31,    35,    39,    43,    46,    49,    51,    53,    58,
    62,    66,    70,    75,    82,    91,   102,   115,   118,   122,
   124,   126,   128,   133,   135,   137,   141,   145,   149,   153,
   157,   161,   164,   167,   171,   175,   179,   185,   191,   197,
   200,   204,   208,   212,   216,   222,   228,   238,   243,   250,
   259,   270,   283,   286,   289,   292,   295,   297,   299,   304,
   308,   312,   316,   320,   324,   328,   332,   336,   340,   344,
   348,   352,   356,   360,   364,   368,   372,   376,   380,   384,
   388,   392,   396,   402,   408,   412,   416,   420,   426,   434,
   446,   462,   465,   469,   475,   485,   489,   497,   507,   512,
   519,   528,   539,   552,   555,   559,   561,   563,   568,   570,
   574,   578,   584,   590
};

static const short ffrhs[] = {    -1,
    54,    55,     0,    50,     0,    58,    50,     0,    59,    50,
     0,    61,    50,     0,    60,    50,     0,     1,    50,     0,
    23,    59,     0,    56,    20,    59,     0,    23,    58,     0,
    57,    20,    58,     0,    57,    20,    59,     0,    56,    20,
    58,     0,    57,    24,     0,    56,    24,     0,     7,     0,
    16,     0,    16,    23,    58,    24,     0,    60,    41,    60,
     0,    60,    40,    60,     0,    60,    35,    60,     0,    60,
    47,    58,    51,     0,    60,    47,    58,    20,    58,    51,
     0,    60,    47,    58,    20,    58,    20,    58,    51,     0,
    60,    47,    58,    20,    58,    20,    58,    20,    58,    51,
     0,    60,    47,    58,    20,    58,    20,    58,    20,    58,
    20,    58,    51,     0,    43,    60,     0,    52,    60,    53,
     0,     4,     0,     5,     0,    13,     0,    13,    23,    58,
    24,     0,    17,     0,    18,     0,    58,    37,    58,     0,
    58,    35,    58,     0,    58,    36,    58,     0,    58,    38,
    58,     0,    58,    39,    58,     0,    58,    42,    58,     0,
    35,    58,     0,    36,    58,     0,    52,    58,    53,     0,
    58,    38,    59,     0,    59,    38,    58,     0,    59,    25,
    58,    22,    58,     0,    59,    25,    59,    22,    58,     0,
    59,    25,    58,    22,    59,     0,     8,    53,     0,     8,
    59,    53,     0,     8,    61,    53,     0,     8,    60,    53,
     0,     8,    58,    53,     0,    10,    61,    20,    61,    53,
     0,     8,    58,    20,    58,    53,     0,     8,    58,    20,
    58,    20,    58,    20,    58,    53,     0,    58,    47,    58,
    51,     0,    58,    47,    58,    20,    58,    51,     0,    58,
    47,    58,    20,    58,    20,    58,    51,     0,    58,    47,
    58,    20,    58,    20,    58,    20,    58,    51,     0,    58,
    47,    58,    20,    58,    20,    58,    20,    58,    20,    58,
    51,     0,    44,    58,     0,    44,    59,     0,    45,    58,
     0,    45,    59,     0,     3,     0,    14,     0,    14,    23,
    58,    24,     0,    60,    28,    60,     0,    60,    29,    60,
     0,    60,    32,    60,     0,    60,    33,    60,     0,    60,
    31,    60,     0,    60,    34,    60,     0,    58,    31,    58,
     0,    58,    32,    58,     0,    58,    34,    58,     0,    58,
    33,    58,     0,    58,    30,    58,     0,    58,    28,    58,
     0,    58,    29,    58,     0,    61,    28,    61,     0,    61,
    29,    61,     0,    61,    31,    61,     0,    61,    34,    61,
     0,    61,    32,    61,     0,    61,    33,    61,     0,    59,
    27,    59,     0,    59,    26,    59,     0,    59,    28,    59,
     0,    59,    29,    59,     0,    58,    21,    58,    22,    58,
     0,    59,    25,    59,    22,    59,     0,     9,    58,    53,
     0,     9,    59,    53,     0,     9,    61,    53,     0,     8,
    59,    20,    59,    53,     0,     9,    58,    20,    58,    20,
    58,    53,     0,     9,    58,    20,    58,    20,    58,    20,
    58,    20,    58,    53,     0,     9,    58,    20,    58,    20,
    58,    20,    58,    20,    58,    20,    58,    20,    58,    53,
     0,    11,    53,     0,    11,     6,    53,     0,    11,     6,
    20,    58,    53,     0,    11,     6,    20,    58,    20,     6,
    20,     6,    53,     0,    12,     6,    53,     0,    12,     6,
    20,    58,    20,    58,    53,     0,    12,     6,    20,    58,
    20,    58,    20,     6,    53,     0,    59,    47,    58,    51,
     0,    59,    47,    58,    20,    58,    51,     0,    59,    47,
    58,    20,    58,    20,    58,    51,     0,    59,    47,    58,
    20,    58,    20,    58,    20,    58,    51,     0,    59,    47,
    58,    20,    58,    20,    58,    20,    58,    20,    58,    51,
     0,    43,    59,     0,    52,    59,    53,     0,     6,     0,
    15,     0,    15,    23,    58,    24,     0,    19,     0,    52,
    61,    53,     0,    61,    35,    61,     0,    59,    25,    61,
    22,    61,     0,     8,    61,    20,    61,    53,     0,     8,
    61,    20,    58,    20,    58,    53,     0
};

#endif

#if FFDEBUG != 0
static const short ffrline[] = { 0,
   241,   242,   245,   246,   252,   258,   264,   270,   273,   275,
   288,   290,   303,   314,   328,   332,   336,   340,   342,   351,
   354,   357,   366,   368,   370,   372,   374,   376,   379,   383,
   385,   387,   389,   398,   400,   402,   405,   408,   411,   414,
   417,   420,   422,   424,   426,   430,   434,   453,   472,   491,
   504,   518,   530,   561,   659,   667,   729,   753,   755,   757,
   759,   761,   763,   765,   767,   769,   773,   775,   777,   786,
   789,   792,   795,   798,   801,   804,   807,   810,   813,   816,
   819,   822,   825,   828,   831,   834,   837,   840,   843,   845,
   847,   849,   852,   859,   876,   889,   902,   913,   929,   953,
   981,  1018,  1022,  1026,  1029,  1033,  1037,  1040,  1044,  1046,
  1048,  1050,  1052,  1054,  1056,  1060,  1063,  1065,  1074,  1076,
  1078,  1087,  1106,  1125
};
#endif


#if FFDEBUG != 0 || defined (FFERROR_VERBOSE)

static const char * const fftname[] = {   "$","error","$undefined.","BOOLEAN",
"LONG","DOUBLE","STRING","BITSTR","FUNCTION","BFUNCTION","IFUNCTION","GTIFILTER",
"REGFILTER","COLUMN","BCOLUMN","SCOLUMN","BITCOL","ROWREF","NULLREF","SNULLREF",
"','","'='","':'","'{'","'}'","'?'","OR","AND","EQ","NE","'~'","GT","LT","LTE",
"GTE","'+'","'-'","'%'","'*'","'/'","'|'","'&'","POWER","NOT","INTCAST","FLTCAST",
"UMINUS","'['","ACCUM","DIFF","'\\n'","']'","'('","')'","lines","line","bvector",
"vector","expr","bexpr","bits","sexpr", NULL
};
#endif

static const short ffr1[] = {     0,
    54,    54,    55,    55,    55,    55,    55,    55,    56,    56,
    57,    57,    57,    57,    58,    59,    60,    60,    60,    60,
    60,    60,    60,    60,    60,    60,    60,    60,    60,    58,
    58,    58,    58,    58,    58,    58,    58,    58,    58,    58,
    58,    58,    58,    58,    58,    58,    58,    58,    58,    58,
    58,    58,    58,    58,    58,    58,    58,    58,    58,    58,
    58,    58,    58,    58,    58,    58,    59,    59,    59,    59,
    59,    59,    59,    59,    59,    59,    59,    59,    59,    59,
    59,    59,    59,    59,    59,    59,    59,    59,    59,    59,
    59,    59,    59,    59,    59,    59,    59,    59,    59,    59,
    59,    59,    59,    59,    59,    59,    59,    59,    59,    59,
    59,    59,    59,    59,    59,    61,    61,    61,    61,    61,
    61,    61,    61,    61
};

static const short ffr2[] = {     0,
     0,     2,     1,     2,     2,     2,     2,     2,     2,     3,
     2,     3,     3,     3,     2,     2,     1,     1,     4,     3,
     3,     3,     4,     6,     8,    10,    12,     2,     3,     1,
     1,     1,     4,     1,     1,     3,     3,     3,     3,     3,
     3,     2,     2,     3,     3,     3,     5,     5,     5,     2,
     3,     3,     3,     3,     5,     5,     9,     4,     6,     8,
    10,    12,     2,     2,     2,     2,     1,     1,     4,     3,
     3,     3,     3,     3,     3,     3,     3,     3,     3,     3,
     3,     3,     3,     3,     3,     3,     3,     3,     3,     3,
     3,     3,     5,     5,     3,     3,     3,     5,     7,    11,
    15,     2,     3,     5,     9,     3,     7,     9,     4,     6,
     8,    10,    12,     2,     3,     1,     1,     4,     1,     3,
     3,     5,     5,     7
};

static const short ffdefact[] = {     1,
     0,     0,    67,    30,    31,   116,    17,     0,     0,     0,
     0,     0,    32,    68,   117,    18,    34,    35,   119,     0,
     0,     0,     0,     0,     0,     3,     0,     2,     0,     0,
     0,     0,     0,     0,     8,    50,     0,     0,     0,     0,
     0,     0,     0,     0,     0,     0,     0,     0,   102,     0,
     0,     0,     0,     0,    11,     9,     0,    42,    43,   114,
    28,    63,    64,    65,    66,     0,     0,     0,     0,     0,
    16,     0,    15,     0,     0,     0,     0,     0,     0,     0,
     0,     0,     0,     0,     0,     0,     0,     0,     4,     0,
     0,     0,     0,     0,     0,     0,     5,     0,     0,     0,
     0,     0,     0,     0,     0,     0,     0,     7,     0,     0,
     0,     0,     0,     0,     0,     6,     0,    54,     0,    51,
    53,     0,    52,     0,    95,    96,    97,     0,     0,   103,
     0,   106,     0,     0,     0,     0,    44,   115,    29,   120,
    14,    10,    12,    13,     0,    81,    82,    80,    76,    77,
    79,    78,    37,    38,    36,    39,    45,    40,    41,     0,
     0,     0,     0,    90,    89,    91,    92,    46,     0,     0,
     0,    70,    71,    74,    72,    73,    75,    22,    21,    20,
     0,    83,    84,    85,    87,    88,    86,   121,     0,     0,
     0,     0,     0,     0,     0,     0,    33,    69,   118,    19,
     0,     0,    58,     0,     0,     0,     0,   109,    28,     0,
     0,    23,     0,    56,    98,     0,   123,     0,    55,     0,
   104,     0,    93,     0,    47,    49,    48,    94,   122,     0,
     0,     0,     0,     0,     0,     0,     0,    59,     0,   110,
     0,    24,     0,   124,     0,    99,     0,     0,   107,     0,
     0,     0,     0,     0,     0,     0,     0,    60,     0,   111,
     0,    25,    57,     0,   105,   108,     0,     0,     0,     0,
     0,    61,     0,   112,     0,    26,     0,   100,     0,     0,
     0,     0,    62,   113,    27,     0,     0,   101,     0,     0
};

static const short ffdefgoto[] = {     1,
    28,    29,    30,    45,    46,    43,    57
};

static const short ffpact[] = {-32768,
   301,   -41,-32768,-32768,-32768,-32768,-32768,   351,   402,   402,
    -5,    12,     8,    33,    34,    41,-32768,-32768,-32768,   402,
   402,   402,   402,   402,   402,-32768,   402,-32768,   -18,     9,
  1092,   403,  1438,    79,-32768,-32768,   428,   143,   294,    10,
   456,   224,  1478,   125,  1390,  1436,  1523,    -6,-32768,     2,
   402,   402,   402,   402,  1390,  1436,  1129,    19,    19,    20,
    21,    19,    20,    19,    20,   623,   240,   344,  1120,   402,
-32768,   402,-32768,   402,   402,   402,   402,   402,   402,   402,
   402,   402,   402,   402,   402,   402,   402,   402,-32768,   402,
   402,   402,   402,   402,   402,   402,-32768,    -3,    -3,    -3,
    -3,    -3,    -3,    -3,    -3,    -3,   402,-32768,   402,   402,
   402,   402,   402,   402,   402,-32768,   402,-32768,   402,-32768,
-32768,   402,-32768,   402,-32768,-32768,-32768,   402,   402,-32768,
   402,-32768,  1266,  1286,  1306,  1326,-32768,-32768,-32768,-32768,
  1390,  1436,  1390,  1436,  1348,  1503,  1503,  1503,    23,    23,
    23,    23,   160,   160,   160,   -15,    20,   -15,   -15,   732,
  1370,  1413,  1531,   146,   -13,   -35,   -35,   -15,   756,    -3,
    -3,   -30,   -30,   -30,   -30,   -30,   -30,    50,    21,    21,
   780,    67,    67,    11,    11,    11,    11,-32768,   484,  1118,
  1146,  1415,  1166,  1424,   512,  1186,-32768,-32768,-32768,-32768,
   402,   402,-32768,   402,   402,   402,   402,-32768,    21,  1480,
   402,-32768,   402,-32768,-32768,   402,-32768,   402,-32768,    66,
-32768,   402,  1461,   804,  1461,  1436,  1461,  1436,  1129,   828,
   852,  1206,   650,   540,    68,   568,   402,-32768,   402,-32768,
   402,-32768,   402,-32768,   402,-32768,    86,    87,-32768,   876,
   900,   924,   677,  1226,    52,    56,   402,-32768,   402,-32768,
   402,-32768,-32768,   402,-32768,-32768,   948,   972,   996,   596,
   402,-32768,   402,-32768,   402,-32768,   402,-32768,  1020,  1044,
  1068,  1246,-32768,-32768,-32768,   402,   704,-32768,   126,-32768
};

static const short ffpgoto[] = {-32768,
-32768,-32768,-32768,    -1,    95,   124,    27
};


#define	FFLAST		1566


static const short fftable[] = {    31,
    48,    70,    95,     7,   104,    71,    37,    41,    35,   105,
   106,    96,    16,   129,    93,    94,   107,    50,    55,    58,
    59,   131,    62,    64,    95,    66,    87,    34,    72,   122,
    51,    88,    73,    96,    40,    44,    47,   109,   110,   170,
   111,   112,   113,   114,   115,   115,   130,    49,   171,   133,
   134,   135,   136,    69,   132,    52,    53,    82,    83,    84,
    85,    86,   123,    54,    87,    88,    96,   107,   141,    88,
   143,   235,   145,   146,   147,   148,   149,   150,   151,   152,
   153,   154,   155,   156,   158,   159,   160,   247,   161,   105,
   106,   255,   256,   168,   169,    32,   107,   111,   112,   113,
   114,   115,    38,    42,   265,   181,   109,   110,   266,   111,
   112,   113,   114,   115,    56,   189,   163,    60,    63,    65,
   191,    67,   193,     0,    33,   290,     0,   195,   116,   196,
     0,    39,     0,     0,     0,   182,   183,   184,   185,   186,
   187,   188,     0,     0,     0,     0,    61,     0,   192,     0,
    68,     0,   109,   110,   194,   111,   112,   113,   114,   115,
     0,     0,   119,     0,   142,     0,   144,    90,    91,    92,
    93,    94,    92,    93,    94,     0,     0,   127,     0,   157,
    95,     0,     0,    95,   162,   164,   165,   166,   167,    96,
     0,     0,    96,     0,     0,   120,     0,    85,    86,   223,
   224,    87,   225,   227,     0,   230,    88,     0,     0,   231,
     0,   232,     0,   190,   233,     0,   234,     0,     0,     0,
   236,   172,   173,   174,   175,   176,   177,   178,   179,   180,
     0,     0,   229,     0,     0,   250,     0,   251,     0,   252,
     0,   253,     0,   254,     0,     0,     0,     0,    90,    91,
    92,    93,    94,     0,     0,   267,     0,   268,     0,   269,
     0,    95,   270,     0,    90,    91,    92,    93,    94,   279,
    96,   280,     0,   281,     0,   282,   126,    95,     0,     0,
     0,     0,     0,     0,   287,     0,    96,     0,     0,     0,
     0,     0,   138,   209,   210,     0,     0,     0,   226,   228,
   289,     2,     0,     3,     4,     5,     6,     7,     8,     9,
    10,    11,    12,    13,    14,    15,    16,    17,    18,    19,
     0,    98,    99,    20,   100,   101,   102,   103,   104,     0,
     0,     0,     0,   105,   106,    21,    22,     0,     0,     0,
   107,     0,     0,    23,    24,    25,   121,     0,     0,     0,
    26,     0,    27,     3,     4,     5,     6,     7,     8,     9,
    10,    11,    12,    13,    14,    15,    16,    17,    18,    19,
     0,    98,    99,    20,   100,   101,   102,   103,   104,     0,
     0,     0,     0,   105,   106,    21,    22,     0,     0,     0,
   107,     0,     0,    23,    24,    25,   139,     0,     0,     0,
     0,     0,    27,    36,     3,     4,     5,     6,     7,     8,
     9,    10,    11,    12,    13,    14,    15,    16,    17,    18,
    19,     0,     0,     0,    20,     0,     0,    90,    91,    92,
    93,    94,     0,     0,     0,     0,    21,    22,     0,     0,
    95,     0,     0,     0,    23,    24,    25,   117,    74,    96,
     0,     0,    97,    27,     0,    75,    76,    77,    78,    79,
    80,    81,    82,    83,    84,    85,    86,     0,     0,    87,
     0,     0,     0,     0,    88,   124,    74,     0,     0,     0,
   118,     0,     0,    75,    76,    77,    78,    79,    80,    81,
    82,    83,    84,    85,    86,     0,     0,    87,     0,     0,
     0,     0,    88,   213,    74,     0,     0,     0,   125,     0,
     0,    75,    76,    77,    78,    79,    80,    81,    82,    83,
    84,    85,    86,     0,     0,    87,     0,     0,     0,     0,
    88,   220,    74,     0,     0,     0,   214,     0,     0,    75,
    76,    77,    78,    79,    80,    81,    82,    83,    84,    85,
    86,     0,     0,    87,     0,     0,     0,     0,    88,   245,
    74,     0,     0,     0,   221,     0,     0,    75,    76,    77,
    78,    79,    80,    81,    82,    83,    84,    85,    86,     0,
     0,    87,     0,     0,     0,     0,    88,   248,    74,     0,
     0,     0,   246,     0,     0,    75,    76,    77,    78,    79,
    80,    81,    82,    83,    84,    85,    86,     0,     0,    87,
     0,     0,     0,     0,    88,   277,    74,     0,     0,     0,
   249,     0,     0,    75,    76,    77,    78,    79,    80,    81,
    82,    83,    84,    85,    86,     0,     0,    87,     0,     0,
     0,     0,    88,    74,     0,     0,     0,     0,   278,     0,
    75,    76,    77,    78,    79,    80,    81,    82,    83,    84,
    85,    86,     0,     0,    87,     0,     0,     0,     0,    88,
    74,     0,     0,     0,     0,   137,     0,    75,    76,    77,
    78,    79,    80,    81,    82,    83,    84,    85,    86,     0,
     0,    87,     0,     0,     0,     0,    88,    74,     0,     0,
     0,     0,   244,     0,    75,    76,    77,    78,    79,    80,
    81,    82,    83,    84,    85,    86,     0,     0,    87,     0,
     0,     0,     0,    88,    74,     0,     0,     0,     0,   263,
     0,    75,    76,    77,    78,    79,    80,    81,    82,    83,
    84,    85,    86,     0,     0,    87,     0,     0,     0,     0,
    88,   202,    74,     0,     0,     0,   288,     0,     0,    75,
    76,    77,    78,    79,    80,    81,    82,    83,    84,    85,
    86,     0,     0,    87,     0,   207,    74,     0,    88,     0,
     0,     0,   203,    75,    76,    77,    78,    79,    80,    81,
    82,    83,    84,    85,    86,     0,     0,    87,     0,   211,
    74,     0,    88,     0,     0,     0,   208,    75,    76,    77,
    78,    79,    80,    81,    82,    83,    84,    85,    86,     0,
     0,    87,     0,   237,    74,     0,    88,     0,     0,     0,
   212,    75,    76,    77,    78,    79,    80,    81,    82,    83,
    84,    85,    86,     0,     0,    87,     0,   239,    74,     0,
    88,     0,     0,     0,   238,    75,    76,    77,    78,    79,
    80,    81,    82,    83,    84,    85,    86,     0,     0,    87,
     0,   241,    74,     0,    88,     0,     0,     0,   240,    75,
    76,    77,    78,    79,    80,    81,    82,    83,    84,    85,
    86,     0,     0,    87,     0,   257,    74,     0,    88,     0,
     0,     0,   242,    75,    76,    77,    78,    79,    80,    81,
    82,    83,    84,    85,    86,     0,     0,    87,     0,   259,
    74,     0,    88,     0,     0,     0,   258,    75,    76,    77,
    78,    79,    80,    81,    82,    83,    84,    85,    86,     0,
     0,    87,     0,   261,    74,     0,    88,     0,     0,     0,
   260,    75,    76,    77,    78,    79,    80,    81,    82,    83,
    84,    85,    86,     0,     0,    87,     0,   271,    74,     0,
    88,     0,     0,     0,   262,    75,    76,    77,    78,    79,
    80,    81,    82,    83,    84,    85,    86,     0,     0,    87,
     0,   273,    74,     0,    88,     0,     0,     0,   272,    75,
    76,    77,    78,    79,    80,    81,    82,    83,    84,    85,
    86,     0,     0,    87,     0,   275,    74,     0,    88,     0,
     0,     0,   274,    75,    76,    77,    78,    79,    80,    81,
    82,    83,    84,    85,    86,     0,     0,    87,     0,     0,
    74,     0,    88,     0,     0,     0,   276,    75,    76,    77,
    78,    79,    80,    81,    82,    83,    84,    85,    86,     0,
     0,    87,     0,     0,    74,     0,    88,     0,     0,     0,
   283,    75,    76,    77,    78,    79,    80,    81,    82,    83,
    84,    85,    86,     0,     0,    87,     0,     0,    74,     0,
    88,     0,     0,     0,   284,    75,    76,    77,    78,    79,
    80,    81,    82,    83,    84,    85,    86,     0,     0,    87,
     0,     0,    74,     0,    88,     0,     0,     0,   285,    75,
    76,    77,    78,    79,    80,    81,    82,    83,    84,    85,
    86,     0,     0,    87,     0,     0,     0,     0,    88,     0,
     0,    89,    90,    91,    92,    93,    94,   109,   110,     0,
   111,   112,   113,   114,   115,    95,   109,   110,     0,   111,
   112,   113,   114,   115,    96,   216,    74,     0,     0,     0,
   215,     0,   140,    75,    76,    77,    78,    79,    80,    81,
    82,    83,    84,    85,    86,   218,    74,    87,     0,     0,
     0,     0,    88,    75,    76,    77,    78,    79,    80,    81,
    82,    83,    84,    85,    86,   222,    74,    87,     0,     0,
     0,     0,    88,    75,    76,    77,    78,    79,    80,    81,
    82,    83,    84,    85,    86,   243,    74,    87,     0,     0,
     0,     0,    88,    75,    76,    77,    78,    79,    80,    81,
    82,    83,    84,    85,    86,   264,    74,    87,     0,     0,
     0,     0,    88,    75,    76,    77,    78,    79,    80,    81,
    82,    83,    84,    85,    86,   286,    74,    87,     0,     0,
     0,     0,    88,    75,    76,    77,    78,    79,    80,    81,
    82,    83,    84,    85,    86,     0,    74,    87,     0,   197,
     0,     0,    88,    75,    76,    77,    78,    79,    80,    81,
    82,    83,    84,    85,    86,     0,    74,    87,     0,   198,
     0,     0,    88,    75,    76,    77,    78,    79,    80,    81,
    82,    83,    84,    85,    86,     0,    74,    87,     0,   199,
     0,     0,    88,    75,    76,    77,    78,    79,    80,    81,
    82,    83,    84,    85,    86,     0,    74,    87,     0,   200,
     0,     0,    88,    75,    76,    77,    78,    79,    80,    81,
    82,    83,    84,    85,    86,     0,     0,    87,    74,   201,
     0,     0,    88,     0,     0,    75,    76,    77,    78,    79,
    80,    81,    82,    83,    84,    85,    86,     0,     0,    87,
    74,   204,     0,     0,    88,     0,     0,    75,    76,    77,
    78,    79,    80,    81,    82,    83,    84,    85,    86,     0,
    74,    87,     0,     0,     0,     0,    88,    75,    76,    77,
    78,    79,    80,    81,    82,    83,    84,    85,    86,     0,
     0,    87,     0,     0,   205,     0,    88,    90,    91,    92,
    93,    94,   109,   110,     0,   111,   112,   113,   114,   115,
    95,   109,   110,     0,   111,   112,   113,   114,   115,    96,
    90,    91,    92,    93,    94,    98,    99,   217,   100,   101,
   102,   103,   104,    95,     0,     0,   219,   105,   106,     0,
     0,     0,    96,     0,   107,     0,     0,   108,    75,    76,
    77,    78,    79,    80,    81,    82,    83,    84,    85,    86,
     0,     0,    87,     0,     0,    98,    99,    88,   100,   101,
   102,   103,   104,     0,   104,     0,     0,   105,   106,   105,
   106,     0,     0,     0,   107,     0,   107,     0,     0,     0,
     0,     0,   139,    78,    79,    80,    81,    82,    83,    84,
    85,    86,   128,     0,    87,     0,     0,     0,     0,    88,
   109,   110,   206,   111,   112,   113,   114,   115,   109,   110,
     0,   111,   112,   113,   114,   115
};

static const short ffcheck[] = {     1,
     6,    20,    38,     7,    35,    24,     8,     9,    50,    40,
    41,    47,    16,    20,    28,    29,    47,     6,    20,    21,
    22,    20,    24,    25,    38,    27,    42,     1,    20,    20,
    23,    47,    24,    47,     8,     9,    10,    28,    29,    43,
    31,    32,    33,    34,    35,    35,    53,    53,    52,    51,
    52,    53,    54,    27,    53,    23,    23,    35,    36,    37,
    38,    39,    53,    23,    42,    47,    47,    47,    70,    47,
    72,     6,    74,    75,    76,    77,    78,    79,    80,    81,
    82,    83,    84,    85,    86,    87,    88,    20,    90,    40,
    41,     6,     6,    95,    96,     1,    47,    31,    32,    33,
    34,    35,     8,     9,    53,   107,    28,    29,    53,    31,
    32,    33,    34,    35,    20,   117,    90,    23,    24,    25,
   122,    27,   124,    -1,     1,     0,    -1,   129,    50,   131,
    -1,     8,    -1,    -1,    -1,   109,   110,   111,   112,   113,
   114,   115,    -1,    -1,    -1,    -1,    23,    -1,   122,    -1,
    27,    -1,    28,    29,   128,    31,    32,    33,    34,    35,
    -1,    -1,    20,    -1,    70,    -1,    72,    25,    26,    27,
    28,    29,    27,    28,    29,    -1,    -1,    53,    -1,    85,
    38,    -1,    -1,    38,    90,    91,    92,    93,    94,    47,
    -1,    -1,    47,    -1,    -1,    53,    -1,    38,    39,   201,
   202,    42,   204,   205,    -1,   207,    47,    -1,    -1,   211,
    -1,   213,    -1,   119,   216,    -1,   218,    -1,    -1,    -1,
   222,    98,    99,   100,   101,   102,   103,   104,   105,   106,
    -1,    -1,   206,    -1,    -1,   237,    -1,   239,    -1,   241,
    -1,   243,    -1,   245,    -1,    -1,    -1,    -1,    25,    26,
    27,    28,    29,    -1,    -1,   257,    -1,   259,    -1,   261,
    -1,    38,   264,    -1,    25,    26,    27,    28,    29,   271,
    47,   273,    -1,   275,    -1,   277,    53,    38,    -1,    -1,
    -1,    -1,    -1,    -1,   286,    -1,    47,    -1,    -1,    -1,
    -1,    -1,    53,   170,   171,    -1,    -1,    -1,   204,   205,
     0,     1,    -1,     3,     4,     5,     6,     7,     8,     9,
    10,    11,    12,    13,    14,    15,    16,    17,    18,    19,
    -1,    28,    29,    23,    31,    32,    33,    34,    35,    -1,
    -1,    -1,    -1,    40,    41,    35,    36,    -1,    -1,    -1,
    47,    -1,    -1,    43,    44,    45,    53,    -1,    -1,    -1,
    50,    -1,    52,     3,     4,     5,     6,     7,     8,     9,
    10,    11,    12,    13,    14,    15,    16,    17,    18,    19,
    -1,    28,    29,    23,    31,    32,    33,    34,    35,    -1,
    -1,    -1,    -1,    40,    41,    35,    36,    -1,    -1,    -1,
    47,    -1,    -1,    43,    44,    45,    53,    -1,    -1,    -1,
    -1,    -1,    52,    53,     3,     4,     5,     6,     7,     8,
     9,    10,    11,    12,    13,    14,    15,    16,    17,    18,
    19,    -1,    -1,    -1,    23,    -1,    -1,    25,    26,    27,
    28,    29,    -1,    -1,    -1,    -1,    35,    36,    -1,    -1,
    38,    -1,    -1,    -1,    43,    44,    45,    20,    21,    47,
    -1,    -1,    50,    52,    -1,    28,    29,    30,    31,    32,
    33,    34,    35,    36,    37,    38,    39,    -1,    -1,    42,
    -1,    -1,    -1,    -1,    47,    20,    21,    -1,    -1,    -1,
    53,    -1,    -1,    28,    29,    30,    31,    32,    33,    34,
    35,    36,    37,    38,    39,    -1,    -1,    42,    -1,    -1,
    -1,    -1,    47,    20,    21,    -1,    -1,    -1,    53,    -1,
    -1,    28,    29,    30,    31,    32,    33,    34,    35,    36,
    37,    38,    39,    -1,    -1,    42,    -1,    -1,    -1,    -1,
    47,    20,    21,    -1,    -1,    -1,    53,    -1,    -1,    28,
    29,    30,    31,    32,    33,    34,    35,    36,    37,    38,
    39,    -1,    -1,    42,    -1,    -1,    -1,    -1,    47,    20,
    21,    -1,    -1,    -1,    53,    -1,    -1,    28,    29,    30,
    31,    32,    33,    34,    35,    36,    37,    38,    39,    -1,
    -1,    42,    -1,    -1,    -1,    -1,    47,    20,    21,    -1,
    -1,    -1,    53,    -1,    -1,    28,    29,    30,    31,    32,
    33,    34,    35,    36,    37,    38,    39,    -1,    -1,    42,
    -1,    -1,    -1,    -1,    47,    20,    21,    -1,    -1,    -1,
    53,    -1,    -1,    28,    29,    30,    31,    32,    33,    34,
    35,    36,    37,    38,    39,    -1,    -1,    42,    -1,    -1,
    -1,    -1,    47,    21,    -1,    -1,    -1,    -1,    53,    -1,
    28,    29,    30,    31,    32,    33,    34,    35,    36,    37,
    38,    39,    -1,    -1,    42,    -1,    -1,    -1,    -1,    47,
    21,    -1,    -1,    -1,    -1,    53,    -1,    28,    29,    30,
    31,    32,    33,    34,    35,    36,    37,    38,    39,    -1,
    -1,    42,    -1,    -1,    -1,    -1,    47,    21,    -1,    -1,
    -1,    -1,    53,    -1,    28,    29,    30,    31,    32,    33,
    34,    35,    36,    37,    38,    39,    -1,    -1,    42,    -1,
    -1,    -1,    -1,    47,    21,    -1,    -1,    -1,    -1,    53,
    -1,    28,    29,    30,    31,    32,    33,    34,    35,    36,
    37,    38,    39,    -1,    -1,    42,    -1,    -1,    -1,    -1,
    47,    20,    21,    -1,    -1,    -1,    53,    -1,    -1,    28,
    29,    30,    31,    32,    33,    34,    35,    36,    37,    38,
    39,    -1,    -1,    42,    -1,    20,    21,    -1,    47,    -1,
    -1,    -1,    51,    28,    29,    30,    31,    32,    33,    34,
    35,    36,    37,    38,    39,    -1,    -1,    42,    -1,    20,
    21,    -1,    47,    -1,    -1,    -1,    51,    28,    29,    30,
    31,    32,    33,    34,    35,    36,    37,    38,    39,    -1,
    -1,    42,    -1,    20,    21,    -1,    47,    -1,    -1,    -1,
    51,    28,    29,    30,    31,    32,    33,    34,    35,    36,
    37,    38,    39,    -1,    -1,    42,    -1,    20,    21,    -1,
    47,    -1,    -1,    -1,    51,    28,    29,    30,    31,    32,
    33,    34,    35,    36,    37,    38,    39,    -1,    -1,    42,
    -1,    20,    21,    -1,    47,    -1,    -1,    -1,    51,    28,
    29,    30,    31,    32,    33,    34,    35,    36,    37,    38,
    39,    -1,    -1,    42,    -1,    20,    21,    -1,    47,    -1,
    -1,    -1,    51,    28,    29,    30,    31,    32,    33,    34,
    35,    36,    37,    38,    39,    -1,    -1,    42,    -1,    20,
    21,    -1,    47,    -1,    -1,    -1,    51,    28,    29,    30,
    31,    32,    33,    34,    35,    36,    37,    38,    39,    -1,
    -1,    42,    -1,    20,    21,    -1,    47,    -1,    -1,    -1,
    51,    28,    29,    30,    31,    32,    33,    34,    35,    36,
    37,    38,    39,    -1,    -1,    42,    -1,    20,    21,    -1,
    47,    -1,    -1,    -1,    51,    28,    29,    30,    31,    32,
    33,    34,    35,    36,    37,    38,    39,    -1,    -1,    42,
    -1,    20,    21,    -1,    47,    -1,    -1,    -1,    51,    28,
    29,    30,    31,    32,    33,    34,    35,    36,    37,    38,
    39,    -1,    -1,    42,    -1,    20,    21,    -1,    47,    -1,
    -1,    -1,    51,    28,    29,    30,    31,    32,    33,    34,
    35,    36,    37,    38,    39,    -1,    -1,    42,    -1,    -1,
    21,    -1,    47,    -1,    -1,    -1,    51,    28,    29,    30,
    31,    32,    33,    34,    35,    36,    37,    38,    39,    -1,
    -1,    42,    -1,    -1,    21,    -1,    47,    -1,    -1,    -1,
    51,    28,    29,    30,    31,    32,    33,    34,    35,    36,
    37,    38,    39,    -1,    -1,    42,    -1,    -1,    21,    -1,
    47,    -1,    -1,    -1,    51,    28,    29,    30,    31,    32,
    33,    34,    35,    36,    37,    38,    39,    -1,    -1,    42,
    -1,    -1,    21,    -1,    47,    -1,    -1,    -1,    51,    28,
    29,    30,    31,    32,    33,    34,    35,    36,    37,    38,
    39,    -1,    -1,    42,    -1,    -1,    -1,    -1,    47,    -1,
    -1,    50,    25,    26,    27,    28,    29,    28,    29,    -1,
    31,    32,    33,    34,    35,    38,    28,    29,    -1,    31,
    32,    33,    34,    35,    47,    20,    21,    -1,    -1,    -1,
    53,    -1,    53,    28,    29,    30,    31,    32,    33,    34,
    35,    36,    37,    38,    39,    20,    21,    42,    -1,    -1,
    -1,    -1,    47,    28,    29,    30,    31,    32,    33,    34,
    35,    36,    37,    38,    39,    20,    21,    42,    -1,    -1,
    -1,    -1,    47,    28,    29,    30,    31,    32,    33,    34,
    35,    36,    37,    38,    39,    20,    21,    42,    -1,    -1,
    -1,    -1,    47,    28,    29,    30,    31,    32,    33,    34,
    35,    36,    37,    38,    39,    20,    21,    42,    -1,    -1,
    -1,    -1,    47,    28,    29,    30,    31,    32,    33,    34,
    35,    36,    37,    38,    39,    20,    21,    42,    -1,    -1,
    -1,    -1,    47,    28,    29,    30,    31,    32,    33,    34,
    35,    36,    37,    38,    39,    -1,    21,    42,    -1,    24,
    -1,    -1,    47,    28,    29,    30,    31,    32,    33,    34,
    35,    36,    37,    38,    39,    -1,    21,    42,    -1,    24,
    -1,    -1,    47,    28,    29,    30,    31,    32,    33,    34,
    35,    36,    37,    38,    39,    -1,    21,    42,    -1,    24,
    -1,    -1,    47,    28,    29,    30,    31,    32,    33,    34,
    35,    36,    37,    38,    39,    -1,    21,    42,    -1,    24,
    -1,    -1,    47,    28,    29,    30,    31,    32,    33,    34,
    35,    36,    37,    38,    39,    -1,    -1,    42,    21,    22,
    -1,    -1,    47,    -1,    -1,    28,    29,    30,    31,    32,
    33,    34,    35,    36,    37,    38,    39,    -1,    -1,    42,
    21,    22,    -1,    -1,    47,    -1,    -1,    28,    29,    30,
    31,    32,    33,    34,    35,    36,    37,    38,    39,    -1,
    21,    42,    -1,    -1,    -1,    -1,    47,    28,    29,    30,
    31,    32,    33,    34,    35,    36,    37,    38,    39,    -1,
    -1,    42,    -1,    -1,    22,    -1,    47,    25,    26,    27,
    28,    29,    28,    29,    -1,    31,    32,    33,    34,    35,
    38,    28,    29,    -1,    31,    32,    33,    34,    35,    47,
    25,    26,    27,    28,    29,    28,    29,    53,    31,    32,
    33,    34,    35,    38,    -1,    -1,    53,    40,    41,    -1,
    -1,    -1,    47,    -1,    47,    -1,    -1,    50,    28,    29,
    30,    31,    32,    33,    34,    35,    36,    37,    38,    39,
    -1,    -1,    42,    -1,    -1,    28,    29,    47,    31,    32,
    33,    34,    35,    -1,    35,    -1,    -1,    40,    41,    40,
    41,    -1,    -1,    -1,    47,    -1,    47,    -1,    -1,    -1,
    -1,    -1,    53,    31,    32,    33,    34,    35,    36,    37,
    38,    39,    20,    -1,    42,    -1,    -1,    -1,    -1,    47,
    28,    29,    22,    31,    32,    33,    34,    35,    28,    29,
    -1,    31,    32,    33,    34,    35
};
/* -*-C-*-  Note some compilers choke on comments on `#line' lines.  */
#line 3 "/usr1/local/share/bison.simple"

/* Skeleton output parser for bison,
   Copyright (C) 1984, 1989, 1990 Free Software Foundation, Inc.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.  */

/* As a special exception, when this file is copied by Bison into a
   Bison output file, you may use that output file without restriction.
   This special exception was added by the Free Software Foundation
   in version 1.24 of Bison.  */

#ifndef alloca
#ifdef __GNUC__
#define alloca __builtin_alloca
#else /* not GNU C.  */
#if (!defined (__STDC__) && defined (sparc)) || defined (__sparc__) || defined (__sparc) || defined (__sgi)
#include <alloca.h>
#else /* not sparc */
#if defined (MSDOS) && !defined (__TURBOC__)
#include <malloc.h>
#else /* not MSDOS, or __TURBOC__ */
#if defined(_AIX)
#include <malloc.h>
 #pragma alloca
#else /* not MSDOS, __TURBOC__, or _AIX */
#ifdef __hpux
#ifdef __cplusplus
extern "C" {
void *alloca (unsigned int);
};
#else /* not __cplusplus */
void *alloca ();
#endif /* not __cplusplus */
#endif /* __hpux */
#endif /* not _AIX */
#endif /* not MSDOS, or __TURBOC__ */
#endif /* not sparc.  */
#endif /* not GNU C.  */
#endif /* alloca not defined.  */

/* This is the parser code that is written into each bison parser
  when the %semantic_parser declaration is not specified in the grammar.
  It was written by Richard Stallman by simplifying the hairy parser
  used when %semantic_parser is specified.  */

/* Note: there must be only one dollar sign in this file.
   It is replaced by the list of actions, each action
   as one case of the switch.  */

#define fferrok		(fferrstatus = 0)
#define ffclearin	(ffchar = FFEMPTY)
#define FFEMPTY		-2
#define FFEOF		0
#define FFACCEPT	return(0)
#define FFABORT 	return(1)
#define FFERROR		goto fferrlab1
/* Like FFERROR except do call fferror.
   This remains here temporarily to ease the
   transition to the new meaning of FFERROR, for GCC.
   Once GCC version 2 has supplanted version 1, this can go.  */
#define FFFAIL		goto fferrlab
#define FFRECOVERING()  (!!fferrstatus)
#define FFBACKUP(token, value) \
do								\
  if (ffchar == FFEMPTY && fflen == 1)				\
    { ffchar = (token), fflval = (value);			\
      ffchar1 = FFTRANSLATE (ffchar);				\
      FFPOPSTACK;						\
      goto ffbackup;						\
    }								\
  else								\
    { fferror ("syntax error: cannot back up"); FFERROR; }	\
while (0)

#define FFTERROR	1
#define FFERRCODE	256

#ifndef FFPURE
#define FFLEX		fflex()
#endif

#ifdef FFPURE
#ifdef FFLSP_NEEDED
#ifdef FFLEX_PARAM
#define FFLEX		fflex(&fflval, &fflloc, FFLEX_PARAM)
#else
#define FFLEX		fflex(&fflval, &fflloc)
#endif
#else /* not FFLSP_NEEDED */
#ifdef FFLEX_PARAM
#define FFLEX		fflex(&fflval, FFLEX_PARAM)
#else
#define FFLEX		fflex(&fflval)
#endif
#endif /* not FFLSP_NEEDED */
#endif

/* If nonreentrant, generate the variables here */

#ifndef FFPURE

int	ffchar;			/*  the lookahead symbol		*/
FFSTYPE	fflval;			/*  the semantic value of the		*/
				/*  lookahead symbol			*/

#ifdef FFLSP_NEEDED
FFLTYPE fflloc;			/*  location data for the lookahead	*/
				/*  symbol				*/
#endif

int ffnerrs;			/*  number of parse errors so far       */
#endif  /* not FFPURE */

#if FFDEBUG != 0
int ffdebug;			/*  nonzero means print parse trace	*/
/* Since this is uninitialized, it does not stop multiple parsers
   from coexisting.  */
#endif

/*  FFINITDEPTH indicates the initial size of the parser's stacks	*/

#ifndef	FFINITDEPTH
#define FFINITDEPTH 200
#endif

/*  FFMAXDEPTH is the maximum size the stacks can grow to
    (effective only if the built-in stack extension method is used).  */

#if FFMAXDEPTH == 0
#undef FFMAXDEPTH
#endif

#ifndef FFMAXDEPTH
#define FFMAXDEPTH 10000
#endif

/* Prevent warning if -Wstrict-prototypes.  */
#ifdef __GNUC__
int ffparse (void);
#endif

#if __GNUC__ > 1		/* GNU C and GNU C++ define this.  */
#define __ff_memcpy(TO,FROM,COUNT)	__builtin_memcpy(TO,FROM,COUNT)
#else				/* not GNU C or C++ */
#ifndef __cplusplus

/* This is the most reliable way to avoid incompatibilities
   in available built-in functions on various systems.  */
static void
__ff_memcpy (to, from, count)
     char *to;
     char *from;
     int count;
{
  register char *f = from;
  register char *t = to;
  register int i = count;

  while (i-- > 0)
    *t++ = *f++;
}

#else /* __cplusplus */

/* This is the most reliable way to avoid incompatibilities
   in available built-in functions on various systems.  */
static void
__ff_memcpy (char *to, char *from, int count)
{
  register char *f = from;
  register char *t = to;
  register int i = count;

  while (i-- > 0)
    *t++ = *f++;
}

#endif
#endif

#line 196 "/usr1/local/share/bison.simple"

/* The user can define FFPARSE_PARAM as the name of an argument to be passed
   into ffparse.  The argument should have type void *.
   It should actually point to an object.
   Grammar actions can access the variable by casting it
   to the proper pointer type.  */

#ifdef FFPARSE_PARAM
#ifdef __cplusplus
#define FFPARSE_PARAM_ARG void *FFPARSE_PARAM
#define FFPARSE_PARAM_DECL
#else /* not __cplusplus */
#define FFPARSE_PARAM_ARG FFPARSE_PARAM
#define FFPARSE_PARAM_DECL void *FFPARSE_PARAM;
#endif /* not __cplusplus */
#else /* not FFPARSE_PARAM */
#define FFPARSE_PARAM_ARG
#define FFPARSE_PARAM_DECL
#endif /* not FFPARSE_PARAM */

int
ffparse(FFPARSE_PARAM_ARG)
     FFPARSE_PARAM_DECL
{
  register int ffstate;
  register int ffn;
  register short *ffssp;
  register FFSTYPE *ffvsp;
  int fferrstatus;	/*  number of tokens to shift before error messages enabled */
  int ffchar1 = 0;		/*  lookahead token as an internal (translated) token number */

  short	ffssa[FFINITDEPTH];	/*  the state stack			*/
  FFSTYPE ffvsa[FFINITDEPTH];	/*  the semantic value stack		*/

  short *ffss = ffssa;		/*  refer to the stacks thru separate pointers */
  FFSTYPE *ffvs = ffvsa;	/*  to allow ffoverflow to reallocate them elsewhere */

#ifdef FFLSP_NEEDED
  FFLTYPE fflsa[FFINITDEPTH];	/*  the location stack			*/
  FFLTYPE *ffls = fflsa;
  FFLTYPE *fflsp;

#define FFPOPSTACK   (ffvsp--, ffssp--, fflsp--)
#else
#define FFPOPSTACK   (ffvsp--, ffssp--)
#endif

  int ffstacksize = FFINITDEPTH;

#ifdef FFPURE
  int ffchar;
  FFSTYPE fflval;
  int ffnerrs;
#ifdef FFLSP_NEEDED
  FFLTYPE fflloc;
#endif
#endif

  FFSTYPE ffval;		/*  the variable used to return		*/
				/*  semantic values from the action	*/
				/*  routines				*/

  int fflen;

#if FFDEBUG != 0
  if (ffdebug)
    fprintf(stderr, "Starting parse\n");
#endif

  ffstate = 0;
  fferrstatus = 0;
  ffnerrs = 0;
  ffchar = FFEMPTY;		/* Cause a token to be read.  */

  /* Initialize stack pointers.
     Waste one element of value and location stack
     so that they stay on the same level as the state stack.
     The wasted elements are never initialized.  */

  ffssp = ffss - 1;
  ffvsp = ffvs;
#ifdef FFLSP_NEEDED
  fflsp = ffls;
#endif

/* Push a new state, which is found in  ffstate  .  */
/* In all cases, when you get here, the value and location stacks
   have just been pushed. so pushing a state here evens the stacks.  */
ffnewstate:

  *++ffssp = ffstate;

  if (ffssp >= ffss + ffstacksize - 1)
    {
      /* Give user a chance to reallocate the stack */
      /* Use copies of these so that the &'s don't force the real ones into memory. */
      FFSTYPE *ffvs1 = ffvs;
      short *ffss1 = ffss;
#ifdef FFLSP_NEEDED
      FFLTYPE *ffls1 = ffls;
#endif

      /* Get the current used size of the three stacks, in elements.  */
      int size = ffssp - ffss + 1;

#ifdef ffoverflow
      /* Each stack pointer address is followed by the size of
	 the data in use in that stack, in bytes.  */
#ifdef FFLSP_NEEDED
      /* This used to be a conditional around just the two extra args,
	 but that might be undefined if ffoverflow is a macro.  */
      ffoverflow("parser stack overflow",
		 &ffss1, size * sizeof (*ffssp),
		 &ffvs1, size * sizeof (*ffvsp),
		 &ffls1, size * sizeof (*fflsp),
		 &ffstacksize);
#else
      ffoverflow("parser stack overflow",
		 &ffss1, size * sizeof (*ffssp),
		 &ffvs1, size * sizeof (*ffvsp),
		 &ffstacksize);
#endif

      ffss = ffss1; ffvs = ffvs1;
#ifdef FFLSP_NEEDED
      ffls = ffls1;
#endif
#else /* no ffoverflow */
      /* Extend the stack our own way.  */
      if (ffstacksize >= FFMAXDEPTH)
	{
	  fferror("parser stack overflow");
	  return 2;
	}
      ffstacksize *= 2;
      if (ffstacksize > FFMAXDEPTH)
	ffstacksize = FFMAXDEPTH;
      ffss = (short *) alloca (ffstacksize * sizeof (*ffssp));
      __ff_memcpy ((char *)ffss, (char *)ffss1, size * sizeof (*ffssp));
      ffvs = (FFSTYPE *) alloca (ffstacksize * sizeof (*ffvsp));
      __ff_memcpy ((char *)ffvs, (char *)ffvs1, size * sizeof (*ffvsp));
#ifdef FFLSP_NEEDED
      ffls = (FFLTYPE *) alloca (ffstacksize * sizeof (*fflsp));
      __ff_memcpy ((char *)ffls, (char *)ffls1, size * sizeof (*fflsp));
#endif
#endif /* no ffoverflow */

      ffssp = ffss + size - 1;
      ffvsp = ffvs + size - 1;
#ifdef FFLSP_NEEDED
      fflsp = ffls + size - 1;
#endif

#if FFDEBUG != 0
      if (ffdebug)
	fprintf(stderr, "Stack size increased to %d\n", ffstacksize);
#endif

      if (ffssp >= ffss + ffstacksize - 1)
	FFABORT;
    }

#if FFDEBUG != 0
  if (ffdebug)
    fprintf(stderr, "Entering state %d\n", ffstate);
#endif

  goto ffbackup;
 ffbackup:

/* Do appropriate processing given the current state.  */
/* Read a lookahead token if we need one and don't already have one.  */
/* ffresume: */

  /* First try to decide what to do without reference to lookahead token.  */

  ffn = ffpact[ffstate];
  if (ffn == FFFLAG)
    goto ffdefault;

  /* Not known => get a lookahead token if don't already have one.  */

  /* ffchar is either FFEMPTY or FFEOF
     or a valid token in external form.  */

  if (ffchar == FFEMPTY)
    {
#if FFDEBUG != 0
      if (ffdebug)
	fprintf(stderr, "Reading a token: ");
#endif
      ffchar = FFLEX;
    }

  /* Convert token to internal form (in ffchar1) for indexing tables with */

  if (ffchar <= 0)		/* This means end of input. */
    {
      ffchar1 = 0;
      ffchar = FFEOF;		/* Don't call FFLEX any more */

#if FFDEBUG != 0
      if (ffdebug)
	fprintf(stderr, "Now at end of input.\n");
#endif
    }
  else
    {
      ffchar1 = FFTRANSLATE(ffchar);

#if FFDEBUG != 0
      if (ffdebug)
	{
	  fprintf (stderr, "Next token is %d (%s", ffchar, fftname[ffchar1]);
	  /* Give the individual parser a way to print the precise meaning
	     of a token, for further debugging info.  */
#ifdef FFPRINT
	  FFPRINT (stderr, ffchar, fflval);
#endif
	  fprintf (stderr, ")\n");
	}
#endif
    }

  ffn += ffchar1;
  if (ffn < 0 || ffn > FFLAST || ffcheck[ffn] != ffchar1)
    goto ffdefault;

  ffn = fftable[ffn];

  /* ffn is what to do for this token type in this state.
     Negative => reduce, -ffn is rule number.
     Positive => shift, ffn is new state.
       New state is final state => don't bother to shift,
       just return success.
     0, or most negative number => error.  */

  if (ffn < 0)
    {
      if (ffn == FFFLAG)
	goto fferrlab;
      ffn = -ffn;
      goto ffreduce;
    }
  else if (ffn == 0)
    goto fferrlab;

  if (ffn == FFFINAL)
    FFACCEPT;

  /* Shift the lookahead token.  */

#if FFDEBUG != 0
  if (ffdebug)
    fprintf(stderr, "Shifting token %d (%s), ", ffchar, fftname[ffchar1]);
#endif

  /* Discard the token being shifted unless it is eof.  */
  if (ffchar != FFEOF)
    ffchar = FFEMPTY;

  *++ffvsp = fflval;
#ifdef FFLSP_NEEDED
  *++fflsp = fflloc;
#endif

  /* count tokens shifted since error; after three, turn off error status.  */
  if (fferrstatus) fferrstatus--;

  ffstate = ffn;
  goto ffnewstate;

/* Do the default action for the current state.  */
ffdefault:

  ffn = ffdefact[ffstate];
  if (ffn == 0)
    goto fferrlab;

/* Do a reduction.  ffn is the number of a rule to reduce with.  */
ffreduce:
  fflen = ffr2[ffn];
  if (fflen > 0)
    ffval = ffvsp[1-fflen]; /* implement default value of the action */

#if FFDEBUG != 0
  if (ffdebug)
    {
      int i;

      fprintf (stderr, "Reducing via rule %d (line %d), ",
	       ffn, ffrline[ffn]);

      /* Print the symbols being reduced, and their result.  */
      for (i = ffprhs[ffn]; ffrhs[i] > 0; i++)
	fprintf (stderr, "%s ", fftname[ffrhs[i]]);
      fprintf (stderr, " -> %s\n", fftname[ffr1[ffn]]);
    }
#endif


  switch (ffn) {

case 3:
#line 245 "eval.y"
{;
    break;}
case 4:
#line 247 "eval.y"
{ if( ffvsp[-1].Node<0 ) {
		     fferror("Couldn't build node structure: out of memory?");
		     FFERROR;  }
                  gParse.resultNode = ffvsp[-1].Node;
		;
    break;}
case 5:
#line 253 "eval.y"
{ if( ffvsp[-1].Node<0 ) {
		     fferror("Couldn't build node structure: out of memory?");
		     FFERROR;  }
                  gParse.resultNode = ffvsp[-1].Node;
		;
    break;}
case 6:
#line 259 "eval.y"
{ if( ffvsp[-1].Node<0 ) {
		     fferror("Couldn't build node structure: out of memory?");
		     FFERROR;  } 
                  gParse.resultNode = ffvsp[-1].Node;
		;
    break;}
case 7:
#line 265 "eval.y"
{ if( ffvsp[-1].Node<0 ) {
		     fferror("Couldn't build node structure: out of memory?");
		     FFERROR;  }
                  gParse.resultNode = ffvsp[-1].Node;
		;
    break;}
case 8:
#line 270 "eval.y"
{  fferrok;  ;
    break;}
case 9:
#line 274 "eval.y"
{ ffval.Node = New_Vector( ffvsp[0].Node ); TEST(ffval.Node); ;
    break;}
case 10:
#line 276 "eval.y"
{
                  if( gParse.Nodes[ffvsp[-2].Node].nSubNodes >= MAXSUBS ) {
		     ffvsp[-2].Node = Close_Vec( ffvsp[-2].Node ); TEST(ffvsp[-2].Node);
		     ffval.Node = New_Vector( ffvsp[-2].Node ); TEST(ffval.Node);
                  } else {
                     ffval.Node = ffvsp[-2].Node;
                  }
		  gParse.Nodes[ffval.Node].SubNodes[ gParse.Nodes[ffval.Node].nSubNodes++ ]
		     = ffvsp[0].Node;
                ;
    break;}
case 11:
#line 289 "eval.y"
{ ffval.Node = New_Vector( ffvsp[0].Node ); TEST(ffval.Node); ;
    break;}
case 12:
#line 291 "eval.y"
{
                  if( TYPE(ffvsp[-2].Node) < TYPE(ffvsp[0].Node) )
                     TYPE(ffvsp[-2].Node) = TYPE(ffvsp[0].Node);
                  if( gParse.Nodes[ffvsp[-2].Node].nSubNodes >= MAXSUBS ) {
		     ffvsp[-2].Node = Close_Vec( ffvsp[-2].Node ); TEST(ffvsp[-2].Node);
		     ffval.Node = New_Vector( ffvsp[-2].Node ); TEST(ffval.Node);
                  } else {
                     ffval.Node = ffvsp[-2].Node;
                  }
		  gParse.Nodes[ffval.Node].SubNodes[ gParse.Nodes[ffval.Node].nSubNodes++ ]
		     = ffvsp[0].Node;
                ;
    break;}
case 13:
#line 304 "eval.y"
{
                  if( gParse.Nodes[ffvsp[-2].Node].nSubNodes >= MAXSUBS ) {
		     ffvsp[-2].Node = Close_Vec( ffvsp[-2].Node ); TEST(ffvsp[-2].Node);
		     ffval.Node = New_Vector( ffvsp[-2].Node ); TEST(ffval.Node);
                  } else {
                     ffval.Node = ffvsp[-2].Node;
                  }
		  gParse.Nodes[ffval.Node].SubNodes[ gParse.Nodes[ffval.Node].nSubNodes++ ]
		     = ffvsp[0].Node;
                ;
    break;}
case 14:
#line 315 "eval.y"
{
                  TYPE(ffvsp[-2].Node) = TYPE(ffvsp[0].Node);
                  if( gParse.Nodes[ffvsp[-2].Node].nSubNodes >= MAXSUBS ) {
		     ffvsp[-2].Node = Close_Vec( ffvsp[-2].Node ); TEST(ffvsp[-2].Node);
		     ffval.Node = New_Vector( ffvsp[-2].Node ); TEST(ffval.Node);
                  } else {
                     ffval.Node = ffvsp[-2].Node;
                  }
		  gParse.Nodes[ffval.Node].SubNodes[ gParse.Nodes[ffval.Node].nSubNodes++ ]
		     = ffvsp[0].Node;
                ;
    break;}
case 15:
#line 329 "eval.y"
{ ffval.Node = Close_Vec( ffvsp[-1].Node ); TEST(ffval.Node); ;
    break;}
case 16:
#line 333 "eval.y"
{ ffval.Node = Close_Vec( ffvsp[-1].Node ); TEST(ffval.Node); ;
    break;}
case 17:
#line 337 "eval.y"
{
                  ffval.Node = New_Const( BITSTR, ffvsp[0].str, strlen(ffvsp[0].str)+1 ); TEST(ffval.Node);
		  SIZE(ffval.Node) = strlen(ffvsp[0].str); ;
    break;}
case 18:
#line 341 "eval.y"
{ ffval.Node = New_Column( ffvsp[0].lng ); TEST(ffval.Node); ;
    break;}
case 19:
#line 343 "eval.y"
{
                  if( TYPE(ffvsp[-1].Node) != LONG
		      || OPER(ffvsp[-1].Node) != CONST_OP ) {
		     fferror("Offset argument must be a constant integer");
		     FFERROR;
		  }
                  ffval.Node = New_Offset( ffvsp[-3].lng, ffvsp[-1].Node ); TEST(ffval.Node);
                ;
    break;}
case 20:
#line 352 "eval.y"
{ ffval.Node = New_BinOp( BITSTR, ffvsp[-2].Node, '&', ffvsp[0].Node ); TEST(ffval.Node);
                  SIZE(ffval.Node) = ( SIZE(ffvsp[-2].Node)>SIZE(ffvsp[0].Node) ? SIZE(ffvsp[-2].Node) : SIZE(ffvsp[0].Node) );  ;
    break;}
case 21:
#line 355 "eval.y"
{ ffval.Node = New_BinOp( BITSTR, ffvsp[-2].Node, '|', ffvsp[0].Node ); TEST(ffval.Node);
                  SIZE(ffval.Node) = ( SIZE(ffvsp[-2].Node)>SIZE(ffvsp[0].Node) ? SIZE(ffvsp[-2].Node) : SIZE(ffvsp[0].Node) );  ;
    break;}
case 22:
#line 358 "eval.y"
{ 
		  if (SIZE(ffvsp[-2].Node)+SIZE(ffvsp[0].Node) >= MAX_STRLEN) {
		    fferror("Combined bit string size exceeds " MAX_STRLEN_S " bits");
		    FFERROR;
		  }
		  ffval.Node = New_BinOp( BITSTR, ffvsp[-2].Node, '+', ffvsp[0].Node ); TEST(ffval.Node);
                  SIZE(ffval.Node) = SIZE(ffvsp[-2].Node) + SIZE(ffvsp[0].Node); 
		;
    break;}
case 23:
#line 367 "eval.y"
{ ffval.Node = New_Deref( ffvsp[-3].Node, 1, ffvsp[-1].Node,  0,  0,  0,   0 ); TEST(ffval.Node); ;
    break;}
case 24:
#line 369 "eval.y"
{ ffval.Node = New_Deref( ffvsp[-5].Node, 2, ffvsp[-3].Node, ffvsp[-1].Node,  0,  0,   0 ); TEST(ffval.Node); ;
    break;}
case 25:
#line 371 "eval.y"
{ ffval.Node = New_Deref( ffvsp[-7].Node, 3, ffvsp[-5].Node, ffvsp[-3].Node, ffvsp[-1].Node,  0,   0 ); TEST(ffval.Node); ;
    break;}
case 26:
#line 373 "eval.y"
{ ffval.Node = New_Deref( ffvsp[-9].Node, 4, ffvsp[-7].Node, ffvsp[-5].Node, ffvsp[-3].Node, ffvsp[-1].Node,   0 ); TEST(ffval.Node); ;
    break;}
case 27:
#line 375 "eval.y"
{ ffval.Node = New_Deref( ffvsp[-11].Node, 5, ffvsp[-9].Node, ffvsp[-7].Node, ffvsp[-5].Node, ffvsp[-3].Node, ffvsp[-1].Node ); TEST(ffval.Node); ;
    break;}
case 28:
#line 377 "eval.y"
{ ffval.Node = New_Unary( BITSTR, NOT, ffvsp[0].Node ); TEST(ffval.Node);     ;
    break;}
case 29:
#line 380 "eval.y"
{ ffval.Node = ffvsp[-1].Node; ;
    break;}
case 30:
#line 384 "eval.y"
{ ffval.Node = New_Const( LONG,   &(ffvsp[0].lng), sizeof(long)   ); TEST(ffval.Node); ;
    break;}
case 31:
#line 386 "eval.y"
{ ffval.Node = New_Const( DOUBLE, &(ffvsp[0].dbl), sizeof(double) ); TEST(ffval.Node); ;
    break;}
case 32:
#line 388 "eval.y"
{ ffval.Node = New_Column( ffvsp[0].lng ); TEST(ffval.Node); ;
    break;}
case 33:
#line 390 "eval.y"
{
                  if( TYPE(ffvsp[-1].Node) != LONG
		      || OPER(ffvsp[-1].Node) != CONST_OP ) {
		     fferror("Offset argument must be a constant integer");
		     FFERROR;
		  }
                  ffval.Node = New_Offset( ffvsp[-3].lng, ffvsp[-1].Node ); TEST(ffval.Node);
                ;
    break;}
case 34:
#line 399 "eval.y"
{ ffval.Node = New_Func( LONG, row_fct,  0, 0, 0, 0, 0, 0, 0, 0 ); ;
    break;}
case 35:
#line 401 "eval.y"
{ ffval.Node = New_Func( LONG, null_fct, 0, 0, 0, 0, 0, 0, 0, 0 ); ;
    break;}
case 36:
#line 403 "eval.y"
{ PROMOTE(ffvsp[-2].Node,ffvsp[0].Node); ffval.Node = New_BinOp( TYPE(ffvsp[-2].Node), ffvsp[-2].Node, '%', ffvsp[0].Node );
		  TEST(ffval.Node);                                                ;
    break;}
case 37:
#line 406 "eval.y"
{ PROMOTE(ffvsp[-2].Node,ffvsp[0].Node); ffval.Node = New_BinOp( TYPE(ffvsp[-2].Node), ffvsp[-2].Node, '+', ffvsp[0].Node );
		  TEST(ffval.Node);                                                ;
    break;}
case 38:
#line 409 "eval.y"
{ PROMOTE(ffvsp[-2].Node,ffvsp[0].Node); ffval.Node = New_BinOp( TYPE(ffvsp[-2].Node), ffvsp[-2].Node, '-', ffvsp[0].Node ); 
		  TEST(ffval.Node);                                                ;
    break;}
case 39:
#line 412 "eval.y"
{ PROMOTE(ffvsp[-2].Node,ffvsp[0].Node); ffval.Node = New_BinOp( TYPE(ffvsp[-2].Node), ffvsp[-2].Node, '*', ffvsp[0].Node ); 
		  TEST(ffval.Node);                                                ;
    break;}
case 40:
#line 415 "eval.y"
{ PROMOTE(ffvsp[-2].Node,ffvsp[0].Node); ffval.Node = New_BinOp( TYPE(ffvsp[-2].Node), ffvsp[-2].Node, '/', ffvsp[0].Node ); 
		  TEST(ffval.Node);                                                ;
    break;}
case 41:
#line 418 "eval.y"
{ PROMOTE(ffvsp[-2].Node,ffvsp[0].Node); ffval.Node = New_BinOp( TYPE(ffvsp[-2].Node), ffvsp[-2].Node, POWER, ffvsp[0].Node );
		  TEST(ffval.Node);                                                ;
    break;}
case 42:
#line 421 "eval.y"
{ ffval.Node = ffvsp[0].Node; ;
    break;}
case 43:
#line 423 "eval.y"
{ ffval.Node = New_Unary( TYPE(ffvsp[0].Node), UMINUS, ffvsp[0].Node ); TEST(ffval.Node); ;
    break;}
case 44:
#line 425 "eval.y"
{ ffval.Node = ffvsp[-1].Node; ;
    break;}
case 45:
#line 427 "eval.y"
{ ffvsp[0].Node = New_Unary( TYPE(ffvsp[-2].Node), 0, ffvsp[0].Node );
                  ffval.Node = New_BinOp( TYPE(ffvsp[-2].Node), ffvsp[-2].Node, '*', ffvsp[0].Node ); 
		  TEST(ffval.Node);                                ;
    break;}
case 46:
#line 431 "eval.y"
{ ffvsp[-2].Node = New_Unary( TYPE(ffvsp[0].Node), 0, ffvsp[-2].Node );
                  ffval.Node = New_BinOp( TYPE(ffvsp[0].Node), ffvsp[-2].Node, '*', ffvsp[0].Node );
                  TEST(ffval.Node);                                ;
    break;}
case 47:
#line 435 "eval.y"
{
                  PROMOTE(ffvsp[-2].Node,ffvsp[0].Node);
                  if( ! Test_Dims(ffvsp[-2].Node,ffvsp[0].Node) ) {
                     fferror("Incompatible dimensions in '?:' arguments");
		     FFERROR;
                  }
                  ffval.Node = New_Func( 0, ifthenelse_fct, 3, ffvsp[-2].Node, ffvsp[0].Node, ffvsp[-4].Node,
                                 0, 0, 0, 0 );
                  TEST(ffval.Node);
                  if( SIZE(ffvsp[-2].Node)<SIZE(ffvsp[0].Node) )  Copy_Dims(ffval.Node, ffvsp[0].Node);
                  TYPE(ffvsp[-4].Node) = TYPE(ffvsp[-2].Node);
                  if( ! Test_Dims(ffvsp[-4].Node,ffval.Node) ) {
                     fferror("Incompatible dimensions in '?:' condition");
		     FFERROR;
                  }
                  TYPE(ffvsp[-4].Node) = BOOLEAN;
                  if( SIZE(ffval.Node)<SIZE(ffvsp[-4].Node) )  Copy_Dims(ffval.Node, ffvsp[-4].Node);
                ;
    break;}
case 48:
#line 454 "eval.y"
{
                  PROMOTE(ffvsp[-2].Node,ffvsp[0].Node);
                  if( ! Test_Dims(ffvsp[-2].Node,ffvsp[0].Node) ) {
                     fferror("Incompatible dimensions in '?:' arguments");
		     FFERROR;
                  }
                  ffval.Node = New_Func( 0, ifthenelse_fct, 3, ffvsp[-2].Node, ffvsp[0].Node, ffvsp[-4].Node,
                                 0, 0, 0, 0 );
                  TEST(ffval.Node);
                  if( SIZE(ffvsp[-2].Node)<SIZE(ffvsp[0].Node) )  Copy_Dims(ffval.Node, ffvsp[0].Node);
                  TYPE(ffvsp[-4].Node) = TYPE(ffvsp[-2].Node);
                  if( ! Test_Dims(ffvsp[-4].Node,ffval.Node) ) {
                     fferror("Incompatible dimensions in '?:' condition");
		     FFERROR;
                  }
                  TYPE(ffvsp[-4].Node) = BOOLEAN;
                  if( SIZE(ffval.Node)<SIZE(ffvsp[-4].Node) )  Copy_Dims(ffval.Node, ffvsp[-4].Node);
                ;
    break;}
case 49:
#line 473 "eval.y"
{
                  PROMOTE(ffvsp[-2].Node,ffvsp[0].Node);
                  if( ! Test_Dims(ffvsp[-2].Node,ffvsp[0].Node) ) {
                     fferror("Incompatible dimensions in '?:' arguments");
		     FFERROR;
                  }
                  ffval.Node = New_Func( 0, ifthenelse_fct, 3, ffvsp[-2].Node, ffvsp[0].Node, ffvsp[-4].Node,
                                 0, 0, 0, 0 );
                  TEST(ffval.Node);
                  if( SIZE(ffvsp[-2].Node)<SIZE(ffvsp[0].Node) )  Copy_Dims(ffval.Node, ffvsp[0].Node);
                  TYPE(ffvsp[-4].Node) = TYPE(ffvsp[-2].Node);
                  if( ! Test_Dims(ffvsp[-4].Node,ffval.Node) ) {
                     fferror("Incompatible dimensions in '?:' condition");
		     FFERROR;
                  }
                  TYPE(ffvsp[-4].Node) = BOOLEAN;
                  if( SIZE(ffval.Node)<SIZE(ffvsp[-4].Node) )  Copy_Dims(ffval.Node, ffvsp[-4].Node);
                ;
    break;}
case 50:
#line 492 "eval.y"
{ if (FSTRCMP(ffvsp[-1].str,"RANDOM(") == 0) {  /* Scalar RANDOM() */
                     srand( (unsigned int) time(NULL) );
                     ffval.Node = New_Func( DOUBLE, rnd_fct, 0, 0, 0, 0, 0, 0, 0, 0 );
		  } else if (FSTRCMP(ffvsp[-1].str,"RANDOMN(") == 0) {/*Scalar RANDOMN()*/
		     srand( (unsigned int) time(NULL) );
		     ffval.Node = New_Func( DOUBLE, gasrnd_fct, 0, 0, 0, 0, 0, 0, 0, 0 );
                  } else {
                     fferror("Function() not supported");
		     FFERROR;
		  }
                  TEST(ffval.Node); 
                ;
    break;}
case 51:
#line 505 "eval.y"
{ if (FSTRCMP(ffvsp[-2].str,"SUM(") == 0) {
		     ffval.Node = New_Func( LONG, sum_fct, 1, ffvsp[-1].Node, 0, 0, 0, 0, 0, 0 );
                  } else if (FSTRCMP(ffvsp[-2].str,"NELEM(") == 0) {
                     ffval.Node = New_Const( LONG, &( SIZE(ffvsp[-1].Node) ), sizeof(long) );
                  } else if (FSTRCMP(ffvsp[-2].str,"ACCUM(") == 0) {
		    long zero = 0;
		    ffval.Node = New_BinOp( LONG , ffvsp[-1].Node, ACCUM, New_Const( LONG, &zero, sizeof(zero) ));
		  } else {
                     fferror("Function(bool) not supported");
		     FFERROR;
		  }
                  TEST(ffval.Node); 
		;
    break;}
case 52:
#line 519 "eval.y"
{ if (FSTRCMP(ffvsp[-2].str,"NELEM(") == 0) {
                     ffval.Node = New_Const( LONG, &( SIZE(ffvsp[-1].Node) ), sizeof(long) );
		  } else if (FSTRCMP(ffvsp[-2].str,"NVALID(") == 0) {
		     ffval.Node = New_Func( LONG, nonnull_fct, 1, ffvsp[-1].Node,
				    0, 0, 0, 0, 0, 0 );
		  } else {
                     fferror("Function(str) not supported");
		     FFERROR;
		  }
                  TEST(ffval.Node); 
		;
    break;}
case 53:
#line 531 "eval.y"
{ if (FSTRCMP(ffvsp[-2].str,"NELEM(") == 0) {
                     ffval.Node = New_Const( LONG, &( SIZE(ffvsp[-1].Node) ), sizeof(long) );
		} else if (FSTRCMP(ffvsp[-2].str,"NVALID(") == 0) { /* Bit arrays do not have NULL */
                     ffval.Node = New_Const( LONG, &( SIZE(ffvsp[-1].Node) ), sizeof(long) );
		} else if (FSTRCMP(ffvsp[-2].str,"SUM(") == 0) {
		     ffval.Node = New_Func( LONG, sum_fct, 1, ffvsp[-1].Node,
				    0, 0, 0, 0, 0, 0 );
		} else if (FSTRCMP(ffvsp[-2].str,"MIN(") == 0) {
		     ffval.Node = New_Func( TYPE(ffvsp[-1].Node),  /* Force 1D result */
				    min1_fct, 1, ffvsp[-1].Node, 0, 0, 0, 0, 0, 0 );
		     /* Note: $2 is a vector so the result can never
		        be a constant.  Therefore it will never be set
		        inside New_Func(), and it is safe to set SIZE() */
		     SIZE(ffval.Node) = 1;
		} else if (FSTRCMP(ffvsp[-2].str,"ACCUM(") == 0) {
		    long zero = 0;
		    ffval.Node = New_BinOp( LONG , ffvsp[-1].Node, ACCUM, New_Const( LONG, &zero, sizeof(zero) ));
		} else if (FSTRCMP(ffvsp[-2].str,"MAX(") == 0) {
		     ffval.Node = New_Func( TYPE(ffvsp[-1].Node),  /* Force 1D result */
				    max1_fct, 1, ffvsp[-1].Node, 0, 0, 0, 0, 0, 0 );
		     /* Note: $2 is a vector so the result can never
		        be a constant.  Therefore it will never be set
		        inside New_Func(), and it is safe to set SIZE() */
		     SIZE(ffval.Node) = 1;
		} else {
                     fferror("Function(bits) not supported");
		     FFERROR;
		  }
                  TEST(ffval.Node); 
		;
    break;}
case 54:
#line 562 "eval.y"
{ if (FSTRCMP(ffvsp[-2].str,"SUM(") == 0)
		     ffval.Node = New_Func( TYPE(ffvsp[-1].Node), sum_fct, 1, ffvsp[-1].Node,
				    0, 0, 0, 0, 0, 0 );
		  else if (FSTRCMP(ffvsp[-2].str,"AVERAGE(") == 0)
		     ffval.Node = New_Func( DOUBLE, average_fct, 1, ffvsp[-1].Node,
				    0, 0, 0, 0, 0, 0 );
		  else if (FSTRCMP(ffvsp[-2].str,"STDDEV(") == 0)
		     ffval.Node = New_Func( DOUBLE, stddev_fct, 1, ffvsp[-1].Node,
				    0, 0, 0, 0, 0, 0 );
		  else if (FSTRCMP(ffvsp[-2].str,"MEDIAN(") == 0)
		     ffval.Node = New_Func( TYPE(ffvsp[-1].Node), median_fct, 1, ffvsp[-1].Node,
				    0, 0, 0, 0, 0, 0 );
		  else if (FSTRCMP(ffvsp[-2].str,"NELEM(") == 0)
                     ffval.Node = New_Const( LONG, &( SIZE(ffvsp[-1].Node) ), sizeof(long) );
		  else if (FSTRCMP(ffvsp[-2].str,"NVALID(") == 0)
		     ffval.Node = New_Func( LONG, nonnull_fct, 1, ffvsp[-1].Node,
				    0, 0, 0, 0, 0, 0 );
		  else if   ((FSTRCMP(ffvsp[-2].str,"ACCUM(") == 0) && (TYPE(ffvsp[-1].Node) == LONG)) {
		    long zero = 0;
		    ffval.Node = New_BinOp( LONG ,   ffvsp[-1].Node, ACCUM, New_Const( LONG,   &zero, sizeof(zero) ));
		  } else if ((FSTRCMP(ffvsp[-2].str,"ACCUM(") == 0) && (TYPE(ffvsp[-1].Node) == DOUBLE)) {
		    double zero = 0;
		    ffval.Node = New_BinOp( DOUBLE , ffvsp[-1].Node, ACCUM, New_Const( DOUBLE, &zero, sizeof(zero) ));
		  } else if ((FSTRCMP(ffvsp[-2].str,"SEQDIFF(") == 0) && (TYPE(ffvsp[-1].Node) == LONG)) {
		    long zero = 0;
		    ffval.Node = New_BinOp( LONG ,   ffvsp[-1].Node, DIFF, New_Const( LONG,   &zero, sizeof(zero) ));
		  } else if ((FSTRCMP(ffvsp[-2].str,"SEQDIFF(") == 0) && (TYPE(ffvsp[-1].Node) == DOUBLE)) {
		    double zero = 0;
		    ffval.Node = New_BinOp( DOUBLE , ffvsp[-1].Node, DIFF, New_Const( DOUBLE, &zero, sizeof(zero) ));
		  } else if (FSTRCMP(ffvsp[-2].str,"ABS(") == 0)
		     ffval.Node = New_Func( 0, abs_fct, 1, ffvsp[-1].Node, 0, 0, 0, 0, 0, 0 );
 		  else if (FSTRCMP(ffvsp[-2].str,"MIN(") == 0)
		     ffval.Node = New_Func( TYPE(ffvsp[-1].Node),  /* Force 1D result */
				    min1_fct, 1, ffvsp[-1].Node, 0, 0, 0, 0, 0, 0 );
		  else if (FSTRCMP(ffvsp[-2].str,"MAX(") == 0)
		     ffval.Node = New_Func( TYPE(ffvsp[-1].Node),  /* Force 1D result */
				    max1_fct, 1, ffvsp[-1].Node, 0, 0, 0, 0, 0, 0 );
		  else if (FSTRCMP(ffvsp[-2].str,"RANDOM(") == 0) { /* Vector RANDOM() */
                     srand( (unsigned int) time(NULL) );
                     ffval.Node = New_Func( 0, rnd_fct, 1, ffvsp[-1].Node, 0, 0, 0, 0, 0, 0 );
		     TEST(ffval.Node);
		     TYPE(ffval.Node) = DOUBLE;
		  } else if (FSTRCMP(ffvsp[-2].str,"RANDOMN(") == 0) {
		     srand( (unsigned int) time(NULL) ); /* Vector RANDOMN() */
		     ffval.Node = New_Func( 0, gasrnd_fct, 1, ffvsp[-1].Node, 0, 0, 0, 0, 0, 0 );
		     TEST(ffval.Node);
		     TYPE(ffval.Node) = DOUBLE;
                  } 
  		  else {  /*  These all take DOUBLE arguments  */
		     if( TYPE(ffvsp[-1].Node) != DOUBLE ) ffvsp[-1].Node = New_Unary( DOUBLE, 0, ffvsp[-1].Node );
                     if (FSTRCMP(ffvsp[-2].str,"SIN(") == 0)
			ffval.Node = New_Func( 0, sin_fct,  1, ffvsp[-1].Node, 0, 0, 0, 0, 0, 0 );
		     else if (FSTRCMP(ffvsp[-2].str,"COS(") == 0)
			ffval.Node = New_Func( 0, cos_fct,  1, ffvsp[-1].Node, 0, 0, 0, 0, 0, 0 );
		     else if (FSTRCMP(ffvsp[-2].str,"TAN(") == 0)
			ffval.Node = New_Func( 0, tan_fct,  1, ffvsp[-1].Node, 0, 0, 0, 0, 0, 0 );
		     else if (FSTRCMP(ffvsp[-2].str,"ARCSIN(") == 0
			      || FSTRCMP(ffvsp[-2].str,"ASIN(") == 0)
			ffval.Node = New_Func( 0, asin_fct, 1, ffvsp[-1].Node, 0, 0, 0, 0, 0, 0 );
		     else if (FSTRCMP(ffvsp[-2].str,"ARCCOS(") == 0
			      || FSTRCMP(ffvsp[-2].str,"ACOS(") == 0)
			ffval.Node = New_Func( 0, acos_fct, 1, ffvsp[-1].Node, 0, 0, 0, 0, 0, 0 );
		     else if (FSTRCMP(ffvsp[-2].str,"ARCTAN(") == 0
			      || FSTRCMP(ffvsp[-2].str,"ATAN(") == 0)
			ffval.Node = New_Func( 0, atan_fct, 1, ffvsp[-1].Node, 0, 0, 0, 0, 0, 0 );
		     else if (FSTRCMP(ffvsp[-2].str,"SINH(") == 0)
			ffval.Node = New_Func( 0, sinh_fct,  1, ffvsp[-1].Node, 0, 0, 0, 0, 0, 0 );
		     else if (FSTRCMP(ffvsp[-2].str,"COSH(") == 0)
			ffval.Node = New_Func( 0, cosh_fct,  1, ffvsp[-1].Node, 0, 0, 0, 0, 0, 0 );
		     else if (FSTRCMP(ffvsp[-2].str,"TANH(") == 0)
			ffval.Node = New_Func( 0, tanh_fct,  1, ffvsp[-1].Node, 0, 0, 0, 0, 0, 0 );
		     else if (FSTRCMP(ffvsp[-2].str,"EXP(") == 0)
			ffval.Node = New_Func( 0, exp_fct,  1, ffvsp[-1].Node, 0, 0, 0, 0, 0, 0 );
		     else if (FSTRCMP(ffvsp[-2].str,"LOG(") == 0)
			ffval.Node = New_Func( 0, log_fct,  1, ffvsp[-1].Node, 0, 0, 0, 0, 0, 0 );
		     else if (FSTRCMP(ffvsp[-2].str,"LOG10(") == 0)
			ffval.Node = New_Func( 0, log10_fct, 1, ffvsp[-1].Node, 0, 0, 0, 0, 0, 0 );
		     else if (FSTRCMP(ffvsp[-2].str,"SQRT(") == 0)
			ffval.Node = New_Func( 0, sqrt_fct, 1, ffvsp[-1].Node, 0, 0, 0, 0, 0, 0 );
		     else if (FSTRCMP(ffvsp[-2].str,"ROUND(") == 0)
			ffval.Node = New_Func( 0, round_fct, 1, ffvsp[-1].Node, 0, 0, 0, 0, 0, 0 );
		     else if (FSTRCMP(ffvsp[-2].str,"FLOOR(") == 0)
			ffval.Node = New_Func( 0, floor_fct, 1, ffvsp[-1].Node, 0, 0, 0, 0, 0, 0 );
		     else if (FSTRCMP(ffvsp[-2].str,"CEIL(") == 0)
			ffval.Node = New_Func( 0, ceil_fct, 1, ffvsp[-1].Node, 0, 0, 0, 0, 0, 0 );
		     else if (FSTRCMP(ffvsp[-2].str,"RANDOMP(") == 0) {
		       srand( (unsigned int) time(NULL) );
		       ffval.Node = New_Func( 0, poirnd_fct, 1, ffvsp[-1].Node, 
				      0, 0, 0, 0, 0, 0 );
		       TYPE(ffval.Node) = LONG;
		     } else {
			fferror("Function(expr) not supported");
			FFERROR;
		     }
		  }
                  TEST(ffval.Node); 
                ;
    break;}
case 55:
#line 660 "eval.y"
{ 
		  if (FSTRCMP(ffvsp[-4].str,"STRSTR(") == 0) {
		    ffval.Node = New_Func( LONG, strpos_fct, 2, ffvsp[-3].Node, ffvsp[-1].Node, 0, 
				   0, 0, 0, 0 );
		    TEST(ffval.Node);
		  }
                ;
    break;}
case 56:
#line 668 "eval.y"
{ 
		   if (FSTRCMP(ffvsp[-4].str,"DEFNULL(") == 0) {
		      if( SIZE(ffvsp[-3].Node)>=SIZE(ffvsp[-1].Node) && Test_Dims( ffvsp[-3].Node, ffvsp[-1].Node ) ) {
			 PROMOTE(ffvsp[-3].Node,ffvsp[-1].Node);
			 ffval.Node = New_Func( 0, defnull_fct, 2, ffvsp[-3].Node, ffvsp[-1].Node, 0,
					0, 0, 0, 0 );
			 TEST(ffval.Node); 
		      } else {
			 fferror("Dimensions of DEFNULL arguments "
				 "are not compatible");
			 FFERROR;
		      }
		   } else if (FSTRCMP(ffvsp[-4].str,"ARCTAN2(") == 0) {
		     if( TYPE(ffvsp[-3].Node) != DOUBLE ) ffvsp[-3].Node = New_Unary( DOUBLE, 0, ffvsp[-3].Node );
		     if( TYPE(ffvsp[-1].Node) != DOUBLE ) ffvsp[-1].Node = New_Unary( DOUBLE, 0, ffvsp[-1].Node );
		     if( Test_Dims( ffvsp[-3].Node, ffvsp[-1].Node ) ) {
			ffval.Node = New_Func( 0, atan2_fct, 2, ffvsp[-3].Node, ffvsp[-1].Node, 0, 0, 0, 0, 0 );
			TEST(ffval.Node); 
			if( SIZE(ffvsp[-3].Node)<SIZE(ffvsp[-1].Node) ) Copy_Dims(ffval.Node, ffvsp[-1].Node);
		     } else {
			fferror("Dimensions of arctan2 arguments "
				"are not compatible");
			FFERROR;
		     }
		   } else if (FSTRCMP(ffvsp[-4].str,"MIN(") == 0) {
		      PROMOTE( ffvsp[-3].Node, ffvsp[-1].Node );
		      if( Test_Dims( ffvsp[-3].Node, ffvsp[-1].Node ) ) {
			ffval.Node = New_Func( 0, min2_fct, 2, ffvsp[-3].Node, ffvsp[-1].Node, 0, 0, 0, 0, 0 );
			TEST(ffval.Node);
			if( SIZE(ffvsp[-3].Node)<SIZE(ffvsp[-1].Node) ) Copy_Dims(ffval.Node, ffvsp[-1].Node);
		      } else {
			fferror("Dimensions of min(a,b) arguments "
				"are not compatible");
			FFERROR;
		      }
		   } else if (FSTRCMP(ffvsp[-4].str,"MAX(") == 0) {
		      PROMOTE( ffvsp[-3].Node, ffvsp[-1].Node );
		      if( Test_Dims( ffvsp[-3].Node, ffvsp[-1].Node ) ) {
			ffval.Node = New_Func( 0, max2_fct, 2, ffvsp[-3].Node, ffvsp[-1].Node, 0, 0, 0, 0, 0 );
			TEST(ffval.Node);
			if( SIZE(ffvsp[-3].Node)<SIZE(ffvsp[-1].Node) ) Copy_Dims(ffval.Node, ffvsp[-1].Node);
		      } else {
			fferror("Dimensions of max(a,b) arguments "
				"are not compatible");
			FFERROR;
		      }
#if 0
		   } else if (FSTRCMP(ffvsp[-4].str,"STRSTR(") == 0) {
		     if( TYPE(ffvsp[-3].Node) != STRING || TYPE(ffvsp[-1].Node) != STRING) {
		       fferror("Arguments to strstr(s,r) must be strings");
		       FFERROR;
		     }
		     ffval.Node = New_Func( LONG, strpos_fct, 2, ffvsp[-3].Node, ffvsp[-1].Node, 0, 
				    0, 0, 0, 0 );
		     TEST(ffval.Node);
#endif
		   } else {
		      fferror("Function(expr,expr) not supported");
		      FFERROR;
		   }
                ;
    break;}
case 57:
#line 730 "eval.y"
{ 
		  if (FSTRCMP(ffvsp[-8].str,"ANGSEP(") == 0) {
		    if( TYPE(ffvsp[-7].Node) != DOUBLE ) ffvsp[-7].Node = New_Unary( DOUBLE, 0, ffvsp[-7].Node );
		    if( TYPE(ffvsp[-5].Node) != DOUBLE ) ffvsp[-5].Node = New_Unary( DOUBLE, 0, ffvsp[-5].Node );
		    if( TYPE(ffvsp[-3].Node) != DOUBLE ) ffvsp[-3].Node = New_Unary( DOUBLE, 0, ffvsp[-3].Node );
		    if( TYPE(ffvsp[-1].Node) != DOUBLE ) ffvsp[-1].Node = New_Unary( DOUBLE, 0, ffvsp[-1].Node );
		    if( Test_Dims( ffvsp[-7].Node, ffvsp[-5].Node ) && Test_Dims( ffvsp[-5].Node, ffvsp[-3].Node ) && 
			Test_Dims( ffvsp[-3].Node, ffvsp[-1].Node ) ) {
		      ffval.Node = New_Func( 0, angsep_fct, 4, ffvsp[-7].Node, ffvsp[-5].Node, ffvsp[-3].Node, ffvsp[-1].Node,0,0,0 );
		      TEST(ffval.Node); 
		      if( SIZE(ffvsp[-7].Node)<SIZE(ffvsp[-5].Node) ) Copy_Dims(ffval.Node, ffvsp[-5].Node);
		      if( SIZE(ffvsp[-5].Node)<SIZE(ffvsp[-3].Node) ) Copy_Dims(ffval.Node, ffvsp[-3].Node);
		      if( SIZE(ffvsp[-3].Node)<SIZE(ffvsp[-1].Node) ) Copy_Dims(ffval.Node, ffvsp[-1].Node);
		    } else {
		      fferror("Dimensions of ANGSEP arguments "
			      "are not compatible");
		      FFERROR;
		    }
		   } else {
		      fferror("Function(expr,expr,expr,expr) not supported");
		      FFERROR;
		   }
                ;
    break;}
case 58:
#line 754 "eval.y"
{ ffval.Node = New_Deref( ffvsp[-3].Node, 1, ffvsp[-1].Node,  0,  0,  0,   0 ); TEST(ffval.Node); ;
    break;}
case 59:
#line 756 "eval.y"
{ ffval.Node = New_Deref( ffvsp[-5].Node, 2, ffvsp[-3].Node, ffvsp[-1].Node,  0,  0,   0 ); TEST(ffval.Node); ;
    break;}
case 60:
#line 758 "eval.y"
{ ffval.Node = New_Deref( ffvsp[-7].Node, 3, ffvsp[-5].Node, ffvsp[-3].Node, ffvsp[-1].Node,  0,   0 ); TEST(ffval.Node); ;
    break;}
case 61:
#line 760 "eval.y"
{ ffval.Node = New_Deref( ffvsp[-9].Node, 4, ffvsp[-7].Node, ffvsp[-5].Node, ffvsp[-3].Node, ffvsp[-1].Node,   0 ); TEST(ffval.Node); ;
    break;}
case 62:
#line 762 "eval.y"
{ ffval.Node = New_Deref( ffvsp[-11].Node, 5, ffvsp[-9].Node, ffvsp[-7].Node, ffvsp[-5].Node, ffvsp[-3].Node, ffvsp[-1].Node ); TEST(ffval.Node); ;
    break;}
case 63:
#line 764 "eval.y"
{ ffval.Node = New_Unary( LONG,   INTCAST, ffvsp[0].Node );  TEST(ffval.Node);  ;
    break;}
case 64:
#line 766 "eval.y"
{ ffval.Node = New_Unary( LONG,   INTCAST, ffvsp[0].Node );  TEST(ffval.Node);  ;
    break;}
case 65:
#line 768 "eval.y"
{ ffval.Node = New_Unary( DOUBLE, FLTCAST, ffvsp[0].Node );  TEST(ffval.Node);  ;
    break;}
case 66:
#line 770 "eval.y"
{ ffval.Node = New_Unary( DOUBLE, FLTCAST, ffvsp[0].Node );  TEST(ffval.Node);  ;
    break;}
case 67:
#line 774 "eval.y"
{ ffval.Node = New_Const( BOOLEAN, &(ffvsp[0].log), sizeof(char) ); TEST(ffval.Node); ;
    break;}
case 68:
#line 776 "eval.y"
{ ffval.Node = New_Column( ffvsp[0].lng ); TEST(ffval.Node); ;
    break;}
case 69:
#line 778 "eval.y"
{
                  if( TYPE(ffvsp[-1].Node) != LONG
		      || OPER(ffvsp[-1].Node) != CONST_OP ) {
		     fferror("Offset argument must be a constant integer");
		     FFERROR;
		  }
                  ffval.Node = New_Offset( ffvsp[-3].lng, ffvsp[-1].Node ); TEST(ffval.Node);
                ;
    break;}
case 70:
#line 787 "eval.y"
{ ffval.Node = New_BinOp( BOOLEAN, ffvsp[-2].Node, EQ,  ffvsp[0].Node ); TEST(ffval.Node);
		  SIZE(ffval.Node) = 1;                                     ;
    break;}
case 71:
#line 790 "eval.y"
{ ffval.Node = New_BinOp( BOOLEAN, ffvsp[-2].Node, NE,  ffvsp[0].Node ); TEST(ffval.Node); 
		  SIZE(ffval.Node) = 1;                                     ;
    break;}
case 72:
#line 793 "eval.y"
{ ffval.Node = New_BinOp( BOOLEAN, ffvsp[-2].Node, LT,  ffvsp[0].Node ); TEST(ffval.Node); 
		  SIZE(ffval.Node) = 1;                                     ;
    break;}
case 73:
#line 796 "eval.y"
{ ffval.Node = New_BinOp( BOOLEAN, ffvsp[-2].Node, LTE, ffvsp[0].Node ); TEST(ffval.Node); 
		  SIZE(ffval.Node) = 1;                                     ;
    break;}
case 74:
#line 799 "eval.y"
{ ffval.Node = New_BinOp( BOOLEAN, ffvsp[-2].Node, GT,  ffvsp[0].Node ); TEST(ffval.Node); 
		  SIZE(ffval.Node) = 1;                                     ;
    break;}
case 75:
#line 802 "eval.y"
{ ffval.Node = New_BinOp( BOOLEAN, ffvsp[-2].Node, GTE, ffvsp[0].Node ); TEST(ffval.Node); 
		  SIZE(ffval.Node) = 1;                                     ;
    break;}
case 76:
#line 805 "eval.y"
{ PROMOTE(ffvsp[-2].Node,ffvsp[0].Node); ffval.Node = New_BinOp( BOOLEAN, ffvsp[-2].Node, GT,  ffvsp[0].Node );
                  TEST(ffval.Node);                                               ;
    break;}
case 77:
#line 808 "eval.y"
{ PROMOTE(ffvsp[-2].Node,ffvsp[0].Node); ffval.Node = New_BinOp( BOOLEAN, ffvsp[-2].Node, LT,  ffvsp[0].Node );
                  TEST(ffval.Node);                                               ;
    break;}
case 78:
#line 811 "eval.y"
{ PROMOTE(ffvsp[-2].Node,ffvsp[0].Node); ffval.Node = New_BinOp( BOOLEAN, ffvsp[-2].Node, GTE, ffvsp[0].Node );
                  TEST(ffval.Node);                                               ;
    break;}
case 79:
#line 814 "eval.y"
{ PROMOTE(ffvsp[-2].Node,ffvsp[0].Node); ffval.Node = New_BinOp( BOOLEAN, ffvsp[-2].Node, LTE, ffvsp[0].Node );
                  TEST(ffval.Node);                                               ;
    break;}
case 80:
#line 817 "eval.y"
{ PROMOTE(ffvsp[-2].Node,ffvsp[0].Node); ffval.Node = New_BinOp( BOOLEAN, ffvsp[-2].Node, '~', ffvsp[0].Node );
                  TEST(ffval.Node);                                               ;
    break;}
case 81:
#line 820 "eval.y"
{ PROMOTE(ffvsp[-2].Node,ffvsp[0].Node); ffval.Node = New_BinOp( BOOLEAN, ffvsp[-2].Node, EQ,  ffvsp[0].Node );
                  TEST(ffval.Node);                                               ;
    break;}
case 82:
#line 823 "eval.y"
{ PROMOTE(ffvsp[-2].Node,ffvsp[0].Node); ffval.Node = New_BinOp( BOOLEAN, ffvsp[-2].Node, NE,  ffvsp[0].Node );
                  TEST(ffval.Node);                                               ;
    break;}
case 83:
#line 826 "eval.y"
{ ffval.Node = New_BinOp( BOOLEAN, ffvsp[-2].Node, EQ,  ffvsp[0].Node ); TEST(ffval.Node);
                  SIZE(ffval.Node) = 1; ;
    break;}
case 84:
#line 829 "eval.y"
{ ffval.Node = New_BinOp( BOOLEAN, ffvsp[-2].Node, NE,  ffvsp[0].Node ); TEST(ffval.Node);
                  SIZE(ffval.Node) = 1; ;
    break;}
case 85:
#line 832 "eval.y"
{ ffval.Node = New_BinOp( BOOLEAN, ffvsp[-2].Node, GT,  ffvsp[0].Node ); TEST(ffval.Node);
                  SIZE(ffval.Node) = 1; ;
    break;}
case 86:
#line 835 "eval.y"
{ ffval.Node = New_BinOp( BOOLEAN, ffvsp[-2].Node, GTE, ffvsp[0].Node ); TEST(ffval.Node);
                  SIZE(ffval.Node) = 1; ;
    break;}
case 87:
#line 838 "eval.y"
{ ffval.Node = New_BinOp( BOOLEAN, ffvsp[-2].Node, LT,  ffvsp[0].Node ); TEST(ffval.Node);
                  SIZE(ffval.Node) = 1; ;
    break;}
case 88:
#line 841 "eval.y"
{ ffval.Node = New_BinOp( BOOLEAN, ffvsp[-2].Node, LTE, ffvsp[0].Node ); TEST(ffval.Node);
                  SIZE(ffval.Node) = 1; ;
    break;}
case 89:
#line 844 "eval.y"
{ ffval.Node = New_BinOp( BOOLEAN, ffvsp[-2].Node, AND, ffvsp[0].Node ); TEST(ffval.Node); ;
    break;}
case 90:
#line 846 "eval.y"
{ ffval.Node = New_BinOp( BOOLEAN, ffvsp[-2].Node, OR,  ffvsp[0].Node ); TEST(ffval.Node); ;
    break;}
case 91:
#line 848 "eval.y"
{ ffval.Node = New_BinOp( BOOLEAN, ffvsp[-2].Node, EQ,  ffvsp[0].Node ); TEST(ffval.Node); ;
    break;}
case 92:
#line 850 "eval.y"
{ ffval.Node = New_BinOp( BOOLEAN, ffvsp[-2].Node, NE,  ffvsp[0].Node ); TEST(ffval.Node); ;
    break;}
case 93:
#line 853 "eval.y"
{ PROMOTE(ffvsp[-4].Node,ffvsp[-2].Node); PROMOTE(ffvsp[-4].Node,ffvsp[0].Node); PROMOTE(ffvsp[-2].Node,ffvsp[0].Node);
		  ffvsp[-2].Node = New_BinOp( BOOLEAN, ffvsp[-2].Node, LTE, ffvsp[-4].Node );
                  ffvsp[0].Node = New_BinOp( BOOLEAN, ffvsp[-4].Node, LTE, ffvsp[0].Node );
                  ffval.Node = New_BinOp( BOOLEAN, ffvsp[-2].Node, AND, ffvsp[0].Node );
                  TEST(ffval.Node);                                         ;
    break;}
case 94:
#line 860 "eval.y"
{
                  if( ! Test_Dims(ffvsp[-2].Node,ffvsp[0].Node) ) {
                     fferror("Incompatible dimensions in '?:' arguments");
		     FFERROR;
                  }
                  ffval.Node = New_Func( 0, ifthenelse_fct, 3, ffvsp[-2].Node, ffvsp[0].Node, ffvsp[-4].Node,
                                 0, 0, 0, 0 );
                  TEST(ffval.Node);
                  if( SIZE(ffvsp[-2].Node)<SIZE(ffvsp[0].Node) )  Copy_Dims(ffval.Node, ffvsp[0].Node);
                  if( ! Test_Dims(ffvsp[-4].Node,ffval.Node) ) {
                     fferror("Incompatible dimensions in '?:' condition");
		     FFERROR;
                  }
                  if( SIZE(ffval.Node)<SIZE(ffvsp[-4].Node) )  Copy_Dims(ffval.Node, ffvsp[-4].Node);
                ;
    break;}
case 95:
#line 877 "eval.y"
{
		   if (FSTRCMP(ffvsp[-2].str,"ISNULL(") == 0) {
		      ffval.Node = New_Func( 0, isnull_fct, 1, ffvsp[-1].Node, 0, 0,
				     0, 0, 0, 0 );
		      TEST(ffval.Node); 
                      /* Use expression's size, but return BOOLEAN */
		      TYPE(ffval.Node) = BOOLEAN;
		   } else {
		      fferror("Boolean Function(expr) not supported");
		      FFERROR;
		   }
		;
    break;}
case 96:
#line 890 "eval.y"
{
		   if (FSTRCMP(ffvsp[-2].str,"ISNULL(") == 0) {
		      ffval.Node = New_Func( 0, isnull_fct, 1, ffvsp[-1].Node, 0, 0,
				     0, 0, 0, 0 );
		      TEST(ffval.Node); 
                      /* Use expression's size, but return BOOLEAN */
		      TYPE(ffval.Node) = BOOLEAN;
		   } else {
		      fferror("Boolean Function(expr) not supported");
		      FFERROR;
		   }
		;
    break;}
case 97:
#line 903 "eval.y"
{
		   if (FSTRCMP(ffvsp[-2].str,"ISNULL(") == 0) {
		      ffval.Node = New_Func( BOOLEAN, isnull_fct, 1, ffvsp[-1].Node, 0, 0,
				     0, 0, 0, 0 );
		      TEST(ffval.Node); 
		   } else {
		      fferror("Boolean Function(expr) not supported");
		      FFERROR;
		   }
		;
    break;}
case 98:
#line 914 "eval.y"
{
		   if (FSTRCMP(ffvsp[-4].str,"DEFNULL(") == 0) {
		      if( SIZE(ffvsp[-3].Node)>=SIZE(ffvsp[-1].Node) && Test_Dims( ffvsp[-3].Node, ffvsp[-1].Node ) ) {
			 ffval.Node = New_Func( 0, defnull_fct, 2, ffvsp[-3].Node, ffvsp[-1].Node, 0,
					0, 0, 0, 0 );
			 TEST(ffval.Node); 
		      } else {
			 fferror("Dimensions of DEFNULL arguments are not compatible");
			 FFERROR;
		      }
		   } else {
		      fferror("Boolean Function(expr,expr) not supported");
		      FFERROR;
		   }
		;
    break;}
case 99:
#line 930 "eval.y"
{
		   if( TYPE(ffvsp[-5].Node) != DOUBLE ) ffvsp[-5].Node = New_Unary( DOUBLE, 0, ffvsp[-5].Node );
		   if( TYPE(ffvsp[-3].Node) != DOUBLE ) ffvsp[-3].Node = New_Unary( DOUBLE, 0, ffvsp[-3].Node );
		   if( TYPE(ffvsp[-1].Node) != DOUBLE ) ffvsp[-1].Node = New_Unary( DOUBLE, 0, ffvsp[-1].Node );
		   if( ! (Test_Dims( ffvsp[-5].Node, ffvsp[-3].Node ) && Test_Dims( ffvsp[-3].Node, ffvsp[-1].Node ) ) ) {
		       fferror("Dimensions of NEAR arguments "
			       "are not compatible");
		       FFERROR;
		   } else {
		     if (FSTRCMP(ffvsp[-6].str,"NEAR(") == 0) {
		       ffval.Node = New_Func( BOOLEAN, near_fct, 3, ffvsp[-5].Node, ffvsp[-3].Node, ffvsp[-1].Node,
				      0, 0, 0, 0 );
		     } else {
		       fferror("Boolean Function not supported");
		       FFERROR;
		     }
		     TEST(ffval.Node); 

		     if( SIZE(ffval.Node)<SIZE(ffvsp[-5].Node) )  Copy_Dims(ffval.Node, ffvsp[-5].Node);
		     if( SIZE(ffvsp[-5].Node)<SIZE(ffvsp[-3].Node) )  Copy_Dims(ffval.Node, ffvsp[-3].Node);
		     if( SIZE(ffvsp[-3].Node)<SIZE(ffvsp[-1].Node) )  Copy_Dims(ffval.Node, ffvsp[-1].Node);
		   }
		;
    break;}
case 100:
#line 954 "eval.y"
{
		   if( TYPE(ffvsp[-9].Node) != DOUBLE ) ffvsp[-9].Node = New_Unary( DOUBLE, 0, ffvsp[-9].Node );
		   if( TYPE(ffvsp[-7].Node) != DOUBLE ) ffvsp[-7].Node = New_Unary( DOUBLE, 0, ffvsp[-7].Node );
		   if( TYPE(ffvsp[-5].Node) != DOUBLE ) ffvsp[-5].Node = New_Unary( DOUBLE, 0, ffvsp[-5].Node );
		   if( TYPE(ffvsp[-3].Node) != DOUBLE ) ffvsp[-3].Node = New_Unary( DOUBLE, 0, ffvsp[-3].Node );
		   if( TYPE(ffvsp[-1].Node)!= DOUBLE ) ffvsp[-1].Node= New_Unary( DOUBLE, 0, ffvsp[-1].Node);
		   if( ! (Test_Dims( ffvsp[-9].Node, ffvsp[-7].Node ) && Test_Dims( ffvsp[-7].Node, ffvsp[-5].Node ) && 
			  Test_Dims( ffvsp[-5].Node, ffvsp[-3].Node ) && Test_Dims( ffvsp[-3].Node, ffvsp[-1].Node )) ) {
		     fferror("Dimensions of CIRCLE arguments "
			     "are not compatible");
		     FFERROR;
		   } else {
		     if (FSTRCMP(ffvsp[-10].str,"CIRCLE(") == 0) {
		       ffval.Node = New_Func( BOOLEAN, circle_fct, 5, ffvsp[-9].Node, ffvsp[-7].Node, ffvsp[-5].Node, ffvsp[-3].Node,
				      ffvsp[-1].Node, 0, 0 );
		     } else {
		       fferror("Boolean Function not supported");
		       FFERROR;
		     }
		     TEST(ffval.Node); 
		     if( SIZE(ffval.Node)<SIZE(ffvsp[-9].Node) )  Copy_Dims(ffval.Node, ffvsp[-9].Node);
		     if( SIZE(ffvsp[-9].Node)<SIZE(ffvsp[-7].Node) )  Copy_Dims(ffval.Node, ffvsp[-7].Node);
		     if( SIZE(ffvsp[-7].Node)<SIZE(ffvsp[-5].Node) )  Copy_Dims(ffval.Node, ffvsp[-5].Node);
		     if( SIZE(ffvsp[-5].Node)<SIZE(ffvsp[-3].Node) )  Copy_Dims(ffval.Node, ffvsp[-3].Node);
		     if( SIZE(ffvsp[-3].Node)<SIZE(ffvsp[-1].Node) ) Copy_Dims(ffval.Node, ffvsp[-1].Node);
		   }
		;
    break;}
case 101:
#line 982 "eval.y"
{
		   if( TYPE(ffvsp[-13].Node) != DOUBLE ) ffvsp[-13].Node = New_Unary( DOUBLE, 0, ffvsp[-13].Node );
		   if( TYPE(ffvsp[-11].Node) != DOUBLE ) ffvsp[-11].Node = New_Unary( DOUBLE, 0, ffvsp[-11].Node );
		   if( TYPE(ffvsp[-9].Node) != DOUBLE ) ffvsp[-9].Node = New_Unary( DOUBLE, 0, ffvsp[-9].Node );
		   if( TYPE(ffvsp[-7].Node) != DOUBLE ) ffvsp[-7].Node = New_Unary( DOUBLE, 0, ffvsp[-7].Node );
		   if( TYPE(ffvsp[-5].Node)!= DOUBLE ) ffvsp[-5].Node= New_Unary( DOUBLE, 0, ffvsp[-5].Node);
		   if( TYPE(ffvsp[-3].Node)!= DOUBLE ) ffvsp[-3].Node= New_Unary( DOUBLE, 0, ffvsp[-3].Node);
		   if( TYPE(ffvsp[-1].Node)!= DOUBLE ) ffvsp[-1].Node= New_Unary( DOUBLE, 0, ffvsp[-1].Node);
		   if( ! (Test_Dims( ffvsp[-13].Node, ffvsp[-11].Node ) && Test_Dims( ffvsp[-11].Node, ffvsp[-9].Node ) && 
			  Test_Dims( ffvsp[-9].Node, ffvsp[-7].Node ) && Test_Dims( ffvsp[-7].Node, ffvsp[-5].Node ) &&
			  Test_Dims(ffvsp[-5].Node,ffvsp[-3].Node ) && Test_Dims(ffvsp[-3].Node, ffvsp[-1].Node ) ) ) {
		     fferror("Dimensions of BOX or ELLIPSE arguments "
			     "are not compatible");
		     FFERROR;
		   } else {
		     if (FSTRCMP(ffvsp[-14].str,"BOX(") == 0) {
		       ffval.Node = New_Func( BOOLEAN, box_fct, 7, ffvsp[-13].Node, ffvsp[-11].Node, ffvsp[-9].Node, ffvsp[-7].Node,
				      ffvsp[-5].Node, ffvsp[-3].Node, ffvsp[-1].Node );
		     } else if (FSTRCMP(ffvsp[-14].str,"ELLIPSE(") == 0) {
		       ffval.Node = New_Func( BOOLEAN, elps_fct, 7, ffvsp[-13].Node, ffvsp[-11].Node, ffvsp[-9].Node, ffvsp[-7].Node,
				      ffvsp[-5].Node, ffvsp[-3].Node, ffvsp[-1].Node );
		     } else {
		       fferror("SAO Image Function not supported");
		       FFERROR;
		     }
		     TEST(ffval.Node); 
		     if( SIZE(ffval.Node)<SIZE(ffvsp[-13].Node) )  Copy_Dims(ffval.Node, ffvsp[-13].Node);
		     if( SIZE(ffvsp[-13].Node)<SIZE(ffvsp[-11].Node) )  Copy_Dims(ffval.Node, ffvsp[-11].Node);
		     if( SIZE(ffvsp[-11].Node)<SIZE(ffvsp[-9].Node) )  Copy_Dims(ffval.Node, ffvsp[-9].Node);
		     if( SIZE(ffvsp[-9].Node)<SIZE(ffvsp[-7].Node) )  Copy_Dims(ffval.Node, ffvsp[-7].Node);
		     if( SIZE(ffvsp[-7].Node)<SIZE(ffvsp[-5].Node) ) Copy_Dims(ffval.Node, ffvsp[-5].Node);
		     if( SIZE(ffvsp[-5].Node)<SIZE(ffvsp[-3].Node) ) Copy_Dims(ffval.Node, ffvsp[-3].Node);
		     if( SIZE(ffvsp[-3].Node)<SIZE(ffvsp[-1].Node) ) Copy_Dims(ffval.Node, ffvsp[-1].Node);
		   }
		;
    break;}
case 102:
#line 1019 "eval.y"
{ /* Use defaults for all elements */
                   ffval.Node = New_GTI( "", -99, "*START*", "*STOP*" );
                   TEST(ffval.Node);                                        ;
    break;}
case 103:
#line 1023 "eval.y"
{ /* Use defaults for all except filename */
                   ffval.Node = New_GTI( ffvsp[-1].str, -99, "*START*", "*STOP*" );
                   TEST(ffval.Node);                                        ;
    break;}
case 104:
#line 1027 "eval.y"
{  ffval.Node = New_GTI( ffvsp[-3].str, ffvsp[-1].Node, "*START*", "*STOP*" );
                   TEST(ffval.Node);                                        ;
    break;}
case 105:
#line 1030 "eval.y"
{  ffval.Node = New_GTI( ffvsp[-7].str, ffvsp[-5].Node, ffvsp[-3].str, ffvsp[-1].str );
                   TEST(ffval.Node);                                        ;
    break;}
case 106:
#line 1034 "eval.y"
{ /* Use defaults for all except filename */
                   ffval.Node = New_REG( ffvsp[-1].str, -99, -99, "" );
                   TEST(ffval.Node);                                        ;
    break;}
case 107:
#line 1038 "eval.y"
{  ffval.Node = New_REG( ffvsp[-5].str, ffvsp[-3].Node, ffvsp[-1].Node, "" );
                   TEST(ffval.Node);                                        ;
    break;}
case 108:
#line 1041 "eval.y"
{  ffval.Node = New_REG( ffvsp[-7].str, ffvsp[-5].Node, ffvsp[-3].Node, ffvsp[-1].str );
                   TEST(ffval.Node);                                        ;
    break;}
case 109:
#line 1045 "eval.y"
{ ffval.Node = New_Deref( ffvsp[-3].Node, 1, ffvsp[-1].Node,  0,  0,  0,   0 ); TEST(ffval.Node); ;
    break;}
case 110:
#line 1047 "eval.y"
{ ffval.Node = New_Deref( ffvsp[-5].Node, 2, ffvsp[-3].Node, ffvsp[-1].Node,  0,  0,   0 ); TEST(ffval.Node); ;
    break;}
case 111:
#line 1049 "eval.y"
{ ffval.Node = New_Deref( ffvsp[-7].Node, 3, ffvsp[-5].Node, ffvsp[-3].Node, ffvsp[-1].Node,  0,   0 ); TEST(ffval.Node); ;
    break;}
case 112:
#line 1051 "eval.y"
{ ffval.Node = New_Deref( ffvsp[-9].Node, 4, ffvsp[-7].Node, ffvsp[-5].Node, ffvsp[-3].Node, ffvsp[-1].Node,   0 ); TEST(ffval.Node); ;
    break;}
case 113:
#line 1053 "eval.y"
{ ffval.Node = New_Deref( ffvsp[-11].Node, 5, ffvsp[-9].Node, ffvsp[-7].Node, ffvsp[-5].Node, ffvsp[-3].Node, ffvsp[-1].Node ); TEST(ffval.Node); ;
    break;}
case 114:
#line 1055 "eval.y"
{ ffval.Node = New_Unary( BOOLEAN, NOT, ffvsp[0].Node ); TEST(ffval.Node); ;
    break;}
case 115:
#line 1057 "eval.y"
{ ffval.Node = ffvsp[-1].Node; ;
    break;}
case 116:
#line 1061 "eval.y"
{ ffval.Node = New_Const( STRING, ffvsp[0].str, strlen(ffvsp[0].str)+1 ); TEST(ffval.Node);
                  SIZE(ffval.Node) = strlen(ffvsp[0].str); ;
    break;}
case 117:
#line 1064 "eval.y"
{ ffval.Node = New_Column( ffvsp[0].lng ); TEST(ffval.Node); ;
    break;}
case 118:
#line 1066 "eval.y"
{
                  if( TYPE(ffvsp[-1].Node) != LONG
		      || OPER(ffvsp[-1].Node) != CONST_OP ) {
		     fferror("Offset argument must be a constant integer");
		     FFERROR;
		  }
                  ffval.Node = New_Offset( ffvsp[-3].lng, ffvsp[-1].Node ); TEST(ffval.Node);
                ;
    break;}
case 119:
#line 1075 "eval.y"
{ ffval.Node = New_Func( STRING, null_fct, 0, 0, 0, 0, 0, 0, 0, 0 ); ;
    break;}
case 120:
#line 1077 "eval.y"
{ ffval.Node = ffvsp[-1].Node; ;
    break;}
case 121:
#line 1079 "eval.y"
{ 
		  if (SIZE(ffvsp[-2].Node)+SIZE(ffvsp[0].Node) >= MAX_STRLEN) {
		    fferror("Combined string size exceeds " MAX_STRLEN_S " characters");
		    FFERROR;
		  }
		  ffval.Node = New_BinOp( STRING, ffvsp[-2].Node, '+', ffvsp[0].Node );  TEST(ffval.Node);
		  SIZE(ffval.Node) = SIZE(ffvsp[-2].Node) + SIZE(ffvsp[0].Node);
		;
    break;}
case 122:
#line 1088 "eval.y"
{
		  int outSize;
                  if( SIZE(ffvsp[-4].Node)!=1 ) {
                     fferror("Cannot have a vector string column");
		     FFERROR;
                  }
		  /* Since the output can be calculated now, as a constant
		     scalar, we must precalculate the output size, in
		     order to avoid an overflow. */
		  outSize = SIZE(ffvsp[-2].Node);
		  if (SIZE(ffvsp[0].Node) > outSize) outSize = SIZE(ffvsp[0].Node);
                  ffval.Node = New_FuncSize( 0, ifthenelse_fct, 3, ffvsp[-2].Node, ffvsp[0].Node, ffvsp[-4].Node,
				     0, 0, 0, 0, outSize);
		  
                  TEST(ffval.Node);
                  if( SIZE(ffvsp[-2].Node)<SIZE(ffvsp[0].Node) )  Copy_Dims(ffval.Node, ffvsp[0].Node);
                ;
    break;}
case 123:
#line 1107 "eval.y"
{ 
		  if (FSTRCMP(ffvsp[-4].str,"DEFNULL(") == 0) {
		     int outSize;
		     /* Since the output can be calculated now, as a constant
			scalar, we must precalculate the output size, in
			order to avoid an overflow. */
		     outSize = SIZE(ffvsp[-3].Node);
		     if (SIZE(ffvsp[-1].Node) > outSize) outSize = SIZE(ffvsp[-1].Node);
		     
		     ffval.Node = New_FuncSize( 0, defnull_fct, 2, ffvsp[-3].Node, ffvsp[-1].Node, 0,
					0, 0, 0, 0, outSize );
		     TEST(ffval.Node); 
		     if( SIZE(ffvsp[-1].Node)>SIZE(ffvsp[-3].Node) ) SIZE(ffval.Node) = SIZE(ffvsp[-1].Node);
		  } else {
		     fferror("Function(string,string) not supported");
		     FFERROR;
		  }
		;
    break;}
case 124:
#line 1126 "eval.y"
{ 
		  if (FSTRCMP(ffvsp[-6].str,"STRMID(") == 0) {
		    int len;
		    if( TYPE(ffvsp[-3].Node) != LONG || SIZE(ffvsp[-3].Node) != 1 ||
			TYPE(ffvsp[-1].Node) != LONG || SIZE(ffvsp[-1].Node) != 1) {
		      fferror("When using STRMID(S,P,N), P and N must be integers (and not vector columns)");
		      FFERROR;
		    }
		    if (OPER(ffvsp[-1].Node) == CONST_OP) {
		      /* Constant value: use that directly */
		      len = (gParse.Nodes[ffvsp[-1].Node].value.data.lng);
		    } else {
		      /* Variable value: use the maximum possible (from $2) */
		      len = SIZE(ffvsp[-5].Node);
		    }
		    if (len <= 0 || len >= MAX_STRLEN) {
		      fferror("STRMID(S,P,N), N must be 1-" MAX_STRLEN_S);
		      FFERROR;
		    }
		    ffval.Node = New_FuncSize( 0, strmid_fct, 3, ffvsp[-5].Node, ffvsp[-3].Node,ffvsp[-1].Node,0,0,0,0,len);
		    TEST(ffval.Node);
		  } else {
		     fferror("Function(string,expr,expr) not supported");
		     FFERROR;
		  }
		;
    break;}
}
   /* the action file gets copied in in place of this dollarsign */
#line 498 "/usr1/local/share/bison.simple"

  ffvsp -= fflen;
  ffssp -= fflen;
#ifdef FFLSP_NEEDED
  fflsp -= fflen;
#endif

#if FFDEBUG != 0
  if (ffdebug)
    {
      short *ssp1 = ffss - 1;
      fprintf (stderr, "state stack now");
      while (ssp1 != ffssp)
	fprintf (stderr, " %d", *++ssp1);
      fprintf (stderr, "\n");
    }
#endif

  *++ffvsp = ffval;

#ifdef FFLSP_NEEDED
  fflsp++;
  if (fflen == 0)
    {
      fflsp->first_line = fflloc.first_line;
      fflsp->first_column = fflloc.first_column;
      fflsp->last_line = (fflsp-1)->last_line;
      fflsp->last_column = (fflsp-1)->last_column;
      fflsp->text = 0;
    }
  else
    {
      fflsp->last_line = (fflsp+fflen-1)->last_line;
      fflsp->last_column = (fflsp+fflen-1)->last_column;
    }
#endif

  /* Now "shift" the result of the reduction.
     Determine what state that goes to,
     based on the state we popped back to
     and the rule number reduced by.  */

  ffn = ffr1[ffn];

  ffstate = ffpgoto[ffn - FFNTBASE] + *ffssp;
  if (ffstate >= 0 && ffstate <= FFLAST && ffcheck[ffstate] == *ffssp)
    ffstate = fftable[ffstate];
  else
    ffstate = ffdefgoto[ffn - FFNTBASE];

  goto ffnewstate;

fferrlab:   /* here on detecting error */

  if (! fferrstatus)
    /* If not already recovering from an error, report this error.  */
    {
      ++ffnerrs;

#ifdef FFERROR_VERBOSE
      ffn = ffpact[ffstate];

      if (ffn > FFFLAG && ffn < FFLAST)
	{
	  int size = 0;
	  char *msg;
	  int x, count;

	  count = 0;
	  /* Start X at -ffn if nec to avoid negative indexes in ffcheck.  */
	  for (x = (ffn < 0 ? -ffn : 0);
	       x < (sizeof(fftname) / sizeof(char *)); x++)
	    if (ffcheck[x + ffn] == x)
	      size += strlen(fftname[x]) + 15, count++;
	  msg = (char *) malloc(size + 15);
	  if (msg != 0)
	    {
	      strcpy(msg, "parse error");

	      if (count < 5)
		{
		  count = 0;
		  for (x = (ffn < 0 ? -ffn : 0);
		       x < (sizeof(fftname) / sizeof(char *)); x++)
		    if (ffcheck[x + ffn] == x)
		      {
			strcat(msg, count == 0 ? ", expecting `" : " or `");
			strcat(msg, fftname[x]);
			strcat(msg, "'");
			count++;
		      }
		}
	      fferror(msg);
	      free(msg);
	    }
	  else
	    fferror ("parse error; also virtual memory exceeded");
	}
      else
#endif /* FFERROR_VERBOSE */
	fferror("parse error");
    }

  goto fferrlab1;
fferrlab1:   /* here on error raised explicitly by an action */

  if (fferrstatus == 3)
    {
      /* if just tried and failed to reuse lookahead token after an error, discard it.  */

      /* return failure if at end of input */
      if (ffchar == FFEOF)
	FFABORT;

#if FFDEBUG != 0
      if (ffdebug)
	fprintf(stderr, "Discarding token %d (%s).\n", ffchar, fftname[ffchar1]);
#endif

      ffchar = FFEMPTY;
    }

  /* Else will try to reuse lookahead token
     after shifting the error token.  */

  fferrstatus = 3;		/* Each real token shifted decrements this */

  goto fferrhandle;

fferrdefault:  /* current state does not do anything special for the error token. */

#if 0
  /* This is wrong; only states that explicitly want error tokens
     should shift them.  */
  ffn = ffdefact[ffstate];  /* If its default is to accept any token, ok.  Otherwise pop it.*/
  if (ffn) goto ffdefault;
#endif

fferrpop:   /* pop the current state because it cannot handle the error token */

  if (ffssp == ffss) FFABORT;
  ffvsp--;
  ffstate = *--ffssp;
#ifdef FFLSP_NEEDED
  fflsp--;
#endif

#if FFDEBUG != 0
  if (ffdebug)
    {
      short *ssp1 = ffss - 1;
      fprintf (stderr, "Error: state stack now");
      while (ssp1 != ffssp)
	fprintf (stderr, " %d", *++ssp1);
      fprintf (stderr, "\n");
    }
#endif

fferrhandle:

  ffn = ffpact[ffstate];
  if (ffn == FFFLAG)
    goto fferrdefault;

  ffn += FFTERROR;
  if (ffn < 0 || ffn > FFLAST || ffcheck[ffn] != FFTERROR)
    goto fferrdefault;

  ffn = fftable[ffn];
  if (ffn < 0)
    {
      if (ffn == FFFLAG)
	goto fferrpop;
      ffn = -ffn;
      goto ffreduce;
    }
  else if (ffn == 0)
    goto fferrpop;

  if (ffn == FFFINAL)
    FFACCEPT;

#if FFDEBUG != 0
  if (ffdebug)
    fprintf(stderr, "Shifting error token, ");
#endif

  *++ffvsp = fflval;
#ifdef FFLSP_NEEDED
  *++fflsp = fflloc;
#endif

  ffstate = ffn;
  goto ffnewstate;
}
#line 1155 "eval.y"


/*************************************************************************/
/*  Start of "New" routines which build the expression Nodal structure   */
/*************************************************************************/

static int Alloc_Node( void )
{
                      /* Use this for allocation to guarantee *Nodes */
   Node *newNodePtr;  /* survives on failure, making it still valid  */
                      /* while working our way out of this error     */

   if( gParse.nNodes == gParse.nNodesAlloc ) {
      if( gParse.Nodes ) {
	 gParse.nNodesAlloc += gParse.nNodesAlloc;
	 newNodePtr = (Node *)realloc( gParse.Nodes,
				       sizeof(Node)*gParse.nNodesAlloc );
      } else {
	 gParse.nNodesAlloc = 100;
	 newNodePtr = (Node *)malloc ( sizeof(Node)*gParse.nNodesAlloc );
      }	 

      if( newNodePtr ) {
	 gParse.Nodes = newNodePtr;
      } else {
	 gParse.status = MEMORY_ALLOCATION;
	 return( -1 );
      }
   }

   return ( gParse.nNodes++ );
}

static void Free_Last_Node( void )
{
   if( gParse.nNodes ) gParse.nNodes--;
}

static int New_Const( int returnType, void *value, long len )
{
   Node *this;
   int n;

   n = Alloc_Node();
   if( n>=0 ) {
      this             = gParse.Nodes + n;
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

static int New_Column( int ColNum )
{
   Node *this;
   int  n, i;

   n = Alloc_Node();
   if( n>=0 ) {
      this              = gParse.Nodes + n;
      this->operation   = -ColNum;
      this->DoOp        = NULL;
      this->nSubNodes   = 0;
      this->type        = gParse.varData[ColNum].type;
      this->value.nelem = gParse.varData[ColNum].nelem;
      this->value.naxis = gParse.varData[ColNum].naxis;
      for( i=0; i<gParse.varData[ColNum].naxis; i++ )
	 this->value.naxes[i] = gParse.varData[ColNum].naxes[i];
   }
   return(n);
}

static int New_Offset( int ColNum, int offsetNode )
{
   Node *this;
   int  n, i, colNode;

   colNode = New_Column( ColNum );
   if( colNode<0 ) return(-1);

   n = Alloc_Node();
   if( n>=0 ) {
      this              = gParse.Nodes + n;
      this->operation   = '{';
      this->DoOp        = Do_Offset;
      this->nSubNodes   = 2;
      this->SubNodes[0] = colNode;
      this->SubNodes[1] = offsetNode;
      this->type        = gParse.varData[ColNum].type;
      this->value.nelem = gParse.varData[ColNum].nelem;
      this->value.naxis = gParse.varData[ColNum].naxis;
      for( i=0; i<gParse.varData[ColNum].naxis; i++ )
	 this->value.naxes[i] = gParse.varData[ColNum].naxes[i];
   }
   return(n);
}

static int New_Unary( int returnType, int Op, int Node1 )
{
   Node *this, *that;
   int  i,n;

   if( Node1<0 ) return(-1);
   that = gParse.Nodes + Node1;

   if( !Op ) Op = returnType;

   if( (Op==DOUBLE || Op==FLTCAST) && that->type==DOUBLE  ) return( Node1 );
   if( (Op==LONG   || Op==INTCAST) && that->type==LONG    ) return( Node1 );
   if( (Op==BOOLEAN              ) && that->type==BOOLEAN ) return( Node1 );
   
   n = Alloc_Node();
   if( n>=0 ) {
      this              = gParse.Nodes + n;
      this->operation   = Op;
      this->DoOp        = Do_Unary;
      this->nSubNodes   = 1;
      this->SubNodes[0] = Node1;
      this->type        = returnType;

      that              = gParse.Nodes + Node1; /* Reset in case .Nodes mv'd */
      this->value.nelem = that->value.nelem;
      this->value.naxis = that->value.naxis;
      for( i=0; i<that->value.naxis; i++ )
	 this->value.naxes[i] = that->value.naxes[i];

      if( that->operation==CONST_OP ) this->DoOp( this );
   }
   return( n );
}

static int New_BinOp( int returnType, int Node1, int Op, int Node2 )
{
   Node *this,*that1,*that2;
   int  n,i,constant;

   if( Node1<0 || Node2<0 ) return(-1);

   n = Alloc_Node();
   if( n>=0 ) {
      this             = gParse.Nodes + n;
      this->operation  = Op;
      this->nSubNodes  = 2;
      this->SubNodes[0]= Node1;
      this->SubNodes[1]= Node2;
      this->type       = returnType;

      that1            = gParse.Nodes + Node1;
      that2            = gParse.Nodes + Node2;
      constant         = (that1->operation==CONST_OP
                          && that2->operation==CONST_OP);
      if( that1->type!=STRING && that1->type!=BITSTR )
	 if( !Test_Dims( Node1, Node2 ) ) {
	    Free_Last_Node();
	    fferror("Array sizes/dims do not match for binary operator");
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
      if( constant ) this->DoOp( this );
   }
   return( n );
}

static int New_Func( int returnType, funcOp Op, int nNodes,
		     int Node1, int Node2, int Node3, int Node4, 
		     int Node5, int Node6, int Node7 )
{
  return New_FuncSize(returnType, Op, nNodes,
		      Node1, Node2, Node3, Node4, 
		      Node5, Node6, Node7, 0);
}

static int New_FuncSize( int returnType, funcOp Op, int nNodes,
		     int Node1, int Node2, int Node3, int Node4, 
			 int Node5, int Node6, int Node7, int Size )
/* If returnType==0 , use Node1's type and vector sizes as returnType, */
/* else return a single value of type returnType                       */
{
   Node *this, *that;
   int  i,n,constant;

   if( Node1<0 || Node2<0 || Node3<0 || Node4<0 || 
       Node5<0 || Node6<0 || Node7<0 ) return(-1);

   n = Alloc_Node();
   if( n>=0 ) {
      this              = gParse.Nodes + n;
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
	 that              = gParse.Nodes + Node1;
	 this->type        = that->type;
	 this->value.nelem = that->value.nelem;
	 this->value.naxis = that->value.naxis;
	 for( i=0; i<that->value.naxis; i++ )
	    this->value.naxes[i] = that->value.naxes[i];
      }
      /* Force explicit size before evaluating */
      if (Size > 0) this->value.nelem = Size;

      if( constant ) this->DoOp( this );
   }
   return( n );
}

static int New_Deref( int Var,  int nDim,
		      int Dim1, int Dim2, int Dim3, int Dim4, int Dim5 )
{
   int n, idx, constant;
   long elem=0;
   Node *this, *theVar, *theDim[MAXDIMS];

   if( Var<0 || Dim1<0 || Dim2<0 || Dim3<0 || Dim4<0 || Dim5<0 ) return(-1);

   theVar = gParse.Nodes + Var;
   if( theVar->operation==CONST_OP || theVar->value.nelem==1 ) {
      fferror("Cannot index a scalar value");
      return(-1);
   }

   n = Alloc_Node();
   if( n>=0 ) {
      this              = gParse.Nodes + n;
      this->nSubNodes   = nDim+1;
      theVar            = gParse.Nodes + (this->SubNodes[0]=Var);
      theDim[0]         = gParse.Nodes + (this->SubNodes[1]=Dim1);
      theDim[1]         = gParse.Nodes + (this->SubNodes[2]=Dim2);
      theDim[2]         = gParse.Nodes + (this->SubNodes[3]=Dim3);
      theDim[3]         = gParse.Nodes + (this->SubNodes[4]=Dim4);
      theDim[4]         = gParse.Nodes + (this->SubNodes[5]=Dim5);
      constant          = theVar->operation==CONST_OP;
      for( idx=0; idx<nDim; idx++ )
	 constant = (constant && theDim[idx]->operation==CONST_OP);

      for( idx=0; idx<nDim; idx++ )
	 if( theDim[idx]->value.nelem>1 ) {
	    Free_Last_Node();
	    fferror("Cannot use an array as an index value");
	    return(-1);
	 } else if( theDim[idx]->type!=LONG ) {
	    Free_Last_Node();
	    fferror("Index value must be an integer type");
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
	 Free_Last_Node();
	 fferror("Must specify just one or all indices for vector");
	 return(-1);
      }
      if( constant ) this->DoOp( this );
   }
   return(n);
}

extern int ffGetVariable( char *varName, FFSTYPE *varVal );

static int New_GTI( char *fname, int Node1, char *start, char *stop )
{
   fitsfile *fptr;
   Node *this, *that0, *that1;
   int  type,i,n, startCol, stopCol, Node0;
   int  hdutype, hdunum, evthdu, samefile, extvers, movetotype, tstat;
   char extname[100];
   long nrows;
   double timeZeroI[2], timeZeroF[2], dt, timeSpan;
   char xcol[20], xexpr[20];
   FFSTYPE colVal;

   if( Node1==-99 ) {
      type = ffGetVariable( "TIME", &colVal );
      if( type==COLUMN ) {
	 Node1 = New_Column( (int)colVal.lng );
      } else {
	 fferror("Could not build TIME column for GTIFILTER");
	 return(-1);
      }
   }
   Node1 = New_Unary( DOUBLE, 0, Node1 );
   Node0 = Alloc_Node(); /* This will hold the START/STOP times */
   if( Node1<0 || Node0<0 ) return(-1);

   /*  Record current HDU number in case we need to move within this file  */

   fptr = gParse.def_fptr;
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
		 xcol, xexpr, &gParse.status );
         if( *extname ) {
	    ffmnhd( fptr, movetotype, extname, extvers, &gParse.status );
	    ffghdn( fptr, &hdunum );
	 } else if( hdunum ) {
	    ffmahd( fptr, ++hdunum, &hdutype, &gParse.status );
	 } else if( !gParse.status ) {
	    fferror("Cannot use primary array for GTI filter");
	    return( -1 );
	 }
      } else {
	 fferror("File extension specifier lacks closing ']'");
	 return( -1 );
      }
      break;
   case '+':
      samefile = 1;
      hdunum = atoi( fname ) + 1;
      if( hdunum>1 )
	 ffmahd( fptr, hdunum, &hdutype, &gParse.status );
      else {
	 fferror("Cannot use primary array for GTI filter");
	 return( -1 );
      }
      break;
   default:
      samefile = 0;
      if( ! ffopen( &fptr, fname, READONLY, &gParse.status ) )
	 ffghdn( fptr, &hdunum );
      break;
   }
   if( gParse.status ) return(-1);

   /*  If at primary, search for GTI extension  */

   if( hdunum==1 ) {
      while( 1 ) {
	 hdunum++;
	 if( ffmahd( fptr, hdunum, &hdutype, &gParse.status ) ) break;
	 if( hdutype==IMAGE_HDU ) continue;
	 tstat = 0;
	 if( ffgkys( fptr, "EXTNAME", extname, NULL, &tstat ) ) continue;
	 ffupch( extname );
	 if( strstr( extname, "GTI" ) ) break;
      }
      if( gParse.status ) {
	 if( gParse.status==END_OF_FILE )
	    fferror("GTI extension not found in this file");
	 return(-1);
      }
   }

   /*  Locate START/STOP Columns  */

   ffgcno( fptr, CASEINSEN, start, &startCol, &gParse.status );
   ffgcno( fptr, CASEINSEN, stop,  &stopCol,  &gParse.status );
   if( gParse.status ) return(-1);

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

   n = Alloc_Node();
   if( n >= 0 ) {
      this                 = gParse.Nodes + n;
      this->nSubNodes      = 2;
      this->SubNodes[1]    = Node1;
      this->operation      = (int)gtifilt_fct;
      this->DoOp           = Do_GTI;
      this->type           = BOOLEAN;
      that1                = gParse.Nodes + Node1;
      this->value.nelem    = that1->value.nelem;
      this->value.naxis    = that1->value.naxis;
      for( i=0; i < that1->value.naxis; i++ )
	 this->value.naxes[i] = that1->value.naxes[i];

      /* Init START/STOP node to be treated as a "constant" */

      this->SubNodes[0]    = Node0;
      that0                = gParse.Nodes + Node0;
      that0->operation     = CONST_OP;
      that0->DoOp          = NULL;
      that0->value.data.ptr= NULL;

      /*  Read in START/STOP times  */

      if( ffgkyj( fptr, "NAXIS2", &nrows, NULL, &gParse.status ) )
	 return(-1);
      that0->value.nelem = nrows;
      if( nrows ) {

	 that0->value.data.dblptr = (double*)malloc( 2*nrows*sizeof(double) );
	 if( !that0->value.data.dblptr ) {
	    gParse.status = MEMORY_ALLOCATION;
	    return(-1);
	 }
	 
	 ffgcvd( fptr, startCol, 1L, 1L, nrows, 0.0,
		 that0->value.data.dblptr, &i, &gParse.status );
	 ffgcvd( fptr, stopCol, 1L, 1L, nrows, 0.0,
		 that0->value.data.dblptr+nrows, &i, &gParse.status );
	 if( gParse.status ) {
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
	 
	 /*  Handle TIMEZERO offset, if any  */
	 
	 dt = (timeZeroI[1] - timeZeroI[0]) + (timeZeroF[1] - timeZeroF[0]);
	 timeSpan = that0->value.data.dblptr[nrows+nrows-1]
	    - that0->value.data.dblptr[0];
	 
	 if( fabs( dt / timeSpan ) > 1e-12 ) {
	    for( i=0; i<(nrows+nrows); i++ )
	       that0->value.data.dblptr[i] += dt;
	 }
      }
      if( OPER(Node1)==CONST_OP )
	 this->DoOp( this );
   }

   if( samefile )
      ffmahd( fptr, evthdu, &hdutype, &gParse.status );
   else
      ffclos( fptr, &gParse.status );

   return( n );
}

static int New_REG( char *fname, int NodeX, int NodeY, char *colNames )
{
   Node *this, *that0;
   int  type, n, Node0;
   int  Xcol, Ycol, tstat;
   WCSdata wcs;
   SAORegion *Rgn;
   char *cX, *cY;
   FFSTYPE colVal;

   if( NodeX==-99 ) {
      type = ffGetVariable( "X", &colVal );
      if( type==COLUMN ) {
	 NodeX = New_Column( (int)colVal.lng );
      } else {
	 fferror("Could not build X column for REGFILTER");
	 return(-1);
      }
   }
   if( NodeY==-99 ) {
      type = ffGetVariable( "Y", &colVal );
      if( type==COLUMN ) {
	 NodeY = New_Column( (int)colVal.lng );
      } else {
	 fferror("Could not build Y column for REGFILTER");
	 return(-1);
      }
   }
   NodeX = New_Unary( DOUBLE, 0, NodeX );
   NodeY = New_Unary( DOUBLE, 0, NodeY );
   Node0 = Alloc_Node(); /* This will hold the Region Data */
   if( NodeX<0 || NodeY<0 || Node0<0 ) return(-1);

   if( ! (Test_Dims( NodeX, NodeY ) ) ) {
     fferror("Dimensions of REGFILTER arguments are not compatible");
     return (-1);
   }

   n = Alloc_Node();
   if( n >= 0 ) {
      this                 = gParse.Nodes + n;
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
      
      Copy_Dims(n, NodeX);
      if( SIZE(NodeX)<SIZE(NodeY) )  Copy_Dims(n, NodeY);

      /* Init Region node to be treated as a "constant" */

      that0                = gParse.Nodes + Node0;
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
	    fferror("Could not extract valid pair of column names from REGFILTER");
	    Free_Last_Node();
	    return( -1 );
	 }
	 fits_get_colnum( gParse.def_fptr, CASEINSEN, cX, &Xcol,
			  &gParse.status );
	 fits_get_colnum( gParse.def_fptr, CASEINSEN, cY, &Ycol,
			  &gParse.status );
	 if( gParse.status ) {
	    fferror("Could not locate columns indicated for WCS info");
	    Free_Last_Node();
	    return( -1 );
	 }

      } else {
	 /*  Try to find columns used in X/Y expressions  */
	 Xcol = Locate_Col( gParse.Nodes + NodeX );
	 Ycol = Locate_Col( gParse.Nodes + NodeY );
	 if( Xcol<0 || Ycol<0 ) {
	    fferror("Found multiple X/Y column references in REGFILTER");
	    Free_Last_Node();
	    return( -1 );
	 }
      }

      /*  Now, get the WCS info, if it exists, from the indicated columns  */
      wcs.exists = 0;
      if( Xcol>0 && Ycol>0 ) {
	 tstat = 0;
	 ffgtcs( gParse.def_fptr, Xcol, Ycol,
		 &wcs.xrefval, &wcs.yrefval,
		 &wcs.xrefpix, &wcs.yrefpix,
		 &wcs.xinc,    &wcs.yinc,
		 &wcs.rot,      wcs.type,
		 &tstat );
	 if( tstat==NO_WCS_KEY ) {
	    wcs.exists = 0;
	 } else if( tstat ) {
	    gParse.status = tstat;
	    Free_Last_Node();
	    return( -1 );
	 } else {
	    wcs.exists = 1;
	 }
      }

      /*  Read in Region file  */

      fits_read_rgnfile( fname, &wcs, &Rgn, &gParse.status );
      if( gParse.status ) {
	 Free_Last_Node();
	 return( -1 );
      }

      that0->value.data.ptr = Rgn;

      if( OPER(NodeX)==CONST_OP && OPER(NodeY)==CONST_OP )
	 this->DoOp( this );
   }

   return( n );
}

static int New_Vector( int subNode )
{
   Node *this, *that;
   int n;

   n = Alloc_Node();
   if( n >= 0 ) {
      this              = gParse.Nodes + n;
      that              = gParse.Nodes + subNode;
      this->type        = that->type;
      this->nSubNodes   = 1;
      this->SubNodes[0] = subNode;
      this->operation   = '{';
      this->DoOp        = Do_Vector;
   }

   return( n );
}

static int Close_Vec( int vecNode )
{
   Node *this;
   int n, nelem=0;

   this = gParse.Nodes + vecNode;
   for( n=0; n < this->nSubNodes; n++ ) {
      if( TYPE( this->SubNodes[n] ) != this->type ) {
	 this->SubNodes[n] = New_Unary( this->type, 0, this->SubNodes[n] );
	 if( this->SubNodes[n]<0 ) return(-1);
      }
      nelem += SIZE(this->SubNodes[n]);
   }
   this->value.naxis    = 1;
   this->value.nelem    = nelem;
   this->value.naxes[0] = nelem;

   return( vecNode );
}

static int Locate_Col( Node *this )
/*  Locate the TABLE column number of any columns in "this" calculation.  */
/*  Return ZERO if none found, or negative if more than 1 found.          */
{
   Node *that;
   int  i, col=0, newCol, nfound=0;
   
   if( this->nSubNodes==0
       && this->operation<=0 && this->operation!=CONST_OP )
      return gParse.colData[ - this->operation].colnum;

   for( i=0; i<this->nSubNodes; i++ ) {
      that = gParse.Nodes + this->SubNodes[i];
      if( that->operation>0 ) {
	 newCol = Locate_Col( that );
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
	 newCol = gParse.colData[- that->operation].colnum;
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

static int Test_Dims( int Node1, int Node2 )
{
   Node *that1, *that2;
   int valid, i;

   if( Node1<0 || Node2<0 ) return(0);

   that1 = gParse.Nodes + Node1;
   that2 = gParse.Nodes + Node2;

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

static void Copy_Dims( int Node1, int Node2 )
{
   Node *that1, *that2;
   int i;

   if( Node1<0 || Node2<0 ) return;

   that1 = gParse.Nodes + Node1;
   that2 = gParse.Nodes + Node2;

   that1->value.nelem = that2->value.nelem;
   that1->value.naxis = that2->value.naxis;
   for( i=0; i<that2->value.naxis; i++ )
      that1->value.naxes[i] = that2->value.naxes[i];
}

/********************************************************************/
/*    Routines for actually evaluating the expression start here    */
/********************************************************************/

void Evaluate_Parser( long firstRow, long nRows )
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

   gParse.firstRow = firstRow;
   gParse.nRows    = nRows;

   /*  Reset Column Nodes' pointers to point to right data and UNDEF arrays  */

   rowOffset = firstRow - gParse.firstDataRow;
   for( i=0; i<gParse.nNodes; i++ ) {
     if(    OPER(i) >  0 || OPER(i) == CONST_OP ) continue;

      column = -OPER(i);
      offset = gParse.varData[column].nelem * rowOffset;

      gParse.Nodes[i].value.undef = gParse.varData[column].undef + offset;

      switch( gParse.Nodes[i].type ) {
      case BITSTR:
	 gParse.Nodes[i].value.data.strptr =
	    (char**)gParse.varData[column].data + rowOffset;
	 gParse.Nodes[i].value.undef       = NULL;
	 break;
      case STRING:
	 gParse.Nodes[i].value.data.strptr = 
	    (char**)gParse.varData[column].data + rowOffset;
	 gParse.Nodes[i].value.undef = gParse.varData[column].undef + rowOffset;
	 break;
      case BOOLEAN:
	 gParse.Nodes[i].value.data.logptr = 
	    (char*)gParse.varData[column].data + offset;
	 break;
      case LONG:
	 gParse.Nodes[i].value.data.lngptr = 
	    (long*)gParse.varData[column].data + offset;
	 break;
      case DOUBLE:
	 gParse.Nodes[i].value.data.dblptr = 
	    (double*)gParse.varData[column].data + offset;
	 break;
      }
   }

   Evaluate_Node( gParse.resultNode );
}

static void Evaluate_Node( int thisNode )
    /**********************************************************************/
    /*  Recursively evaluate thisNode's subNodes, then call one of the    */
    /*  Do_<Action> functions pointed to by thisNode's DoOp element.      */
    /**********************************************************************/
{
   Node *this;
   int i;
   
   if( gParse.status ) return;

   this = gParse.Nodes + thisNode;
   if( this->operation>0 ) {  /* <=0 indicate constants and columns */
      i = this->nSubNodes;
      while( i-- ) {
	 Evaluate_Node( this->SubNodes[i] );
	 if( gParse.status ) return;
      }
      this->DoOp( this );
   }
}

static void Allocate_Ptrs( Node *this )
{
   long elem, row, size;

   if( this->type==BITSTR || this->type==STRING ) {

      this->value.data.strptr = (char**)malloc( gParse.nRows
						* sizeof(char*) );
      if( this->value.data.strptr ) {
	 this->value.data.strptr[0] = (char*)malloc( gParse.nRows
						     * (this->value.nelem+2)
						     * sizeof(char) );
	 if( this->value.data.strptr[0] ) {
	    row = 0;
	    while( (++row)<gParse.nRows ) {
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
	    gParse.status = MEMORY_ALLOCATION;
	    free( this->value.data.strptr );
	 }
      } else {
	 gParse.status = MEMORY_ALLOCATION;
      }

   } else {

      elem = this->value.nelem * gParse.nRows;
      switch( this->type ) {
      case DOUBLE:  size = sizeof( double ); break;
      case LONG:    size = sizeof( long   ); break;
      case BOOLEAN: size = sizeof( char   ); break;
      default:      size = 1;                break;
      }

      this->value.data.ptr = calloc(size+1, elem);

      if( this->value.data.ptr==NULL ) {
	 gParse.status = MEMORY_ALLOCATION;
      } else {
	 this->value.undef = (char *)this->value.data.ptr + elem*size;
      }
   }
}

static void Do_Unary( Node *this )
{
   Node *that;
   long elem;

   that = gParse.Nodes + this->SubNodes[0];

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

      Allocate_Ptrs( this );

      if( !gParse.status ) {

	 if( this->type!=BITSTR ) {
	    elem = gParse.nRows;
	    if( this->type!=STRING )
	       elem *= this->value.nelem;
	    while( elem-- )
	       this->value.undef[elem] = that->value.undef[elem];
	 }

	 elem = gParse.nRows * this->value.nelem;

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
	       elem = gParse.nRows;
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

static void Do_Offset( Node *this )
{
   Node *col;
   long fRow, nRowOverlap, nRowReload, rowOffset;
   long nelem, elem, offset, nRealElem;
   int status;

   col       = gParse.Nodes + this->SubNodes[0];
   rowOffset = gParse.Nodes[  this->SubNodes[1] ].value.data.lng;

   Allocate_Ptrs( this );

   fRow   = gParse.firstRow + rowOffset;
   if( this->type==STRING || this->type==BITSTR )
      nRealElem = 1;
   else
      nRealElem = this->value.nelem;

   nelem = nRealElem;

   if( fRow < gParse.firstDataRow ) {

      /* Must fill in data at start of array */

      nRowReload = gParse.firstDataRow - fRow;
      if( nRowReload > gParse.nRows ) nRowReload = gParse.nRows;
      nRowOverlap = gParse.nRows - nRowReload;

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

   } else if( fRow + gParse.nRows > gParse.firstDataRow + gParse.nDataRows ) {

      /* Must fill in data at end of array */

      nRowReload = (fRow+gParse.nRows) - (gParse.firstDataRow+gParse.nDataRows);
      if( nRowReload>gParse.nRows ) {
	 nRowReload = gParse.nRows;
      } else {
	 fRow = gParse.firstDataRow + gParse.nDataRows;
      }
      nRowOverlap = gParse.nRows - nRowReload;

      offset = nRowOverlap * nelem;

      /*  NULLify any values falling out of bounds  */

      elem = gParse.nRows * nelem;
      while( fRow+nRowReload>gParse.totalRows && nRowReload>0 ) {
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
      nRowOverlap = gParse.nRows;
      offset      = 0;

   }

   if( nRowReload>0 ) {
      switch( this->type ) {
      case BITSTR:
      case STRING:
	 status = (*gParse.loadData)( -col->operation, fRow, nRowReload,
				      this->value.data.strptr+offset,
				      this->value.undef+offset );
	 break;
      case BOOLEAN:
	 status = (*gParse.loadData)( -col->operation, fRow, nRowReload,
				      this->value.data.logptr+offset,
				      this->value.undef+offset );
	 break;
      case LONG:
	 status = (*gParse.loadData)( -col->operation, fRow, nRowReload,
				      this->value.data.lngptr+offset,
				      this->value.undef+offset );
	 break;
      case DOUBLE:
	 status = (*gParse.loadData)( -col->operation, fRow, nRowReload,
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
      elem = gParse.nRows * nelem;

   offset = nelem * rowOffset;
   while( nRowOverlap-- && !gParse.status ) {
      while( nelem-- && !gParse.status ) {
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

static void Do_BinOp_bit( Node *this )
{
   Node *that1, *that2;
   char *sptr1=NULL, *sptr2=NULL;
   int  const1, const2;
   long rows;

   that1 = gParse.Nodes + this->SubNodes[0];
   that2 = gParse.Nodes + this->SubNodes[1];

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

      Allocate_Ptrs( this );

      if( !gParse.status ) {
	 rows  = gParse.nRows;
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

static void Do_BinOp_str( Node *this )
{
   Node *that1, *that2;
   char *sptr1, *sptr2, null1=0, null2=0;
   int const1, const2, val;
   long rows;

   that1 = gParse.Nodes + this->SubNodes[0];
   that2 = gParse.Nodes + this->SubNodes[1];

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

      Allocate_Ptrs( this );

      if( !gParse.status ) {

	 rows = gParse.nRows;
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

static void Do_BinOp_log( Node *this )
{
   Node *that1, *that2;
   int vector1, vector2;
   char val1=0, val2=0, null1=0, null2=0;
   long rows, nelem, elem;

   that1 = gParse.Nodes + this->SubNodes[0];
   that2 = gParse.Nodes + this->SubNodes[1];

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
      rows  = gParse.nRows;
      nelem = this->value.nelem;
      elem  = this->value.nelem * rows;
      
      Allocate_Ptrs( this );
      
      if( !gParse.status ) {
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
      rows  = gParse.nRows;
      nelem = this->value.nelem;
      elem  = this->value.nelem * rows;

      Allocate_Ptrs( this );

      if( !gParse.status ) {
	
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

static void Do_BinOp_lng( Node *this )
{
   Node *that1, *that2;
   int  vector1, vector2;
   long val1=0, val2=0;
   char null1=0, null2=0;
   long rows, nelem, elem;

   that1 = gParse.Nodes + this->SubNodes[0];
   that2 = gParse.Nodes + this->SubNodes[1];

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

      case '%':
	 if( val2 ) this->value.data.lng = (val1 % val2);
	 else       fferror("Divide by Zero");
	 break;
      case '/': 
	 if( val2 ) this->value.data.lng = (val1 / val2); 
	 else       fferror("Divide by Zero");
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
      rows  = gParse.nRows;
      nelem = this->value.nelem;
      elem  = this->value.nelem * rows;
      
      Allocate_Ptrs( this );
      
      if( !gParse.status ) {
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

      rows  = gParse.nRows;
      nelem = this->value.nelem;
      elem  = this->value.nelem * rows;

      Allocate_Ptrs( this );

      while( rows-- && !gParse.status ) {
	 while( nelem-- && !gParse.status ) {
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

static void Do_BinOp_dbl( Node *this )
{
   Node   *that1, *that2;
   int    vector1, vector2;
   double val1=0.0, val2=0.0;
   char   null1=0, null2=0;
   long   rows, nelem, elem;

   that1 = gParse.Nodes + this->SubNodes[0];
   that2 = gParse.Nodes + this->SubNodes[1];

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
	 else       fferror("Divide by Zero");
	 break;
      case '/': 
	 if( val2 ) this->value.data.dbl = (val1 / val2); 
	 else       fferror("Divide by Zero");
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
      rows  = gParse.nRows;
      nelem = this->value.nelem;
      elem  = this->value.nelem * rows;
      
      Allocate_Ptrs( this );
      
      if( !gParse.status ) {
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

      rows  = gParse.nRows;
      nelem = this->value.nelem;
      elem  = this->value.nelem * rows;

      Allocate_Ptrs( this );

      while( rows-- && !gParse.status ) {
	 while( nelem-- && !gParse.status ) {
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
  double cd;
  static double deg = 0;
  double a, sdec, sra;
  
  if (deg == 0) deg = ((double)4)*atan((double)1)/((double)180);
  /* deg = 1.0; **** UNCOMMENT IF YOU WANT RADIANS */


  
/*
This (commented out) algorithm uses the Low of Cosines, which becomes
 unstable for angles less than 0.1 arcsec. 
 
  cd = sin(dec1*deg)*sin(dec2*deg) 
    + cos(dec1*deg)*cos(dec2*deg)*cos((ra1-ra2)*deg);
  if (cd < (-1)) cd = -1;
  if (cd > (+1)) cd = +1;
  return acos(cd)/deg;
*/

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






static double ran1()
{
  static double dval = 0.0;
  double rndVal;

  if (dval == 0.0) {
    if( rand()<32768 && rand()<32768 )
      dval =      32768.0;
    else
      dval = 2147483648.0;
  }

  rndVal = (double)rand();
  while( rndVal > dval ) dval *= 2.0;
  return rndVal/dval;
}

/* Gaussian deviate routine from Numerical Recipes */
static double gasdev()
{
  static int iset = 0;
  static double gset;
  double fac, rsq, v1, v2;

  if (iset == 0) {
    do {
      v1 = 2.0*ran1()-1.0;
      v2 = 2.0*ran1()-1.0;
      rsq = v1*v1 + v2*v2;
    } while (rsq >= 1.0 || rsq == 0.0);
    fac = sqrt(-2.0*log(rsq)/rsq);
    gset = v1*fac;
    iset = 1;
    return v2*fac;
  } else {
    iset = 0;
    return gset;
  }

}

/* lgamma function - from Numerical Recipes */

float gammaln(float xx)
     /* Returns the value ln Gamma[(xx)] for xx > 0. */
{
  /* 
     Internal arithmetic will be done in double precision, a nicety
     that you can omit if five-figure accuracy is good enough. */
  double x,y,tmp,ser;
  static double cof[6]={76.18009172947146,-86.50532032941677,
			24.01409824083091,-1.231739572450155,
			0.1208650973866179e-2,-0.5395239384953e-5};
  int j;
  y=x=xx;
  tmp=x+5.5;
  tmp -= (x+0.5)*log(tmp);
  ser=1.000000000190015;
  for (j=0;j<=5;j++) ser += cof[j]/++y;
  return (float) -tmp+log(2.5066282746310005*ser/x);
}

/* Poisson deviate - derived from Numerical Recipes */
static long poidev(double xm)
{
  static double sq, alxm, g, oldm = -1.0;
  static double pi = 0;
  double em, t, y;

  if (pi == 0) pi = ((double)4)*atan((double)1);

  if (xm < 20.0) {
    if (xm != oldm) {
      oldm = xm;
      g = exp(-xm);
    }
    em = -1;
    t = 1.0;
    do {
      em += 1;
      t *= ran1();
    } while (t > g);
  } else {
    if (xm != oldm) {
      oldm = xm;
      sq = sqrt(2.0*xm);
      alxm = log(xm);
      g = xm*alxm-gammaln( (float) (xm+1.0));
    }
    do {
      do {
	y = tan(pi*ran1());
	em = sq*y+xm;
      } while (em < 0.0);
      em = floor(em);
      t = 0.9*(1.0+y*y)*exp(em*alxm-gammaln( (float) (em+1.0) )-g);
    } while (ran1() > t);
  }

  /* Return integer version */
  return (long int) floor(em+0.5);
}

static void Do_Func( Node *this )
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
      theParams[i] = gParse.Nodes + this->SubNodes[i];
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
	      this->value.data.lng = poidev(pVals[0].data.dbl);
	    else
	      this->value.data.lng = poidev(pVals[0].data.lng);
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
	       fferror("Out of range argument to arcsin");
	    else
	       this->value.data.dbl = asin( dval );
	    break;
	 case acos_fct:
	    dval = pVals[0].data.dbl;
	    if( dval<-1.0 || dval>1.0 )
	       fferror("Out of range argument to arccos");
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
	       fferror("Out of range argument to log");
	    else
	       this->value.data.dbl = log( dval );
	    break;
	 case log10_fct:
	    dval = pVals[0].data.dbl;
	    if( dval<=0.0 )
	       fferror("Out of range argument to log10");
	    else
	       this->value.data.dbl = log10( dval );
	    break;
	 case sqrt_fct:
	    dval = pVals[0].data.dbl;
	    if( dval<0.0 )
	       fferror("Out of range argument to sqrt");
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
	   cstrmid(this->value.data.str, this->value.nelem, 
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

      Allocate_Ptrs( this );

      row  = gParse.nRows;
      elem = row * this->value.nelem;

      if( !gParse.status ) {
	 switch( this->operation ) {

	    /* Special functions with no arguments */

	 case row_fct:
	    while( row-- ) {
	       this->value.data.lngptr[row] = gParse.firstRow + row;
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
	 case rnd_fct:
	   while( elem-- ) {
	     this->value.data.dblptr[elem] = ran1();
	     this->value.undef[elem] = 0;
	    }
	    break;

	 case gasrnd_fct:
	    while( elem-- ) {
	       this->value.data.dblptr[elem] = gasdev();
	       this->value.undef[elem] = 0;
	    }
	    break;

	 case poirnd_fct:
	   if( theParams[0]->type==DOUBLE ) {
	      if (theParams[0]->operation == CONST_OP) {
		while( elem-- ) {
		  this->value.undef[elem] = (pVals[0].data.dbl < 0);
		  if (! this->value.undef[elem]) {
		    this->value.data.lngptr[elem] = poidev(pVals[0].data.dbl);
		  }
		} 
	      } else {
		while( elem-- ) {
		  this->value.undef[elem] = theParams[0]->value.undef[elem];
		  if (theParams[0]->value.data.dblptr[elem] < 0) 
		    this->value.undef[elem] = 1;
		  if (! this->value.undef[elem]) {
		    this->value.data.lngptr[elem] = 
		      poidev(theParams[0]->value.data.dblptr[elem]);
		  }
		} /* while */
	      } /* ! CONST_OP */
	   } else {
	     /* LONG */
	      if (theParams[0]->operation == CONST_OP) {
		while( elem-- ) {
		  this->value.undef[elem] = (pVals[0].data.lng < 0);
		  if (! this->value.undef[elem]) {
		    this->value.data.lngptr[elem] = poidev(pVals[0].data.lng);
		  }
		} 
	      } else {
		while( elem-- ) {
		  this->value.undef[elem] = theParams[0]->value.undef[elem];
		  if (theParams[0]->value.data.lngptr[elem] < 0) 
		    this->value.undef[elem] = 1;
		  if (! this->value.undef[elem]) {
		    this->value.data.lngptr[elem] = 
		      poidev(theParams[0]->value.data.lngptr[elem]);
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
		 fferror("Could not allocate temporary memory in median function");
		 free( this->value.data.ptr );
		 break;
	       }

	       for (irow=0; irow<row; irow++) {
		  long *p = mptr;
		  int nelem1 = nelem;
		  int count = 0;

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
		 fferror("Could not allocate temporary memory in median function");
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
		    if (cstrmid(this->value.data.strptr[row], len,
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
      } /* End if (!gParse.status) */
   } /* End non-constant operations */

   i = this->nSubNodes;
   while( i-- ) {
      if( theParams[i]->operation>0 ) {
	 /*  Currently only numeric params allowed  */
	 free( theParams[i]->value.data.ptr );
      }
   }
}

static void Do_Deref( Node *this )
{
   Node *theVar, *theDims[MAXDIMS];
   int  isConst[MAXDIMS], allConst;
   long dimVals[MAXDIMS];
   int  i, nDims;
   long row, elem, dsize;

   theVar = gParse.Nodes + this->SubNodes[0];

   i = nDims = this->nSubNodes-1;
   allConst = 1;
   while( i-- ) {
      theDims[i] = gParse.Nodes + this->SubNodes[i+1];
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

   Allocate_Ptrs( this );

   if( !gParse.status ) {

      if( allConst && theVar->value.naxis==nDims ) {

	 /* Dereference completely using constant indices */

	 elem = 0;
	 i    = nDims;
	 while( i-- ) {
	    if( dimVals[i]<1 || dimVals[i]>theVar->value.naxes[i] ) break;
	    elem = theVar->value.naxes[i]*elem + dimVals[i]-1;
	 }
	 if( i<0 ) {
	    for( row=0; row<gParse.nRows; row++ ) {
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
	    fferror("Index out of range");
	    free( this->value.data.ptr );
	 }
	 
      } else if( allConst && nDims==1 ) {
	 
	 /* Reduce dimensions by 1, using a constant index */
	 
	 if( dimVals[0] < 1 ||
	     dimVals[0] > theVar->value.naxes[ theVar->value.naxis-1 ] ) {
	    fferror("Index out of range");
	    free( this->value.data.ptr );
	 } else if ( this->type == BITSTR || this->type == STRING ) {
	    elem = this->value.nelem * (dimVals[0]-1);
	    for( row=0; row<gParse.nRows; row++ ) {
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
	    for( row=0; row<gParse.nRows; row++ ) {
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

	 for( row=0; row<gParse.nRows; row++ ) {

	    for( i=0; i<nDims; i++ ) {
	       if( !isConst[i] ) {
		  if( theDims[i]->value.undef[row] ) {
		     fferror("Null encountered as vector index");
		     free( this->value.data.ptr );
		     break;
		  } else
		     dimVals[i] = theDims[i]->value.data.lngptr[row];
	       }
	    }
	    if( gParse.status ) break;

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
	       fferror("Index out of range");
	       free( this->value.data.ptr );
	    }
	 }

      } else {

	 /* Reduce dimensions by 1, using a nonconstant expression */

	 for( row=0; row<gParse.nRows; row++ ) {

	    /* Index cannot be a constant */

	    if( theDims[0]->value.undef[row] ) {
	       fferror("Null encountered as vector index");
	       free( this->value.data.ptr );
	       break;
	    } else
	       dimVals[0] = theDims[0]->value.data.lngptr[row];

	    if( dimVals[0] < 1 ||
		dimVals[0] > theVar->value.naxes[ theVar->value.naxis-1 ] ) {
	       fferror("Index out of range");
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

static void Do_GTI( Node *this )
{
   Node *theExpr, *theTimes;
   double *start, *stop, *times;
   long elem, nGTI, gti;
   int ordered;

   theTimes = gParse.Nodes + this->SubNodes[0];
   theExpr  = gParse.Nodes + this->SubNodes[1];

   nGTI    = theTimes->value.nelem;
   start   = theTimes->value.data.dblptr;
   stop    = theTimes->value.data.dblptr + nGTI;
   ordered = theTimes->type;

   if( theExpr->operation==CONST_OP ) {

      this->value.data.log = 
	 (Search_GTI( theExpr->value.data.dbl, nGTI, start, stop, ordered )>=0);
      this->operation      = CONST_OP;

   } else {

      Allocate_Ptrs( this );

      times = theExpr->value.data.dblptr;
      if( !gParse.status ) {

	 elem = gParse.nRows * this->value.nelem;
	 if( nGTI ) {
	    gti = -1;
	    while( elem-- ) {
	       if( (this->value.undef[elem] = theExpr->value.undef[elem]) )
		  continue;

            /*  Before searching entire GTI, check the GTI found last time  */
	       if( gti<0 || times[elem]<start[gti] || times[elem]>stop[gti] ) {
		  gti = Search_GTI( times[elem], nGTI, start, stop, ordered );
	       }
	       this->value.data.logptr[elem] = ( gti>=0 );
	    }
	 } else
	    while( elem-- ) {
	       this->value.data.logptr[elem] = 0;
	       this->value.undef[elem]       = 0;
	    }
      }
   }

   if( theExpr->operation>0 )
      free( theExpr->value.data.ptr );
}

static long Search_GTI( double evtTime, long nGTI, double *start,
			double *stop, int ordered )
{
   long gti, step;
                             
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
		  gti = -1L;
		  break;
	       }
	    } else if( evtTime<start[gti] ) {
	       if( evtTime<=stop[gti-1] )
		  gti -= step;
	       else {
		  gti = -1L;
		  break;
	       }
	    } else {
	       break;
	    }
	 }
      } else
	 gti = -1L;
      
   } else { /*  Use "SLOW" linear search  */
      gti = nGTI;
      while( gti-- )
	 if( evtTime>=start[gti] && evtTime<=stop[gti] )
	    break;
   }
   return( gti );
}

static void Do_REG( Node *this )
{
   Node *theRegion, *theX, *theY;
   double Xval=0.0, Yval=0.0;
   char   Xnull=0, Ynull=0;
   int    Xvector, Yvector;
   long   nelem, elem, rows;

   theRegion = gParse.Nodes + this->SubNodes[0];
   theX      = gParse.Nodes + this->SubNodes[1];
   theY      = gParse.Nodes + this->SubNodes[2];

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

      Allocate_Ptrs( this );

      if( !gParse.status ) {

	 rows  = gParse.nRows;
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

static void Do_Vector( Node *this )
{
   Node *that;
   long row, elem, idx, jdx, offset=0;
   int node;

   Allocate_Ptrs( this );

   if( !gParse.status ) {

      for( node=0; node<this->nSubNodes; node++ ) {

	 that = gParse.Nodes + this->SubNodes[node];

	 if( that->operation == CONST_OP ) {

	    idx = gParse.nRows*this->value.nelem + offset;
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
	       
	    row  = gParse.nRows;
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
       free( gParse.Nodes[this->SubNodes[node]].value.data.ptr );
}

/*****************************************************************************/
/*  Utility routines which perform the calculations on bits and SAO regions  */
/*****************************************************************************/

static char bitlgte(char *bits1, int oper, char *bits2)
{
 int val1, val2, nextbit;
 char result;
 int i, l1, l2, length, ldiff;
 char stream[256];
 char chr1, chr2;

 l1 = strlen(bits1);
 l2 = strlen(bits2);
 if (l1 < l2)
   {
    length = l2;
    ldiff = l2 - l1;
    i=0;
    while( ldiff-- ) stream[i++] = '0';
    while( l1--    ) stream[i++] = *(bits1++);
    stream[i] = '\0';
    bits1 = stream;
   }
 else if (l2 < l1)
   {
    length = l1;
    ldiff = l1 - l2;
    i=0;
    while( ldiff-- ) stream[i++] = '0';
    while( l2--    ) stream[i++] = *(bits2++);
    stream[i] = '\0';
    bits2 = stream;
   }
 else
    length = l1;

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
 return (result);
}

static void bitand(char *result,char *bitstrm1,char *bitstrm2)
{
 int i, l1, l2, ldiff;
 char stream[256];
 char chr1, chr2;

 l1 = strlen(bitstrm1);
 l2 = strlen(bitstrm2);
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
 *result = '\0';
}

static void bitor(char *result,char *bitstrm1,char *bitstrm2)
{
 int i, l1, l2, ldiff;
 char stream[256];
 char chr1, chr2;

 l1 = strlen(bitstrm1);
 l2 = strlen(bitstrm2);
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
 int i, l1, l2, ldiff;
 char stream[256];
 char chr1, chr2;

 l1 = strlen(bitstrm1);
 l2 = strlen(bitstrm2);
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
	  return( 0 );
    }
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
int cstrmid(char *dest_str, int dest_len,
	    char *src_str,  int src_len,
	    int pos)
{
  /* char fill_char = ' '; */
  char fill_char = '\0';
  if (src_len == 0) { src_len = strlen(src_str); } /* .. if constant */

  /* Fill destination with blanks */
  if (pos < 0) { 
    fferror("STRMID(S,P,N) P must be 0 or greater");
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


static void fferror(char *s)
{
    char msg[80];

    if( !gParse.status ) gParse.status = PARSE_SYNTAX_ERR;

    strncpy(msg, s, 80);
    msg[79] = '\0';
    ffpmsg(msg);
}
