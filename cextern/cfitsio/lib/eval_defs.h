#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#if defined(__sgi) || defined(__hpux)
#include <alloca.h>
#endif
#ifdef sparc
#include <malloc.h>
#endif
#include "fitsio2.h"

#define MAXDIMS       5
#define MAXSUBS      10
#define MAXVARNAME   80
#define CONST_OP  -1000
#define pERROR       -1
#define MAX_STRLEN  256
#define MAX_STRLEN_S "255"

typedef struct ParseData_struct ParseData;
typedef void* yyscan_t;
#ifndef FFBISON
#include "eval_tab.h"
#endif


typedef struct {
                  char   name[MAXVARNAME+1];
                  int    type;
                  long   nelem;
                  int    naxis;
                  long   naxes[MAXDIMS];
                  char   *undef;
                  void   *data;
                                } DataInfo;

typedef struct {
                  long   nelem;
                  int    naxis;
                  long   naxes[MAXDIMS];
                  char   *undef;
                  union {
                         double dbl;
                         long   lng;
                         char   log;
                         char   str[MAX_STRLEN];
                         double *dblptr;
                         long   *lngptr;
                         char   *logptr;
                         char   **strptr;
                         void   *ptr;
		  } data;
                                } lval;

typedef struct Node {
                  int    operation;
  		  void   (*DoOp)(ParseData *, struct Node *this);
                  int    nSubNodes;
                  int    SubNodes[MAXSUBS];
                  int    type;
                  lval   value;
                                } Node;

struct ParseData_struct {
                  fitsfile    *def_fptr;
  		  int         (*getData)( ParseData *, char *dataName, void *dataValue );
  		  int         (*loadData)( ParseData *, int varNum, long fRow, long nRows,
					   void *data, char *undef );

                  int         compressed;
                  int         timeCol;
                  int         parCol;
                  int         valCol;

                  char        *expr;
                  int         index;
                  int         is_eobuf;

                  Node        *Nodes;
                  int         nNodes;
                  int         nNodesAlloc;
                  int         resultNode;
                  
                  long        firstRow;
                  long        nRows;

                  int         nCols;
                  long 	      nElements;
                  int         nAxis;
                  long        nAxes[MAXDIMS];
                  iteratorCol *colData;
                  DataInfo    *varData;
                  PixelFilter *pixFilter;

                  long        firstDataRow;
                  long        nDataRows;
                  long        totalRows;
                  long        nPrevDataRows;

                  int         datatype;
                  int         hdutype;

                  int         status;
};

typedef enum {
                  rnd_fct = 1001,
                  sum_fct,
                  nelem_fct,
                  sin_fct,
                  cos_fct,
                  tan_fct,
                  asin_fct,
                  acos_fct,
                  atan_fct,
                  sinh_fct,
                  cosh_fct,
                  tanh_fct,
                  exp_fct,
                  log_fct,
                  log10_fct,
                  sqrt_fct,
                  abs_fct,
                  atan2_fct,
                  ceil_fct,
                  floor_fct,
                  round_fct,
		  min1_fct,
		  min2_fct,
		  max1_fct,
		  max2_fct,
                  near_fct,
                  circle_fct,
                  box_fct,
                  elps_fct,
                  isnull_fct,
                  defnull_fct,
                  gtifilt_fct,
                  regfilt_fct,
                  ifthenelse_fct,
                  row_fct,
                  null_fct,
		  median_fct,
		  average_fct,
		  stddev_fct,
		  nonnull_fct,
		  angsep_fct,
		  gasrnd_fct,
		  poirnd_fct,
		  strmid_fct,
		  strpos_fct,
		  setnull_fct,
		  gtiover_fct,
		  gtifind_fct,
		  elemnum_fct,
		  axiselem_fct,
		  array_fct
                                } funcOp;


typedef struct parseInfo_struct parseInfo;

struct ParseStatusVariables { /* These variables were 'static' in fits_parse_workfn() */
  void *Data, *Null;
  int  datasize;
  long lastRow, repeat, resDataSize;
  LONGLONG jnull;
  parseInfo *userInfo;
  long zeros[4];
};

struct parseInfo_struct {
     int  datatype;   /* Data type to cast parse results into for user       */
     void *dataPtr;   /* Pointer to array of results, NULL if to use iterCol */
     void *nullPtr;   /* Pointer to nulval, use zero if NULL                 */
     long maxRows;    /* Max No. of rows to process, -1=all, 0=1 iteration   */
     int  anyNull;    /* Flag indicating at least 1 undef value encountered  */
     ParseData *parseData; /* Pointer to parser configuration */
     struct ParseStatusVariables parseVariables;
};

#ifdef __cplusplus
extern "C" {
#endif

/* Not sure why this is needed but it is */
#define YYSTYPE FITS_PARSER_YYSTYPE
/* How ParseData is accessed from the lexer, i.e. by yyextra */
#define YY_EXTRA_TYPE ParseData *

   int  fits_parser_yyparse(yyscan_t yyscaner, ParseData *lParse);
   int  fits_parser_yylex(FITS_PARSER_YYSTYPE *, yyscan_t yyscanner);
   void fits_parser_yyrestart(FILE*, yyscan_t yyscanner);
   int  fits_parser_yylex_init_extra ( YY_EXTRA_TYPE user_defined, yyscan_t* scanner);
   int  fits_parser_yylex_destroy (yyscan_t scanner);

   void Evaluate_Parser( ParseData *lParse, long firstRow, long nRows );
   int  fits_parser_allocateCol( ParseData *lParse, int nCol, int *status );
   int fits_parser_set_temporary_col(ParseData *lParse, parseInfo *Info,
				     long int nrows, void *nulval, int *status);

#ifdef __cplusplus
    }
#endif
