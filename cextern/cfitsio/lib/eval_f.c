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
/*   Peter D Wilson   Aug 1999  Add row-offset capability               */
/*   Peter D Wilson   Sep 1999  Add row-range capability to ffcalc_rng  */
/*                                                                      */
/************************************************************************/

#include <limits.h>
#include <ctype.h>
#include "eval_defs.h"
#include "region.h"


/*  Internal routines needed to allow the evaluator to operate on FITS data  */

static void Setup_DataArrays( ParseData *lParse, int nCols, iteratorCol *cols,
                              long fRow, long nRows );
static int  find_column( ParseData *lParse, char *colName, void *itslval );
static int  find_keywd ( ParseData *lParse, char *key,     void *itslval );
static int  load_column( ParseData *lParse, int varNum, long fRow, long nRows,
                         void *data, char *undef );

static int DEBUG_PIXFILTER;

#define FREE(x) { if (x) free(x); else printf("invalid free(" #x ") at %s:%d\n", __FILE__, __LINE__); }

/*---------------------------------------------------------------------------*/
int fffrow( fitsfile *fptr,         /* I - Input FITS file                   */
            char     *expr,         /* I - Boolean expression                */
            long     firstrow,      /* I - First row of table to eval        */
            long     nrows,         /* I - Number of rows to evaluate        */
            long     *n_good_rows,  /* O - Number of rows eval to True       */
            char     *row_status,   /* O - Array of boolean results          */
            int      *status )      /* O - Error status                      */
/*                                                                           */
/* Evaluate a boolean expression using the indicated rows, returning an      */
/* array of flags indicating which rows evaluated to TRUE/FALSE              */
/*---------------------------------------------------------------------------*/
{
   parseInfo Info;
   int naxis, constant;
   long nelem, naxes[MAXDIMS], elem;
   char result;
   ParseData lParse;

   if( *status ) return( *status );
   memset(&Info, 0, sizeof(Info));   

   if( ffiprs( fptr, 0, expr, MAXDIMS, &Info.datatype, &nelem, &naxis,
               naxes, &lParse, status ) ) {
      ffcprs(&lParse);
      return( *status );
   }
   if( nelem<0 ) {
      constant = 1;
      nelem = -nelem;
   } else
      constant = 0;

   if( Info.datatype!=TLOGICAL || nelem!=1 ) {
      ffcprs(&lParse);
      ffpmsg("Expression does not evaluate to a logical scalar.");
      return( *status = PARSE_BAD_TYPE );
   }

   if( constant ) { /* No need to call parser... have result from ffiprs */
      result = lParse.Nodes[lParse.resultNode].value.data.log;
      *n_good_rows = nrows;
      for( elem=0; elem<nrows; elem++ )
         row_status[elem] = result;
   } else {
      firstrow     = (firstrow>1 ? firstrow : 1);
      Info.dataPtr = row_status;
      Info.nullPtr = NULL;
      Info.maxRows = nrows;
      Info.parseData = &lParse;

      if( ffiter( lParse.nCols, lParse.colData, firstrow-1, 0,
                  fits_parser_workfn, (void*)&Info, status ) == -1 )
         *status = 0;  /* -1 indicates exitted without error before end... OK */

      if( *status ) {

         /***********************/
         /* Error... Do nothing */
         /***********************/

      } else {

         /***********************************/
         /* Count number of good rows found */
         /***********************************/

         *n_good_rows = 0L;
         for( elem=0; elem<Info.maxRows; elem++ ) {
            if( row_status[elem]==1 ) ++*n_good_rows;
         }
      }
   }

   ffcprs(&lParse);
   return(*status);
}

/*--------------------------------------------------------------------------*/
int ffsrow( fitsfile *infptr,   /* I - Input FITS file                      */
            fitsfile *outfptr,  /* I - Output FITS file                     */
            char     *expr,     /* I - Boolean expression                   */
            int      *status )  /* O - Error status                         */
/*                                                                          */
/* Evaluate an expression on all rows of a table.  If the input and output  */
/* files are not the same, copy the TRUE rows to the output file.  If the   */
/* files are the same, delete the FALSE rows (preserve the TRUE rows).      */
/* Can copy rows between extensions of the same file, *BUT* if output       */
/* extension is before the input extension, the second extension *MUST* be  */
/* opened using ffreopen, so that CFITSIO can handle changing file lengths. */
/*--------------------------------------------------------------------------*/
{
   parseInfo Info;
   int naxis, constant;
   long nelem, rdlen, naxes[MAXDIMS], maxrows, nbuff, nGood, inloc, outloc;
   LONGLONG ntodo, inbyteloc, outbyteloc, hsize;
   long freespace;
   unsigned char *buffer, result;
   struct {
      LONGLONG rowLength, numRows, heapSize;
      LONGLONG dataStart, heapStart;
   } inExt, outExt;
   ParseData lParse;

   if( *status ) return( *status );

   memset(&Info, 0, sizeof(Info));   
   memset(&inExt, 0, sizeof(inExt));
   memset(&outExt, 0, sizeof(outExt));

   if( ffiprs( infptr, 0, expr, MAXDIMS, &Info.datatype, &nelem, &naxis,
               naxes, &lParse, status ) ) {
      ffcprs(&lParse);
      return( *status );
   }

   if( nelem<0 ) {
      constant = 1;
      nelem = -nelem;
   } else
      constant = 0;

   /**********************************************************************/
   /* Make sure expression evaluates to the right type... logical scalar */
   /**********************************************************************/

   if( Info.datatype!=TLOGICAL || nelem!=1 ) {
      ffcprs(&lParse);
      ffpmsg("Expression does not evaluate to a logical scalar.");
      return( *status = PARSE_BAD_TYPE );
   }

   /***********************************************************/
   /*  Extract various table information from each extension  */
   /***********************************************************/

   if( infptr->HDUposition != (infptr->Fptr)->curhdu )
      ffmahd( infptr, (infptr->HDUposition) + 1, NULL, status );
   if( *status ) {
      ffcprs(&lParse);
      return( *status );
   }
   inExt.rowLength = (long) (infptr->Fptr)->rowlength;
   inExt.numRows   = (infptr->Fptr)->numrows;
   inExt.heapSize  = (infptr->Fptr)->heapsize;
   if( inExt.numRows == 0 ) { /* Nothing to copy */
      ffcprs(&lParse);
      return( *status );
   }

   if( outfptr->HDUposition != (outfptr->Fptr)->curhdu )
      ffmahd( outfptr, (outfptr->HDUposition) + 1, NULL, status );
   if( (outfptr->Fptr)->datastart < 0 )
      ffrdef( outfptr, status );
   if( *status ) {
      ffcprs(&lParse);
      return( *status );
   }
   outExt.rowLength = (long) (outfptr->Fptr)->rowlength;
   outExt.numRows   = (outfptr->Fptr)->numrows;
   if( !outExt.numRows )
      (outfptr->Fptr)->heapsize = 0L;
   outExt.heapSize  = (outfptr->Fptr)->heapsize;

   if( inExt.rowLength != outExt.rowLength ) {
      ffpmsg("Output table has different row length from input");
      ffcprs(&lParse);
      return( *status = PARSE_BAD_OUTPUT );
   }

   /***********************************/
   /*  Fill out Info data for parser  */
   /***********************************/

   Info.dataPtr = (char *)malloc( (size_t) ((inExt.numRows + 1) * sizeof(char)) );
   Info.nullPtr = NULL;
   Info.maxRows = (long) inExt.numRows;
   Info.parseData = &lParse;
   if( !Info.dataPtr ) {
      ffpmsg("Unable to allocate memory for row selection");
      ffcprs(&lParse);
      return( *status = MEMORY_ALLOCATION );
   }
   
   /* make sure array is zero terminated */
   ((char*)Info.dataPtr)[inExt.numRows] = 0;

   if( constant ) { /*  Set all rows to the same value from constant result  */

      result = lParse.Nodes[lParse.resultNode].value.data.log;
      for( ntodo = 0; ntodo<inExt.numRows; ntodo++ )
         ((char*)Info.dataPtr)[ntodo] = result;
      nGood = (long) (result ? inExt.numRows : 0);

   } else {

      ffiter( lParse.nCols, lParse.colData, 0L, 0L,
              fits_parser_workfn, (void*)&Info, status );

      nGood = 0;
      for( ntodo = 0; ntodo<inExt.numRows; ntodo++ )
         if( ((char*)Info.dataPtr)[ntodo] ) nGood++;
   }

   if( *status ) {
      /* Error... Do nothing */
   } else {
      rdlen  = (long) inExt.rowLength;
      buffer = (unsigned char *)malloc(maxvalue(500000,rdlen) * sizeof(char) );
      if( buffer==NULL ) {
         ffcprs(&lParse);
         return( *status=MEMORY_ALLOCATION );
      }
      maxrows = maxvalue( (500000L/rdlen), 1);
      nbuff = 0;
      inloc = 1;
      if( infptr==outfptr ) { /* Skip initial good rows if input==output file */
         while( ((char*)Info.dataPtr)[inloc-1] ) inloc++;
         outloc = inloc;
      } else {
         outloc = (long) (outExt.numRows + 1);
         if (outloc > 1) 
            ffirow( outfptr, outExt.numRows, nGood, status );
      }

      do {
         if( ((char*)Info.dataPtr)[inloc-1] ) {
            ffgtbb( infptr, inloc, 1L, rdlen, buffer+rdlen*nbuff, status );
            nbuff++;
            if( nbuff==maxrows ) {
               ffptbb( outfptr, outloc, 1L, rdlen*nbuff, buffer,  status );
               outloc += nbuff;
               nbuff = 0;
            }
         }
         inloc++;
      } while( !*status && inloc<=inExt.numRows );

      if( nbuff ) {
         ffptbb( outfptr, outloc, 1L, rdlen*nbuff, buffer,  status );
         outloc += nbuff;
      }

      if( infptr==outfptr ) {

         if( outloc<=inExt.numRows )
            ffdrow( infptr, outloc, inExt.numRows-outloc+1, status );

      } else if( inExt.heapSize && nGood ) {

         /* Copy heap, if it exists and at least one row copied */

         /********************************************************/
         /*  Get location information from the output extension  */
         /********************************************************/

         if( outfptr->HDUposition != (outfptr->Fptr)->curhdu )
            ffmahd( outfptr, (outfptr->HDUposition) + 1, NULL, status );
         outExt.dataStart = (outfptr->Fptr)->datastart;
         outExt.heapStart = (outfptr->Fptr)->heapstart;

         /*************************************************/
         /*  Insert more space into outfptr if necessary  */
         /*************************************************/

         hsize     = outExt.heapStart + outExt.heapSize;
         freespace = (long) (( ( (hsize + 2879) / 2880) * 2880) - hsize);
         ntodo     = inExt.heapSize;

         if ( (freespace - ntodo) < 0) {       /* not enough existing space? */
            ntodo = (ntodo - freespace + 2879) / 2880;  /* number of blocks  */
            ffiblk(outfptr, (long) ntodo, 1, status);   /* insert the blocks */
         }
         ffukyj( outfptr, "PCOUNT", inExt.heapSize+outExt.heapSize,
                 NULL, status );

         /*******************************************************/
         /*  Get location information from the input extension  */
         /*******************************************************/

         if( infptr->HDUposition != (infptr->Fptr)->curhdu )
            ffmahd( infptr, (infptr->HDUposition) + 1, NULL, status );
         inExt.dataStart = (infptr->Fptr)->datastart;
         inExt.heapStart = (infptr->Fptr)->heapstart;

         /**********************************/
         /*  Finally copy heap to outfptr  */
         /**********************************/

         ntodo  =  inExt.heapSize;
         inbyteloc  =  inExt.heapStart +  inExt.dataStart;
         outbyteloc = outExt.heapStart + outExt.dataStart + outExt.heapSize;

         while ( ntodo && !*status ) {
            rdlen = (long) minvalue(ntodo,500000);
            ffmbyt( infptr,  inbyteloc,  REPORT_EOF, status );
            ffgbyt( infptr,  rdlen,  buffer,     status );
            ffmbyt( outfptr, outbyteloc, IGNORE_EOF, status );
            ffpbyt( outfptr, rdlen,  buffer,     status );
            inbyteloc  += rdlen;
            outbyteloc += rdlen;
            ntodo  -= rdlen;
         }

         /***********************************************************/
         /*  But must update DES if data is being appended to a     */
         /*  pre-existing heap space.  Edit each new entry in file  */
         /***********************************************************/

         if( outExt.heapSize ) {
            LONGLONG repeat, offset, j;
            int i;
            for( i=1; i<=(outfptr->Fptr)->tfield; i++ ) {
               if( (outfptr->Fptr)->tableptr[i-1].tdatatype<0 ) {
                  for( j=outExt.numRows+1; j<=outExt.numRows+nGood; j++ ) {
                     ffgdesll( outfptr, i, j, &repeat, &offset, status );
                     offset += outExt.heapSize;
                     ffpdes( outfptr, i, j, repeat, offset, status );
                  }
               }
            }
         }

      } /*  End of HEAP copy  */

      FREE(buffer);
   }

   FREE(Info.dataPtr);
   ffcprs(&lParse);

   ffcmph(outfptr, status);  /* compress heap, deleting any orphaned data */
   return(*status);
}

/*---------------------------------------------------------------------------*/
int ffcrow( fitsfile *fptr,      /* I - Input FITS file                      */
            int      datatype,   /* I - Datatype to return results as        */
            char     *expr,      /* I - Arithmetic expression                */
            long     firstrow,   /* I - First row to evaluate                */
            long     nelements,  /* I - Number of elements to return         */
            void     *nulval,    /* I - Ptr to value to use as UNDEF         */
            void     *array,     /* O - Array of results                     */
            int      *anynul,    /* O - Were any UNDEFs encountered?         */
            int      *status )   /* O - Error status                         */
/*                                                                           */
/* Calculate an expression for the indicated rows of a table, returning      */
/* the results, cast as datatype (TSHORT, TDOUBLE, etc), in array.  If       */
/* nulval==NULL, UNDEFs will be zeroed out.  For vector results, the number  */
/* of elements returned may be less than nelements if nelements is not an    */
/* even multiple of the result dimension.  Call fftexp to obtain the         */
/* dimensions of the results.                                                */
/*---------------------------------------------------------------------------*/
{
   parseInfo Info;
   int naxis;
   long nelem1, naxes[MAXDIMS];
   ParseData lParse;

   if( *status ) return( *status );

   memset(&Info, 0, sizeof(Info));   

   if( ffiprs( fptr, 0, expr, MAXDIMS, &Info.datatype, &nelem1, &naxis,
               naxes, &lParse, status ) ) {
      ffcprs(&lParse);
      return( *status );
   }
   if( nelem1<0 ) nelem1 = - nelem1;

   if( nelements<nelem1 ) {
      ffcprs(&lParse);
      ffpmsg("Array not large enough to hold at least one row of data.");
      return( *status = PARSE_LRG_VECTOR );
   }

   firstrow = (firstrow>1 ? firstrow : 1);

   if( datatype ) Info.datatype = datatype;

   Info.dataPtr = array;
   Info.nullPtr = nulval;
   Info.maxRows = nelements / nelem1;
   Info.parseData = &lParse;
   
   if( ffiter( lParse.nCols, lParse.colData, firstrow-1, 0,
               fits_parser_workfn, (void*)&Info, status ) == -1 )
      *status=0;  /* -1 indicates exitted without error before end... OK */

   *anynul = Info.anyNull;
   ffcprs(&lParse);
   return( *status );
}

/*--------------------------------------------------------------------------*/
int ffcalc( fitsfile *infptr,   /* I - Input FITS file                      */
            char     *expr,     /* I - Arithmetic expression                */
            fitsfile *outfptr,  /* I - Output fits file                     */
            char     *parName,  /* I - Name of output parameter             */
            char     *parInfo,  /* I - Extra information on parameter       */
            int      *status )  /* O - Error status                         */
/*                                                                          */
/* Evaluate an expression for all rows of a table.  Call ffcalc_rng with    */
/* a row range of 1-MAX.                                                    */
{
   long start=1, end=LONG_MAX;

   return ffcalc_rng( infptr, expr, outfptr, parName, parInfo,
                      1, &start, &end, status );
}

/*--------------------------------------------------------------------------*/
int ffcalc_rng( fitsfile *infptr,   /* I - Input FITS file                  */
                char     *expr,     /* I - Arithmetic expression            */
                fitsfile *outfptr,  /* I - Output fits file                 */
                char     *parName,  /* I - Name of output parameter         */
                char     *parInfo,  /* I - Extra information on parameter   */
                int      nRngs,     /* I - Row range info                   */
                long     *start,    /* I - Row range info                   */
                long     *end,      /* I - Row range info                   */
                int      *status )  /* O - Error status                     */
/*                                                                          */
/* Evaluate an expression using the data in the input FITS file and place   */
/* the results into either a column or keyword in the output fits file,     */
/* depending on the value of parName (keywords normally prefixed with '#')  */
/* and whether the expression evaluates to a constant or a table column.    */
/* The logic is as follows:                                                 */
/*    (1) If a column exists with name, parName, put results there.         */
/*    (2) If parName starts with '#', as in #NAXIS, put result there,       */
/*        with parInfo used as the comment. If expression does not evaluate */
/*        to a constant, flag an error.                                     */
/*    (3) If a keyword exists with name, parName, and expression is a       */
/*        constant, put result there, using parInfo as the new comment.     */
/*    (4) Else, create a new column with name parName and TFORM parInfo.    */
/*        If parInfo is NULL, use a default data type for the column.       */
/*--------------------------------------------------------------------------*/
{
   parseInfo Info;
   int naxis, constant, typecode, newNullKwd=0;
   long nelem, naxes[MAXDIMS], repeat, width;
   int col_cnt, colNo;
   Node *result;
   char card[81], tform[16], nullKwd[9], tdimKwd[9];
   ParseData lParse;

   if( *status ) return( *status );

   memset(&Info, 0, sizeof(Info));   

   if( ffiprs( infptr, 0, expr, MAXDIMS, &Info.datatype, &nelem, &naxis,
               naxes, &lParse, status ) ) {

      ffcprs(&lParse);
      return( *status );
   }
   if( nelem<0 ) {
      constant = 1;
      nelem = -nelem;
   } else
      constant = 0;

   Info.parseData = &lParse;
   /*  Case (1): If column exists put it there  */

   colNo = 0;
   if( ffgcno( outfptr, CASEINSEN, parName, &colNo, status )==COL_NOT_FOUND ) {

      /*  Output column doesn't exist.  Test for keyword. */

      /* Case (2): Does parName indicate result should be put into keyword */

      *status = 0;
      if( parName[0]=='#' ) {
         if( ! constant ) {
            ffcprs(&lParse);
            ffpmsg( "Cannot put tabular result into keyword (ffcalc)" );
            return( *status = PARSE_BAD_TYPE );
         }
         parName++;  /* Advance past '#' */
	 if ( (fits_strcasecmp(parName,"HISTORY") == 0 || fits_strcasecmp(parName,"COMMENT") == 0) &&
	      Info.datatype != TSTRING ) {
            ffcprs(&lParse);
            ffpmsg( "HISTORY and COMMENT values must be strings (ffcalc)" );
	    return( *status = PARSE_BAD_TYPE );
	 }

      } else if( constant ) {

         /* Case (3): Does a keyword named parName already exist */

         if( ffgcrd( outfptr, parName, card, status )==KEY_NO_EXIST ) {
            colNo = -1;
         } else if( *status ) {
            ffcprs(&lParse);
            return( *status );
         }

      } else
         colNo = -1;

      if( colNo<0 ) {

         /* Case (4): Create new column */

         *status = 0;
         ffgncl( outfptr, &colNo, status );
         colNo++;
         if( parInfo==NULL || *parInfo=='\0' ) {
            /*  Figure out best default column type  */
            if( lParse.hdutype==BINARY_TBL ) {
               snprintf(tform,15,"%ld",nelem);
               switch( Info.datatype ) {
               case TLOGICAL:  strcat(tform,"L");  break;
               case TLONG:     strcat(tform,"J");  break;
               case TDOUBLE:   strcat(tform,"D");  break;
               case TSTRING:   strcat(tform,"A");  break;
               case TBIT:      strcat(tform,"X");  break;
               case TLONGLONG: strcat(tform,"K");  break;
               }
            } else {
               switch( Info.datatype ) {
               case TLOGICAL:
                  ffcprs(&lParse);
                  ffpmsg("Cannot create LOGICAL column in ASCII table");
                  return( *status = NOT_BTABLE );
               case TLONG:     strcpy(tform,"I11");     break;
               case TDOUBLE:   strcpy(tform,"D23.15");  break;
               case TSTRING:   
               case TBIT:      snprintf(tform,16,"A%ld",nelem);  break;
               }
            }
            parInfo = tform;
         } else if( !(isdigit((int) *parInfo)) && lParse.hdutype==BINARY_TBL ) {
            if( Info.datatype==TBIT && *parInfo=='B' )
               nelem = (nelem+7)/8;
            snprintf(tform,16,"%ld%s",nelem,parInfo);
            parInfo = tform;
         }
         fficol( outfptr, colNo, parName, parInfo, status );
         if( naxis>1 )
            ffptdm( outfptr, colNo, naxis, naxes, status );

         /*  Setup TNULLn keyword in case NULLs are encountered  */

         ffkeyn("TNULL", colNo, nullKwd, status);
         if( ffgcrd( outfptr, nullKwd, card, status )==KEY_NO_EXIST ) {
            *status = 0;
            if( lParse.hdutype==BINARY_TBL ) {
	       LONGLONG nullVal=0;
               fits_binary_tform( parInfo, &typecode, &repeat, &width, status );
               if( typecode==TBYTE )
                  nullVal = UCHAR_MAX;
               else if( typecode==TSHORT )
                  nullVal = SHRT_MIN;
               else if( typecode==TINT )
                  nullVal = INT_MIN;
               else if( typecode==TLONG ) {
                  if (sizeof(long) == 8 && sizeof(int) == 4)
                     nullVal = INT_MIN;
                  else
                     nullVal = LONG_MIN;
               }
               else if( typecode==TLONGLONG )
                  nullVal = LONGLONG_MIN;
		  
               if( nullVal ) {
                  ffpkyj( outfptr, nullKwd, nullVal, "Null value", status );
                  fits_set_btblnull( outfptr, colNo, nullVal, status );
                  newNullKwd = 1;
               }
            } else if( lParse.hdutype==ASCII_TBL ) {
               ffpkys( outfptr, nullKwd, "NULL", "Null value string", status );
               fits_set_atblnull( outfptr, colNo, "NULL", status );
               newNullKwd = 1;
            }
         }

      }

   } else if( *status ) {
      ffcprs(&lParse);
      return( *status );
   } else {

      /********************************************************/
      /*  Check if a TDIM keyword should be written/updated.  */
      /********************************************************/

      ffkeyn("TDIM", colNo, tdimKwd, status);
      ffgcrd( outfptr, tdimKwd, card, status );
      if( *status==0 ) {
         /*  TDIM exists, so update it with result's dimension  */
         ffptdm( outfptr, colNo, naxis, naxes, status );
      } else if( *status==KEY_NO_EXIST ) {
         /*  TDIM does not exist, so clear error stack and     */
         /*  write a TDIM only if result is multi-dimensional  */
         *status = 0;
         ffcmsg();
         if( naxis>1 )
            ffptdm( outfptr, colNo, naxis, naxes, status );
      }
      if( *status ) {
         /*  Either some other error happened in ffgcrd   */
         /*  or one happened in ffptdm                    */
         ffcprs(&lParse);
         return( *status );
      }

   }

   if( colNo>0 ) {

      /*  Output column exists (now)... put results into it  */

      int anyNull = 0;
      int nPerLp, i;
      long totaln;

      ffgkyj(infptr, "NAXIS2", &totaln, 0, status);

      /*************************************/
      /* Create new iterator Output Column */
      /*************************************/

      col_cnt = lParse.nCols;
      if( fits_parser_allocateCol( &lParse, col_cnt, status ) ) {
         ffcprs(&lParse);
         return( *status );
      }

      fits_iter_set_by_num( lParse.colData+col_cnt, outfptr,
                            colNo, 0, OutputCol );
      lParse.nCols++;

      for( i=0; i<nRngs; i++ ) {
         Info.dataPtr = NULL;
         Info.maxRows = end[i]-start[i]+1;

          /*
            If there is only 1 range, and it includes all the rows,
            and there are 10 or more rows, then set nPerLp = 0 so
            that the iterator function will dynamically choose the
            most efficient number of rows to process in each loop.
            Otherwise, set nPerLp to the number of rows in this range.
         */

         if( (Info.maxRows >= 10) && (nRngs == 1) &&
             (start[0] == 1) && (end[0] == totaln))
              nPerLp = 0;
         else
              nPerLp = Info.maxRows;

         if( ffiter( lParse.nCols, lParse.colData, start[i]-1,
                     nPerLp, fits_parser_workfn, (void*)&Info, status ) == -1 )
            *status = 0;
         else if( *status ) {
            ffcprs(&lParse);
            return( *status );
         }
         if( Info.anyNull ) anyNull = 1;
      }

      if( newNullKwd && !anyNull ) {
         ffdkey( outfptr, nullKwd, status );
      }

   } else {

      /* Put constant result into keyword */

      result  = lParse.Nodes + lParse.resultNode;
      switch( Info.datatype ) {
      case TDOUBLE:
         ffukyd( outfptr, parName, result->value.data.dbl, 15,
                 parInfo, status );
         break;
      case TLONG:
         ffukyj( outfptr, parName, result->value.data.lng, parInfo, status );
         break;
      case TLOGICAL:
         ffukyl( outfptr, parName, result->value.data.log, parInfo, status );
         break;
      case TBIT:
      case TSTRING:
	 if (fits_strcasecmp(parName,"HISTORY") == 0) {
	   ffphis( outfptr, result->value.data.str, status);
	 } else if (fits_strcasecmp(parName,"COMMENT") == 0) {
	   ffpcom( outfptr, result->value.data.str, status);
	 } else {
	   ffukys( outfptr, parName, result->value.data.str, parInfo, status );
	 }
         break;
      }
   }

   ffcprs(&lParse);
   return( *status );
}

/*--------------------------------------------------------------------------*/
int fftexp( fitsfile *fptr,      /* I - Input FITS file                     */
            char     *expr,      /* I - Arithmetic expression               */
            int      maxdim,     /* I - Max Dimension of naxes              */
            int      *datatype,  /* O - Data type of result                 */
            long     *nelem,     /* O - Vector length of result             */
            int      *naxis,     /* O - # of dimensions of result           */
            long     *naxes,     /* O - Size of each dimension              */
            int      *status )   /* O - Error status                        */
/*                                                                          */
/* Evaluate the given expression and return information on the result.      */
/*--------------------------------------------------------------------------*/
{
   ParseData lParse;

   ffiprs( fptr, 0, expr, maxdim, datatype, nelem, naxis, naxes, &lParse, status );
   ffcprs(&lParse);
   return( *status );
}


/*--------------------------------------------------------------------------*/
int ffiprs( fitsfile *fptr,      /* I - Input FITS file                     */
            int      compressed, /* I - Is FITS file hkunexpanded?          */
            char     *expr,      /* I - Arithmetic expression               */
            int      maxdim,     /* I - Max Dimension of naxes              */
            int      *datatype,  /* O - Data type of result                 */
            long     *nelem,     /* O - Vector length of result             */
            int      *naxis,     /* O - # of dimensions of result           */
            long     *naxes,     /* O - Size of each dimension              */
	    ParseData *lParse,   /* O - parser status                       */
            int      *status )   /* O - Error status                        */
/*                                                                          */
/* Initialize the parser and determine what type of result the expression   */
/* produces.                                                                */
/*--------------------------------------------------------------------------*/
{
   Node *result;
   int  i,lexpr, tstatus = 0;
   int xaxis, bitpix;
   long xaxes[9];
   yyscan_t yylex_scanner; /* Used internally by FLEX lexer */
   PixelFilter *pixFilter = 0;

   if( *status ) return( *status );

   /* make sure all internal structures for this HDU are current */
   if ( ffrdef(fptr, status) ) return(*status);

   /*  Initialize the Parser structure  */

   /* Unfortunately we need to preserve the pixFilter value since it
      is pre-set when ffiprs() is called */
   pixFilter = lParse->pixFilter;
   memset(lParse, 0, sizeof(*lParse));
   lParse->pixFilter = pixFilter;

   lParse->def_fptr   = fptr;
   lParse->compressed = compressed;
   lParse->nCols      = 0;
   lParse->colData    = NULL;
   lParse->varData    = NULL;
   lParse->getData    = find_column;
   lParse->loadData   = load_column;
   lParse->Nodes      = NULL;
   lParse->nNodesAlloc= 0;
   lParse->nNodes     = 0;
   lParse->hdutype    = 0;
   lParse->status     = 0;

   fits_get_hdu_type(fptr, &(lParse->hdutype), status );

   if (lParse->hdutype == IMAGE_HDU) {

      fits_get_img_param(fptr, 9, &bitpix, &xaxis, xaxes, status);
      if (*status) {
         ffpmsg("ffiprs: unable to get image dimensions");
         return( *status );
      }
      lParse->totalRows = xaxis > 0 ? 1 : 0;
      for (i = 0; i < xaxis; ++i)
         lParse->totalRows *= xaxes[i];
      if (DEBUG_PIXFILTER)
         printf("naxis=%d, lParse->totalRows=%ld\n", xaxis, lParse->totalRows);
   }
   else if( ffgkyj(fptr, "NAXIS2", &lParse->totalRows, 0, &tstatus) )
   {
      /* this might be a 1D or null image with no NAXIS2 keyword */
      lParse->totalRows = 0;
   } 
   

   /*  Copy expression into parser... read from file if necessary  */


   if( expr[0]=='@' ) {
      if( ffimport_file( expr+1, &lParse->expr, status ) ) return( *status );
      lexpr = strlen(lParse->expr);
   } else {
      lexpr = strlen(expr);
      lParse->expr = (char*)malloc( (2+lexpr)*sizeof(char));
      strcpy(lParse->expr,expr);
   }
   strcat(lParse->expr + lexpr,"\n");
   lParse->index    = 0;
   lParse->is_eobuf = 0;

   /*  Parse the expression, building the Nodes and determing  */
   /*  which columns are needed and what data type is returned  */
   
   fits_parser_yylex_init_extra(lParse, &yylex_scanner);
   fits_parser_yyrestart(NULL, yylex_scanner);
   *status = fits_parser_yyparse(yylex_scanner, lParse);
   fits_parser_yylex_destroy(yylex_scanner);

   if( *status  ) return( *status = PARSE_SYNTAX_ERR );

   /*  Check results  */
   *status = lParse->status;
   if( *status ) return(*status);

   if( !lParse->nNodes ) {
      ffpmsg("Blank expression");
      return( *status = PARSE_SYNTAX_ERR );
   }
   if( !lParse->nCols ) {
     lParse->colData = (iteratorCol *) malloc(sizeof(iteratorCol));
     if (lParse->colData == 0) {
       ffpmsg("memory allocation failed (ffiprs)");
       return( *status = MEMORY_ALLOCATION );
     }
     /* This allows iterator to know value of */ 
     /* fptr when no columns are referenced   */
     memset(lParse->colData, 0, sizeof(iteratorCol));
     lParse->colData[0].fptr = fptr;
   }

   result = lParse->Nodes + lParse->resultNode;

   *naxis = lParse->nAxis     = result->value.naxis;
   *nelem = lParse->nElements = result->value.nelem;
   for( i=0; i<*naxis && i<maxdim; i++ )
      naxes[i] = lParse->nAxes[i] = result->value.naxes[i];

   switch( result->type ) {
   case BOOLEAN:
      *datatype = TLOGICAL;
      break;
   case LONG:
      *datatype = TLONG;
      break;
   case DOUBLE:
      *datatype = TDOUBLE;
      break;
   case BITSTR:
      *datatype = TBIT;
      break;
   case STRING:
      *datatype = TSTRING;
      break;
   default:
      *datatype = 0;
      ffpmsg("Bad return data type");
      *status = lParse->status = PARSE_BAD_TYPE;
      break;
   }
   lParse->datatype = *datatype;
   FREE(lParse->expr);

   if( result->operation==CONST_OP ) *nelem = - *nelem;
   return(*status);
}

/*--------------------------------------------------------------------------*/
void ffcprs( ParseData *lParse )
/*                                                                          */
/* Clear the parser, making it ready to accept a new expression.            */
/*--------------------------------------------------------------------------*/
{
   int col, node, i;

   if( lParse->nCols > 0 ) {
      FREE( lParse->colData  );
      for( col=0; col<lParse->nCols; col++ ) {
         if( lParse->varData[col].undef == NULL ) continue;
         if( lParse->varData[col].type  == BITSTR )
           FREE( ((char**)lParse->varData[col].data)[0] );
         free( lParse->varData[col].undef );
      }
      FREE( lParse->varData );
      lParse->nCols = 0;
   } else if ( lParse->colData ) {
     /* Special case if colData needed to be created with no columns */
     FREE( lParse->colData );
   }

   if( lParse->nNodes > 0 ) {
      node = lParse->nNodes;
      while( node-- ) {
         if( lParse->Nodes[node].operation==gtifilt_fct ) {
            i = lParse->Nodes[node].SubNodes[0];
            if (lParse->Nodes[ i ].value.data.ptr)
	        FREE( lParse->Nodes[ i ].value.data.ptr );
         }
         else if( lParse->Nodes[node].operation==regfilt_fct ) {
            i = lParse->Nodes[node].SubNodes[0];
            fits_free_region( (SAORegion *)lParse->Nodes[ i ].value.data.ptr );
         }
      }
      lParse->nNodes = 0;
   }
   if( lParse->Nodes ) free( lParse->Nodes );
   lParse->Nodes = NULL;

   lParse->hdutype = ANY_HDU;
   lParse->pixFilter = 0;
   lParse->nDataRows = lParse->nPrevDataRows = 0;
}

/*---------------------------------------------------------------------------*/
int fits_parser_workfn( long    totalrows,     /* I - Total rows to be processed     */
                long    offset,        /* I - Number of rows skipped at start*/
                long    firstrow,      /* I - First row of this iteration    */
                long    nrows,         /* I - Number of rows in this iter    */
                int      nCols,        /* I - Number of columns in use       */
                iteratorCol *colData,  /* IO- Column information/data        */
                void    *userPtr )     /* I - Data handling instructions     */
/*                                                                           */
/* Iterator work function which calls the parser and copies the results      */
/* into either an OutputCol or a data pointer supplied in the userPtr        */
/* structure.                                                                */
/*---------------------------------------------------------------------------*/
{
    int status, constant=0, anyNullThisTime=0;
    long jj, kk, idx, remain, ntodo;
    Node *result;
    iteratorCol * outcol;
    ParseData *lParse = ((parseInfo*)userPtr)->parseData;
    struct ParseStatusVariables *pv = &( ((parseInfo*)userPtr)->parseVariables );
    void *Data0 = 0;

    /* declare variables static to preserve their values between calls */
    long zeros[4] = {0,0,0,0};

    if (DEBUG_PIXFILTER)
       printf("fits_parser_workfn(total=%ld, offset=%ld, first=%ld, rows=%ld, cols=%d)\n",
                totalrows, offset, firstrow, nrows, nCols);
    /*--------------------------------------------------------*/
    /*  Initialization procedures: execute on the first call  */
    /*--------------------------------------------------------*/
    outcol = colData + (nCols - 1);
    if (firstrow == offset+1)
    {
       (pv->userInfo) = (parseInfo*)userPtr;
       (pv->userInfo)->anyNull = 0;

       /* Unfortunately there are two copies of the iterator columns,
	  one inside the parser and one outside maintained by the
	  higher level.  (This could happen if the histogramming
	  routines are binning multiple columns, and so there are
	  multiple parsers being managed at one time.) Upon the first
	  call we make sure they match */
       for (jj = 0; jj<nCols; jj++) {
	 lParse->colData[jj].repeat = colData[jj].repeat;
       }

       if( (pv->userInfo)->maxRows>0 )
          (pv->userInfo)->maxRows = minvalue(totalrows,(pv->userInfo)->maxRows);
       else if( (pv->userInfo)->maxRows<0 )
          (pv->userInfo)->maxRows = totalrows;
       else
          (pv->userInfo)->maxRows = nrows;

       (pv->lastRow) = firstrow + (pv->userInfo)->maxRows - 1;

       /* dataPtr == NULL indicates an iterator-derived column, which
	  means that the first value will be a null value and the remaining
	  values will be the where the outputs are placed */
       if( (pv->userInfo)->dataPtr==NULL ) {

          if( outcol->iotype == InputCol ) {
             ffpmsg("Output column for parser results not found!");
             return( PARSE_NO_OUTPUT );
          }
          /* Data gets set later */
          (pv->Null) = outcol->array;
          (pv->userInfo)->datatype = outcol->datatype;

          /* Check for a TNULL/BLANK keyword for output column/image */

          status = 0;
          (pv->jnull) = 0;
          if (lParse->hdutype == IMAGE_HDU) {
             if (lParse->pixFilter->blank)
                (pv->jnull) = (LONGLONG) lParse->pixFilter->blank;
          }
          else {
	    if (outcol->iotype != TemporaryCol) {
	      ffgknjj( outcol->fptr, "TNULL", outcol->colnum,
		       1, &(pv->jnull), (int*)&jj, &status );
	    }

             if( status==BAD_INTKEY || outcol->iotype == TemporaryCol) {
                /*  Probably ASCII table with text TNULL keyword  */
                switch( (pv->userInfo)->datatype ) {
                   case TSHORT:  (pv->jnull) = (LONGLONG) SHRT_MIN;      break;
                   case TINT:    (pv->jnull) = (LONGLONG) INT_MIN;       break;
                   case TLONG:   (pv->jnull) = (LONGLONG) LONG_MIN;      break;
                }
             }
          }
          (pv->repeat) = outcol->repeat;
/*
          if (DEBUG_PIXFILTER)
            printf("fits_parser_workfn: using null value %ld\n", (pv->jnull));
*/
       } else {

	  /* This clause applies if the user is passing user-allocated 
	     data arrays, which is where the data will be placed.  This 
	     means they should also be passing null values */
          (pv->Data) = (pv->userInfo)->dataPtr;
          (pv->Null) = ((pv->userInfo)->nullPtr ? (pv->userInfo)->nullPtr : zeros);
          (pv->repeat) = lParse->Nodes[lParse->resultNode].value.nelem;

       }

       /* Determine the size of each element of the returned result */

       switch( (pv->userInfo)->datatype ) {
       case TBIT:       /*  Fall through to TBYTE  */
       case TLOGICAL:   /*  Fall through to TBYTE  */
       case TBYTE:     (pv->datasize) = sizeof(char);     break;
       case TSHORT:    (pv->datasize) = sizeof(short);    break;
       case TINT:      (pv->datasize) = sizeof(int);      break;
       case TLONG:     (pv->datasize) = sizeof(long);     break;
       case TLONGLONG: (pv->datasize) = sizeof(LONGLONG); break;
       case TFLOAT:    (pv->datasize) = sizeof(float);    break;
       case TDOUBLE:   (pv->datasize) = sizeof(double);   break;
       case TSTRING:   (pv->datasize) = sizeof(char*);    break;
       }

       /* Determine the size of each element of the calculated result */
       /*   (only matters for numeric/logical data)                   */

       switch( lParse->Nodes[lParse->resultNode].type ) {
       case BOOLEAN:   (pv->resDataSize) = sizeof(char);    break;
       case LONG:      (pv->resDataSize) = sizeof(long);    break;
       case DOUBLE:    (pv->resDataSize) = sizeof(double);  break;
       }
    }

    /*-------------------------------------------*/
    /*  Main loop: process all the rows of data  */
    /*-------------------------------------------*/

    /*  If writing to output column, set first element to appropriate  */
    /*  null value.  If no NULLs encounter, zero out before returning. */
/*
          if (DEBUG_PIXFILTER)
            printf("fits_parser_workfn: using null value %ld\n", (pv->jnull));
*/

    if( (pv->userInfo)->dataPtr == NULL ) {
       /* First, reset Data pointer to start of output array, plus 1
	  because the 0th element is the null value (cute undocumented
	  feature of the iterator!) */
       (pv->Data) = (char*) outcol->array + (pv->datasize);

       /* A TemporaryCol with null value specified explicitly */
       if (outcol->iotype == TemporaryCol && (pv->userInfo)->nullPtr) {

	 pv->Null = (pv->userInfo)->nullPtr;

       } else {

	 /* ... or an OutputCol or TemporaryCol with no explicit null */
	 switch( (pv->userInfo)->datatype ) {
	 case TLOGICAL: *(char  *)(pv->Null) = 'U';             break;
	 case TBYTE:    *(char  *)(pv->Null) = (char )(pv->jnull);    break;
	 case TSHORT:   *(short *)(pv->Null) = (short)(pv->jnull);    break;
	 case TINT:     *(int   *)(pv->Null) = (int  )(pv->jnull);    break;
	 case TLONG:    *(long  *)(pv->Null) = (long )(pv->jnull);    break;
	 case TLONGLONG: *(LONGLONG  *)(pv->Null) = (LONGLONG )(pv->jnull);    break;
	 case TFLOAT:   *(float *)(pv->Null) = FLOATNULLVALUE;  break;
	 case TDOUBLE:  *(double*)(pv->Null) = DOUBLENULLVALUE; break;
	 case TSTRING: (*(char **)(pv->Null))[0] = '\1';
	               (*(char **)(pv->Null))[1] = '\0';        break;
	 }
       }
    }

    /* Alter nrows in case calling routine didn't want to do all rows */

    Data0 = pv->Data; /* Record starting point */
    nrows = minvalue(nrows,(pv->lastRow)-firstrow+1);

    Setup_DataArrays( lParse, nCols, colData, firstrow, nrows );

    /* Parser allocates arrays for each column and calculation it performs. */
    /* Limit number of rows processed during each pass to reduce memory     */
    /* requirements... In most cases, iterator will limit rows to less      */
    /* than 10000 rows per iteration, so this is really only relevant for    */
    /* hk-compressed files which must be decompressed in memory and sent    */
    /* whole to fits_parser_workfn in a single iteration.                           */

    remain = nrows;
    while( remain ) {
       ntodo = minvalue(remain,10000);
       Evaluate_Parser ( lParse, firstrow, ntodo );
       if( lParse->status ) break;

       firstrow += ntodo;
       remain   -= ntodo;

       /*  Copy results into data array  */

       result = lParse->Nodes + lParse->resultNode;
       if( result->operation==CONST_OP ) constant = 1;

       switch( result->type ) {

       case BOOLEAN:
       case LONG:
       case DOUBLE:
          if( constant ) {
             char undef=0;
             for( kk=0; kk<ntodo; kk++ )
                for( jj=0; jj<(pv->repeat); jj++ )
                   ffcvtn( lParse->datatype,
                           &(result->value.data),
                           &undef, result->value.nelem /* 1 */,
                           (pv->userInfo)->datatype, (pv->Null),
                           (char*)(pv->Data) + (kk*(pv->repeat)+jj)*(pv->datasize),
                           &anyNullThisTime, &lParse->status );
          } else {
             if ( (pv->repeat) == result->value.nelem ) {
                ffcvtn( lParse->datatype,
                        result->value.data.ptr,
                        result->value.undef,
                        result->value.nelem*ntodo,
                        (pv->userInfo)->datatype, (pv->Null), (pv->Data),
                        &anyNullThisTime, &lParse->status );
             } else if( result->value.nelem == 1 ) {
                for( kk=0; kk<ntodo; kk++ )
                   for( jj=0; jj<(pv->repeat); jj++ ) {
                      ffcvtn( lParse->datatype,
                              (char*)result->value.data.ptr + kk*(pv->resDataSize),
                              (char*)result->value.undef + kk,
                              1, (pv->userInfo)->datatype, (pv->Null),
                              (char*)(pv->Data) + (kk*(pv->repeat)+jj)*(pv->datasize),
                              &anyNullThisTime, &lParse->status );
                   }
             } else {
                int nCopy;
                nCopy = minvalue( (pv->repeat), result->value.nelem );
                for( kk=0; kk<ntodo; kk++ ) {
                   ffcvtn( lParse->datatype,
                           (char*)result->value.data.ptr
                                  + kk*result->value.nelem*(pv->resDataSize),
                           (char*)result->value.undef
                                  + kk*result->value.nelem,
                           nCopy, (pv->userInfo)->datatype, (pv->Null),
                           (char*)(pv->Data) + (kk*(pv->repeat))*(pv->datasize),
                           &anyNullThisTime, &lParse->status );
                   if( nCopy < (pv->repeat) ) {
                      memset( (char*)(pv->Data) + (kk*(pv->repeat)+nCopy)*(pv->datasize),
                              0, ((pv->repeat)-nCopy)*(pv->datasize));
                   }
                }

             }
             if( result->operation>0 ) {
                FREE( result->value.data.ptr );
             }
          }
          if( lParse->status==OVERFLOW_ERR ) {
             lParse->status = NUM_OVERFLOW;
             ffpmsg("Numerical overflow while converting expression to necessary datatype");
          }
          break;

       case BITSTR:
          switch( (pv->userInfo)->datatype ) {
          case TBYTE:
             idx = -1;
             for( kk=0; kk<ntodo; kk++ ) {
                for( jj=0; jj<result->value.nelem; jj++ ) {
                   if( jj%8 == 0 )
                      ((char*)(pv->Data))[++idx] = 0;
                   if( constant ) {
                      if( result->value.data.str[jj]=='1' )
                         ((char*)(pv->Data))[idx] |= 128>>(jj%8);
                   } else {
                      if( result->value.data.strptr[kk][jj]=='1' )
                         ((char*)(pv->Data))[idx] |= 128>>(jj%8);
                   }
                }
             }
             break;
          case TBIT:
          case TLOGICAL:
             if( constant ) {
                for( kk=0; kk<ntodo; kk++ )
                   for( jj=0; jj<result->value.nelem; jj++ ) {
                      ((char*)(pv->Data))[ jj+kk*result->value.nelem ] =
                         ( result->value.data.str[jj]=='1' );
                   }
             } else {
                for( kk=0; kk<ntodo; kk++ )
                   for( jj=0; jj<result->value.nelem; jj++ ) {
                      ((char*)(pv->Data))[ jj+kk*result->value.nelem ] =
                         ( result->value.data.strptr[kk][jj]=='1' );
                   }
             }
             break; 
          case TSTRING:
             if( constant ) {
                for( jj=0; jj<ntodo; jj++ ) {
                   strcpy( ((char**)(pv->Data))[jj], result->value.data.str );
                }
             } else {
                for( jj=0; jj<ntodo; jj++ ) {
                   strcpy( ((char**)(pv->Data))[jj], result->value.data.strptr[jj] );
                }
             }
             break;
          default:
             ffpmsg("Cannot convert bit expression to desired type.");
             lParse->status = PARSE_BAD_TYPE;
             break;
          }
          if( result->operation>0 ) {
             FREE( result->value.data.strptr[0] );
             FREE( result->value.data.strptr );
          }
          break;

       case STRING:
          if( (pv->userInfo)->datatype==TSTRING ) {
             if( constant ) {
                for( jj=0; jj<ntodo; jj++ )
                   strcpy( ((char**)(pv->Data))[jj], result->value.data.str );
             } else {
                for( jj=0; jj<ntodo; jj++ )
                   if( result->value.undef[jj] ) {
                      anyNullThisTime = 1;
                      strcpy( ((char**)(pv->Data))[jj],
                              *(char **)(pv->Null) );
                   } else {
                      strcpy( ((char**)(pv->Data))[jj],
                              result->value.data.strptr[jj] );
                   }
             }
          } else {
             ffpmsg("Cannot convert string expression to desired type.");
             lParse->status = PARSE_BAD_TYPE;
          }
          if( result->operation>0 ) {
             FREE( result->value.data.strptr[0] );
             FREE( result->value.data.strptr );
          }
          break;
       }

       if( lParse->status ) break;

       /*  Increment Data to point to where the next block should go  */

       if( result->type==BITSTR && (pv->userInfo)->datatype==TBYTE )
          (pv->Data) = (char*)(pv->Data)
                    + (pv->datasize) * ( (result->value.nelem+7)/8 ) * ntodo;
       else if( result->type==STRING )
          (pv->Data) = (char*)(pv->Data) + (pv->datasize) * ntodo;
       else
          (pv->Data) = (char*)(pv->Data) + (pv->datasize) * ntodo * (pv->repeat);
    }

    /* If a TemporaryCol output is used, we want to inform the caller
       what the null value is expected to be */
    if (pv->Null != outcol->array && 
	(Data0) == (char*) outcol->array + (pv->datasize)) {
      if( (pv->userInfo)->datatype == TSTRING )
	memcpy( outcol->array, *(char **)(pv->Null), 2 );
      else 
	memcpy( outcol->array, (pv->Null), (pv->datasize) );
    }

    /* If no NULLs encountered during this pass, set Null value to */
    /* zero to make the writing of the output column data faster   */

    if( anyNullThisTime )
       (pv->userInfo)->anyNull = 1;
    else if( pv->Null == outcol->array ) {
       if( (pv->userInfo)->datatype == TSTRING )
          memcpy( *(char **)(pv->Null), zeros, 2 );
       else 
          memcpy( (pv->Null), zeros, (pv->datasize) );
    }

    /*-------------------------------------------------------*/
    /*  Clean up procedures:  after processing all the rows  */
    /*-------------------------------------------------------*/

    /*  if the calling routine specified that only a limited number    */
    /*  of rows in the table should be processed, return a value of -1 */
    /*  once all the rows have been done, if no other error occurred.  */

    if (lParse->hdutype != IMAGE_HDU && firstrow - 1 == (pv->lastRow)) {
           if (!lParse->status && (pv->userInfo)->maxRows<totalrows) {
                  return (-1);
           }
    }

    return(lParse->status);  /* return successful status */
}

static void Setup_DataArrays( ParseData *lParse, int nCols, iteratorCol *cols,
                              long fRow, long nRows )
    /***********************************************************************/
    /*  Setup the varData array in gParse to contain the fits column data. */
    /*  Then, allocate and initialize the necessary UNDEF arrays for each  */
    /*  column used by the parser.                                         */
    /***********************************************************************/
{
   int     i;
   long    nelem, len, row, idx;
   char  **bitStrs;
   char  **sptr;
   char   *barray;
   long   *iarray;
   double *rarray;
   char msg[80];
   int do_realloc = 0;

   lParse->firstDataRow = fRow;
   lParse->nDataRows    = nRows;
   /* Only perform reallocations if the number of rows changed */
   if (lParse->nPrevDataRows != nRows) do_realloc = 1;

   /*  Resize and fill in UNDEF arrays for each column  */

   for( i=0; i<nCols; i++ ) {

      iteratorCol *icol = cols + i;
      DataInfo *varData = lParse->varData + i;

      if( icol->iotype == OutputCol || icol->iotype == TemporaryCol ) continue;

      nelem  = varData->nelem;
      len    = nelem * nRows;

      switch ( varData->type ) {

      case BITSTR:
      /* No need for UNDEF array, but must make string DATA array */
         len = (nelem+1)*nRows;   /* Count '\0' */
         bitStrs = (char**)varData->data;
	 if (do_realloc) {
	   if( bitStrs ) FREE( bitStrs[0] );
	   free( bitStrs );
	   bitStrs = (char**)malloc( nRows*sizeof(char*) );
	   if( bitStrs==NULL ) {
	     varData->data = varData->undef = NULL;
	     lParse->status = MEMORY_ALLOCATION;
	     break;
	   }
	   bitStrs[0] = (char*)malloc( len*sizeof(char) );
	   if( bitStrs[0]==NULL ) {
	     free( bitStrs );
	     varData->data = varData->undef = NULL;
	     lParse->status = MEMORY_ALLOCATION;
	     break;
	   }
	 }

         for( row=0; row<nRows; row++ ) {
            bitStrs[row] = bitStrs[0] + row*(nelem+1);
            idx = (row)*( (nelem+7)/8 ) + 1;
            for(len=0; len<nelem; len++) {
               if( ((char*)icol->array)[idx] & (1<<(7-len%8)) )
                  bitStrs[row][len] = '1';
               else
                  bitStrs[row][len] = '0';
               if( len%8==7 ) idx++;
            }
            bitStrs[row][len] = '\0';
         }
         varData->undef = (char*)bitStrs;
         varData->data  = (char*)bitStrs;
         break;

      case STRING:
         sptr = (char**)icol->array;
	 if (do_realloc) {
	   if (varData->undef)
	     free( varData->undef );
	   varData->undef = (char*)malloc( nRows*sizeof(char) );
	   if( varData->undef==NULL ) {
	     lParse->status = MEMORY_ALLOCATION;
	     break;
	   }
	 }
         row = nRows;
         while( row-- )
            varData->undef[row] =
               ( **sptr != '\0' && FSTRCMP( sptr[0], sptr[row+1] )==0 );
         varData->data  = sptr + 1;
         break;

      case BOOLEAN:
         barray = (char*)icol->array;
	 if (do_realloc) {
	   if (varData->undef)
	     free( varData->undef );
	   varData->undef = (char*)malloc( len*sizeof(char) );
	   if( varData->undef==NULL ) {
	     lParse->status = MEMORY_ALLOCATION;
	     break;
	   }
	 }
         while( len-- ) {
            varData->undef[len] = 
               ( barray[0]!=0 && barray[0]==barray[len+1] );
         }
         varData->data  = barray + 1;
         break;

      case LONG:
         iarray = (long*)icol->array;
	 if (do_realloc) {
	   if (varData->undef)
	     free( varData->undef );
	   varData->undef = (char*)malloc( len*sizeof(char) );
	   if( varData->undef==NULL ) {
	     lParse->status = MEMORY_ALLOCATION;
	     break;
	   }
	 }
         while( len-- ) {
            varData->undef[len] = 
               ( iarray[0]!=0L && iarray[0]==iarray[len+1] );
         }
         varData->data  = iarray + 1;
         break;

      case DOUBLE:
         rarray = (double*)icol->array;
	 if (do_realloc) {
	   if (varData->undef)
	     free( varData->undef );
	   varData->undef = (char*)malloc( len*sizeof(char) );
	   if( varData->undef==NULL ) {
	     lParse->status = MEMORY_ALLOCATION;
	     break;
	   }
	 }
         while( len-- ) {
            varData->undef[len] = 
               ( rarray[0]!=0.0 && rarray[0]==rarray[len+1]);
         }
         varData->data  = rarray + 1;
         break;

      default:
         snprintf(msg, 80, "SetupDataArrays, unhandled type %d\n",
                varData->type);
         ffpmsg(msg);
      }

      if( lParse->status ) {  /*  Deallocate NULL arrays of previous columns */
         while( i-- ) {
            varData = lParse->varData + i;
            if( varData->type==BITSTR )
               FREE( ((char**)varData->data)[0] );
            FREE( varData->undef );
            varData->undef = NULL;
         }
	 lParse->nPrevDataRows = 0;
         return;
      }
   }

   lParse->nPrevDataRows = nRows;
}

/*--------------------------------------------------------------------------*/
int ffcvtn( int   inputType,  /* I - Data type of input array               */
            void  *input,     /* I - Input array of type inputType          */
            char  *undef,     /* I - Array of flags indicating UNDEF elems  */
            long  ntodo,      /* I - Number of elements to process          */
            int   outputType, /* I - Data type of output array              */
            void  *nulval,    /* I - Ptr to value to use for UNDEF elements */
            void  *output,    /* O - Output array of type outputType        */
            int   *anynull,   /* O - Any nulls flagged?                     */
            int   *status )   /* O - Error status                           */
/*                                                                          */
/* Convert an array of any input data type to an array of any output        */
/* data type, using an array of UNDEF flags to assign nulvals to            */
/*--------------------------------------------------------------------------*/
{
   long i;

   switch( outputType ) {

   case TLOGICAL:
      switch( inputType ) {
      case TLOGICAL:
      case TBYTE:
         for( i=0; i<ntodo; i++ )
            if( ((unsigned char*)input)[i] )
                ((unsigned char*)output)[i] = 1;
            else
                ((unsigned char*)output)[i] = 0;
         break;
      case TSHORT:
         for( i=0; i<ntodo; i++ )
            if( ((short*)input)[i] )
                ((unsigned char*)output)[i] = 1;
            else
                ((unsigned char*)output)[i] = 0;
         break;
      case TLONG:
         for( i=0; i<ntodo; i++ )
            if( ((long*)input)[i] )
                ((unsigned char*)output)[i] = 1;
            else
                ((unsigned char*)output)[i] = 0;
         break;
      case TFLOAT:
         for( i=0; i<ntodo; i++ )
            if( ((float*)input)[i] )
                ((unsigned char*)output)[i] = 1;
            else
                ((unsigned char*)output)[i] = 0;
         break;
      case TDOUBLE:
         for( i=0; i<ntodo; i++ )
            if( ((double*)input)[i] )
                ((unsigned char*)output)[i] = 1;
            else
                ((unsigned char*)output)[i] = 0;
         break;
      default:
         *status = BAD_DATATYPE;
         break;
      }
      for(i=0;i<ntodo;i++) {
         if( undef[i] ) {
            ((unsigned char*)output)[i] = *(unsigned char*)nulval;
            *anynull = 1;
         }
      }
      break;

   case TBYTE:
      switch( inputType ) {
      case TLOGICAL:
      case TBYTE:
         for( i=0; i<ntodo; i++ )
            ((unsigned char*)output)[i] = ((unsigned char*)input)[i];
         break;
      case TSHORT:
         fffi2i1((short*)input,ntodo,1.,0.,0,0,0,NULL,NULL,(unsigned char*)output,status);
         break;
      case TLONG:
         for (i = 0; i < ntodo; i++) {
            if( undef[i] ) {
               ((unsigned char*)output)[i] = *(unsigned char*)nulval;
               *anynull = 1;
            } else {
               if( ((long*)input)[i] < 0 ) {
                  *status = OVERFLOW_ERR;
                  ((unsigned char*)output)[i] = 0;
               } else if( ((long*)input)[i] > UCHAR_MAX ) {
                  *status = OVERFLOW_ERR;
                  ((unsigned char*)output)[i] = UCHAR_MAX;
               } else
                  ((unsigned char*)output)[i] = 
                     (unsigned char) ((long*)input)[i];
            }
         }
         return( *status );
      case TFLOAT:
         fffr4i1((float*)input,ntodo,1.,0.,0,0,NULL,NULL,
                 (unsigned char*)output,status);
         break;
      case TDOUBLE:
         fffr8i1((double*)input,ntodo,1.,0.,0,0,NULL,NULL,
                 (unsigned char*)output,status);
         break;
      default:
         *status = BAD_DATATYPE;
         break;
      }
      for(i=0;i<ntodo;i++) {
         if( undef[i] ) {
            ((unsigned char*)output)[i] = *(unsigned char*)nulval;
            *anynull = 1;
         }
      }
      break;

   case TSHORT:
      switch( inputType ) {
      case TLOGICAL:
      case TBYTE:
         for( i=0; i<ntodo; i++ )
            ((short*)output)[i] = ((unsigned char*)input)[i];
         break;
      case TSHORT:
         for( i=0; i<ntodo; i++ )
            ((short*)output)[i] = ((short*)input)[i];
         break;
      case TLONG:
         for (i = 0; i < ntodo; i++) {
            if( undef[i] ) {
               ((short*)output)[i] = *(short*)nulval;
               *anynull = 1;
            } else {
               if( ((long*)input)[i] < SHRT_MIN ) {
                  *status = OVERFLOW_ERR;
                  ((short*)output)[i] = SHRT_MIN;
               } else if ( ((long*)input)[i] > SHRT_MAX ) {
                  *status = OVERFLOW_ERR;
                  ((short*)output)[i] = SHRT_MAX;
               } else
                  ((short*)output)[i] = (short) ((long*)input)[i];
            }
         }
         return( *status );
      case TFLOAT:
         fffr4i2((float*)input,ntodo,1.,0.,0,0,NULL,NULL,
                 (short*)output,status);
         break;
      case TDOUBLE:
         fffr8i2((double*)input,ntodo,1.,0.,0,0,NULL,NULL,
                 (short*)output,status);
         break;
      default:
         *status = BAD_DATATYPE;
         break;
      }
      for(i=0;i<ntodo;i++) {
         if( undef[i] ) {
            ((short*)output)[i] = *(short*)nulval;
            *anynull = 1;
         }
      }
      break;

   case TINT:
      switch( inputType ) {
      case TLOGICAL:
      case TBYTE:
         for( i=0; i<ntodo; i++ )
            ((int*)output)[i] = ((unsigned char*)input)[i];
         break;
      case TSHORT:
         for( i=0; i<ntodo; i++ )
            ((int*)output)[i] = ((short*)input)[i];
         break;
      case TLONG:
         for( i=0; i<ntodo; i++ )
            ((int*)output)[i] = ((long*)input)[i];
         break;
      case TFLOAT:
         fffr4int((float*)input,ntodo,1.,0.,0,0,NULL,NULL,
                  (int*)output,status);
         break;
      case TDOUBLE:
         fffr8int((double*)input,ntodo,1.,0.,0,0,NULL,NULL,
                  (int*)output,status);
         break;
      default:
         *status = BAD_DATATYPE;
         break;
      }
      for(i=0;i<ntodo;i++) {
         if( undef[i] ) {
            ((int*)output)[i] = *(int*)nulval;
            *anynull = 1;
         }
      }
      break;

   case TLONG:
      switch( inputType ) {
      case TLOGICAL:
      case TBYTE:
         for( i=0; i<ntodo; i++ )
            ((long*)output)[i] = ((unsigned char*)input)[i];
         break;
      case TSHORT:
         for( i=0; i<ntodo; i++ )
            ((long*)output)[i] = ((short*)input)[i];
         break;
      case TLONG:
         for( i=0; i<ntodo; i++ )
            ((long*)output)[i] = ((long*)input)[i];
         break;
      case TFLOAT:
         fffr4i4((float*)input,ntodo,1.,0.,0,0,NULL,NULL,
                 (long*)output,status);
         break;
      case TDOUBLE:
         fffr8i4((double*)input,ntodo,1.,0.,0,0,NULL,NULL,
                 (long*)output,status);
         break;
      default:
         *status = BAD_DATATYPE;
         break;
      }
      for(i=0;i<ntodo;i++) {
         if( undef[i] ) {
            ((long*)output)[i] = *(long*)nulval;
            *anynull = 1;
         }
      }
      break;

   case TLONGLONG:
      switch( inputType ) {
      case TLOGICAL:
      case TBYTE:
         for( i=0; i<ntodo; i++ )
            ((LONGLONG*)output)[i] = ((unsigned char*)input)[i];
         break;
      case TSHORT:
         for( i=0; i<ntodo; i++ )
            ((LONGLONG*)output)[i] = ((short*)input)[i];
         break;
      case TLONG:
         for( i=0; i<ntodo; i++ )
            ((LONGLONG*)output)[i] = ((long*)input)[i];
         break;
      case TFLOAT:
         fffr4i8((float*)input,ntodo,1.,0.,0,0,NULL,NULL,
                 (LONGLONG*)output,status);
         break;
      case TDOUBLE:
         fffr8i8((double*)input,ntodo,1.,0.,0,0,NULL,NULL,
                 (LONGLONG*)output,status);

         break;
      default:
         *status = BAD_DATATYPE;
         break;
      }
      for(i=0;i<ntodo;i++) {
         if( undef[i] ) {
            ((LONGLONG*)output)[i] = *(LONGLONG*)nulval;
            *anynull = 1;
         }
      }
      break;

   case TFLOAT:
      switch( inputType ) {
      case TLOGICAL:
      case TBYTE:
         for( i=0; i<ntodo; i++ )
            ((float*)output)[i] = ((unsigned char*)input)[i];
         break;
      case TSHORT:
         for( i=0; i<ntodo; i++ )
            ((float*)output)[i] = ((short*)input)[i];
         break;
      case TLONG:
         for( i=0; i<ntodo; i++ )
            ((float*)output)[i] = (float) ((long*)input)[i];
         break;
      case TFLOAT:
         for( i=0; i<ntodo; i++ )
            ((float*)output)[i] = ((float*)input)[i];
         break;
      case TDOUBLE:
         fffr8r4((double*)input,ntodo,1.,0.,0,0,NULL,NULL,
                 (float*)output,status);
         break;
      default:
         *status = BAD_DATATYPE;
         break;
      }
      for(i=0;i<ntodo;i++) {
         if( undef[i] ) {
            ((float*)output)[i] = *(float*)nulval;
            *anynull = 1;
         }
      }
      break;

   case TDOUBLE:
      switch( inputType ) {
      case TLOGICAL:
      case TBYTE:
         for( i=0; i<ntodo; i++ )
            ((double*)output)[i] = ((unsigned char*)input)[i];
         break;
      case TSHORT:
         for( i=0; i<ntodo; i++ )
            ((double*)output)[i] = ((short*)input)[i];
         break;
      case TLONG:
         for( i=0; i<ntodo; i++ )
            ((double*)output)[i] = ((long*)input)[i];
         break;
      case TFLOAT:
         for( i=0; i<ntodo; i++ )
            ((double*)output)[i] = ((float*)input)[i];
         break;
      case TDOUBLE:
         for( i=0; i<ntodo; i++ )
            ((double*)output)[i] = ((double*)input)[i];
         break;
      default:
         *status = BAD_DATATYPE;
         break;
      }
      for(i=0;i<ntodo;i++) {
         if( undef[i] ) {
            ((double*)output)[i] = *(double*)nulval;
            *anynull = 1;
         }
      }
      break;

   default:
      *status = BAD_DATATYPE;
      break;
   }

   return ( *status );
}

/*---------------------------------------------------------------------------*/
int fffrwc( fitsfile *fptr,        /* I - Input FITS file                    */
            char     *expr,        /* I - Boolean expression                 */
            char     *timeCol,     /* I - Name of time column                */
            char     *parCol,      /* I - Name of parameter column           */
            char     *valCol,      /* I - Name of value column               */
            long     ntimes,       /* I - Number of distinct times in file   */
            double   *times,       /* O - Array of times in file             */
            char     *time_status, /* O - Array of boolean results           */
            int      *status )     /* O - Error status                       */
/*                                                                           */
/* Evaluate a boolean expression for each time in a compressed file,         */
/* returning an array of flags indicating which times evaluated to TRUE/FALSE*/
/*---------------------------------------------------------------------------*/
{
   parseInfo Info;
   long alen, width;
   int parNo, typecode;
   int naxis, constant, nCol=0;
   long nelem, naxes[MAXDIMS], elem;
   char result;
   ParseData lParse;

   if( *status ) return( *status );

   memset(&Info, 0, sizeof(Info));   

   if( ffiprs( fptr, 1, expr, MAXDIMS, &Info.datatype, &nelem,
               &naxis, naxes, &lParse, status ) ) {
      ffcprs(&lParse);
      return( *status );
   }

   fits_get_colnum( fptr, CASEINSEN, timeCol, &lParse.timeCol, status );
   fits_get_colnum( fptr, CASEINSEN, parCol,  &lParse.parCol , status );
   fits_get_colnum( fptr, CASEINSEN, valCol,  &lParse.valCol, status );
   if( *status ) return( *status );
   
   if( nelem<0 ) {
      constant = 1;
      nelem = -nelem;
      nCol = lParse.nCols;
      lParse.nCols = 0;    /*  Ignore all column references  */
   } else
      constant = 0;

   if( Info.datatype!=TLOGICAL || nelem!=1 ) {
      ffcprs(&lParse);
      ffpmsg("Expression does not evaluate to a logical scalar.");
      return( *status = PARSE_BAD_TYPE );
   }

   /*******************************************/
   /* Allocate data arrays for each parameter */
   /*******************************************/
   
   parNo = lParse.nCols;
   while( parNo-- ) {
      switch( lParse.colData[parNo].datatype ) {
      case TLONG:
         if( (lParse.colData[parNo].array =
              (long *)malloc( (ntimes+1)*sizeof(long) )) )
            ((long*)lParse.colData[parNo].array)[0] = 1234554321;
         else
            *status = MEMORY_ALLOCATION;
         break;
      case TDOUBLE:
         if( (lParse.colData[parNo].array =
              (double *)malloc( (ntimes+1)*sizeof(double) )) )
            ((double*)lParse.colData[parNo].array)[0] = DOUBLENULLVALUE;
         else
            *status = MEMORY_ALLOCATION;
         break;
      case TSTRING:
         if( !fits_get_coltype( fptr, lParse.valCol, &typecode,
                                &alen, &width, status ) ) {
            alen++;
            if( (lParse.colData[parNo].array =
                 (char **)malloc( (ntimes+1)*sizeof(char*) )) ) {
               if( (((char **)lParse.colData[parNo].array)[0] =
                    (char *)malloc( (ntimes+1)*sizeof(char)*alen )) ) {
                  for( elem=1; elem<=ntimes; elem++ )
                     ((char **)lParse.colData[parNo].array)[elem] =
                        ((char **)lParse.colData[parNo].array)[elem-1]+alen;
                  ((char **)lParse.colData[parNo].array)[0][0] = '\0';
               } else {
                  free( lParse.colData[parNo].array );
                  *status = MEMORY_ALLOCATION;
               }
            } else {
               *status = MEMORY_ALLOCATION;
            }
         }
         break;
      }
      if( *status ) {
         while( parNo-- ) {
            if( lParse.colData[parNo].datatype==TSTRING )
               FREE( ((char **)lParse.colData[parNo].array)[0] );
            FREE( lParse.colData[parNo].array );
         }
         return( *status );
      }
   }
   
   /**********************************************************************/
   /* Read data from columns needed for the expression and then parse it */
   /**********************************************************************/
   
   if( !fits_uncompress_hkdata( &lParse, fptr, ntimes, times, status ) ) {
      if( constant ) {
         result = lParse.Nodes[lParse.resultNode].value.data.log;
         elem = ntimes;
         while( elem-- ) time_status[elem] = result;
      } else {
         Info.dataPtr  = time_status;
         Info.nullPtr  = NULL;
         Info.maxRows  = ntimes;
         *status       = fits_parser_workfn( ntimes, 0, 1, ntimes, lParse.nCols,
                                     lParse.colData, (void*)&Info );
      }
   }
   
   /************/
   /* Clean up */
   /************/
   
   parNo = lParse.nCols;
   while ( parNo-- ) {
      if( lParse.colData[parNo].datatype==TSTRING )
         FREE( ((char **)lParse.colData[parNo].array)[0] );
      FREE( lParse.colData[parNo].array );
   }
   
   if( constant ) lParse.nCols = nCol;

   ffcprs(&lParse);
   return(*status);
}

/*---------------------------------------------------------------------------*/
int fits_parser_set_temporary_col(ParseData *lParse,
				  parseInfo *Info,
				  long int nrows,
				  void *nulval,
				  int *status)
{
  int col_cnt;
  /* Setup iterator column and parser information to be ready to compute 
     temporary calculator expression */

  if (*status) return *status;

  col_cnt = lParse->nCols;

  if( fits_parser_allocateCol( lParse, col_cnt, status ) ) return *status;
	
  /* Set important variables for TemporaryCol where calculated results end up */
  fits_iter_set_by_num( &(lParse->colData[col_cnt]), 0, 0, TDOUBLE, TemporaryCol);
  lParse->colData[col_cnt].repeat = lParse->nElements;
  Info->dataPtr = NULL;
  Info->nullPtr = nulval;
  Info->maxRows = nrows;
  Info->parseData = lParse;
  lParse->nCols ++;

  return 0;
}

/*---------------------------------------------------------------------------*/
int fits_uncompress_hkdata( ParseData *lParse,
		       fitsfile *fptr,
                       long     ntimes,
                       double   *times,
                       int      *status )
/*                                                                           */
/* description                                                               */
/*---------------------------------------------------------------------------*/
{
   char parName[256], *sPtr[1], found[1000];
   int parNo, anynul;
   long naxis2, row, currelem;
   double currtime, newtime;

   sPtr[0] = parName;
   currelem = 0;
   currtime = -1e38;

   parNo=lParse->nCols;
   while( parNo-- ) found[parNo] = 0;

   if( ffgkyj( fptr, "NAXIS2", &naxis2, NULL, status ) ) return( *status );

   for( row=1; row<=naxis2; row++ ) {
      if( ffgcvd( fptr, lParse->timeCol, row, 1L, 1L, 0.0,
                  &newtime, &anynul, status ) ) return( *status );
      if( newtime != currtime ) {
         /*  New time encountered... propogate parameters to next row  */
         if( currelem==ntimes ) {
            ffpmsg("Found more unique time stamps than caller indicated");
            return( *status = PARSE_BAD_COL );
         }
         times[currelem++] = currtime = newtime;
         parNo = lParse->nCols;
         while( parNo-- ) {
            switch( lParse->colData[parNo].datatype ) {
            case TLONG:
               ((long*)lParse->colData[parNo].array)[currelem] =
                  ((long*)lParse->colData[parNo].array)[currelem-1];
               break;
            case TDOUBLE:
               ((double*)lParse->colData[parNo].array)[currelem] =
                  ((double*)lParse->colData[parNo].array)[currelem-1];
               break;
            case TSTRING:
               strcpy( ((char **)lParse->colData[parNo].array)[currelem],
                       ((char **)lParse->colData[parNo].array)[currelem-1] );
               break;
            }
         }
      }

      if( ffgcvs( fptr, lParse->parCol, row, 1L, 1L, "",
                  sPtr, &anynul, status ) ) return( *status );
      parNo = lParse->nCols;
      while( parNo-- )
         if( !fits_strcasecmp( parName, lParse->varData[parNo].name ) ) break;

      if( parNo>=0 ) {
         found[parNo] = 1; /* Flag this parameter as found */
         switch( lParse->colData[parNo].datatype ) {
         case TLONG:
            ffgcvj( fptr, lParse->valCol, row, 1L, 1L,
                    ((long*)lParse->colData[parNo].array)[0],
                    ((long*)lParse->colData[parNo].array)+currelem,
                    &anynul, status );
            break;
         case TDOUBLE:
            ffgcvd( fptr, lParse->valCol, row, 1L, 1L,
                    ((double*)lParse->colData[parNo].array)[0],
                    ((double*)lParse->colData[parNo].array)+currelem,
                    &anynul, status );
            break;
         case TSTRING:
            ffgcvs( fptr, lParse->valCol, row, 1L, 1L,
                    ((char**)lParse->colData[parNo].array)[0],
                    ((char**)lParse->colData[parNo].array)+currelem,
                    &anynul, status );
            break;
         }
         if( *status ) return( *status );
      }
   }

   if( currelem<ntimes ) {
      ffpmsg("Found fewer unique time stamps than caller indicated");
      return( *status = PARSE_BAD_COL );
   }

   /*  Check for any parameters which were not located in the table  */
   parNo = lParse->nCols;
   while( parNo-- )
      if( !found[parNo] ) {
         snprintf( parName, 256, "Parameter not found: %-30s", 
                  lParse->varData[parNo].name );
         ffpmsg( parName );
         *status = PARSE_SYNTAX_ERR;
      }
   return( *status );
}

typedef struct {
  long *prownum;
  ParseData *lParse;
} ffffrw_workdata;

/*---------------------------------------------------------------------------*/
int ffffrw( fitsfile *fptr,         /* I - Input FITS file                   */
            char     *expr,         /* I - Boolean expression                */
            long     *rownum,       /* O - First row of table to eval to T   */
            int      *status )      /* O - Error status                      */
/*                                                                           */
/* Evaluate a boolean expression, returning the row number of the first      */
/* row which evaluates to TRUE                                               */
/*---------------------------------------------------------------------------*/
{
   int naxis, constant, dtype;
   long nelem, naxes[MAXDIMS];
   char result;
   ParseData lParse;

   if( *status ) return( *status );

   if( ffiprs( fptr, 0, expr, MAXDIMS, &dtype, &nelem, &naxis,
               naxes, &lParse, status ) ) {
      ffcprs(&lParse);
      return( *status );
   }
   if( nelem<0 ) {
      constant = 1;
      nelem = -nelem;
   } else
      constant = 0;

   if( dtype!=TLOGICAL || nelem!=1 ) {
      ffcprs(&lParse);
      ffpmsg("Expression does not evaluate to a logical scalar.");
      return( *status = PARSE_BAD_TYPE );
   }

   *rownum = 0;
   if( constant ) { /* No need to call parser... have result from ffiprs */
      result = lParse.Nodes[lParse.resultNode].value.data.log;
      if( result ) {
         /*  Make sure there is at least 1 row in table  */
         ffgnrw( fptr, &nelem, status );
         if( nelem )
            *rownum = 1;
      }
   } else {
      ffffrw_workdata workData;
      workData.prownum = rownum;
      workData.lParse = &lParse;
      if( ffiter( lParse.nCols, lParse.colData, 0, 0,
                  ffffrw_work, (void*)&workData, status ) == -1 )
         *status = 0;  /* -1 indicates exitted without error before end... OK */
   }

   ffcprs(&lParse);
   return(*status);
}

/*---------------------------------------------------------------------------*/
int ffffrw_work(long        totalrows, /* I - Total rows to be processed     */
                long        offset,    /* I - Number of rows skipped at start*/
                long        firstrow,  /* I - First row of this iteration    */
                long        nrows,     /* I - Number of rows in this iter    */
                int         nCols,     /* I - Number of columns in use       */
                iteratorCol *colData,  /* IO- Column information/data        */
                void        *userPtr ) /* I - Data handling instructions     */
/*                                                                           */
/* Iterator work function which calls the parser and searches for the        */
/* first row which evaluates to TRUE.                                        */
/*---------------------------------------------------------------------------*/
{
    long idx;
    Node *result;
    ffffrw_workdata *workData = userPtr;
    ParseData *lParse = workData->lParse;

    Evaluate_Parser( lParse, firstrow, nrows );

    if( !lParse->status ) {

       result = lParse->Nodes + lParse->resultNode;
       if( result->operation==CONST_OP ) {

          if( result->value.data.log ) {
	     *(workData->prownum) = firstrow;
             return( -1 );
          }

       } else {

          for( idx=0; idx<nrows; idx++ )
             if( result->value.data.logptr[idx] && !result->value.undef[idx] ) {
	        *(workData->prownum) = firstrow + idx;;
                return( -1 );
             }
       }
    }

    return( lParse->status );
}


static int set_image_col_types (ParseData *lParse,
				fitsfile * fptr, const char * name, int bitpix,
				DataInfo * varInfo, iteratorCol *colIter) {

   int istatus;
   double tscale, tzero;
   char temp[80];

   switch (bitpix) {
      case BYTE_IMG:
      case SHORT_IMG:
      case LONG_IMG:
         istatus = 0;
         if (fits_read_key(fptr, TDOUBLE, "BZERO", &tzero, NULL, &istatus))
            tzero = 0.0;

         istatus = 0;
         if (fits_read_key(fptr, TDOUBLE, "BSCALE", &tscale, NULL, &istatus))
            tscale = 1.0;

         if (tscale == 1.0 && (tzero == 0.0 || tzero == 32768.0 )) {
            varInfo->type     = LONG;
            colIter->datatype = TLONG;
         }
         else {
            varInfo->type     = DOUBLE;
            colIter->datatype = TDOUBLE;
            if (DEBUG_PIXFILTER)
                printf("use DOUBLE for %s with BSCALE=%g/BZERO=%g\n",
                        name, tscale, tzero);
         }
         break;

      case LONGLONG_IMG:
      case FLOAT_IMG:
      case DOUBLE_IMG:
         varInfo->type     = DOUBLE;
         colIter->datatype = TDOUBLE;
         break;
      default:
         snprintf(temp, 80,"set_image_col_types: unrecognized image bitpix [%d]\n",
                bitpix);
         ffpmsg(temp);
         return lParse->status = PARSE_BAD_TYPE;
   }
   return 0;
}


/*************************************************************************

        Functions used by the evaluator to access FITS data
            (find_column, find_keywd, fits_parser_allocateCol, load_column)

 *************************************************************************/

static int find_column( ParseData *lParse, char *colName, void *itslval )
{
   FITS_PARSER_YYSTYPE *thelval = (FITS_PARSER_YYSTYPE*)itslval;
   int col_cnt, status;
   int colnum, typecode, type;
   long repeat, width;
   fitsfile *fptr;
   char temp[80];
   double tzero,tscale;
   int istatus;
   DataInfo *varInfo;
   iteratorCol *colIter;

if (DEBUG_PIXFILTER)
   printf("find_column(%s)\n", colName);

   if( *colName == '#' )
     return( find_keywd( lParse, colName + 1, itslval ) );

   fptr = lParse->def_fptr;

   status = 0;
   col_cnt = lParse->nCols;

if (lParse->hdutype == IMAGE_HDU) {
   int i;
   if (!lParse->pixFilter) {
      lParse->status = COL_NOT_FOUND;
      ffpmsg("find_column: IMAGE_HDU but no PixelFilter");
      return pERROR;
   }

   colnum = -1;
   for (i = 0; i < lParse->pixFilter->count; ++i) {
      if (!fits_strcasecmp(colName, lParse->pixFilter->tag[i]))
         colnum = i;
   }
   if (colnum < 0) {
      snprintf(temp, 80, "find_column: PixelFilter tag %s not found", colName);
      ffpmsg(temp);
      lParse->status = COL_NOT_FOUND;
      return pERROR;
   }

   if( fits_parser_allocateCol( lParse, col_cnt, &lParse->status ) ) return pERROR;

   varInfo = lParse->varData + col_cnt;
   colIter = lParse->colData + col_cnt;

   fptr = lParse->pixFilter->ifptr[colnum];
   fits_get_img_param(fptr,
                MAXDIMS,
                &typecode, /* actually bitpix */
                &varInfo->naxis,
                &varInfo->naxes[0],
                &status);
   varInfo->nelem = 1;
   type = COLUMN;
   if (set_image_col_types(lParse, fptr, colName, typecode, varInfo, colIter))
      return pERROR;
   colIter->fptr = fptr;
   colIter->iotype = InputCol;
}
else { /* HDU holds a table */
   if( lParse->compressed )
      colnum = lParse->valCol;
   else
      if( fits_get_colnum( fptr, CASEINSEN, colName, &colnum, &status ) ) {
         if( status == COL_NOT_FOUND ) {
	   type = find_keywd( lParse, colName, itslval );
            if( type != pERROR ) ffcmsg();
            return( type );
         }
         lParse->status = status;
         return pERROR;
      }
   
   if( fits_get_coltype( fptr, colnum, &typecode,
                         &repeat, &width, &status ) ) {
      lParse->status = status;
      return pERROR;
   }

   if( fits_parser_allocateCol( lParse, col_cnt, &lParse->status ) ) return pERROR;

   varInfo = lParse->varData + col_cnt;
   colIter = lParse->colData + col_cnt;

   fits_iter_set_by_num( colIter, fptr, colnum, 0, InputCol );
}

   /*  Make sure we don't overflow variable name array  */
   strncpy(varInfo->name,colName,MAXVARNAME);
   varInfo->name[MAXVARNAME] = '\0';

if (lParse->hdutype != IMAGE_HDU) {
   switch( typecode ) {
   case TBIT:
      varInfo->type     = BITSTR;
      colIter->datatype = TBYTE;
      type = BITCOL;
      break;
   case TBYTE:
   case TSHORT:
   case TLONG:
      /* The datatype of column with TZERO and TSCALE keywords might be 
         float or double. 
      */
      snprintf(temp,80,"TZERO%d",colnum);
      istatus = 0;
      if(fits_read_key(fptr,TDOUBLE,temp,&tzero,NULL,&istatus)) {
          tzero = 0.0;
      } 
      snprintf(temp,80,"TSCAL%d",colnum);
      istatus = 0;
      if(fits_read_key(fptr,TDOUBLE,temp,&tscale,NULL,&istatus)) {
          tscale = 1.0;
      } 
      if (tscale == 1.0 && (tzero == 0.0 || tzero == 32768.0 )) {
          varInfo->type     = LONG;
          colIter->datatype = TLONG;
/*    Reading an unsigned long column as a long can cause overflow errors.
      Treat the column as a double instead.
      } else if (tscale == 1.0 &&  tzero == 2147483648.0 ) {
          varInfo->type     = LONG;
          colIter->datatype = TULONG;
 */

      }
      else {
          varInfo->type     = DOUBLE;
          colIter->datatype = TDOUBLE;
      }
      type = COLUMN;
      break;
/* 
  For now, treat 8-byte integer columns as type double.
  This can lose precision, so the better long term solution
  will be to add support for TLONGLONG as a separate datatype.
*/
   case TLONGLONG:
   case TFLOAT:
   case TDOUBLE:
      varInfo->type     = DOUBLE;
      colIter->datatype = TDOUBLE;
      type = COLUMN;
      break;
   case TLOGICAL:
      varInfo->type     = BOOLEAN;
      colIter->datatype = TLOGICAL;
      type = BCOLUMN;
      break;
   case TSTRING:
      varInfo->type     = STRING;
      colIter->datatype = TSTRING;
      type = SCOLUMN;
      if ( width >= MAX_STRLEN ) {
	snprintf(temp, 80, "column %d is wider than maximum %d characters",
		colnum, MAX_STRLEN-1);
        ffpmsg(temp);
	lParse->status = PARSE_LRG_VECTOR;
	return pERROR;
      }
      if( lParse->hdutype == ASCII_TBL ) repeat = width;
      break;
   default:
      if (typecode < 0) {
        snprintf(temp, 80,"variable-length array columns are not supported. typecode = %d", typecode);
        ffpmsg(temp);
      }
      lParse->status = PARSE_BAD_TYPE;
      return pERROR;
   }
   varInfo->nelem = repeat;
   colIter->repeat = 0; /* ffiter() will fill in this value */
   if( repeat>1 && typecode!=TSTRING ) {
      if( fits_read_tdim( fptr, colnum, MAXDIMS,
                          &varInfo->naxis,
                          &varInfo->naxes[0], &status )
          ) {
         lParse->status = status;
         return pERROR;
      }
   } else {
      varInfo->naxis = 1;
      varInfo->naxes[0] = 1;
   }
}
   lParse->nCols++;
   thelval->lng = col_cnt;

   return( type );
}

static int find_keywd(ParseData *lParse, char *keyname, void *itslval )
{
   FITS_PARSER_YYSTYPE *thelval = (FITS_PARSER_YYSTYPE*)itslval;
   int status, type;
   char keyvalue[FLEN_VALUE], dtype;
   fitsfile *fptr;
   double rval;
   int bval;
   long ival;

   status = 0;
   fptr = lParse->def_fptr;
   if( fits_read_keyword( fptr, keyname, keyvalue, NULL, &status ) ) {
      if( status == KEY_NO_EXIST ) {
         /*  Do this since ffgkey doesn't put an error message on stack  */
         snprintf(keyvalue,FLEN_VALUE, "ffgkey could not find keyword: %s",keyname);
         ffpmsg(keyvalue);
      }
      lParse->status = status;
      return( pERROR );
   }
      
   if( fits_get_keytype( keyvalue, &dtype, &status ) ) {
      lParse->status = status;
      return( pERROR );
   }

   /* Read appropriate value type and set to CONST_OP */
   switch( dtype ) {
   case 'C':
      fits_read_key_str( fptr, keyname, keyvalue, NULL, &status );
      type = STRING;
      strcpy( thelval->str , keyvalue );
      break;
   case 'L':
      fits_read_key_log( fptr, keyname, &bval, NULL, &status );
      type = BOOLEAN;
      thelval->log = bval;
      break;
   case 'I':
      fits_read_key_lng( fptr, keyname, &ival, NULL, &status );
      type = LONG;
      thelval->lng = ival;
      break;
   case 'F':
      fits_read_key_dbl( fptr, keyname, &rval, NULL, &status );
      type = DOUBLE;
      thelval->dbl = rval;
      break;
   default:
      type = pERROR;
      break;
   }

   if( status ) {
      lParse->status=status;
      return pERROR;
   }

   return( type );
}

/* Allocates parser iterator column storage for 25 columns *more* than
   nCols */
int fits_parser_allocateCol( ParseData *lParse, int nCol, int *status )
{
   if( (nCol%25)==0 ) {
     lParse->colData = (iteratorCol*) fits_recalloc( lParse->colData,
						     nCol, nCol+25,
						     sizeof(iteratorCol) );
     lParse->varData = (DataInfo   *) fits_recalloc( lParse->varData,
						     nCol, nCol+25,
						     sizeof(DataInfo) );

     memset( (lParse->colData + nCol), 0x00, sizeof(iteratorCol)*25 );
     memset( (lParse->varData + nCol), 0x00, sizeof(DataInfo)*25    );

     if(    lParse->colData  == NULL
	    || lParse->varData  == NULL    ) {
         if( lParse->colData  ) free(lParse->colData);
         if( lParse->varData  ) free(lParse->varData);
         lParse->colData = NULL;
         lParse->varData = NULL;
         return( *status = MEMORY_ALLOCATION );
     }
   }
   lParse->varData[nCol].data  = NULL;
   lParse->varData[nCol].undef = NULL;
   return 0;
}

static int load_column( ParseData *lParse, int varNum, long fRow, long nRows,
                        void *data, char *undef )
{
   iteratorCol *var;
   long nelem,nbytes,row,len,idx;
   char **bitStrs, msg[80];
   unsigned char *bytes;
   int status = 0, anynul;

   var = lParse->colData+varNum;
   if (lParse->hdutype == IMAGE_HDU) {
    /* This test would need to be on a per varNum basis to support
     * cross HDU operations */
    fits_read_imgnull(var->fptr, var->datatype, fRow, nRows,
                data, undef, &anynul, &status);
    if (DEBUG_PIXFILTER)
        printf("load_column: IMAGE_HDU fRow=%ld, nRows=%ld => %d\n",
                        fRow, nRows, status);
  } else { 

   nelem = nRows * var->repeat;

   switch( var->datatype ) {
   case TBYTE:
      nbytes = ((var->repeat+7)/8) * nRows;
      bytes = (unsigned char *)malloc( nbytes * sizeof(char) );

      ffgcvb(var->fptr, var->colnum, fRow, 1L, nbytes,
             0, bytes, &anynul, &status);

      nelem = var->repeat;
      bitStrs = (char **)data;
      for( row=0; row<nRows; row++ ) {
         idx = (row)*( (nelem+7)/8 ) + 1;
         for(len=0; len<nelem; len++) {
            if( bytes[idx] & (1<<(7-len%8)) )
               bitStrs[row][len] = '1';
            else
               bitStrs[row][len] = '0';
            if( len%8==7 ) idx++;
         }
         bitStrs[row][len] = '\0';
      }

      FREE( (char *)bytes );
      break;
   case TSTRING:
      ffgcfs(var->fptr, var->colnum, fRow, 1L, nRows,
             (char **)data, undef, &anynul, &status);
      break;
   case TLOGICAL:
      ffgcfl(var->fptr, var->colnum, fRow, 1L, nelem,
             (char *)data, undef, &anynul, &status);
      break;
   case TLONG:
      ffgcfj(var->fptr, var->colnum, fRow, 1L, nelem,
             (long *)data, undef, &anynul, &status);
      break;
   case TDOUBLE:
      ffgcfd(var->fptr, var->colnum, fRow, 1L, nelem,
             (double *)data, undef, &anynul, &status);
      break;
   default:
      snprintf(msg,80,"load_column: unexpected datatype %d", var->datatype);
      ffpmsg(msg);
   }
  }
   if( status ) {
      lParse->status = status;
      return pERROR;
   }

   return 0;
}


/*--------------------------------------------------------------------------*/
int fits_pixel_filter (PixelFilter * filter, int * status)
/* Evaluate an expression using the data in the input FITS file(s)          */
/*--------------------------------------------------------------------------*/
{
   parseInfo Info = { 0 };
   int naxis, bitpix;
   long nelem, naxes[MAXDIMS];
   int col_cnt;
   Node *result;
   int datatype;
   fitsfile * infptr;
   fitsfile * outfptr;
   char * DEFAULT_TAGS[] = { "X" };
   char msg[256];
   int writeBlankKwd = 0;   /* write BLANK if any output nulls? */
   ParseData lParse;

   DEBUG_PIXFILTER = getenv("DEBUG_PIXFILTER") ? 1 : 0;

   memset(&Info, 0, sizeof(Info));   

   if (*status)
      return (*status);

   if (!filter->tag || !filter->tag[0] || !filter->tag[0][0]) {
      filter->tag = DEFAULT_TAGS;
      if (DEBUG_PIXFILTER)
         printf("using default tag '%s'\n", filter->tag[0]);
   }

   infptr = filter->ifptr[0];
   outfptr = filter->ofptr;
   lParse.pixFilter = filter;

   if (ffiprs(infptr, 0, filter->expression, MAXDIMS,
	      &Info.datatype, &nelem, &naxis, naxes, &lParse, status)) {
      goto CLEANUP;
   }


   if (nelem < 0) {
      nelem = -nelem;
   }

   {
      /* validate result type */
      const char * type = 0;
      switch (Info.datatype) {
         case TLOGICAL:  type = "LOGICAL"; break;
         case TLONG:     type = "LONG"; break;
         case TDOUBLE:   type = "DOUBLE"; break;
         case TSTRING:   type = "STRING";
                         *status = pERROR;
                         ffpmsg("pixel_filter: cannot have string image");
         case TBIT:      type = "BIT";
                         if (DEBUG_PIXFILTER)
                            printf("hmm, image from bits?\n");
                         break;
         default:       type = "UNKNOWN?!";
                        *status = pERROR;
                        ffpmsg("pixel_filter: unexpected result datatype");
      }
      if (DEBUG_PIXFILTER)
         printf("result type is %s [%d]\n", type, Info.datatype);
      if (*status)
         goto CLEANUP;
   }

   if (fits_get_img_param(infptr, MAXDIMS,
            &bitpix, &naxis, &naxes[0], status)) {
      ffpmsg("pixel_filter: unable to read input image parameters");
      goto CLEANUP;
   }

   if (DEBUG_PIXFILTER)
      printf("input bitpix %d\n", bitpix);

   if (Info.datatype == TDOUBLE) {
       /*  for floating point expressions, set the default output image to
           bitpix = -32 (float) unless the default is already a double */
       if (bitpix != DOUBLE_IMG)
           bitpix = FLOAT_IMG;
   }

   /* override output image bitpix if specified by caller */
   if (filter->bitpix)
      bitpix = filter->bitpix;
   if (DEBUG_PIXFILTER)
      printf("output bitpix %d\n", bitpix);

   if (fits_create_img(outfptr, bitpix, naxis, naxes, status)) {
      ffpmsg("pixel_filter: unable to create output image");
      goto CLEANUP;
   }

   /* transfer keycards */
   {
      int i, ncards, more;
      if (fits_get_hdrspace(infptr, &ncards, &more, status)) {
         ffpmsg("pixel_filter: unable to determine number of keycards");
         goto CLEANUP;
      }

      for (i = 1; i <= ncards; ++i) {

         int keyclass;
         char card[FLEN_CARD];

         if (fits_read_record(infptr, i, card, status)) {
            snprintf(msg, 256,"pixel_filter: unable to read keycard %d", i);
            ffpmsg(msg);
            goto CLEANUP;
         }

         keyclass = fits_get_keyclass(card);
         if (keyclass == TYP_STRUC_KEY) {
            /* output structure defined by fits_create_img */
         }
         else if (keyclass == TYP_COMM_KEY && i < 12) {
            /* assume this is one of the FITS standard comments */
         }
         else if (keyclass == TYP_NULL_KEY && bitpix < 0) {
            /* do not transfer BLANK to real output image */
         }
         else if (keyclass == TYP_SCAL_KEY && bitpix < 0) {
            /* do not transfer BZERO, BSCALE to real output image */
         }
         else if (fits_write_record(outfptr, card, status)) {
            snprintf(msg,256, "pixel_filter: unable to write keycard '%s' [%d]\n",
                        card, *status);
            ffpmsg(msg);
            goto CLEANUP;
         }
      }
   }

   switch (bitpix) {
      case BYTE_IMG: datatype = TLONG; Info.datatype = TBYTE; break;
      case SHORT_IMG: datatype = TLONG; Info.datatype = TSHORT; break;
      case LONG_IMG: datatype = TLONG; Info.datatype = TLONG; break;
      case FLOAT_IMG: datatype = TDOUBLE; Info.datatype = TFLOAT; break;
      case DOUBLE_IMG: datatype = TDOUBLE; Info.datatype = TDOUBLE; break;

      default:
           snprintf(msg, 256,"pixel_filter: unexpected output bitpix %d\n", bitpix);
           ffpmsg(msg);
           *status = pERROR;
           goto CLEANUP;
   }

   if (bitpix > 0) { /* arrange for NULLs in output */
      long nullVal = filter->blank;
      if (!filter->blank) {
         int tstatus = 0;
         if (fits_read_key_lng(infptr, "BLANK", &nullVal, 0, &tstatus)) {

            writeBlankKwd = 1;

            if (bitpix == BYTE_IMG)
                nullVal = UCHAR_MAX;
            else if (bitpix == SHORT_IMG)
                nullVal = SHRT_MIN;
            else if (bitpix == LONG_IMG) {
               if (sizeof(long) == 8 && sizeof(int) == 4)
                  nullVal = INT_MIN;
               else
                  nullVal = LONG_MIN;
            }
            else
                printf("unhandled positive output BITPIX %d\n", bitpix);
         }

         filter->blank = nullVal;
      }

      fits_set_imgnull(outfptr, filter->blank, status);
      if (DEBUG_PIXFILTER)
         printf("using blank %ld\n", nullVal);

   }

   if (!filter->keyword[0]) {
      iteratorCol * colIter;
      DataInfo * varInfo;

      /*************************************/
      /* Create new iterator Output Column */
      /*************************************/
      col_cnt = lParse.nCols;
      if (fits_parser_allocateCol(&lParse, col_cnt, status))
         goto CLEANUP;
      lParse.nCols++;

      colIter = &lParse.colData[col_cnt];
      colIter->fptr = filter->ofptr;
      colIter->iotype = OutputCol;
      varInfo = &lParse.varData[col_cnt];
      set_image_col_types(&lParse, colIter->fptr, "CREATED", bitpix, varInfo, colIter);

      Info.maxRows = -1;
      Info.parseData = &lParse;

      if (ffiter(lParse.nCols, lParse.colData, 0,
                     0, fits_parser_workfn, &Info, status) == -1)
            *status = 0;
      else if (*status)
         goto CLEANUP;

      if (Info.anyNull) {
         if (writeBlankKwd) {
            fits_update_key_lng(outfptr, "BLANK", filter->blank, "NULL pixel value", status);
            if (*status)
                ffpmsg("pixel_filter: unable to write BLANK keyword");
            if (DEBUG_PIXFILTER) {
                printf("output has NULLs\n");
                printf("wrote blank [%d]\n", *status);
            }
         }
      }
      else if (bitpix > 0) /* never used a null */
         if (fits_set_imgnull(outfptr, -1234554321, status))
            ffpmsg("pixel_filter: unable to reset imgnull");
   }
   else {

      /* Put constant result into keyword */
      char * parName = filter->keyword;
      char * parInfo = filter->comment;

      result  = lParse.Nodes + lParse.resultNode;
      switch (Info.datatype) {
      case TDOUBLE:
         ffukyd(outfptr, parName, result->value.data.dbl, 15, parInfo, status);
         break;
      case TLONG:
         ffukyj(outfptr, parName, result->value.data.lng, parInfo, status);
         break;
      case TLOGICAL:
         ffukyl(outfptr, parName, result->value.data.log, parInfo, status);
         break;
      case TBIT:
      case TSTRING:
         ffukys(outfptr, parName, result->value.data.str, parInfo, status);
         break;
      default:
         snprintf(msg, 256,"pixel_filter: unexpected constant result type [%d]\n",
                Info.datatype);
         ffpmsg(msg);
      }
   }

CLEANUP:
   ffcprs(&lParse);
   return (*status);
}
