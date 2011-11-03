/* $Id$ 
*/

/* "compression module */

/*****************************************************************************/
/*                                                                           */
/* The compression software is a python module implemented in C that, when    */
/* accessed through the pyfits module, supports the storage of compressed    */
/* images in FITS binary tables.  An n-dimensional image is divided into a   */
/* rectabgular grid of subimages or 'tiles'.  Each tile is then compressed   */
/* as a continuous block of data, and the resulting compressed byte stream   */
/* is stored in a row of a variable length column in a FITS binary table.    */
/* The default tiling pattern treates each row of a 2-dimensional image      */
/* (or higher dimensional cube) as a tile, such that each tile contains      */
/* NAXIS1 pixels.                                                            */
/*                                                                           */
/* This module contains two functions that are callable from python.  The    */
/* first is compressData.  This function takes a numpy.ndarray object        */
/* containing the uncompressed image data and returns a list of byte streams */
/* in which each element of the list holds the compressed data for a single  */
/* tile of the image.  The second function is decompressData.  It takes a    */
/* list of byte streams that hold the compressed data for each tile in the   */
/* image.  It returns a list containing the decompressed data for each tile  */
/* in the image.                                                             */
/*                                                                           */
/* Copyright (C) 2004 Association of Universities for Research in Astronomy  */
/* (AURA)                                                                    */
/*                                                                           */
/* Redistribution and use in source and binary forms, with or without        */
/* modification, are permitted provided that the following conditions are    */
/* met:                                                                      */
/*                                                                           */
/*    1. Redistributions of source code must retain the above copyright      */
/*      notice, this list of conditions and the following disclaimer.        */
/*                                                                           */
/*    2. Redistributions in binary form must reproduce the above             */
/*      copyright notice, this list of conditions and the following          */
/*      disclaimer in the documentation and/or other materials provided      */
/*      with the distribution.                                               */
/*                                                                           */
/*    3. The name of AURA and its representatives may not be used to         */
/*      endorse or promote products derived from this software without       */
/*      specific prior written permission.                                   */
/*                                                                           */
/* THIS SOFTWARE IS PROVIDED BY AURA ``AS IS'' AND ANY EXPRESS OR IMPLIED    */
/* WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF      */
/* MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE                  */
/* DISCLAIMED. IN NO EVENT SHALL AURA BE LIABLE FOR ANY DIRECT, INDIRECT,    */
/* INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,      */
/* BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS     */
/* OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND    */
/* ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR     */
/* TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE    */
/* USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH          */
/* DAMAGE.                                                                   */
/*                                                                           */
/* Some of the source code used by this module was copied and modified from  */
/* the FITSIO software that was written by William Pence at the High Energy  */
/* Astrophysic Science Archive Research Center (HEASARC) at the NASA Goddard */
/* Space Flight Center.  That software contained the following copyright and */
/* warranty notices:                                                         */
/*                                                                           */
/* Copyright (Unpublished--all rights reserved under the copyright laws of   */
/* the United States), U.S. Government as represented by the Administrator   */
/* of the National Aeronautics and Space Administration.  No copyright is    */
/* claimed in the United States under Title 17, U.S. Code.                   */
/*                                                                           */
/* Permission to freely use, copy, modify, and distribute this software      */
/* and its documentation without fee is hereby granted, provided that this   */
/* copyright notice and disclaimer of warranty appears in all copies.        */
/*                                                                           */
/* DISCLAIMER:                                                               */
/*                                                                           */
/* THE SOFTWARE IS PROVIDED 'AS IS' WITHOUT ANY WARRANTY OF ANY KIND,        */
/* EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT LIMITED TO,   */
/* ANY WARRANTY THAT THE SOFTWARE WILL CONFORM TO SPECIFICATIONS, ANY        */
/* IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR           */
/* PURPOSE, AND FREEDOM FROM INFRINGEMENT, AND ANY WARRANTY THAT THE         */
/* DOCUMENTATION WILL CONFORM TO THE SOFTWARE, OR ANY WARRANTY THAT THE      */
/* SOFTWARE WILL BE ERROR FREE.  IN NO EVENT SHALL NASA BE LIABLE FOR ANY    */
/* DAMAGES, INCLUDING, BUT NOT LIMITED TO, DIRECT, INDIRECT, SPECIAL OR      */
/* CONSEQUENTIAL DAMAGES, ARISING OUT OF, RESULTING FROM, OR IN ANY WAY      */
/* CONNECTED WITH THIS SOFTWARE, WHETHER OR NOT BASED UPON WARRANTY,         */
/* CONTRACT, TORT , OR OTHERWISE, WHETHER OR NOT INJURY WAS SUSTAINED BY     */
/* PERSONS OR PROPERTY OR OTHERWISE, AND WHETHER OR NOT LOSS WAS SUSTAINED   */
/* FROM, OR AROSE OUT OF THE RESULTS OF, OR USE OF, THE SOFTWARE OR          */
/* SERVICES PROVIDED HEREUNDER."                                             */
/*                                                                           */
/*****************************************************************************/

/* Include the Python C API */

#include "Python.h"
#include <numpy/arrayobject.h>
#include "fitsio.h"
#include "string.h"

/* Some defines for Python3 support--bytes objects should be used where */
/* strings were previously used                                         */
#if PY_MAJOR_VERSION >= 3
#define PyString_AsString PyBytes_AsString
#define PyString_FromStringAndSize PyBytes_FromStringAndSize
#define PyString_Size PyBytes_Size
#endif

/* Function to get the input long values from the input list */

static long* get_long_array(PyObject* data, const char* description,
                            int* data_size)
{
   int    i;
   int    size;
   long*  out;
   int    seq;
   char   errMsg[80];

   seq = PyList_Check(data);

   if (!seq)
   {
      strncpy(errMsg,description,79);
      strncat(errMsg," argument must be a list.",79-strlen(errMsg));
      PyErr_SetString(PyExc_TypeError,errMsg);
      return NULL;
   }

   size = PyList_Size(data);

   if (size < 0)
   {
      strncpy(errMsg,description,79);
      strncat(errMsg," list has invalid size.",79-strlen(errMsg));
      PyErr_SetString(PyExc_ValueError,errMsg);
      return NULL;
   }

   if (data_size)
   {
      *data_size = size;
   }

   out = (long*) PyMem_Malloc(size * sizeof(long));

   if(!out)
   {
      PyErr_NoMemory();
      return NULL;
   }

   for (i = 0; i < size; i++)
   {
      out[i] = PyLong_AsLong(PyList_GetItem(data, i));
   }

   if ( PyErr_Occurred())
   {
      PyMem_Free(out);
      out = NULL;
   }

   return out;
}


/* Function to get the input character arrays from the input list */

static unsigned char** get_char_array(PyObject* data, const char* description,
                                      int* data_size, int** dataLen)
{
   int             i;
   int             size;
   unsigned char** out;
   int             seq;
   char            errMsg[80];

   seq = PyList_Check(data);

   if (!seq)
   {
      strncpy(errMsg,description,79);
      strncat(errMsg," argument must be a list.",79-strlen(errMsg));
      PyErr_SetString(PyExc_TypeError,errMsg);
      return NULL;
   }

   size = PyList_Size(data);

   if ( size < 0)
   {
      strncpy(errMsg,description,79);
      strncat(errMsg," list has invalid size.",79-strlen(errMsg));
      PyErr_SetString(PyExc_ValueError,errMsg);
      return NULL;
   }

   if (data_size)
   {
      *data_size = size;
   }

   out = (unsigned char**) PyMem_Malloc(size * sizeof(char*));

   if(!out)
   {
      PyErr_NoMemory();
      return NULL;
   }

   *dataLen = (int*) PyMem_Malloc(size * sizeof(int));

   if(!(*dataLen))
   {
      PyMem_Free(out);
      PyErr_NoMemory();
      return NULL;
   }

   for (i = 0; i < size; i++)
   {
      out[i] = (unsigned char*)PyString_AsString(PyList_GetItem(data,i));
      (*dataLen)[i] = PyString_Size(PyList_GetItem(data,i));
   }

   if (PyErr_Occurred())
   {
      PyMem_Free(out);
      PyMem_Free(*dataLen);
      out = NULL;
   }

   return out;
}

/* Function to get the input float values from the input list */

static float* get_float_array(PyObject* data, const char* description,
                              int* data_size)
{
   int    i;
   int    size;
   float* out;
   int    seq;
   char   errMsg[80];

   seq = PyList_Check(data);

   if (!seq)
   {
      strncpy(errMsg,description,79);
      strncat(errMsg," argument must be a list.",79-strlen(errMsg));
      PyErr_SetString(PyExc_TypeError,errMsg);
      return NULL;
   }

   size = PyList_Size(data);

   if (size < 0)
   {
      strncpy(errMsg,description,79);
      strncat(errMsg," list has invalid size.",79-strlen(errMsg));
      PyErr_SetString(PyExc_ValueError,errMsg);
      return NULL;
   }

   if (data_size)
   {
      *data_size = size;
   }

   out = (float*) PyMem_Malloc(size * sizeof(float));

   if(!out)
   {
      PyErr_NoMemory();
      return NULL;
   }

   for (i = 0; i < size; i++)
   {
      out[i] = PyFloat_AsDouble(PyList_GetItem(data, i));
   }

   if ( PyErr_Occurred())
   {
      PyMem_Free(out);
      out = NULL;
   }

   return out;
}

/* Function to get the input double values from the input list */

static double* get_double_array(PyObject* data, const char* description,
                                int* data_size)
{
   int     i;
   int     size;
   double* out;
   int     seq;
   char    errMsg[80];

   seq = PyList_Check(data);

   if (!seq)
   {
      strncpy(errMsg,description,79);
      strncat(errMsg," argument must be a list.",79-strlen(errMsg));
      PyErr_SetString(PyExc_TypeError,errMsg);
      return NULL;
   }

   size = PyList_Size(data);

   if (size < 0)
   {
      strncpy(errMsg,description,79);
      strncat(errMsg," list has invalid size.",79-strlen(errMsg));
      PyErr_SetString(PyExc_ValueError,errMsg);
      return NULL;
   }

   if (data_size)
   {
      *data_size = size;
   }

   out = (double*) PyMem_Malloc(size * sizeof(double));

   if(!out)
   {
      PyErr_NoMemory();
      return NULL;
   }

   for (i = 0; i < size; i++)
   {
      out[i] = PyFloat_AsDouble(PyList_GetItem(data, i));
   }

   if ( PyErr_Occurred())
   {
      PyMem_Free(out);
      out = NULL;
   }

   return out;
}

/* Report any error based on the status returned from cfitsio. */

void processStatusErr(int status)
{
   PyObject* exceptType;
   char      errMsg[81];
   char      defErrMsg[81];

   errMsg[0] = '\0';
   defErrMsg[0] = '\0';

   switch (status)
   {
      case MEMORY_ALLOCATION:
         exceptType = PyExc_MemoryError;
         break;
      case OVERFLOW_ERR:
         exceptType = PyExc_OverflowError;
         break;
      case BAD_COL_NUM:
         strcpy(defErrMsg,"bad column number");
         exceptType = PyExc_ValueError;
         break;
      case BAD_PIX_NUM:
         strcpy(defErrMsg,"bad pixel number");
         exceptType = PyExc_ValueError;
         break;
      case NEG_AXIS:
         strcpy(defErrMsg,"negative axis number");
         exceptType = PyExc_ValueError;
         break;
      case BAD_DATATYPE:
         strcpy(defErrMsg,"bad data type");
         exceptType = PyExc_TypeError;
         break;
      case NO_COMPRESSED_TILE:
         strcpy(defErrMsg,"no compressed or uncompressed data for tile.");
         exceptType = PyExc_ValueError;
         break;
      default:
         exceptType = PyExc_RuntimeError;
   }

   if (_pyfits_ffgmsg(errMsg))
   {
      PyErr_SetString(exceptType,errMsg);
   }
   else if (*defErrMsg)
   {
      PyErr_SetString(exceptType,defErrMsg);
   }
   else
   {
      PyErr_SetString(exceptType, "unknown error.");
   }
}

/* Wrapper for the _pyfits_fits_write_img() function */

PyObject* compression_compressData(PyObject* self, PyObject* args)
{
   int             status;
   PyObject*       naxesObj;
   PyObject*       tileSizeObj;
   PyObject*       zvalObj;
   PyObject*       outList;
   PyObject*       outScale;
   PyObject*       outZero;
   PyObject*       outUncompressed;
   PyObject*       uncompressedTileDataList;
   PyObject*       returnTuple = NULL;
   long*           tileSize = 0;
   long*           zval = 0;
   int             i;
   int             loop;
   int             numzVals;
   char*           compressTypeStr = "";
   int             datatype;
   int             bitpix;
   int             firstelem;
   int             naxis;
   int             ntiles;
   int             ii;
   int             zblank;
   int             cn_zblank;
   int             cn_zscale;
   int             cn_zzero;
   int             cn_uncompressed;
   double          cn_bscale;
   double          cn_bzero;
   double          quantize_level;
   double          hcomp_scale;
   long*           naxes = 0;
   long            nelem;

   FITSfile        fileParms;
   fitsfile        theFile;

   PyArrayObject*  array;

   status = 0;

   /* Get Python arguments */

   if (!PyArg_ParseTuple(args, "O!iOOiiddiiiddOsiil:compression.compressData",
                         &PyArray_Type, &array, &naxis, &naxesObj,
                         &tileSizeObj, &cn_zblank, &zblank, &cn_bscale, 
                         &cn_bzero, &cn_zscale, &cn_zzero, &cn_uncompressed,
                         &quantize_level, &hcomp_scale, &zvalObj,
                         &compressTypeStr, &bitpix, &firstelem, &nelem))
   {
      PyErr_SetString(PyExc_TypeError,"Couldn't parse agruments");
      return NULL;
   }

   /* Determine the data type based on the bitpix value from the header */

   switch (bitpix)
   {
      case BYTE_IMG:
         datatype = TBYTE;
         break;
      case SHORT_IMG:
         datatype = TSHORT;
         break;
      case LONG_IMG:
         datatype = TINT;
         break;
      case LONGLONG_IMG:
         datatype = TLONGLONG;
         break;
      case FLOAT_IMG:
         datatype = TFLOAT;
         break;
      case DOUBLE_IMG:
         datatype = TDOUBLE;
         break;
      default:
         PyErr_SetString(PyExc_ValueError,"Invalid value for BITPIX");
         return NULL;
   }

   /* Initialize allocated array pointers to zero so we can free them */
   /* without allocating memory for them.                             */

   theFile.Fptr = &fileParms;
   (theFile.Fptr)->c_zscale = 0;
   (theFile.Fptr)->c_zzero = 0;
   (theFile.Fptr)->ucDataLen = 0;
   (theFile.Fptr)->ucData = 0;
   (theFile.Fptr)->dataLen = 0;
   (theFile.Fptr)->data = 0;

   /* The loop allows you to break out if there is an error */

   for (loop = 0; loop == 0; loop++)
   {
      /* Convert the NAXISn, ZTILEn, and ZVALn lists into a C type arrays */

      naxes = get_long_array(naxesObj, "ZNAXISn", NULL);

      if (!naxes)
      {
         break;
      }

      tileSize = get_long_array(tileSizeObj, "ZTILEn", NULL);

      if (!tileSize)
      {
         break;
      }

      zval = get_long_array(zvalObj, "ZVALn", &numzVals);

      if (!zval)
      {
         break;
      }

      /* Set up the fitsfile object */
      /* Store the compression type and compression parameters */

      (theFile.Fptr)->hcomp_smooth = 0;
      (theFile.Fptr)->rice_blocksize = 32;
      (theFile.Fptr)->rice_bytepix = 4;

      if (strcmp(compressTypeStr, "RICE_1") == 0)
      {
         (theFile.Fptr)->compress_type = RICE_1;

         if (numzVals > 0)
         {
            (theFile.Fptr)->rice_blocksize = zval[0];

            if (numzVals > 1)
            {
               (theFile.Fptr)->rice_bytepix = zval[1];
 
            }
         }
      }
      else if (strcmp(compressTypeStr, "GZIP_1") == 0)
      {
         (theFile.Fptr)->compress_type = GZIP_1;

      }
      else if (strcmp(compressTypeStr, "HCOMPRESS_1") == 0)
      {
         (theFile.Fptr)->compress_type = HCOMPRESS_1;

         if (numzVals > 0)
         {
           (theFile.Fptr)->hcomp_smooth = zval[1];
         }
      }
      else if (strcmp(compressTypeStr, "PLIO_1") == 0)
      {
         (theFile.Fptr)->compress_type = PLIO_1;

      }
      else
      {
         (theFile.Fptr)->compress_type = 0;
      }

      (theFile.Fptr)->zndim = naxis;
      (theFile.Fptr)->maxtilelen = 1;
      (theFile.Fptr)->zbitpix = bitpix;
      (theFile.Fptr)->cn_zblank = cn_zblank;
      (theFile.Fptr)->zblank = zblank;
      (theFile.Fptr)->cn_bscale = cn_bscale;
      (theFile.Fptr)->cn_bzero = cn_bzero;
      (theFile.Fptr)->cn_zscale = cn_zscale;
      (theFile.Fptr)->cn_zzero = cn_zzero;
      (theFile.Fptr)->quantize_level = quantize_level;
      (theFile.Fptr)->hcomp_scale = hcomp_scale;

      /* Initialize arrays */

      for (ii = 0; ii < MAX_COMPRESS_DIM; ii++)
      {
         ((theFile.Fptr)->tilesize)[ii] = 1;
         ((theFile.Fptr)->znaxis)[ii] = 1;
      }

      ntiles = 1;

      for (ii = 0; ii < naxis; ii++)
      {
         ((theFile.Fptr)->znaxis)[ii] = naxes[ii];
         ((theFile.Fptr)->tilesize)[ii] = tileSize[ii];
         (theFile.Fptr)->maxtilelen *= tileSize[ii];
         ntiles *= (naxes[ii] - 1) / tileSize[ii] + 1;
      }

      (theFile.Fptr)->maxelem = _pyfits_imcomp_calc_max_elem(
                                 (theFile.Fptr)->compress_type,
                                 (theFile.Fptr)->maxtilelen,
                                 (theFile.Fptr)->zbitpix,
                                 (theFile.Fptr)->rice_blocksize);

      if (cn_zscale > 0)
      {
         (theFile.Fptr)->c_zscale = 
                         (double*)PyMem_Malloc(ntiles * sizeof(double));
         (theFile.Fptr)->c_zzero = 
                         (double*)PyMem_Malloc(ntiles * sizeof(double));

         if(!(theFile.Fptr)->c_zzero)
         {
            PyErr_NoMemory();
            break;
         }
      }

      (theFile.Fptr)->cn_uncompressed = cn_uncompressed;

      if (cn_uncompressed > 0)
      {
         (theFile.Fptr)->ucDataLen = 
                          (int*)PyMem_Malloc(ntiles * sizeof(int));
         (theFile.Fptr)->ucData =
                          (void**)PyMem_Malloc(ntiles * sizeof(double*));

         if(!(theFile.Fptr)->ucData)
         {
            PyErr_NoMemory();
            break;
         }

         for (i = 0; i < ntiles; i++)
         {
            (theFile.Fptr)->ucDataLen[i] = 0;
            (theFile.Fptr)->ucData[i] = 0;
         }
      }

      (theFile.Fptr)->dataLen = 
                         (int*) PyMem_Malloc(ntiles * sizeof(int));
      (theFile.Fptr)->data =
                         (unsigned char**) PyMem_Malloc(ntiles*sizeof(char*));

      if(!(theFile.Fptr)->data)
      {
         PyErr_NoMemory();
         break;
      }

      for (i = 0; i < ntiles; i++)
      {
         (theFile.Fptr)->dataLen[i] = 0;
         (theFile.Fptr)->data[i] = 0;
      }

      status = _pyfits_fits_write_img(&theFile, datatype, firstelem,
                                      nelem, (void*)array->data, &status);

      if (status == 0)
      {
         outList = PyList_New(0);
         outScale = PyList_New(0);
         outZero = PyList_New(0);
         outUncompressed = PyList_New(0);

         for ( i = 0; i < ntiles; i++)
         {
            PyList_Append(outList, PyString_FromStringAndSize(
                          (const char*)((theFile.Fptr)->data[i]),
                          (theFile.Fptr)->dataLen[i]));
            free((theFile.Fptr)->data[i]);

            if (cn_zscale > 0)
            {
               PyList_Append(outScale, 
                             PyFloat_FromDouble((theFile.Fptr)->c_zscale[i]));
               PyList_Append(outZero, 
                             PyFloat_FromDouble((theFile.Fptr)->c_zzero[i]));
            }

            if (cn_uncompressed > 0)
            {
               uncompressedTileDataList = PyList_New(0);
   
               for (ii = 0; ii < (theFile.Fptr)->ucDataLen[i]; ii++)
               {
                   PyList_Append(uncompressedTileDataList, 
                   PyFloat_FromDouble(
                               ((double**)((theFile.Fptr)->ucData))[i][ii]));
               }

               free((theFile.Fptr)->ucData[i]);
               PyList_Append(outUncompressed, uncompressedTileDataList);
            }
         }
      }
      else
      {
         processStatusErr(status);
         break;
      }

      returnTuple = PyTuple_New(5);
      PyTuple_SetItem(returnTuple, 0, Py_BuildValue("i",status));
      PyTuple_SetItem(returnTuple, 1, outList);
      PyTuple_SetItem(returnTuple, 2, outScale);
      PyTuple_SetItem(returnTuple, 3, outZero);
      PyTuple_SetItem(returnTuple, 4, outUncompressed);
   }
   
   /* Free any allocated memory */

   PyMem_Free((theFile.Fptr)->dataLen);
   PyMem_Free((theFile.Fptr)->data);
   PyMem_Free((theFile.Fptr)->c_zscale);
   PyMem_Free((theFile.Fptr)->c_zzero);
   PyMem_Free((theFile.Fptr)->ucData);
   PyMem_Free((theFile.Fptr)->ucDataLen);
   PyMem_Free(naxes);
   PyMem_Free(tileSize);
   PyMem_Free(zval);

   if (loop == 0)
   {
      /* Error has occurred */

      return NULL;
   }
   else
   {
      return returnTuple;
   }
}

/* Wrapper for the _pyfits_fits_read_img() function */

PyObject* compression_decompressData(PyObject* self, PyObject* args)
{
   int             status;
   int             nUcTiles = 0;
   int*            inDataLen = 0;
   int             naxis;
   int             numzVals;
   int*            numUncompressedVals = 0;
   int             cn_zblank;
   int             cn_zscale;
   int             cn_zzero;
   int             cn_uncompressed;
   int             bitpix;
   int             datatype;
   int             firstelem;
   int             anynul;
   long            nelem;
   long*           naxes = 0;
   long*           tileSize = 0;
   long*           zval = 0;
   double          nulval;
   double*         bscale;
   double*         bzero;
   double          quantize_level;
   double          hcomp_scale;
   long*           nullDVals;
   unsigned char** inData = 0;
   void**          uncompressedData = 0;
   int             i;
   int             ii;
   char*           compressTypeStr = "";

   PyObject*       inDataObj;
   PyObject*       naxesObj;
   PyObject*       tileSizeObj;
   PyObject*       uncompressedDataObj;
   PyObject*       zvalObj;

   FITSfile        fileParms;
   fitsfile        theFile;

   PyArrayObject*  bscaleArray;
   PyArrayObject*  bzeroArray;
   PyArrayObject*  nullDvalsArray;
   PyArrayObject*  decompDataArray;

   PyArrayObject*  bscaleArray1 = 0;
   PyArrayObject*  bzeroArray1 = 0;
   PyArrayObject*  nullDvalsArray1 = 0;

   /* Get Python arguments */

   if (!PyArg_ParseTuple(args, 
                         "OiOOO!iO!iO!iOiddOsiildO!:compression.decompressData",
                         &inDataObj, 
                         &naxis, &naxesObj, &tileSizeObj, &PyArray_Type, 
                         &bscaleArray, &cn_zscale, &PyArray_Type, &bzeroArray,
                         &cn_zzero, &PyArray_Type, &nullDvalsArray, 
                         &cn_zblank, &uncompressedDataObj,
                         &cn_uncompressed, &quantize_level, &hcomp_scale,
                         &zvalObj, &compressTypeStr, &bitpix, &firstelem,
                         &nelem, &nulval, &PyArray_Type, &decompDataArray))
   {
      PyErr_SetString(PyExc_TypeError,"Couldn't parse arguments");
      return NULL;
   }

   /* Convert the input lists into C type arrays */

   inData = get_char_array(inDataObj, "Compressed Data", NULL, &inDataLen);

   if (!inData)
   {
      goto error;
   }

   naxes = get_long_array(naxesObj, "ZNAXISn", NULL);

   if (!naxes)
   {
      goto error;
   }

   tileSize = get_long_array(tileSizeObj, "ZTILEn", NULL);

   if (!tileSize)
   {
      goto error;
   }

   if (cn_zzero != 1)
   {
      bzero = (double*)bzeroArray->data;
   }
   else
   {
      bzeroArray1 = (PyArrayObject*)PyArray_ContiguousFromObject(
                     (PyObject*)bzeroArray, PyArray_DOUBLE, 1, 1);
      bzero = (double*)bzeroArray1->data;
   }

   if (cn_zscale != 1)
   {
      bscale = (double*)bscaleArray->data;
   }
   else
   {
      bscaleArray1 = (PyArrayObject*)PyArray_ContiguousFromObject(
                      (PyObject*)bscaleArray, PyArray_DOUBLE, 1, 1);
      bscale = (double*)bscaleArray1->data;
   }

   if (cn_zblank != 1)
   {
      nullDVals = (long*)nullDvalsArray->data;
   }
   else
   {
      nullDvalsArray1 = (PyArrayObject*)PyArray_ContiguousFromObject(
                         (PyObject*)nullDvalsArray, PyArray_LONG, 1, 1);
      nullDVals = (long*)nullDvalsArray1->data;
   }

   zval = get_long_array(zvalObj, "ZVALn", &numzVals);

   if (!zval)
   {
      goto error;
   }

   switch (bitpix)
   {
      case BYTE_IMG:
         datatype = TBYTE;
         break;
      case SHORT_IMG:
         datatype = TSHORT;
         break;
      case LONG_IMG:
         datatype = TINT;
         break;
      case LONGLONG_IMG:
         datatype = TLONGLONG;
         break;
      case FLOAT_IMG:
         datatype = TFLOAT;

         if (cn_uncompressed == 1)
         {
            nUcTiles = PyList_Size(uncompressedDataObj);
            uncompressedData = (void**) PyMem_Malloc(nUcTiles*sizeof(float*));

            if (!uncompressedData)
            {
               goto error;
            }

            numUncompressedVals = PyMem_Malloc(nUcTiles*sizeof(int));

            if (!numUncompressedVals)
            {
               goto error;
            }

            for (i = 0; i < nUcTiles; i++)
            {
                uncompressedData[i] = 
                       get_float_array(PyList_GetItem(uncompressedDataObj, i),
                                       "Uncompressed Data",
                                       &numUncompressedVals[i]);

                if (!uncompressedData[i])
                {
                   goto error;
                }
            }
         }
         break;
      case DOUBLE_IMG:
         datatype = TDOUBLE;

         if (cn_uncompressed == 1)
         {
            nUcTiles = PyList_Size(uncompressedDataObj);
            uncompressedData = (void**) PyMem_Malloc(nUcTiles*sizeof(double*));

            if (!uncompressedData)
            {
               goto error;
            }

            numUncompressedVals = PyMem_Malloc(nUcTiles*sizeof(int));

            if (!numUncompressedVals)
            {
               goto error;
            }

            for (i = 0; i < nUcTiles; i++)
            {
                uncompressedData[i] = 
                       get_double_array(PyList_GetItem(uncompressedDataObj, i),
                                        "Uncompressed Data",
                                        &numUncompressedVals[i]);

                if (!uncompressedData[i])
                {
                   goto error;
                }
            }
         }
         break;
      default:
         PyErr_SetString(PyExc_ValueError,"Invalid value for BITPIX");
         return NULL;
   }

   /* Set up the fitsfile object */

   theFile.Fptr = &fileParms;

   (theFile.Fptr)->rice_blocksize = 32;
   (theFile.Fptr)->hcomp_smooth = 0;
   (theFile.Fptr)->rice_bytepix = 4;

   if (strcmp(compressTypeStr, "RICE_1") == 0)
   {
      (theFile.Fptr)->compress_type = RICE_1;

      if (numzVals > 0)
      {
         (theFile.Fptr)->rice_blocksize = zval[0];

         if (numzVals > 1)
         {
            (theFile.Fptr)->rice_bytepix = zval[1];
         }
      }
   }
   else if (strcmp(compressTypeStr, "GZIP_1") == 0)
   {
      (theFile.Fptr)->compress_type = GZIP_1;
   }
   else if (strcmp(compressTypeStr, "HCOMPRESS_1") == 0)
   {
      (theFile.Fptr)->compress_type = HCOMPRESS_1;

      if (numzVals > 0)
      {
        (theFile.Fptr)->hcomp_smooth = zval[1];
      }
   }
   else if (strcmp(compressTypeStr, "PLIO_1") == 0)
   {
      (theFile.Fptr)->compress_type = PLIO_1;
   }
   else
   {
      (theFile.Fptr)->compress_type = 0;
   }

   (theFile.Fptr)->zndim = naxis;
   (theFile.Fptr)->maxtilelen = 1;
   (theFile.Fptr)->zbitpix = bitpix;
   (theFile.Fptr)->data = inData;
   (theFile.Fptr)->dataLen = inDataLen;

   (theFile.Fptr)->bscale = bscale;
   (theFile.Fptr)->cn_zscale = cn_zscale;
   (theFile.Fptr)->quantize_level = quantize_level;
   (theFile.Fptr)->hcomp_scale = hcomp_scale;

   if (cn_zscale == -1)
   {
      (theFile.Fptr)->zscale = bscale[0];
      (theFile.Fptr)->cn_bscale = bscale[0];
   }
   else
   {
      (theFile.Fptr)->zscale = 1.0;
      (theFile.Fptr)->cn_bscale = 1.0;
   }

   (theFile.Fptr)->bzero = bzero;
   (theFile.Fptr)->cn_zzero = cn_zzero;

   if (cn_zzero == -1)
   {
      (theFile.Fptr)->zzero = bzero[0];
      (theFile.Fptr)->cn_bzero = bzero[0];
   }
   else
   {
      (theFile.Fptr)->zzero = 0.0;
      (theFile.Fptr)->cn_bzero = 0.0;
   }

   (theFile.Fptr)->blank = nullDVals;
   (theFile.Fptr)->cn_zblank = cn_zblank;

   if (cn_zblank == -1)
   {
      (theFile.Fptr)->zblank = nullDVals[0];
   }
   else
   {
      (theFile.Fptr)->zblank = 0;
   }

   /* Initialize arrays */

   for (ii = 0; ii < MAX_COMPRESS_DIM; ii++)
   {
      ((theFile.Fptr)->tilesize)[ii] = 1;
      ((theFile.Fptr)->znaxis)[ii] = 1;
   }

   for (ii = 0; ii < naxis; ii++)
   {
      ((theFile.Fptr)->znaxis)[ii] = naxes[ii];
      ((theFile.Fptr)->tilesize)[ii] = tileSize[ii];
      (theFile.Fptr)->maxtilelen *= tileSize[ii];
   }

   (theFile.Fptr)->cn_uncompressed = cn_uncompressed;

   if (cn_uncompressed == 1)
   {
       (theFile.Fptr)->ucData = uncompressedData;
       (theFile.Fptr)->ucDataLen = numUncompressedVals;
   }

   /* Call the C function */

   status = 0;
   status = _pyfits_fits_read_img(&theFile, datatype, firstelem,
                                  nelem, &nulval, decompDataArray->data,
                                  &anynul, &status);

   if (status != 0)
   {
      processStatusErr(status);
   }

   error:
      PyMem_Free(inData);
      PyMem_Free(inDataLen);
      PyMem_Free(naxes);
      PyMem_Free(tileSize);
      PyMem_Free(zval);

      if (cn_uncompressed == 1)
      {
         for (i = 0; i < nUcTiles; i++)
         {
             PyMem_Free(uncompressedData[i]);
         }

         PyMem_Free(uncompressedData);
         PyMem_Free(numUncompressedVals);
      }

      if (bscaleArray1 != 0)
      {
         Py_DECREF(bscaleArray1);
      }

      if (bzeroArray1 != 0)
      {
         Py_DECREF(bzeroArray1);
      }

      if (nullDvalsArray1 != 0)
      {
         Py_DECREF(nullDvalsArray1);
      }
   
      if (status != 0)
      {
         return NULL;
      }
      else
      {
         return Py_BuildValue("i",status);
      }
}


/* Method table mapping names to wrappers */
static PyMethodDef compression_methods[] =
{
   {"decompressData", compression_decompressData, METH_VARARGS},
   {"compressData", compression_compressData, METH_VARARGS},
   {NULL,NULL}
};

#if PY_MAJOR_VERSION >=3
static struct PyModuleDef compressionmodule = {
    PyModuleDef_HEAD_INIT,
    "compression",
    "pyfits.compression module",
    -1, /* No global state */
    compression_methods
};

PyObject *
PyInit_compression(void)
{
    PyObject *module = PyModule_Create(&compressionmodule);
    import_array();
    return module;
}
#else
PyMODINIT_FUNC initcompression(void)
{
   Py_InitModule4("compression", compression_methods, "compression module",
                  NULL, PYTHON_API_VERSION);
   import_array();
}
#endif
