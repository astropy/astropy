/* "compression module */

/*****************************************************************************/
/*                                                                           */
/* The compression software is a python module implemented in C that, when   */
/* accessed through the astropy module, supports the storage of compressed   */
/* images in FITS binary tables.  An n-dimensional image is divided into a   */
/* rectangular grid of subimages or 'tiles'.  Each tile is then compressed   */
/* as a continuous block of data, and the resulting compressed byte stream   */
/* is stored in a row of a variable length column in a FITS binary table.    */
/* The default tiling pattern treates each row of a 2-dimensional image      */
/* (or higher dimensional cube) as a tile, such that each tile contains      */
/* NAXIS1 pixels.                                                            */
/*                                                                           */
/* This module contains three functions that are callable from python.  The  */
/* first is compress_hdu.  This function takes an                            */
/* astropy.io.fits.CompImageHDU object containing the uncompressed image     */
/* data and returns the compressed data for all tiles into the               */
/* .compressed_data attribute of that HDU.                                   */
/*                                                                           */
/* The second function is decompress_hdu.  It takes an                       */
/* astropy.io.fits.CompImageHDU object that already has compressed data in   */
/* its .compressed_data attribute.  It returns the decompressed image data   */
/* into the HDU's .data attribute.                                           */
/*                                                                           */
/* Copyright (C) 2013 Association of Universities for Research in Astronomy  */
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

#include <float.h>
#include <limits.h>
#include <math.h>
#include <string.h>

#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>
#include <fitsio2.h>
#include "compressionmodule.h"


/* These defaults mirror the defaults in astropy.io.fits.hdu.compressed */
#define DEFAULT_COMPRESSION_TYPE "RICE_1"
#define DEFAULT_QUANTIZE_LEVEL 16.0
#define DEFAULT_HCOMP_SCALE 0
#define DEFAULT_HCOMP_SMOOTH 0
#define DEFAULT_BLOCK_SIZE 32
#define DEFAULT_BYTE_PIX 4

/* Flags to pass to get_header_* functions to control error messages. */
typedef enum {
    HDR_NOFLAG = 0,
    HDR_FAIL_KEY_MISSING = 1 << 0,
    HDR_FAIL_VAL_NEGATIVE = 1 << 1,
} HeaderGetFlags;


/* Report any error based on the status returned from cfitsio. */
void process_status_err(int status)
{
   PyObject* except_type;
   char      err_msg[81];
   char      def_err_msg[81];

   err_msg[0] = '\0';
   def_err_msg[0] = '\0';

   switch (status) {
      case MEMORY_ALLOCATION:
         except_type = PyExc_MemoryError;
         break;
      case OVERFLOW_ERR:
         except_type = PyExc_OverflowError;
         break;
      case BAD_COL_NUM:
         strcpy(def_err_msg, "bad column number");
         except_type = PyExc_ValueError;
         break;
      case BAD_PIX_NUM:
         strcpy(def_err_msg, "bad pixel number");
         except_type = PyExc_ValueError;
         break;
      case NEG_AXIS:
         strcpy(def_err_msg, "negative axis number");
         except_type = PyExc_ValueError;
         break;
      case BAD_DATATYPE:
         strcpy(def_err_msg, "bad data type");
         except_type = PyExc_TypeError;
         break;
      case NO_COMPRESSED_TILE:
         strcpy(def_err_msg, "no compressed or uncompressed data for tile.");
         except_type = PyExc_ValueError;
         break;
      default:
         except_type = PyExc_RuntimeError;
         break;
   }

   if (fits_read_errmsg(err_msg)) {
      PyErr_SetString(except_type, err_msg);
   } else if (*def_err_msg) {
      PyErr_SetString(except_type, def_err_msg);
   } else {
      PyErr_Format(except_type, "unknown error %i.", status);
   }
}


void bitpix_to_datatypes(int bitpix, int* datatype, int* npdatatype) {
    /* Given a FITS BITPIX value, returns the appropriate CFITSIO type code and
       Numpy type code for that BITPIX into datatype and npdatatype
       respectively.
     */
    switch (bitpix) {
        case BYTE_IMG:
            *datatype = TBYTE;
            *npdatatype = NPY_UINT8;
            break;
        case SHORT_IMG:
            *datatype = TSHORT;
            *npdatatype = NPY_INT16;
            break;
        case LONG_IMG:
            *datatype = TINT;
            *npdatatype = NPY_INT32;
            break;
        case LONGLONG_IMG:
            *datatype = TLONGLONG;
            *npdatatype = NPY_LONGLONG;
            break;
        case FLOAT_IMG:
            *datatype = TFLOAT;
            *npdatatype = NPY_FLOAT;
            break;
        case DOUBLE_IMG:
            *datatype = TDOUBLE;
            *npdatatype = NPY_DOUBLE;
            break;
        default:
            PyErr_Format(PyExc_ValueError, "Invalid value for BITPIX: %d",
                         bitpix);
            break;
   }

   return;
}



int compress_type_from_string(char* zcmptype) {
    if (0 == strcmp(zcmptype, "RICE_1")) {
        return RICE_1;
    } else if (0 == strcmp(zcmptype, "GZIP_1")) {
        return GZIP_1;
    } else if (0 == strcmp(zcmptype, "GZIP_2")) {
        return GZIP_2;
    } else if (0 == strcmp(zcmptype, "PLIO_1")) {
        return PLIO_1;
    } else if (0 == strcmp(zcmptype, "HCOMPRESS_1")) {
        return HCOMPRESS_1;
    }
#ifdef CFITSIO_SUPPORTS_SUBTRACTIVE_DITHER_2
    /* CFITSIO adds a compression type alias for RICE_1 compression
       as a flag for using subtractive_dither_2 */
    else if (0 == strcmp(zcmptype, "RICE_ONE")) {
        return RICE_1;
    }
#endif
    else {
        PyErr_Format(PyExc_ValueError, "Unrecognized compression type: %s",
                     zcmptype);
        return -1;
    }
}


PyObject *
get_header_value(PyObject* header, const char* key, HeaderGetFlags flags) {
    PyObject* hdrkey;
    PyObject* hdrval;
    hdrkey = PyUnicode_FromString(key);
    if (hdrkey == NULL) {
        return NULL;
    }
    hdrval = PyObject_GetItem(header, hdrkey);
    Py_DECREF(hdrkey);
    if ((flags & HDR_FAIL_KEY_MISSING) == 0) {
        /* Normally we have a default so we want to ignore the exception in
           any case. But if the flag was given this step must be skipped. */
        PyErr_Clear();
    }
    return hdrval;
}


// TODO: It might be possible to simplify these further by making the
// conversion function (eg. PyString_AsString) an argument to a macro or
// something, but I'm not sure yet how easy it is to generalize the error
// handling
/* The get_header_* functions resemble "Header.get" where "def" is the default
   value, "keyword" is a string representing the header-key and the result is
   stored in "val".
   The function returns 0 on success, 1 if the header didn't have the keyword
   and the default was applied and -1 (with an exception set) if an Exception
   happened (like a MemoryError or Overflow).
*/
#define GET_HEADER_SUCCESS 0
#define GET_HEADER_DEFAULT_USED 1
#define GET_HEADER_FAILED -1
int get_header_string(PyObject* header, const char* keyword, char* val,
                      const char* def, HeaderGetFlags flags) {
    /* nonnegative doesn't make sense for strings*/
    assert(!(flags & HDR_FAIL_VAL_NEGATIVE));
    PyObject* keyval = get_header_value(header, keyword, flags);

    if (keyval == NULL) {
        strncpy(val, def, 72);
        return PyErr_Occurred() ? GET_HEADER_FAILED : GET_HEADER_DEFAULT_USED;
    }
    PyObject* tmp = PyUnicode_AsLatin1String(keyval);
    // FITS header values should always be ASCII, but Latin1 is on the
    // safe side
    Py_DECREF(keyval);
    if (tmp == NULL) {
        /* could always fail to allocate the memory or such like. */
        return GET_HEADER_FAILED;
    }
    strncpy(val, PyBytes_AsString(tmp), 72);
    Py_DECREF(tmp);
    return GET_HEADER_SUCCESS;
}


int get_header_long(PyObject* header, const char* keyword, long* val, long def,
                    HeaderGetFlags flags) {
    PyObject* keyval = get_header_value(header, keyword, flags);

    if (keyval == NULL) {
        *val = def;
        return PyErr_Occurred() ? GET_HEADER_FAILED : GET_HEADER_DEFAULT_USED;
    }
    long tmp = PyLong_AsLong(keyval);
    Py_DECREF(keyval);
    if (PyErr_Occurred()) {
        return GET_HEADER_FAILED;
    }
    if ((flags & HDR_FAIL_VAL_NEGATIVE) && (tmp < 0)) {
        PyErr_Format(PyExc_ValueError, "%s should not be negative.", keyword);
        return GET_HEADER_FAILED;
    }
    *val = tmp;
    return GET_HEADER_SUCCESS;
}


int get_header_int(PyObject* header, const char* keyword, int* val, int def,
                   HeaderGetFlags flags) {
    long tmp;
    int ret = get_header_long(header, keyword, &tmp, def, flags);
    if (ret == GET_HEADER_SUCCESS) {
        if (tmp >= INT_MIN && tmp <= INT_MAX) {
            *val = (int) tmp;
        } else {
            PyErr_Format(PyExc_OverflowError, "Cannot convert %ld to C 'int'", tmp);
            ret = GET_HEADER_FAILED;
        }
    }
    return ret;
}


int get_header_double(PyObject* header, const char* keyword, double* val,
                      double def, HeaderGetFlags flags) {
    /* nonnegative isn't currently used for doubles/floats. But if needed one
       could simply remove the assert again and implement the negative check. */
    assert(!(flags & HDR_FAIL_VAL_NEGATIVE));
    PyObject* keyval = get_header_value(header, keyword, flags);

    if (keyval == NULL) {
        *val = def;
        return PyErr_Occurred() ? GET_HEADER_FAILED : GET_HEADER_DEFAULT_USED;
    }
    double tmp = PyFloat_AsDouble(keyval);
    Py_DECREF(keyval);
    if (PyErr_Occurred()) {
        return GET_HEADER_FAILED;
    }
    *val = tmp;
    return GET_HEADER_SUCCESS;
}


int get_header_float(PyObject* header, const char* keyword, float* val,
                     float def, HeaderGetFlags flags) {
    double tmp;
    int ret = get_header_double(header, keyword, &tmp, def, flags);
    if (ret == GET_HEADER_SUCCESS) {
        if (tmp == 0.0 || (fabs(tmp) >= FLT_MIN && fabs(tmp) <= FLT_MAX)) {
            *val = (float) tmp;
        } else {
            PyErr_SetString(PyExc_OverflowError,
                            "Cannot convert 'double' to 'float'");
            ret = GET_HEADER_FAILED;
        }
    }
    return ret;
}


int get_header_longlong(PyObject* header, const char* keyword, long long* val,
                        long long def, HeaderGetFlags flags) {
    PyObject* keyval = get_header_value(header, keyword, flags);

    if (keyval == NULL) {
        *val = def;
        return PyErr_Occurred() ? GET_HEADER_FAILED : GET_HEADER_DEFAULT_USED;
    }
    long long tmp = PyLong_AsLongLong(keyval);
    Py_DECREF(keyval);
    if (PyErr_Occurred()) {
        return GET_HEADER_FAILED;
    }
    if ((flags & HDR_FAIL_VAL_NEGATIVE) && (tmp < 0)) {
        PyErr_Format(PyExc_ValueError, "%s should not be negative.", keyword);
        return GET_HEADER_FAILED;
    }
    *val = tmp;
    return GET_HEADER_SUCCESS;
}


void tcolumns_from_header(fitsfile* fileptr, PyObject* header,
                          tcolumn** columns) {
    // Creates the array of tcolumn structures from the table column keywords
    // read from the astropy.io.fits.Header object; caller is responsible for
    // freeing the memory allocated for this array

    tcolumn* column;
    char tkw[9];

    int tfields;
    char ttype[72];
    char tform[72];
    int dtcode;
    long trepeat;
    long twidth;
    long long totalwidth;
    int status = 0;
    int idx;

    if (get_header_int(header, "TFIELDS", &tfields, 0, HDR_FAIL_VAL_NEGATIVE) == GET_HEADER_FAILED) {
        return;
    }
    /* To avoid issues in the loop we need to limit the number of TFIELDs to
       999. Otherwise we would exceed the maximum length of the keyword name of
       8. This could lead to multiple accesses of the same header keyword with
       snprintf because we limit it to 8 characters + null-termination. */
    if (tfields > 999) {
        PyErr_SetString(PyExc_ValueError, "The TFIELDS value exceeds 999.");
        return;
    }

    // This used to use PyMem_New, but don't do that; CFITSIO will later
    // free() this object when the file is closed, so just use malloc here
    // *columns = column = PyMem_New(tcolumn, (size_t) tfields);
    *columns = column = calloc((size_t) tfields, sizeof(tcolumn));
    if (column == NULL) {
        PyErr_SetString(PyExc_MemoryError,
                        "Couldn't allocate memory for columns.");
        return;
    }


    for (idx = 1; idx <= tfields; idx++, column++) {
        /* set some invalid defaults */
        column->ttype[0] = '\0';
        column->tbcol = 0;
        column->tdatatype = -9999; /* this default used by cfitsio */
        column->trepeat = 1;
        column->strnull[0] = '\0';
        column->tform[0] = '\0';
        column->twidth = 0;

        snprintf(tkw, 9, "TTYPE%u", idx);
        if (get_header_string(header, tkw, ttype, "", HDR_NOFLAG) == GET_HEADER_FAILED) {
            return;
        }
        strncpy(column->ttype, ttype, 69);
        column->ttype[69] = '\0';

        snprintf(tkw, 9, "TFORM%u", idx);
        if (get_header_string(header, tkw, tform, "", HDR_NOFLAG) == GET_HEADER_FAILED) {
            return;
        }
        strncpy(column->tform, tform, 9);
        column->tform[9] = '\0';
        fits_binary_tform(tform, &dtcode, &trepeat, &twidth, &status);
        if (status != 0) {
            process_status_err(status);
            return;
        }

        column->tdatatype = dtcode;
        column->trepeat = trepeat;
        column->twidth = twidth;

        snprintf(tkw, 9, "TSCAL%u", idx);
        if (get_header_double(header, tkw, &(column->tscale), 1.0, HDR_NOFLAG) == GET_HEADER_FAILED) {
            return;
        }

        snprintf(tkw, 9, "TZERO%u", idx);
        if (get_header_double(header, tkw, &(column->tzero), 0.0, HDR_NOFLAG) == GET_HEADER_FAILED) {
            return;
        }

        snprintf(tkw, 9, "TNULL%u", idx);
        if (get_header_longlong(header, tkw, &(column->tnull), NULL_UNDEFINED, HDR_NOFLAG) == GET_HEADER_FAILED) {
            return;
        }
    }

    fileptr->Fptr->tableptr = *columns;
    fileptr->Fptr->tfield = tfields;

    // This routine from CFITSIO calculates the byte offset of each column
    // and stores it in the column->tbcol field
    ffgtbc(fileptr, &totalwidth, &status);
    if (status != 0) {
        process_status_err(status);
    }

    return;
}



void configure_compression(fitsfile* fileptr, PyObject* header) {
    /* Configure the compression-related elements in the fitsfile struct
       using values in the FITS header. */

    FITSfile* Fptr;

    int tfields;
    tcolumn* columns;

    char keyword[9];
    char zname[72];
    int znaxis;
    char tmp[72];
    float version;

    int idx;

    Fptr = fileptr->Fptr;
    tfields = Fptr->tfield;
    columns = Fptr->tableptr;

    int tmp_retval;

    // Get the ZBITPIX header value; if this is missing we're in trouble
    if (get_header_int(header, "ZBITPIX", &(Fptr->zbitpix), 0, HDR_FAIL_KEY_MISSING) != GET_HEADER_SUCCESS) {
        return;
    }

    // By default assume there is no ZBLANK column and check for ZBLANK or
    // BLANK in the header
    Fptr->cn_zblank = Fptr->cn_zzero = Fptr->cn_zscale = -1;
    Fptr->cn_uncompressed = 0;
#ifdef CFITSIO_SUPPORTS_GZIPDATA
    Fptr->cn_gzip_data = 0;
#endif

    // Check for a ZBLANK, ZZERO, ZSCALE, and
    // UNCOMPRESSED_DATA/GZIP_COMPRESSED_DATA columns in the compressed data
    // table
    for (idx = 0; idx < tfields; idx++) {
        if (0 == strncmp(columns[idx].ttype, "UNCOMPRESSED_DATA", 18)) {
            Fptr->cn_uncompressed = idx + 1;
#ifdef CFITSIO_SUPPORTS_GZIPDATA
        } else if (0 == strncmp(columns[idx].ttype,
                                "GZIP_COMPRESSED_DATA", 21)) {
            Fptr->cn_gzip_data = idx + 1;
#endif
        } else if (0 == strncmp(columns[idx].ttype, "ZSCALE", 7)) {
            Fptr->cn_zscale = idx + 1;
        } else if (0 == strncmp(columns[idx].ttype, "ZZERO", 6)) {
            Fptr->cn_zzero = idx + 1;
        } else if (0 == strncmp(columns[idx].ttype, "ZBLANK", 7)) {
            Fptr->cn_zblank = idx + 1;
        }
    }

    Fptr->zblank = 0;
    if (Fptr->cn_zblank < 1) {
        // No ZBLANK column--check the ZBLANK and BLANK heard keywords
        switch (get_header_int(header, "ZBLANK", &(Fptr->zblank), 0, HDR_NOFLAG)) {
          case GET_HEADER_FAILED:
            return;
          case GET_HEADER_DEFAULT_USED:
            // ZBLANK keyword not found
            if (get_header_int(header, "BLANK", &(Fptr->zblank), 0, HDR_NOFLAG) == GET_HEADER_FAILED) {
              return;
            }
            break;
          default:
            break;
        }
    }

    Fptr->zscale = 1.0;
    if (Fptr->cn_zscale < 1) {
        switch (get_header_double(header, "ZSCALE", &(Fptr->zscale), 1.0, HDR_NOFLAG)) {
          case GET_HEADER_FAILED:
            return;
          case GET_HEADER_DEFAULT_USED:
            Fptr->cn_zscale = 0;
            break;
          default:
            break;
        }
    }
    Fptr->cn_bscale = Fptr->zscale;

    Fptr->zzero = 0.0;
    if (Fptr->cn_zzero < 1) {
        switch (get_header_double(header, "ZZERO", &(Fptr->zzero), 0.0, HDR_NOFLAG)) {
          case GET_HEADER_FAILED:
            return;
          case GET_HEADER_DEFAULT_USED:
            Fptr->cn_zzero = 0;
            break;
          default:
            break;
        }
    }
    Fptr->cn_bzero = Fptr->zzero;

    if (get_header_string(header, "ZCMPTYPE", tmp, DEFAULT_COMPRESSION_TYPE, HDR_NOFLAG) == GET_HEADER_FAILED) {
        return;
    }
    strncpy(Fptr->zcmptype, tmp, 11);
    Fptr->zcmptype[strlen(tmp)] = '\0';

    Fptr->compress_type = compress_type_from_string(Fptr->zcmptype);
    if (PyErr_Occurred()) {
        return;
    }

    if (get_header_int(header, "ZNAXIS", &znaxis, 0, HDR_NOFLAG) == GET_HEADER_FAILED) {
        return;
    }
    Fptr->zndim = znaxis;

    if (znaxis > MAX_COMPRESS_DIM) {
        // The CFITSIO compression code currently only supports up to 6
        // dimensions by default.
        znaxis = MAX_COMPRESS_DIM;
    }

    Fptr->tilerow = NULL;
    Fptr->maxtilelen = 1;
    for (idx = 1; idx <= znaxis; idx++) {
        snprintf(keyword, 9, "ZNAXIS%u", idx);
        if (get_header_long(header, keyword, Fptr->znaxis + idx - 1, 0, HDR_NOFLAG) == GET_HEADER_FAILED) {
            return;
        }
        snprintf(keyword, 9, "ZTILE%u", idx);
        if (get_header_long(header, keyword, Fptr->tilesize + idx - 1, 0, HDR_NOFLAG) == GET_HEADER_FAILED) {
            return;
        }
        Fptr->maxtilelen *= Fptr->tilesize[idx - 1];
    }

    // Set some more default compression options
    Fptr->rice_blocksize = DEFAULT_BLOCK_SIZE;
    Fptr->rice_bytepix = DEFAULT_BYTE_PIX;
    Fptr->quantize_level = DEFAULT_QUANTIZE_LEVEL;
    Fptr->hcomp_smooth = DEFAULT_HCOMP_SMOOTH;
    Fptr->hcomp_scale = DEFAULT_HCOMP_SCALE;

    // Now process the ZVALn keywords
    idx = 1;
    while (1) {
        snprintf(keyword, 9, "ZNAME%u", idx);
        // Assumes there are no gaps in the ZNAMEn keywords; this same
        // assumption was made in the Python code.  This could be done slightly
        // more flexibly by using a wildcard slice of the header
        tmp_retval = get_header_string(header, keyword, zname, "", HDR_NOFLAG);
        if (tmp_retval == GET_HEADER_FAILED) {
            return;
        } else if (tmp_retval == 1) {
            break;
        }

        snprintf(keyword, 9, "ZVAL%u", idx);
        if (Fptr->compress_type == RICE_1) {
            if (0 == strcmp(zname, "BLOCKSIZE")) {
                if (get_header_int(header, keyword, &(Fptr->rice_blocksize),
                                   DEFAULT_BLOCK_SIZE, HDR_NOFLAG) == GET_HEADER_FAILED) {
                    return;
                }
            } else if (0 == strcmp(zname, "BYTEPIX")) {
                if (get_header_int(header, keyword, &(Fptr->rice_bytepix),
                                   DEFAULT_BYTE_PIX, HDR_NOFLAG) == GET_HEADER_FAILED) {
                    return;
                }
            }
        } else if (Fptr->compress_type == HCOMPRESS_1) {
            if (0 == strcmp(zname, "SMOOTH")) {
                if (get_header_int(header, keyword, &(Fptr->hcomp_smooth),
                                   DEFAULT_HCOMP_SMOOTH, HDR_NOFLAG) == GET_HEADER_FAILED) {
                    return;
                }
            } else if (0 == strcmp(zname, "SCALE")) {
                if (get_header_float(header, keyword, &(Fptr->hcomp_scale),
                                     DEFAULT_HCOMP_SCALE, HDR_NOFLAG) == GET_HEADER_FAILED) {
                    return;
                }
            }
        }
        if (Fptr->zbitpix < 0 && 0 == strcmp(zname, "NOISEBIT")) {
            if (get_header_float(header, keyword, &(Fptr->quantize_level),
                                 DEFAULT_QUANTIZE_LEVEL, HDR_NOFLAG) == GET_HEADER_FAILED) {
                return;
            }
            if (Fptr->quantize_level == 0.0) {
                /* NOISEBIT == 0 is equivalent to no quantize */
                Fptr->quantize_level = NO_QUANTIZE;
            }
        }

        idx++;
    }

    /* The ZQUANTIZ keyword determines the quantization algorithm; NO_QUANTIZE
       implies lossless compression */
    tmp_retval = get_header_string(header, "ZQUANTIZ", tmp, "", HDR_NOFLAG);
    if (tmp_retval == GET_HEADER_FAILED) {
        return;
    } else if (tmp_retval == GET_HEADER_SUCCESS) {
        /* Ugh; the fact that cfitsio defines its version as a float makes
           preprocessor comparison impossible */
        fits_get_version(&version);
        if ((version >= CFITSIO_LOSSLESS_COMP_SUPPORTED_VERS) &&
                (0 == strcmp(tmp, "NONE"))) {
            Fptr->quantize_level = NO_QUANTIZE;
        } else if (0 == strcmp(tmp, "SUBTRACTIVE_DITHER_1")) {
#ifdef CFITSIO_SUPPORTS_SUBTRACTIVE_DITHER_2
            // Added in CFITSIO 3.35, this also changed the name of the
            // quantize_dither struct member to quantize_method
            Fptr->quantize_method = SUBTRACTIVE_DITHER_1;
        } else if (0 == strcmp(tmp, "SUBTRACTIVE_DITHER_2")) {
            Fptr->quantize_method = SUBTRACTIVE_DITHER_2;
        } else {
            Fptr->quantize_method = NO_DITHER;
        }
    } else {
        Fptr->quantize_method = NO_DITHER;
    }

    if (Fptr->quantize_method != NO_DITHER) {
        switch (get_header_int(header, "ZDITHER0", &(Fptr->dither_seed), 0, HDR_NOFLAG)) {
          case GET_HEADER_FAILED:
            return;
          case GET_HEADER_DEFAULT_USED: // ZDITHER0 keyword not found
            Fptr->dither_seed = 0;
            Fptr->request_dither_seed = 0;
            break;
          default:
            break;
        }
    }
#else
            Fptr->quantize_dither = SUBTRACTIVE_DITHER_1;
        } else {
            Fptr->quantize_dither = NO_DITHER;
        }
    } else {
        Fptr->quantize_dither = NO_DITHER;
    }

    if (Fptr->quantize_dither != NO_DITHER) {
        switch (get_header_int(header, "ZDITHER0", &(Fptr->dither_offset), 0, HDR_NOFLAG)) {
          case GET_HEADER_FAILED:
            return;
          case GET_HEADER_DEFAULT_USED: // ZDITHER0 keyword no found
            /* TODO: Find out if that's actually working and not invalid... */
            Fptr->dither_offset = 0;
            Fptr->request_dither_offset = 0;
            break;
          default:
            break;
        }
    }
#endif

    Fptr->compressimg = 1;
    Fptr->maxelem = imcomp_calc_max_elem(Fptr->compress_type,
                                         Fptr->maxtilelen,
                                         Fptr->zbitpix,
                                         Fptr->rice_blocksize);
    Fptr->cn_compressed = 1;
    return;
}


void init_output_buffer(PyObject* hdu, void** buf, size_t* bufsize) {
    // Determines a good size for the output data buffer and allocates
    // memory for it, returning the address and size of the allocated
    // memory into **buf and *bufsize respectively.

    PyObject* header = NULL;
    char keyword[9];
    char tmp[72];
    int znaxis;
    int compress_type;
    int zbitpix;
    int rice_blocksize = 0;
    long long rowlen;
    long long nrows;
    long maxelem;
    long tilelen;
    unsigned long maxtilelen = 1;
    int idx;

    header = PyObject_GetAttrString(hdu, "_header");
    if (header == NULL) {
        return;
    }

    if (get_header_int(header, "ZNAXIS", &znaxis, 0,
                       HDR_FAIL_KEY_MISSING | HDR_FAIL_VAL_NEGATIVE) != GET_HEADER_SUCCESS) {
        goto fail;
    }

    if (znaxis > 999) {
        PyErr_SetString(PyExc_ValueError, "ZNAXIS is greater than 999.");
        goto fail;
    }

    for (idx = 1; idx <= znaxis; idx++) {
        snprintf(keyword, 9, "ZTILE%u", idx);
        if (get_header_long(header, keyword, &tilelen, 1, HDR_NOFLAG) == GET_HEADER_FAILED) {
            goto fail;
        }
        maxtilelen *= tilelen;
    }

    if (get_header_string(header, "ZCMPTYPE", tmp, DEFAULT_COMPRESSION_TYPE, HDR_NOFLAG) == GET_HEADER_FAILED) {
        goto fail;
    }
    compress_type = compress_type_from_string(tmp);
    if (PyErr_Occurred()) {
        goto fail;
    }
    if (compress_type == RICE_1) {
        if (get_header_int(header, "ZVAL1", &rice_blocksize, 0, HDR_NOFLAG) == GET_HEADER_FAILED) {
            goto fail;
        }
    }

    /* Because we calculate the size of the buffer based on these values they
       must not be negative. Otherwise it would wrap around during the casting
       to size_t and give huge values. */
    if (get_header_longlong(header, "NAXIS1", &rowlen, 0, HDR_FAIL_VAL_NEGATIVE) == GET_HEADER_FAILED) {
        goto fail;
    }
    if (get_header_longlong(header, "NAXIS2", &nrows, 0, HDR_FAIL_VAL_NEGATIVE) == GET_HEADER_FAILED) {
        goto fail;
    }

    // Get the ZBITPIX header value; if this is missing we're in trouble
    if (get_header_int(header, "ZBITPIX", &zbitpix, 0, HDR_FAIL_KEY_MISSING) != GET_HEADER_SUCCESS) {
        goto fail;
    }

    maxelem = imcomp_calc_max_elem(compress_type, maxtilelen, zbitpix,
                                   rice_blocksize);

    *bufsize = ((size_t) (rowlen * nrows) + (nrows * maxelem));

    if (*bufsize < IOBUFLEN) {
        // We must have a full FITS block at a minimum
        *bufsize = IOBUFLEN;
    } else if (*bufsize % IOBUFLEN != 0) {
        // Still make sure to pad out to a multiple of 2880 byte blocks
        // otherwise CFITSIO can get read errors when it tries to read
        // a partial block that goes past the end of the file
        *bufsize += ((size_t) (IOBUFLEN - (*bufsize % IOBUFLEN)));
    }

    *buf = calloc(*bufsize, sizeof(char));
    if (*buf == NULL) {
        // Checking if calloc failed.
        PyErr_SetString(PyExc_MemoryError,
                        "Failed to allocate memory for output data buffer.");
        goto fail;
    }

fail:
    Py_DECREF(header);
    return;
}


void get_hdu_data_base(PyObject* hdu, void** buf, size_t* bufsize) {
    // Given a pointer to an HDU object, returns a pointer to the deepest base
    // array of that HDU's data array into **buf, and the size of that array
    // into *bufsize.

    PyArrayObject* data = NULL;
    PyArrayObject* base;
    PyArrayObject* tmp;

    data = (PyArrayObject*) PyObject_GetAttrString(hdu, "compressed_data");
    if (data == NULL) {
        goto fail;
    }

    // Walk the array data bases until we find the lowest ndarray base; for
    // CompImageHDUs there should always be at least one contiguous byte array
    // allocated for the table and its heap
    if (!PyObject_TypeCheck(data, &PyArray_Type)) {
        PyErr_SetString(PyExc_TypeError,
                        "CompImageHDU.compressed_data must be a numpy.ndarray");
        goto fail;
    }

    tmp = base = data;
    while (PyObject_TypeCheck((PyObject*) tmp, &PyArray_Type)) {
        base = tmp;
        *bufsize = (size_t) PyArray_NBYTES(base);
        tmp = (PyArrayObject*) PyArray_BASE(base);
        if (tmp == NULL) {
            break;
        }
    }

    *buf = PyArray_DATA(base);
fail:
    Py_XDECREF(data);
    return;
}


void open_from_hdu(fitsfile** fileptr, void** buf, size_t* bufsize,
                   PyObject* hdu, tcolumn** columns, int mode) {

    PyObject* header = NULL;
    FITSfile* Fptr;

    int status = 0;
    long long rowlen;
    long long nrows;
    long long heapsize;
    long long theap;

    header = PyObject_GetAttrString(hdu, "_header");
    if (header == NULL) {
        goto fail;
    }

    if (get_header_longlong(header, "NAXIS1", &rowlen, 0, HDR_NOFLAG) == GET_HEADER_FAILED) {
        goto fail;
    }
    if (get_header_longlong(header, "NAXIS2", &nrows, 0, HDR_NOFLAG) == GET_HEADER_FAILED) {
        goto fail;
    }

    // The PCOUNT keyword contains the number of bytes in the table heap
    if (get_header_longlong(header, "PCOUNT", &heapsize, 0, HDR_FAIL_VAL_NEGATIVE) == GET_HEADER_FAILED) {
        goto fail;
    }

    // The THEAP keyword gives the offset of the heap from the beginning of
    // the HDU data portion; normally this offset is 0 but it can be set
    // to something else with THEAP
    if (get_header_longlong(header, "THEAP", &theap, 0, HDR_NOFLAG) == GET_HEADER_FAILED) {
        goto fail;
    }

    fits_create_memfile(fileptr, buf, bufsize, 0, realloc, &status);
    if (status != 0) {
        process_status_err(status);
        goto fail;
    }

    Fptr = (*fileptr)->Fptr;

    // Now we have some fun munging some of the elements in the fitsfile struct
    Fptr->writemode = mode;
    Fptr->open_count = 1;
    Fptr->hdutype = BINARY_TBL;  /* This is a binary table HDU */
    Fptr->lasthdu = 1;
    Fptr->headstart[0] = 0;
    Fptr->headend = 0;
    Fptr->datastart = 0;  /* There is no header, data starts at 0 */
    Fptr->origrows = Fptr->numrows = nrows;
    Fptr->rowlength = rowlen;
    if (theap != 0) {
        Fptr->heapstart = theap;
    } else {
        Fptr->heapstart = rowlen * nrows;
    }

    Fptr->heapsize = heapsize;

    // Configure the array of table column structs from the Astropy header
    // instead of allowing CFITSIO to try to read from the header
    tcolumns_from_header(*fileptr, header, columns);
    if (PyErr_Occurred()) {
        goto fail;
    }

    // If any errors occur in this function they'll bubble up from here to
    // compression_decompress_hdu
    configure_compression(*fileptr, header);

fail:
    Py_XDECREF(header);
    return;
}


PyObject* compression_compress_hdu(PyObject* self, PyObject* args)
{
    PyObject* hdu;
    PyObject* retval = NULL;
    tcolumn* columns = NULL;

    void* outbuf = NULL;
    size_t outbufsize;

    PyObject* tmp_indata;
    PyArrayObject* indata = NULL;
    PyArrayObject* tmp;
    npy_intp znaxis;
    int datatype;
    int npdatatype;
    unsigned long long heapsize;

    fitsfile* fileptr = NULL;
    FITSfile* Fptr = NULL;
    int status = 0;

    if (!PyArg_ParseTuple(args, "O:compression.compress_hdu", &hdu)) {
        return NULL;
    }

    // For HDU compression never use CFITSIO to write directly to the file;
    // although there's nothing wrong with CFITSIO, right now that would cause
    // too much confusion to Astropy's internal book keeping.
    // We just need to get the compressed bytes and Astropy will handle the
    // writing of them.
    init_output_buffer(hdu, &outbuf, &outbufsize);
    if (outbuf == NULL) {
        return NULL;
    }

    open_from_hdu(&fileptr, &outbuf, &outbufsize, hdu, &columns, READWRITE);
    if (PyErr_Occurred()) {
        goto fail;
    }

    Fptr = fileptr->Fptr;

    bitpix_to_datatypes(Fptr->zbitpix, &datatype, &npdatatype);
    if (PyErr_Occurred()) {
        goto fail;
    }

    /* The data attribute could be something different from an array, i.e. None */
    tmp_indata = PyObject_GetAttrString(hdu, "data");
    if (tmp_indata == NULL) {
        goto fail;
    }

    if (!PyObject_TypeCheck(tmp_indata, &PyArray_Type)) {
        PyErr_SetString(PyExc_TypeError,
                        "CompImageHDU.data must be a numpy.ndarray");
        Py_DECREF(tmp_indata);
        goto fail;
    }

    indata = (PyArrayObject*) tmp_indata;

    fits_write_img(fileptr, datatype, 1, PyArray_SIZE(indata),
                   PyArray_DATA(indata), &status);
    if (status != 0) {
        process_status_err(status);
        goto fail;
    }

    fits_flush_buffer(fileptr, 1, &status);
    if (status != 0) {
        process_status_err(status);
        goto fail;
    }

    // Previously this used outbufsize as the size to use for the new Numpy
    // byte array. However outbufsize is usually larger than necessary to
    // store all the compressed data exactly; instead use the exact size
    // of the compressed data from the heapsize plus the size of the table
    // itself
    heapsize = (unsigned long long) Fptr->heapsize;
    znaxis = (npy_intp) (Fptr->heapstart + heapsize);

    if (znaxis < outbufsize) {
        void* tmp_outbuf = NULL;
        // Go ahead and truncate to the size in znaxis to free the
        // redundant allocation
        if (znaxis == 0) {
            /* This really shouldn't happen, but if it did, we would have a
               problem because realloc would deallocate outbuf AND return NULL.
               */
            PyErr_SetString(PyExc_ValueError,
                            "Calculated array size is zero. This shouldn't happen!");
            goto fail;
        }
        tmp_outbuf = realloc(outbuf, (size_t) znaxis);
        if (tmp_outbuf == NULL) {
            PyErr_SetString(PyExc_MemoryError,
                            "Couldn't resize the output-buffer.");
            goto fail;
        }
        outbuf = tmp_outbuf;
    }

    tmp = (PyArrayObject*) PyArray_SimpleNewFromData(1, &znaxis, NPY_UBYTE,
                                                     outbuf);
    if (tmp == NULL) {
        /* Really not sure if it's always safe to free outbuf when
           PyArray_SimpleNewFromData failed (which is unlikely but could happen)
           but it seems like if it fails then the outbuf NEEDS to be freed... */
        goto fail;
    }
    PyArray_ENABLEFLAGS(tmp, NPY_ARRAY_OWNDATA);
    /* From this point on outbuf MUST NOT BE FREED! */

    // Leaves refcount of tmp untouched, so its refcount should remain as 1
    retval = Py_BuildValue("KN", heapsize, tmp);
    if (retval == NULL) {
        Py_DECREF(tmp);
        goto cleanup;
    }

    goto cleanup;

fail:
    if (outbuf != NULL) {
        // At this point outbuf should never not be NULL, but in principle
        // buggy code somewhere in CFITSIO or Numpy could set it to NULL
        free(outbuf);
    }
cleanup:
    if (columns != NULL) {
        free(columns);
        /* See https://github.com/astropy/astropy/pull/4489
           We can only set the tableptr to NULL if Fptr is actually not NULL.
           */
        if (fileptr != NULL && fileptr->Fptr != NULL) {
            fileptr->Fptr->tableptr = NULL;
        }
    }

    if (fileptr != NULL) {
        status = 1; // Disable header-related errors
        fits_close_file(fileptr, &status);
        if (status != 1) {
            process_status_err(status);
            retval = NULL;
        }
    }

    Py_XDECREF(indata);

    // Clear any messages remaining in CFITSIO's error stack
    fits_clear_errmsg();

    return retval;
}


PyObject* compression_decompress_hdu(PyObject* self, PyObject* args)
{

    PyObject* hdu;
    tcolumn* columns = NULL;

    void* inbuf;
    size_t inbufsize;

    PyArrayObject* outdata = NULL;
    int datatype;
    int npdatatype;
    npy_intp zndim;
    npy_intp* znaxis = NULL;
    long arrsize;

    fitsfile* fileptr = NULL;
    int anynul = 0;
    int status = 0;
    int idx;

    int free_columns_manually = 1;

    if (!PyArg_ParseTuple(args, "O:compression.decompress_hdu", &hdu)) {
        return NULL;
    }

    // Grab a pointer to the input data from the HDU's compressed_data
    // attribute
    get_hdu_data_base(hdu, &inbuf, &inbufsize);
    if (PyErr_Occurred()) {
        return NULL;
    } else if (inbufsize == 0) {
        // The compressed data buffer is empty (probably zero rows, for an
        // empty "compressed" image.  Just return None in this case.
        Py_RETURN_NONE;
    }

    open_from_hdu(&fileptr, &inbuf, &inbufsize, hdu, &columns, READONLY);
    if (PyErr_Occurred()) {
        goto fail;
    }

    bitpix_to_datatypes(fileptr->Fptr->zbitpix, &datatype, &npdatatype);
    if (PyErr_Occurred()) {
        goto fail;
    }

    zndim = (npy_intp)fileptr->Fptr->zndim;
    znaxis = PyMem_Malloc(sizeof(npy_intp) * zndim);
    if (znaxis == NULL) {
        goto fail;
    }

    arrsize = 1;
    for (idx = 0; idx < zndim; idx++) {
        znaxis[zndim - idx - 1] = fileptr->Fptr->znaxis[idx];
        arrsize *= fileptr->Fptr->znaxis[idx];
    }

    /* Create and allocate a new array for the decompressed data */
    outdata = (PyArrayObject*) PyArray_SimpleNew(zndim, znaxis, npdatatype);
    if (outdata == NULL) {
        goto fail;
    }

    fits_read_img(fileptr, datatype, 1, arrsize, NULL, PyArray_DATA(outdata),
                  &anynul, &status);
    /* At this point we need to let CFITSIO clean up the tableptr and the
       compressed tile cache. */
    free_columns_manually = 0;
    if (status != 0) {
        process_status_err(status);
        Py_DECREF(outdata);
        outdata = NULL;
    }

fail:
    // CFITSIO will free this object in the ffchdu function by way of
    // fits_close_file; we need to let CFITSIO handle this so that it also
    // cleans up the compressed tile cache - but that's only necessary in case
    // we called "fits_read_img"...
    if (free_columns_manually && columns != NULL) {
        free(columns);
        if (fileptr != NULL && fileptr->Fptr != NULL) {
            fileptr->Fptr->tableptr = NULL;
        }
    }

    if (fileptr != NULL) {
        status = 1;// Disable header-related errors
        fits_close_file(fileptr, &status);
        if (status != 1) {
            process_status_err(status);
            outdata = NULL;
        }
    }

    if (znaxis != NULL) {
        PyMem_Free(znaxis);
    }

    // Clear any messages remaining in CFITSIO's error stack
    fits_clear_errmsg();

    return (PyObject*) outdata;
}


/* CFITSIO version float as returned by fits_get_version() */
static double cfitsio_version;


int compression_module_init(PyObject* module) {
    /* Python version-independent initialization routine for the
       compression module. Returns 0 on success and -1 (with exception set)
       on failure. */
    PyObject* tmp;
    float version_tmp;
    int ret;

    fits_get_version(&version_tmp);
    cfitsio_version = (double) version_tmp;
    /* The conversion to double can lead to some rounding errors; round to the
       nearest 3 decimal places, which should be accurate for any past or
       current CFITSIO version. This is why relying on floats for version
       comparison isn't generally a bright idea... */
    cfitsio_version = floor((1000 * version_tmp + 0.5)) / 1000;

    tmp = PyFloat_FromDouble(cfitsio_version);
    if (tmp == NULL) {
        return -1;
    }
    ret = PyObject_SetAttrString(module, "CFITSIO_VERSION", tmp);
    Py_DECREF(tmp);
    return ret;
}


/* Method table mapping names to wrappers */
static PyMethodDef compression_methods[] =
{
   {"compress_hdu", compression_compress_hdu, METH_VARARGS},
   {"decompress_hdu", compression_decompress_hdu, METH_VARARGS},
   {NULL, NULL}
};

static struct PyModuleDef compressionmodule = {
    PyModuleDef_HEAD_INIT,
    "compression",
    "astropy.compression module",
    -1, /* No global state */
    compression_methods
};

PyObject *
PyInit_compression(void)
{
    PyObject* module = PyModule_Create(&compressionmodule);
    if (module == NULL) {
        return NULL;
    }
    if (compression_module_init(module)) {
        Py_DECREF(module);
        return NULL;
    }

    /* Needed to use Numpy routines */
    /* Note -- import_array() is a macro that behaves differently in Python2.x
     * vs. Python 3. See the discussion at:
     * https://groups.google.com/d/topic/astropy-dev/6_AesAsCauM/discussion
     */
    import_array();
    return module;
}
