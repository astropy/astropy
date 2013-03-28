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
/* data and returns the compressed data for all tiles into the .compData     */
/* attribute of that HDU.                                                    */
/*                                                                           */
/* The second function is decompress_hdu.  It takes an                       */
/* astropy.io.fits.CompImageHDU object that already has compressed data in   */
/* its .compData attribute. It returns the decompressed image data into the  */
/* HDU's .data attribute.                                                    */
/*                                                                           */
/* Finally, the utility function calc_max_elm is used internally by the      */
/* CompImageHDU class to estimate how much memory to allocate for the        */
/* compressed data when compressing an image.                                */
/*                                                                           */
/* Copyright (C) 2012 Association of Universities for Research in Astronomy  */
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

#include <math.h>

#include <Python.h>
#include <numpy/arrayobject.h>
#include <fitsio2.h>
#include <string.h>
#include "compressionmodule.h"

/* Some defines for Python3 support--bytes objects should be used where */
/* strings were previously used                                         */
#if PY_MAJOR_VERSION >= 3
#define IS_PY3K
#endif

#ifdef IS_PY3K
#define PyString_FromString PyUnicode_FromString
#define PyInt_AsLong PyLong_AsLong
#endif


/* These defaults mirror the defaults in astropy.io.fits.hdu.compressed */
#define DEFAULT_COMPRESSION_TYPE "RICE_1"
#define DEFAULT_QUANTIZE_LEVEL 16.0
#define DEFAULT_HCOMP_SCALE 0
#define DEFAULT_HCOMP_SMOOTH 0
#define DEFAULT_BLOCK_SIZE 32
#define DEFAULT_BYTE_PIX 4


/* Report any error based on the status returned from cfitsio. */
void process_status_err(int status)
{
   PyObject* except_type;
   char      err_msg[81];
   char      def_err_msg[81];

   err_msg[0] = '\0';
   def_err_msg[0] = '\0';

   switch (status)
   {
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
   }

   if (fits_read_errmsg(err_msg))
   {
      PyErr_SetString(except_type, err_msg);
   }
   else if (*def_err_msg)
   {
      PyErr_SetString(except_type, def_err_msg);
   }
   else
   {
      PyErr_SetString(except_type, "unknown error.");
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
            *npdatatype = NPY_INT8;
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
   }

   return;
}



int compress_type_from_string(char* zcmptype) {
    if (0 == strcmp(zcmptype, "RICE_1")) {
        return RICE_1;
    } else if (0 == strcmp(zcmptype, "GZIP_1")) {
        return GZIP_1;
    } else if (0 == strcmp(zcmptype, "PLIO_1")) {
        return PLIO_1;
    } else if (0 == strcmp(zcmptype, "HCOMPRESS_1")) {
        return HCOMPRESS_1;
    } else {
        PyErr_Format(PyExc_ValueError, "Unrecognized compression type: %s",
                     zcmptype);
        return -1;
    }
}


// TODO: It might be possible to simplify these further by making the
// conversion function (eg. PyString_AsString) an argument to a macro or
// something, but I'm not sure yet how easy it is to generalize the error
// handling
int get_header_string(PyObject* header, char* keyword, char** val, char* def) {
    PyObject* keystr;
    PyObject* keyval;
#ifdef IS_PY3K
    PyObject* tmp;  // Temp pointer to decoded bytes object
#endif
    int retval;

    keystr = PyString_FromString(keyword);
    keyval = PyObject_GetItem(header, keystr);

    if (keyval != NULL) {
#ifdef IS_PY3K
        // FITS header values should always be ASCII, but Latin1 is on the
        // safe side
        tmp = PyUnicode_AsLatin1String(keyval);
        *val = PyBytes_AsString(tmp);
        Py_DECREF(tmp);
#else
        *val = PyString_AsString(keyval);
#endif
        retval = 0;
    }
    else {
        PyErr_Clear();
        *val = def;
        retval = 1;
    }

    Py_DECREF(keystr);
    Py_XDECREF(keyval);
    return retval;
}


int get_header_int(PyObject* header, char* keyword, int* val, int def) {
    PyObject* keystr;
    PyObject* keyval;
    int retval;

    keystr = PyString_FromString(keyword);
    keyval = PyObject_GetItem(header, keystr);

    if (keyval != NULL) {
        *val = (int) PyInt_AsLong(keyval);
        retval = 0;
    }
    else {
        PyErr_Clear();
        *val = def;
        retval = 1;
    }

    Py_DECREF(keystr);
    Py_XDECREF(keyval);
    return retval;
}


int get_header_long(PyObject* header, char* keyword, long* val, long def) {
    PyObject* keystr;
    PyObject* keyval;
    int retval;

    keystr = PyString_FromString(keyword);
    keyval = PyObject_GetItem(header, keystr);

    if (keyval != NULL) {
        *val = PyLong_AsLong(keyval);
        retval = 0;
    }
    else {
        PyErr_Clear();
        *val = def;
        retval = 1;
    }

    Py_DECREF(keystr);
    Py_XDECREF(keyval);
    return retval;
}


int get_header_float(PyObject* header, char* keyword, float* val,
                     float def) {
    PyObject* keystr;
    PyObject* keyval;
    int retval;

    keystr = PyString_FromString(keyword);
    keyval = PyObject_GetItem(header, keystr);

    if (keyval != NULL) {
        *val = (float) PyFloat_AsDouble(keyval);
        retval = 0;
    }
    else {
        PyErr_Clear();
        *val = def;
        retval = 1;
    }

    Py_DECREF(keystr);
    Py_XDECREF(keyval);
    return retval;
}


int get_header_double(PyObject* header, char* keyword, double* val,
                      double def) {
    PyObject* keystr;
    PyObject* keyval;
    int retval;

    keystr = PyString_FromString(keyword);
    keyval = PyObject_GetItem(header, keystr);

    if (keyval != NULL) {
        *val = PyFloat_AsDouble(keyval);
        retval = 0;
    }
    else {
        PyErr_Clear();
        *val = def;
        retval = 1;
    }

    Py_DECREF(keystr);
    Py_XDECREF(keyval);
    return retval;
}


int get_header_longlong(PyObject* header, char* keyword, long long* val,
                        long long def) {
    PyObject* keystr;
    PyObject* keyval;
    int retval;

    keystr = PyString_FromString(keyword);
    keyval = PyObject_GetItem(header, keystr);

    if (keyval != NULL) {
        *val = PyLong_AsLongLong(keyval);
        retval = 0;
    }
    else {
        PyErr_Clear();
        *val = def;
        retval = 1;
    }

    Py_DECREF(keystr);
    Py_XDECREF(keyval);
    return retval;
}


void* compression_realloc(void* ptr, size_t size) {
    // This realloc()-like function actually just mallocs the requested
    // size and copies from the original memory address into the new one, and
    // returns the newly malloc'd address.
    // This is generally less efficient than an actual realloc(), but the
    // problem with using realloc in this case is that when it succeeds it will
    // free() the original memory, which may still be in use by the ndarray
    // using that memory as its data buffer.  This seems like the least hacky
    // way around that for now.
    // I'm open to other ideas though.
    void* tmp;
    tmp = malloc(size);
    memcpy(tmp, ptr, size);
    return tmp;
}


void tcolumns_from_header(fitsfile* fileptr, PyObject* header,
                          tcolumn** columns) {
    // Creates the array of tcolumn structures from the table column keywords
    // read from the astropy.io.fits.Header object; caller is responsible for
    // freeing the memory allocated for this array

    tcolumn* column;
    char tkw[9];
    unsigned int idx;

    int tfields;
    char* ttype;
    char* tform;
    int dtcode;
    long trepeat;
    long twidth;
    long long totalwidth;
    int status = 0;

    get_header_int(header, "TFIELDS", &tfields, 0);

    *columns = column = PyMem_New(tcolumn, (size_t) tfields);
    if (column == NULL) {
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
        get_header_string(header, tkw, &ttype, "");
        strncpy(column->ttype, ttype, 69);
        column->ttype[69] = '\0';

        snprintf(tkw, 9, "TFORM%u", idx);
        get_header_string(header, tkw, &tform, "");
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
        get_header_double(header, tkw, &(column->tscale), 1.0);

        snprintf(tkw, 9, "TZERO%u", idx);
        get_header_double(header, tkw, &(column->tzero), 0.0);

        snprintf(tkw, 9, "TNULL%u", idx);
        get_header_longlong(header, tkw, &(column->tnull), NULL_UNDEFINED);
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
    char* zname;
    int znaxis;
    char* tmp;
    float version;

    unsigned int idx;

    Fptr = fileptr->Fptr;
    tfields = Fptr->tfield;
    columns = Fptr->tableptr;

    // Get the ZBITPIX header value; if this is missing we're in trouble
    if (0 != get_header_int(header, "ZBITPIX", &(Fptr->zbitpix), 0)) {
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
        if(0 != get_header_int(header, "ZBLANK", &(Fptr->zblank), 0)) {
            // ZBLANK keyword not found
            get_header_int(header, "BLANK", &(Fptr->zblank), 0);
        }
    }

    Fptr->zscale = 1.0;
    if (Fptr->cn_zscale < 1) {
        if (0 != get_header_double(header, "ZSCALE", &(Fptr->zscale), 1.0)) {
            Fptr->cn_zscale = 0;
        }
    }
    Fptr->cn_bscale = Fptr->zscale;

    Fptr->zzero = 0.0;
    if (Fptr->cn_zzero < 1) {
        if (0 != get_header_double(header, "ZZERO", &(Fptr->zzero), 0.0)) {
            Fptr->cn_zzero = 0;
        }
    }
    Fptr->cn_bzero = Fptr->zzero;

    get_header_string(header, "ZCMPTYPE", &tmp, DEFAULT_COMPRESSION_TYPE);
    strncpy(Fptr->zcmptype, tmp, 11);
    Fptr->zcmptype[strlen(tmp)] = '\0';

    Fptr->compress_type = compress_type_from_string(Fptr->zcmptype);
    if (PyErr_Occurred()) {
        return;
    }

    get_header_int(header, "ZNAXIS", &znaxis, 0);
    Fptr->zndim = znaxis;

    if (znaxis > MAX_COMPRESS_DIM) {
        // The CFITSIO compression code currently only supports up to 6
        // dimensions by default.
        znaxis = MAX_COMPRESS_DIM;
    }

    Fptr->maxtilelen = 1;
    for (idx = 1; idx <= znaxis; idx++) {
        snprintf(keyword, 9, "ZNAXIS%u", idx);
        get_header_long(header, keyword, Fptr->znaxis + idx - 1, 0);
        snprintf(keyword, 9, "ZTILE%u", idx);
        get_header_long(header, keyword, Fptr->tilesize + idx - 1, 0);
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
        if (0 != get_header_string(header, keyword, &zname, "")) {
            break;
        }
        snprintf(keyword, 9, "ZVAL%u", idx);
        if (Fptr->compress_type == RICE_1) {
            if (0 == strcmp(zname, "BLOCKSIZE")) {
                get_header_int(header, keyword, &(Fptr->rice_blocksize),
                               DEFAULT_BLOCK_SIZE);
            } else if (0 == strcmp(zname, "BYTEPIX")) {
                get_header_int(header, keyword, &(Fptr->rice_bytepix),
                               DEFAULT_BYTE_PIX);
            }
        } else if (Fptr->compress_type == HCOMPRESS_1) {
            if (0 == strcmp(zname, "SMOOTH")) {
                get_header_int(header, keyword, &(Fptr->hcomp_smooth),
                               DEFAULT_HCOMP_SMOOTH);
            } else if (0 == strcmp(zname, "SCALE")) {
                get_header_float(header, keyword, &(Fptr->hcomp_scale),
                                 DEFAULT_HCOMP_SCALE);
            }
        } else if (Fptr->zbitpix < 0 && 0 == strcmp(zname, "NOISEBIT")) {
             get_header_float(header, keyword, &(Fptr->quantize_level),
                              DEFAULT_QUANTIZE_LEVEL);
             if (Fptr->quantize_level == 0.0) {
                 /* NOISEBIT == 0 is equivalent to no quantize */
                 Fptr->quantize_level = NO_QUANTIZE;
             }
        }

        idx++;
    }

    /* The ZQUANTIZ keyword determines the quantization algorithm; NO_QUANTIZE
       implies lossless compression */
    if (0 == get_header_string(header, "ZQUANTIZ", &tmp, "")) {
        /* Ugh; the fact that cfitsio defines its version as a float makes
           preprocessor comparison impossible */
        fits_get_version(&version);
        if ((version >= CFITSIO_LOSSLESS_COMP_SUPPORTED_VERS) &&
                (0 == strcmp(tmp, "NONE"))) {
            Fptr->quantize_level = NO_QUANTIZE;
        } else if (0 == strcmp(tmp, "SUBTRACTIVE_DITHER_1")) {
            Fptr->quantize_dither = SUBTRACTIVE_DITHER_1;
        } else {
            Fptr->quantize_dither = 0;
        }
    } else {
        Fptr->quantize_dither = 0;
    }

    Fptr->compressimg = 1;
    Fptr->maxelem = imcomp_calc_max_elem(Fptr->compress_type,
                                         Fptr->maxtilelen,
                                         Fptr->zbitpix,
                                         Fptr->rice_blocksize);
    Fptr->cn_compressed = 1;
    return;
}


void open_from_hdu(fitsfile** fileptr, void** buf, size_t* bufsize,
                   PyObject* hdu, tcolumn* columns) {
    PyObject* header = NULL;
    PyArrayObject* data = NULL;
    PyArrayObject* base;
    PyArrayObject* tmp;
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

    data = (PyArrayObject*) PyObject_GetAttrString(hdu, "compData");
    if (data == NULL) {
        goto fail;
    }


    get_header_longlong(header, "NAXIS1", &rowlen, 0);
    get_header_longlong(header, "NAXIS2", &nrows, 0);

    // The PCOUNT keyword contains the number of bytes in the table heap
    get_header_longlong(header, "PCOUNT", &heapsize, 0);

    // The THEAP keyword gives the offset of the heap from the beginning of
    // the HDU data portion; normally this offset is 0 but it can be set
    // to something else with THEAP
    get_header_longlong(header, "THEAP", &theap, 0);


    // Walk the array data bases until we find the lowest ndarray base; for
    // CompImageHDUs there should always be at least one contiguous byte array
    // allocated for the table and its heap
    if (!PyObject_TypeCheck(data, &PyArray_Type)) {
        PyErr_SetString(PyExc_TypeError,
                        "CompImageHDU.compData must be a numpy.ndarray");
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

    fits_create_memfile(fileptr, buf, bufsize, 0, compression_realloc,
                        &status);
    if (status != 0) {
        process_status_err(status);
        goto fail;
    }

    Fptr = (*fileptr)->Fptr;

    // Now we have some fun munging some of the elements in the fitsfile struct
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
    tcolumns_from_header(*fileptr, header, &columns);
    if (PyErr_Occurred()) {
        goto fail;
    }

    // If any errors occur in this function they'll bubble up from here to
    // compression_decompress_hdu
    configure_compression(*fileptr, header);

fail:
    Py_XDECREF(header);
    Py_XDECREF(data);
    return;
}



PyObject* compression_compress_hdu(PyObject* self, PyObject* args)
{
    PyObject* hdu;
    PyObject* retval = NULL;
    tcolumn* columns = NULL;

    void* orig_outbuf;
    size_t orig_outbufsize;
    void* outbuf;
    size_t outbufsize;

    PyArrayObject* indata;
    PyArrayObject* tmp;
    long znaxis;
    int datatype;
    int npdatatype;
    unsigned long long heapsize;

    fitsfile* fileptr;
    int status = 0;

    if (!PyArg_ParseTuple(args, "O:compression.compress_hdu", &hdu))
    {
        PyErr_SetString(PyExc_TypeError, "Couldn't parse arguments");
        return NULL;
    }

    // For HDU compression never use CFITSIO to write directly to the file;
    // although there's nothing wrong with CFITSIO, right now that would cause
    // too much confusion to Astropy's internal book keeping.
    // We just need to get the compressed bytes and Astropy will handle the
    // writing of them.
    open_from_hdu(&fileptr, &outbuf, &outbufsize, hdu, columns);
    if (PyErr_Occurred()) {
        return NULL;
    }
    orig_outbuf = outbuf;
    orig_outbufsize = outbufsize;

    bitpix_to_datatypes(fileptr->Fptr->zbitpix, &datatype, &npdatatype);
    if (PyErr_Occurred()) {
        return NULL;
    }

    indata = (PyArrayObject*) PyObject_GetAttrString(hdu, "data");

    fits_write_img(fileptr, datatype, 1, PyArray_SIZE(indata), indata->data,
                   &status);
    if (status != 0) {
        process_status_err(status);
        goto fail;
    }

    fits_flush_buffer(fileptr, 1, &status);
    if (status != 0) {
        process_status_err(status);
        goto fail;
    }

    // It's possible that between calls to fits_write_img and fits_flush_buffer
    // the size of the output buffer was insufficient and the buffer had to be
    // reallocated.  The address and size of the new buffer should be in
    // outbuf and outbufsize.  We need to create a new PyArrayObject using the
    // new buffer and size
    if (orig_outbuf != outbuf || orig_outbufsize != outbufsize) {
        // It's possible, albeit unlikely, that realloc can return a block of
        // memory with the same address but different size.
        znaxis = (long) outbufsize;  // The output array is just one dimension.
        tmp = (PyArrayObject*) PyArray_SimpleNewFromData(1, &znaxis, NPY_UBYTE,
                                                         outbuf);
        PyObject_SetAttrString(hdu, "compData", (PyObject*) tmp);
        Py_DECREF(tmp);
    }

    heapsize = (unsigned long long) fileptr->Fptr->heapsize;

    retval = PyLong_FromUnsignedLongLong(heapsize);

fail:
    if (fileptr != NULL) {
        status = 1; // Disable header-related errors
        fits_close_file(fileptr, &status);
        if (status != 1) {
            process_status_err(status);
            retval = NULL;
        }
    }

    if (columns != NULL) {
        PyMem_Free(columns);
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

    PyArrayObject* outdata;
    int datatype;
    int npdatatype;
    int zndim;
    long* znaxis;
    long arrsize;
    unsigned int idx;

    fitsfile* fileptr;
    int anynul = 0;
    int status = 0;

    if (!PyArg_ParseTuple(args, "O:compression.decompress_hdu", &hdu))
    {
        PyErr_SetString(PyExc_TypeError, "Couldn't parse arguments");
        return NULL;
    }

    open_from_hdu(&fileptr, &inbuf, &inbufsize, hdu, columns);
    if (PyErr_Occurred()) {
        return NULL;
    }

    bitpix_to_datatypes(fileptr->Fptr->zbitpix, &datatype, &npdatatype);
    if (PyErr_Occurred()) {
        return NULL;
    }

    zndim = fileptr->Fptr->zndim;
    znaxis = (long*) PyMem_Malloc(sizeof(long) * zndim);
    arrsize = 1;
    for (idx = 0; idx < zndim; idx++) {
        znaxis[zndim - idx - 1] = fileptr->Fptr->znaxis[idx];
        arrsize *= fileptr->Fptr->znaxis[idx];
    }

    /* Create and allocate a new array for the decompressed data */
    outdata = (PyArrayObject*) PyArray_SimpleNew(zndim, znaxis, npdatatype);

    fits_read_img(fileptr, datatype, 1, arrsize, NULL, outdata->data, &anynul,
                  &status);
    if (status != 0) {
        process_status_err(status);
        outdata = NULL;
        goto fail;
    }

fail:
    if (fileptr != NULL) {
        status = 1;// Disable header-related errors
        fits_close_file(fileptr, &status);
        if (status != 1) {
            process_status_err(status);
            outdata = NULL;
        }
    }

    if (columns != NULL) {
        PyMem_Free(columns);
    }
    PyMem_Free(znaxis);

    // Clear any messages remaining in CFITSIO's error stack
    fits_clear_errmsg();

    return (PyObject*) outdata;
}


PyObject* compression_calc_max_elem(PyObject* self, PyObject* args) {
    char* zcmptype;
    long maxtilelen;
    int zbitpix;
    int rice_blocksize;

    if (!PyArg_ParseTuple(args, "slii:compression.calc_max_elem", &zcmptype,
                          &maxtilelen, &zbitpix, &rice_blocksize)) {
        PyErr_SetString(PyExc_TypeError, "Couldn't parse arguments");
        return NULL;
    }

    return PyLong_FromLong(
        (long) imcomp_calc_max_elem(compress_type_from_string(zcmptype),
                                    maxtilelen, zbitpix, rice_blocksize));
}


/* CFITSIO version float as returned by fits_get_version() */
static double cfitsio_version;


void compression_module_init(PyObject* module) {
    /* Python version-indendependent initialization routine for the
       compression module */
    PyObject* tmp;
    float version_tmp;

    fits_get_version(&version_tmp);
    cfitsio_version = (double) version_tmp;
    /* The conversion to double can lead to some rounding errors; round to the
       nearest 3 decimal places, which should be accurate for any past or
       current CFITSIO version. This is why relying on floats for version
       comparison isn't generally a bright idea... */
    cfitsio_version = floor((1000 * version_tmp + 0.5)) / 1000;

    tmp = PyFloat_FromDouble(cfitsio_version);
    PyObject_SetAttrString(module, "CFITSIO_VERSION", tmp);
    Py_XDECREF(tmp);

    return;
}


/* Method table mapping names to wrappers */
static PyMethodDef compression_methods[] =
{
   {"compress_hdu", compression_compress_hdu, METH_VARARGS},
   {"decompress_hdu", compression_decompress_hdu, METH_VARARGS},
   {"calc_max_elem", compression_calc_max_elem, METH_VARARGS},
   {NULL, NULL}
};

#ifdef IS_PY3K
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
    compression_module_init(module);

    /* Needed to use Numpy routines */
    /* Note -- import_array() is a macro that behaves differently in Python2.x
     * vs. Python 3. See the discussion at:
     * https://groups.google.com/d/topic/astropy-dev/6_AesAsCauM/discussion
     */
    import_array();
    return module;
}
#else
PyMODINIT_FUNC initcompression(void)
{
   PyObject* module = Py_InitModule4("compression", compression_methods,
                                     "astropy.io.fits.compression module",
                                     NULL, PYTHON_API_VERSION);
   compression_module_init(module);
   import_array();
}
#endif
