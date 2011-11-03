/* $Id$
*/

/*****************************************************************************/
/*                                                                           */
/* This file, compress.c, contains the code required to compress and         */
/* uncompress data using the GZIP_1 compression format.                      */
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
/* This file was copied and heavily modified from the FITSIO software that   */
/* was written by William Pence at the High Energy Astrophysic Science       */
/* Archive Research Center (HEASARC) at the NASA Goddard Space Flight Center.*/
/* That software contained the following copyright and warranty notices:     */
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
/* This code calls routines from and links to the ZLIB compression library   */
/* that was written by Jean-loup Gailly and Mark Adler.  This package is     */
/* normally destributed with Python 2.5.  That software containes the        */
/* following copyright and warranty notices:     */
/*                                                                           */
/*  Copyright (C) 1995-2010 Jean-loup Gailly and Mark Adler                  */
/*                                                                           */
/*  This software is provided 'as-is', without any express or implied        */
/*  warranty.  In no event will the authors be held liable for any damages   */
/*  arising from the use of this software.                                   */
/*                                                                           */
/*  Permission is granted to anyone to use this software for any purpose,    */
/*  including commercial applications, and to alter it and redistribute it   */
/*  freely, subject to the following restrictions:                           */
/*                                                                           */
/*  1. The origin of this software must not be misrepresented; you must not  */
/*     claim that you wrote the original software. If you use this software  */
/*     in a product, an acknowledgment in the product documentation would be */
/*     appreciated but is not required.                                      */
/*  2. Altered source versions must be plainly marked as such, and must not  */
/*     be misrepresented as being the original software.                     */
/*  3. This notice may not be removed or altered from any source             */
/*     distribution.                                                         */
/*                                                                           */
/*  Jean-loup Gailly        Mark Adler                                       */
/*  jloup@gzip.org          madler@alumni.caltech.edu                        */
/*                                                                           */
/*                                                                           */
/*  The data format used by the zlib library is described by RFCs (Request   */
/*  for Comments) 1950 to 1952 in the files                                  */
/*  http://www.ietf.org/rfc/rfc1950.txt (zlib format),                       */
/*  rfc1951.txt (deflate format) and rfc1952.txt (gzip format).              */
/*                                                                           */
/*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "zlib.h"  

int _pyfits_uncompress2mem_from_mem(
             char *inmemptr,     
             size_t inmemsize, 
             char **buffptr,  
             size_t *buffsize,  
             void *(*mem_realloc)(void *p, size_t newsize), 
             size_t *filesize,  
             int *status);

int _pyfits_compress2mem_from_mem(
             char *inmemptr,     
             size_t inmemsize, 
             char **buffptr,  
             size_t *buffsize,  
             void *(*mem_realloc)(void *p, size_t newsize), 
             size_t *filesize,  
             int *status);

/*--------------------------------------------------------------------------*/
int _pyfits_uncompress2mem_from_mem(
             char *inmemptr,     /* I - memory pointer to compressed bytes */
             size_t inmemsize,   /* I - size of input compressed file      */
             char **buffptr,   /* IO - memory pointer                      */
             size_t *buffsize,   /* IO - size of buffer, in bytes           */
             void *(*mem_realloc)(void *p, size_t newsize), /* function     */
             size_t *filesize,   /* O - size of file, in bytes              */
             int *status)        /* IO - error status                       */

/*
  Uncompress the file into memory.  Fill whatever amount of memory has
  already been allocated, then realloc more memory, using the supplied
  input function, if necessary.
*/
{
    int err; 
    uLong uncomprLen;
    Byte *uncompr;
    z_stream d_stream;   /* decompression stream */
    uLong bytes_out_so_far = 0;  /* Keeps track of the number of bytes put in
                                    the output buffer so far */


    if (*status > 0) 
        return(*status); 

    /* Allocate memory as a temporary buffer in which to uncompress. */
    uncomprLen = *buffsize;
    uncompr = (Byte*)malloc(*buffsize);

    d_stream.zalloc = (alloc_func)0;
    d_stream.zfree = (free_func)0;
    d_stream.opaque = (voidpf)0;

    d_stream.next_in = (unsigned char*)inmemptr;
    d_stream.avail_in = inmemsize;

    /* Initialize the decompression.  The argument (15+16) tells the
       decompressor that we are to use the gzip algorithm */
    err = inflateInit2(&d_stream, (15+16));

    if (err != Z_OK)
    {
        /* free temporary output data buffer */
        free(uncompr);
        return(*status = 414);
    }

    for (;;)
    {
        /* Output to the temporary buffer.  This will overwrite the
           previous data each time. */
        d_stream.next_out = uncompr;
        d_stream.avail_out = uncomprLen;

        err = _pyfits_inflate(&d_stream, Z_NO_FLUSH);

        if (err != Z_OK && err != Z_STREAM_END)
        {
            /* free temporary output data buffer */
            free(uncompr);
            return(*status = 414);
        }

        if (d_stream.total_out > *buffsize)
        {
            /* OK, we need more memory for the output so reallocate it */
            *buffsize = d_stream.total_out;
            *buffptr = mem_realloc(*buffptr,*buffsize);

            if (*buffptr == NULL)
            {
                /* free temporary output data buffer */
                free(uncompr);
                return(*status = 414);
            }
        }

        /* copy from the temporary buffer into the output memory buffer */
        memcpy((char *) *buffptr + bytes_out_so_far, (char *) uncompr,
               d_stream.total_out-bytes_out_so_far);
        bytes_out_so_far = d_stream.total_out;

        if (err == Z_STREAM_END) break;  /* We reached the end of the input */
    }

    /* Set the output file size to be the total output data */
    *filesize = d_stream.total_out;

    /* End the decompression */
    err = _pyfits_inflateEnd(&d_stream);

    /* free temporary output data buffer */
    free(uncompr);

    if (err != Z_OK)
    {
        return(*status = 414);
    }
    
    return(*status);
}
/*--------------------------------------------------------------------------*/
int _pyfits_compress2mem_from_mem(
             char *inmemptr,     /* I - memory pointer to uncompressed bytes */
             size_t inmemsize,   /* I - size of input uncompressed file      */
             char **buffptr,   /* IO - memory pointer for compressed file    */
             size_t *buffsize,   /* IO - size of buffer, in bytes           */
             void *(*mem_realloc)(void *p, size_t newsize), /* function     */
             size_t *filesize,   /* O - size of file, in bytes              */
             int *status)        /* IO - error status                       */

/*
  Compress the file into memory.  Fill whatever amount of memory has
  already been allocated, then realloc more memory, using the supplied
  input function, if necessary.
*/
{
    int err;
    uLong comprLen;
    Byte *compr;

    z_stream c_stream;  /* compression stream */

    uLong bytes_out_so_far = 0;  /* Keeps track of the number of bytes put in
                                    the output buffer so far */

    if (*status > 0)
        return(*status);

    /* Allocate memory as a temporary buffer in which to compress. */
    comprLen = *buffsize;
    compr = (Byte*)malloc(*buffsize);

    c_stream.zalloc = (alloc_func)0;
    c_stream.zfree = (free_func)0;
    c_stream.opaque = (voidpf)0;

    /* Initialize the compression.  The argument (15+16) tells the 
       compressor that we are to use the gzip algorythm */
    err = deflateInit2(&c_stream, Z_DEFAULT_COMPRESSION, Z_DEFLATED,
                       (15+16), 8, Z_DEFAULT_STRATEGY);

    if (err != Z_OK)
    {
        return(*status = 413);
    }

    c_stream.next_in = (unsigned char*)inmemptr;
    c_stream.avail_in = inmemsize;

    for (;;)
    {
        /* Output to the temporary buffer.  This will overwrite the
           previous data each time. */
        c_stream.next_out = compr;
        c_stream.avail_out = comprLen;

        err = _pyfits_deflate(&c_stream, Z_FINISH);

        if (err != Z_OK && err != Z_STREAM_END)
        {
            /* free temporary output data buffer */
            free(compr);
            return(*status = 413);
        }

        if (c_stream.total_out > *buffsize)
        {
            /* OK, we need more memory for the output so reallocate it */
            *buffsize = c_stream.total_out;
            *buffptr = mem_realloc(*buffptr,*buffsize);

            if (*buffptr == NULL)
            {
                /* free temporary output data buffer */
                free(compr);
                return(*status = 413);
            }
        }

        /* copy from the temporary buffer into the output memory buffer */
        memcpy((char *) *buffptr + bytes_out_so_far, (char *) compr,
               c_stream.total_out-bytes_out_so_far);
        bytes_out_so_far = c_stream.total_out;

        if (err == Z_STREAM_END) break;  /* We reached the end of the input */
    }

    /* Set the output file size to be the total output data */
    *filesize = c_stream.total_out;

    /* End the compression */
    err = _pyfits_deflateEnd(&c_stream);

    /* free temporary output data buffer */
    free(compr);

    if (err != Z_OK)
    {
        return(*status = 413);
    }
     
    return(*status);
}

