/*
** Copyright (C) 2016-2017, NumFOCUS Foundation. 
**
** Licensed under a 3-clause BSD style license - see LICENSE
**
** This file is NOT derived from SOFA sources
*/


#ifndef _ERFA_EXTRA_H
#define _ERFA_EXTRA_H

#ifdef __cplusplus
extern "C" {
#endif


/* 
** Returns the package version
** as defined in configure.ac
** in string format
*/
const char* eraVersion();

/* 
** Returns the package major version
** as defined in configure.ac
** as integer
*/
int eraVersionMajor();

/* 
** Returns the package minor version
** as defined in configure.ac
** as integer
*/
int eraVersionMinor();

/* 
** Returns the package micro version
** as defined in configure.ac
** as integer
*/
int eraVersionMicro();

/* 
** Returns the orresponding SOFA version
** as defined in configure.ac
** in string format
*/
const char* eraSofaVersion();


#ifdef __cplusplus
}
#endif

#endif /* _ERFA_EXTRA_H */

