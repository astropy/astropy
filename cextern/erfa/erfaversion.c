/*
** Copyright (C) 2016-2017, NumFOCUS Foundation.
**
** Licensed under a 3-clause BSD style license - see LICENSE
**
** This file is NOT derived from SOFA sources
*/

/* Define to the version of this package. */
#define PACKAGE_VERSION "1.6.0"

/* Define to the major version of this package. */
#define PACKAGE_VERSION_MAJOR 1

/* Define to the micro version of this package. */
#define PACKAGE_VERSION_MICRO 0

/* Define to the minor version of this package. */
#define PACKAGE_VERSION_MINOR 6

/* Define to the version of SOFA */
#define SOFA_VERSION "20190722"

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */


const char* eraVersion(void) {
  return PACKAGE_VERSION;
}


int eraVersionMajor(void) {
  return PACKAGE_VERSION_MAJOR;
}


int eraVersionMinor(void) {
  return PACKAGE_VERSION_MINOR;
}


int eraVersionMicro(void) {
  return PACKAGE_VERSION_MICRO;
}


const char* eraSofaVersion(void) {
  return SOFA_VERSION;
}
