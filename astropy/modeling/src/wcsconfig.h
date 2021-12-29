
    /* The bundled version has WCSLIB_VERSION */
    #define HAVE_WCSLIB_VERSION 1

    /* WCSLIB library version number. */
    #define WCSLIB_VERSION 7.7

    /* 64-bit integer data type. */
    #define WCSLIB_INT64 long long int

    /* Windows needs some other defines to prevent inclusion of wcsset()
       which conflicts with wcslib's wcsset().  These need to be set
       on code that *uses* astropy.wcs, in addition to astropy.wcs itself.
       */
    #if defined(_WIN32) || defined(_MSC_VER) || defined(__MINGW32__) || defined (__MINGW64__)

    #ifndef YY_NO_UNISTD_H
    #define YY_NO_UNISTD_H
    #endif

    #ifndef _CRT_SECURE_NO_WARNINGS
    #define _CRT_SECURE_NO_WARNINGS
    #endif

    #ifndef _NO_OLDNAMES
    #define _NO_OLDNAMES
    #endif

    #ifndef NO_OLDNAMES
    #define NO_OLDNAMES
    #endif

    #ifndef __STDC__
    #define __STDC__ 1
    #endif

    #endif
    