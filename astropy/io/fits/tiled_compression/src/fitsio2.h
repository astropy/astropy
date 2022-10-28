#ifndef LONGLONG_TYPE   /* this may have been previously defined */
#if defined(_MSC_VER)   /* Microsoft Visual C++ */

#if (_MSC_VER < 1300)   /* versions earlier than V7.0 do not have 'long long' */
    typedef __int64 LONGLONG;
    typedef unsigned __int64 ULONGLONG;

#else                   /* newer versions do support 'long long' */
    typedef long long LONGLONG;
    typedef unsigned long long ULONGLONG;

#endif

#elif defined( __BORLANDC__)  /* for the Borland 5.5 compiler, in particular */
    typedef __int64 LONGLONG;
    typedef unsigned __int64 ULONGLONG;
#else
    typedef long long LONGLONG;
    typedef unsigned long long ULONGLONG;
#endif

#define LONGLONG_TYPE
#endif

# define DATA_COMPRESSION_ERR 413  /* error in imcompress routines */
# define DATA_DECOMPRESSION_ERR 414 /* error in imcompress routines */
# define NO_COMPRESSED_TILE  415 /* compressed tile doesn't exist */

void ffpmsg(const char *err_message);

#ifdef _REENTRANT
#include <pthread.h>
/*  #include <assert.h>  not needed any more */
extern pthread_mutex_t Fitsio_Lock;
extern int Fitsio_Pthread_Status;

#define FFLOCK1(lockname)   (Fitsio_Pthread_Status = pthread_mutex_lock(&lockname))
#define FFUNLOCK1(lockname) (Fitsio_Pthread_Status = pthread_mutex_unlock(&lockname))
#define FFLOCK   FFLOCK1(Fitsio_Lock)
#define FFUNLOCK FFUNLOCK1(Fitsio_Lock)
#define ffstrtok(str, tok, save) strtok_r(str, tok, save)

#else
#define FFLOCK
#define FFUNLOCK
#define ffstrtok(str, tok, save) strtok(str, tok)
#endif
