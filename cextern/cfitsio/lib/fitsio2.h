#ifndef LONGLONG_TYPE
typedef long long LONGLONG;
typedef unsigned long long ULONGLONG;
#define LONGLONG_TYPE
#endif

# define DATA_COMPRESSION_ERR 413
# define DATA_DECOMPRESSION_ERR 414

void ffpmsg(const char *err_message);

#define FFLOCK
#define FFUNLOCK

#define N_RANDOM 10000

#define MEMORY_ALLOCATION 113

#define NO_DITHER -1
#define SUBTRACTIVE_DITHER_1 1
#define SUBTRACTIVE_DITHER_2 2

int fits_init_randoms(void);
