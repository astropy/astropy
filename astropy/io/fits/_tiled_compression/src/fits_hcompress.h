typedef long long LONGLONG;

int fits_hcompress(int *a, int ny, int nx, int scale, char *output,
                   long *nbytes, int *status);
int fits_hcompress64(LONGLONG *a, int ny, int nx, int scale, char *output,
                  long *nbytes, int *status);
