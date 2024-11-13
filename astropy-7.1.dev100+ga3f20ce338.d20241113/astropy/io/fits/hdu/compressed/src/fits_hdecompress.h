typedef long long LONGLONG;

int fits_hdecompress(unsigned char *input, int smooth, int *a, int *ny, int *nx,
                     int *scale, int *status);
int fits_hdecompress64(unsigned char *input, int smooth, LONGLONG *a, int *ny, int *nx,
                     int *scale, int *status);
