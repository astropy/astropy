int fits_rcomp_byte(signed char a[], int nx, unsigned char *c, int clen, int nblock);
int fits_rcomp_short(short a[], int nx, unsigned char *c, int clen, int nblock);
int fits_rcomp(int a[], int nx, unsigned char *c, int clen, int nblock);
int fits_rdecomp_byte(unsigned char *c,	int clen, unsigned char array[], int nx, int nblock);
int fits_rdecomp_short(unsigned char *c, int clen, unsigned short array[], int nx, int nblock);
int fits_rdecomp(unsigned char *c, int clen, unsigned int array[],	int nx, int nblock);
