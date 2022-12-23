int fits_quantize_float (long row, float fdata[], long nxpix, long nypix, int nullcheck,
	float in_null_value, float qlevel, int dither_method, int idata[], double *bscale,
	double *bzero, int *iminval, int *imaxval);

int fits_quantize_double (long row, double fdata[], long nxpix, long nypix, int nullcheck,
	double in_null_value, float qlevel, int dither_method, int idata[], double *bscale,
	double *bzero, int *iminval, int *imaxval);
