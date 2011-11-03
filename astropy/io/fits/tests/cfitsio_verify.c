/* This script verifies .fits checksums using CFITSIO to demonstrate
compatibility with PyFITS.   Since running it requires compiling and
linking against cfitsio,  the script is included as a maintenance 
asset but not automatically compiled and run.

After installing cfitsio to ~/include and ~/lib,  I built cfitsio_verify
like this:

% gcc cfitsio_verify.c -I~/include -L~/lib -lcfitsio -lm -o cfitsio_verify

Run cfitsio_verify like this:

% cfitsio_verify tmp.fits

*/

#include <fitsio.h>

char * verify_status(int status)
{
	if (status == 1) { 
		return "ok";
	} else if (status == 0) {
		return "missing";
	} else if (status == -1) {
		return "error";
	}
}

int main(int argc, char *argv[])
{
	fitsfile *fptr;
	int i, j, status, dataok, hduok, hdunum, hdutype;
	char *hdustr, *datastr;

	for (i=1; i<argc; i++) {

		fits_open_file(&fptr, argv[i], READONLY, &status);
		if (status) { 
			fits_report_error(stderr, status);
			exit(-1); 
		}

		fits_get_num_hdus(fptr, &hdunum, &status);
		if (status) { 
			fprintf(stderr, "Bad get_num_hdus status for '%s' = %d", 
				argv[i], status);
			exit(-1); 
		}

		for (j=0; j<hdunum; j++) {
			fits_movabs_hdu(fptr, hdunum, &hdutype, &status);
			if (status) { 
				fprintf(stderr, "Bad movabs status for '%s[%d]' = %d.", 
					argv[i], j, status);
				exit(-1);
			}
			fits_verify_chksum(fptr, &dataok, &hduok, &status);
			if (status) { 
				fprintf(stderr, "Bad verify status for '%s[%d]' = %d.", 
					argv[i], j, status);
				exit(-1);
			}
			datastr = verify_status(dataok);
			hdustr = verify_status(hduok);
			printf("Verifying '%s[%d]'  data='%s'   hdu='%s'.\n",
			       argv[i], j, datastr, hdustr);
		}
	}
}

