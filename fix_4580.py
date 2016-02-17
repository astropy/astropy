# fix for the FITS I/O bug #4580
import warnings
from astropy.io import fits
import numpy as np
from astropy.io.fits.hdu import PrimaryHDU,BinTableHDU,HDUList
dum = np.empty((5,100), dtype=[('w','>f8'),('f','>f4')])
prihdu = fits.PrimaryHDU()
tblhdu = fits.BinTableHDU(dum, name='DATA')
hdulist = fits.HDUList([prihdu, tblhdu])
hdulist.writeto('file.fits', clobber=True) # Does not fail
warnings.filterwarnings('ignore', category=UserWarning, append=True)
tmp = fits.open('jhfaj.fits',ignore_missing_end=True)


