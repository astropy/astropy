"""Downloads the FITS files that are used in image testing and for building documentation.


"""

from astropy.utils.data import download_file
from astropy.io import fits


def msx_hdu(cache=True):
    filename = download_file("http://astrofrog.github.io/wcsaxes-datasets/msx.fits", cache=cache)
    return fits.open(filename)[0]


def rosat_hdu(cache=True):
    filename = download_file("http://astrofrog.github.io/wcsaxes-datasets/rosat.fits", cache=cache)
    return fits.open(filename)[0]


def twoMASS_k_hdu(cache=True):
    filename = download_file("http://astrofrog.github.io/wcsaxes-datasets/2MASS_k.fits", cache=cache)
    return fits.open(filename)[0]


def l1448_co_hdu(cache=True):
    filename = download_file("http://astrofrog.github.io/wcsaxes-datasets/L1448_13CO_subset.fits", cache=cache)
    return fits.open(filename)[0]


def bolocam_hdu(cache=True):
    filename = download_file("http://astrofrog.github.io/wcsaxes-datasets/bolocam_v2.0.fits", cache=cache)
    return fits.open(filename)[0]
