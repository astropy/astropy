"""Downloads the FITS files that are used in image testing and for building documentation.
"""

from astropy.utils.data import download_file
from astropy.io import fits

URL = 'http://astrofrog.github.io/wcsaxes-datasets/'


def get_hdu(filename, cache=True):
    path = download_file(URL + filename, cache=cache)
    return fits.open(path)[0]


def msx_hdu(cache=True):
    return get_hdu('msx.fits', cache=cache)


def rosat_hdu(cache=True):
    return get_hdu('rosat.fits', cache=cache)


def twoMASS_k_hdu(cache=True):
    return get_hdu('2MASS_k.fits', cache=cache)


def l1448_co_hdu(cache=True):
    return get_hdu('L1448_13CO_subset.fits', cache=cache)


def bolocam_hdu(cache=True):
    return get_hdu('bolocam_v2.0.fits', cache=cache)
