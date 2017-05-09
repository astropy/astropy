# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Downloads the FITS files that are used in image testing and for building documentation.
"""

import time

from ....extern.six.moves import urllib
from ....utils.data import download_file
from ....io import fits

__all__ = ['fetch_msx_hdu',
           'fetch_rosat_hdu',
           'fetch_twoMASS_k_hdu',
           'fetch_l1448_co_hdu',
           'fetch_bolocam_hdu',
           ]

MAX_RETRIES = 10
TIME_BETWEEN_RETRIES = 5
URL = 'http://data.astropy.org/'


def fetch_hdu(filename, cache=True):
    """Download a FITS file to the cache and open HDU 0.
    """
    for retry in range(MAX_RETRIES):
        try:
            path = download_file(URL + filename, cache=cache, timeout=30)
        except urllib.error.URLError:
            if retry == MAX_RETRIES - 1:
                raise
            else:
                time.sleep(TIME_BETWEEN_RETRIES)
        else:
            break
    else:
        raise Exception("Failed to download file {0}".format(filename))
    return fits.open(path)[0]


def fetch_msx_hdu(cache=True):
    """Fetch the MSX example dataset HDU.

    Returns
    -------
    hdu : `~astropy.io.fits.ImageHDU`
        Image HDU
    """
    return fetch_hdu('galactic_center/gc_msx_e.fits', cache=cache)


def fetch_rosat_hdu(cache=True):
    return fetch_hdu('allsky/allsky_rosat.fits', cache=cache)


def fetch_twoMASS_k_hdu(cache=True):
    return fetch_hdu('galactic_center/gc_2mass_k.fits', cache=cache)


def fetch_l1448_co_hdu(cache=True):
    return fetch_hdu('l1448/l1448_13co.fits', cache=cache)


def fetch_bolocam_hdu(cache=True):
    return fetch_hdu('galactic_center/gc_bolocam_gps.fits', cache=cache)
