# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Downloads the FITS files that are used in image testing and for building documentation.
"""

from ....utils.data import get_pkg_data_filename
from ....io import fits

__all__ = ['fetch_msx_hdu',
           'fetch_rosat_hdu',
           'fetch_twoMASS_k_hdu',
           'fetch_l1448_co_hdu',
           'fetch_bolocam_hdu',
           ]


def fetch_hdu(filename):
    """
    Download a FITS file to the cache and open HDU 0.
    """
    path = get_pkg_data_filename(filename)
    return fits.open(path)[0]


def fetch_msx_hdu():
    """Fetch the MSX example dataset HDU.

    Returns
    -------
    hdu : `~astropy.io.fits.ImageHDU`
        Image HDU
    """
    return fetch_hdu('galactic_center/gc_msx_e.fits')


def fetch_rosat_hdu():
    return fetch_hdu('allsky/allsky_rosat.fits')


def fetch_twoMASS_k_hdu():
    return fetch_hdu('galactic_center/gc_2mass_k.fits')


def fetch_l1448_co_hdu():
    return fetch_hdu('l1448/l1448_13co.fits')


def fetch_bolocam_hdu():
    return fetch_hdu('galactic_center/gc_bolocam_gps.fits')
