# Licensed under a 3-clause BSD style license - see LICENSE.rst

import json

from ..io import fits
from ..utils.data import download_file, get_readable_fileobj

__all__ = ['dataset']

URL = 'http://data.astropy.org/'


DATASETS = {
    'galcen_2mass_k': 'galactic_center/gc_2mass_k.fits',
    'galcen_bolocam_gps': 'galactic_center/gc_bolocam_gps.fits',
    'galcen_msx_e': 'galactic_center/gc_msx_e.fits',
    'allsky_rosat': 'allsky/allsky_rosat.fits',
    'l1448_13co': 'l1448/l1448_13co.fits'
}


def dataset(name, cache=True):
    """
    Return the path to a local version of dataset ``name``
    """
    return download_file(URL + DATASETS[name], cache=cache)
