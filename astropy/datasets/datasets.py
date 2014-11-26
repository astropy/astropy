# Licensed under a 3-clause BSD style license - see LICENSE.rst

import json

from ..io import fits
from ..utils.data import download_file, get_readable_fileobj

__all__ = ['fetch_hdu', 'list_available_datasets']

URL = 'http://data.astropy.org/'


def fetch_hdu(filename, cache=True):
    """
    Download a FITS file to the cache and open HDU 0.
    """
    path = download_file(URL + filename, cache=cache)
    return fits.open(path)[0]


def list_available_datasets():

    with get_readable_fileobj(URL + 'file_index.json', cache=False, show_progress=False) as f:
        datasets = json.load(f)

    datasets.remove('')  # last entry in json file

    datasets = [d for d in datasets if not d.startswith('tutorials')]

    return datasets
