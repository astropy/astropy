# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import argparse
import astropy.io.fits as fits
from astropy import log


def fitsinfo(filename):
    """
    Print a summary of the HDUs in a FITS file.

    Parameters
    ----------
    filename : str
        The path to a FITS file.
    """

    try:
        hdulist = fits.open(filename)
        hdulist.info()
        hdulist.close()
    except IOError as e:
        log.error(str(e))
    return


def main(args=None):
    """The main function called by the `fitsinfo` script."""
    parser = argparse.ArgumentParser(
        description=('Print a summary of the HDUs in a FITS file(s).'))
    parser.add_argument('filename', nargs='+',
                        help='Path to one or more FITS files. '
                             'Wildcards are supported.')
    args = parser.parse_args(args)

    for idx, filename in enumerate(args.filename):
        if idx > 0:
            print()
        fitsinfo(filename)
