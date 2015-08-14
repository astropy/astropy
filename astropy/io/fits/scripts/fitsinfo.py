# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
``fitsinfo`` is a command-line script based on astropy.io.fits for
printing a summary of the HDUs in one or more FITS files(s) to the
standard output.

Example usage of ``fitsinfo``:

1. Print a summary of the HDUs in a FITS file::

    $ fitsinfo filename.fits

    Filename: filename.fits
    No.    Name         Type      Cards   Dimensions   Format
    0    PRIMARY     PrimaryHDU     138   ()
    1    SCI         ImageHDU        61   (800, 800)   int16
    2    SCI         ImageHDU        61   (800, 800)   int16
    3    SCI         ImageHDU        61   (800, 800)   int16
    4    SCI         ImageHDU        61   (800, 800)   int16

2. Print a summary of HDUs of all the FITS files in the current directory::

    $ fitsinfo *.fits
"""

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
        fits.info(filename)
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
