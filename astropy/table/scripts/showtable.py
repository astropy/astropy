# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""

``showtable`` is a command-line script based on astropy.io and astropy.table
for printing ASCII, FITS, HDF5(?) or VOTABLE files(s) to the standard output.

Example usage of ``showtable``:

1. Print a summary of the HDUs in a FITS file::

    $ showtable filename.fits

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

import argparse
from astropy.table import Table
from astropy import log


def showtable(filename, **kwargs):
    """
    Print a summary of the HDUs in a FITS file.

    Parameters
    ----------
    filename : str
        The path to a FITS file.
    """

    try:
        table = Table.read(filename)
        table.pprint(**kwargs)
    except IOError as e:
        log.error(str(e))


def main(args=None):
    """The main function called by the `fitsinfo` script."""
    parser = argparse.ArgumentParser(
        description=('Print tables from ASCII, FITS, HDF5, VOTABLE file(s).'))
    parser.add_argument(
        'filename', nargs='+',
        help='Path to one or more FITS files. Wildcards are supported.')
    args = parser.parse_args(args)

    for idx, filename in enumerate(args.filename):
        if idx > 0:
            print()
        showtable(filename)
