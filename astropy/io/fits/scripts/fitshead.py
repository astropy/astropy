# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
``fitshead`` is a command line script based on astropy.io.fits for printing
the header(s) of a FITS file to the standard output.

Example uses of fitshead:

1. Print the header of all the HDUs of a single .fits file:

    $ fitshead filename.fits

2. Print the header of the third HDU extension:

    $ fitshead -e 3 filename.fits

3. Print the headers of all fits files in a directory:

    $ fitshead *.fits
"""
from __future__ import absolute_import, division, print_function, unicode_literals

from ... import fits
from .... import log


def print_header(filename, extno=None):
    """Prints the header(s) of a FITS file to the standard output.

    Parameters
    ----------
    filename : str
        Path to the FITS file for which the header should be printed.

    extno : int, optional
        Number of the HDU extension to print. If unspecified, then all
        extensions will be printed.
    """
    try:
        f = fits.open(filename)
    except IOError as e:
        # Display a user-friendly error message rather than an exception trace
        log.error(e)
        return

    if extno is None:  # If no extension was specified, display all
        extno_list = range(len(f))
    else:
        if extno > len(f)-1 or extno < 0:  # Make sure the extension exists
            log.error('Extension #{0} is not present in {1}.'.format(extno,
                                                                     filename))
            return
        extno_list = [extno]

    for number in extno_list:
        if number != extno_list[0]:
            print('')  # Separate multiple HDUs by newlines
        print('# HDU {0} in {1}:'.format(number, filename))
        print(repr(f[number].header))
        print('END')


def main(args=None):
    from astropy.utils.compat import argparse

    parser = argparse.ArgumentParser(
        description=("Print the header(s) of a FITS file. "
                     "By default, all HDU extensions are shown."))
    parser.add_argument('-e', '--ext', metavar='hdu', type=int,
                        help='display the specified hdu extension number')
    parser.add_argument('filename', nargs='+',
                        help='path to one or more FITS files to display')
    args = parser.parse_args(args)

    for filename in args.filename:
        print_header(filename, args.ext)
