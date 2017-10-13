# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
``showtable`` is a command-line script based on astropy.io and astropy.table
for printing ASCII, FITS, HDF5(?) or VOTABLE files(s) to the standard output.

Example usage of ``showtable``:

1. Print a summary of the HDUs in a FITS file::

    $ showtable astropy/io/fits/tests/data/table.fits

    target V_mag
    ------- -----
    NGC1001  11.1
    NGC1002  12.3
    NGC1003  15.2

3. ASCII::

    $ showtable astropy/io/ascii/tests/t/simple_csv.csv

     a   b   c
    --- --- ---
      1   2   3
      4   5   6

3. XML::

    $ showtable astropy/io/votable/tests/data/names.xml
               col1             col2     col3    col4     col5   ... col13 col14 col15 col16 col17
               ---              deg      deg     deg      deg    ...  mag   mag   mag   mag   ---
    ------------------------- -------- ------- -------- -------- ... ----- ----- ----- ----- -----
    SSTGLMC G000.0000+00.1611   0.0000  0.1611 266.2480 -28.8521 ...  9.13  8.17    --    --    AA



2. Print a summary of HDUs of all the FITS files in the current directory::

    $ showtable *.fits

"""

import argparse
from astropy.table import Table
from astropy import log


def showtable(filename, **kwargs):
    """
    Read a table and print to the standard output.

    Parameters
    ----------
    filename : str
        The path to a FITS file.

    """
    # print(kwargs)
    try:
        table = Table.read(filename)
        table.pprint(**kwargs)
    except IOError as e:
        log.error(str(e))


def main(args=None):
    """The main function called by the `fitsinfo` script."""
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=('Print tables from ASCII, FITS, HDF5, VOTABLE file(s).'))
    parser.add_argument('--max-lines', type=int,
                        help='Maximum number of lines in table output.')
    parser.add_argument('--max-width', type=int,
                        help='Maximum character width of output.')
    parser.add_argument('--show-name', default=True,
                        help='Include a header row for column names.')
    parser.add_argument('--show-unit',
                        help='Include a header row for unit.  Default is to '
                        'show a row for units only if one or more columns')
    parser.add_argument('--show-dtype', default=False,
                        help='Include a header row for column dtypes.')
    parser.add_argument('filename', nargs='+',
                        help='Path to one or more FITS files.')

    # TODO: fix bool args
    # TODO: add `read` kwargs (to choose format, hdu, etc.)

    args = parser.parse_args(args)
    args = vars(args)
    filenames = args.pop('filename')

    for idx, filename in enumerate(filenames):
        if idx > 0:
            print()
        showtable(filename, **args)
