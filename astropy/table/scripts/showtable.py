# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
``showtable`` is a command-line script based on astropy.io and astropy.table
for printing ASCII, FITS, HDF5(?) or VOTable files(s) to the standard output.

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


def showtable(filename, args):
    """
    Read a table and print to the standard output.

    Parameters
    ----------
    filename : str
        The path to a FITS file.

    """
    print(args)
    read_kwargs = {k: v for k, v in vars(args).items()
                   if k in ('hdu', 'format', 'table_id') and k is not None}
    try:
        table = Table.read(filename, **read_kwargs)
        table.pprint(max_lines=args.max_lines, max_width=args.max_width,
                     show_unit=not args.hide_unit, show_dtype=args.show_dtype)
    except IOError as e:
        log.error(str(e))


def main(args=None):
    """The main function called by the `fitsinfo` script."""
    parser = argparse.ArgumentParser(
        description=('Print tables from ASCII, FITS, HDF5, VOTable file(s).'))
    addarg = parser.add_argument

    # pprint arguments
    addarg('--max-lines', type=int,
           help='Maximum number of lines in table output.')
    addarg('--max-width', type=int,
           help='Maximum character width of output.')
    addarg('--hide-unit', action='store_true',
           help='Hide the header row for unit (which is shown '
           'only if one or more columns has a unit).')
    addarg('--show-dtype', action='store_true',
           help='Include a header row for column dtypes.')

    # ASCII-specific arguments
    # FIXME: add more args ? (delimiter, guess ?)
    addarg('--format', help='Input table format (only for ASCII files).')

    # FITS-specific arguments
    addarg('--hdu', help='Name of the HDU to show (only for FITS files).')

    # HDF5-specific arguments
    addarg('--path', help='The path from which to read the table (only '
           'for HDF5 files).')

    # VOTable-specific arguments
    addarg('--table_id', help='The table to read in (only for VOTable files).')

    addarg('filename', nargs='+', help='Path to one or more files.')

    args = parser.parse_args(args)

    for idx, filename in enumerate(args.filename):
        if idx > 0:
            print()
        showtable(filename, args)
