# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
``showtable`` is a command-line script based on ``astropy.io`` and
``astropy.table`` for printing ASCII, FITS, HDF5 or VOTable files(s) to the
standard output.

Example usage of ``showtable``:

1. FITS::

    $ showtable astropy/io/fits/tests/data/table.fits

    target V_mag
    ------- -----
    NGC1001  11.1
    NGC1002  12.3
    NGC1003  15.2

2. ASCII::

    $ showtable astropy/io/ascii/tests/t/simple_csv.csv

     a   b   c
    --- --- ---
      1   2   3
      4   5   6

3. XML::

    $ showtable astropy/io/votable/tests/data/names.xml --max-width 70

               col1             col2     col3  ... col15 col16 col17
               ---              deg      deg   ...  mag   mag   ---
    ------------------------- -------- ------- ... ----- ----- -----
    SSTGLMC G000.0000+00.1611   0.0000  0.1611 ...    --    --    AA



4. Print all the FITS tables in the current directory::

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
    read_kwargs = {k: v for k, v in vars(args).items()
                   if k in ('hdu', 'format', 'table_id') and v is not None}
    try:
        table = Table.read(filename, **read_kwargs)
        formatter = table.more if args.more else table.pprint
        formatter(max_lines=args.max_lines, max_width=args.max_width,
                  show_unit=not args.hide_unit, show_dtype=args.show_dtype)
    except IOError as e:
        log.error(str(e))


def main(args=None):
    """The main function called by the `showtable` script."""
    parser = argparse.ArgumentParser(
        description=(
            'Print tables from ASCII, FITS, HDF5, VOTable file(s).'
            'The default behavior is make the table output fit onto a single '
            'screen page.  For a long and wide table this will mean cutting '
            'out inner rows and columns.  To print **all** the rows or columns'
            ' use ``--max-lines=-1`` or ``max-width=-1``, respectively.'
        ))

    addarg = parser.add_argument

    # pprint arguments
    addarg('--more', action='store_true',
           help='Use the pager mode from Table.more.')
    addarg('--max-lines', type=int,
           help='Maximum number of lines in table output (default=screen '
           'length, -1 for no limit).')
    addarg('--max-width', type=int,
           help='Maximum width in table output (default=screen width, '
           '-1 for no limit).')
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
