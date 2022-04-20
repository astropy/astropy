# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
``fitsheader`` is a command line script based on astropy.io.fits for printing
the header(s) of one or more FITS file(s) to the standard output in a human-
readable format.

Example uses of fitsheader:

1. Print the header of all the HDUs of a .fits file::

    $ fitsheader filename.fits

2. Print the header of the third and fifth HDU extension::

    $ fitsheader --extension 3 --extension 5 filename.fits

3. Print the header of a named extension, e.g. select the HDU containing
   keywords EXTNAME='SCI' and EXTVER='2'::

    $ fitsheader --extension "SCI,2" filename.fits

4. Print only specific keywords::

    $ fitsheader --keyword BITPIX --keyword NAXIS filename.fits

5. Print keywords NAXIS, NAXIS1, NAXIS2, etc using a wildcard::

    $ fitsheader --keyword NAXIS* filename.fits

6. Dump the header keywords of all the files in the current directory into a
   machine-readable csv file::

    $ fitsheader --table ascii.csv *.fits > keywords.csv

7. Specify hierarchical keywords with the dotted or spaced notation::

    $ fitsheader --keyword ESO.INS.ID filename.fits
    $ fitsheader --keyword "ESO INS ID" filename.fits

8. Compare the headers of different fits files, following ESO's ``fitsort``
   format::

    $ fitsheader --fitsort --extension 0 --keyword ESO.INS.ID *.fits

9. Same as above, sorting the output along a specified keyword::

    $ fitsheader -f -s DATE-OBS -e 0 -k DATE-OBS -k ESO.INS.ID *.fits

10. Sort first by OBJECT, then DATE-OBS::

    $ fitsheader -f -s OBJECT -s DATE-OBS *.fits

Note that compressed images (HDUs of type
:class:`~astropy.io.fits.CompImageHDU`) really have two headers: a real
BINTABLE header to describe the compressed data, and a fake IMAGE header
representing the image that was compressed. Astropy returns the latter by
default. You must supply the ``--compressed`` option if you require the real
header that describes the compression.

With Astropy installed, please run ``fitsheader --help`` to see the full usage
documentation.
"""

import sys
import argparse

import numpy as np

from astropy.io import fits
from astropy import log, __version__


DESCRIPTION = """
Print the header(s) of a FITS file. Optional arguments allow the desired
extension(s), keyword(s), and output format to be specified.
Note that in the case of a compressed image, the decompressed header is
shown by default.

This script is part of the Astropy package. See
https://docs.astropy.org/en/latest/io/fits/usage/scripts.html#module-astropy.io.fits.scripts.fitsheader
for further documentation.
""".strip()


class ExtensionNotFoundException(Exception):
    """Raised if an HDU extension requested by the user does not exist."""
    pass


class HeaderFormatter:
    """Class to format the header(s) of a FITS file for display by the
    `fitsheader` tool; essentially a wrapper around a `HDUList` object.

    Example usage:
    fmt = HeaderFormatter('/path/to/file.fits')
    print(fmt.parse(extensions=[0, 3], keywords=['NAXIS', 'BITPIX']))

    Parameters
    ----------
    filename : str
        Path to a single FITS file.
    verbose : bool
        Verbose flag, to show more information about missing extensions,
        keywords, etc.

    Raises
    ------
    OSError
        If `filename` does not exist or cannot be read.
    """

    def __init__(self, filename, verbose=True):
        self.filename = filename
        self.verbose = verbose
        self._hdulist = fits.open(filename)

    def parse(self, extensions=None, keywords=None, compressed=False):
        """Returns the FITS file header(s) in a readable format.

        Parameters
        ----------
        extensions : list of int or str, optional
            Format only specific HDU(s), identified by number or name.
            The name can be composed of the "EXTNAME" or "EXTNAME,EXTVER"
            keywords.

        keywords : list of str, optional
            Keywords for which the value(s) should be returned.
            If not specified, then the entire header is returned.

        compressed : bool, optional
            If True, shows the header describing the compression, rather than
            the header obtained after decompression. (Affects FITS files
            containing `CompImageHDU` extensions only.)

        Returns
        -------
        formatted_header : str or astropy.table.Table
            Traditional 80-char wide format in the case of `HeaderFormatter`;
            an Astropy Table object in the case of `TableHeaderFormatter`.
        """
        # `hdukeys` will hold the keys of the HDUList items to display
        if extensions is None:
            hdukeys = range(len(self._hdulist))  # Display all by default
        else:
            hdukeys = []
            for ext in extensions:
                try:
                    # HDU may be specified by number
                    hdukeys.append(int(ext))
                except ValueError:
                    # The user can specify "EXTNAME" or "EXTNAME,EXTVER"
                    parts = ext.split(',')
                    if len(parts) > 1:
                        extname = ','.join(parts[0:-1])
                        extver = int(parts[-1])
                        hdukeys.append((extname, extver))
                    else:
                        hdukeys.append(ext)

        # Having established which HDUs the user wants, we now format these:
        return self._parse_internal(hdukeys, keywords, compressed)

    def _parse_internal(self, hdukeys, keywords, compressed):
        """The meat of the formatting; in a separate method to allow overriding.
        """
        result = []
        for idx, hdu in enumerate(hdukeys):
            try:
                cards = self._get_cards(hdu, keywords, compressed)
            except ExtensionNotFoundException:
                continue

            if idx > 0:  # Separate HDUs by a blank line
                result.append('\n')
            result.append(f'# HDU {hdu} in {self.filename}:\n')
            for c in cards:
                result.append(f'{c}\n')
        return ''.join(result)

    def _get_cards(self, hdukey, keywords, compressed):
        """Returns a list of `astropy.io.fits.card.Card` objects.

        This function will return the desired header cards, taking into
        account the user's preference to see the compressed or uncompressed
        version.

        Parameters
        ----------
        hdukey : int or str
            Key of a single HDU in the HDUList.

        keywords : list of str, optional
            Keywords for which the cards should be returned.

        compressed : bool, optional
            If True, shows the header describing the compression.

        Raises
        ------
        ExtensionNotFoundException
            If the hdukey does not correspond to an extension.
        """
        # First we obtain the desired header
        try:
            if compressed:
                # In the case of a compressed image, return the header before
                # decompression (not the default behavior)
                header = self._hdulist[hdukey]._header
            else:
                header = self._hdulist[hdukey].header
        except (IndexError, KeyError):
            message = f'{self.filename}: Extension {hdukey} not found.'
            if self.verbose:
                log.warning(message)
            raise ExtensionNotFoundException(message)

        if not keywords:  # return all cards
            cards = header.cards
        else:  # specific keywords are requested
            cards = []
            for kw in keywords:
                try:
                    crd = header.cards[kw]
                    if isinstance(crd, fits.card.Card):  # Single card
                        cards.append(crd)
                    else:  # Allow for wildcard access
                        cards.extend(crd)
                except KeyError:  # Keyword does not exist
                    if self.verbose:
                        log.warning('{filename} (HDU {hdukey}): '
                                    'Keyword {kw} not found.'.format(
                                        filename=self.filename,
                                        hdukey=hdukey,
                                        kw=kw))
        return cards

    def close(self):
        self._hdulist.close()


class TableHeaderFormatter(HeaderFormatter):
    """Class to convert the header(s) of a FITS file into a Table object.
    The table returned by the `parse` method will contain four columns:
    filename, hdu, keyword, and value.

    Subclassed from HeaderFormatter, which contains the meat of the formatting.
    """

    def _parse_internal(self, hdukeys, keywords, compressed):
        """Method called by the parse method in the parent class."""
        tablerows = []
        for hdu in hdukeys:
            try:
                for card in self._get_cards(hdu, keywords, compressed):
                    tablerows.append({'filename': self.filename,
                                      'hdu': hdu,
                                      'keyword': card.keyword,
                                      'value': str(card.value)})
            except ExtensionNotFoundException:
                pass

        if tablerows:
            from astropy import table
            return table.Table(tablerows)
        return None


def print_headers_traditional(args):
    """Prints FITS header(s) using the traditional 80-char format.

    Parameters
    ----------
    args : argparse.Namespace
        Arguments passed from the command-line as defined below.
    """
    for idx, filename in enumerate(args.filename):  # support wildcards
        if idx > 0 and not args.keyword:
            print()  # print a newline between different files

        formatter = None
        try:
            formatter = HeaderFormatter(filename)
            print(formatter.parse(args.extensions,
                                  args.keyword,
                                  args.compressed), end='')
        except OSError as e:
            log.error(str(e))
        finally:
            if formatter:
                formatter.close()


def print_headers_as_table(args):
    """Prints FITS header(s) in a machine-readable table format.

    Parameters
    ----------
    args : argparse.Namespace
        Arguments passed from the command-line as defined below.
    """
    tables = []
    # Create a Table object for each file
    for filename in args.filename:  # Support wildcards
        formatter = None
        try:
            formatter = TableHeaderFormatter(filename)
            tbl = formatter.parse(args.extensions,
                                  args.keyword,
                                  args.compressed)
            if tbl:
                tables.append(tbl)
        except OSError as e:
            log.error(str(e))  # file not found or unreadable
        finally:
            if formatter:
                formatter.close()

    # Concatenate the tables
    if len(tables) == 0:
        return False
    elif len(tables) == 1:
        resulting_table = tables[0]
    else:
        from astropy import table
        resulting_table = table.vstack(tables)
    # Print the string representation of the concatenated table
    resulting_table.write(sys.stdout, format=args.table)


def print_headers_as_comparison(args):
    """Prints FITS header(s) with keywords as columns.

    This follows the dfits+fitsort format.

    Parameters
    ----------
    args : argparse.Namespace
        Arguments passed from the command-line as defined below.
    """
    from astropy import table
    tables = []
    # Create a Table object for each file
    for filename in args.filename:  # Support wildcards
        formatter = None
        try:
            formatter = TableHeaderFormatter(filename, verbose=False)
            tbl = formatter.parse(args.extensions,
                                  args.keyword,
                                  args.compressed)
            if tbl:
                # Remove empty keywords
                tbl = tbl[np.where(tbl['keyword'] != '')]
            else:
                tbl = table.Table([[filename]], names=('filename',))
            tables.append(tbl)
        except OSError as e:
            log.error(str(e))  # file not found or unreadable
        finally:
            if formatter:
                formatter.close()

    # Concatenate the tables
    if len(tables) == 0:
        return False
    elif len(tables) == 1:
        resulting_table = tables[0]
    else:
        resulting_table = table.vstack(tables)

    # If we obtained more than one hdu, merge hdu and keywords columns
    hdus = resulting_table['hdu']
    if np.ma.isMaskedArray(hdus):
        hdus = hdus.compressed()
    if len(np.unique(hdus)) > 1:
        for tab in tables:
            new_column = table.Column(
                [f"{row['hdu']}:{row['keyword']}" for row in tab])
            tab.add_column(new_column, name='hdu+keyword')
        keyword_column_name = 'hdu+keyword'
    else:
        keyword_column_name = 'keyword'

    # Check how many hdus we are processing
    final_tables = []
    for tab in tables:
        final_table = [table.Column([tab['filename'][0]], name='filename')]
        if 'value' in tab.colnames:
            for row in tab:
                if row['keyword'] in ('COMMENT', 'HISTORY'):
                    continue
                final_table.append(table.Column([row['value']],
                                                name=row[keyword_column_name]))
        final_tables.append(table.Table(final_table))
    final_table = table.vstack(final_tables)
    # Sort if requested
    if args.sort:
        final_table.sort(args.sort)
    # Reorganise to keyword by columns
    final_table.pprint(max_lines=-1, max_width=-1)


def main(args=None):
    """This is the main function called by the `fitsheader` script."""

    parser = argparse.ArgumentParser(
        description=DESCRIPTION,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument(
        '--version', action='version',
        version=f'%(prog)s {__version__}')

    parser.add_argument('-e', '--extension', metavar='HDU',
                        action='append', dest='extensions',
                        help='specify the extension by name or number; '
                             'this argument can be repeated '
                             'to select multiple extensions')
    parser.add_argument('-k', '--keyword', metavar='KEYWORD',
                        action='append', type=str,
                        help='specify a keyword; this argument can be '
                             'repeated to select multiple keywords; '
                             'also supports wildcards')
    mode_group = parser.add_mutually_exclusive_group()
    mode_group.add_argument('-t', '--table',
                            nargs='?', default=False, metavar='FORMAT',
                            help='print the header(s) in machine-readable table '
                                 'format; the default format is '
                                 '"ascii.fixed_width" (can be "ascii.csv", '
                                 '"ascii.html", "ascii.latex", "fits", etc)')
    mode_group.add_argument('-f', '--fitsort', action='store_true',
                            help='print the headers as a table with each unique '
                                 'keyword in a given column (fitsort format) ')
    parser.add_argument('-s', '--sort', metavar='SORT_KEYWORD',
                        action='append', type=str,
                        help='sort output by the specified header keywords, '
                             'can be repeated to sort by multiple keywords; '
                             'Only supported with -f/--fitsort')
    parser.add_argument('-c', '--compressed', action='store_true',
                        help='for compressed image data, '
                             'show the true header which describes '
                             'the compression rather than the data')
    parser.add_argument('filename', nargs='+',
                        help='path to one or more files; '
                             'wildcards are supported')
    args = parser.parse_args(args)

    # If `--table` was used but no format specified,
    # then use ascii.fixed_width by default
    if args.table is None:
        args.table = 'ascii.fixed_width'

    if args.sort:
        args.sort = [key.replace('.', ' ') for key in args.sort]
        if not args.fitsort:
            log.error('Sorting with -s/--sort is only supported in conjunction with -f/--fitsort')
            # 2: Unix error convention for command line syntax
            sys.exit(2)

    if args.keyword:
        args.keyword = [key.replace('.', ' ') for key in args.keyword]

    # Now print the desired headers
    try:
        if args.table:
            print_headers_as_table(args)
        elif args.fitsort:
            print_headers_as_comparison(args)
        else:
            print_headers_traditional(args)
    except OSError:
        # A 'Broken pipe' OSError may occur when stdout is closed prematurely,
        # eg. when calling `fitsheader file.fits | head`. We let this pass.
        pass
