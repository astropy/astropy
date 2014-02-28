# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
``fitsheader`` is a command line script based on astropy.io.fits for printing
the header(s) of one ore more FITS file(s) to the standard output.

Example uses of fitsheader:

1. Print the header of all the HDUs of a .fits file:

    $ fitsheader filename.fits

2. Print the header of the third HDU extension:

    $ fitsheader --ext 3 filename.fits

3. Print the header of a named extension having EXTNAME='SCI' and EXTVER='2'
in the header:

    $ fitsheader --ext "SCI,2" filename.fits

4. Print the headers of all fits files in a directory:

    $ fitsheader *.fits

Note that compressed images (HDUs of type `CompImageHDU`) really have two
headers: a real BINTABLE header to describe the compressed data, and a fake
IMAGE header representing the image that was compressed. Astropy returns the
latter by default. You must supply the "--compressed" option if you require the
header that describes the compression.

With Astropy installed, please run ``fitscheck --help`` to see the full program
usage documentation.
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from ... import fits
from .... import log


class FormattingException(Exception):
    pass


class HeaderFormatter(object):
    """
    Base class to format the header(s) of a FITS file for terminal display.

    Parameters
    ----------
    filename : str
        Path to the FITS file.
    """
    def __init__(self, filename, compressed=False):
        try:
            self.hdulist = fits.open(filename)
        except IOError as e:
            raise FormattingException(e.message)
        self.compressed = compressed

    def parse(self, extension=None):
        """Returns the FITS file header(s) in a readable format.

        Parameters
        ----------
        extension : int or str, optional
            Format only a specific HDU, identified by its number or its name.
            The name is the "EXTNAME" or "EXTNAME,EXTVER" string.

        Returns
        -------
        formatted_header : str
            Nicely formatted header information.
        """
        # `hdukeys` will hold the keys of the HDUList items to display
        if extension is None:
            hdukeys = range(len(self.hdulist))  # Display all by default
        else:
            try:
                hdukeys = [int(extension)]  # HDU may be specified by number
            except ValueError:
                # The user can specify "EXTNAME" or "EXTNAME,EXTVER" as well
                parts = extension.split(',')
                if len(parts) > 1:
                    extname = str(','.join(parts[0:-1]))
                    extver = int(parts[-1])
                    hdukeys = [(extname, extver)]
                else:
                    hdukeys = [extension]

        # Having established which HDUs are wanted, we can now format these
        return self._format_hdulist(hdukeys)

    def _get_header(self, hdukey):
        """Returns the `astropy.io.fits.header.Header` object for the HDU."""
        try:
            if self.compressed:
                # In the case of a compressed image, return the header before
                # decompression (not the default behavior)
                return self.hdulist[hdukey]._header
            else:
                return self.hdulist[hdukey].header
        except IndexError:
            raise FormattingException('{0}: Extension #{1} not found.'.
                                      format(self.hdulist.filename(), hdukey))
        except KeyError as e:
            raise FormattingException('{0}: {1}'.format(
                                      self.hdulist.filename(), e.message))

    def _format_hdulist(self, hdukeys):
        """Returns the formatted version of the header; the important bit."""
        text = ''
        for i, key in enumerate(hdukeys):
            if i > 0:
                prefix = '\n\n'  # Separate HDUs by a blank line
            else:
                prefix = ''
            text += '{prefix}# HDU {key} in {filename}:\n{header}'.format(
                    prefix=prefix,
                    key=key,
                    filename=self.hdulist.filename(),
                    header=self._get_header(key).tostring(sep='\n',
                                                          padding=False))
        return text


def main(args=None):
    from astropy.utils.compat import argparse

    parser = argparse.ArgumentParser(
        description=('Print the header(s) of a FITS file. '
                     'All HDU extensions are shown by default. '
                     'In the case of a compressed image, '
                     'the decompressed header is shown.'))
    parser.add_argument('-e', '--ext', metavar='hdu',
                        help='specify the HDU extension number or name')
    parser.add_argument('-c', '--compressed', action='store_true',
                        help='for compressed image data, '
                             'show the true header which describes '
                             'the compression rather than the data')
    parser.add_argument('filename', nargs='+',
                        help='path to one or more FITS files to display')
    args = parser.parse_args(args)

    try:
        for filename in args.filename:
            print(HeaderFormatter(filename, args.compressed).parse(args.ext))
    except FormattingException as e:
        log.error(e)
    except IOError as e:
        # A 'Broken pipe' IOError may occur when stdout is closed prematurely,
        # eg when using `fitsheader file.fits | head`. We let this pass.
        pass
