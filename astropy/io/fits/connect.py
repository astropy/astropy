# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import print_function

import os

from ...utils import OrderedDict
from ...table import io_registry, Table
from . import HDUList, TableHDU, BinTableHDU
from . import open as fits_open


def is_fits(origin, args, kwargs):
    """
    Determine whether `origin` is a FITS file.

    Parameters
    ----------
    origin : str or readable file-like object
        Path or file object containing a potential FITS file.

    Returns
    -------
    is_fits : bool
        Returns `True` if the given file is a FITS file.
    """
    if isinstance(origin, basestring):
        if origin.lower().endswith(('.fits', '.fits.gz', '.fit', '.fit.gz')):
            return True
        else:
            return False
    else:
        return False

def read_table_fits(input, hdu_id=None):
    """
    Read a Table object from an FITS file

    Parameters
    ----------
    input : str or `~astropy.io.fits.hdu.table.TableHDU` or `~astropy.io.fits.hdu.table.BinTableHDU` or `~astropy.io.fits.hdu.hdulist.HDUList`
        If a string, the filename to read the table from. If a
        :class:`~astropy.io.fits.hdu.table.TableHDU` or
        :class:`~astropy.io.fits.hdu.table.BinTableHDU` or
        :class:`~astropy.io.fits.hdu.hdulist.HDUList`, the object to extract
        the table from.
    hdu_id : str, optional
        The HDU to read the table from
    """
    if isinstance(input, basestring):
        input = fits_open(input)

    # Parse all table objects
    tables = OrderedDict()
    if isinstance(input, HDUList):
        for ihdu, hdu in enumerate(input):
            if isinstance(hdu, (TableHDU, BinTableHDU)):
                tables[ihdu] = hdu

        if len(tables) > 1:
            if hdu_id is None:
                raise ValueError(
                    "Multiple tables found: HDU id should be set via "
                    "the hdu= argument. The available tables HDUs are " +
                    ', '.join([str(x) for x in tables.keys()]))
            else:
                if hdu_id in tables:
                    table = tables[hdu_id]
                else:
                    raise ValueError(
                        "No tables with hdu_id={0} found".format(hdu_id))
        elif len(tables) == 1:
            table = tables[tables.keys()[0]]
        else:
            raise ValueError("No table found")

    # Convert to an astropy.table.Table object
    return Table(table.data)


def write_table_fits(input, output, overwrite=False):
    """
    Write a Table object to a FITS file

    Parameters
    ----------
    input : Table
        The table to write out.
    output : str
        The filename to write the table to.
    overwrite : bool
        Whether to overwrite any existing file without warning.
    """

    # Check if output file already exists
    if isinstance(output, basestring) and os.path.exists(output):
        if overwrite:
            os.remove(output)
        else:
            raise IOError("File exists: {0}".format(output))

    # Create a new HDU object
    table_hdu = BinTableHDU(input._data)

    # Write out file
    table_hdu.writeto(output)


io_registry.register_reader('fits', read_table_fits)
io_registry.register_writer('fits', write_table_fits)
io_registry.register_identifier('fits', is_fits)
