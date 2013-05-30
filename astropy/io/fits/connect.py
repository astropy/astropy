# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import print_function

import os
import re
import warnings

import numpy as np

from ...utils import OrderedDict
from .. import registry as io_registry
from ...table import Table
from ... import log
from ... import units as u

from . import HDUList, TableHDU, BinTableHDU, GroupsHDU
from .hdu.hdulist import fitsopen as fits_open


# FITS file signature as per RFC 4047
FITS_SIGNATURE = (b"\x53\x49\x4d\x50\x4c\x45\x20\x20\x3d\x20\x20\x20\x20\x20"
                  b"\x20\x20\x20\x20\x20\x20\x20\x20\x20\x20\x20\x20\x20\x20"
                  b"\x20\x54")

# Keywords to remove for all tables that are read in
REMOVE_KEYWORDS = ['XTENSION', 'BITPIX', 'NAXIS', 'NAXIS1', 'NAXIS2',
                   'PCOUNT', 'GCOUNT', 'TFIELDS']

# Column-specific keywords
COLUMN_KEYWORDS = ['TFORM[0-9]+',
                   'TBCOL[0-9]+',
                   'TSCAL[0-9]+',
                   'TZERO[0-9]+',
                   'TNULL[0-9]+',
                   'TTYPE[0-9]+',
                   'TUNIT[0-9]+',
                   'TDISP[0-9]+',
                   'TDIM[0-9]+',
                   'THEAP']


def is_column_keyword(keyword):
    for c in COLUMN_KEYWORDS:
        if re.match(c, keyword) is not None:
            return True
    return False


def is_fits(origin, filepath, fileobj, *args, **kwargs):
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
    if fileobj is not None:
        pos = fileobj.tell()
        sig = fileobj.read(30)
        fileobj.seek(pos)
        return sig == FITS_SIGNATURE
    elif filepath is not None:
        if filepath.lower().endswith(('.fits', '.fits.gz', '.fit', '.fit.gz')):
            return True
    elif isinstance(args[0], (HDUList, TableHDU, BinTableHDU, GroupsHDU)):
        return True
    else:
        return False


def read_table_fits(input, hdu=None):
    """
    Read a Table object from an FITS file

    Parameters
    ----------
    input : str or file-like object or compatible `astropy.io.fits` HDU object
        If a string, the filename to read the table from. If a file object, or
        a compatible HDU object, the object to extract the table from. The
        following `astropy.io.fits` HDU objects can be used as input:
        - :class:`~astropy.io.fits.hdu.table.TableHDU`
        - :class:`~astropy.io.fits.hdu.table.BinTableHDU`
        - :class:`~astropy.io.fits.hdu.table.GroupsHDU`
        - :class:`~astropy.io.fits.hdu.hdulist.HDUList`
    hdu : int or str, optional
        The HDU to read the table from.
    """

    if isinstance(input, HDUList):

        # Parse all table objects
        tables = OrderedDict()
        for ihdu, hdu_item in enumerate(input):
            if isinstance(hdu_item, (TableHDU, BinTableHDU, GroupsHDU)):
                tables[ihdu] = hdu_item

        if len(tables) > 1:

            if hdu is None:
                warnings.warn("hdu= was not specified but multiple tables"
                              " are present, reading in first available"
                              " table (hdu={0})".format(tables.keys()[0]))
                hdu = tables.keys()[0]

            # hdu might not be an integer, so we first need to convert it
            # to the correct HDU index
            hdu = input.index_of(hdu)

            if hdu in tables:
                table = tables[hdu]
            else:
                raise ValueError("No table found in hdu={0}".format(hdu))

        elif len(tables) == 1:
            table = tables[tables.keys()[0]]
        else:
            raise ValueError("No table found")

    elif isinstance(input, (TableHDU, BinTableHDU, GroupsHDU)):

        table = input

    else:

        hdulist = fits_open(input)

        try:
            return read_table_fits(hdulist, hdu=hdu)
        finally:
            hdulist.close()

    # Check if table is masked
    masked = False
    for col in table.columns:
        if col.null is not None:
            masked = True
            break

    # Convert to an astropy.table.Table object
    t = Table(table.data, masked=masked)

    # Copy over null values if needed
    if masked:
        for col in table.columns:
            t[col.name].set_fill_value(col.null)
            t[col.name].mask[t[col.name] == col.null] = True

    # Copy over units
    for col in table.columns:
        if col.unit is not None:
            t[col.name].units = u.Unit(
                col.unit, format='fits', parse_strict='warn')

    # TODO: deal properly with unsigned integers

    for key, value, comment in table.header.cards:

        if key in ['COMMENT', 'HISTORY']:
            if key in t.meta:
                t.meta[key].append(value)
            else:
                t.meta[key] = [value]

        elif key in t.meta:  # key is duplicate

            if isinstance(t.meta[key], list):
                t.meta[key].append(value)
            else:
                t.meta[key] = [t.meta[key], value]

        elif (is_column_keyword(key.upper()) or
              key.upper() in REMOVE_KEYWORDS):

            pass

        else:

            t.meta[key] = value

    # TODO: implement masking

    return t


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
    if input.masked:
        table_hdu = BinTableHDU(np.array(input.filled()))
        for col in table_hdu.columns:
            # The astype is necessary because if the string column is less
            # than one character, the fill value will be N/A by default which
            # is too long, and so no values will get masked.
            fill_value = input[col.name].get_fill_value()
            col.null = fill_value.astype(input[col.name].dtype)
    else:
        table_hdu = BinTableHDU(np.array(input))

    # Set units for output HDU
    for col in table_hdu.columns:
        if input[col.name].units is not None:
            col.unit = input[col.name].units.to_string(format='fits')

    for key, value in input.meta.items():

        if is_column_keyword(key.upper()) or key.upper() in REMOVE_KEYWORDS:

            log.warn("Meta-data keyword {0} will be ignored since it "
                     "conflicts with a FITS reserved keyword".format(key))

        if isinstance(value, list):
            for item in value:
                try:
                    table_hdu.header.append((key, item))
                except ValueError:
                    log.warn("Attribute `{0}` of type {1} cannot be written "
                             "to FITS files - skipping".format(key,
                                                               type(value)))
        else:
            try:
                table_hdu.header[key] = value
            except ValueError:
                log.warn("Attribute `{0}` of type {1} cannot be written to "
                         "FITS files - skipping".format(key, type(value)))

    # Write out file
    table_hdu.writeto(output)

io_registry.register_reader('fits', Table, read_table_fits)
io_registry.register_writer('fits', Table, write_table_fits)
io_registry.register_identifier('fits', Table, is_fits)
