# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import print_function

import os
import re
import warnings

import numpy as np

from .. import registry as io_registry
from ... import log
from ... import units as u
from ...extern import six
from ...extern.six import string_types
from ...table import Table
from ...utils import OrderedDict
from ...utils.exceptions import AstropyUserWarning
from astropy.units.format.fits import UnitScaleError

from . import HDUList, TableHDU, BinTableHDU, GroupsHDU
from . import FITS_rec
from .hdu.hdulist import fitsopen as fits_open
from .util import first


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
                              " table (hdu={0})".format(first(tables)),
                              AstropyUserWarning)
                hdu = first(tables)

            # hdu might not be an integer, so we first need to convert it
            # to the correct HDU index
            hdu = input.index_of(hdu)

            if hdu in tables:
                table = tables[hdu]
            else:
                raise ValueError("No table found in hdu={0}".format(hdu))

        elif len(tables) == 1:
            table = tables[first(tables)]
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
            if col.null is not None:
                t[col.name].set_fill_value(col.null)
                t[col.name].mask[t[col.name] == col.null] = True

    # Copy over units
    for col in table.columns:
        if col.unit is not None:
            t[col.name].unit = u.Unit(
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

    # Tables with mixin columns are not supported
    if input.has_mixin_columns:
        mixin_names = [name for name, col in input.columns.items()
                       if not isinstance(col, input.ColumnClass)]
        raise ValueError('cannot write table with mixin column(s) {0} to FITS'
                         .format(mixin_names))

    # Check if output file already exists
    if isinstance(output, string_types) and os.path.exists(output):
        if overwrite:
            os.remove(output)
        else:
            raise IOError("File exists: {0}".format(output))

    # Create a new HDU object
    if input.masked:
        #float column's default mask value needs to be Nan
        for column in six.itervalues(input.columns):
            fill_value = column.get_fill_value()
            if column.dtype.kind == 'f' and np.allclose(fill_value, 1e20):
                column.set_fill_value(np.nan)

        fits_rec = FITS_rec.from_columns(np.array(input.filled()))
        table_hdu = BinTableHDU(fits_rec)
        for col in table_hdu.columns:
            # Binary FITS tables support TNULL *only* for integer data columns
            # TODO: Determine a schema for handling non-integer masked columns
            # in FITS (if at all possible)
            int_formats = ('B', 'I', 'J', 'K')
            if not (col.format in int_formats or
                    col.format.p_format in int_formats):
                continue

            # The astype is necessary because if the string column is less
            # than one character, the fill value will be N/A by default which
            # is too long, and so no values will get masked.
            fill_value = input[col.name].get_fill_value()

            col.null = fill_value.astype(input[col.name].dtype)
    else:
        fits_rec = FITS_rec.from_columns(np.array(input.filled()))
        table_hdu = BinTableHDU(fits_rec)

    # Set units for output HDU
    for col in table_hdu.columns:
        unit = input[col.name].unit
        if unit is not None:
            try:
                col.unit = unit.to_string(format='fits')
            except UnitScaleError:
                scale = unit.scale
                raise UnitScaleError(
                    "The column '{0}' could not be stored in FITS format "
                    "because it has a scale '({1})' that "
                    "is not recognized by the FITS standard. Either scale "
                    "the data or change the units.".format(col.name, str(scale)))
            except ValueError:
                warnings.warn(
                    "The unit '{0}' could not be saved to FITS format".format(
                        unit.to_string()), AstropyUserWarning)

    for key, value in input.meta.items():

        if is_column_keyword(key.upper()) or key.upper() in REMOVE_KEYWORDS:

            warnings.warn(
                "Meta-data keyword {0} will be ignored since it conflicts "
                "with a FITS reserved keyword".format(key), AstropyUserWarning)

        if isinstance(value, list):
            for item in value:
                try:
                    table_hdu.header.append((key, item))
                except ValueError:
                    warnings.warn(
                        "Attribute `{0}` of type {1} cannot be written to "
                        "FITS files - skipping".format(key, type(value)),
                        AstropyUserWarning)
        else:
            try:
                table_hdu.header[key] = value
            except ValueError:
                warnings.warn(
                    "Attribute `{0}` of type {1} cannot be written to FITS "
                    "files - skipping".format(key, type(value)),
                    AstropyUserWarning)

    # Write out file
    table_hdu.writeto(output)

io_registry.register_reader('fits', Table, read_table_fits)
io_registry.register_writer('fits', Table, write_table_fits)
io_registry.register_identifier('fits', Table, is_fits)
