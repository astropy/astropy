# Licensed under a 3-clause BSD style license - see LICENSE.rst


import os
import re
import warnings
from collections import OrderedDict

from .. import registry as io_registry
from ... import units as u
from ...table import Table, Column, MaskedColumn
from ...utils.exceptions import AstropyUserWarning
from . import HDUList, TableHDU, BinTableHDU, GroupsHDU
from .column import KEYWORD_NAMES
from .convenience import table_to_hdu
from .hdu.hdulist import fitsopen as fits_open
from .util import first


# FITS file signature as per RFC 4047
FITS_SIGNATURE = (b"\x53\x49\x4d\x50\x4c\x45\x20\x20\x3d\x20\x20\x20\x20\x20"
                  b"\x20\x20\x20\x20\x20\x20\x20\x20\x20\x20\x20\x20\x20\x20"
                  b"\x20\x54")

# Keywords to remove for all tables that are read in
REMOVE_KEYWORDS = ['XTENSION', 'BITPIX', 'NAXIS', 'NAXIS1', 'NAXIS2',
                   'PCOUNT', 'GCOUNT', 'TFIELDS', 'THEAP']

# Column-specific keywords regex
COLUMN_KEYWORD_REGEXP = '(' + '|'.join(KEYWORD_NAMES) + ')[0-9]+'


def is_column_keyword(keyword):
    return re.match(COLUMN_KEYWORD_REGEXP, keyword) is not None


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
        if filepath.lower().endswith(('.fits', '.fits.gz', '.fit', '.fit.gz',
                                      '.fts', '.fts.gz')):
            return True
    elif isinstance(args[0], (HDUList, TableHDU, BinTableHDU, GroupsHDU)):
        return True
    else:
        return False


def read_table_fits(input, hdu=None, astropy_native=False, memmap=False,
                    character_as_bytes=True):
    """
    Read a Table object from an FITS file

    If the ``astropy_native`` argument is ``True``, then input FITS columns
    which are representations of an astropy core object will be converted to
    that class and stored in the ``Table`` as "mixin columns".  Currently this
    is limited to FITS columns which adhere to the FITS Time standard, in which
    case they will be converted to a `~astropy.time.Time` column in the output
    table.

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
    astropy_native : bool, optional
        Read in FITS columns as native astropy objects where possible instead
        of standard Table Column objects. Default is False.
    memmap : bool, optional
        Whether to use memory mapping, which accesses data on disk as needed. If
        you are only accessing part of the data, this is often more efficient.
        If you want to access all the values in the table, and you are able to
        fit the table in memory, you may be better off leaving memory mapping
        off. However, if your table would not fit in memory, you should set this
        to `True`.
    character_as_bytes : bool, optional
        If `True`, string columns are stored as Numpy byte arrays (dtype ``S``)
        and are converted on-the-fly to unicode strings when accessing
        individual elements. If you need to use Numpy unicode arrays (dtype
        ``U``) internally, you should set this to `False`, but note that this
        will use more memory. If set to `False`, string columns will not be
        memory-mapped even if ``memmap`` is `True`.
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

        hdulist = fits_open(input, character_as_bytes=character_as_bytes,
                            memmap=memmap)

        try:
            return read_table_fits(hdulist, hdu=hdu,
                                   astropy_native=astropy_native)
        finally:
            hdulist.close()

    # Check if table is masked
    masked = any(col.null is not None for col in table.columns)

    # TODO: in future, it may make more sense to do this column-by-column,
    # rather than via the structured array.

    # In the loop below we access the data using data[col.name] rather than
    # col.array to make sure that the data is scaled correctly if needed.
    data = table.data

    columns = []
    for col in data.columns:

        # Set column data
        if masked:
            column = MaskedColumn(data=data[col.name], name=col.name, copy=False)
            if col.null is not None:
                column.set_fill_value(col.null)
                column.mask[column.data == col.null] = True
        else:
            column = Column(data=data[col.name], name=col.name, copy=False)

        # Copy over units
        if col.unit is not None:
            column.unit = u.Unit(col.unit, format='fits', parse_strict='silent')

        columns.append(column)

    # Create Table object
    t = Table(columns, masked=masked, copy=False)

    # TODO: deal properly with unsigned integers

    hdr = table.header
    if astropy_native:
        # Avoid circular imports, and also only import if necessary.
        from .fitstime import fits_to_time
        hdr = fits_to_time(hdr, t)

    for key, value, comment in hdr.cards:

        if key in ['COMMENT', 'HISTORY']:
            # Convert to io.ascii format
            if key == 'COMMENT':
                key = 'comments'

            if key in t.meta:
                t.meta[key].append(value)
            else:
                t.meta[key] = [value]

        elif key in t.meta:  # key is duplicate

            if isinstance(t.meta[key], list):
                t.meta[key].append(value)
            else:
                t.meta[key] = [t.meta[key], value]

        elif is_column_keyword(key) or key in REMOVE_KEYWORDS:

            pass

        else:

            t.meta[key] = value

    # TODO: implement masking

    return t


def write_table_fits(input, output, overwrite=False, character_as_bytes=True):
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

    table_hdu = table_to_hdu(input, character_as_bytes=character_as_bytes)

    # Check if output file already exists
    if isinstance(output, str) and os.path.exists(output):
        if overwrite:
            os.remove(output)
        else:
            raise IOError("File exists: {0}".format(output))

    table_hdu.writeto(output)


io_registry.register_reader('fits', Table, read_table_fits)
io_registry.register_writer('fits', Table, write_table_fits)
io_registry.register_identifier('fits', Table, is_fits)
