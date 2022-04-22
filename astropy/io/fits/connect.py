# Licensed under a 3-clause BSD style license - see LICENSE.rst


import os
import re
import warnings
from copy import deepcopy

import numpy as np

from astropy.io import registry as io_registry
from astropy import units as u
from astropy.table import Table, serialize, meta, Column, MaskedColumn
from astropy.time import Time
from astropy.utils.data_info import serialize_context_as
from astropy.utils.exceptions import (AstropyUserWarning,
                                      AstropyDeprecationWarning)
from astropy.utils.misc import NOT_OVERWRITING_MSG
from . import HDUList, TableHDU, BinTableHDU, GroupsHDU, append as fits_append
from .column import KEYWORD_NAMES, _fortran_to_python_format
from .convenience import table_to_hdu
from .hdu.hdulist import fitsopen as fits_open, FITS_SIGNATURE
from .util import first


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
    origin : str or readable file-like
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


def _decode_mixins(tbl):
    """Decode a Table ``tbl`` that has astropy Columns + appropriate meta-data into
    the corresponding table with mixin columns (as appropriate).
    """
    # If available read in __serialized_columns__ meta info which is stored
    # in FITS COMMENTS between two sentinels.
    try:
        i0 = tbl.meta['comments'].index('--BEGIN-ASTROPY-SERIALIZED-COLUMNS--')
        i1 = tbl.meta['comments'].index('--END-ASTROPY-SERIALIZED-COLUMNS--')
    except (ValueError, KeyError):
        return tbl

    # The YAML data are split into COMMENT cards, with lines longer than 70
    # characters being split with a continuation character \ (backslash).
    # Strip the backslashes and join together.
    continuation_line = False
    lines = []
    for line in tbl.meta['comments'][i0 + 1:i1]:
        if continuation_line:
            lines[-1] = lines[-1] + line[:70]
        else:
            lines.append(line[:70])
        continuation_line = len(line) == 71

    del tbl.meta['comments'][i0:i1 + 1]
    if not tbl.meta['comments']:
        del tbl.meta['comments']

    info = meta.get_header_from_yaml(lines)

    # Add serialized column information to table meta for use in constructing mixins
    tbl.meta['__serialized_columns__'] = info['meta']['__serialized_columns__']

    # Use the `datatype` attribute info to update column attributes that are
    # NOT already handled via standard FITS column keys (name, dtype, unit).
    for col in info['datatype']:
        for attr in ['description', 'meta']:
            if attr in col:
                setattr(tbl[col['name']].info, attr, col[attr])

    # Construct new table with mixins, using tbl.meta['__serialized_columns__']
    # as guidance.
    tbl = serialize._construct_mixins_from_columns(tbl)

    return tbl


def read_table_fits(input, hdu=None, astropy_native=False, memmap=False,
                    character_as_bytes=True, unit_parse_strict='warn',
                    mask_invalid=True):
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
    input : str or file-like or compatible `astropy.io.fits` HDU object
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
        When set to `True` then ``mask_invalid`` is set to `False` since the
        masking would cause loading the full data array.
    character_as_bytes : bool, optional
        If `True`, string columns are stored as Numpy byte arrays (dtype ``S``)
        and are converted on-the-fly to unicode strings when accessing
        individual elements. If you need to use Numpy unicode arrays (dtype
        ``U``) internally, you should set this to `False`, but note that this
        will use more memory. If set to `False`, string columns will not be
        memory-mapped even if ``memmap`` is `True`.
    unit_parse_strict : str, optional
        Behaviour when encountering invalid column units in the FITS header.
        Default is "warn", which will emit a ``UnitsWarning`` and create a
        :class:`~astropy.units.core.UnrecognizedUnit`.
        Values are the ones allowed by the ``parse_strict`` argument of
        :class:`~astropy.units.core.Unit`: ``raise``, ``warn`` and ``silent``.
    mask_invalid : bool, optional
        By default the code masks NaNs in float columns and empty strings in
        string columns. Set this parameter to `False` to avoid the performance
        penalty of doing this masking step. The masking is always deactivated
        when using ``memmap=True`` (see above).

    """

    if isinstance(input, HDUList):

        # Parse all table objects
        tables = dict()
        for ihdu, hdu_item in enumerate(input):
            if isinstance(hdu_item, (TableHDU, BinTableHDU, GroupsHDU)):
                tables[ihdu] = hdu_item

        if len(tables) > 1:
            if hdu is None:
                warnings.warn("hdu= was not specified but multiple tables"
                              " are present, reading in first available"
                              f" table (hdu={first(tables)})",
                              AstropyUserWarning)
                hdu = first(tables)

            # hdu might not be an integer, so we first need to convert it
            # to the correct HDU index
            hdu = input.index_of(hdu)

            if hdu in tables:
                table = tables[hdu]
            else:
                raise ValueError(f"No table found in hdu={hdu}")

        elif len(tables) == 1:
            if hdu is not None:
                msg = None
                try:
                    hdi = input.index_of(hdu)
                except KeyError:
                    msg = f"Specified hdu={hdu} not found"
                else:
                    if hdi >= len(input):
                        msg = f"Specified hdu={hdu} not found"
                    elif hdi not in tables:
                        msg = f"No table found in specified hdu={hdu}"
                if msg is not None:
                    warnings.warn(f"{msg}, reading in first available table "
                                  f"(hdu={first(tables)}) instead. This will"
                                  " result in an error in future versions!",
                                  AstropyDeprecationWarning)
            table = tables[first(tables)]

        else:
            raise ValueError("No table found")

    elif isinstance(input, (TableHDU, BinTableHDU, GroupsHDU)):

        table = input

    else:

        if memmap:
            # using memmap is not compatible with masking invalid value by
            # default so we deactivate the masking
            mask_invalid = False

        hdulist = fits_open(input, character_as_bytes=character_as_bytes,
                            memmap=memmap)

        try:
            return read_table_fits(
                hdulist, hdu=hdu,
                astropy_native=astropy_native,
                unit_parse_strict=unit_parse_strict,
                mask_invalid=mask_invalid,
            )
        finally:
            hdulist.close()

    # In the loop below we access the data using data[col.name] rather than
    # col.array to make sure that the data is scaled correctly if needed.
    data = table.data

    columns = []
    for col in data.columns:
        # Check if column is masked. Here, we make a guess based on the
        # presence of FITS mask values. For integer columns, this is simply
        # the null header, for float and complex, the presence of NaN, and for
        # string, empty strings.
        # Since Multi-element columns with dtypes such as '2f8' have a subdtype,
        # we should look up the type of column on that.
        masked = mask = False
        coltype = (col.dtype.subdtype[0].type if col.dtype.subdtype
                   else col.dtype.type)
        if col.null is not None:
            mask = data[col.name] == col.null
            # Return a MaskedColumn even if no elements are masked so
            # we roundtrip better.
            masked = True
        elif mask_invalid and issubclass(coltype, np.inexact):
            mask = np.isnan(data[col.name])
        elif mask_invalid and issubclass(coltype, np.character):
            mask = col.array == b''

        if masked or np.any(mask):
            column = MaskedColumn(data=data[col.name], name=col.name,
                                  mask=mask, copy=False)
        else:
            column = Column(data=data[col.name], name=col.name, copy=False)

        # Copy over units
        if col.unit is not None:
            column.unit = u.Unit(col.unit, format='fits', parse_strict=unit_parse_strict)

        # Copy over display format
        if col.disp is not None:
            column.format = _fortran_to_python_format(col.disp)

        columns.append(column)

    # Create Table object
    t = Table(columns, copy=False)

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

    # Decode any mixin columns that have been stored as standard Columns.
    t = _decode_mixins(t)

    return t


def _encode_mixins(tbl):
    """Encode a Table ``tbl`` that may have mixin columns to a Table with only
    astropy Columns + appropriate meta-data to allow subsequent decoding.
    """
    # Determine if information will be lost without serializing meta.  This is hardcoded
    # to the set difference between column info attributes and what FITS can store
    # natively (name, dtype, unit).  See _get_col_attributes() in table/meta.py for where
    # this comes from.
    info_lost = any(any(getattr(col.info, attr, None) not in (None, {})
                        for attr in ('description', 'meta'))
                    for col in tbl.itercols())

    # Convert the table to one with no mixins, only Column objects.  This adds
    # meta data which is extracted with meta.get_yaml_from_table.  This ignores
    # Time-subclass columns and leave them in the table so that the downstream
    # FITS Time handling does the right thing.

    with serialize_context_as('fits'):
        encode_tbl = serialize.represent_mixins_as_columns(
            tbl, exclude_classes=(Time,))

    # If the encoded table is unchanged then there were no mixins.  But if there
    # is column metadata (format, description, meta) that would be lost, then
    # still go through the serialized columns machinery.
    if encode_tbl is tbl and not info_lost:
        return tbl

    # Copy the meta dict if it was not copied by represent_mixins_as_columns.
    # We will modify .meta['comments'] below and we do not want to see these
    # comments in the input table.
    if encode_tbl is tbl:
        meta_copy = deepcopy(tbl.meta)
        encode_tbl = Table(tbl.columns, meta=meta_copy, copy=False)

    # Get the YAML serialization of information describing the table columns.
    # This is re-using ECSV code that combined existing table.meta with with
    # the extra __serialized_columns__ key.  For FITS the table.meta is handled
    # by the native FITS connect code, so don't include that in the YAML
    # output.
    ser_col = '__serialized_columns__'

    # encode_tbl might not have a __serialized_columns__ key if there were no mixins,
    # but machinery below expects it to be available, so just make an empty dict.
    encode_tbl.meta.setdefault(ser_col, {})

    tbl_meta_copy = encode_tbl.meta.copy()
    try:
        encode_tbl.meta = {ser_col: encode_tbl.meta[ser_col]}
        meta_yaml_lines = meta.get_yaml_from_table(encode_tbl)
    finally:
        encode_tbl.meta = tbl_meta_copy
    del encode_tbl.meta[ser_col]

    if 'comments' not in encode_tbl.meta:
        encode_tbl.meta['comments'] = []
    encode_tbl.meta['comments'].append('--BEGIN-ASTROPY-SERIALIZED-COLUMNS--')

    for line in meta_yaml_lines:
        if len(line) == 0:
            lines = ['']
        else:
            # Split line into 70 character chunks for COMMENT cards
            idxs = list(range(0, len(line) + 70, 70))
            lines = [line[i0:i1] + '\\' for i0, i1 in zip(idxs[:-1], idxs[1:])]
            lines[-1] = lines[-1][:-1]
        encode_tbl.meta['comments'].extend(lines)

    encode_tbl.meta['comments'].append('--END-ASTROPY-SERIALIZED-COLUMNS--')

    return encode_tbl


def write_table_fits(input, output, overwrite=False, append=False):
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
    append : bool
        Whether to append the table to an existing file
    """

    # Encode any mixin columns into standard Columns.
    input = _encode_mixins(input)

    table_hdu = table_to_hdu(input, character_as_bytes=True)

    # Check if output file already exists
    if isinstance(output, str) and os.path.exists(output):
        if overwrite:
            os.remove(output)
        elif not append:
            raise OSError(NOT_OVERWRITING_MSG.format(output))

    if append:
        # verify=False stops it reading and checking the existing file.
        fits_append(output, table_hdu.data, table_hdu.header, verify=False)
    else:
        table_hdu.writeto(output)


io_registry.register_reader('fits', Table, read_table_fits)
io_registry.register_writer('fits', Table, write_table_fits)
io_registry.register_identifier('fits', Table, is_fits)
