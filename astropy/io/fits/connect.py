# Licensed under a 3-clause BSD style license - see LICENSE.rst


import os
import re
import warnings
from collections import OrderedDict

from .. import registry as io_registry
from ... import units as u
from ...table import Table, serialize, meta, Column, MaskedColumn
from ...table.table import has_info_class
from ...time import Time
from ...utils.exceptions import AstropyUserWarning
from ...utils.data_info import MixinInfo, serialize_context_as
from . import HDUList, TableHDU, BinTableHDU, GroupsHDU
from .column import KEYWORD_NAMES, ASCII_DEFAULT_WIDTHS
from .convenience import table_to_hdu
from .hdu.hdulist import fitsopen as fits_open
from .util import first
from .verify import VerifyError, VerifyWarning


# FITS file signature as per RFC 4047
FITS_SIGNATURE = (b"\x53\x49\x4d\x50\x4c\x45\x20\x20\x3d\x20\x20\x20\x20\x20"
                  b"\x20\x20\x20\x20\x20\x20\x20\x20\x20\x20\x20\x20\x20\x20"
                  b"\x20\x54")

# Keywords to remove for all tables that are read in
REMOVE_KEYWORDS = ['XTENSION', 'BITPIX', 'NAXIS', 'NAXIS1', 'NAXIS2',
                   'PCOUNT', 'GCOUNT', 'TFIELDS', 'THEAP']

# Column-specific keywords regex
COLUMN_KEYWORD_REGEXP = '(' + '|'.join(KEYWORD_NAMES) + ')[0-9]+'

# TDISPn for both ASCII and Binary tables
TDISP_RE_DICT = {}
TDISP_RE_DICT['F'] = re.compile(r'(?:(?P<formatc>[F])(?:(?P<width>[0-9]+)\.{1}'
                                r'(?P<precision>[0-9])+)+)|')
TDISP_RE_DICT['A'] = TDISP_RE_DICT['L'] = \
    re.compile(r'(?:(?P<formatc>[AL])(?P<width>[0-9]+)+)|')
TDISP_RE_DICT['I'] = TDISP_RE_DICT['B'] = \
    TDISP_RE_DICT['O'] = TDISP_RE_DICT['Z'] =  \
    re.compile(r'(?:(?P<formatc>[IBOZ])(?:(?P<width>[0-9]+)'
               r'(?:\.{0,1}(?P<precision>[0-9]+))?))|')
TDISP_RE_DICT['E'] = TDISP_RE_DICT['G'] = \
    TDISP_RE_DICT['D'] = \
    re.compile(r'(?:(?P<formatc>[EGD])(?:(?P<width>[0-9]+)\.'
               r'(?P<precision>[0-9]+))+)'
               r'(?:E{0,1}(?P<exponential>[0-9]+)?)|')
TDISP_RE_DICT['EN'] = TDISP_RE_DICT['ES'] = \
    re.compile(r'(?:(?P<formatc>E[NS])(?:(?P<width>[0-9]+)\.{1}'
               r'(?P<precision>[0-9])+)+)')

# mapping from TDISP format to python format
# A: Character
# L: Logical (Boolean)
# I: 16-bit Integer
#    Can't predefine zero padding and space padding before hand without
#    knowing the value being formatted, so grabbing precision and using that
#    to zero pad, ignoring width. Same with B, O, and Z
# B: Binary Integer
# O: Octal Integer
# Z: Hexadecimal Integer
# F: Float (64-bit; fixed decimal notation)
# EN: Float (engineering fortran format, exponential multiple of thee
# ES: Float (scientific, same as EN but non-zero leading digit
# E: Float, exponential notation
#    Can't get exponential restriction to work without knowing value
#    before hand, so just using width and precision, same with D, G, EN, and
#    ES formats
# D: Double-precision Floating Point with exponential
#    (E but for double precision)
# G: Double-precision Floating Point, may or may not show exponent
TDISP_FMT_DICT = {'I' : '{{:{width}d}}',
                  'B' : '{{:{width}b}}',
                  'O' : '{{:{width}o}}',
                  'Z' : '{{:{width}x}}',
                  'F' : '{{:{width}.{precision}f}}',
                  'G' : '{{:{width}.{precision}g}}'}
TDISP_FMT_DICT['A'] = TDISP_FMT_DICT['L'] = '{{:>{width}}}'
TDISP_FMT_DICT['E'] = TDISP_FMT_DICT['D'] =  \
    TDISP_FMT_DICT['EN'] = TDISP_FMT_DICT['ES'] ='{{:{width}.{precision}e}}'

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
        for attr in ['format', 'description', 'meta']:
            if attr in col:
                setattr(tbl[col['name']].info, attr, col[attr])

    # Construct new table with mixins, using tbl.meta['__serialized_columns__']
    # as guidance.
    tbl = serialize._construct_mixins_from_columns(tbl)

    return tbl


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

        # Copy over display format
        if col.disp is not None:
            column.format = _fortran_to_python_format(col.disp)

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
                        for attr in ('format', 'description', 'meta'))
                    for col in tbl.itercols())

    # If PyYAML is not available then check to see if there are any mixin cols
    # that *require* YAML serialization.  FITS already has support for Time,
    # Quantity, so if those are the only mixins the proceed without doing the
    # YAML bit, for backward compatibility (i.e. not requiring YAML to write
    # Time or Quantity).  In this case other mixin column meta (e.g.
    # description or meta) will be silently dropped, consistent with astropy <=
    # 2.0 behavior.
    try:
        import yaml
    except ImportError:
        for col in tbl.itercols():
            if (has_info_class(col, MixinInfo) and
                    col.__class__ not in (u.Quantity, Time)):
                raise TypeError("cannot write type {} column '{}' "
                                "to FITS without PyYAML installed."
                                .format(col.__class__.__name__, col.info.name))
        else:
            if info_lost:
                warnings.warn("table contains column(s) with defined 'format',"
                              " 'description', or 'meta' info attributes. These"
                              " will be dropped unless you install PyYAML.",
                              AstropyUserWarning)
            return tbl

    # Convert the table to one with no mixins, only Column objects.  This adds
    # meta data which is extracted with meta.get_yaml_from_table.  This ignores
    # Time-subclass columns and leave them in the table so that the downstream
    # FITS Time handling does the right thing.

    with serialize_context_as('fits'):
        encode_tbl = serialize._represent_mixins_as_columns(
            tbl, exclude_classes=(Time,))

    # If the encoded table is unchanged then there were no mixins.  But if there
    # is column metadata (format, description, meta) that would be lost, then
    # still go through the serialized columns machinery.
    if encode_tbl is tbl and not info_lost:
        return tbl

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
        # Split line into 70 character chunks for COMMENT cards
        idxs = list(range(0, len(line) + 70, 70))
        lines = [line[i0:i1] + '\\' for i0, i1 in zip(idxs[:-1], idxs[1:])]
        lines[-1] = lines[-1][:-1]
        encode_tbl.meta['comments'].extend(lines)

    encode_tbl.meta['comments'].append('--END-ASTROPY-SERIALIZED-COLUMNS--')

    return encode_tbl


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

    # Encode any mixin columns into standard Columns.
    input = _encode_mixins(input)

    table_hdu = table_to_hdu(input, character_as_bytes=True)

    # Check if output file already exists
    if isinstance(output, str) and os.path.exists(output):
        if overwrite:
            os.remove(output)
        else:
            raise OSError("File exists: {0}".format(output))

    table_hdu.writeto(output)


io_registry.register_reader('fits', Table, read_table_fits)
io_registry.register_writer('fits', Table, write_table_fits)
io_registry.register_identifier('fits', Table, is_fits)


def _parse_tdisp_format(tdisp):
    """
    Parse the ``TDISPn`` keywords for ASCII and binary tables into a
    ``(format, width, precision, exponential)`` tuple (the TDISP values
    for ASCII and binary are identical except for 'Lw',
    which is only present in BINTABLE extensions

    Parameters
    ----------
    tdisp: str
        TDISPn FITS Header keyword.  Used to specify display formatting.

    Returns
    -------
    formatc: str
        The format characters from TDISPn
    width: str
        The width int value from TDISPn
    precision: str
        The precision int value from TIDISPn
    exponential: str
        The exponential int value from TDISPn

    """

    # Use appropriate regex for format type
    tdisp = tdisp.strip()
    fmt_key = tdisp[0] if tdisp[0] !='E' or tdisp[1] not in 'NS' else tdisp[:2]
    try:
        tdisp_re = TDISP_RE_DICT[fmt_key]
    except KeyError:
        raise VerifyError('Format {} is not recognized.'.format(tdisp))


    match = tdisp_re.match(tdisp.strip())
    if not match or match.group('formatc') is None:
        raise VerifyError('Format {} is not recognized.'.format(tdisp))

    formatc = match.group('formatc')
    width = match.group('width')
    precision = None
    exponential = None

    # Some formats have precision and exponential
    if tdisp[0] in ('I', 'B', 'O', 'Z', 'F', 'E', 'G', 'D'):
        precision = match.group('precision')
        if precision is None:
            precision = 1
    if tdisp[0] in ('E', 'D', 'G') and tdisp[1] not in ('N', 'S'):
        exponential = match.group('exponential')
        if exponential is None:
            exponential = 1

    # Once parsed, check format dict to do conversion to a formatting string
    return formatc, width, precision, exponential


def _fortran_to_python_format(tdisp):
    """
    Turn the TDISPn fortran format pieces into a final Python format string.
    See the format_type definitions above the TDISP_FMT_DICT. If codes is
    changed to take advantage of the exponential specification, will need to
    add it as another input parameter.

    Parameters
    ----------
    tdisp: str
        TDISPn FITS Header keyword.  Used to specify display formatting.

    Returns
    -------
    format_string: str
        The TDISPn keyword string translated into a Python format string.
    """
    format_type, width, precision, exponential = _parse_tdisp_format(tdisp)

    try:
        fmt = TDISP_FMT_DICT[format_type]
        return fmt.format(width=width, precision=precision)

    except KeyError:
        raise VerifyError('Format {} is not recognized.'.format(format_type))


def python_to_tdisp(format_string, logical_dtype = False):
    """
    Turn the Python format string to a TDISP FITS compliant format string. Not
    all formats convert. these will cause a Warning and return None.

    Parameters
    ----------
    format_string: str
        TDISPn FITS Header keyword.  Used to specify display formatting.
    logical_dtype: bool
        True is this format type should be a logical type, 'L'. Needs special
        handeling.

    Returns
    -------
    tdsip_string: str
        The TDISPn keyword string translated into a Python format string.
    """

    fmt_to_tdisp = {'a': 'A', 's': 'A', 'd': 'I', 'b': 'B', 'o': 'O', 'x': 'Z',
                    'X': 'Z', 'f': 'F', 'F': 'F', 'g': 'G', 'G': 'G', 'e': 'E',
                    'E': 'E'}

    if format_string in [None, "", "{}"]:
        return None

    # Strip out extra format characters that aren't a type or a width/precision
    if format_string[0] == '{' and format_string != "{}":
        fmt_str = format_string.lstrip("{:").rstrip('}')
    elif format_string[0] == '%':
            fmt_str = format_string.lstrip("%")
    else:
        fmt_str = format_string

    precision, sep = '', ''

    # Character format, only translate right aligned, and don't take zero fills
    if fmt_str[-1].isdigit() and fmt_str[0] == '>' and fmt_str[1] != '0':
        ftype = fmt_to_tdisp['a']
        width = fmt_str[1:]

    elif fmt_str[-1] == 's' and fmt_str != 's':
        ftype = fmt_to_tdisp['a']
        width = fmt_str[:-1].lstrip('0')

    # Number formats, don't take zero fills
    elif fmt_str[-1].isalpha() and len(fmt_str) > 1 and fmt_str[0] != '0':
        ftype = fmt_to_tdisp[fmt_str[-1]]
        fmt_str = fmt_str[:-1]

        # If format has a "." split out the width and precision
        if '.' in fmt_str:
            width, precision = fmt_str.split('.')
            sep = '.'
            if width == "":
                ascii_key = ftype if ftype != 'G' else 'F'
                width = str(int(precision) + (ASCII_DEFAULT_WIDTHS[ascii_key][0] -
                                     ASCII_DEFAULT_WIDTHS[ascii_key][1]))
        # Otherwise we just have a width
        else:
            width = fmt_str

    else:
        warnings.warn('Format {} cannot be mapped to the accepted '
                      'TDISPn keyword values.  Format will not be '
                      'moved into TDISPn keyword.'.format(format_string),
                      AstropyUserWarning)
        return None

    # Catch logical data type, set the format type back to L in this case
    if logical_dtype:
        ftype = 'L'

    return ftype + width + sep + precision
