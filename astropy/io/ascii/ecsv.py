# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Define the Enhanced Character-Separated-Values (ECSV) which allows for reading and
writing all the meta data associated with an astropy Table object.
"""

import re
from collections import OrderedDict
import warnings
import json

import numpy as np

from . import core, basic
from astropy.table import meta, serialize
from astropy.utils.data_info import serialize_context_as
from astropy.utils.exceptions import AstropyUserWarning
from astropy.io.ascii.core import convert_numpy

ECSV_VERSION = '1.0'
DELIMITERS = (' ', ',')
ECSV_DATATYPES = (
    'bool', 'int8', 'int16', 'int32', 'int64', 'uint8', 'uint16',
    'uint32', 'uint64', 'float16', 'float32', 'float64',
    'float128', 'string')  # Raise warning if not one of these standard dtypes


class InvalidEcsvDatatypeWarning(AstropyUserWarning):
    """
    ECSV specific Astropy warning class.
    """


class EcsvHeader(basic.BasicHeader):
    """Header class for which the column definition line starts with the
    comment character.  See the :class:`CommentedHeader` class  for an example.
    """

    def process_lines(self, lines):
        """Return only non-blank lines that start with the comment regexp.  For these
        lines strip out the matching characters and leading/trailing whitespace."""
        re_comment = re.compile(self.comment)
        for line in lines:
            line = line.strip()
            if not line:
                continue
            match = re_comment.match(line)
            if match:
                out = line[match.end():]
                if out:
                    yield out
            else:
                # Stop iterating on first failed match for a non-blank line
                return

    def write(self, lines):
        """
        Write header information in the ECSV ASCII format.

        This function is called at the point when preprocessing has been done to
        convert the input table columns to `self.cols` which is a list of
        `astropy.io.ascii.core.Column` objects. In particular `col.str_vals`
        is available for each column with the string representation of each
        column item for output.

        This format starts with a delimiter separated list of the column names
        in order to make this format readable by humans and simple csv-type
        readers. It then encodes the full table meta and column attributes and
        meta as YAML and pretty-prints this in the header.  Finally the
        delimited column names are repeated again, for humans and readers that
        look for the *last* comment line as defining the column names.
        """
        if self.splitter.delimiter not in DELIMITERS:
            raise ValueError('only space and comma are allowed for delimiter in ECSV format')

        # Now assemble the header dict that will be serialized by the YAML dumper
        header = {'cols': self.cols, 'schema': 'astropy-2.0'}

        if self.table_meta:
            header['meta'] = self.table_meta

        # Set the delimiter only for the non-default option(s)
        if self.splitter.delimiter != ' ':
            header['delimiter'] = self.splitter.delimiter

        header_yaml_lines = ([f'%ECSV {ECSV_VERSION}',
                              '---']
                             + meta.get_yaml_from_header(header))

        lines.extend([self.write_comment + line for line in header_yaml_lines])
        lines.append(self.splitter.join([x.info.name for x in self.cols]))

    def write_comments(self, lines, meta):
        """
        WRITE: Override the default write_comments to do nothing since this is handled
        in the custom write method.
        """
        pass

    def update_meta(self, lines, meta):
        """
        READ: Override the default update_meta to do nothing.  This process is done
        in get_cols() for this reader.
        """
        pass

    def get_cols(self, lines):
        """
        READ: Initialize the header Column objects from the table ``lines``.

        Parameters
        ----------
        lines : list
            List of table lines

        """
        # Cache a copy of the original input lines before processing below
        raw_lines = lines

        # Extract non-blank comment (header) lines with comment character stripped
        lines = list(self.process_lines(lines))

        # Validate that this is a ECSV file
        ecsv_header_re = r"""%ECSV [ ]
                             (?P<major> \d+)
                             \. (?P<minor> \d+)
                             \.? (?P<bugfix> \d+)? $"""

        no_header_msg = ('ECSV header line like "# %ECSV <version>" not found as first line.'
                         '  This is required for a ECSV file.')

        if not lines:
            raise core.InconsistentTableError(no_header_msg)

        match = re.match(ecsv_header_re, lines[0].strip(), re.VERBOSE)
        if not match:
            raise core.InconsistentTableError(no_header_msg)

        try:
            header = meta.get_header_from_yaml(lines)
        except meta.YamlParseError:
            raise core.InconsistentTableError('unable to parse yaml in meta header')

        if 'meta' in header:
            self.table_meta = header['meta']

        if 'delimiter' in header:
            delimiter = header['delimiter']
            if delimiter not in DELIMITERS:
                raise ValueError('only space and comma are allowed for delimiter in ECSV format')
            self.splitter.delimiter = delimiter
            self.data.splitter.delimiter = delimiter

        # Create the list of io.ascii column objects from `header`
        header_cols = OrderedDict((x['name'], x) for x in header['datatype'])
        self.names = [x['name'] for x in header['datatype']]

        # Read the first non-commented line of table and split to get the CSV
        # header column names.  This is essentially what the Basic reader does.
        header_line = next(super().process_lines(raw_lines))
        header_names = next(self.splitter([header_line]))

        # Check for consistency of the ECSV vs. CSV header column names
        if header_names != self.names:
            raise core.InconsistentTableError('column names from ECSV header {} do not '
                                              'match names from header line of CSV data {}'
                                              .format(self.names, header_names))

        # BaseHeader method to create self.cols, which is a list of
        # io.ascii.core.Column objects (*not* Table Column objects).
        self._set_cols_from_names()

        # Transfer attributes from the column descriptor stored in the input
        # header YAML metadata to the new columns to create this table.
        for col in self.cols:
            for attr in ('description', 'format', 'unit', 'meta', 'subtype'):
                if attr in header_cols[col.name]:
                    setattr(col, attr, header_cols[col.name][attr])

            col.dtype = header_cols[col.name]['datatype']
            # Warn if col dtype is not a valid ECSV datatype, but allow reading for
            # back-compatibility with existing older files that have numpy datatypes
            # like datetime64 or object or python str, which are not in the ECSV standard.
            if col.dtype not in ECSV_DATATYPES:
                msg = (f'unexpected datatype {col.dtype!r} of column {col.name!r} '
                       f'is not in allowed ECSV datatypes {ECSV_DATATYPES}. '
                       'Using anyway as a numpy dtype but beware since unexpected '
                       'results are possible.')
                warnings.warn(msg, category=InvalidEcsvDatatypeWarning)

            # Subtype is written like "int64[2,null]" and we want to split this
            # out to "int64" and [2, None].
            subtype = col.subtype
            if subtype and '[' in subtype:
                idx = subtype.index('[')
                col.subtype = subtype[:idx]
                col.shape = json.loads(subtype[idx:])

            # Convert ECSV "string" to numpy "str"
            for attr in ('dtype', 'subtype'):
                if getattr(col, attr) == 'string':
                    setattr(col, attr, 'str')

            # ECSV subtype of 'json' maps to numpy 'object' dtype
            if col.subtype == 'json':
                col.subtype = 'object'


def _check_dtype_is_str(col):
    if col.dtype != 'str':
        raise ValueError(f'datatype of column {col.name!r} must be "string"')


class EcsvOutputter(core.TableOutputter):
    """
    After reading the input lines and processing, convert the Reader columns
    and metadata to an astropy.table.Table object.  This overrides the default
    converters to be an empty list because there is no "guessing" of the
    conversion function.
    """
    default_converters = []

    def __call__(self, cols, meta):
        # Convert to a Table with all plain Column subclass columns
        out = super().__call__(cols, meta)

        # If mixin columns exist (based on the special '__mixin_columns__'
        # key in the table ``meta``), then use that information to construct
        # appropriate mixin columns and remove the original data columns.
        # If no __mixin_columns__ exists then this function just passes back
        # the input table.
        out = serialize._construct_mixins_from_columns(out)

        return out

    def _convert_vals(self, cols):
        """READ: Convert str_vals in `cols` to final arrays with correct dtypes.

        This is adapted from ``BaseOutputter._convert_vals``. In the case of ECSV
        there is no guessing and all types are known in advance. A big change
        is handling the possibility of JSON-encoded values, both unstructured
        object data and structured values that may contain masked data.
        """
        for col in cols:
            try:
                # 1-d or N-d object columns are serialized as JSON.
                if col.subtype == 'object':
                    _check_dtype_is_str(col)
                    col_vals = [json.loads(val) for val in col.str_vals]
                    col.data = np.empty([len(col_vals)] + col.shape, dtype=object)
                    col.data[...] = col_vals

                # Variable length arrays with shape (n, m, ..., *) for fixed
                # n, m, .. and variable in last axis. Masked values here are
                # not currently supported.
                elif col.shape and col.shape[-1] is None:
                    _check_dtype_is_str(col)

                    # Empty (blank) values in original ECSV are changed to "0"
                    # in str_vals with corresponding col.mask being created and
                    # set accordingly. Instead use an empty list here.
                    if hasattr(col, 'mask'):
                        for idx in np.nonzero(col.mask)[0]:
                            col.str_vals[idx] = '[]'

                    # Remake as a 1-d object column of numpy ndarrays or
                    # MaskedArray using the datatype specified in the ECSV file.
                    col_vals = []
                    for str_val in col.str_vals:
                        obj_val = json.loads(str_val)  # list or nested lists
                        try:
                            arr_val = np.array(obj_val, dtype=col.subtype)
                        except TypeError:
                            # obj_val has entries that are inconsistent with
                            # dtype. For a valid ECSV file the only possibility
                            # is None values (indicating missing values).
                            data = np.array(obj_val, dtype=object)
                            # Replace all the None with an appropriate fill value
                            mask = (data == None)  # noqa: E711
                            kind = np.dtype(col.subtype).kind
                            data[mask] = {'U': '', 'S': b''}.get(kind, 0)
                            arr_val = np.ma.array(data.astype(col.subtype), mask=mask)

                        col_vals.append(arr_val)

                    col.shape = ()
                    col.dtype = np.dtype(object)
                    # np.array(col_vals_arr, dtype=object) fails ?? so this workaround:
                    col.data = np.empty(len(col_vals), dtype=object)
                    col.data[:] = col_vals

                # Multidim columns with consistent shape (n, m, ...). These
                # might be masked.
                elif col.shape:
                    _check_dtype_is_str(col)

                    # Change empty (blank) values in original ECSV to something
                    # like "[[null, null],[null,null]]" so subsequent JSON
                    # decoding works. Delete `col.mask` so that later code in
                    # core TableOutputter.__call__() that deals with col.mask
                    # does not run (since handling is done here already).
                    if hasattr(col, 'mask'):
                        all_none_arr = np.full(shape=col.shape, fill_value=None, dtype=object)
                        all_none_json = json.dumps(all_none_arr.tolist())
                        for idx in np.nonzero(col.mask)[0]:
                            col.str_vals[idx] = all_none_json
                        del col.mask

                    col_vals = [json.loads(val) for val in col.str_vals]
                    # Make a numpy object array of col_vals to look for None
                    # (masked values)
                    data = np.array(col_vals, dtype=object)
                    mask = (data == None)  # noqa: E711
                    if not np.any(mask):
                        # No None's, just convert to required dtype
                        col.data = data.astype(col.subtype)
                    else:
                        # Replace all the None with an appropriate fill value
                        kind = np.dtype(col.subtype).kind
                        data[mask] = {'U': '', 'S': b''}.get(kind, 0)
                        # Finally make a MaskedArray with the filled data + mask
                        col.data = np.ma.array(data.astype(col.subtype), mask=mask)

                # Regular scalar value column
                else:
                    if col.subtype:
                        warnings.warn(f'unexpected subtype {col.subtype!r} set for column '
                                      f'{col.name!r}, using dtype={col.dtype!r} instead.',
                                      category=InvalidEcsvDatatypeWarning)
                    converter_func, _ = convert_numpy(col.dtype)
                    col.data = converter_func(col.str_vals)

                if col.data.shape[1:] != tuple(col.shape):
                    raise ValueError('shape mismatch between value and column specifier')

            except json.JSONDecodeError:
                raise ValueError(f'column {col.name!r} failed to convert: '
                                 'column value is not valid JSON')
            except Exception as exc:
                raise ValueError(f'column {col.name!r} failed to convert: {exc}')


class EcsvData(basic.BasicData):
    def _set_fill_values(self, cols):
        """READ: Set the fill values of the individual cols based on fill_values of BaseData

        For ECSV handle the corner case of data that has been serialized using
        the serialize_method='data_mask' option, which writes the full data and
        mask directly, AND where that table includes a string column with zero-length
        string entries ("") which are valid data.

        Normally the super() method will set col.fill_value=('', '0') to replace
        blanks with a '0'.  But for that corner case subset, instead do not do
        any filling.
        """
        super()._set_fill_values(cols)

        # Get the serialized columns spec.  It might not exist and there might
        # not even be any table meta, so punt in those cases.
        try:
            scs = self.header.table_meta['__serialized_columns__']
        except (AttributeError, KeyError):
            return

        # Got some serialized columns, so check for string type and serialized
        # as a MaskedColumn.  Without 'data_mask', MaskedColumn objects are
        # stored to ECSV as normal columns.
        for col in cols:
            if (col.dtype == 'str' and col.name in scs
                    and scs[col.name]['__class__'] == 'astropy.table.column.MaskedColumn'):
                col.fill_values = {}  # No data value replacement

    def str_vals(self):
        """WRITE: convert all values in table to a list of lists of strings

        This version considerably simplifies the base method:
        - No need to set fill values and column formats
        - No per-item formatting, just use repr()
        - Use JSON for object-type or multidim values
        - Only Column or MaskedColumn can end up as cols here.
        - Only replace masked values with "", not the generalized filling
        """
        for col in self.cols:
            if len(col.shape) > 1 or col.info.dtype.kind == 'O':
                def format_col_item(idx):
                    obj = col[idx]
                    try:
                        obj = obj.tolist()
                    except AttributeError:
                        pass
                    return json.dumps(obj, separators=(',', ':'))
            else:
                def format_col_item(idx):
                    return str(col[idx])

            try:
                col.str_vals = [format_col_item(idx) for idx in range(len(col))]
            except TypeError as exc:
                raise TypeError(f'could not convert column {col.info.name!r}'
                                f' to string: {exc}') from exc

            # Replace every masked value in a 1-d column with an empty string.
            # For multi-dim columns this gets done by JSON via "null".
            if hasattr(col, 'mask') and col.ndim == 1:
                for idx in col.mask.nonzero()[0]:
                    col.str_vals[idx] = ""

        out = [col.str_vals for col in self.cols]
        return out


class Ecsv(basic.Basic):
    """ECSV (Enhanced Character Separated Values) format table.

    Th ECSV format allows for specification of key table and column meta-data, in
    particular the data type and unit.

    See: https://github.com/astropy/astropy-APEs/blob/main/APE6.rst

    Examples
    --------

    >>> from astropy.table import Table
    >>> ecsv_content = '''# %ECSV 0.9
    ... # ---
    ... # datatype:
    ... # - {name: a, unit: m / s, datatype: int64, format: '%03d'}
    ... # - {name: b, unit: km, datatype: int64, description: This is column b}
    ... a b
    ... 001 2
    ... 004 3
    ... '''

    >>> Table.read(ecsv_content, format='ascii.ecsv')
    <Table length=2>
      a     b
    m / s   km
    int64 int64
    ----- -----
      001     2
      004     3

    """
    _format_name = 'ecsv'
    _description = 'Enhanced CSV'
    _io_registry_suffix = '.ecsv'

    header_class = EcsvHeader
    data_class = EcsvData
    outputter_class = EcsvOutputter

    max_ndim = None  # No limit on column dimensionality

    def update_table_data(self, table):
        """
        Update table columns in place if mixin columns are present.

        This is a hook to allow updating the table columns after name
        filtering but before setting up to write the data.  This is currently
        only used by ECSV and is otherwise just a pass-through.

        Parameters
        ----------
        table : `astropy.table.Table`
            Input table for writing

        Returns
        -------
        table : `astropy.table.Table`
            Output table for writing
        """
        with serialize_context_as('ecsv'):
            out = serialize.represent_mixins_as_columns(table)
        return out
