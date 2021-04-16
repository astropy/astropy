# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Define the Enhanced Character-Separated-Values (ECSV) which allows for reading and
writing all the meta data associated with an astropy Table object.
"""

import re
from collections import OrderedDict
import warnings

import numpy as np

from . import core, basic
from astropy.table import meta, serialize
from astropy.utils.data_info import serialize_context_as
from astropy.utils.exceptions import AstropyWarning
from astropy.io.ascii.core import convert_numpy

__doctest_requires__ = {'Ecsv': ['yaml']}

ECSV_VERSION = '1.0'
DELIMITERS = (' ', ',')


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
        Write header information in the ECSV ASCII format.  This format
        starts with a delimiter separated list of the column names in order
        to make this format readable by humans and simple csv-type readers.
        It then encodes the full table meta and column attributes and meta
        as YAML and pretty-prints this in the header.  Finally the delimited
        column names are repeated again, for humans and readers that look
        for the *last* comment line as defining the column names.
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
        Override the default write_comments to do nothing since this is handled
        in the custom write method.
        """
        pass

    def update_meta(self, lines, meta):
        """
        Override the default update_meta to do nothing.  This process is done
        in get_cols() for this reader.
        """
        pass

    def get_cols(self, lines):
        """
        Initialize the header Column objects from the table ``lines``.

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
        # ecsv_version could be constructed here, but it is not currently used.

        try:
            header = meta.get_header_from_yaml(lines)
        except ImportError as exc:
            if 'PyYAML package is required' in str(exc):
                warnings.warn("file looks like ECSV format but PyYAML is not installed "
                              "so it cannot be parsed as ECSV",
                              AstropyWarning)
            raise core.InconsistentTableError('unable to parse yaml in meta header'
                                              ' (PyYAML package is required)')
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
            for attr in ('description', 'format', 'unit', 'meta', 'shape'):
                if attr in header_cols[col.name]:
                    setattr(col, attr, header_cols[col.name][attr])
            col.dtype = header_cols[col.name]['datatype']
            # ECSV "string" means numpy dtype.kind == 'U' AKA str in Python 3
            if col.dtype == 'string':
                col.dtype = 'str'
            if col.dtype.startswith('complex'):
                raise TypeError('ecsv reader does not support complex number types')


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

        This is adapted (and shortened) from BaseOutputter._convert_vals. In
        the case of ECSV there is no guessing and all types are known in
        advance.
        """
        for col in cols:
            converter_func, _ = convert_numpy(col.dtype)
            try:
                # Object columns are serialized as JSON.
                if np.dtype(col.dtype).kind == 'O':
                    import json
                    col_vals = [json.loads(val) for val in col.str_vals]
                    col.data = np.array(col_vals, dtype=col.dtype)

                # Multidim columns are serialized as JSON
                elif col.shape:
                    import json
                    col_vals = [json.loads(val) for val in col.str_vals]
                    # Make a numpy object array of col_vals to look for None
                    data = np.array(col_vals, dtype=object)
                    mask = (data == None)  # noqa
                    if not np.any(mask):
                        # No None's, just convert to required dtype
                        col.data = data.astype(col.dtype)
                    else:
                        # Replace all the None with an appropriate fill value
                        kind = np.dtype(col.dtype).kind
                        if kind == 'U':
                            fill = ''
                        elif kind == 'S':
                            fill = b''
                        else:
                            fill = 0
                        data[mask] = fill
                        # Finally make a MaskedArray with the filled data + mask
                        col.data = np.ma.array(data.astype(col.dtype), mask=mask)
                else:
                    col.data = converter_func(col.str_vals)

                if col.data.shape[1:] != tuple(col.shape):
                    raise ValueError('shape mismatch')
            except Exception as exc:
                raise ValueError(f'Column {col.name} failed to convert: {exc}')


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
            if len(col.shape) > 1:
                import json

                def format_col_item(idx):
                    return json.dumps(col[idx].tolist(), separators=(',', ':'))
            elif col.info.dtype.kind == 'O':
                import json

                def format_col_item(idx):
                    return json.dumps(col[idx], separators=(',', ':'))
            else:
                def format_col_item(idx):
                    return str(col[idx])

            col.str_vals = [format_col_item(idx) for idx in range(len(col))]

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

    def _check_multidim_table(self, table, max_ndim=1):
        """Check that ``table`` has only 1-d columns.

        This is overriding a base method that raises an exception if non-allowed
        columns (dim > max_ndim) are present. The ECSV reader class supports N-d
        columns so just pass.

        Parameters
        ----------
        table : `~astropy.table.Table` Input table. max_ndim : int Max allowed
            number of dimensions (default=1)
        max_ndim : int
            Max allowed number of dimensions (default=1)
        """
        pass
