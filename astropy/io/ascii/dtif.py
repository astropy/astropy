# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Define the Data-table Text Interchange Format DTIF which allows for reading and
writing all the meta data associated with an astropy Table object.
"""

import re
import json

from ...utils import OrderedDict

from . import core


def _decode_list(data):
    rv = []
    for item in data:
        if isinstance(item, unicode):
            item = item.encode('utf-8')
        elif isinstance(item, list):
            item = _decode_list(item)
        elif isinstance(item, dict):
            item = _decode_dict(item)
        rv.append(item)
    return rv


def _decode_dict(data):
    rv = OrderedDict()
    for key, value in data.iteritems():
        if isinstance(key, unicode):
            key = key.encode('utf-8')
        if isinstance(value, unicode):
            value = value.encode('utf-8')
        elif isinstance(value, list):
            value = _decode_list(value)
        elif isinstance(value, dict):
            value = _decode_dict(value)
        rv[key] = value
    return rv


def _get_col_attributes(col):
    """
    Extract information from a column (apart from the values) that is required
    to fully serialize the column.
    """
    attrs = OrderedDict()
    attrs['name'] = col.name
    attrs['unit'] = str(col.unit)
    attrs['format'] = col.format
    attrs['description'] = col.description
    attrs['type'] = col.dtype.type.__name__
    attrs['meta'] = col.meta

    return attrs


class DtifHeader(core.BaseHeader):
    """Header class for which the column definition line starts with the
    comment character.  See the :class:`CommentedHeader` class  for an example.
    """
    def process_lines(self, lines):
        """Return only non-blank lines that start with the comment regexp.  For these
        lines strip out the matching characters and leading/trailing whitespace."""
        re_comment = re.compile(self.comment)
        for line in lines:
            match = re_comment.match(line)
            if match:
                out = line[match.end():].strip()
                if out:
                    yield out

    def write(self, lines):
        """
        Write header information in the DTIF ASCII format.  This format
        starts with a delimiter separated list of the column names in order
        to make this format readable by humans and simple csv-type readers.
        It then encodes the full table meta and column attributes and meta
        as JSON and pretty-prints this in the header.  Finally the delimited
        column names are repeated again, for humans and readers that look
        for the *last* comment line as defining the column names.
        """
        meta = OrderedDict()
        meta['version'] = 1.0
        meta['schema'] = 'astropy.table'
        meta['table_meta'] = self.table_meta
        meta['columns'] = [_get_col_attributes(col) for col in self.cols]

        outs = ['<DTIF encoding=ascii>']
        meta_json = json.dumps(meta, indent=2, separators=(',', ': '))
        outs.extend(meta_json.splitlines())

        lines.extend([self.write_comment + line for line in outs])
        lines.append(self.splitter.join([x.name for x in self.cols]))

    def get_cols(self, lines):
        """Initialize the header Column objects from the table ``lines``.

        :param lines: list of table lines
        :returns: None (but sets self.cols)
        """
        # Extract non-blank comment (header) lines with comment character striped
        lines = list(self.process_lines(lines))

        m = re.match(r'< \s* DTIF \s+ encoding=(\S+) \s* >', lines[0], re.VERBOSE)
        if not m:
            raise core.InconsistentTableError('DTIF file must start with "# <DTIF encoding=[encoding]>')
        encoding = m.group(1)  # Nothing done with encoding at the moment

        # Now actually load the JSON data structure into `meta`
        meta_json = '\n'.join(lines[1:])
        meta = json.loads(meta_json, object_pairs_hook=OrderedDict)
        meta = _decode_dict(meta)
        self.table_meta = meta['table_meta']

        # Create the list of io.ascii column objects from `meta`
        meta_cols = OrderedDict((x['name'], x) for x in meta['columns'])
        self.names = meta_cols.keys()
        self._set_cols_from_names()  # BaseHeader method to create self.cols

        # Transfer attributes from the column descriptor stored in the input
        # header JSON metadata to the new columns to create this table.
        for col in self.cols:
            for attr in ('description', 'format', 'unit', 'meta'):
                setattr(col, attr, meta_cols[col.name][attr])
            col.dtype = meta_cols[col.name]['type']


class Dtif(core.BaseReader):
    """Read a file where the column names are given in a line that begins with
    the header comment character. `header_start` can be used to specify the
    line index of column names, and it can be a negative index (for example -1
    for the last commented line).  The default delimiter is the <space>
    character.::

      # col1 col2 col3
      # Comment line
      1 2 3
      4 5 6
    """
    _format_name = 'dtif'
    _description = 'Data-table Text Interchange Format'

    def __init__(self):
        core.BaseReader.__init__(self)
        self.header = DtifHeader()
        self.header.data = self.data
        self.data.header = self.header
        self.header.splitter.delimiter = ' '
        self.data.splitter.delimiter = ' '
        self.header.start_line = 0
        self.data.start_line = 1
        self.header.comment = r'\s*#'
        self.header.write_comment = '# '
        self.data.comment = r'\s*#'
        self.data.write_comment = '# '
