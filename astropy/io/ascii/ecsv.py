# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Define the Enhanced Character-Separated-Values (ECSV) which allows for reading and
writing all the meta data associated with an astropy Table object.
"""

import re

from ...utils import OrderedDict
from ...extern import six

from . import core, basic
from ...table import meta

ECSV_VERSION = '0.9'
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
            raise ValueError('only space and comma are allowed for delimiter in ECVS format')

        for col in self.cols:
            if len(getattr(col, 'shape', ())) > 1:
                raise ValueError("ECSV format does not support multidimensional column '{0}'"
                                 .format(col.info.name))

        # Now assemble the header dict that will be serialized by the YAML dumper
        header = {'cols': self.cols}

        if self.table_meta:
            header['meta'] = self.table_meta

        # Set the delimiter only for the non-default option(s)
        if self.splitter.delimiter != ' ':
            header['delimiter'] = self.splitter.delimiter

        header_yaml_lines = (['%ECSV {0}'.format(ECSV_VERSION),
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
        except meta.YamlParseError:
            raise core.InconsistentTableError('unable to parse yaml in meta header')

        if 'meta' in header:
            self.table_meta = header['meta']

        if 'delimiter' in header:
            delimiter = header['delimiter']
            if delimiter not in DELIMITERS:
                raise ValueError('only space and comma are allowed for delimiter in ECVS format')
            self.splitter.delimiter = delimiter
            self.data.splitter.delimiter = delimiter

        # Create the list of io.ascii column objects from `header`
        header_cols = OrderedDict((x['name'], x) for x in header['datatype'])
        self.names = [x['name'] for x in header['datatype']]
        self._set_cols_from_names()  # BaseHeader method to create self.cols

        # Transfer attributes from the column descriptor stored in the input
        # header YAML metadata to the new columns to create this table.
        for col in self.cols:
            for attr in ('description', 'format', 'unit', 'meta'):
                if attr in header_cols[col.name]:
                    setattr(col, attr, header_cols[col.name][attr])
            col.dtype = header_cols[col.name]['datatype']
            # ECSV "string" means numpy dtype.kind == 'U' AKA str in Python 3
            if six.PY3 and col.dtype == 'string':
                col.dtype = 'str'
            if col.dtype.startswith('complex'):
                raise TypeError('ecsv reader does not support complex number types')


class EcsvOutputter(core.TableOutputter):
    """
    Output the table as an astropy.table.Table object.  This overrides the
    default converters to be an empty list because there is no "guessing"
    of the conversion function.
    """
    default_converters = []


class Ecsv(basic.Basic):
    """
    Read a file which conforms to the ECSV (Enhanced Character Separated
    Values) format.  This format allows for specification of key table
    and column meta-data, in particular the data type and unit.  For details
    see: https://github.com/astropy/astropy-APEs/blob/master/APE6.rst.

    For example::

      # %ECSV 0.9
      # ---
      # columns:
      # - {name: a, unit: m / s, type: int64, format: '%03d'}
      # - {name: b, unit: km, type: int64, description: This is column b}
      a b
      001 2
      004 3
    """
    _format_name = 'ecsv'
    _description = 'Enhanced CSV'

    header_class = EcsvHeader
    outputter_class = EcsvOutputter
