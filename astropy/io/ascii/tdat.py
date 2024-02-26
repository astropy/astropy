# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""This package contains functions for reading and writing QDP tables that are
not meant to be used directly, but instead are available as readers/writers in
`astropy.table`. See :ref:`astropy:table_io` for more details.

:Author: Abdu Zoghbi
"""
import re
from collections import OrderedDict
from warnings import warn
from astropy.utils.exceptions import AstropyWarning

from . import core

_STD_MSG = 'See details in https://heasarc.gsfc.nasa.gov/docs/software/dbdocs/tdat.html'

class TdatFormatError(Exception):
    """Tdat Format Error"""
    def __str__(self):
        return "{}\n{}".format(super().__str__(),_STD_MSG)

class TdatFormatWarning(AstropyWarning):
    """Tdat Format Warning"""


class TdatMeta:

    # comments: # or //
    _comment = r'\s*(#|//)(.*)$'

    # keywords in the header: name = value
    _keys = r'(?P<key>\w+)\s*=\s*([\'"`])?(?P<value>.*?)(?(2)\2|\s*$)'
    # keywords in the header: name[text] = some_other_text; name: relate|line
    _extra_keys = r'\s*(relate|line)\[(\w+)\]\s*=\s*([\w\s]+)(?:\((\w+)\))?'

    # descriptors that shouldn't be registered as comments
    _desc_comments = ['Table Parameters', 'Virtual Parameters',
                      'Relationship Definitions', 'Data Format Specification']

    _comments = None
    _keywords = None
    _col_lines = None
    _line_fields = None
    _delimiter = '|'

    _deprecated_keywords = ('record_delimiter', 'field_delimiter')
    _required_keywords = ('table_name')
    
    
    def _process_lines(self, lines, ltype):
        """Process the lines to separate header|data|comment"""

        is_header = False
        is_data = False
        comment_parser = re.compile(self._comment)

        desc_comments = self._desc_comments

        for line in lines:

            # handle comments first
            cmatch = comment_parser.match(line)

            if cmatch:
                if ltype != 'comment':
                    continue
                add = True
                for cc in desc_comments:
                    if cc in cmatch.group(2):
                        # ignore keywords from desc_comments
                        add = False
                        break
                # exclude empty comment lines
                if cmatch.group(2) == '':
                    add = False
                if not add:
                    continue
                yield cmatch.group(2)

            if '<header>' in line.lower():
                is_header = True
                # we don't want <header> itself
                continue
            
            if '<data>' in line.lower():
                is_header = False
                is_data = True
                continue

            if '<end>' in line.lower():
                is_data = False
                continue

            if ltype == 'header' and is_header and cmatch is None:
                yield line
            
            if ltype == 'data' and is_data and cmatch is None:
                yield line


    def _parse_header(self, lines):
        """Parse header into keywords and comments and column descriptors"""

        # if already parsed, do not repeat
        if (self._comments is not None and 
            self._keywords is not None and
            self._col_lines is not None):
            return

        keywords = OrderedDict()
        col_lines = []
        line_fields = OrderedDict()

        keys_parser = re.compile(self._keys)
        extra_keys_parser = re.compile(self._extra_keys)

        for line in self._process_lines(lines, 'header'):
            kmatch = keys_parser.match(line)
            ematch = extra_keys_parser.match(line)
            if kmatch:
                key = kmatch.group('key')
                if key in self._deprecated_keywords:
                    warn(
                        f'"{key}" keyword is deprecated from the tdat standard. '
                        f'It will be ignored. {_STD_MSG}',
                        TdatFormatWarning
                    )
                else:
                    keywords[kmatch.group('key')] = kmatch.group('value')
            elif ematch:
                # match extra keywords
                if ematch.group(1) == 'relate':
                    warn(
                        '"relate" keyword is deprecated from the tdat standard. '
                        f'It will be ignored. {_STD_MSG}',
                        TdatFormatWarning
                    )
                    continue
                if ematch.group(1) == 'line':
                    # fields in line; typically only 1, but can be more
                    line_fields[ematch.group(2)] = ematch.group(3).split()
            else:
                if 'field' in line:
                    col_lines.append(line)
        # check we have the required keywords

        self._comments = [c for c in self._process_lines(lines, 'comment')]
        self._keywords = keywords
        self._col_lines = col_lines

        # raise an error if no 'line[...] is found in the header'
        if len(line_fields) == 0:
            raise TdatFormatError(
                '"line[..]" keyword is required but not found in the header.\n'
            )
        self._line_fields = line_fields


class TdatHeader(core.BaseHeader, TdatMeta):
    """Header Reader for TDAT format"""

    start_line = 0
    write_comment = "# "

    def update_meta(self, lines, meta):
        """Extract meta information: comments and key/values in the header"""
        
        self._parse_header(lines)
        meta['table']['comments'] = self._comments
        meta['table']['keywords'] = self._keywords
        

    def process_lines(self, lines):
        """Select lines between <HEADER> and <DATA>"""

        # if already processed, yield the lines
        self._parse_header(lines)
        for line in self._col_lines:
            yield line
        

    def get_cols(self, lines):
        """Identify the columns and table description"""
        # generator returning valid header lines

        col_parser = re.compile(r"""
            \s*field
            \[(?P<name>\w+)\]\s*=\s*
            (?P<ctype>[^\W\s_]+)
            (?:\:(?P<fmt>[^\s_]+))?
            (?:_(?P<unit>[^\s_]+))?
            \s*
            (?:\[(?P<ucd>[\w\.]+)\])?
            \s+
            (?:\((?P<idx>\w+)\))?
            \s+
            (?:[//|#]+\s*(?P<desc>[^/#]*))?
            (?:[//|#]+\s*(?P<comment>[^/#]*))?
            \s*
            """, re.VERBOSE)
        
        cols = []
        for line in self.process_lines(lines):
            
            # look for field[..]= ... column definitions
            cmatch = col_parser.match(line)
            if cmatch:
                col = core.Column(name=cmatch.group('name'))
                col.ctype = cmatch.group('ctype')
                col.fmt = cmatch.group('fmt')
                col.unit = cmatch.group('unit')
                col.ucd = cmatch.group('ucd')
                col.idx = cmatch.group('idx')
                col.description = cmatch.group('desc')
                col.comment = cmatch.group('comment')
                cols.append(col)

        self.names = [col.name for col in cols]
        self.cols = cols

        # check that cols and _line_fields are consistent or throw an error
        if sum([len(val) for val in self._line_fields.values()]) != len(cols):
            lsummary = [v for val in self._line_fields.values() for v in val]
            raise TdatFormatError(
                'The columns "field" descriptors are not consistent with '
                'the line[..] keyword.\n'
                f'"field" values: {self.names}\n'
                f'line[..] values: {lsummary}'
            )

                
class TdatDataSplitter(core.BaseSplitter):
    """Splitter for tdat data."""

    delimiter = '|'

    def __call__(self, lines, field_delimiter='|', nlines=1):
        """Handle the case of multiple delimiter.

        The delimiter is specified in the header, and can have multiple values
        e.g. field_delimiter = "|!"

        nlines gives the number of lines to read at once to give one data row.

        """
        if self.process_line:
            lines = (self.process_line(x) for x in lines)
        iline = 0
        _lines = []
        for line in lines:
            _lines.append(line)
            iline += 1
            if iline == nlines:
                vals = [val for line in _lines
                        for val in re.split(fr'[{field_delimiter}]', line)[:-1]]
                iline = 0
                _lines = []
            else:
                continue
            if self.process_val:
                yield [self.process_val(x) for x in vals]
            else:
                yield vals


class TdatData(core.BaseData, TdatMeta):
    """Data Reader for TDAT format
    """

    write_comment = "# "
    splitter_class = TdatDataSplitter

    def process_lines(self, lines):
        """Select lines between <DATA> and <END>"""

        data_lines = [line for line in self._process_lines(lines, 'data')]
        return data_lines


    def get_str_vals(self):
        """Return a generator that returns a list of column values (as strings)
        for each data line.
        """
        field_delimiter = self._delimiter
        nlines = len(self.header._line_fields)
        # if we have  self._delimiter from the header, user it.
        if hasattr(self.header, '_delimiter'):
            field_delimiter = self.header._delimiter
        return self.splitter(self.data_lines, field_delimiter, nlines)


class Tdat(core.BaseReader):
    r"""Read TDAT format
    """

    _format_name = "tdat"
    _description = "HEASARC tdat format"
    _io_registry_can_write = True
    _io_registry_suffix = ".tdat"

    header_class = TdatHeader
    data_class = TdatData

    def inconsistent_handler(self, str_vals, ncols):
        """Remove the last field separator if it exists"""
        if len(str_vals) == ncols + 1 and str_vals[-1] == '':
            str_vals = str_vals[:-1]
        return str_vals
