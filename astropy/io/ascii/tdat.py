# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""This package contains functions for reading and writing QDP tables that are
not meant to be used directly, but instead are available as readers/writers in
`astropy.table`. See :ref:`astropy:table_io` for more details.

:Author: Abdu Zoghbi
"""
import re
from collections import OrderedDict

from . import core

class TdatMeta:

    _comment = r'\s*(#|//)(.*)$'
    _keys = r'(?P<key>\w+)\s*=\s*([\'"`])?(?P<value>.*?)(?(2)\2|\s*$)'
    _desc_comments = ['Table Parameters', 'Virtual Parameters',
                      'Relationship Definitions', 'Data Format Specification']
    _comments = None
    _keywords = None
    _col_lines = None
    _delimiter = '|'
    
    
    def _process_lines(self, lines, ltype):
        """Process the lines to separate header|data"""

        is_header = False
        is_data = False

        for line in lines:
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

            if ltype == 'header' and is_header:
                yield line
            
            if ltype == 'data' and is_data:
                yield line


    def _parse_header(self, lines):
        """Parse header into keywords and comments and column descriptors"""

        # if already parsed, do not repeat
        if (self._comments is not None and 
            self._keywords is not None and
            self._col_lines is not None):
            return

        comments = []
        keywords = OrderedDict()
        col_lines = []

        keys_parser = re.compile(self._keys)
        comment_parser = re.compile(self._comment)
        desc_comments = self._desc_comments

        for line in self._process_lines(lines, 'header'):
            # get comments
            match = comment_parser.match(line)
            if match:
                add = True
                for cc in desc_comments:
                    if cc in match.group(2):
                        add = False
                        break
                if add and match.group(2) != '':
                    comments.append(match.group(2))
            
            kmatch = keys_parser.match(line)
            if kmatch:
                keywords[kmatch.group('key')] = kmatch.group('value')
            else:
                if 'field' in line:
                    col_lines.append(line)
        self._comments = comments
        self._keywords = keywords
        self._col_lines = col_lines

        self._delimiter = keywords.get('field_delimiter', '|')
            
    

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

                
class TdatDataSplitter(core.BaseSplitter):
    """Splitter for tdat data."""

    delimiter = '|'


class TdatData(core.BaseData, TdatMeta):
    """Data Reader for TDAT format
    """

    start_line = 1
    comment = r"\s*(#|//).*$"
    write_comment = "# "
    splitter_class = TdatDataSplitter

    def process_lines(self, lines):
        """Select lines between <DATA> and <END>"""
        is_data = False
        re_comment = re.compile(self.comment)
        data_lines = []
        for line in lines:
            # skip comments
            if re_comment.match(line):
                continue
            if '<data>' in line.lower():
                is_data = True
                # we don't want <data> itself
                continue
            if '<end>' in line.lower():
                # end of data
                break
            if is_data:
                data_lines.append(line)
        return data_lines


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
