# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""This package contains functions for reading and writing QDP tables that are
not meant to be used directly, but instead are available as readers/writers in
`astropy.table`. See :ref:`astropy:table_io` for more details.

:Author: Abdu Zoghbi
"""
import re
from collections import OrderedDict

from . import core


class TdatHeader(core.BaseHeader):
    """Header Reader for TDAT format"""

    start_line = 0
    comment = r"\s*(#|//)(.*)$"
    write_comment = "# "

    def update_meta(self, lines, meta):
        """Extract meta information: comments and key/values in the header"""
        comments = []
        keywords = OrderedDict()
        col_lines = []

        key_parser = re.compile(r'(?P<key>\w+)\s*=\s*([\'"`])?(?P<value>.*?)(?(2)\2|\s*$)')

        desc_comments = ['Table Parameters', 'Virtual Parameters',
                         'Relationship Definitions', 'Data Format Specification']
        
        is_header = False
        re_comment = re.compile(self.comment)
        for line in lines:
            if '<header>' in line.lower():
                is_header = True
                # we don't want <header> itself
                continue
            if '<data>' in line.lower():
                # end of header
                is_header = False
            
            if is_header:
                # get comments
                match = re_comment.match(line)
                kmatch = key_parser.match(line)
                if match:
                    add = True
                    for cc in desc_comments:
                        if cc in match.group(2):
                            add = False
                    if add and match.group(2) != '':
                        comments.append(match.group(2))
                elif kmatch:
                    keywords[kmatch.group('key')] = kmatch.group('value')
                else:
                    if 'field' in line:
                        col_lines.append(line)
        
        meta['table']['comments'] = comments
        meta['table']['keywords'] = keywords
        self.col_lines = col_lines
        

    def process_lines(self, lines):
        """Select lines between <HEADER> and <DATA>"""

        # if already processed, yield the lines
        if hasattr(self, 'col_lines'):
            for line in self.col_lines:
                yield line
        else:
            is_header = False
            re_comment = re.compile(self.comment)
            for line in lines:
                # skip comments
                if re_comment.match(line):
                    continue
                if '<header>' in line.lower():
                    is_header = True
                    # we don't want <header> itself
                    continue
                if '<data>' in line.lower():
                    # end of header
                    return
                if is_header:
                    yield line

    def get_cols(self, lines):
        """Identify the columns and table description"""
        # generator returning valid header lines
        header_lines = self.process_lines(lines)

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
        for line in header_lines:
            
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
    def __call__(self, lines):
        for line in lines:
            vals = line.strip().split(self.delimiter)[:-1]
            if self.process_val:
                yield [self.process_val(x) for x in vals]
            else:
                yield vals


class TdatData(core.BaseData):
    """Data Reader for TDAT format
    """

    start_line = 1
    comment = r"\s*(#|//).*$"
    write_comment = "# "
    splitter_class = TdatDataSplitter

    def process_lines(self, lines):
        """Select lines between <DATA> and <END>"""
        ll =  [line.strip() for line in lines if '|' in line]
        return ll
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
                return
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
