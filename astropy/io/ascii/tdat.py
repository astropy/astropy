# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""This package contains functions for reading and writing QDP tables that are
not meant to be used directly, but instead are available as readers/writers in
`astropy.table`. See :ref:`astropy:table_io` for more details.

:Author: Abdu Zoghbi
"""
import re

from . import core


class TdatHeader(core.BaseHeader):
    """Header Reader for TDAT format"""

    start_line = 0
    comment = r"\s*(#|//).*$"
    write_comment = "# "

    def process_lines(self, lines):
        """Select lines between <HEADER> and <DATA>"""
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
            match = col_parser.match(line)
            if match:
                col = core.Column(name=match.group('name'))
                col.ctype = match.group('ctype')
                col.fmt = match.group('fmt')
                col.unit = match.group('unit')
                col.ucd = match.group('ucd')
                col.idx = match.group('idx')
                col.desc = match.group('desc')
                col.comment = match.group('comment')
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
