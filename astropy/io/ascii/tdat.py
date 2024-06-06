# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""This package contains functions for reading and writing QDP tables that are
not meant to be used directly, but instead are available as readers/writers in
`astropy.table`. See :ref:`astropy:table_io` for more details.

:Author: Abdu Zoghbi
"""

import json
import re
from collections import OrderedDict
from warnings import warn

from astropy.utils.exceptions import AstropyWarning

from . import core

_STD_MSG = "See details in https://heasarc.gsfc.nasa.gov/docs/software/dbdocs/tdat.html"


class TdatFormatError(Exception):
    """Tdat Format Error"""

    def __str__(self):
        return f"{super().__str__()}\n{_STD_MSG}"


class TdatFormatWarning(AstropyWarning):
    """Tdat Format Warning"""


class TdatMeta:
    # comments: # or //
    _comment = r"\s*(#|//)(.*)$"

    # keywords in the header: name = value
    _keys = r"(?P<key>\w+)\s*=\s*([\'])?(?P<value>.*?)(?(2)\2|\s*$)"
    # keywords in the header: name[text] = some_other_text; name: relate|line
    _extra_keys = r"\s*(relate|line)\[(\w+)\]\s*=\s*([\w\s]+)(?:\((\w+)\))?"

    # descriptors that shouldn't be registered as comments
    _desc_comments = [
        "Table Parameters",
        "Virtual Parameters",
        "Relationship Definitions",
        "Data Format Specification",
    ]

    _comments = None
    _keywords = None
    _col_lines = None
    _line_fields = None
    _delimiter = "|"

    _deprecated_keywords = ("record_delimiter", "field_delimiter")
    _required_keywords = ("table_name",)

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
                if ltype != "comment":
                    continue
                add = True
                for cc in desc_comments:
                    if cc in cmatch.group(2):
                        # ignore keywords from desc_comments
                        add = False
                        break
                # exclude empty comment lines
                if cmatch.group(2) == "":
                    add = False
                if not add:
                    continue
                yield cmatch.group(2)

            if "<header>" in line.lower():
                is_header = True
                # we don't want <header> itself
                continue

            if "<data>" in line.lower():
                is_header = False
                is_data = True
                continue

            if "<end>" in line.lower():
                is_data = False
                continue

            if ltype == "header" and is_header and cmatch is None:
                yield line

            if ltype == "data" and is_data and cmatch is None:
                yield line

    def _parse_header(self, lines):
        """Parse header into keywords and comments and column descriptors"""
        # if already parsed, do not repeat
        if (
            self._comments is not None
            and self._keywords is not None
            and self._col_lines is not None
        ):
            return

        keywords = OrderedDict()
        col_lines = []
        line_fields = OrderedDict()

        keys_parser = re.compile(self._keys)
        extra_keys_parser = re.compile(self._extra_keys)

        for line in self._process_lines(lines, "header"):
            kmatch = keys_parser.match(line)
            ematch = extra_keys_parser.match(line)
            if kmatch:
                key = kmatch.group("key")
                if key in self._deprecated_keywords:
                    warn(
                        f'"{key}" keyword is obsolete and will be ignored. {_STD_MSG}',
                        TdatFormatWarning,
                    )
                else:
                    keywords[kmatch.group("key")] = kmatch.group("value")
            elif ematch:
                # match extra keywords
                if ematch.group(1) == "relate":
                    warn(
                        f'"relate" keyword is obsolete and will be ignored. {_STD_MSG}',
                        TdatFormatWarning,
                    )
                    continue
                if ematch.group(1) == "line":
                    # fields in line; typically only 1, but can be more
                    line_fields[ematch.group(2)] = ematch.group(3).split()
            else:
                if "field" in line:
                    col_lines.append(line)
        # check we have the required keywords

        self._comments = [c for c in self._process_lines(lines, "comment")]
        self._col_lines = col_lines

        # raise an error if no 'line[...] is found in the header'
        if len(line_fields) == 0:
            raise TdatFormatError(
                '"line[..]" keyword is required but not found in the header.\n'
            )
        self._line_fields = line_fields

        # raise and error if a required keyword is not present
        for key in self._required_keywords:
            if key not in keywords:
                raise TdatFormatError(
                    f'"{key}" keyword is required but not found in the header.\n'
                )
        self._keywords = keywords


class TdatHeader(core.BaseHeader, TdatMeta):
    """Header Reader for TDAT format"""

    start_line = 0
    write_comment = "# "

    def update_meta(self, lines, meta):
        """Extract meta information: comments and key/values in the header
        READ: Override the default update_meta
        """
        self._parse_header(lines)
        meta["table"]["comments"] = self._comments
        meta["table"]["keywords"] = self._keywords
        meta["table"]["field_lines"] = self._col_lines

    def process_lines(self, lines):
        """Select lines between <HEADER> and <DATA>"""
        # if already processed, yield the lines
        self._parse_header(lines)
        yield from self._col_lines

    def get_cols(self, lines):
        """Identify the columns and table description"""
        # generator returning valid header lines

        col_parser = re.compile(
            r"""
            \s*field
            \[(?P<name>\w+)\]\s*=\s*
            (?P<ctype>[^\W\s_]+)
            (?:\:(?P<fmt>[^\s_]+))?
            (?:_(?P<unit>[^\s_]+))?
            \s*
            (?:\[(?P<ucd>[\w\.\;]+)\])?
            \s*
            (?:\((?P<idx>\w+)\))?
            \s*
            (?:[//|#]+\s*(?P<desc>[^/#]*))?
            (?:[//|#]+\s*(?P<comment>[^/#]*))?
            \s*
            """,
            re.VERBOSE,
        )

        cols = []
        for line in self.process_lines(lines):
            # look for field[..]= ... column definitions
            cmatch = col_parser.match(line)
            if cmatch:
                col = core.Column(name=cmatch.group("name"))
                ctype = cmatch.group("ctype")
                if "int" in ctype:
                    col.dtype = int
                elif "char" in ctype:
                    col.dtype = str
                else:
                    col.dtype = float
                col.unit = cmatch.group("unit")
                col.format = cmatch.group("fmt")
                col.description = f'{cmatch.group("desc")}{cmatch.group("comment")}'
                col.comment = cmatch.group("comment")
                cols.append(col)

        self.names = [col.name for col in cols]
        self.cols = cols
        # check that cols and _line_fields are consistent or throw an error
        if sum([len(val) for val in self._line_fields.values()]) != len(cols):
            lsummary = [v for val in self._line_fields.values() for v in val]
            raise TdatFormatError(
                'The columns "field" descriptors are not consistent with '
                "the line[..] keyword.\n"
                f'"field" values: {self.names}\n'
                f"line[..] values: {lsummary}"
            )

    def write_comments(self, lines, meta):
        """
        WRITE: Override the default write_comments to include <HEADER> as first line
        """
        lines.append("<HEADER>")
        if self.write_comment not in (False, None):
            for comment in meta.get("comments", []):
                lines.append(self.write_comment + comment)

    def write(self, lines):
        """Write the keywords and column descriptors"""
        keywords = self.table_meta.get("keywords", None)
        col_lines = self.table_meta.get("field_lines", None)
        if keywords is not None:
            if "table_name" not in keywords:
                raise TdatFormatError(
                    '"table_name" keyword is required but not found in the header.\n'
                )
            lines.append(f'table_name = {keywords["table_name"]}')

            # loop through option table keywords
            for key in ["table_description", "table_document_url", "table_security"]:
                if key in keywords:
                    lines.append(f"{key} = {keywords[key]}")
        else:
            lines.append("table_name = astropy_table")

        # add table columns as fields
        lines.append("#")
        lines.append("# Table Parameters")
        lines.append("#")
        if col_lines is not None:
            for col in col_lines:
                lines.append(col)
        else:
            for col in self.cols:
                line = f"field[{col.name}] = {col.dtype}:{col.format}_{col.unit}"
                lines.append(line)
        lines.append("#")
        lines.append("# Data Format Specification")
        lines.append("#")
        lines.append(f"line[1] = {' '.join([col.name for col in self.cols])}")
        lines.append("#")

class TdatDataSplitter(core.BaseSplitter):
    """Splitter for tdat data."""

    delimiter = "|"

    def __call__(self, lines, field_delimiter="|", nlines=1):
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
                vals = [
                    val
                    for line in _lines
                    for val in re.split(rf"[{field_delimiter}]", line)[:-1]
                ]
                iline = 0
                _lines = []
            else:
                continue
            if self.process_val:
                yield [self.process_val(x) for x in vals]
            else:
                yield vals


class TdatData(core.BaseData, TdatMeta):
    """Data Reader for TDAT format"""

    write_comment = "# "
    splitter_class = TdatDataSplitter

    def process_lines(self, lines):
        """Select lines between <DATA> and <END>"""
        data_lines = list(self._process_lines(lines, "data"))
        return data_lines

    def get_str_vals(self):
        """Return a generator that returns a list of column values (as strings)
        for each data line.
        """
        field_delimiter = self._delimiter
        nlines = len(self.header._line_fields)
        # if we have  self._delimiter from the header, user it.
        if hasattr(self.header, "_delimiter"):
            field_delimiter = self.header._delimiter
        return self.splitter(self.data_lines, field_delimiter, nlines)

    def write(self, lines):
        """Write ``self.cols`` in place to ``lines``.

        Parameters
        ----------
        lines : list
            List for collecting output of writing self.cols.
        """
        lines.append("<DATA>")
        if callable(self.start_line):
            raise TypeError("Start_line attribute cannot be callable for write()")
        else:
            data_start_line = self.start_line or 0
        while len(lines) < data_start_line:
            lines.append(itertools.cycle(self.write_spacer_lines))

        col_str_iters = self.str_vals()
        for vals in zip(*col_str_iters):
            lines.append(self.splitter.join(vals)+"|")
        lines.append("<END>")


class Tdat(core.BaseReader):
    r"""Read TDAT format"""

    _format_name = "tdat"
    _description = "HEASARC tdat format"
    _io_registry_can_write = True
    _io_registry_suffix = ".tdat"

    header_class = TdatHeader
    data_class = TdatData

    def inconsistent_handler(self, str_vals, ncols):
        """Remove the last field separator if it exists"""
        if len(str_vals) == ncols + 1 and str_vals[-1] == "":
            str_vals = str_vals[:-1]
        return str_vals
