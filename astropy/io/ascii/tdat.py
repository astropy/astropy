# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""This package contains functions for reading and writing QDP tables that are
not meant to be used directly, but instead are available as readers/writers in
`astropy.table`. See :ref:`astropy:table_io` for more details.

:Author: Abdu Zoghbi
"""

import re
from collections import OrderedDict
from copy import deepcopy
from warnings import warn
import numpy as np
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
    _extra_keys = (
        r"\s*(relate|line|parameter_defaults)\[(\w+)\]\s*=\s*([\w\s]+)(?:\((\w+)\))?"
    )
    _deprecated_keys = r"\s*(record_delimiter|field_delimiter)\s*=\s\"(.+)\""
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
    _idx_lines = None
    _delimiter = "|"

    # _deprecated_keywords = ("record_delimiter", "field_delimiter")
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
        idx_lines = {"key": None, "index": []}

        keys_parser = re.compile(self._keys)
        extra_keys_parser = re.compile(self._extra_keys)
        deprecated_keys_parser = re.compile(self._deprecated_keys)

        for line in self._process_lines(lines, "header"):
            kmatch = keys_parser.match(line)
            ematch = extra_keys_parser.match(line)
            dmatch = deprecated_keys_parser.match(line)
            if dmatch:
                warn(
                    f'"{dmatch.group(1)}" keyword is deprecated. {_STD_MSG}',
                    TdatFormatWarning,
                )
                if dmatch.group(1) == "record_delimiter":
                    self._record_delimiter = dmatch.group(2)
                elif dmatch.group(1) == "field_delimiter":
                    self._delimiter = dmatch.group(2)
            elif kmatch:
                key = kmatch.group("key")
                # if key in self._deprecated_keywords:
                #     warn(
                #         f'"{key}" keyword is obsolete and will be ignored. {_STD_MSG}',
                #         TdatFormatWarning,
                #     )
                # else:
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
                    if "(key)" in line:
                        idx_lines["key"] = line.split("]")[0].split("[")[1]
                    elif "(index)" in line:
                        idx_lines["index"] += [line.split("]")[0].split("[")[1]]
        # check we have the required keywords

        self._comments = [c for c in self._process_lines(lines, "comment")]
        self._col_lines = col_lines
        self._idx_lines = idx_lines

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
    write_comment = "#"
    cols = None

    def update_meta(self, lines, meta):
        """Extract meta information: comments and key/values in the header
        READ: Override the default update_meta
        """
        self._parse_header(lines)
        meta["table"]["comments"] = self._comments
        meta["table"]["keywords"] = self._keywords
        meta["table"]["field_lines"] = self._col_lines
        meta["table"]["index_lines"] = self._idx_lines

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
                col.meta = OrderedDict()

                ctype = cmatch.group("ctype")
                if "int" in ctype:
                    col.dtype = int
                elif "char" in ctype:
                    col.dtype = str
                else:
                    col.dtype = float
                col.unit = cmatch.group("unit")
                col.format = cmatch.group("fmt")
                col.description = f'{cmatch.group("desc")}'
                col.meta['comment'] = cmatch.group("comment")
                col.meta['ucd'] = cmatch.group('ucd')
                col.meta['index'] = cmatch.group('idx')
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
        """Write the Table out to a TDAT formatted file.

        The header will be auto-populated by information in the Table,
        prioritizing information given in the Table.meta. The keys in the meta
        are as follows:
            comments : list or string, (optional)
                Table information which provide context. This information is
                included in the header preceding all other lines and commented
                out.
            keywords : OrderedDict, (optional, recommended)
                Header keywords which will appear in the file as "name=value" lines.
                Of particular importance are `table_name`, `table_description`,
                and `table_document_url`.
                `table_name` is a required keyword for the TDAT format and will
                be autopopulated with "astropy_table" if not specified in
                Table.meta["keywords"]
            field_lines : list, (optional)
                Specifications for each data column in the TDAT field format.
                See https://heasarc.gsfc.nasa.gov/docs/software/dbdocs/tdat.html
                for details. If this is not specified, field_lines will be
                inferred fromt the Column information. The `name`, `unit`, `format`
                and `description` are inferred directly from Column attributes,
                and the `index` tag, `ucd` and `comment` specifications can be
                included in Column.meta.
                For `ucd`, see https://cdsweb.u-strasbg.fr/UCD/ for acceptable entries.

        If there is no Table.meta, this writer will attempt to automatically
        generate the appropriate header information based on the table and
        column properties and the recommendations for the TDAT fromat by HEASARC.

        Parameters
        ----------
        lines : _type_
            _description_

        Raises
        ------
        ValueError
            _description_
        TdatFormatError
            _description_
        """
        if self.splitter.delimiter not in [" ", "|"]:
            raise ValueError("only pipe and space delimitter is allowed in tdat format")
        # Write the keywords and column descriptors
        keywords = deepcopy(self.table_meta.get("keywords", None))
        col_lines = self.table_meta.get("field_lines", None)

        if keywords is not None:
            if "table_name" not in keywords:
                raise TdatFormatError(
                    '"table_name" keyword is required but not found in the header.\n'
                )
            # Table names are limited to 20 characters
            lines.append(f'table_name = {keywords["table_name"][:20]}')

            # loop through option table keywords
            for key in ["table_description", "table_document_url", "table_security"]:
                if key == "table_description":
                    lines.append(f"{key} = {keywords.pop(key)[:80]}")
                elif key in keywords:
                    lines.append(f"{key} = {keywords.pop(key)}")
        else:
            lines.append("table_name = astropy_table")
            lines.append('table_description = "A table created via astropy"')

        # add table columns as fields
        lines.append("#")
        lines.append("# Table Parameters")
        lines.append("#")
        if col_lines is not None:
            for col in col_lines:
                lines.append(col)
        else:
            for col in self.cols:
                if col.dtype == int:
                    ctype = "integer"
                elif col.dtype == float:
                    ctype = "float"
                else:
                    ctype = f"char{str(col.dtype).split('<U')[-1]}"
                field_line = f"field[{col.name}] = {ctype}"
                if col.format is not None:
                    field_line += f":{col.format}"
                if col.unit is not None:
                    field_line += f"_{col.unit:vounit}"
                if "ucd" in col.meta:
                    field_line += f" [{col.meta['ucd']}]"
                if "index" in col.meta:
                    field_line += f" ({col.meta['index']})"
                if col.description is not None:
                    field_line += f" // {col.description}"
                if "comment" in col.meta:
                    field_line += f" // {col.meta['comment']}"
                lines.append(field_line)

        if keywords is not None and len(keywords) != 0:
            if "parameter_defaults" in keywords:
                lines.append("#")
                lines.append(
                    f"{'parameter_defaults'} = {keywords.pop('parameter_defaults')}"
                )
        if keywords is not None and len(keywords) != 0:
            lines.append("#")
            lines.append("# Virtual Parameters")
            lines.append("#")
            for key in keywords:
                lines.append(f"{key} = {keywords[key]}")
        lines.append("#")
        lines.append("# Data Format Specification")
        lines.append("#")
        lines.append(f"line[1] = {' '.join([col.name for col in self.cols])}")
        lines.append("#")


class TdatDataSplitter(core.BaseSplitter):
    """Splitter for tdat data."""

    delimiter = "|"

    def preprocess_data_lines(self, lines, record_delimiter):
        """Handle the case of a record delimiter.

        The record_delimiter (deprecated) can be specified in the header. By
        default there is no record delimiter and new records should be set on
        new lines.
        The following list standard escaped character sequences and their
        equivalent meanings can be used: \t (tab), \b (backspace), \r (carriage
        return), \f (form feed), \v (vertical tab), \a (audible alert/bell),
        and \### (where ### is a number between 1 and 127 and represents the
        ASCII character with that numerical code). Note: Specifying a record
        delimiter value of "" is interpreted as a single blank line between
        records.

        Parameters
        ----------
        lines : list
        record_delimiter : str

        Returns
        -------
        data_lines : list
        """
        data_lines = []
        for line in lines:
            data_lines += re.split(rf"[{record_delimiter}]", line)
        return data_lines

    def __call__(self, lines, field_delimiter="|", record_delimiter=None, nlines=1):
        """Handle the case of multiple delimiter.

        The delimiter is specified in the header, and can have multiple values
        e.g. field_delimiter = "|!"

        nlines gives the number of lines to read at once to give one data row.

        """
        if record_delimiter is not None:
            lines = self.preprocess_data_lines(lines, record_delimiter)
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
        record_delimiter = None
        nlines = len(self.header._line_fields)
        # if we have  self._delimiter from the header, user it.
        if hasattr(self.header, "_delimiter"):
            field_delimiter = self.header._delimiter
        if hasattr(self.header, "_record_delimiter"):
            record_delimiter = self.header._record_delimiter
        return self.splitter(self.data_lines, field_delimiter, record_delimiter, nlines)

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
            # TDAT specification requires a delimiter at the end of each record
            lines.append(self.splitter.join(vals) + "|")
        lines.append("<END>")


class TdatOutputter(core.TableOutputter):
    """
    Output the table as an astropy.table.Table object.
    """

    def __call__(self, cols, meta):
        """
        READ: Override the default outputter.
        TDAT files may (optionally) specify which field lines should be used as
        the primary index and secondary indices Astropy tables support adding
        indices after creation. This overwrite adds labeled indices on read.
        """
        # Sets col.data to numpy array and col.type to io.ascii Type class (e.g.
        # FloatType) for each col.
        self._convert_vals(cols)

        t_cols = [
            np.ma.MaskedArray(x.data, mask=x.mask)
            if hasattr(x, "mask") and np.any(x.mask)
            else x.data
            for x in cols
        ]
        out = core.Table(t_cols, names=[x.name for x in cols], meta=meta["table"])

        for col, out_col in zip(cols, out.columns.values()):
            for attr in ("format", "unit", "description"):
                if hasattr(col, attr):
                    setattr(out_col, attr, getattr(col, attr))
            if hasattr(col, "meta"):
                out_col.meta.update(col.meta)
        # Add indices, if specified
        if meta["table"]["index_lines"]["key"] is not None:
            out.add_index(meta["table"]["index_lines"]["key"])
        if len(meta["table"]["index_lines"]["index"]) > 0:
            for idx in meta["table"]["index_lines"]["index"]:
                out.add_index(idx)
        return out


class Tdat(core.BaseReader):
    r"""Read TDAT format"""

    _format_name = "tdat"
    _description = "HEASARC tdat format"
    _io_registry_can_write = True
    _io_registry_suffix = ".tdat"

    header_class = TdatHeader
    data_class = TdatData
    outputter_class = TdatOutputter

    def inconsistent_handler(self, str_vals, ncols):
        """Remove the last field separator if it exists"""
        if len(str_vals) == ncols + 1 and str_vals[-1] == "":
            str_vals = str_vals[:-1]
        return str_vals
