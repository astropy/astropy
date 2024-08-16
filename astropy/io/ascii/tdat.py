# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""This package contains functions for reading and writing TDAT tables that are
not meant to be used directly, but instead are available as readers/writers in
`astropy.table`. See :ref:`astropy:table_io` for more details.

:Author: Daniel Giles (daniel.k.giles@gmail.com)
         Abdu Zoghbi

"""

import re
from collections import OrderedDict
from copy import deepcopy
from warnings import warn

import numpy as np

from astropy.utils.exceptions import AstropyWarning

from . import basic, core

_STD_MSG = "See details in https://heasarc.gsfc.nasa.gov/docs/software/dbdocs/tdat.html"


def make_example_data():
    example_lines = [
        "<HEADER>",
        "# # and // are comments",
        "table_name = example_table",
        'table_description = "Example table"',
        "#",
        "# Table Parameters",
        "#",
        "field[id] = integer [meta.id] (key) // Unique ID",
        "field[ra] = float:.4f_degree [pos.eq.ra] (index) // Right Ascension",
        "field[name] = char12 [meta.id] // Name",
        "#",
        "# Virtual Parameters",
        "#",
        "table_author = Example et al.",
        "#",
        "# Data Format Specification",
        "#",
        "line[1] = id name ra",
        "<DATA>",
        "1|TargetOne|1.0|",
        "2|TargetTwo|2.0|",
        "<END>",
    ]
    return example_lines


class TdatFormatError(Exception):
    """Tdat Format Error"""

    def __str__(self):
        return f"{super().__str__()}\n{_STD_MSG}"


class TdatFormatWarning(AstropyWarning):
    """Tdat Format Warning"""


class TdatHeader(basic.BasicHeader):
    """TDAT table header"""

    cols = None
    comment = r"\s*(#|//)"
    # comments: # or //
    _comment = r"\s*(#|//)(.*)$"
    # descriptors that shouldn't be registered as comments
    _desc_comments = [
        "Table Parameters",
        "Virtual Parameters",
        "Relationship Definitions",
        "Data Format Specification",
    ]
    # keywords in the header: name = value
    _keys = r"(?P<key>\w+)\s*=\s*([\'])?(?P<value>.*?)(?(2)\2|\s*$)"
    # keywords in the header: name[text] = some_other_text;
    # names: relate|line
    _extra_keys = r"\s*(relate|line)\[(\w+)\]\s*=\s*([\w\s]+)(?:\((\w+)\))?"
    _field_line = r"\s*field\[\w+\]\s*=\.*"

    _deprecated_keys = r"\s*(record_delimiter|field_delimiter)\s*=\s\"(.+)\""
    _required_keywords = ("table_name",)
    _dtype_dict = {
        "int": "int4",
        "int32": "int4",
        "int16": "int2",
        "int8": "int1",
        "float": "float8",
        "float64": "float8",
        "float32": "float4",
        "float8": float,
        "float4": np.float32,
        "real": np.float32,
        "int4": int,
        "int2": np.int16,
        "integer2": np.int16,
        "smallint": np.int16,
        "int1": np.int8,
        "integer1": np.int8,
        "tinyint": np.int8,
    }

    def _validate_comment(self, line) -> bool:
        """Check if line is a valid comment, return comment if so"""
        comment_parser = re.compile(self._comment)
        cmatch = comment_parser.match(line)
        if cmatch:
            # Ignore common section headers
            for cc in self._desc_comments:
                if cc in cmatch.group(2):
                    return False
            # Ignore empty comments
            return cmatch.group(2) != ""
        return False

    def _process_comment(self, line: str) -> str:
        line = re.sub("^" + self.comment, "", line).strip()
        return line

    def _process_keywords(self, lines):
        line_fields = OrderedDict()
        keywords = getattr(self, "_keywords", OrderedDict())
        keys_parser = re.compile(self._keys)
        extra_keys_parser = re.compile(self._extra_keys)
        deprecated_keys_parser = re.compile(self._deprecated_keys)

        for line in lines:
            kmatch = keys_parser.match(line)
            ematch = extra_keys_parser.match(line)
            dmatch = deprecated_keys_parser.match(line)
            if dmatch:
                warn(
                    f'"{dmatch.group(1)}" keyword is deprecated. {_STD_MSG}',
                    TdatFormatWarning,
                )
                if dmatch.group(1) == "record_delimiter":
                    self.record_delimiter = dmatch.group(2)
                elif dmatch.group(1) == "field_delimiter":
                    self.delimiter = dmatch.group(2)
            elif ematch:
                # match extra keywords
                if ematch.group(1) == "relate":
                    warn(
                        f'"relate" keyword is obsolete and will be ignored. {_STD_MSG}',
                        TdatFormatWarning,
                    )
                elif ematch.group(1) == "line":
                    # fields in line; typically only 1, but can be more
                    line_fields[ematch.group(2)] = ematch.group(3).split()

            elif kmatch:
                key = kmatch.group("key")
                keywords[key] = kmatch.group("value")

        self._line_fields = line_fields

        if len(keywords) > len(getattr(self, "_keywords", OrderedDict())):
            self._keywords = keywords

    def update_meta(self, lines, meta):
        """Extract meta information: comments and key/values in the header
        READ: Overrides the default update_meta
        """
        try:
            start_line = (
                min(i for i, line in enumerate(lines) if line.strip() == "<HEADER>") + 1
            )
        except ValueError:
            raise TdatFormatError("<HEADER> not found in file." + _STD_MSG)
        try:
            end_line = min(
                i for i, line in enumerate(lines) if line.strip() == "<DATA>"
            )
        except ValueError:
            raise TdatFormatError("<DATA> not found in file." + _STD_MSG)

        meta["table"]["comments"] = [
            self._process_comment(line)
            for line in lines[start_line:end_line]
            if self._validate_comment(line)
        ]
        self._process_keywords(lines)
        meta["table"]["keywords"] = self._keywords

    def process_lines(self, lines):
        """Select and process lines between <HEADER> and <DATA> lines
        READ: Override default process_lines
        """
        fl_parser = re.compile(self._field_line)
        try:
            start_line = (
                min(i for i, line in enumerate(lines) if line.strip() == "<HEADER>") + 1
            )
        except ValueError:
            raise TdatFormatError("<HEADER> not found in file." + _STD_MSG)
        try:
            end_line = min(
                i for i, line in enumerate(lines) if line.strip() == "<DATA>"
            )
        except ValueError:
            raise TdatFormatError("<DATA> not found in file." + _STD_MSG)
        for line in lines[start_line:end_line]:
            if fl_parser.match(line):
                yield line

    def get_cols(self, lines):
        """Initialize the header Column objects from TDAT field lines
        READ: Overrides default get_cols
        """
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
            (?:\((?P<index>\w+)\))?
            \s*
            (?:[//|#]+\s*(?P<desc>[^/#]*))?
            (?:[//|#]+\s*(?P<comment>[^/#]*))?
            \s*
            """,
            re.VERBOSE,
        )

        cols = {}
        keywords = getattr(self, "_keywords", OrderedDict())
        for line in self.process_lines(lines):
            # look for field[..]= ... column definitions
            cmatch = col_parser.match(line)
            if cmatch:
                col = core.Column(name=cmatch.group("name"))
                col.meta = OrderedDict()

                ctype = cmatch.group("ctype")
                if ctype in self._dtype_dict:
                    col.dtype = self._dtype_dict[ctype]
                elif "int" in ctype:
                    col.dtype = int
                elif "float" in ctype:
                    col.dtype = float
                elif "char" in ctype:
                    col.dtype = str
                else:
                    raise TdatFormatError(
                        f"Unrecognized data type {ctype} for {col.name}."
                    )
                col.unit = cmatch.group("unit")
                col.format = cmatch.group("fmt")
                col.description = f'{cmatch.group("desc")}'.strip()
                for val in ["comment", "ucd", "index"]:
                    if cmatch.group(val) is not None:
                        col.meta[val] = cmatch.group(val).strip()
                cols[col.name] = col
        if len(keywords) > len(getattr(self, "_keywords", OrderedDict())):
            self._keywords = keywords

        self.names = [
            val for line_field in self._line_fields.values() for val in line_field
        ]

        # check that cols and _line_fields are consistent or throw an error
        if len(self.names) != len(cols):
            colnames = [col.name for col in cols]
            raise TdatFormatError(
                'The columns "field" descriptors are not consistent with '
                "the line[..] keyword.\n"
                f'"field" values: {colnames}\n'
                f"line[..] values: {self.names}"
            )
        self.cols = []
        for name in self.names:
            self.cols.append(cols[name])

    def write_comments(self, lines, meta):
        """Write comment lines in the header
        WRITE: Override the default write_comments to include <HEADER> as first line
        """
        lines.append("<HEADER>")
        if self.write_comment not in (False, None):
            for comment in meta.get("comments", []):
                lines.append(self.write_comment + comment.lstrip())

    def write(self, lines):
        """Write the Table out to a TDAT formatted file."""
        if self.splitter.delimiter not in [" ", "|"]:
            raise ValueError("only pipe and space delimiter is allowed in tdat format")

        # Write the keywords and column descriptors
        keywords = deepcopy(self.table_meta.get("keywords", {}))
        meta_keys = ["keywords", "comments"]
        # In case a user puts a keyword directly in meta, instead of meta.keywords
        for key in [
            key.lower()
            for key in self.table_meta.keys()
            if key.lower() not in meta_keys
        ]:
            keywords[key] = self.table_meta.get(key)

        indices = [col.name for col in self.cols if col.indices != []]

        if "table_name" in keywords:
            if len(keywords["table_name"]) > 20:
                warn(
                    "'table_name' is too long, truncating to 20 characters",
                    TdatFormatWarning,
                )
            lines.append(f'table_name = {keywords.pop("table_name")[:20]}')
        else:
            warn(
                "'table_name' must be specified\n"  # noqa: ISC003
                + f"{_STD_MSG}\n"
                + "This should be specified in the Table.meta.\n"
                + "default value of 'astropy_table' being assigned.",
                TdatFormatWarning,
            )
            lines.append("table_name = astropy_table")

        # loop through optional table keywords
        table_keywords = [
            "table_description",
            "table_document_url",
            "table_security",
        ]
        table_keywords = [kw for kw in table_keywords if kw in keywords.keys()]
        for key in table_keywords:
            if key == "table_description":
                if len(keywords[key]) > 80:
                    warn(
                        "'table_description' is too long, truncating to 80 characters",
                        TdatFormatWarning,
                    )
                new_desc = keywords.pop(key)[:80].replace("'", "")
                new_desc = new_desc.replace('"', "")
                lines.append(f'{key} = "{new_desc}"')
            elif key in keywords:
                lines.append(f"{key} = {keywords.pop(key)}")

        # add table columns as fields
        lines.append("#")
        lines.append("# Table Parameters")
        lines.append("#")
        for col in self.cols:
            if str(col.dtype) in self._dtype_dict:
                ctype = self._dtype_dict[str(col.dtype)]
            elif "int" in str(col.dtype):
                ctype = "int4"
            elif "float" in str(col.dtype):
                ctype = "float8"
            elif "<U" in str(col.dtype):
                ctype = f"char{str(col.dtype).rsplit('<U', maxsplit=1)[-1]}"
            else:
                raise TdatFormatError(
                    f"Unrecognized data type {col.dtype} for column {col.name}."
                )
            field_line = f"field[{col.name}] = {ctype}"

            if col.format is not None:
                field_line += f":{col.format}"
            if col.unit is not None:
                field_line += f"_{col.unit:vounit}"
            if "ucd" in col.meta:
                field_line += f" [{col.meta['ucd']}]"
            if "index" in col.meta:
                field_line += f" ({col.meta['index']})"
            elif (indices != []) and (col.name == indices[0]):
                field_line += " (key)"
            elif col.name in indices:
                field_line += " (index)"
            if col.description is not None:
                field_line += f" // {col.description}"
            elif "comment" in col.meta:
                field_line += " //"
            if "comment" in col.meta:
                field_line += f" // {col.meta['comment']}"
            lines.append(field_line)

        if len(keywords) != 0:
            if "parameter_defaults" in keywords:
                lines.append("#")
                lines.append(
                    f"{'parameter_defaults'} = {keywords.pop('parameter_defaults')}"
                )
        if len(keywords) != 0:
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
    """Splitter for tdat data.

    Handles the (deprecated) cases of multiple data delimiters, record
    delimiters, and multi-line records.
    Multiple data delimiters - Multiple delimiters can be specified in the
        header, e.g. field_delimiter = "|!" would treat both | and !
        as individual delimiters. Default: "|"
    Record Delimiters - The record_delimiter can be specified in the header. By
        default there is no record delimiter and new records should be set on
        new lines. The following list standard escaped character sequences and their
        equivalent meanings can be used:
        * \t (tab)
        * \b (backspace)
        * \r (carriage return)
        * \f (form feed)
        * \v (vertical tab)
        * \a (audible alert/bell),
        * \\### (where ### is a number between 1 and 127 and represents the
        ASCII character with that numerical code).
        Note: Specifying a record delimiter value of "" is interpreted as a
            single blank line between records.
    Multi-line records - A single record may take more than one line, indicated in the header

    """

    delimiter = "|"
    record_delimiter = None

    def preprocess_data_lines(self, lines, record_delimiter):
        """Split lines into multiple lines if a record_delimiter is specified."""
        data_lines = []
        for line in lines:
            data_lines += re.split(rf"[{record_delimiter}]", line)
        return data_lines

    def __call__(self, lines, field_delimiter="|", record_delimiter=None, nlines=1):
        """ """
        if record_delimiter is not None:
            lines = self.preprocess_data_lines(lines, record_delimiter)
        if hasattr(self, "process_line"):
            lines = (self.process_line(x) for x in lines)

        iline = 0
        _lines = []
        for line in lines:
            _lines.append(line)
            iline += 1
            if iline == nlines:
                vals = [
                    val
                    for _line in _lines
                    for val in re.split(rf"[{field_delimiter}]", _line)[:-1]
                ]
                iline = 0
                _lines = []
                if hasattr(self, "process_val"):
                    yield [self.process_val(x) for x in vals]
                else:
                    yield vals
            else:
                continue

    def join(self, vals):
        delimiter = getattr(self, "delimiter", "|")
        # TDAT specification requires a delimiter at the end of each record
        return delimiter.join(str(x) for x in vals) + delimiter


class TdatData(core.BaseData):
    """Data Reader for TDAT format"""

    comment = r"\s*(#|//)"
    write_comment = "# "
    splitter_class = TdatDataSplitter

    def get_data_lines(self, lines):
        """
        READ: Override the default get_data_lines to find start and end lines.
        """
        # Select lines between <DATA> and <END> in file
        try:
            start_line = (
                min(i for i, line in enumerate(lines) if line.strip() == "<DATA>") + 1
            )
        except ValueError:
            raise TdatFormatError("<DATA> not found in file." + _STD_MSG)
        try:
            end_line = min(i for i, line in enumerate(lines) if line.strip() == "<END>")
        except ValueError:
            raise TdatFormatError("<END> not found in file." + _STD_MSG)

        self.data_lines = self.process_lines(lines[start_line:end_line])

    def get_str_vals(self):
        """
        READ: Override the default get_str_vals to handle TDAT delimiters.
        """
        # if we have  delimiter from the header, use it.
        field_delimiter = getattr(self.header, "delimiter", self.splitter.delimiter)
        record_delimiter = getattr(
            self.header, "record_delimiter", self.splitter.record_delimiter
        )
        nlines = len(self.header._line_fields)
        return self.splitter(self.data_lines, field_delimiter, record_delimiter, nlines)

    def write(self, lines):
        """
        WRITE: Override the default write to include <DATA> and <END> lines.
        """
        lines.append("<DATA>")
        super().write(lines)
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
        indices = []
        for col, out_col in zip(cols, out.columns.values()):
            for attr in ("format", "unit", "description"):
                if hasattr(col, attr):
                    setattr(out_col, attr, getattr(col, attr))
            if hasattr(col, "meta"):
                out_col.meta.update(col.meta)
                if "index" in col.meta:
                    if col.meta["index"] == "key":
                        indices.insert(0, col.name)
                    else:
                        indices.append(col.name)
        # Add indices, if specified
        if len(indices) > 0:
            for name in indices:
                out.add_index(name)
        return out


class Tdat(core.BaseReader):
    """TDAT format
    See: https://heasarc.gsfc.nasa.gov/docs/software/dbdocs/tdat.html"

    Example::
      <HEADER>
      # # and // are comments
      table_name = example_table
      table_description = "Example table"
      #
      # Table Parameters
      #
      field[id] = integer [meta.id] (key) // Unique ID
      field[ra] = float:.4f_degree [pos.eq.ra] (index) // Right Ascension
      field[name] = char12 [meta.id] // Name
      #
      # Virtual Parameters
      #
      table_author = Example et al.
      #
      # Data Format Specification
      #
      line[1] = id name ra
      <DATA>
      1|TargetOne|1.0|
      2|TargetTwo|2.0|
      <END>

    The comments and keywords defined in the header, excepting common header
    section titles and blank comments, are available via the output table
    ``meta`` attribute::

      >>> from astropy.io import ascii
      >>> lines = ascii.tdat.make_example_data()
      >>> data = ascii.read(lines, format='tdat')
      >>> print(data.meta['comments'])
      ['# and // are comments']
      >>> for name, keyword in data.meta['keywords'].items():
      ...   print(name, keyword['value'])
      table_name example_table
      table_description "Example table"
      table_author Example et al.


    When writing to the TDAT format, the header will be auto-populated by
    information in the Table, prioritizing information given in the Table.meta:
        comments : list or string, (optional)
            Table information which provide context. This information is
            included in the header preceding all other lines and commented
            out with #
        keywords : OrderedDict, (optional, recommended)
            Header keywords which will appear in the file as "name=value" lines.
            Of particular importance are `table_name`, `table_description`,
            and `table_document_url`.
            `table_name` is a required keyword for the TDAT format and will
            be autopopulated with "astropy_table" if not specified in
            Table.meta["keywords"]
    If there is no Table.meta, this writer will attempt to automatically
    generate the appropriate header information based on the table and
    column properties and the recommendations for the TDAT format by HEASARC.

    Example::
    >>> from astropy.table import Table
    >>> from io import StringIO
    >>> t = Table(names=('reference_id', 'RA', 'Name'),
    ...           data=[[1, 2, 3], [1.0, 2.0, 3.0], ['c', 'd', 'e']])
    >>> out = StringIO()
    >>> t.write(out, format="ascii.tdat")
    >>> out.getvalue().splitlines()
    ['<HEADER>',
     'table_name = astropy_table',
     'table_description = "A table created via astropy"',
     '#',
     '# Table Parameters',
     '#',
     'field[reference_id] = integer',
     'field[RA] = float',
     'field[Name] = char1',
     '#',
     '# Data Format Specification',
     '#',
     'line[1] = reference_id RA Name',
     '#',
     '<DATA>',
     '1|1.0|c|',
     '2|2.0|d|',
     '3|3.0|e|',
     '<END>']

    Including relevant metadata for the table and columns separately
    is possible with a mixture of attribute assignment and additions to the
    metadata::
    >>> from astropy.table import Table
    >>> from io import StringIO
    >>> t = Table(names=('reference_id', 'RA', 'Name'),
    ...           data=[[1, 2, 3], [1.0, 2.0, 3.0], ['c', 'd', 'e']])
    >>> t.meta["table_name"] = "example_table"
    >>> t.meta["table_description"] = "An example table for the tdat writer."
    >>> t.add_index('reference_id')
    >>> t.columns['reference_id'].meta['comment'] = "For internal reference only"
    >>> t.add_index('RA')
    >>> t.columns['RA'].unit = "degree"
    >>> t.columns['RA'].format = ".4f"
    >>> t.columns['RA'].meta['ucd'] = "pos.eq.ra"
    >>> t.columns['Name'].description = "The name of the source (if available)"
    >>> out = StringIO()
    >>> t.write(out, format="ascii.tdat")
    >>> out.getvalue().splitlines()
    ['<HEADER>',
     'table_name = example_table',
     'table_description = An example table for the tdat writer.',
     '#',
     '# Table Parameters',
     '#',
     'field[reference_id] = integer (key) // // For internal reference only',
     'field[RA] = float:.4f_deg [pos.eq.ra] (index)',
     'field[Name] = char1 // The name of the source (if available)',
     '#',
     '# Virtual Parameters',
     '#',
     'table_name = example_table',
     '#',
     '# Data Format Specification',
     '#',
     'line[1] = reference_id RA Name',
     '#',
     '<DATA>',
     '1|1.0000|c|',
     '2|2.0000|d|',
     '3|3.0000|e|',
     '<END>']
    """

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
