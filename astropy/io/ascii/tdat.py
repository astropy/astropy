# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""This package contains functions for reading and writing TDAT tables that are
not meant to be used directly, but instead are available as readers/writers in
`astropy.table`. See :ref:`astropy:table_io` for more details.

:Author: Daniel Giles (daniel.k.giles@gmail.com)
         Abdu Zoghbi

"""

import re
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

    @property
    def literals_dict(self):
        """Return a dictionary of placeholders to be used in place of
        backslashed delimiter characters.
        """
        return {c: f"literal{n}" for n, c in enumerate(self.delimiter)}

    def preprocess_data_lines(self, lines, record_delimiter):
        """Split lines into multiple lines if a record_delimiter is specified."""
        data_lines = []
        for line in lines:
            data_lines += re.split(rf"{record_delimiter}", line)
        return data_lines

    def process_line(self, line: str) -> str:
        """Remove whitespace at the beginning or end of line.  This is especially useful for
        whitespace-delimited files to prevent spurious columns at the beginning or end.

        READ: override default to handle backslashed delimiter characters.
        """
        for c in self.delimiter:
            line = line.replace(f"\\{c}", self.literals_dict[c])
        line = re.sub(rf"[{self.delimiter}]", "|", line)
        return line.strip()

    def __call__(self, lines, field_delimiter="|", record_delimiter=None, nlines=1):
        """ """
        self.delimiter = field_delimiter

        def replace_placeholders(line):
            for c in field_delimiter:
                line = line.replace(self.literals_dict[c], c)
            return line

        if " " in field_delimiter:
            warn(
                TdatFormatWarning("Double check your data when using space delimiters.")
            )

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
                    replace_placeholders(val)
                    for _line in _lines
                    for val in re.split(r"\|", _line)[:-1]
                ]
                # Reset
                iline = 0
                _lines = []
                if hasattr(self, "process_val"):
                    yield [self.process_val(x) for x in vals]
                else:  # pragma: no cover
                    yield vals
            else:  # pragma: no cover
                continue

    def join(self, vals):
        delimiter = getattr(self, "delimiter", "|")
        # TDAT specification requires a delimiter at the end of each record
        return delimiter.join(str(x) for x in vals) + delimiter


class TdatHeader(basic.BasicHeader):
    """TDAT table header"""

    splitter_class = TdatDataSplitter
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
    _keys = r"(?P<key>\w+)\s*=\s*(?P<value>.*?)(\s*(#|//)(?P<comment>.*$)|(\s*$))"
    # keywords in the header: name[text] = some_other_text;
    # names: relate|line
    _extra_keys = r"\s*(relate|line)\[(\w+)\]\s*=\s*([\w\s]+)(?:\((\w+)\))?"
    _field_line = r"\s*field\[\w+\]\s*=\.*"

    _deprecated_keys = r"\s*(record_delimiter|field_delimiter)\s*=\s\"(.+)\""
    _required_keywords = ("table_name",)
    _dtype_dict_in = {
        "int1": np.int8,
        "integer1": np.int8,
        "tinyint": np.int8,
        "int2": np.int16,
        "integer2": np.int16,
        "smallint": np.int16,
        "int4": np.int32,
        "integer4": np.int32,
        "integer": np.int32,
        "float4": np.float32,
        "real": np.float32,
        "float": float,
        "float8": float,
    }

    _dtype_dict_out = {
        "int": "int4",
        "int32": "int4",
        "int16": "int2",
        "int8": "int1",
        "float": "float8",
        "float64": "float8",
        "float32": "float4",
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
        line_fields = {}
        keywords = getattr(self, "_keywords", {})
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
                _quotes = r"^(?P<quote>[\"\'\`]*)(?P<value>.*?)(?P<endquote>[\"\'\`]*$)"
                quote_parse = re.compile(_quotes)
                quotes_match = quote_parse.match(kmatch.group("value"))
                if quotes_match.group("quote") != quotes_match.group("endquote")[::-1]:
                    raise TdatFormatError("Mismatched quotes for value of " + key)

                value = kmatch.group("value").strip("\"'`")
                if key == "table_name":
                    if len(value) > 20:
                        warn(
                            "The table_name has a maximum length of 20 characters, truncating.",
                            TdatFormatWarning,
                        )
                        value = value[:20]
                elif key == "table_description":
                    if len(value) > 80:
                        warn(
                            "The table_description has a maximum length of 80 characters, truncating.",
                            TdatFormatWarning,
                        )
                        value = value[:80]
                elif key == "table_security":
                    if value.lower() not in ("public", "private"):
                        warn(
                            "Value for table_security not recognized, should be public|private.",
                            TdatFormatWarning,
                        )
                keywords[key] = value

        self._line_fields = line_fields

        if len(keywords) > len(getattr(self, "_keywords", {})):
            self._keywords = keywords

    def update_meta(self, lines, meta):
        """Extract meta information: comments and key/values in the header
        READ: Overrides the default update_meta
        """
        try:
            start_line = (
                min(
                    i
                    for i, line in enumerate(lines)
                    if line.strip().upper() == "<HEADER>"
                )
                + 1
            )
        except ValueError:
            raise TdatFormatError("<HEADER> not found in file." + _STD_MSG)
        try:
            end_line = min(
                i for i, line in enumerate(lines) if line.strip().upper() == "<DATA>"
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
        start_line = (
            min(i for i, line in enumerate(lines) if line.strip().upper() == "<HEADER>")
            + 1
        )
        end_line = min(
            i for i, line in enumerate(lines) if line.strip().upper() == "<DATA>"
        )

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
            (?P<ctype>[^\W\s_]+(\(\d+\))?)
            (?:\:(?P<fmt>[^_\s\[\(/#]+))?
            (?:_(?P<unit>[^\[\(#]+))?\s*
            (?:\[(?P<ucd>[\w\.\;]+)\])?
            \s*
            (?:\((?P<index>\w+)\))?
            \s*
            (?:(//|\#)\s*(?P<desc>[^/#]*))?
            (?:(//|\#)\s*(?P<comment>.*))?
            \s*
            """,
            re.VERBOSE,
        )

        cols = {}
        keywords = getattr(self, "_keywords", {})
        for line in self.process_lines(lines):
            # look for field[..]= ... column definitions
            cmatch = col_parser.match(line)
            if cmatch:
                name = cmatch.group("name")
                if len(name) > 24:
                    warn(
                        "The field name must be shorter than 24 characters.",
                        TdatFormatWarning,
                    )
                col = core.Column(name=name)
                col.meta = {}

                ctype = cmatch.group("ctype")
                if ctype in self._dtype_dict_in:
                    col.dtype = self._dtype_dict_in[ctype]
                elif "char" in ctype:
                    col.dtype = str
                else:
                    raise TdatFormatError(
                        f"Unrecognized or unsupported data type {ctype} for {col.name}."
                    )
                col.unit = cmatch.group("unit")
                col.format = cmatch.group("fmt")
                if len(str(ctype) + str(col.format)) > 23:
                    # max 24 characters with ":" separator
                    raise TdatFormatError(
                        "The type:fmt specifier has a\
                        maximum length of 24 characters. The offending line is:\
                            \n{line}\n"
                        + _STD_MSG
                    )
                col.description = f"{cmatch.group('desc')}".strip()
                for val in ["comment", "ucd", "index"]:
                    if cmatch.group(val) is not None:
                        text = cmatch.group(val).strip()
                        if (val == "comment") and (len(text) > 80):
                            warn(
                                TdatFormatWarning(
                                    "Comments are limited to 80 characters or less, truncating."
                                )
                            )
                        col.meta[val] = text[:80]
                cols[col.name] = col

        self.names = [
            val for line_field in self._line_fields.values() for val in line_field
        ]

        # check that cols and _line_fields are consistent or throw an error
        if len(self.names) != len(cols):
            colnames = [col.name for col in cols.values()]
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
        # the reader can handle alternate delimiters, but the writer should
        # comply with the current convention
        if self.splitter.delimiter.strip() != "|":
            warn(
                TdatFormatWarning(
                    "Delimiters other than the pipe character, '|', are deprecated. Using '|'."
                )
            )
        self.splitter.delimiter = "|"

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

        indices = [col.info.name for col in self.cols if col.info.indices != []]

        if "table_name" in keywords:
            if len(keywords["table_name"]) > 20:
                warn(
                    "'table_name' is too long, truncating to 20 characters",
                    TdatFormatWarning,
                )
            lines.append(f"table_name = {keywords.pop('table_name')[:20]}")
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
        for kw in ["record_delimiter", "field_delimiter"]:
            val = keywords.pop(kw, None)
            if val is not None:
                warn(TdatFormatWarning(f"Skipping deprecated keyword: {kw}."))

        for key in table_keywords:
            if key == "table_description":
                if len(keywords[key]) > 80:
                    warn(
                        "'table_description' is too long, truncating to 80 characters",
                        TdatFormatWarning,
                    )
                new_desc = keywords.pop(key)
                new_desc = new_desc[:80]
                lines.append(f"{key} = {new_desc}")
            elif key == "table_security":
                if keywords.get(key).lower() not in ("public", "private"):
                    warn(
                        "Value for table_security not recognized, should be public|private.",
                        TdatFormatWarning,
                    )
                lines.append(f"{key} = {keywords.pop(key)}")
            elif key in keywords:
                lines.append(f"{key} = {keywords.pop(key)}")

        # add table columns as fields
        lines.append("#")
        lines.append("# Table Parameters")
        lines.append("#")
        for col in self.cols:
            if str(col_type := col.info.dtype) in self._dtype_dict_out:
                ctype = self._dtype_dict_out[str(col_type)]
            elif col_type.kind == "i":
                ctype = "int4"
            elif col_type.kind == "f":
                ctype = "float8"
            elif col_type.kind == "U":
                ctype = f"char{col_type.itemsize // 4}"
            else:
                raise TdatFormatError(
                    f'Unrecognized data type `{col_type}` for column "{col.info.name}".'
                )
            col_name = col.info.name
            if len(col_name) >= 24:
                warn(
                    "The field name must be shorter than 24 characters, truncating.",
                    TdatFormatWarning,
                )
                col_name = col_name[:23]
            field_line = f"field[{col_name}] = {ctype}"

            col_info_meta = col.info.meta or {}
            if col.info.format is not None:
                field_line += f":{col.info.format}"
            if col.info.unit is not None:
                field_line += f"_{col.info.unit:cds}"
            if "ucd" in col_info_meta:
                field_line += f" [{col_info_meta['ucd']}]"
            if "index" in col_info_meta:
                field_line += f" ({col_info_meta['index']})"
            elif (indices != []) and (col.info.name == indices[0]):
                field_line += " (key)"
            elif col.info.name in indices:
                field_line += " (index)"
            if col.info.description is not None and col.info.description != "None":
                field_line += f" // {col.info.description}"
            elif "comment" in col_info_meta:
                field_line += " //"
            if "comment" in col_info_meta:
                field_line += f" // {col_info_meta['comment']}"
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
        lines.append(f"line[1] = {' '.join([col.info.name for col in self.cols])}")
        lines.append("#")


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
                min(
                    i
                    for i, line in enumerate(lines)
                    if line.strip().upper() == "<DATA>"
                )
                + 1
            )
        except ValueError:
            raise TdatFormatError("<DATA> not found in file." + _STD_MSG)
        try:
            end_line = min(
                i for i, line in enumerate(lines) if line.strip().upper() == "<END>"
            )
        except ValueError:
            # <END> is an optional keyword to demarcate the end of the data.
            # If not present the end of the document is assumed to be the end
            # end of the data.
            end_line = len(lines) - 1

        self.data_lines = self.process_lines(lines[start_line:end_line])

    def str_vals(self):
        """WRITE: convert all values in table to a list of lists of strings.

        This sets the fill values and possibly column formats from the input
        formats={} keyword, then ends up calling table.pprint._pformat_col_iter()
        by a circuitous path. That function does the real work of formatting.
        Finally replace anything matching the fill_values.

        Returns
        -------
        values : list of list of str
        """
        self._set_fill_values(self.cols)
        self._set_col_formats()
        for col in self.cols:
            if np.issubdtype(col.dtype, np.integer):
                if np.any(col.data > 2**31 - 1) or np.any(col.data < -(2**31)):
                    raise TdatFormatError(
                        "Values cannot be converted to a TDAT compatible integer."
                    )
            col.str_vals = []
            for val in col.info.iter_str_vals():
                col.str_vals.append(val.replace("|", "\\|"))

        self._replace_vals(self.cols)
        return [col.str_vals for col in self.cols]

    def get_str_vals(self):
        """
        READ: Override the default get_str_vals to handle TDAT delimiters.
        """
        # if we have  delimiter from the header, use it.
        field_delimiter = getattr(self.header, "delimiter", self.splitter.delimiter)
        self.splitter.delimiter = field_delimiter
        record_delimiter = getattr(
            self.header, "record_delimiter", self.splitter.record_delimiter
        )
        allowed_delimiters = list(map(chr, range(1, 128)))
        if record_delimiter is not None:
            if record_delimiter not in allowed_delimiters:
                raise TdatFormatError(
                    f"Unsupported record delimiter: {record_delimiter}."
                )
            if record_delimiter == "|":
                raise TdatFormatError("'|' character reserved for field delimiter.")
        for c in field_delimiter:
            if c not in allowed_delimiters:
                raise TdatFormatError(f"Unsupported field delimiter: {c}.")
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

    See: https://heasarc.gsfc.nasa.gov/docs/software/dbdocs/tdat.html

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
      ...   print(name, keyword)
      table_name example_table
      table_description Example table
      table_author Example et al.

    When writing to the TDAT format, the header will be auto-populated by
    information in the Table, prioritizing information given in the Table.meta:

    **comments** : list or string, (optional)
      Table information which provide context. This information is
      included in the header preceding all other lines and commented
      out with #

    **keywords** : dict, (optional, recommended)
      Header keywords which will appear in the file as "name=value" lines.
      Of particular importance are table_name, table_description,
      and table_document_url.

    If there is no Table.meta, this writer will attempt to automatically
    generate the appropriate header information based on the table and
    column properties and the recommendations for the TDAT format by HEASARC.
    Column ``units`` are written using the CDS format.

    Example::

      >>> from astropy.table import Table
      >>> import sys
      >>> t = Table(names=('reference_id', 'RA', 'Name'),
      ...           data=[[1, 2, 3], [1.0, 2.0, 3.0], ['c', 'd', 'e']])
      >>> t.meta['table_name'] = "astropy_table"
      >>> t.write(sys.stdout, format="ascii.tdat")
      <HEADER>
      table_name = astropy_table
      #
      # Table Parameters
      #
      field[reference_id] = int4
      field[RA] = float8
      field[Name] = char1
      #
      # Data Format Specification
      #
      line[1] = reference_id RA Name
      #
      <DATA>
      1|1.0|c|
      2|2.0|d|
      3|3.0|e|
      <END>

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
      >>> t.write(sys.stdout, format="ascii.tdat")
      <HEADER>
      table_name = example_table
      table_description = An example table for the tdat writer.
      #
      # Table Parameters
      #
      field[reference_id] = int4 (key) // // For internal reference only
      field[RA] = float8:.4f_deg [pos.eq.ra] (index)
      field[Name] = char1 // The name of the source (if available)
      #
      # Data Format Specification
      #
      line[1] = reference_id RA Name
      #
      <DATA>
      1|1.0000|c|
      2|2.0000|d|
      3|3.0000|e|
      <END>
    """

    _format_name = "tdat"
    _description = "HEASARC tdat format"
    _io_registry_can_write = True
    _io_registry_suffix = ".tdat"

    header_class = TdatHeader
    data_class = TdatData
    outputter_class = TdatOutputter
