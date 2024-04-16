# Licensed under a 3-clause BSD style license - see LICENSE.rst

import copy
import re

from astropy.table import Table
from astropy.utils.misc import _set_locale

from . import core, cparser


class FastBasic(metaclass=core.MetaBaseReader):
    """
    This class is intended to handle the same format addressed by the
    ordinary :class:`Basic` writer, but it acts as a wrapper for underlying C
    code and is therefore much faster. Unlike the other ASCII readers and
    writers, this class is not very extensible and is restricted
    by optimization requirements.
    """

    _format_name = "fast_basic"
    _description = "Basic table with custom delimiter using the fast C engine"
    _fast = True
    fill_extra_cols = False
    guessing = False
    strict_names = False

    def __init__(self, default_kwargs={}, **user_kwargs):
        # Make sure user does not set header_start to None for a reader
        # that expects a non-None value (i.e. a number >= 0).  This mimics
        # what happens in the Basic reader.
        if (
            default_kwargs.get("header_start", 0) is not None
            and user_kwargs.get("header_start", 0) is None
        ):
            raise ValueError("header_start cannot be set to None for this Reader")

        # Set up kwargs and copy any user kwargs.  Use deepcopy user kwargs
        # since they may contain a dict item which would end up as a ref to the
        # original and get munged later (e.g. in cparser.pyx validation of
        # fast_reader dict).
        kwargs = copy.deepcopy(default_kwargs)
        kwargs.update(copy.deepcopy(user_kwargs))

        delimiter = kwargs.pop("delimiter", " ")
        self.delimiter = str(delimiter) if delimiter is not None else None
        self.write_comment = kwargs.get("comment", "# ")
        self.comment = kwargs.pop("comment", "#")
        if self.comment is not None:
            self.comment = str(self.comment)
        self.quotechar = str(kwargs.pop("quotechar", '"'))
        self.header_start = kwargs.pop("header_start", 0)
        # If data_start is not specified, start reading
        # data right after the header line
        data_start_default = user_kwargs.get(
            "data_start", self.header_start + 1 if self.header_start is not None else 1
        )
        self.data_start = kwargs.pop("data_start", data_start_default)
        self.kwargs = kwargs
        self.strip_whitespace_lines = True
        self.strip_whitespace_fields = True

    def _read_header(self):
        # Use the tokenizer by default -- this method
        # can be overridden for specialized headers
        self.engine.read_header()

    def read(self, table):
        """
        Read input data (file-like object, filename, list of strings, or
        single string) into a Table and return the result.
        """
        if self.comment is not None and len(self.comment) != 1:
            raise core.ParameterError("The C reader does not support a comment regex")
        elif self.data_start is None:
            raise core.ParameterError(
                "The C reader does not allow data_start to be None"
            )
        elif (
            self.header_start is not None
            and self.header_start < 0
            and not isinstance(self, FastCommentedHeader)
        ):
            raise core.ParameterError(
                "The C reader does not allow header_start to be "
                "negative except for commented-header files"
            )
        elif self.data_start < 0:
            raise core.ParameterError(
                "The C reader does not allow data_start to be negative"
            )
        elif len(self.delimiter) != 1:
            raise core.ParameterError("The C reader only supports 1-char delimiters")
        elif len(self.quotechar) != 1:
            raise core.ParameterError(
                "The C reader only supports a length-1 quote character"
            )
        elif "converters" in self.kwargs:
            raise core.ParameterError(
                "The C reader does not support passing specialized converters"
            )
        elif "encoding" in self.kwargs:
            raise core.ParameterError(
                "The C reader does not use the encoding parameter"
            )
        elif "outputter_cls" in self.kwargs:
            raise core.ParameterError(
                "The C reader does not use the outputter_cls parameter"
            )
        elif "inputter_cls" in self.kwargs:
            raise core.ParameterError(
                "The C reader does not use the inputter_cls parameter"
            )
        elif "data_splitter_cls" in self.kwargs or "header_splitter_cls" in self.kwargs:
            raise core.ParameterError("The C reader does not use a Splitter class")

        self.strict_names = self.kwargs.pop("strict_names", False)

        # Process fast_reader kwarg, which may or may not exist (though ui.py will always
        # pass this as a dict with at least 'enable' set).
        fast_reader = self.kwargs.get("fast_reader", True)
        if not isinstance(fast_reader, dict):
            fast_reader = {}

        fast_reader.pop("enable", None)
        self.return_header_chars = fast_reader.pop("return_header_chars", False)
        # Put fast_reader dict back into kwargs.
        self.kwargs["fast_reader"] = fast_reader

        self.engine = cparser.CParser(
            table,
            self.strip_whitespace_lines,
            self.strip_whitespace_fields,
            delimiter=self.delimiter,
            header_start=self.header_start,
            comment=self.comment,
            quotechar=self.quotechar,
            data_start=self.data_start,
            fill_extra_cols=self.fill_extra_cols,
            **self.kwargs,
        )
        conversion_info = self._read_header()
        self.check_header()
        if conversion_info is not None:
            try_int, try_float, try_string = conversion_info
        else:
            try_int = {}
            try_float = {}
            try_string = {}

        with _set_locale("C"):
            data, comments = self.engine.read(try_int, try_float, try_string)
        out = self.make_table(data, comments)

        if self.return_header_chars:
            out.meta["__ascii_fast_reader_header_chars__"] = self.engine.header_chars

        return out

    def make_table(self, data, comments):
        """Actually make the output table give the data and comments."""
        meta = {}
        if comments:
            meta["comments"] = comments

        names = core._deduplicate_names(self.engine.get_names())
        return Table(data, names=names, meta=meta)

    def check_header(self):
        names = self.engine.get_header_names() or self.engine.get_names()
        if self.strict_names:
            # Impose strict requirements on column names (normally used in guessing)
            bads = [" ", ",", "|", "\t", "'", '"']
            for name in names:
                if (
                    core._is_number(name)
                    or len(name) == 0
                    or name[0] in bads
                    or name[-1] in bads
                ):
                    raise ValueError(
                        f"Column name {name!r} does not meet strict name requirements"
                    )
        # When guessing require at least two columns
        if self.guessing and len(names) <= 1:
            raise ValueError(
                f"Table format guessing requires at least two columns, got {names}"
            )

    def write(self, table, output):
        """
        Use a fast Cython method to write table data to output,
        where output is a filename or file-like object.
        """
        self._write(table, output, {})

    def _write(
        self, table, output, default_kwargs, header_output=True, output_types=False
    ):
        # Fast writer supports only 1-d columns
        core._check_multidim_table(table, max_ndim=1)

        write_kwargs = {
            "delimiter": self.delimiter,
            "quotechar": self.quotechar,
            "strip_whitespace": self.strip_whitespace_fields,
            "comment": self.write_comment,
        }
        write_kwargs.update(default_kwargs)
        # user kwargs take precedence over default kwargs
        write_kwargs.update(self.kwargs)
        writer = cparser.FastWriter(table, **write_kwargs)
        writer.write(output, header_output, output_types)


class FastCsv(FastBasic):
    """
    A faster version of the ordinary :class:`Csv` writer that uses the
    optimized C parsing engine. Note that this reader will append empty
    field values to the end of any row with not enough columns, while
    :class:`FastBasic` simply raises an error.
    """

    _format_name = "fast_csv"
    _description = "Comma-separated values table using the fast C engine"
    _fast = True
    fill_extra_cols = True

    def __init__(self, **kwargs):
        super().__init__({"delimiter": ",", "comment": None}, **kwargs)

    def write(self, table, output):
        """
        Override the default write method of `FastBasic` to
        output masked values as empty fields.
        """
        self._write(table, output, {"fill_values": [(core.masked, "")]})


class FastTab(FastBasic):
    """
    A faster version of the ordinary :class:`Tab` reader that uses
    the optimized C parsing engine.
    """

    _format_name = "fast_tab"
    _description = "Tab-separated values table using the fast C engine"
    _fast = True

    def __init__(self, **kwargs):
        super().__init__({"delimiter": "\t"}, **kwargs)
        self.strip_whitespace_lines = False
        self.strip_whitespace_fields = False


class FastNoHeader(FastBasic):
    """
    This class uses the fast C engine to read tables with no header line. If
    the names parameter is unspecified, the columns will be autonamed with
    "col{}".
    """

    _format_name = "fast_no_header"
    _description = "Basic table with no headers using the fast C engine"
    _fast = True

    def __init__(self, **kwargs):
        super().__init__({"header_start": None, "data_start": 0}, **kwargs)

    def write(self, table, output):
        """
        Override the default writing behavior in `FastBasic` so
        that columns names are not included in output.
        """
        self._write(table, output, {}, header_output=None)


class FastCommentedHeader(FastBasic):
    """
    A faster version of the :class:`CommentedHeader` reader, which looks for
    column names in a commented line. ``header_start`` denotes the index of
    the header line among all commented lines and is 0 by default.
    """

    _format_name = "fast_commented_header"
    _description = "Columns name in a commented line using the fast C engine"
    _fast = True

    def __init__(self, **kwargs):
        super().__init__({}, **kwargs)
        # Mimic CommentedHeader's behavior in which data_start
        # is relative to header_start if unspecified; see #2692
        if "data_start" not in kwargs:
            self.data_start = 0

    def make_table(self, data, comments):
        """
        Actually make the output table give the data and comments.  This is
        slightly different from the base FastBasic method in the way comments
        are handled.
        """
        meta = {}
        if comments:
            idx = self.header_start
            if idx < 0:
                idx = len(comments) + idx
            meta["comments"] = comments[:idx] + comments[idx + 1 :]
            if not meta["comments"]:
                del meta["comments"]

        names = core._deduplicate_names(self.engine.get_names())
        return Table(data, names=names, meta=meta)

    def _read_header(self):
        tmp = self.engine.source
        commented_lines = []

        for line in tmp.splitlines():
            line = line.lstrip()
            if line and line[0] == self.comment:  # line begins with a comment
                commented_lines.append(line[1:])
                if len(commented_lines) == self.header_start + 1:
                    break

        if len(commented_lines) <= self.header_start:
            raise cparser.CParserError("not enough commented lines")

        self.engine.setup_tokenizer([commented_lines[self.header_start]])
        self.engine.header_start = 0
        self.engine.read_header()
        self.engine.setup_tokenizer(tmp)

    def write(self, table, output):
        """
        Override the default writing behavior in `FastBasic` so
        that column names are commented.
        """
        self._write(table, output, {}, header_output="comment")


class FastRdb(FastBasic):
    """
    A faster version of the :class:`Rdb` reader. This format is similar to
    tab-delimited, but it also contains a header line after the column
    name line denoting the type of each column (N for numeric, S for string).
    """

    _format_name = "fast_rdb"
    _description = "Tab-separated with a type definition header line"
    _fast = True

    def __init__(self, **kwargs):
        super().__init__({"delimiter": "\t", "data_start": 2}, **kwargs)
        self.strip_whitespace_lines = False
        self.strip_whitespace_fields = False

    def _read_header(self):
        tmp = self.engine.source
        line1 = ""
        line2 = ""
        for line in tmp.splitlines():
            # valid non-comment line
            if not line1 and line.strip() and line.lstrip()[0] != self.comment:
                line1 = line
            elif not line2 and line.strip() and line.lstrip()[0] != self.comment:
                line2 = line
                break
        else:  # less than 2 lines in table
            raise ValueError("RDB header requires 2 lines")

        # Tokenize the two header lines separately.
        # Each call to self.engine.read_header by default
        #  - calls _deduplicate_names to ensure unique header_names
        #  - sets self.names from self.header_names if not provided as kwarg
        #  - applies self.include_names/exclude_names to self.names.
        # For parsing the types disable 1+3, but self.names needs to be set.
        self.engine.setup_tokenizer([line2])
        self.engine.header_start = 0
        self.engine.read_header(deduplicate=False, filter_names=False)
        types = self.engine.get_header_names()

        # If no kwarg names have been passed, reset to have column names read from header line 1.
        if types == self.engine.get_names():
            self.engine.set_names([])
        self.engine.setup_tokenizer([line1])
        # Get full list of column names prior to applying include/exclude_names,
        # which have to be applied to the unique name set after deduplicate.
        self.engine.read_header(deduplicate=True, filter_names=False)
        col_names = self.engine.get_names()
        self.engine.read_header(deduplicate=False)
        if len(col_names) != len(types):
            raise core.InconsistentTableError(
                "RDB header mismatch between number of column names and column types"
            )
        # If columns have been removed via include/exclude_names, extract matching types.
        if len(self.engine.get_names()) != len(types):
            types = [types[col_names.index(n)] for n in self.engine.get_names()]

        if any(not re.match(r"\d*(N|S)$", x, re.IGNORECASE) for x in types):
            raise core.InconsistentTableError(
                f"RDB type definitions do not all match [num](N|S): {types}"
            )

        try_int = {}
        try_float = {}
        try_string = {}

        for name, col_type in zip(self.engine.get_names(), types):
            if col_type[-1].lower() == "s":
                try_int[name] = 0
                try_float[name] = 0
                try_string[name] = 1
            else:
                try_int[name] = 1
                try_float[name] = 1
                try_string[name] = 0

        self.engine.setup_tokenizer(tmp)
        return (try_int, try_float, try_string)

    def write(self, table, output):
        """
        Override the default writing behavior in `FastBasic` to
        output a line with column types after the column name line.
        """
        self._write(table, output, {}, output_types=True)
