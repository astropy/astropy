# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""An extensible HTML table reader and writer.

html.py:
  Classes to read and write HTML tables

`BeautifulSoup <http://www.crummy.com/software/BeautifulSoup/>`_
must be installed to read HTML tables.
"""

import warnings
from copy import deepcopy

from astropy.table import Column
from astropy.utils.xml import writer

from . import core


class SoupString(str):
    """
    Allows for strings to hold BeautifulSoup data.
    """

    def __new__(cls, *args, **kwargs):
        return str.__new__(cls, *args, **kwargs)

    def __init__(self, val):
        self.soup = val


class ListWriter:
    """
    Allows for XMLWriter to write to a list instead of a file.
    """

    def __init__(self, out):
        self.out = out

    def write(self, data):
        self.out.append(data)


def identify_table(soup, htmldict, numtable):
    """
    Checks whether the given BeautifulSoup tag is the table
    the user intends to process.
    """
    if soup is None or soup.name != "table":
        return False  # Tag is not a <table>

    elif "table_id" not in htmldict:
        return numtable == 1
    table_id = htmldict["table_id"]

    if isinstance(table_id, str):
        return "id" in soup.attrs and soup["id"] == table_id
    elif isinstance(table_id, int):
        return table_id == numtable

    # Return False if an invalid parameter is given
    return False


class HTMLInputter(core.BaseInputter):
    """
    Input lines of HTML in a valid form.

    This requires `BeautifulSoup
    <http://www.crummy.com/software/BeautifulSoup/>`_ to be installed.
    """

    def process_lines(self, lines):
        """
        Convert the given input into a list of SoupString rows
        for further processing.
        """
        try:
            from bs4 import BeautifulSoup
        except ImportError:
            raise core.OptionalTableImportError(
                "BeautifulSoup must be installed to read HTML tables"
            )

        if "parser" not in self.html:
            with warnings.catch_warnings():
                # Ignore bs4 parser warning #4550.
                warnings.filterwarnings(
                    "ignore", ".*no parser was explicitly specified.*"
                )
                soup = BeautifulSoup("\n".join(lines))
        else:  # use a custom backend parser
            soup = BeautifulSoup("\n".join(lines), self.html["parser"])
        tables = soup.find_all("table")
        for i, possible_table in enumerate(tables):
            if identify_table(possible_table, self.html, i + 1):
                table = possible_table  # Find the correct table
                break
        else:
            if isinstance(self.html["table_id"], int):
                err_descr = f"number {self.html['table_id']}"
            else:
                err_descr = f"id '{self.html['table_id']}'"
            raise core.InconsistentTableError(
                f"ERROR: HTML table {err_descr} not found"
            )

        # Get all table rows
        soup_list = [SoupString(x) for x in table.find_all("tr")]

        return soup_list


class HTMLSplitter(core.BaseSplitter):
    """
    Split HTML table data.
    """

    def __call__(self, lines):
        """
        Return HTML data from lines as a generator.
        """
        for line in lines:
            if not isinstance(line, SoupString):
                raise TypeError("HTML lines should be of type SoupString")
            soup = line.soup
            header_elements = soup.find_all("th")
            if header_elements:
                # Return multicolumns as tuples for HTMLHeader handling
                yield [
                    (el.text.strip(), el["colspan"])
                    if el.has_attr("colspan")
                    else el.text.strip()
                    for el in header_elements
                ]
            data_elements = soup.find_all("td")
            if data_elements:
                yield [el.text.strip() for el in data_elements]
        if len(lines) == 0:
            raise core.InconsistentTableError(
                "HTML tables must contain data in a <table> tag"
            )


class HTMLOutputter(core.TableOutputter):
    """
    Output the HTML data as an ``astropy.table.Table`` object.

    This subclass allows for the final table to contain
    multidimensional columns (defined using the colspan attribute
    of <th>).
    """

    default_converters = [
        core.convert_numpy(int),
        core.convert_numpy(float),
        core.convert_numpy(str),
    ]

    def __call__(self, cols, meta):
        """
        Process the data in multidimensional columns.
        """
        new_cols = []
        col_num = 0

        while col_num < len(cols):
            col = cols[col_num]
            if hasattr(col, "colspan"):
                # Join elements of spanned columns together into list of tuples
                span_cols = cols[col_num : col_num + col.colspan]
                new_col = core.Column(col.name)
                new_col.str_vals = list(zip(*[x.str_vals for x in span_cols]))
                new_cols.append(new_col)
                col_num += col.colspan
            else:
                new_cols.append(col)
                col_num += 1

        return super().__call__(new_cols, meta)


class HTMLHeader(core.BaseHeader):
    splitter_class = HTMLSplitter

    def start_line(self, lines):
        """
        Return the line number at which header data begins.
        """
        for i, line in enumerate(lines):
            if not isinstance(line, SoupString):
                raise TypeError("HTML lines should be of type SoupString")
            soup = line.soup
            if soup.th is not None:
                return i

        return None

    def _set_cols_from_names(self):
        """
        Set columns from header names, handling multicolumns appropriately.
        """
        self.cols = []
        new_names = []

        for name in self.names:
            if isinstance(name, tuple):
                col = core.Column(name=name[0])
                col.colspan = int(name[1])
                self.cols.append(col)
                new_names.append(name[0])
                for i in range(1, int(name[1])):
                    # Add dummy columns
                    self.cols.append(core.Column(""))
                    new_names.append("")
            else:
                self.cols.append(core.Column(name=name))
                new_names.append(name)

        self.names = new_names


class HTMLData(core.BaseData):
    splitter_class = HTMLSplitter

    def start_line(self, lines):
        """
        Return the line number at which table data begins.
        """
        for i, line in enumerate(lines):
            if not isinstance(line, SoupString):
                raise TypeError("HTML lines should be of type SoupString")
            soup = line.soup

            if soup.td is not None:
                if soup.th is not None:
                    raise core.InconsistentTableError(
                        "HTML tables cannot have headings and data in the same row"
                    )
                return i

        raise core.InconsistentTableError("No start line found for HTML data")

    def end_line(self, lines):
        """
        Return the line number at which table data ends.
        """
        last_index = -1

        for i, line in enumerate(lines):
            if not isinstance(line, SoupString):
                raise TypeError("HTML lines should be of type SoupString")
            soup = line.soup
            if soup.td is not None:
                last_index = i

        if last_index == -1:
            return None
        return last_index + 1


class HTML(core.BaseReader):
    """HTML format table.

    In order to customize input and output, a dict of parameters may
    be passed to this class holding specific customizations.

    **htmldict** : Dictionary of parameters for HTML input/output.

        * css : Customized styling
            If present, this parameter will be included in a <style>
            tag and will define stylistic attributes of the output.

        * table_id : ID for the input table
            If a string, this defines the HTML id of the table to be processed.
            If an integer, this specifies the index of the input table in the
            available tables. Unless this parameter is given, the reader will
            use the first table found in the input file.

        * multicol : Use multi-dimensional columns for output
            The writer will output tuples as elements of multi-dimensional
            columns if this parameter is true, and if not then it will
            use the syntax 1.36583e-13 .. 1.36583e-13 for output. If not
            present, this parameter will be true by default.

        * raw_html_cols : column name or list of names with raw HTML content
            This allows one to include raw HTML content in the column output,
            for instance to include link references in a table.  This option
            requires that the bleach package be installed.  Only whitelisted
            tags are allowed through for security reasons (see the
            raw_html_clean_kwargs arg).

        * raw_html_clean_kwargs : dict of keyword args controlling HTML cleaning
            Raw HTML will be cleaned to prevent unsafe HTML from ending up in
            the table output.  This is done by calling ``bleach.clean(data,
            **raw_html_clean_kwargs)``.  For details on the available options
            (e.g. tag whitelist) see:
            https://bleach.readthedocs.io/en/latest/clean.html

        * parser : Specific HTML parsing library to use
            If specified, this specifies which HTML parsing library
            BeautifulSoup should use as a backend. The options to choose
            from are 'html.parser' (the standard library parser), 'lxml'
            (the recommended parser), 'xml' (lxml's XML parser), and
            'html5lib'. html5lib is a highly lenient parser and therefore
            might work correctly for unusual input if a different parser
            fails.

        * jsfiles : list of js files to include when writing table.

        * cssfiles : list of css files to include when writing table.

        * js : js script to include in the body when writing table.

        * table_class : css class for the table

    """

    _format_name = "html"
    _io_registry_format_aliases = ["html"]
    _io_registry_suffix = ".html"
    _description = "HTML table"

    header_class = HTMLHeader
    data_class = HTMLData
    inputter_class = HTMLInputter

    max_ndim = 2  # HTML supports writing 2-d columns with shape (n, m)

    def __init__(self, htmldict={}):
        """
        Initialize classes for HTML reading and writing.
        """
        super().__init__()
        self.html = deepcopy(htmldict)
        if "multicol" not in htmldict:
            self.html["multicol"] = True
        if "table_id" not in htmldict:
            self.html["table_id"] = 1
        self.inputter.html = self.html

    def read(self, table):
        """
        Read the ``table`` in HTML format and return a resulting ``Table``.
        """
        self.outputter = HTMLOutputter()
        return super().read(table)

    def write(self, table):
        """
        Return data in ``table`` converted to HTML as a list of strings.
        """
        # Check that table has only 1-d or 2-d columns. Above that fails.
        self._check_multidim_table(table)

        cols = list(table.columns.values())

        self.data.header.cols = cols
        self.data.cols = cols

        if isinstance(self.data.fill_values, tuple):
            self.data.fill_values = [self.data.fill_values]

        self.data._set_fill_values(cols)
        self.data._set_col_formats()

        lines = []

        # Set HTML escaping to False for any column in the raw_html_cols input
        raw_html_cols = self.html.get("raw_html_cols", [])
        if isinstance(raw_html_cols, str):
            raw_html_cols = [raw_html_cols]  # Allow for a single string as input
        cols_escaped = [col.info.name not in raw_html_cols for col in cols]

        # Kwargs that get passed on to bleach.clean() if that is available.
        raw_html_clean_kwargs = self.html.get("raw_html_clean_kwargs", {})

        # Use XMLWriter to output HTML to lines
        w = writer.XMLWriter(ListWriter(lines))

        with w.tag("html"):
            with w.tag("head"):
                # Declare encoding and set CSS style for table
                with w.tag("meta", attrib={"charset": "utf-8"}):
                    pass
                with w.tag(
                    "meta",
                    attrib={
                        "http-equiv": "Content-type",
                        "content": "text/html;charset=UTF-8",
                    },
                ):
                    pass
                if "css" in self.html:
                    with w.tag("style"):
                        w.data(self.html["css"])
                if "cssfiles" in self.html:
                    for filename in self.html["cssfiles"]:
                        with w.tag(
                            "link", rel="stylesheet", href=filename, type="text/css"
                        ):
                            pass
                if "jsfiles" in self.html:
                    for filename in self.html["jsfiles"]:
                        with w.tag("script", src=filename):
                            # need this instead of pass to get <script></script>
                            w.data("")
            with w.tag("body"):
                if "js" in self.html:
                    with w.xml_cleaning_method("none"):
                        with w.tag("script"):
                            w.data(self.html["js"])
                if isinstance(self.html["table_id"], str):
                    html_table_id = self.html["table_id"]
                else:
                    html_table_id = None
                if "table_class" in self.html:
                    html_table_class = self.html["table_class"]
                    attrib = {"class": html_table_class}
                else:
                    attrib = {}
                with w.tag("table", id=html_table_id, attrib=attrib):
                    with w.tag("thead"):
                        with w.tag("tr"):
                            for col in cols:
                                if len(col.shape) > 1 and self.html["multicol"]:
                                    # Set colspan attribute for multicolumns
                                    w.start("th", colspan=col.shape[1])
                                else:
                                    w.start("th")
                                w.data(col.info.name.strip())
                                w.end(indent=False)
                        col_str_iters = []
                        new_cols_escaped = []

                        # Make a container to hold any new_col objects created
                        # below for multicolumn elements.  This is purely to
                        # maintain a reference for these objects during
                        # subsequent iteration to format column values.  This
                        # requires that the weakref info._parent be maintained.
                        new_cols = []

                        for col, col_escaped in zip(cols, cols_escaped):
                            if len(col.shape) > 1 and self.html["multicol"]:
                                span = col.shape[1]
                                for i in range(span):
                                    # Split up multicolumns into separate columns
                                    new_col = Column([el[i] for el in col])

                                    new_col_iter_str_vals = self.fill_values(
                                        col, new_col.info.iter_str_vals()
                                    )
                                    col_str_iters.append(new_col_iter_str_vals)
                                    new_cols_escaped.append(col_escaped)
                                    new_cols.append(new_col)
                            else:
                                col_iter_str_vals = self.fill_values(
                                    col, col.info.iter_str_vals()
                                )
                                col_str_iters.append(col_iter_str_vals)

                                new_cols_escaped.append(col_escaped)

                    for row in zip(*col_str_iters):
                        with w.tag("tr"):
                            for el, col_escaped in zip(row, new_cols_escaped):
                                # Potentially disable HTML escaping for column
                                method = "escape_xml" if col_escaped else "bleach_clean"
                                with w.xml_cleaning_method(
                                    method, **raw_html_clean_kwargs
                                ):
                                    w.start("td")
                                    w.data(el.strip())
                                    w.end(indent=False)

        # Fixes XMLWriter's insertion of unwanted line breaks
        return ["".join(lines)]

    def fill_values(self, col, col_str_iters):
        """
        Return an iterator of the values with replacements based on fill_values.
        """
        # check if the col is a masked column and has fill values
        is_masked_column = hasattr(col, "mask")
        has_fill_values = hasattr(col, "fill_values")

        for idx, col_str in enumerate(col_str_iters):
            if is_masked_column and has_fill_values:
                if col.mask[idx]:
                    yield col.fill_values[core.masked]
                    continue

            if has_fill_values:
                if col_str in col.fill_values:
                    yield col.fill_values[col_str]
                    continue

            yield col_str
