# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""An extensible ASCII table reader and writer.

latex.py:
  Classes to read and write LaTeX tables

:Copyright: Smithsonian Astrophysical Observatory (2011)
:Author: Tom Aldcroft (aldcroft@head.cfa.harvard.edu)
"""

from __future__ import annotations

import re
from typing import TYPE_CHECKING

from . import core

if TYPE_CHECKING:
    from collections.abc import Generator
    from re import Pattern
    from typing import ClassVar, Final

latexdicts = {
    "AA": {
        "tabletype": "table",
        "header_start": r"\hline \hline",
        "header_end": r"\hline",
        "data_end": r"\hline",
    },
    "doublelines": {
        "tabletype": "table",
        "header_start": r"\hline \hline",
        "header_end": r"\hline\hline",
        "data_end": r"\hline\hline",
    },
    "template": {
        "tabletype": "tabletype",
        "caption": "caption",
        "tablealign": "tablealign",
        "col_align": "col_align",
        "preamble": "preamble",
        "header_start": "header_start",
        "header_end": "header_end",
        "data_start": "data_start",
        "data_end": "data_end",
        "tablefoot": "tablefoot",
        "units": {"col1": "unit of col1", "col2": "unit of col2"},
    },
}


RE_COMMENT: Final[Pattern[str]] = re.compile(r"(?<!\\)%")  # % character but not \%


def add_dictval_to_list(adict, key, alist):
    """
    Add a value from a dictionary to a list.

    Parameters
    ----------
    adict : dictionary
    key : hashable
    alist : list
        List where value should be added
    """
    if key in adict:
        if isinstance(adict[key], str):
            alist.append(adict[key])
        else:
            alist.extend(adict[key])


def find_latex_line(lines: list[str], latex: str) -> int | None:
    """
    Find the first line which matches a pattern.

    Parameters
    ----------
    lines : list
        List of strings
    latex : str
        Search pattern

    Returns
    -------
    line_num : int, None
        Line number. Returns None, if no match was found

    """
    re_string = re.compile(r"\s*" + latex.replace("\\", "\\\\"))
    for i, line in enumerate(lines):
        if re_string.match(line):
            return i
    return None


class LatexInputter(core.BaseInputter):
    def process_lines(self, lines: list[str]) -> list[str]:
        return [lin.strip() for lin in lines]


class LatexSplitter(core.BaseSplitter):
    """Split LaTeX table data. Default delimiter is `&`."""

    delimiter = "&"

    def __call__(self, lines: list[str]) -> Generator[list[str], None, None]:
        last_line = RE_COMMENT.split(lines[-1])[0].strip()
        if not last_line.endswith(r"\\"):
            lines[-1] = last_line + r"\\"

        return super().__call__(lines)

    def process_line(self, line: str) -> str:
        """Remove whitespace at the beginning or end of line. Also remove
        \\ at end of line.
        """
        line = RE_COMMENT.split(line)[0].strip()
        if not line.endswith(r"\\"):
            raise core.InconsistentTableError(
                r"Lines in LaTeX table have to end with \\"
            )
        return line.removesuffix(r"\\")

    def process_val(self, val: str) -> str:
        """Remove whitespace and {} at the beginning or end of value."""
        val = val.strip()
        if val and (val[0] == "{") and (val[-1] == "}"):
            val = val[1:-1]
        return val

    def join(self, vals: list[str]) -> str:
        """Join values together and add a few extra spaces for readability."""
        delimiter = " " + self.delimiter + " "
        return delimiter.join(x.strip() for x in vals) + r" \\"


class LatexHeader(core.BaseHeader):
    """Class to read the header of Latex Tables."""

    header_start = r"\begin{tabular}"
    splitter_class = LatexSplitter

    def start_line(self, lines):
        line = find_latex_line(lines, self.header_start)
        if line is not None:
            return line + 1
        else:
            return None

    def _get_units(self) -> dict[str, str]:
        units = {}
        col_units = [col.info.unit for col in self.cols]
        for name, unit in zip(self.colnames, col_units):
            if unit:
                try:
                    units[name] = unit.to_string(format="latex_inline")
                except AttributeError:
                    units[name] = unit
        return units

    def write(self, lines):
        if "col_align" not in self.latex:
            self.latex["col_align"] = len(self.cols) * "c"
        if "tablealign" in self.latex:
            align = "[" + self.latex["tablealign"] + "]"
        else:
            align = ""
        if self.latex["tabletype"] is not None:
            lines.append(r"\begin{" + self.latex["tabletype"] + r"}" + align)
        add_dictval_to_list(self.latex, "preamble", lines)
        if "caption" in self.latex:
            lines.append(r"\caption{" + self.latex["caption"] + "}")
        lines.append(self.header_start + r"{" + self.latex["col_align"] + r"}")
        add_dictval_to_list(self.latex, "header_start", lines)
        lines.append(self.splitter.join(self.colnames))
        units = self._get_units()
        if "units" in self.latex:
            units.update(self.latex["units"])
        if units:
            lines.append(
                self.splitter.join([units.get(name, " ") for name in self.colnames])
            )
        add_dictval_to_list(self.latex, "header_end", lines)


class LatexData(core.BaseData):
    """Class to read the data in LaTeX tables."""

    data_start: ClassVar[str | None] = None
    data_end = r"\end{tabular}"
    splitter_class = LatexSplitter

    def start_line(self, lines):
        if self.data_start:
            return find_latex_line(lines, self.data_start)
        else:
            start = self.header.start_line(lines)
            if start is None:
                raise core.InconsistentTableError(r"Could not find table start")
            return start + 1

    def end_line(self, lines):
        if self.data_end:
            return find_latex_line(lines, self.data_end)
        else:
            return None

    def write(self, lines):
        add_dictval_to_list(self.latex, "data_start", lines)
        core.BaseData.write(self, lines)
        add_dictval_to_list(self.latex, "data_end", lines)
        lines.append(self.data_end)
        add_dictval_to_list(self.latex, "tablefoot", lines)
        if self.latex["tabletype"] is not None:
            lines.append(r"\end{" + self.latex["tabletype"] + "}")


class Latex(core.BaseReader):
    r"""LaTeX format table.

    This class implements some LaTeX specific commands.  Its main
    purpose is to write out a table in a form that LaTeX can compile. It
    is beyond the scope of this class to implement every possible LaTeX
    command, instead the focus is to generate a syntactically valid
    LaTeX tables.

    This class can also read simple LaTeX tables (one line per table
    row, no ``\multicolumn`` or similar constructs), specifically, it
    can read the tables that it writes.
    When reading, it will look for the Latex commands to start and end tabular
    data (``\begin{tabular}`` and ``\end{tabular}``). That means that
    those lines have to be present in the input file; the benefit is that this
    reader can be used on a LaTeX file with text, tables, and figures and it
    will read the first valid table.

    .. note:: **Units in LaTeX tables**

        The LaTeX writer will output units in the table if they are present in the
        column info::

            >>> import io
            >>> out = io.StringIO()
            >>> import sys
            >>> import astropy.units as u
            >>> from astropy.table import Table
            >>> t = Table({'v': [1, 2] * u.km/u.s, 'class': ['star', 'jet']})
            >>> t.write(out, format='ascii.latex')
            >>> print(out.getvalue())
            \begin{table}
            \begin{tabular}{cc}
            v & class \\
            $\mathrm{km\,s^{-1}}$ &  \\
            1.0 & star \\
            2.0 & jet \\
            \end{tabular}
            \end{table}

        However, it will fail to read a table with units. There are so
        many ways to write units in LaTeX (enclosed in parenthesis or square brackets,
        as a separate row are as part of the column headers, using plain text, LaTeX
        symbols etc. ) that it is not feasible to implement a
        general reader for this. If you need to read a table with units, you can
        skip reading the lines with units to just read the numerical values using the
        ``data_start`` parameter to set the first line where numerical data values appear::

            >>> Table.read(out.getvalue(), format='ascii.latex', data_start=4)
            <Table length=2>
               v    class
            float64  str4
            ------- -----
                1.0  star
                2.0   jet

        Alternatively, you can write a custom reader using your knowledge of the exact
        format of the units in that case, by extending this class.

    Reading a LaTeX table, the following keywords are accepted:

    **ignore_latex_commands** :
        Lines starting with these LaTeX commands will be treated as comments (i.e. ignored).

    When writing a LaTeX table, the some keywords can customize the
    format.  Care has to be taken here, because python interprets ``\\``
    in a string as an escape character.  In order to pass this to the
    output either format your strings as raw strings with the ``r``
    specifier or use a double ``\\\\``.

    Examples::

        caption = r'My table \label{mytable}'
        caption = 'My table \\\\label{mytable}'

    **latexdict** : Dictionary of extra parameters for the LaTeX output

        * tabletype : used for first and last line of table.
            The default is ``\\begin{table}``.  The following would generate a table,
            which spans the whole page in a two-column document::

                ascii.write(data, sys.stdout, format="latex",
                            latexdict={'tabletype': 'table*'})

            If ``None``, the table environment will be dropped, keeping only
            the ``tabular`` environment.

        * tablealign : positioning of table in text.
            The default is not to specify a position preference in the text.
            If, e.g. the alignment is ``ht``, then the LaTeX will be ``\\begin{table}[ht]``.

        * col_align : Alignment of columns
            If not present all columns will be centered.

        * caption : Table caption (string or list of strings)
            This will appear above the table as it is the standard in
            many scientific publications.  If you prefer a caption below
            the table, just write the full LaTeX command as
            ``latexdict['tablefoot'] = r'\caption{My table}'``

        * preamble, header_start, header_end, data_start, data_end, tablefoot: Pure LaTeX
            Each one can be a string or a list of strings. These strings
            will be inserted into the table without any further
            processing. See the examples below.

        * units : dictionary of strings
            Keys in this dictionary should be names of columns. If
            present, a line in the LaTeX table directly below the column
            names is added, which contains the values of the
            dictionary. Example::

              from astropy.io import ascii
              data = {'name': ['bike', 'car'], 'mass': [75,1200], 'speed': [10, 130]}
              ascii.write(data, format="latex",
                          latexdict={'units': {'mass': 'kg', 'speed': 'km/h'}})

            If the column has no entry in the ``units`` dictionary, it defaults
            to the **unit** attribute of the column. If this attribute is not
            specified (i.e. it is None), the unit will be written as ``' '``.

        Run the following code to see where each element of the
        dictionary is inserted in the LaTeX table::

            from astropy.io import ascii
            data = {'cola': [1,2], 'colb': [3,4]}
            ascii.write(data, format="latex", latexdict=ascii.latex.latexdicts['template'])

        Some table styles are predefined in the dictionary
        ``ascii.latex.latexdicts``. The following generates in table in
        style preferred by A&A and some other journals::

            ascii.write(data, format="latex", latexdict=ascii.latex.latexdicts['AA'])

        As an example, this generates a table, which spans all columns
        and is centered on the page::

            ascii.write(data, format="latex", col_align='|lr|',
                        latexdict={'preamble': r'\begin{center}',
                                   'tablefoot': r'\end{center}',
                                   'tabletype': 'table*'})

    **caption** : Set table caption
        Shorthand for::

            latexdict['caption'] = caption

    **col_align** : Set the column alignment.
        If not present this will be auto-generated for centered
        columns. Shorthand for::

            latexdict['col_align'] = col_align

    """

    _format_name = "latex"
    _io_registry_format_aliases = ["latex"]
    _io_registry_suffix = ".tex"
    _description = "LaTeX table"

    header_class = LatexHeader
    data_class = LatexData
    inputter_class = LatexInputter

    # Strictly speaking latex only supports 1-d columns so this should inherit
    # the base max_ndim = 1. But as reported in #11695 this causes a strange
    # problem with Jupyter notebook, which displays a table by first calling
    # _repr_latex_. For a multidimensional table this issues a stack traceback
    # before moving on to _repr_html_. Here we prioritize fixing the issue with
    # Jupyter displaying a Table with multidimensional columns.
    max_ndim = None

    def __init__(
        self,
        ignore_latex_commands=[
            "hline",
            "vspace",
            "tableline",
            "toprule",
            "midrule",
            "bottomrule",
        ],
        latexdict={},
        caption="",
        col_align=None,
    ):
        super().__init__()

        self.latex = {}
        # The latex dict drives the format of the table and needs to be shared
        # with data and header
        self.header.latex = self.latex
        self.data.latex = self.latex
        self.latex["tabletype"] = "table"
        self.latex.update(latexdict)
        if caption:
            self.latex["caption"] = caption
        if col_align:
            self.latex["col_align"] = col_align

        self.ignore_latex_commands = ignore_latex_commands
        self.header.comment = "%|" + "|".join(
            [r"\\" + command for command in self.ignore_latex_commands]
        )
        self.data.comment = self.header.comment

    def write(self, table=None):
        self.header.start_line = None
        self.data.start_line = None
        return core.BaseReader.write(self, table=table)


class AASTexHeaderSplitter(LatexSplitter):
    r"""Extract column names from a `deluxetable`_.

    This splitter expects the following LaTeX code **in a single line**:

        \tablehead{\colhead{col1} & ... & \colhead{coln}}
    """

    def __call__(self, lines: list[str]) -> Generator[list[str], None, None]:
        return super(LatexSplitter, self).__call__(lines)

    def process_line(self, line: str) -> str:
        """extract column names from tablehead."""
        line = line.split("%")[0]
        line = line.replace(r"\tablehead", "")
        line = line.strip()
        if (line[0] == "{") and (line[-1] == "}"):
            line = line[1:-1]
        else:
            raise core.InconsistentTableError(r"\tablehead is missing {}")
        return line.replace(r"\colhead", "")

    def join(self, vals: list[str]) -> str:
        return " & ".join([r"\colhead{" + str(x) + "}" for x in vals])


class AASTexHeader(LatexHeader):
    r"""In a `deluxetable
    <http://fits.gsfc.nasa.gov/standard30/deluxetable.sty>`_ some header
    keywords differ from standard LaTeX.

    This header is modified to take that into account.
    """

    header_start = r"\tablehead"
    splitter_class = AASTexHeaderSplitter

    def start_line(self, lines):
        return find_latex_line(lines, r"\tablehead")

    def write(self, lines):
        if "col_align" not in self.latex:
            self.latex["col_align"] = len(self.cols) * "c"
        if "tablealign" in self.latex:
            align = "[" + self.latex["tablealign"] + "]"
        else:
            align = ""
        lines.append(
            r"\begin{"
            + self.latex["tabletype"]
            + r"}{"
            + self.latex["col_align"]
            + r"}"
            + align
        )
        add_dictval_to_list(self.latex, "preamble", lines)
        if "caption" in self.latex:
            lines.append(r"\tablecaption{" + self.latex["caption"] + "}")
        tablehead = " & ".join([r"\colhead{" + name + "}" for name in self.colnames])
        units = self._get_units()
        if "units" in self.latex:
            units.update(self.latex["units"])
        if units:
            tablehead += r"\\ " + self.splitter.join(
                [units.get(name, " ") for name in self.colnames]
            )
        lines.append(r"\tablehead{" + tablehead + "}")


class AASTexData(LatexData):
    r"""In a `deluxetable`_ the data is enclosed in `\startdata` and `\enddata`."""

    data_start = r"\startdata"
    data_end = r"\enddata"

    def start_line(self, lines):
        return find_latex_line(lines, self.data_start) + 1

    def write(self, lines):
        lines.append(self.data_start)
        lines_length_initial = len(lines)
        core.BaseData.write(self, lines)
        # To remove extra space(s) and // appended which creates an extra new line
        # in the end.
        if len(lines) > lines_length_initial:
            lines[-1] = re.sub(r"\s* \\ \\ \s* $", "", lines[-1], flags=re.VERBOSE)
        lines.append(self.data_end)
        add_dictval_to_list(self.latex, "tablefoot", lines)
        lines.append(r"\end{" + self.latex["tabletype"] + r"}")


class AASTex(Latex):
    """AASTeX format table.

    This class implements some AASTeX specific commands.
    AASTeX is used for the AAS (American Astronomical Society)
    publications like ApJ, ApJL and AJ.

    It derives from the ``Latex`` reader and accepts the same
    keywords.  However, the keywords ``header_start``, ``header_end``,
    ``data_start`` and ``data_end`` in ``latexdict`` have no effect.
    """

    _format_name = "aastex"
    _io_registry_format_aliases = ["aastex"]
    _io_registry_suffix = ""  # AASTex inherits from Latex, so override this class attr
    _description = "AASTeX deluxetable used for AAS journals"

    header_class = AASTexHeader
    data_class = AASTexData

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        # check if tabletype was explicitly set by the user
        if not (("latexdict" in kwargs) and ("tabletype" in kwargs["latexdict"])):
            self.latex["tabletype"] = "deluxetable"
