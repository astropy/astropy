# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This package contains functions for reading and writing QDP tables that are
not meant to be used directly, but instead are available as readers/writers in
`astropy.table`. See :ref:`astropy:table_io` for more details.
"""
import re
import copy
from collections.abc import Iterable
import numpy as np
import warnings
from astropy.utils.exceptions import AstropyUserWarning
from astropy.table import Table

from . import core, basic


def _line_type(line, delimiter=None):
    """Interpret a QDP file line

    Parameters
    ----------
    line : str
        a single line of the file

    Returns
    -------
    type : str
        Line type: "comment", "command", or "data"

    Examples
    --------
    >>> _line_type("READ SERR 3")
    'command'
    >>> _line_type(" \\n    !some gibberish")
    'comment'
    >>> _line_type("   ")
    'comment'
    >>> _line_type(" 21345.45")
    'data,1'
    >>> _line_type(" 21345.45 1.53e-3 1e-3 .04 NO nan")
    'data,6'
    >>> _line_type(" 21345.45,1.53e-3,1e-3,.04,NO,nan", delimiter=',')
    'data,6'
    >>> _line_type(" 21345.45 ! a comment to disturb")
    'data,1'
    >>> _line_type("NO NO NO NO NO")
    'new'
    >>> _line_type("NO,NO,NO,NO,NO", delimiter=',')
    'new'
    >>> _line_type("N O N NOON OON O")
    Traceback (most recent call last):
        ...
    ValueError: Unrecognized QDP line...
    >>> _line_type(" some non-comment gibberish")
    Traceback (most recent call last):
        ...
    ValueError: Unrecognized QDP line...
    """
    _decimal_re = r'[+-]?(\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?'
    _command_re = r'READ [TS]ERR(\s+[0-9]+)+'

    sep = delimiter
    if delimiter is None:
        sep = r'\s+'
    _new_re = rf'NO({sep}NO)+'
    _data_re = rf'({_decimal_re}|NO|[-+]?nan)({sep}({_decimal_re}|NO|[-+]?nan))*)'
    _type_re = rf'^\s*((?P<command>{_command_re})|(?P<new>{_new_re})|(?P<data>{_data_re})?\s*(\!(?P<comment>.*))?\s*$'
    _line_type_re = re.compile(_type_re)
    line = line.strip()
    if not line:
        return 'comment'
    match = _line_type_re.match(line)

    if match is None:
        raise ValueError(f'Unrecognized QDP line: {line}')
    for type_, val in match.groupdict().items():
        if val is None:
            continue
        if type_ == 'data':
            return f'data,{len(val.split(sep=delimiter))}'
        else:
            return type_


def _get_type_from_list_of_lines(lines, delimiter=None):
    """Read through the list of QDP file lines and label each line by type

    Parameters
    ----------
    lines : list
        List containing one file line in each entry

    Returns
    -------
    contents : list
        List containing the type for each line (see `line_type_and_data`)
    ncol : int
        The number of columns in the data lines. Must be the same throughout
        the file

    Examples
    --------
    >>> line0 = "! A comment"
    >>> line1 = "543 12 456.0"
    >>> lines = [line0, line1]
    >>> types, ncol = _get_type_from_list_of_lines(lines)
    >>> types[0]
    'comment'
    >>> types[1]
    'data,3'
    >>> ncol
    3
    >>> lines.append("23")
    >>> _get_type_from_list_of_lines(lines)
    Traceback (most recent call last):
        ...
    ValueError: Inconsistent number of columns
    """

    types = [_line_type(line, delimiter=delimiter) for line in lines]
    current_ncol = None
    for type_ in types:
        if type_.startswith('data', ):
            ncol = int(type_[5:])
            if current_ncol is None:
                current_ncol = ncol
            elif ncol != current_ncol:
                raise ValueError('Inconsistent number of columns')

    return types, current_ncol


def _get_lines_from_file(qdp_file):
    if "\n" in qdp_file:
        lines = qdp_file.split("\n")
    elif isinstance(qdp_file, str):
        with open(qdp_file) as fobj:
            lines = [line.strip() for line in fobj.readlines()]
    elif isinstance(qdp_file, Iterable):
        lines = qdp_file
    else:
        raise ValueError('invalid value of qdb_file')

    return lines


def _interpret_err_lines(err_specs, ncols, names=None):
    """Give list of column names from the READ SERR and TERR commands

    Parameters
    ----------
    err_specs : dict
        ``{'serr': [n0, n1, ...], 'terr': [n2, n3, ...]}``
        Error specifications for symmetric and two-sided errors
    ncols : int
        Number of data columns

    Other Parameters
    ----------------
    names : list of str
        Name of data columns (defaults to ['col1', 'col2', ...]), _not_
        including error columns.

    Returns
    -------
    colnames : list
        List containing the column names. Error columns will have the name
        of the main column plus ``_err`` for symmetric errors, and ``_perr``
        and ``_nerr`` for positive and negative errors respectively

    Examples
    --------
    >>> col_in = ['MJD', 'Rate']
    >>> cols = _interpret_err_lines(None, 2, names=col_in)
    >>> cols[0]
    'MJD'
    >>> err_specs = {'terr': [1], 'serr': [2]}
    >>> ncols = 5
    >>> cols = _interpret_err_lines(err_specs, ncols, names=col_in)
    >>> cols[0]
    'MJD'
    >>> cols[2]
    'MJD_nerr'
    >>> cols[4]
    'Rate_err'
    >>> _interpret_err_lines(err_specs, 6, names=col_in)
    Traceback (most recent call last):
        ...
    ValueError: Inconsistent number of input colnames
    """

    colnames = ["" for i in range(ncols)]
    if err_specs is None:
        serr_cols = terr_cols = []

    else:
        # I don't want to empty the original one when using `pop` below
        err_specs = copy.deepcopy(err_specs)

        serr_cols = err_specs.pop("serr", [])
        terr_cols = err_specs.pop("terr", [])

    if names is not None:
        all_error_cols = len(serr_cols) + len(terr_cols) * 2
        if all_error_cols + len(names) != ncols:
            raise ValueError("Inconsistent number of input colnames")

    shift = 0
    for i in range(ncols):
        col_num = i + 1 - shift
        if colnames[i] != "":
            continue

        colname_root = f"col{col_num}"

        if names is not None:
            colname_root = names[col_num - 1]

        colnames[i] = f"{colname_root}"
        if col_num in serr_cols:
            colnames[i + 1] = f"{colname_root}_err"
            shift += 1
            continue

        if col_num in terr_cols:
            colnames[i + 1] = f"{colname_root}_perr"
            colnames[i + 2] = f"{colname_root}_nerr"
            shift += 2
            continue

    assert not np.any([c == "" for c in colnames])

    return colnames


def _get_tables_from_qdp_file(qdp_file, input_colnames=None, delimiter=None):
    """Get all tables from a QDP file

    Parameters
    ----------
    qdp_file : str
        Input QDP file name

    Other Parameters
    ----------------
    input_colnames : list of str
        Name of data columns (defaults to ['col1', 'col2', ...]), _not_
        including error columns.
    delimiter : str
        Delimiter for the values in the table.

    Returns
    -------
    list of `~astropy.table.Table`
        List containing all the tables present inside the QDP file
    """

    lines = _get_lines_from_file(qdp_file)
    contents, ncol = _get_type_from_list_of_lines(lines, delimiter=delimiter)

    table_list = []
    err_specs = {}
    colnames = None

    comment_text = ""
    initial_comments = ""
    command_lines = ""
    current_rows = None

    for line, datatype in zip(lines, contents):
        line = line.strip().lstrip('!')
        # Is this a comment?
        if datatype == "comment":
            comment_text += line + '\n'
            continue

        if datatype == "command":
            # The first time I find commands, I save whatever comments into
            # The initial comments.
            if command_lines == "":
                initial_comments = comment_text
                comment_text = ""

            if err_specs != {}:
                warnings.warn(
                    "This file contains multiple command blocks. Please verify",
                    AstropyUserWarning
                )
            command_lines += line + '\n'
            continue

        if datatype.startswith("data"):
            # The first time I find data, I define err_specs
            if err_specs == {} and command_lines != "":
                for cline in command_lines.strip().split('\n'):
                    command = cline.strip().split()
                    # This should never happen, but just in case.
                    if len(command) < 3:
                        continue
                    err_specs[command[1].lower()] = [int(c) for c in
                                                     command[2:]]
            if colnames is None:
                colnames = _interpret_err_lines(
                    err_specs, ncol, names=input_colnames
                )

            if current_rows is None:
                current_rows = []

            values = []
            for v in line.split(delimiter):
                if v == "NO":
                    values.append(np.ma.masked)
                else:
                    # Understand if number is int or float
                    try:
                        values.append(int(v))
                    except ValueError:
                        values.append(float(v))
            current_rows.append(values)
            continue

        if datatype == "new":
            # Save table to table_list and reset
            if current_rows is not None:
                new_table = Table(names=colnames, rows=current_rows)
                new_table.meta["initial_comments"] = initial_comments.strip().split("\n")
                new_table.meta["comments"] = comment_text.strip().split("\n")
                # Reset comments
                comment_text = ""
                table_list.append(new_table)
                current_rows = None
            continue

    # At the very end, if there is still a table being written, let's save
    # it to the table_list
    if current_rows is not None:
        new_table = Table(names=colnames, rows=current_rows)
        new_table.meta["initial_comments"] = initial_comments.strip().split("\n")
        new_table.meta["comments"] = comment_text.strip().split("\n")
        table_list.append(new_table)

    return table_list


def _understand_err_col(colnames):
    """Get which column names are error columns

    Examples
    --------
    >>> colnames = ['a', 'a_err', 'b', 'b_perr', 'b_nerr']
    >>> serr, terr = _understand_err_col(colnames)
    >>> np.allclose(serr, [1])
    True
    >>> np.allclose(terr, [2])
    True
    >>> serr, terr = _understand_err_col(['a', 'a_nerr'])
    Traceback (most recent call last):
    ...
    ValueError: Missing positive error...
    >>> serr, terr = _understand_err_col(['a', 'a_perr'])
    Traceback (most recent call last):
    ...
    ValueError: Missing negative error...
    """
    shift = 0
    serr = []
    terr = []

    for i, col in enumerate(colnames):
        if col.endswith("_err"):
            # The previous column, but they're numbered from 1!
            # Plus, take shift into account
            serr.append(i - shift)
            shift += 1
        elif col.endswith("_perr"):
            terr.append(i - shift)
            if len(colnames) == i + 1 or not colnames[i + 1].endswith('_nerr'):
                raise ValueError("Missing negative error")
            shift += 2
        elif col.endswith("_nerr") and not colnames[i - 1].endswith('_perr'):
            raise ValueError("Missing positive error")
    return serr, terr


def _read_table_qdp(qdp_file, names=None, table_id=None, delimiter=None):
    """Read a table from a QDP file

    Parameters
    ----------
    qdp_file : str
        Input QDP file name

    Other Parameters
    ----------------
    names : list of str
        Name of data columns (defaults to ['col1', 'col2', ...]), _not_
        including error columns.

    table_id : int, default 0
        Number of the table to be read from the QDP file. This is useful
        when multiple tables present in the file. By default, the first is read.

    delimiter : str
        Any delimiter accepted by the `sep` argument of str.split()

    Returns
    -------
    tables : list of `~astropy.table.Table`
        List containing all the tables present inside the QDP file
    """
    if table_id is None:
        warnings.warn("table_id not specified. Reading the first available "
                      "table", AstropyUserWarning)
        table_id = 0

    tables = _get_tables_from_qdp_file(qdp_file, input_colnames=names, delimiter=delimiter)

    return tables[table_id]


def _write_table_qdp(table, filename=None, err_specs=None):
    """Write a table to a QDP file

    Parameters
    ----------
    table : :class:`~astropy.table.Table`
        Input table to be written
    filename : str
        Output QDP file name

    Other Parameters
    ----------------
    err_specs : dict
        Dictionary of the format {'serr': [1], 'terr': [2, 3]}, specifying
        which columns have symmetric and two-sided errors (see QDP format
        specification)
    """
    import io
    fobj = io.StringIO()

    if 'initial_comments' in table.meta and table.meta['initial_comments'] != []:
        for line in table.meta['initial_comments']:
            line = line.strip()
            if not line.startswith("!"):
                line = "!" + line
            print(line, file=fobj)

    if err_specs is None:
        serr_cols, terr_cols = _understand_err_col(table.colnames)
    else:
        serr_cols = err_specs.pop("serr", [])
        terr_cols = err_specs.pop("terr", [])
    if serr_cols != []:
        col_string = " ".join([str(val) for val in serr_cols])
        print(f"READ SERR {col_string}", file=fobj)
    if terr_cols != []:
        col_string = " ".join([str(val) for val in terr_cols])
        print(f"READ TERR {col_string}", file=fobj)

    if 'comments' in table.meta and table.meta['comments'] != []:
        for line in table.meta['comments']:
            line = line.strip()
            if not line.startswith("!"):
                line = "!" + line
            print(line, file=fobj)

    colnames = table.colnames
    print("!" + " ".join(colnames), file=fobj)
    for row in table:
        values = []
        for val in row:
            if not np.ma.is_masked(val):
                rep = str(val)
            else:
                rep = "NO"
            values.append(rep)
        print(" ".join(values), file=fobj)

    full_string = fobj.getvalue()
    fobj.close()

    if filename is not None:
        with open(filename, 'w') as fobj:
            print(full_string, file=fobj)

    return full_string.split("\n")


class QDPSplitter(core.DefaultSplitter):
    """
    Split on space for QDP tables
    """
    delimiter = ' '


class QDPHeader(basic.CommentedHeaderHeader):
    """
    Header that uses the :class:`astropy.io.ascii.basic.QDPSplitter`
    """
    splitter_class = QDPSplitter
    comment = "!"
    write_comment = "!"


class QDPData(basic.BasicData):
    """
    Data that uses the :class:`astropy.io.ascii.basic.CsvSplitter`
    """
    splitter_class = QDPSplitter
    fill_values = [(core.masked, 'NO')]
    comment = "!"
    write_comment = None


class QDP(basic.Basic):
    """Quick and Dandy Plot table.

    Example::

        ! Initial comment line 1
        ! Initial comment line 2
        READ TERR 1
        READ SERR 3
        ! Table 0 comment
        !a a(pos) a(neg) b be c d
        53000.5   0.25  -0.5   1  1.5  3.5 2
        54000.5   1.25  -1.5   2  2.5  4.5 3
        NO NO NO NO NO
        ! Table 1 comment
        !a a(pos) a(neg) b be c d
        54000.5   2.25  -2.5   NO  3.5  5.5 5
        55000.5   3.25  -3.5   4  4.5  6.5 nan

    The input table above contains some initial comments, the error commands,
    then two tables.
    This file format can contain multiple tables, separated by a line full
    of ``NO``s. Comments are exclamation marks, and missing values are single
    ``NO`` entries. The delimiter is usually whitespace, more rarely a comma.
    The QDP format differentiates between data and error columns. The table
    above has commands::

        READ TERR 1
        READ SERR 3

    which mean that after data column 1 there will be two error columns
    containing its positive and engative error bars, then data column 2 without
    error bars, then column 3, then a column with the symmetric error of column
    3, then the remaining data columns.

    As explained below, table headers are highly inconsistent. Possible
    comments containing column names will be ignored and columns will be called
    ``col1``, ``col2``, etc. unless the user specifies their names with the
    ``names=`` keyword argument,
    When passing column names, pass **only the names of the data columns, not
    the error columns.**
    Error information will be encoded in the names of the table columns.
    (e.g. ``a_perr`` and ``a_nerr`` for the positive and negative error of
    column ``a``, ``b_err`` the symmetric error of column ``b``.)

    When writing tables to this format, users can pass an ``err_specs`` keyword
    passing a dictionary ``{'serr': [3], 'terr': [1, 2]}``, meaning that data
    columns 1 and two will have two additional columns each with their positive
    and negative errors, and data column 3 will have an additional column with
    a symmetric error (just like the ``READ SERR`` and ``READ TERR`` commands
    above)

    Headers are just comments, and tables distributed by various missions
    can differ greatly in their use of conventions. For example, light curves
    distributed by the Swift-Gehrels mission have an extra space in one header
    entry that makes the number of labels inconsistent with the number of cols.
    For this reason, we ignore the comments that might encode the column names
    and leave the name specification to the user.

    Example::

        >               Extra space
        >                   |
        >                   v
        >!     MJD       Err (pos)       Err(neg)        Rate            Error
        >53000.123456   2.378e-05     -2.378472e-05     NO             0.212439

    These readers and writer classes will strive to understand which of the
    comments belong to all the tables, and which ones to each single table.
    General comments will be stored in the ``initial_comments`` meta of each
    table. The comments of each table will be stored in the ``comments`` meta.

    Example::

        t = Table.read(example_qdp, format='ascii.qdp', table_id=1, names=['a', 'b', 'c', 'd'])

    reads the second table (``table_id=1``) in file ``example.qdp`` containing
    the table above. There are four column names but seven data columns, why?
    Because the ``READ SERR`` and ``READ TERR`` commands say that there are
    three error columns.
    ``t.meta['initial_comments']`` will contain the initial two comment lines
    in the file, while ``t.meta['comments']`` will contain ``Table 1 comment``

    The table can be written to another file, preserving the same information,
    as::

        t.write(test_file, err_specs={'terr': [1], 'serr': [3]})

    Note how the ``terr`` and ``serr`` commands are passed to the writer.

    """
    _format_name = 'qdp'
    _io_registry_can_write = True
    _io_registry_suffix = '.qdp'
    _description = 'Quick and Dandy Plotter'

    header_class = QDPHeader
    data_class = QDPData

    def __init__(self, table_id=None, names=None, err_specs=None, sep=None):
        super().__init__()
        self.table_id = table_id
        self.names = names
        self.err_specs = err_specs
        self.delimiter = sep

    def read(self, table):
        self.lines = self.inputter.get_lines(table, newline="\n")
        return _read_table_qdp(self.lines, table_id=self.table_id,
                               names=self.names, delimiter=self.delimiter)

    def write(self, table):
        self._check_multidim_table(table)
        lines = _write_table_qdp(table, err_specs=self.err_specs)
        return lines
