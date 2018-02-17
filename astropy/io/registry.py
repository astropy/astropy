# Licensed under a 3-clause BSD style license - see LICENSE.rst


import contextlib
import pathlib
import re
import sys

from collections import OrderedDict
from operator import itemgetter

import numpy as np


__all__ = ['register_reader', 'register_writer', 'register_identifier',
           'identify_format', 'get_reader', 'get_writer', 'read', 'write',
           'get_formats', 'IORegistryError', 'delay_doc_updates']


__doctest_skip__ = ['register_identifier']


_readers = OrderedDict()
_writers = OrderedDict()
_identifiers = OrderedDict()

PATH_TYPES = (str, pathlib.Path)


class IORegistryError(Exception):
    """Custom error for registry clashes.
    """
    pass


# If multiple formats are added to one class the update of the docs is quite
# expensive. Classes for which the doc update is temporarly delayed are added
# to this set.
_delayed_docs_classes = set()


@contextlib.contextmanager
def delay_doc_updates(cls):
    """Contextmanager to disable documentation updates when registering
    reader and writer. The documentation is only built once when the
    contextmanager exits.

    .. versionadded:: 1.3

    Parameters
    ----------
    cls : class
        Class for which the documentation updates should be delayed.

    Notes
    -----
    Registering mutliple readers and writers can cause significant overhead
    because the documentation of the corresponding ``read`` and ``write``
    methods are build every time.

    .. warning::
        This contextmanager is experimental and may be replaced by a more
        general approach.

    Examples
    --------
    see for example the source code of ``astropy.table.__init__``.
    """
    _delayed_docs_classes.add(cls)

    yield

    _delayed_docs_classes.discard(cls)
    _update__doc__(cls, 'read')
    _update__doc__(cls, 'write')


def get_formats(data_class=None, readwrite=None):
    """
    Get the list of registered I/O formats as a Table.

    Parameters
    ----------
    data_class : classobj, optional
        Filter readers/writer to match data class (default = all classes).

    readwrite : str or None, optional
        Search only for readers (``"Read"``) or writers (``"Write"``). If None
        search for both.  Default is None.

        .. versionadded:: 1.3

    Returns
    -------
    format_table : Table
        Table of available I/O formats.
    """
    from ..table import Table

    format_classes = sorted(set(_readers) | set(_writers), key=itemgetter(0))
    rows = []

    for format_class in format_classes:
        if (data_class is not None and not _is_best_match(
                data_class, format_class[1], format_classes)):
            continue

        has_read = 'Yes' if format_class in _readers else 'No'
        has_write = 'Yes' if format_class in _writers else 'No'
        has_identify = 'Yes' if format_class in _identifiers else 'No'

        # Check if this is a short name (e.g. 'rdb') which is deprecated in
        # favor of the full 'ascii.rdb'.
        ascii_format_class = ('ascii.' + format_class[0], format_class[1])

        deprecated = 'Yes' if ascii_format_class in format_classes else ''

        rows.append((format_class[1].__name__, format_class[0], has_read,
                     has_write, has_identify, deprecated))

    if readwrite is not None:
        if readwrite == 'Read':
            rows = [row for row in rows if row[2] == 'Yes']
        elif readwrite == 'Write':
            rows = [row for row in rows if row[3] == 'Yes']
        else:
            raise ValueError('unrecognized value for "readwrite": {0}.\n'
                             'Allowed are "Read" and "Write" and None.')

    # Sorting the list of tuples is much faster than sorting it after the table
    # is created. (#5262)
    if rows:
        # Indices represent "Data Class", "Deprecated" and "Format".
        data = list(zip(*sorted(rows, key=itemgetter(0, 5, 1))))
    else:
        data = None
    format_table = Table(data, names=('Data class', 'Format', 'Read', 'Write',
                                      'Auto-identify', 'Deprecated'))

    if not np.any(format_table['Deprecated'] == 'Yes'):
        format_table.remove_column('Deprecated')

    return format_table


def _update__doc__(data_class, readwrite):
    """
    Update the docstring to include all the available readers / writers for the
    ``data_class.read`` or ``data_class.write`` functions (respectively).
    """
    FORMATS_TEXT = 'The available built-in formats are:'

    # Get the existing read or write method and its docstring
    class_readwrite_func = getattr(data_class, readwrite)

    if not isinstance(class_readwrite_func.__doc__, str):
        # No docstring--could just be test code, or possibly code compiled
        # without docstrings
        return

    lines = class_readwrite_func.__doc__.splitlines()

    # Find the location of the existing formats table if it exists
    sep_indices = [ii for ii, line in enumerate(lines) if FORMATS_TEXT in line]
    if sep_indices:
        # Chop off the existing formats table, including the initial blank line
        chop_index = sep_indices[0]
        lines = lines[:chop_index]

    # Find the minimum indent, skipping the first line because it might be odd
    matches = [re.search(r'(\S)', line) for line in lines[1:]]
    left_indent = ' ' * min(match.start() for match in matches if match)

    # Get the available unified I/O formats for this class
    # Include only formats that have a reader, and drop the 'Data class' column
    format_table = get_formats(data_class, readwrite.capitalize())
    format_table.remove_column('Data class')

    # Get the available formats as a table, then munge the output of pformat()
    # a bit and put it into the docstring.
    new_lines = format_table.pformat(max_lines=-1, max_width=80)
    table_rst_sep = re.sub('-', '=', new_lines[1])
    new_lines[1] = table_rst_sep
    new_lines.insert(0, table_rst_sep)
    new_lines.append(table_rst_sep)

    # Check for deprecated names and include a warning at the end.
    if 'Deprecated' in format_table.colnames:
        new_lines.extend(['',
                          'Deprecated format names like ``aastex`` will be '
                          'removed in a future version. Use the full ',
                          'name (e.g. ``ascii.aastex``) instead.'])

    new_lines = [FORMATS_TEXT, ''] + new_lines
    lines.extend([left_indent + line for line in new_lines])

    # Depending on Python version and whether class_readwrite_func is
    # an instancemethod or classmethod, one of the following will work.
    try:
        class_readwrite_func.__doc__ = '\n'.join(lines)
    except AttributeError:
        class_readwrite_func.__func__.__doc__ = '\n'.join(lines)


def register_reader(data_format, data_class, function, force=False):
    """
    Register a reader function.

    Parameters
    ----------
    data_format : str
        The data format identifier. This is the string that will be used to
        specify the data type when reading.
    data_class : classobj
        The class of the object that the reader produces.
    function : function
        The function to read in a data object.
    force : bool, optional
        Whether to override any existing function if already present.
        Default is ``False``.
    """

    if not (data_format, data_class) in _readers or force:
        _readers[(data_format, data_class)] = function
    else:
        raise IORegistryError("Reader for format '{0}' and class '{1}' is "
                              'already defined'
                              ''.format(data_format, data_class.__name__))

    if data_class not in _delayed_docs_classes:
        _update__doc__(data_class, 'read')


def unregister_reader(data_format, data_class):
    """
    Unregister a reader function

    Parameters
    ----------
    data_format : str
        The data format identifier.
    data_class : classobj
        The class of the object that the reader produces.
    """

    if (data_format, data_class) in _readers:
        _readers.pop((data_format, data_class))
    else:
        raise IORegistryError("No reader defined for format '{0}' and class '{1}'"
                              ''.format(data_format, data_class.__name__))

    if data_class not in _delayed_docs_classes:
        _update__doc__(data_class, 'read')


def register_writer(data_format, data_class, function, force=False):
    """
    Register a table writer function.

    Parameters
    ----------
    data_format : str
        The data format identifier. This is the string that will be used to
        specify the data type when writing.
    data_class : classobj
        The class of the object that can be written.
    function : function
        The function to write out a data object.
    force : bool, optional
        Whether to override any existing function if already present.
        Default is ``False``.
    """

    if not (data_format, data_class) in _writers or force:
        _writers[(data_format, data_class)] = function
    else:
        raise IORegistryError("Writer for format '{0}' and class '{1}' is "
                              'already defined'
                              ''.format(data_format, data_class.__name__))

    if data_class not in _delayed_docs_classes:
        _update__doc__(data_class, 'write')


def unregister_writer(data_format, data_class):
    """
    Unregister a writer function

    Parameters
    ----------
    data_format : str
        The data format identifier.
    data_class : classobj
        The class of the object that can be written.
    """

    if (data_format, data_class) in _writers:
        _writers.pop((data_format, data_class))
    else:
        raise IORegistryError("No writer defined for format '{0}' and class '{1}'"
                              ''.format(data_format, data_class.__name__))

    if data_class not in _delayed_docs_classes:
        _update__doc__(data_class, 'write')


def register_identifier(data_format, data_class, identifier, force=False):
    """
    Associate an identifier function with a specific data type.

    Parameters
    ----------
    data_format : str
        The data format identifier. This is the string that is used to
        specify the data type when reading/writing.
    data_class : classobj
        The class of the object that can be written.
    identifier : function
        A function that checks the argument specified to `read` or `write` to
        determine whether the input can be interpreted as a table of type
        ``data_format``. This function should take the following arguments:

           - ``origin``: A string ``"read"`` or ``"write"`` identifying whether
             the file is to be opened for reading or writing.
           - ``path``: The path to the file.
           - ``fileobj``: An open file object to read the file's contents, or
             `None` if the file could not be opened.
           - ``*args``: Positional arguments for the `read` or `write`
             function.
           - ``**kwargs``: Keyword arguments for the `read` or `write`
             function.

        One or both of ``path`` or ``fileobj`` may be `None`.  If they are
        both `None`, the identifier will need to work from ``args[0]``.

        The function should return True if the input can be identified
        as being of format ``data_format``, and False otherwise.
    force : bool, optional
        Whether to override any existing function if already present.
        Default is ``False``.

    Examples
    --------
    To set the identifier based on extensions, for formats that take a
    filename as a first argument, you can do for example::

        >>> def my_identifier(*args, **kwargs):
        ...     return isinstance(args[0], str) and args[0].endswith('.tbl')
        >>> register_identifier('ipac', Table, my_identifier)
    """

    if not (data_format, data_class) in _identifiers or force:
        _identifiers[(data_format, data_class)] = identifier
    else:
        raise IORegistryError("Identifier for format '{0}' and class '{1}' is "
                              'already defined'.format(data_format,
                                                       data_class.__name__))


def unregister_identifier(data_format, data_class):
    """
    Unregister an identifier function

    Parameters
    ----------
    data_format : str
        The data format identifier.
    data_class : classobj
        The class of the object that can be read/written.
    """

    if (data_format, data_class) in _identifiers:
        _identifiers.pop((data_format, data_class))
    else:
        raise IORegistryError("No identifier defined for format '{0}' and class"
                              " '{1}'".format(data_format, data_class.__name__))


def identify_format(origin, data_class_required, path, fileobj, args, kwargs):
    """Loop through identifiers to see which formats match.

    Parameters
    ----------
    origin : str
        A string ``"read`` or ``"write"`` identifying whether the file is to be
        opened for reading or writing.
    data_class_required : object
        The specified class for the result of `read` or the class that is to be
        written.
    path : str, other path object or None
        The path to the file or None.
    fileobj : File object or None.
        An open file object to read the file's contents, or ``None`` if the
        file could not be opened.
    args : sequence
        Positional arguments for the `read` or `write` function. Note that
        these must be provided as sequence.
    kwargs : dict-like
        Keyword arguments for the `read` or `write` function. Note that this
        parameter must be `dict`-like.

    Returns
    -------
    valid_formats : list
        List of matching formats.
    """
    valid_formats = []
    for data_format, data_class in _identifiers:
        if _is_best_match(data_class_required, data_class, _identifiers):
            if _identifiers[(data_format, data_class)](
                    origin, path, fileobj, *args, **kwargs):
                valid_formats.append(data_format)

    return valid_formats


def _get_format_table_str(data_class, readwrite):
    format_table = get_formats(data_class, readwrite=readwrite)
    format_table.remove_column('Data class')
    format_table_str = '\n'.join(format_table.pformat(max_lines=-1))
    return format_table_str


def get_reader(data_format, data_class):
    """Get reader for ``data_format``.

    Parameters
    ----------
    data_format : str
        The data format identifier. This is the string that is used to
        specify the data type when reading/writing.
    data_class : classobj
        The class of the object that can be written.

    Returns
    -------
    reader : callable
        The registered reader function for this format and class.
    """
    readers = [(fmt, cls) for fmt, cls in _readers if fmt == data_format]
    for reader_format, reader_class in readers:
        if _is_best_match(data_class, reader_class, readers):
            return _readers[(reader_format, reader_class)]
    else:
        format_table_str = _get_format_table_str(data_class, 'Read')
        raise IORegistryError(
            "No reader defined for format '{0}' and class '{1}'.\nThe "
            "available formats are:\n{2}".format(
                data_format, data_class.__name__, format_table_str))


def get_writer(data_format, data_class):
    """Get writer for ``data_format``.

    Parameters
    ----------
    data_format : str
        The data format identifier. This is the string that is used to
        specify the data type when reading/writing.
    data_class : classobj
        The class of the object that can be written.

    Returns
    -------
    writer : callable
        The registered writer function for this format and class.
    """
    writers = [(fmt, cls) for fmt, cls in _writers if fmt == data_format]
    for writer_format, writer_class in writers:
        if _is_best_match(data_class, writer_class, writers):
            return _writers[(writer_format, writer_class)]
    else:
        format_table_str = _get_format_table_str(data_class, 'Write')
        raise IORegistryError(
            "No writer defined for format '{0}' and class '{1}'.\nThe "
            "available formats are:\n{2}".format(
                data_format, data_class.__name__, format_table_str))


def read(cls, *args, format=None, **kwargs):
    """
    Read in data.

    The arguments passed to this method depend on the format.
    """

    ctx = None
    try:
        if format is None:
            path = None
            fileobj = None

            if len(args):
                if isinstance(args[0], PATH_TYPES):
                    from ..utils.data import get_readable_fileobj
                    # path might be a pathlib.Path object
                    if isinstance(args[0], pathlib.Path):
                        args = (str(args[0]),) + args[1:]
                    path = args[0]
                    try:
                        ctx = get_readable_fileobj(args[0], encoding='binary')
                        fileobj = ctx.__enter__()
                    except OSError:
                        raise
                    except Exception:
                        fileobj = None
                    else:
                        args = [fileobj] + list(args[1:])
                elif hasattr(args[0], 'read'):
                    path = None
                    fileobj = args[0]

            format = _get_valid_format(
                'read', cls, path, fileobj, args, kwargs)

        reader = get_reader(format, cls)
        data = reader(*args, **kwargs)

        if not isinstance(data, cls):
            # User has read with a subclass where only the parent class is
            # registered.  This returns the parent class, so try coercing
            # to desired subclass.
            try:
                data = cls(data)
            except Exception:
                raise TypeError('could not convert reader output to {0} '
                                'class.'.format(cls.__name__))
    finally:
        if ctx is not None:
            ctx.__exit__(*sys.exc_info())

    return data


def write(data, *args, format=None, **kwargs):
    """
    Write out data.

    The arguments passed to this method depend on the format.
    """

    if format is None:
        path = None
        fileobj = None
        if len(args):
            if isinstance(args[0], PATH_TYPES):
                # path might be a pathlib.Path object
                if isinstance(args[0], pathlib.Path):
                    args = (str(args[0]),) + args[1:]
                path = args[0]
                fileobj = None
            elif hasattr(args[0], 'read'):
                path = None
                fileobj = args[0]

        format = _get_valid_format(
            'write', data.__class__, path, fileobj, args, kwargs)

    writer = get_writer(format, data.__class__)
    writer(data, *args, **kwargs)


def _is_best_match(class1, class2, format_classes):
    """
    Determine if class2 is the "best" match for class1 in the list
    of classes.  It is assumed that (class2 in classes) is True.
    class2 is the the best match if:

    - ``class1`` is a subclass of ``class2`` AND
    - ``class2`` is the nearest ancestor of ``class1`` that is in classes
      (which includes the case that ``class1 is class2``)
    """
    if issubclass(class1, class2):
        classes = {cls for fmt, cls in format_classes}
        for parent in class1.__mro__:
            if parent is class2:  # class2 is closest registered ancestor
                return True
            if parent in classes:  # class2 was superceded
                return False
    return False


def _get_valid_format(mode, cls, path, fileobj, args, kwargs):
    """
    Returns the first valid format that can be used to read/write the data in
    question.  Mode can be either 'read' or 'write'.
    """

    valid_formats = identify_format(mode, cls, path, fileobj, args, kwargs)

    if len(valid_formats) == 0:
        format_table_str = _get_format_table_str(cls, mode.capitalize())
        raise IORegistryError("Format could not be identified.\n"
                              "The available formats are:\n"
                              "{0}".format(format_table_str))
    elif len(valid_formats) > 1:
        raise IORegistryError(
            "Format is ambiguous - options are: {0}".format(
                ', '.join(sorted(valid_formats, key=itemgetter(0)))))

    return valid_formats[0]
