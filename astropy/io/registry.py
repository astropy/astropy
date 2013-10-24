# Licensed under a 3-clause BSD style license - see LICENSE.rst

import re
import io
import sys

from ..utils import OrderedDict

__all__ = ['register_reader', 'register_writer', 'register_identifier',
           'identify_format', 'get_reader', 'get_writer', 'read', 'write']

__doctest_skip__ = ['register_identifier']

_readers = OrderedDict()
_writers = OrderedDict()
_identifiers = OrderedDict()


def get_formats(data_class=None):
    """
    Get the list of registered I/O formats as a Table.

    Parameters
    ----------
    data_class : classobj
        Filter readers/writer to match data class (default = all classes)

    Returns
    -------
    format_table: Table
        Table of available I/O formats
    """
    from ..table import Table
    format_classes = sorted(set(_readers) | set(_writers))
    rows = []

    for format_class in format_classes:
        if data_class is not None and format_class[1] is not data_class:
            continue

        has_read = 'Yes' if format_class in _readers else 'No'
        has_write = 'Yes' if format_class in _writers else 'No'
        has_identify = 'Yes' if format_class in _identifiers else 'No'

        rows.append((format_class[1].__name__, format_class[0], has_read, has_write, has_identify))

    if not rows:
        raise ValueError('No formats have data class {!r:}'.format(data_class.__name__))

    format_table = Table(zip(*rows),
                         names=('Data class', 'Format', 'Read', 'Write', 'Auto-identify'))
    format_table.sort(['Data class', 'Format'])

    return format_table


def _update__doc__(data_class, readwrite):
    """
    Update the docstring to include all the available readers / writers for the
    ``data_class.read`` or ``data_class.write`` functions (respectively).
    """
    class_readwrite_func = getattr(data_class, readwrite)
    lines = class_readwrite_func.__doc__.splitlines()
    sep_indices = [ii for ii, line in enumerate(lines) if '=======' in line]
    left_indent = lines[sep_indices[0]].index('=')
    lines = lines[:sep_indices[0]]

    format_table = get_formats(data_class)

    # Include only formats that have a reader, and drop the 'Data class' column
    has_read = format_table[readwrite.capitalize()] == 'Yes'
    format_table = format_table[has_read]
    format_table.remove_column('Data class')

    new_lines = format_table.pformat(max_lines=-1)
    table_rst_sep = re.sub('-', '=', new_lines[1])
    new_lines[1] = table_rst_sep
    new_lines.insert(0, table_rst_sep)
    new_lines.append(table_rst_sep)
    lines.extend([' ' * left_indent + line for line in new_lines])

    class_readwrite_func.__func__.__doc__ = '\n'.join(lines)


def register_reader(data_format, data_class, function, force=False):
    """
    Register a reader function.

    Parameters
    ----------
    data_format : string
        The data type identifier. This is the string that will be used to
        specify the data type when reading.
    data_class : classobj
        The class of the object that the reader produces
    function : function
        The function to read in a data object.
    force : bool
        Whether to override any existing function if already present.
    """

    if not (data_format, data_class) in _readers or force:
        _readers[(data_format, data_class)] = function
    else:
        raise Exception('Reader for format {0!r} and class {1!r} is '
                        'already defined'.format(data_format,
                                                 data_class.__name__))

    _update__doc__(data_class, 'read')


def register_writer(data_format, data_class, function, force=False):
    """
    Register a table writer function.

    Parameters
    ----------
    data_format : string
        The data type identifier. This is the string that will be used to
        specify the data type when writing.
    data_class : classobj
        The class of the object that can be written
    function : function
        The function to write out a data object.
    force : bool
        Whether to override any existing function if already present.
    """

    if not (data_format, data_class) in _writers or force:
        _writers[(data_format, data_class)] = function
    else:
        raise Exception('Writer for format {0!r} and class {1!r} is '
                        'already defined'.format(data_format,
                                                 data_class.__name__))

    _update__doc__(data_class, 'write')


def register_identifier(data_format, data_class, identifier, force=False):
    """
    Associate an identifier function with a specific data type.

    Parameters
    ----------
    data_format : str
        The data type identifier. This is the string that is used to
        specify the data type when reading/writing.
    data_class : classobj
        The class of the object that can be written
    identifier : function
        A function that checks the argument specified to `read` or `write` to
        determine whether the input can be interpreted as a table of type
        `data_format`. This function should take the following arguments:

           - `origin`: A string `read` or `write` identifying whether
             the file is to be opened for reading or writing.
           - `path`: The path to the file.
           - `fileobj`: An open file object to read the file's contents, or
             `None` if the file could not be opened.
           - `*args`: A list of positional arguments to the `read` or
             `write` function.
           - `**kwargs`: A list of keyword arguments to the `read` or
             `write` function.

        One or both of `path` or `fileobj` may be `None`.  If they are
        both `None`, the identifier will need to work from `args[0]`.

        The function should return True if the input can be identified
        as being of format `data_format`, and False otherwise.
    force : bool
        Whether to override any existing function if already present.

    Examples
    --------

    To set the identifier based on extensions, for formats that take a
    filename as a first argument, you can do for example::

        >>> def my_identifier(*args, **kwargs):
        ...     return (isinstance(args[0], basestring) and
        ...             args[0].endswith('.tbl'))
        >>> register_identifier('ipac', Table, my_identifier)
    """

    if not (data_format, data_class) in _identifiers or force:
        _identifiers[(data_format, data_class)] = identifier
    else:
        raise Exception('Identifier for format {0!r} and class {1!r} is '
                        'already defined'.format(data_format,
                                                 data_class.__name__))


def identify_format(origin, data_class_required, path, fileobj, args, kwargs):
    # Loop through identifiers to see which formats match
    valid_formats = []
    for data_format, data_class in _identifiers:
        if data_class is data_class_required:
            if _identifiers[(data_format, data_class)](
                origin, path, fileobj, *args, **kwargs):
                valid_formats.append(data_format)

    return valid_formats


def get_reader(data_format, data_class):
    if (data_format, data_class) in _readers:
        return _readers[(data_format, data_class)]
    else:
        raise Exception('No reader defined for format {0!r} and class '
                        '{1!r}'.format(data_format, data_class.__name__))


def get_writer(data_format, data_class):
    if (data_format, data_class) in _writers:
        return _writers[(data_format, data_class)]
    else:
        raise Exception('No writer defined for format {0!r} and class '
                        '{1!r}'.format(data_format, data_class.__name__))


def read(cls, *args, **kwargs):
    """
    Read in data

    The arguments passed to this method depend on the format
    """

    if 'format' in kwargs:
        format = kwargs.pop('format')
    else:
        format = None

    ctx = None
    try:
        if format is None:
            path = None
            fileobj = None

            if len(args):
                if isinstance(args[0], basestring):
                    from ..utils.data import get_readable_fileobj
                    path = args[0]
                    try:
                        ctx = get_readable_fileobj(args[0], encoding='binary')
                        fileobj = ctx.__enter__()
                    except Exception as e:
                        fileobj = None
                    else:
                        args = [fileobj] + list(args[1:])
                elif hasattr(args[0], 'read'):
                    path = None
                    fileobj = args[0]

            format = _get_valid_format(
                'read', cls, path, fileobj, args, kwargs)

        reader = get_reader(format, cls)
        table = reader(*args, **kwargs)

        if not isinstance(table, cls):
            raise TypeError(
                "reader should return a {0} instance".format(cls.__name__))
    finally:
        if ctx is not None:
            ctx.__exit__(*sys.exc_info())

    return table


def write(data, *args, **kwargs):
    """
    Write out data

    The arguments passed to this method depend on the format
    """

    if 'format' in kwargs:
        format = kwargs.pop('format')
    else:
        format = None

    if format is None:
        path = None
        fileobj = None
        if len(args):
            if isinstance(args[0], basestring):
                path = args[0]
                fileobj = None
            elif hasattr(args[0], 'read'):
                path = None
                fileobj = args[0]

        format = _get_valid_format(
            'write', data.__class__, path, fileobj, args, kwargs)

    writer = get_writer(format, data.__class__)
    writer(data, *args, **kwargs)


def _get_valid_format(mode, cls, path, fileobj, args, kwargs):
    """
    Returns the first valid format that can be used to read/write the data in
    question.  Mode can be either 'read' or 'write'.
    """

    if mode == 'read':
        funcs = _readers
    elif mode == 'write':
        funcs = _writers

    valid_formats = identify_format(mode, cls, path, fileobj, args, kwargs)

    if len(valid_formats) == 0:
        raise Exception(
            "Format could not be identified. ",
            "Valid formats are {0}".format(
                ', '.join(sorted(x[0] for x in funcs))))
    elif len(valid_formats) > 1:
        raise Exception(
            "Format is ambiguous - options are: {0}".format(
                ', '.join(sorted(valid_formats))))

    return valid_formats[0]
