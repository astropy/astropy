# Licensed under a 3-clause BSD style license - see LICENSE.rst

import sys

from ..utils import OrderedDict

__all__ = ['register_reader', 'register_writer', 'register_identifier',
           'identify_format', 'get_reader', 'get_writer', 'read', 'write']

_readers = OrderedDict()
_writers = OrderedDict()
_identifiers = OrderedDict()


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
        `data_format`. This function should take two arguments, which will be
        set to the list of arguments and a dictionary of keyword arguments
        passed to `read` or `write`. The function should return True if the
        input can be identified as being of format `data_format`, and False
        otherwise.
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


def identify_format(origin, data_class_required, fileobj, *args, **kwargs):
    # Loop through identifiers to see which formats match
    valid_formats = []
    for data_format, data_class in _identifiers:
        if data_class is data_class_required:
            # We try first with the fileobj, and failing that with the
            # original arguments
            if (fileobj is not None and
                _identifiers[(data_format, data_class)](
                    origin, [fileobj] + list(args[1:]), kwargs)):
                valid_formats.append((data_format, [fileobj] + list(args[1:])))
            elif _identifiers[(data_format, data_class)](
                    origin, *args, **kwargs):
                valid_formats.append((data_format, args))

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

    fileobj = None
    ctx = None
    try:
        if format is None:
            if len(args) and isinstance(args[0], (bytes, unicode)):
                from ..utils.data import get_readable_fileobj
                try:
                    ctx = get_readable_fileobj(args[0], encoding='binary')
                    fileobj = ctx.__enter__()
                except:
                    ctx = None
                    fileobj = None

            format, args = _get_valid_format('read', cls, fileobj, *args,
                                             **kwargs)

        reader = get_reader(format, cls)
        if fileobj is not None:
            fileobj.seek(0)
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
        format, args = _get_valid_format('write', data.__class__, None, *args,
                                         **kwargs)

    writer = get_writer(format, data.__class__)
    writer(data, *args, **kwargs)


def _get_valid_format(mode, cls, fileobj, *args, **kwargs):
    """
    Returns the first valid format that can be used to read/write the data in
    question.  Mode can be either 'read' or 'write'.
    """

    if mode == 'read':
        funcs = _readers
    elif mode == 'write':
        funcs = _writers

    valid_formats = identify_format(mode, cls, fileobj, *args, **kwargs)

    if len(valid_formats) == 0:
        raise Exception(
            "Format could not be identified. ",
            "Valid formats are {0}".format(
                ', '.join(sorted(r[0] for r in funcs))))
    elif len(valid_formats) > 1:
        raise Exception(
            "Format is ambiguous - options are: {0}".format(
                ', '.join(sorted(x[0] for x in valid_formats))))

    return valid_formats[0]
