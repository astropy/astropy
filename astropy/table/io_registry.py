# Licensed under a 3-clause BSD style license - see LICENSE.rst

__all__ = ['register_reader', 'register_writer', 'register_identifier',
           'identify_format', 'get_reader', 'get_writer']

_readers = {}
_writers = {}
_identifiers = {}


def register_reader(table_format, function, force=False):
    '''
    Register a table reader function.

    Parameters
    ----------
    table_format : string
        The table type identifier. This is the string that will be used to
        specify the table type when reading.
    function : function
        The function to read in a table.
    force : bool
        Whether to override any existing function if already present.
    '''

    if not table_format in _readers or force:
        _readers[table_format] = function
    else:
        raise Exception("Reader for format {0:s} is already defined".format(table_format))


def register_writer(table_format, function, force=False):
    '''
    Register a table writer function.

    Parameters
    ----------
    table_format : string
        The table type identifier. This is the string that will be used to
        specify the table type when writing.
    function : function
        The function to write in a table.
    force : bool
        Whether to override any existing function if already present.
    '''

    if not table_format in _writers or force:
        _writers[table_format] = function
    else:
        raise Exception("Writer for format {0:s} is already defined".format(table_format))


def register_identifier(table_format, identifier, force=False):
    '''
    Associate an identifier function with a specific table type.

    Parameters
    ----------
    table_format : str
        The table type identifier. This is the string that is used to
        specify the table type when reading/writing.
    identifier : function
        A function that checks the argument specified to `Table.read` or
        `Table.write` to determine whether the input can be interpreted as a
        table of type `table_format`. This function should take two arguments,
        which will be set to the list of arguments and a dictionary of keyword
        arguments passed to `Table.read` or `Table.write`. The function should
        return True if the input can be identified as being of format
        `table_format`, and False otherwise.
    force : bool
        Whether to override any existing function if already present.

    Examples
    --------

    To set the identifier based on extensions, for formats that take a
    filename as a first argument, you can do for example::

    >>> register_identifier('ipac', lambda args, kwargs: isinstance(args[0], basestring) and args[0].endswith('.tbl'))
    '''

    if not table_format in _identifiers or force:
        _identifiers[table_format] = identifier
    else:
        raise Exception("Identifier for format %s is already defined" % table_format)


def identify_format(origin, args, kwargs):
    # Loop through identifiers to see which formats match
    valid_formats = []
    for table_format in _identifiers:
        if _identifiers[table_format](origin, args, kwargs):
            valid_formats.append(table_format)

    return valid_formats


def get_reader(table_format):
    if table_format in _readers:
        return _readers[table_format]
    else:
        raise Exception("No reader defined for format '{0}'".format(table_format))


def get_writer(table_format):
    if table_format in _writers:
        return _writers[table_format]
    else:
        raise Exception("No writer defined for format '{0}'".format(table_format))
