# Licensed under a 3-clause BSD style license - see LICENSE.rst

from ..utils import OrderedDict

__all__ = ['register_reader', 'register_writer', 'register_identifier',
           'identify_format', 'get_reader', 'get_writer']

_readers = OrderedDict()
_writers = OrderedDict()
_identifiers = OrderedDict()


def register_reader(data_format, data_class, function, force=False):
    '''
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
    '''

    if not issubclass(data_class, DataIO):
        raise Exception("data_class should be a sub-class of DataIO")

    if not (data_format, data_class) in _readers or force:
        _readers[(data_format, data_class)] = function
    else:
        raise Exception("Reader for format '{0:s}' and class '{1:s}' is already defined".format(data_format, data_class.__name__))


def register_writer(data_format, data_class, function, force=False):
    '''
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
    '''

    if not issubclass(data_class, DataIO):
        raise Exception("data_class should be a sub-class of DataIO")

    if not (data_format, data_class) in _writers or force:
        _writers[(data_format, data_class)] = function
    else:
        raise Exception("Writer for format '{0:s}' and class '{1:s}' is already defined".format(data_format, data_class.__name__))


def register_identifier(data_format, data_class, identifier, force=False):
    '''
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

    >>> register_identifier('ipac', lambda args, kwargs: isinstance(args[0], basestring) and args[0].endswith('.tbl'))
    '''

    if not issubclass(data_class, DataIO):
        raise Exception("data_class should be a sub-class of DataIO")

    if not (data_format, data_class) in _identifiers or force:
        _identifiers[(data_format, data_class)] = identifier
    else:
        raise Exception("Identifier for format '{0:s}' and class '{1:s}' is already defined".format(data_format, data_class.__name__))


def identify_format(origin, data_class_required, args, kwargs):
    # Loop through identifiers to see which formats match
    valid_formats = []
    for data_format, data_class in _identifiers:
        if data_class is data_class_required:
            if _identifiers[(data_format, data_class)](origin, args, kwargs):
                valid_formats.append(data_format)

    return valid_formats


def get_reader(data_format, data_class):
    if (data_format, data_class) in _readers:
        return _readers[(data_format, data_class)]
    else:
        raise Exception("No reader defined for format '{0}' and class '{1}'".format(data_format, data_class.__name__))


def get_writer(data_format, data_class):
    if (data_format, data_class) in _writers:
        return _writers[(data_format, data_class)]
    else:
        raise Exception("No writer defined for format '{0}' and class '{1}'".format(data_format, data_class.__name__))


class DataIO(object):

    @classmethod
    def read(cls, *args, **kwargs):
        '''
        Read in data

        The arguments passed to this method depend on the format
        '''

        if 'format' in kwargs:
            format = kwargs.pop('format')
        else:
            format = None

        if format is None:

            valid_formats = identify_format('read', cls, args, kwargs)

            if len(valid_formats) == 0:
                raise Exception("Format could not be identified")
            elif len(valid_formats) > 1:
                raise Exception("Format is ambiguous - options are: {0:s}".format(', '.join(valid_formats)))
            else:
                format = valid_formats[0]

        reader = get_reader(format, cls)
        table = reader(*args, **kwargs)
        if not isinstance(table, cls):
            raise TypeError("reader should return a {0:s} instance".format(cls.__name__))
        return table

    def write(self, *args, **kwargs):
        '''
        Write out data

        The arguments passed to this method depend on the format
        '''

        if 'format' in kwargs:
            format = kwargs.pop('format')
        else:
            format = None

        if format is None:

            valid_formats = identify_format('write', self.__class__, args, kwargs)

            if len(valid_formats) == 0:
                raise Exception("Format could not be identified")
            elif len(valid_formats) > 1:
                raise Exception("Format is ambiguous - options are: {0:s}".format(', '.join(valid_formats)))
            else:
                format = valid_formats[0]

        writer = get_writer(format, self.__class__)
        writer(self, *args, **kwargs)
