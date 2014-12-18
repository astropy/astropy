# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import os
import re
import sys
import json

import numpy as np

from ...utils import OrderedDict
from ...extern import six
from ...extern.six.moves import zip

__all__ = ['identify_format', 'get_reader', 'get_writer', 'read', 'write',
           'get_formats', 'BaseIO']

__doctest_skip__ = ['register_identifier']


_registry = OrderedDict()

def _load_builtins():
    global _registry
    with open(os.path.join(os.path.abspath(os.path.dirname(__file__)), "registry.json")) as f_reg:
        _registry.update(json.load(f_reg))

# Read in registry of built-in formats
_load_builtins()



def _get_object_by_string(path_to_function):
    """
    Given a path, e.g. 'astropy.table.table.Table', return the actual object at
    that path.
    """
    module, function = path_to_function.rsplit('.', 1)
    if sys.version_info[:2] < (2, 7):
        __import__(module)
        return getattr(sys.modules[module], function)
    else:
        from importlib import import_module
        return getattr(import_module(module), function)


def _get_object_name(object_class):
    """
    Given either a class or a string giving the path to a class, return the
    name of the class.
    """
    if isinstance(object_class, six.string_types):
        return object_class.split('.')[-1]
    else:
        return object_class.__name__


def get_registry_subset(data_class, identify=None, read=None, write=None,
                        data_format=None, resolve=False):
    """
    Return a subset of the registry.

    Parameters
    ----------
    data_class : class
        The class to return a subset of the registry for. If this class matches
        both a class in the registry as well as another class that it is a
        subclass of, the formats that match the class directly take precedence
        over the parent classes.
    identify : bool
        If set, only include the subset of the registry that have the identify
        flag set to this value.
    read : bool
        If set, only include the subset of the registry that have the read
        flag set to this value.
    write : bool
        If set, only include the subset of the registry that have the write
        flag set to this value.
    data_format : str
        If set, only include the subset of the registry for this format.
    resolve : bool
        If set, resolve the module and io_class in the registry.
    """

    # We now need to resolve the data classes in the repository since we are
    # doing subclass checks. This only imports e.g. Table, not all the IO
    # classes.
    resolve_data_classes()

    # We start off by selecting the subset of the registry that matches the
    # data class exactly.
    if data_class in _registry:
        registry_subset = _registry[data_class]
    else:
        registry_subset = {}

    for dc in _registry:

        if dc is data_class:
            continue  # we already did this

        if not issubclass(data_class, dc):
            continue

        for fmt_name in _registry[dc]:
            if fmt_name not in registry_subset:
                registry_subset[fmt_name] = _registry[dc][fmt_name]

    # Now remove non-matching formats
    filtered_registry_subset = {}
    for fmt_name in registry_subset:

        if data_format is not None and fmt_name != data_format:
            continue

        fmt = registry_subset[fmt_name]

        if identify is not None and identify is not fmt['identify']:
            continue

        if read is not None and read is not fmt['read']:
            continue

        if write is not None and write is not fmt['write']:
            continue

        filtered_registry_subset[fmt_name] = fmt

    # Finally we resolve IO classes if requested
    resolve_io_classes(filtered_registry_subset, subset=True)

    return filtered_registry_subset


def resolve_data_classes(registry_class=None):
    """
    Resolve the primary keys of the registry.
    """

    if registry_class is None:
        registry_class = _registry

    for data_class in registry_class:
        if isinstance(data_class, six.string_types):
            registry_class[_get_object_by_string(data_class)] = registry_class.pop(data_class)


def resolve_io_classes(registry_class=None, subset=False):
    """
    Resolve the I/O classes in the registry

    Parameters
    ----------
    registry_class : dict, optional
        If specified, only resolve I/O classes in this registry class
    subset : bool, optional
        If specified, the registry class provided is a subset with no primary
        data class key.
    """

    if registry_class is None:
        if subset:
            raise ValueError("Cannot use subset mode without specifying a registry class")
        registry_class = _registry

    if subset:
        for fmt_name in registry_class:
            fmt = registry_class[fmt_name]
            if isinstance(fmt['io_class'], six.string_types):
                fmt['io_class'] = _get_object_by_string(fmt['module'] + '.' + fmt['io_class'])
    else:
        for data_class in registry_class:
            resolve_io_classes(registry_class[data_class], subset=True)


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

    from ...table import Table

    rows = []

    if data_class is None:
        formats = _registry
    else:
        formats = {data_class: get_registry_subset(data_class)}

    for dc in formats:

        for fmt in formats[dc]:

            has_read = 'Yes' if formats[dc][fmt]['read'] else 'No'
            has_write = 'Yes' if formats[dc][fmt]['write'] else 'No'
            has_identify = 'Yes' if formats[dc][fmt]['identify'] else 'No'

            # TODO: add back support for deprecated ASCII formats
            deprecated = ''

            rows.append((_get_object_name(dc),
                         fmt, has_read, has_write, has_identify, deprecated))

    data = list(zip(*rows)) if rows else None
    format_table = Table(data, names=('Data class', 'Format', 'Read', 'Write',
                                      'Auto-identify', 'Deprecated'))
    format_table.sort(['Data class', 'Deprecated', 'Format'])

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
    lines = class_readwrite_func.__doc__.splitlines()

    # Find the location of the existing formats table if it exists
    sep_indices = [ii for ii, line in enumerate(lines) if FORMATS_TEXT in line]
    if sep_indices:
        # Chop off the existing formats table, including the initial blank line.
        chop_index = sep_indices[0]
        lines = lines[:chop_index]

    # Find the minimum indent, skipping the first line because it might be odd
    matches = [re.search('(\S)', line) for line in lines[1:]]
    left_indent = min(match.start() for match in matches if match)

    # Get the available unified I/O formats for this class
    format_table = get_formats(data_class)

    # Include only formats that have a reader, and drop the 'Data class' column
    has_readwrite = format_table[readwrite.capitalize()] == 'Yes'
    format_table = format_table[has_readwrite]
    format_table.remove_column('Data class')

    # Get the available formats as a table, then munge the output of pformat() a bit and
    # put it into the docstring.
    new_lines = format_table.pformat(max_lines=-1, max_width=80)
    table_rst_sep = re.sub('-', '=', new_lines[1])
    new_lines[1] = table_rst_sep
    new_lines.insert(0, table_rst_sep)
    new_lines.append(table_rst_sep)

    # Check for deprecated names and include a warning at the end.
    if 'Deprecated' in format_table.colnames:
        new_lines.extend(['',
                          'Deprecated format names like ``aastex`` will be removed in a '
                          'future version.',
                          'Use the full name (e.g. ``ascii.aastex``) instead.'])

    new_lines = [FORMATS_TEXT, ''] + new_lines
    lines.extend([' ' * left_indent + line for line in new_lines])

    # Depending on Python version and whether class_readwrite_func is
    # an instancemethod or classmethod, one of the following will work.
    try:
        class_readwrite_func.__doc__ = '\n'.join(lines)
    except AttributeError:
        class_readwrite_func.__func__.__doc__ = '\n'.join(lines)


def identify_format(origin, data_class_required, path, fileobj, args, kwargs):

    # Initialize list of valid formats
    valid_formats = []

    # Get only the entries in the registry that match the data class required.
    # This returns a dictionary with unique entries for each format.
    registry_subset = get_registry_subset(data_class_required, identify=True, resolve=True)

    # We now loop over the supported formats and try all matching identifiers
    for fmt_name in registry_subset:
        fmt = registry_subset[fmt_name]
        if fmt["io_class"] is None:
            identifier = fmt["identifier"]
        else:
            identifier = fmt["io_class"]().identify
        if identifier(origin, path, fileobj, *args, **kwargs):
            valid_formats.append(fmt_name)

    return valid_formats


def _get_format_table_str(data_class, readwrite):
    format_table = get_formats(data_class)
    if len(format_table) > 0:
        has_readwrite = format_table[readwrite] == 'Yes'
        format_table = format_table[has_readwrite]
    format_table.remove_column('Data class')
    format_table_str = '\n'.join(format_table.pformat(max_lines=-1))
    return format_table_str


def get_reader(data_format, data_class):

    # Get only the entries in the registry that match the data class required.
    # This returns a dictionary with unique entries for each format.
    registry_subset = get_registry_subset(data_class, read=True,
                                          data_format=data_format, resolve=True)

    if len(registry_subset) == 0:

        format_table_str = _get_format_table_str(data_class, 'Read')
        raise Exception("No reader defined for format '{0}' and class '{1}'.\n"
                        'The available formats are:\n'
                        '{2}'
                        .format(data_format, data_class.__name__, format_table_str))

    else:

        fmt = registry_subset[data_format]

        if fmt['io_class'] is None:  # backward compatibility
            return fmt['reader']
        else:
            return fmt['io_class']().read


def get_writer(data_format, data_class):

    # Get only the entries in the registry that match the data class required.
    # This returns a dictionary with unique entries for each format.
    registry_subset = get_registry_subset(data_class, write=True,
                                          data_format=data_format, resolve=True)

    if len(registry_subset) == 0:

        format_table_str = _get_format_table_str(data_class, 'Read')
        raise Exception("No writer defined for format '{0}' and class '{1}'.\n"
                        'The available formats are:\n'
                        '{2}'
                        .format(data_format, data_class.__name__, format_table_str))

    else:

        fmt = registry_subset[data_format]

        if fmt['io_class'] is None:  # backward compatibility
            return fmt['writer']
        else:
            return fmt['io_class']().write


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
                if isinstance(args[0], six.string_types):
                    from ...utils.data import get_readable_fileobj
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
        data = reader(*args, **kwargs)

        if not isinstance(data, cls):
            if issubclass(cls, data.__class__):
                # User has read with a subclass where only the parent class is
                # registered.  This returns the parent class, so try coercing to
                # desired subclass.
                try:
                    data = cls(data)
                except:
                    raise TypeError('could not convert reader output to {0} class'
                                    .format(cls.__name__))
            else:
                raise TypeError("reader should return a {0} instance".format(cls.__name__))
    finally:
        if ctx is not None:
            ctx.__exit__(*sys.exc_info())

    return data


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
            if isinstance(args[0], six.string_types):
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

    valid_formats = identify_format(mode, cls, path, fileobj, args, kwargs)

    if len(valid_formats) == 0:
        format_table_str = _get_format_table_str(cls, mode.capitalize())
        raise Exception("Format could not be identified.\n"
                        "The available formats are:\n"
                        "{0}".format(format_table_str))
    elif len(valid_formats) > 1:
        raise Exception(
            "Format is ambiguous - options are: {0}".format(
                ', '.join(sorted(valid_formats))))

    return valid_formats[0]


class MetaRegisterBaseIO(type):

    def __init__(cls, name, bases, members):

        super(MetaRegisterBaseIO, cls).__init__(name, bases, members)

        format_abbreviation = members.get('_format_name')
        if format_abbreviation is None:
            if cls.__name__ == 'BaseIO':
                return
            else:
                raise ValueError("_format_name is not defined")

        supported_class = members.get('_supported_class')
        if supported_class is None:
            raise ValueError("_supported_class is not defined")

        import inspect
        frm = inspect.stack()[1]
        module = inspect.getmodule(frm[0])

        resolve_data_classes()

        if not supported_class in _registry:
            _registry[supported_class] = {}

        _registry[supported_class][format_abbreviation] = {
            'module': module,
            'io_class': cls,
            'read': members.get('read') is not None,
            'write':  members.get('write') is not None,
            'identify':  members.get('identify') is not None
        }


@six.add_metaclass(MetaRegisterBaseIO)
class BaseIO(object):

    _format_name = None
    _supported_class = None
    _flexible = False


# Backward-compatibility functions
# TODO: add deprecation warnings


def _get_flexible_io_dict(data_class, data_format):

    resolve_data_classes()

    if data_class not in _registry:
        _registry[data_class] = {}

    if data_format in _registry[data_class]:
        if _registry[data_class][data_format]["io_class"] is not None:
            raise Exception("I/O class already defined for format "
                            "'{0}' and class '{1}'".format(data_format, data_class.__name__))
        io_dict = _registry[data_class][data_format]
    else:
        io_dict = {}
        io_dict['module'] = None
        io_dict['io_class'] = None
        io_dict['read'] = False
        io_dict['reader'] = None
        io_dict['write'] = False
        io_dict['writer'] = None
        io_dict['identify'] = False
        io_dict['identifier'] = None
        _registry[data_class][data_format] = io_dict

    return io_dict


def register_reader(data_format, data_class, function, force=False):
    """
    Register a reader function (deprecated).
    """

    io_dict = _get_flexible_io_dict(data_class, data_format)

    if io_dict['reader'] is None or force:
        io_dict['read'] = True
        io_dict['reader'] = function
    else:
        raise Exception("Reader for format '{0}' and class '{1}' is "
                        "already defined".format(data_format, data_class.__name__))


def register_writer(data_format, data_class, function, force=False):
    """
    Register a writer function (deprecated).
    """

    io_dict = _get_flexible_io_dict(data_class, data_format)

    if io_dict['writer'] is None or force:
        io_dict['write'] = True
        io_dict['writer'] = function
    else:
        raise Exception("Writer for format '{0}' and class '{1}' is "
                        "already defined".format(data_format, data_class.__name__))

def register_identifier(data_format, data_class, function, force=False):
    """
    Register an identifier function (deprecated).
    """

    io_dict = _get_flexible_io_dict(data_class, data_format)

    if io_dict['identifier'] is None or force:
        io_dict['identify'] = True
        io_dict['identifier'] = function
    else:
        raise Exception("Identifier for format '{0}' and class '{1}' is "
                        "already defined".format(data_format, data_class.__name__))
