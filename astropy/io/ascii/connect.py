# Licensed under a 3-clause BSD style license - see LICENSE.rst
# This file connects the readers/writers to the astropy.table.Table class


import re
import functools

from astropy.io import registry as io_registry
from astropy.table import Table

__all__ = []


def io_read(format, filename, **kwargs):
    from .ui import read
    if format != 'ascii':
        format = re.sub(r'^ascii\.', '', format)
        kwargs['format'] = format
    return read(filename, **kwargs)


def io_write(format, table, filename, **kwargs):
    from .ui import write
    if format != 'ascii':
        format = re.sub(r'^ascii\.', '', format)
        kwargs['format'] = format
    return write(table, filename, **kwargs)


def io_identify(suffix, origin, filepath, fileobj, *args, **kwargs):
    return filepath is not None and filepath.endswith(suffix)


def _get_connectors_table():
    from .core import FORMAT_CLASSES

    rows = []
    rows.append(('ascii', '', 'Yes', 'ASCII table in any supported format (uses guessing)'))
    for format in sorted(FORMAT_CLASSES):
        cls = FORMAT_CLASSES[format]

        io_format = 'ascii.' + cls._format_name
        description = getattr(cls, '_description', '')
        class_link = ':class:`~{0}.{1}`'.format(cls.__module__, cls.__name__)
        suffix = getattr(cls, '_io_registry_suffix', '')
        can_write = 'Yes' if getattr(cls, '_io_registry_can_write', True) else ''

        rows.append((io_format, suffix, can_write,
                     '{0}: {1}'.format(class_link, description)))
    out = Table(list(zip(*rows)), names=('Format', 'Suffix', 'Write', 'Description'))
    for colname in ('Format', 'Description'):
        width = max(len(x) for x in out[colname])
        out[colname].format = '%-{0}s'.format(width)

    return out
