# Licensed under a 3-clause BSD style license - see LICENSE.rst
# This file connects the readers/writers to the astropy.table.Table class

import re

from .. import registry as io_registry
from ...table import Table

__all__ = []


# Generic
# =======


def read_asciitable(filename, **kwargs):
    from .ui import read
    return read(filename, **kwargs)

io_registry.register_reader('ascii', Table, read_asciitable)


def write_asciitable(table, filename, **kwargs):
    from .ui import write
    return write(table, filename, **kwargs)

io_registry.register_writer('ascii', Table, write_asciitable)


def io_read(format, filename, **kwargs):
    from .ui import read
    if 'guess' not in kwargs:
        kwargs['guess'] = False
    format = re.sub(r'^ascii\.', '', format)
    return read(filename, format=format, **kwargs)


def io_write(format, table, filename, **kwargs):
    from .ui import write
    format = re.sub(r'^ascii\.', '', format)
    return write(table, filename, format=format, **kwargs)


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
        class_link = ':class:`~{}.{}`'.format(cls.__module__, cls.__name__)
        suffix = getattr(cls, '_io_registry_suffix', '')
        can_write = 'Yes' if getattr(cls, '_io_registry_can_write', True) else ''

        rows.append((io_format, suffix, can_write,
                     '{}: {}'.format(class_link, description)))
    out = Table(zip(*rows), names=('Format', 'Suffix', 'Write', 'Description'))
    for colname in ('Format', 'Description'):
        width = max(len(x) for x in out[colname])
        out[colname].format = '%-{}s'.format(width)

    return out
