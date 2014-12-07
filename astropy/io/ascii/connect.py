# Licensed under a 3-clause BSD style license - see LICENSE.rst
# This file connects the readers/writers to the astropy.table.Table class

from __future__ import absolute_import, division, print_function

import re

from ..registry import BaseIO
from ...table import Table

__all__ = []


# Generic
# =======


def read_asciitable(filename, **kwargs):
    from .ui import read
    return read(filename, **kwargs)


def write_asciitable(table, filename, **kwargs):
    from .ui import write
    return write(table, filename, **kwargs)


# Generic ASCII I/O class

class ASCIITableIO(BaseIO):

    _format_name = 'ascii'
    _supported_class = Table

    @staticmethod
    def read(filename, **kwargs):
        from .ui import read
        return read(filename, **kwargs)

    @staticmethod
    def write(table, filename, **kwargs):
        from .ui import write
        return write(table, filename, **kwargs)

# Format-specific classes

BUILTIN_ASCII_FORMATS = ['ascii.aastex',
                         'ascii.basic',
                         'ascii.cds',
                         'ascii.commented_header',
                         'ascii.csv',
                         'ascii.daophot',
                         'ascii.fast_basic',
                         'ascii.fast_commented_header',
                         'ascii.fast_csv',
                         'ascii.fast_no_header',
                         'ascii.fast_rdb',
                         'ascii.fast_tab',
                         'ascii.fixed_width',
                         'ascii.fixed_width_no_header',
                         'ascii.fixed_width_two_line',
                         'ascii.html',
                         'ascii.ipac',
                         'ascii.latex',
                         'ascii.no_header',
                         'ascii.rdb',
                         'ascii.sextractor',
                         'ascii.tab']

EXTENSIONS = {
    'ascii.csv': ('csv'),
    'ascii.rdb': ('rdb'),
    'ascii.html': ('html'),
    'ascii.latex': ('tex')
}

ASCII_CLASSES = {}

for ascii_format in BUILTIN_ASCII_FORMATS:

    _format_name = 'ascii'
    _supported_class = Table

    def read(self, filename, **kwargs):
        from .ui import read
        fmt = re.sub(r'^ascii\.', '', self._format_name)
        return read(filename, format=fmt, **kwargs)

    def write(self, table, filename, **kwargs):
        from .ui import write
        fmt = re.sub(r'^ascii\.', '', self._format_name)
        return write(table, filename, format=fmt, **kwargs)

    def identify(self, origin, filepath, fileobj, *args, **kwargs):
        return filepath is not None and filepath.endswith(self._extensions)

    name = "".join([x.capitalize() for x in ascii_format.split(".")]) + "TableIO"

    attributes = {
        "_format_name": ascii_format,
        "_supported_class": Table,
        "read": read,
        "write": write
    }

    if ascii_format in EXTENSIONS:
        attributes["_extensions"] = EXTENSIONS[ascii_format],
        attributes["identify"] = identify

    globals()[name] = type(name, (BaseIO,), attributes)

del read, write, identify, name

# TODO: add back deprecated formats?