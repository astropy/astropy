# Licensed under a 3-clause BSD style license - see LICENSE.rst
# This file connects the readers/writers to the astropy.table.Table class

import functools

from astropy.table import Table
import astropy.io.registry as io_registry

__all__ = ['PANDAS_FMTS']

# By default pandas ``to_<fmt>`` methods want to write the index row, which
# for ``Table.to_pandas()`` is going to be just the integer row number.
# Normally astropy users don't expect that in their output table.
# But JSON is different for some reason.
PANDAS_FMTS = {'csv': {'read': {},
                       'write': {'index': False}},
               'fwf': {'read': {}},  # No writer
               'html': {'read': {'flavor': 'bs4'},
                        'write': {'index': False}},
               'json': {'read': {'orient': 'columns'},
                        'write': {'orient': 'columns', 'index': True}}}
PANDAS_PREFIX = 'pandas.'


def _pandas_read(fmt, filespec, **kwargs):
    """Provide io Table connector to read table using pandas.

    """
    try:
        import pandas
    except ImportError:
        raise ImportError('pandas must be installed to use pandas table reader')

    pandas_fmt = fmt[len(PANDAS_PREFIX):]  # chop the 'pandas.' in front
    read_func = getattr(pandas, 'read_' + pandas_fmt)

    # Get defaults and then override with user-supplied values
    read_kwargs = PANDAS_FMTS[pandas_fmt]['read'].copy()
    read_kwargs.update(kwargs)

    df = read_func(filespec, **read_kwargs)

    # Special case for HTML
    if pandas_fmt == 'html':
        df = df[0]

    return Table.from_pandas(df)


def _pandas_write(fmt, tbl, filespec, **kwargs):
    """Provide io Table connector to write table using pandas.

    """
    pandas_fmt = fmt[len(PANDAS_PREFIX):]  # chop the 'pandas.' in front

    # Get defaults and then override with user-supplied values
    write_kwargs = PANDAS_FMTS[pandas_fmt]['write'].copy()
    write_kwargs.update(kwargs)

    df = tbl.to_pandas()
    write_method = getattr(df, 'to_' + pandas_fmt)

    return write_method(filespec, **write_kwargs)


for pandas_fmt, defaults in PANDAS_FMTS.items():
    fmt = PANDAS_PREFIX + pandas_fmt  # Full format specifier

    if 'read' in defaults:
        func = functools.partial(_pandas_read, fmt)
        io_registry.register_reader(fmt, Table, func)

    if 'write' in defaults:
        func = functools.partial(_pandas_write, fmt)
        io_registry.register_writer(fmt, Table, func)
