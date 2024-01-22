# Licensed under a 3-clause BSD style license - see LICENSE.rst
# This file connects the readers/writers to the astropy.table.Table class

import functools
import os.path

import astropy.io.registry as io_registry
from astropy.table import Table
from astropy.utils.misc import NOT_OVERWRITING_MSG

__all__ = ["PANDAS_FMTS"]

# Astropy users normally expect to not have an index, so default to turn
# off writing the index.  This structure allows for astropy-specific
# customization of all options.
PANDAS_FMTS = {
    "csv": {"read": {}, "write": {"index": False}},
    "fwf": {"read": {}},  # No writer
    "html": {"read": {}, "write": {"index": False}},
    "json": {"read": {}, "write": {}},
}

PANDAS_PREFIX = "pandas."

# Imports for reading HTML
_IMPORTS = False
_HAS_BS4 = False
_HAS_LXML = False
_HAS_HTML5LIB = False


def import_html_libs():
    """Try importing dependencies for reading HTML.

    This is copied from pandas.io.html
    """
    # import things we need
    # but make this done on a first use basis

    global _IMPORTS
    if _IMPORTS:
        return

    global _HAS_BS4, _HAS_LXML, _HAS_HTML5LIB

    from astropy.utils.compat.optional_deps import HAS_BS4 as _HAS_BS4
    from astropy.utils.compat.optional_deps import HAS_HTML5LIB as _HAS_HTML5LIB
    from astropy.utils.compat.optional_deps import HAS_LXML as _HAS_LXML

    _IMPORTS = True


def _pandas_read(fmt, filespec, **kwargs):
    """Provide io Table connector to read table using pandas."""
    try:
        import pandas as pd
    except ImportError:
        raise ImportError("pandas must be installed to use pandas table reader")

    pandas_fmt = fmt[len(PANDAS_PREFIX) :]  # chop the 'pandas.' in front
    read_func = getattr(pd, "read_" + pandas_fmt)

    # Get defaults and then override with user-supplied values
    read_kwargs = PANDAS_FMTS[pandas_fmt]["read"].copy()
    read_kwargs.update(kwargs)

    # Special case: pandas defaults to HTML lxml for reading, but does not attempt
    # to fall back to bs4 + html5lib.  So do that now for convenience if user has
    # not specifically selected a flavor.  If things go wrong the pandas exception
    # with instruction to install a library will come up.
    if pandas_fmt == "html" and "flavor" not in kwargs:
        import_html_libs()
        if not _HAS_LXML and _HAS_HTML5LIB and _HAS_BS4:
            read_kwargs["flavor"] = "bs4"

    df = read_func(filespec, **read_kwargs)

    # Special case for HTML
    if pandas_fmt == "html":
        df = df[0]

    return Table.from_pandas(df)


def _pandas_write(fmt, tbl, filespec, overwrite=False, **kwargs):
    """Provide io Table connector to write table using pandas."""
    pandas_fmt = fmt[len(PANDAS_PREFIX) :]  # chop the 'pandas.' in front

    # Get defaults and then override with user-supplied values
    write_kwargs = PANDAS_FMTS[pandas_fmt]["write"].copy()
    write_kwargs.update(kwargs)

    df = tbl.to_pandas()
    write_method = getattr(df, "to_" + pandas_fmt)

    if not overwrite:
        try:  # filespec is not always a path-like
            exists = os.path.exists(filespec)
        except TypeError:  # skip invalid arguments
            pass
        else:
            if exists:  # only error if file already exists
                raise OSError(NOT_OVERWRITING_MSG.format(filespec))

    return write_method(filespec, **write_kwargs)


for pandas_fmt, defaults in PANDAS_FMTS.items():
    fmt = PANDAS_PREFIX + pandas_fmt  # Full format specifier

    if "read" in defaults:
        func = functools.partial(_pandas_read, fmt)
        io_registry.register_reader(fmt, Table, func)

    if "write" in defaults:
        func = functools.partial(_pandas_write, fmt)
        io_registry.register_writer(fmt, Table, func)
