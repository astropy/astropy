"""
Table property for providing information about table.
"""
import os

# Licensed under a 3-clause BSD style license - see LICENSE.rst
import sys
from contextlib import contextmanager
from inspect import isclass

import numpy as np

from astropy.utils.data_info import DataInfo

__all__ = ["table_info", "TableInfo", "serialize_method_as"]


def table_info(tbl, option="attributes", out=""):
    """
    Write summary information about column to the ``out`` filehandle.
    By default this prints to standard output via sys.stdout.

    The ``option`` argument specifies what type of information
    to include.  This can be a string, a function, or a list of
    strings or functions.  Built-in options are:

    - ``attributes``: basic column meta data like ``dtype`` or ``format``
    - ``stats``: basic statistics: minimum, mean, and maximum

    If a function is specified then that function will be called with the
    column as its single argument.  The function must return an OrderedDict
    containing the information attributes.

    If a list is provided then the information attributes will be
    appended for each of the options, in order.

    Examples
    --------
    >>> from astropy.table.table_helpers import simple_table
    >>> t = simple_table(size=2, kinds='if')
    >>> t['a'].unit = 'm'
    >>> t.info()
    <Table length=2>
    name  dtype  unit
    ---- ------- ----
       a   int64    m
       b float64

    >>> t.info('stats')
    <Table length=2>
    name mean std min max
    ---- ---- --- --- ---
       a  1.5 0.5   1   2
       b  1.5 0.5   1   2

    Parameters
    ----------
    option : str, callable, list of (str or callable)
        Info option, defaults to 'attributes'.
    out : file-like, None
        Output destination, default is sys.stdout.  If None then a
        Table with information attributes is returned

    Returns
    -------
    info : `~astropy.table.Table` if out==None else None
    """
    from .table import Table

    if out == "":
        out = sys.stdout

    descr_vals = [tbl.__class__.__name__]
    if tbl.masked:
        descr_vals.append("masked=True")
    descr_vals.append(f"length={len(tbl)}")

    outlines = ["<" + " ".join(descr_vals) + ">"]

    cols = list(tbl.columns.values())
    if tbl.colnames:
        infos = []
        for col in cols:
            infos.append(col.info(option, out=None))

        info = Table(infos, names=list(infos[0]))
    else:
        info = Table()

    if out is None:
        return info

    # Since info is going to a filehandle for viewing then remove uninteresting
    # columns.
    if "class" in info.colnames:
        # Remove 'class' info column if all table columns are the same class
        # and they are the default column class for that table.
        uniq_types = {type(col) for col in cols}
        if len(uniq_types) == 1 and isinstance(cols[0], tbl.ColumnClass):
            del info["class"]

    if "n_bad" in info.colnames and np.all(info["n_bad"] == 0):
        del info["n_bad"]

    # Standard attributes has 'length' but this is typically redundant
    if "length" in info.colnames and np.all(info["length"] == len(tbl)):
        del info["length"]

    for name in info.colnames:
        if info[name].dtype.kind in "SU" and np.all(info[name] == ""):
            del info[name]

    if tbl.colnames:
        outlines.extend(info.pformat(max_width=-1, max_lines=-1, show_unit=False))
    else:
        outlines.append("<No columns>")

    out.writelines(outline + os.linesep for outline in outlines)


class TableInfo(DataInfo):
    def __call__(self, option="attributes", out=""):
        return table_info(self._parent, option, out)

    __call__.__doc__ = table_info.__doc__


@contextmanager
def serialize_method_as(tbl, serialize_method):
    """Context manager to temporarily override individual
    column info.serialize_method dict values.  The serialize_method
    attribute is an optional dict which might look like ``{'fits':
    'jd1_jd2', 'ecsv': 'formatted_value', ..}``.

    ``serialize_method`` is a str or dict.  If str then it the the value
    is the ``serialize_method`` that will be used for all formats.
    If dict then the key values can be either:

    - Column name.  This has higher precedence than the second option of
      matching class.
    - Class (matches any column which is an instance of the class)

    This context manager is expected to be used only within ``Table.write``.
    It could have been a private method on Table but prefer not to add
    clutter to that class.

    Parameters
    ----------
    tbl : Table object
        Input table
    serialize_method : dict, str
        Dict with key values of column names or types, or str

    Returns
    -------
    None (context manager)
    """

    def get_override_sm(col):
        """
        Determine if the ``serialize_method`` str or dict specifies an
        override of column presets for ``col``.  Returns the matching
        serialize_method value or ``None``.
        """
        # If a string then all columns match
        if isinstance(serialize_method, str):
            return serialize_method

        # If column name then return that serialize_method
        if col.info.name in serialize_method:
            return serialize_method[col.info.name]

        # Otherwise look for subclass matches
        for key in serialize_method:
            if isclass(key) and isinstance(col, key):
                return serialize_method[key]

        return None

    # Setup for the context block.  Set individual column.info.serialize_method
    # values as appropriate and keep a backup copy.  If ``serialize_method``
    # is None or empty then don't do anything.

    # Original serialize_method dict, keyed by column name.  This only
    # gets used and set if there is an override.
    original_sms = {}

    if serialize_method:
        # Go through every column and if it has a serialize_method info
        # attribute then potentially update it for the duration of the write.
        for col in tbl.itercols():
            if hasattr(col.info, "serialize_method"):
                override_sm = get_override_sm(col)
                if override_sm:
                    # Make a reference copy of the column serialize_method
                    # dict which maps format (e.g. 'fits') to the
                    # appropriate method (e.g. 'data_mask').
                    original_sms[col.info.name] = col.info.serialize_method

                    # Set serialize method for *every* available format.  This is
                    # brute force, but at this point the format ('fits', 'ecsv', etc)
                    # is not actually known (this gets determined by the write function
                    # in registry.py).  Note this creates a new temporary dict object
                    # so that the restored version is the same original object.
                    col.info.serialize_method = {
                        fmt: override_sm for fmt in col.info.serialize_method
                    }

    # Finally yield for the context block
    try:
        yield
    finally:
        # Teardown (restore) for the context block.  Be sure to do this even
        # if an exception occurred.
        if serialize_method:
            for name, original_sm in original_sms.items():
                tbl[name].info.serialize_method = original_sm
