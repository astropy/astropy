"""
Table property for providing information about table.
"""
# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import sys
import os

import numpy as np
from ..extern import six
from ..utils.data_info import ParentDtypeInfo

__all__ = ['table_info', 'TableInfo']

def table_info(tbl, option='attributes', out=''):
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
       a   int32    m
       b float32

    >>> t.info('stats')
    <Table length=2>
    name mean std min max
    ---- ---- --- --- ---
       a  1.5 0.5   1   2
       b  1.5 0.5 1.0 2.0

    Parameters
    ----------
    option: str, function, list of (str or function)
        Info option (default='attributes')
    out: file-like object, None
        Output destination (default=sys.stdout).  If None then a
        Table with information attributes is returned

    Returns
    -------
    info: `~astropy.table.Table` if out==None else None
    """
    from .table import Table

    if out == '':
        out = sys.stdout

    descr_vals = [tbl.__class__.__name__]
    if tbl.masked:
        descr_vals.append('masked=True')
    descr_vals.append('length={0}'.format(len(tbl)))

    outlines = ['<' + ' '.join(descr_vals) + '>']

    cols = tbl.columns.values()
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
    for name in info.colnames:
        if np.all(info[name] == ''):
            del info[name]

    if 'class' in info.colnames:
        # Remove 'class' info column if all table columns are the same class
        # and they are the default column class for that table.
        uniq_types = set(type(col) for col in cols)
        if len(uniq_types) == 1 and isinstance(cols[0], tbl.ColumnClass):
            del info['class']

    if 'n_bad' in info.colnames and np.all(info['n_bad'] == 0):
        del info['n_bad']

    # Standard attributes has 'length' but this is typically redundant
    if 'length' in info.colnames and np.all(info['length'] == len(tbl)):
        del info['length']

    if tbl.colnames:
        outlines.extend(info.pformat(max_width=-1, max_lines=-1, show_unit=False))
    else:
        outlines.append('<No columns>')

    out.writelines(outline + os.linesep for outline in outlines)

class TableInfo(ParentDtypeInfo):
    _parent = None

    @staticmethod
    def default_format(val):
        return str(val.as_void())

    @property
    def unit(self):
        """Use the unit row to show column names"""
        unit = '(' + ', '.join(self._parent.colnames) + ')'
        return unit

    def __call__(self, option='attributes', out=''):
        return table_info(self._parent, option, out)

    __call__.__doc__ = table_info.__doc__

    def __repr__(self):
        out = six.moves.cStringIO()
        self.__call__(out=out)
        return out.getvalue()
