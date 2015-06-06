"""
Table method and functions related to providing information about
columns and tables.
"""
# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import sys
import os

import numpy as np

from ..extern import six
from ..utils import OrderedDict

__all__ = ['column_info_factory']

OPTIONS = ('meta', 'stats')

def column_info_factory(names, funcs):
    """
    Factory to create a function that can be used as an ``option``
    for outputting table or column summary information.

    Examples
    --------
    >>> from astropy.table.info import column_info_factory
    >>> from astropy.table import Column
    >>> c = Column([4., 3., 2., 1.])
    >>> mystats = column_info_factory(names=['min', 'median', 'max'],
    ...                               funcs=[np.min, np.median, np.max])
    >>> c.info(option=mystats)
    min = 1.0
    median = 2.5
    max = 4.0

    Parameters
    ----------
    names: list
        List of information attribute names
    funcs: list
        List of functions that compute the corresponding information attribute

    Returns
    -------
    func: function
        Function that can be used as a table or column info option
    """
    def func(col):
        outs = []
        for name, func in zip(names, funcs):
            try:
                out = func(col)
            except:
                outs.append('--')
            else:
                outs.append(str(out))

        return OrderedDict(zip(names, outs))
    return func


def meta(col):
    """
    Info function that returns basic metadata.
    """
    from .column import BaseColumn, col_getattr

    attrs = ('dtype', 'unit', 'format', 'description', 'class')
    info = OrderedDict()
    for attr in attrs:
        if attr == 'class':
            val = '' if isinstance(col, BaseColumn) else col.__class__.__name__
        elif attr == 'dtype':
            val = col_getattr(col, attr).name
        else:
            val = col_getattr(col, attr)
        if val is None:
            val = ''
        info[attr] = str(val)

    return info


stats = column_info_factory(names=['min', 'mean', 'max'],
                            funcs=[np.min, np.mean, np.max])

def column_info(self, option='meta', out=''):
    """
    Write summary information about column to the ``out`` filehandle.
    By default this prints to standard output via sys.stdout.

    The ``option` argument specifies what type of information
    to include.  This can be a string, a function, or a list of
    strings or functions.  Built-in options are:

    - ``meta``: basic column meta data like ``dtype`` or ``format``
    - ``stats``: basic statistics: minimum, mean, and maximum

    If a function is specified then that function will be called with the
    column as its single argument.  The function must return an OrderedDict
    containing the information attributes.

    If a list is provided then the information attributes will be
    appended for each of the options, in order.

    Examples
    --------

    >>> from astropy.table import Column
    >>> c = Column([1, 2, 3], unit='m', dtype='int32')
    >>> c.info()
    dtype = int64
    unit = m
    >>> c.info(['meta', 'stats'])
    dtype = int32
    unit = m
    min = 1
    mean = 2.0
    max = 3

    Parameters
    ----------
    option: str, function, list of (str or function)
        Info option (default='meta')
    out: file-like object, None
        Output destination (default=sys.stdout).  If None then the
        OrderedDict with information attributes is returned

    Returns
    -------
    info: OrderedDict if out==None else None
    """
    from .column import col_getattr

    if out == '':
        out = sys.stdout

    info = OrderedDict()
    name=col_getattr(self, 'name')
    if name is not None:
        info['name'] = name

    options = option if isinstance(option, (list, tuple)) else [option]
    for option in options:
        if isinstance(option, six.string_types):
            if option in OPTIONS:
                option = globals()[option]
            else:
                raise ValueError('option={0} is not in allowed values for option: {1}',
                                 option, ', '.join(OPTIONS))
        info.update(option(self))

    if out is None:
        return info

    for key, val in info.items():
        if val != '':
            out.write('{0} = {1}\n'.format(key, val))


def table_info(self, option='meta', out=''):
    """
    Write summary information about column to the ``out`` filehandle.
    By default this prints to standard output via sys.stdout.

    The ``option` argument specifies what type of information
    to include.  This can be a string, a function, or a list of
    strings or functions.  Built-in options are:

    - ``meta``: basic column meta data like ``dtype`` or ``format``
    - ``stats``: basic statistics: minimum, mean, and maximum

    If a function is specified then that function will be called with the
    column as its single argument.  The function must return an OrderedDict
    containing the information attributes.

    If a list is provided then the information attributes will be
    appended for each of the options, in order.

    Examples
    --------
    >>> from astropy.table.table_helpers import simple_table
    >>> t = simple_table()
    >>> t['a'].unit = 'm'
    >>> t.info()
    <Table length=3>
    name  dtype  unit
    ---- ------- ----
       a   int32    m
       b float32
       c string8

    >>> t.info('stats')
    <Table length=3>
    name min mean max
    ---- --- ---- ---
       a   1  2.0   3
       b 1.0  2.0 3.0
       c  --   --  --

    Parameters
    ----------
    option: str, function, list of (str or function)
        Info option (default='meta')
    out: file-like object, None
        Output destination (default=sys.stdout).  If None then the
        OrderedDict with information attributes is returned
    """
    from .table import Table

    if out == '':
        out = sys.stdout

    descr_vals = [self.__class__.__name__]
    if self.masked:
        descr_vals.append('masked=True')
    descr_vals.append('length={0}'.format(len(self)))

    outlines = ['<' + ' '.join(descr_vals) + '>']

    infos = []
    for col in self.columns.values():
        # SHOULD BE as follows after #3731
        # infos.append(col.info(option, out='object'))
        infos.append(column_info(col, option, out=None))

    info = Table(infos, names=list(infos[0].keys()))
    for name in info.colnames:
        if np.all(info[name] == ''):
            del info[name]

    outlines.extend(info.pformat(max_width=-1, max_lines=-1, show_unit=False))
    out.writelines(outline + os.linesep for outline in outlines)
