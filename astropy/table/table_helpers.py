# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Helper functions for table development, mostly creating useful
tables for testing.
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from itertools import cycle
import string
import numpy as np

from .table import Table, Column
from ..extern.six.moves import zip, range
from ..utils.data_info import DataInfo, BaseInfo

class TimingTables(object):
    """
    Object which contains two tables and various other attributes that
    are useful for timing and other API tests.
    """
    def __init__(self, size=1000, masked=False):
        self.masked = masked

        # Initialize table
        self.table = Table(masked=self.masked)

        # Create column with mixed types
        np.random.seed(12345)
        self.table['i'] = np.arange(size)
        self.table['a'] = np.random.random(size)  # float
        self.table['b'] = np.random.random(size) > 0.5  # bool
        self.table['c'] = np.random.random((size,10))  # 2d column
        self.table['d'] = np.random.choice(np.array(list(string.ascii_letters)),size)

        self.extra_row = {'a':1.2, 'b':True, 'c':np.repeat(1, 10), 'd':'Z'}
        self.extra_column = np.random.randint(0, 100, size)
        self.row_indices = np.where(self.table['a'] > 0.9)[0]
        self.table_grouped = self.table.group_by('d')

        # Another table for testing joining
        self.other_table = Table(masked=self.masked)
        self.other_table['i'] = np.arange(1,size,3)
        self.other_table['f'] = np.random.random()
        self.other_table.sort('f')

        # Another table for testing hstack
        self.other_table_2 = Table(masked=self.masked)
        self.other_table_2['g'] = np.random.random(size)
        self.other_table_2['h'] = np.random.random((size, 10))

        self.bool_mask = self.table['a'] > 0.6


def simple_table(size=3, cols=None, kinds='ifS', masked=False):
    """
    Return a simple table for testing.

    Example
    --------
    ::

      >>> from astropy.table.table_helpers import simple_table
      >>> print(simple_table(3, 6, masked=True, kinds='ifOS'))
       a   b     c      d   e   f
      --- --- -------- --- --- ---
       -- 1.0 {'c': 2}  --   5 5.0
        2 2.0       --   e   6  --
        3  -- {'e': 4}   f  -- 7.0

    Parameters
    ----------
    size : int
        Number of table rows
    cols : int, default=number of kinds
        Number of table columns
    kinds : str
        String consisting of the column dtype.kinds.  This string
        will be cycled through to generate the column dtype.
        The allowed values are 'i', 'f', 'S', 'O'.

    Returns
    -------
    out : `Table`
        New table with appropriate characteristics
    """
    if cols is None:
        cols = len(kinds)
    if cols > 26:
        raise ValueError("Max 26 columns in SimpleTable")

    columns = []
    names = [chr(ord('a') + ii) for ii in range(cols)]
    letters = np.array([c for c in string.ascii_letters])
    for jj, kind in zip(range(cols), cycle(kinds)):
        if kind == 'i':
            data = np.arange(1, size + 1, dtype=int) + jj
        elif kind == 'f':
            data = np.arange(size, dtype=float) + jj
        elif kind == 'S':
            indices = (np.arange(size) + jj) % len(letters)
            data = letters[indices]
        elif kind == 'O':
            indices = (np.arange(size) + jj) % len(letters)
            vals = letters[indices]
            data = [{val: index} for val, index in zip(vals, indices)]
        else:
            raise ValueError('Unknown data kind')
        columns.append(Column(data, dtype=kind))

    table = Table(columns, names=names, masked=masked)
    if masked:
        for ii, col in enumerate(table.columns.values()):
            mask = np.array((np.arange(size) + ii) % 3, dtype=bool)
            col.mask = ~mask

    return table


def complex_table():
    """
    Return a masked table from the io.votable test set that has a wide variety
    of stressing types.
    """
    from ..utils.data import get_pkg_data_filename
    from ..io.votable.table import parse
    import warnings


    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        votable = parse(get_pkg_data_filename('../io/votable/tests/data/regression.xml'),
                        pedantic=False)
    first_table = votable.get_first_table()
    table = first_table.to_table()

    return table

class ArrayWrapper(object):
    """
    Minimal mixin using a simple wrapper around a numpy array
    """
    def __init__(self, data):
        self.data = np.array(data)
        self.info.dtype = self.data.dtype

    def __getitem__(self, item):
        if isinstance(item, (int, np.integer)):
            out = self.data[item]
        else:
            out = self.__class__(self.data[item])
        return out

    def __setitem__(self, item, value):
        self.data[item] = value

    def __len__(self):
        return len(self.data)

    info = DataInfo(BaseInfo)

    @property
    def shape(self):
        return self.data.shape

    def __repr__(self):
        return ("<{0} name='{1}' data={2}>"
                .format(self.__class__.__name__, self.info.name, self.data))
