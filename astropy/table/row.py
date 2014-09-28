# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import collections
from distutils import version

import numpy as np
from numpy import ma

from ..extern import six

class Row(object):
    """A class to represent one row of a Table object.

    A Row object is returned when a Table object is indexed with an integer
    or when iterating over a table::

      >>> from astropy.table import Table
      >>> table = Table([(1, 2), (3, 4)], names=('a', 'b'),
      ...               dtype=('int32', 'int32'))
      >>> row = table[1]
      >>> row
      <Row 1 of table
       values=(2, 4)
       dtype=[('a', '<i4'), ('b', '<i4')]>
      >>> row['a']
      2
      >>> row[1]
      4
    """

    def __init__(self, table, index):
        self._table = table
        self._index = index

        # TODO: Make sure these corner cases are still handled.
        #
        # try:
        #     self._data = self._table._data[0]
        #     # MaskedArray __getitem__ has a strange behavior where if a
        #     # row mask is all False then it returns a np.void which
        #     # has no mask attribute. This makes it impossible to then set
        #     # the mask. Here we recast back to mvoid. This was fixed in
        #     # Numpy following issue numpy/numpy#483, and the fix should be
        #     # included in Numpy 1.8.0.
        #     if self._table.masked and isinstance(self._data, np.void):
        #         self._data = ma.core.mvoid(self._data,
        #                                    mask=self._table._mask[index])
        # except ValueError as err:
        #     # Another bug (or maybe same?) that is fixed in 1.8 prevents
        #     # accessing a row in masked array if it has object-type members.
        #     # >>> x = np.ma.empty(1, dtype=[('a', 'O')])
        #     # >>> x['a'] = 1
        #     # >>> x['a'].mask = True
        #     # >>> x[0]
        #     # ValueError: Setting void-array with object members using buffer. [numpy.ma.core]
        #     #
        #     # All we do here is re-raise with a more informative message
        #     if (six.text_type(err).startswith('Setting void-array with object members')
        #             and version.LooseVersion(np.__version__) < version.LooseVersion('1.8')):
        #         raise ValueError('Cannot access table row with Object type columns, due to '
        #                          'a bug in numpy {0}.  Please upgrade to numpy 1.8 or newer.'
        #                          .format(np.__version__))
        #     else:
        #         raise


    def __getitem__(self, item):
        return self._table.columns[item][self._index]

    def __setitem__(self, item, val):
        self._table.columns[item][self._index] = val

    def __eq__(self, other):
        if self._table.masked:
            # Sent bug report to numpy-discussion group on 2012-Oct-21, subject:
            # "Comparing rows in a structured masked array raises exception"
            # No response, so this is still unresolved.
            raise ValueError('Unable to compare rows for masked table due to numpy.ma bug')
        return self.data == other

    def __ne__(self, other):
        if self._table.masked:
            raise ValueError('Unable to compare rows for masked table due to numpy.ma bug')
        return self.data != other

    @property
    def _mask(self):
        return self.data.mask

    def __array__(self, dtype=None):
        """Support converting Row to np.array via np.array(table).

        Coercion to a different dtype via np.array(table, dtype) is not
        supported and will raise a ValueError.
        """
        if dtype is not None:
            raise ValueError('Datatype coercion is not allowed')

        return np.array(self.data, dtype=self.dtype)

    def __len__(self):
        return len(self._table.columns)

    @property
    def table(self):
        return self._table

    @property
    def index(self):
        return self._index

    @property
    def data(self):
        if not hasattr(self, '_data'):
            index = self._index
            cols = self._table.columns.values()
            vals = tuple(col[index] for col in cols)
            if self._table.masked:
                mask = tuple(col.mask[index] for col in cols)
                self._data = np.ma.array([vals], mask=[mask], dtype=self.dtype)[0]
            else:
                self._data = np.array([vals], dtype=self.dtype)[0]
        return self._data

    @property
    def meta(self):
        return self._table.meta

    @property
    def columns(self):
        return self._table.columns

    @property
    def colnames(self):
        return self._table.colnames

    @property
    def dtype(self):
        return self._table.dtype

    def __repr__(self):
        return "<{3} {0} of table\n values={1!r}\n dtype={2}>".format(
            self.index, self.data, self.dtype, self.__class__.__name__)


collections.Sequence.register(Row)
