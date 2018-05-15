# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import collections
import operator

import numpy as np

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
      <Row index=1>
        a     b
      int32 int32
      ----- -----
          2     4
      >>> row['a']
      2
      >>> row[1]
      4
    """

    def __init__(self, table, index):
        self._table = table
        self._index = operator.index(index)

        n = len(table)
        if index < -n or index >= n:
            raise IndexError('index {0} out of range for table with length {1}'
                             .format(index, len(table)))

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
        return self.as_void() == other

    def __ne__(self, other):
        if self._table.masked:
            raise ValueError('Unable to compare rows for masked table due to numpy.ma bug')
        return self.as_void() != other

    def __array__(self, dtype=None):
        """Support converting Row to np.array via np.array(table).

        Coercion to a different dtype via np.array(table, dtype) is not
        supported and will raise a ValueError.

        If the parent table is masked then the mask information is dropped.
        """
        if dtype is not None:
            raise ValueError('Datatype coercion is not allowed')

        return np.asarray(self.as_void())

    def __len__(self):
        return len(self._table.columns)

    def __iter__(self):
        index = self._index
        for col in six.itervalues(self._table.columns):
            yield col[index]

    @property
    def table(self):
        return self._table

    @property
    def index(self):
        return self._index

    def as_void(self):
        """
        Returns a *read-only* copy of the row values in the form of np.void or
        np.ma.mvoid objects.  This corresponds to the object types returned for
        row indexing of a pure numpy structured array or masked array. This
        method is slow and its use is discouraged when possible.

        Returns
        -------
        void_row : np.void (unmasked) or np.ma.mvoid (masked)
            Copy of row values
        """
        index = self._index
        cols = self._table.columns.values()
        vals = tuple(np.asarray(col)[index] for col in cols)
        if self._table.masked:
            # The logic here is a little complicated to work around
            # bug in numpy < 1.8 (numpy/numpy#483).  Need to build up
            # a np.ma.mvoid object by hand.
            from .table import descr

            # Make np.void version of masks.  Use the table dtype but
            # substitute bool for data type
            masks = tuple(col.mask[index] if hasattr(col, 'mask') else False
                          for col in cols)
            descrs = (descr(col) for col in cols)
            mask_dtypes = [(name, np.bool, shape) for name, type_, shape in descrs]
            row_mask = np.array([masks], dtype=mask_dtypes)[0]

            # Make np.void version of values, and then the final mvoid row
            row_vals = np.array([vals], dtype=self.dtype)[0]
            void_row = np.ma.mvoid(data=row_vals, mask=row_mask)
        else:
            void_row = np.array([vals], dtype=self.dtype)[0]
        return void_row

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

    def _base_repr_(self, html=False):
        """
        Display row as a single-line table but with appropriate header line.
        """
        index = self.index if (self.index >= 0) else self.index + len(self._table)
        table = self._table[index:index + 1]
        descr_vals = [self.__class__.__name__,
                      'index={0}'.format(self.index)]
        if table.masked:
            descr_vals.append('masked=True')

        return table._base_repr_(html, descr_vals, max_width=-1,
                                 tableid='table{0}'.format(id(self._table)))

    def _repr_html_(self):
        return self._base_repr_(html=True)

    def __repr__(self):
        return self._base_repr_(html=False)

    def __unicode__(self):
        index = self.index if (self.index >= 0) else self.index + len(self._table)
        return '\n'.join(self.table[index:index + 1].pformat(max_width=-1))
    if not six.PY2:
        __str__ = __unicode__

    def __bytes__(self):
        return six.text_type(self).encode('utf-8')
    if six.PY2:
        __str__ = __bytes__


collections.Sequence.register(Row)
