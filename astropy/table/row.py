# Licensed under a 3-clause BSD style license - see LICENSE.rst

import collections
from collections import OrderedDict
from operator import index as operator_index

import numpy as np

from astropy.utils.compat import COPY_IF_NEEDED


class Row:
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
        # Ensure that the row index is a valid index (int)
        index = operator_index(index)

        n = len(table)

        if index < -n or index >= n:
            raise IndexError(
                f"index {index} out of range for table with length {len(table)}"
            )

        # Finally, ensure the index is positive [#8422] and set Row attributes
        self._index = index % n
        self._table = table

    def __getitem__(self, item):
        try:
            # Try the most common use case of accessing a single column in the Row.
            # Bypass the TableColumns __getitem__ since that does more testing
            # and allows a list of tuple or str, which is not the right thing here.
            out = OrderedDict.__getitem__(self._table.columns, item)[self._index]
        except (KeyError, TypeError):
            if self._table._is_list_or_tuple_of_str(item):
                cols = [self._table[name] for name in item]
                out = self._table.__class__(cols, copy=False)[self._index]
            elif isinstance(item, slice):
                # https://github.com/astropy/astropy/issues/14007
                out = tuple(self.values())[item]
            else:
                # This is only to raise an exception
                out = self._table.columns[item][self._index]
        return out

    def __setitem__(self, item, val):
        if self._table._is_list_or_tuple_of_str(item):
            self._table._set_row(self._index, colnames=item, vals=val)
        else:
            self._table.columns[item][self._index] = val

    def _ipython_key_completions_(self):
        return self.colnames

    def __eq__(self, other):
        if self._table.masked:
            # Sent bug report to numpy-discussion group on 2012-Oct-21, subject:
            # "Comparing rows in a structured masked array raises exception"
            # No response, so this is still unresolved.
            raise ValueError(
                "Unable to compare rows for masked table due to numpy.ma bug"
            )
        return self.as_void() == other

    def __ne__(self, other):
        if self._table.masked:
            raise ValueError(
                "Unable to compare rows for masked table due to numpy.ma bug"
            )
        return self.as_void() != other

    def __array__(self, dtype=None, copy=COPY_IF_NEEDED):
        """Support converting Row to np.array via np.array(table).

        Coercion to a different dtype via np.array(table, dtype) is not
        supported and will raise a ValueError.

        If the parent table is masked then the mask information is dropped.
        """
        if dtype is not None:
            raise ValueError("Datatype coercion is not allowed")

        return np.array(self.as_void(), copy=copy)

    def __len__(self):
        return len(self._table.columns)

    def __iter__(self):
        index = self._index
        for col in self._table.columns.values():
            yield col[index]

    def get(self, key, default=None, /):
        """Return the value for key if key is in the columns, else default.

        Parameters
        ----------
        key : `str`, positional-only
            The name of the column to look for.
        default : `object`, optional, positional-only
            The value to return if the ``key`` is not among the columns.

        Returns
        -------
        `object`
            The value in the ``key`` column of the row if present,
            ``default`` otherwise.

        Examples
        --------
        >>> from astropy.table import Table
        >>> t = Table({"a": [2, 3, 5], "b": [7, 11, 13]})
        >>> t[0].get("a")
        2
        >>> t[1].get("b", 0)
        11
        >>> t[2].get("c", 0)
        0
        """
        return self[key] if key in self._table.columns else default

    def keys(self):
        return self._table.columns.keys()

    def values(self):
        return self.__iter__()

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
        void_row : ``numpy.void`` or ``numpy.ma.mvoid``
            Copy of row values.
            ``numpy.void`` if unmasked, ``numpy.ma.mvoid`` else.
        """
        index = self._index
        cols = self._table.columns.values()
        vals = tuple(np.asarray(col)[index] for col in cols)
        if self._table.masked:
            mask = tuple(
                col.mask[index] if hasattr(col, "mask") else False for col in cols
            )
            void_row = np.ma.array([vals], mask=[mask], dtype=self.dtype)[0]
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
        table = self._table[index : index + 1]
        descr_vals = [self.__class__.__name__, f"index={self.index}"]
        if table.masked:
            descr_vals.append("masked=True")

        return table._base_repr_(
            html, descr_vals, max_width=-1, tableid=f"table{id(self._table)}"
        )

    def _repr_html_(self):
        return self._base_repr_(html=True)

    def __repr__(self):
        return self._base_repr_(html=False)

    def __str__(self):
        index = self.index if (self.index >= 0) else self.index + len(self._table)
        return "\n".join(self.table[index : index + 1].pformat(max_width=-1))

    def __bytes__(self):
        return str(self).encode("utf-8")


collections.abc.Sequence.register(Row)
