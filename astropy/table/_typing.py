# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Type aliases shared by the modules of `astropy.table`.

This module is private.  The aliases defined here are an implementation detail
of the annotations in `astropy.table` and are not part of the public API.  They
are meant to be imported inside an ``if TYPE_CHECKING:`` block.
"""

__all__ = ["ColumnLike", "DataLike", "SortKind", "TableLike"]

from typing import TYPE_CHECKING, Any, Literal

if TYPE_CHECKING:
    from .column import Column, MaskedColumn

type ColumnLike = Column | MaskedColumn | Any
"""An actual table column.

This is a `~astropy.table.Column`, a `~astropy.table.MaskedColumn`, or a mixin
column such as `~astropy.units.Quantity`, `~astropy.time.Time` or
`~astropy.coordinates.SkyCoord`.  Mixin columns share no common base class (the
only requirement is a working ``info`` attribute), hence the ``Any``.
"""

type DataLike = Any
"""Anything that can be turned into a table column.

This includes a `ColumnLike` object, a `~numpy.ndarray`, a plain sequence, or a
scalar / length-1 object that gets broadcast to the table length.
"""

type TableLike = Any
"""A `~astropy.table.Table` or anything that can be used to initialize one.

For example a `dict` of columns, a list of rows, a structured
`~numpy.ndarray`, or an object implementing ``__astropy_table__``.
"""

type SortKind = Literal["quicksort", "mergesort", "heapsort", "stable"]
"""Sorting algorithm accepted by `numpy.argsort`."""
