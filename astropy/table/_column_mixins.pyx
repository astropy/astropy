#cython: language_level=3
"""
This module provides mixin bases classes for the Column and MaskedColumn
classes to provide those classes with their custom __getitem__ implementations.

The reason for this is that implementing a __getitem__ in pure Python actually
significantly slows down the array subscript operation, especially if it needs
to call the subclass's __getitem__ (i.e. ndarray.__getitem__ in this case).  By
providing __getitem__ through a base type implemented in C, the __getitem__
implementation will go straight into the class's tp_as_mapping->mp_subscript
slot, rather than going through a class __dict__ and calling a pure Python
method.  Furthermore, the C implementation of __getitem__ can easily directly
call the base class's implementation (as seen in _ColumnGetitemShim, which
directly calls to ndarray->tp_as_mapping->mp_subscript).

The main reason for overriding __getitem__ in the Column class is for
returning elements out of a multi-dimensional column.  That is, if the
elements of a Column are themselves arrays, the default ndarray.__getitem__
applies the subclass to those arrays, so they are returned as Column instances
(when really they're just an array that was in a Column).  This overrides that
behavior in the case where the element returned from a single row of the
Column is itself an array.

This module also provides helpers that scan sequences of Python objects to pick a
string dtype for ``Column`` and ``Table.from_pandas`` / ``Table.from_df``.  These
scans visit every element, so they are written in Cython to keep them fast.
"""

import sys

import numpy as np
cimport numpy as np
from cpython.ref cimport PyTypeObject
from libc.math cimport isnan

np.import_array()


cdef tuple INTEGER_TYPES = (int, np.integer)
cdef tuple STRING_TYPES = (str, np.str_)


# Annoying boilerplate that we shouldn't have to write; Cython should
# have this built in (some versions do, but the ctypedefs are still lacking,
# or what is available is Cython version dependent)
ctypedef object (*binaryfunc)(object, object)


cdef extern from "Python.h":
    void *PyType_GetSlot(PyTypeObject *type, int slot)
    int Py_mp_subscript

ctypedef object (*item_getter)(object, object)


cdef inline object base_getitem(object self, object item, item_getter getitem):
    if (<np.ndarray>self).ndim > 1 and isinstance(item, INTEGER_TYPES):
        return self.data[item]

    dtype_kind = (<np.ndarray>self).dtype.kind
    if dtype_kind == 'V' and isinstance(item, STRING_TYPES):
        return self.data[item]

    value = getitem(self, item)

    try:
        if dtype_kind == 'S' and not value.shape:
            value = value.decode('utf-8', errors='replace')
    except AttributeError:
        pass

    return value


cdef inline object column_getitem(object self, object item):
    return (<binaryfunc>PyType_GetSlot(<PyTypeObject *>np.ndarray, Py_mp_subscript))(self, item)


cdef class _ColumnGetitemShim:
    def __getitem__(self, item):
        return base_getitem(self, item, column_getitem)


MaskedArray = np.ma.MaskedArray


cdef inline object masked_column_getitem(object self, object item):
    value = MaskedArray.__getitem__(self, item)
    return self._copy_attrs_slice(value)


cdef class _MaskedColumnGetitemShim(_ColumnGetitemShim):
    def __getitem__(self, item):
        return base_getitem(self, item, masked_column_getitem)


def _get_string_dtype(value):
    """
    Choose a string dtype for a flat sequence without NumPy bytes/text promotion.

    This is used by ``Column._encode_str`` when a value is assigned to, compared
    with, or used to create a bytes (``S``) column.  Letting NumPy infer the dtype
    of a sequence that mixes ``bytes`` and ``str`` relies on implicit bytes/text
    promotion, which NumPy is deprecating (numpy/numpy#32757).  Picking the dtype
    here keeps the existing results without relying on that promotion.

    Parameters
    ----------
    value : object
        Value to inspect.  Only a ``list`` or ``tuple`` whose elements are all
        ``bytes``, ``str`` or numeric scalars is classified.

    Returns
    -------
    dtype : {"S", "U", None}
        ``"S"`` if every element is ``bytes``, ``"U"`` if at least one element is
        ``str``, and None otherwise.  None means the caller should let NumPy
        infer the dtype as before.  This includes nested sequences, arrays,
        sequences of only bytes and numbers, and values that are not a list or
        tuple.
    """
    cdef bint has_text = False
    cdef bint all_bytes = True
    # Type of the last element that was classified as text or a number.  Long
    # sequences usually repeat one type, so later elements of that exact type
    # skip the isinstance checks below.
    cdef type scalar_type = None

    if not isinstance(value, (list, tuple)):
        return None
    for item in value:
        if isinstance(item, bytes) or type(item) is scalar_type:
            continue
        if isinstance(item, str):
            has_text = True
        elif not isinstance(item, (int, float, complex, np.number, np.bool_)):
            # Let NumPy handle nested sequences and other types.
            return None
        scalar_type = type(item)
        all_bytes = False
    # Numbers are allowed alongside text because NumPy converts them to strings
    # when it creates a "U" array.
    if has_text:
        return "U"
    return "S" if all_bytes else None


def _get_string_fill_value(data, *, bint legacy_pandas=False):
    """
    Choose a placeholder for missing values in a string-like object column.

    ``Table.from_pandas`` and ``Table.from_df`` convert an object column whose
    values are all strings or missing values into a native string column.  The
    caller builds a list with the placeholder in place of each missing value, and
    NumPy then infers a ``U`` or ``S`` dtype from that list.  Building a new list
    also avoids modifying a view of the dataframe's array.

    The placeholder matches the type of the first string found.  Using ``b""`` for
    a text column would rely on NumPy's implicit bytes/text promotion, which NumPy
    is deprecating (numpy/numpy#32757).

    Parameters
    ----------
    data : ndarray
        Object array of column values.
    legacy_pandas : bool, optional
        If True, apply the rules of ``Table.from_pandas``, where only the
        ``numpy.nan`` object itself counts as a missing value.  If False, apply the
        rules of ``Table.from_df``, where None and any NaN value count as missing.
        Default is False.

    Returns
    -------
    fill : str, bytes or None
        ``""`` if the first string is ``str``, ``b""`` if it is ``bytes`` or if
        every value is missing, and None if any value is neither a string nor
        missing.  None means the column should stay an object column.
    """
    cdef object fill = None
    cdef object nan = np.nan

    for value in data:
        if isinstance(value, (str, bytes)):
            if fill is None:
                fill = "" if isinstance(value, str) else b""
        elif legacy_pandas:
            if value is not nan:
                return None
        elif value is not None and not isnan(value):
            return None
    return b"" if fill is None else fill
