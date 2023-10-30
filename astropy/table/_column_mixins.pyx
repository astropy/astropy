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
"""

import sys

import numpy as np
cimport numpy as np

np.import_array()


cdef tuple INTEGER_TYPES = (int, np.integer)
cdef tuple STRING_TYPES = (str, np.str_)


# Annoying boilerplate that we shouldn't have to write; Cython should
# have this built in (some versions do, but the ctypedefs are still lacking,
# or what is available is Cython version dependent)
ctypedef object (*binaryfunc)(object, object)


cdef extern from "Python.h":
    ctypedef struct PyMappingMethods:
        binaryfunc mp_subscript

    ctypedef struct PyTypeObject:
        PyMappingMethods* tp_as_mapping


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
    return (<PyTypeObject *>np.ndarray).tp_as_mapping.mp_subscript(self, item)


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
