import sys
import numpy as np


if sys.version_info[0] == 3:
    INTEGER_TYPES = (int, np.integer)
else:
    INTEGER_TYPES = (int, long, np.integer)


# Annoying boilerplate that we shouldn't have to write; Cython should
# have this built in (some versions do, but the ctypedefs are still lacking,
# or what is available is Cython version dependent)
ctypedef object (*binaryfunc)(object, object)


cdef extern from "Python.h":
    ctypedef struct PyMappingMethods:
        binaryfunc mp_subscript

    ctypedef struct PyTypeObject:
        PyMappingMethods* tp_as_mapping


cdef extern from "numpy/arrayobject.h":
    ctypedef class numpy.ndarray [object PyArrayObject]:
        cdef int ndim "nd"


ctypedef object (*item_getter)(object, object)


cdef inline object base_getitem(object self, object item, item_getter getitem):
    if (<ndarray>self).ndim > 1 and isinstance(item, INTEGER_TYPES):
        return self.data[item]

    value = getitem(self, item)

    if type(value) is type(self):
        value = self.info.slice_indices(value, item, len(self))

    return value


cdef inline object column_getitem(object self, object item):
    return (<PyTypeObject *>ndarray).tp_as_mapping.mp_subscript(self, item)


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
