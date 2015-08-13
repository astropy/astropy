import sys
import numpy as np
cimport numpy as cnp


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


ndarray = np.ndarray
ctypedef cnp.ndarray ndarray_t


cdef class _ColumnGetitemShim:
    def __getitem__(self, item):
        if (<ndarray_t>self).ndim > 1 and isinstance(item, INTEGER_TYPES):
            return self.data[item]
        else:
            return (<PyTypeObject *>ndarray).tp_as_mapping.mp_subscript(self, item)


MaskedArray = np.ma.MaskedArray


cdef class _MaskedColumnGetitemShim(_ColumnGetitemShim):
    def __getitem__(self, item):
        if (<ndarray_t>self).ndim > 1 and isinstance(item, INTEGER_TYPES):
            return self.data[item]
        else:
            return MaskedArray.__getitem__(self, item)
