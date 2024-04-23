#cython: language_level=3
cimport cython
cimport numpy as np
from libc cimport math
from libc cimport stdlib

np.import_array()

ctypedef fused dtype:
    short
    int
    long
    float
    double


cdef bint is_invalid_lat(dtype angle, dtype limit):
    if dtype is float or dtype is double:
        angle = math.fabs(angle)
    else:
        angle = stdlib.abs(angle)
    return angle > limit

@cython.boundscheck(False)
@cython.wraparound(False)
cpdef bint any_invalid_lat(const dtype[:] angles, dtype limit):
    cdef ssize_t size
    size = angles.size
    for i in range(size):
        if is_invalid_lat(angles[i], limit):
            return True
    return False
