# cython: language_level = 3
import numpy as np
cimport numpy as np
cimport cython
from libc.math cimport isfinite as _isfinite, M_PI
from libc.stdint cimport int32_t, int64_t


ctypedef fused dtype_t:
    int32_t
    int64_t
    float
    double


cdef bint isfinite(dtype_t val):
    if dtype_t is int32_t or dtype_t is int64_t:
        return True
    else:
        return _isfinite(val)

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
@cython.cdivision(True)
cpdef _wrap_at(np.ndarray[dtype_t] angle, dtype_t wrap_angle, double full_circle):
    cdef np.npy_intp size
    cdef np.npy_intp i
    cdef dtype_t wrap_angle_floor
    cdef dtype_t wrap_by

    wrap_angle_floor = <dtype_t>(wrap_angle - full_circle)

    size = angle.shape[0]
    for i in range(size):
        wrap_by =  <dtype_t>(full_circle * ((angle[i] - wrap_angle_floor) // full_circle))

        if wrap_by != 0 and isfinite(wrap_by):
            angle[i] -= wrap_by

        if angle[i] < wrap_angle_floor:
            angle[i] = <dtype_t>(angle[i] + full_circle)

        if angle[i] >= wrap_angle:
            angle[i] = <dtype_t>(angle[i] - full_circle)


cpdef bint _needs_wrapping(np.ndarray[dtype_t] angle, dtype_t wrap_angle, double full_circle):
    cdef np.npy_intp size
    cdef np.npy_intp i
    cdef dtype_t wrap_angle_floor

    size = angle.shape[0]
    wrap_angle_floor = <dtype_t>(wrap_angle - full_circle)
    for i in range(size):
        if not isfinite(angle[i]):
            continue

        if angle[i] >= wrap_angle or angle[i] < wrap_angle_floor:
            return True

    return False
