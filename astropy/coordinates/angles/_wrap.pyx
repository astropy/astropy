#cython: language_level=3
cimport cython
cimport numpy as np
from libc cimport math

np.import_array()

ctypedef fused dtype:
    short
    int
    long
    float
    double


@cython.ufunc
@cython.cdivision(True)
cdef dtype wrap_at(dtype angle, dtype wrap_angle, dtype period):
    cdef dtype wrap_angle_floor
    wrap_angle_floor = wrap_angle - period

    if dtype is float or dtype is double:
        if math.isnan(angle):
            return angle

    if (angle >= wrap_angle_floor) and (angle < wrap_angle):
        return angle

    cdef int wraps
    wraps = int((angle - wrap_angle_floor) // period)

    angle -= wraps * period

    # Rounding errors can cause problems.
    if angle >= wrap_angle:
        angle -= period
    elif angle < wrap_angle_floor:
        angle += period

    return angle


@cython.ufunc
@cython.cdivision(True)
cdef dtype needs_wrapping(dtype angle, dtype wrap_angle, dtype period):
    cdef dtype wrap_angle_floor
    wrap_angle_floor = wrap_angle - period

    if dtype is float or dtype is double:
        if math.isnan(angle):
            return False

    if (angle >= wrap_angle_floor) and (angle < wrap_angle):
        return False
    return True
