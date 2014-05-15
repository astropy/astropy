"""
Cython utilities for numpy structured arrays.

join_inner():  Do the inner-loop cartesian product for np_utils.join() processing.
               (The "inner" is about the inner loop, not inner join).
"""

import numpy as np
import numpy.ma as ma
from numpy.lib.recfunctions import drop_fields

cimport cython
cimport numpy as np
DTYPE = np.int
ctypedef np.intp_t DTYPE_t

@cython.wraparound(False)
@cython.boundscheck(False)
def join_inner(np.ndarray[DTYPE_t, ndim=1] idxs,
               np.ndarray[DTYPE_t, ndim=1] idx_sort,
               int len_left,
               int jointype):
    """
    Do the inner-loop cartesian product for np_utils.join() processing.
    (The "inner" is about the inner loop, not inner join).
    """
    cdef int n_out = 0
    cdef int max_key_idxs = 0
    cdef DTYPE_t ii, key_idxs, n_left, n_right, idx0, idx1, idx, i
    cdef DTYPE_t i_left, i_right, i_out
    cdef int masked

    # First count the final number of rows and max number of indexes
    # for a single key
    masked = 0
    for ii in range(idxs.shape[0] - 1):
        idx0 = idxs[ii]
        idx1 = idxs[ii + 1]

        # Number of indexes for this key
        key_idxs = idx1 - idx0
        if key_idxs > max_key_idxs:
            max_key_idxs = key_idxs

        # Number of rows for this key
        n_left = 0
        n_right = 0
        for idx in range(idx0, idx1):
            i = idx_sort[idx]
            if i < len_left:
                n_left += 1
            else:
                n_right += 1

        # Fix n_left and n_right for different join types
        if jointype == 0:
            pass
        elif jointype == 1:
            if n_left == 0:
                masked = 1
                n_left = 1
            if n_right == 0:
                masked = 1
                n_right = 1
        elif jointype == 2:
            if n_right == 0:
                masked = 1
                n_right = 1
        elif jointype == 3:
            if n_left == 0:
                masked = 1
                n_left = 1

        n_out += n_left * n_right

    cdef np.ndarray left_out = np.empty(n_out, dtype=DTYPE)
    cdef np.ndarray right_out = np.empty(n_out, dtype=DTYPE)
    cdef np.ndarray left_mask = np.zeros(n_out, dtype=np.bool)
    cdef np.ndarray right_mask = np.zeros(n_out, dtype=np.bool)
    cdef np.ndarray left_idxs = np.empty(max_key_idxs, dtype=DTYPE)
    cdef np.ndarray right_idxs = np.empty(max_key_idxs, dtype=DTYPE)

    i_out = 0
    for ii in range(idxs.shape[0] - 1):
        idx0 = idxs[ii]
        idx1 = idxs[ii + 1]

        # Number of rows for this key
        n_left = 0
        n_right = 0
        for idx in range(idx0, idx1):
            i = idx_sort[idx]
            if i < len_left:
                left_idxs[n_left] = i
                n_left += 1
            else:
                right_idxs[n_right] = i - len_left
                n_right += 1

        if jointype == 0:
            pass
        elif jointype == 1:
            if n_left == 0:
                left_idxs[0] = -1
                n_left = 1
            if n_right == 0:
                right_idxs[0] = -1
                n_right = 1
        elif jointype == 2:
            if n_right == 0:
                right_idxs[0] = -1
                n_right = 1
        elif jointype == 3:
            if n_left == 0:
                left_idxs[0] = -1
                n_left = 1

        for i_left in range(n_left):
            for i_right in range(n_right):
                idx = left_idxs[i_left]
                if idx < 0:
                    idx = 0
                    left_mask[i_out] = 1
                left_out[i_out] = idx

                idx = right_idxs[i_right]
                if idx < 0:
                    idx = 0
                    right_mask[i_out] = 1
                right_out[i_out] = idx

                i_out += 1

    return masked, n_out, left_out, left_mask, right_out, right_mask
