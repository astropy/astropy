#cython: language_level=3

# Fast implementation of common cases for sigma clipping

import numpy as np
cimport numpy as np

DTYPE = np.float
ctypedef np.float_t DTYPE_FLOAT
ctypedef np.uint8_t DTYPE_BOOL

cimport cython


@cython.boundscheck(False)
@cython.cdivision(True)
def sigma_clip_fast_mean(np.ndarray[DTYPE_FLOAT, ndim=2] array,
                         np.ndarray[DTYPE_BOOL, ndim=2] mask,
                        int maxiters, double sigma_lower, double sigma_upper):

    # Given a 2-d array, iterate over the second dimension and compute sigma clipping
    # for each resulting 1-d array. Returns a mask indicating which values should be removed
    # by the sigma clipping.

    cdef int m = array.shape[1]
    cdef int n = array.shape[0]

    cdef int iter, i, j
    cdef double count, count_prev, std, mean, lower, upper
    cdef double inf = np.inf

    for j in range(m):

        lower = -inf
        upper = +inf
        iter = 0

        while True:

            mean = 0
            count = 0

            for i in range(n):
                if not mask[i, j]:
                    mean += array[i, j]
                    count += 1

            if count > 0:
                mean /= count

            std = 0

            for i in range(n):
                if not mask[i, j]:
                    std += (mean - array[i, j]) ** 2

            if count > 0:
                std = (std / count) ** 0.5

            lower = mean - sigma_lower * std
            upper = mean + sigma_upper * std

            mask_changed = False

            for i in range(n):
                if not mask[i, j] and (array[i, j] < lower or array[i, j] > upper):
                    mask[i, j] = True
                    if not mask_changed:
                        mask_changed = True

            if not mask_changed:
                break

            iter += 1

            if maxiters != -1 and iter >= maxiters:
                break

    return lower, upper
