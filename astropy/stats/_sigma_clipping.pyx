#cython: language_level=3

# Fast implementation of common cases for sigma clipping

import numpy as np
cimport numpy as np
from libc.stdlib cimport qsort

DTYPE = np.float
ctypedef np.float_t DTYPE_FLOAT
ctypedef np.uint8_t DTYPE_BOOL
ctypedef np.uint8_t DTYPE_INT8

cimport cython


@cython.boundscheck(False)
@cython.cdivision(True)
cdef int cmp_func(const void* a, const void* b) nogil:
    """
    Comparison function for use in qsort
    """
    cdef double a_v = (<double*>a)[0]
    cdef double b_v = (<double*>b)[0]
    if a_v < b_v:
        return -1
    else:
        return 1


@cython.boundscheck(False)
@cython.cdivision(True)
cpdef sort_inplace(double[:] a, size_t count):
    """
    Sort a memoryview inplace using the native C quicksort.
    This assumes the array is contiguous in memory.
    """
    cdef int size = 8
    qsort(&a[0], count, size, &cmp_func)


@cython.boundscheck(False)
@cython.cdivision(True)
cpdef size_t locate(double[:] array, double target, size_t start, size_t end):
    """
    Find the position that the value 'target' would have if it was inserted
    into 'array', assuming 'array' is sorted.
    """
    cdef size_t mid
    while True:
        if end == start + 1:
            return start
        mid = (start + end) // 2
        if array[mid] < target:
            start = mid
        else:
            end = mid


@cython.boundscheck(False)
@cython.cdivision(True)
cpdef double calc_mean(double[:] arr, size_t start, size_t end):
    """
    Compute the mean of values in a memory view
    """
    cdef double total = 0
    cdef size_t i
    for i in range(start, end):
        total += arr[i]
    total /= (end - start)
    return total


@cython.boundscheck(False)
@cython.cdivision(True)
cpdef double calc_std(double[:] arr, size_t start, size_t end, double mean):
    """
    Compute the standard deviation of values in a memory view
    """
    cdef double total = 0
    cdef size_t i
    for i in range(start, end):
        total += (mean - arr[i]) ** 2
    total /= (end - start)
    total = total ** 0.5
    return total


@cython.boundscheck(False)
@cython.cdivision(True)
def sigma_clip_fast(np.ndarray[DTYPE_FLOAT, ndim=2] array,
                    np.ndarray[DTYPE_BOOL, ndim=2] mask,
                    int use_median,
                    int maxiters, double sigma_lower, double sigma_upper):

    cdef size_t m = array.shape[1]
    cdef size_t n = array.shape[0]
    cdef size_t i, j
    cdef np.ndarray[DTYPE_FLOAT, ndim=1] array_sorted = np.zeros(n, dtype=float)

    cdef int iteration
    cdef size_t median_index, count
    cdef size_t start_prev, end_prev, start, end
    cdef double std, mean, median, lower, upper
    cdef double inf = np.inf

    # This function is constructed to take a 2-d array of values and assumes
    # that each 1-d array when looping over the last dimension should be
    # treated separately.

    for j in range(m):

        iteration = 0

        # We copy all finite values from array into the array_sorted buffer,
        # sort it, then keep track of the range of values that are
        # not rejected by the sigma clipping. Having the values in a
        # sorted array means that the sigma clipping is removing values
        # from either or both ends of the array - [start:end] gives the
        # range of values that have not been rejected. This has the
        # advantage that we don't then need to update a data or mask
        # array in each iteration - just the start and end indices.

        start = 0
        end = 0

        for i in range(n):
            if array[i, j] == array[i, j]:
                array_sorted[end] = array[i, j]
                end += 1

        # If end == 0, no values have been copied over (this can happen
        # for example if all the values are NaN). In this case, we just
        # proceed to the next array.
        if end == 0:
            continue

        # We now sort the values in the array up to the end index.
        sort_inplace(array_sorted, end)

        while True:

            count = end - start

            # Calculate the mean and standard deviation of values so far.
            mean = calc_mean(array_sorted, start, end)
            std = calc_std(array_sorted, start, end, mean)

            # If needed, we compute the median
            if use_median == 1:

                if count % 2 == 0:
                    median_index = start + count // 2 - 1
                    median = 0.5 * (array_sorted[median_index] + array_sorted[median_index + 1])
                else:
                    median_index = start + (count - 1) // 2
                    median = array_sorted[median_index]

                lower = median - sigma_lower * std
                upper = median + sigma_upper * std

            else:

                lower = mean - sigma_lower * std
                upper = mean + sigma_upper * std

            # If all array values in the [start:end] range are still inside
            # (lower, upper) then the process has converged and we can exit
            # the loop over iterations.
            if array_sorted[start] > lower and array_sorted[end - 1] < upper:
                break

            # We need to keep track of the previous start/end values as
            # we need the original values for both locate calls.
            start_prev = start
            end_prev = end

            # Update the start/end values based on the new lower/upper values
            if array_sorted[start] < lower:
                start = locate(array_sorted, lower, start_prev, end_prev) + 1
            if array_sorted[end - 1] > upper:
                end = locate(array_sorted, upper, start_prev, end_prev) + 1

            iteration += 1

            if maxiters != -1 and iteration >= maxiters:
                break

        # Populate the final (unsorted) mask
        for i in range(n):
            mask[i, j] = array[i, j] < lower or array[i, j] > upper

    # FIXME: the following should actually be arrays, currently just
    # return the last used lower/upper values.
    return lower, upper