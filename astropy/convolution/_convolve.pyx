# cython: language_level=3

cimport numpy as np
from libcpp cimport bool

np.import_array()


cdef extern from "src/convolve.h":
    void convolveNd_c(np.float64_t * const result,
            const np.float64_t * const f,
            const unsigned n_dim,
            const size_t * const image_shape,
            const np.float64_t * const g,
            const size_t * const kernel_shape,
            const bool nan_interpolate,
            const bool embed_result_within_padded_region,
            const unsigned n_threads) nogil


def _convolveNd_c(np.ndarray result,
                  np.ndarray array_to_convolve,
                  np.ndarray kernel,
                  bool nan_interpolate,
                  bool embed_result_within_padded_region,
                  int n_threads):
    # Cache pointer/shape attributes under the GIL so we can release it
    # around the C call. ``convolveNd_c`` is declared ``nogil`` in the C
    # header, so this is safe; releasing the GIL lets threaded executors
    # (dask, joblib threading backend, ...) run multiple convolutions in
    # parallel instead of being serialised at the Python layer.
    cdef:
        np.float64_t * result_ptr = <np.float64_t*>np.PyArray_DATA(result)
        np.float64_t * array_ptr = <np.float64_t*>np.PyArray_DATA(array_to_convolve)
        unsigned ndim = array_to_convolve.ndim
        size_t * array_shape = <size_t*>array_to_convolve.shape
        np.float64_t * kernel_ptr = <np.float64_t*>np.PyArray_DATA(kernel)
        size_t * kernel_shape = <size_t*>kernel.shape
    with nogil:
        convolveNd_c(
            result_ptr,
            array_ptr,
            ndim,
            array_shape,
            kernel_ptr,
            kernel_shape,
            nan_interpolate,
            embed_result_within_padded_region,
            n_threads,
        )
