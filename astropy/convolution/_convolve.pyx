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
    # Cache attributes we need with the GIL before releasing it.
    cdef:
        np.float64_t * result_ptr = <np.float64_t*>np.PyArray_DATA(result)
        const np.float64_t * array_ptr = <np.float64_t*>np.PyArray_DATA(array_to_convolve)
        const unsigned ndim = array_to_convolve.ndim
        const size_t * array_shape = <size_t*>array_to_convolve.shape
        const np.float64_t * kernel_ptr = <np.float64_t*>np.PyArray_DATA(kernel)
        const size_t * kernel_shape = <size_t*>kernel.shape
    # convolveNd_c is declared nogil in the C header; release the GIL
    # around the call so threaded executors (e.g. dask, joblib threading
    # backend) can run multiple convolutions in parallel instead of being
    # serialised at the Python layer.
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
