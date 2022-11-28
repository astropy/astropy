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
    convolveNd_c(
        <np.float64_t*>np.PyArray_DATA(result),
        <np.float64_t*>np.PyArray_DATA(array_to_convolve),
        array_to_convolve.ndim,
        <size_t*>array_to_convolve.shape,
        <np.float64_t*>np.PyArray_DATA(kernel),
        <size_t*>kernel.shape,
        nan_interpolate,
        embed_result_within_padded_region,
        n_threads,
    )
