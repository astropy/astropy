# Licensed under a 3-clause BSD style license - see LICENSE.rst
#cython: language_level=3

import numpy as np

cimport cython
cimport numpy as np
from libc.math cimport sqrt
from libc.stdlib cimport free, malloc

np.import_array()

DTYPE = np.float64
ctypedef np.float64_t DTYPE_t

IDTYPE = np.int64
ctypedef np.int64_t IDTYPE_t

cdef extern int run_bls (
    int N,                   # Length of the time array
    double* t,               # The list of timestamps
    double* y,               # The y measured at ``t``
    double* ivar,            # The inverse variance of the y array

    int n_periods,           #
    double* periods,         # The period to test in units of ``t``

    int n_durations,         # Length of the durations array
    double* durations,       # The durations to test in units of ``bin_duration``
    int oversample,          # The number of ``bin_duration`` bins in the maximum duration

    int obj_flag,            # A flag indicating the periodogram type
                             # 0 - depth signal-to-noise
                             # 1 - log likelihood

    # Outputs
    double* best_objective,  # The value of the periodogram at maximum
    double* best_depth,      # The estimated depth at maximum
    double* best_depth_std,  # The uncertainty on ``best_depth``
    double* best_duration,   # The best fitting duration in units of ``t``
    double* best_phase,      # The phase of the mid-transit time in units of
                             # ``t``
    double* best_depth_snr,  # The signal-to-noise ratio of the depth estimate
    double* best_log_like    # The log likelihood at maximum
) nogil


@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
def bls_impl(
    np.ndarray[DTYPE_t, mode='c'] t_array,
    np.ndarray[DTYPE_t, mode='c'] y_array,
    np.ndarray[DTYPE_t, mode='c'] ivar_array,
    np.ndarray[DTYPE_t, mode='c'] period_array,
    np.ndarray[DTYPE_t, mode='c'] duration_array,
    int oversample,
    int obj_flag
):

    cdef np.ndarray[DTYPE_t, mode='c'] out_objective = np.empty_like(period_array, dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, mode='c'] out_depth     = np.empty_like(period_array, dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, mode='c'] out_depth_err = np.empty_like(period_array, dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, mode='c'] out_duration  = np.empty_like(period_array, dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, mode='c'] out_phase     = np.empty_like(period_array, dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, mode='c'] out_depth_snr = np.empty_like(period_array, dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, mode='c'] out_log_like  = np.empty_like(period_array, dtype=DTYPE)
    cdef int flag, N = len(t_array), n_periods = len(period_array), n_durations = len(duration_array)

    with nogil:
        flag = run_bls(
            N,
            <double*>t_array.data,
            <double*>y_array.data,
            <double*>ivar_array.data,
            n_periods,
            <double*>period_array.data,
            n_durations,
            <double*>duration_array.data,
            oversample,
            obj_flag,
            <double*>out_objective.data,
            <double*>out_depth.data,
            <double*>out_depth_err.data,
            <double*>out_duration.data,
            <double*>out_phase.data,
            <double*>out_depth_snr.data,
            <double*>out_log_like.data
        )

    if flag < 0:
        raise MemoryError()
    if flag > 0:
        raise ValueError("Invalid inputs for period and/or duration")

    return (out_objective, out_depth, out_depth_err, out_duration, out_phase,
            out_depth_snr, out_log_like)
