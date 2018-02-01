# Licensed under a 3-clause BSD style license - see LICENSE.rst

import warnings
from multiprocessing import cpu_count
from ..openmp_enabled import is_openmp_enabled
from .exceptions import AstropyUserWarning

def handle_n_threads_usage(n_threads):
    '''
    Helper function to check n_thread usage and correct where possible.
    '''

    if not isinstance(n_threads, int) or n_threads < 0:
        raise ValueError("n_threads must be a positive integer")

    total_cpus = cpu_count()
    if n_threads > 1 and not is_openmp_enabled():
        warnings.warn("n_threads={0} used yet Astropy was NOT built "
                      "with OpenMP support. "
                      "Running single threaded only.".format(n_threads),
                      AstropyUserWarning)
        n_threads = 1
    elif n_threads > total_cpus:
        warnings.warn("n_threads is greater than the total number "
                      "of CPUs: {0}. Over commiting threads is unlikely "
                      "to boost performance and may lock-up your machine. "
                      "However, if mainly IO intensive the rule of thumb "
                      "is to use 1.5*total. NOTE: Check hyperthreading status, "
                      "if ON, using the total number of CPUs will already be over commiting.".format(total_cpus),
                      AstropyUserWarning)

    if n_threads == 0:
        n_threads = total_cpus

    return n_threads
