import time

import numpy as np

from astropy.utils.misc import NumpyRNGContext


def func(i):
    """An identity function that jitters its execution time by a
    pseudo-random amount.

    FIXME: This function should be defined in test_console.py, but Astropy's
    `python setup.py test` interacts strangely with Python's `multiprocessing`
    module. I was getting a mysterious PicklingError until I moved this
    function into a separate module. (It worked fine in a standalone pytest
    script.)"""
    with NumpyRNGContext(i):
        time.sleep(np.random.uniform(0, 0.01))
    return i
