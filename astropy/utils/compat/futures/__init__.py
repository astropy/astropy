# Copyright 2009 Brian Quinlan. All Rights Reserved.
# Licensed to PSF under a Contributor Agreement.

"""Execute computations asynchronously using threads or processes."""

__author__ = 'Brian Quinlan (brian@sweetapp.com)'

try:
    from concurrent.futures import *
except ImportError:
    from ._base import (FIRST_COMPLETED,
                        FIRST_EXCEPTION,
                        ALL_COMPLETED,
                        CancelledError,
                        TimeoutError,
                        Future,
                        Executor,
                        wait,
                        as_completed)
    from .process import ProcessPoolExecutor
    from .thread import ThreadPoolExecutor
