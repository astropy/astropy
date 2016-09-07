# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import with_statement


import os
import signal

import numpy as np

from ....tests.helper import pytest, catch_warnings
from ..util import ignore_sigint
from .._numpy_hacks import realign_dtype

from . import FitsTestCase


class TestUtils(FitsTestCase):
    @pytest.mark.skipif("sys.platform.startswith('win')")
    def test_ignore_sigint(self):
        @ignore_sigint
        def test():
            with catch_warnings(UserWarning) as w:
                pid = os.getpid()
                os.kill(pid, signal.SIGINT)
                # One more time, for good measure
                os.kill(pid, signal.SIGINT)
            assert len(w) == 2
            assert (str(w[0].message) ==
                    'KeyboardInterrupt ignored until test is complete!')

        pytest.raises(KeyboardInterrupt, test)

    def test_realign_dtype(self):
        """
        Tests a few corner-cases for the realign_dtype hack.

        These are unlikely to come in practice given how this is currently
        used in astropy.io.fits, but nonetheless tests for bugs that were
        present in earlier versions of the function.
        """

        dt = np.dtype([('a', np.int32), ('b', np.int16)])
        dt2 = realign_dtype(dt, [0, 0])
        assert dt2.itemsize == 4

        dt2 = realign_dtype(dt, [0, 1])
        assert dt2.itemsize == 4

        dt2 = realign_dtype(dt, [1, 0])
        assert dt2.itemsize == 5

        dt = np.dtype([('a', np.float64), ('b', np.int8), ('c', np.int8)])
        dt2 = realign_dtype(dt, [0, 0, 0])
        assert dt2.itemsize == 8

        dt2 = realign_dtype(dt, [0, 0, 1])
        assert dt2.itemsize == 8

        dt2 = realign_dtype(dt, [0, 0, 27])
        assert dt2.itemsize == 28
