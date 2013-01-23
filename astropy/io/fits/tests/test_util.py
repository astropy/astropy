# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import with_statement


import os
import signal
import sys
import warnings

from ....tests.helper import pytest
from ..util import ignore_sigint

from . import FitsTestCase


class TestUtils(FitsTestCase):
    @pytest.mark.skipif("sys.platform.startswith('win')")
    def test_ignore_sigint(self):
        @ignore_sigint
        def test():
            with warnings.catch_warnings(record=True) as w:
                pid = os.getpid()
                os.kill(pid, signal.SIGINT)
                # One more time, for good measure
                os.kill(pid, signal.SIGINT)
                assert len(w) == 2
                assert (str(w[0].message) ==
                        'KeyboardInterrupt ignored until test is complete!')

        pytest.raises(KeyboardInterrupt, test)
