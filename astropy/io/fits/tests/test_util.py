import os
import signal
import warnings

from ....tests.helper import pytest
from ..util import ignore_sigint

from . import FitsTestCase


class TestUtils(FitsTestCase):
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
