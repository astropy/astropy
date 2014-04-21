# Licensed under a 3-clause BSD style license - see LICENSE.rst
# this test doesn't actually use any online data, it should just be skipped
# by run_tests because it has the remote_data decorator.
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from ..helper import remote_data, slow
from ..helper import pytest


@remote_data
def test_skip_remote_data(pytestconfig):
    # this test was called from the command line and it should behave as if
    # astropy.test() has remote_data=True
    if not hasattr(pytestconfig.option, 'remotedata'):
        assert True

    # astropy.test() has remote_data=False but we still got here somehow,
    # so fail with a helpful message
    elif not getattr(pytestconfig.option, 'remotedata'):
        pytest.fail('@remote_data was not skipped with remote_data=False')

    # astropy.test() has remote_data=True, so pass
    elif getattr(pytestconfig.option, 'remotedata'):
        assert True

@slow
def test_skip_slow(pytestconfig):
    # this test was called from the command line and it should behave as if
    # astropy.test() has slow=True
    if not hasattr(pytestconfig.option, 'slow'):
        assert True

    # astropy.test() has slow=False but we still got here somehow,
    # so fail with a helpful message
    elif not getattr(pytestconfig.option, 'slow'):
        pytest.fail('@slow was not skipped with slow=False')

    # astropy.test() has slow=True, so pass
    elif getattr(pytestconfig.option, 'slow'):
        assert True
