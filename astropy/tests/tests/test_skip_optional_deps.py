# Licensed under a 3-clause BSD style license - see LICENSE.rst
# This test doesn't actually have any optional dependencies, it should
# just be skipped by run_tests because it has the optional_deps decorator.

from ..helper import optional_deps
from ..helper import pytest


@optional_deps
def test_skip_optional_deps(pytestconfig):
    # this test was called from the command line and it should behave as if
    # astropy.test() has optional_deps=True
    if not hasattr(pytestconfig.option, 'optional_deps'):
        assert True

    # astropy.test() has optional_deps=False but we still got here somehow,
    # so fail with a helpful message
    elif not getattr(pytestconfig.option, 'optional_deps'):
        pytest.fail('@optional_deps was not skipped with optional_deps=False')

    # astropy.test() has optional_deps=True, so pass
    elif getattr(pytestconfig.option, 'optional_deps'):
        assert True
