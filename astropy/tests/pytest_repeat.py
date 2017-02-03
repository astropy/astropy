# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
These plugins modify the behavior of py.test and are meant to be imported
into conftest.py in the root directory.
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from ..extern.six.moves import range


def pytest_addoption(parser):

    parser.addoption('--repeat', action='store',
                     help='Number of times to repeat each test')


def pytest_generate_tests(metafunc):

    # If the repeat option is set, we add a fixture for the repeat count and
    # parametrize the tests over the repeats. Solution adapted from:
    # http://stackoverflow.com/q/21764473/180783

    if metafunc.config.option.repeat is not None:
        count = int(metafunc.config.option.repeat)
        metafunc.fixturenames.append('tmp_ct')
        metafunc.parametrize('tmp_ct', range(count))
