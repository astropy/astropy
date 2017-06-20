# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Tests for blackbody functions."""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

# LOCAL
from ..blackbody import blackbody_nu, blackbody_lambda
from ... import units as u
from ...tests.helper import catch_warnings
from ...utils.exceptions import AstropyDeprecationWarning

__doctest_skip__ = ['*']


def test_deprecated_blackbodies():
    with catch_warnings(AstropyDeprecationWarning) as w:
        blackbody_nu(5000 * u.AA, 6000 * u.K)
    assert len(w) == 1

    with catch_warnings(AstropyDeprecationWarning) as w:
        blackbody_lambda(5000 * u.AA, 6000 * u.K)
    assert len(w) == 1
