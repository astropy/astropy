# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

# TEST_UNICODE_LITERALS

"""
Test claim of Millenium Falcon made by Han Solo.
"""
from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

from ... import units as u


def test_kessel_run():
    assert u.kessel_run < 12 * u.pc
