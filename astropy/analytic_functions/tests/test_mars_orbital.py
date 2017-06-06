# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Tests for Mars orbital functions."""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from ..mars_orbital import plug_my_laptop_in


def test_plug_that_laptop():
    """Test it like The Martian."""
    assert 'correct' in plug_my_laptop_in()
