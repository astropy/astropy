# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np

from ...tests.helper import pytest
from ...table import Column, Table
from ... import units as u


class TestQuantityColumn():
    """Test that quantities are properly represented in Columns"""

    def test_examples_from_eteq(self):
        """From make tables work consistently with quantites [#2000]"""
        t = Table()
        t['length'] = np.arange(10.) * u.km
        t['time'] = np.arange(1.,11.) * u.s
        v = t['length']/t['time']
        assert isinstance(v, u.Quantity)
        assert v.unit == u.Unit('km/s')
        check = t['length']*u.km / (t['time']*u.second)
        assert isinstance(check, u.Quantity)
        assert check.unit == u.Unit('km2/s2')
        v2 = u.Quantity(t['length'])/u.Quantity(t['time'])
        assert isinstance(v2, u.Quantity)
        assert v2.unit == u.Unit('km/s')
