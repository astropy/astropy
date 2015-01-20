# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Tests for miscellaneous functionality in the `funcs` module
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np

from ...tests.helper import pytest

from ... import units as u
from ...time import Time


def test_sun():
    """
    Test that `get_sun` works and it behaves roughly as it should (in GCRS)
    """
    from ..funcs import get_sun

    northern_summer_solstice = Time('2010-6-21')
    northern_winter_solstice = Time('2010-12-21')
    equinox_1 = Time('2010-3-21')
    equinox_2 = Time('2010-9-21')

    gcrs1 = get_sun(equinox_1)
    assert np.abs(gcrs1.dec.deg) < 1

    gcrs2 = get_sun(Time([northern_summer_solstice, equinox_2, northern_winter_solstice]))
    assert np.all(np.abs(gcrs2.dec - [23.5, 0, -23.5]*u.deg) < 1*u.deg)

def test_concatenate():
    from .. import FK5, SkyCoord
    from ..funcs import concatenate

    fk5 = FK5(1*u.deg, 2*u.deg)
    sc = SkyCoord(3*u.deg, 4*u.deg, frame='fk5')

    res = concatenate([fk5, sc])
    np.testing.assert_allclose(res.ra, [1, 3]*u.deg)
    np.testing.assert_allclose(res.dec, [2, 4]*u.deg)

    with pytest.raises(TypeError):
        concatenate(fk5)

    with pytest.raises(TypeError):
        concatenate(1*u.deg)
