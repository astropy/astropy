# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Tests for miscellaneous functionality in the `funcs` module
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import pytest
import numpy as np
from numpy import testing as npt

from ...extern import six

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


def test_constellations():
    from .. import ICRS, FK5, SkyCoord
    from ..funcs import get_constellation

    inuma = ICRS(9*u.hour, 65*u.deg)
    res = get_constellation(inuma)
    res_short = get_constellation(inuma, short_name=True)
    assert res == 'Ursa Major'
    assert res_short == 'UMa'
    assert isinstance(res, six.string_types) or getattr(res, 'shape', None) == tuple()

    # these are taken from the ReadMe for Roman 1987
    ras = [9, 23.5, 5.12, 9.4555, 12.8888, 15.6687, 19, 6.2222]
    decs = [65, -20, 9.12, -19.9, 22, -12.1234, -40, -81.1234]
    shortnames = ['UMa', 'Aqr', 'Ori', 'Hya', 'Com', 'Lib', 'CrA', 'Men']

    testcoos = FK5(ras*u.hour, decs*u.deg, equinox='B1950')
    npt.assert_equal(get_constellation(testcoos, short_name=True), shortnames)

    # test on a SkyCoord, *and* test Boötes, which is special in that it has a
    # non-ASCII character
    bootest = SkyCoord(15*u.hour, 30*u.deg, frame='icrs')
    boores = get_constellation(bootest)
    assert boores == u'Boötes'
    assert isinstance(boores, six.string_types) or getattr(boores, 'shape', None) == tuple()
