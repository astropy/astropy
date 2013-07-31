# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import print_function, division

from ...tests.helper import pytest

import numpy as np
from numpy import testing as npt

from ... import units as u


def test_angle_arrays():
    """
    Test arrays values with Angle objects.
    """
    from .. import Angle

    # Tests incomplete
    a1 = Angle([0, 45, 90, 180, 270, 360, 720.], unit=u.degree)
    npt.assert_almost_equal([0., 45., 90., 180., 270., 360., 360.], a1.value)

    a2 = Angle(np.array([-90, -45, 0, 45, 90, 180, 270, 360]), unit=u.degree,
               bounds=(0, 360))
    npt.assert_almost_equal([270., 315., 0., 45., 90., 180., 270., 360.],
                            a2.value)

    a3 = Angle(["12 degrees", "3 hours", "5 deg", "4rad"])
    npt.assert_almost_equal([0.20943951, 0.78539816, 0.08726646, 4], a3.value)
    assert a3.unit == u.radian

    a4 = Angle(["12 degrees", "3 hours", "5 deg", "4rad"], u.degree)
    npt.assert_almost_equal(a4.radian, a3.value)
    assert a4.unit == u.degree

    a5 = Angle([0, 45, 90, 180, 270, 360], unit=u.degree, bounds=(0, 360))
    a6 = a5.sum()
    npt.assert_almost_equal(a6.value, 225.0)
    assert a6.unit is u.degree
    assert a6.bounds == (Angle(0, u.degree), Angle(360, u.degree))

    with pytest.raises(TypeError):
        # Arrays of Angle objects are not supported -- that's really
        # tricky to do correctly, if at all, due to the possibility of
        # nesting.
        a7 = Angle([a1, a2, a3], unit=u.degree)
