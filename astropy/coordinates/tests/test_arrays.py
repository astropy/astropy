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
    with pytest.raises(NotImplementedError):
        a1 = Angle([0, 45, 90, 180, 270, 360], unit=u.degree)

    with pytest.raises(NotImplementedError):
        a2 = Angle(np.array([0, 45, 90, 180, 270, 360]), unit=u.degree)

    with pytest.raises(NotImplementedError):
        a3 = Angle(["12 degrees", "3 hours", "5 deg", "4rad"])
