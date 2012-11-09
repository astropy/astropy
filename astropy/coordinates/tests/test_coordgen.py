# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import print_function

from ...tests.helper import pytest

from .. import Coordinates, ICRSCoordinates, GalacticCoordinates
from .. import HorizontalCoordinates, Distance
from ... import units as u


def test_coordinate_generator():
    """
    Tests to make sure the Coordinate generator class works correctly
    """
    c1 = Coordinates(ra=1, dec=2, unit=u.radian)
    assert isinstance(c1, ICRSCoordinates)
    c2 = Coordinates(l=1, b=2, unit=u.radian)
    assert isinstance(c2, GalacticCoordinates)
    c3 = Coordinates(az=1, el=2, unit=u.radian)
    assert isinstance(c3, HorizontalCoordinates)

    c4 = Coordinates(ra=1, dec=2, unit=u.radian, distance=Distance(1, u.kpc))

    with pytest.raises(TypeError):
        c = Coordinates(ra=1, dec=2,unit=u.radian, extra='foo')
