# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Test importability of common Masked classes."""

import pytest

from astropy.coordinates.angles.core import Angle, Latitude, Longitude
from astropy.units.quantity import Quantity
from astropy.utils.masked import core


@pytest.mark.parametrize("base_class", [Quantity, Angle, Latitude, Longitude])
def test_importable(base_class):
    assert getattr(core, f"Masked{base_class.__name__}") is core.Masked(base_class)
