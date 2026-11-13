# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pickle

import numpy as np
import pytest

import astropy.units as u
from astropy.coordinates import Angle, Latitude, Longitude
from astropy.units import Quantity
from astropy.utils.masked import Masked


@pytest.mark.parametrize(
    "data",
    [
        np.array([1, 2, 3]),
        np.array([(1, 2), (3, 4), (5, 6)], dtype=[("a", "i4"), ("b", "f8")]),
        np.rec.array([(1, 2), (3, 4), (5, 6)], names=["a", "b"]),
        Quantity([1, 2, 3], u.m),
        Angle([1, 2, 3], u.deg),
        Latitude([1, 2, 3], u.deg),
        Longitude([1, 2, 3], u.deg, wrap_angle=180 * u.deg),
    ],
)
@pytest.mark.parametrize("protocol", range(pickle.HIGHEST_PROTOCOL + 1))
def test_masked_pickle(data, protocol):
    mask = [True, False, False]
    m = Masked(data, mask=mask)

    # Force creation of the info object (see #19142)
    m.info.name = "masked_data"

    m2 = pickle.loads(pickle.dumps(m, protocol=protocol))

    np.testing.assert_array_equal(m.unmasked, m2.unmasked)
    np.testing.assert_array_equal(m.mask, m2.mask)

    assert getattr(m, "unit", None) == getattr(m2, "unit", None)
    if isinstance(m, Longitude):
        assert m.wrap_angle == m2.wrap_angle
    assert type(m) is type(m2)
    assert type(m.info) is type(m2.info)
    assert m.info.name == m2.info.name
