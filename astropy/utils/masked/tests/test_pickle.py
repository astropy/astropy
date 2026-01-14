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
        Quantity([1, 2, 3], u.m),
        Angle([1, 2, 3], u.deg),
        Latitude([1, 2, 3], u.deg),
        Longitude([1, 2, 3], u.deg),
    ],
)
def test_masked_pickle(data):
    mask = [True, False, False]
    m = Masked(data, mask=mask)

    # Force creation of the info object (see #19142)
    assert m.info is not None

    m2 = pickle.loads(pickle.dumps(m))

    np.testing.assert_array_equal(m.unmasked, m2.unmasked)
    np.testing.assert_array_equal(m.mask, m2.mask)

    assert m.unit == m2.unit
    assert type(m) is type(m2)
    assert type(m.info) is type(m2.info)
