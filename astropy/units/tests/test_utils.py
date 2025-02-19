# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Test utilities for `astropy.units`.
"""

import numpy as np
import pytest

from astropy.units.quantity import Quantity
from astropy.units.utils import quantity_asanyarray


def test_quantity_asanyarray():
    array_of_quantities = [Quantity(1), Quantity(2), Quantity(3)]
    quantity_array = quantity_asanyarray(array_of_quantities)
    assert isinstance(quantity_array, Quantity)
    assert np.issubdtype(quantity_array.dtype, np.inexact)

    array_of_integers = [1, 2, 3]
    np_array = quantity_asanyarray(array_of_integers)
    assert isinstance(np_array, np.ndarray)
    assert np.issubdtype(np_array.dtype, np.integer)

    np_array = quantity_asanyarray(array_of_integers, dtype=np.inexact)
    assert np.issubdtype(np_array.dtype, np.inexact)


def test_composite_unit_definition():
    # Regression test for #17781 - error message was not informative
    with pytest.raises(
        TypeError, match=r"^'km_per_h' must be defined with 'def_unit\(\)'$"
    ):
        from astropy.units.tests.data import bad_module  # noqa: F401
