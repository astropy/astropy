# Licensed under a 3-clause BSD style license - see LICENSE.rst

# TEST_UNICODE_LITERALS

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np
from numpy.testing import assert_array_equal

from ..nduncertainty import (StdDevUncertainty, NDUncertainty,
                             IncompatibleUncertaintiesException)
from ...tests.helper import pytest
from ... import units as u

# Regarding setter tests:
# No need to test setters since the uncertainty is considered immutable after
# creation except of the parent_nddata attribute and this accepts just
# everything.
# Additionally they should be covered by NDData, NDArithmeticMixin which rely
# on it

# Regarding propagate, _convert_uncert, _propagate_* tests:
# They should be covered by NDArithmeticMixin since there is generally no need
# to test them without this mixin.

# Regarding __getitem__ tests:
# Should be covered by NDSlicingMixin.

# Regarding StdDevUncertainty tests:
# This subclass only overrides the methods for propagation so the same
# they should be covered in NDArithmeticMixin.

# Not really fake but the minimum an uncertainty has to override not to be
# abstract.
class FakeUncertainty(NDUncertainty):

    @property
    def uncertainty_type(self):
        return 'fake'

    def _propagate_add(self, data, final_data):
        pass

    def _propagate_subtract(self, data, final_data):
        pass

    def _propagate_multiply(self, data, final_data):
        pass

    def _propagate_divide(self, data, final_data):
        pass

# Test the fake

def test_init_fake_with_list():
    fake_uncert = FakeUncertainty([1,2,3])
    assert_array_equal(fake_uncert.array, np.array([1,2,3]))
    # Copy makes no difference since casting a list to an np.ndarray always
    # makes a copy.
    # But let's give the uncertainty a unit too
    fake_uncert = FakeUncertainty([1,2,3], unit=u.adu)
    assert_array_equal(fake_uncert.array, np.array([1,2,3]))
    assert fake_uncert.unit is u.adu


def test_init_fake_with_ndarray():
    uncert = np.arange(100).reshape(10,10)
    fake_uncert = FakeUncertainty(uncert)
    # Numpy Arrays are copied by default
    assert_array_equal(fake_uncert.array, uncert)
    assert fake_uncert.array is not uncert
    # Now try it without copy
    fake_uncert = FakeUncertainty(uncert, copy=False)
    assert fake_uncert.array is uncert
    # let's provide a unit
    fake_uncert = FakeUncertainty(uncert, unit=u.adu)
    assert_array_equal(fake_uncert.array, uncert)
    assert fake_uncert.array is not uncert
    assert fake_uncert.unit is u.adu


def test_init_fake_with_quantity():
    uncert = np.arange(10).reshape(2,5) * u.adu
    fake_uncert = FakeUncertainty(uncert)
    # Numpy Arrays are copied by default
    assert_array_equal(fake_uncert.array, uncert.value)
    assert fake_uncert.array is not uncert.value
    assert fake_uncert.unit is u.adu
    # Try without copy (should not work, quantity.value always returns a copy)
    fake_uncert = FakeUncertainty(uncert, copy=False)
    assert fake_uncert.array is not uncert.value
    assert fake_uncert.unit is u.adu
    # Now try with an explicit unit parameter too
    fake_uncert = FakeUncertainty(uncert, unit=u.m)
    assert_array_equal(fake_uncert.array, uncert.value) # No conversion is done
    assert fake_uncert.array is not uncert.value
    assert fake_uncert.unit is u.m # It took the explicit one


def test_init_fake_with_fake():
    uncert = np.arange(5).reshape(5,1)
    fake_uncert1 = FakeUncertainty(uncert)
    fake_uncert2 = FakeUncertainty(fake_uncert1)
    assert_array_equal(fake_uncert2.array, uncert)
    assert fake_uncert2.array is not uncert
    # Without making copies
    fake_uncert1 = FakeUncertainty(uncert, copy=False)
    fake_uncert2 = FakeUncertainty(fake_uncert1, copy=False)
    assert_array_equal(fake_uncert2.array, fake_uncert1.array)
    assert fake_uncert2.array is fake_uncert1.array
    # With a unit
    uncert = np.arange(5).reshape(5,1) * u.adu
    fake_uncert1 = FakeUncertainty(uncert)
    fake_uncert2 = FakeUncertainty(fake_uncert1)
    assert_array_equal(fake_uncert2.array, uncert.value)
    assert fake_uncert2.array is not uncert.value
    assert fake_uncert2.unit is u.adu
    # With a unit and an explicit unit-parameter
    fake_uncert2 = FakeUncertainty(fake_uncert1, unit=u.cm)
    assert_array_equal(fake_uncert2.array, uncert.value)
    assert fake_uncert2.array is not uncert.value
    assert fake_uncert2.unit is u.cm


def test_init_fake_with_somethingElse():
    # What about a dict?
    uncert = {'rdnoise': 2.9, 'gain': 0.6}
    fake_uncert = FakeUncertainty(uncert)
    assert fake_uncert.array == uncert
    # We can pass a unit too but since we cannot do uncertainty propagation
    # the interpretation is up to the user
    fake_uncert = FakeUncertainty(uncert, unit=u.s)
    assert fake_uncert.array == uncert
    assert fake_uncert.unit is u.s
    # So, now check what happens if copy is False
    fake_uncert = FakeUncertainty(uncert, copy=False)
    assert fake_uncert.array == uncert
    assert id(fake_uncert) != id(uncert)
    # dicts cannot be referenced without copy
    # TODO : Find something that can be referenced without copy :-)


def test_init_fake_with_StdDevUncertainty():
    # Different instances of uncertainties are not directly convertable so this
    # should fail
    uncert = np.arange(5).reshape(5,1)
    std_uncert = StdDevUncertainty(uncert)
    with pytest.raises(IncompatibleUncertaintiesException):
        FakeUncertainty(std_uncert)
    # Ok try it the other way around
    fake_uncert = FakeUncertainty(uncert)
    with pytest.raises(IncompatibleUncertaintiesException):
        StdDevUncertainty(fake_uncert)


def test_uncertainty_type():
    fake_uncert = FakeUncertainty([10,2])
    assert fake_uncert.uncertainty_type == 'fake'
    std_uncert = StdDevUncertainty([10,2])
    assert std_uncert.uncertainty_type == 'std'


def test_uncertainty_correlated():
    fake_uncert = FakeUncertainty([10,2])
    assert not fake_uncert.supports_correlated
    std_uncert = StdDevUncertainty([10,2])
    assert std_uncert.supports_correlated
