# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pickle

import pytest
import numpy as np
from numpy.testing import assert_array_equal, assert_allclose

from astropy.nddata.nduncertainty import (StdDevUncertainty,
                             VarianceUncertainty,
                             InverseVariance,
                             NDUncertainty,
                             IncompatibleUncertaintiesException,
                             MissingDataAssociationException,
                             UnknownUncertainty)
from astropy.nddata.nddata import NDData
from astropy.nddata.compat import NDDataArray
from astropy.nddata.ccddata import CCDData
from astropy import units as u

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

    def _data_unit_to_uncertainty_unit(self, value):
        return None

    def _propagate_add(self, data, final_data):
        pass

    def _propagate_subtract(self, data, final_data):
        pass

    def _propagate_multiply(self, data, final_data):
        pass

    def _propagate_divide(self, data, final_data):
        pass


# Test the fake (added also StdDevUncertainty which should behave identical)

# the list of classes used for parametrization in tests below
uncertainty_types_to_be_tested = [
    FakeUncertainty,
    StdDevUncertainty,
    VarianceUncertainty,
    InverseVariance,
    UnknownUncertainty
]

uncertainty_types_with_conversion_support = (
    StdDevUncertainty, VarianceUncertainty, InverseVariance)
uncertainty_types_without_conversion_support = (
    FakeUncertainty, UnknownUncertainty)


@pytest.mark.parametrize(('UncertClass'), uncertainty_types_to_be_tested)
def test_init_fake_with_list(UncertClass):
    fake_uncert = UncertClass([1, 2, 3])
    assert_array_equal(fake_uncert.array, np.array([1, 2, 3]))
    # Copy makes no difference since casting a list to an np.ndarray always
    # makes a copy.
    # But let's give the uncertainty a unit too
    fake_uncert = UncertClass([1, 2, 3], unit=u.adu)
    assert_array_equal(fake_uncert.array, np.array([1, 2, 3]))
    assert fake_uncert.unit is u.adu


@pytest.mark.parametrize(('UncertClass'), uncertainty_types_to_be_tested)
def test_init_fake_with_ndarray(UncertClass):
    uncert = np.arange(100).reshape(10, 10)
    fake_uncert = UncertClass(uncert)
    # Numpy Arrays are copied by default
    assert_array_equal(fake_uncert.array, uncert)
    assert fake_uncert.array is not uncert
    # Now try it without copy
    fake_uncert = UncertClass(uncert, copy=False)
    assert fake_uncert.array is uncert
    # let's provide a unit
    fake_uncert = UncertClass(uncert, unit=u.adu)
    assert_array_equal(fake_uncert.array, uncert)
    assert fake_uncert.array is not uncert
    assert fake_uncert.unit is u.adu


@pytest.mark.parametrize(('UncertClass'), uncertainty_types_to_be_tested)
def test_init_fake_with_quantity(UncertClass):
    uncert = np.arange(10).reshape(2, 5) * u.adu
    fake_uncert = UncertClass(uncert)
    # Numpy Arrays are copied by default
    assert_array_equal(fake_uncert.array, uncert.value)
    assert fake_uncert.array is not uncert.value
    assert fake_uncert.unit is u.adu
    # Try without copy (should not work, quantity.value always returns a copy)
    fake_uncert = UncertClass(uncert, copy=False)
    assert fake_uncert.array is not uncert.value
    assert fake_uncert.unit is u.adu
    # Now try with an explicit unit parameter too
    fake_uncert = UncertClass(uncert, unit=u.m)
    assert_array_equal(fake_uncert.array, uncert.value)  # No conversion done
    assert fake_uncert.array is not uncert.value
    assert fake_uncert.unit is u.m  # It took the explicit one


@pytest.mark.parametrize(('UncertClass'), uncertainty_types_to_be_tested)
def test_init_fake_with_fake(UncertClass):
    uncert = np.arange(5).reshape(5, 1)
    fake_uncert1 = UncertClass(uncert)
    fake_uncert2 = UncertClass(fake_uncert1)
    assert_array_equal(fake_uncert2.array, uncert)
    assert fake_uncert2.array is not uncert
    # Without making copies
    fake_uncert1 = UncertClass(uncert, copy=False)
    fake_uncert2 = UncertClass(fake_uncert1, copy=False)
    assert_array_equal(fake_uncert2.array, fake_uncert1.array)
    assert fake_uncert2.array is fake_uncert1.array
    # With a unit
    uncert = np.arange(5).reshape(5, 1) * u.adu
    fake_uncert1 = UncertClass(uncert)
    fake_uncert2 = UncertClass(fake_uncert1)
    assert_array_equal(fake_uncert2.array, uncert.value)
    assert fake_uncert2.array is not uncert.value
    assert fake_uncert2.unit is u.adu
    # With a unit and an explicit unit-parameter
    fake_uncert2 = UncertClass(fake_uncert1, unit=u.cm)
    assert_array_equal(fake_uncert2.array, uncert.value)
    assert fake_uncert2.array is not uncert.value
    assert fake_uncert2.unit is u.cm


@pytest.mark.parametrize(('UncertClass'), uncertainty_types_to_be_tested)
def test_init_fake_with_somethingElse(UncertClass):
    # What about a dict?
    uncert = {'rdnoise': 2.9, 'gain': 0.6}
    fake_uncert = UncertClass(uncert)
    assert fake_uncert.array == uncert
    # We can pass a unit too but since we cannot do uncertainty propagation
    # the interpretation is up to the user
    fake_uncert = UncertClass(uncert, unit=u.s)
    assert fake_uncert.array == uncert
    assert fake_uncert.unit is u.s
    # So, now check what happens if copy is False
    fake_uncert = UncertClass(uncert, copy=False)
    assert fake_uncert.array == uncert
    assert id(fake_uncert) != id(uncert)
    # dicts cannot be referenced without copy
    # TODO : Find something that can be referenced without copy :-)


def test_init_fake_with_StdDevUncertainty():
    # Different instances of uncertainties are not directly convertible so this
    # should fail
    uncert = np.arange(5).reshape(5, 1)
    std_uncert = StdDevUncertainty(uncert)
    with pytest.raises(IncompatibleUncertaintiesException):
        FakeUncertainty(std_uncert)
    # Ok try it the other way around
    fake_uncert = FakeUncertainty(uncert)
    with pytest.raises(IncompatibleUncertaintiesException):
        StdDevUncertainty(fake_uncert)


def test_uncertainty_type():
    fake_uncert = FakeUncertainty([10, 2])
    assert fake_uncert.uncertainty_type == 'fake'
    std_uncert = StdDevUncertainty([10, 2])
    assert std_uncert.uncertainty_type == 'std'
    var_uncert = VarianceUncertainty([10, 2])
    assert var_uncert.uncertainty_type == 'var'
    ivar_uncert = InverseVariance([10, 2])
    assert ivar_uncert.uncertainty_type == 'ivar'


def test_uncertainty_correlated():
    fake_uncert = FakeUncertainty([10, 2])
    assert not fake_uncert.supports_correlated
    std_uncert = StdDevUncertainty([10, 2])
    assert std_uncert.supports_correlated


def test_for_leak_with_uncertainty():
    # Regression test for memory leak because of cyclic references between
    # NDData and uncertainty
    from collections import defaultdict
    from gc import get_objects

    def test_leak(func, specific_objects=None):
        """Function based on gc.get_objects to determine if any object or
        a specific object leaks.

        It requires a function to be given and if any objects survive the
        function scope it's considered a leak (so don't return anything).
        """
        before = defaultdict(int)
        for i in get_objects():
            before[type(i)] += 1

        func()

        after = defaultdict(int)
        for i in get_objects():
            after[type(i)] += 1

        if specific_objects is None:
            assert all(after[k] - before[k] == 0 for k in after)
        else:
            assert after[specific_objects] - before[specific_objects] == 0

    def non_leaker_nddata():
        # Without uncertainty there is no reason to assume that there is a
        # memory leak but test it nevertheless.
        NDData(np.ones(100))

    def leaker_nddata():
        # With uncertainty there was a memory leak!
        NDData(np.ones(100), uncertainty=StdDevUncertainty(np.ones(100)))

    test_leak(non_leaker_nddata, NDData)
    test_leak(leaker_nddata, NDData)

    # Same for NDDataArray:

    from astropy.nddata.compat import NDDataArray

    def non_leaker_nddataarray():
        NDDataArray(np.ones(100))

    def leaker_nddataarray():
        NDDataArray(np.ones(100), uncertainty=StdDevUncertainty(np.ones(100)))

    test_leak(non_leaker_nddataarray, NDDataArray)
    test_leak(leaker_nddataarray, NDDataArray)


def test_for_stolen_uncertainty():
    # Sharing uncertainties should not overwrite the parent_nddata attribute
    ndd1 = NDData(1, uncertainty=1)
    ndd2 = NDData(2, uncertainty=ndd1.uncertainty)
    # uncertainty.parent_nddata.data should be the original data!
    assert ndd1.uncertainty.parent_nddata.data == ndd1.data
    assert ndd2.uncertainty.parent_nddata.data == ndd2.data


def test_stddevuncertainty_pickle():
    uncertainty = StdDevUncertainty(np.ones(3), unit=u.m)
    uncertainty_restored = pickle.loads(pickle.dumps(uncertainty))
    np.testing.assert_array_equal(uncertainty.array, uncertainty_restored.array)
    assert uncertainty.unit == uncertainty_restored.unit
    with pytest.raises(MissingDataAssociationException):
        uncertainty_restored.parent_nddata


@pytest.mark.parametrize(('UncertClass'), uncertainty_types_to_be_tested)
def test_quantity(UncertClass):
    fake_uncert = UncertClass([1, 2, 3], unit=u.adu)
    assert isinstance(fake_uncert.quantity, u.Quantity)
    assert fake_uncert.quantity.unit.is_equivalent(u.adu)

    fake_uncert_nounit = UncertClass([1, 2, 3])
    assert isinstance(fake_uncert_nounit.quantity, u.Quantity)
    assert fake_uncert_nounit.quantity.unit.is_equivalent(u.dimensionless_unscaled)


@pytest.mark.parametrize(('UncertClass'),
                         [VarianceUncertainty,
                          StdDevUncertainty,
                          InverseVariance])
def test_setting_uncertainty_unit_results_in_unit_object(UncertClass):
    v = UncertClass([1, 1])
    v.unit = 'electron'
    assert isinstance(v.unit, u.UnitBase)


@pytest.mark.parametrize('NDClass', [NDData, NDDataArray, CCDData])
@pytest.mark.parametrize(('UncertClass'),
                         [VarianceUncertainty,
                          StdDevUncertainty,
                          InverseVariance])
def test_changing_unit_to_value_inconsistent_with_parent_fails(NDClass,
                                                               UncertClass):
    ndd1 = NDClass(1, unit='adu')
    v = UncertClass(1)
    # Sets the uncertainty unit to whatever makes sense with this data.
    ndd1.uncertainty = v

    with pytest.raises(u.UnitConversionError):
        # Nothing special about 15 except no one would ever use that unit
        v.unit = ndd1.unit ** 15


@pytest.mark.parametrize('NDClass', [NDData, NDDataArray, CCDData])
@pytest.mark.parametrize(('UncertClass, expected_unit'),
                         [(VarianceUncertainty, u.adu ** 2),
                          (StdDevUncertainty, u.adu),
                          (InverseVariance, 1 / u.adu ** 2)])
def test_assigning_uncertainty_to_parent_gives_correct_unit(NDClass,
                                                            UncertClass,
                                                            expected_unit):
    # Does assigning a unitless uncertainty to an NDData result in the
    # expected unit?
    ndd = NDClass([1, 1], unit=u.adu)
    v = UncertClass([1, 1])
    ndd.uncertainty = v
    assert v.unit == expected_unit


@pytest.mark.parametrize('NDClass', [NDData, NDDataArray, CCDData])
@pytest.mark.parametrize(('UncertClass, expected_unit'),
                         [(VarianceUncertainty, u.adu ** 2),
                          (StdDevUncertainty, u.adu),
                          (InverseVariance, 1 / u.adu ** 2)])
def test_assigning_uncertainty_with_unit_to_parent_with_unit(NDClass,
                                                             UncertClass,
                                                             expected_unit):
    # Does assigning an uncertainty with an appropriate unit to an NDData
    # with a unit work?
    ndd = NDClass([1, 1], unit=u.adu)
    v = UncertClass([1, 1], unit=expected_unit)
    ndd.uncertainty = v
    assert v.unit == expected_unit


@pytest.mark.parametrize('NDClass', [NDData, NDDataArray, CCDData])
@pytest.mark.parametrize(('UncertClass'),
                         [(VarianceUncertainty),
                          (StdDevUncertainty),
                          (InverseVariance)])
def test_assigning_uncertainty_with_bad_unit_to_parent_fails(NDClass,
                                                             UncertClass):
    # Does assigning an uncertainty with a non-matching unit to an NDData
    # with a unit work?
    ndd = NDClass([1, 1], unit=u.adu)
    # Set the unit to something inconsistent with ndd's unit
    v = UncertClass([1, 1], unit=u.second)
    with pytest.raises(u.UnitConversionError):
        ndd.uncertainty = v


@pytest.mark.parametrize('UncertClass', uncertainty_types_with_conversion_support)
def test_self_conversion_via_variance_supported(UncertClass):
    uncert = np.arange(1, 11).reshape(2, 5) * u.adu
    start_uncert = UncertClass(uncert)
    final_uncert = start_uncert.represent_as(UncertClass)
    assert_array_equal(start_uncert.array, final_uncert.array)
    assert start_uncert.unit == final_uncert.unit


@pytest.mark.parametrize(
    'UncertClass,to_variance_func',
    zip(uncertainty_types_with_conversion_support,
    (lambda x: x ** 2, lambda x: x, lambda x: 1 / x))
)
def test_conversion_to_from_variance_supported(UncertClass, to_variance_func):
    uncert = np.arange(1, 11).reshape(2, 5) * u.adu
    start_uncert = UncertClass(uncert)
    var_uncert = start_uncert.represent_as(VarianceUncertainty)
    final_uncert = var_uncert.represent_as(UncertClass)
    assert_allclose(to_variance_func(start_uncert.array), var_uncert.array)
    assert_array_equal(start_uncert.array, final_uncert.array)
    assert start_uncert.unit == final_uncert.unit


@pytest.mark.parametrize('UncertClass', uncertainty_types_without_conversion_support)
def test_self_conversion_via_variance_not_supported(UncertClass):
    uncert = np.arange(1, 11).reshape(2, 5) * u.adu
    start_uncert = UncertClass(uncert)
    with pytest.raises(TypeError):
        final_uncert = start_uncert.represent_as(UncertClass)
