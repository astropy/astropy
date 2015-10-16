# Licensed under a 3-clause BSD style license - see LICENSE.rst

# TEST_UNICODE_LITERALS

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np
from numpy.testing import assert_array_equal, assert_array_almost_equal

from ... import NDData, NDArithmeticMixin
from ...nduncertainty import NDUncertainty, StdDevUncertainty
from ....units import UnitsError, Quantity
from ....tests.helper import pytest
from .... import units as u

# TODO: Tests with UnknownUncertainty


# Just add the Mixin to NDData
class NDDataArithmetic(NDArithmeticMixin, NDData):

    pass


# Test with Data covers:
# scalars, 1D, 2D and 3D
# broadcasting between them
@pytest.mark.parametrize(('data1', 'data2'), [
                         (np.array(5), np.array(10)),
                         (np.array(5), np.arange(10)),
                         (np.array(5), np.arange(10).reshape(2, 5)),
                         (np.arange(10), np.ones(10) * 2),
                         (np.arange(10), np.ones((10, 10)) * 2),
                         (np.arange(10).reshape(2, 5), np.ones((2, 5)) * 3),
                         (np.arange(1000).reshape(20, 5, 10),
                          np.ones((20, 5, 10)) * 3)
                         ])
def test_arithmetics_data(data1, data2):

    nd1 = NDDataArithmetic(data1)
    nd2 = NDDataArithmetic(data2)

    # Addition
    nd3 = nd1.add(nd2)
    assert_array_equal(data1+data2, nd3.data)
    # Subtraction
    nd4 = nd1.subtract(nd2)
    assert_array_equal(data1-data2, nd4.data)
    # Multiplication
    nd5 = nd1.multiply(nd2)
    assert_array_equal(data1*data2, nd5.data)
    # Division
    nd6 = nd1.divide(nd2)
    assert_array_equal(data1/data2, nd6.data)
    for nd in [nd3, nd4, nd5, nd6]:
        # Check that broadcasting worked as expected
        if data1.ndim > data2.ndim:
            assert data1.shape == nd.data.shape
        else:
            assert data2.shape == nd.data.shape
        # Check all other attributes are not set
        assert nd.unit is None
        assert nd.uncertainty is None
        assert nd.mask is None
        assert len(nd.meta) == 0
        assert nd.wcs is None


# Invalid arithmetic operations for data covering:
# not broadcastable data
def test_arithmetics_data_invalid():
    nd1 = NDDataArithmetic([1, 2, 3])
    nd2 = NDDataArithmetic([1, 2])
    with pytest.raises(ValueError):
        nd1.add(nd2)


# Test with Data and unit and covers:
# identical units (even dimensionless unscaled vs. no unit),
# equivalent units (such as meter and kilometer)
# equivalent composite units (such as m/s and km/h)
@pytest.mark.parametrize(('data1', 'data2'), [
    (np.array(5) * u.s, np.array(10) * u.s),
    (np.array(5) * u.s, np.arange(10) * u.h),
    (np.array(5) * u.s, np.arange(10).reshape(2, 5) * u.min),
    (np.arange(10) * u.m / u.s, np.ones(10) * 2 * u.km / u.s),
    (np.arange(10) * u.m / u.s, np.ones((10, 10)) * 2 * u.m / u.h),
    (np.arange(10).reshape(2, 5) * u.m / u.s,
     np.ones((2, 5)) * 3 * u.km / u.h),
    (np.arange(1000).reshape(20, 5, 10),
     np.ones((20, 5, 10)) * 3 * u.dimensionless_unscaled),
    (np.array(5), np.array(10) * u.s / u.h),
    ])
def test_arithmetics_data_unit_identical(data1, data2):

    nd1 = NDDataArithmetic(data1)
    nd2 = NDDataArithmetic(data2)

    # Addition
    nd3 = nd1.add(nd2)
    ref = data1 + data2
    ref_unit, ref_data = ref.unit, ref.value
    assert_array_equal(ref_data, nd3.data)
    assert nd3.unit == ref_unit
    # Subtraction
    nd4 = nd1.subtract(nd2)
    ref = data1 - data2
    ref_unit, ref_data = ref.unit, ref.value
    assert_array_equal(ref_data, nd4.data)
    assert nd4.unit == ref_unit
    # Multiplication
    nd5 = nd1.multiply(nd2)
    ref = data1 * data2
    ref_unit, ref_data = ref.unit, ref.value
    assert_array_equal(ref_data, nd5.data)
    assert nd5.unit == ref_unit
    # Division
    nd6 = nd1.divide(nd2)
    ref = data1 / data2
    ref_unit, ref_data = ref.unit, ref.value
    assert_array_equal(ref_data, nd6.data)
    assert nd6.unit == ref_unit
    for nd in [nd3, nd4, nd5, nd6]:
        # Check that broadcasting worked as expected
        if data1.ndim > data2.ndim:
            assert data1.shape == nd.data.shape
        else:
            assert data2.shape == nd.data.shape
        # Check all other attributes are not set
        assert nd.uncertainty is None
        assert nd.mask is None
        assert len(nd.meta) == 0
        assert nd.wcs is None

# Test with Data and unit and covers:
# not identical not convertable units
# one with unit (which is not dimensionless) and one without
@pytest.mark.parametrize(('data1','data2'), [
    (np.array(5)*u.s, np.array(10)*u.m),
    (np.array(5)*u.Mpc, np.array(10)*u.km / u.s),
    (np.array(5)*u.Mpc, np.array(10)),
    (np.array(5), np.array(10)*u.s),
    ])
def test_arithmetics_data_unit_NOT_identical(data1, data2):

    nd1 = NDDataArithmetic(data1)
    nd2 = NDDataArithmetic(data2)

    # Addition should not be possible
    with pytest.raises(UnitsError):
        nd1.add(nd2)
    # Subtraction should not be possible
    with pytest.raises(UnitsError):
        nd1.subtract(nd2)
    # Multiplication is possible
    nd3 = nd1.multiply(nd2)
    ref = data1 * data2
    ref_unit, ref_data = ref.unit, ref.value
    assert_array_equal(ref_data, nd3.data)
    assert nd3.unit == ref_unit
    # Division is possible
    nd4 = nd1.divide(nd2)
    ref = data1 / data2
    ref_unit, ref_data = ref.unit, ref.value
    assert_array_equal(ref_data, nd4.data)
    assert nd4.unit == ref_unit
    for nd in [nd3, nd4]:
        # Check all other attributes are not set
        assert nd.uncertainty is None
        assert nd.mask is None
        assert len(nd.meta) == 0
        assert nd.wcs is None


# Tests with wcs (not very senseable because there is no operation between them
# covering:
# both set and identical/not identical
# one set
# None set
@pytest.mark.parametrize(('wcs1','wcs2'), [
    (None, None),
    (None, 5),
    (5, None),
    (5, 5),
    (7, 5),
    ])
def test_arithmetics_data_wcs(wcs1, wcs2):

    nd1 = NDDataArithmetic(1, wcs=wcs1)
    nd2 = NDDataArithmetic(1, wcs=wcs2)

    if wcs1 is None and wcs2 is None:
        ref_wcs = None
    elif wcs1 is None:
        ref_wcs = wcs2
    elif wcs2 is None:
        ref_wcs = wcs1
    else:
        ref_wcs = wcs1

    # Addition
    nd3 = nd1.add(nd2)
    assert ref_wcs == nd3.wcs
    # Subtraction
    nd4 = nd1.subtract(nd2)
    assert ref_wcs == nd3.wcs
    # Multiplication
    nd5 = nd1.multiply(nd2)
    assert ref_wcs == nd3.wcs
    # Division
    nd6 = nd1.divide(nd2)
    assert ref_wcs == nd3.wcs
    for nd in [nd3, nd4, nd5, nd6]:
        # Check all other attributes are not set
        assert nd.unit is None
        assert nd.uncertainty is None
        assert len(nd.meta) == 0
        assert nd.mask is None

# Masks are completely seperated in the NDArithmetics from the data so we need
# no correlated tests but covering:
# masks 1D, 2D and mixed cases with broadcasting
@pytest.mark.parametrize(('mask1','mask2'), [
    (None, None),
    (None, False),
    (True, None),
    (False, False),
    (True, False),
    (False, True),
    (True, True),
    (np.array(False), np.array(True)),
    (np.array(False), np.array([0,1,0,1,1], dtype=np.bool_)),
    (np.array(True), np.array([[0,1,0,1,1],[1,1,0,1,1]], dtype=np.bool_)),
    (np.array([0,1,0,1,1], dtype=np.bool_), np.array([1,1,0,0,1], dtype=np.bool_)),
    (np.array([0,1,0,1,1], dtype=np.bool_), np.array([[0,1,0,1,1],[1,0,0,1,1]], dtype=np.bool_)),
    (np.array([[0,1,0,1,1],[1,0,0,1,1]], dtype=np.bool_), np.array([[0,1,0,1,1],[1,1,0,1,1]], dtype=np.bool_)),
    ])
def test_arithmetics_data_masks(mask1, mask2):

    nd1 = NDDataArithmetic(1, mask=mask1)
    nd2 = NDDataArithmetic(1, mask=mask2)

    if mask1 is None and mask2 is None:
        ref_mask = None
    elif mask1 is None:
        ref_mask = mask2
    elif mask2 is None:
        ref_mask = mask1
    else:
        ref_mask = mask1 | mask2

    # Addition
    nd3 = nd1.add(nd2)
    assert_array_equal(ref_mask, nd3.mask)
    # Subtraction
    nd4 = nd1.subtract(nd2)
    assert_array_equal(ref_mask, nd4.mask)
    # Multiplication
    nd5 = nd1.multiply(nd2)
    assert_array_equal(ref_mask, nd5.mask)
    # Division
    nd6 = nd1.divide(nd2)
    assert_array_equal(ref_mask, nd6.mask)
    for nd in [nd3, nd4, nd5, nd6]:
        # Check all other attributes are not set
        assert nd.unit is None
        assert nd.uncertainty is None
        assert len(nd.meta) == 0
        assert nd.wcs is None

# Masks can only be used in that way if they are booleans or np.ndarrays of
# booleans so everything else should result in an TypeError
# Covered cases:
# both: strings, integers, floats, classes
# only one is such a type (which just copies it no exception)
@pytest.mark.parametrize(('mask1','mask2'), [
    ('String', 'String'),
    (100, 100),
    (2.7, 2.7),
    ])
def test_arithmetics_data_masks_invalid(mask1, mask2):

    # Only one mask is not the expected type should work

    nd1 = NDDataArithmetic(1, mask=mask1)
    nd2 = NDDataArithmetic(1, mask=None)

    assert nd1.add(nd2).mask == mask1
    assert nd1.subtract(nd2).mask == mask1
    assert nd1.multiply(nd2).mask == mask1
    assert nd1.divide(nd2).mask == mask1

    nd1 = NDDataArithmetic(1, mask=None)
    nd2 = NDDataArithmetic(1, mask=mask2)

    assert nd1.add(nd2).mask == mask2
    assert nd1.subtract(nd2).mask == mask2
    assert nd1.multiply(nd2).mask == mask2
    assert nd1.divide(nd2).mask == mask2

    # Operations should not be possible if both are of some unexpected type

    nd1 = NDDataArithmetic(1, mask=mask1)
    nd2 = NDDataArithmetic(1, mask=mask2)

    with pytest.raises(TypeError):
        nd1.add(nd2)
    with pytest.raises(TypeError):
        nd1.multiply(nd2)
    with pytest.raises(TypeError):
        nd1.subtract(nd2)
    with pytest.raises(TypeError):
        nd1.divide(nd2)

# One additional case which can not be easily incorporated in the test above
# what happens if the masks are numpy ndarrays are not broadcastable
def test_arithmetics_data_masks_invalid_2():

    nd1 = NDDataArithmetic(1, mask=np.array([1,0], dtype=np.bool_))
    nd2 = NDDataArithmetic(1, mask=np.array([1,0,1], dtype=np.bool_))

    with pytest.raises(ValueError):
        nd1.add(nd2)
    with pytest.raises(ValueError):
        nd1.multiply(nd2)
    with pytest.raises(ValueError):
        nd1.subtract(nd2)
    with pytest.raises(ValueError):
        nd1.divide(nd2)

# Covering:
# both have meta with and without meta_keyword_operations
def test_arithmetics_meta_both():
    meta1 = {'exposure': 100, 'dummy': 5}
    meta2 = {'exposure': 50, 'dummy': 10}
    nd1 = NDDataArithmetic(1, meta=meta1)
    nd2 = NDDataArithmetic(1, meta=meta2)

    # No keywords specified to be part of arithmetics
    nd3 = nd1.add(nd2)
    assert nd3.meta['exposure'] == 100
    assert nd3.meta['dummy'] == 5

    # keyword IS specified to be part of arithmetics
    nd4 = nd1.add(nd2, meta_kwds_operate=['exposure'])
    assert nd4.meta['exposure'] == 150
    assert nd4.meta['dummy'] == 5

    nd4 = nd1.subtract(nd2, meta_kwds_operate=['exposure'])
    assert nd4.meta['exposure'] == 50
    assert nd4.meta['dummy'] == 5

    nd4 = nd1.multiply(nd2, meta_kwds_operate=['exposure'])
    assert nd4.meta['exposure'] == 5000
    assert nd4.meta['dummy'] == 5

    nd4 = nd1.divide(nd2, meta_kwds_operate=['exposure'])
    assert nd4.meta['exposure'] == 2
    assert nd4.meta['dummy'] == 5

# Covering:
# both have meta but the keyword is missing in one of them
def test_arithmetics_meta_both_invalid():
    meta1 = {'exposure': 100, 'dummy': 5}
    meta2 = {'exposures': 50, 'dummy': 10}
    nd1 = NDDataArithmetic(1, meta=meta1)
    nd2 = NDDataArithmetic(1, meta=meta2)

    # keyword IS specified to be part of arithmetics but exists only in one
    with pytest.raises(KeyError):
        nd1.add(nd2, meta_kwds_operate=['exposure'])
    with pytest.raises(KeyError):
        nd1.subtract(nd2, meta_kwds_operate=['exposure'])
    with pytest.raises(KeyError):
        nd1.multiply(nd2, meta_kwds_operate=['exposure'])
    with pytest.raises(KeyError):
        nd1.divide(nd2, meta_kwds_operate=['exposure'])

# Covering:
# only first has meta and the other is a scalar
def test_arithmetics_meta_one():
    meta1 = {'exposure': 100, 'dummy': 5}
    nd1 = NDDataArithmetic(1, meta=meta1)
    nd2 = NDDataArithmetic(50)

    # No keywords specified to be part of arithmetics
    nd3 = nd1.add(nd2)
    assert nd3.meta['exposure'] == 100
    assert nd3.meta['dummy'] == 5

    # keyword IS specified to be part of arithmetics
    nd4 = nd1.add(nd2, meta_kwds_operate=['exposure'])
    assert nd4.meta['exposure'] == 150
    assert nd4.meta['dummy'] == 5

    nd4 = nd1.subtract(nd2, meta_kwds_operate=['exposure'])
    assert nd4.meta['exposure'] == 50
    assert nd4.meta['dummy'] == 5

    nd4 = nd1.multiply(nd2, meta_kwds_operate=['exposure'])
    assert nd4.meta['exposure'] == 5000
    assert nd4.meta['dummy'] == 5

    nd4 = nd1.divide(nd2, meta_kwds_operate=['exposure'])
    assert nd4.meta['exposure'] == 2
    assert nd4.meta['dummy'] == 5

    # Short check that units are ignored:
    nd2 = NDDataArithmetic(50, unit=u.s)
    nd4 = nd1.multiply(nd2, meta_kwds_operate=['exposure'])
    assert nd4.meta['exposure'] == 5000 # would fail if it was a quantity

# Covering:
# only second has meta and the other is a scalar, the inverse case
def test_arithmetics_meta_one_inverse():
    meta1 = {'exposure': 100, 'dummy': 5}
    nd1 = NDDataArithmetic(1, meta=meta1)
    nd2 = NDDataArithmetic(50)

    # Now also the inverse operations are of interest
    nd4 = nd2.add(nd1, meta_kwds_operate=['exposure'])
    assert nd4.meta['exposure'] == 150
    assert nd4.meta['dummy'] == 5

    nd4 = nd2.subtract(nd1, meta_kwds_operate=['exposure'])
    assert nd4.meta['exposure'] == -50
    assert nd4.meta['dummy'] == 5

    nd4 = nd2.multiply(nd1, meta_kwds_operate=['exposure'])
    assert nd4.meta['exposure'] == 5000
    assert nd4.meta['dummy'] == 5

    nd4 = nd2.divide(nd1, meta_kwds_operate=['exposure'])
    assert nd4.meta['exposure'] == 0.5
    assert nd4.meta['dummy'] == 5

# Covering:
# only one has meta and the other is NOT a scalar
def test_arithmetics_meta_one_invalid():
    meta1 = {'exposure': 100, 'dummy': 5}
    nd1 = NDDataArithmetic(1, meta=meta1)
    nd2 = NDDataArithmetic([1,2,3])

    # Considered invalid would be if only one has a meta and the other
    # is an ndarray but this case only raises a warning
    nd3 = nd1.add(nd2, meta_kwds_operate=['exposure'])
    assert nd3.meta['exposure'] == 100

    # Check also for reverse
    nd3 = nd2.add(nd1, meta_kwds_operate=['exposure'])
    assert nd3.meta['exposure'] == 100

# Covering:
# both have uncertainties (data and uncertainty without unit)
# tested against manually determined resulting uncertainties to verify the
# implemented formulas
# this test only works as long as data1 and data2 do not contain any 0
def test_arithmetics_stddevuncertainty_basic():
    nd1 = NDDataArithmetic([1,2,3], uncertainty=StdDevUncertainty([1,1,3]))
    nd2 = NDDataArithmetic([2,2,2], uncertainty=StdDevUncertainty([2,2,2]))
    nd3 = nd1.add(nd2)
    nd4 = nd2.add(nd1)
    # Inverse operation should result in the same uncertainty
    assert_array_equal(nd3.uncertainty.array, nd4.uncertainty.array)
    # Compare it to the theoretical uncertainty
    ref_uncertainty = np.sqrt(np.array([1,1,3])**2+np.array([2,2,2])**2)
    assert_array_equal(nd3.uncertainty.array, ref_uncertainty)

    nd3 = nd1.subtract(nd2)
    nd4 = nd2.subtract(nd1)
    # Inverse operation should result in the same uncertainty
    assert_array_equal(nd3.uncertainty.array, nd4.uncertainty.array)
    # Compare it to the theoretical uncertainty (same as for add)
    assert_array_equal(nd3.uncertainty.array, ref_uncertainty)

    # Multiplication and Division only work with almost equal array comparisons
    # since the formula implemented and the formula used as reference are
    # slightly different.
    nd3 = nd1.multiply(nd2)
    nd4 = nd2.multiply(nd1)
    # Inverse operation should result in the same uncertainty
    assert_array_almost_equal(nd3.uncertainty.array, nd4.uncertainty.array)
    # Compare it to the theoretical uncertainty
    ref_uncertainty = np.abs(np.array([2,4,6])) * np.sqrt(
        (np.array([1,1,3]) / np.array([1,2,3]))**2 +
        (np.array([2,2,2]) / np.array([2,2,2]))**2)
    assert_array_almost_equal(nd3.uncertainty.array, ref_uncertainty)

    nd3 = nd1.divide(nd2)
    nd4 = nd2.divide(nd1)
    # Inverse operation gives a different uncertainty!
    # Compare it to the theoretical uncertainty
    ref_uncertainty_1 = np.abs(np.array([1/2,2/2,3/2])) * np.sqrt(
        (np.array([1,1,3]) / np.array([1,2,3]))**2 +
        (np.array([2,2,2]) / np.array([2,2,2]))**2)
    assert_array_almost_equal(nd3.uncertainty.array, ref_uncertainty_1)
    ref_uncertainty_2 = np.abs(np.array([2,1,2/3])) * np.sqrt(
        (np.array([1,1,3]) / np.array([1,2,3]))**2 +
        (np.array([2,2,2]) / np.array([2,2,2]))**2)
    assert_array_almost_equal(nd4.uncertainty.array, ref_uncertainty_2)

# Tests for correlation, covering
# correlation between -1 and 1 with correlation term being positive / negative
# also with one data being once positive and once completly negative
# The point of this test is to compare the used formula to the theoretical one.
# TODO: Maybe covering units too but I think that should work because of
# the next tests. Also this may be reduced somehow.
@pytest.mark.parametrize(('cor','uncert1','data2'), [
    (-1, [1,1,3], [2,2,7]),
    (-0.5, [1,1,3], [2,2,7]),
    (-0.25, [1,1,3], [2,2,7]),
    (0, [1,1,3], [2,2,7]),
    (0.25, [1,1,3], [2,2,7]),
    (0.5, [1,1,3], [2,2,7]),
    (1, [1,1,3], [2,2,7]),
    (-1, [-1,-1,-3], [2,2,7]),
    (-0.5, [-1,-1,-3], [2,2,7]),
    (-0.25, [-1,-1,-3], [2,2,7]),
    (0, [-1,-1,-3], [2,2,7]),
    (0.25, [-1,-1,-3], [2,2,7]),
    (0.5, [-1,-1,-3], [2,2,7]),
    (1, [-1,-1,-3], [2,2,7]),
    (-1, [1,1,3], [-2,-3,-2]),
    (-0.5, [1,1,3], [-2,-3,-2]),
    (-0.25, [1,1,3], [-2,-3,-2]),
    (0, [1,1,3], [-2,-3,-2]),
    (0.25, [1,1,3], [-2,-3,-2]),
    (0.5, [1,1,3], [-2,-3,-2]),
    (1, [1,1,3], [-2,-3,-2]),
    (-1, [-1,-1,-3], [-2,-3,-2]),
    (-0.5, [-1,-1,-3], [-2,-3,-2]),
    (-0.25, [-1,-1,-3], [-2,-3,-2]),
    (0, [-1,-1,-3], [-2,-3,-2]),
    (0.25, [-1,-1,-3], [-2,-3,-2]),
    (0.5, [-1,-1,-3], [-2,-3,-2]),
    (1, [-1,-1,-3], [-2,-3,-2]),
    ])
def test_arithmetics_stddevuncertainty_basic_with_correlation(cor, uncert1, data2):
    data1 = np.array([1,2,3])
    data2 = np.array(data2)
    uncert1 = np.array(uncert1)
    uncert2 = np.array([2,2,2])
    nd1 = NDDataArithmetic(data1, uncertainty=StdDevUncertainty(uncert1))
    nd2 = NDDataArithmetic(data2, uncertainty=StdDevUncertainty(uncert2))
    nd3 = nd1.add(nd2, uncertainty_correlation=cor)
    nd4 = nd2.add(nd1, uncertainty_correlation=cor)
    # Inverse operation should result in the same uncertainty
    assert_array_equal(nd3.uncertainty.array, nd4.uncertainty.array)
    # Compare it to the theoretical uncertainty
    ref_uncertainty = np.sqrt(uncert1**2 + uncert2**2 + 2 * cor * uncert1 * uncert2)
    assert_array_equal(nd3.uncertainty.array, ref_uncertainty)

    nd3 = nd1.subtract(nd2, uncertainty_correlation=cor)
    nd4 = nd2.subtract(nd1, uncertainty_correlation=cor)
    # Inverse operation should result in the same uncertainty
    assert_array_equal(nd3.uncertainty.array, nd4.uncertainty.array)
    # Compare it to the theoretical uncertainty
    ref_uncertainty = np.sqrt(uncert1**2 + uncert2**2 - 2 * cor * uncert1 * uncert2)
    assert_array_equal(nd3.uncertainty.array, ref_uncertainty)

    # Multiplication and Division only work with almost equal array comparisons
    # since the formula implemented and the formula used as reference are
    # slightly different.
    nd3 = nd1.multiply(nd2, uncertainty_correlation=cor)
    nd4 = nd2.multiply(nd1, uncertainty_correlation=cor)
    # Inverse operation should result in the same uncertainty
    assert_array_almost_equal(nd3.uncertainty.array, nd4.uncertainty.array)
    # Compare it to the theoretical uncertainty
    ref_uncertainty = (np.abs(data1 * data2)) * np.sqrt(
        (uncert1 / data1)**2 + (uncert2 / data2)**2 +
        (2 * cor * uncert1 * uncert2 / (data1 * data2)))
    assert_array_almost_equal(nd3.uncertainty.array, ref_uncertainty)

    nd3 = nd1.divide(nd2, uncertainty_correlation=cor)
    nd4 = nd2.divide(nd1, uncertainty_correlation=cor)
    # Inverse operation gives a different uncertainty!
    # Compare it to the theoretical uncertainty
    ref_uncertainty_1 = (np.abs(data1 / data2)) * np.sqrt(
        (uncert1 / data1)**2 + (uncert2 / data2)**2 -
        (2 * cor * uncert1 * uncert2 / (data1 * data2)))
    assert_array_almost_equal(nd3.uncertainty.array, ref_uncertainty_1)
    ref_uncertainty_2 = (np.abs(data2 / data1)) * np.sqrt(
        (uncert1 / data1)**2 + (uncert2 / data2)**2 -
        (2 * cor * uncert1 * uncert2 / (data1 * data2)))
    assert_array_almost_equal(nd4.uncertainty.array, ref_uncertainty_2)

# Covering:
# only one has an uncertainty (data and uncertainty without unit)
# tested against the case where the other one has zero uncertainty. (this case
# must be correct because we tested it in the last case)
# Also verify that if the result of the data has negative values the resulting
# uncertainty has no negative values.
def test_arithmetics_stddevuncertainty_one_missing():
    nd1 = NDDataArithmetic([1,-2,3])
    nd1_ref = NDDataArithmetic([1,-2,3], uncertainty=StdDevUncertainty([0,0,0]))
    nd2 = NDDataArithmetic([2,2,-2], uncertainty=StdDevUncertainty([2,2,2]))

    # Addition
    nd3 = nd1.add(nd2)
    nd3_ref = nd1_ref.add(nd2)
    assert_array_equal(nd3.uncertainty.array, nd3_ref.uncertainty.array)
    assert_array_equal(np.abs(nd3.uncertainty.array), nd3.uncertainty.array)

    nd3 = nd2.add(nd1)
    nd3_ref = nd2.add(nd1_ref)
    assert_array_equal(nd3.uncertainty.array, nd3_ref.uncertainty.array)
    assert_array_equal(np.abs(nd3.uncertainty.array), nd3.uncertainty.array)

    # Subtraction
    nd3 = nd1.subtract(nd2)
    nd3_ref = nd1_ref.subtract(nd2)
    assert_array_equal(nd3.uncertainty.array, nd3_ref.uncertainty.array)
    assert_array_equal(np.abs(nd3.uncertainty.array), nd3.uncertainty.array)

    nd3 = nd2.subtract(nd1)
    nd3_ref = nd2.subtract(nd1_ref)
    assert_array_equal(nd3.uncertainty.array, nd3_ref.uncertainty.array)
    assert_array_equal(np.abs(nd3.uncertainty.array), nd3.uncertainty.array)

    # Multiplication
    nd3 = nd1.multiply(nd2)
    nd3_ref = nd1_ref.multiply(nd2)
    assert_array_equal(nd3.uncertainty.array, nd3_ref.uncertainty.array)
    assert_array_equal(np.abs(nd3.uncertainty.array), nd3.uncertainty.array)

    nd3 = nd2.multiply(nd1)
    nd3_ref = nd2.multiply(nd1_ref)
    assert_array_equal(nd3.uncertainty.array, nd3_ref.uncertainty.array)
    assert_array_equal(np.abs(nd3.uncertainty.array), nd3.uncertainty.array)

    # Division
    nd3 = nd1.divide(nd2)
    nd3_ref = nd1_ref.divide(nd2)
    assert_array_equal(nd3.uncertainty.array, nd3_ref.uncertainty.array)
    assert_array_equal(np.abs(nd3.uncertainty.array), nd3.uncertainty.array)

    nd3 = nd2.divide(nd1)
    nd3_ref = nd2.divide(nd1_ref)
    assert_array_equal(nd3.uncertainty.array, nd3_ref.uncertainty.array)
    assert_array_equal(np.abs(nd3.uncertainty.array), nd3.uncertainty.array)

# Covering:
# data with unit and uncertainty with unit (but equivalent units)
# compared against correctly scaled NDDatas
@pytest.mark.parametrize(('uncert1', 'uncert2'), [
    (np.array([1,2,3])*u.m, None),
    (np.array([1,2,3])*u.cm, None),
    (None, np.array([1,2,3])*u.m),
    (None, np.array([1,2,3])*u.cm),
    (np.array([1,2,3]), np.array([2,3,4])),
    (np.array([1,2,3])*u.m, np.array([2,3,4])),
    (np.array([1,2,3]), np.array([2,3,4]))*u.m,
    (np.array([1,2,3])*u.m, np.array([2,3,4]))*u.m,
    (np.array([1,2,3])*u.cm, np.array([2,3,4])),
    (np.array([1,2,3]), np.array([2,3,4]))*u.cm,
    (np.array([1,2,3])*u.cm, np.array([2,3,4]))*u.cm,
    (np.array([1,2,3])*u.km, np.array([2,3,4]))*u.cm,
    ])
def test_arithmetics_stddevuncertainty_with_units(uncert1, uncert2):
    # Data has same units
    data1 = np.array([1,2,3]) * u.m
    data2 = np.array([-4,7,0]) * u.m
    if uncert1 is not None:
        uncert1 = StdDevUncertainty(uncert1)
        if isinstance(uncert1, Quantity):
            uncert1_ref = uncert1.to(data1.unit).value
        else:
            uncert1_ref = uncert1
        uncert_ref1 = StdDevUncertainty(uncert1_ref, copy=True)
    else:
        uncert1 = None
        uncert_ref1 = None

    if uncert2 is not None:
        uncert2 = StdDevUncertainty(uncert2)
        if isinstance(uncert2, Quantity):
            uncert2_ref = uncert2.to(data2.unit).value
        else:
            uncert2_ref = uncert2
        uncert_ref2 = StdDevUncertainty(uncert2_ref, copy=True)
    else:
        uncert2 = None
        uncert_ref2 = None

    nd1 = NDDataArithmetic(data1, uncertainty = uncert1)
    nd2 = NDDataArithmetic(data2, uncertainty = uncert2)

    nd1_ref = NDDataArithmetic(data1, uncertainty = uncert_ref1)
    nd2_ref = NDDataArithmetic(data2, uncertainty = uncert_ref2)

    # Let's start the tests
    # Addition
    nd3 = nd1.add(nd2)
    nd3_ref = nd1_ref.add(nd2_ref)
    assert nd3.unit == nd3_ref.unit
    assert nd3.uncertainty.unit == nd3_ref.uncertainty.unit
    assert_array_equal(nd3.uncertainty.array, nd3.uncertainty.array)

    nd3 = nd2.add(nd1)
    nd3_ref = nd2_ref.add(nd1_ref)
    assert nd3.unit == nd3_ref.unit
    assert nd3.uncertainty.unit == nd3_ref.uncertainty.unit
    assert_array_equal(nd3.uncertainty.array, nd3.uncertainty.array)

    # Subtraction
    nd3 = nd1.subtract(nd2)
    nd3_ref = nd1_ref.subtract(nd2_ref)
    assert nd3.unit == nd3_ref.unit
    assert nd3.uncertainty.unit == nd3_ref.uncertainty.unit
    assert_array_equal(nd3.uncertainty.array, nd3.uncertainty.array)

    nd3 = nd2.subtract(nd1)
    nd3_ref = nd2_ref.subtract(nd1_ref)
    assert nd3.unit == nd3_ref.unit
    assert nd3.uncertainty.unit == nd3_ref.uncertainty.unit
    assert_array_equal(nd3.uncertainty.array, nd3.uncertainty.array)

    # Multiplication
    nd3 = nd1.multiply(nd2)
    nd3_ref = nd1_ref.multiply(nd2_ref)
    assert nd3.unit == nd3_ref.unit
    assert nd3.uncertainty.unit == nd3_ref.uncertainty.unit
    assert_array_equal(nd3.uncertainty.array, nd3.uncertainty.array)

    nd3 = nd2.multiply(nd1)
    nd3_ref = nd2_ref.multiply(nd1_ref)
    assert nd3.unit == nd3_ref.unit
    assert nd3.uncertainty.unit == nd3_ref.uncertainty.unit
    assert_array_equal(nd3.uncertainty.array, nd3.uncertainty.array)

    # Division
    nd3 = nd1.divide(nd2)
    nd3_ref = nd1_ref.divide(nd2_ref)
    assert nd3.unit == nd3_ref.unit
    assert nd3.uncertainty.unit == nd3_ref.uncertainty.unit
    assert_array_equal(nd3.uncertainty.array, nd3.uncertainty.array)

    nd3 = nd2.divide(nd1)
    nd3_ref = nd2_ref.divide(nd1_ref)
    assert nd3.unit == nd3_ref.unit
    assert nd3.uncertainty.unit == nd3_ref.uncertainty.unit
    assert_array_equal(nd3.uncertainty.array, nd3.uncertainty.array)


def test_arithmetics_handle_switches():
    meta1 = {'a': 1}
    meta2 = {'b': 2}
    mask1 = True
    mask2 = False
    uncertainty1 = StdDevUncertainty([1,2,3])
    uncertainty2 = StdDevUncertainty([1,2,3])
    wcs1 = 5
    wcs2 = 100
    data1 = [1,1,1]
    data2 = [1,1,1]

    nd1 = NDDataArithmetic(data1, meta=meta1, mask=mask1, wcs=wcs1,
                           uncertainty=uncertainty1)
    nd2 = NDDataArithmetic(data2, meta=meta2, mask=mask2, wcs=wcs2,
                           uncertainty=uncertainty2)
    nd3 = NDDataArithmetic(data1)

    # Both have the attributes but option None is chosen
    nd_ = nd1.add(nd2, propagate_uncertainties=None, handle_meta=None,
                  handle_mask=None, compare_wcs=None)
    assert nd_.wcs is None
    assert len(nd_.meta) == 0
    assert nd_.mask is None
    assert nd_.uncertainty is None

    # Only second has attributes and False is chosen
    nd_ = nd3.add(nd2, propagate_uncertainties=False, handle_meta=False,
                  handle_mask=False, compare_wcs=False)
    assert nd_.wcs == wcs2
    assert nd_.meta == meta2
    assert nd_.mask == mask2
    assert_array_equal(nd_.uncertainty.array, uncertainty2.array)

    # Only first has attributes and False is chosen
    nd_ = nd1.add(nd3, propagate_uncertainties=False, handle_meta=False,
                  handle_mask=False, compare_wcs=False)
    assert nd_.wcs == wcs1
    assert nd_.meta == meta1
    assert nd_.mask == mask1
    assert_array_equal(nd_.uncertainty.array, uncertainty1.array)




# Old tests belonging to NDDataArray and only partly touch NDArithmetic itself



from ...compat import NDDataArray
from ...nduncertainty import IncompatibleUncertaintiesException
from ....utils import NumpyRNGContext


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


# Uncertainty tests
@pytest.mark.xfail # Shape checks for uncertainty are deleted
def test_nddata_uncertainty_init_invalid_shape_1():
    u = StdDevUncertainty(array=np.ones((6, 6)))
    with pytest.raises(ValueError) as exc:
        NDDataArray(np.ones((5, 5)), uncertainty=u)
    assert exc.value.args[0] == 'parent shape does not match array data shape'

@pytest.mark.xfail # Shape checks for uncertainty are deleted
def test_nddata_uncertainty_init_invalid_shape_2():
    u = StdDevUncertainty()
    NDDataArray(np.ones((5, 5)), uncertainty=u)
    with pytest.raises(ValueError) as exc:
        u.array = np.ones((6, 6))
    assert exc.value.args[0] == 'array shape does not match parent data shape'


@pytest.mark.parametrize(('uncertainty'), [1., 'spam', np.ones((5, 5))])
def test_nddata_uncertainty_invalid_type(uncertainty):
    with pytest.raises(TypeError) as exc:
        NDDataArray(np.ones((5, 5)), uncertainty=uncertainty)
    assert exc.value.args[0] == ('Uncertainty must be an instance of '
                                 'a NDUncertainty object')


# slicing tests
def test_simple_slicing():
    u1 = StdDevUncertainty(array=np.ones((5, 5)) * 3)
    d1 = NDDataArray(np.ones((5, 5)), uncertainty=u1)
    assert d1.data.shape == (5, 5)
    d2 = d1[2:3, 2:3]
    assert d2.data.shape == (1, 1)
    d3 = d1[2, 2]
    assert d3.data.shape == ()


def test_slicing_reference():
    u1 = StdDevUncertainty(array=np.ones((5, 5)) * 3)
    d1 = NDDataArray(np.ones((5, 5)), uncertainty=u1)
    d2 = d1[2:3, 2:3]
    # asserting that the new nddata contains references to the original nddata
    assert d2.data.base is d1.data
    assert d2.uncertainty.array.base is d1.uncertainty.array


def test_slicing_with_mask_or_flag():
    # Regression test for #2170
    ndd = NDDataArray(np.array([1, 2, 3]),
                           mask=np.array([False, False, False]))
    assert ndd[0].data.shape == ()
    assert not ndd[0].mask


def test_nddata_add():
    d1 = NDDataArray(np.ones((5, 5)))
    d2 = NDDataArray(np.ones((5, 5)))
    d3 = d1.add(d2)
    assert np.all(d3.data == 2.)

@pytest.mark.xfail # Just a warning now not an exception
def test_nddata_add_mismatch_wcs():
    d1 = NDDataArray(np.ones((5, 5)), wcs=1.)
    d2 = NDDataArray(np.ones((5, 5)), wcs=2.)
    with pytest.raises(ValueError) as exc:
        d1.add(d2)
    assert exc.value.args[0] == "WCS properties do not match"

@pytest.mark.xfail # There is no unit test anymore
def test_nddata_add_mismatch_units():
    d1 = NDDataArray(np.ones((5, 5)), unit='Jy')
    d2 = NDDataArray(np.ones((5, 5)), unit='erg/s')
    with pytest.raises(ValueError) as exc:
        d1.add(d2)
    assert exc.value.args[0] == "operand units do not match"

@pytest.mark.xfail # There is no shape check anymore
def test_nddata_add_mismatch_shape():
    d1 = NDDataArray(np.ones((5, 5)))
    d2 = NDDataArray(np.ones((6, 6)))
    with pytest.raises(ValueError) as exc:
        d1.add(d2)
    assert exc.value.args[0] == "operand shapes do not match"


def test_nddata_add_with_masks():
    # numpy masked arrays mask the result of binary operations if the
    # mask of either operand is set.
    # Does NDData?
    ndd1 = NDDataArray(np.array([1, 2]), mask=np.array([True, False]))
    other_mask = ~ ndd1.mask
    ndd2 = NDDataArray(np.array([1, 2]), mask=other_mask)
    result = ndd1.add(ndd2)
    # The result should have all entries masked...
    assert result.mask.all()


def test_nddata_add_uncertainties():
    u1 = StdDevUncertainty(array=np.ones((5, 5)) * 3)
    u2 = StdDevUncertainty(array=np.ones((5, 5)))
    d1 = NDDataArray(np.ones((5, 5)), uncertainty=u1)
    d2 = NDDataArray(np.ones((5, 5)), uncertainty=u2)
    d3 = d1.add(d2)
    assert np.all(d3.data == 2.)
    assert_array_equal(d3.uncertainty.array, np.sqrt(10.))

@pytest.mark.xfail # Altered the way in which propagation checks for compatibility
def test_nddata_add_uncertainties_mismatch():
    u1 = StdDevUncertainty(array=np.ones((5, 5)) * 3)
    u2 = FakeUncertainty()
    d1 = NDDataArray(np.ones((5, 5)), uncertainty=u1)
    d2 = NDDataArray(np.ones((5, 5)), uncertainty=u2)
    with pytest.raises(IncompatibleUncertaintiesException) as exc:
        d3 = d1.add(d2)
    assert exc.value.args[0] == ('Cannot propagate uncertainties of type '
                                 'StdDevUncertainty with uncertainties of '
                                 'type FakeUncertainty for addition')

@pytest.mark.xfail # Altered the way in which uncertainty init checks for compatibility
def test_initializing_nduncertainty_from_quantity():
    # Until nddata and quantity are integrated initializing with a quantity
    # should raise an error.
    input_ndd_unit = u.kg
    ndd = NDDataArray(np.array([1, 2, 3]), unit=input_ndd_unit)
    std_data = np.array([1, 2, 3])

    # Unit of the uncertainty not convertible to unit of ndd, should raise
    # an error.
    std_error = StdDevUncertainty(u.adu * std_data)
    assert std_error._unit is u.adu
    with pytest.raises(u.UnitsError):
        ndd.uncertainty = std_error

    # Uncertainty should be settable without any change in its values
    # because uncertainty unit is same as data unit.
    std_error = StdDevUncertainty(u.kg * std_data)
    ndd.uncertainty = std_error
    assert_array_equal(std_data, ndd.uncertainty.array)

    # If the uncertainty unit is grams there should be no error, but the
    # values of the uncertainty should be scaled.
    std_error = StdDevUncertainty(u.g * std_data)
    ndd.uncertainty = std_error
    assert_array_equal(std_data, 1000 * ndd.uncertainty.array)

    # If ndd has no unit but the uncertainty does an error should be raised.
    ndd = NDDataArray(np.array([1, 2, 3]), unit=None)
    with pytest.raises(ValueError):
        ndd.uncertainty = std_error


# Mask tests
def test_unmasked_masked_array_input():
    # Test for #2784
    marr = np.ma.array([1, 2, 5])  # Masked array with no masked entries
    nd = NDDataArray(marr)  # Before fix this raised a ValueError

    # Check that masks are correct
    assert marr.mask is np.ma.nomask
    # Internal representation is np.ma.nomask but getter returns None.
    assert nd.mask is None


def test_nddata_unmasked_in_operation_with_masked_numpy_array():
    # test for #2417
    ndd = NDDataArray(np.array([1, 2, 3]))
    np_data = -np.ones_like(ndd)
    np_mask = np.array([True, False, True])
    np_arr_masked = np.ma.masked_array(np_data, mask=np_mask, copy=True)
    # check multiplication in both orders as in test above
    result1 = ndd * np_arr_masked
    result2 = np_arr_masked * ndd
    for result in [result1, result2]:
        # multiplying by a masked numpy array should return a masked array
        assert isinstance(result, np.ma.MaskedArray)
        assert np.all(result.mask == np_mask)
        assert np.all(result[~result.mask] == -ndd.data[~np_mask])


@pytest.mark.parametrize(('shape'), [(10,), (5, 5), (3, 10, 10)])
def test_nddata_mask_invalid_shape(shape):
    with pytest.raises(ValueError) as exc:
        with NumpyRNGContext(789):
            NDDataArray(np.random.random((10, 10)),
                             mask=np.random.random(shape) > 0.5)
    assert exc.value.args[0] == 'dimensions of mask do not match data'


@pytest.mark.parametrize('mask_in', [
                         np.array([True, False]),
                         np.array([1, 0]),
                         [True, False],
                         [1, 0]])
def test_nddata_mask_init_without_np_array(mask_in):
    ndd = NDDataArray(np.array([1, 1]), mask=mask_in)
    assert (ndd.mask == mask_in).all()


def test_ndddata_with_mask_acts_like_masked_array():
    # test for #2414
    input_mask = np.array([True, False, False])
    input_data = np.array([1, 2, 3])
    ndd_masked = NDDataArray(input_data.copy(), mask=input_mask.copy())
    other = - np.ones_like(input_data)
    result1 = ndd_masked * other
    result2 = other * ndd_masked
    # Test for both orders of multiplication -- if multiplication is
    # eventually overridden for NDData the result can depend on order.
    for result in [result1, result2]:
        # Result should be a masked array because input NDData was masked
        assert isinstance(result, np.ma.MaskedArray)
        # Result mask should match input mask because other has no mask
        assert np.all(result.mask == input_mask)
        assert np.all(result[~result.mask] == - input_data[~input_mask])


# Arithmetic tests

def test_nddata_subtract():
    d1 = NDDataArray(np.ones((5, 5)))
    d2 = NDDataArray(np.ones((5, 5)) * 2.)
    d3 = d1.subtract(d2)
    assert np.all(d3.data == -1.)

@pytest.mark.xfail # WCS inequality only raises a warning now
def test_nddata_subtract_mismatch_wcs():
    d1 = NDDataArray(np.ones((5, 5)), wcs=1.)
    d2 = NDDataArray(np.ones((5, 5)) * 2., wcs=2.)
    with pytest.raises(ValueError) as exc:
        d1.subtract(d2)
    assert exc.value.args[0] == "WCS properties do not match"

@pytest.mark.xfail # Unit mismatch is not explicitly checked and catched
def test_nddata_subtract_mismatch_units():
    d1 = NDDataArray(np.ones((5, 5)), unit='Jy')
    d2 = NDDataArray(np.ones((5, 5)) * 2., unit='erg/s')
    with pytest.raises(ValueError) as exc:
        d1.subtract(d2)
    assert exc.value.args[0] == "operand units do not match"

@pytest.mark.xfail # Shape mismatch is not explicitly checked and catched
def test_nddata_subtract_mismatch_shape():
    d1 = NDDataArray(np.ones((5, 5)))
    d2 = NDDataArray(np.ones((6, 6)) * 2.)
    with pytest.raises(ValueError) as exc:
        d1.subtract(d2)
    assert exc.value.args[0] == "operand shapes do not match"


def test_nddata_subtract_uncertainties():
    u1 = StdDevUncertainty(array=np.ones((5, 5)) * 3)
    u2 = StdDevUncertainty(array=np.ones((5, 5)))
    d1 = NDDataArray(np.ones((5, 5)), uncertainty=u1)
    d2 = NDDataArray(np.ones((5, 5)) * 2., uncertainty=u2)
    d3 = d1.subtract(d2)
    assert np.all(d3.data == -1.)
    assert_array_equal(d3.uncertainty.array, np.sqrt(10.))


def test_nddata_multiply_uncertainties():
    u1 = StdDevUncertainty(array=np.ones((5, 5)) * 3)
    u2 = StdDevUncertainty(array=np.ones((5, 5)))
    d1 = NDDataArray(np.ones((5, 5)), uncertainty=u1)
    d2 = NDDataArray(np.ones((5, 5)) * 2., uncertainty=u2)
    d3 = d1.multiply(d2)
    assert np.all(d3.data == 2.)
    assert_array_equal(d3.uncertainty.array, 2 * np.sqrt(9.25))


def test_nddata_divide_uncertainties():
    u1 = StdDevUncertainty(array=np.ones((5, 5)) * 3)
    u2 = StdDevUncertainty(array=np.ones((5, 5)))
    d1 = NDDataArray(np.ones((5, 5)), uncertainty=u1)
    d2 = NDDataArray(np.ones((5, 5)) * 2., uncertainty=u2)
    d3 = d1.divide(d2)
    assert np.all(d3.data == 0.5)
    assert_array_equal(d3.uncertainty.array, 0.5 * np.sqrt(9.25))

@pytest.mark.xfail # Altered the check which uncertainties are compatible
def test_nddata_subtract_uncertainties_mismatch():
    u1 = StdDevUncertainty(array=np.ones((5, 5)) * 3)
    u2 = FakeUncertainty()
    d1 = NDDataArray(np.ones((5, 5)), uncertainty=u1)
    d2 = NDDataArray(np.ones((5, 5)) * 2., uncertainty=u2)
    with pytest.raises(IncompatibleUncertaintiesException) as exc:
        d3 = d1.subtract(d2)
    assert exc.value.args[0] == ('Cannot propagate uncertainties of type '
                                 'StdDevUncertainty with uncertainties of '
                                 'type FakeUncertainty for subtraction')


@pytest.mark.parametrize(('op1_unc','op2_unc'), [
                         (None, None),
                         (StdDevUncertainty([1]), None),
                         (None, StdDevUncertainty([1])),
                         (StdDevUncertainty([1]), StdDevUncertainty([1]))
                         ])
def test_arithmetic_result_not_tied_to_operands_uncertainty(op1_unc, op2_unc):
    # Expectation is that the result of an arithmetic operation should be a
    # new object whose members are not tied to the members of the operand.
    # The only reliable test of this is to change elements of the result and
    # see if the corresponding elements of the operands change.
    # All four of the cases parametrized in this test do need to be checked
    # because each of the four cases is handled separately in the code (well,
    # except for the None, None case).
    # Only one of the arithmetic operations need to be checked because the
    # logic for propagating the uncertainties is common to all of the
    # operations.
    op1 = NDDataArray(np.array([1]), uncertainty=op1_unc)
    op2 = NDDataArray(np.array([1]), uncertainty=op2_unc)

    result = op1.add(op2)
    if result.uncertainty:
        result.uncertainty.array[0] = 0

    if op1_unc:
        assert op1.uncertainty.array[0] == 1
    if op2_unc:
        assert op2.uncertainty.array[0] == 1

    result.data[0] = np.pi
    assert op1.data[0] == 1
    assert op2.data[0] == 1


@pytest.mark.parametrize(('op1_mask','op2_mask'), [
                         (None, None),
                         (None, np.array([False])),
                         (np.array([False]), None),
                         (np.array([False]), np.array([False]))])
def test_arithmetic_result_not_tied_to_operands_mask(op1_mask, op2_mask):
    # See test_arithmetic_result_not_tied_to_operands_uncertainty for comments
    op1 = NDDataArray(np.array([1]), mask=op1_mask)
    op2 = NDDataArray(np.array([1]), mask=op2_mask)
    result = op1.add(op2)

    if result.mask is not None:
        result.mask[0] = True

    if op1_mask is not None:
        assert op1.mask[0] == (not result.mask[0])

    if op2_mask is not None:
        assert op2.mask[0] == (not result.mask[0])


def test_arithmetic_result_not_tied_to_operands_wcs():
    # unit is no longer settable, so test that result unit is different object
    # than operands was removed.

    # Unlike the previous two tests, we only need to check a case where both
    # operands have the same wcs because operands with different wcs is not
    # supported
    op1 = NDDataArray(np.array([1]), wcs=np.array([1]), unit='m')
    op2 = NDDataArray(np.array([1]), wcs=np.array([1]), unit='m')
    result = op1.add(op2)
    result.wcs[0] = 12345
    assert op1.wcs[0] != result.wcs[0]
    assert op2.wcs[0] != result.wcs[0]


# first operand has unit km, second has unit m
@pytest.mark.xfail # TODO First two tests may fail, but I have not checked why
@pytest.mark.parametrize(('operation','result_unit'), [
                         ('add', u.km),
                         ('subtract', u.km),
                         ('multiply', u.km * u.m),
                         ('divide', u.km / u.m)])
def test_uncertainty_unit_conversion_add_subtract(operation, result_unit):
    in_km = NDDataArray(np.array([1, 1]), unit=u.km,
                             uncertainty=StdDevUncertainty([.1, .1]))
    in_m = NDDataArray(in_km.data * 1000, unit=u.m)
    in_m.uncertainty = StdDevUncertainty(in_km.uncertainty.array * 1000)
    operator_km = in_km.__getattribute__(operation)
    combined = operator_km(in_m)
    assert combined.unit == result_unit
    if operation in ['add', 'subtract']:
        # uncertainty is not scaled by result values
        assert_array_equal(combined.uncertainty.array,
                           np.sqrt(2) * in_km.uncertainty.array)
    else:
        # uncertainty is scaled by result
        assert_array_equal(combined.uncertainty.array,
            np.sqrt(2) * in_km.uncertainty.array * combined.data)


@pytest.mark.parametrize(('unit1','unit2','op','result_unit'), [
                         (None, None, 'add', None),
                         (None, None, 'multiply', None),
                         (None, u.m, 'multiply', u.m),
                         (u.dimensionless_unscaled, None, 'multiply',
                          u.dimensionless_unscaled),
                         (u.adu, u.adu, 'add', u.adu),
                         (u.adu, u.adu, 'subtract', u.adu),
                         (u.adu, u.adu, 'divide', u.dimensionless_unscaled),
                         (u.adu, u.m, 'multiply', u.m * u.adu)
                         ])
def test_arithmetic_unit_calculation(unit1, unit2, op, result_unit):
    # Test for #2413
    ndd1 = NDDataArray(np.array([1]), unit=unit1)
    ndd2 = NDDataArray(np.array([1]), unit=unit2)
    ndd1_method = ndd1.__getattribute__(op)
    result = ndd1_method(ndd2)
    assert result.unit == result_unit


# check that subclasses can require wcs and/or unit to be present and use
# _arithmetic and convert_unit_to
class SubNDData(NDDataArray):
    """
    Subclass for test initialization of subclasses in NDData._arithmetic and
    NDData.convert_unit_to
    """
    def __init__(self, *arg, **kwd):
        super(SubNDData, self).__init__(*arg, **kwd)
        if self.unit is None:
            raise ValueError("Unit for subclass must be specified")
        if self.wcs is None:
            raise ValueError("WCS for subclass must be specified")


def test_init_of_subclasses_in_arithmetic():
    with NumpyRNGContext(12345):
        data = np.ones([10, 10])
    # The wcs only needs to be not None for this test to succeed
    arr1 = SubNDData(data, unit='adu', wcs=5)
    arr2 = SubNDData(data, unit='adu', wcs=5)
    result = arr1.add(arr2)
    assert result.unit == arr1.unit
    assert result.wcs == arr1.wcs
