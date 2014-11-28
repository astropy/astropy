# Licensed under a 3-clause BSD style license - see LICENSE.rst

# TEST_UNICODE_LITERALS

from __future__ import absolute_import, division, print_function, unicode_literals

import numpy as np
from numpy.testing import assert_array_equal

from ..arithmetic import NDDataArithmetic
from ..nduncertainty import (StdDevUncertainty,
                             IncompatibleUncertaintiesException,
                             NDUncertainty)
from ...tests.helper import pytest
from ... import units as u
from ...utils import NumpyRNGContext


class FakeUncertainty(NDUncertainty):

    def __init__(self, *arg, **kwd):
        self._unit = None
        pass

    def propagate_add(self, data, final_data):
        pass

    def propagate_subtract(self, data, final_data):
        pass

    def propagate_multiply(self, data, final_data):
        pass

    def propagate_divide(self, data, final_data):
        pass

    def array(self):
        pass


# Uncertainty tests
def test_nddata_uncertainty_init_invalid_shape_1():
    u = StdDevUncertainty(array=np.ones((6, 6)))
    with pytest.raises(ValueError) as exc:
        NDDataArithmetic(np.ones((5, 5)), uncertainty=u)
    assert exc.value.args[0] == 'parent shape does not match array data shape'


def test_nddata_uncertainty_init_invalid_shape_2():
    u = StdDevUncertainty()
    NDDataArithmetic(np.ones((5, 5)), uncertainty=u)
    with pytest.raises(ValueError) as exc:
        u.array = np.ones((6, 6))
    assert exc.value.args[0] == 'array shape does not match parent data shape'


@pytest.mark.parametrize(('uncertainty'), [1., 'spam', np.ones((5, 5))])
def test_nddata_uncertainty_invalid_type(uncertainty):
    with pytest.raises(TypeError) as exc:
        NDDataArithmetic(np.ones((5, 5)), uncertainty=uncertainty)
    assert exc.value.args[0] == 'Uncertainty must be an instance of a NDUncertainty object'


def test_nddata_add():
    d1 = NDDataArithmetic(np.ones((5, 5)))
    d2 = NDDataArithmetic(np.ones((5, 5)))
    d3 = d1.add(d2)
    assert np.all(d3.data == 2.)


def test_nddata_add_mismatch_wcs():
    d1 = NDDataArithmetic(np.ones((5, 5)), wcs=1.)
    d2 = NDDataArithmetic(np.ones((5, 5)), wcs=2.)
    with pytest.raises(ValueError) as exc:
        d1.add(d2)
    assert exc.value.args[0] == "WCS properties do not match"


def test_nddata_add_mismatch_units():
    d1 = NDDataArithmetic(np.ones((5, 5)), unit='Jy')
    d2 = NDDataArithmetic(np.ones((5, 5)), unit='erg/s')
    with pytest.raises(ValueError) as exc:
        d1.add(d2)
    assert exc.value.args[0] == "operand units do not match"


def test_nddata_add_mismatch_shape():
    d1 = NDDataArithmetic(np.ones((5, 5)))
    d2 = NDDataArithmetic(np.ones((6, 6)))
    with pytest.raises(ValueError) as exc:
        d1.add(d2)
    assert exc.value.args[0] == "operand shapes do not match"


def test_nddata_add_with_masks():
    # numpy masked arrays mask the result of binary operations if the
    # mask of either operand is set.
    # Does NDData?
    ndd1 = NDDataArithmetic([1, 2], mask=np.array([True, False]))
    other_mask = ~ ndd1.mask
    ndd2 = NDDataArithmetic([1, 2], mask=other_mask)
    result = ndd1.add(ndd2)
    # The result should have all entries masked...
    assert result.mask.all()


def test_nddata_add_uncertainties():
    u1 = StdDevUncertainty(array=np.ones((5, 5)) * 3)
    u2 = StdDevUncertainty(array=np.ones((5, 5)))
    d1 = NDDataArithmetic(np.ones((5, 5)), uncertainty=u1)
    d2 = NDDataArithmetic(np.ones((5, 5)), uncertainty=u2)
    d3 = d1.add(d2)
    assert np.all(d3.data == 2.)
    assert_array_equal(d3.uncertainty.array, np.sqrt(10.))


def test_nddata_add_uncertainties_mismatch():
    u1 = StdDevUncertainty(array=np.ones((5, 5)) * 3)
    u2 = FakeUncertainty()
    d1 = NDDataArithmetic(np.ones((5, 5)), uncertainty=u1)
    d2 = NDDataArithmetic(np.ones((5, 5)), uncertainty=u2)
    with pytest.raises(IncompatibleUncertaintiesException) as exc:
        d3 = d1.add(d2)
    assert exc.value.args[0] == 'Cannot propagate uncertainties of type StdDevUncertainty with uncertainties of type FakeUncertainty for addition'


def test_initializing_nduncertainty_from_quantity():
    # Until nddata and quantity are integrated initializing with a quantity
    # should raise an error.
    input_ndd_unit = u.kg
    ndd = NDDataArithmetic(np.array([1, 2, 3]), unit=input_ndd_unit)
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
    ndd = NDDataArithmetic(np.array([1, 2, 3]), unit=None)
    with pytest.raises(ValueError):
        ndd.uncertainty = std_error


# Mask tests
def test_unmasked_masked_array_input():
    # Test for #2784
    marr = np.ma.array([1, 2, 5])  # Masked array with no masked entries
    nd = NDDataArithmetic(marr)  # Before fix this raised a ValueError

    # Check that masks are correct
    assert marr.mask is np.ma.nomask
    assert nd.mask is None  # Internal representation is np.ma.nomask but getter returns None


def test_nddata_unmasked_in_operation_with_masked_numpy_array():
    # test for #2417
    ndd = NDDataArithmetic([1, 2, 3])
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
            NDDataArithmetic(np.random.random((10, 10)), mask=np.random.random(shape) > 0.5)
    assert exc.value.args[0] == 'dimensions of mask do not match data'


@pytest.mark.parametrize('mask_in', [
                         np.array([True, False]),
                         np.array([1, 0]),
                         [True, False],
                         [1, 0]])
def test_nddata_mask_init_without_np_array(mask_in):
    ndd = NDDataArithmetic([1, 1], mask=mask_in)
    assert (ndd.mask == mask_in).all()


def test_ndddata_with_mask_acts_like_masked_array():
    # test for #2414
    input_mask = np.array([True, False, False])
    input_data = np.array([1, 2, 3])
    ndd_masked = NDDataArithmetic(input_data.copy(), mask=input_mask.copy())
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
    d1 = NDDataArithmetic(np.ones((5, 5)))
    d2 = NDDataArithmetic(np.ones((5, 5)) * 2.)
    d3 = d1.subtract(d2)
    assert np.all(d3.data == -1.)


def test_nddata_subtract_mismatch_wcs():
    d1 = NDDataArithmetic(np.ones((5, 5)), wcs=1.)
    d2 = NDDataArithmetic(np.ones((5, 5)) * 2., wcs=2.)
    with pytest.raises(ValueError) as exc:
        d1.subtract(d2)
    assert exc.value.args[0] == "WCS properties do not match"


def test_nddata_subtract_mismatch_units():
    d1 = NDDataArithmetic(np.ones((5, 5)), unit='Jy')
    d2 = NDDataArithmetic(np.ones((5, 5)) * 2., unit='erg/s')
    with pytest.raises(ValueError) as exc:
        d1.subtract(d2)
    assert exc.value.args[0] == "operand units do not match"


def test_nddata_subtract_mismatch_shape():
    d1 = NDDataArithmetic(np.ones((5, 5)))
    d2 = NDDataArithmetic(np.ones((6, 6)) * 2.)
    with pytest.raises(ValueError) as exc:
        d1.subtract(d2)
    assert exc.value.args[0] == "operand shapes do not match"


def test_nddata_subtract_uncertainties():
    u1 = StdDevUncertainty(array=np.ones((5, 5)) * 3)
    u2 = StdDevUncertainty(array=np.ones((5, 5)))
    d1 = NDDataArithmetic(np.ones((5, 5)), uncertainty=u1)
    d2 = NDDataArithmetic(np.ones((5, 5)) * 2., uncertainty=u2)
    d3 = d1.subtract(d2)
    assert np.all(d3.data == -1.)
    assert_array_equal(d3.uncertainty.array, np.sqrt(10.))


def test_nddata_multiply_uncertainties():
    u1 = StdDevUncertainty(array=np.ones((5, 5)) * 3)
    u2 = StdDevUncertainty(array=np.ones((5, 5)))
    d1 = NDDataArithmetic(np.ones((5, 5)), uncertainty=u1)
    d2 = NDDataArithmetic(np.ones((5, 5)) * 2., uncertainty=u2)
    d3 = d1.multiply(d2)
    assert np.all(d3.data == 2.)
    assert_array_equal(d3.uncertainty.array, 2 * np.sqrt(9.25))


def test_nddata_divide_uncertainties():
    u1 = StdDevUncertainty(array=np.ones((5, 5)) * 3)
    u2 = StdDevUncertainty(array=np.ones((5, 5)))
    d1 = NDDataArithmetic(np.ones((5, 5)), uncertainty=u1)
    d2 = NDDataArithmetic(np.ones((5, 5)) * 2., uncertainty=u2)
    d3 = d1.divide(d2)
    assert np.all(d3.data == 0.5)
    assert_array_equal(d3.uncertainty.array, 0.5 * np.sqrt(9.25))


def test_nddata_subtract_uncertainties_mismatch():
    u1 = StdDevUncertainty(array=np.ones((5, 5)) * 3)
    u2 = FakeUncertainty()
    d1 = NDDataArithmetic(np.ones((5, 5)), uncertainty=u1)
    d2 = NDDataArithmetic(np.ones((5, 5)) * 2., uncertainty=u2)
    with pytest.raises(IncompatibleUncertaintiesException) as exc:
        d3 = d1.subtract(d2)
    assert exc.value.args[0] == 'Cannot propagate uncertainties of type StdDevUncertainty with uncertainties of type FakeUncertainty for subtraction'


@pytest.mark.parametrize('op1_unc,op2_unc', [
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
    op1 = NDDataArithmetic([1], uncertainty=op1_unc)
    op2 = NDDataArithmetic([1], uncertainty=op2_unc)

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


@pytest.mark.parametrize('op1_mask,op2_mask', [
                         (None, None),
                         (None, np.array([False])),
                         (np.array([False]), None),
                         (np.array([False]), np.array([False]))])
def test_arithmetic_result_not_tied_to_operands_mask(op1_mask, op2_mask):
    # See test_arithmetic_result_not_tied_to_operands_uncertainty for comments
    op1 = NDDataArithmetic([1], mask=op1_mask)
    op2 = NDDataArithmetic([1], mask=op2_mask)
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
    op1 = NDDataArithmetic([1], wcs=np.array([1]), unit='m')
    op2 = NDDataArithmetic([1], wcs=np.array([1]), unit='m')
    result = op1.add(op2)
    result.wcs[0] = 12345
    assert op1.wcs[0] != result.wcs[0]
    assert op2.wcs[0] != result.wcs[0]


# first operand has unit km, second has unit m
@pytest.mark.parametrize('operation,result_unit', [
                         ('add', u.km),
                         ('subtract', u.km),
                         ('multiply', u.km * u.m),
                         ('divide', u.km / u.m)])
def test_uncertainty_unit_conversion_add_subtract(operation, result_unit):
    in_km = NDDataArithmetic([1, 1], unit=u.km, uncertainty=StdDevUncertainty([.1, .1]))
    in_m = NDDataArithmetic(in_km.data * 1000, unit=u.m)
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


@pytest.mark.parametrize('unit1,unit2,op,result_unit', [
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
    ndd1 = NDDataArithmetic([1], unit=unit1)
    ndd2 = NDDataArithmetic([1], unit=unit2)
    ndd1_method = ndd1.__getattribute__(op)
    result = ndd1_method(ndd2)
    assert result.unit == result_unit


# check that subclasses can require wcs and/or unit to be present and use
# _arithmetic and convert_unit_to
class SubNDData(NDDataArithmetic):
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
