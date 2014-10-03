# Licensed under a 3-clause BSD style license - see LICENSE.rst

# TEST_UNICODE_LITERALS

from __future__ import absolute_import, division, print_function, unicode_literals

import textwrap

import numpy as np
from numpy.testing import assert_array_equal

from ..nddata import NDData
from ..nduncertainty import StdDevUncertainty, NDUncertainty
from ...tests.helper import pytest, raises
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


def test_nddata_empty():
    with pytest.raises(TypeError):
        NDData()  # empty initializer should fail


def test_nddata_simple():
    with NumpyRNGContext(123):
        nd = NDData(np.random.random((10, 10)))
    assert nd.shape == (10, 10)
    assert nd.size == 100
    assert nd.dtype == np.dtype(float)


def test_nddata_str():
    arr1d = NDData([1, 2, 3])
    assert str(arr1d) == '[1 2 3]'

    arr2d = NDData([[1, 2], [3, 4]])
    assert str(arr2d) == textwrap.dedent("""
        [[1 2]
         [3 4]]"""[1:])

    arr3d = NDData([[[1, 2], [3, 4]], [[5, 6], [7, 8]]])
    assert str(arr3d) == textwrap.dedent("""
        [[[1 2]
          [3 4]]

         [[5 6]
          [7 8]]]"""[1:])


def test_nddata_repr():
    arr1d = NDData([1, 2, 3])
    assert repr(arr1d) == 'NDData([1, 2, 3])'

    arr2d = NDData([[1, 2], [3, 4]])
    assert repr(arr2d) == textwrap.dedent("""
        NDData([[1, 2],
                [3, 4]])"""[1:])

    arr3d = NDData([[[1, 2], [3, 4]], [[5, 6], [7, 8]]])
    assert repr(arr3d) == textwrap.dedent("""
        NDData([[[1, 2],
                 [3, 4]],

                [[5, 6],
                 [7, 8]]])"""[1:])


def test_nddata_mask_valid():
    with NumpyRNGContext(456):
        NDData(np.random.random((10, 10)), mask=np.random.random((10, 10)) > 0.5)


@pytest.mark.parametrize('mask_in', [
                         np.array([True, False]),
                         np.array([1, 0]),
                         [True, False],
                         [1, 0]])
def test_nddata_mask_init_without_np_array(mask_in):
    ndd = NDData([1, 1], mask=mask_in)
    assert (ndd.mask == mask_in).all()


@pytest.mark.parametrize(('shape'), [(10,), (5, 5), (3, 10, 10)])
def test_nddata_mask_invalid_shape(shape):
    with pytest.raises(ValueError) as exc:
        with NumpyRNGContext(789):
            NDData(np.random.random((10, 10)), mask=np.random.random(shape) > 0.5)
    assert exc.value.args[0] == 'dimensions of mask do not match data'


def test_nddata_uncertainty_init():
    u = StdDevUncertainty(array=np.ones((5, 5)))
    d = NDData(np.ones((5, 5)), uncertainty=u)


def test_nddata_init_from_nddata_data_argument_only():
    ndd1 = NDData([1])
    ndd2 = NDData(ndd1)
    assert ndd2.wcs == ndd1.wcs
    assert ndd2.uncertainty == ndd1.uncertainty
    assert ndd2.mask == ndd1.mask
    assert ndd2.unit == ndd1.unit
    assert ndd2.meta == ndd1.meta


def test_nddata_copy_ref():
    """
    Tests to ensure that creating a new NDData object copies by *reference*.
    """
    a = np.ones((10, 10))
    nd_ref = NDData(a)
    a[0, 0] = 0
    assert nd_ref.data[0, 0] == 0


def test_nddata_conversion():
    nd = NDData([[1, 2, 3], [4, 5, 6]])
    assert nd.size == 6
    assert nd.dtype == np.dtype(int)


def test_ndddata_with_mask_acts_like_masked_array():
    # test for #2414
    input_mask = np.array([True, False, False])
    input_data = np.array([1, 2, 3])
    ndd_masked = NDData(input_data.copy(), mask=input_mask.copy())
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


def test_nddata_unmasked_in_operation_with_masked_numpy_array():
    # test for #2417
    ndd = NDData([1, 2, 3])
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


def test_convert_unit_to():
    # convert_unit_to should return a copy of its input
    d = NDData(np.ones((5, 5)))
    d.unit = 'km'
    d.uncertainty = StdDevUncertainty(0.1 + np.zeros_like(d))
    # workaround because zeros_like does not support dtype arg until v1.6
    # and NDData accepts only bool ndarray as mask
    tmp = np.zeros_like(d)
    d.mask = np.array(tmp, dtype=np.bool)
    d1 = d.convert_unit_to('m')
    assert np.all(d1 == np.array(1000.0))
    assert np.all(d1.uncertainty.array == 1000.0 * d.uncertainty.array)
    assert d1.unit == u.m
    # changing the output mask should not change the original
    d1.mask[0, 0] = True
    assert d.mask[0, 0] != d1.mask[0, 0]
    d.flags = np.zeros_like(d.data)
    d1 = d.convert_unit_to('m')


@raises(ValueError)
def test_invalid_unit():
    d = NDData(np.ones((5, 5)), unit="NotAValidUnit")


def test_simple_slicing():
    u1 = StdDevUncertainty(array=np.ones((5, 5)) * 3)
    d1 = NDData(np.ones((5, 5)), uncertainty=u1)
    assert d1.shape == (5, 5)
    d2 = d1[2:3, 2:3]
    assert d2.shape == (1, 1)
    d3 = d1[2, 2]
    assert d3.shape == ()


def test_slicing_reference():
    u1 = StdDevUncertainty(array=np.ones((5, 5)) * 3)
    d1 = NDData(np.ones((5, 5)), uncertainty=u1)
    d2 = d1[2:3, 2:3]
    # asserting that the new nddata contains references to the original nddata
    assert d2.data.base is d1.data
    assert d2.uncertainty.array.base is d1.uncertainty.array


def test_slicing_with_mask_or_flag():
    # Regression test for #2170
    ndd = NDData([1, 2, 3], mask=np.array([False, False, False]))
    assert ndd[0].shape == ()
    assert not ndd[0].mask


def test_initializing_from_nddata():
    d1 = NDData(np.ones((5, 5)))
    d2 = NDData(d1)

    assert d1.data is d2.data


def test_initializing_from_nduncertainty():
    u1 = StdDevUncertainty(np.ones((5, 5)) * 3)
    u2 = StdDevUncertainty(u1, copy=False)

    assert u1.array is u2.array


# Test an array and a scalar because a scalar Quantity does not always
# behaves the same way as an array.
@pytest.mark.parametrize('data', [np.array([1, 2, 3]), 5])
def test_initializing_nddata_from_quantity(data):
    # Until nddata and quantity are integrated initializing with a quantity
    # should raise an error.
    unit = u.adu
    ndd = NDData(data * unit)
    assert ndd.unit == unit
    np.testing.assert_array_equal(ndd.data, np.array(data))


def test_initializing_nddata_from_quantity_and_unit_raises_error():
    # Should raise an error if a Quantity is provided for the data and
    # an explicit unit is given.
    with pytest.raises(ValueError):
        NDData([1, 2, 3] * u.adu, unit=u.adu)


def test_masked_array_input():

    with NumpyRNGContext(12345):
        a = np.random.randn(100)
        marr = np.ma.masked_where(a > 0, a)

    nd = NDData(marr)

    # check that masks and data match
    assert_array_equal(nd.mask, marr.mask)
    assert_array_equal(nd.data, marr.data)

    # check that they are both by reference
    marr.mask[10] = ~marr.mask[10]
    marr.data[11] = 123456789

    assert_array_equal(nd.mask, marr.mask)
    assert_array_equal(nd.data, marr.data)


def test_unmasked_masked_array_input():
    # Test for #2784
    marr = np.ma.array([1,2,5]) # Masked array with no masked entries
    nd = NDData(marr) # Before fix this raised a ValueError

    # Check that masks are correct
    assert marr.mask is np.ma.nomask
    assert nd.mask is None  # Internal representation is np.ma.nomask but getter returns None


# Check that the meta descriptor is working as expected. The MetaBaseTest class
# takes care of defining all the tests, and we simply have to define the class
# and any minimal set of args to pass.

from ...utils.tests.test_metadata import MetaBaseTest


class TestMetaNDData(MetaBaseTest):
    test_class = NDData
    args = np.array([[1.]])


# check that subclasses can require wcs and/or unit to be present and use
# _arithmetic and convert_unit_to
class SubNDData(NDData):
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


def test_init_of_subclass_in_convert_unit_to():
    with NumpyRNGContext(12345):
        data = np.ones([10, 10])
    arr1 = SubNDData(data, unit='m', wcs=5)
    result = arr1.convert_unit_to('km')
    assert_array_equal(arr1.data, 1000 * result.data)
