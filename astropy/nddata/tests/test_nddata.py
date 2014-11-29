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


@raises(ValueError)
def test_invalid_unit():
    d = NDData(np.ones((5, 5)), unit="NotAValidUnit")



def test_slicing_not_supported():
    ndd =  NDData(np.ones((5, 5)))
    with pytest.raises(TypeError):
        ndd[0]


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
