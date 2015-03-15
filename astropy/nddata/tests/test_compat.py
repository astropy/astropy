# Licensed under a 3-clause BSD style license - see LICENSE.rst
# This module contains tests of a class equivalent to pre-1.0 NDData.

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np

from ...tests.helper import pytest
from ..compat import NDDataArray
from ..nduncertainty import StdDevUncertainty
from ... import units as u


NDDATA_ATTRIBUTES = ['mask', 'flags', 'uncertainty', 'unit', 'shape', 'size',
                     'dtype', 'ndim', 'wcs', 'convert_unit_to']


def test_nddataarray_has_attributes_of_old_nddata():
    ndd = NDDataArray([1, 2, 3])
    for attr in NDDATA_ATTRIBUTES:
        assert hasattr(ndd, attr)


def test_nddata_simple():
    nd = NDDataArray(np.zeros((10, 10)))
    assert nd.shape == (10, 10)
    assert nd.size == 100
    assert nd.dtype == np.dtype(float)


def test_nddata_conversion():
    nd = NDDataArray(np.array([[1, 2, 3], [4, 5, 6]]))
    assert nd.size == 6
    assert nd.dtype == np.dtype(int)


@pytest.mark.parametrize('flags_in', [
                         np.array([True, False]),
                         np.array([1, 0]),
                         [True, False],
                         [1, 0],
                         np.array(['a', 'b']),
                         ['a', 'b']])
def test_nddata_flags_init_without_np_array(flags_in):
    ndd = NDDataArray([1, 1], flags=flags_in)
    assert (ndd.flags == flags_in).all()


@pytest.mark.parametrize(('shape'), [(10,), (5, 5), (3, 10, 10)])
def test_nddata_flags_invalid_shape(shape):
    with pytest.raises(ValueError) as exc:
        NDDataArray(np.zeros((10, 10)), flags=np.ones(shape))
    assert exc.value.args[0] == 'dimensions of flags do not match data'


def test_convert_unit_to():
    # convert_unit_to should return a copy of its input
    d = NDDataArray(np.ones((5, 5)))
    d.unit = 'km'
    d.uncertainty = StdDevUncertainty(0.1 + np.zeros_like(d))
    # workaround because zeros_like does not support dtype arg until v1.6
    # and NDData accepts only bool ndarray as mask
    tmp = np.zeros_like(d.data)
    d.mask = np.array(tmp, dtype=np.bool)
    d1 = d.convert_unit_to('m')
    assert np.all(d1.data == np.array(1000.0))
    assert np.all(d1.uncertainty.array == 1000.0 * d.uncertainty.array)
    assert d1.unit == u.m
    # changing the output mask should not change the original
    d1.mask[0, 0] = True
    assert d.mask[0, 0] != d1.mask[0, 0]
    d.flags = np.zeros_like(d.data)
    d1 = d.convert_unit_to('m')


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


def test_init_of_subclass_in_convert_unit_to():
    data = np.ones([10, 10])
    arr1 = SubNDData(data, unit='m', wcs=5)
    result = arr1.convert_unit_to('km')
    np.testing.assert_array_equal(arr1.data, 1000 * result.data)
