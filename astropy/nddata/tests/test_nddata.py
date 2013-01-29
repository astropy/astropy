# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np
from numpy.testing import assert_array_equal

from ..nddata import NDData
from ..nduncertainty import StdDevUncertainty, IncompatibleUncertaintiesException, NDUncertainty
from ...tests.helper import pytest, raises
from ...io import fits
from ...utils import NumpyRNGContext


class FakeUncertainty(NDUncertainty):

    def propagate_add(self, data, final_data):
        pass

    def propagate_subtract(self, data, final_data):
        pass

    def propagate_multiply(self, data, final_data):
        pass

    def propagate_divide(self, data, final_data):
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


def test_nddata_mask_valid():
    with NumpyRNGContext(456):
        NDData(np.random.random((10, 10)), mask=np.random.random((10, 10)) > 0.5)


@pytest.mark.parametrize(('shape'), [(10,), (5, 5), (3, 10, 10)])
def test_nddata_mask_invalid_shape(shape):
    with pytest.raises(ValueError) as exc:
        with NumpyRNGContext(789):
            NDData(np.random.random((10, 10)), mask=np.random.random(shape) > 0.5)
    assert exc.value.args[0] == 'dimensions of mask do not match data'


def test_nddata_uncertainty_init():
    u =  StdDevUncertainty(array=np.ones((5, 5)))
    d = NDData(np.ones((5, 5)), uncertainty=u)


def test_nddata_uncertainty_init_invalid_shape_1():
    u =  StdDevUncertainty(array=np.ones((6, 6)))
    with pytest.raises(ValueError) as exc:
        NDData(np.ones((5, 5)), uncertainty=u)
    assert exc.value.args[0] == 'parent shape does not match array data shape'


def test_nddata_uncertainty_init_invalid_shape_2():
    u =  StdDevUncertainty()
    NDData(np.ones((5, 5)), uncertainty=u)
    with pytest.raises(ValueError) as exc:
        u.array = np.ones((6, 6))
    assert exc.value.args[0] == 'array shape does not match parent data shape'


@pytest.mark.parametrize(('uncertainty'), [1., 'spam', np.ones((5, 5))])
def test_nddata_uncertainty_invalid_type(uncertainty):
    with pytest.raises(TypeError) as exc:
        NDData(np.ones((5, 5)), uncertainty=uncertainty)
    assert exc.value.args[0] == 'Uncertainty must be an instance of a NDUncertainty object'


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


def test_nddata_add():
    d1 = NDData(np.ones((5, 5)))
    d2 = NDData(np.ones((5, 5)))
    d3 = d1.add(d2)
    assert np.all(d3.data == 2.)


def test_nddata_add_mismatch_wcs():
    d1 = NDData(np.ones((5, 5)), wcs=1.)
    d2 = NDData(np.ones((5, 5)), wcs=2.)
    with pytest.raises(ValueError) as exc:
        d1.add(d2)
    assert exc.value.args[0] == "WCS properties do not match"


def test_nddata_add_mismatch_units():
    d1 = NDData(np.ones((5, 5)), unit='Jy')
    d2 = NDData(np.ones((5, 5)), unit='erg/s')
    with pytest.raises(ValueError) as exc:
        d1.add(d2)
    assert exc.value.args[0] == "operand units do not match"


def test_nddata_add_mismatch_shape():
    d1 = NDData(np.ones((5, 5)))
    d2 = NDData(np.ones((6, 6)))
    with pytest.raises(ValueError) as exc:
        d1.add(d2)
    assert exc.value.args[0] == "operand shapes do not match"


def test_nddata_add_uncertainties():
    u1 = StdDevUncertainty(array=np.ones((5, 5)) * 3)
    u2 = StdDevUncertainty(array=np.ones((5, 5)))
    d1 = NDData(np.ones((5, 5)), uncertainty=u1)
    d2 = NDData(np.ones((5, 5)), uncertainty=u2)
    d3 = d1.add(d2)
    assert np.all(d3.data == 2.)
    assert np.all(d3.uncertainty.array == np.sqrt(10.))


def test_nddata_add_uncertainties_mismatch():
    u1 = StdDevUncertainty(array=np.ones((5, 5)) * 3)
    u2 = FakeUncertainty()
    print u2.__class__
    d1 = NDData(np.ones((5, 5)), uncertainty=u1)
    d2 = NDData(np.ones((5, 5)), uncertainty=u2)
    with pytest.raises(IncompatibleUncertaintiesException) as exc:
        d3 = d1.add(d2)
    assert exc.value.args[0] == 'Cannot propagate uncertainties of type StdDevUncertainty with uncertainties of type FakeUncertainty for addition'


def test_nddata_subtract():
    d1 = NDData(np.ones((5, 5)))
    d2 = NDData(np.ones((5, 5)) * 2.)
    d3 = d1.subtract(d2)
    assert np.all(d3.data == -1.)


def test_nddata_subtract_mismatch_wcs():
    d1 = NDData(np.ones((5, 5)), wcs=1.)
    d2 = NDData(np.ones((5, 5)) * 2., wcs=2.)
    with pytest.raises(ValueError) as exc:
        d1.subtract(d2)
    assert exc.value.args[0] == "WCS properties do not match"


def test_nddata_subtract_mismatch_units():
    d1 = NDData(np.ones((5, 5)), unit='Jy')
    d2 = NDData(np.ones((5, 5)) * 2., unit='erg/s')
    with pytest.raises(ValueError) as exc:
        d1.subtract(d2)
    assert exc.value.args[0] == "operand units do not match"


def test_nddata_subtract_mismatch_shape():
    d1 = NDData(np.ones((5, 5)))
    d2 = NDData(np.ones((6, 6)) * 2.)
    with pytest.raises(ValueError) as exc:
        d1.subtract(d2)
    assert exc.value.args[0] == "operand shapes do not match"


def test_nddata_subtract_uncertainties():
    u1 = StdDevUncertainty(array=np.ones((5, 5)) * 3)
    u2 = StdDevUncertainty(array=np.ones((5, 5)))
    d1 = NDData(np.ones((5, 5)), uncertainty=u1)
    d2 = NDData(np.ones((5, 5)) * 2., uncertainty=u2)
    d3 = d1.subtract(d2)
    assert np.all(d3.data == -1.)
    assert np.all(d3.uncertainty.array == np.sqrt(10.))


def test_nddata_subtract_uncertainties_mismatch():
    u1 = StdDevUncertainty(array=np.ones((5, 5)) * 3)
    u2 = FakeUncertainty()
    print u2.__class__
    d1 = NDData(np.ones((5, 5)), uncertainty=u1)
    d2 = NDData(np.ones((5, 5)) * 2., uncertainty=u2)
    with pytest.raises(IncompatibleUncertaintiesException) as exc:
        d3 = d1.subtract(d2)
    assert exc.value.args[0] == 'Cannot propagate uncertainties of type StdDevUncertainty with uncertainties of type FakeUncertainty for subtraction'


def test_convert_unit_to():
    d = NDData(np.ones((5, 5)))
    d.unit = 'km'
    d1 = d.convert_unit_to('m')
    assert np.all(d1 == np.array(1000.0))


@raises(ValueError)
def test_invalid_unit():
    d = NDData(np.ones((5, 5)), unit="NotAValidUnit")

def test_simple_slicing():
    u1 = StdDevUncertainty(array=np.ones((5, 5)) * 3)
    d1 = NDData(np.ones((5, 5)), uncertainty=u1)
    assert d1.shape == (5,5)
    d2 = d1[2:3, 2:3]
    assert d2.shape == (1,1)

def test_slicing_reference():
    u1 = StdDevUncertainty(array=np.ones((5, 5)) * 3)
    d1 = NDData(np.ones((5, 5)), uncertainty=u1)
    d2 = d1[2:3, 2:3]
    #asserting that the new nddata contains references to the original nddata
    assert d2.data.base is d1.data
    assert d2.uncertainty.array.base is d1.uncertainty.array


def test_initializing_from_nddata():
    d1 = NDData(np.ones((5, 5)))
    d2 = NDData(d1)

    assert d1.data is d2.data


def test_initializing_from_nduncertainty():
    u1 = StdDevUncertainty(np.ones((5, 5)) * 3)
    u2 = StdDevUncertainty(u1, copy=False)

    assert u1.array is u2.array

def test_meta2ordered_dict():
    hdr = fits.header.Header()
    hdr.set('observer', 'Edwin Hubble')
    hdr.set('exptime', '3600')

    d1 = NDData(np.ones((5, 5)), meta=hdr)
    assert d1.meta['OBSERVER'] == 'Edwin Hubble'

@raises(TypeError)
def test_meta2ordered_dict_fail():
    hdr = 'this is not a valid header'
    d1 = NDData(np.ones((5, 5)), meta=hdr)

def test_masked_array_input():

    with NumpyRNGContext(12345):
        a = np.random.randn(100)
        marr = np.ma.masked_where(a > 0, a)

    nd = NDData(marr)

    #check that masks and data match
    assert_array_equal(nd.mask, marr.mask)
    assert_array_equal(nd.data, marr.data)

    #check that they are both by reference
    marr.mask[10] = ~marr.mask[10]
    marr.data[11] = 123456789

    assert_array_equal(nd.mask, marr.mask)
    assert_array_equal(nd.data, marr.data)
