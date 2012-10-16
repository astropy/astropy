# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest
import numpy as np

from ..nddata import NDData
from ..nderror import StandardDeviationError, IncompatibleErrorsException, NDError
from ...tests.helper import raises
from ...io import fits

np.random.seed(12345)


class FakeErrors(NDError):

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
    nd = NDData(np.random.random((10, 10)))
    assert nd.shape == (10, 10)
    assert nd.size == 100
    assert nd.dtype == np.dtype(float)


def test_nddata_mask_valid():
    NDData(np.random.random((10, 10)), mask=np.random.random((10, 10)) > 0.5)


@pytest.mark.parametrize(('shape'), [(10,), (5, 5), (3, 10, 10)])
def test_nddata_mask_invalid_shape(shape):
    with pytest.raises(ValueError) as exc:
        NDData(np.random.random((10, 10)), mask=np.random.random(shape) > 0.5)
    assert exc.value.args[0] == 'dimensions of mask do not match data'


def test_nddata_error_init():
    e = StandardDeviationError(array=np.ones((5, 5)))
    d = NDData(np.ones((5, 5)), error=e)


def test_nddata_error_init_invalid_shape_1():
    e = StandardDeviationError(array=np.ones((6, 6)))
    with pytest.raises(ValueError) as exc:
        NDData(np.ones((5, 5)), error=e)
    assert exc.value.args[0] == 'parent shape does not match array data shape'


def test_nddata_error_init_invalid_shape_2():
    e = StandardDeviationError()
    NDData(np.ones((5, 5)), error=e)
    with pytest.raises(ValueError) as exc:
        e.array = np.ones((6, 6))
    assert exc.value.args[0] == 'array shape does not match parent data shape'


@pytest.mark.parametrize(('error'), [1., 'spam', np.ones((5, 5))])
def test_nddata_error_invalid_type(error):
    with pytest.raises(TypeError) as exc:
        NDData(np.ones((5, 5)), error=error)
    assert exc.value.args[0] == 'error must be an instance of a NDError object'


def test_nddata_copy():
    a = np.ones((10, 10))
    nd_copy = NDData(a, copy=True)
    nd_ref = NDData(a, copy=False)
    a[0, 0] = 0
    assert nd_ref.data[0, 0] == 0
    assert nd_copy.data[0, 0] == 1


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
    d1 = NDData(np.ones((5, 5)), units='Jy')
    d2 = NDData(np.ones((5, 5)), units='erg/s')
    with pytest.raises(ValueError) as exc:
        d1.add(d2)
    assert exc.value.args[0] == "operand units do not match"


def test_nddata_add_mismatch_shape():
    d1 = NDData(np.ones((5, 5)))
    d2 = NDData(np.ones((6, 6)))
    with pytest.raises(ValueError) as exc:
        d1.add(d2)
    assert exc.value.args[0] == "operand shapes do not match"


def test_nddata_add_errors():
    e1 = StandardDeviationError(array=np.ones((5, 5)) * 3)
    e2 = StandardDeviationError(array=np.ones((5, 5)))
    d1 = NDData(np.ones((5, 5)), error=e1)
    d2 = NDData(np.ones((5, 5)), error=e2)
    d3 = d1.add(d2)
    assert np.all(d3.data == 2.)
    assert np.all(d3.error.array == np.sqrt(10.))


def test_nddata_add_errors_mismatch():
    e1 = StandardDeviationError(array=np.ones((5, 5)) * 3)
    e2 = FakeErrors()
    print e2.__class__
    d1 = NDData(np.ones((5, 5)), error=e1)
    d2 = NDData(np.ones((5, 5)), error=e2)
    with pytest.raises(IncompatibleErrorsException) as exc:
        d3 = d1.add(d2)
    assert exc.value.args[0] == 'Cannot propagate errors of type StandardDeviationError with errors of type FakeErrors for addition'


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
    d1 = NDData(np.ones((5, 5)), units='Jy')
    d2 = NDData(np.ones((5, 5)) * 2., units='erg/s')
    with pytest.raises(ValueError) as exc:
        d1.subtract(d2)
    assert exc.value.args[0] == "operand units do not match"


def test_nddata_subtract_mismatch_shape():
    d1 = NDData(np.ones((5, 5)))
    d2 = NDData(np.ones((6, 6)) * 2.)
    with pytest.raises(ValueError) as exc:
        d1.subtract(d2)
    assert exc.value.args[0] == "operand shapes do not match"


def test_nddata_subtract_errors():
    e1 = StandardDeviationError(array=np.ones((5, 5)) * 3)
    e2 = StandardDeviationError(array=np.ones((5, 5)))
    d1 = NDData(np.ones((5, 5)), error=e1)
    d2 = NDData(np.ones((5, 5)) * 2., error=e2)
    d3 = d1.subtract(d2)
    assert np.all(d3.data == -1.)
    assert np.all(d3.error.array == np.sqrt(10.))


def test_nddata_subtract_errors_mismatch():
    e1 = StandardDeviationError(array=np.ones((5, 5)) * 3)
    e2 = FakeErrors()
    print e2.__class__
    d1 = NDData(np.ones((5, 5)), error=e1)
    d2 = NDData(np.ones((5, 5)) * 2., error=e2)
    with pytest.raises(IncompatibleErrorsException) as exc:
        d3 = d1.subtract(d2)
    assert exc.value.args[0] == 'Cannot propagate errors of type StandardDeviationError with errors of type FakeErrors for subtraction'


def test_convert_units_to():
    d = NDData(np.ones((5, 5)))
    d.units = 'km'
    d1 = d.convert_units_to('m')
    assert np.all(d1 == np.array(1000.0))


@raises(ValueError)
def test_invalid_unit():
    d = NDData(np.ones((5, 5)), units="NotAValidUnit")

def test_simple_slicing():
    e1 = StandardDeviationError(array=np.ones((5, 5)) * 3)
    d1 = NDData(np.ones((5, 5)), error=e1)
    assert d1.shape == (5,5)
    d2 = d1[2:3, 2:3]
    assert d2.shape == (1,1)

def test_slicing_reference():
    e1 = StandardDeviationError(array=np.ones((5, 5)) * 3)
    d1 = NDData(np.ones((5, 5)), error=e1)
    d2 = d1[2:3, 2:3]
    #asserting that the new nddata contains references to the original nddata
    assert d2.data.base is d1.data
    assert d2.error.array.base is d1.error.array


def test_initializing_from_nddata():
    d1 = NDData(np.ones((5, 5)))
    d2 = NDData(d1, copy=False)

    assert d1.data is d2.data


def test_initializing_from_nderror():
    e1 = StandardDeviationError(np.ones((5, 5)) * 3)
    e2 = StandardDeviationError(e1, copy=False)

    assert e1.array is e2.array

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
