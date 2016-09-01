# Licensed under a 3-clause BSD style license - see LICENSE.rst

# TEST_UNICODE_LITERALS

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np
from numpy.testing import assert_array_equal

from ... import NDData, NDSlicingMixin
from ...nduncertainty import NDUncertainty, StdDevUncertainty
from ....tests.helper import pytest
from .... import units as u


# Just add the Mixin to NDData
# TODO: Make this use NDDataRef instead!
class NDDataSliceable(NDSlicingMixin, NDData):
    pass


# Just some uncertainty (following the StdDevUncertainty implementation of
# storing the uncertainty in a property 'array') with slicing.
class SomeUncertainty(NDUncertainty):

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


def test_slicing_only_data():
    data = np.arange(10)
    nd = NDDataSliceable(data)
    nd2 = nd[2:5]
    assert_array_equal(data[2:5], nd2.data)


def test_slicing_data_scalar_fail():
    data = np.array(10)
    nd = NDDataSliceable(data)
    with pytest.raises(TypeError):  # as exc
        nd[:]
    # assert exc.value.args[0] == 'Scalars cannot be sliced.'


def test_slicing_1ddata_ndslice():
    data = np.array([10, 20])
    nd = NDDataSliceable(data)
    # Standard numpy warning here:
    with pytest.raises(IndexError):
        nd[:, :]


@pytest.mark.parametrize('prop_name', ['mask', 'wcs', 'uncertainty'])
def test_slicing_1dmask_ndslice(prop_name):
    # Data is 2d but mask only 1d so this should let the IndexError when
    # slicing the mask rise to the user.
    data = np.ones((3, 3))
    kwarg = {prop_name: np.ones(3)}
    nd = NDDataSliceable(data, **kwarg)
    # Standard numpy warning here:
    with pytest.raises(IndexError):
        nd[:, :]


def test_slicing_all_npndarray_1d():
    data = np.arange(10)
    mask = data > 3
    uncertainty = np.linspace(10, 20, 10)
    wcs = np.linspace(1, 1000, 10)
    # Just to have them too
    unit = u.s
    meta = {'observer': 'Brian'}

    nd = NDDataSliceable(data, mask=mask, uncertainty=uncertainty, wcs=wcs,
                         unit=unit, meta=meta)
    nd2 = nd[2:5]
    assert_array_equal(data[2:5], nd2.data)
    assert_array_equal(mask[2:5], nd2.mask)
    assert_array_equal(uncertainty[2:5], nd2.uncertainty.array)
    assert_array_equal(wcs[2:5], nd2.wcs)
    assert unit is nd2.unit
    assert meta == nd.meta


def test_slicing_all_npndarray_nd():
    # See what happens for multidimensional properties
    data = np.arange(1000).reshape(10, 10, 10)
    mask = data > 3
    uncertainty = np.linspace(10, 20, 1000).reshape(10, 10, 10)
    wcs = np.linspace(1, 1000, 1000).reshape(10, 10, 10)

    nd = NDDataSliceable(data, mask=mask, uncertainty=uncertainty, wcs=wcs)
    # Slice only 1D
    nd2 = nd[2:5]
    assert_array_equal(data[2:5], nd2.data)
    assert_array_equal(mask[2:5], nd2.mask)
    assert_array_equal(uncertainty[2:5], nd2.uncertainty.array)
    assert_array_equal(wcs[2:5], nd2.wcs)
    # Slice 3D
    nd2 = nd[2:5, :, 4:7]
    assert_array_equal(data[2:5, :, 4:7], nd2.data)
    assert_array_equal(mask[2:5, :, 4:7], nd2.mask)
    assert_array_equal(uncertainty[2:5, :, 4:7], nd2.uncertainty.array)
    assert_array_equal(wcs[2:5, :, 4:7], nd2.wcs)


def test_slicing_all_npndarray_shape_diff():
    data = np.arange(10)
    mask = (data > 3)[0:9]
    uncertainty = np.linspace(10, 20, 15)
    wcs = np.linspace(1, 1000, 12)

    nd = NDDataSliceable(data, mask=mask, uncertainty=uncertainty, wcs=wcs)
    nd2 = nd[2:5]
    assert_array_equal(data[2:5], nd2.data)
    # All are sliced even if the shapes differ (no Info)
    assert_array_equal(mask[2:5], nd2.mask)
    assert_array_equal(uncertainty[2:5], nd2.uncertainty.array)
    assert_array_equal(wcs[2:5], nd2.wcs)


def test_slicing_all_something_wrong():
    data = np.arange(10)
    mask = [False]*10
    uncertainty = {'rdnoise': 2.9, 'gain': 1.4}
    wcs = 145 * u.degree

    nd = NDDataSliceable(data, mask=mask, uncertainty=uncertainty, wcs=wcs)
    nd2 = nd[2:5]
    # Sliced properties:
    assert_array_equal(data[2:5], nd2.data)
    assert_array_equal(mask[2:5], nd2.mask)
    # Not sliced attributes (they will raise a Info nevertheless)
    uncertainty is nd2.uncertainty
    assert_array_equal(wcs, nd2.wcs)


def test_boolean_slicing():
    data = np.arange(10)
    mask = data.copy()
    uncertainty = StdDevUncertainty(data.copy())
    wcs = data.copy()
    nd = NDDataSliceable(data, mask=mask, uncertainty=uncertainty, wcs=wcs)

    nd2 = nd[(nd.data >= 3) & (nd.data < 8)]
    assert_array_equal(data[3:8], nd2.data)
    assert_array_equal(mask[3:8], nd2.mask)
    assert_array_equal(wcs[3:8], nd2.wcs)
    assert_array_equal(uncertainty.array[3:8], nd2.uncertainty.array)
