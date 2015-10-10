# Licensed under a 3-clause BSD style license - see LICENSE.rst

# TEST_UNICODE_LITERALS

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np
from numpy.testing import assert_array_equal

from ... import NDData, NDSlicingMixin
from ...nduncertainty import NDUncertainty
from ....tests.helper import pytest
from .... import units as u

# Just add the Mixin to NDData
class NDDataSliceable(NDSlicingMixin, NDData):

    pass

# Just some uncertainty (following the StdDevUncertainty implementation of
# storing the uncertainty in a propery 'array') with slicing.
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


def test_slicing_data():
    data = np.arange(10)
    nd = NDDataSliceable(data)
    nd2 = nd[2:5]
    assert_array_equal(data[2:5], nd2.data)


def test_slicing_data_scalar_fail():
    data = np.array(10)
    nd = NDDataSliceable(data)
    with pytest.raises(TypeError): # as exc
        nd[:]
    # assert exc.value.args[0] == 'Scalars cannot be sliced.'


def test_slicing_all_npndarray():
    data = np.arange(10)
    mask = data > 3
    uncertainty = np.linspace(10,20,10)
    wcs = np.linspace(1,1000,10)
    # Just to have them too
    unit = u.s
    meta = {'observer': 'Brian'}

    nd = NDDataSliceable(data, mask=mask, uncertainty=uncertainty, wcs=wcs,
                         unit=unit, meta=meta)
    nd2 = nd[2:5]
    assert_array_equal(data[2:5], nd2.data)
    assert_array_equal(mask[2:5], nd2.mask)
    assert_array_equal(uncertainty[2:5], nd2.uncertainty)
    assert_array_equal(wcs[2:5], nd2.wcs)
    assert unit is nd2.unit
    assert meta == nd.meta

    # See what happens for multidimensional properties
    data = np.arange(1000).reshape(10,10,10)
    mask = data > 3
    uncertainty = np.linspace(10,20,1000).reshape(10,10,10)
    wcs = np.linspace(1,1000,1000).reshape(10,10,10)
    # Just to have them too
    unit = u.s
    meta = {'observer': 'Brian'}

    nd = NDDataSliceable(data, mask=mask, uncertainty=uncertainty, wcs=wcs,
                         unit=unit, meta=meta)
    # Slice only 1D
    nd2 = nd[2:5]
    assert_array_equal(data[2:5], nd2.data)
    assert_array_equal(mask[2:5], nd2.mask)
    assert_array_equal(uncertainty[2:5], nd2.uncertainty)
    assert_array_equal(wcs[2:5], nd2.wcs)
    assert unit is nd2.unit
    assert meta == nd.meta
    # Slice 3D
    nd2 = nd[2:5,:,4:7]
    assert_array_equal(data[2:5,:,4:7], nd2.data)
    assert_array_equal(mask[2:5,:,4:7], nd2.mask)
    assert_array_equal(uncertainty[2:5,:,4:7], nd2.uncertainty)
    assert_array_equal(wcs[2:5,:,4:7], nd2.wcs)
    assert unit is nd2.unit
    assert meta == nd.meta


def test_slicing_all_npndarray_shape_diff():
    data = np.arange(10)
    mask = (data > 3)[0:9]
    uncertainty = np.linspace(10,20,15)
    wcs = np.linspace(1,1000,12)
    unit = u.s
    meta = {'observer': 'Brian'}

    nd = NDDataSliceable(data, mask=mask, uncertainty=uncertainty, wcs=wcs,
                         unit=unit, meta=meta)
    nd2 = nd[2:5]
    assert_array_equal(data[2:5], nd2.data)
    # All other properties should not be sliced since their shape differs!
    assert mask is nd2.mask
    assert uncertainty is nd2.uncertainty
    assert wcs is nd2.wcs
    # But unit and meta remain the same
    assert unit is nd2.unit
    assert meta == nd.meta

    # See what happens for multidimensional properties
    data = np.arange(1000).reshape(10,10,10)
    mask = (data > 3)[0:9,:,:]
    uncertainty = np.linspace(10,20,1200).reshape(12,10,10)
    wcs = np.linspace(1,1000,700).reshape(10,7,10)
    # Just to have them too
    unit = u.s
    meta = {'observer': 'Brian'}

    nd = NDDataSliceable(data, mask=mask, uncertainty=uncertainty, wcs=wcs,
                         unit=unit, meta=meta)
    # Slice only 1D
    nd2 = nd[2:5]
    assert_array_equal(data[2:5], nd2.data)
    # All other properties should not be sliced since their shape differs!
    assert mask is nd2.mask
    assert uncertainty is nd2.uncertainty
    assert wcs is nd2.wcs
    # But unit and meta remain the same
    assert unit is nd2.unit
    assert meta == nd.meta
    # Slice 3D
    nd2 = nd[2:5,:,4:7]
    assert_array_equal(data[2:5,:,4:7], nd2.data)
    # All other properties should not be sliced since their shape differs!
    assert mask is nd2.mask
    assert uncertainty is nd2.uncertainty
    assert wcs is nd2.wcs
    # But unit and meta remain the same
    assert unit is nd2.unit
    assert meta == nd.meta


def test_slicing_uncertainty_is_nduncertainty():
    data = np.arange(10)
    uncertainty = SomeUncertainty(np.linspace(10,20,10))

    nd = NDDataSliceable(data, uncertainty=uncertainty)
    nd2 = nd[2:5]
    assert_array_equal(data[2:5], nd2.data)
    assert_array_equal(uncertainty.array[2:5], nd2.uncertainty.array)

    # Check the NDUncertainty with a wrong shape is not considered sliceable
    uncertainty_wrong = SomeUncertainty(np.linspace(10,20,15))
    nd = NDDataSliceable(data, uncertainty=uncertainty_wrong)
    nd2 = nd[2:5]
    assert_array_equal(data[2:5], nd2.data)
    # not sliced
    assert_array_equal(uncertainty_wrong.array, nd2.uncertainty.array)


def test_slicing_wcs_is_WCS():
    # TODO: it's not that easy to make such a test, I think it's best to
    # rely on astropy.wcs.WCS tests.
    pass


def test_slicing_all_something_wrong():
    data = np.arange(10)
    mask = [False]*10
    uncertainty = {'rdnoise': 2.9, 'gain': 1.4}
    wcs = 145 * u.degree
    unit = u.s
    meta = {'observer': 'Brian'}

    nd = NDDataSliceable(data, mask=mask, uncertainty=uncertainty, wcs=wcs,
                         unit=unit, meta=meta)
    nd2 = nd[2:5]
    assert_array_equal(data[2:5], nd2.data)
    # All other properties should not be sliced since their types are
    # are considered not sliceable
    assert_array_equal(mask, nd2.mask)
    assert_array_equal(uncertainty, nd2.uncertainty)
    assert_array_equal(wcs, nd2.wcs)
    # But unit and meta remain the same
    assert unit is nd2.unit
    assert meta == nd.meta
