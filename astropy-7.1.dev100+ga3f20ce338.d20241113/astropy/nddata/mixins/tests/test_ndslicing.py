# Licensed under a 3-clause BSD style license - see LICENSE.rst


import numpy as np
import pytest
from numpy.testing import assert_array_equal

from astropy import units as u
from astropy.nddata import NDData, NDSlicingMixin
from astropy.nddata import _testing as nd_testing
from astropy.nddata.nduncertainty import (
    NDUncertainty,
    StdDevUncertainty,
    UnknownUncertainty,
)


# Just add the Mixin to NDData
# TODO: Make this use NDDataRef instead!
class NDDataSliceable(NDSlicingMixin, NDData):
    pass


# Just some uncertainty (following the StdDevUncertainty implementation of
# storing the uncertainty in a property 'array') with slicing.
class SomeUncertainty(NDUncertainty):
    @property
    def uncertainty_type(self):
        return "fake"

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


@pytest.mark.parametrize("prop_name", ["mask", "uncertainty"])
def test_slicing_1dmask_ndslice(prop_name):
    # Data is 2d but mask/uncertainty only 1d so this should let the IndexError when
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
    uncertainty = StdDevUncertainty(np.linspace(10, 20, 10))
    naxis = 1
    wcs = nd_testing._create_wcs_simple(
        naxis=naxis,
        ctype=["deg"] * naxis,
        crpix=[3] * naxis,
        crval=[10] * naxis,
        cdelt=[1] * naxis,
    )
    # Just to have them too
    unit = u.s
    meta = {"observer": "Brian"}

    nd = NDDataSliceable(
        data, mask=mask, uncertainty=uncertainty, wcs=wcs, unit=unit, meta=meta
    )
    nd2 = nd[2:5]
    assert_array_equal(data[2:5], nd2.data)
    assert_array_equal(mask[2:5], nd2.mask)
    assert_array_equal(uncertainty[2:5].array, nd2.uncertainty.array)
    assert nd2.wcs.pixel_to_world(1) == nd.wcs.pixel_to_world(3)
    assert unit is nd2.unit
    assert meta == nd.meta


def test_slicing_all_npndarray_nd():
    # See what happens for multidimensional properties
    data = np.arange(1000).reshape(10, 10, 10)
    mask = data > 3
    uncertainty = np.linspace(10, 20, 1000).reshape(10, 10, 10)
    naxis = 3
    wcs = nd_testing._create_wcs_simple(
        naxis=naxis,
        ctype=["deg"] * naxis,
        crpix=[3] * naxis,
        crval=[10] * naxis,
        cdelt=[1] * naxis,
    )

    nd = NDDataSliceable(data, mask=mask, uncertainty=uncertainty, wcs=wcs)
    # Slice only 1D
    nd2 = nd[2:5]
    assert_array_equal(data[2:5], nd2.data)
    assert_array_equal(mask[2:5], nd2.mask)
    assert_array_equal(uncertainty[2:5], nd2.uncertainty.array)
    # Slice 3D
    nd2 = nd[2:5, :, 4:7]
    assert_array_equal(data[2:5, :, 4:7], nd2.data)
    assert_array_equal(mask[2:5, :, 4:7], nd2.mask)
    assert_array_equal(uncertainty[2:5, :, 4:7], nd2.uncertainty.array)
    assert nd2.wcs.pixel_to_world(1, 5, 1) == nd.wcs.pixel_to_world(5, 5, 3)


def test_slicing_all_npndarray_shape_diff():
    data = np.arange(10)
    mask = (data > 3)[0:9]
    uncertainty = np.linspace(10, 20, 15)
    naxis = 1
    wcs = nd_testing._create_wcs_simple(
        naxis=naxis,
        ctype=["deg"] * naxis,
        crpix=[3] * naxis,
        crval=[10] * naxis,
        cdelt=[1] * naxis,
    )

    nd = NDDataSliceable(data, mask=mask, uncertainty=uncertainty, wcs=wcs)
    nd2 = nd[2:5]
    assert_array_equal(data[2:5], nd2.data)
    # All are sliced even if the shapes differ (no Info)
    assert_array_equal(mask[2:5], nd2.mask)
    assert_array_equal(uncertainty[2:5], nd2.uncertainty.array)
    assert nd2.wcs.pixel_to_world(1) == nd.wcs.pixel_to_world(3)


def test_slicing_all_something_wrong():
    data = np.arange(10)
    mask = [False] * 10
    uncertainty = UnknownUncertainty({"rdnoise": 2.9, "gain": 1.4})
    naxis = 1
    wcs = nd_testing._create_wcs_simple(
        naxis=naxis,
        ctype=["deg"] * naxis,
        crpix=[3] * naxis,
        crval=[10] * naxis,
        cdelt=[1] * naxis,
    )

    nd = NDDataSliceable(data, mask=mask, uncertainty=uncertainty, wcs=wcs)
    nd2 = nd[2:5]
    # Sliced properties:
    assert_array_equal(data[2:5], nd2.data)
    assert_array_equal(mask[2:5], nd2.mask)
    # Not sliced attributes (they will raise a Info nevertheless)
    assert uncertainty.array == nd2.uncertainty.array
    assert uncertainty.uncertainty_type == nd2.uncertainty.uncertainty_type
    assert uncertainty.unit == nd2.uncertainty.unit
    assert nd2.wcs.pixel_to_world(1) == nd.wcs.pixel_to_world(3)


def test_boolean_slicing():
    data = np.arange(10)
    mask = data.copy()
    uncertainty = StdDevUncertainty(data.copy())
    naxis = 1
    wcs = nd_testing._create_wcs_simple(
        naxis=naxis,
        ctype=["deg"] * naxis,
        crpix=[3] * naxis,
        crval=[10] * naxis,
        cdelt=[1] * naxis,
    )

    nd = NDDataSliceable(data, mask=mask, uncertainty=uncertainty, wcs=wcs)

    with pytest.raises(ValueError):
        nd2 = nd[(nd.data >= 3) & (nd.data < 8)]

    nd.wcs = None
    nd2 = nd[(nd.data >= 3) & (nd.data < 8)]
    assert_array_equal(data[3:8], nd2.data)
    assert_array_equal(mask[3:8], nd2.mask)
