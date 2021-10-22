# Licensed under a 3-clause BSD style license - see LICENSE.rst


import pickle
import textwrap
from collections import OrderedDict

import pytest
import numpy as np
from numpy.testing import assert_array_equal

from astropy.nddata.nddata import NDData
from astropy.nddata.nduncertainty import StdDevUncertainty
from astropy import units as u
from astropy.utils import NumpyRNGContext
from astropy.wcs import WCS
from astropy.wcs.wcsapi import HighLevelWCSWrapper, SlicedLowLevelWCS, \
                               BaseHighLevelWCS

from .test_nduncertainty import FakeUncertainty
from astropy.nddata import _testing as nd_testing


class FakeNumpyArray:
    """
    Class that has a few of the attributes of a numpy array.

    These attributes are checked for by NDData.
    """
    def __init__(self):
        super().__init__()

    def shape(self):
        pass

    def __getitem__(self):
        pass

    def __array__(self):
        pass

    @property
    def dtype(self):
        return 'fake'


class MinimalUncertainty:
    """
    Define the minimum attributes acceptable as an uncertainty object.
    """
    def __init__(self, value):
        self._uncertainty = value

    @property
    def uncertainty_type(self):
        return "totally and completely fake"


class BadNDDataSubclass(NDData):

    def __init__(self, data, uncertainty=None, mask=None, wcs=None,
                 meta=None, unit=None):
        self._data = data
        self._uncertainty = uncertainty
        self._mask = mask
        self._wcs = wcs
        self._unit = unit
        self._meta = meta


# Setter tests
def test_uncertainty_setter():
    nd = NDData([1, 2, 3])
    good_uncertainty = MinimalUncertainty(5)
    nd.uncertainty = good_uncertainty
    assert nd.uncertainty is good_uncertainty
    # Check the fake uncertainty (minimal does not work since it has no
    # parent_nddata attribute from NDUncertainty)
    nd.uncertainty = FakeUncertainty(5)
    assert nd.uncertainty.parent_nddata is nd
    # Check that it works if the uncertainty was set during init
    nd = NDData(nd)
    assert isinstance(nd.uncertainty, FakeUncertainty)
    nd.uncertainty = 10
    assert not isinstance(nd.uncertainty, FakeUncertainty)
    assert nd.uncertainty.array == 10


def test_mask_setter():
    # Since it just changes the _mask attribute everything should work
    nd = NDData([1, 2, 3])
    nd.mask = True
    assert nd.mask
    nd.mask = False
    assert not nd.mask
    # Check that it replaces a mask from init
    nd = NDData(nd, mask=True)
    assert nd.mask
    nd.mask = False
    assert not nd.mask


# Init tests
def test_nddata_empty():
    with pytest.raises(TypeError):
        NDData()  # empty initializer should fail


def test_nddata_init_data_nonarray():
    inp = [1, 2, 3]
    nd = NDData(inp)
    assert (np.array(inp) == nd.data).all()


def test_nddata_init_data_ndarray():
    # random floats
    with NumpyRNGContext(123):
        nd = NDData(np.random.random((10, 10)))
    assert nd.data.shape == (10, 10)
    assert nd.data.size == 100
    assert nd.data.dtype == np.dtype(float)

    # specific integers
    nd = NDData(np.array([[1, 2, 3], [4, 5, 6]]))
    assert nd.data.size == 6
    assert nd.data.dtype == np.dtype(int)

    # Tests to ensure that creating a new NDData object copies by *reference*.
    a = np.ones((10, 10))
    nd_ref = NDData(a)
    a[0, 0] = 0
    assert nd_ref.data[0, 0] == 0

    # Except we choose copy=True
    a = np.ones((10, 10))
    nd_ref = NDData(a, copy=True)
    a[0, 0] = 0
    assert nd_ref.data[0, 0] != 0


def test_nddata_init_data_maskedarray():
    with NumpyRNGContext(456):
        NDData(np.random.random((10, 10)),
               mask=np.random.random((10, 10)) > 0.5)

    # Another test (just copied here)
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

    # or not if we choose copy=True
    nd = NDData(marr, copy=True)
    marr.mask[10] = ~marr.mask[10]
    marr.data[11] = 0
    assert nd.mask[10] != marr.mask[10]
    assert nd.data[11] != marr.data[11]


@pytest.mark.parametrize('data', [np.array([1, 2, 3]), 5])
def test_nddata_init_data_quantity(data):
    # Test an array and a scalar because a scalar Quantity does not always
    # behaves the same way as an array.
    quantity = data * u.adu
    ndd = NDData(quantity)
    assert ndd.unit == quantity.unit
    assert_array_equal(ndd.data, np.array(quantity.value))
    if ndd.data.size > 1:
        # check that if it is an array it is not copied
        quantity.value[1] = 100
        assert ndd.data[1] == quantity.value[1]

        # or is copied if we choose copy=True
        ndd = NDData(quantity, copy=True)
        quantity.value[1] = 5
        assert ndd.data[1] != quantity.value[1]


def test_nddata_init_data_masked_quantity():
    a = np.array([2, 3])
    q = a * u.m
    m = False
    mq = np.ma.array(q, mask=m)
    nd = NDData(mq)
    assert_array_equal(nd.data, a)
    # This test failed before the change in nddata init because the masked
    # arrays data (which in fact was a quantity was directly saved)
    assert nd.unit == u.m
    assert not isinstance(nd.data, u.Quantity)
    np.testing.assert_array_equal(nd.mask, np.array(m))


def test_nddata_init_data_nddata():
    nd1 = NDData(np.array([1]))
    nd2 = NDData(nd1)
    assert nd2.wcs == nd1.wcs
    assert nd2.uncertainty == nd1.uncertainty
    assert nd2.mask == nd1.mask
    assert nd2.unit == nd1.unit
    assert nd2.meta == nd1.meta

    # Check that it is copied by reference
    nd1 = NDData(np.ones((5, 5)))
    nd2 = NDData(nd1)
    assert nd1.data is nd2.data

    # Check that it is really copied if copy=True
    nd2 = NDData(nd1, copy=True)
    nd1.data[2, 3] = 10
    assert nd1.data[2, 3] != nd2.data[2, 3]

    # Now let's see what happens if we have all explicitly set
    nd1 = NDData(np.array([1]), mask=False, uncertainty=StdDevUncertainty(10), unit=u.s,
                 meta={'dest': 'mordor'}, wcs=WCS(naxis=1))
    nd2 = NDData(nd1)
    assert nd2.data is nd1.data
    assert nd2.wcs is nd1.wcs
    assert nd2.uncertainty.array == nd1.uncertainty.array
    assert nd2.mask == nd1.mask
    assert nd2.unit == nd1.unit
    assert nd2.meta == nd1.meta

    # now what happens if we overwrite them all too
    nd3 = NDData(nd1, mask=True, uncertainty=StdDevUncertainty(200), unit=u.km,
                 meta={'observer': 'ME'}, wcs=WCS(naxis=1))
    assert nd3.data is nd1.data
    assert nd3.wcs is not nd1.wcs
    assert nd3.uncertainty.array != nd1.uncertainty.array
    assert nd3.mask != nd1.mask
    assert nd3.unit != nd1.unit
    assert nd3.meta != nd1.meta


def test_nddata_init_data_nddata_subclass():
    uncert = StdDevUncertainty(3)
    # There might be some incompatible subclasses of NDData around.
    bnd = BadNDDataSubclass(False, True, 3, 2, 'gollum', 100)
    # Before changing the NDData init this would not have raised an error but
    # would have lead to a compromised nddata instance
    with pytest.raises(TypeError):
        NDData(bnd)

    # but if it has no actual incompatible attributes it passes
    bnd_good = BadNDDataSubclass(np.array([1, 2]), uncert, 3, HighLevelWCSWrapper(WCS(naxis=1)),
                                 {'enemy': 'black knight'}, u.km)
    nd = NDData(bnd_good)
    assert nd.unit == bnd_good.unit
    assert nd.meta == bnd_good.meta
    assert nd.uncertainty == bnd_good.uncertainty
    assert nd.mask == bnd_good.mask
    assert nd.wcs is bnd_good.wcs
    assert nd.data is bnd_good.data


def test_nddata_init_data_fail():
    # First one is sliceable but has no shape, so should fail.
    with pytest.raises(TypeError):
        NDData({'a': 'dict'})

    # This has a shape but is not sliceable
    class Shape:
        def __init__(self):
            self.shape = 5

        def __repr__(self):
            return '7'

    with pytest.raises(TypeError):
        NDData(Shape())


def test_nddata_init_data_fakes():
    ndd1 = NDData(FakeNumpyArray())
    # First make sure that NDData isn't converting its data to a numpy array.
    assert isinstance(ndd1.data, FakeNumpyArray)
    # Make a new NDData initialized from an NDData
    ndd2 = NDData(ndd1)
    # Check that the data wasn't converted to numpy
    assert isinstance(ndd2.data, FakeNumpyArray)


# Specific parameters
def test_param_uncertainty():
    u = StdDevUncertainty(array=np.ones((5, 5)))
    d = NDData(np.ones((5, 5)), uncertainty=u)
    # Test that the parent_nddata is set.
    assert d.uncertainty.parent_nddata is d
    # Test conflicting uncertainties (other NDData)
    u2 = StdDevUncertainty(array=np.ones((5, 5))*2)
    d2 = NDData(d, uncertainty=u2)
    assert d2.uncertainty is u2
    assert d2.uncertainty.parent_nddata is d2


def test_param_wcs():
    # Since everything is allowed we only need to test something
    nd = NDData([1], wcs=WCS(naxis=1))
    assert nd.wcs is not None
    # Test conflicting wcs (other NDData)
    nd2 = NDData(nd, wcs=WCS(naxis=1))
    assert nd2.wcs is not None and nd2.wcs is not nd.wcs


def test_param_meta():
    # everything dict-like is allowed
    with pytest.raises(TypeError):
        NDData([1], meta=3)
    nd = NDData([1, 2, 3], meta={})
    assert len(nd.meta) == 0
    nd = NDData([1, 2, 3])
    assert isinstance(nd.meta, OrderedDict)
    assert len(nd.meta) == 0
    # Test conflicting meta (other NDData)
    nd2 = NDData(nd, meta={'image': 'sun'})
    assert len(nd2.meta) == 1
    nd3 = NDData(nd2, meta={'image': 'moon'})
    assert len(nd3.meta) == 1
    assert nd3.meta['image'] == 'moon'


def test_param_mask():
    # Since everything is allowed we only need to test something
    nd = NDData([1], mask=False)
    assert not nd.mask
    # Test conflicting mask (other NDData)
    nd2 = NDData(nd, mask=True)
    assert nd2.mask
    # (masked array)
    nd3 = NDData(np.ma.array([1], mask=False), mask=True)
    assert nd3.mask
    # (masked quantity)
    mq = np.ma.array(np.array([2, 3])*u.m, mask=False)
    nd4 = NDData(mq, mask=True)
    assert nd4.mask


def test_param_unit():
    with pytest.raises(ValueError):
        NDData(np.ones((5, 5)), unit="NotAValidUnit")
    NDData([1, 2, 3], unit='meter')
    # Test conflicting units (quantity as data)
    q = np.array([1, 2, 3]) * u.m
    nd = NDData(q, unit='cm')
    assert nd.unit != q.unit
    assert nd.unit == u.cm
    # (masked quantity)
    mq = np.ma.array(np.array([2, 3])*u.m, mask=False)
    nd2 = NDData(mq, unit=u.s)
    assert nd2.unit == u.s
    # (another NDData as data)
    nd3 = NDData(nd, unit='km')
    assert nd3.unit == u.km


def test_pickle_nddata_with_uncertainty():
    ndd = NDData(np.ones(3),
                 uncertainty=StdDevUncertainty(np.ones(5), unit=u.m),
                 unit=u.m)
    ndd_dumped = pickle.dumps(ndd)
    ndd_restored = pickle.loads(ndd_dumped)
    assert type(ndd_restored.uncertainty) is StdDevUncertainty
    assert ndd_restored.uncertainty.parent_nddata is ndd_restored
    assert ndd_restored.uncertainty.unit == u.m


def test_pickle_uncertainty_only():
    ndd = NDData(np.ones(3),
                 uncertainty=StdDevUncertainty(np.ones(5), unit=u.m),
                 unit=u.m)
    uncertainty_dumped = pickle.dumps(ndd.uncertainty)
    uncertainty_restored = pickle.loads(uncertainty_dumped)
    np.testing.assert_array_equal(ndd.uncertainty.array,
                                  uncertainty_restored.array)
    assert ndd.uncertainty.unit == uncertainty_restored.unit
    # Even though it has a parent there is no one that references the parent
    # after unpickling so the weakref "dies" immediately after unpickling
    # finishes.
    assert uncertainty_restored.parent_nddata is None


def test_pickle_nddata_without_uncertainty():
    ndd = NDData(np.ones(3), unit=u.m)
    dumped = pickle.dumps(ndd)
    ndd_restored = pickle.loads(dumped)
    np.testing.assert_array_equal(ndd.data, ndd_restored.data)


# Check that the meta descriptor is working as expected. The MetaBaseTest class
# takes care of defining all the tests, and we simply have to define the class
# and any minimal set of args to pass.
from astropy.utils.tests.test_metadata import MetaBaseTest


class TestMetaNDData(MetaBaseTest):
    test_class = NDData
    args = np.array([[1.]])


# Representation tests
def test_nddata_str():
    arr1d = NDData(np.array([1, 2, 3]))
    assert str(arr1d) == '[1 2 3]'

    arr2d = NDData(np.array([[1, 2], [3, 4]]))
    assert str(arr2d) == textwrap.dedent("""
        [[1 2]
         [3 4]]"""[1:])

    arr3d = NDData(np.array([[[1, 2], [3, 4]], [[5, 6], [7, 8]]]))
    assert str(arr3d) == textwrap.dedent("""
        [[[1 2]
          [3 4]]

         [[5 6]
          [7 8]]]"""[1:])

    # let's add units!
    arr = NDData(np.array([1, 2, 3]), unit="km")
    assert str(arr) == '[1 2 3] km'

    # what if it had these units?
    arr = NDData(np.array([1, 2, 3]), unit="erg cm^-2 s^-1 A^-1")
    assert str(arr) == '[1 2 3] erg / (A cm2 s)'


def test_nddata_repr():
    # The big test is eval(repr()) should be equal to the original!

    arr1d = NDData(np.array([1, 2, 3]))
    s = repr(arr1d)
    assert s == 'NDData([1, 2, 3])'
    got = eval(s)
    assert np.all(got.data == arr1d.data)
    assert got.unit == arr1d.unit

    arr2d = NDData(np.array([[1, 2], [3, 4]]))
    s = repr(arr2d)
    assert s == textwrap.dedent("""
        NDData([[1, 2],
                [3, 4]])"""[1:])
    got = eval(s)
    assert np.all(got.data == arr2d.data)
    assert got.unit == arr2d.unit

    arr3d = NDData(np.array([[[1, 2], [3, 4]], [[5, 6], [7, 8]]]))
    s = repr(arr3d)
    assert s == textwrap.dedent("""
        NDData([[[1, 2],
                 [3, 4]],

                [[5, 6],
                 [7, 8]]])"""[1:])
    got = eval(s)
    assert np.all(got.data == arr3d.data)
    assert got.unit == arr3d.unit

    # let's add units!
    arr = NDData(np.array([1, 2, 3]), unit="km")
    s = repr(arr)
    assert s == "NDData([1, 2, 3], unit='km')"
    got = eval(s)
    assert np.all(got.data == arr.data)
    assert got.unit == arr.unit


# Not supported features
def test_slicing_not_supported():
    ndd = NDData(np.ones((5, 5)))
    with pytest.raises(TypeError):
        ndd[0]


def test_arithmetic_not_supported():
    ndd = NDData(np.ones((5, 5)))
    with pytest.raises(TypeError):
        ndd + ndd


def test_nddata_wcs_setter_error_cases():
    ndd = NDData(np.ones((5, 5)))

    # Setting with a non-WCS should raise an error
    with pytest.raises(TypeError):
        ndd.wcs = "I am not a WCS"

    naxis = 2
    # This should succeed since the WCS is currently None
    ndd.wcs = nd_testing._create_wcs_simple(naxis=naxis,
                                            ctype=['deg'] * naxis,
                                            crpix=[0] * naxis,
                                            crval=[10] * naxis,
                                            cdelt=[1] * naxis)
    with pytest.raises(ValueError):
        # This should fail since the WCS is not None
        ndd.wcs = nd_testing._create_wcs_simple(naxis=naxis,
                                                ctype=['deg'] * naxis,
                                                crpix=[0] * naxis,
                                                crval=[10] * naxis,
                                                cdelt=[1] * naxis)


def test_nddata_wcs_setter_with_low_level_wcs():
    ndd = NDData(np.ones((5, 5)))
    wcs = WCS()
    # If the wcs property is set with a low level WCS it should get
    # wrapped to high level.
    low_level = SlicedLowLevelWCS(wcs, 5)
    assert not isinstance(low_level, BaseHighLevelWCS)

    ndd.wcs = low_level

    assert isinstance(ndd.wcs, BaseHighLevelWCS)


def test_nddata_init_with_low_level_wcs():
    wcs = WCS()
    low_level = SlicedLowLevelWCS(wcs, 5)
    ndd = NDData(np.ones((5, 5)), wcs=low_level)
    assert isinstance(ndd.wcs, BaseHighLevelWCS)


class NDDataCustomWCS(NDData):
    @property
    def wcs(self):
        return WCS()


def test_overriden_wcs():
    # Check that a sub-class that overrides `.wcs` without providing a setter
    # works
    NDDataCustomWCS(np.ones((5, 5)))
