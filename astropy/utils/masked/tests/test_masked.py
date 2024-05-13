# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Test masked class initialization, methods, and operators.

Functions, including ufuncs, are tested in test_functions.py
"""

import operator
import sys

import numpy as np
import pytest
from numpy.testing import assert_array_equal

from astropy import units as u
from astropy.coordinates import Longitude
from astropy.units import Quantity
from astropy.utils.compat import NUMPY_LT_2_0
from astropy.utils.compat.optional_deps import HAS_PLT
from astropy.utils.masked import Masked, MaskedNDArray


def assert_masked_equal(a, b):
    assert_array_equal(a.unmasked, b.unmasked)
    assert_array_equal(a.mask, b.mask)


VARIOUS_ITEMS = [(1, 1), slice(None, 1), (), 1]


class ArraySetup:
    _data_cls = np.ndarray

    @classmethod
    def setup_class(self):
        self.a = np.arange(6.0).reshape(2, 3)
        self.mask_a = np.array([[True, False, False], [False, True, False]])
        self.b = np.array([-3.0, -2.0, -1.0])
        self.mask_b = np.array([False, True, False])
        self.c = np.array([[0.25], [0.5]])
        self.mask_c = np.array([[False], [True]])
        self.sdt = np.dtype([("a", "f8"), ("b", "f8")])
        self.mask_sdt = np.dtype([("a", "?"), ("b", "?")])
        self.sa = np.array(
            [
                [(1.0, 2.0), (3.0, 4.0)],
                [(11.0, 12.0), (13.0, 14.0)],
            ],
            dtype=self.sdt,
        )
        self.mask_sa = np.array(
            [
                [(True, True), (False, False)],
                [(False, True), (True, False)],
            ],
            dtype=self.mask_sdt,
        )
        self.sb = np.array([(1.0, 2.0), (-3.0, 4.0)], dtype=self.sdt)
        self.mask_sb = np.array([(True, False), (False, False)], dtype=self.mask_sdt)
        self.scdt = np.dtype([("sa", "2f8"), ("sb", "i8", (2, 2))])
        self.sc = np.array(
            [
                ([1.0, 2.0], [[1, 2], [3, 4]]),
                ([-1.0, -2.0], [[-1, -2], [-3, -4]]),
            ],
            dtype=self.scdt,
        )
        self.mask_scdt = np.dtype([("sa", "2?"), ("sb", "?", (2, 2))])
        self.mask_sc = np.array(
            [
                ([True, False], [[False, False], [True, True]]),
                ([False, True], [[True, False], [False, True]]),
            ],
            dtype=self.mask_scdt,
        )


class QuantitySetup(ArraySetup):
    _data_cls = Quantity

    @classmethod
    def setup_class(self):
        super().setup_class()
        self.a = Quantity(self.a, u.m)
        self.b = Quantity(self.b, u.cm)
        self.c = Quantity(self.c, u.km)
        self.sa = Quantity(self.sa, u.m, dtype=self.sdt)
        self.sb = Quantity(self.sb, u.cm, dtype=self.sdt)


class LongitudeSetup(ArraySetup):
    _data_cls = Longitude

    @classmethod
    def setup_class(self):
        super().setup_class()
        self.a = Longitude(self.a, u.deg)
        self.b = Longitude(self.b, u.deg)
        self.c = Longitude(self.c, u.deg)
        # Note: Longitude does not work on structured arrays, so
        # leaving it as regular array (which just reruns some tests).


class TestMaskedArrayInitialization(ArraySetup):
    def test_simple(self):
        ma = Masked(self.a, mask=self.mask_a)
        assert isinstance(ma, np.ndarray)
        assert isinstance(ma, type(self.a))
        assert isinstance(ma, Masked)
        assert_array_equal(ma.unmasked, self.a)
        assert_array_equal(ma.mask, self.mask_a)
        assert ma.mask is not self.mask_a
        assert np.may_share_memory(ma.mask, self.mask_a)

    def test_structured(self):
        ma = Masked(self.sa, mask=self.mask_sa)
        assert isinstance(ma, np.ndarray)
        assert isinstance(ma, type(self.sa))
        assert isinstance(ma, Masked)
        assert_array_equal(ma.unmasked, self.sa)
        assert_array_equal(ma.mask, self.mask_sa)
        assert ma.mask is not self.mask_sa
        assert np.may_share_memory(ma.mask, self.mask_sa)


def test_masked_ndarray_init():
    # Note: as a straight ndarray subclass, MaskedNDArray passes on
    # the arguments relevant for np.ndarray, not np.array.
    a_in = np.arange(3, dtype=int)
    m_in = np.array([True, False, False])
    buff = a_in.tobytes()
    # Check we're doing things correctly using regular ndarray.
    a = np.ndarray(shape=(3,), dtype=int, buffer=buff)
    assert_array_equal(a, a_in)
    # Check with and without mask.
    ma = MaskedNDArray((3,), dtype=int, mask=m_in, buffer=buff)
    assert_array_equal(ma.unmasked, a_in)
    assert_array_equal(ma.mask, m_in)
    ma = MaskedNDArray((3,), dtype=int, buffer=buff)
    assert_array_equal(ma.unmasked, a_in)
    assert_array_equal(ma.mask, np.zeros(3, bool))


def test_cannot_initialize_with_masked():
    with pytest.raises(ValueError, match="cannot handle np.ma.masked"):
        Masked(np.ma.masked)


def test_cannot_just_use_anything_with_a_mask_attribute():
    class my_array(np.ndarray):
        mask = True

    a = np.array([1.0, 2.0]).view(my_array)
    with pytest.raises(AttributeError, match="unmasked"):
        Masked(a)


class TestMaskedClassCreation:
    """Try creating a MaskedList and subclasses.

    By no means meant to be realistic, just to check that the basic
    machinery allows it.
    """

    @classmethod
    def setup_class(self):
        self._base_classes_orig = Masked._base_classes.copy()
        self._masked_classes_orig = Masked._masked_classes.copy()

        class MaskedList(Masked, list, base_cls=list, data_cls=list):
            def __new__(cls, *args, mask=None, copy=False, **kwargs):
                self = super().__new__(cls)
                self._unmasked = self._data_cls(*args, **kwargs)
                self.mask = mask
                return self

            # Need to have shape for basics to work.
            @property
            def shape(self):
                return (len(self._unmasked),)

        self.MaskedList = MaskedList

    def teardown_class(self):
        Masked._base_classes = self._base_classes_orig
        Masked._masked_classes = self._masked_classes_orig

    def test_setup(self):
        assert issubclass(self.MaskedList, Masked)
        assert issubclass(self.MaskedList, list)
        assert Masked(list) is self.MaskedList

    def test_masked_list(self):
        ml = self.MaskedList(range(3), mask=[True, False, False])
        assert ml.unmasked == [0, 1, 2]
        assert_array_equal(ml.mask, np.array([True, False, False]))
        ml01 = ml[:2]
        assert ml01.unmasked == [0, 1]
        assert_array_equal(ml01.mask, np.array([True, False]))

    def test_from_list(self):
        ml = Masked([1, 2, 3], mask=[True, False, False])
        assert ml.unmasked == [1, 2, 3]
        assert_array_equal(ml.mask, np.array([True, False, False]))

    def test_masked_list_subclass(self):
        class MyList(list):
            pass

        ml = MyList(range(3))
        mml = Masked(ml, mask=[False, True, False])
        assert isinstance(mml, Masked)
        assert isinstance(mml, MyList)
        assert isinstance(mml.unmasked, MyList)
        assert mml.unmasked == [0, 1, 2]
        assert_array_equal(mml.mask, np.array([False, True, False]))

        assert Masked(MyList) is type(mml)


class TestMaskedNDArraySubclassCreation:
    """Test that masked subclasses can be created directly and indirectly."""

    @classmethod
    def setup_class(self):
        class MyArray(np.ndarray):
            def __new__(cls, *args, **kwargs):
                return np.asanyarray(*args, **kwargs).view(cls)

        self.MyArray = MyArray
        self.a = np.array([1.0, 2.0]).view(self.MyArray)
        self.m = np.array([True, False], dtype=bool)

    def teardown_method(self, method):
        Masked._masked_classes.pop(self.MyArray, None)

    def test_direct_creation(self):
        assert self.MyArray not in Masked._masked_classes
        mcls = Masked(self.MyArray)
        assert issubclass(mcls, Masked)
        assert issubclass(mcls, self.MyArray)
        assert mcls.__name__ == "MaskedMyArray"
        assert mcls.__doc__.startswith("Masked version of MyArray")
        mms = mcls(self.a, mask=self.m)
        assert isinstance(mms, mcls)
        assert_array_equal(mms.unmasked, self.a)
        assert_array_equal(mms.mask, self.m)

    def test_initialization_without_mask(self):
        # Default for not giving a mask should be False.
        mcls = Masked(self.MyArray)
        mms = mcls(self.a)
        assert isinstance(mms, mcls)
        assert_array_equal(mms.unmasked, self.a)
        assert_array_equal(mms.mask, np.zeros(mms.shape, bool))

    @pytest.mark.parametrize("masked_array", [Masked, np.ma.MaskedArray])
    def test_initialization_with_masked_values(self, masked_array):
        mcls = Masked(self.MyArray)
        ma = masked_array(np.asarray(self.a), mask=self.m)
        mms = mcls(ma)
        assert isinstance(mms, Masked)
        assert isinstance(mms, self.MyArray)
        assert_array_equal(mms.unmasked, self.a)
        assert_array_equal(mms.mask, self.m)

    def test_indirect_creation(self):
        assert self.MyArray not in Masked._masked_classes
        mms = Masked(self.a, mask=self.m)
        assert isinstance(mms, Masked)
        assert isinstance(mms, self.MyArray)
        assert_array_equal(mms.unmasked, self.a)
        assert_array_equal(mms.mask, self.m)
        assert self.MyArray in Masked._masked_classes
        assert Masked(self.MyArray) is type(mms)

    def test_can_initialize_with_masked_values(self):
        mcls = Masked(self.MyArray)
        mms = mcls(Masked(np.asarray(self.a), mask=self.m))
        assert isinstance(mms, Masked)
        assert isinstance(mms, self.MyArray)
        assert_array_equal(mms.unmasked, self.a)
        assert_array_equal(mms.mask, self.m)

    def test_viewing(self):
        mms = Masked(self.a, mask=self.m)
        mms2 = mms.view()
        assert type(mms2) is mms.__class__
        assert_masked_equal(mms2, mms)

        ma = mms.view(np.ndarray)
        assert type(ma) is MaskedNDArray
        assert_array_equal(ma.unmasked, self.a.view(np.ndarray))
        assert_array_equal(ma.mask, self.m)

    def test_viewing_independent_shape(self):
        mms = Masked(self.a, mask=self.m)
        mms2 = mms.view()
        mms2.shape = mms2.shape[::-1]
        assert mms2.shape == mms.shape[::-1]
        assert mms2.mask.shape == mms.shape[::-1]
        # This should not affect the original array!
        assert mms.shape == self.a.shape
        assert mms.mask.shape == self.a.shape


class TestMaskedQuantityInitialization(TestMaskedArrayInitialization, QuantitySetup):
    @classmethod
    def setup_class(self):
        super().setup_class()
        # Ensure we have used MaskedQuantity before - just in case a single test gets
        # called; see gh-15316.
        self.MQ = Masked(Quantity)

    def test_masked_quantity_getting(self):
        # First check setup_class (or previous use) defined a cache entry.
        mcls = Masked._masked_classes[type(self.a)]
        # Next check this is what one gets now.
        MQ = Masked(Quantity)
        assert MQ is mcls
        assert issubclass(MQ, Quantity)
        assert issubclass(MQ, Masked)
        # And also what one gets implicitly.
        mq = Masked([1.0, 2.0] * u.m)
        assert isinstance(mq, MQ)

    def test_masked_quantity_class_init(self):
        # This is not a very careful test.
        mq = self.MQ([1.0, 2.0], mask=[True, False], unit=u.s)
        assert mq.unit == u.s
        assert np.all(mq.value.unmasked == [1.0, 2.0])
        assert np.all(mq.value.mask == [True, False])
        assert np.all(mq.mask == [True, False])

    def test_initialization_without_mask(self):
        # Default for not giving a mask should be False.
        mq = self.MQ([1.0, 2.0], u.s)
        assert mq.unit == u.s
        assert np.all(mq.value.unmasked == [1.0, 2.0])
        assert np.all(mq.mask == [False, False])

    @pytest.mark.parametrize("masked_array", [Masked, np.ma.MaskedArray])
    def test_initialization_with_masked_values(self, masked_array):
        a = np.array([1.0, 2.0])
        m = np.array([True, False])
        ma = masked_array(a, m)
        mq = self.MQ(ma)
        assert isinstance(mq, self.MQ)
        assert_array_equal(mq.value.unmasked, a)
        assert_array_equal(mq.mask, m)


class TestMaskSetting(ArraySetup):
    def test_whole_mask_setting_simple(self):
        ma = Masked(self.a)
        assert ma.mask.shape == ma.shape
        assert not ma.mask.any()
        ma.mask = True
        assert ma.mask.shape == ma.shape
        assert ma.mask.all()
        ma.mask = [[True], [False]]
        assert ma.mask.shape == ma.shape
        assert_array_equal(ma.mask, np.array([[True] * 3, [False] * 3]))
        ma.mask = self.mask_a
        assert ma.mask.shape == ma.shape
        assert_array_equal(ma.mask, self.mask_a)
        assert ma.mask is not self.mask_a
        assert np.may_share_memory(ma.mask, self.mask_a)

    def test_whole_mask_setting_structured(self):
        ma = Masked(self.sa)
        assert ma.mask.shape == ma.shape
        assert not ma.mask["a"].any() and not ma.mask["b"].any()
        ma.mask = True
        assert ma.mask.shape == ma.shape
        assert ma.mask["a"].all() and ma.mask["b"].all()
        ma.mask = [[True], [False]]
        assert ma.mask.shape == ma.shape
        assert_array_equal(
            ma.mask,
            np.array([[(True, True)] * 2, [(False, False)] * 2], dtype=self.mask_sdt),
        )
        ma.mask = self.mask_sa
        assert ma.mask.shape == ma.shape
        assert_array_equal(ma.mask, self.mask_sa)
        assert ma.mask is not self.mask_sa
        assert np.may_share_memory(ma.mask, self.mask_sa)

    @pytest.mark.parametrize("item", VARIOUS_ITEMS)
    def test_part_mask_setting(self, item):
        ma = Masked(self.a)
        ma.mask[item] = True
        expected = np.zeros(ma.shape, bool)
        expected[item] = True
        assert_array_equal(ma.mask, expected)
        ma.mask[item] = False
        assert_array_equal(ma.mask, np.zeros(ma.shape, bool))
        # Mask propagation
        mask = np.zeros(self.a.shape, bool)
        ma = Masked(self.a, mask)
        ma.mask[item] = True
        assert np.may_share_memory(ma.mask, mask)
        assert_array_equal(ma.mask, mask)

    @pytest.mark.parametrize("item", ["a"] + VARIOUS_ITEMS)
    def test_part_mask_setting_structured(self, item):
        ma = Masked(self.sa)
        ma.mask[item] = True
        expected = np.zeros(ma.shape, self.mask_sdt)
        expected[item] = True
        assert_array_equal(ma.mask, expected)
        ma.mask[item] = False
        assert_array_equal(ma.mask, np.zeros(ma.shape, self.mask_sdt))
        # Mask propagation
        mask = np.zeros(self.sa.shape, self.mask_sdt)
        ma = Masked(self.sa, mask)
        ma.mask[item] = True
        assert np.may_share_memory(ma.mask, mask)
        assert_array_equal(ma.mask, mask)


# Following are tests where we trust the initializer works.


class MaskedArraySetup(ArraySetup):
    @classmethod
    def setup_class(self):
        super().setup_class()
        self.ma = Masked(self.a, mask=self.mask_a)
        self.mb = Masked(self.b, mask=self.mask_b)
        self.mc = Masked(self.c, mask=self.mask_c)
        self.msa = Masked(self.sa, mask=self.mask_sa)
        self.msb = Masked(self.sb, mask=self.mask_sb)
        self.msc = Masked(self.sc, mask=self.mask_sc)


class TestViewing(MaskedArraySetup):
    def test_viewing_as_new_type(self):
        ma2 = self.ma.view(type(self.ma))
        assert_masked_equal(ma2, self.ma)

        ma3 = self.ma.view()
        assert_masked_equal(ma3, self.ma)

    def test_viewing_as_new_dtype(self):
        # Not very meaningful, but possible...
        ma2 = self.ma.view("c8")
        assert_array_equal(ma2.unmasked, self.a.view("c8"))
        assert_array_equal(ma2.mask, self.mask_a)

    def test_viewing_as_new_structured_dtype(self):
        ma2 = self.ma.view("f8,f8,f8")
        assert_array_equal(ma2.unmasked, self.a.view("f8,f8,f8"))
        assert_array_equal(ma2.mask, self.mask_a.view("?,?,?"))
        # Check round-trip
        ma3 = ma2.view(self.ma.dtype)
        assert_array_equal(ma3.unmasked, self.ma.unmasked)
        assert_array_equal(ma3.mask, self.mask_a)

    @pytest.mark.parametrize("new_dtype", ["f4", "2f4"])
    def test_viewing_as_new_dtype_not_implemented(self, new_dtype):
        # But cannot (yet) view in way that would need to create a new mask,
        # even though that view is possible for a regular array.
        check = self.a.view(new_dtype)
        with pytest.raises(NotImplementedError, match="different.*size"):
            self.ma.view(new_dtype)

    def test_viewing_as_something_impossible(self):
        with pytest.raises(TypeError):
            # Use intp to ensure have the same size as object,
            # otherwise we get a different error message
            Masked(np.array([1, 2], dtype=np.intp)).view(Masked)


class TestMaskedArrayCopyFilled(MaskedArraySetup):
    def test_copy(self):
        ma_copy = self.ma.copy()
        assert type(ma_copy) is type(self.ma)
        assert_array_equal(ma_copy.unmasked, self.ma.unmasked)
        assert_array_equal(ma_copy.mask, self.ma.mask)
        assert not np.may_share_memory(ma_copy.unmasked, self.ma.unmasked)
        assert not np.may_share_memory(ma_copy.mask, self.ma.mask)

    @pytest.mark.parametrize("fill_value", (0, 1))
    def test_filled(self, fill_value):
        fill_value = fill_value * getattr(self.a, "unit", 1)
        expected = self.a.copy()
        expected[self.ma.mask] = fill_value
        result = self.ma.filled(fill_value)
        assert_array_equal(expected, result)

    def test_filled_no_fill_value(self):
        with pytest.raises(TypeError, match="missing 1 required"):
            self.ma.filled()

    @pytest.mark.parametrize("fill_value", [(0, 1), (-1, -1)])
    def test_filled_structured(self, fill_value):
        fill_value = np.array(fill_value, dtype=self.sdt)
        if hasattr(self.sa, "unit"):
            fill_value = fill_value << self.sa.unit
        expected = self.sa.copy()
        expected["a"][self.msa.mask["a"]] = fill_value["a"]
        expected["b"][self.msa.mask["b"]] = fill_value["b"]
        result = self.msa.filled(fill_value)
        assert_array_equal(expected, result)

    def test_flat(self):
        ma_copy = self.ma.copy()
        ma_flat = ma_copy.flat
        # Check that single item keeps class and mask
        ma_flat1 = ma_flat[1]
        assert ma_flat1.unmasked == self.a.flat[1]
        assert ma_flat1.mask == self.mask_a.flat[1]
        # As well as getting items via iteration.
        assert all(
            (ma.unmasked == a and ma.mask == m)
            for (ma, a, m) in zip(self.ma.flat, self.a.flat, self.mask_a.flat)
        )

        # check that flat works like a view of the real array
        ma_flat[1] = self.b[1]
        assert ma_flat[1] == self.b[1]
        assert ma_copy[0, 1] == self.b[1]


class TestMaskedQuantityCopyFilled(TestMaskedArrayCopyFilled, QuantitySetup):
    pass


class TestMaskedLongitudeCopyFilled(TestMaskedArrayCopyFilled, LongitudeSetup):
    pass


class TestMaskedArrayShaping(MaskedArraySetup):
    def test_reshape(self):
        ma_reshape = self.ma.reshape((6,))
        expected_data = self.a.reshape((6,))
        expected_mask = self.mask_a.reshape((6,))
        assert ma_reshape.shape == expected_data.shape
        assert_array_equal(ma_reshape.unmasked, expected_data)
        assert_array_equal(ma_reshape.mask, expected_mask)

    def test_shape_setting(self):
        ma_reshape = self.ma.copy()
        ma_reshape.shape = (6,)
        expected_data = self.a.reshape((6,))
        expected_mask = self.mask_a.reshape((6,))
        assert ma_reshape.shape == expected_data.shape
        assert_array_equal(ma_reshape.unmasked, expected_data)
        assert_array_equal(ma_reshape.mask, expected_mask)

    def test_shape_setting_failure(self):
        ma = self.ma.copy()
        with pytest.raises(ValueError, match="cannot reshape"):
            ma.shape = (5,)

        assert ma.shape == self.ma.shape
        assert ma.mask.shape == self.ma.shape

        # Here, mask can be reshaped but array cannot.
        ma2 = Masked(np.broadcast_to([[1.0], [2.0]], self.a.shape), mask=self.mask_a)
        with pytest.raises(AttributeError, match="ncompatible shape"):
            ma2.shape = (6,)

        assert ma2.shape == self.ma.shape
        assert ma2.mask.shape == self.ma.shape

        # Here, array can be reshaped but mask cannot.
        ma3 = Masked(
            self.a.copy(), mask=np.broadcast_to([[True], [False]], self.mask_a.shape)
        )
        with pytest.raises(AttributeError, match="ncompatible shape"):
            ma3.shape = (6,)

        assert ma3.shape == self.ma.shape
        assert ma3.mask.shape == self.ma.shape

    def test_ravel(self):
        ma_ravel = self.ma.ravel()
        expected_data = self.a.ravel()
        expected_mask = self.mask_a.ravel()
        assert ma_ravel.shape == expected_data.shape
        assert_array_equal(ma_ravel.unmasked, expected_data)
        assert_array_equal(ma_ravel.mask, expected_mask)

    def test_transpose(self):
        ma_transpose = self.ma.transpose()
        expected_data = self.a.transpose()
        expected_mask = self.mask_a.transpose()
        assert ma_transpose.shape == expected_data.shape
        assert_array_equal(ma_transpose.unmasked, expected_data)
        assert_array_equal(ma_transpose.mask, expected_mask)

    def test_iter(self):
        for ma, d, m in zip(self.ma, self.a, self.mask_a):
            assert_array_equal(ma.unmasked, d)
            assert_array_equal(ma.mask, m)


class MaskedItemTests(MaskedArraySetup):
    @pytest.mark.parametrize("item", VARIOUS_ITEMS)
    def test_getitem(self, item):
        ma_part = self.ma[item]
        expected_data = self.a[item]
        expected_mask = self.mask_a[item]
        assert_array_equal(ma_part.unmasked, expected_data)
        assert_array_equal(ma_part.mask, expected_mask)

    @pytest.mark.parametrize("item", ["a"] + VARIOUS_ITEMS)
    def test_getitem_structured(self, item):
        ma_part = self.msa[item]
        expected_data = self.sa[item]
        expected_mask = self.mask_sa[item]
        assert_array_equal(ma_part.unmasked, expected_data)
        assert_array_equal(ma_part.mask, expected_mask)

    @pytest.mark.parametrize(
        "indices,axis",
        [([0, 1], 1), ([0, 1], 0), ([0, 1], None), ([[0, 1], [2, 3]], None)],
    )
    def test_take(self, indices, axis):
        ma_take = self.ma.take(indices, axis=axis)
        expected_data = self.a.take(indices, axis=axis)
        expected_mask = self.mask_a.take(indices, axis=axis)
        assert_array_equal(ma_take.unmasked, expected_data)
        assert_array_equal(ma_take.mask, expected_mask)
        ma_take2 = np.take(self.ma, indices, axis=axis)
        assert_masked_equal(ma_take2, ma_take)

    @pytest.mark.parametrize("item", VARIOUS_ITEMS)
    @pytest.mark.parametrize("mask", [None, True, False])
    def test_setitem(self, item, mask):
        base = self.ma.copy()
        expected_data = self.a.copy()
        expected_mask = self.mask_a.copy()
        value = self.a[0, 0] if mask is None else Masked(self.a[0, 0], mask)
        base[item] = value
        expected_data[item] = value if mask is None else value.unmasked
        expected_mask[item] = False if mask is None else value.mask
        assert_array_equal(base.unmasked, expected_data)
        assert_array_equal(base.mask, expected_mask)

    @pytest.mark.parametrize("item", ["a"] + VARIOUS_ITEMS)
    @pytest.mark.parametrize("mask", [None, True, False])
    def test_setitem_structured(self, item, mask):
        base = self.msa.copy()
        expected_data = self.sa.copy()
        expected_mask = self.mask_sa.copy()
        value = self.sa["b"] if item == "a" else self.sa[0, 0]
        if mask is not None:
            value = Masked(value, mask)
        base[item] = value
        expected_data[item] = value if mask is None else value.unmasked
        expected_mask[item] = False if mask is None else value.mask
        assert_array_equal(base.unmasked, expected_data)
        assert_array_equal(base.mask, expected_mask)

    @pytest.mark.parametrize("item", VARIOUS_ITEMS)
    def test_setitem_np_ma_masked(self, item):
        base = self.ma.copy()
        expected_mask = self.mask_a.copy()
        base[item] = np.ma.masked
        expected_mask[item] = True
        assert_array_equal(base.unmasked, self.a)
        assert_array_equal(base.mask, expected_mask)

    @pytest.mark.parametrize("item", VARIOUS_ITEMS)
    def test_hash(self, item):
        ma_part = self.ma[item]
        if ma_part.ndim > 0 or ma_part.mask.any():
            with pytest.raises(TypeError, match="unhashable"):
                hash(ma_part)
        else:
            assert hash(ma_part) == hash(ma_part.unmasked)


class TestMaskedArrayItems(MaskedItemTests):
    @classmethod
    def setup_class(self):
        super().setup_class()
        self.d = np.array(["aa", "bb"])
        self.mask_d = np.array([True, False])
        self.md = Masked(self.d, self.mask_d)

    # Quantity, Longitude cannot hold strings.
    def test_getitem_strings(self):
        md = self.md.copy()
        md0 = md[0]
        assert md0.unmasked == self.d[0]
        assert md0.mask
        md_all = md[:]
        assert_masked_equal(md_all, md)

    def test_setitem_strings_np_ma_masked(self):
        md = self.md.copy()
        md[1] = np.ma.masked
        assert_array_equal(md.unmasked, self.d)
        assert_array_equal(md.mask, np.ones(2, bool))


class TestMaskedQuantityItems(MaskedItemTests, QuantitySetup):
    pass


class TestMaskedLongitudeItems(MaskedItemTests, LongitudeSetup):
    pass


class MaskedOperatorTests(MaskedArraySetup):
    @pytest.mark.parametrize("op", (operator.add, operator.sub))
    def test_add_subtract(self, op):
        mapmb = op(self.ma, self.mb)
        expected_data = op(self.a, self.b)
        expected_mask = self.ma.mask | self.mb.mask
        # Note: assert_array_equal also checks type, i.e., that, e.g.,
        # Longitude decays into an Angle.
        assert_array_equal(mapmb.unmasked, expected_data)
        assert_array_equal(mapmb.mask, expected_mask)

    @pytest.mark.parametrize("op", (operator.eq, operator.ne))
    def test_equality(self, op):
        mapmb = op(self.ma, self.mb)
        expected_data = op(self.a, self.b)
        expected_mask = self.ma.mask | self.mb.mask
        # Note: assert_array_equal also checks type, i.e., that boolean
        # output is represented as plain Masked ndarray.
        assert_array_equal(mapmb.unmasked, expected_data)
        assert_array_equal(mapmb.mask, expected_mask)

    def test_not_implemented(self):
        with pytest.raises(TypeError):
            self.ma > "abc"  # noqa: B015

    @pytest.mark.parametrize("different_names", [False, True])
    @pytest.mark.parametrize("op", (operator.eq, operator.ne))
    def test_structured_equality(self, op, different_names):
        msb = self.msb
        if different_names:
            msb = msb.astype(
                [(f"different_{name}", dt) for name, dt in msb.dtype.fields.items()]
            )
        mapmb = op(self.msa, self.msb)
        # Expected is a bit tricky here: only unmasked fields count
        expected_data = np.ones(mapmb.shape, bool)
        expected_mask = np.ones(mapmb.shape, bool)
        for field in self.sdt.names:
            fa, mfa = self.sa[field], self.mask_sa[field]
            fb, mfb = self.sb[field], self.mask_sb[field]
            mfequal = mfa | mfb
            fequal = (fa == fb) | mfequal
            expected_data &= fequal
            expected_mask &= mfequal

        if op is operator.ne:
            expected_data = ~expected_data

        # Note: assert_array_equal also checks type, i.e., that boolean
        # output is represented as plain Masked ndarray.
        assert_array_equal(mapmb.unmasked, expected_data)
        assert_array_equal(mapmb.mask, expected_mask)

    def test_matmul(self):
        result = self.ma.T @ self.ma
        assert_array_equal(result.unmasked, self.a.T @ self.a)
        mask1 = np.any(self.mask_a, axis=0)
        expected_mask = np.logical_or.outer(mask1, mask1)
        assert_array_equal(result.mask, expected_mask)
        result2 = self.ma.T @ self.a
        assert_array_equal(result2.unmasked, self.a.T @ self.a)
        expected_mask2 = np.logical_or.outer(mask1, np.zeros(3, bool))
        assert_array_equal(result2.mask, expected_mask2)
        result3 = self.a.T @ self.ma
        assert_array_equal(result3.unmasked, self.a.T @ self.a)
        expected_mask3 = np.logical_or.outer(np.zeros(3, bool), mask1)
        assert_array_equal(result3.mask, expected_mask3)

    def test_matmul_axes(self):
        m1 = Masked(np.arange(27.0).reshape(3, 3, 3))
        m2 = Masked(np.arange(-27.0, 0.0).reshape(3, 3, 3))
        mxm1 = np.matmul(m1, m2)
        exp = np.matmul(m1.unmasked, m2.unmasked)
        assert_array_equal(mxm1.unmasked, exp)
        assert_array_equal(mxm1.mask, False)
        m1.mask[0, 1, 2] = True
        m2.mask[0, 2, 0] = True
        axes = [(0, 2), (-2, -1), (0, 1)]
        mxm2 = np.matmul(m1, m2, axes=axes)
        exp2 = np.matmul(m1.unmasked, m2.unmasked, axes=axes)
        # Any unmasked result will have all elements contributing unity,
        # while masked entries mean the total will be lower.
        mask2 = (
            np.matmul(
                (~m1.mask).astype(int),
                (~m2.mask).astype(int),
                axes=axes,
            )
            != m1.shape[axes[0][1]]
        )
        assert_array_equal(mxm2.unmasked, exp2)
        assert_array_equal(mxm2.mask, mask2)

    def test_matvec(self):
        result = self.ma @ self.mb
        assert np.all(result.mask)
        assert_array_equal(result.unmasked, self.a @ self.b)
        # Just using the masked vector still has all elements masked.
        result2 = self.a @ self.mb
        assert np.all(result2.mask)
        assert_array_equal(result2.unmasked, self.a @ self.b)
        new_ma = self.ma.copy()
        new_ma.mask[0, 0] = False
        result3 = new_ma @ self.b
        assert_array_equal(result3.unmasked, self.a @ self.b)
        assert_array_equal(result3.mask, new_ma.mask.any(-1))

    def test_vecmat(self):
        result = self.mb @ self.ma.T
        assert np.all(result.mask)
        assert_array_equal(result.unmasked, self.b @ self.a.T)
        result2 = self.b @ self.ma.T
        assert np.all(result2.mask)
        assert_array_equal(result2.unmasked, self.b @ self.a.T)
        new_ma = self.ma.T.copy()
        new_ma.mask[0, 0] = False
        result3 = self.b @ new_ma
        assert_array_equal(result3.unmasked, self.b @ self.a.T)
        assert_array_equal(result3.mask, new_ma.mask.any(0))

    def test_vecvec(self):
        result = self.mb @ self.mb
        assert result.shape == ()
        assert result.mask
        assert result.unmasked == self.b @ self.b
        mb_no_mask = Masked(self.b, False)
        result2 = mb_no_mask @ mb_no_mask
        assert not result2.mask


class TestMaskedArrayOperators(MaskedOperatorTests):
    # Some further tests that use strings, which are not useful for Quantity.
    @pytest.mark.parametrize("op", (operator.eq, operator.ne))
    def test_equality_strings(self, op):
        m1 = Masked(np.array(["a", "b", "c"]), mask=[True, False, False])
        m2 = Masked(np.array(["a", "b", "d"]), mask=[False, False, False])
        result = op(m1, m2)
        assert_array_equal(result.unmasked, op(m1.unmasked, m2.unmasked))
        assert_array_equal(result.mask, m1.mask | m2.mask)

        result2 = op(m1, m2.unmasked)
        assert_masked_equal(result2, result)

    def test_not_implemented(self):
        with pytest.raises(TypeError):
            Masked(["a", "b"]) > object()  # noqa: B015


class TestMaskedQuantityOperators(MaskedOperatorTests, QuantitySetup):
    pass


class TestMaskedLongitudeOperators(MaskedOperatorTests, LongitudeSetup):
    pass


class TestMaskedArrayMethods(MaskedArraySetup):
    def test_round(self):
        # Goes via ufunc, hence easy.
        mrc = self.mc.round()
        expected = Masked(self.c.round(), self.mask_c)
        assert_masked_equal(mrc, expected)

    @pytest.mark.parametrize("axis", (0, 1, None))
    def test_sum(self, axis):
        ma_sum = self.ma.sum(axis)
        expected_data = self.a.sum(axis)
        expected_mask = self.ma.mask.any(axis)
        assert_array_equal(ma_sum.unmasked, expected_data)
        assert_array_equal(ma_sum.mask, expected_mask)

    @pytest.mark.parametrize("axis", (0, 1, None))
    def test_sum_where(self, axis):
        where = np.array(
            [
                [True, False, False],
                [True, True, True],
            ]
        )
        where_final = ~self.ma.mask & where
        ma_sum = self.ma.sum(axis, where=where_final)
        expected_data = self.ma.unmasked.sum(axis, where=where_final)
        expected_mask = np.logical_or.reduce(
            self.ma.mask, axis=axis, where=where_final
        ) | (~where_final).all(axis)
        assert_array_equal(ma_sum.unmasked, expected_data)
        assert_array_equal(ma_sum.mask, expected_mask)

    def test_sum_hash(self):
        ma_sum = self.ma.sum()
        assert ma_sum.mask
        # Masked scalars cannot be hashed.
        with pytest.raises(TypeError, match="unhashable"):
            hash(ma_sum)
        ma_sum2 = Masked(self.a).sum()
        # But an unmasked scalar can.
        assert hash(ma_sum2) == hash(self.a.sum())

    @pytest.mark.parametrize("axis", (0, 1, None))
    def test_cumsum(self, axis):
        ma_sum = self.ma.cumsum(axis)
        expected_data = self.a.cumsum(axis)
        mask = self.mask_a
        if axis is None:
            mask = mask.ravel()
        expected_mask = np.logical_or.accumulate(mask, axis=axis)
        assert_array_equal(ma_sum.unmasked, expected_data)
        assert_array_equal(ma_sum.mask, expected_mask)

    @pytest.mark.parametrize("axis", (0, 1, None))
    def test_mean(self, axis):
        ma_mean = self.ma.mean(axis)
        filled = self.a.copy()
        filled[self.mask_a] = 0.0
        count = 1 - self.ma.mask.astype(int)
        expected_data = filled.sum(axis) / count.sum(axis)
        expected_mask = self.ma.mask.all(axis)
        assert_array_equal(ma_mean.unmasked, expected_data)
        assert_array_equal(ma_mean.mask, expected_mask)

    @pytest.mark.parametrize("axis", (0, 1, None))
    def test_mean_all_masked(self, axis):
        # test corner case when all values are masked
        md = Masked(self.a, np.ones(self.a.shape, dtype=bool))
        md_mean = md.mean(axis)
        assert np.all(np.isnan(md_mean.unmasked))
        assert np.all(md_mean.mask)

    def test_mean_int16(self):
        ma = self.ma.astype("i2")
        ma_mean = ma.mean()
        assert ma_mean.dtype == "f8"
        expected = ma.astype("f8").mean()
        assert_masked_equal(ma_mean, expected)

    def test_mean_float16(self):
        ma = self.ma.astype("f2")
        ma_mean = ma.mean()
        assert ma_mean.dtype == "f2"
        expected = self.ma.mean().astype("f2")
        assert_masked_equal(ma_mean, expected)

    def test_mean_inplace(self):
        expected = self.ma.mean(1)
        out = Masked(np.zeros_like(expected.unmasked))
        result = self.ma.mean(1, out=out)
        assert result is out
        assert_masked_equal(out, expected)

    @pytest.mark.filterwarnings("ignore:.*encountered in.*divide")
    @pytest.mark.filterwarnings("ignore:Mean of empty slice")
    @pytest.mark.parametrize("axis", (0, 1, None))
    def test_mean_where(self, axis):
        where = np.array(
            [
                [True, False, False],
                [True, True, True],
            ]
        )
        where_final = ~self.ma.mask & where
        ma_mean = self.ma.mean(axis, where=where)
        expected_data = self.ma.unmasked.mean(axis, where=where_final)
        expected_mask = np.logical_or.reduce(
            self.ma.mask, axis=axis, where=where_final
        ) | (~where_final).all(axis)
        assert_array_equal(ma_mean.unmasked, expected_data)
        assert_array_equal(ma_mean.mask, expected_mask)

    def test_mean_hash(self):
        ma_mean = self.ma.mean()
        assert hash(ma_mean) == hash(ma_mean.unmasked[()])

    @pytest.mark.filterwarnings("ignore:.*encountered in.*divide")
    @pytest.mark.parametrize("axis", (0, 1, None))
    def test_var(self, axis):
        ma_var = self.ma.var(axis)
        filled = (self.a - self.ma.mean(axis, keepdims=True)) ** 2
        filled[self.mask_a] = 0.0
        count = (1 - self.ma.mask.astype(int)).sum(axis)
        expected_data = filled.sum(axis) / count
        expected_mask = self.ma.mask.all(axis)
        assert_array_equal(ma_var.unmasked, expected_data)
        assert_array_equal(ma_var.mask, expected_mask)
        ma_var1 = self.ma.var(axis, ddof=1)
        expected_data1 = filled.sum(axis) / (count - 1)
        expected_mask1 = self.ma.mask.all(axis) | (count <= 1)
        assert_array_equal(ma_var1.unmasked, expected_data1)
        assert_array_equal(ma_var1.mask, expected_mask1)
        ma_var5 = self.ma.var(axis, ddof=5)
        assert np.all(~np.isfinite(ma_var5.unmasked))
        assert ma_var5.mask.all()

    def test_var_int16(self):
        ma = self.ma.astype("i2")
        ma_var = ma.var()
        assert ma_var.dtype == "f8"
        expected = ma.astype("f8").var()
        assert_masked_equal(ma_var, expected)

    @pytest.mark.filterwarnings("ignore:.*encountered in.*divide")
    @pytest.mark.filterwarnings("ignore:Degrees of freedom <= 0 for slice")
    @pytest.mark.parametrize("axis", (0, 1, None))
    def test_var_where(self, axis):
        where = np.array(
            [
                [True, False, False],
                [True, True, True],
            ]
        )
        where_final = ~self.ma.mask & where
        ma_var = self.ma.var(axis, where=where)
        expected_data = self.ma.unmasked.var(axis, where=where_final)
        expected_mask = np.logical_or.reduce(
            self.ma.mask, axis=axis, where=where_final
        ) | (~where_final).all(axis)
        assert_array_equal(ma_var.unmasked, expected_data)
        assert_array_equal(ma_var.mask, expected_mask)

    def test_std(self):
        ma_std = self.ma.std(1, ddof=1)
        ma_var1 = self.ma.var(1, ddof=1)
        expected = np.sqrt(ma_var1)
        assert_masked_equal(ma_std, expected)

    def test_std_inplace(self):
        expected = self.ma.std(1, ddof=1)
        out = Masked(np.zeros_like(expected.unmasked))
        result = self.ma.std(1, ddof=1, out=out)
        assert result is out
        assert_masked_equal(result, expected)

    @pytest.mark.filterwarnings("ignore:.*encountered in.*divide")
    @pytest.mark.filterwarnings("ignore:Degrees of freedom <= 0 for slice")
    @pytest.mark.parametrize("axis", (0, 1, None))
    def test_std_where(self, axis):
        where = np.array(
            [
                [True, False, False],
                [True, True, True],
            ]
        )
        where_final = ~self.ma.mask & where
        ma_std = self.ma.std(axis, where=where)
        expected_data = self.ma.unmasked.std(axis, where=where_final)
        expected_mask = np.logical_or.reduce(
            self.ma.mask, axis=axis, where=where_final
        ) | (~where_final).all(axis)
        assert_array_equal(ma_std.unmasked, expected_data)
        assert_array_equal(ma_std.mask, expected_mask)

    @pytest.mark.parametrize("axis", (0, 1, None))
    def test_min(self, axis):
        ma_min = self.ma.min(axis)
        filled = self.a.copy()
        filled[self.mask_a] = self.a.max()
        expected_data = filled.min(axis)
        assert_array_equal(ma_min.unmasked, expected_data)
        assert not np.any(ma_min.mask)

    def test_min_with_masked_nan(self):
        ma = Masked([3.0, np.nan, 2.0], mask=[False, True, False])
        ma_min = ma.min()
        assert_array_equal(ma_min.unmasked, np.array(2.0))
        assert not ma_min.mask

    @pytest.mark.parametrize("axis", (0, 1, None))
    def test_min_where(self, axis):
        where = np.array(
            [
                [True, False, False],
                [True, True, True],
            ]
        )
        where_final = ~self.ma.mask & where
        ma_min = self.ma.min(axis, where=where_final, initial=np.inf)
        expected_data = self.ma.unmasked.min(axis, where=where_final, initial=np.inf)
        expected_mask = np.logical_or.reduce(
            self.ma.mask, axis=axis, where=where_final
        ) | (~where_final).all(axis)
        assert_array_equal(ma_min.unmasked, expected_data)
        assert_array_equal(ma_min.mask, expected_mask)

    @pytest.mark.parametrize("axis", (0, 1, None))
    def test_max(self, axis):
        ma_max = self.ma.max(axis)
        filled = self.a.copy()
        filled[self.mask_a] = self.a.min()
        expected_data = filled.max(axis)
        assert_array_equal(ma_max.unmasked, expected_data)
        assert not np.any(ma_max.mask)

    @pytest.mark.parametrize("axis", (0, 1, None))
    def test_max_where(self, axis):
        where = np.array(
            [
                [True, False, False],
                [True, True, True],
            ]
        )
        where_final = ~self.ma.mask & where
        ma_max = self.ma.max(axis, where=where_final, initial=-np.inf)
        expected_data = self.ma.unmasked.max(axis, where=where_final, initial=-np.inf)
        expected_mask = np.logical_or.reduce(
            self.ma.mask, axis=axis, where=where_final
        ) | (~where_final).all(axis)
        assert_array_equal(ma_max.unmasked, expected_data)
        assert_array_equal(ma_max.mask, expected_mask)

    @pytest.mark.parametrize("axis", (0, 1, None))
    def test_argmin(self, axis):
        ma_argmin = self.ma.argmin(axis)
        filled = self.a.copy()
        filled[self.mask_a] = self.a.max()
        expected_data = filled.argmin(axis)
        assert_array_equal(ma_argmin, expected_data)

    def test_argmin_only_one_unmasked_element(self):
        # Regression test for example from @taldcroft at
        # https://github.com/astropy/astropy/pull/11127#discussion_r600864559
        ma = Masked(data=[1, 2], mask=[True, False])
        assert ma.argmin() == 1

    def test_argmin_keepdims(self):
        ma = Masked(data=[[1, 2], [3, 4]], mask=[[True, False], [False, False]])
        assert_array_equal(ma.argmin(axis=0, keepdims=True), np.array([[1, 0]]))

    @pytest.mark.parametrize("axis", (0, 1, None))
    def test_argmax(self, axis):
        ma_argmax = self.ma.argmax(axis)
        filled = self.a.copy()
        filled[self.mask_a] = self.a.min()
        expected_data = filled.argmax(axis)
        assert_array_equal(ma_argmax, expected_data)

    def test_argmax_keepdims(self):
        ma = Masked(data=[[1, 2], [3, 4]], mask=[[True, False], [False, False]])
        assert_array_equal(ma.argmax(axis=1, keepdims=True), np.array([[1], [1]]))

    @pytest.mark.parametrize("axis", (0, 1, None))
    def test_argsort(self, axis):
        ma_argsort = self.ma.argsort(axis)
        filled = self.a.copy()
        filled[self.mask_a] = self.a.max() * 1.1
        expected_data = filled.argsort(axis)
        assert_array_equal(ma_argsort, expected_data)

    @pytest.mark.parametrize("order", [None, "a", ("a", "b"), ("b", "a")])
    @pytest.mark.parametrize("axis", [0, 1])
    def test_structured_argsort(self, axis, order):
        ma_argsort = self.msa.argsort(axis, order=order)
        filled = self.msa.filled(fill_value=np.array((np.inf, np.inf), dtype=self.sdt))
        expected_data = filled.argsort(axis, order=order)
        assert_array_equal(ma_argsort, expected_data)

    def test_argsort_error(self):
        with pytest.raises(ValueError, match="when the array has no fields"):
            self.ma.argsort(axis=0, order="a")

    @pytest.mark.parametrize("axis", (0, 1))
    def test_sort(self, axis):
        ma_sort = self.ma.copy()
        ma_sort.sort(axis)
        indices = self.ma.argsort(axis)
        expected_data = np.take_along_axis(self.ma.unmasked, indices, axis)
        expected_mask = np.take_along_axis(self.ma.mask, indices, axis)
        assert_array_equal(ma_sort.unmasked, expected_data)
        assert_array_equal(ma_sort.mask, expected_mask)

    @pytest.mark.parametrize("kth", [1, 3])
    def test_argpartition(self, kth):
        ma = self.ma.ravel()
        ma_argpartition = ma.argpartition(kth)
        partitioned = ma[ma_argpartition]
        assert (partitioned[:kth] < partitioned[kth]).all()
        assert (partitioned[kth:] >= partitioned[kth]).all()
        if partitioned[kth].mask:
            assert all(partitioned.mask[kth:])
        else:
            assert not any(partitioned.mask[:kth])

    @pytest.mark.parametrize("kth", [1, 3])
    def test_partition(self, kth):
        partitioned = self.ma.flatten()
        partitioned.partition(kth)
        assert (partitioned[:kth] < partitioned[kth]).all()
        assert (partitioned[kth:] >= partitioned[kth]).all()
        if partitioned[kth].mask:
            assert all(partitioned.mask[kth:])
        else:
            assert not any(partitioned.mask[:kth])

    def test_all_explicit(self):
        a1 = np.array(
            [
                [1.0, 2.0],
                [3.0, 4.0],
            ]
        )
        a2 = np.array(
            [
                [1.0, 0.0],
                [3.0, 4.0],
            ]
        )
        if self._data_cls is not np.ndarray:
            a1 = self._data_cls(a1, self.a.unit)
            a2 = self._data_cls(a2, self.a.unit)
        ma1 = Masked(
            a1,
            mask=[
                [False, False],
                [True, True],
            ],
        )
        ma2 = Masked(
            a2,
            mask=[
                [False, True],
                [False, True],
            ],
        )
        ma1_eq_ma2 = ma1 == ma2
        assert_array_equal(
            ma1_eq_ma2.unmasked,
            np.array(
                [
                    [True, False],
                    [True, True],
                ]
            ),
        )
        assert_array_equal(
            ma1_eq_ma2.mask,
            np.array(
                [
                    [False, True],
                    [True, True],
                ]
            ),
        )
        assert ma1_eq_ma2.all()
        assert not (ma1 != ma2).all()
        ma_eq1 = ma1_eq_ma2.all(1)
        assert_array_equal(ma_eq1.mask, np.array([False, True]))
        assert bool(ma_eq1[0]) is True
        assert bool(ma_eq1[1]) is False
        ma_eq0 = ma1_eq_ma2.all(0)
        assert_array_equal(ma_eq0.mask, np.array([False, True]))
        assert bool(ma_eq1[0]) is True
        assert bool(ma_eq1[1]) is False

    @pytest.mark.parametrize("method", ["any", "all"])
    @pytest.mark.parametrize(
        "array,axis",
        [("a", 0), ("a", 1), ("a", None), ("b", None), ("c", 0), ("c", 1), ("c", None)],
    )
    def test_all_and_any(self, array, axis, method):
        ma = getattr(self, "m" + array)
        ma_eq = ma == ma
        ma_all_or_any = getattr(ma_eq, method)(axis=axis)
        filled = ma_eq.unmasked.copy()
        filled[ma_eq.mask] = method == "all"
        a_all_or_any = getattr(filled, method)(axis=axis)
        all_masked = ma.mask.all(axis)
        assert_array_equal(ma_all_or_any.mask, all_masked)
        assert_array_equal(ma_all_or_any.unmasked, a_all_or_any)
        # interpretation as bool
        as_bool = [bool(a) for a in ma_all_or_any.ravel()]
        expected = [bool(a) for a in (a_all_or_any & ~all_masked).ravel()]
        assert as_bool == expected

    def test_any_inplace(self):
        ma_eq = self.ma == self.ma
        expected = ma_eq.any(1)
        out = Masked(np.zeros_like(expected.unmasked))
        result = ma_eq.any(1, out=out)
        assert result is out
        assert_masked_equal(result, expected)

    @pytest.mark.parametrize("method", ("all", "any"))
    @pytest.mark.parametrize("axis", (0, 1, None))
    def test_all_and_any_where(self, method, axis):
        where = np.array(
            [
                [True, False, False],
                [True, True, True],
            ]
        )
        where_final = ~self.ma.mask & where
        ma_eq = self.ma == self.ma
        ma_any = getattr(ma_eq, method)(axis, where=where)
        expected_data = getattr(ma_eq.unmasked, method)(axis, where=where_final)
        expected_mask = np.logical_or.reduce(
            self.ma.mask, axis=axis, where=where_final
        ) | (~where_final).all(axis)
        assert_array_equal(ma_any.unmasked, expected_data)
        assert_array_equal(ma_any.mask, expected_mask)

    @pytest.mark.parametrize("offset", (0, 1))
    def test_diagonal(self, offset):
        mda = self.ma.diagonal(offset=offset)
        expected = Masked(
            self.a.diagonal(offset=offset), self.mask_a.diagonal(offset=offset)
        )
        assert_masked_equal(mda, expected)

    @pytest.mark.parametrize("offset", (0, 1))
    def test_trace(self, offset):
        mta = self.ma.trace(offset=offset)
        expected = Masked(
            self.a.trace(offset=offset), self.mask_a.trace(offset=offset, dtype=bool)
        )
        assert_masked_equal(mta, expected)

    def test_clip(self):
        maclip = self.ma.clip(self.b, self.c)
        expected = Masked(self.a.clip(self.b, self.c), self.mask_a)
        assert_masked_equal(maclip, expected)

    def test_clip_masked_min_max(self):
        maclip = self.ma.clip(self.mb, self.mc)
        # Need to be careful with min, max because of Longitude, which wraps.
        dmax = np.maximum(np.maximum(self.a, self.b), self.c).max()
        dmin = np.minimum(np.minimum(self.a, self.b), self.c).min()
        expected = Masked(
            self.a.clip(self.mb.filled(dmin), self.mc.filled(dmax)), mask=self.mask_a
        )
        assert_masked_equal(maclip, expected)


class TestMaskedQuantityMethods(TestMaskedArrayMethods, QuantitySetup):
    pass


class TestMaskedLongitudeMethods(TestMaskedArrayMethods, LongitudeSetup):
    pass


class TestMaskedArrayProductMethods(MaskedArraySetup):
    # These cannot work on Quantity, so done separately
    @pytest.mark.parametrize("axis", (0, 1, None))
    def test_prod(self, axis):
        ma_sum = self.ma.prod(axis)
        expected_data = self.a.prod(axis)
        expected_mask = self.ma.mask.any(axis)
        assert_array_equal(ma_sum.unmasked, expected_data)
        assert_array_equal(ma_sum.mask, expected_mask)

    @pytest.mark.parametrize("axis", (0, 1, None))
    def test_cumprod(self, axis):
        ma_sum = self.ma.cumprod(axis)
        expected_data = self.a.cumprod(axis)
        mask = self.mask_a
        if axis is None:
            mask = mask.ravel()
        expected_mask = np.logical_or.accumulate(mask, axis=axis)
        assert_array_equal(ma_sum.unmasked, expected_data)
        assert_array_equal(ma_sum.mask, expected_mask)


def test_masked_str_repr_explicit_float():
    sa = np.array([np.pi, 2 * np.pi])
    msa = Masked(sa, [False, True])
    # Test masking  the array works as expected, including truncating digits
    assert str(msa) == "[3.14159265        ———]"
    # Test the digits are kept for scalars.
    assert str(msa[0]) == "3.141592653589793" == str(sa[0])
    # And that the masked string has the same length.
    assert str(msa[1]) == "              ———"
    # Test temporary precision change (which does not affect scalars).
    with np.printoptions(precision=3, floatmode="fixed"):
        assert str(msa) == "[3.142   ———]"
        assert str(msa[0]) == "3.141592653589793" == str(sa[0])
    assert repr(msa) == "MaskedNDArray([3.14159265,        ———])"


def test_masked_str_explicit_string():
    sa = np.array(["2001-02-03", "2002-03-04"])
    msa = Masked(sa, [False, True])
    assert str(msa) == "['2001-02-03'          ———]"
    assert str(msa[0]) == "2001-02-03" == str(sa[0])
    assert str(msa[1]) == "       ———"
    byteorder = "<" if sys.byteorder == "little" else ">"
    repr_ = f"MaskedNDArray(['2001-02-03',          ———], dtype='{byteorder}U10')"
    assert repr(msa) == repr_


@pytest.mark.usefixtures("without_legacy_printoptions")
def test_masked_str_explicit_structured():
    sa = np.array([(1.0, 2.0), (3.0, 4.0)], dtype="f8,f8")
    msa = Masked(sa, [(False, True), (False, False)])
    assert str(msa) == "[(1., ——) (3., 4.)]"
    assert str(msa[0]) == ("(1., ——)" if NUMPY_LT_2_0 else "(1.0, ———)")
    assert str(msa[1]) == str(sa[1]) == ("(3., 4.)" if NUMPY_LT_2_0 else "(3.0, 4.0)")
    with np.printoptions(precision=3, floatmode="fixed"):
        assert str(msa) == "[(1.000,   ———) (3.000, 4.000)]"
        assert str(msa[0]) == ("(1.000,   ———)" if NUMPY_LT_2_0 else "(1.0, ———)")
        assert (
            str(msa[1])
            == str(sa[1])
            == ("(3.000, 4.000)" if NUMPY_LT_2_0 else "(3.0, 4.0)")
        )


def test_masked_repr_explicit_structured():
    # Use explicit endianness to ensure tests pass on all architectures
    sa = np.array([(1.0, 2.0), (3.0, 4.0)], dtype=">f8,>f8")
    msa = Masked(sa, [(False, True), (False, False)])
    assert (
        repr(msa)
        == "MaskedNDArray([(1., ——), (3., 4.)], dtype=[('f0', '>f8'), ('f1', '>f8')])"
    )
    assert (
        repr(msa[0]) == "MaskedNDArray((1., ——), dtype=[('f0', '>f8'), ('f1', '>f8')])"
    )
    assert (
        repr(msa[1]) == "MaskedNDArray((3., 4.), dtype=[('f0', '>f8'), ('f1', '>f8')])"
    )


def test_masked_repr_summary():
    ma = Masked(np.arange(15.0), mask=[True] + [False] * 14)
    with np.printoptions(threshold=2):
        assert repr(ma) == "MaskedNDArray([———,  1.,  2., ..., 12., 13., 14.])"


def test_masked_repr_nodata():
    assert repr(Masked([])) == "MaskedNDArray([], dtype=float64)"


class TestMaskedArrayRepr(MaskedArraySetup):
    def test_array_str(self):
        # very blunt check they work at all.
        str(self.ma)
        str(self.mb)
        str(self.mc)
        str(self.msa)
        str(self.msb)
        str(self.msc)

    def test_scalar_str(self):
        assert self.mb[0].shape == ()
        str(self.mb[0])
        assert self.msb[0].shape == ()
        str(self.msb[0])
        assert self.msc[0].shape == ()
        str(self.msc[0])

    def test_array_repr(self):
        repr(self.ma)
        repr(self.mb)
        repr(self.mc)
        repr(self.msa)
        repr(self.msb)
        repr(self.msc)

    def test_scalar_repr(self):
        repr(self.mb[0])
        repr(self.msb[0])
        repr(self.msc[0])


class TestMaskedQuantityRepr(TestMaskedArrayRepr, QuantitySetup):
    pass


class TestMaskedRecarray(MaskedArraySetup):
    @classmethod
    def setup_class(self):
        super().setup_class()
        self.ra = self.sa.view(np.recarray)
        self.mra = Masked(self.ra, mask=self.mask_sa)

    def test_recarray_setup(self):
        assert isinstance(self.mra, Masked)
        assert isinstance(self.mra, np.recarray)
        assert np.all(self.mra.unmasked == self.ra)
        assert np.all(self.mra.mask == self.mask_sa)
        assert_array_equal(self.mra.view(np.ndarray), self.sa)
        assert isinstance(self.mra.a, Masked)
        assert_array_equal(self.mra.a.unmasked, self.sa["a"])
        assert_array_equal(self.mra.a.mask, self.mask_sa["a"])

    def test_recarray_setting(self):
        mra = self.mra.copy()
        mra.a = self.msa["b"]
        assert_array_equal(mra.a.unmasked, self.msa["b"].unmasked)
        assert_array_equal(mra.a.mask, self.msa["b"].mask)

    @pytest.mark.parametrize("attr", [0, "a"])
    def test_recarray_field_getting(self, attr):
        mra_a = self.mra.field(attr)
        assert isinstance(mra_a, Masked)
        assert_array_equal(mra_a.unmasked, self.sa["a"])
        assert_array_equal(mra_a.mask, self.mask_sa["a"])

    @pytest.mark.parametrize("attr", [0, "a"])
    def test_recarray_field_setting(self, attr):
        mra = self.mra.copy()
        mra.field(attr, self.msa["b"])
        assert_array_equal(mra.a.unmasked, self.msa["b"].unmasked)
        assert_array_equal(mra.a.mask, self.msa["b"].mask)

    def test_recarray_repr(self):
        # Omit dtype part with endian-dependence.
        assert repr(self.mra).startswith(
            "MaskedRecarray([[(———, ———), ( 3.,  4.)],\n"
            "                [(11., ———), (———, 14.)]],\n"
        )

    def test_recarray_represent_as_dict(self):
        rasd = self.mra.info._represent_as_dict()
        assert type(rasd["data"]) is np.ma.MaskedArray
        assert type(rasd["data"].base) is np.ndarray
        mra2 = type(self.mra).info._construct_from_dict(rasd)
        assert type(mra2) is type(self.mra)
        assert_array_equal(mra2.unmasked, self.ra)
        assert_array_equal(mra2.mask, self.mra.mask)


class TestMaskedArrayInteractionWithNumpyMA(MaskedArraySetup):
    def test_masked_array_from_masked(self):
        """Check that we can initialize a MaskedArray properly."""
        np_ma = np.ma.MaskedArray(self.ma)
        assert type(np_ma) is np.ma.MaskedArray
        assert type(np_ma.data) is self._data_cls
        assert type(np_ma.mask) is np.ndarray
        assert_array_equal(np_ma.data, self.a)
        assert_array_equal(np_ma.mask, self.mask_a)

    def test_view_as_masked_array(self):
        """Test that we can be viewed as a MaskedArray."""
        np_ma = self.ma.view(np.ma.MaskedArray)
        assert type(np_ma) is np.ma.MaskedArray
        assert type(np_ma.data) is self._data_cls
        assert type(np_ma.mask) is np.ndarray
        assert_array_equal(np_ma.data, self.a)
        assert_array_equal(np_ma.mask, self.mask_a)


class TestMaskedQuantityInteractionWithNumpyMA(
    TestMaskedArrayInteractionWithNumpyMA, QuantitySetup
):
    pass


@pytest.mark.skipif(not HAS_PLT, reason="requires matplotlib.pyplot")
def test_plt_scatter_masked():
    import matplotlib as mpl

    from astropy.utils import minversion

    MPL_LT_3_8 = not minversion(mpl, "3.8")

    if MPL_LT_3_8:
        pytest.skip("never worked before matplotlib 3.8")

    # check that plotting Masked data doesn't raise an exception
    # see https://github.com/astropy/astropy/issues/12481
    import matplotlib.pyplot as plt

    _, ax = plt.subplots()

    # no mask
    x = Masked([1, 2, 3])
    ax.scatter(x, x, c=x)

    # all masked
    x = Masked([1, 2, 3], mask=True)
    ax.scatter(x, x, c=x)

    # *some* masked
    x = Masked([1, 2, 3], mask=[False, True, False])
    ax.scatter(x, x, c=x)
