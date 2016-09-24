# Licensed under a 3-clause BSD style license - see LICENSE.rst

# TEST_UNICODE_LITERALS

import operator

import numpy as np

from ...tests.helper import pytest, assert_follows_unicode_guidelines
from ... import table
from ... import units as u
from ...extern import six
from ...utils.compat import NUMPY_LT_1_8


class TestColumn():

    def test_subclass(self, Column):
        c = Column(name='a')
        assert isinstance(c, np.ndarray)
        c2 = c * 2
        assert isinstance(c2, Column)
        assert isinstance(c2, np.ndarray)

    def test_numpy_ops(self, Column):
        """Show that basic numpy operations with Column behave sensibly"""

        arr = np.array([1, 2, 3])
        c = Column(arr, name='a')

        for op, test_equal in ((operator.eq, True),
                               (operator.ne, False),
                               (operator.ge, True),
                               (operator.gt, False),
                               (operator.le, True),
                               (operator.lt, False)):
            for eq in (op(c, arr), op(arr, c)):

                assert np.all(eq) if test_equal else not np.any(eq)
                assert len(eq) == 3
                if Column is table.Column:
                    assert type(eq) == np.ndarray
                else:
                    assert type(eq) == np.ma.core.MaskedArray
                assert eq.dtype.str == '|b1'

        lt = c - 1 < arr
        assert np.all(lt)

    def test_numpy_boolean_ufuncs(self, Column):
        """Show that basic numpy operations with Column behave sensibly"""

        arr = np.array([1, 2, 3])
        c = Column(arr, name='a')

        for ufunc, test_true in ((np.isfinite, True),
                                 (np.isinf, False),
                                 (np.isnan, False),
                                 (np.sign, True),
                                 (np.signbit, False)):
            result = ufunc(c)
            assert len(result) == len(c)
            assert np.all(result) if test_true else not np.any(result)
            if Column is table.Column:
                assert type(result) == np.ndarray
            else:
                assert type(result) == np.ma.core.MaskedArray
                if ufunc is not np.sign:
                    assert result.dtype.str == '|b1'

    def test_view(self, Column):
        c = np.array([1, 2, 3], dtype=np.int64).view(Column)
        assert repr(c) == "<{0} dtype='int64' length=3>\n1\n2\n3".format(Column.__name__)

    def test_format(self, Column):
        """Show that the formatted output from str() works"""
        from ... import conf
        with conf.set_temp('max_lines', 8):
            c1 = Column(np.arange(2000), name='a', dtype=float,
                        format='%6.2f')
            assert str(c1).splitlines() == ['   a   ',
                                            '-------',
                                            '   0.00',
                                            '   1.00',
                                            '    ...',
                                            '1998.00',
                                            '1999.00',
                                            'Length = 2000 rows']

    def test_convert_numpy_array(self, Column):
        d = Column([1, 2, 3], name='a', dtype='i8')

        np_data = np.array(d)
        assert np.all(np_data == d)
        np_data = np.array(d, copy=False)
        assert np.all(np_data == d)
        np_data = np.array(d, dtype='i4')
        assert np.all(np_data == d)

    def test_convert_unit(self, Column):
        d = Column([1, 2, 3], name='a', dtype="f8", unit="m")
        d.convert_unit_to("km")
        assert np.all(d.data == [0.001, 0.002, 0.003])

    def test_array_wrap(self):
        """Test that the __array_wrap__ method converts a reduction ufunc
        output that has a different shape into an ndarray view.  Without this a
        method call like c.mean() returns a Column array object with length=1."""
        # Mean and sum for a 1-d float column
        c = table.Column(name='a', data=[1., 2., 3.])
        assert np.allclose(c.mean(), 2.0)
        assert isinstance(c.mean(), (np.floating, float))
        assert np.allclose(c.sum(), 6.)
        assert isinstance(c.sum(), (np.floating, float))

        # Non-reduction ufunc preserves Column class
        assert isinstance(np.cos(c), table.Column)

        # Sum for a 1-d int column
        c = table.Column(name='a', data=[1, 2, 3])
        assert np.allclose(c.sum(), 6)
        assert isinstance(c.sum(), (np.integer, int))

        # Sum for a 2-d int column
        c = table.Column(name='a', data=[[1, 2, 3],
                                         [4, 5, 6]])
        assert c.sum() == 21
        assert isinstance(c.sum(), (np.integer, int))
        assert np.all(c.sum(axis=0) == [5, 7, 9])
        assert c.sum(axis=0).shape == (3,)
        assert isinstance(c.sum(axis=0), np.ndarray)

        # Sum and mean for a 1-d masked column
        c = table.MaskedColumn(name='a', data=[1., 2., 3.], mask=[0, 0, 1])
        assert np.allclose(c.mean(), 1.5)
        assert isinstance(c.mean(), (np.floating, float))
        assert np.allclose(c.sum(), 3.)
        assert isinstance(c.sum(), (np.floating, float))

    def test_name_none(self, Column):
        """Can create a column without supplying name, which defaults to None"""
        c = Column([1, 2])
        assert c.name is None
        assert np.all(c == np.array([1, 2]))

    def test_quantity_init(self, Column):

        c = Column(data=np.array([1, 2, 3]) * u.m)
        assert np.all(c.data == np.array([1, 2, 3]))
        assert np.all(c.unit == u.m)

        c = Column(data=np.array([1, 2, 3]) * u.m, unit=u.cm)
        assert np.all(c.data == np.array([100, 200, 300]))
        assert np.all(c.unit == u.cm)

    def test_attrs_survive_getitem_after_change(self, Column):
        """
        Test for issue #3023: when calling getitem with a MaskedArray subclass
        the original object attributes are not copied.
        """
        c1 = Column([1, 2, 3], name='a', unit='m', format='i',
                    description='aa', meta={'a': 1})
        c1.name = 'b'
        c1.unit = 'km'
        c1.format = 'i2'
        c1.description = 'bb'
        c1.meta = {'bbb': 2}

        for item in (slice(None, None), slice(None, 1), np.array([0, 2]),
                     np.array([False, True, False])):
            c2 = c1[item]
            assert c2.name == 'b'
            assert c2.unit is u.km
            assert c2.format == 'i2'
            assert c2.description == 'bb'
            assert c2.meta == {'bbb': 2}

        # Make sure that calling getitem resulting in a scalar does
        # not copy attributes.
        val = c1[1]
        for attr in ('name', 'unit', 'format', 'description', 'meta'):
            assert not hasattr(val, attr)

    def test_to_quantity(self, Column):
        d = Column([1, 2, 3], name='a', dtype="f8", unit="m")

        assert np.all(d.quantity == ([1, 2, 3.] * u.m))
        assert np.all(d.quantity.value == ([1, 2, 3.] * u.m).value)
        assert np.all(d.quantity == d.to('m'))
        assert np.all(d.quantity.value == d.to('m').value)

        np.testing.assert_allclose(d.to(u.km).value, ([.001, .002, .003] * u.km).value)
        np.testing.assert_allclose(d.to('km').value, ([.001, .002, .003] * u.km).value)

        np.testing.assert_allclose(d.to(u.MHz,u.equivalencies.spectral()).value,
                                   [299.792458, 149.896229,  99.93081933])

        d_nounit = Column([1, 2, 3], name='a', dtype="f8", unit=None)
        with pytest.raises(u.UnitsError):
            d_nounit.to(u.km)
        assert np.all(d_nounit.to(u.dimensionless_unscaled) == np.array([1, 2, 3]))

        #make sure the correct copy/no copy behavior is happening
        q = [1, 3, 5]*u.km

        # to should always make a copy
        d.to(u.km)[:] = q
        np.testing.assert_allclose(d, [1, 2, 3])

        # explcit copying of the quantity should not change the column
        d.quantity.copy()[:] = q
        np.testing.assert_allclose(d, [1, 2, 3])

        # but quantity directly is a "view", accessing the underlying column
        d.quantity[:] = q
        np.testing.assert_allclose(d, [1000, 3000, 5000])

        #view should also work for integers
        d2 = Column([1, 2, 3], name='a', dtype=int, unit="m")
        d2.quantity[:] = q
        np.testing.assert_allclose(d2, [1000, 3000, 5000])

        #but it should fail for strings or other non-numeric tables
        d3 = Column(['arg', 'name', 'stuff'], name='a', unit="m")
        with pytest.raises(TypeError):
            d3.quantity

    def test_item_access_type(self, Column):
        """
        Tests for #3095, which forces integer item access to always return a plain
        ndarray or MaskedArray, even in the case of a multi-dim column.
        """
        integer_types = (int, long, np.int) if six.PY2 else (int, np.int)

        for int_type in integer_types:
            c = Column([[1, 2], [3, 4]])
            i0 = int_type(0)
            i1 = int_type(1)
            assert np.all(c[i0] == [1, 2])
            assert type(c[i0]) == (np.ma.MaskedArray if hasattr(Column, 'mask') else np.ndarray)
            assert c[i0].shape == (2,)

            c01 = c[i0:i1]
            assert np.all(c01 == [[1, 2]])
            assert isinstance(c01, Column)
            assert c01.shape == (1, 2)

            c = Column([1, 2])
            assert np.all(c[i0] == 1)
            assert isinstance(c[i0], np.integer)
            assert c[i0].shape == ()

            c01 = c[i0:i1]
            assert np.all(c01 == [1])
            assert isinstance(c01, Column)
            assert c01.shape == (1,)

    def test_insert_basic(self, Column):
        c = Column([0, 1, 2], name='a', dtype=int, unit='mJy', format='%i',
                   description='test column', meta={'c': 8, 'd': 12})

        # Basic insert
        c1 = c.insert(1, 100)
        assert np.all(c1 == [0, 100, 1, 2])
        assert c1.attrs_equal(c)
        assert type(c) is type(c1)
        if hasattr(c1, 'mask'):
            assert c1.data.shape == c1.mask.shape

        c1 = c.insert(-1, 100)
        assert np.all(c1 == [0, 1, 100, 2])

        c1 = c.insert(3, 100)
        assert np.all(c1 == [0, 1, 2, 100])

        c1 = c.insert(-3, 100)
        assert np.all(c1 == [100, 0, 1, 2])

        c1 = c.insert(1, [100, 200, 300])
        if hasattr(c1, 'mask'):
            assert c1.data.shape == c1.mask.shape

        # Out of bounds index
        with pytest.raises((ValueError, IndexError)):
            c1 = c.insert(-4, 100)
        with pytest.raises((ValueError,IndexError)):
            c1 = c.insert(4, 100)

    def test_insert_multidim(self, Column):
        c = Column([[1, 2],
                    [3, 4]], name='a', dtype=int)

        # Basic insert
        c1 = c.insert(1, [100, 200])
        assert np.all(c1 == [[1, 2], [100, 200], [3, 4]])

        # Broadcast
        c1 = c.insert(1, 100)
        assert np.all(c1 == [[1, 2], [100, 100], [3, 4]])

        # Wrong shape
        with pytest.raises(ValueError):
            c1 = c.insert(1, [100, 200, 300])

    def test_insert_object(self, Column):
        c = Column(['a', 1, None], name='a', dtype=object)

        # Basic insert
        c1 = c.insert(1, [100, 200])
        assert np.all(c1 == ['a', [100, 200], 1, None])

    def test_insert_masked(self):
        c = table.MaskedColumn([0, 1, 2], name='a', mask=[False, True, False])

        # Basic insert
        c1 = c.insert(1, 100)
        assert np.all(c1.data.data == [0, 100, 1, 2])
        assert np.all(c1.data.mask == [False, False, True, False])
        assert type(c) is type(c1)

        for mask in (False, True):
            c1 = c.insert(1, 100, mask=mask)
            assert np.all(c1.data.data == [0, 100, 1, 2])
            assert np.all(c1.data.mask == [False, mask, True, False])

    def test_insert_masked_multidim(self):
        c = table.MaskedColumn([[1, 2],
                                [3, 4]], name='a', dtype=int)

        c1 = c.insert(1, [100, 200], mask=True)
        assert np.all(c1.data.data == [[1, 2], [100, 200], [3, 4]])
        assert np.all(c1.data.mask == [[False, False], [True, True], [False, False]])

        c1 = c.insert(1, [100, 200], mask=[True, False])
        assert np.all(c1.data.data == [[1, 2], [100, 200], [3, 4]])
        assert np.all(c1.data.mask == [[False, False], [True, False], [False, False]])

        with pytest.raises(ValueError):
            c1 = c.insert(1, [100, 200], mask=[True, False, True])

    def test_mask_on_non_masked_table(self):
        """
        When table is not masked and trying to set mask on column then
        it's Raise AttributeError.
        """

        t = table.Table([[1, 2], [3, 4]], names=('a', 'b'), dtype=('i4', 'f8'))

        with pytest.raises(AttributeError):
            t['a'].mask = [True, False]


class TestAttrEqual():
    """Bunch of tests originally from ATpy that test the attrs_equal method."""

    def test_5(self, Column):
        c1 = Column(name='a', dtype=int, unit='mJy')
        c2 = Column(name='a', dtype=int, unit='mJy')
        assert c1.attrs_equal(c2)

    def test_6(self, Column):
        c1 = Column(name='a', dtype=int, unit='mJy', format='%i',
                    description='test column', meta={'c': 8, 'd': 12})
        c2 = Column(name='a', dtype=int, unit='mJy', format='%i',
                    description='test column', meta={'c': 8, 'd': 12})
        assert c1.attrs_equal(c2)

    def test_7(self, Column):
        c1 = Column(name='a', dtype=int, unit='mJy', format='%i',
                    description='test column', meta={'c': 8, 'd': 12})
        c2 = Column(name='b', dtype=int, unit='mJy', format='%i',
                    description='test column', meta={'c': 8, 'd': 12})
        assert not c1.attrs_equal(c2)

    def test_8(self, Column):
        c1 = Column(name='a', dtype=int, unit='mJy', format='%i',
                    description='test column', meta={'c': 8, 'd': 12})
        c2 = Column(name='a', dtype=float, unit='mJy', format='%i',
                    description='test column', meta={'c': 8, 'd': 12})
        assert not c1.attrs_equal(c2)

    def test_9(self, Column):
        c1 = Column(name='a', dtype=int, unit='mJy', format='%i',
                    description='test column', meta={'c': 8, 'd': 12})
        c2 = Column(name='a', dtype=int, unit='erg.cm-2.s-1.Hz-1', format='%i',
                    description='test column', meta={'c': 8, 'd': 12})
        assert not c1.attrs_equal(c2)

    def test_10(self, Column):
        c1 = Column(name='a', dtype=int, unit='mJy', format='%i',
                    description='test column', meta={'c': 8, 'd': 12})
        c2 = Column(name='a', dtype=int, unit='mJy', format='%g',
                    description='test column', meta={'c': 8, 'd': 12})
        assert not c1.attrs_equal(c2)

    def test_11(self, Column):
        c1 = Column(name='a', dtype=int, unit='mJy', format='%i',
                    description='test column', meta={'c': 8, 'd': 12})
        c2 = Column(name='a', dtype=int, unit='mJy', format='%i',
                    description='another test column', meta={'c': 8, 'd': 12})
        assert not c1.attrs_equal(c2)

    def test_12(self, Column):
        c1 = Column(name='a', dtype=int, unit='mJy', format='%i',
                    description='test column', meta={'c': 8, 'd': 12})
        c2 = Column(name='a', dtype=int, unit='mJy', format='%i',
                    description='test column', meta={'e': 8, 'd': 12})
        assert not c1.attrs_equal(c2)

    def test_13(self, Column):
        c1 = Column(name='a', dtype=int, unit='mJy', format='%i',
                    description='test column', meta={'c': 8, 'd': 12})
        c2 = Column(name='a', dtype=int, unit='mJy', format='%i',
                    description='test column', meta={'c': 9, 'd': 12})
        assert not c1.attrs_equal(c2)

    def test_col_and_masked_col(self):
        c1 = table.Column(name='a', dtype=int, unit='mJy', format='%i',
                          description='test column', meta={'c': 8, 'd': 12})
        c2 = table.MaskedColumn(name='a', dtype=int, unit='mJy', format='%i',
                                description='test column', meta={'c': 8, 'd': 12})
        assert c1.attrs_equal(c2)
        assert c2.attrs_equal(c1)

# Check that the meta descriptor is working as expected. The MetaBaseTest class
# takes care of defining all the tests, and we simply have to define the class
# and any minimal set of args to pass.

from ...utils.tests.test_metadata import MetaBaseTest


class TestMetaColumn(MetaBaseTest):
    test_class = table.Column
    args = ()


class TestMetaMaskedColumn(MetaBaseTest):
    test_class = table.MaskedColumn
    args = ()


def test_getitem_metadata_regression():
    """
    Regression test for #1471: MaskedArray does not call __array_finalize__ so
    the meta-data was not getting copied over. By overloading _update_from we
    are able to work around this bug.
    """

    # Make sure that meta-data gets propagated with __getitem__

    c = table.Column(data=[1,2], name='a', description='b', unit='m', format="%i", meta={'c': 8})
    assert c[1:2].name == 'a'
    assert c[1:2].description == 'b'
    assert c[1:2].unit == 'm'
    assert c[1:2].format == '%i'
    assert c[1:2].meta['c'] == 8

    c = table.MaskedColumn(data=[1,2], name='a', description='b', unit='m', format="%i", meta={'c': 8})
    assert c[1:2].name == 'a'
    assert c[1:2].description == 'b'
    assert c[1:2].unit == 'm'
    assert c[1:2].format == '%i'
    assert c[1:2].meta['c'] == 8

    # As above, but with take() - check the method and the function

    c = table.Column(data=[1,2,3], name='a', description='b', unit='m', format="%i", meta={'c': 8})
    for subset in [c.take([0, 1]), np.take(c, [0, 1])]:
        assert subset.name == 'a'
        assert subset.description == 'b'
        assert subset.unit == 'm'
        assert subset.format == '%i'
        assert subset.meta['c'] == 8

    # Metadata isn't copied for scalar values
    if NUMPY_LT_1_8:
        with pytest.raises(ValueError):
            c.take(0)
        with pytest.raises(ValueError):
            np.take(c, 0)
    else:
        for subset in [c.take(0), np.take(c, 0)]:
            assert subset == 1
            assert subset.shape == ()
            assert not isinstance(subset, table.Column)

    c = table.MaskedColumn(data=[1,2,3], name='a', description='b', unit='m', format="%i", meta={'c': 8})
    for subset in [c.take([0, 1]), np.take(c, [0, 1])]:
        assert subset.name == 'a'
        assert subset.description == 'b'
        assert subset.unit == 'm'
        assert subset.format == '%i'
        assert subset.meta['c'] == 8

    # Metadata isn't copied for scalar values
    if NUMPY_LT_1_8:
        with pytest.raises(ValueError):
            c.take(0)
        with pytest.raises(ValueError):
            np.take(c, 0)
    else:
        for subset in [c.take(0), np.take(c, 0)]:
            assert subset == 1
            assert subset.shape == ()
            assert not isinstance(subset, table.MaskedColumn)


def test_unicode_guidelines():
    arr = np.array([1, 2, 3])
    c = table.Column(arr, name='a')

    assert_follows_unicode_guidelines(c)


def test_scalar_column():
    """
    Column is not designed to hold scalars, but for numpy 1.6 this can happen:

      >> type(np.std(table.Column([1, 2])))
      astropy.table.column.Column
    """
    c = table.Column(1.5)
    assert repr(c) == '1.5'
    assert str(c) == '1.5'


def test_qtable_column_conversion():
    """
    Ensures that a QTable that gets assigned a unit switches to be Quantity-y
    """
    qtab = table.QTable([[1, 2], [3, 4.2]], names=['i', 'f'])

    assert isinstance(qtab['i'], table.column.Column)
    assert isinstance(qtab['f'], table.column.Column)

    qtab['i'].unit = 'km/s'
    assert isinstance(qtab['i'], u.Quantity)
    assert isinstance(qtab['f'], table.column.Column)

    # should follow from the above, but good to make sure as a #4497 regression test
    assert isinstance(qtab['i'][0], u.Quantity)
    assert isinstance(qtab[0]['i'], u.Quantity)
    assert not isinstance(qtab['f'][0], u.Quantity)
    assert not isinstance(qtab[0]['f'], u.Quantity)

    # Regression test for #5342: if a function unit is assigned, the column
    # should become the appropriate FunctionQuantity subclass.
    qtab['f'].unit = u.dex(u.cm/u.s**2)
    assert isinstance(qtab['f'], u.Dex)
