# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

import operator

import pytest
import numpy as np
from numpy.testing import assert_array_equal

from astropy.tests.helper import assert_follows_unicode_guidelines, catch_warnings
from astropy import table
from astropy import units as u


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
        assert repr(c) == f"<{Column.__name__} dtype='int64' length=3>\n1\n2\n3"

    def test_format(self, Column):
        """Show that the formatted output from str() works"""
        from astropy import conf
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

    def test_quantity_comparison(self, Column):
        # regression test for gh-6532
        c = Column([1, 2100, 3], unit='Hz')
        q = 2 * u.kHz
        check = c < q
        assert np.all(check == [True, False, True])
        # This already worked, but just in case.
        check = q >= c
        assert np.all(check == [True, False, True])

    def test_attrs_survive_getitem_after_change(self, Column):
        """
        Test for issue #3023: when calling getitem with a MaskedArray subclass
        the original object attributes are not copied.
        """
        c1 = Column([1, 2, 3], name='a', unit='m', format='%i',
                    description='aa', meta={'a': 1})
        c1.name = 'b'
        c1.unit = 'km'
        c1.format = '%d'
        c1.description = 'bb'
        c1.meta = {'bbb': 2}

        for item in (slice(None, None), slice(None, 1), np.array([0, 2]),
                     np.array([False, True, False])):
            c2 = c1[item]
            assert c2.name == 'b'
            assert c2.unit is u.km
            assert c2.format == '%d'
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

        np.testing.assert_allclose(d.to(u.MHz, u.equivalencies.spectral()).value,
                                   [299.792458, 149.896229, 99.93081933])

        d_nounit = Column([1, 2, 3], name='a', dtype="f8", unit=None)
        with pytest.raises(u.UnitsError):
            d_nounit.to(u.km)
        assert np.all(d_nounit.to(u.dimensionless_unscaled) == np.array([1, 2, 3]))

        # make sure the correct copy/no copy behavior is happening
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

        # view should also work for integers
        d2 = Column([1, 2, 3], name='a', dtype=int, unit="m")
        d2.quantity[:] = q
        np.testing.assert_allclose(d2, [1000, 3000, 5000])

        # but it should fail for strings or other non-numeric tables
        d3 = Column(['arg', 'name', 'stuff'], name='a', unit="m")
        with pytest.raises(TypeError):
            d3.quantity

    def test_to_funcunit_quantity(self, Column):
        """
        Tests for #8424, check if function-unit can be retrieved from column.
        """
        d = Column([1, 2, 3], name='a', dtype="f8", unit="dex(AA)")

        assert np.all(d.quantity == ([1, 2, 3] * u.dex(u.AA)))
        assert np.all(d.quantity.value == ([1, 2, 3] * u.dex(u.AA)).value)
        assert np.all(d.quantity == d.to("dex(AA)"))
        assert np.all(d.quantity.value == d.to("dex(AA)").value)

        # make sure, casting to linear unit works
        q = [10, 100, 1000] * u.AA
        np.testing.assert_allclose(d.to(u.AA), q)

    def test_item_access_type(self, Column):
        """
        Tests for #3095, which forces integer item access to always return a plain
        ndarray or MaskedArray, even in the case of a multi-dim column.
        """
        integer_types = (int, np.int_)

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
        with pytest.raises((ValueError, IndexError)):
            c1 = c.insert(4, 100)

    def test_insert_axis(self, Column):
        """Insert with non-default axis kwarg"""
        c = Column([[1, 2], [3, 4]])

        c1 = c.insert(1, [5, 6], axis=None)
        assert np.all(c1 == [1, 5, 6, 2, 3, 4])

        c1 = c.insert(1, [5, 6], axis=1)
        assert np.all(c1 == [[1, 5, 2], [3, 6, 4]])

    def test_insert_string_expand(self, Column):
        c = Column(['a', 'b'])
        c1 = c.insert(0, 'abc')
        assert np.all(c1 == ['abc', 'a', 'b'])

        c = Column(['a', 'b'])
        c1 = c.insert(0, ['c', 'def'])
        assert np.all(c1 == ['c', 'def', 'a', 'b'])

    def test_insert_string_masked_values(self):
        c = table.MaskedColumn(['a', 'b'])
        c1 = c.insert(0, np.ma.masked)
        assert np.all(c1 == ['', 'a', 'b'])
        assert np.all(c1.mask == [True, False, False])
        assert c1.dtype == 'U1'
        c2 = c.insert(1, np.ma.MaskedArray(['ccc', 'dd'], mask=[True, False]))
        assert np.all(c2 == ['a', 'ccc', 'dd', 'b'])
        assert np.all(c2.mask == [False, True, False, False])
        assert c2.dtype == 'U3'

    def test_insert_string_type_error(self, Column):
        c = Column([1, 2])
        with pytest.raises(ValueError, match='invalid literal for int'):
            c.insert(0, 'string')

        c = Column(['a', 'b'])
        with pytest.raises(TypeError, match='string operation on non-string array'):
            c.insert(0, 1)

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
        assert np.all(c1 == np.array(['a', [100, 200], 1, None],
                                     dtype=object))

    def test_insert_masked(self):
        c = table.MaskedColumn([0, 1, 2], name='a', fill_value=9999,
                               mask=[False, True, False])

        # Basic insert
        c1 = c.insert(1, 100)
        assert np.all(c1.data.data == [0, 100, 1, 2])
        assert c1.fill_value == 9999
        assert np.all(c1.data.mask == [False, False, True, False])
        assert type(c) is type(c1)

        for mask in (False, True):
            c1 = c.insert(1, 100, mask=mask)
            assert np.all(c1.data.data == [0, 100, 1, 2])
            assert np.all(c1.data.mask == [False, mask, True, False])

    def test_masked_multidim_as_list(self):
        data = np.ma.MaskedArray([1, 2], mask=[True, False])
        c = table.MaskedColumn([data])
        assert c.shape == (1, 2)
        assert np.all(c[0].mask == [True, False])

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


from astropy.utils.tests.test_metadata import MetaBaseTest


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

    c = table.Column(data=[1, 2], name='a', description='b', unit='m', format="%i", meta={'c': 8})
    assert c[1:2].name == 'a'
    assert c[1:2].description == 'b'
    assert c[1:2].unit == 'm'
    assert c[1:2].format == '%i'
    assert c[1:2].meta['c'] == 8

    c = table.MaskedColumn(data=[1, 2], name='a', description='b', unit='m', format="%i", meta={'c': 8})
    assert c[1:2].name == 'a'
    assert c[1:2].description == 'b'
    assert c[1:2].unit == 'm'
    assert c[1:2].format == '%i'
    assert c[1:2].meta['c'] == 8

    # As above, but with take() - check the method and the function

    c = table.Column(data=[1, 2, 3], name='a', description='b', unit='m', format="%i", meta={'c': 8})
    for subset in [c.take([0, 1]), np.take(c, [0, 1])]:
        assert subset.name == 'a'
        assert subset.description == 'b'
        assert subset.unit == 'm'
        assert subset.format == '%i'
        assert subset.meta['c'] == 8

    # Metadata isn't copied for scalar values
    for subset in [c.take(0), np.take(c, 0)]:
        assert subset == 1
        assert subset.shape == ()
        assert not isinstance(subset, table.Column)

    c = table.MaskedColumn(data=[1, 2, 3], name='a', description='b', unit='m', format="%i", meta={'c': 8})
    for subset in [c.take([0, 1]), np.take(c, [0, 1])]:
        assert subset.name == 'a'
        assert subset.description == 'b'
        assert subset.unit == 'm'
        assert subset.format == '%i'
        assert subset.meta['c'] == 8

    # Metadata isn't copied for scalar values
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


@pytest.mark.parametrize('masked', [True, False])
def test_string_truncation_warning(masked):
    """
    Test warnings associated with in-place assignment to a string
    column that results in truncation of the right hand side.
    """
    t = table.Table([['aa', 'bb']], names=['a'], masked=masked)

    with catch_warnings() as w:
        from inspect import currentframe, getframeinfo
        t['a'][1] = 'cc'
        assert len(w) == 0

        t['a'][:] = 'dd'
        assert len(w) == 0

    with catch_warnings() as w:
        frameinfo = getframeinfo(currentframe())
        t['a'][0] = 'eee'  # replace item with string that gets truncated
        assert t['a'][0] == 'ee'
        assert len(w) == 1
        assert ('truncated right side string(s) longer than 2 character(s)'
                in str(w[0].message))

        # Make sure the warning points back to the user code line
        assert w[0].lineno == frameinfo.lineno + 1
        assert w[0].category is table.StringTruncateWarning
        assert 'test_column' in w[0].filename

    with catch_warnings() as w:
        t['a'][:] = ['ff', 'ggg']  # replace item with string that gets truncated
        assert np.all(t['a'] == ['ff', 'gg'])
        assert len(w) == 1
        assert ('truncated right side string(s) longer than 2 character(s)'
                in str(w[0].message))

    with catch_warnings() as w:
        # Test the obscure case of assigning from an array that was originally
        # wider than any of the current elements (i.e. dtype is U4 but actual
        # elements are U1 at the time of assignment).
        val = np.array(['ffff', 'gggg'])
        val[:] = ['f', 'g']
        t['a'][:] = val
        assert np.all(t['a'] == ['f', 'g'])
        assert len(w) == 0


def test_string_truncation_warning_masked():
    """
    Test warnings associated with in-place assignment to a string
    to a masked column, specifically where the right hand side
    contains np.ma.masked.
    """

    # Test for strings, but also cover assignment of np.ma.masked to
    # int and float masked column setting.  This was previously only
    # covered in an unrelated io.ascii test (test_line_endings) which
    # showed an unexpected difference between handling of str and numeric
    # masked arrays.
    for values in (['a', 'b'], [1, 2], [1.0, 2.0]):
        mc = table.MaskedColumn(values)

        with catch_warnings() as w:
            mc[1] = np.ma.masked
            assert len(w) == 0
            assert np.all(mc.mask == [False, True])

            mc[:] = np.ma.masked
            assert len(w) == 0
            assert np.all(mc.mask == [True, True])

    mc = table.MaskedColumn(['aa', 'bb'])

    with catch_warnings() as w:
        mc[:] = [np.ma.masked, 'ggg']  # replace item with string that gets truncated
        assert mc[1] == 'gg'
        assert np.all(mc.mask == [True, False])
        assert len(w) == 1
        assert ('truncated right side string(s) longer than 2 character(s)'
                in str(w[0].message))


@pytest.mark.parametrize('Column', (table.Column, table.MaskedColumn))
def test_col_unicode_sandwich_create_from_str(Column):
    """
    Create a bytestring Column from strings (including unicode) in Py3.
    """
    # a-umlaut is a 2-byte character in utf-8, test fails with ascii encoding.
    # Stress the system by injecting non-ASCII characters.
    uba = 'bä'
    c = Column([uba, 'def'], dtype='S')
    assert c.dtype.char == 'S'
    assert c[0] == uba
    assert isinstance(c[0], str)
    assert isinstance(c[:0], table.Column)
    assert np.all(c[:2] == np.array([uba, 'def']))


@pytest.mark.parametrize('Column', (table.Column, table.MaskedColumn))
def test_col_unicode_sandwich_bytes_obj(Column):
    """
    Create a Column of dtype object with bytestring in it and make sure
    it keeps the bytestring and not convert to str with accessed.
    """
    c = Column([None, b'def'])
    assert c.dtype.char == 'O'
    assert not c[0]
    assert c[1] == b'def'
    assert isinstance(c[1], bytes)
    assert not isinstance(c[1], str)
    assert isinstance(c[:0], table.Column)
    assert np.all(c[:2] == np.array([None, b'def']))
    assert not np.all(c[:2] == np.array([None, 'def']))


@pytest.mark.parametrize('Column', (table.Column, table.MaskedColumn))
def test_col_unicode_sandwich_bytes(Column):
    """
    Create a bytestring Column from bytes and ensure that it works in Python 3 in
    a convenient way like in Python 2.
    """
    # a-umlaut is a 2-byte character in utf-8, test fails with ascii encoding.
    # Stress the system by injecting non-ASCII characters.
    uba = 'bä'
    uba8 = uba.encode('utf-8')
    c = Column([uba8, b'def'])
    assert c.dtype.char == 'S'
    assert c[0] == uba
    assert isinstance(c[0], str)
    assert isinstance(c[:0], table.Column)
    assert np.all(c[:2] == np.array([uba, 'def']))

    assert isinstance(c[:], table.Column)
    assert c[:].dtype.char == 'S'

    # Array / list comparisons
    assert np.all(c == [uba, 'def'])

    ok = c == [uba8, b'def']
    assert type(ok) is type(c.data)
    assert ok.dtype.char == '?'
    assert np.all(ok)

    assert np.all(c == np.array([uba, 'def']))
    assert np.all(c == np.array([uba8, b'def']))

    # Scalar compare
    cmps = (uba, uba8)
    for cmp in cmps:
        ok = c == cmp
        assert type(ok) is type(c.data)
        assert np.all(ok == [True, False])


def test_col_unicode_sandwich_unicode():
    """
    Sanity check that Unicode Column behaves normally.
    """
    # On Py2 the unicode must be ASCII-compatible, else the final test fails.
    uba = 'bä'
    uba8 = uba.encode('utf-8')

    c = table.Column([uba, 'def'], dtype='U')
    assert c[0] == uba
    assert isinstance(c[:0], table.Column)
    assert isinstance(c[0], str)
    assert np.all(c[:2] == np.array([uba, 'def']))

    assert isinstance(c[:], table.Column)
    assert c[:].dtype.char == 'U'

    ok = c == [uba, 'def']
    assert type(ok) == np.ndarray
    assert ok.dtype.char == '?'
    assert np.all(ok)

    assert np.all(c != [uba8, b'def'])


def test_masked_col_unicode_sandwich():
    """
    Create a bytestring MaskedColumn and ensure that it works in Python 3 in
    a convenient way like in Python 2.
    """
    c = table.MaskedColumn([b'abc', b'def'])
    c[1] = np.ma.masked
    assert isinstance(c[:0], table.MaskedColumn)
    assert isinstance(c[0], str)

    assert c[0] == 'abc'
    assert c[1] is np.ma.masked

    assert isinstance(c[:], table.MaskedColumn)
    assert c[:].dtype.char == 'S'

    ok = c == ['abc', 'def']
    assert ok[0] == True
    assert ok[1] is np.ma.masked
    assert np.all(c == [b'abc', b'def'])
    assert np.all(c == np.array(['abc', 'def']))
    assert np.all(c == np.array([b'abc', b'def']))

    for cmp in ('abc', b'abc'):
        ok = c == cmp
        assert type(ok) is np.ma.MaskedArray
        assert ok[0] == True
        assert ok[1] is np.ma.masked


@pytest.mark.parametrize('Column', (table.Column, table.MaskedColumn))
def test_unicode_sandwich_set(Column):
    """
    Test setting
    """
    uba = 'bä'

    c = Column([b'abc', b'def'])

    c[0] = b'aa'
    assert np.all(c == ['aa', 'def'])

    c[0] = uba  # a-umlaut is a 2-byte character in utf-8, test fails with ascii encoding
    assert np.all(c == [uba, 'def'])
    assert c.pformat() == ['None', '----', '  ' + uba, ' def']

    c[:] = b'cc'
    assert np.all(c == ['cc', 'cc'])

    c[:] = uba
    assert np.all(c == [uba, uba])

    c[:] = ''
    c[:] = [uba, b'def']
    assert np.all(c == [uba, b'def'])


@pytest.mark.parametrize('class1', [table.MaskedColumn, table.Column])
@pytest.mark.parametrize('class2', [table.MaskedColumn, table.Column, str, list])
def test_unicode_sandwich_compare(class1, class2):
    """Test that comparing a bytestring Column/MaskedColumn with various
    str (unicode) object types gives the expected result.  Tests #6838.
    """
    obj1 = class1([b'a', b'c'])
    if class2 is str:
        obj2 = 'a'
    elif class2 is list:
        obj2 = ['a', 'b']
    else:
        obj2 = class2(['a', 'b'])

    assert np.all((obj1 == obj2) == [True, False])
    assert np.all((obj2 == obj1) == [True, False])

    assert np.all((obj1 != obj2) == [False, True])
    assert np.all((obj2 != obj1) == [False, True])

    assert np.all((obj1 > obj2) == [False, True])
    assert np.all((obj2 > obj1) == [False, False])

    assert np.all((obj1 <= obj2) == [True, False])
    assert np.all((obj2 <= obj1) == [True, True])

    assert np.all((obj1 < obj2) == [False, False])
    assert np.all((obj2 < obj1) == [False, True])

    assert np.all((obj1 >= obj2) == [True, True])
    assert np.all((obj2 >= obj1) == [True, False])


def test_unicode_sandwich_masked_compare():
    """Test the fix for #6839 from #6899."""
    c1 = table.MaskedColumn(['a', 'b', 'c', 'd'],
                            mask=[True, False, True, False])
    c2 = table.MaskedColumn([b'a', b'b', b'c', b'd'],
                            mask=[True, True, False, False])

    for cmp in ((c1 == c2), (c2 == c1)):
        assert cmp[0] is np.ma.masked
        assert cmp[1] is np.ma.masked
        assert cmp[2] is np.ma.masked
        assert cmp[3]

    for cmp in ((c1 != c2), (c2 != c1)):
        assert cmp[0] is np.ma.masked
        assert cmp[1] is np.ma.masked
        assert cmp[2] is np.ma.masked
        assert not cmp[3]

    # Note: comparisons <, >, >=, <= fail to return a masked array entirely,
    # see https://github.com/numpy/numpy/issues/10092.


def test_structured_masked_column_roundtrip():
    mc = table.MaskedColumn([(1., 2.), (3., 4.)],
                            mask=[(False, False), (False, False)], dtype='f8,f8')
    assert len(mc.dtype.fields) == 2
    mc2 = table.MaskedColumn(mc)
    assert_array_equal(mc2, mc)
