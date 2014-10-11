# Licensed under a 3-clause BSD style license - see LICENSE.rst

# TEST_UNICODE_LITERALS

import operator

from distutils import version

import numpy as np

from ...tests.helper import pytest, catch_warnings, assert_follows_unicode_guidelines
from ...utils.exceptions import AstropyDeprecationWarning
from ... import table
from ... import units as u

NUMPY_LT_1P8 = [int(x) for x in np.__version__.split('.')[:2]] < [1, 8]


class TestColumn():

    def test_1(self, Column):
        Column(name='a')

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
        c = np.array([1, 2, 3]).view(Column)
        if Column == table.MaskedColumn:
            assert repr(c) == ('<MaskedColumn name=None unit=None format=None description=None>\n'
                               'masked_array(data = [1 2 3],\n'
                               '             mask = False,\n'
                               '       fill_value = 999999)\n')
        else:
            assert repr(c) == ('<Column name=None unit=None format=None description=None>\n'
                               'array([1, 2, 3])')

    def test_format(self, Column):
        """Show that the formatted output from str() works"""
        from ... import conf
        with conf.set_temp('max_lines', 7):
            c1 = Column(np.arange(2000), name='a', dtype=float,
                        format='%6.2f')
            assert str(c1) == ('   a   \n-------\n   0.00\n'
                               '   1.00\n    ...\n1998.00\n1999.00')

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

        c = Column(data=np.array([1,2,3]) * u.m)
        assert np.all(c.data == np.array([1,2,3]))
        assert np.all(c.unit == u.m)

        c = Column(data=np.array([1,2,3]) * u.m, unit=u.cm)
        assert np.all(c.data == np.array([100,200,300]))
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
    if NUMPY_LT_1P8:
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
    if NUMPY_LT_1P8:
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
