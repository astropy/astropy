# Licensed under a 3-clause BSD style license - see LICENSE.rst

import operator

from distutils import version

import numpy as np

from ...tests.helper import pytest, catch_warnings
from ...utils.exceptions import AstropyDeprecationWarning
from ... import table
from ... import units as u

@pytest.fixture(params=[table.Column, table.MaskedColumn])
def Column(request):
    # Fixture to run all the Column tests for both an unmasked (ndarray)
    # and masked (MaskedArray) column.
    return request.param


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
        MAX_LINES_val = table.pprint.MAX_LINES()
        table.pprint.MAX_LINES.set(7)
        c1 = Column(np.arange(2000), name='a', dtype=float,
                    format='%6.2f')
        assert str(c1) == ('   a   \n-------\n   0.00\n'
                           '   1.00\n    ...\n1998.00\n1999.00')
        table.pprint.MAX_LINES.set(MAX_LINES_val)

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

    def test_deprecated_attributes(self, Column, recwarn):
        d = Column([1, 2, 3], name='a', dtype="f8", unit="m")

        with catch_warnings(AstropyDeprecationWarning) as warning_lines:
            d.units
            assert warning_lines[0].category == AstropyDeprecationWarning

        with catch_warnings(AstropyDeprecationWarning) as warning_lines:
            c = Column([1,2,3], name='a', dtypes="f8", unit="m")
            assert warning_lines[0].category == AstropyDeprecationWarning

        with catch_warnings(AstropyDeprecationWarning) as warning_lines:
            c = Column([1,2,3], name='a', dtype="f8", units="m")
            assert warning_lines[0].category == AstropyDeprecationWarning

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
