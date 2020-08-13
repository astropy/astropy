# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Test behavior related to masked tables"""

import pytest
import numpy as np
import numpy.ma as ma

from astropy.table import Column, MaskedColumn, Table, QTable
from astropy.table.column import BaseColumn
from astropy.time import Time
import astropy.units as u


class SetupData:
    def setup_method(self, method):
        self.a = MaskedColumn(name='a', data=[1, 2, 3], fill_value=1)
        self.b = MaskedColumn(name='b', data=[4, 5, 6], mask=True)
        self.c = MaskedColumn(name='c', data=[7, 8, 9], mask=False)
        self.d_mask = np.array([False, True, False])
        self.d = MaskedColumn(name='d', data=[7, 8, 7], mask=self.d_mask)
        self.t = Table([self.a, self.b], masked=True)
        self.ca = Column(name='ca', data=[1, 2, 3])
        self.sc = MaskedColumn(name='sc', data=[(1, 1.), (2, 2.), (3, 3.)],
                               dtype='i8,f8', fill_value=(0, -1.))


class TestPprint(SetupData):
    def test_pformat(self):
        assert self.t.pformat() == [' a   b ', '--- ---', '  1  --', '  2  --', '  3  --']


class TestFilled:
    """Test the filled method in MaskedColumn and Table"""

    def setup_method(self, method):
        mask = [True, False, False]
        self.meta = {'a': 1, 'b': [2, 3]}
        self.a = MaskedColumn(name='a', data=[1, 2, 3], fill_value=10, mask=mask, meta={'a': 1})
        self.b = MaskedColumn(name='b', data=[4.0, 5.0, 6.0], fill_value=10.0, mask=mask)
        self.c = MaskedColumn(name='c', data=['7', '8', '9'], fill_value='1', mask=mask)

    def test_filled_column(self):
        f = self.a.filled()
        assert np.all(f == [10, 2, 3])
        assert isinstance(f, Column)
        assert not isinstance(f, MaskedColumn)

        # Confirm copy, not ref
        assert f.meta['a'] == 1
        f.meta['a'] = 2
        f[1] = 100
        assert self.a[1] == 2
        assert self.a.meta['a'] == 1

        # Fill with arg fill_value not column fill_value
        f = self.a.filled(20)
        assert np.all(f == [20, 2, 3])

        f = self.b.filled()
        assert np.all(f == [10.0, 5.0, 6.0])
        assert isinstance(f, Column)

        f = self.c.filled()
        assert np.all(f == ['1', '8', '9'])
        assert isinstance(f, Column)

    def test_filled_masked_table(self, tableclass):
        t = tableclass([self.a, self.b, self.c], meta=self.meta)

        f = t.filled()
        assert isinstance(f, Table)
        assert f.masked is False
        assert np.all(f['a'] == [10, 2, 3])
        assert np.allclose(f['b'], [10.0, 5.0, 6.0])
        assert np.all(f['c'] == ['1', '8', '9'])

        # Confirm copy, not ref
        assert f.meta['b'] == [2, 3]
        f.meta['b'][0] = 20
        assert t.meta['b'] == [2, 3]
        f['a'][2] = 100
        assert t['a'][2] == 3

    def test_filled_unmasked_table(self, tableclass):
        t = tableclass([(1, 2), ('3', '4')], names=('a', 'b'), meta=self.meta)
        f = t.filled()
        assert isinstance(f, Table)
        assert f.masked is False
        assert np.all(f['a'] == t['a'])
        assert np.all(f['b'] == t['b'])

        # Confirm copy, not ref
        assert f.meta['b'] == [2, 3]
        f.meta['b'][0] = 20
        assert t.meta['b'] == [2, 3]
        f['a'][1] = 100
        assert t['a'][1] == 2


class TestFillValue(SetupData):
    """Test setting and getting fill value in MaskedColumn and Table"""

    def test_init_set_fill_value(self):
        """Check that setting fill_value in the MaskedColumn init works"""
        assert self.a.fill_value == 1
        c = MaskedColumn(name='c', data=['xxxx', 'yyyy'], fill_value='none')
        assert c.fill_value == 'none'

    def test_set_get_fill_value_for_bare_column(self):
        """Check set and get of fill value works for bare Column"""
        self.d.fill_value = -999
        assert self.d.fill_value == -999
        assert np.all(self.d.filled() == [7, -999, 7])

    def test_set_get_fill_value_for_str_column(self):
        c = MaskedColumn(name='c', data=['xxxx', 'yyyy'], mask=[True, False])
        # assert np.all(c.filled() == ['N/A', 'yyyy'])
        c.fill_value = 'ABCDEF'
        assert c.fill_value == 'ABCD'  # string truncated to dtype length
        assert np.all(c.filled() == ['ABCD', 'yyyy'])
        assert np.all(c.filled('XY') == ['XY', 'yyyy'])

    def test_set_get_fill_value_for_structured_column(self):
        assert self.sc.fill_value == np.array((0, -1.), self.sc.dtype)
        sc = self.sc.copy()
        assert sc.fill_value.item() == (0, -1.)
        sc.fill_value = (-1, np.inf)
        assert sc.fill_value == np.array((-1, np.inf), self.sc.dtype)
        sc2 = MaskedColumn(sc, fill_value=(-2, -np.inf))
        assert sc2.fill_value == np.array((-2, -np.inf), sc2.dtype)

    def test_table_column_mask_not_ref(self):
        """Table column mask is not ref of original column mask"""
        self.b.fill_value = -999
        assert self.t['b'].fill_value != -999

    def test_set_get_fill_value_for_table_column(self):
        """Check set and get of fill value works for Column in a Table"""
        self.t['b'].fill_value = 1
        assert self.t['b'].fill_value == 1
        assert np.all(self.t['b'].filled() == [1, 1, 1])

    def test_data_attribute_fill_and_mask(self):
        """Check that .data attribute preserves fill_value and mask"""
        self.t['b'].fill_value = 1
        self.t['b'].mask = [True, False, True]
        assert self.t['b'].data.fill_value == 1
        assert np.all(self.t['b'].data.mask == [True, False, True])


class TestMaskedColumnInit(SetupData):
    """Initialization of a masked column"""

    def test_set_mask_and_not_ref(self):
        """Check that mask gets set properly and that it is a copy, not ref"""
        assert np.all(~self.a.mask)
        assert np.all(self.b.mask)
        assert np.all(~self.c.mask)
        assert np.all(self.d.mask == self.d_mask)
        self.d.mask[0] = True
        assert not np.all(self.d.mask == self.d_mask)

    def test_set_mask_from_list(self):
        """Set mask from a list"""
        mask_list = [False, True, False]
        a = MaskedColumn(name='a', data=[1, 2, 3], mask=mask_list)
        assert np.all(a.mask == mask_list)

    def test_override_existing_mask(self):
        """Override existing mask values"""
        mask_list = [False, True, False]
        b = MaskedColumn(name='b', data=self.b, mask=mask_list)
        assert np.all(b.mask == mask_list)

    def test_incomplete_mask_spec(self):
        """Incomplete mask specification raises MaskError"""
        mask_list = [False, True]
        with pytest.raises(ma.MaskError):
            MaskedColumn(name='b', length=4, mask=mask_list)


class TestTableInit(SetupData):
    """Initializing a table"""

    @pytest.mark.parametrize('type_str', ('?', 'b', 'i2', 'f4', 'c8', 'S', 'U', 'O'))
    @pytest.mark.parametrize('shape', ((8,), (4, 2), (2, 2, 2)))
    def test_init_from_sequence_data_numeric_typed(self, type_str, shape):
        """Test init from list or list of lists with dtype specified, optionally
        including an np.ma.masked element.
        """
        # Make data of correct dtype and shape, then turn into a list,
        # then use that to init Table with spec'd type_str.
        data = list(range(8))
        np_data = np.array(data, dtype=type_str).reshape(shape)
        np_data_list = np_data.tolist()
        t = Table([np_data_list], dtype=[type_str])
        col = t['col0']
        assert col.dtype == np_data.dtype
        assert np.all(col == np_data)
        assert type(col) is Column

        # Introduce np.ma.masked in the list input and confirm dtype still OK.
        if len(shape) == 1:
            np_data_list[-1] = np.ma.masked
        elif len(shape) == 2:
            np_data_list[-1][-1] = np.ma.masked
        else:
            np_data_list[-1][-1][-1] = np.ma.masked
        last_idx = tuple(-1 for _ in shape)
        t = Table([np_data_list], dtype=[type_str])
        col = t['col0']
        assert col.dtype == np_data.dtype
        assert np.all(col == np_data)
        assert col.mask[last_idx]
        assert type(col) is MaskedColumn

    @pytest.mark.parametrize('type_str', ('?', 'b', 'i2', 'f4', 'c8', 'S', 'U', 'O'))
    @pytest.mark.parametrize('shape', ((8,), (4, 2), (2, 2, 2)))
    def test_init_from_sequence_data_numeric_untyped(self, type_str, shape):
        """Test init from list or list of lists with dtype NOT specified,
        optionally including an np.ma.masked element.
        """
        data = list(range(8))
        np_data = np.array(data, dtype=type_str).reshape(shape)
        np_data_list = np_data.tolist()
        t = Table([np_data_list])
        # Grab the dtype that numpy assigns for the Python list inputs
        dtype_expected = t['col0'].dtype

        # Introduce np.ma.masked in the list input and confirm dtype still OK.
        if len(shape) == 1:
            np_data_list[-1] = np.ma.masked
        elif len(shape) == 2:
            np_data_list[-1][-1] = np.ma.masked
        else:
            np_data_list[-1][-1][-1] = np.ma.masked
        last_idx = tuple(-1 for _ in shape)
        t = Table([np_data_list])
        col = t['col0']

        # Confirm dtype is same as for untype list input w/ no mask
        assert col.dtype == dtype_expected
        assert np.all(col == np_data)
        assert col.mask[last_idx]
        assert type(col) is MaskedColumn

    def test_initialization_with_all_columns(self):
        t1 = Table([self.a, self.b, self.c, self.d, self.ca, self.sc])
        assert t1.colnames == ['a', 'b', 'c', 'd', 'ca', 'sc']
        # Check we get the same result by passing in as list of dict.
        # (Regression test for error uncovered by scintillometry package.)
        lofd = [{k: row[k] for k in t1.colnames} for row in t1]
        t2 = Table(lofd)
        for k in t1.colnames:
            assert t1[k].dtype == t2[k].dtype
            assert np.all(t1[k] == t2[k]) in (True, np.ma.masked)
            assert np.all(getattr(t1[k], 'mask', False)
                          == getattr(t2[k], 'mask', False))

    def test_mask_false_if_input_mask_not_true(self):
        """Masking is always False if initial masked arg is not True"""
        t = Table([self.ca, self.a])
        assert t.masked is False  # True before astropy 4.0
        t = Table([self.ca])
        assert t.masked is False
        t = Table([self.ca, ma.array([1, 2, 3])])
        assert t.masked is False  # True before astropy 4.0

    def test_mask_false_if_no_input_masked(self):
        """Masking not true if not (requested or input requires mask)"""
        t0 = Table([[3, 4]], masked=False)
        t1 = Table(t0, masked=True)
        t2 = Table(t1, masked=False)
        assert not t0.masked
        assert t1.masked
        assert not t2.masked

    def test_mask_property(self):
        t = self.t
        # Access table mask (boolean structured array) by column name
        assert np.all(t.mask['a'] == np.array([False, False, False]))
        assert np.all(t.mask['b'] == np.array([True, True, True]))
        # Check that setting mask from table mask has the desired effect on column
        t.mask['b'] = np.array([False, True, False])
        assert np.all(t['b'].mask == np.array([False, True, False]))
        # Non-masked table returns None for mask attribute
        t2 = Table([self.ca], masked=False)
        assert t2.mask is None
        # Set mask property globally and verify local correctness
        for mask in (True, False):
            t.mask = mask
            for name in ('a', 'b'):
                assert np.all(t[name].mask == mask)


class TestAddColumn:

    def test_add_masked_column_to_masked_table(self):
        t = Table(masked=True)
        assert t.masked
        t.add_column(MaskedColumn(name='a', data=[1, 2, 3], mask=[0, 1, 0]))
        assert t.masked
        t.add_column(MaskedColumn(name='b', data=[4, 5, 6], mask=[1, 0, 1]))
        assert t.masked
        assert isinstance(t['a'], MaskedColumn)
        assert isinstance(t['b'], MaskedColumn)
        assert np.all(t['a'] == np.array([1, 2, 3]))
        assert np.all(t['a'].mask == np.array([0, 1, 0], bool))
        assert np.all(t['b'] == np.array([4, 5, 6]))
        assert np.all(t['b'].mask == np.array([1, 0, 1], bool))

    def test_add_masked_column_to_non_masked_table(self):
        t = Table(masked=False)
        assert not t.masked
        t.add_column(Column(name='a', data=[1, 2, 3]))
        assert not t.masked
        t.add_column(MaskedColumn(name='b', data=[4, 5, 6], mask=[1, 0, 1]))
        assert not t.masked  # Changed in 4.0, table no longer auto-upgrades
        assert isinstance(t['a'], Column)  # Was MaskedColumn before 4.0
        assert isinstance(t['b'], MaskedColumn)
        assert np.all(t['a'] == np.array([1, 2, 3]))
        assert not hasattr(t['a'], 'mask')
        assert np.all(t['b'] == np.array([4, 5, 6]))
        assert np.all(t['b'].mask == np.array([1, 0, 1], bool))

    def test_add_non_masked_column_to_masked_table(self):
        t = Table(masked=True)
        assert t.masked
        t.add_column(Column(name='a', data=[1, 2, 3]))
        assert t.masked
        t.add_column(MaskedColumn(name='b', data=[4, 5, 6], mask=[1, 0, 1]))
        assert t.masked
        assert isinstance(t['a'], MaskedColumn)
        assert isinstance(t['b'], MaskedColumn)
        assert np.all(t['a'] == np.array([1, 2, 3]))
        assert np.all(t['a'].mask == np.array([0, 0, 0], bool))
        assert np.all(t['b'] == np.array([4, 5, 6]))
        assert np.all(t['b'].mask == np.array([1, 0, 1], bool))

    def test_convert_to_masked_table_only_if_necessary(self):
        # Do not convert to masked table, if new column has no masked value.
        # See #1185 for details.
        t = Table(masked=False)
        assert not t.masked
        t.add_column(Column(name='a', data=[1, 2, 3]))
        assert not t.masked
        t.add_column(MaskedColumn(name='b', data=[4, 5, 6], mask=[0, 0, 0]))
        assert not t.masked
        assert np.all(t['a'] == np.array([1, 2, 3]))
        assert np.all(t['b'] == np.array([4, 5, 6]))


class TestRenameColumn:

    def test_rename_masked_column(self):
        t = Table(masked=True)
        t.add_column(MaskedColumn(name='a', data=[1, 2, 3], mask=[0, 1, 0]))
        t['a'].fill_value = 42
        t.rename_column('a', 'b')
        assert t.masked
        assert np.all(t['b'] == np.array([1, 2, 3]))
        assert np.all(t['b'].mask == np.array([0, 1, 0], bool))
        assert t['b'].fill_value == 42
        assert t.colnames == ['b']


class TestRemoveColumn:

    def test_remove_masked_column(self):
        t = Table(masked=True)
        t.add_column(MaskedColumn(name='a', data=[1, 2, 3], mask=[0, 1, 0]))
        t['a'].fill_value = 42
        t.add_column(MaskedColumn(name='b', data=[4, 5, 6], mask=[1, 0, 1]))
        t.remove_column('b')
        assert t.masked
        assert np.all(t['a'] == np.array([1, 2, 3]))
        assert np.all(t['a'].mask == np.array([0, 1, 0], bool))
        assert t['a'].fill_value == 42
        assert t.colnames == ['a']


class TestAddRow:

    def test_add_masked_row_to_masked_table_iterable(self):
        t = Table(masked=True)
        t.add_column(MaskedColumn(name='a', data=[1], mask=[0]))
        t.add_column(MaskedColumn(name='b', data=[4], mask=[1]))
        t.add_row([2, 5], mask=[1, 0])
        t.add_row([3, 6], mask=[0, 1])
        assert t.masked
        assert np.all(np.array(t['a']) == np.array([1, 2, 3]))
        assert np.all(t['a'].mask == np.array([0, 1, 0], bool))
        assert np.all(np.array(t['b']) == np.array([4, 5, 6]))
        assert np.all(t['b'].mask == np.array([1, 0, 1], bool))

    def test_add_masked_row_to_masked_table_mapping1(self):
        t = Table(masked=True)
        t.add_column(MaskedColumn(name='a', data=[1], mask=[0]))
        t.add_column(MaskedColumn(name='b', data=[4], mask=[1]))
        t.add_row({'b': 5, 'a': 2}, mask={'a': 1, 'b': 0})
        t.add_row({'a': 3, 'b': 6}, mask={'b': 1, 'a': 0})
        assert t.masked
        assert np.all(np.array(t['a']) == np.array([1, 2, 3]))
        assert np.all(t['a'].mask == np.array([0, 1, 0], bool))
        assert np.all(np.array(t['b']) == np.array([4, 5, 6]))
        assert np.all(t['b'].mask == np.array([1, 0, 1], bool))

    def test_add_masked_row_to_masked_table_mapping2(self):
        # When adding values to a masked table, if the mask is specified as a
        # dict, then values not specified will have mask values set to True
        t = Table(masked=True)
        t.add_column(MaskedColumn(name='a', data=[1], mask=[0]))
        t.add_column(MaskedColumn(name='b', data=[4], mask=[1]))
        t.add_row({'b': 5}, mask={'b': 0})
        t.add_row({'a': 3}, mask={'a': 0})
        assert t.masked
        assert t['a'][0] == 1 and t['a'][2] == 3
        assert np.all(t['a'].mask == np.array([0, 1, 0], bool))
        assert t['b'][1] == 5
        assert np.all(t['b'].mask == np.array([1, 0, 1], bool))

    def test_add_masked_row_to_masked_table_mapping3(self):
        # When adding values to a masked table, if mask is not passed to
        # add_row, then the mask should be set to False if values are present
        # and True if not.
        t = Table(masked=True)
        t.add_column(MaskedColumn(name='a', data=[1], mask=[0]))
        t.add_column(MaskedColumn(name='b', data=[4], mask=[1]))
        t.add_row({'b': 5})
        t.add_row({'a': 3})
        assert t.masked
        assert t['a'][0] == 1 and t['a'][2] == 3
        assert np.all(t['a'].mask == np.array([0, 1, 0], bool))
        assert t['b'][1] == 5
        assert np.all(t['b'].mask == np.array([1, 0, 1], bool))

    def test_add_masked_row_to_masked_table_mapping4(self):
        # When adding values to a masked table, if the mask is specified as a
        # dict, then keys in values should match keys in mask
        t = Table(masked=True)
        t.add_column(MaskedColumn(name='a', data=[1], mask=[0]))
        t.add_column(MaskedColumn(name='b', data=[4], mask=[1]))
        with pytest.raises(ValueError) as exc:
            t.add_row({'b': 5}, mask={'a': True})
        assert exc.value.args[0] == 'keys in mask should match keys in vals'

    def test_add_masked_row_to_masked_table_mismatch(self):
        t = Table(masked=True)
        t.add_column(MaskedColumn(name='a', data=[1], mask=[0]))
        t.add_column(MaskedColumn(name='b', data=[4], mask=[1]))
        with pytest.raises(TypeError) as exc:
            t.add_row([2, 5], mask={'a': 1, 'b': 0})
        assert exc.value.args[0] == "Mismatch between type of vals and mask"
        with pytest.raises(TypeError) as exc:
            t.add_row({'b': 5, 'a': 2}, mask=[1, 0])
        assert exc.value.args[0] == "Mismatch between type of vals and mask"

    def test_add_masked_row_to_non_masked_table_iterable(self):
        t = Table(masked=False)
        t['a'] = [1]
        t['b'] = [4]
        t['c'] = Time([1], format='cxcsec')

        tm = Time(2, format='cxcsec')
        assert not t.masked
        t.add_row([2, 5, tm])
        assert not t.masked
        t.add_row([3, 6, tm], mask=[0, 1, 1])
        assert not t.masked

        assert type(t['a']) is Column
        assert type(t['b']) is MaskedColumn
        assert type(t['c']) is Time

        assert np.all(t['a'] == [1, 2, 3])
        assert np.all(t['b'].data == [4, 5, 6])
        assert np.all(t['b'].mask == [False, False, True])
        assert np.all(t['c'][:2] == Time([1, 2], format='cxcsec'))
        assert np.all(t['c'].mask == [False, False, True])

    def test_add_row_cannot_mask_column_raises_typeerror(self):
        t = QTable()
        t['a'] = [1, 2] * u.m
        t.add_row((3 * u.m,))  # No problem
        with pytest.raises(ValueError) as exc:
            t.add_row((3 * u.m,), mask=(True,))
        assert (exc.value.args[0].splitlines()
                == ["Unable to insert row because of exception in column 'a':",
                    "mask was supplied for column 'a' but it does not support masked values"])


def test_setting_from_masked_column():
    """Test issue in #2997"""
    mask_b = np.array([True, True, False, False])
    for select in (mask_b, slice(0, 2)):
        t = Table(masked=True)
        t['a'] = Column([1, 2, 3, 4])
        t['b'] = MaskedColumn([11, 22, 33, 44], mask=mask_b)
        t['c'] = MaskedColumn([111, 222, 333, 444], mask=[True, False, True, False])

        t['b'][select] = t['c'][select]
        assert t['b'][1] == t[1]['b']
        assert t['b'][0] is np.ma.masked  # Original state since t['c'][0] is masked
        assert t['b'][1] == 222  # New from t['c'] since t['c'][1] is unmasked
        assert t['b'][2] == 33
        assert t['b'][3] == 44
        assert np.all(t['b'].mask == t.mask['b'])  # Avoid t.mask in general, this is for testing

        mask_before_add = t.mask.copy()
        t['d'] = np.arange(len(t))
        assert np.all(t.mask['b'] == mask_before_add['b'])


def test_coercing_fill_value_type():
    """
    Test that masked column fill_value is coerced into the correct column type.
    """
    # This is the original example posted on the astropy@scipy mailing list
    t = Table({'a': ['1']}, masked=True)
    t['a'].set_fill_value('0')
    t2 = Table(t, names=['a'], dtype=[np.int32])
    assert isinstance(t2['a'].fill_value, np.int32)

    # Unit test the same thing.
    c = MaskedColumn(['1'])
    c.set_fill_value('0')
    c2 = MaskedColumn(c, dtype=np.int32)
    assert isinstance(c2.fill_value, np.int32)


def test_mask_copy():
    """Test that the mask is copied when copying a table (issue #7362)."""

    c = MaskedColumn([1, 2], mask=[False, True])
    c2 = MaskedColumn(c, copy=True)
    c2.mask[0] = True
    assert np.all(c.mask == [False, True])
    assert np.all(c2.mask == [True, True])


def test_masked_as_array_with_mixin():
    """Test that as_array() and Table.mask attr work with masked mixin columns"""
    t = Table()
    t['a'] = Time([1, 2], format='cxcsec')
    t['b'] = [3, 4]
    t['c'] = [5, 6] * u.m

    # With no mask, the output should be ndarray
    ta = t.as_array()
    assert isinstance(ta, np.ndarray) and not isinstance(ta, np.ma.MaskedArray)

    # With a mask, output is MaskedArray
    t['a'][1] = np.ma.masked
    ta = t.as_array()
    assert isinstance(ta, np.ma.MaskedArray)
    assert np.all(ta['a'].mask == [False, True])
    assert np.isclose(ta['a'][0].cxcsec, 1.0)
    assert np.all(ta['b'].mask == False)  # noqa
    assert np.all(ta['c'].mask == False)  # noqa

    # Check table ``mask`` property
    tm = t.mask
    assert np.all(tm['a'] == [False, True])
    assert np.all(tm['b'] == False)  # noqa
    assert np.all(tm['c'] == False)  # noqa


def test_masked_column_with_unit_in_qtable():
    """Test that adding a MaskedColumn with a unit to QTable issues warning"""
    t = QTable()
    t['a'] = MaskedColumn([1, 2])
    assert isinstance(t['a'], MaskedColumn)

    t['b'] = MaskedColumn([1, 2], unit=u.m)
    assert isinstance(t['b'], u.Quantity)

    with pytest.warns(UserWarning, match="dropping mask in Quantity column 'c'") as w:
        t['c'] = MaskedColumn([1, 2], unit=u.m, mask=[True, False])
    assert len(w) == 1
    assert isinstance(t['b'], u.Quantity)


def test_masked_column_data_attribute_is_plain_masked_array():
    c = MaskedColumn([1, 2], mask=[False, True])
    c_data = c.data
    assert type(c_data) is np.ma.MaskedArray
    assert type(c_data.data) is np.ndarray


def test_mask_slicing_count_array_finalize():
    """Check that we don't finalize MaskedColumn too often.

    Regression test for gh-6721.
    """
    # Create a new BaseColumn class that counts how often
    # ``__array_finalize__`` is called.
    class MyBaseColumn(BaseColumn):
        counter = 0

        def __array_finalize__(self, obj):
            super().__array_finalize__(obj)
            MyBaseColumn.counter += 1

    # Base a new MaskedColumn class on it.  The normal MaskedColumn
    # hardcodes the initialization to BaseColumn, so we exchange that.
    class MyMaskedColumn(MaskedColumn, Column, MyBaseColumn):
        def __new__(cls, *args, **kwargs):
            self = super().__new__(cls, *args, **kwargs)
            self._baseclass = MyBaseColumn
            return self

    # Creation really needs 2 finalizations (once for the BaseColumn
    # call inside ``__new__`` and once when the view as a MaskedColumn
    # is taken), but since the first is hardcoded, we do not capture it
    # and thus the count is only 1.
    c = MyMaskedColumn([1, 2], mask=[False, True])
    assert MyBaseColumn.counter == 1
    # slicing should need only one ``__array_finalize__`` (used to be 3).
    c0 = c[:]
    assert MyBaseColumn.counter == 2
    # repr should need none (used to be 2!!)
    repr(c0)
    assert MyBaseColumn.counter == 2
