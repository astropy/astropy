# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Test behavior related to masked tables"""

from distutils import version
import numpy as np
import numpy.ma as ma

from ...tests.helper import pytest
from ...table import Column, MaskedColumn, Table

numpy_lt_1p5 = version.LooseVersion(np.__version__) < version.LooseVersion('1.5')


class SetupData(object):
    def setup_method(self, method):
        self.a = MaskedColumn(name='a', data=[1, 2, 3], fill_value=1)
        self.b = MaskedColumn(name='b', data=[4, 5, 6], mask=True)
        self.c = MaskedColumn(name='c', data=[7, 8, 9], mask=False)
        self.d_mask = np.array([False, True, False])
        self.d = MaskedColumn(name='d', data=[7, 8, 7], mask=self.d_mask)
        self.t = Table([self.a, self.b], masked=True)
        self.ca = Column(name='ca', data=[1, 2, 3])


@pytest.mark.xfail('numpy_lt_1p5')
class TestPprint(SetupData):
    def test_pformat(self):
        assert self.t.pformat() == [' a   b ', '--- ---', '  1  --', '  2  --', '  3  --']


@pytest.mark.xfail('numpy_lt_1p5')
class TestFilled(object):
    """Test the filled method in MaskedColumn and Table"""
    def setup_method(self, method):
        mask = [True, False, False]
        self.meta = {'a': 1, 'b': [2, 3]}
        a = self.a = MaskedColumn(name='a', data=[1, 2, 3], fill_value=10, mask=mask, meta={'a': 1})
        b = self.b = MaskedColumn(name='b', data=[4.0, 5.0, 6.0], fill_value=10.0, mask=mask)
        c = self.c = MaskedColumn(name='c', data=['7', '8', '9'], fill_value='1', mask=mask)

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

    def test_filled_masked_table(self):
        t = Table([self.a, self.b, self.c], meta=self.meta)

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

    def test_filled_unmasked_table(self):
        t = Table([(1, 2), ('3', '4')], names=('a', 'b'), meta=self.meta)
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


@pytest.mark.xfail('numpy_lt_1p5')
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

    def test_table_column_mask_not_ref(self):
        """Table column mask is not ref of original column mask"""
        self.b.fill_value = -999
        assert self.t['b'].fill_value != -999

    def test_set_get_fill_value_for_table_column(self):
        """Check set and get of fill value works for Column in a Table"""
        self.t['b'].fill_value = 1
        assert self.t['b'].fill_value == 1
        assert np.all(self.t['b'].filled() == [1, 1, 1])
        assert self.t._data['b'].fill_value == 1

    def test_data_attribute_fill_and_mask(self):
        """Check that .data attribute preserves fill_value and mask"""
        self.t['b'].fill_value = 1
        self.t['b'].mask = [True, False, True]
        assert self.t['b'].data.fill_value == 1
        assert np.all(self.t['b'].data.mask == [True, False, True])


@pytest.mark.xfail('numpy_lt_1p5')
class TestMaskedColumnInit(SetupData):
    """Initialization of a masked column"""

    def test_set_mask_and_not_ref(self):
        """Check that mask gets set properly and that it is a copy, not ref"""
        assert np.all(self.a.mask == False)
        assert np.all(self.b.mask == True)
        assert np.all(self.c.mask == False)
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
        """Incomplete mask specification (mask values cycle through available)"""
        mask_list = [False, True]
        b = MaskedColumn(name='b', length=4, mask=mask_list)
        assert np.all(b.mask == mask_list + mask_list)


@pytest.mark.xfail('numpy_lt_1p5')
class TestTableInit(SetupData):
    """Initializing a table"""

    def test_mask_true_if_any_input_masked(self):
        """Masking is True if any input is masked"""
        t = Table([self.ca, self.a])
        assert t.masked is True
        t = Table([self.ca])
        assert t.masked is False
        t = Table([self.ca, ma.array([1, 2, 3])])
        assert t.masked is True

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


@pytest.mark.xfail('numpy_lt_1p5')
class TestAddColumn(object):

    def test_add_masked_column_to_masked_table(self):
        t = Table(masked=True)
        assert t.masked
        t.add_column(MaskedColumn(name='a', data=[1,2,3], mask=[0,1,0]))
        assert t.masked
        t.add_column(MaskedColumn(name='b', data=[4,5,6], mask=[1,0,1]))
        assert t.masked
        assert np.all(t['a'] == np.array([1,2,3]))
        assert np.all(t['a'].mask == np.array([0,1,0], bool))
        assert np.all(t['b'] == np.array([4,5,6]))
        assert np.all(t['b'].mask == np.array([1,0,1], bool))

    def test_add_masked_column_to_non_masked_table(self):
        t = Table(masked=False)
        assert not t.masked
        t.add_column(Column(name='a', data=[1,2,3]))
        assert not t.masked
        t.add_column(MaskedColumn(name='b', data=[4,5,6], mask=[1,0,1]))
        assert t.masked
        assert np.all(t['a'] == np.array([1,2,3]))
        assert np.all(t['a'].mask == np.array([0,0,0], bool))
        assert np.all(t['b'] == np.array([4,5,6]))
        assert np.all(t['b'].mask == np.array([1,0,1], bool))

    def test_add_non_masked_column_to_masked_table(self):
        t = Table(masked=True)
        assert t.masked
        t.add_column(Column(name='a', data=[1,2,3]))
        assert t.masked
        t.add_column(MaskedColumn(name='b', data=[4,5,6], mask=[1,0,1]))
        assert t.masked
        assert np.all(t['a'] == np.array([1,2,3]))
        assert np.all(t['a'].mask == np.array([0,0,0], bool))
        assert np.all(t['b'] == np.array([4,5,6]))
        assert np.all(t['b'].mask == np.array([1,0,1], bool))


@pytest.mark.xfail('numpy_lt_1p5')
class TestRenameColumn(object):

    def test_rename_masked_column(self):
        t = Table(masked=True)
        t.add_column(MaskedColumn(name='a', data=[1,2,3], mask=[0,1,0]))
        t['a'].fill_value = 42
        t.rename_column('a', 'b')
        assert t.masked
        assert np.all(t['b'] == np.array([1,2,3]))
        assert np.all(t['b'].mask == np.array([0,1,0], bool))
        assert t['a'].fill_value == 42
        assert t.colnames == ['b']

@pytest.mark.xfail('numpy_lt_1p5')
class TestRemoveColumn(object):

    def test_remove_masked_column(self):
        t = Table(masked=True)
        t.add_column(MaskedColumn(name='a', data=[1,2,3], mask=[0,1,0]))
        t['a'].fill_value = 42
        t.add_column(MaskedColumn(name='b', data=[4,5,6], mask=[1,0,1]))
        t.remove_column('b')
        assert t.masked
        assert np.all(t['a'] == np.array([1,2,3]))
        assert np.all(t['a'].mask == np.array([0,1,0], bool))
        assert t['a'].fill_value == 42
        assert t.colnames == ['a']


@pytest.mark.xfail('numpy_lt_1p5')
class TestAddRow(object):

    def test_add_masked_row_to_masked_table_iterable(self):
        t = Table(masked=True)
        t.add_column(MaskedColumn(name='a', data=[1], mask=[0]))
        t.add_column(MaskedColumn(name='b', data=[4], mask=[1]))
        t.add_row([2,5], mask=[1,0])
        t.add_row([3,6], mask=[0,1])
        assert t.masked
        assert np.all(np.array(t['a']) == np.array([1,2,3]))
        assert np.all(t['a'].mask == np.array([0,1,0], bool))
        assert np.all(np.array(t['b']) == np.array([4,5,6]))
        assert np.all(t['b'].mask == np.array([1,0,1], bool))

    def test_add_masked_row_to_masked_table_mapping1(self):
        t = Table(masked=True)
        t.add_column(MaskedColumn(name='a', data=[1], mask=[0]))
        t.add_column(MaskedColumn(name='b', data=[4], mask=[1]))
        t.add_row({'b':5, 'a':2}, mask={'a':1, 'b':0})
        t.add_row({'a':3, 'b':6}, mask={'b':1, 'a':0})
        assert t.masked
        assert np.all(np.array(t['a']) == np.array([1,2,3]))
        assert np.all(t['a'].mask == np.array([0,1,0], bool))
        assert np.all(np.array(t['b']) == np.array([4,5,6]))
        assert np.all(t['b'].mask == np.array([1,0,1], bool))

    def test_add_masked_row_to_masked_table_mapping2(self):
        # When adding values to a masked table, if the mask is specified as a
        # dict, then values not specified will have mask values set to True
        t = Table(masked=True)
        t.add_column(MaskedColumn(name='a', data=[1], mask=[0]))
        t.add_column(MaskedColumn(name='b', data=[4], mask=[1]))
        t.add_row({'b':5}, mask={'b':0})
        t.add_row({'a':3}, mask={'a':0})
        assert t.masked
        assert t['a'][0] == 1 and t['a'][2] == 3
        assert np.all(t['a'].mask == np.array([0,1,0], bool))
        assert t['b'][1] == 5
        assert np.all(t['b'].mask == np.array([1,0,1], bool))

    def test_add_masked_row_to_masked_table_mapping3(self):
        # When adding values to a masked table, if mask is not passed to
        # add_row, then the mask should be set to False if values are present
        # and True if not.
        t = Table(masked=True)
        t.add_column(MaskedColumn(name='a', data=[1], mask=[0]))
        t.add_column(MaskedColumn(name='b', data=[4], mask=[1]))
        t.add_row({'b':5})
        t.add_row({'a':3})
        assert t.masked
        assert t['a'][0] == 1 and t['a'][2] == 3
        assert np.all(t['a'].mask == np.array([0,1,0], bool))
        assert t['b'][1] == 5
        assert np.all(t['b'].mask == np.array([1,0,1], bool))

    def test_add_masked_row_to_masked_table_mapping4(self):
        # When adding values to a masked table, if the mask is specified as a
        # dict, then keys in values should match keys in mask
        t = Table(masked=True)
        t.add_column(MaskedColumn(name='a', data=[1], mask=[0]))
        t.add_column(MaskedColumn(name='b', data=[4], mask=[1]))
        with pytest.raises(ValueError) as exc:
            t.add_row({'b':5}, mask={'a':True})
        assert exc.value.args[0] == 'keys in mask should match keys in vals'

    def test_add_masked_row_to_masked_table_mismatch(self):
        t = Table(masked=True)
        t.add_column(MaskedColumn(name='a', data=[1], mask=[0]))
        t.add_column(MaskedColumn(name='b', data=[4], mask=[1]))
        with pytest.raises(TypeError) as exc:
            t.add_row([2,5], mask={'a':1, 'b':0})
        assert exc.value.args[0] == "Mismatch between type of vals and mask"
        with pytest.raises(TypeError) as exc:
            t.add_row({'b':5, 'a':2}, mask=[1,0])
        assert exc.value.args[0] == "Mismatch between type of vals and mask"

    def test_add_masked_row_to_non_masked_table_iterable(self):
        t = Table(masked=False)
        t.add_column(Column(name='a', data=[1]))
        t.add_column(Column(name='b', data=[4]))
        assert not t.masked
        t.add_row([2,5])
        assert not t.masked
        t.add_row([3,6], mask=[0,1])
        assert t.masked
        assert np.all(np.array(t['a']) == np.array([1,2,3]))
        assert np.all(t['a'].mask == np.array([0,0,0], bool))
        assert np.all(np.array(t['b']) == np.array([4,5,6]))
        assert np.all(t['b'].mask == np.array([0,0,1], bool))
