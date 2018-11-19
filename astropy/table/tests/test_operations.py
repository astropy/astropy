# Licensed under a 3-clause BSD style license - see LICENSE.rst


from collections import OrderedDict

import pytest
import numpy as np

from ...tests.helper import catch_warnings
from ...table import Table, QTable, TableMergeError
from ...table.operations import _get_out_class
from ... import units as u
from ...utils import metadata
from ...utils.metadata import MergeConflictError
from ... import table
from ...time import Time
from ...coordinates import SkyCoord


def sort_eq(list1, list2):
    return sorted(list1) == sorted(list2)


def skycoord_equal(sc1, sc2):
    if not sc1.is_equivalent_frame(sc2):
        return False
    if sc1.representation_type is not sc2.representation_type:
        return False
    if sc1.shape != sc2.shape:
        return False  # Maybe raise ValueError corresponding to future numpy behavior
    eq = np.ones(shape=sc1.shape, dtype=bool)
    for comp in sc1.data.components:
        eq &= getattr(sc1.data, comp) == getattr(sc2.data, comp)
    return np.all(eq)


class TestJoin():

    def _setup(self, t_cls=Table):
        lines1 = [' a   b   c ',
                  '  0 foo  L1',
                  '  1 foo  L2',
                  '  1 bar  L3',
                  '  2 bar  L4']
        lines2 = [' a   b   d ',
                  '  1 foo  R1',
                  '  1 foo  R2',
                  '  2 bar  R3',
                  '  4 bar  R4']
        self.t1 = t_cls.read(lines1, format='ascii')
        self.t2 = t_cls.read(lines2, format='ascii')
        self.t3 = t_cls(self.t2, copy=True)

        self.t1.meta.update(OrderedDict([('b', [1, 2]), ('c', {'a': 1}), ('d', 1)]))
        self.t2.meta.update(OrderedDict([('b', [3, 4]), ('c', {'b': 1}), ('a', 1)]))
        self.t3.meta.update(OrderedDict([('b', 3), ('c', [1, 2]), ('d', 2), ('a', 1)]))

        self.meta_merge = OrderedDict([('b', [1, 2, 3, 4]),
                                       ('c', {'a': 1, 'b': 1}),
                                       ('d', 1),
                                       ('a', 1)])

    def test_table_meta_merge(self, operation_table_type):
        self._setup(operation_table_type)
        out = table.join(self.t1, self.t2, join_type='inner')
        assert out.meta == self.meta_merge

    def test_table_meta_merge_conflict(self, operation_table_type):
        self._setup(operation_table_type)

        with catch_warnings() as w:
            out = table.join(self.t1, self.t3, join_type='inner')
        assert len(w) == 3

        assert out.meta == self.t3.meta

        with catch_warnings() as w:
            out = table.join(self.t1, self.t3, join_type='inner', metadata_conflicts='warn')
        assert len(w) == 3

        assert out.meta == self.t3.meta

        with catch_warnings() as w:
            out = table.join(self.t1, self.t3, join_type='inner', metadata_conflicts='silent')
        assert len(w) == 0

        assert out.meta == self.t3.meta

        with pytest.raises(MergeConflictError):
            out = table.join(self.t1, self.t3, join_type='inner', metadata_conflicts='error')

        with pytest.raises(ValueError):
            out = table.join(self.t1, self.t3, join_type='inner', metadata_conflicts='nonsense')

    def test_both_unmasked_inner(self, operation_table_type):
        self._setup(operation_table_type)
        t1 = self.t1
        t2 = self.t2

        # Basic join with default parameters (inner join on common keys)
        t12 = table.join(t1, t2)
        assert type(t12) is operation_table_type
        assert type(t12['a']) is type(t1['a'])
        assert type(t12['b']) is type(t1['b'])
        assert type(t12['c']) is type(t1['c'])
        assert type(t12['d']) is type(t2['d'])
        assert t12.masked is False
        assert sort_eq(t12.pformat(), [' a   b   c   d ',
                                       '--- --- --- ---',
                                       '  1 foo  L2  R1',
                                       '  1 foo  L2  R2',
                                       '  2 bar  L4  R3'])
        # Table meta merged properly
        assert t12.meta == self.meta_merge

    def test_both_unmasked_left_right_outer(self, operation_table_type):
        if operation_table_type is QTable:
            pytest.xfail('Quantity columns do not support masking.')
        self._setup(operation_table_type)
        t1 = self.t1
        t2 = self.t2

        # Left join
        t12 = table.join(t1, t2, join_type='left')
        assert t12.masked is True
        assert sort_eq(t12.pformat(), [' a   b   c   d ',
                                       '--- --- --- ---',
                                       '  0 foo  L1  --',
                                       '  1 bar  L3  --',
                                       '  1 foo  L2  R1',
                                       '  1 foo  L2  R2',
                                       '  2 bar  L4  R3'])

        # Right join
        t12 = table.join(t1, t2, join_type='right')
        assert t12.masked is True
        assert sort_eq(t12.pformat(), [' a   b   c   d ',
                                       '--- --- --- ---',
                                       '  1 foo  L2  R1',
                                       '  1 foo  L2  R2',
                                       '  2 bar  L4  R3',
                                       '  4 bar  --  R4'])

        # Outer join
        t12 = table.join(t1, t2, join_type='outer')
        assert t12.masked is True
        assert sort_eq(t12.pformat(), [' a   b   c   d ',
                                       '--- --- --- ---',
                                       '  0 foo  L1  --',
                                       '  1 bar  L3  --',
                                       '  1 foo  L2  R1',
                                       '  1 foo  L2  R2',
                                       '  2 bar  L4  R3',
                                       '  4 bar  --  R4'])

        # Check that the common keys are 'a', 'b'
        t12a = table.join(t1, t2, join_type='outer')
        t12b = table.join(t1, t2, join_type='outer', keys=['a', 'b'])
        assert np.all(t12a.as_array() == t12b.as_array())

    def test_both_unmasked_single_key_inner(self, operation_table_type):
        self._setup(operation_table_type)
        t1 = self.t1
        t2 = self.t2

        # Inner join on 'a' column
        t12 = table.join(t1, t2, keys='a')
        assert type(t12) is operation_table_type
        assert type(t12['a']) is type(t1['a'])
        assert type(t12['b_1']) is type(t1['b'])
        assert type(t12['c']) is type(t1['c'])
        assert type(t12['b_2']) is type(t2['b'])
        assert type(t12['d']) is type(t2['d'])
        assert t12.masked is False
        assert sort_eq(t12.pformat(), [' a  b_1  c  b_2  d ',
                                       '--- --- --- --- ---',
                                       '  1 foo  L2 foo  R1',
                                       '  1 foo  L2 foo  R2',
                                       '  1 bar  L3 foo  R1',
                                       '  1 bar  L3 foo  R2',
                                       '  2 bar  L4 bar  R3'])

    def test_both_unmasked_single_key_left_right_outer(self, operation_table_type):
        if operation_table_type is QTable:
            pytest.xfail('Quantity columns do not support masking.')
        self._setup(operation_table_type)
        t1 = self.t1
        t2 = self.t2

        # Left join
        t12 = table.join(t1, t2, join_type='left', keys='a')
        assert t12.masked is True
        assert sort_eq(t12.pformat(), [' a  b_1  c  b_2  d ',
                                       '--- --- --- --- ---',
                                       '  0 foo  L1  --  --',
                                       '  1 foo  L2 foo  R1',
                                       '  1 foo  L2 foo  R2',
                                       '  1 bar  L3 foo  R1',
                                       '  1 bar  L3 foo  R2',
                                       '  2 bar  L4 bar  R3'])

        # Right join
        t12 = table.join(t1, t2, join_type='right', keys='a')
        assert t12.masked is True
        assert sort_eq(t12.pformat(), [' a  b_1  c  b_2  d ',
                                       '--- --- --- --- ---',
                                       '  1 foo  L2 foo  R1',
                                       '  1 foo  L2 foo  R2',
                                       '  1 bar  L3 foo  R1',
                                       '  1 bar  L3 foo  R2',
                                       '  2 bar  L4 bar  R3',
                                       '  4  --  -- bar  R4'])

        # Outer join
        t12 = table.join(t1, t2, join_type='outer', keys='a')
        assert t12.masked is True
        assert sort_eq(t12.pformat(), [' a  b_1  c  b_2  d ',
                                       '--- --- --- --- ---',
                                       '  0 foo  L1  --  --',
                                       '  1 foo  L2 foo  R1',
                                       '  1 foo  L2 foo  R2',
                                       '  1 bar  L3 foo  R1',
                                       '  1 bar  L3 foo  R2',
                                       '  2 bar  L4 bar  R3',
                                       '  4  --  -- bar  R4'])

    def test_masked_unmasked(self, operation_table_type):
        if operation_table_type is QTable:
            pytest.xfail('Quantity columns do not support masking.')
        self._setup(operation_table_type)
        t1 = self.t1
        t1m = operation_table_type(self.t1, masked=True)
        t2 = self.t2

        # Result should be masked even though not req'd by inner join
        t1m2 = table.join(t1m, t2, join_type='inner')
        assert t1m2.masked is True

        # Result should match non-masked result
        t12 = table.join(t1, t2)
        assert np.all(t12.as_array() == np.array(t1m2))

        # Mask out some values in left table and make sure they propagate
        t1m['b'].mask[1] = True
        t1m['c'].mask[2] = True
        t1m2 = table.join(t1m, t2, join_type='inner', keys='a')
        assert sort_eq(t1m2.pformat(), [' a  b_1  c  b_2  d ',
                                        '--- --- --- --- ---',
                                        '  1  --  L2 foo  R1',
                                        '  1  --  L2 foo  R2',
                                        '  1 bar  -- foo  R1',
                                        '  1 bar  -- foo  R2',
                                        '  2 bar  L4 bar  R3'])

        t21m = table.join(t2, t1m, join_type='inner', keys='a')
        assert sort_eq(t21m.pformat(), [' a  b_1  d  b_2  c ',
                                        '--- --- --- --- ---',
                                        '  1 foo  R2  --  L2',
                                        '  1 foo  R2 bar  --',
                                        '  1 foo  R1  --  L2',
                                        '  1 foo  R1 bar  --',
                                        '  2 bar  R3 bar  L4'])

    def test_masked_masked(self, operation_table_type):
        self._setup(operation_table_type)
        """Two masked tables"""
        if operation_table_type is QTable:
            pytest.xfail('Quantity columns do not support masking.')
        t1 = self.t1
        t1m = operation_table_type(self.t1, masked=True)
        t2 = self.t2
        t2m = operation_table_type(self.t2, masked=True)

        # Result should be masked even though not req'd by inner join
        t1m2m = table.join(t1m, t2m, join_type='inner')
        assert t1m2m.masked is True

        # Result should match non-masked result
        t12 = table.join(t1, t2)
        assert np.all(t12.as_array() == np.array(t1m2m))

        # Mask out some values in both tables and make sure they propagate
        t1m['b'].mask[1] = True
        t1m['c'].mask[2] = True
        t2m['d'].mask[2] = True
        t1m2m = table.join(t1m, t2m, join_type='inner', keys='a')
        assert sort_eq(t1m2m.pformat(), [' a  b_1  c  b_2  d ',
                                         '--- --- --- --- ---',
                                         '  1  --  L2 foo  R1',
                                         '  1  --  L2 foo  R2',
                                         '  1 bar  -- foo  R1',
                                         '  1 bar  -- foo  R2',
                                         '  2 bar  L4 bar  --'])

    def test_col_rename(self, operation_table_type):
        self._setup(operation_table_type)
        """
        Test auto col renaming when there is a conflict.  Use
        non-default values of uniq_col_name and table_names.
        """
        t1 = self.t1
        t2 = self.t2
        t12 = table.join(t1, t2, uniq_col_name='x_{table_name}_{col_name}_y',
                         table_names=['L', 'R'], keys='a')
        assert t12.colnames == ['a', 'x_L_b_y', 'c', 'x_R_b_y', 'd']

    def test_rename_conflict(self, operation_table_type):
        self._setup(operation_table_type)
        """
        Test that auto-column rename fails because of a conflict
        with an existing column
        """
        t1 = self.t1
        t2 = self.t2
        t1['b_1'] = 1  # Add a new column b_1 that will conflict with auto-rename
        with pytest.raises(TableMergeError):
            table.join(t1, t2, keys='a')

    def test_missing_keys(self, operation_table_type):
        self._setup(operation_table_type)
        """Merge on a key column that doesn't exist"""
        t1 = self.t1
        t2 = self.t2
        with pytest.raises(TableMergeError):
            table.join(t1, t2, keys=['a', 'not there'])

    def test_bad_join_type(self, operation_table_type):
        self._setup(operation_table_type)
        """Bad join_type input"""
        t1 = self.t1
        t2 = self.t2
        with pytest.raises(ValueError):
            table.join(t1, t2, join_type='illegal value')

    def test_no_common_keys(self, operation_table_type):
        self._setup(operation_table_type)
        """Merge tables with no common keys"""
        t1 = self.t1
        t2 = self.t2
        del t1['a']
        del t1['b']
        del t2['a']
        del t2['b']
        with pytest.raises(TableMergeError):
            table.join(t1, t2)

    def test_masked_key_column(self, operation_table_type):
        self._setup(operation_table_type)
        """Merge on a key column that has a masked element"""
        if operation_table_type is QTable:
            pytest.xfail('Quantity columns do not support masking.')
        t1 = self.t1
        t2 = operation_table_type(self.t2, masked=True)
        table.join(t1, t2)  # OK
        t2['a'].mask[0] = True
        with pytest.raises(TableMergeError):
            table.join(t1, t2)

    def test_col_meta_merge(self, operation_table_type):
        self._setup(operation_table_type)
        t1 = self.t1
        t2 = self.t2
        t2.rename_column('d', 'c')  # force col conflict and renaming
        meta1 = OrderedDict([('b', [1, 2]), ('c', {'a': 1}), ('d', 1)])
        meta2 = OrderedDict([('b', [3, 4]), ('c', {'b': 1}), ('a', 1)])

        # Key col 'a', should first value ('cm')
        t1['a'].unit = 'cm'
        t2['a'].unit = 'm'
        # Key col 'b', take first value 't1_b'
        t1['b'].info.description = 't1_b'
        # Key col 'b', take first non-empty value 't1_b'
        t2['b'].info.format = '%6s'
        # Key col 'a', should be merged meta
        t1['a'].info.meta = meta1
        t2['a'].info.meta = meta2
        # Key col 'b', should be meta2
        t2['b'].info.meta = meta2

        # All these should pass through
        t1['c'].info.format = '%3s'
        t1['c'].info.description = 't1_c'

        t2['c'].info.format = '%6s'
        t2['c'].info.description = 't2_c'

        with catch_warnings(metadata.MergeConflictWarning) as warning_lines:
            t12 = table.join(t1, t2, keys=['a', 'b'])

        if operation_table_type is Table:
            assert warning_lines[0].category == metadata.MergeConflictWarning
            assert ("In merged column 'a' the 'unit' attribute does not match (cm != m)"
                    in str(warning_lines[0].message))
        else:
            assert len(warning_lines) == 0

        assert t12['a'].unit == 'm'
        assert t12['b'].info.description == 't1_b'
        assert t12['b'].info.format == '%6s'
        assert t12['a'].info.meta == self.meta_merge
        assert t12['b'].info.meta == meta2
        assert t12['c_1'].info.format == '%3s'
        assert t12['c_1'].info.description == 't1_c'
        assert t12['c_2'].info.format == '%6s'
        assert t12['c_2'].info.description == 't2_c'

    def test_join_multidimensional(self, operation_table_type):
        self._setup(operation_table_type)

        # Regression test for #2984, which was an issue where join did not work
        # on multi-dimensional columns.

        t1 = operation_table_type()
        t1['a'] = [1, 2, 3]
        t1['b'] = np.ones((3, 4))

        t2 = operation_table_type()
        t2['a'] = [1, 2, 3]
        t2['c'] = [4, 5, 6]

        t3 = table.join(t1, t2)

        np.testing.assert_allclose(t3['a'], t1['a'])
        np.testing.assert_allclose(t3['b'], t1['b'])
        np.testing.assert_allclose(t3['c'], t2['c'])

    def test_join_multidimensional_masked(self, operation_table_type):
        self._setup(operation_table_type)
        """
        Test for outer join with multidimensional columns where masking is required.
        (Issue #4059).
        """
        if operation_table_type is QTable:
            pytest.xfail('Quantity columns do not support masking.')

        a = table.MaskedColumn([1, 2, 3], name='a')
        a2 = table.Column([1, 3, 4], name='a')
        b = table.MaskedColumn([[1, 2],
                                [3, 4],
                                [5, 6]],
                               name='b',
                               mask=[[1, 0],
                                     [0, 1],
                                     [0, 0]])
        c = table.Column([[1, 1],
                          [2, 2],
                          [3, 3]],
                         name='c')
        t1 = operation_table_type([a, b])
        t2 = operation_table_type([a2, c])
        t12 = table.join(t1, t2, join_type='inner')

        assert np.all(t12['b'].mask == [[True, False],
                                        [False, False]])
        assert np.all(t12['c'].mask == [[False, False],
                                        [False, False]])

        t12 = table.join(t1, t2, join_type='outer')
        assert np.all(t12['b'].mask == [[True, False],
                                        [False, True],
                                        [False, False],
                                        [True, True]])
        assert np.all(t12['c'].mask == [[False, False],
                                        [True, True],
                                        [False, False],
                                        [False, False]])

    def test_mixin_functionality(self, mixin_cols):
        col = mixin_cols['m']
        cls_name = type(col).__name__
        len_col = len(col)
        idx = np.arange(len_col)
        t1 = table.QTable([idx, col], names=['idx', 'm1'])
        t2 = table.QTable([idx, col], names=['idx', 'm2'])
        # Set up join mismatches for different join_type cases
        t1 = t1[[0, 1, 3]]
        t2 = t2[[0, 2, 3]]

        # Test inner join, which works for all mixin_cols
        out = table.join(t1, t2, join_type='inner')
        assert len(out) == 2
        assert out['m2'].__class__ is col.__class__
        assert np.all(out['idx'] == [0, 3])
        if cls_name == 'SkyCoord':
            # SkyCoord doesn't support __eq__ so use our own
            assert skycoord_equal(out['m1'], col[[0, 3]])
            assert skycoord_equal(out['m2'], col[[0, 3]])
        else:
            assert np.all(out['m1'] == col[[0, 3]])
            assert np.all(out['m2'] == col[[0, 3]])

        # Check for left, right, outer join which requires masking.  Only Time
        # supports this currently.
        if cls_name == 'Time':
            out = table.join(t1, t2, join_type='left')
            assert len(out) == 3
            assert np.all(out['idx'] == [0, 1, 3])
            assert np.all(out['m1'] == t1['m1'])
            assert np.all(out['m2'] == t2['m2'])
            assert np.all(out['m1'].mask == [False, False, False])
            assert np.all(out['m2'].mask == [False, True, False])

            out = table.join(t1, t2, join_type='right')
            assert len(out) == 3
            assert np.all(out['idx'] == [0, 2, 3])
            assert np.all(out['m1'] == t1['m1'])
            assert np.all(out['m2'] == t2['m2'])
            assert np.all(out['m1'].mask == [False, True, False])
            assert np.all(out['m2'].mask == [False, False, False])

            out = table.join(t1, t2, join_type='outer')
            assert len(out) == 4
            assert np.all(out['idx'] == [0, 1, 2, 3])
            assert np.all(out['m1'] == col)
            assert np.all(out['m2'] == col)
            assert np.all(out['m1'].mask == [False, False, True, False])
            assert np.all(out['m2'].mask == [False, True, False, False])
        else:
            # Otherwise make sure it fails with the right exception message
            for join_type in ('outer', 'left', 'right'):
                with pytest.raises(NotImplementedError) as err:
                    table.join(t1, t2, join_type='outer')
                assert ('join requires masking' in str(err) or
                        'join unavailable' in str(err))


class TestSetdiff():

    def _setup(self, t_cls=Table):
        lines1 = [' a   b ',
                 '  0 foo ',
                 '  1 foo ',
                 '  1 bar ',
                 '  2 bar ']
        lines2 = [' a   b ',
                  '  0 foo ',
                  '  3 foo ',
                  '  4 bar ',
                  '  2 bar ']
        lines3 = [' a   b   d ',
                  '  0 foo  R1',
                  '  8 foo  R2',
                  '  1 bar  R3',
                  '  4 bar  R4']
        self.t1 = t_cls.read(lines1, format='ascii')
        self.t2 = t_cls.read(lines2, format='ascii')
        self.t3 = t_cls.read(lines3, format='ascii')

    def test_default_same_columns(self, operation_table_type):
        self._setup(operation_table_type)
        out = table.setdiff(self.t1, self.t2)
        assert type(out['a']) is type(self.t1['a'])
        assert type(out['b']) is type(self.t1['b'])
        assert out.pformat() == [' a   b ',
                                 '--- ---',
                                 '  1 bar',
                                 '  1 foo']

    def test_default_same_tables(self, operation_table_type):
        self._setup(operation_table_type)
        out = table.setdiff(self.t1, self.t1)

        assert type(out['a']) is type(self.t1['a'])
        assert type(out['b']) is type(self.t1['b'])
        assert out.pformat() == [' a   b ',
                                 '--- ---']

    def test_extra_col_left_table(self, operation_table_type):
        self._setup(operation_table_type)

        with pytest.raises(ValueError):
            out = table.setdiff(self.t3, self.t1)

    def test_extra_col_right_table(self, operation_table_type):
        self._setup(operation_table_type)
        out = table.setdiff(self.t1, self.t3)

        assert type(out['a']) is type(self.t1['a'])
        assert type(out['b']) is type(self.t1['b'])
        assert out.pformat() == [' a   b ',
                                 '--- ---',
                                 '  1 foo',
                                 '  2 bar']

    def test_keys(self, operation_table_type):
        self._setup(operation_table_type)
        out = table.setdiff(self.t3, self.t1, keys=['a', 'b'])

        assert type(out['a']) is type(self.t1['a'])
        assert type(out['b']) is type(self.t1['b'])
        assert out.pformat() == [' a   b   d ',
                                 '--- --- ---',
                                 '  4 bar  R4',
                                 '  8 foo  R2']

    def test_missing_key(self, operation_table_type):
        self._setup(operation_table_type)

        with pytest.raises(ValueError):
            out = table.setdiff(self.t3, self.t1, keys=['a', 'd'])


class TestVStack():

    def _setup(self, t_cls=Table):
        self.t1 = t_cls.read([' a   b',
                              ' 0. foo',
                              ' 1. bar'], format='ascii')

        self.t2 = t_cls.read([' a    b   c',
                              ' 2.  pez  4',
                              ' 3.  sez  5'], format='ascii')

        self.t3 = t_cls.read([' a    b',
                              ' 4.   7',
                              ' 5.   8',
                              ' 6.   9'], format='ascii')
        self.t4 = t_cls(self.t1, copy=True, masked=t_cls is Table)

        # The following table has meta-data that conflicts with t1
        self.t5 = t_cls(self.t1, copy=True)

        self.t1.meta.update(OrderedDict([('b', [1, 2]), ('c', {'a': 1}), ('d', 1)]))
        self.t2.meta.update(OrderedDict([('b', [3, 4]), ('c', {'b': 1}), ('a', 1)]))
        self.t4.meta.update(OrderedDict([('b', [5, 6]), ('c', {'c': 1}), ('e', 1)]))
        self.t5.meta.update(OrderedDict([('b', 3), ('c', 'k'), ('d', 1)]))
        self.meta_merge = OrderedDict([('b', [1, 2, 3, 4, 5, 6]),
                                       ('c', {'a': 1, 'b': 1, 'c': 1}),
                                       ('d', 1),
                                       ('a', 1),
                                       ('e', 1)])

    def test_stack_rows(self, operation_table_type):
        self._setup(operation_table_type)
        t2 = self.t1.copy()
        t2.meta.clear()
        out = table.vstack([self.t1, t2[1]])
        assert type(out['a']) is type(self.t1['a'])
        assert type(out['b']) is type(self.t1['b'])
        assert out.pformat() == [' a   b ',
                                 '--- ---',
                                 '0.0 foo',
                                 '1.0 bar',
                                 '1.0 bar']

    def test_stack_table_column(self, operation_table_type):
        self._setup(operation_table_type)
        t2 = self.t1.copy()
        t2.meta.clear()
        out = table.vstack([self.t1, t2['a']])
        assert out.pformat() == [' a   b ',
                                 '--- ---',
                                 '0.0 foo',
                                 '1.0 bar',
                                 '0.0  --',
                                 '1.0  --']

    def test_table_meta_merge(self, operation_table_type):
        self._setup(operation_table_type)
        out = table.vstack([self.t1, self.t2, self.t4], join_type='inner')
        assert out.meta == self.meta_merge

    def test_table_meta_merge_conflict(self, operation_table_type):
        self._setup(operation_table_type)

        with catch_warnings() as w:
            out = table.vstack([self.t1, self.t5], join_type='inner')
        assert len(w) == 2

        assert out.meta == self.t5.meta

        with catch_warnings() as w:
            out = table.vstack([self.t1, self.t5], join_type='inner', metadata_conflicts='warn')
        assert len(w) == 2

        assert out.meta == self.t5.meta

        with catch_warnings() as w:
            out = table.vstack([self.t1, self.t5], join_type='inner', metadata_conflicts='silent')
        assert len(w) == 0

        assert out.meta == self.t5.meta

        with pytest.raises(MergeConflictError):
            out = table.vstack([self.t1, self.t5], join_type='inner', metadata_conflicts='error')

        with pytest.raises(ValueError):
            out = table.vstack([self.t1, self.t5], join_type='inner', metadata_conflicts='nonsense')

    def test_bad_input_type(self, operation_table_type):
        self._setup(operation_table_type)
        with pytest.raises(ValueError):
            table.vstack([])
        with pytest.raises(TypeError):
            table.vstack(1)
        with pytest.raises(TypeError):
            table.vstack([self.t2, 1])
        with pytest.raises(ValueError):
            table.vstack([self.t1, self.t2], join_type='invalid join type')

    def test_stack_basic_inner(self, operation_table_type):
        self._setup(operation_table_type)
        t1 = self.t1
        t2 = self.t2
        t4 = self.t4

        t12 = table.vstack([t1, t2], join_type='inner')
        assert t12.masked is False
        assert type(t12) is operation_table_type
        assert type(t12['a']) is type(t1['a'])
        assert type(t12['b']) is type(t1['b'])
        assert t12.pformat() == [' a   b ',
                                 '--- ---',
                                 '0.0 foo',
                                 '1.0 bar',
                                 '2.0 pez',
                                 '3.0 sez']

        t124 = table.vstack([t1, t2, t4], join_type='inner')
        assert type(t124) is operation_table_type
        assert type(t12['a']) is type(t1['a'])
        assert type(t12['b']) is type(t1['b'])
        assert t124.pformat() == [' a   b ',
                                  '--- ---',
                                  '0.0 foo',
                                  '1.0 bar',
                                  '2.0 pez',
                                  '3.0 sez',
                                  '0.0 foo',
                                  '1.0 bar']

    def test_stack_basic_outer(self, operation_table_type):
        if operation_table_type is QTable:
            pytest.xfail('Quantity columns do not support masking.')
        self._setup(operation_table_type)
        t1 = self.t1
        t2 = self.t2
        t4 = self.t4
        t12 = table.vstack([t1, t2], join_type='outer')
        assert t12.pformat() == [' a   b   c ',
                                 '--- --- ---',
                                 '0.0 foo  --',
                                 '1.0 bar  --',
                                 '2.0 pez   4',
                                 '3.0 sez   5']

        t124 = table.vstack([t1, t2, t4], join_type='outer')
        assert t124.pformat() == [' a   b   c ',
                                  '--- --- ---',
                                  '0.0 foo  --',
                                  '1.0 bar  --',
                                  '2.0 pez   4',
                                  '3.0 sez   5',
                                  '0.0 foo  --',
                                  '1.0 bar  --']

    def test_stack_incompatible(self, operation_table_type):
        self._setup(operation_table_type)
        with pytest.raises(TableMergeError) as excinfo:
            table.vstack([self.t1, self.t3], join_type='inner')
        assert ("The 'b' columns have incompatible types: {0}"
                .format([self.t1['b'].dtype.name, self.t3['b'].dtype.name])
                in str(excinfo))

        with pytest.raises(TableMergeError) as excinfo:
            table.vstack([self.t1, self.t3], join_type='outer')
        assert "The 'b' columns have incompatible types:" in str(excinfo)

        with pytest.raises(TableMergeError):
            table.vstack([self.t1, self.t2], join_type='exact')

        t1_reshape = self.t1.copy()
        t1_reshape['b'].shape = [2, 1]
        with pytest.raises(TableMergeError) as excinfo:
            table.vstack([self.t1, t1_reshape])
        assert "have different shape" in str(excinfo)

    def test_vstack_one_masked(self, operation_table_type):
        if operation_table_type is QTable:
            pytest.xfail('Quantity columns do not support masking.')
        self._setup(operation_table_type)
        t1 = self.t1
        t4 = self.t4
        t4['b'].mask[1] = True
        assert table.vstack([t1, t4]).pformat() == [' a   b ',
                                                    '--- ---',
                                                    '0.0 foo',
                                                    '1.0 bar',
                                                    '0.0 foo',
                                                    '1.0  --']

    def test_col_meta_merge_inner(self, operation_table_type):
        self._setup(operation_table_type)
        t1 = self.t1
        t2 = self.t2
        t4 = self.t4

        # Key col 'a', should last value ('km')
        t1['a'].info.unit = 'cm'
        t2['a'].info.unit = 'm'
        t4['a'].info.unit = 'km'

        # Key col 'a' format should take last when all match
        t1['a'].info.format = '%f'
        t2['a'].info.format = '%f'
        t4['a'].info.format = '%f'

        # Key col 'b', take first value 't1_b'
        t1['b'].info.description = 't1_b'

        # Key col 'b', take first non-empty value '%6s'
        t4['b'].info.format = '%6s'

        # Key col 'a', should be merged meta
        t1['a'].info.meta.update(OrderedDict([('b', [1, 2]), ('c', {'a': 1}), ('d', 1)]))
        t2['a'].info.meta.update(OrderedDict([('b', [3, 4]), ('c', {'b': 1}), ('a', 1)]))
        t4['a'].info.meta.update(OrderedDict([('b', [5, 6]), ('c', {'c': 1}), ('e', 1)]))

        # Key col 'b', should be meta2
        t2['b'].info.meta.update(OrderedDict([('b', [3, 4]), ('c', {'b': 1}), ('a', 1)]))

        with catch_warnings(metadata.MergeConflictWarning) as warning_lines:
            out = table.vstack([t1, t2, t4], join_type='inner')

        if operation_table_type is Table:
            assert warning_lines[0].category == metadata.MergeConflictWarning
            assert ("In merged column 'a' the 'unit' attribute does not match (cm != m)"
                    in str(warning_lines[0].message))
            assert warning_lines[1].category == metadata.MergeConflictWarning
            assert ("In merged column 'a' the 'unit' attribute does not match (m != km)"
                    in str(warning_lines[1].message))
            # Check units are suitably ignored for a regular Table
            assert out.pformat() == ['   a       b   ',
                                     '   km          ',
                                     '-------- ------',
                                     '0.000000    foo',
                                     '1.000000    bar',
                                     '2.000000    pez',
                                     '3.000000    sez',
                                     '0.000000    foo',
                                     '1.000000    bar']
        else:
            assert len(warning_lines) == 0
            # Check QTable correctly dealt with units.
            assert out.pformat() == ['   a       b   ',
                                     '   km          ',
                                     '-------- ------',
                                     '0.000000    foo',
                                     '0.000010    bar',
                                     '0.002000    pez',
                                     '0.003000    sez',
                                     '0.000000    foo',
                                     '1.000000    bar']
        assert out['a'].info.unit == 'km'
        assert out['a'].info.format == '%f'
        assert out['b'].info.description == 't1_b'
        assert out['b'].info.format == '%6s'
        assert out['a'].info.meta == self.meta_merge
        assert out['b'].info.meta == OrderedDict([('b', [3, 4]), ('c', {'b': 1}), ('a', 1)])

    def test_col_meta_merge_outer(self, operation_table_type):
        if operation_table_type is QTable:
            pytest.xfail('Quantity columns do not support masking.')
        self._setup(operation_table_type)
        t1 = self.t1
        t2 = self.t2
        t4 = self.t4

        # Key col 'a', should last value ('km')
        t1['a'].unit = 'cm'
        t2['a'].unit = 'm'
        t4['a'].unit = 'km'

        # Key col 'a' format should take last when all match
        t1['a'].info.format = '%0d'
        t2['a'].info.format = '%0d'
        t4['a'].info.format = '%0d'

        # Key col 'b', take first value 't1_b'
        t1['b'].info.description = 't1_b'

        # Key col 'b', take first non-empty value '%6s'
        t4['b'].info.format = '%6s'

        # Key col 'a', should be merged meta
        t1['a'].info.meta.update(OrderedDict([('b', [1, 2]), ('c', {'a': 1}), ('d', 1)]))
        t2['a'].info.meta.update(OrderedDict([('b', [3, 4]), ('c', {'b': 1}), ('a', 1)]))
        t4['a'].info.meta.update(OrderedDict([('b', [5, 6]), ('c', {'c': 1}), ('e', 1)]))

        # Key col 'b', should be meta2
        t2['b'].info.meta.update(OrderedDict([('b', [3, 4]), ('c', {'b': 1}), ('a', 1)]))

        # All these should pass through
        t2['c'].unit = 'm'
        t2['c'].info.format = '%6s'
        t2['c'].info.description = 't2_c'

        with catch_warnings(metadata.MergeConflictWarning) as warning_lines:
            out = table.vstack([t1, t2, t4], join_type='outer')

        assert warning_lines[0].category == metadata.MergeConflictWarning
        assert ("In merged column 'a' the 'unit' attribute does not match (cm != m)"
                in str(warning_lines[0].message))
        assert warning_lines[1].category == metadata.MergeConflictWarning
        assert ("In merged column 'a' the 'unit' attribute does not match (m != km)"
                in str(warning_lines[1].message))
        assert out['a'].unit == 'km'
        assert out['a'].info.format == '%0d'
        assert out['b'].info.description == 't1_b'
        assert out['b'].info.format == '%6s'
        assert out['a'].info.meta == self.meta_merge
        assert out['b'].info.meta == OrderedDict([('b', [3, 4]), ('c', {'b': 1}), ('a', 1)])
        assert out['c'].info.unit == 'm'
        assert out['c'].info.format == '%6s'
        assert out['c'].info.description == 't2_c'

    def test_vstack_one_table(self, operation_table_type):
        self._setup(operation_table_type)
        """Regression test for issue #3313"""
        assert (self.t1 == table.vstack(self.t1)).all()
        assert (self.t1 == table.vstack([self.t1])).all()

    def test_mixin_functionality(self, mixin_cols):
        col = mixin_cols['m']
        len_col = len(col)
        t = table.QTable([col], names=['a'])
        cls_name = type(col).__name__

        # Vstack works for these classes:
        implemented_mixin_classes = ['Quantity', 'Angle', 'Time',
                                     'Latitude', 'Longitude',
                                     'EarthLocation']
        if cls_name in implemented_mixin_classes:
            out = table.vstack([t, t])
            assert len(out) == len_col * 2
            assert np.all(out['a'][:len_col] == col)
            assert np.all(out['a'][len_col:] == col)
        else:
            with pytest.raises(NotImplementedError) as err:
                table.vstack([t, t])
            assert ('vstack unavailable for mixin column type(s): {}'
                    .format(cls_name) in str(err))

        # Check for outer stack which requires masking.  Only Time supports
        # this currently.
        t2 = table.QTable([col], names=['b'])  # different from col name for t
        if cls_name == 'Time':
            out = table.vstack([t, t2], join_type='outer')
            assert len(out) == len_col * 2
            assert np.all(out['a'][:len_col] == col)
            assert np.all(out['b'][len_col:] == col)
            assert np.all(out['a'].mask == [False] * len_col + [True] * len_col)
            assert np.all(out['b'].mask == [True] * len_col + [False] * len_col)
            # check directly stacking mixin columns:
            out2 = table.vstack([t, t2['b']])
            assert np.all(out['a'] == out2['a'])
            assert np.all(out['b'] == out2['b'])
        else:
            with pytest.raises(NotImplementedError) as err:
                table.vstack([t, t2], join_type='outer')
            assert ('vstack requires masking' in str(err) or
                    'vstack unavailable' in str(err))


class TestHStack():

    def _setup(self, t_cls=Table):
        self.t1 = t_cls.read([' a    b',
                              ' 0. foo',
                              ' 1. bar'], format='ascii')

        self.t2 = t_cls.read([' a    b   c',
                              ' 2.  pez  4',
                              ' 3.  sez  5'], format='ascii')

        self.t3 = t_cls.read([' d    e',
                              ' 4.   7',
                              ' 5.   8',
                              ' 6.   9'], format='ascii')
        self.t4 = t_cls(self.t1, copy=True, masked=True)
        self.t4['a'].name = 'f'
        self.t4['b'].name = 'g'

        # The following table has meta-data that conflicts with t1
        self.t5 = t_cls(self.t1, copy=True)

        self.t1.meta.update(OrderedDict([('b', [1, 2]), ('c', {'a': 1}), ('d', 1)]))
        self.t2.meta.update(OrderedDict([('b', [3, 4]), ('c', {'b': 1}), ('a', 1)]))
        self.t4.meta.update(OrderedDict([('b', [5, 6]), ('c', {'c': 1}), ('e', 1)]))
        self.t5.meta.update(OrderedDict([('b', 3), ('c', 'k'), ('d', 1)]))
        self.meta_merge = OrderedDict([('b', [1, 2, 3, 4, 5, 6]),
                                       ('c', {'a': 1, 'b': 1, 'c': 1}),
                                       ('d', 1),
                                       ('a', 1),
                                       ('e', 1)])

    def test_stack_same_table(self, operation_table_type):
        """
        From #2995, test that hstack'ing references to the same table has the
        expected output.
        """
        self._setup(operation_table_type)
        out = table.hstack([self.t1, self.t1])
        assert out.pformat() == ['a_1 b_1 a_2 b_2',
                                 '--- --- --- ---',
                                 '0.0 foo 0.0 foo',
                                 '1.0 bar 1.0 bar']

    def test_stack_rows(self, operation_table_type):
        self._setup(operation_table_type)
        out = table.hstack([self.t1[0], self.t2[1]])
        assert out.pformat() == ['a_1 b_1 a_2 b_2  c ',
                                 '--- --- --- --- ---',
                                 '0.0 foo 3.0 sez   5']

    def test_stack_columns(self, operation_table_type):
        self._setup(operation_table_type)
        out = table.hstack([self.t1, self.t2['c']])
        assert type(out['a']) is type(self.t1['a'])
        assert type(out['b']) is type(self.t1['b'])
        assert type(out['c']) is type(self.t2['c'])
        assert out.pformat() == [' a   b   c ',
                                 '--- --- ---',
                                 '0.0 foo   4',
                                 '1.0 bar   5']

    def test_table_meta_merge(self, operation_table_type):
        self._setup(operation_table_type)
        out = table.hstack([self.t1, self.t2, self.t4], join_type='inner')
        assert out.meta == self.meta_merge

    def test_table_meta_merge_conflict(self, operation_table_type):
        self._setup(operation_table_type)

        with catch_warnings() as w:
            out = table.hstack([self.t1, self.t5], join_type='inner')
        assert len(w) == 2

        assert out.meta == self.t5.meta

        with catch_warnings() as w:
            out = table.hstack([self.t1, self.t5], join_type='inner', metadata_conflicts='warn')
        assert len(w) == 2

        assert out.meta == self.t5.meta

        with catch_warnings() as w:
            out = table.hstack([self.t1, self.t5], join_type='inner', metadata_conflicts='silent')
        assert len(w) == 0

        assert out.meta == self.t5.meta

        with pytest.raises(MergeConflictError):
            out = table.hstack([self.t1, self.t5], join_type='inner', metadata_conflicts='error')

        with pytest.raises(ValueError):
            out = table.hstack([self.t1, self.t5], join_type='inner', metadata_conflicts='nonsense')

    def test_bad_input_type(self, operation_table_type):
        self._setup(operation_table_type)
        with pytest.raises(ValueError):
            table.hstack([])
        with pytest.raises(TypeError):
            table.hstack(1)
        with pytest.raises(TypeError):
            table.hstack([self.t2, 1])
        with pytest.raises(ValueError):
            table.hstack([self.t1, self.t2], join_type='invalid join type')

    def test_stack_basic(self, operation_table_type):
        self._setup(operation_table_type)
        t1 = self.t1
        t2 = self.t2
        t3 = self.t3
        t4 = self.t4

        out = table.hstack([t1, t2], join_type='inner')
        assert out.masked is False
        assert type(out) is operation_table_type
        assert type(out['a_1']) is type(t1['a'])
        assert type(out['b_1']) is type(t1['b'])
        assert type(out['a_2']) is type(t2['a'])
        assert type(out['b_2']) is type(t2['b'])
        assert out.pformat() == ['a_1 b_1 a_2 b_2  c ',
                                 '--- --- --- --- ---',
                                 '0.0 foo 2.0 pez   4',
                                 '1.0 bar 3.0 sez   5']

        # stacking as a list gives same result
        out_list = table.hstack([t1, t2], join_type='inner')
        assert out.pformat() == out_list.pformat()

        out = table.hstack([t1, t2], join_type='outer')
        assert out.pformat() == out_list.pformat()

        out = table.hstack([t1, t2, t3, t4], join_type='outer')
        assert out.pformat() == ['a_1 b_1 a_2 b_2  c   d   e   f   g ',
                                 '--- --- --- --- --- --- --- --- ---',
                                 '0.0 foo 2.0 pez   4 4.0   7 0.0 foo',
                                 '1.0 bar 3.0 sez   5 5.0   8 1.0 bar',
                                 ' --  --  --  --  -- 6.0   9  --  --']

        out = table.hstack([t1, t2, t3, t4], join_type='inner')
        assert out.pformat() == ['a_1 b_1 a_2 b_2  c   d   e   f   g ',
                                 '--- --- --- --- --- --- --- --- ---',
                                 '0.0 foo 2.0 pez   4 4.0   7 0.0 foo',
                                 '1.0 bar 3.0 sez   5 5.0   8 1.0 bar']

    def test_stack_incompatible(self, operation_table_type):
        self._setup(operation_table_type)
        # For join_type exact, which will fail here because n_rows
        # does not match
        with pytest.raises(TableMergeError):
            table.hstack([self.t1, self.t3], join_type='exact')

    def test_hstack_one_masked(self, operation_table_type):
        if operation_table_type is QTable:
            pytest.xfail()
        self._setup(operation_table_type)
        t1 = self.t1
        t2 = operation_table_type(t1, copy=True, masked=True)
        t2.meta.clear()
        t2['b'].mask[1] = True
        assert table.hstack([t1, t2]).pformat() == ['a_1 b_1 a_2 b_2',
                                                    '--- --- --- ---',
                                                    '0.0 foo 0.0 foo',
                                                    '1.0 bar 1.0  --']

    def test_table_col_rename(self, operation_table_type):
        self._setup(operation_table_type)
        out = table.hstack([self.t1, self.t2], join_type='inner',
                           uniq_col_name='{table_name}_{col_name}',
                           table_names=('left', 'right'))
        assert out.masked is False
        assert out.pformat() == ['left_a left_b right_a right_b  c ',
                                 '------ ------ ------- ------- ---',
                                 '   0.0    foo     2.0     pez   4',
                                 '   1.0    bar     3.0     sez   5']

    def test_col_meta_merge(self, operation_table_type):
        self._setup(operation_table_type)
        t1 = self.t1
        t3 = self.t3[:2]
        t4 = self.t4

        # Just set a bunch of meta and make sure it is the same in output
        meta1 = OrderedDict([('b', [1, 2]), ('c', {'a': 1}), ('d', 1)])
        t1['a'].unit = 'cm'
        t1['b'].info.description = 't1_b'
        t4['f'].info.format = '%6s'
        t1['b'].info.meta.update(meta1)
        t3['d'].info.meta.update(OrderedDict([('b', [3, 4]), ('c', {'b': 1}), ('a', 1)]))
        t4['g'].info.meta.update(OrderedDict([('b', [5, 6]), ('c', {'c': 1}), ('e', 1)]))
        t3['e'].info.meta.update(OrderedDict([('b', [3, 4]), ('c', {'b': 1}), ('a', 1)]))
        t3['d'].unit = 'm'
        t3['d'].info.format = '%6s'
        t3['d'].info.description = 't3_c'

        with catch_warnings(metadata.MergeConflictWarning) as warning_lines:
            out = table.hstack([t1, t3, t4], join_type='exact')

        assert len(warning_lines) == 0

        for t in [t1, t3, t4]:
            for name in t.colnames:
                for attr in ('meta', 'unit', 'format', 'description'):
                    assert getattr(out[name].info, attr) == getattr(t[name].info, attr)

        # Make sure we got a copy of meta, not ref
        t1['b'].info.meta['b'] = None
        assert out['b'].info.meta['b'] == [1, 2]

    def test_hstack_one_table(self, operation_table_type):
        self._setup(operation_table_type)
        """Regression test for issue #3313"""
        assert (self.t1 == table.hstack(self.t1)).all()
        assert (self.t1 == table.hstack([self.t1])).all()

    def test_mixin_functionality(self, mixin_cols):
        col1 = mixin_cols['m']
        col2 = col1[2:4]  # Shorter version of col1
        t1 = table.QTable([col1])
        t2 = table.QTable([col2])

        cls_name = type(col1).__name__

        out = table.hstack([t1, t2], join_type='inner')
        assert type(out['col0_1']) is type(out['col0_2'])
        assert len(out) == len(col2)

        # Check that columns are as expected.
        if cls_name == 'SkyCoord':
            assert skycoord_equal(out['col0_1'], col1[:len(col2)])
            assert skycoord_equal(out['col0_2'], col2)
        else:
            assert np.all(out['col0_1'] == col1[:len(col2)])
            assert np.all(out['col0_2'] == col2)

        # Time class supports masking, all other mixins do not
        if cls_name == 'Time':
            out = table.hstack([t1, t2], join_type='outer')
            assert len(out) == len(t1)
            assert np.all(out['col0_1'] == col1)
            assert np.all(out['col0_2'][:len(col2)] == col2)
            assert np.all(out['col0_2'].mask == [False, False, True, True])

            # check directly stacking mixin columns:
            out2 = table.hstack([t1, t2['col0']], join_type='outer')
            assert np.all(out['col0_1'] == out2['col0_1'])
            assert np.all(out['col0_2'] == out2['col0_2'])
        else:
            with pytest.raises(NotImplementedError) as err:
                table.hstack([t1, t2], join_type='outer')
            assert 'hstack requires masking' in str(err)


def test_unique(operation_table_type):
    t = operation_table_type.read(
        [' a b  c  d',
         ' 2 b 7.0 0',
         ' 1 c 3.0 5',
         ' 2 b 6.0 2',
         ' 2 a 4.0 3',
         ' 1 a 1.0 7',
         ' 2 b 5.0 1',
         ' 0 a 0.0 4',
         ' 1 a 2.0 6',
         ' 1 c 3.0 5',
        ], format='ascii')

    tu = operation_table_type(np.sort(t[:-1]))

    t_all = table.unique(t)
    assert sort_eq(t_all.pformat(), tu.pformat())
    t_s = t.copy()
    del t_s['b', 'c', 'd']
    t_all = table.unique(t_s)
    assert sort_eq(t_all.pformat(), [' a ',
                                     '---',
                                     '  0',
                                     '  1',
                                     '  2'])

    key1 = 'a'
    t1a = table.unique(t, key1)
    assert sort_eq(t1a.pformat(), [' a   b   c   d ',
                                   '--- --- --- ---',
                                   '  0   a 0.0   4',
                                   '  1   c 3.0   5',
                                   '  2   b 7.0   0'])
    t1b = table.unique(t, key1, keep='last')
    assert sort_eq(t1b.pformat(), [' a   b   c   d ',
                                   '--- --- --- ---',
                                   '  0   a 0.0   4',
                                   '  1   c 3.0   5',
                                   '  2   b 5.0   1'])
    t1c = table.unique(t, key1, keep='none')
    assert sort_eq(t1c.pformat(), [' a   b   c   d ',
                                   '--- --- --- ---',
                                   '  0   a 0.0   4'])

    key2 = ['a', 'b']
    t2a = table.unique(t, key2)
    assert sort_eq(t2a.pformat(), [' a   b   c   d ',
                                   '--- --- --- ---',
                                   '  0   a 0.0   4',
                                   '  1   a 1.0   7',
                                   '  1   c 3.0   5',
                                   '  2   a 4.0   3',
                                   '  2   b 7.0   0'])

    t2b = table.unique(t, key2, keep='last')
    assert sort_eq(t2b.pformat(), [' a   b   c   d ',
                                   '--- --- --- ---',
                                   '  0   a 0.0   4',
                                   '  1   a 2.0   6',
                                   '  1   c 3.0   5',
                                   '  2   a 4.0   3',
                                   '  2   b 5.0   1'])
    t2c = table.unique(t, key2, keep='none')
    assert sort_eq(t2c.pformat(), [' a   b   c   d ',
                                   '--- --- --- ---',
                                   '  0   a 0.0   4',
                                   '  2   a 4.0   3'])

    key2 = ['a', 'a']
    with pytest.raises(ValueError) as exc:
        t2a = table.unique(t, key2)
    assert exc.value.args[0] == "duplicate key names"

    with pytest.raises(ValueError) as exc:
        table.unique(t, key2, keep=True)
    assert exc.value.args[0] == (
        "'keep' should be one of 'first', 'last', 'none'")

    t1_m = operation_table_type(t1a, masked=True)
    t1_m['a'].mask[1] = True

    with pytest.raises(ValueError) as exc:
        t1_mu = table.unique(t1_m)
    assert exc.value.args[0] == (
        "cannot use columns with masked values as keys; "
        "remove column 'a' from keys and rerun unique()")

    t1_mu = table.unique(t1_m, silent=True)
    assert t1_mu.pformat() == [' a   b   c   d ',
                               '--- --- --- ---',
                               '  0   a 0.0   4',
                               '  2   b 7.0   0',
                               ' --   c 3.0   5']

    with pytest.raises(ValueError) as e:
        t1_mu = table.unique(t1_m, silent=True, keys='a')

    t1_m = operation_table_type(t, masked=True)
    t1_m['a'].mask[1] = True
    t1_m['d'].mask[3] = True

    # Test that multiple masked key columns get removed in the correct
    # order
    t1_mu = table.unique(t1_m, keys=['d', 'a', 'b'], silent=True)
    assert t1_mu.pformat() == [' a   b   c   d ',
                               '--- --- --- ---',
                               '  2   a 4.0  --',
                               '  2   b 7.0   0',
                               ' --   c 3.0   5']


def test_vstack_bytes(operation_table_type):
    """
    Test for issue #5617 when vstack'ing bytes columns in Py3.
    This is really an upsteam numpy issue numpy/numpy/#8403.
    """
    t = operation_table_type([[b'a']], names=['a'])
    assert t['a'].itemsize == 1

    t2 = table.vstack([t, t])
    assert len(t2) == 2
    assert t2['a'].itemsize == 1


def test_vstack_unicode():
    """
    Test for problem related to issue #5617 when vstack'ing *unicode*
    columns.  In this case the character size gets multiplied by 4.
    """
    t = table.Table([[u'a']], names=['a'])
    assert t['a'].itemsize == 4  # 4-byte / char for U dtype

    t2 = table.vstack([t, t])
    assert len(t2) == 2
    assert t2['a'].itemsize == 4


def test_get_out_class():
    c = table.Column([1, 2])
    mc = table.MaskedColumn([1, 2])
    q = [1, 2] * u.m

    assert _get_out_class([c, mc]) is mc.__class__
    assert _get_out_class([mc, c]) is mc.__class__
    assert _get_out_class([c, c]) is c.__class__
    assert _get_out_class([c]) is c.__class__

    with pytest.raises(ValueError):
        _get_out_class([c, q])

    with pytest.raises(ValueError):
        _get_out_class([q, c])


def test_masking_required_exception():
    """
    Test that outer join, hstack and vstack fail for a mixin column which
    does not support masking.
    """
    col = [1, 2, 3, 4] * u.m
    t1 = table.QTable([[1, 2, 3, 4], col], names=['a', 'b'])
    t2 = table.QTable([[1, 2], col[:2]], names=['a', 'c'])

    with pytest.raises(NotImplementedError) as err:
        table.vstack([t1, t2], join_type='outer')
    assert 'vstack requires masking' in str(err)

    with pytest.raises(NotImplementedError) as err:
        table.hstack([t1, t2], join_type='outer')
    assert 'hstack requires masking' in str(err)

    with pytest.raises(NotImplementedError) as err:
        table.join(t1, t2, join_type='outer')
    assert 'join requires masking' in str(err)


def test_stack_columns():
    c = table.Column([1, 2])
    mc = table.MaskedColumn([1, 2])
    q = [1, 2] * u.m
    time = Time(['2001-01-02T12:34:56', '2001-02-03T00:01:02'])
    sc = SkyCoord([1, 2], [3, 4], unit='deg')
    cq = table.Column([11, 22], unit=u.m)

    t = table.hstack([c, q])
    assert t.__class__ is table.QTable
    assert t.masked is False
    t = table.hstack([q, c])
    assert t.__class__ is table.QTable
    assert t.masked is False

    t = table.hstack([mc, q])
    assert t.__class__ is table.QTable
    assert t.masked is True

    t = table.hstack([c, mc])
    assert t.__class__ is table.Table
    assert t.masked is True

    t = table.vstack([q, q])
    assert t.__class__ is table.QTable

    t = table.vstack([c, c])
    assert t.__class__ is table.Table

    t = table.hstack([c, time])
    assert t.__class__ is table.Table
    t = table.hstack([c, sc])
    assert t.__class__ is table.Table
    t = table.hstack([q, time, sc])
    assert t.__class__ is table.QTable

    with pytest.raises(ValueError):
        table.vstack([c, q])

    with pytest.raises(ValueError):
        t = table.vstack([q, cq])
