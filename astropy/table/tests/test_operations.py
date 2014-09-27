# Licensed under a 3-clause BSD style license - see LICENSE.rst

# TEST_UNICODE_LITERALS

from distutils import version
import warnings

import numpy as np

from ...tests.helper import pytest, catch_warnings
from ...table import Table
from ...utils import OrderedDict, metadata
from ...utils.metadata import MergeConflictError
from .. import np_utils
from ... import table
from ... import log

def sort_eq(list1, list2):
    return sorted(list1) == sorted(list2)


class TestJoin():

    def setup_method(self, method):
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
        self.t1 = Table.read(lines1, format='ascii')
        self.t2 = Table.read(lines2, format='ascii')
        self.t3 = Table(self.t2, copy=True)

        self.t1.meta.update(OrderedDict([('b', [1, 2]), ('c', {'a': 1}), ('d', 1)]))
        self.t2.meta.update(OrderedDict([('b', [3, 4]), ('c', {'b': 1}), ('a', 1)]))
        self.t3.meta.update(OrderedDict([('b', 3), ('c', [1, 2]), ('d', 2), ('a', 1)]))

        self.meta_merge = OrderedDict([('b', [1, 2, 3, 4]),
                                       ('c', {'a': 1, 'b': 1}),
                                       ('d', 1),
                                       ('a', 1)])

    def test_table_meta_merge(self):
        out = table.join(self.t1, self.t2, join_type='inner')
        assert out.meta == self.meta_merge

    def test_table_meta_merge_conflict(self):

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

    def test_both_unmasked_inner(self):
        t1 = self.t1
        t2 = self.t2

        # Basic join with default parameters (inner join on common keys)
        t12 = table.join(t1, t2)
        assert t12.masked is False
        assert sort_eq(t12.pformat(), [' a   b   c   d ',
                                       '--- --- --- ---',
                                       '  1 foo  L2  R1',
                                       '  1 foo  L2  R2',
                                       '  2 bar  L4  R3'])

        # Table meta merged properly
        assert t12.meta == self.meta_merge

    def test_both_unmasked_left_right_outer(self):
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
        assert np.all(t12a._data == t12b._data)

    def test_both_unmasked_single_key_inner(self):
        t1 = self.t1
        t2 = self.t2

        # Inner join on 'a' column
        t12 = table.join(t1, t2, keys='a')
        assert t12.masked is False
        assert sort_eq(t12.pformat(), [' a  b_1  c  b_2  d ',
                                       '--- --- --- --- ---',
                                       '  1 foo  L2 foo  R1',
                                       '  1 foo  L2 foo  R2',
                                       '  1 bar  L3 foo  R1',
                                       '  1 bar  L3 foo  R2',
                                       '  2 bar  L4 bar  R3'])

    def test_both_unmasked_single_key_left_right_outer(self):
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

    def test_masked_unmasked(self):
        t1 = self.t1
        t1m = Table(self.t1, masked=True)
        t2 = self.t2

        # Result should be masked even though not req'd by inner join
        t1m2 = table.join(t1m, t2, join_type='inner')
        assert t1m2.masked is True

        # Result should match non-masked result
        t12 = table.join(t1, t2)
        assert np.all(t12._data == np.array(t1m2._data))

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

    def test_masked_masked(self):
        """Two masked tables"""
        t1 = self.t1
        t1m = Table(self.t1, masked=True)
        t2 = self.t2
        t2m = Table(self.t2, masked=True)

        # Result should be masked even though not req'd by inner join
        t1m2m = table.join(t1m, t2m, join_type='inner')
        assert t1m2m.masked is True

        # Result should match non-masked result
        t12 = table.join(t1, t2)
        assert np.all(t12._data == np.array(t1m2m._data))

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

    def test_col_rename(self):
        """
        Test auto col renaming when there is a conflict.  Use
        non-default values of uniq_col_name and table_names.
        """
        t1 = self.t1
        t2 = self.t2
        t12 = table.join(t1, t2, uniq_col_name='x_{table_name}_{col_name}_y',
                         table_names=['L', 'R'], keys='a')
        assert t12.colnames == ['a', 'x_L_b_y', 'c', 'x_R_b_y', 'd']

    def test_rename_conflict(self):
        """
        Test that auto-column rename fails because of a conflict
        with an existing column
        """
        t1 = self.t1
        t2 = self.t2
        t1['b_1'] = 1  # Add a new column b_1 that will conflict with auto-rename
        with pytest.raises(np_utils.TableMergeError):
            table.join(t1, t2, keys='a')

    def test_missing_keys(self):
        """Merge on a key column that doesn't exist"""
        t1 = self.t1
        t2 = self.t2
        with pytest.raises(np_utils.TableMergeError):
            table.join(t1, t2, keys=['a', 'not there'])

    def test_bad_join_type(self):
        """Bad join_type input"""
        t1 = self.t1
        t2 = self.t2
        with pytest.raises(ValueError):
            table.join(t1, t2, join_type='illegal value')

    def test_no_common_keys(self):
        """Merge tables with no common keys"""
        t1 = self.t1
        t2 = self.t2
        del t1['a']
        del t1['b']
        del t2['a']
        del t2['b']
        with pytest.raises(np_utils.TableMergeError):
            table.join(t1, t2)

    def test_masked_key_column(self):
        """Merge on a key column that has a masked element"""
        t1 = self.t1
        t2 = Table(self.t2, masked=True)
        table.join(t1, t2)  # OK
        t2['a'].mask[0] = True
        with pytest.raises(np_utils.TableMergeError):
            table.join(t1, t2)

    def test_col_meta_merge(self):
        t1 = self.t1
        t2 = self.t2
        t2.rename_column('d', 'c')  # force col conflict and renaming
        meta1 = OrderedDict([('b', [1, 2]), ('c', {'a': 1}), ('d', 1)])
        meta2 = OrderedDict([('b', [3, 4]), ('c', {'b': 1}), ('a', 1)])

        # Key col 'a', should first value ('cm')
        t1['a'].unit = 'cm'
        t2['a'].unit = 'm'
        # Key col 'b', take first value 't1_b'
        t1['b'].description = 't1_b'
        # Key col 'b', take first non-empty value 't1_b'
        t2['b'].format = '%6s'
        # Key col 'a', should be merged meta
        t1['a'].meta = meta1
        t2['a'].meta = meta2
        # Key col 'b', should be meta2
        t2['b'].meta = meta2

        # All these should pass through
        t1['c'].unit = 'cm'
        t1['c'].format = '%3s'
        t1['c'].description = 't1_c'

        t2['c'].unit = 'm'
        t2['c'].format = '%6s'
        t2['c'].description = 't2_c'

        with catch_warnings(metadata.MergeConflictWarning) as warning_lines:
            t12 = table.join(t1, t2, keys=['a', 'b'])

            assert t12['a'].unit == 'm'
            assert t12['b'].description == 't1_b'
            assert t12['b'].format == '%6s'
            assert t12['a'].meta == self.meta_merge
            assert t12['b'].meta == meta2
            assert t12['c_1'].unit == 'cm'
            assert t12['c_1'].format == '%3s'
            assert t12['c_1'].description == 't1_c'
            assert t12['c_2'].unit == 'm'
            assert t12['c_2'].format == '%6s'
            assert t12['c_2'].description == 't2_c'

            assert warning_lines[0].category == metadata.MergeConflictWarning
            assert ("In merged column 'a' the 'unit' attribute does not match (cm != m)"
                    in str(warning_lines[0].message))

    def test_join_multidimensional(self):

        # Regression test for #2984, which was an issue where join did not work
        # on multi-dimensional columns.

        t1 = Table()
        t1['a'] = [1,2,3]
        t1['b'] = np.ones((3,4))

        t2 = Table()
        t2['a'] = [1,2,3]
        t2['c'] = [4,5,6]

        t3 = table.join(t1, t2)

        np.testing.assert_allclose(t3['a'], t1['a'])
        np.testing.assert_allclose(t3['b'], t1['b'])
        np.testing.assert_allclose(t3['c'], t2['c'])


class TestVStack():

    def setup_method(self, method):

        self.t1 = Table.read([' a   b',
                              ' 0 foo',
                              ' 1 bar'], format='ascii')

        self.t2 = Table.read([' a   b   c',
                              ' 2  pez  4',
                              ' 3  sez  5'], format='ascii')

        self.t3 = Table.read([' a   b',
                              ' 4   7',
                              ' 5   8',
                              ' 6   9'], format='ascii')
        self.t4 = Table(self.t1, copy=True, masked=True)

        # The following table has meta-data that conflicts with t1
        self.t5 = Table(self.t1, copy=True)

        self.t1.meta.update(OrderedDict([('b', [1, 2]), ('c', {'a': 1}), ('d', 1)]))
        self.t2.meta.update(OrderedDict([('b', [3, 4]), ('c', {'b': 1}), ('a', 1)]))
        self.t4.meta.update(OrderedDict([('b', [5, 6]), ('c', {'c': 1}), ('e', 1)]))
        self.t5.meta.update(OrderedDict([('b', 3), ('c', 'k'), ('d', 1)]))
        self.meta_merge = OrderedDict([('b', [1, 2, 3, 4, 5, 6]),
                                       ('c', {'a': 1, 'b': 1, 'c': 1}),
                                       ('d', 1),
                                       ('a', 1),
                                       ('e', 1)])

    def test_stack_rows(self):
        t2 = self.t1.copy()
        t2.meta.clear()
        out = table.vstack([self.t1, t2[1]])
        assert out.pformat() == [' a   b ',
                                 '--- ---',
                                 '  0 foo',
                                 '  1 bar',
                                 '  1 bar']

    def test_table_meta_merge(self):
        out = table.vstack([self.t1, self.t2, self.t4], join_type='inner')
        assert out.meta == self.meta_merge

    def test_table_meta_merge_conflict(self):

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


    def test_bad_input_type(self):
        with pytest.raises(TypeError):
            table.vstack(1)
        with pytest.raises(TypeError):
            table.vstack([self.t2, 1])
        with pytest.raises(ValueError):
            table.vstack([self.t1, self.t2], join_type='invalid join type')

    def test_stack_basic(self):
        t1 = self.t1
        t2 = self.t2
        t4 = self.t4

        t12 = table.vstack([t1, t2], join_type='inner')
        assert t12.masked is False
        assert t12.pformat() == [' a   b ',
                                 '--- ---',
                                 '  0 foo',
                                 '  1 bar',
                                 '  2 pez',
                                 '  3 sez']


        t12 = table.vstack([t1, t2], join_type='outer')
        assert t12.pformat() == [' a   b   c ',
                                 '--- --- ---',
                                 '  0 foo  --',
                                 '  1 bar  --',
                                 '  2 pez   4',
                                 '  3 sez   5']

        t124 = table.vstack([t1, t2, t4], join_type='outer')
        assert t124.pformat() == [' a   b   c ',
                                  '--- --- ---',
                                  '  0 foo  --',
                                  '  1 bar  --',
                                  '  2 pez   4',
                                  '  3 sez   5',
                                  '  0 foo  --',
                                  '  1 bar  --']

        t124 = table.vstack([t1, t2, t4], join_type='inner')
        assert t124.pformat() == [' a   b ',
                                  '--- ---',
                                  '  0 foo',
                                  '  1 bar',
                                  '  2 pez',
                                  '  3 sez',
                                  '  0 foo',
                                  '  1 bar']

    def test_stack_incompatible(self):
        with pytest.raises(np_utils.TableMergeError) as excinfo:
            table.vstack([self.t1, self.t3], join_type='inner')
        assert "The 'b' columns have incompatible types:" in str(excinfo)

        with pytest.raises(np_utils.TableMergeError) as excinfo:
            table.vstack([self.t1, self.t3], join_type='outer')
        assert "The 'b' columns have incompatible types:" in str(excinfo)

        with pytest.raises(np_utils.TableMergeError):
            table.vstack([self.t1, self.t2], join_type='exact')

    def test_vstack_one_masked(self):
        t1 = self.t1
        t4 = self.t4
        t4['b'].mask[1] = True
        assert table.vstack([t1, t4]).pformat() == [' a   b ',
                                                    '--- ---',
                                                    '  0 foo',
                                                    '  1 bar',
                                                    '  0 foo',
                                                    '  1  --']

    def test_col_meta_merge(self):
        t1 = self.t1
        t2 = self.t2
        t4 = self.t4

        # Key col 'a', should last value ('km')
        t1['a'].unit = 'cm'
        t2['a'].unit = 'm'
        t4['a'].unit = 'km'

        # Key col 'a' format should take last when all match
        t1['a'].format = '%0d'
        t2['a'].format = '%0d'
        t4['a'].format = '%0d'

        # Key col 'b', take first value 't1_b'
        t1['b'].description = 't1_b'

        # Key col 'b', take first non-empty value '%6s'
        t4['b'].format = '%6s'

        # Key col 'a', should be merged meta
        t1['a'].meta.update(OrderedDict([('b', [1, 2]), ('c', {'a': 1}), ('d', 1)]))
        t2['a'].meta.update(OrderedDict([('b', [3, 4]), ('c', {'b': 1}), ('a', 1)]))
        t4['a'].meta.update(OrderedDict([('b', [5, 6]), ('c', {'c': 1}), ('e', 1)]))

        # Key col 'b', should be meta2
        t2['b'].meta.update(OrderedDict([('b', [3, 4]), ('c', {'b': 1}), ('a', 1)]))

        # All these should pass through
        t2['c'].unit = 'm'
        t2['c'].format = '%6s'
        t2['c'].description = 't2_c'

        with catch_warnings(metadata.MergeConflictWarning) as warning_lines:
            out = table.vstack([t1, t2, t4], join_type='outer')

            assert out['a'].unit == 'km'
            assert out['a'].format == '%0d'
            assert out['b'].description == 't1_b'
            assert out['b'].format == '%6s'
            assert out['a'].meta == self.meta_merge
            assert out['b'].meta == OrderedDict([('b', [3, 4]), ('c', {'b': 1}), ('a', 1)])
            assert out['c'].unit == 'm'
            assert out['c'].format == '%6s'
            assert out['c'].description == 't2_c'

            assert warning_lines[0].category == metadata.MergeConflictWarning
            assert ("In merged column 'a' the 'unit' attribute does not match (cm != m)"
                    in str(warning_lines[0].message))
            assert warning_lines[1].category == metadata.MergeConflictWarning
            assert ("In merged column 'a' the 'unit' attribute does not match (m != km)"
                    in str(warning_lines[1].message))


class TestHStack():

    def setup_method(self, method):
        self.t1 = Table.read([' a   b',
                              ' 0 foo',
                              ' 1 bar'], format='ascii')

        self.t2 = Table.read([' a   b   c',
                              ' 2  pez  4',
                              ' 3  sez  5'], format='ascii')

        self.t3 = Table.read([' d   e',
                              ' 4   7',
                              ' 5   8',
                              ' 6   9'], format='ascii')
        self.t4 = Table(self.t1, copy=True, masked=True)
        self.t4['a'].name = 'f'
        self.t4['b'].name = 'g'

        # The following table has meta-data that conflicts with t1
        self.t5 = Table(self.t1, copy=True)

        self.t1.meta.update(OrderedDict([('b', [1, 2]), ('c', {'a': 1}), ('d', 1)]))
        self.t2.meta.update(OrderedDict([('b', [3, 4]), ('c', {'b': 1}), ('a', 1)]))
        self.t4.meta.update(OrderedDict([('b', [5, 6]), ('c', {'c': 1}), ('e', 1)]))
        self.t5.meta.update(OrderedDict([('b', 3), ('c', 'k'), ('d', 1)]))
        self.meta_merge = OrderedDict([('b', [1, 2, 3, 4, 5, 6]),
                                       ('c', {'a': 1, 'b': 1, 'c': 1}),
                                       ('d', 1),
                                       ('a', 1),
                                       ('e', 1)])

    def test_stack_rows(self):
        out = table.hstack([self.t1[0], self.t2[1]])
        assert out.pformat() == ['a_1 b_1 a_2 b_2  c ',
                                 '--- --- --- --- ---',
                                 '  0 foo   3 sez   5']

    def test_table_meta_merge(self):
        out = table.hstack([self.t1, self.t2, self.t4], join_type='inner')
        assert out.meta == self.meta_merge

    def test_table_meta_merge_conflict(self):

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

    def test_bad_input_type(self):
        with pytest.raises(TypeError):
            table.hstack(1)
        with pytest.raises(TypeError):
            table.hstack([self.t2, 1])
        with pytest.raises(ValueError):
            table.hstack([self.t1, self.t2], join_type='invalid join type')

    def test_stack_basic(self):
        t1 = self.t1
        t2 = self.t2
        t3 = self.t3
        t4 = self.t4

        out = table.hstack([t1, t2], join_type='inner')
        assert out.masked is False
        assert out.pformat() == ['a_1 b_1 a_2 b_2  c ',
                                 '--- --- --- --- ---',
                                 '  0 foo   2 pez   4',
                                 '  1 bar   3 sez   5']

        # stacking as a list gives same result
        out_list = table.hstack([t1, t2], join_type='inner')
        assert out.pformat() == out_list.pformat()

        out = table.hstack([t1, t2], join_type='outer')
        assert out.pformat() == out_list.pformat()

        out = table.hstack([t1, t2, t3, t4], join_type='outer')
        assert out.pformat() == ['a_1 b_1 a_2 b_2  c   d   e   f   g ',
                                 '--- --- --- --- --- --- --- --- ---',
                                 '  0 foo   2 pez   4   4   7   0 foo',
                                 '  1 bar   3 sez   5   5   8   1 bar',
                                 ' --  --  --  --  --   6   9  --  --']

        out = table.hstack([t1, t2, t3, t4], join_type='inner')
        assert out.pformat() == ['a_1 b_1 a_2 b_2  c   d   e   f   g ',
                                 '--- --- --- --- --- --- --- --- ---',
                                 '  0 foo   2 pez   4   4   7   0 foo',
                                 '  1 bar   3 sez   5   5   8   1 bar']

    def test_stack_incompatible(self):
        # For join_type exact, which will fail here because n_rows
        # does not match
        with pytest.raises(np_utils.TableMergeError):
            table.hstack([self.t1, self.t3], join_type='exact')

    def test_hstack_one_masked(self):
        t1 = self.t1
        t2 = Table(t1, copy=True, masked=True)
        t2.meta.clear()
        t2['b'].mask[1] = True
        assert table.hstack([t1, t2]).pformat() == ['a_1 b_1 a_2 b_2',
                                                    '--- --- --- ---',
                                                    '  0 foo   0 foo',
                                                    '  1 bar   1  --']

    def test_table_col_rename(self):
        out = table.hstack([self.t1, self.t2], join_type='inner',
                           uniq_col_name='{table_name}_{col_name}',
                           table_names=('left', 'right'))
        assert out.masked is False
        assert out.pformat() == ['left_a left_b right_a right_b  c ',
                                 '------ ------ ------- ------- ---',
                                 '     0    foo       2     pez   4',
                                 '     1    bar       3     sez   5']

    def test_col_meta_merge(self):
        t1 = self.t1
        t3 = self.t3
        t4 = self.t4

        # Just set a bunch of meta and make sure it is the same in output
        meta1 = OrderedDict([('b', [1, 2]), ('c', {'a': 1}), ('d', 1)])
        t1['a'].unit = 'cm'
        t1['b'].description = 't1_b'
        t4['f'].format = '%6s'
        t1['b'].meta.update(meta1)
        t3['d'].meta.update(OrderedDict([('b', [3, 4]), ('c', {'b': 1}), ('a', 1)]))
        t4['g'].meta.update(OrderedDict([('b', [5, 6]), ('c', {'c': 1}), ('e', 1)]))
        t3['e'].meta.update(OrderedDict([('b', [3, 4]), ('c', {'b': 1}), ('a', 1)]))
        t3['d'].unit = 'm'
        t3['d'].format = '%6s'
        t3['d'].description = 't3_c'

        with catch_warnings(metadata.MergeConflictWarning) as warning_lines:
            out = table.vstack([t1, t3, t4], join_type='outer')

            for t in [t1, t3, t4]:
                for name in t.colnames:
                    for attr in ('meta', 'unit', 'format', 'description'):
                        assert getattr(out[name], attr) == getattr(t[name], attr)

            assert len(warning_lines) == 0

            # Make sure we got a copy of meta, not ref
            t1['b'].meta['b'] = None
            assert out['b'].meta['b'] == [1, 2]
