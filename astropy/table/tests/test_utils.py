# Licensed under a 3-clause BSD style license - see LICENSE.rst
from distutils import version
import numpy as np

from ...tests.helper import pytest
from ... import table
from ...io import ascii
from ...utils import OrderedDict
from .. import np_utils

NUMPY_LT_1P5 = version.LooseVersion(np.__version__) < version.LooseVersion('1.5')

class TestJoin():

    def setup_method(self, method):
        lines1 = [' a   b   c ',
                  '--- --- ---',
                  '  0 foo  L1',
                  '  1 foo  L2',
                  '  1 bar  L3',
                  '  2 bar  L4']
        lines2 = [' a   b   d ',
                  '--- --- ---',
                  '  1 foo  R1',
                  '  1 foo  R2',
                  '  2 bar  R3',
                  '  4 bar  R4']
        self.t1 = ascii.read(lines1, Reader=ascii.FixedWidthTwoLine)
        self.t2 = ascii.read(lines2, Reader=ascii.FixedWidthTwoLine)
        self.t1.meta.update(OrderedDict([('b', [1, 2]), ('c', {'a': 1}), ('d', 1)]))
        self.t2.meta.update(OrderedDict([('b', [3, 4]), ('c', {'b': 1}), ('a', 1)]))
        self.meta_merge = OrderedDict([('b', [1, 2, 3, 4]),
                                       ('c', {'a': 1, 'b': 1}),
                                       ('d', 1),
                                       ('a', 1)])
        
    def test_both_unmasked_inner(self):
        t1 = self.t1
        t2 = self.t2

        # Basic join with default parameters (inner join on common keys)
        t12 = t1.join(t2)
        assert t12.masked is False
        assert t12.pformat() == [' a   b   c   d ',
                                 '--- --- --- ---',
                                 '  1 foo  L2  R1',
                                 '  1 foo  L2  R2',
                                 '  2 bar  L4  R3']

        # Table meta merged properly
        assert t12.meta == self.meta_merge

    @pytest.mark.xfail('NUMPY_LT_1P5')
    def test_both_unmasked_left_right_outer(self):
        t1 = self.t1
        t2 = self.t2

        # Left join
        t12 = t1.join(t2, join_type='left')
        assert t12.masked is True
        assert t12.pformat() == [' a   b   c   d ',
                                 '--- --- --- ---',
                                 '  0 foo  L1  --',
                                 '  1 bar  L3  --',
                                 '  1 foo  L2  R1',
                                 '  1 foo  L2  R2',
                                 '  2 bar  L4  R3']

        # Right join
        t12 = t1.join(t2, join_type='right')
        assert t12.masked is True
        assert t12.pformat() == [' a   b   c   d ',
                                 '--- --- --- ---',
                                 '  1 foo  L2  R1',
                                 '  1 foo  L2  R2',
                                 '  2 bar  L4  R3',
                                 '  4 bar  --  R4']

        # Outer join
        t12 = t1.join(t2, join_type='outer')
        assert t12.masked is True
        assert t12.pformat() == [' a   b   c   d ',
                                 '--- --- --- ---',
                                 '  0 foo  L1  --',
                                 '  1 bar  L3  --',
                                 '  1 foo  L2  R1',
                                 '  1 foo  L2  R2',
                                 '  2 bar  L4  R3',
                                 '  4 bar  --  R4']

        # Check that the common keys are 'a', 'b'
        t12a = t1.join(t2, join_type='outer')
        t12b = t1.join(t2, join_type='outer', keys=['a', 'b'])
        assert np.all(t12a._data == t12b._data)

    def test_both_unmasked_single_key_inner(self):
        t1 = self.t1
        t2 = self.t2

        # Inner join on 'a' column
        t12 = t1.join(t2, keys='a')
        assert t12.masked is False
        assert t12.pformat() == [' a  b_1  c  b_2  d ',
                                 '--- --- --- --- ---',
                                 '  1 foo  L2 foo  R1',
                                 '  1 foo  L2 foo  R2',
                                 '  1 bar  L3 foo  R1',
                                 '  1 bar  L3 foo  R2',
                                 '  2 bar  L4 bar  R3']

    @pytest.mark.xfail('NUMPY_LT_1P5')
    def test_both_unmasked_single_key_left_right_outer(self):
        t1 = self.t1
        t2 = self.t2

        # Left join
        t12 = t1.join(t2, join_type='left', keys='a')
        assert t12.masked is True
        assert t12.pformat() == [' a  b_1  c  b_2  d ',
                                 '--- --- --- --- ---',
                                 '  0 foo  L1  --  --',
                                 '  1 foo  L2 foo  R1',
                                 '  1 foo  L2 foo  R2',
                                 '  1 bar  L3 foo  R1',
                                 '  1 bar  L3 foo  R2',
                                 '  2 bar  L4 bar  R3']

        # Right join
        t12 = t1.join(t2, join_type='right', keys='a')
        assert t12.masked is True
        assert t12.pformat() == [' a  b_1  c  b_2  d ',
                                 '--- --- --- --- ---',
                                 '  1 foo  L2 foo  R1',
                                 '  1 foo  L2 foo  R2',
                                 '  1 bar  L3 foo  R1',
                                 '  1 bar  L3 foo  R2',
                                 '  2 bar  L4 bar  R3',
                                 '  4  --  -- bar  R4']

        # Outer join
        t12 = t1.join(t2, join_type='outer', keys='a')
        assert t12.masked is True
        assert t12.pformat() == [' a  b_1  c  b_2  d ',
                                 '--- --- --- --- ---',
                                 '  0 foo  L1  --  --',
                                 '  1 foo  L2 foo  R1',
                                 '  1 foo  L2 foo  R2',
                                 '  1 bar  L3 foo  R1',
                                 '  1 bar  L3 foo  R2',
                                 '  2 bar  L4 bar  R3',
                                 '  4  --  -- bar  R4']


    @pytest.mark.xfail('NUMPY_LT_1P5')
    def test_masked_unmasked(self):
        t1 = self.t1
        t1m = table.Table(self.t1, masked=True)
        t2 = self.t2

        # Result should be masked even though not req'd by inner join
        t1m2 = t1m.join(t2, join_type='inner')
        assert t1m2.masked is True
        
        # Result should match non-masked result
        t12 = t1.join(t2)
        assert np.all(t12._data == np.array(t1m2._data))

        # Mask out some values in left table and make sure they propagate
        t1m['b'].mask[1] = True
        t1m['c'].mask[2] = True        
        t1m2 = t1m.join(t2, join_type='inner', keys='a')
        assert t1m2.pformat() == [' a  b_1  c  b_2  d ',
                                  '--- --- --- --- ---',
                                  '  1  --  L2 foo  R1',
                                  '  1  --  L2 foo  R2',
                                  '  1 bar  -- foo  R1',
                                  '  1 bar  -- foo  R2',
                                  '  2 bar  L4 bar  R3']

        t21m = t2.join(t1m, join_type='inner', keys='a')
        assert t21m.pformat() == [' a  b_1  d  b_2  c ',
                                  '--- --- --- --- ---',
                                  '  1 foo  R2  --  L2',
                                  '  1 foo  R2 bar  --',
                                  '  1 foo  R1  --  L2',
                                  '  1 foo  R1 bar  --',
                                  '  2 bar  R3 bar  L4']

    @pytest.mark.xfail('NUMPY_LT_1P5')
    def test_masked_masked(self):
        """Two masked tables"""
        t1 = self.t1
        t1m = table.Table(self.t1, masked=True)
        t2 = self.t2
        t2m = table.Table(self.t2, masked=True)

        # Result should be masked even though not req'd by inner join
        t1m2m = t1m.join(t2m, join_type='inner')
        assert t1m2m.masked is True
        
        # Result should match non-masked result
        t12 = t1.join(t2)
        assert np.all(t12._data == np.array(t1m2m._data))

        # Mask out some values in both tables and make sure they propagate
        t1m['b'].mask[1] = True
        t1m['c'].mask[2] = True        
        t2m['d'].mask[2] = True        
        t1m2m = t1m.join(t2m, join_type='inner', keys='a')
        assert t1m2m.pformat() == [' a  b_1  c  b_2  d ',
                                   '--- --- --- --- ---',
                                   '  1  --  L2 foo  R1',
                                   '  1  --  L2 foo  R2',
                                   '  1 bar  -- foo  R1',
                                   '  1 bar  -- foo  R2',
                                   '  2 bar  L4 bar  --']

    def test_col_rename(self):
        """
        Test auto col renaming when there is a conflict.  Use
        non-default values of uniq_col_name and table_names.
        """
        t1 = self.t1
        t2 = self.t2
        t12 = t1.join(t2, uniq_col_name='x_{table_name}_{col_name}_y',
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
            t1.join(t2, keys='a')

    def test_col_meta_merge(self):
        t1 = self.t1
        t2 = self.t2
        t2.rename_column('d', 'c')  # force col conflict and renaming
        meta1 = OrderedDict([('b', [1, 2]), ('c', {'a': 1}), ('d', 1)])
        meta2 = OrderedDict([('b', [3, 4]), ('c', {'b': 1}), ('a', 1)])

        # Key col 'a', should first value ('cm')
        t1['a'].units = 'cm'
        t2['a'].units = 'm'
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
        t1['c'].units = 'cm'
        t1['c'].format = '%3s'
        t1['c'].description = 't1_c'

        t2['c'].units = 'm'
        t2['c'].format = '%6s'
        t2['c'].description = 't2_c'

        t12 = t1.join(t2, keys=['a', 'b'])
        assert t12['a'].units == 'cm'
        assert t12['b'].description == 't1_b'
        assert t12['b'].format == '%6s'
        assert t12['a'].meta == self.meta_merge
        assert t12['b'].meta == meta2
        assert t12['c_1'].units == 'cm'
        assert t12['c_1'].format == '%3s'
        assert t12['c_1'].description == 't1_c'
        assert t12['c_2'].units == 'm'
        assert t12['c_2'].format == '%6s'
        assert t12['c_2'].description == 't2_c'
