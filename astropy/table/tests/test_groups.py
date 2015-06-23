# Licensed under a 3-clause BSD style license - see LICENSE.rst

# TEST_UNICODE_LITERALS

import numpy as np

from ...tests.helper import pytest, catch_warnings
from ...table import Table, Column
from ...utils.exceptions import AstropyUserWarning


def sort_eq(list1, list2):
    return sorted(list1) == sorted(list2)


T1 = Table.read([' a b c d',
                 ' 2 c 7.0 0',
                 ' 2 b 5.0 1',
                 ' 2 b 6.0 2',
                 ' 2 a 4.0 3',
                 ' 0 a 0.0 4',
                 ' 1 b 3.0 5',
                 ' 1 a 2.0 6',
                 ' 1 a 1.0 7',
                 ], format='ascii')
T1.meta.update({'ta': 1})
T1['c'].meta.update({'a': 1})
T1['c'].description = 'column c'


def test_column_group_by():
    for masked in (False, True):
        t1 = Table(T1, masked=masked)
        t1a = t1['a'].copy()

        # Group by a Column (i.e. numpy array)
        t1ag = t1a.group_by(t1['a'])
        assert np.all(t1ag.groups.indices == np.array([0, 1, 4, 8]))

        # Group by a Table
        t1ag = t1a.group_by(t1['a', 'b'])
        assert np.all(t1ag.groups.indices == np.array([0, 1, 3, 4, 5, 7, 8]))

        # Group by a numpy structured array
        t1ag = t1a.group_by(t1['a', 'b'].as_array())
        assert np.all(t1ag.groups.indices == np.array([0, 1, 3, 4, 5, 7, 8]))

##TODO: parametrize other tests as well
@pytest.mark.parametrize("index", [True, False])
def test_table_group_by(index):
    """
    Test basic table group_by functionality for possible key types and for
    masked/unmasked tables.
    """
    for masked in (False, True):
        t1 = Table(T1, masked=masked)
        if index:
            t1.add_index('a')
        # Group by a single column key specified by name
        tg = t1.group_by('a')
        assert np.all(tg.groups.indices == np.array([0, 1, 4, 8]))
        assert str(tg.groups) == "<TableGroups indices=[0 1 4 8]>"
        assert str(tg['a'].groups) == "<ColumnGroups indices=[0 1 4 8]>"

        # Sorted by 'a' and in original order for rest
        assert tg.pformat() == [' a   b   c   d ',
                                '--- --- --- ---',
                                '  0   a 0.0   4',
                                '  1   b 3.0   5',
                                '  1   a 2.0   6',
                                '  1   a 1.0   7',
                                '  2   c 7.0   0',
                                '  2   b 5.0   1',
                                '  2   b 6.0   2',
                                '  2   a 4.0   3']
        assert tg.meta['ta'] == 1
        assert tg['c'].meta['a'] == 1
        assert tg['c'].description == 'column c'

        # Group by a table column
        tg2 = t1.group_by(t1['a'])
        assert tg.pformat() == tg2.pformat()

        # Group by two columns spec'd by name
        for keys in (['a', 'b'], ('a', 'b')):
            tg = t1.group_by(keys)
            assert np.all(tg.groups.indices == np.array([0, 1, 3, 4, 5, 7, 8]))
            # Sorted by 'a', 'b' and in original order for rest
            assert tg.pformat() == [' a   b   c   d ',
                                    '--- --- --- ---',
                                    '  0   a 0.0   4',
                                    '  1   a 2.0   6',
                                    '  1   a 1.0   7',
                                    '  1   b 3.0   5',
                                    '  2   a 4.0   3',
                                    '  2   b 5.0   1',
                                    '  2   b 6.0   2',
                                    '  2   c 7.0   0']

        # Group by a Table
        tg2 = t1.group_by(t1['a', 'b'])
        assert tg.pformat() == tg2.pformat()

        # Group by a structured array
        tg2 = t1.group_by(t1['a', 'b'].as_array())
        assert tg.pformat() == tg2.pformat()

        # Group by a simple ndarray
        tg = t1.group_by(np.array([0, 1, 0, 1, 2, 1, 0, 0]))
        assert np.all(tg.groups.indices == np.array([0, 4, 7, 8]))
        assert tg.pformat() == [' a   b   c   d ',
                                '--- --- --- ---',
                                '  2   c 7.0   0',
                                '  2   b 6.0   2',
                                '  1   a 2.0   6',
                                '  1   a 1.0   7',
                                '  2   b 5.0   1',
                                '  2   a 4.0   3',
                                '  1   b 3.0   5',
                                '  0   a 0.0   4']


def test_groups_keys():
    tg = T1.group_by('a')
    keys = tg.groups.keys
    assert keys.dtype.names == ('a',)
    assert np.all(keys['a'] == np.array([0, 1, 2]))

    tg = T1.group_by(['a', 'b'])
    keys = tg.groups.keys
    assert keys.dtype.names == ('a', 'b')
    assert np.all(keys['a'] == np.array([0, 1, 1, 2, 2, 2]))
    assert np.all(keys['b'] == np.array(['a', 'a', 'b', 'a', 'b', 'c']))

    # Grouping by Column ignores column name
    tg = T1.group_by(T1['b'])
    keys = tg.groups.keys
    assert keys.dtype.names is None


def test_groups_iterator():
    tg = T1.group_by('a')
    for ii, group in enumerate(tg.groups):
        assert group.pformat() == tg.groups[ii].pformat()
        assert group['a'][0] == tg['a'][tg.groups.indices[ii]]


def test_grouped_copy():
    """
    Test that copying a table or column copies the groups properly
    """
    for masked in (False, True):
        t1 = Table(T1, masked=masked)
        tg = t1.group_by('a')
        tgc = tg.copy()
        assert np.all(tgc.groups.indices == tg.groups.indices)
        assert np.all(tgc.groups.keys == tg.groups.keys)

        tac = tg['a'].copy()
        assert np.all(tac.groups.indices == tg['a'].groups.indices)

        c1 = t1['a'].copy()
        gc1 = c1.group_by(t1['a'])
        gc1c = gc1.copy()
        assert np.all(gc1c.groups.indices == np.array([0, 1, 4, 8]))


def test_grouped_slicing():
    """
    Test that slicing a table removes previous grouping
    """

    for masked in (False, True):
        t1 = Table(T1, masked=masked)

        # Regular slice of a table
        tg = t1.group_by('a')
        tg2 = tg[3:5]
        assert np.all(tg2.groups.indices == np.array([0, len(tg2)]))
        assert tg2.groups.keys is None


def test_group_column_from_table():
    """
    Group a column that is part of a table
    """
    cg = T1['c'].group_by(np.array(T1['a']))
    assert np.all(cg.groups.keys == np.array([0, 1, 2]))
    assert np.all(cg.groups.indices == np.array([0, 1, 4, 8]))


def test_table_groups_mask_index():
    """
    Use boolean mask as item in __getitem__ for groups
    """
    for masked in (False, True):
        t1 = Table(T1, masked=masked).group_by('a')

        t2 = t1.groups[np.array([True, False, True])]
        assert len(t2.groups) == 2
        assert t2.groups[0].pformat() == t1.groups[0].pformat()
        assert t2.groups[1].pformat() == t1.groups[2].pformat()
        assert np.all(t2.groups.keys['a'] == np.array([0, 2]))


def test_table_groups_array_index():
    """
    Use numpy array as item in __getitem__ for groups
    """
    for masked in (False, True):
        t1 = Table(T1, masked=masked).group_by('a')

        t2 = t1.groups[np.array([0, 2])]
        assert len(t2.groups) == 2
        assert t2.groups[0].pformat() == t1.groups[0].pformat()
        assert t2.groups[1].pformat() == t1.groups[2].pformat()
        assert np.all(t2.groups.keys['a'] == np.array([0, 2]))


def test_table_groups_slicing():
    """
    Test that slicing table groups works
    """

    for masked in (False, True):
        t1 = Table(T1, masked=masked).group_by('a')

        # slice(0, 2)
        t2 = t1.groups[0:2]
        assert len(t2.groups) == 2
        assert t2.groups[0].pformat() == t1.groups[0].pformat()
        assert t2.groups[1].pformat() == t1.groups[1].pformat()
        assert np.all(t2.groups.keys['a'] == np.array([0, 1]))

        # slice(1, 2)
        t2 = t1.groups[1:2]
        assert len(t2.groups) == 1
        assert t2.groups[0].pformat() == t1.groups[1].pformat()
        assert np.all(t2.groups.keys['a'] == np.array([1]))

        # slice(0, 3, 2)
        t2 = t1.groups[0:3:2]
        assert len(t2.groups) == 2
        assert t2.groups[0].pformat() == t1.groups[0].pformat()
        assert t2.groups[1].pformat() == t1.groups[2].pformat()
        assert np.all(t2.groups.keys['a'] == np.array([0, 2]))


def test_grouped_item_access():
    """
    Test that column slicing preserves grouping
    """
    for masked in (False, True):
        t1 = Table(T1, masked=masked)

        # Regular slice of a table
        tg = t1.group_by('a')
        tgs = tg['a', 'c', 'd']
        assert np.all(tgs.groups.keys == tg.groups.keys)
        assert np.all(tgs.groups.indices == tg.groups.indices)
        tgsa = tgs.groups.aggregate(np.sum)
        assert tgsa.pformat() == [' a   c    d ',
                                  '--- ---- ---',
                                  '  0  0.0   4',
                                  '  1  6.0  18',
                                  '  2 22.0   6']

        tgs = tg['c', 'd']
        assert np.all(tgs.groups.keys == tg.groups.keys)
        assert np.all(tgs.groups.indices == tg.groups.indices)
        tgsa = tgs.groups.aggregate(np.sum)
        assert tgsa.pformat() == [' c    d ',
                                  '---- ---',
                                  ' 0.0   4',
                                  ' 6.0  18',
                                  '22.0   6']


def test_mutable_operations():
    """
    Operations like adding or deleting a row should removing grouping,
    but adding or removing or renaming a column should retain grouping.
    """
    for masked in (False, True):
        t1 = Table(T1, masked=masked)

        # add row
        tg = t1.group_by('a')
        tg.add_row((0, 'a', 3.0, 4))
        assert np.all(tg.groups.indices == np.array([0, len(tg)]))
        assert tg.groups.keys is None

        # remove row
        tg = t1.group_by('a')
        tg.remove_row(4)
        assert np.all(tg.groups.indices == np.array([0, len(tg)]))
        assert tg.groups.keys is None

        # add column
        tg = t1.group_by('a')
        indices = tg.groups.indices.copy()
        tg.add_column(Column(name='e', data=np.arange(len(tg))))
        assert np.all(tg.groups.indices == indices)
        assert np.all(tg['e'].groups.indices == indices)
        assert np.all(tg['e'].groups.keys == tg.groups.keys)

        # remove column (not key column)
        tg = t1.group_by('a')
        tg.remove_column('b')
        assert np.all(tg.groups.indices == indices)
        # Still has original key col names
        assert tg.groups.keys.dtype.names == ('a',)
        assert np.all(tg['a'].groups.indices == indices)

        # remove key column
        tg = t1.group_by('a')
        tg.remove_column('a')
        assert np.all(tg.groups.indices == indices)
        assert tg.groups.keys.dtype.names == ('a',)
        assert np.all(tg['b'].groups.indices == indices)

        # rename key column
        tg = t1.group_by('a')
        tg.rename_column('a', 'aa')
        assert np.all(tg.groups.indices == indices)
        assert tg.groups.keys.dtype.names == ('a',)
        assert np.all(tg['aa'].groups.indices == indices)


def test_group_by_masked():
    t1m = Table(T1, masked=True)
    t1m['c'].mask[4] = True
    t1m['d'].mask[5] = True
    assert t1m.group_by('a').pformat() == [' a   b   c   d ',
                                           '--- --- --- ---',
                                           '  0   a  --   4',
                                           '  1   b 3.0  --',
                                           '  1   a 2.0   6',
                                           '  1   a 1.0   7',
                                           '  2   c 7.0   0',
                                           '  2   b 5.0   1',
                                           '  2   b 6.0   2',
                                           '  2   a 4.0   3']


def test_group_by_errors():
    """
    Appropriate errors get raised.
    """
    # Bad column name as string
    with pytest.raises(ValueError):
        T1.group_by('f')

    # Bad column names in list
    with pytest.raises(ValueError):
        T1.group_by(['f', 'g'])

    # Wrong length array
    with pytest.raises(ValueError):
        T1.group_by(np.array([1, 2]))

    # Wrong type
    with pytest.raises(TypeError):
        T1.group_by(None)

    # Masked key column
    t1 = Table(T1, masked=True)
    t1['a'].mask[4] = True
    with pytest.raises(ValueError):
        t1.group_by('a')


def test_groups_keys_meta():
    """
    Make sure the keys meta['grouped_by_table_cols'] is working.
    """
    # Group by column in this table
    tg = T1.group_by('a')
    assert tg.groups.keys.meta['grouped_by_table_cols'] is True
    assert tg['c'].groups.keys.meta['grouped_by_table_cols'] is True
    assert tg.groups[1].groups.keys.meta['grouped_by_table_cols'] is True
    assert (tg['d'].groups[np.array([False, True, True])]
            .groups.keys.meta['grouped_by_table_cols'] is True)

    # Group by external Table
    tg = T1.group_by(T1['a', 'b'])
    assert tg.groups.keys.meta['grouped_by_table_cols'] is False
    assert tg['c'].groups.keys.meta['grouped_by_table_cols'] is False
    assert tg.groups[1].groups.keys.meta['grouped_by_table_cols'] is False

    # Group by external numpy array
    tg = T1.group_by(T1['a', 'b'].as_array())
    assert not hasattr(tg.groups.keys, 'meta')
    assert not hasattr(tg['c'].groups.keys, 'meta')

    # Group by Column
    tg = T1.group_by(T1['a'])
    assert 'grouped_by_table_cols' not in tg.groups.keys.meta
    assert 'grouped_by_table_cols' not in tg['c'].groups.keys.meta


def test_table_aggregate():
    """
    Aggregate a table
    """
    # Table with only summable cols
    t1 = T1['a', 'c', 'd']
    tg = t1.group_by('a')
    tga = tg.groups.aggregate(np.sum)
    assert tga.pformat() == [' a   c    d ',
                             '--- ---- ---',
                             '  0  0.0   4',
                             '  1  6.0  18',
                             '  2 22.0   6']
    # Reverts to default groups
    assert np.all(tga.groups.indices == np.array([0, 3]))
    assert tga.groups.keys is None

    # metadata survives
    assert tga.meta['ta'] == 1
    assert tga['c'].meta['a'] == 1
    assert tga['c'].description == 'column c'

    # Aggregate with np.sum with masked elements.  This results
    # in one group with no elements, hence a nan result and conversion
    # to float for the 'd' column.
    t1m = Table(t1, masked=True)
    t1m['c'].mask[4:6] = True
    t1m['d'].mask[4:6] = True
    tg = t1m.group_by('a')
    with catch_warnings(Warning) as warning_lines:
        tga = tg.groups.aggregate(np.sum)
        assert warning_lines[0].category == UserWarning
        assert "converting a masked element to nan" in str(warning_lines[0].message)

    assert tga.pformat() == [' a   c    d  ',
                             '--- ---- ----',
                             '  0  nan  nan',
                             '  1  3.0 13.0',
                             '  2 22.0  6.0']

    # Aggregrate with np.sum with masked elements, but where every
    # group has at least one remaining (unmasked) element.  Then
    # the int column stays as an int.
    t1m = Table(t1, masked=True)
    t1m['c'].mask[5] = True
    t1m['d'].mask[5] = True
    tg = t1m.group_by('a')
    tga = tg.groups.aggregate(np.sum)
    assert tga.pformat() == [' a   c    d ',
                             '--- ---- ---',
                             '  0  0.0   4',
                             '  1  3.0  13',
                             '  2 22.0   6']

    # Aggregate with a column type that cannot by supplied to the aggregating
    # function.  This raises a warning but still works.
    tg = T1.group_by('a')
    with catch_warnings(Warning) as warning_lines:
        tga = tg.groups.aggregate(np.sum)
        assert warning_lines[0].category == AstropyUserWarning
        assert "Cannot aggregate column" in str(warning_lines[0].message)
    assert tga.pformat() == [' a   c    d ',
                             '--- ---- ---',
                             '  0  0.0   4',
                             '  1  6.0  18',
                             '  2 22.0   6']


def test_table_aggregate_reduceat():
    """
    Aggregate table with functions which have a reduceat method
    """
    # Comparison functions without reduceat
    def np_mean(x):
        return np.mean(x)
    def np_sum(x):
        return np.sum(x)
    def np_add(x):
        return np.add(x)

    # Table with only summable cols
    t1 = T1['a', 'c', 'd']
    tg = t1.group_by('a')
    # Comparison
    tga_r = tg.groups.aggregate(np.sum)
    tga_a = tg.groups.aggregate(np.add)
    tga_n = tg.groups.aggregate(np_sum)

    assert np.all(tga_r == tga_n)
    assert np.all(tga_a == tga_n)
    assert tga_n.pformat() == [' a   c    d ',
                               '--- ---- ---',
                               '  0  0.0   4',
                               '  1  6.0  18',
                               '  2 22.0   6']

    tga_r = tg.groups.aggregate(np.mean)
    tga_n = tg.groups.aggregate(np_mean)
    assert np.all(tga_r == tga_n)
    assert tga_n.pformat() == [' a   c   d ',
                               '--- --- ---',
                               '  0 0.0 4.0',
                               '  1 2.0 6.0',
                               '  2 5.5 1.5']

    # Binary ufunc np_add should raise warning without reduceat
    t2 = T1['a', 'c']
    tg = t2.group_by('a')

    with catch_warnings(Warning) as warning_lines:
        tga = tg.groups.aggregate(np_add)
        assert warning_lines[0].category == AstropyUserWarning
        assert "Cannot aggregate column" in str(warning_lines[0].message)
    assert tga.pformat() == [' a ',
                             '---',
                             '  0',
                             '  1',
                             '  2']


def test_column_aggregate():
    """
    Aggregate a single table column
    """
    for masked in (False, True):
        tg = Table(T1, masked=masked).group_by('a')
        tga = tg['c'].groups.aggregate(np.sum)
        assert tga.pformat() == [' c  ',
                                 '----',
                                 ' 0.0',
                                 ' 6.0',
                                 '22.0']


def test_table_filter():
    """
    Table groups filtering
    """
    def all_positive(table, key_colnames):
        colnames = [name for name in table.colnames if name not in key_colnames]
        for colname in colnames:
            if np.any(table[colname] < 0):
                return False
        return True

    # Negative value in 'a' column should not filter because it is a key col
    t = Table.read([' a c d',
                    ' -2 7.0 0',
                    ' -2 5.0 1',
                    ' 0 0.0 4',
                    ' 1 3.0 5',
                    ' 1 2.0 -6',
                    ' 1 1.0 7',
                    ' 3 3.0 5',
                    ' 3 -2.0 6',
                    ' 3 1.0 7',
                    ], format='ascii')
    tg = t.group_by('a')
    t2 = tg.groups.filter(all_positive)
    assert t2.groups[0].pformat() == [' a   c   d ',
                                      '--- --- ---',
                                      ' -2 7.0   0',
                                      ' -2 5.0   1']
    assert t2.groups[1].pformat() == [' a   c   d ',
                                      '--- --- ---',
                                      '  0 0.0   4']


def test_column_filter():
    """
    Table groups filtering
    """
    def all_positive(column):
        if np.any(column < 0):
            return False
        return True

    # Negative value in 'a' column should not filter because it is a key col
    t = Table.read([' a c d',
                    ' -2 7.0 0',
                    ' -2 5.0 1',
                    ' 0 0.0 4',
                    ' 1 3.0 5',
                    ' 1 2.0 -6',
                    ' 1 1.0 7',
                    ' 3 3.0 5',
                    ' 3 -2.0 6',
                    ' 3 1.0 7',
                    ], format='ascii')
    tg = t.group_by('a')
    c2 = tg['c'].groups.filter(all_positive)
    assert len(c2.groups) == 3
    assert c2.groups[0].pformat() == [' c ', '---', '7.0', '5.0']
    assert c2.groups[1].pformat() == [' c ', '---', '0.0']
    assert c2.groups[2].pformat() == [' c ', '---', '3.0', '2.0', '1.0']
