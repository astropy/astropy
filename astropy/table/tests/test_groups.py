# Licensed under a 3-clause BSD style license - see LICENSE.rst
import numpy as np

from ...tests.helper import pytest, catch_warnings
from ...table import Table, Column


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
        t1ag = t1a.group_by(t1['a', 'b']._data)
        assert np.all(t1ag.groups.indices == np.array([0, 1, 3, 4, 5, 7, 8]))


def test_table_group_by():
    """
    Test basic table group_by functionality for possible key types and for
    masked/unmasked tables.
    """

    for masked in (False, True):
        t1 = Table(T1, masked=masked)

        # Group by a single column key specified by name
        tg = t1.group_by('a')
        assert np.all(tg.groups.indices == np.array([0, 1, 4, 8]))
        assert str(tg.groups) == "<TableGroups group_keys=('a',) indices=[0 1 4 8]>"
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
        tg2 = t1.group_by(t1['a', 'b']._data)
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
    keys = list(tg.groups.keys())
    assert keys == [0, 1, 2]

    tg = T1.group_by(['a', 'b'])
    keys = list(tg.groups.keys())
    assert keys == [(0, 'a'), (1, 'a'), (1, 'b'), (2, 'a'), (2, 'b'), (2, 'c')]

    tg = T1.group_by(T1['b'])
    keys = list(tg.groups.keys())
    assert keys == [0, 1, 2]


def test_groups_values():
    tg = T1.group_by('a')
    values = list(tg.groups.values())
    assert values[0].pformat() == [' a   b   c   d ',
                                   '--- --- --- ---',
                                   '  0   a 0.0   4']
    assert values[1].pformat() == [' a   b   c   d ',
                                   '--- --- --- ---',
                                   '  1   b 3.0   5',
                                   '  1   a 2.0   6',
                                   '  1   a 1.0   7',
                                   ]
    assert values[2].pformat() == [' a   b   c   d ',
                                   '--- --- --- ---',
                                   '  2   c 7.0   0',
                                   '  2   b 5.0   1',
                                   '  2   b 6.0   2',
                                   '  2   a 4.0   3',
                                   ]


def test_grouped_copy():
    """
    Test that copying a table or column copies the groups properly
    """
    for masked in (False, True):
        t1 = Table(T1, masked=masked)
        tg = t1.group_by('a')
        tgc = tg.copy()
        assert np.all(tgc.groups.indices == tg.groups.indices)
        assert tgc.groups.group_keys == tg.groups.group_keys

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
        assert tg2.groups.group_keys == ()


def test_grouped_item_access():
    """
    Test that column slicing preserves grouping
    """
    for masked in (False, True):
        t1 = Table(T1, masked=masked)

        # Regular slice of a table
        tg = t1.group_by('a')
        tgs = tg['a', 'c', 'd']
        assert tgs.groups.group_keys == ('a',)
        assert np.all(tgs.groups.indices == np.array([0, 1, 4, 8]))
        tgsa = tgs.groups.aggregate(np.sum)
        assert tgsa.pformat() == [' a   c    d ',
                                  '--- ---- ---',
                                  '  0  0.0   4',
                                  '  1  6.0  18',
                                  '  2 22.0   6']

        tgs = tg['c', 'd']
        assert tgs.groups.group_keys == ()
        assert np.all(tgs.groups.indices == np.array([0, 1, 4, 8]))
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

        # remove row
        tg = t1.group_by('a')
        tg.remove_row(4)
        assert np.all(tg.groups.indices == np.array([0, len(tg)]))

        # add column
        tg = t1.group_by('a')
        indices = tg.groups.indices.copy()
        tg.add_column(Column(name='e', data=np.arange(len(tg))))
        assert np.all(tg.groups.indices == indices)
        assert np.all(tg['e'].groups.indices == indices)

        # remove column (not key column)
        tg = t1.group_by('a')
        tg.remove_column('b')
        assert np.all(tg.groups.indices == indices)
        assert tg.groups.group_keys == ('a',)
        assert np.all(tg['a'].groups.indices == indices)

        # remove key column (removes key from group_keys)
        tg = t1.group_by('a')
        tg.remove_column('a')
        assert np.all(tg.groups.indices == indices)
        assert tg.groups.group_keys == ()
        assert np.all(tg['b'].groups.indices == indices)

        # remove key column (removes key from group_keys)
        tg = t1.group_by('a')
        tg.rename_column('a', 'aa')
        assert np.all(tg.groups.indices == indices)
        assert tg.groups.group_keys == ('aa',)
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
    assert tga.groups.group_keys == ()

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
        assert warning_lines[0].category == UserWarning
        assert "Cannot aggregate column" in str(warning_lines[0].message)
    assert tga.pformat() == [' a   c    d ',
                             '--- ---- ---',
                             '  0  0.0   4',
                             '  1  6.0  18',
                             '  2 22.0   6']


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
