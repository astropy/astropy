"""
High-level table operations:

- join()
- setdiff()
- hstack()
- vstack()
- dstack()
"""
# Licensed under a 3-clause BSD style license - see LICENSE.rst

from copy import deepcopy
import collections
import itertools
from collections import OrderedDict, Counter
from collections.abc import Mapping, Sequence

import numpy as np

from astropy.utils import metadata
from astropy.utils.masked import Masked
from .table import Table, QTable, Row, Column, MaskedColumn
from astropy.units import Quantity

from . import _np_utils
from .np_utils import TableMergeError

__all__ = ['join', 'setdiff', 'hstack', 'vstack', 'unique',
           'join_skycoord', 'join_distance']

__doctest_requires__ = {'join_skycoord': ['scipy'], 'join_distance': ['scipy']}


def _merge_table_meta(out, tables, metadata_conflicts='warn'):
    out_meta = deepcopy(tables[0].meta)
    for table in tables[1:]:
        out_meta = metadata.merge(out_meta, table.meta, metadata_conflicts=metadata_conflicts)
    out.meta.update(out_meta)


def _get_list_of_tables(tables):
    """
    Check that tables is a Table or sequence of Tables.  Returns the
    corresponding list of Tables.
    """

    # Make sure we have a list of things
    if not isinstance(tables, Sequence):
        tables = [tables]

    # Make sure there is something to stack
    if len(tables) == 0:
        raise ValueError('no values provided to stack.')

    # Convert inputs (Table, Row, or anything column-like) to Tables.
    # Special case that Quantity converts to a QTable.
    for ii, val in enumerate(tables):
        if isinstance(val, Table):
            pass
        elif isinstance(val, Row):
            tables[ii] = Table(val)
        elif isinstance(val, Quantity):
            tables[ii] = QTable([val])
        else:
            try:
                tables[ii] = Table([val])
            except (ValueError, TypeError) as err:
                raise TypeError(f'Cannot convert {val} to table column.') from err

    return tables


def _get_out_class(objs):
    """
    From a list of input objects ``objs`` get merged output object class.

    This is just taken as the deepest subclass. This doesn't handle complicated
    inheritance schemes, but as a special case, classes which share ``info``
    are taken to be compatible.
    """
    out_class = objs[0].__class__
    for obj in objs[1:]:
        if issubclass(obj.__class__, out_class):
            out_class = obj.__class__

    if any(not (issubclass(out_class, obj.__class__)
                or out_class.info is obj.__class__.info) for obj in objs):
        raise ValueError('unmergeable object classes {}'
                         .format([obj.__class__.__name__ for obj in objs]))

    return out_class


def join_skycoord(distance, distance_func='search_around_sky'):
    """Helper function to join on SkyCoord columns using distance matching.

    This function is intended for use in ``table.join()`` to allow performing a
    table join where the key columns are both ``SkyCoord`` objects, matched by
    computing the distance between points and accepting values below
    ``distance``.

    The distance cross-matching is done using either
    `~astropy.coordinates.search_around_sky` or
    `~astropy.coordinates.search_around_3d`, depending on the value of
    ``distance_func``.  The default is ``'search_around_sky'``.

    One can also provide a function object for ``distance_func``, in which case
    it must be a function that follows the same input and output API as
    `~astropy.coordinates.search_around_sky`. In this case the function will
    be called with ``(skycoord1, skycoord2, distance)`` as arguments.

    Parameters
    ----------
    distance : `~astropy.units.Quantity` ['angle', 'length']
        Maximum distance between points to be considered a join match.
        Must have angular or distance units.
    distance_func : str or function
        Specifies the function for performing the cross-match based on
        ``distance``. If supplied as a string this specifies the name of a
        function in `astropy.coordinates`. If supplied as a function then that
        function is called directly.

    Returns
    -------
    join_func : function
        Function that accepts two ``SkyCoord`` columns (col1, col2) and returns
        the tuple (ids1, ids2) of pair-matched unique identifiers.

    Examples
    --------
    This example shows an inner join of two ``SkyCoord`` columns, taking any
    sources within 0.2 deg to be a match.  Note the new ``sc_id`` column which
    is added and provides a unique source identifier for the matches.

      >>> from astropy.coordinates import SkyCoord
      >>> import astropy.units as u
      >>> from astropy.table import Table, join_skycoord
      >>> from astropy import table

      >>> sc1 = SkyCoord([0, 1, 1.1, 2], [0, 0, 0, 0], unit='deg')
      >>> sc2 = SkyCoord([0.5, 1.05, 2.1], [0, 0, 0], unit='deg')

      >>> join_func = join_skycoord(0.2 * u.deg)
      >>> join_func(sc1, sc2)  # Associate each coordinate with unique source ID
      (array([3, 1, 1, 2]), array([4, 1, 2]))

      >>> t1 = Table([sc1], names=['sc'])
      >>> t2 = Table([sc2], names=['sc'])
      >>> t12 = table.join(t1, t2, join_funcs={'sc': join_skycoord(0.2 * u.deg)})
      >>> print(t12)  # Note new `sc_id` column with the IDs from join_func()
      sc_id   sc_1    sc_2
            deg,deg deg,deg
      ----- ------- --------
          1 1.0,0.0 1.05,0.0
          1 1.1,0.0 1.05,0.0
          2 2.0,0.0  2.1,0.0

    """
    if isinstance(distance_func, str):
        import astropy.coordinates as coords
        try:
            distance_func = getattr(coords, distance_func)
        except AttributeError as err:
            raise ValueError('distance_func must be a function in astropy.coordinates') from err
    else:
        from inspect import isfunction
        if not isfunction(distance_func):
            raise ValueError('distance_func must be a str or function')

    def join_func(sc1, sc2):

        # Call the appropriate SkyCoord method to find pairs within distance
        idxs1, idxs2, d2d, d3d = distance_func(sc1, sc2, distance)

        # Now convert that into unique identifiers for each near-pair. This is
        # taken to be transitive, so that if points 1 and 2 are "near" and points
        # 1 and 3 are "near", then 1, 2, and 3 are all given the same identifier.
        # This identifier will then be used in the table join matching.

        # Identifiers for each column, initialized to all zero.
        ids1 = np.zeros(len(sc1), dtype=int)
        ids2 = np.zeros(len(sc2), dtype=int)

        # Start the identifier count at 1
        id_ = 1
        for idx1, idx2 in zip(idxs1, idxs2):
            # If this col1 point is previously identified then set corresponding
            # col2 point to same identifier.  Likewise for col2 and col1.
            if ids1[idx1] > 0:
                ids2[idx2] = ids1[idx1]
            elif ids2[idx2] > 0:
                ids1[idx1] = ids2[idx2]
            else:
                # Not yet seen so set identifier for col1 and col2
                ids1[idx1] = id_
                ids2[idx2] = id_
                id_ += 1

        # Fill in unique identifiers for points with no near neighbor
        for ids in (ids1, ids2):
            for idx in np.flatnonzero(ids == 0):
                ids[idx] = id_
                id_ += 1

        # End of enclosure join_func()
        return ids1, ids2

    return join_func


def join_distance(distance, kdtree_args=None, query_args=None):
    """Helper function to join table columns using distance matching.

    This function is intended for use in ``table.join()`` to allow performing
    a table join where the key columns are matched by computing the distance
    between points and accepting values below ``distance``. This numerical
    "fuzzy" match can apply to 1-D or 2-D columns, where in the latter case
    the distance is a vector distance.

    The distance cross-matching is done using `scipy.spatial.cKDTree`. If
    necessary you can tweak the default behavior by providing ``dict`` values
    for the ``kdtree_args`` or ``query_args``.

    Parameters
    ----------
    distance : float or `~astropy.units.Quantity` ['length']
        Maximum distance between points to be considered a join match
    kdtree_args : dict, None
        Optional extra args for `~scipy.spatial.cKDTree`
    query_args : dict, None
        Optional extra args for `~scipy.spatial.cKDTree.query_ball_tree`

    Returns
    -------
    join_func : function
        Function that accepts (skycoord1, skycoord2) and returns the tuple
        (ids1, ids2) of pair-matched unique identifiers.

    Examples
    --------

      >>> from astropy.table import Table, join_distance
      >>> from astropy import table

      >>> c1 = [0, 1, 1.1, 2]
      >>> c2 = [0.5, 1.05, 2.1]

      >>> t1 = Table([c1], names=['col'])
      >>> t2 = Table([c2], names=['col'])
      >>> t12 = table.join(t1, t2, join_type='outer', join_funcs={'col': join_distance(0.2)})
      >>> print(t12)
      col_id col_1 col_2
      ------ ----- -----
           1   1.0  1.05
           1   1.1  1.05
           2   2.0   2.1
           3   0.0    --
           4    --   0.5

    """
    try:
        from scipy.spatial import cKDTree
    except ImportError as exc:
        raise ImportError('scipy is required to use join_distance()') from exc

    if kdtree_args is None:
        kdtree_args = {}
    if query_args is None:
        query_args = {}

    def join_func(col1, col2):
        if col1.ndim > 2 or col2.ndim > 2:
            raise ValueError('columns for isclose_join must be 1- or 2-dimensional')

        if isinstance(distance, Quantity):
            # Convert to np.array with common unit
            col1 = col1.to_value(distance.unit)
            col2 = col2.to_value(distance.unit)
            dist = distance.value
        else:
            # Convert to np.array to allow later in-place shape changing
            col1 = np.asarray(col1)
            col2 = np.asarray(col2)
            dist = distance

        # Ensure columns are pure np.array and are 2-D for use with KDTree
        if col1.ndim == 1:
            col1.shape = col1.shape + (1,)
        if col2.ndim == 1:
            col2.shape = col2.shape + (1,)

        # Cross-match col1 and col2 within dist using KDTree
        kd1 = cKDTree(col1, **kdtree_args)
        kd2 = cKDTree(col2, **kdtree_args)
        nears = kd1.query_ball_tree(kd2, r=dist, **query_args)

        # Output of above is nears which is a list of lists, where the outer
        # list corresponds to each item in col1, and where the inner lists are
        # indexes into col2 of elements within the distance tolerance.  This
        # identifies col1 / col2 near pairs.

        # Now convert that into unique identifiers for each near-pair. This is
        # taken to be transitive, so that if points 1 and 2 are "near" and points
        # 1 and 3 are "near", then 1, 2, and 3 are all given the same identifier.
        # This identifier will then be used in the table join matching.

        # Identifiers for each column, initialized to all zero.
        ids1 = np.zeros(len(col1), dtype=int)
        ids2 = np.zeros(len(col2), dtype=int)

        # Start the identifier count at 1
        id_ = 1
        for idx1, idxs2 in enumerate(nears):
            for idx2 in idxs2:
                # If this col1 point is previously identified then set corresponding
                # col2 point to same identifier.  Likewise for col2 and col1.
                if ids1[idx1] > 0:
                    ids2[idx2] = ids1[idx1]
                elif ids2[idx2] > 0:
                    ids1[idx1] = ids2[idx2]
                else:
                    # Not yet seen so set identifier for col1 and col2
                    ids1[idx1] = id_
                    ids2[idx2] = id_
                    id_ += 1

        # Fill in unique identifiers for points with no near neighbor
        for ids in (ids1, ids2):
            for idx in np.flatnonzero(ids == 0):
                ids[idx] = id_
                id_ += 1

        # End of enclosure join_func()
        return ids1, ids2

    return join_func


def join(left, right, keys=None, join_type='inner', *,
         keys_left=None, keys_right=None,
         uniq_col_name='{col_name}_{table_name}',
         table_names=['1', '2'], metadata_conflicts='warn',
         join_funcs=None):
    """
    Perform a join of the left table with the right table on specified keys.

    Parameters
    ----------
    left : `~astropy.table.Table`-like object
        Left side table in the join. If not a Table, will call ``Table(left)``
    right : `~astropy.table.Table`-like object
        Right side table in the join. If not a Table, will call ``Table(right)``
    keys : str or list of str
        Name(s) of column(s) used to match rows of left and right tables.
        Default is to use all columns which are common to both tables.
    join_type : str
        Join type ('inner' | 'outer' | 'left' | 'right' | 'cartesian'), default is 'inner'
    keys_left : str or list of str or list of column-like, optional
        Left column(s) used to match rows instead of ``keys`` arg. This can be
        be a single left table column name or list of column names, or a list of
        column-like values with the same lengths as the left table.
    keys_right : str or list of str or list of column-like, optional
        Same as ``keys_left``, but for the right side of the join.
    uniq_col_name : str or None
        String generate a unique output column name in case of a conflict.
        The default is '{col_name}_{table_name}'.
    table_names : list of str or None
        Two-element list of table names used when generating unique output
        column names.  The default is ['1', '2'].
    metadata_conflicts : str
        How to proceed with metadata conflicts. This should be one of:
            * ``'silent'``: silently pick the last conflicting meta-data value
            * ``'warn'``: pick the last conflicting meta-data value, but emit a warning (default)
            * ``'error'``: raise an exception.
    join_funcs : dict, None
        Dict of functions to use for matching the corresponding key column(s).
        See `~astropy.table.join_skycoord` for an example and details.

    Returns
    -------
    joined_table : `~astropy.table.Table` object
        New table containing the result of the join operation.
    """

    # Try converting inputs to Table as needed
    if not isinstance(left, Table):
        left = Table(left)
    if not isinstance(right, Table):
        right = Table(right)

    col_name_map = OrderedDict()
    out = _join(left, right, keys, join_type,
                uniq_col_name, table_names, col_name_map, metadata_conflicts,
                join_funcs,
                keys_left=keys_left, keys_right=keys_right)

    # Merge the column and table meta data. Table subclasses might override
    # these methods for custom merge behavior.
    _merge_table_meta(out, [left, right], metadata_conflicts=metadata_conflicts)

    return out


def setdiff(table1, table2, keys=None):
    """
    Take a set difference of table rows.

    The row set difference will contain all rows in ``table1`` that are not
    present in ``table2``. If the keys parameter is not defined, all columns in
    ``table1`` will be included in the output table.

    Parameters
    ----------
    table1 : `~astropy.table.Table`
        ``table1`` is on the left side of the set difference.
    table2 : `~astropy.table.Table`
        ``table2`` is on the right side of the set difference.
    keys : str or list of str
        Name(s) of column(s) used to match rows of left and right tables.
        Default is to use all columns in ``table1``.

    Returns
    -------
    diff_table : `~astropy.table.Table`
        New table containing the set difference between tables. If the set
        difference is none, an empty table will be returned.

    Examples
    --------
    To get a set difference between two tables::

      >>> from astropy.table import setdiff, Table
      >>> t1 = Table({'a': [1, 4, 9], 'b': ['c', 'd', 'f']}, names=('a', 'b'))
      >>> t2 = Table({'a': [1, 5, 9], 'b': ['c', 'b', 'f']}, names=('a', 'b'))
      >>> print(t1)
       a   b
      --- ---
        1   c
        4   d
        9   f
      >>> print(t2)
       a   b
      --- ---
        1   c
        5   b
        9   f
      >>> print(setdiff(t1, t2))
       a   b
      --- ---
        4   d

      >>> print(setdiff(t2, t1))
       a   b
      --- ---
        5   b
    """
    if keys is None:
        keys = table1.colnames

    # Check that all keys are in table1 and table2
    for tbl, tbl_str in ((table1, 'table1'), (table2, 'table2')):
        diff_keys = np.setdiff1d(keys, tbl.colnames)
        if len(diff_keys) != 0:
            raise ValueError("The {} columns are missing from {}, cannot take "
                             "a set difference.".format(diff_keys, tbl_str))

    # Make a light internal copy of both tables
    t1 = table1.copy(copy_data=False)
    t1.meta = {}
    t1.keep_columns(keys)
    t1['__index1__'] = np.arange(len(table1))  # Keep track of rows indices

    # Make a light internal copy to avoid touching table2
    t2 = table2.copy(copy_data=False)
    t2.meta = {}
    t2.keep_columns(keys)
    # Dummy column to recover rows after join
    t2['__index2__'] = np.zeros(len(t2), dtype=np.uint8)  # dummy column

    t12 = _join(t1, t2, join_type='left', keys=keys,
                metadata_conflicts='silent')

    # If t12 index2 is masked then that means some rows were in table1 but not table2.
    if hasattr(t12['__index2__'], 'mask'):
        # Define bool mask of table1 rows not in table2
        diff = t12['__index2__'].mask
        # Get the row indices of table1 for those rows
        idx = t12['__index1__'][diff]
        # Select corresponding table1 rows straight from table1 to ensure
        # correct table and column types.
        t12_diff = table1[idx]
    else:
        t12_diff = table1[[]]

    return t12_diff


def dstack(tables, join_type='outer', metadata_conflicts='warn'):
    """
    Stack columns within tables depth-wise

    A ``join_type`` of 'exact' means that the tables must all have exactly
    the same column names (though the order can vary).  If ``join_type``
    is 'inner' then the intersection of common columns will be the output.
    A value of 'outer' (default) means the output will have the union of
    all columns, with table values being masked where no common values are
    available.

    Parameters
    ----------
    tables : `~astropy.table.Table` or `~astropy.table.Row` or list thereof
        Table(s) to stack along depth-wise with the current table
        Table columns should have same shape and name for depth-wise stacking
    join_type : str
        Join type ('inner' | 'exact' | 'outer'), default is 'outer'
    metadata_conflicts : str
        How to proceed with metadata conflicts. This should be one of:
            * ``'silent'``: silently pick the last conflicting meta-data value
            * ``'warn'``: pick the last conflicting meta-data value, but emit a warning (default)
            * ``'error'``: raise an exception.

    Returns
    -------
    stacked_table : `~astropy.table.Table` object
        New table containing the stacked data from the input tables.

    Examples
    --------
    To stack two tables along rows do::

      >>> from astropy.table import vstack, Table
      >>> t1 = Table({'a': [1, 2], 'b': [3, 4]}, names=('a', 'b'))
      >>> t2 = Table({'a': [5, 6], 'b': [7, 8]}, names=('a', 'b'))
      >>> print(t1)
       a   b
      --- ---
        1   3
        2   4
      >>> print(t2)
       a   b
      --- ---
        5   7
        6   8
      >>> print(dstack([t1, t2]))
         a        b
      int64[2] int64[2]
      -------- --------
        1 .. 5   3 .. 7
        2 .. 6   4 .. 8
    """
    _check_join_type(join_type, 'dstack')

    tables = _get_list_of_tables(tables)
    if len(tables) == 1:
        return tables[0]  # no point in stacking a single table

    n_rows = set(len(table) for table in tables)
    if len(n_rows) != 1:
        raise ValueError('Table lengths must all match for dstack')
    n_row = n_rows.pop()

    out = vstack(tables, join_type, metadata_conflicts)

    for name, col in out.columns.items():
        col = out[name]

        # Reshape to so each original column is now in a row.
        # If entries are not 0-dim then those additional shape dims
        # are just carried along.
        # [x x x y y y] => [[x x x],
        #                   [y y y]]
        new_shape = (len(tables), n_row) + col.shape[1:]
        try:
            col.shape = (len(tables), n_row) + col.shape[1:]
        except AttributeError:
            col = col.reshape(new_shape)

        # Transpose the table and row axes to get to
        # [[x, y],
        #  [x, y]
        #  [x, y]]
        axes = np.arange(len(col.shape))
        axes[:2] = [1, 0]

        # This temporarily makes `out` be corrupted (columns of different
        # length) but it all works out in the end.
        out.columns.__setitem__(name, col.transpose(axes), validated=True)

    return out


def vstack(tables, join_type='outer', metadata_conflicts='warn'):
    """
    Stack tables vertically (along rows)

    A ``join_type`` of 'exact' means that the tables must all have exactly
    the same column names (though the order can vary).  If ``join_type``
    is 'inner' then the intersection of common columns will be the output.
    A value of 'outer' (default) means the output will have the union of
    all columns, with table values being masked where no common values are
    available.

    Parameters
    ----------
    tables : `~astropy.table.Table` or `~astropy.table.Row` or list thereof
        Table(s) to stack along rows (vertically) with the current table
    join_type : str
        Join type ('inner' | 'exact' | 'outer'), default is 'outer'
    metadata_conflicts : str
        How to proceed with metadata conflicts. This should be one of:
            * ``'silent'``: silently pick the last conflicting meta-data value
            * ``'warn'``: pick the last conflicting meta-data value, but emit a warning (default)
            * ``'error'``: raise an exception.

    Returns
    -------
    stacked_table : `~astropy.table.Table` object
        New table containing the stacked data from the input tables.

    Examples
    --------
    To stack two tables along rows do::

      >>> from astropy.table import vstack, Table
      >>> t1 = Table({'a': [1, 2], 'b': [3, 4]}, names=('a', 'b'))
      >>> t2 = Table({'a': [5, 6], 'b': [7, 8]}, names=('a', 'b'))
      >>> print(t1)
       a   b
      --- ---
        1   3
        2   4
      >>> print(t2)
       a   b
      --- ---
        5   7
        6   8
      >>> print(vstack([t1, t2]))
       a   b
      --- ---
        1   3
        2   4
        5   7
        6   8
    """
    _check_join_type(join_type, 'vstack')

    tables = _get_list_of_tables(tables)  # validates input
    if len(tables) == 1:
        return tables[0]  # no point in stacking a single table
    col_name_map = OrderedDict()

    out = _vstack(tables, join_type, col_name_map, metadata_conflicts)

    # Merge table metadata
    _merge_table_meta(out, tables, metadata_conflicts=metadata_conflicts)

    return out


def hstack(tables, join_type='outer',
           uniq_col_name='{col_name}_{table_name}', table_names=None,
           metadata_conflicts='warn'):
    """
    Stack tables along columns (horizontally)

    A ``join_type`` of 'exact' means that the tables must all
    have exactly the same number of rows.  If ``join_type`` is 'inner' then
    the intersection of rows will be the output.  A value of 'outer' (default)
    means the output will have the union of all rows, with table values being
    masked where no common values are available.

    Parameters
    ----------
    tables : `~astropy.table.Table` or `~astropy.table.Row` or list thereof
        Tables to stack along columns (horizontally) with the current table
    join_type : str
        Join type ('inner' | 'exact' | 'outer'), default is 'outer'
    uniq_col_name : str or None
        String generate a unique output column name in case of a conflict.
        The default is '{col_name}_{table_name}'.
    table_names : list of str or None
        Two-element list of table names used when generating unique output
        column names.  The default is ['1', '2', ..].
    metadata_conflicts : str
        How to proceed with metadata conflicts. This should be one of:
            * ``'silent'``: silently pick the last conflicting meta-data value
            * ``'warn'``: pick the last conflicting meta-data value,
              but emit a warning (default)
            * ``'error'``: raise an exception.

    Returns
    -------
    stacked_table : `~astropy.table.Table` object
        New table containing the stacked data from the input tables.

    See Also
    --------
    Table.add_columns, Table.replace_column, Table.update

    Examples
    --------
    To stack two tables horizontally (along columns) do::

      >>> from astropy.table import Table, hstack
      >>> t1 = Table({'a': [1, 2], 'b': [3, 4]}, names=('a', 'b'))
      >>> t2 = Table({'c': [5, 6], 'd': [7, 8]}, names=('c', 'd'))
      >>> print(t1)
       a   b
      --- ---
        1   3
        2   4
      >>> print(t2)
       c   d
      --- ---
        5   7
        6   8
      >>> print(hstack([t1, t2]))
       a   b   c   d
      --- --- --- ---
        1   3   5   7
        2   4   6   8
    """
    _check_join_type(join_type, 'hstack')

    tables = _get_list_of_tables(tables)  # validates input
    if len(tables) == 1:
        return tables[0]  # no point in stacking a single table
    col_name_map = OrderedDict()

    out = _hstack(tables, join_type, uniq_col_name, table_names,
                  col_name_map)

    _merge_table_meta(out, tables, metadata_conflicts=metadata_conflicts)

    return out


def unique(input_table, keys=None, silent=False, keep='first'):
    """
    Returns the unique rows of a table.

    Parameters
    ----------
    input_table : table-like
    keys : str or list of str
        Name(s) of column(s) used to create unique rows.
        Default is to use all columns.
    keep : {'first', 'last', 'none'}
        Whether to keep the first or last row for each set of
        duplicates. If 'none', all rows that are duplicate are
        removed, leaving only rows that are already unique in
        the input.
        Default is 'first'.
    silent : bool
        If `True`, masked value column(s) are silently removed from
        ``keys``. If `False`, an exception is raised when ``keys``
        contains masked value column(s).
        Default is `False`.

    Returns
    -------
    unique_table : `~astropy.table.Table` object
        New table containing only the unique rows of ``input_table``.

    Examples
    --------
    >>> from astropy.table import unique, Table
    >>> import numpy as np
    >>> table = Table(data=[[1,2,3,2,3,3],
    ... [2,3,4,5,4,6],
    ... [3,4,5,6,7,8]],
    ... names=['col1', 'col2', 'col3'],
    ... dtype=[np.int32, np.int32, np.int32])
    >>> table
    <Table length=6>
     col1  col2  col3
    int32 int32 int32
    ----- ----- -----
        1     2     3
        2     3     4
        3     4     5
        2     5     6
        3     4     7
        3     6     8
    >>> unique(table, keys='col1')
    <Table length=3>
     col1  col2  col3
    int32 int32 int32
    ----- ----- -----
        1     2     3
        2     3     4
        3     4     5
    >>> unique(table, keys=['col1'], keep='last')
    <Table length=3>
     col1  col2  col3
    int32 int32 int32
    ----- ----- -----
        1     2     3
        2     5     6
        3     6     8
    >>> unique(table, keys=['col1', 'col2'])
    <Table length=5>
     col1  col2  col3
    int32 int32 int32
    ----- ----- -----
        1     2     3
        2     3     4
        2     5     6
        3     4     5
        3     6     8
    >>> unique(table, keys=['col1', 'col2'], keep='none')
    <Table length=4>
     col1  col2  col3
    int32 int32 int32
    ----- ----- -----
        1     2     3
        2     3     4
        2     5     6
        3     6     8
    >>> unique(table, keys=['col1'], keep='none')
    <Table length=1>
     col1  col2  col3
    int32 int32 int32
    ----- ----- -----
        1     2     3

    """

    if keep not in ('first', 'last', 'none'):
        raise ValueError("'keep' should be one of 'first', 'last', 'none'")

    if isinstance(keys, str):
        keys = [keys]
    if keys is None:
        keys = input_table.colnames
    else:
        if len(set(keys)) != len(keys):
            raise ValueError("duplicate key names")

    # Check for columns with masked values
    for key in keys[:]:
        col = input_table[key]
        if hasattr(col, 'mask') and np.any(col.mask):
            if not silent:
                raise ValueError(
                    "cannot use columns with masked values as keys; "
                    "remove column '{}' from keys and rerun "
                    "unique()".format(key))
            del keys[keys.index(key)]
    if len(keys) == 0:
        raise ValueError("no column remained in ``keys``; "
                         "unique() cannot work with masked value "
                         "key columns")

    grouped_table = input_table.group_by(keys)
    indices = grouped_table.groups.indices
    if keep == 'first':
        indices = indices[:-1]
    elif keep == 'last':
        indices = indices[1:] - 1
    else:
        indices = indices[:-1][np.diff(indices) == 1]

    return grouped_table[indices]


def get_col_name_map(arrays, common_names, uniq_col_name='{col_name}_{table_name}',
                     table_names=None):
    """
    Find the column names mapping when merging the list of tables
    ``arrays``.  It is assumed that col names in ``common_names`` are to be
    merged into a single column while the rest will be uniquely represented
    in the output.  The args ``uniq_col_name`` and ``table_names`` specify
    how to rename columns in case of conflicts.

    Returns a dict mapping each output column name to the input(s).  This takes the form
    {outname : (col_name_0, col_name_1, ...), ... }.  For key columns all of input names
    will be present, while for the other non-key columns the value will be (col_name_0,
    None, ..) or (None, col_name_1, ..) etc.
    """

    col_name_map = collections.defaultdict(lambda: [None] * len(arrays))
    col_name_list = []

    if table_names is None:
        table_names = [str(ii + 1) for ii in range(len(arrays))]

    for idx, array in enumerate(arrays):
        table_name = table_names[idx]
        for name in array.colnames:
            out_name = name

            if name in common_names:
                # If name is in the list of common_names then insert into
                # the column name list, but just once.
                if name not in col_name_list:
                    col_name_list.append(name)
            else:
                # If name is not one of the common column outputs, and it collides
                # with the names in one of the other arrays, then rename
                others = list(arrays)
                others.pop(idx)
                if any(name in other.colnames for other in others):
                    out_name = uniq_col_name.format(table_name=table_name, col_name=name)
                col_name_list.append(out_name)

            col_name_map[out_name][idx] = name

    # Check for duplicate output column names
    col_name_count = Counter(col_name_list)
    repeated_names = [name for name, count in col_name_count.items() if count > 1]
    if repeated_names:
        raise TableMergeError('Merging column names resulted in duplicates: {}.  '
                              'Change uniq_col_name or table_names args to fix this.'
                              .format(repeated_names))

    # Convert col_name_map to a regular dict with tuple (immutable) values
    col_name_map = OrderedDict((name, col_name_map[name]) for name in col_name_list)

    return col_name_map


def get_descrs(arrays, col_name_map):
    """
    Find the dtypes descrs resulting from merging the list of arrays' dtypes,
    using the column name mapping ``col_name_map``.

    Return a list of descrs for the output.
    """

    out_descrs = []

    for out_name, in_names in col_name_map.items():
        # List of input arrays that contribute to this output column
        in_cols = [arr[name] for arr, name in zip(arrays, in_names) if name is not None]

        # List of names of the columns that contribute to this output column.
        names = [name for name in in_names if name is not None]

        # Output dtype is the superset of all dtypes in in_arrays
        try:
            dtype = common_dtype(in_cols)
        except TableMergeError as tme:
            # Beautify the error message when we are trying to merge columns with incompatible
            # types by including the name of the columns that originated the error.
            raise TableMergeError("The '{}' columns have incompatible types: {}"
                                  .format(names[0], tme._incompat_types)) from tme

        # Make sure all input shapes are the same
        uniq_shapes = set(col.shape[1:] for col in in_cols)
        if len(uniq_shapes) != 1:
            raise TableMergeError(f'Key columns {names!r} have different shape')
        shape = uniq_shapes.pop()

        if out_name is not None:
            out_name = str(out_name)
        out_descrs.append((out_name, dtype, shape))

    return out_descrs


def common_dtype(cols):
    """
    Use numpy to find the common dtype for a list of columns.

    Only allow columns within the following fundamental numpy data types:
    np.bool_, np.object_, np.number, np.character, np.void
    """
    try:
        return metadata.common_dtype(cols)
    except metadata.MergeConflictError as err:
        tme = TableMergeError(f'Columns have incompatible types {err._incompat_types}')
        tme._incompat_types = err._incompat_types
        raise tme from err


def _get_join_sort_idxs(keys, left, right):
    # Go through each of the key columns in order and make columns for
    # a new structured array that represents the lexical ordering of those
    # key columns. This structured array is then argsort'ed. The trick here
    # is that some columns (e.g. Time) may need to be expanded into multiple
    # columns for ordering here.

    ii = 0  # Index for uniquely naming the sort columns
    sort_keys_dtypes = []  # sortable_table dtypes as list of (name, dtype_str, shape) tuples
    sort_keys = []  # sortable_table (structured ndarray) column names
    sort_left = {}  # sortable ndarrays from left table
    sort_right = {}  # sortable ndarray from right table

    for key in keys:
        # get_sortable_arrays() returns a list of ndarrays that can be lexically
        # sorted to represent the order of the column. In most cases this is just
        # a single element of the column itself.
        left_sort_cols = left[key].info.get_sortable_arrays()
        right_sort_cols = right[key].info.get_sortable_arrays()

        if len(left_sort_cols) != len(right_sort_cols):
            # Should never happen because cols are screened beforehand for compatibility
            raise RuntimeError('mismatch in sort cols lengths')

        for left_sort_col, right_sort_col in zip(left_sort_cols, right_sort_cols):
            # Check for consistency of shapes. Mismatch should never happen.
            shape = left_sort_col.shape[1:]
            if shape != right_sort_col.shape[1:]:
                raise RuntimeError('mismatch in shape of left vs. right sort array')

            if shape != ():
                raise ValueError(f'sort key column {key!r} must be 1-d')

            sort_key = str(ii)
            sort_keys.append(sort_key)
            sort_left[sort_key] = left_sort_col
            sort_right[sort_key] = right_sort_col

            # Build up dtypes for the structured array that gets sorted.
            dtype_str = common_dtype([left_sort_col, right_sort_col])
            sort_keys_dtypes.append((sort_key, dtype_str))
            ii += 1

    # Make the empty sortable table and fill it
    len_left = len(left)
    sortable_table = np.empty(len_left + len(right), dtype=sort_keys_dtypes)
    for key in sort_keys:
        sortable_table[key][:len_left] = sort_left[key]
        sortable_table[key][len_left:] = sort_right[key]

    # Finally do the (lexical) argsort and make a new sorted version
    idx_sort = sortable_table.argsort(order=sort_keys)
    sorted_table = sortable_table[idx_sort]

    # Get indexes of unique elements (i.e. the group boundaries)
    diffs = np.concatenate(([True], sorted_table[1:] != sorted_table[:-1], [True]))
    idxs = np.flatnonzero(diffs)

    return idxs, idx_sort


def _apply_join_funcs(left, right, keys, join_funcs):
    """Apply join_funcs
    """
    # Make light copies of left and right, then add new index columns.
    left = left.copy(copy_data=False)
    right = right.copy(copy_data=False)
    for key, join_func in join_funcs.items():
        ids1, ids2 = join_func(left[key], right[key])
        # Define a unique id_key name, and keep adding underscores until we have
        # a name not yet present.
        id_key = key + '_id'
        while id_key in left.columns or id_key in right.columns:
            id_key = id_key[:-2] + '_id'

        keys = tuple(id_key if orig_key == key else orig_key for orig_key in keys)
        left.add_column(ids1, index=0, name=id_key)  # [id_key] = ids1
        right.add_column(ids2, index=0, name=id_key)  # [id_key] = ids2

    return left, right, keys


def _join(left, right, keys=None, join_type='inner',
          uniq_col_name='{col_name}_{table_name}',
          table_names=['1', '2'],
          col_name_map=None, metadata_conflicts='warn',
          join_funcs=None,
          keys_left=None, keys_right=None):
    """
    Perform a join of the left and right Tables on specified keys.

    Parameters
    ----------
    left : Table
        Left side table in the join
    right : Table
        Right side table in the join
    keys : str or list of str
        Name(s) of column(s) used to match rows of left and right tables.
        Default is to use all columns which are common to both tables.
    join_type : str
        Join type ('inner' | 'outer' | 'left' | 'right' | 'cartesian'), default is 'inner'
    uniq_col_name : str or None
        String generate a unique output column name in case of a conflict.
        The default is '{col_name}_{table_name}'.
    table_names : list of str or None
        Two-element list of table names used when generating unique output
        column names.  The default is ['1', '2'].
    col_name_map : empty dict or None
        If passed as a dict then it will be updated in-place with the
        mapping of output to input column names.
    metadata_conflicts : str
        How to proceed with metadata conflicts. This should be one of:
            * ``'silent'``: silently pick the last conflicting meta-data value
            * ``'warn'``: pick the last conflicting meta-data value, but emit a warning (default)
            * ``'error'``: raise an exception.
    join_funcs : dict, None
        Dict of functions to use for matching the corresponding key column(s).
        See `~astropy.table.join_skycoord` for an example and details.

    Returns
    -------
    joined_table : `~astropy.table.Table` object
        New table containing the result of the join operation.
    """
    # Store user-provided col_name_map until the end
    _col_name_map = col_name_map

    # Special column name for cartesian join, should never collide with real column
    cartesian_index_name = '__table_cartesian_join_temp_index__'

    if join_type not in ('inner', 'outer', 'left', 'right', 'cartesian'):
        raise ValueError("The 'join_type' argument should be in 'inner', "
                         "'outer', 'left', 'right', or 'cartesian' "
                         "(got '{}' instead)".
                         format(join_type))

    if join_type == 'cartesian':
        if keys:
            raise ValueError('cannot supply keys for a cartesian join')

        if join_funcs:
            raise ValueError('cannot supply join_funcs for a cartesian join')

        # Make light copies of left and right, then add temporary index columns
        # with all the same value so later an outer join turns into a cartesian join.
        left = left.copy(copy_data=False)
        right = right.copy(copy_data=False)
        left[cartesian_index_name] = np.uint8(0)
        right[cartesian_index_name] = np.uint8(0)
        keys = (cartesian_index_name, )

    # Handle the case of join key columns that are different between left and
    # right via keys_left/keys_right args. This is done by saving the original
    # input tables and making new left and right tables that contain only the
    # key cols but with common column names ['0', '1', etc]. This sets `keys` to
    # those fake key names in the left and right tables
    if keys_left is not None or keys_right is not None:
        left_orig = left
        right_orig = right
        left, right, keys = _join_keys_left_right(
            left, right, keys, keys_left, keys_right, join_funcs)

    if keys is None:
        keys = tuple(name for name in left.colnames if name in right.colnames)
        if len(keys) == 0:
            raise TableMergeError('No keys in common between left and right tables')
    elif isinstance(keys, str):
        # If we have a single key, put it in a tuple
        keys = (keys,)

    # Check the key columns
    for arr, arr_label in ((left, 'Left'), (right, 'Right')):
        for name in keys:
            if name not in arr.colnames:
                raise TableMergeError('{} table does not have key column {!r}'
                                      .format(arr_label, name))
            if hasattr(arr[name], 'mask') and np.any(arr[name].mask):
                raise TableMergeError('{} key column {!r} has missing values'
                                      .format(arr_label, name))

    if join_funcs is not None:
        if not all(key in keys for key in join_funcs):
            raise ValueError(f'join_funcs keys {join_funcs.keys()} must be a '
                             f'subset of join keys {keys}')
        left, right, keys = _apply_join_funcs(left, right, keys, join_funcs)

    len_left, len_right = len(left), len(right)

    if len_left == 0 or len_right == 0:
        raise ValueError('input tables for join must both have at least one row')

    try:
        idxs, idx_sort = _get_join_sort_idxs(keys, left, right)
    except NotImplementedError:
        raise TypeError('one or more key columns are not sortable')

    # Now that we have idxs and idx_sort, revert to the original table args to
    # carry on with making the output joined table. `keys` is set to to an empty
    # list so that all original left and right columns are included in the
    # output table.
    if keys_left is not None or keys_right is not None:
        keys = []
        left = left_orig
        right = right_orig

    # Joined array dtype as a list of descr (name, type_str, shape) tuples
    col_name_map = get_col_name_map([left, right], keys, uniq_col_name, table_names)
    out_descrs = get_descrs([left, right], col_name_map)

    # Main inner loop in Cython to compute the cartesian product
    # indices for the given join type
    int_join_type = {'inner': 0, 'outer': 1, 'left': 2, 'right': 3,
                     'cartesian': 1}[join_type]
    masked, n_out, left_out, left_mask, right_out, right_mask = \
        _np_utils.join_inner(idxs, idx_sort, len_left, int_join_type)

    out = _get_out_class([left, right])()

    for out_name, dtype, shape in out_descrs:
        if out_name == cartesian_index_name:
            continue

        left_name, right_name = col_name_map[out_name]
        if left_name and right_name:  # this is a key which comes from left and right
            cols = [left[left_name], right[right_name]]

            col_cls = _get_out_class(cols)
            if not hasattr(col_cls.info, 'new_like'):
                raise NotImplementedError('join unavailable for mixin column type(s): {}'
                                          .format(col_cls.__name__))

            out[out_name] = col_cls.info.new_like(cols, n_out, metadata_conflicts, out_name)
            out[out_name][:] = np.where(right_mask,
                                        left[left_name].take(left_out),
                                        right[right_name].take(right_out))
            continue
        elif left_name:  # out_name came from the left table
            name, array, array_out, array_mask = left_name, left, left_out, left_mask
        elif right_name:
            name, array, array_out, array_mask = right_name, right, right_out, right_mask
        else:
            raise TableMergeError('Unexpected column names (maybe one is ""?)')

        # Select the correct elements from the original table
        col = array[name][array_out]

        # If the output column is masked then set the output column masking
        # accordingly.  Check for columns that don't support a mask attribute.
        if masked and np.any(array_mask):
            # If col is a Column but not MaskedColumn then upgrade at this point
            # because masking is required.
            if isinstance(col, Column) and not isinstance(col, MaskedColumn):
                col = out.MaskedColumn(col, copy=False)

            if isinstance(col, Quantity) and not isinstance(col, Masked):
                col = Masked(col, copy=False)

            # array_mask is 1-d corresponding to length of output column.  We need
            # make it have the correct shape for broadcasting, i.e. (length, 1, 1, ..).
            # Mixin columns might not have ndim attribute so use len(col.shape).
            array_mask.shape = (col.shape[0],) + (1,) * (len(col.shape) - 1)

            # Now broadcast to the correct final shape
            array_mask = np.broadcast_to(array_mask, col.shape)

            try:
                col[array_mask] = col.info.mask_val
            except Exception as err:  # Not clear how different classes will fail here
                raise NotImplementedError(
                    "join requires masking column '{}' but column"
                    " type {} does not support masking"
                    .format(out_name, col.__class__.__name__)) from err

        # Set the output table column to the new joined column
        out[out_name] = col

    # If col_name_map supplied as a dict input, then update.
    if isinstance(_col_name_map, Mapping):
        _col_name_map.update(col_name_map)

    return out


def _join_keys_left_right(left, right, keys, keys_left, keys_right, join_funcs):
    """Do processing to handle keys_left / keys_right args for join.

    This takes the keys_left/right inputs and turns them into a list of left/right
    columns corresponding to those inputs (which can be column names or column
    data values). It also generates the list of fake key column names (strings
    of "1", "2", etc.) that correspond to the input keys.
    """
    def _keys_to_cols(keys, table, label):
        # Process input `keys`, which is a str or list of str column names in
        # `table` or a list of column-like objects. The `label` is just for
        # error reporting.
        if isinstance(keys, str):
            keys = [keys]
        cols = []
        for key in keys:
            if isinstance(key, str):
                try:
                    cols.append(table[key])
                except KeyError:
                    raise ValueError(f'{label} table does not have key column {key!r}')
            else:
                if len(key) != len(table):
                    raise ValueError(f'{label} table has different length from key {key}')
                cols.append(key)
        return cols

    if join_funcs is not None:
        raise ValueError('cannot supply join_funcs arg and keys_left / keys_right')

    if keys_left is None or keys_right is None:
        raise ValueError('keys_left and keys_right must both be provided')

    if keys is not None:
        raise ValueError('keys arg must be None if keys_left and keys_right are supplied')

    cols_left = _keys_to_cols(keys_left, left, 'left')
    cols_right = _keys_to_cols(keys_right, right, 'right')

    if len(cols_left) != len(cols_right):
        raise ValueError('keys_left and keys_right args must have same length')

    # Make two new temp tables for the join with only the join columns and
    # key columns in common.
    keys = [f'{ii}' for ii in range(len(cols_left))]

    left = left.__class__(cols_left, names=keys, copy=False)
    right = right.__class__(cols_right, names=keys, copy=False)

    return left, right, keys


def _check_join_type(join_type, func_name):
    """Check join_type arg in hstack and vstack.

    This specifically checks for the common mistake of call vstack(t1, t2)
    instead of vstack([t1, t2]). The subsequent check of
    ``join_type in ('inner', ..)`` does not raise in this case.
    """
    if not isinstance(join_type, str):
        msg = '`join_type` arg must be a string'
        if isinstance(join_type, Table):
            msg += ('. Did you accidentally '
                    f'call {func_name}(t1, t2, ..) instead of '
                    f'{func_name}([t1, t2], ..)?')
        raise TypeError(msg)

    if join_type not in ('inner', 'exact', 'outer'):
        raise ValueError("`join_type` arg must be one of 'inner', 'exact' or 'outer'")


def _vstack(arrays, join_type='outer', col_name_map=None, metadata_conflicts='warn'):
    """
    Stack Tables vertically (by rows)

    A ``join_type`` of 'exact' (default) means that the arrays must all
    have exactly the same column names (though the order can vary).  If
    ``join_type`` is 'inner' then the intersection of common columns will
    be the output.  A value of 'outer' means the output will have the union of
    all columns, with array values being masked where no common values are
    available.

    Parameters
    ----------
    arrays : list of Tables
        Tables to stack by rows (vertically)
    join_type : str
        Join type ('inner' | 'exact' | 'outer'), default is 'outer'
    col_name_map : empty dict or None
        If passed as a dict then it will be updated in-place with the
        mapping of output to input column names.

    Returns
    -------
    stacked_table : `~astropy.table.Table` object
        New table containing the stacked data from the input tables.
    """
    # Store user-provided col_name_map until the end
    _col_name_map = col_name_map

    # Trivial case of one input array
    if len(arrays) == 1:
        return arrays[0]

    # Start by assuming an outer match where all names go to output
    names = set(itertools.chain(*[arr.colnames for arr in arrays]))
    col_name_map = get_col_name_map(arrays, names)

    # If require_match is True then the output must have exactly the same
    # number of columns as each input array
    if join_type == 'exact':
        for names in col_name_map.values():
            if any(x is None for x in names):
                raise TableMergeError('Inconsistent columns in input arrays '
                                      "(use 'inner' or 'outer' join_type to "
                                      "allow non-matching columns)")
        join_type = 'outer'

    # For an inner join, keep only columns where all input arrays have that column
    if join_type == 'inner':
        col_name_map = OrderedDict((name, in_names) for name, in_names in col_name_map.items()
                                   if all(x is not None for x in in_names))
        if len(col_name_map) == 0:
            raise TableMergeError('Input arrays have no columns in common')

    lens = [len(arr) for arr in arrays]
    n_rows = sum(lens)
    out = _get_out_class(arrays)()

    for out_name, in_names in col_name_map.items():
        # List of input arrays that contribute to this output column
        cols = [arr[name] for arr, name in zip(arrays, in_names) if name is not None]

        col_cls = _get_out_class(cols)
        if not hasattr(col_cls.info, 'new_like'):
            raise NotImplementedError('vstack unavailable for mixin column type(s): {}'
                                      .format(col_cls.__name__))
        try:
            col = col_cls.info.new_like(cols, n_rows, metadata_conflicts, out_name)
        except metadata.MergeConflictError as err:
            # Beautify the error message when we are trying to merge columns with incompatible
            # types by including the name of the columns that originated the error.
            raise TableMergeError("The '{}' columns have incompatible types: {}"
                                  .format(out_name, err._incompat_types)) from err

        idx0 = 0
        for name, array in zip(in_names, arrays):
            idx1 = idx0 + len(array)
            if name in array.colnames:
                col[idx0:idx1] = array[name]
            else:
                # If col is a Column but not MaskedColumn then upgrade at this point
                # because masking is required.
                if isinstance(col, Column) and not isinstance(col, MaskedColumn):
                    col = out.MaskedColumn(col, copy=False)

                if isinstance(col, Quantity) and not isinstance(col, Masked):
                    col = Masked(col, copy=False)

                try:
                    col[idx0:idx1] = col.info.mask_val
                except Exception as err:
                    raise NotImplementedError(
                        "vstack requires masking column '{}' but column"
                        " type {} does not support masking"
                        .format(out_name, col.__class__.__name__)) from err
            idx0 = idx1

        out[out_name] = col

    # If col_name_map supplied as a dict input, then update.
    if isinstance(_col_name_map, Mapping):
        _col_name_map.update(col_name_map)

    return out


def _hstack(arrays, join_type='outer', uniq_col_name='{col_name}_{table_name}',
            table_names=None, col_name_map=None):
    """
    Stack tables horizontally (by columns)

    A ``join_type`` of 'exact' (default) means that the arrays must all
    have exactly the same number of rows.  If ``join_type`` is 'inner' then
    the intersection of rows will be the output.  A value of 'outer' means
    the output will have the union of all rows, with array values being
    masked where no common values are available.

    Parameters
    ----------
    arrays : List of tables
        Tables to stack by columns (horizontally)
    join_type : str
        Join type ('inner' | 'exact' | 'outer'), default is 'outer'
    uniq_col_name : str or None
        String generate a unique output column name in case of a conflict.
        The default is '{col_name}_{table_name}'.
    table_names : list of str or None
        Two-element list of table names used when generating unique output
        column names.  The default is ['1', '2', ..].

    Returns
    -------
    stacked_table : `~astropy.table.Table` object
        New table containing the stacked data from the input tables.
    """

    # Store user-provided col_name_map until the end
    _col_name_map = col_name_map

    if table_names is None:
        table_names = [f'{ii + 1}' for ii in range(len(arrays))]
    if len(arrays) != len(table_names):
        raise ValueError('Number of arrays must match number of table_names')

    # Trivial case of one input arrays
    if len(arrays) == 1:
        return arrays[0]

    col_name_map = get_col_name_map(arrays, [], uniq_col_name, table_names)

    # If require_match is True then all input arrays must have the same length
    arr_lens = [len(arr) for arr in arrays]
    if join_type == 'exact':
        if len(set(arr_lens)) > 1:
            raise TableMergeError("Inconsistent number of rows in input arrays "
                                  "(use 'inner' or 'outer' join_type to allow "
                                  "non-matching rows)")
        join_type = 'outer'

    # For an inner join, keep only the common rows
    if join_type == 'inner':
        min_arr_len = min(arr_lens)
        if len(set(arr_lens)) > 1:
            arrays = [arr[:min_arr_len] for arr in arrays]
        arr_lens = [min_arr_len for arr in arrays]

    # If there are any output rows where one or more input arrays are missing
    # then the output must be masked.  If any input arrays are masked then
    # output is masked.

    n_rows = max(arr_lens)
    out = _get_out_class(arrays)()

    for out_name, in_names in col_name_map.items():
        for name, array, arr_len in zip(in_names, arrays, arr_lens):
            if name is None:
                continue

            if n_rows > arr_len:
                indices = np.arange(n_rows)
                indices[arr_len:] = 0
                col = array[name][indices]

                # If col is a Column but not MaskedColumn then upgrade at this point
                # because masking is required.
                if isinstance(col, Column) and not isinstance(col, MaskedColumn):
                    col = out.MaskedColumn(col, copy=False)

                if isinstance(col, Quantity) and not isinstance(col, Masked):
                    col = Masked(col, copy=False)

                try:
                    col[arr_len:] = col.info.mask_val
                except Exception as err:
                    raise NotImplementedError(
                        "hstack requires masking column '{}' but column"
                        " type {} does not support masking"
                        .format(out_name, col.__class__.__name__)) from err
            else:
                col = array[name][:n_rows]

            out[out_name] = col

    # If col_name_map supplied as a dict input, then update.
    if isinstance(_col_name_map, Mapping):
        _col_name_map.update(col_name_map)

    return out
