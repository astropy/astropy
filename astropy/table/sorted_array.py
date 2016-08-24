# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import numpy as np
from ..extern.six.moves import range, zip

def _searchsorted(array, val, side='left'):
    '''
    Call np.searchsorted or use a custom binary
    search if necessary.
    '''
    if hasattr(array, 'searchsorted'):
        return array.searchsorted(val, side=side)
    # Python binary search
    begin = 0
    end = len(array)
    while begin < end:
        mid = (begin + end) // 2
        if val > array[mid]:
            begin = mid + 1
        elif val < array[mid]:
            end = mid
        elif side == 'right':
            begin = mid + 1
        else:
            end = mid
    return begin

class SortedArray(object):
    '''
    Implements a sorted array container using
    a list of numpy arrays.

    Parameters
    ----------
    data : Table
        Sorted columns of the original table
    row_index : Column object
        Row numbers corresponding to data columns
    unique : bool (defaults to False)
        Whether the values of the index must be unique
    '''
    def __init__(self, data, row_index, unique=False):
        self.data = data
        self.row_index = row_index
        self.num_cols = len(getattr(data, 'colnames', []))
        self.unique = unique

    @property
    def cols(self):
        return self.data.columns.values()

    def add(self, key, row):
        '''
        Add a new entry to the sorted array.

        Parameters
        ----------
        key : tuple
            Column values at the given row
        row : int
            Row number
        '''
        pos = self.find_pos(key, row) # first >= key

        if self.unique and 0 <= pos < len(self.row_index) and \
           all(self.data[pos][i] == key[i] for i in range(len(key))):
            # already exists
            raise ValueError('Cannot add duplicate value "{0}" in a '
                             'unique index'.format(key))
        self.data.insert_row(pos, key)
        self.row_index = self.row_index.insert(pos, row)

    def _get_key_slice(self, i, begin, end):
        '''
        Retrieve the ith slice of the sorted array
        from begin to end.
        '''
        if i < self.num_cols:
            return self.cols[i][begin:end]
        else:
            return self.row_index[begin:end]

    def find_pos(self, key, data, exact=False):
        '''
        Return the index of the largest key in data greater than or
        equal to the given key, data pair.

        Parameters
        ----------
        key : tuple
            Column key
        data : int
            Row number
        exact : bool
            If True, return the index of the given key in data
            or -1 if the key is not present.
        '''
        begin = 0
        end = len(self.row_index)
        num_cols = self.num_cols
        if not self.unique:
            # consider the row value as well
            key = key + (data,)
            num_cols += 1

        # search through keys in lexicographic order
        for i in range(num_cols):
            key_slice = self._get_key_slice(i, begin, end)
            t = _searchsorted(key_slice, key[i])
            # t is the smallest index >= key[i]
            if exact and (t == len(key_slice) or key_slice[t] != key[i]):
                # no match
                return -1
            elif t == len(key_slice) or (t == 0 and len(key_slice) > 0 and
                                         key[i] < key_slice[0]):
                # too small or too large
                return begin + t
            end = begin + _searchsorted(key_slice, key[i], side='right')
            begin += t
            if begin >= len(self.row_index): # greater than all keys
                return begin

        return begin


    def find(self, key):
        '''
        Find all rows matching the given key.

        Parameters
        ----------
        key : tuple
            Column values

        Returns
        -------
        matching_rows : list
            List of rows matching the input key
        '''
        begin = 0
        end = len(self.row_index)

        # search through keys in lexicographic order
        for i in range(self.num_cols):
            key_slice = self._get_key_slice(i, begin, end)
            t = _searchsorted(key_slice, key[i])
            # t is the smallest index >= key[i]
            if t == len(key_slice) or key_slice[t] != key[i]:
                # no match
                return []
            elif t == 0 and len(key_slice) > 0 and key[i] < key_slice[0]:
                # too small or too large
                return []
            end = begin + _searchsorted(key_slice, key[i], side='right')
            begin += t
            if begin >= len(self.row_index): # greater than all keys
                return []

        return self.row_index[begin:end]

    def range(self, lower, upper, bounds):
        '''
        Find values in the given range.

        Parameters
        ----------
        lower : tuple
            Lower search bound
        upper : tuple
            Upper search bound
        bounds : tuple (x, y) of bools
            Indicates whether the search should be inclusive or
            exclusive with respect to the endpoints. The first
            argument x corresponds to an inclusive lower bound,
            and the second argument y to an inclusive upper bound.
        '''
        lower_pos = self.find_pos(lower, 0)
        upper_pos = self.find_pos(upper, 0)
        if lower_pos == len(self.row_index):
            return []

        lower_bound = tuple([col[lower_pos] for col in self.cols])
        if not bounds[0] and lower_bound == lower:
            lower_pos += 1 # data[lower_pos] > lower

        # data[lower_pos] >= lower
        # data[upper_pos] >= upper
        if upper_pos < len(self.row_index):
            upper_bound = tuple([col[upper_pos] for col in self.cols])
            if not bounds[1] and upper_bound == upper:
                upper_pos -= 1 # data[upper_pos] < upper
            elif upper_bound > upper:
                upper_pos -= 1 # data[upper_pos] <= upper
        return self.row_index[lower_pos:upper_pos + 1]

    def remove(self, key, data):
        '''
        Remove the given entry from the sorted array.

        Parameters
        ----------
        key : tuple
            Column values
        data : int
            Row number

        Returns
        -------
        successful : bool
            Whether the entry was successfully removed
        '''
        pos = self.find_pos(key, data, exact=True)
        if pos == -1: # key not found
            return False

        self.data.remove_row(pos)
        keep_mask = np.ones(len(self.row_index), dtype=np.bool)
        keep_mask[pos] = False
        self.row_index = self.row_index[keep_mask]
        return True

    def shift_left(self, row):
        '''
        Decrement all row numbers greater than the input row.

        Parameters
        ----------
        row : int
            Input row number
        '''
        self.row_index[self.row_index > row] -= 1

    def shift_right(self, row):
        '''
        Increment all row numbers greater than or equal to the input row.

        Parameters
        ----------
        row : int
            Input row number
        '''
        self.row_index[self.row_index >= row] += 1

    def replace_rows(self, row_map):
        '''
        Replace all rows with the values they map to in the
        given dictionary. Any rows not present as keys in
        the dictionary will have their entries deleted.

        Parameters
        ----------
        row_map : dict
            Mapping of row numbers to new row numbers
        '''
        num_rows = len(row_map)
        keep_rows = np.zeros(len(self.row_index), dtype=np.bool)
        tagged = 0
        for i, row in enumerate(self.row_index):
            if row in row_map:
                keep_rows[i] = True
                tagged += 1
                if tagged == num_rows:
                    break

        self.data = self.data[keep_rows]
        self.row_index = np.array(
            [row_map[x] for x in self.row_index[keep_rows]])

    def items(self):
        '''
        Retrieve all array items as a list of pairs of the form
        [(key, [row 1, row 2, ...]), ...]
        '''
        array = []
        last_key = None
        for i, key in enumerate(zip(*self.data.columns.values())):
            row = self.row_index[i]
            if key == last_key:
                array[-1][1].append(row)
            else:
                last_key = key
                array.append((key, [row]))
        return array

    def sort(self):
        '''
        Make row order align with key order.
        '''
        self.row_index = np.arange(len(self.row_index))

    def sorted_data(self):
        '''
        Return rows in sorted order.
        '''
        return self.row_index

    def __getitem__(self, item):
        '''
        Return a sliced reference to this sorted array.

        Parameters
        ----------
        item : slice
            Slice to use for referencing
        '''
        return SortedArray(self.data[item], self.row_index[item])

    def __repr__(self):
        t = self.data.copy()
        t['rows'] = self.row_index
        return str(t)

    def __str__(self):
        return repr(self)
