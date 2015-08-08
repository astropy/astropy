# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import numpy as np

class SortedArray(object):
    '''
    Implements a sorted array container using
    a list of numpy arrays.

    Parameters
    ----------
    data : `Table`
        Sorted columns of the original table
    row_index : Column object
        Row numbers corresponding to data columns
    '''
    def __init__(self, data, row_index):
        self._data = data
        self.length = len(data)
        self.num_cols = len(data.colnames)
        self._row_index = row_index

    def set_col(i, val):
        self._data[self._data.colnames[i]] = val

    def resize(self):
        '''
        Regrow the internal numpy column buffers if necessary.
        '''
        buffer_size = len(self._data[0])
        if self.length == buffer_size:
            # resize to 50% larger capacity
            increase = int(0.5 * buffer_size)
            if increase == 0: # started out with buffer_size == 1
                increase = 1
            buffer_size += increase
            for i, col in enumerate(self._data.columns.values()):
                self.set_col(i, np.resize(col, buffer_size))
            for i in range(
            self._row_index[

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
        self.resize()
        pos = self.find_pos(key, row) # first >= key
        key = key + (row,)

        # shift larger keys rightward
        for i, col in enumerate(self._data):
            col[pos:] = np.roll(col[pos:], 1)
            col[pos] = key[i]

        self.length += 1

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
        end = self.length
        key = key + (data,)

        # search through keys in lexicographic order
        for i in range(self.num_cols + 1):
            key_slice = self._data[i][begin:end]
            t = np.searchsorted(key_slice, key[i])
            # t is the smallest index >= key[i]
            if exact and (t == len(key_slice) or key_slice[t] != key[i]):
                # no match
                return -1
            elif t == len(key_slice) or (t == 0 and len(key_slice > 0) and
                                         key_slice[0] > key[i]):
                # too small or too large
                return begin + t
            end = begin + np.searchsorted(key_slice, key[i], side='right')
            begin += t
            if begin >= self.length: # greater than all keys
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
        end = self.length

        # search through keys in lexicographic order
        for i in range(self.num_cols):
            key_slice = self._data[i][begin:end]
            t = np.searchsorted(key_slice, key[i])
            # t is the smallest index >= key[i]
            if t == len(key_slice) or key_slice[t] != key[i]:
                # no match
                return []
            elif t == 0 and len(key_slice > 0) and key_slice[0] > key[i]:
                # too small or too large
                return []
            end = begin + np.searchsorted(key_slice, key[i], side='right')
            begin += t
            if begin >= self.length: # greater than all keys
                return []

        return self._data[-1][begin:end]

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
        if lower_pos == self.length:
            return []

        lower_bound = tuple([col[lower_pos] for col in self._data[:-1]])
        if not bounds[0] and lower_bound == lower:
            lower_pos += 1 # data[lower_pos] > lower

        # data[lower_pos] >= lower
        # data[upper_pos] >= upper
        if upper_pos < self.length:
            upper_bound = tuple([col[upper_pos] for col in self._data[:-1]])
            if not bounds[1] and upper_bound == upper:
                upper_pos -= 1 # data[upper_pos] < upper
            elif upper_bound > upper:
                upper_pos -= 1 # data[upper_pos] <= upper
        return self._data[-1][lower_pos:upper_pos + 1]

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

        # shift larger keys leftward
        for col in self._data:
            col[pos:] = np.roll(col[pos:], -1)
        self.length -= 1
        return True

    def shift_left(self, row):
        '''
        Decrement all row numbers greater than the input row.

        Parameters
        ----------
        row : int
            Input row number
        '''
        self._data[-1][self._data[-1] > row] -= 1

    def shift_right(self, row):
        '''
        Increment all row numbers greater than or equal to the input row.

        Parameters
        ----------
        row : int
            Input row number
        '''
        self._data[-1][self.data[-1] >= row] += 1

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
        for i, col in enumerate(self._data):
            self._data[i] = col[np.array([x in row_map for x in self._data[-1]])]
        self._data[-1] = np.array([row_map[x] for x in self._data[-1]])

    def items(self):
        '''
        Retrieve all array items as a list of pairs of the form
        [(key, [row 1, row 2, ...]), ...]
        '''
        array = []
        last_key = None
        for line in zip(*self.data):
            key, row = line[:-1], line[-1]
            if key == last_key:
                array[-1][1].append(row)
            else:
                last_key = key
                array.append((key, [row]))
        return array

    def sort(self):
        '''
        Return rows in sorted order.
        '''
        return self.data[-1]

    @property
    def data(self):
        return [col[:self.length] for col in self._data]

    def __getitem__(self, item):
        '''
        Return a sliced reference to this sorted array.

        Parameters
        ----------
        item : slice
            Slice to use for referencing
        '''
        ######TODO: test for correct out-of-bounds (data vs _data)
        return SortedArray([col[item] for col in self._data])

    def __repr__(self):
        return '[' + ', '.join([str(x) for x in self.data]) + ']'

    def __str__(self):
        return repr(self)
