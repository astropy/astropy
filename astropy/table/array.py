import numpy as np

class ArraySlot(object):
    __lt__ = lambda x, y: x.key < y.key
    __le__ = lambda x, y: x.key <= y.key
    __eq__ = lambda x, y: x.key == y.key
    __ge__ = lambda x, y: x.key >= y.key
    __gt__ = lambda x, y: x.key > y.key
    __ne__ = lambda x, y: x.key != y.key

    def __init__(self, key, val):
        self.key = key
        self.data = val
    def __repr__(self):
        return str((self.key, self.data))
    def __str__(self):
        return repr(self)

class SortedArray(object):
    '''
    Implements a sorted array container using
    a numpy structured array.
    '''

    def __init__(self, lines, num_cols=None, data=None):
        if data is not None: # sliced reference to data
            self.length = len(data)
            self._data = data
            return
        elif num_cols is not None:
            self.num_cols = num_cols
        else:
            self.num_cols = len(lines.colnames) - 1

        self.length = 0 if lines is None else len(lines)
        key_cols = [('k{0}'.format(i), 'O') for i in range(self.num_cols)]
        self._data = np.zeros(self.length, dtype=key_cols + [('row', 'i')])
        self.cols = ['k{0}'.format(i) for i in range(self.num_cols)] + ['row']

        if len(lines) > 0:
            for i in range(self.num_cols):
                self._data['k{0}'.format(i)] = lines[lines.colnames[i]]
            self._data['row'] = lines[lines.colnames[self.num_cols]]

    def resize(self):
        buffer_size = len(self._data)
        if self.length == buffer_size:
            # resize to 50% larger capacity
            increase = int(0.5 * buffer_size)
            if increase == 0: # started out with buffer_size == 1
                increase = 1
            buffer_size += increase
            self._data = np.resize(self._data, buffer_size)

    def add(self, key, row):
        self.resize()
        pos = self.find_pos(key, row) # first >= key

        # shift larger keys rightward
        self._data[pos:] = np.roll(self._data[pos:], 1)
        self._data[pos] = key + (row,)
        self.length += 1

    def find_pos(self, key, data, exact=False):
        # Return the index of the largest key in data >= given key, data pair.
        # If exact is True, return the index of the given key in
        # data or -1 if the key is not present.
        begin = 0
        end = self.length
        key = key + (data,)

        # search through keys in lexicographic order
        for i in range(self.num_cols + 1):
            key_slice = self._data[self.cols[i]][begin:end]
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
        return self.range(key, key, (True, True))

    def range(self, lower, upper, bounds):
        lower_pos = self.find_pos(lower, 0)
        upper_pos = self.find_pos(upper, 0)
        if lower_pos == self.length:
            return []

        ##TODO: figure out why this tuple thing is necessary
        lower_bound = tuple(self._data[lower_pos])[:-1]
        if not bounds[0] and lower_bound == lower:
            lower_pos += 1 # data[lower_pos] > lower

        # data[lower_pos] >= lower
        # data[upper_pos] >= upper
        if upper_pos < self.length:
            upper_bound = tuple(self._data[upper_pos])[:-1]
            if not bounds[1] and upper_bound == upper:
                upper_pos -= 1 # data[upper_pos] < upper
            elif upper_bound > upper:
                upper_pos -= 1 # data[upper_pos] <= upper
        row_range = self._data['row'][lower_pos:upper_pos + 1]
        return list(np.ravel(list(row_range))) ##TODO: remove list() call

    def remove(self, key, data=None):
        pos = self.find_pos(key, data, exact=True)
        if pos == -1: # key not found
            return False

        # shift larger keys leftward
        self._data[pos:] = np.roll(self._data[pos:], -1)
        self.length -= 1
        return True

    def reorder(self, row):
        self._data['row'][self._data['row'] > row] -= 1

    def replace_rows(self, row_map):
        self._data = self._data[np.array([x in row_map for x in self._data['row']])]
        self._data['row'] = [row_map[x] for x in self._data['row']]

    def items(self):
        array = []
        last_key = None
        for line in self.data:
            line = tuple(line)
            key, row = line[:-1], line[-1]
            if key == last_key:
                array[-1][1].append(row)
            else:
                last_key = key
                array.append((key, [row]))
        return array

    def sort(self):
        return self.data['row']

    @property
    def data(self):
        return self._data[:self.length]

    def __getitem__(self, item):
        # item must be slice
        return SortedArray([], data=self._data[item])

    def __repr__(self):
        return '[' + ', '.join([str(x) for x in self.data]) + ']'

    def __str__(self):
        return repr(self)
