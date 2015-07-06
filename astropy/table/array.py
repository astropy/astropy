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

    def __init__(self, lines, num_cols=None):
        self.root = None
        if num_cols is not None:
            self.num_cols = num_cols
        else:
            for key in lines:
                self.num_cols = len(key) if isinstance(key, tuple) else 1
                break

        self.length = len(lines)
        self.buffer_size = self.length
        key_cols = [('k{0}'.format(i), 'O') for i in range(self.num_cols)]
        self._data = np.zeros(self.length, dtype=key_cols + [('rows', 'O')])
        for i, (key, rows) in enumerate(lines.items()):
            for j in range(self.num_cols):
                self._data['k{0}'.format(j)][i] = key[j]
            self._data['rows'][i] = rows
        self._data = self._data[np.argsort(self._data)]

    @property
    def data(self):
        return self._data[:self.length]

    def resize(self):
        if self.length == self.buffer_size:
            # resize to 50% larger capacity
            increase = int(0.5 * self.buffer_size)
            if increase == 0: # started out with buffer_size == 1
                increase = 1
            self.buffer_size += increase
            self._data = np.resize(self._data, self.buffer_size)
                                                
    def add(self, key, val):
        self.resize()
        begin = 0
        end = self.length
        self.length += 1

    def find_pos(self, key, exact=False):
        # Return the index of the largest key in data <= given key.
        # If exact is True, return the index of the given key in
        # data or -1 if the key is not present.
        begin = 0
        end = self.length

        # search through keys in lexicographic order
        for i in range(self.num_cols):
            key_slice = self._data['k{0}'.format(i)][begin:end]
            t = np.searchsorted(key_slice, key[i])
            begin += t
            if exact and (begin == self.length or key_slice[t] != key[i]):
                # no match
                return -1
            elif begin == self.length: # greater than all keys
                return begin
            end = begin + np.searchsorted(key_slice, key[i], side='right')

        return begin

    def find(self, key):
        pos = self.find_pos(key, exact=True)
        return [] if pos == -1 else tuple(self._data[pos])[-1]

    def range(self, lower, upper, bounds):
        lower_pos = self.find_pos(lower)
        upper_pos = self.find_pos(upper)
        # data[lower_pos] <= lower
        # data[upper_pos] <= upper
        if not bounds[0] and self._data[lower_pos] == lower:
            lower_pos += 1
        if bounds[1] or self._data[upper_pos] == upper:
            upper_pos += 1
        row_range = self._data['rows'][lower_pos:upper_pos]
        return np.ravel(list(row_range)) ##TODO: remove list() call

    def sort(self):
        return self.data

    def __repr__(self):
        return '[' + ', '.join([str(x) for x in self.data]) + ']'

    def __str__(self):
        return repr(self)
