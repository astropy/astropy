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
        key_cols = [('k{0}'.format(i), 'O') for i in range(self.num_cols)]
        self._data = np.zeros(self.length, dtype=key_cols + [('rows', 'O')])
        for i, (key, rows) in enumerate(lines.items()):
            for j in range(self.num_cols):
                self._data['k{0}'.format(j)][i] = key[j]
            self._data['rows'][i] = rows
        self._data = self._data[np.argsort(self._data)]

    def resize(self):
        buffer_size = len(self._data)
        if self.length == buffer_size:
            # resize to 50% larger capacity
            increase = int(0.5 * buffer_size)
            if increase == 0: # started out with buffer_size == 1
                increase = 1
            buffer_size += increase
            self._data = np.resize(self._data, buffer_size)
                                                
    def add(self, key, val):
        pos = self.find_pos(key) # first >= key
        if pos < self.length and tuple(self._data[pos])[:-1] == key:
            # add to existing entry
            rows = self._data['rows'][pos]
            # make sure row lists remain sorted
            rows[:] = np.insert(rows, np.searchsorted(rows, val), val)
        else:
            self.resize()
            self.length += 1
            # shift larger keys rightward
            self._data[pos:] = np.roll(self._data[pos:], 1)
            self._data[pos] = key + ([val],)

    def find_pos(self, key, exact=False):
        # Return the index of the largest key in data >= given key.
        # If exact is True, return the index of the given key in
        # data or -1 if the key is not present.
        begin = 0
        end = self.length

        # search through keys in lexicographic order
        for i in range(self.num_cols):
            key_slice = self._data['k{0}'.format(i)][begin:end]
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
        pos = self.find_pos(key, exact=True)
        return [] if pos == -1 else tuple(self._data[pos])[-1]

    def range(self, lower, upper, bounds):
        lower_pos = self.find_pos(lower)
        upper_pos = self.find_pos(upper)
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
        row_range = self._data['rows'][lower_pos:upper_pos + 1]
        return np.ravel(list(row_range)) ##TODO: remove list() call

    def remove(self, key, data=None):
        pos = self.find_pos(key, exact=True)
        if pos == -1: # key not found
            return False

        if data is not None:
            node = self._data['rows'][pos]
            if data not in node:
                raise ValueError("Data does not belong to correct node")
            elif len(node) > 1: # other rows should be retained
                node.remove(data)
                return True

        # shift larger keys leftward
        self._data[pos:] = np.roll(self._data[pos:], -1)
        self.length -= 1
        return True

    def reorder(self, row):
        for data in self._data['rows']:
            data[:] = [x - 1 if x > row else x for x in data]

    def replace_rows(self, row_map):
        for data in self._data['rows']:
            data[:] = [row_map[x] for x in data if x in row_map]

    def items(self):
        ##TODO: figure this tuple thing out
        return [(tuple(k)[:-1], tuple(k)[-1]) for k in self.data]

    def sort(self):
        rows = [tuple(k)[-1] for k in self.data]
        return [x for sublist in rows for x in sublist]

    @property
    def data(self):
        return self._data[:self.length]

    def __repr__(self):
        return '[' + ', '.join([str(x) for x in self.data]) + ']'

    def __str__(self):
        return repr(self)
