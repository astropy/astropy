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

    def __init__(self, lines):
        self.root = None
        for key in lines:
            self.num_cols = len(key)
            break
        self.length = len(lines)
        self.buffer_size = self.length
        self._data = np.zeros(length, dtype=[
            ('key', '{0}O'.format(num_cols)), ('rows', 'O')])
        for i, (key, rows) in enumerate(lines.items()):
            self._data['key'][i] = key
            self._data['rows'][i] = rows

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
        self.length += 1
        

    def find_pos(self, key):
        # returns index of first data key larger than the supplied key,
        # or len(self.data) if no such key exists
        if key > self.data[-1][0]:
            return len(self.data)
        elif key < self.data[0][0]:
            return 0
        return self._find_pos(key, 0, len(self.data))

    def find(self, key):
        pos = self.find_pos(key)
        if pos < len(self.data) and self.data[pos][0] == key:
            return ArraySlot(*self.data[pos])
        return None

    def _find_pos(self, key, low, high):
        if high == low + 1:
            return low
        elif high == low + 2:
            return low if key <= self.data[low][0] else low + 1
        mid = (low + high) / 2
        if key < self.data[mid][0]:
            return self._find_pos(key, low, mid + 1)
        elif key > self.data[mid][0]:
            return self._find_pos(key, mid + 1, high)
        return mid

    def sort(self):
        return self.data

    def __repr__(self):
        return '[' + ', '.join([str(x) for x in self.data]) + ']'

    def __str__(self):
        return repr(self)
