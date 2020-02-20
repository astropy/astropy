# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
The SCEngine class uses the ``sortedcontainers`` package to implement an
Index engine for Tables.
"""

from collections import OrderedDict
from itertools import starmap

try:
    from sortedcontainers import SortedList
    HAS_SOCO = True
except ImportError:
    HAS_SOCO = False


class Node(object):
    __slots__ = ('key', 'value')

    def __init__(self, key, value):
        self.key = key
        self.value = value

    def __lt__(self, other):
        if other.__class__ is Node:
            return (self.key, self.value) < (other.key, other.value)
        return self.key < other

    def __le__(self, other):
        if other.__class__ is Node:
            return (self.key, self.value) <= (other.key, other.value)
        return self.key <= other

    def __eq__(self, other):
        if other.__class__ is Node:
            return (self.key, self.value) == (other.key, other.value)
        return self.key == other

    def __ne__(self, other):
        if other.__class__ is Node:
            return (self.key, self.value) != (other.key, other.value)
        return self.key != other

    def __gt__(self, other):
        if other.__class__ is Node:
            return (self.key, self.value) > (other.key, other.value)
        return self.key > other

    def __ge__(self, other):
        if other.__class__ is Node:
            return (self.key, self.value) >= (other.key, other.value)
        return self.key >= other

    __hash__ = None

    def __repr__(self):
        return f'Node({self.key!r}, {self.value!r})'


class SCEngine:
    '''
    Fast tree-based implementation for indexing, using the
    ``sortedcontainers`` package.

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
        node_keys = map(tuple, data)
        self._nodes = SortedList(starmap(Node, zip(node_keys, row_index)))
        self._unique = unique

    def add(self, key, value):
        '''
        Add a key, value pair.
        '''
        if self._unique and (key in self._nodes):
            message = f'duplicate {key:!r} in unique index'
            raise ValueError(message)
        self._nodes.add(Node(key, value))

    def find(self, key):
        '''
        Find rows corresponding to the given key.
        '''
        return [node.value for node in self._nodes.irange(key, key)]

    def remove(self, key, data=None):
        '''
        Remove data from the given key.
        '''
        if data is not None:
            item = Node(key, data)
            try:
                self._nodes.remove(item)
            except ValueError:
                return False
            return True
        items = list(self._nodes.irange(key, key))
        for item in items:
            self._nodes.remove(item)
        return bool(items)

    def shift_left(self, row):
        '''
        Decrement rows larger than the given row.
        '''
        for node in self._nodes:
            if node.value > row:
                node.value -= 1

    def shift_right(self, row):
        '''
        Increment rows greater than or equal to the given row.
        '''
        for node in self._nodes:
            if node.value >= row:
                node.value += 1

    def items(self):
        '''
        Return a list of key, data tuples.
        '''
        result = OrderedDict()
        for node in self._nodes:
            if node.key in result:
                result[node.key].append(node.value)
            else:
                result[node.key] = [node.value]
        return result.items()

    def sort(self):
        '''
        Make row order align with key order.
        '''
        for index, node in enumerate(self._nodes):
            node.value = index

    def sorted_data(self):
        '''
        Return a list of rows in order sorted by key.
        '''
        return [node.value for node in self._nodes]

    def range(self, lower, upper, bounds=(True, True)):
        '''
        Return row values in the given range.
        '''
        iterator = self._nodes.irange(lower, upper, bounds)
        return [node.value for node in iterator]

    def replace_rows(self, row_map):
        '''
        Replace rows with the values in row_map.
        '''
        nodes = [node for node in self._nodes if node.value in row_map]
        for node in nodes:
            node.value = row_map[node.value]
        self._nodes.clear()
        self._nodes.update(nodes)

    def __repr__(self):
        return '{!r}'.format(list(self._nodes))
