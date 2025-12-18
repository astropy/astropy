# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
The SCEngine class uses the ``sortedcontainers`` package to implement an
Index engine for Tables.
"""

from collections import OrderedDict
from collections.abc import Hashable, Mapping, Sequence
from numbers import Integral

from astropy.utils.compat.optional_deps import HAS_SORTEDCONTAINERS

if HAS_SORTEDCONTAINERS:
    from sortedcontainers import SortedList


class Node:
    __slots__ = ("key", "value")

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
        return f"Node({self.key!r}, {self.value!r})"


class SCEngine:
    """
    Fast tree-based implementation for indexing, using the
    ``sortedcontainers`` package.

    Parameters
    ----------
    data : Table
        Sorted columns of the original table
    row_index : Column object
        Row numbers corresponding to data columns
    unique : bool
        Whether the values of the index must be unique.
        Defaults to False.
    """

    def __init__(self, data, row_index, unique=False):
        if not HAS_SORTEDCONTAINERS:
            raise ImportError("sortedcontainers is needed for using SCEngine")

        node_keys = map(tuple, data)
        self._nodes = SortedList(map(Node, node_keys, row_index))
        self._unique = unique

    def add(self, key: tuple, row: int) -> None:
        """
        Add a key, value pair.
        """
        if self._unique and (key in self._nodes):
            message = f"duplicate {key!r} in unique index"
            raise ValueError(message)
        self._nodes.add(Node(key, row))

    def find(self, key: tuple) -> Sequence[Integral]:
        """
        Find rows corresponding to the given key.
        """
        return [node.value for node in self._nodes.irange(key, key)]

    def remove(self, key: tuple, data: int) -> bool:
        """
        Remove data from the given key.
        """
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

    def shift_left(self, row: int) -> None:
        """
        Decrement rows larger than the given row.
        """
        for node in self._nodes:
            if node.value > row:
                node.value -= 1

    def shift_right(self, row: int) -> None:
        """
        Increment rows greater than or equal to the given row.
        """
        for node in self._nodes:
            if node.value >= row:
                node.value += 1

    def items(self) -> list[tuple[Hashable, list[Integral]]]:
        """
        Return a list of key, data tuples.
        """
        result = OrderedDict()
        for node in self._nodes:
            if node.key in result:
                result[node.key].append(node.value)
            else:
                result[node.key] = [node.value]
        return list(result.items())

    def sort(self) -> None:
        """
        Make row order align with key order.
        """
        for index, node in enumerate(self._nodes):
            node.value = index

    def sorted_data(self):
        """
        Return a list of rows in order sorted by key.
        """
        return [node.value for node in self._nodes]

    def range(
        self,
        lower: tuple[Hashable, ...] | None,
        upper: tuple[Hashable, ...] | None,
        bounds: tuple[bool, bool],
    ) -> list[int]:
        """
        Return row values in the given range.

        A ``None`` value for ``lower`` or ``upper`` corresponds to no limit just like
        slicing.
        """
        iterator = self._nodes.irange(lower, upper, bounds)
        return [node.value for node in iterator]

    def replace_rows(self, row_map: "Mapping[int, int]") -> None:
        """
        Replace rows with the values in row_map.
        """
        nodes = [node for node in self._nodes if node.value in row_map]
        for node in nodes:
            node.value = row_map[node.value]
        self._nodes.clear()
        self._nodes.update(nodes)

    def __repr__(self):
        if len(self._nodes) > 6:
            nodes = list(self._nodes[:3]) + ["..."] + list(self._nodes[-3:])
        else:
            nodes = self._nodes
        nodes_str = ", ".join(str(node) for node in nodes)
        return f"<{self.__class__.__name__} nodes={nodes_str}>"
