# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
A module containing specialized collection classes.
"""


class HomogeneousList(list):
    """
    A subclass of list that contains only elements of a given type or
    types.  If an item that is not of the specified type is added to
    the list, a `TypeError` is raised.
    """

    def __init__(self, types, values=[]):
        """
        Parameters
        ----------
        types : sequence of types
            The types to accept.

        values : sequence, optional
            An initial set of values.
        """
        self._types = types
        super().__init__()
        self.extend(values)

    def _assert(self, x):
        if not isinstance(x, self._types):
            raise TypeError(
                f"homogeneous list must contain only objects of type '{self._types}'"
            )

    def __iadd__(self, other):
        self.extend(other)
        return self

    def __setitem__(self, idx, value):
        if isinstance(idx, slice):
            value = list(value)
            for item in value:
                self._assert(item)
        else:
            self._assert(value)
        return super().__setitem__(idx, value)

    def append(self, x):
        self._assert(x)
        return super().append(x)

    def insert(self, i, x):
        self._assert(x)
        return super().insert(i, x)

    def extend(self, x):
        for item in x:
            self._assert(item)
            super().append(item)
