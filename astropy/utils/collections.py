# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
A module containing specialized collection classes.
"""


class HomogeneousList(list):
    """
    A subclass of list that contains only elements of a given type or
    types.  If an item that is not of the specified type is added to
    the list, a `~.exceptions.TypeError` is raised.
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
        list.__init__(self, values)

    def _assert(self, x):
        if not isinstance(x, self._types):
            raise TypeError(
                "homogeneous list must contain only objects of type '%s'" %
                (self._types))

    def __iadd__(self, other):
        for x in other:
            self._assert(x)
        return list.__iadd__(self, other)

    def __setitem__(self, x):
        self._assert(x)
        return list.__setitem__(self, x)

    def append(self, x):
        self._assert(x)
        return list.append(self, x)

    def insert(self, i, x):
        self._assert(x)
        return list.insert(self, i, x)

    def extend(self, x):
        for item in x:
            self._assert(item)
        return list.extend(self, x)
