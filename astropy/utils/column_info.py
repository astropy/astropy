from copy import deepcopy
import weakref
import numpy as np

COLUMN_ATTRS = set(['name', 'unit', 'dtype', 'format', 'description', 'meta', 'parent_table'])

class ColumnInfo(object):
    cols_from_parent = set()

    def __init__(self, parent_col=None):
        self._attrs = dict((attr, None) for attr in COLUMN_ATTRS)
        if parent_col is not None:
            self._parent_col = weakref.ref(parent_col)

    def __getstate__(self):
        return self._attrs

    def __setstate__(self, state):
        self._attrs = state

    def __getattr__(self, attr):
        if attr.startswith('_'):
            return super(ColumnInfo, self).__getattribute__(attr)

        if attr in self.cols_from_parent:
            return getattr(self._parent_col(), attr)

        try:
            value = self._attrs[attr]
        except KeyError:
            super(ColumnInfo, self).__getattribute__(attr)  # Generate AttributeError

        # Weak ref for parent table
        if attr == 'parent_table' and callable(value):
            value = value()

        # Mixins have a default dtype of Object if nothing else was set
        if attr == 'dtype' and value is None:
            value = np.dtype('O')

        return value

    def __setattr__(self, attr, value):
        if attr in self.cols_from_parent:
            setattr(self._parent_col(), attr, value)
            return

        if attr.startswith('_'):
            super(ColumnInfo, self).__setattr__(attr, value)
            return

        if attr not in COLUMN_ATTRS:
            raise AttributeError("attribute must be one of {0}".format(COLUMN_ATTRS))

        if attr == 'parent_table':
            value = None if value is None else weakref.ref(value)

        self._attrs[attr] = value

    def copy(self):
        out = self.__class__()
        for attr in COLUMN_ATTRS - self.cols_from_parent:
            setattr(out, attr, deepcopy(getattr(self, attr)))

        return out

