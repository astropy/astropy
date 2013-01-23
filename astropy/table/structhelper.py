# Licensed under a 3-clause BSD style license - see LICENSE.rst
import numpy as np


def _append_field(table, data, dtype=None, position=None):

    newdtype = table.dtype.descr
    if position is None:
        position = len(newdtype)

    newdtype.insert(position, dtype)
    newtable = np.empty(table.shape, dtype=newdtype)

    for field in table.dtype.fields:
        newtable[field] = table[field]

    newtable[dtype[0]] = data

    return newtable


def _drop_fields(table, names):

    names = set(names)
    newdtype = [(name, table.dtype[name]) for name in table.dtype.names
                if name not in names]
    newdtype = np.dtype(newdtype)

    if newdtype:
        newtable = np.empty(table.shape, dtype=newdtype)
    else:
        return None

    for field in newdtype.fields:
        newtable[field] = table[field]

    return newtable
