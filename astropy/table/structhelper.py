import numpy as np
import numpy.ma as ma


def _append_field(sta, data, dtype=None, position=None, masked=False):

    newdtype = sta.dtype.descr
    if np.equal(position, None):
        newdtype.append(dtype)
    else:
        newdtype.insert(position, dtype)
    newdtype = np.dtype(newdtype)

    if masked:
        newsta = ma.empty(sta.shape, dtype=newdtype)
    else:
        newsta = np.empty(sta.shape, dtype=newdtype)

    for field in sta.dtype.fields:
        newsta[field] = sta[field]
        if masked:
            newsta[field].set_fill_value(sta[field].fill_value)

    newsta[dtype[0]] = data
    if masked:
        newsta[dtype[0]].set_fill_value(data.fill_value)

    return newsta


def _drop_fields(sta, names, masked=False):

    names = set(names)

    newdtype = np.dtype([(name, sta.dtype[name]) for name in sta.dtype.names
                       if name not in names])

    if newdtype:
        if masked:
            newsta = ma.empty(sta.shape, dtype=newdtype)
        else:
            newsta = np.empty(sta.shape, dtype=newdtype)
    else:
        return None

    for field in newdtype.fields:
        newsta[field] = sta[field]
        if masked:
            newsta[field].set_fill_value(sta[field].fill_value)

    return newsta
