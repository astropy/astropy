# Licensed under a 3-clause BSD style license - see LICENSE.rst

import re
import itertools

import numpy as np

from astropy.utils.data_info import ParentDtypeInfo


class NdarrayMixinInfo(ParentDtypeInfo):
    _represent_as_dict_primary_data = 'data'

    def _represent_as_dict(self):
        """Represent Column as a dict that can be serialized.

        This handles the special case of serializing in ECSV an N-dimensional
        column (for N > 1) by turning it into the required number of 1-d arrays.

        For instance a Column with shape (2, 2, 20) will turn into four columns
        (each length 20) with "data attribute" names 0_0, 0_1, 1_0, and 1_1.
        These then get constructed back into the same (2, 2, 20) Column.
        """
        col = self._parent

        if col.ndim > 1 and self._serialize_context == 'ecsv':
            # N-d column case for output to ECSV
            out = {}
            # Make cartesian product of all dimensions except length (first axis)
            ranges = [range(ii) for ii in col.shape[1:]]
            for dims in itertools.product(*ranges):
                name = '_'.join(str(ii) for ii in dims)
                out[name] = col[(...,) + dims].view(np.ndarray)
        else:
            out = {'data': col.view(np.ndarray)}

        return out

    def _construct_from_dict(self, map):
        """Construct Column from ``map``.

        Normally this is done by the default constuctor in DataInfo.

        However as a special case if there are any keys in ``map`` that begin
        with a digit, this must correspond to an N-d column that was serialized
        through ECSV. (No actual attribute names can begin with a digit).

        In this case
        """
        idxs_map = {}
        # Look for keys like 0_0 or 2_1 corresponding to serializing an N-d
        # Column. For instance a Column with shape (2, 2, 20) will turn into
        # four columns (each length 20) with "data attribute" names 0_0, 0_1,
        # 1_0, and 1_1.
        for key in map:
            if re.match(r'\d', key):
                idxs = tuple(int(ii) for ii in key.split('_'))
                idxs_map[key] = idxs

        if idxs_map:
            # Found appropriately-named keys, so this is an N-d column. Get the
            # length of the first Column in the set of N-d columns (they must
            # all be the same length).
            key0 = next(iter(map))
            data0 = map[key0]
            length = len(data0)

            # Get the shape by using finding max of the index tuples (tuples do
            # indeed sort correctly).
            max_idxs = max(idxs_map.values())
            shape = (length,) + tuple(ii + 1 for ii in max_idxs)

            # Make the N-d array and fill accordingly
            data = np.empty(shape, dtype=data0.dtype)
            for key, idxs in idxs_map.items():
                data[(...,) + idxs] = map.pop(key)
        else:
            data = map.pop('data')

        out = self._parent_cls(data, **map)

        return out


class NdarrayMixin(np.ndarray):
    """
    Mixin column class to allow storage of arbitrary numpy
    ndarrays within a Table.  This is a subclass of numpy.ndarray
    and has the same initialization options as ndarray().
    """
    __module__ = 'astropy.table.table'
    info = NdarrayMixinInfo()

    def __new__(cls, obj, *args, **kwargs):
        self = np.array(obj, *args, **kwargs).view(cls)
        if 'info' in getattr(obj, '__dict__', ()):
            self.info = obj.info
        return self

    def __array_finalize__(self, obj):
        if obj is None:
            return

        if callable(super().__array_finalize__):
            super().__array_finalize__(obj)

        # Self was created from template (e.g. obj[slice] or (obj * 2))
        # or viewcast e.g. obj.view(Column).  In either case we want to
        # init Column attributes for self from obj if possible.
        if 'info' in getattr(obj, '__dict__', ()):
            self.info = obj.info

    def __reduce__(self):
        # patch to pickle Quantity objects (ndarray subclasses), see
        # http://www.mail-archive.com/numpy-discussion@scipy.org/msg02446.html

        object_state = list(super().__reduce__())
        object_state[2] = (object_state[2], self.__dict__)
        return tuple(object_state)

    def __setstate__(self, state):
        # patch to unpickle NdarrayMixin objects (ndarray subclasses), see
        # http://www.mail-archive.com/numpy-discussion@scipy.org/msg02446.html

        nd_state, own_state = state
        super().__setstate__(nd_state)
        self.__dict__.update(own_state)
