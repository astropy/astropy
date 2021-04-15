# Licensed under a 3-clause BSD style license - see LICENSE.rst

import re
import itertools

import numpy as np

from astropy.utils.data_info import ParentDtypeInfo, MULTIDIM_AS_JSON


class NdarrayMixinInfo(ParentDtypeInfo):
    _represent_as_dict_primary_data = 'data'

    def _represent_as_dict(self):
        """Represent Column as a dict that can be serialized.

        This handles the special case of serializing an N-dimensional column
        (for N > 1) by optionally turning it into the required number of 1-d
        arrays (e.g. for ECSV or None context).

        For instance a Column with shape (2, 2, 20) will turn into four columns
        (each length 20) with "data attribute" names 0_0, 0_1, 1_0, and 1_1.
        These then get constructed back into the same (2, 2, 20) Column.
        """
        col = self._parent
        out = {'data': col.view(np.ndarray)}

        return out

    def _construct_from_dict(self, map):
        """Construct Column from ``map``.

        Normally this is done by the default constuctor in DataInfo.

        However as a special case if there are any keys in ``map`` that begin
        with a digit, this must correspond to an N-d column that was serialized
        through ECSV. (No actual attribute names can begin with a digit).
        """
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
