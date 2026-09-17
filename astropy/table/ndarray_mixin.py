# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import annotations

from typing import TYPE_CHECKING, Any, Self

import numpy as np

from astropy.utils.data_info import ParentDtypeInfo

if TYPE_CHECKING:
    import numpy.typing as npt


class NdarrayMixinInfo(ParentDtypeInfo):
    _represent_as_dict_primary_data = "data"

    def _represent_as_dict(self) -> dict[str, Any]:
        """Represent Column as a dict that can be serialized."""
        col = self._parent
        out = {"data": col.view(np.ndarray)}
        return out

    def _construct_from_dict(self, map: dict[str, Any]) -> NdarrayMixin:
        """Construct Column from ``map``."""
        data = map.pop("data")
        out = self._parent_cls(data, **map)
        return out


class NdarrayMixin(np.ndarray):
    """
    Mixin column class to allow storage of arbitrary numpy
    ndarrays within a Table.  This is a subclass of numpy.ndarray
    and has the same initialization options as ``np.array()``.
    """

    info = NdarrayMixinInfo()

    def __new__(cls, obj: npt.ArrayLike, *args: Any, **kwargs) -> Self:
        self = np.array(obj, *args, **kwargs).view(cls)
        if "info" in getattr(obj, "__dict__", ()):
            self.info = obj.info
        return self

    def __array_finalize__(self, obj: np.ndarray | None) -> None:
        if obj is None:
            return

        super().__array_finalize__(obj)

        # Self was created from template (e.g. obj[slice] or (obj * 2))
        # or viewcast e.g. obj.view(Column).  In either case we want to
        # init Column attributes for self from obj if possible.
        if "info" in getattr(obj, "__dict__", ()):
            self.info = obj.info

    def __reduce__(self) -> tuple[Any, ...]:
        # patch to pickle NdArrayMixin objects (ndarray subclasses), see
        # http://www.mail-archive.com/numpy-discussion@scipy.org/msg02446.html

        object_state = list(super().__reduce__())
        object_state[2] = (object_state[2], self.__dict__)
        return tuple(object_state)

    def __setstate__(self, state: tuple[Any, ...]) -> None:
        # patch to unpickle NdarrayMixin objects (ndarray subclasses), see
        # http://www.mail-archive.com/numpy-discussion@scipy.org/msg02446.html

        nd_state, own_state = state
        super().__setstate__(nd_state)
        self.__dict__.update(own_state)
