"""Typing for `astropy.table`."""

from __future__ import annotations

from typing import (
    TYPE_CHECKING,
    Protocol,
    Self,
    SupportsIndex,
    TypeVar,
    overload,
)

if TYPE_CHECKING:
    import numpy as np
    from numpy.typing import NDArray

    from astropy.utils.data_info import MixinInfo

__all__ = [
    "MixinColumn",
]


T = TypeVar("T")


class MixinColumn(Protocol[T]):
    info: MixinInfo
    shape: tuple[int, ...]

    @overload
    def __getitem__(self, key: SupportsIndex) -> T: ...
    @overload
    def __getitem__(
        self, key: slice | NDArray[np.bool] | NDArray[np.integer]
    ) -> Self: ...
    def __getitem__(self, key): ...
    def __len__(self) -> int:
        """Return the length of the mixin column.

        The MixinColumn provides a default implementation for the length based on the shape attribute.

        Examples
        --------
        >>> from astropy.table.typing import MixinColumn

        >>> class Obj(MixinColumn):
        ...     shape = (5, 2, 3)

        >>> obj = Obj()
        >>> len(obj)
        5

        """
        return self.shape[0]
