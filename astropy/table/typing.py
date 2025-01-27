from typing import Protocol

from astropy.utils.data_info import MixinInfo

__all__ = [
    "MixinColumn",
]


class MixinColumn(Protocol):
    info: MixinInfo
    shape: tuple[int, ...]

    def __getitem__(self, key): ...
    def __len__(self) -> int: ...
