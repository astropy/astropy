from typing import Protocol

from numpy.typing import DTypeLike

from astropy.units import Unit

__all__ = [
    "InfoAttr",
    "AstropyArrayLike",
]


class InfoAttr(Protocol):
    name: str
    format: str
    unit: Unit
    description: str
    meta: dict
    dtype: DTypeLike


class AstropyArrayLike(Protocol):
    # describes astropy ndarray subclasses or ShapedLikeNDArray
    info: InfoAttr
    shape: tuple[int, ...]
