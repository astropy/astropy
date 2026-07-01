import collections.abc
from collections.abc import Callable, Iterator
from typing import TypeAlias

_IterparseEvent: TypeAlias = tuple[
    bool, str, collections.abc.Mapping[str, str] | str, tuple[int, int]
]

class IterParser(Iterator[_IterparseEvent]):
    def __init__(
        self, read_func: Callable[[], bytes], buffersize: int = ...
    ) -> None: ...
    def __iter__(self) -> IterParser: ...
    def __next__(self) -> _IterparseEvent: ...
