import collections.abc
from collections.abc import Callable, Iterator
from typing import AnyStr, Protocol, TypeAlias, final

_IterparseEvent: TypeAlias = tuple[
    bool, str, collections.abc.Mapping[str, str] | str, tuple[int, int]
]

# The C-engine always passes 'buffersize' as a single positional argument to the read function (via Py_BuildValue("(n)", buffersize)).
class _Readable(Protocol):
    def read(self, size: int, /) -> bytes: ...

@final
class IterParser(Iterator[_IterparseEvent]):
    # The Callable fallback explicitly matches the _Readable.read interface, as the C-engine passes exactly one integer (buffersize) to both.
    def __init__(
        self, fd: _Readable | Callable[[int], bytes], buffersize: int = 16384, /
    ) -> None: ...
    def __iter__(self) -> IterParser: ...
    def __next__(self) -> _IterparseEvent: ...

def escape_xml(input: AnyStr, /) -> AnyStr: ...
def escape_xml_cdata(input: AnyStr, /) -> AnyStr: ...
