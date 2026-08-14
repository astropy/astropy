from typing import Any


class CParser:
    def __init__(
        self,
        table: Any,
        strip_whitespace_lines: bool,
        strip_whitespace_fields: bool,
        delimiter: str = ...,
        header_start: int | None = ...,
        data_start: int | None = ...,
        comment: str = ...,
        quotechar: str = ...,
        escapechar: str = ...,
        fill_extra_cols: bool = ...,
        **kwargs: Any,
    ) -> None: ...

    def read(self, nrows: int | None = ...) -> list[tuple[Any, ...]]: ...

    def setup_tokenizer(self) -> None: ...
