from typing import Any


class CParser:
    def __init__(
        self,
        table: Any,
        strip_whitespace_lines: bool,
        strip_whitespace_fields: bool,
        delimiter: str | Ellipsis = ...,
        header_start: int | None | Ellipsis = ...,
        data_start: int | None | Ellipsis = ...,
        comment: str | Ellipsis = ...,
        quotechar: str | Ellipsis = ...,
        escapechar: str | Ellipsis = ...,
        fill_extra_cols: bool | Ellipsis = ...,
        **kwargs: Any,
    ) -> None: ...

    def read(self, nrows: int | None | Ellipsis = ...) -> list[tuple[Any, ...]]: ...

    def setup_tokenizer(self) -> None: ...
