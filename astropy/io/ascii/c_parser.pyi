

class CParser:
    def __init__(
        self,
        table: str | bytes | list[str],
        strip_whitespace_lines: bool,
        strip_whitespace_fields: bool,
        delimiter: str = ",",
        header_start: int | None = None,
        data_start: int | None = None,
        comment: str = "#",
        quotechar: str = '"',
        escapechar: str = "\\",
        fill_extra_cols: bool = False,
        **kwargs: object,
    ) -> None: ...

    def read(self, nrows: int | None = None) -> list[tuple[str | int | float | bool, ...]]: ...
