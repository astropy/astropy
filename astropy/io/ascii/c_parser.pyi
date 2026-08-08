from typing import Any, Dict, List, Optional, Tuple, Union

class CParser:
    def __init__(
        self,
        table: Any,
        strip_whitespace_lines: bool,
        strip_whitespace_fields: bool,
        delimiter: str = ...,
        header_start: Optional[int] = ...,
        data_start: Optional[int] = ...,
        comment: str = ...,
        quotechar: str = ...,
        escapechar: str = ...,
        fill_extra_cols: bool = ...,
        **kwargs: Any,
    ) -> None: ...
    
    def read(self, nrows: Optional[int] = ...) -> List[Tuple[Any, ...]]: ...
    
    def setup_tokenizer(self) -> None: ...
