import contextlib
import io
import tempfile
from pathlib import Path

import numpy as np
import pytest

from astropy.table import Table
from astropy.utils.compat.optional_deps import HAS_PYARROW

if not HAS_PYARROW:
    pytest.skip("pyarrow is not available")


# pytest fixture that returns a sipmle test table like exp
@pytest.fixture
def tbl_simple():
    return Table(
        rows=[
            [1, 2.0, "foo"],
            [4, 5.5, "bår"],
            [7, 8.0, "bazœo"],
        ],
        names=["a", "b", "ç"],
    )


@pytest.fixture
def tbl_simple_masked():
    return Table(
        rows=[
            [1, 2.0, np.ma.masked],
            [4, np.ma.masked, "bår"],
            [np.ma.masked, 8.0, "bazœo"],
        ],
        names=["a", "b", "ç"],
    )


def convert_table_to_text(tbl, delimiter=",", **kwargs):
    """Convert a table to text using the io.ascii CSV format.

    This function is used to convert a table to text for testing purposes.
    It uses the `ascii.csv` format to write the table to a string buffer.
    """
    text_io = io.StringIO()
    tbl.write(text_io, format="ascii.csv", delimiter=delimiter, **kwargs)
    return text_io.getvalue()


@contextlib.contextmanager
def get_input_file(text: str, input_type: str, encoding: str | None = None):
    encode_kwargs = {"encoding": encoding} if encoding else {}
    if input_type in ("str", "path"):
        with tempfile.NamedTemporaryFile(mode="w", **encode_kwargs) as f:
            f.write(text)
            f.flush()
            yield f.name if input_type == "str" else Path(f.name)
    elif input_type == "bytesio":
        yield io.BytesIO(text.encode(**encode_kwargs))
    else:
        raise ValueError(f"Unknown input_type: {input_type}")


def table_read_csv(
    text: str, input_type: str = "bytesio", encoding: str | None = None, **kwargs
):
    """Read ``text`` using ``Table.read`` with format="pyarrow.csv".

    The ``input_type`` parameter determines how the text is passed to ``Table.read``.
    For "str" or "path", a named temporary file is created and that name ("str") or path
    ("path") is passed to ``Table.read``. For "bytesio", a ``BytesIO`` object is created
    and passed to ``Table.read``. The text is encoded using the specified encoding (if
    any) before being passed to ``BytesIO``. The default is "bytesio".

    The ``encoding`` parameter is used only if ``input_type`` is "bytesio". If
    ``encoding`` is None, the text is encoded using the default encoding.

    Parameters
    ----------
    text : str
        The text to read.
    input_type : str
        The type of input for Table.read(). One of "str", "path", or "bytesio".
    **kwargs : dict
        Additional keyword arguments to pass to Table.read().
    """
    if encoding is not None:
        kwargs["encoding"] = encoding

    with get_input_file(text, input_type, encoding) as input_file:
        out = Table.read(input_file, format="pyarrow.csv", **kwargs)
    return out


@pytest.mark.parametrize("input_type", ["str", "path", "bytesio"])
@pytest.mark.parametrize("encoding", [None, "utf-8", "utf-16"])
def test_read_tbl_simple_input_type_encoding(input_type, tbl_simple, encoding):
    """Test reading a simple CSV file with different input types."""
    text = convert_table_to_text(tbl_simple)
    out = table_read_csv(text, input_type, encoding)
    assert np.all(out == tbl_simple)


@pytest.mark.parametrize("encoding", [None, "utf-16"])
def test_read_tbl_masked_input_type_encoding(tbl_simple_masked, encoding):
    """Test reading a simple CSV file with different input types."""
    text = convert_table_to_text(tbl_simple_masked)
    out = table_read_csv(text, "bytesio", encoding)
    assert np.all(out == tbl_simple_masked)
