# Licensed under a 3-clause BSD style license - see LICENSE.rst

import io

from astropy.utils.xml.io import convert_to_writable_filelike


def test_convert_to_writable_filelike_path(tmp_path):
    path = tmp_path / "test.xml"
    with convert_to_writable_filelike(str(path)) as fh:
        fh.write("hello")

    assert path.read_text(encoding="utf8") == "hello"


def test_convert_to_writable_filelike_fileobj():
    buffer = io.BytesIO()
    with convert_to_writable_filelike(buffer) as fh:
        fh.write("hello")

    assert buffer.getvalue() == b"hello"


def test_convert_to_writable_filelike_compressed(tmp_path):
    path = tmp_path / "test.xml.gz"
    with convert_to_writable_filelike(str(path)) as fh:
        fh.write("hello")

    assert path.stat().st_size > 0
    assert path.suffix == ".gz"
