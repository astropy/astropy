import io
from dataclasses import dataclass
from typing import Any

import pytest

from astropy.utils.xml import _iterparser


@dataclass(kw_only=True, slots=True, frozen=True)
class XMLTestData:
    """Strictly typed container for XML parsing validation."""

    xml_bytes: bytes
    expected_events: tuple[tuple[bool, str, Any, tuple[int, int]], ...]


VALID_XML_CASES = [
    pytest.param(
        XMLTestData(
            xml_bytes=b'<?xml version="1.0"?><root><child>text</child></root>',
            expected_events=(
                (True, "xml", {"encoding": "", "version": "1.0"}, (1, 0)),
                (True, "root", {}, (1, 21)),
                (True, "child", {}, (1, 27)),
                (False, "child", "text", (1, 34)),
                (False, "root", "text", (1, 34)),
            ),
        ),
        id="happy_path_standard_xml",
    ),
    pytest.param(
        XMLTestData(
            xml_bytes='<?xml version="1.0" encoding="iso-8859-1"?><root>na\xefve</root>'.encode(
                "iso-8859-1"
            ),
            expected_events=(
                (True, "xml", {"encoding": "iso-8859-1", "version": "1.0"}, (1, 0)),
                (True, "root", {}, (1, 43)),
                (False, "root", "na\xefve", (1, 49)),
            ),
        ),
        id="iso_8859_1_encoded_xml",
    ),
]


def test_c_iterparser_initialization() -> None:
    """Verify that we can spin up the raw C-parser completely on its own."""
    stream = io.BytesIO(b"<dummy></dummy>")
    parser = _iterparser.IterParser(stream.read, buffersize=1024)

    assert type(parser) is _iterparser.IterParser


@pytest.mark.parametrize("test_data", VALID_XML_CASES)
def test_c_iterparser_valid_parsing(test_data: XMLTestData) -> None:
    """Run a 'happy path' check to make sure the C-engine parses valid XML without crashing."""
    stream = io.BytesIO(test_data.xml_bytes)
    parser = _iterparser.IterParser(stream.read, buffersize=1024)

    events = list(parser)

    assert len(events) == len(test_data.expected_events)

    for actual_event, expected_event in zip(
        events, test_data.expected_events, strict=True
    ):
        assert actual_event == expected_event


def test_c_iterparser_malformed_raises_error() -> None:
    """Check that libexpat actively catches the error when we feed it broken XML."""
    # intentionally cutting off the XML stream mid-tag here to make sure Expat raises.
    malformed_bytes = b'<?xml version="1.0"?><root><child>unclosed'
    stream = io.BytesIO(malformed_bytes)
    parser = _iterparser.IterParser(stream.read, buffersize=1024)

    with pytest.raises(ValueError, match="no element found"):
        list(parser)
