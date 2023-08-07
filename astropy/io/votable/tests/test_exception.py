# Licensed under a 3-clause BSD style license - see LICENSE.rst
import pytest

# LOCAL
from astropy.io.votable import converters, exceptions, tree


def test_reraise():
    def fail():
        raise RuntimeError("This failed")

    with pytest.raises(RuntimeError, match="From here"):
        try:
            fail()
        except RuntimeError as e:
            exceptions.vo_reraise(e, additional="From here")


def test_parse_vowarning():
    config = {"verify": "exception", "filename": "foo.xml"}
    pos = (42, 64)
    with pytest.warns(exceptions.W47) as w:
        field = tree.Field(None, name="c", datatype="char", config=config, pos=pos)
        converters.get_converter(field, config=config, pos=pos)

    parts = exceptions.parse_vowarning(str(w[0].message))

    match = {
        "number": 47,
        "is_exception": False,
        "nchar": 64,
        "warning": "W47",
        "is_something": True,
        "message": "Missing arraysize indicates length 1",
        "doc_url": "io/votable/api_exceptions.html#w47",
        "nline": 42,
        "is_warning": True,
    }
    assert parts == match


def test_suppress_warnings():
    cfg = {}
    warn = exceptions.W01("foo")

    with exceptions.conf.set_temp("max_warnings", 2):
        with pytest.warns(exceptions.W01) as record:
            exceptions._suppressed_warning(warn, cfg)
            assert len(record) == 1
            assert "suppressing" not in str(record[0].message)

        with pytest.warns(exceptions.W01, match="suppressing"):
            exceptions._suppressed_warning(warn, cfg)

        exceptions._suppressed_warning(warn, cfg)

    assert cfg["_warning_counts"][exceptions.W01] == 3
    assert exceptions.conf.max_warnings == 10
