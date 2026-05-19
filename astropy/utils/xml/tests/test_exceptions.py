# Licensed under a 3-clause BSD style license - see LICENSE.rst

import warnings

import pytest

from astropy.utils.xml.exceptions import (
    XMLWarning,
    reraise,
    warn,
    warn_or_raise,
)


class XMLTestWarning(XMLWarning):
    message_template = "Test {}"


def test_warn_emits_warning():
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter("always")
        warn(XMLTestWarning, ("message",), config={"verify": "warn"}, pos=(1, 2))

    assert len(rec) == 1
    assert isinstance(rec[0].message, XMLTestWarning)
    assert "Test message" in str(rec[0].message)


def test_warn_or_raise_raises():
    with pytest.raises(XMLTestWarning):
        warn_or_raise(
            XMLTestWarning,
            args=("message",),
            config={"verify": "exception"},
            pos=(1, 2),
        )


def test_reraise_preserves_context():
    try:
        raise RuntimeError("original")
    except RuntimeError as exc:
        with pytest.raises(RuntimeError, match="original From here"):
            reraise(exc, additional="From here")
