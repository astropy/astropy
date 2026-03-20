# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest

from astropy import conf
from astropy.tests.helper import assert_follows_unicode_guidelines


class RoundtripBase:
    def __init__(self, string):
        self.string = string

    def __eq__(self, other):
        return self.string == other.string

    def __repr__(self):
        return f"{type(self).__name__}({self.string!r})"

    def __str__(self):
        return self.string


# The following 7 tests trigger the 7 asserts in assert_follows_unicode_guidelines,
# and they are in the same order as the asserts they check.
def test_non_ascii_format():
    class NonAsciiFormat:
        def __format__(self, fmt):
            return "õ"

    with pytest.raises(AssertionError):
        assert_follows_unicode_guidelines(NonAsciiFormat())


def test_non_ascii_str():
    class NonAsciiStr:
        def __str__(self):
            return "õ"

    with pytest.raises(AssertionError):
        assert_follows_unicode_guidelines(NonAsciiStr())


def test_non_ascii_repr():
    class NonAsciiRepr:
        def __repr__(self):
            return "õ"

    with pytest.raises(AssertionError):
        assert_follows_unicode_guidelines(NonAsciiRepr())


def test_str_does_not_roundtrip():
    class NoRoundtripStr:
        def __init__(self, string):
            self.string = string

    with pytest.raises(AssertionError):
        assert_follows_unicode_guidelines(NoRoundtripStr("bad"), locals())


def test_repr_does_not_roundtrip():
    class NoRoundtripRepr(RoundtripBase):
        def __repr__(self):
            return "NoRoundtripRepr('')"

    with pytest.raises(AssertionError):
        assert_follows_unicode_guidelines(NoRoundtripRepr("bad"), locals())


def test_configurable_repr():
    class ConfigurableRepr:
        def __repr__(self):
            return "õ" if conf.unicode_output else "6"

    with pytest.raises(AssertionError):
        assert_follows_unicode_guidelines(ConfigurableRepr())


def test_unicode_does_not_roundtrip():
    class NoRoundtripConfigurableStr(RoundtripBase):
        def __str__(self):
            return "õ" if conf.unicode_output else self.string

    with pytest.raises(AssertionError):
        assert_follows_unicode_guidelines(NoRoundtripConfigurableStr("bad"), locals())


def test_does_roundtrip():
    assert_follows_unicode_guidelines(RoundtripBase("good"), globals())


class ConfigurableFormat:
    def __format__(self, format_spec):
        return "õ" if conf.unicode_output else "6"


class ConfigurableStr:
    def __str__(self):
        return "õ" if conf.unicode_output else "6"


class FlexibleFormat:
    def __format__(self, format_spec):
        return "õ" if format_spec else "6"


@pytest.mark.parametrize(
    "test_class",
    [
        pytest.param(
            ConfigurableFormat, id="__format__ depends on conf.unicode_output"
        ),
        pytest.param(ConfigurableStr, id="__str__ depends on conf.unicode_output"),
        pytest.param(FlexibleFormat, id="__format__ depends on format_spec"),
    ],
)
def test_allowed_non_ascii(test_class):
    assert_follows_unicode_guidelines(test_class())
