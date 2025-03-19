# Licensed under a 3-clause BSD style license - see LICENSE.rst


import io
import sys

import pytest

from astropy import units as u
from astropy.utils import console

from . import test_progress_bar_func


class FakeTTY(io.StringIO):
    """IOStream that fakes a TTY; provide an encoding to emulate an output
    stream with a specific encoding.
    """

    def __new__(cls, encoding=None):
        # Return a new subclass of FakeTTY with the requested encoding
        if encoding is None:
            return super().__new__(cls)

        cls = type(encoding.title() + cls.__name__, (cls,), {"encoding": encoding})

        return cls.__new__(cls)

    def __init__(self, encoding=None):
        super().__init__()

    def write(self, s):
        if isinstance(s, bytes):
            # Just allow this case to work
            s = s.decode("latin-1")
        elif self.encoding is not None:
            s.encode(self.encoding)

        return super().write(s)

    def isatty(self):
        return True


def test_fake_tty():
    # First test without a specified encoding; we should be able to write
    # arbitrary unicode strings
    f1 = FakeTTY()
    assert f1.isatty()
    f1.write("☃")
    assert f1.getvalue() == "☃"

    # Now test an ASCII-only TTY--it should raise a UnicodeEncodeError when
    # trying to write a string containing non-ASCII characters
    f2 = FakeTTY("ascii")
    assert f2.isatty()
    assert f2.__class__.__name__ == "AsciiFakeTTY"
    assert pytest.raises(UnicodeEncodeError, f2.write, "☃")
    assert f2.getvalue() == ""


@pytest.mark.skipif(sys.platform.startswith("win"), reason="Cannot test on Windows")
def test_color_text():
    assert console._color_text("foo", "green") == "\033[0;32mfoo\033[0m"


def test_color_print():
    # This stuff is hard to test, at least smoke test it
    console.color_print("foo", "green")

    console.color_print("foo", "green", "bar", "red")


def test_color_print2():
    # Test that this automatically detects that io.StringIO is
    # not a tty
    stream = io.StringIO()
    console.color_print("foo", "green", file=stream)
    assert stream.getvalue() == "foo\n"

    stream = io.StringIO()
    console.color_print("foo", "green", "bar", "red", "baz", file=stream)
    assert stream.getvalue() == "foobarbaz\n"


@pytest.mark.skipif(sys.platform.startswith("win"), reason="Cannot test on Windows")
def test_color_print3():
    # Test that this thinks the FakeTTY is a tty and applies colors.

    stream = FakeTTY()
    console.color_print("foo", "green", file=stream)
    assert stream.getvalue() == "\x1b[0;32mfoo\x1b[0m\n"

    stream = FakeTTY()
    console.color_print("foo", "green", "bar", "red", "baz", file=stream)
    assert stream.getvalue() == "\x1b[0;32mfoo\x1b[0m\x1b[0;31mbar\x1b[0mbaz\n"


def test_color_print_unicode():
    console.color_print("überbær", "red")


def test_color_print_invalid_color():
    console.color_print("foo", "unknown")


def test_spinner_non_unicode_console():
    """Regression test for #1760

    Ensures that the spinner can fall go into fallback mode when using the
    unicode spinner on a terminal whose default encoding cannot encode the
    unicode characters.
    """

    stream = FakeTTY("ascii")
    chars = console.Spinner._default_unicode_chars

    with console.Spinner("Reticulating splines", file=stream, chars=chars) as s:
        next(s)


def test_progress_bar():
    # This stuff is hard to test, at least smoke test it
    with console.ProgressBar(50) as bar:
        for _ in range(50):
            bar.update()


def test_progress_bar2():
    for _ in console.ProgressBar(range(50)):
        pass


def test_progress_bar3():
    def do_nothing(*args, **kwargs):
        pass

    console.ProgressBar.map(do_nothing, range(50))


def test_zero_progress_bar():
    with console.ProgressBar(0) as bar:
        pass


def test_progress_bar_as_generator():
    sum = 0
    for x in console.ProgressBar(range(50)):
        sum += x
    assert sum == 1225

    sum = 0
    for x in console.ProgressBar(50):
        sum += x
    assert sum == 1225


def test_progress_bar_map():
    items = list(range(100))
    result = console.ProgressBar.map(
        test_progress_bar_func.func, items, step=10, multiprocess=True
    )
    assert items == result

    result1 = console.ProgressBar.map(
        test_progress_bar_func.func, items, step=10, multiprocess=2
    )

    assert items == result1


@pytest.mark.parametrize(
    ("seconds", "string"),
    [
        (864088, " 1w 3d"),
        (187213, " 2d 4h"),
        (3905, " 1h 5m"),
        (64, " 1m 4s"),
        (15, "   15s"),
        (2, "    2s"),
    ],
)
def test_human_time(seconds, string):
    human_time = console.human_time(seconds)
    assert human_time == string


@pytest.mark.parametrize(
    ("size", "string"),
    [
        (8640882, "8.6M"),
        (187213, "187k"),
        (3905, "3.9k"),
        (64, " 64 "),
        (2, "  2 "),
        (10 * u.GB, " 10G"),
    ],
)
def test_human_file_size(size, string):
    human_time = console.human_file_size(size)
    assert human_time == string


@pytest.mark.parametrize("size", (50 * u.km, 100 * u.g))
def test_bad_human_file_size(size):
    assert pytest.raises(u.UnitConversionError, console.human_file_size, size)
