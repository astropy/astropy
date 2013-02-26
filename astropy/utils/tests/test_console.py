# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-

import io
import sys

from ...tests.helper import raises

from .. import console

def test_color_text():
    assert console._color_text("foo", "green") == u'\033[0;32mfoo\033[0m'

def test_color_print():
    # This stuff is hard to test, at least smoke test it
    console.color_print("foo", "green")

    console.color_print("foo", "green", "bar", "red")


def test_color_print2():
    # Test that this automatically detects that io.StringIO is
    # not a tty
    stream = io.StringIO()
    console.color_print("foo", "green", file=stream)
    assert stream.getvalue() == u'foo\n'

    stream = io.StringIO()
    console.color_print("foo", "green", "bar", "red", "baz", file=stream)
    assert stream.getvalue() == u'foobarbaz\n'


def test_color_print3():
    # Test that this things the FakeTTY is a tty and applies colors.
    class FakeTTY(io.StringIO):
        def isatty(self):
            return True

    stream = FakeTTY()
    console.color_print("foo", "green", file=stream)
    assert stream.getvalue() == u'\x1b[0;32mfoo\x1b[0m\n'

    stream = FakeTTY()
    console.color_print("foo", "green", "bar", "red", "baz", file=stream)
    assert stream.getvalue() == u'\x1b[0;32mfoo\x1b[0m\x1b[0;31mbar\x1b[0mbaz\n'


def test_color_print_unicode():
    console.color_print(u"überbær", "red")


def test_color_print_invalid_color():
    console.color_print("foo", "unknown")


def test_progress_bar():
    # This stuff is hard to test, at least smoke test it
    with console.ProgressBar(50) as bar:
        for i in range(50):
            bar.update()


def test_progress_bar2():
    for x in console.ProgressBar.iterate(range(50)):
        pass


def test_progress_bar3():
    def do_nothing(*args, **kwargs):
        pass

    console.ProgressBar.map(do_nothing, range(50))


def test_zero_progress_bar():
    with console.ProgressBar(0) as bar:
        pass
