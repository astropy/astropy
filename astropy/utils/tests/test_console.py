# -*- coding: utf-8 -*-

import io
import sys

from astropy.tests.helper import raises

from .. import console


def test_color_print():
    # This stuff is hard to test, at least smoke test it
    console.color_print("foo", "green")


def test_color_print2():
    # Test that this automatically detects that io.StringIO is
    # not a tty
    stream = io.StringIO()
    console.color_print("foo", "green", file=stream)
    assert stream.getvalue() == u'foo\n'


def test_color_print3():
    # Test that this things the FakeTTY is a tty and applies colors.
    class FakeTTY(io.StringIO):
        def isatty(self):
            return True

    stream = FakeTTY()
    console.color_print("foo", "green", file=stream)
    assert stream.getvalue() == u'\x1b[0;32mfoo\x1b[0m\n'


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
    for x in console.iterate_with_progress_bar(range(50)):
        pass


def test_progress_bar3():
    def do_nothing(*args, **kwargs):
        pass

    console.map_with_progress_bar(do_nothing, range(50))
