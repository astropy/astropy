# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Utilities for prettifying output to the console.
"""
from __future__ import print_function

import io
import re
import sys
import time


from ..config.configs import ConfigurationItem


__all__ = [
    'Color', 'color_print', 'color_string', 'human_time',
    'ProgressBar', 'iterate_with_progress_bar',
    'map_with_progress_bar', 'Spinner', 'print_code_line']


USE_COLOR = ConfigurationItem(
    'use_color', True,
    'When True, use ANSI color escape sequences when writing to the console. '
    'Changing has no effect after startup time')
USE_UNICODE = ConfigurationItem(
    'use_unicode', True,
    'Use Unicode characters when drawing progress bars etc. at the console.')


if USE_COLOR:
    def _define_color(name, code):
        template = u'\033[{0}{{0}}m{{1}}\033[0m'.format(code)

        def color(s, bold=False):
            if bold:
                styles = u';1'
            else:
                styles = u''
            if isinstance(s, bytes):
                s = s.decode('ascii')
            return template.format(styles, s)

        func = color
        func.__name__ = name
        return func
else:
    def _define_color(name, code):
        def color(s, bold=False):
            if isinstance(s, bytes):
                s = s.decode('ascii')
            return s

        func = color
        func.__name__ = name
        return func


class Color:
    """
    A class for colorizing text.

    It contains a number of static methods, each of which colorizes a
    particular color.  For example, to print out a string in the color
    green, do::

        print(Color.green('go'))
    """
    black        = staticmethod(_define_color('black', '0;30'))
    red          = staticmethod(_define_color('red', '0;31'))
    green        = staticmethod(_define_color('green', '0;32'))
    brown        = staticmethod(_define_color('brown', '0:33'))
    blue         = staticmethod(_define_color('blue', '0;34'))
    magenta      = staticmethod(_define_color('magenta', '0;35'))
    cyan         = staticmethod(_define_color('cyan', '0;36'))
    lightgrey    = staticmethod(_define_color('lightgrey', '0;37'))
    default      = staticmethod(_define_color('default', '0;39'))
    darkgrey     = staticmethod(_define_color('darkgrey', '1;30'))
    lightred     = staticmethod(_define_color('lightred', '1;31'))
    lightgreen   = staticmethod(_define_color('lightgreen', '1;32'))
    yellow       = staticmethod(_define_color('yellow', '1;33'))
    lightblue    = staticmethod(_define_color('lightblue', '1;34'))
    lightmagenta = staticmethod(_define_color('lightmagenta', '1;35'))
    lightcyan    = staticmethod(_define_color('lightcyan', '1;36'))
    white        = staticmethod(_define_color('white', '1;37'))

    @classmethod
    def get_color_func(cls, color):
        """
        Return the function for a particular `color` by name.
        """
        try:
            func = getattr(cls, color)
        except AttributeError:
            func = cls.default
        return func


def color_print(s, color='default', bold=False, file=sys.stdout, end=u'\n'):
    """
    Prints colors and styles to the terminal uses ANSI escape
    sequences.

    Parameters
    ----------
    s : str
        The message to print

    color : str, optional
        An ANSI terminal color name.  Must be one of: black, red,
        green, brown, blue, magenta, cyan, lightgrey, default,
        darkgrey, lightred, lightgreen, yellow, lightblue,
        lightmagenta, lightcyan, white.

    bold : bool, optional
        When `True` use boldface font.

    file : writeable file-like object, optional
        Where to write to.  Defaults to `sys.stdout`.

    end : str, optional
        The ending of the message.  Defaults to ``\\n``.  The end will
        be printed after resetting any color or font state.
    """
    color_func = Color.get_color_func(color)

    print(color_func(s, bold), file=file, end=end)


def color_string(s, color='default', bold=False):
    """
    Returns a string containing ANSI color codes in the given
    color.

    Parameters
    ----------
    s : str
        The message to color

    color : str, optional
        An ANSI terminal color name.  Must be one of: black, red,
        green, brown, blue, magenta, cyan, lightgrey, default,
        darkgrey, lightred, lightgreen, yellow, lightblue,
        lightmagenta, lightcyan, white.

    bold : bool, optional
        When `True` use boldface font.
    """
    color_func = Color.get_color_func(color)

    return color_func(s, bold)


def strip_ansi_codes(s):
    """
    Remove ANSI color codes from the string.
    """
    return re.sub('\033\[([0-9]+)(;[0-9]+)*m', '', s)


def human_time(seconds):
    """
    Returns a human-friendly time string that is always exactly 6
    characters long.

    Depending on the number of seconds given, can be one of::

        1w 3d
        2d 4h
        1h 5m
        1m 4s
          15s

    Will be in color if console coloring is turned on.

    Parameters
    ----------
    seconds : int
        The number of seconds to represent

    Returns
    -------
    time : str
        A human-friendly representation of the given number of seconds
        that is always exactly 6 characters.
    """
    units = [
        ('y', 60 * 60 * 24 * 7 * 52),
        ('w', 60 * 60 * 24 * 7),
        ('d', 60 * 60 * 24),
        ('h', 60 * 60),
        ('m', 60),
        ('s', 1),
        ]

    seconds = int(seconds)

    if seconds < 60:
        return '   {0:02d}s'.format(seconds)
    for i in xrange(len(units) - 1):
        unit1, limit1 = units[i]
        unit2, limit2 = units[i + 1]
        if seconds > limit1:
            return '{0:02d}{1}{2:02d}{3}'.format(
                seconds // limit1, unit1,
                (seconds % limit1) // limit2, unit2)


class ProgressBar:
    """
    A class to display a progress bar in the terminal.

    It is designed to be used with the `with` statement::

        with ProgressBar(len(items)) as bar:
            for item in enumerate(items):
                bar.update()
    """
    def __init__(self, total, file=sys.stdout):
        """
        Parameters
        ----------
        items : int
            The number of increments in the process being tracked.

        file : writable file-like object
            The file to write the progress bar to.  Defaults to
            `sys.stdout`.

        See also
        --------
        map_with_progress_bar

        iterate_with_progress_bar
        """
        self._total = total
        self._file = file
        self._start_time = time.time()
        num_length = len(str(total))
        terminal_width = 78
        if sys.platform.startswith('linux'):
            import subprocess
            p = subprocess.Popen(
                'stty size',
                shell=True,
                stdout=subprocess.PIPE)
            stdout, stderr = p.communicate()
            parts = stdout.split()
            if len(parts) == 2:
                rows, cols = parts
                terminal_width = int(cols)
        self._bar_length = terminal_width - 29 - (num_length * 2)
        self._num_format = '{{0:>{0}}}/{{1:>{0}}}'.format(num_length)
        self.update(0)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        if exc_type is None:
            self.update(self._total)
        self._file.write('\n')
        self._file.flush()

    def update(self, value=None):
        """
        Update the progress bar to the given value (out of the total
        given to the constructor.
        """
        if value is None:
            value = self._current_value = self._current_value + 1
        else:
            self._current_value = value
        if self._total == 0:
            frac = 1.0
        else:
            frac = float(value) / float(self._total)

        write = self._file.write

        bar_fill = int(float(self._bar_length) * frac)
        write(u'\r|')
        write(Color.blue(u'=' * bar_fill))
        if bar_fill < self._bar_length:
            write(Color.green(u'>'))
            write(u'-' * (self._bar_length - bar_fill - 1))
        write(u'|')

        if value >= self._total:
            t = time.time() - self._start_time
            prefix = u'     '
        elif value <= 0:
            t = None
            prefix = u''
        else:
            t = ((time.time() - self._start_time) * (1.0 - frac)) / frac
            prefix = u' ETA '
        write(u' ')
        write(self._num_format.format(value, self._total))
        write(u' ({0:>6s}%)'.format(u'{0:.2f}'.format(frac * 100.0)))
        write(prefix)
        if t is not None:
            write(human_time(t))
        self._file.flush()


def map_with_progress_bar(
        function, items, multiprocess=False, file=sys.stdout):
    """
    Does a `map` operation while displaying a progress bar with
    percentage complete.

    Parameters
    ----------
    function : function
        Function to call for each step

    items : sequence
        Sequence where each element is a tuple of arguments to pass to
        *function*.

    multiprocess : bool, optional
        If `True`, use the `multiprocessing` module to distribute each
        task to a different processor core.

    file : writeable file-like object
        The file to write the progress bar to.  Defaults to
        `sys.stdout`.
    """
    with ProgressBar(len(items), file=file) as bar:
        step_size = max(200, bar._bar_length)
        steps = max(int(float(len(items)) / step_size), 1)
        if not multiprocess:
            for i, item in enumerate(items):
                function(item)
                if (i % steps) == 0:
                    bar.update(i)
        else:
            import multiprocessing
            p = multiprocessing.Pool()
            for i, _ in enumerate(p.imap_unordered(function, items, steps)):
                bar.update(i)


def iterate_with_progress_bar(items, file=sys.stdout):
    """
    Iterate over a sequence while indicating progress with a progress
    bar in the terminal.

    Use as follows::

        for item in iterate_with_progress_bar(items):
            pass

    Parameters
    ----------
    items : sequence
        A sequence of items to iterate over

    file : writeable file-like object
        The file to write the progress bar to.  Defaults to
        `sys.stdout`.

    Returns
    -------
    generator :
        A generator over `items`
    """
    with ProgressBar(len(items), file=file) as bar:
        for item in items:
            yield item
            bar.update()


class Spinner():
    """
    A class to display a spinner in the terminal.

    It is designed to be used with the `with` statement::

        with Spinner("Reticulating splines", "green") as s:
            for item in enumerate(items):
                s.next()
    """
    _default_unicode_chars = u"◓◑◒◐"
    _default_ascii_chars = u"-/|\\"

    def __init__(self, s, file=sys.stdout, step=1, chars=None):
        """
        Parameters
        ----------
        s : str
            The message to print

        file : writeable file-like object, optional
            Where to write to.  Defaults to `sys.stdout`.

        step : int, optional
            Only update the spinner every *step* steps

        chars : str, optional
            The character sequence to use for the spinner
        """
        self._s = s
        self._file = file
        self._step = step
        if chars is None:
            if USE_UNICODE:
                chars = self._default_unicode_chars
            else:
                chars = self._default_ascii_chars
        self._chars = chars

    def _iterator(self):
        chars = self._chars
        index = 0
        write = self._file.write
        flush = self._file.flush

        while True:
            write(u'\r')
            write(self._s)
            write(u' ')
            write(chars[index])
            flush()
            yield

            for i in xrange(self._step):
                yield

            index += 1
            if index == len(chars):
                index = 0

    def __enter__(self):
        return self._iterator()

    def __exit__(self, exc_type, exc_value, traceback):
        write = self._file.write

        write(u'\r')
        write(self._s)
        if exc_type is None:
            write(Color.green(u' Done', True))
        else:
            write(Color.red(u' Failed', True))
        write(u'\n')
        self._file.flush()


def print_code_line(line, col=None, file=sys.stdout, tabwidth=8, width=70):
    u"""
    Prints a line of source code, highlighting a particular character
    position in the line.  Useful for displaying the context of error
    messages.

    If the line is more than `width` characters, the line is truncated
    accordingly and '…' characters are inserted at the front and/or
    end.

    It looks like this::

        there_is_a_syntax_error_here :
                                     ^

    Parameters
    ----------
    line : unicode
        The line of code to display

    col : int, optional
        The character in the line to highlight.  `col` must be less
        than `len(line)`.

    file : writeable file-like object, optional
        Where to write to.  Defaults to `sys.stdout`.

    tabwidth : int, optional
        The number of spaces per tab (``'\t'``) character.  Default is
        8.  All tabs will be converted to spaces to ensure that the
        caret lines up with the correct column.

    width : int, optional
        The width of the display, beyond which the line will be
        truncated.  Defaults to 70 (this matches the default in the
        standard library's `textwrap` module).
    """
    if col is not None:
        assert col < len(line)
        ntabs = line[:col].count(u'\t')
        col += ntabs * (tabwidth - 1)

    line = line.rstrip(u'\n')
    line = line.replace(u'\t', u' ' * tabwidth)

    if col is not None and col > width:
        new_col = min(width / 2, len(line) - col)
        offset = col - new_col
        line = line[offset + 1: ]
        new_col = col
        width = width - 3
        file.write(color.darkgrey(u'…'))

    if len(line) > width:
        file.write(line[:width-1])
        file.write(Color.darkgrey(u'…'))
        file.write(u'\n')
    else:
        file.write(line)
        file.write(u'\n')

    if col is not None:
        file.write(u' ' * col)
        file.write(Color.red(u'^'))
        file.write(u'\n')

