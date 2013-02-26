# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Utilities for console input and output.
"""
from __future__ import division, print_function

import re
import math
import multiprocessing
import sys
import threading
import time

try:
    get_ipython()
except NameError:
    OutStream = None
else:
    try:
        from IPython.zmq.iostream import OutStream
    except ImportError:
        OutStream = None

from ..config import ConfigurationItem


__all__ = [
    'isatty', 'color_print', 'human_time', 'ProgressBar', 'Spinner',
    'print_code_line', 'ProgressBarOrSpinner']


USE_COLOR = ConfigurationItem(
    'use_color', True,
    'When True, use ANSI color escape sequences when writing to the console.')
USE_UNICODE = ConfigurationItem(
    'use_unicode', True,
    'Use Unicode characters when drawing progress bars etc. at the console.')


def isatty(file):
    """
    Returns `True` if `file` is a tty.

    Most built-in Python file-like objects have an `isatty` member,
    but some user-defined types may not, so this assumes those are not
    ttys.
    """
    if (multiprocessing.current_process().name != 'MainProcess' or
        threading.current_thread().getName() != 'MainThread'):
        return False

    if (OutStream is not None and
        isinstance(file, OutStream) and
        file.name == 'stdout'):
        return True
    elif hasattr(file, 'isatty'):
        return file.isatty()
    return False

def _color_text(text, color):
    """
    Returns a string wrapped in ANSI color codes for coloring the
    text in a terminal::

        colored_text = color_text('Here is a message', 'blue')

    This won't actually effect the text until it is printed to the
    terminal.

    Parameters
    ----------
    text : str
        The string to return, bounded by the color codes.
    color : str
        An ANSI terminal color name. Must be one of:
        black, red, green, brown, blue, magenta, cyan, lightgrey,
        default, darkgrey, lightred, lightgreen, yellow, lightblue,
        lightmagenta, lightcyan, white, or '' (the empty string).
    """
    color_mapping = {
        'black': '0;30',
        'red': '0;31',
        'green': '0;32',
        'brown': '0;33',
        'blue': '0;34',
        'magenta': '0;35',
        'cyan': '0;36',
        'lightgrey': '0;37',
        'default': '0;39',
        'darkgrey': '1;30',
        'lightred': '1;31',
        'lightgreen': '1;32',
        'yellow': '1;33',
        'lightblue': '1;34',
        'lightmagenta': '1;35',
        'lightcyan': '1;36',
        'white': '1;37'}

    color_code = color_mapping.get(color, '0;39')
    return u'\033[{0}m{1}\033[0m'.format(color_code, text)

def color_print(*args, **kwargs):
    """
    Prints colors and styles to the terminal uses ANSI escape
    sequences.

    ::

       color_print('This is the color ', 'default', 'GREEN', 'green')

    Parameters
    ----------
    positional args : strings
        The positional arguments come in pairs (*msg*, *color*), where
        *msg* is the string to display and *color* is the color to
        display it in.

        *color* is an ANSI terminal color name.  Must be one of:
        black, red, green, brown, blue, magenta, cyan, lightgrey,
        default, darkgrey, lightred, lightgreen, yellow, lightblue,
        lightmagenta, lightcyan, white, or '' (the empty string).

    file : writeable file-like object, optional
        Where to write to.  Defaults to `sys.stdout`.  If file is not
        a tty (as determined by calling its `isatty` member, if one
        exists), no coloring will be included.

    end : str, optional
        The ending of the message.  Defaults to ``\\n``.  The end will
        be printed after resetting any color or font state.
    """

    file = kwargs.get('file', sys.stdout)
    end = kwargs.get('end', u'\n')

    write = file.write
    if isatty(file) and USE_COLOR():
        for i in xrange(0, len(args), 2):
            msg = args[i]
            if i + 1 == len(args):
                color = ''
            else:
                color = args[i + 1]

            if isinstance(msg, bytes):
                msg = msg.decode('ascii')

            if color == u'' or color is None:
                write(msg)
            else:
                write(_color_text(msg, color))

        write(end)
    else:
        for i in xrange(0, len(args), 2):
            msg = args[i]
            if isinstance(msg, bytes):
                msg = msg.decode('ascii')
            write(msg)
        write(end)

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
        (u'y', 60 * 60 * 24 * 7 * 52),
        (u'w', 60 * 60 * 24 * 7),
        (u'd', 60 * 60 * 24),
        (u'h', 60 * 60),
        (u'm', 60),
        (u's', 1),
        ]

    seconds = int(seconds)

    if seconds < 60:
        return u'   {0:02d}s'.format(seconds)
    for i in xrange(len(units) - 1):
        unit1, limit1 = units[i]
        unit2, limit2 = units[i + 1]
        if seconds >= limit1:
            return u'{0:02d}{1}{2:02d}{3}'.format(
                seconds // limit1, unit1,
                (seconds % limit1) // limit2, unit2)
    return u'  ~inf'


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
        total : int
            The number of increments in the process being tracked.

        file : writable file-like object, optional
            The file to write the progress bar to.  Defaults to
            `sys.stdout`.  If `file` is not a tty (as determined by
            calling its `isatty` member, if any), the scrollbar will
            be completely silent.
        """
        if not isatty(file):
            self.update = self._silent_update
            self._silent = True
        else:
            self._silent = False

        self._total = total
        self._file = file
        self._start_time = time.time()
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
        self._bar_length = terminal_width - 36
        if self._total == 0:
            num_scale = 0
        else:
            num_scale = int(math.floor(math.log(self._total) / math.log(1000)))
        if num_scale > 7:
            self._suffix = '?'
        else:
            suffixes = u' kMGTPEH'
            self._suffix = suffixes[num_scale]
        self._num_scale = int(math.pow(1000, num_scale))
        self.update(0)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        if not self._silent:
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

        file = self._file
        write = file.write

        bar_fill = int(float(self._bar_length) * frac)
        write(u'\r|')
        color_print(u'=' * bar_fill, 'blue', file=file, end=u'')
        if bar_fill < self._bar_length:
            color_print(u'>', 'green', file=file, end=u'')
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
        write(u' {0:>3d}/{1:>3d}{2}'.format(
            value // self._num_scale,
            self._total // self._num_scale,
            self._suffix))
        write(u' ({0:>6s}%)'.format(u'{0:.2f}'.format(frac * 100.0)))
        write(prefix)
        if t is not None:
            write(human_time(t))
        self._file.flush()

    def _silent_update(self, value=None):
        pass

    @classmethod
    def map(cls, function, items, multiprocess=False, file=sys.stdout):
        """
        Does a `map` operation while displaying a progress bar with
        percentage complete.

        ::

            def work(i):
                print(i)

            ProgressBar.map(work, range(50))

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

        file : writeable file-like object, optional
            The file to write the progress bar to.  Defaults to
            `sys.stdout`.  If `file` is not a tty (as determined by
            calling its `isatty` member, if any), the scrollbar will
            be completely silent.
        """
        results = []

        with cls(len(items), file=file) as bar:
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
                for i, result in enumerate(
                    p.imap_unordered(function, items, steps)):
                    bar.update(i)
                    results.append(result)

        return results

    @classmethod
    def iterate(cls, items, file=sys.stdout):
        """
        Iterate over a sequence while indicating progress with a progress
        bar in the terminal.

        ::

            for item in ProgressBar.iterate(items):
                pass

        Parameters
        ----------
        items : sequence
            A sequence of items to iterate over

        file : writeable file-like object, optional
            The file to write the progress bar to.  Defaults to
            `sys.stdout`.  If `file` is not a tty (as determined by
            calling its `isatty` member, if any), the scrollbar will
            be completely silent.

        Returns
        -------
        generator :
            A generator over `items`
        """
        with cls(len(items), file=file) as bar:
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

    def __init__(self, msg, color='default', file=sys.stdout, step=1,
                 chars=None):
        """
        Parameters
        ----------
        msg : str
            The message to print

        color : str, optional
            An ANSI terminal color name.  Must be one of: black, red,
            green, brown, blue, magenta, cyan, lightgrey, default,
            darkgrey, lightred, lightgreen, yellow, lightblue,
            lightmagenta, lightcyan, white.

        file : writeable file-like object, optional
            The file to write the spinner to.  Defaults to
            `sys.stdout`.  If `file` is not a tty (as determined by
            calling its `isatty` member, if any), the scrollbar will
            be completely silent.

        step : int, optional
            Only update the spinner every *step* steps

        chars : str, optional
            The character sequence to use for the spinner
        """
        self._msg = msg
        self._color = color
        self._file = file
        self._step = step
        if chars is None:
            if USE_UNICODE:
                chars = self._default_unicode_chars
            else:
                chars = self._default_ascii_chars
        self._chars = chars

        self._silent = not isatty(file)

    def _iterator(self):
        chars = self._chars
        index = 0
        file = self._file
        write = file.write
        flush = file.flush

        while True:
            write(u'\r')
            color_print(self._msg, self._color, file=file, end=u'')
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
        if self._silent:
            return self._silent_iterator()
        else:
            return self._iterator()

    def __exit__(self, exc_type, exc_value, traceback):
        file = self._file
        write = file.write
        flush = file.flush

        if not self._silent:
            write(u'\r')
            color_print(self._msg, self._color, file=file, end=u'')
        if exc_type is None:
            color_print(u' [Done]', 'green', file=file)
        else:
            color_print(u' [Failed]', 'red', file=file)
        flush()

    def _silent_iterator(self):
        color_print(self._msg, self._color, file=self._file, end=u'')
        self._file.flush()

        while True:
            yield


class ProgressBarOrSpinner:
    """
    A class that displays either a `ProgressBar` or `Spinner`
    depending on whether the total size of the operation is
    known or not.

    It is designed to be used with the `with` statement::

        if file.has_length():
            length = file.get_length()
        else:
            length = None
        bytes_read = 0
        with ProgressBarOrSpinner(length) as bar:
            while file.read(blocksize):
                bytes_read += blocksize
                bar.update(bytes_read)
    """

    def __init__(self, total, msg, color='default', file=sys.stdout):
        """
        Parameters
        ----------
        total : int or None
            If an int, the number of increments in the process being
            tracked and a `ProgressBar` is displayed.  If `None`, a
            `Spinner` is displayed.

        msg : str
            The message to display above the `ProgressBar` or
            alongside the `Spinner`.

        color : str, optional
            The color of `msg`, if any.  Must be an ANSI terminal
            color name.  Must be one of: black, red, green, brown,
            blue, magenta, cyan, lightgrey, default, darkgrey,
            lightred, lightgreen, yellow, lightblue, lightmagenta,
            lightcyan, white.

        file : writable file-like object, optional
            The file to write the to.  Defaults to `sys.stdout`.  If
            `file` is not a tty (as determined by calling its `isatty`
            member, if any), only `msg` will be displayed: the
            `ProgressBar` or `Spinner` will be silent.
        """
        if total is None or not isatty(file):
            self._is_spinner = True
            self._obj = Spinner(msg, color=color, file=file)
        else:
            self._is_spinner = False
            color_print(msg, color, file=file)
            self._obj = ProgressBar(total, file=file)

    def __enter__(self):
        self._iter = self._obj.__enter__()
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        return self._obj.__exit__(exc_type, exc_value, traceback)

    def update(self, value):
        """
        Update the progress bar to the given value (out of the total
        given to the constructor.
        """
        if self._is_spinner:
            self._iter.next()
        else:
            self._obj.update(value)


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
        The number of spaces per tab (``'\\t'``) character.  Default
        is 8.  All tabs will be converted to spaces to ensure that the
        caret lines up with the correct column.

    width : int, optional
        The width of the display, beyond which the line will be
        truncated.  Defaults to 70 (this matches the default in the
        standard library's `textwrap` module).
    """
    write = file.write

    if col is not None:
        assert col < len(line)
        ntabs = line[:col].count(u'\t')
        col += ntabs * (tabwidth - 1)

    line = line.rstrip(u'\n')
    line = line.replace(u'\t', u' ' * tabwidth)

    if col is not None and col > width:
        new_col = min(width // 2, len(line) - col)
        offset = col - new_col
        line = line[offset + 1:]
        new_col = col
        col -= offset
        width = width - 3
        color_print(u'…', 'darkgrey', file=file, end=u'')

    if len(line) > width:
        write(line[:width - 1])
        color_print(u'…', 'darkgrey', file=file)
    else:
        write(line)
        write(u'\n')

    if col is not None:
        write(u' ' * col)
        color_print(u'^', 'red', file=file)


# The following four Getch* classes implement unbuffered character reading from
# stdin on Windows, linux, MacOSX.  This is taken directly from ActiveState
# Code Recipes:
# http://code.activestate.com/recipes/134892-getch-like-unbuffered-character-reading-from-stdin/
#

class Getch(object):
    """Get a single character from standard input without screen echo.

    Returns
    -------
    char : str (one character)
    """
    def __init__(self):
        try:
            self.impl = _GetchWindows()
        except ImportError:
            try:
                self.impl = _GetchMacCarbon()
            except (ImportError, AttributeError):
                self.impl = _GetchUnix()

    def __call__(self):
        return self.impl()


class _GetchUnix(object):
    def __init__(self):
        import tty
        import sys
        import termios  # import termios now or else you'll get the Unix
                        # version on the Mac

    def __call__(self):
        import sys
        import tty
        import termios
        fd = sys.stdin.fileno()
        old_settings = termios.tcgetattr(fd)
        try:
            tty.setraw(sys.stdin.fileno())
            ch = sys.stdin.read(1)
        finally:
            termios.tcsetattr(fd, termios.TCSADRAIN, old_settings)
        return ch


class _GetchWindows(object):
    def __init__(self):
        import msvcrt

    def __call__(self):
        import msvcrt
        return msvcrt.getch()


class _GetchMacCarbon:
    """
    A function which returns the current ASCII key that is down;
    if no ASCII key is down, the null string is returned.  The
    page http://www.mactech.com/macintosh-c/chap02-1.html was
    very helpful in figuring out how to do this.
    """
    def __init__(self):
        import Carbon
        Carbon.Evt  # see if it has this (in Unix, it doesn't)

    def __call__(self):
        import Carbon
        if Carbon.Evt.EventAvail(0x0008)[0] == 0:  # 0x0008 is the keyDownMask
            return ''
        else:
            #
            # The event contains the following info:
            # (what,msg,when,where,mod)=Carbon.Evt.GetNextEvent(0x0008)[1]
            #
            # The message (msg) contains the ASCII char which is
            # extracted with the 0x000000FF charCodeMask; this
            # number is converted to an ASCII character with chr() and
            # returned
            #
            (what, msg, when, where, mod) = Carbon.Evt.GetNextEvent(0x0008)[1]
            return chr(msg & 0x000000FF)
