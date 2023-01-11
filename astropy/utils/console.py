# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Utilities for console input and output.
"""

import codecs
import locale
import math
import multiprocessing
import os
import re
import struct
import sys
import threading
import time

# concurrent.futures imports moved inside functions using them to avoid
# import failure when running in pyodide/Emscripten

try:
    import fcntl
    import signal
    import termios

    _CAN_RESIZE_TERMINAL = True
except ImportError:
    _CAN_RESIZE_TERMINAL = False

from astropy import conf

from .decorators import classproperty
from .misc import isiterable

__all__ = [
    "isatty",
    "color_print",
    "human_time",
    "human_file_size",
    "ProgressBar",
    "Spinner",
    "print_code_line",
    "ProgressBarOrSpinner",
    "terminal_size",
]

_DEFAULT_ENCODING = "utf-8"


class _IPython:
    """Singleton class given access to IPython streams, etc."""

    @classproperty
    def get_ipython(cls):
        try:
            from IPython import get_ipython
        except ImportError:
            pass
        return get_ipython

    @classproperty
    def OutStream(cls):
        if not hasattr(cls, "_OutStream"):
            cls._OutStream = None
            try:
                cls.get_ipython()
            except NameError:
                return None

            try:
                from ipykernel.iostream import OutStream
            except ImportError:
                try:
                    from IPython.zmq.iostream import OutStream
                except ImportError:
                    from IPython import version_info

                    if version_info[0] >= 4:
                        return None

                    try:
                        from IPython.kernel.zmq.iostream import OutStream
                    except ImportError:
                        return None

            cls._OutStream = OutStream

        return cls._OutStream

    @classproperty
    def ipyio(cls):
        if not hasattr(cls, "_ipyio"):
            try:
                from IPython.utils import io
            except ImportError:
                cls._ipyio = None
            else:
                cls._ipyio = io
        return cls._ipyio

    @classmethod
    def get_stream(cls, stream):
        return getattr(cls.ipyio, stream)


def _get_stdout(stderr=False):
    """
    This utility function contains the logic to determine what streams to use
    by default for standard out/err.

    Typically this will just return `sys.stdout`, but it contains additional
    logic for use in IPython on Windows to determine the correct stream to use
    (usually ``IPython.util.io.stdout`` but only if sys.stdout is a TTY).
    """
    if stderr:
        stream = "stderr"
    else:
        stream = "stdout"

    sys_stream = getattr(sys, stream)
    return sys_stream


def isatty(file):
    """
    Returns `True` if ``file`` is a tty.

    Most built-in Python file-like objects have an `isatty` member,
    but some user-defined types may not, so this assumes those are not
    ttys.
    """
    if (
        multiprocessing.current_process().name != "MainProcess"
        or threading.current_thread().name != "MainThread"
    ):
        return False

    if hasattr(file, "isatty"):
        return file.isatty()

    if _IPython.OutStream is None or (not isinstance(file, _IPython.OutStream)):
        return False

    # File is an IPython OutStream. Check whether:
    # - File name is 'stdout'; or
    # - File wraps a Console
    if getattr(file, "name", None) == "stdout":
        return True

    if hasattr(file, "stream"):
        # FIXME: pyreadline has no had new release since 2015, drop it when
        #        IPython minversion is 5.x.
        # On Windows, in IPython 2 the standard I/O streams will wrap
        # pyreadline.Console objects if pyreadline is available; this should
        # be considered a TTY.
        try:
            from pyreadline.console import Console as PyreadlineConsole
        except ImportError:
            return False

        return isinstance(file.stream, PyreadlineConsole)

    return False


def terminal_size(file=None):
    """
    Returns a tuple (height, width) containing the height and width of
    the terminal.

    This function will look for the width in height in multiple areas
    before falling back on the width and height in astropy's
    configuration.
    """
    if file is None:
        file = _get_stdout()

    try:
        s = struct.pack("HHHH", 0, 0, 0, 0)
        x = fcntl.ioctl(file, termios.TIOCGWINSZ, s)
        (lines, width, xpixels, ypixels) = struct.unpack("HHHH", x)
        if lines > 12:
            lines -= 6
        if width > 10:
            width -= 1
        if lines <= 0 or width <= 0:
            raise Exception("unable to get terminal size")
        return (lines, width)
    except Exception:
        try:
            # see if POSIX standard variables will work
            return (int(os.environ.get("LINES")), int(os.environ.get("COLUMNS")))
        except TypeError:
            # fall back on configuration variables, or if not
            # set, (25, 80)
            lines = conf.max_lines
            width = conf.max_width
            if lines is None:
                lines = 25
            if width is None:
                width = 80
            return lines, width


def _color_text(text, color):
    """Returns a string wrapped in ANSI color codes for coloring the text in a terminal.

    ::

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
        "black": "0;30",
        "red": "0;31",
        "green": "0;32",
        "brown": "0;33",
        "blue": "0;34",
        "magenta": "0;35",
        "cyan": "0;36",
        "lightgrey": "0;37",
        "default": "0;39",
        "darkgrey": "1;30",
        "lightred": "1;31",
        "lightgreen": "1;32",
        "yellow": "1;33",
        "lightblue": "1;34",
        "lightmagenta": "1;35",
        "lightcyan": "1;36",
        "white": "1;37",
    }

    if sys.platform == "win32" and _IPython.OutStream is None:
        # On Windows do not colorize text unless in IPython
        return text

    color_code = color_mapping.get(color, "0;39")
    return f"\033[{color_code}m{text}\033[0m"


def _decode_preferred_encoding(s):
    """Decode the supplied byte string using the preferred encoding
    for the locale (`locale.getpreferredencoding`) or, if the default encoding
    is invalid, fall back first on utf-8, then on latin-1 if the message cannot
    be decoded with utf-8.
    """
    enc = locale.getpreferredencoding()
    try:
        try:
            return s.decode(enc)
        except LookupError:
            enc = _DEFAULT_ENCODING
        return s.decode(enc)
    except UnicodeDecodeError:
        return s.decode("latin-1")


def _write_with_fallback(s, write, fileobj):
    """Write the supplied string with the given write function like
    ``write(s)``, but use a writer for the locale's preferred encoding in case
    of a UnicodeEncodeError.  Failing that attempt to write with 'utf-8' or
    'latin-1'.
    """
    try:
        write(s)
        return write
    except UnicodeEncodeError:
        # Let's try the next approach...
        pass

    enc = locale.getpreferredencoding()
    try:
        Writer = codecs.getwriter(enc)
    except LookupError:
        Writer = codecs.getwriter(_DEFAULT_ENCODING)

    f = Writer(fileobj)
    write = f.write

    try:
        write(s)
        return write
    except UnicodeEncodeError:
        Writer = codecs.getwriter("latin-1")
        f = Writer(fileobj)
        write = f.write

    # If this doesn't work let the exception bubble up; I'm out of ideas
    write(s)
    return write


def color_print(*args, end="\n", **kwargs):
    """
    Prints colors and styles to the terminal uses ANSI escape
    sequences.

    ::

       color_print('This is the color ', 'default', 'GREEN', 'green')

    Parameters
    ----------
    positional args : str
        The positional arguments come in pairs (*msg*, *color*), where
        *msg* is the string to display and *color* is the color to
        display it in.

        *color* is an ANSI terminal color name.  Must be one of:
        black, red, green, brown, blue, magenta, cyan, lightgrey,
        default, darkgrey, lightred, lightgreen, yellow, lightblue,
        lightmagenta, lightcyan, white, or '' (the empty string).

    file : writable file-like, optional
        Where to write to.  Defaults to `sys.stdout`.  If file is not
        a tty (as determined by calling its `isatty` member, if one
        exists), no coloring will be included.

    end : str, optional
        The ending of the message.  Defaults to ``\\n``.  The end will
        be printed after resetting any color or font state.
    """
    file = kwargs.get("file", _get_stdout())

    write = file.write
    if isatty(file) and conf.use_color:
        for i in range(0, len(args), 2):
            msg = args[i]
            if i + 1 == len(args):
                color = ""
            else:
                color = args[i + 1]

            if color:
                msg = _color_text(msg, color)

            # Some file objects support writing unicode sensibly on some Python
            # versions; if this fails try creating a writer using the locale's
            # preferred encoding. If that fails too give up.

            write = _write_with_fallback(msg, write, file)

        write(end)
    else:
        for i in range(0, len(args), 2):
            msg = args[i]
            write(msg)
        write(end)


def strip_ansi_codes(s):
    """
    Remove ANSI color codes from the string.
    """
    return re.sub("\033\\[([0-9]+)(;[0-9]+)*m", "", s)


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
        ("y", 60 * 60 * 24 * 7 * 52),
        ("w", 60 * 60 * 24 * 7),
        ("d", 60 * 60 * 24),
        ("h", 60 * 60),
        ("m", 60),
        ("s", 1),
    ]

    seconds = int(seconds)

    if seconds < 60:
        return f"   {seconds:2d}s"
    for i in range(len(units) - 1):
        unit1, limit1 = units[i]
        unit2, limit2 = units[i + 1]
        if seconds >= limit1:
            return "{:2d}{}{:2d}{}".format(
                seconds // limit1, unit1, (seconds % limit1) // limit2, unit2
            )
    return "  ~inf"


def human_file_size(size):
    """
    Returns a human-friendly string representing a file size
    that is 2-4 characters long.

    For example, depending on the number of bytes given, can be one
    of::

        256b
        64k
        1.1G

    Parameters
    ----------
    size : int
        The size of the file (in bytes)

    Returns
    -------
    size : str
        A human-friendly representation of the size of the file
    """
    if hasattr(size, "unit"):
        # Import units only if necessary because the import takes a
        # significant time [#4649]
        from astropy import units as u

        size = u.Quantity(size, u.byte).value

    suffixes = " kMGTPEZY"
    if size == 0:
        num_scale = 0
    else:
        num_scale = int(math.floor(math.log(size) / math.log(1000)))
    if num_scale > 7:
        suffix = "?"
    else:
        suffix = suffixes[num_scale]
    num_scale = int(math.pow(1000, num_scale))
    value = size / num_scale
    str_value = str(value)
    if suffix == " ":
        str_value = str_value[: str_value.index(".")]
    elif str_value[2] == ".":
        str_value = str_value[:2]
    else:
        str_value = str_value[:3]
    return f"{str_value:>3s}{suffix}"


class _mapfunc:
    """
    A function wrapper to support ProgressBar.map().
    """

    def __init__(self, func):
        self._func = func

    def __call__(self, i_arg):
        i, arg = i_arg
        return i, self._func(arg)


class ProgressBar:
    """
    A class to display a progress bar in the terminal.

    It is designed to be used either with the ``with`` statement::

        with ProgressBar(len(items)) as bar:
            for item in enumerate(items):
                bar.update()

    or as a generator::

        for item in ProgressBar(items):
            item.process()
    """

    def __init__(self, total_or_items, ipython_widget=False, file=None):
        """
        Parameters
        ----------
        total_or_items : int or sequence
            If an int, the number of increments in the process being
            tracked.  If a sequence, the items to iterate over.

        ipython_widget : bool, optional
            If `True`, the progress bar will display as an IPython
            notebook widget.

        file : writable file-like, optional
            The file to write the progress bar to.  Defaults to
            `sys.stdout`.  If ``file`` is not a tty (as determined by
            calling its `isatty` member, if any, or special case hacks
            to detect the IPython console), the progress bar will be
            completely silent.
        """
        if file is None:
            file = _get_stdout()

        if not ipython_widget and not isatty(file):
            self.update = self._silent_update
            self._silent = True
        else:
            self._silent = False

        if isiterable(total_or_items):
            self._items = iter(total_or_items)
            self._total = len(total_or_items)
        else:
            try:
                self._total = int(total_or_items)
            except TypeError:
                raise TypeError("First argument must be int or sequence")
            else:
                self._items = iter(range(self._total))

        self._file = file
        self._start_time = time.time()
        self._human_total = human_file_size(self._total)
        self._ipython_widget = ipython_widget

        self._signal_set = False
        if not ipython_widget:
            self._should_handle_resize = _CAN_RESIZE_TERMINAL and self._file.isatty()
            self._handle_resize()
            if self._should_handle_resize:
                signal.signal(signal.SIGWINCH, self._handle_resize)
                self._signal_set = True

        self.update(0)

    def _handle_resize(self, signum=None, frame=None):
        terminal_width = terminal_size(self._file)[1]
        self._bar_length = terminal_width - 37

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        if not self._silent:
            if exc_type is None:
                self.update(self._total)
            self._file.write("\n")
            self._file.flush()
            if self._signal_set:
                signal.signal(signal.SIGWINCH, signal.SIG_DFL)

    def __iter__(self):
        return self

    def __next__(self):
        try:
            rv = next(self._items)
        except StopIteration:
            self.__exit__(None, None, None)
            raise
        else:
            self.update()
            return rv

    def update(self, value=None):
        """
        Update progress bar via the console or notebook accordingly.
        """
        # Update self.value
        if value is None:
            value = self._current_value + 1
        self._current_value = value

        # Choose the appropriate environment
        if self._ipython_widget:
            self._update_ipython_widget(value)
        else:
            self._update_console(value)

    def _update_console(self, value=None):
        """
        Update the progress bar to the given value (out of the total
        given to the constructor).
        """
        if self._total == 0:
            frac = 1.0
        else:
            frac = float(value) / float(self._total)

        file = self._file
        write = file.write

        if frac > 1:
            bar_fill = int(self._bar_length)
        else:
            bar_fill = int(float(self._bar_length) * frac)
        write("\r|")
        color_print("=" * bar_fill, "blue", file=file, end="")
        if bar_fill < self._bar_length:
            color_print(">", "green", file=file, end="")
            write("-" * (self._bar_length - bar_fill - 1))
        write("|")

        if value >= self._total:
            t = time.time() - self._start_time
            prefix = "     "
        elif value <= 0:
            t = None
            prefix = ""
        else:
            t = ((time.time() - self._start_time) * (1.0 - frac)) / frac
            prefix = " ETA "
        write(f" {human_file_size(value):>4s}/{self._human_total:>4s}")
        write(f" ({frac:>6.2%})")
        write(prefix)
        if t is not None:
            write(human_time(t))
        self._file.flush()

    def _update_ipython_widget(self, value=None):
        """
        Update the progress bar to the given value (out of a total
        given to the constructor).

        This method is for use in the IPython notebook 2+.
        """
        # Create and display an empty progress bar widget,
        # if none exists.
        if not hasattr(self, "_widget"):
            # Import only if an IPython widget, i.e., widget in iPython NB
            from IPython import version_info

            if version_info[0] < 4:
                from IPython.html import widgets

                self._widget = widgets.FloatProgressWidget()
            else:
                _IPython.get_ipython()
                from ipywidgets import widgets

                self._widget = widgets.FloatProgress()
            from IPython.display import display

            display(self._widget)
            self._widget.value = 0

        # Calculate percent completion, and update progress bar
        frac = value / self._total
        self._widget.value = frac * 100
        self._widget.description = f" ({frac:>6.2%})"

    def _silent_update(self, value=None):
        pass

    @classmethod
    def map(
        cls,
        function,
        items,
        multiprocess=False,
        file=None,
        step=100,
        ipython_widget=False,
        multiprocessing_start_method=None,
    ):
        """Map function over items while displaying a progress bar with percentage complete.

        The map operation may run in arbitrary order on the items, but the results are
        returned in sequential order.

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

        multiprocess : bool, int, optional
            If `True`, use the `multiprocessing` module to distribute each task
            to a different processor core. If a number greater than 1, then use
            that number of cores.

        ipython_widget : bool, optional
            If `True`, the progress bar will display as an IPython
            notebook widget.

        file : writable file-like, optional
            The file to write the progress bar to.  Defaults to
            `sys.stdout`.  If ``file`` is not a tty (as determined by
            calling its `isatty` member, if any), the scrollbar will
            be completely silent.

        step : int, optional
            Update the progress bar at least every *step* steps (default: 100).
            If ``multiprocess`` is `True`, this will affect the size
            of the chunks of ``items`` that are submitted as separate tasks
            to the process pool.  A large step size may make the job
            complete faster if ``items`` is very long.

        multiprocessing_start_method : str, optional
            Useful primarily for testing; if in doubt leave it as the default.
            When using multiprocessing, certain anomalies occur when starting
            processes with the "spawn" method (the only option on Windows);
            other anomalies occur with the "fork" method (the default on
            Linux).
        """
        if multiprocess:
            function = _mapfunc(function)
            items = list(enumerate(items))

        results = cls.map_unordered(
            function,
            items,
            multiprocess=multiprocess,
            file=file,
            step=step,
            ipython_widget=ipython_widget,
            multiprocessing_start_method=multiprocessing_start_method,
        )

        if multiprocess:
            _, results = zip(*sorted(results))
            results = list(results)

        return results

    @classmethod
    def map_unordered(
        cls,
        function,
        items,
        multiprocess=False,
        file=None,
        step=100,
        ipython_widget=False,
        multiprocessing_start_method=None,
    ):
        """Map function over items, reporting the progress.

        Does a `map` operation while displaying a progress bar with
        percentage complete. The map operation may run on arbitrary order
        on the items, and the results may be returned in arbitrary order.

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

        multiprocess : bool, int, optional
            If `True`, use the `multiprocessing` module to distribute each task
            to a different processor core. If a number greater than 1, then use
            that number of cores.

        ipython_widget : bool, optional
            If `True`, the progress bar will display as an IPython
            notebook widget.

        file : writable file-like, optional
            The file to write the progress bar to.  Defaults to
            `sys.stdout`.  If ``file`` is not a tty (as determined by
            calling its `isatty` member, if any), the scrollbar will
            be completely silent.

        step : int, optional
            Update the progress bar at least every *step* steps (default: 100).
            If ``multiprocess`` is `True`, this will affect the size
            of the chunks of ``items`` that are submitted as separate tasks
            to the process pool.  A large step size may make the job
            complete faster if ``items`` is very long.

        multiprocessing_start_method : str, optional
            Useful primarily for testing; if in doubt leave it as the default.
            When using multiprocessing, certain anomalies occur when starting
            processes with the "spawn" method (the only option on Windows);
            other anomalies occur with the "fork" method (the default on
            Linux).
        """
        # concurrent.futures import here to avoid import failure when running
        # in pyodide/Emscripten
        from concurrent.futures import ProcessPoolExecutor, as_completed

        results = []

        if file is None:
            file = _get_stdout()

        with cls(len(items), ipython_widget=ipython_widget, file=file) as bar:
            if bar._ipython_widget:
                chunksize = step
            else:
                default_step = max(int(float(len(items)) / bar._bar_length), 1)
                chunksize = min(default_step, step)
            if not multiprocess or multiprocess < 1:
                for i, item in enumerate(items):
                    results.append(function(item))
                    if (i % chunksize) == 0:
                        bar.update(i)
            else:
                ctx = multiprocessing.get_context(multiprocessing_start_method)
                kwargs = dict(mp_context=ctx)

                with ProcessPoolExecutor(
                    max_workers=(
                        int(multiprocess) if multiprocess is not True else None
                    ),
                    **kwargs,
                ) as p:
                    for i, f in enumerate(
                        as_completed(p.submit(function, item) for item in items)
                    ):
                        bar.update(i)
                        results.append(f.result())

        return results


class Spinner:
    """
    A class to display a spinner in the terminal.

    It is designed to be used with the ``with`` statement::

        with Spinner("Reticulating splines", "green") as s:
            for item in enumerate(items):
                s.update()
    """

    _default_unicode_chars = "◓◑◒◐"
    _default_ascii_chars = "-/|\\"

    def __init__(self, msg, color="default", file=None, step=1, chars=None):
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

        file : writable file-like, optional
            The file to write the spinner to.  Defaults to
            `sys.stdout`.  If ``file`` is not a tty (as determined by
            calling its `isatty` member, if any, or special case hacks
            to detect the IPython console), the spinner will be
            completely silent.

        step : int, optional
            Only update the spinner every *step* steps

        chars : str, optional
            The character sequence to use for the spinner
        """
        if file is None:
            file = _get_stdout()

        self._msg = msg
        self._color = color
        self._file = file
        self._step = step
        if chars is None:
            if conf.unicode_output:
                chars = self._default_unicode_chars
            else:
                chars = self._default_ascii_chars
        self._chars = chars

        self._silent = not isatty(file)

        if self._silent:
            self._iter = self._silent_iterator()
        else:
            self._iter = self._iterator()

    def _iterator(self):
        chars = self._chars
        index = 0
        file = self._file
        write = file.write
        flush = file.flush
        try_fallback = True

        while True:
            write("\r")
            color_print(self._msg, self._color, file=file, end="")
            write(" ")
            try:
                if try_fallback:
                    write = _write_with_fallback(chars[index], write, file)
                else:
                    write(chars[index])
            except UnicodeError:
                # If even _write_with_fallback failed for any reason just give
                # up on trying to use the unicode characters
                chars = self._default_ascii_chars
                write(chars[index])
                try_fallback = False  # No good will come of using this again
            flush()
            yield

            for i in range(self._step):
                yield

            index = (index + 1) % len(chars)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        file = self._file
        write = file.write
        flush = file.flush

        if not self._silent:
            write("\r")
            color_print(self._msg, self._color, file=file, end="")
        if exc_type is None:
            color_print(" [Done]", "green", file=file)
        else:
            color_print(" [Failed]", "red", file=file)
        flush()

    def __iter__(self):
        return self

    def __next__(self):
        next(self._iter)

    def update(self, value=None):
        """Update the spin wheel in the terminal.

        Parameters
        ----------
        value : int, optional
            Ignored (present just for compatibility with `ProgressBar.update`).

        """
        next(self)

    def _silent_iterator(self):
        color_print(self._msg, self._color, file=self._file, end="")
        self._file.flush()

        while True:
            yield


class ProgressBarOrSpinner:
    """
    A class that displays either a `ProgressBar` or `Spinner`
    depending on whether the total size of the operation is
    known or not.

    It is designed to be used with the ``with`` statement::

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

    def __init__(self, total, msg, color="default", file=None):
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
            The color of ``msg``, if any.  Must be an ANSI terminal
            color name.  Must be one of: black, red, green, brown,
            blue, magenta, cyan, lightgrey, default, darkgrey,
            lightred, lightgreen, yellow, lightblue, lightmagenta,
            lightcyan, white.

        file : writable file-like, optional
            The file to write the to.  Defaults to `sys.stdout`.  If
            ``file`` is not a tty (as determined by calling its `isatty`
            member, if any), only ``msg`` will be displayed: the
            `ProgressBar` or `Spinner` will be silent.
        """
        if file is None:
            file = _get_stdout()

        if total is None or not isatty(file):
            self._is_spinner = True
            self._obj = Spinner(msg, color=color, file=file)
        else:
            self._is_spinner = False
            color_print(msg, color, file=file)
            self._obj = ProgressBar(total, file=file)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        return self._obj.__exit__(exc_type, exc_value, traceback)

    def update(self, value):
        """
        Update the progress bar to the given value (out of the total
        given to the constructor.
        """
        self._obj.update(value)


def print_code_line(line, col=None, file=None, tabwidth=8, width=70):
    """
    Prints a line of source code, highlighting a particular character
    position in the line.  Useful for displaying the context of error
    messages.

    If the line is more than ``width`` characters, the line is truncated
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
        The character in the line to highlight.  ``col`` must be less
        than ``len(line)``.

    file : writable file-like, optional
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
    if file is None:
        file = _get_stdout()

    if conf.unicode_output:
        ellipsis = "…"
    else:
        ellipsis = "..."

    write = file.write

    if col is not None:
        if col >= len(line):
            raise ValueError("col must be less the the line length.")
        ntabs = line[:col].count("\t")
        col += ntabs * (tabwidth - 1)

    line = line.rstrip("\n")
    line = line.replace("\t", " " * tabwidth)

    if col is not None and col > width:
        new_col = min(width // 2, len(line) - col)
        offset = col - new_col
        line = line[offset + len(ellipsis) :]
        width -= len(ellipsis)
        new_col = col
        col -= offset
        color_print(ellipsis, "darkgrey", file=file, end="")

    if len(line) > width:
        write(line[: width - len(ellipsis)])
        color_print(ellipsis, "darkgrey", file=file)
    else:
        write(line)
        write("\n")

    if col is not None:
        write(" " * col)
        color_print("^", "red", file=file)


# The following four Getch* classes implement unbuffered character reading from
# stdin on Windows, linux, MacOSX.  This is taken directly from ActiveState
# Code Recipes:
# http://code.activestate.com/recipes/134892-getch-like-unbuffered-character-reading-from-stdin/
#


class Getch:
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


class _GetchUnix:
    def __init__(self):
        import sys  # noqa: F401

        # import termios now or else you'll get the Unix
        # version on the Mac
        import termios  # noqa: F401
        import tty  # noqa: F401

    def __call__(self):
        import sys
        import termios
        import tty

        fd = sys.stdin.fileno()
        old_settings = termios.tcgetattr(fd)
        try:
            tty.setraw(sys.stdin.fileno())
            ch = sys.stdin.read(1)
        finally:
            termios.tcsetattr(fd, termios.TCSADRAIN, old_settings)
        return ch


class _GetchWindows:
    def __init__(self):
        import msvcrt  # noqa: F401

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
            return ""
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
