"""
Utilities for prettifying output to the console.
"""
import sys
import time


def color_print(s, color='default', bold=False, italic=False, file=sys.stdout,
                end='\n'):
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

    italic : bool, optional
        When `True`, use italic font.

    file : writeable file-like object, optional
        Where to write to.  Defaults to `sys.stdout`.

    end : str, optional
        The ending of the message.  Defaults to ``\\n``.  The end will
        be printed after resetting any color or font state.
    """
    # TODO: Be smart about when to use and not use color.

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

    if bold:
        file.write('\033[1m')
    if italic:
        file.write('\033[3m')
    file.write('\033[%sm' % color_code)
    file.write(s)
    file.write('\033[0m')
    if end is not None:
        file.write(end)


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
            if parts == 2:
                rows, cols = parts
                terminal_width = int(cols)
        self._bar_length = terminal_width - 29 - (num_length * 2)
        self._num_format = '| %%%dd/%%%dd ' % (num_length, num_length)
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
        bar_fill = int(float(self._bar_length) * frac)
        self._file.write('|')
        color_print('=' * bar_fill, color='blue',
                    file=self._file, end='')
        if bar_fill < self._bar_length:
            color_print('>', color='green',
                        file=self._file, end='')
            color_print('-' * (self._bar_length - bar_fill - 1),
                        file=self._file, end='')
        if value >= self._total:
            t = time.time() - self._start_time
            eta = '    %02d:%02d' % (t / 60, t % 60)
        elif value <= 0:
            eta = ''
        else:
            t = ((time.time() - self._start_time) * (1.0 - frac)) / frac
            eta = 'ETA %02d:%02d' % (t / 60, t % 60)
        self._file.write(self._num_format % (value, self._total))
        self._file.write('(%6s%%) ' % ('%.2f' % (frac * 100.0)))
        self._file.write(eta)
        self._file.write('\r')
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
    """
    with ProgressBar(len(items), file=file) as bar:
        for item in items:
            yield item
            bar.update()
