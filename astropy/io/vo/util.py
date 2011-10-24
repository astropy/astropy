# Copyright (C) 2008-2010 Association of Universities for Research in Astronomy (AURA)

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:

#     1. Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.

#     2. Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.

#     3. The name of AURA and its representatives may not be used to
#       endorse or promote products derived from this software without
#       specific prior written permission.

# THIS SOFTWARE IS PROVIDED BY AURA ``AS IS'' AND ANY EXPRESS OR IMPLIED
# WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL AURA BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
# OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
# TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
# USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
# DAMAGE.

"""
Various utilities and cookbook-like things.
"""

from __future__ import division, with_statement, absolute_import

# STDLIB
import collections
from distutils import version
import io
import re
import sys
import time

IS_PY3K = sys.hexversion >= 0x03000000


class HomogeneousList(list):
    """
    A subclass of list that contains only elements of a given type or
    types.
    """
    def __init__(self, types, values=[]):
        """
        *types* is a sequence of acceptable types.

        *values* (optional) is an initial set of values.
        """
        self._types = types
        list.__init__(self, values)

    def _assert(self, x):
        if not isinstance(x, self._types):
            raise TypeError(
                "homogeneous list must contain only objects of type '%s'" %
                (self._types))

    def __iadd__(self, other):
        for x in other:
            self._assert(x)
        return list.__iadd__(self, other)

    def __setitem__(self, x):
        self._assert(x)
        return list.__setitem__(self, x)

    def append(self, x):
        self._assert(x)
        return list.append(self, x)

    def insert(self, i, x):
        self._assert(x)
        return list.insert(self, i, x)


def convert_to_writable_filelike(fd):
    """
    Returns a writable file-like object suitable for streaming output.

    *fd* may be:

      - a file path, in which case it is opened, and the :meth:`write`
        method on the file object is returned.

      - an object with a :meth:`write` method, in which case that
        method is returned.
    """
    if isinstance(fd, basestring):
        if IS_PY3K:
            if fd.endswith('.gz'):
                import gzip
                fd = gzip.GzipFile(fd, 'wb')
                fd = io.TextIOWrapper(fd, encoding='utf8')
            else:
                fd = io.open(fd, 'w', encoding='utf8')
        else:
            if fd.endswith('.gz'):
                import gzip
                fd = gzip.GzipFile(fd, 'wb')
            else:
                fd = open(fd, 'wb')
        write = fd
    elif hasattr(fd, 'write'):
        write = fd
        assert is_callable(fd.write)
    else:
        raise TypeError("Can not be coerced to writable file-like object")

    return write


def convert_to_fd_or_read_function(fd):
    """
    Returns a function suitable for streaming input, or a file object.

    *fd* may be:

      - a file object, in which case it is returned verbatim.

      - a function that reads from a stream, in which case it is
        returned verbatim.

      - a file path, in which case it is opened.  If it ends in `.gz`,
        it is assumed to be a gzipped file, and the :meth:`read`
        method on the file object is returned.  Otherwise, the raw
        file object is returned.

      - an object with a :meth:`read` method, in which case that
        method is returned.
    """
    if IS_PY3K:
        if isinstance(fd, io.IOBase):
            return fd
    else:
        if isinstance(fd, file):
            return fd
    if is_callable(fd):
        return fd
    elif isinstance(fd, basestring):
        if fd.endswith('.gz'):
            import gzip
            fd = gzip.GzipFile(fd, 'rb')
            return fd.read
        else:
            fd = open(fd, 'rb')
            return fd
    elif hasattr(fd, 'read'):
        assert is_callable(fd.read)
        return fd.read
    else:
        raise TypeError("Can not be coerced to read function")


class maxdict(dict):
    """
    A dictionary with a maximum size.
    """
    def __init__(self, maxsize=10):
        dict.__init__(self)
        self._maxsize = maxsize
        self._keys = []

    def __setitem__(self, k, v):
        if k not in self:
            if len(self) >= self._maxsize:
                del self[self._keys[0]]
                del self._keys[0]
            self._keys.append(k)
        dict.__setitem__(self, k, v)

    def __delitem__(self, key):
        if key in self:
            self._keys.remove(key)
        dict.__delitem__(self, key)

    def clear(self):
        self._keys = []
        dict.clear(self)

    def update(self, dict=None, **kwargs):
        if dict is not None:
            for k, v in dict.iteritems():
                self[k] = v
        if len(kwargs):
            self.update(kwargs)

    def setdefault(self, key, failobj=None):
        if key not in self:
            self[key] = failobj
        return self[key]

    def pop(self, key, *args):
        if key in self:
            self._keys.remove(key)
        dict.pop(key, *args)

    def popitem(self):
        key, val = dict.popitem(self)
        self._keys.remove(key)


class memoized(object):
    """
    Decorator that caches a function's return value each time it is called.
    If called later with the same arguments, the cached value is returned, and
    not re-evaluated.
    """
    def __init__(self, func, maxsize=10):
        self.func = func
        self.cache = maxdict(maxsize=maxsize)

    def __call__(self, *args):
        try:
            return self.cache.setdefault(args, self.func(*args))
        except TypeError:
            # uncachable -- for instance, passing a list as an argument.
            # Better to not cache than to blow up entirely.
            return self.func(*args)

    def __repr__(self):
        """Return the function's docstring."""
        return self.func.__doc__


def dict_soft_update(d, u):
    """
    Like dict.update, except if the values in *u* are None, *d* is not updated.
    """
    for key, val in u.iteritems():
        if val is not None:
            d[key] = val


# <http://www.ivoa.net/Documents/REC/DM/STC-20071030.html>
stc_reference_frames = set([
    'FK4', 'FK5', 'ECLIPTIC', 'ICRS', 'GALACTIC', 'GALACTIC_I', 'GALACTIC_II',
    'SUPER_GALACTIC', 'AZ_EL', 'BODY', 'GEO_C', 'GEO_D', 'MAG', 'GSE', 'GSM',
    'SM', 'HGC', 'HGS', 'HEEQ', 'HRTN', 'HPC', 'HPR', 'HCC', 'HGI', 'MERCURY_C',
    'VENUS_C', 'LUNA_C', 'MARS_C', 'JUPITER_C_III', 'SATURN_C_III',
    'URANUS_C_III', 'NEPTUNE_C_III', 'PLUTO_C', 'MERCURY_G', 'VENUS_G',
    'LUNA_G', 'MARS_G', 'JUPITER_G_III', 'SATURN_G_III', 'URANUS_G_III',
    'NEPTUNE_G_III', 'PLUTO_G', 'UNKNOWNFrame'
    ])

def coerce_range_list_param(p, frames=None, numeric=True):
    """
    Coerces and/or verifies the object *p* into a valid
    range-list-format parameter as defined in `Section 8.7.2 of Simple
    Spectral Access Protocol
    <http://www.ivoa.net/Documents/REC/DAL/SSA-20080201.html>`_.

    *p* may be a string as passed verbatim to the service expecting a
    range-list, or a sequence.  If a sequence, each item must be
    either:

      - a numeric value

      - a named value, such as, for example, 'J' for named spectrum
        (if the *numeric* kwarg is False)

      - a 2-tuple indicating a range

      - the last item my be a string indicating the frame of reference

    The result is a tuple:

      - a string suitable for passing to a service as a range-list
        argument

      - an integer counting the number of elements

    *frames*, if provided, should be a sequence of acceptable frame of
    reference keywords.
    """
    def str_or_none(x):
        if x is None:
            return ''
        if numeric:
            x = float(x)
        return str(x)

    def numeric_or_range(x):
        if isinstance(x, tuple) and len(x) == 2:
            return '%s/%s' % (str_or_none(x[0]), str_or_none(x[1]))
        else:
            return str_or_none(x)

    def is_frame_of_reference(x):
        return isinstance(x, basestring)

    if p is None:
        return None, 0

    elif isinstance(p, (tuple, list)):
        has_frame_of_reference = len(p) > 1 and is_frame_of_reference(p[-1])
        if has_frame_of_reference:
            points = p[:-1]
        else:
            points = p[:]

        out = ','.join([numeric_or_range(x) for x in points])
        length = len(points)
        if has_frame_of_reference:
            if frames is not None and p[-1] not in frames:
                raise ValueError("'%s' is not a valid frame of reference" % p[-1])
            out += ';' + p[-1]
            length += 1

        return out, length

    elif isinstance(p, basestring):
        number = r'([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)?'
        if not numeric:
            number = r'(' + number + ')|([A-Z_]+)'
        match = re.match(
            '^' + number + r'([,/]' + number + r')+(;(?P<frame>[<A-Za-z_0-9]+))?$',
            p)

        if match is None:
            raise ValueError("'%s' is not a valid range list" % p)

        frame = match.groupdict()['frame']
        if frames is not None and frame is not None and frame not in frames:
            raise ValueError(
                "'%s' is not a valid frame of reference" % frame)
        return p, p.count(',') + p.count(';') + 1

    try:
        float(p)
        return str(p), 1
    except TypeError:
        raise ValueError("'%s' is not a valid range list" % p)


def version_compare(a, b):
    """
    Compare two version identifiers.
    """
    def version_to_tuple(v):
        if v[0].lower() == 'v':
            v = v[1:]
        return version.StrictVersion(v)
    av = version_to_tuple(a)
    bv = version_to_tuple(b)
    # Can't use cmp because it was removed from Python 3.x
    return (av > bv) - (av < bv)


def color_print(color, s, bold=False, italic=False, stream=sys.stdout,
                newline=True):
    """
    Prints the string *s* in the given *color* where *color* is an
    ANSI terminal color name.  *stream* may be any writable file-like
    object.

    TODO: Be smart about when to use and not use color.
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
        'white': '1;37'
        }

    color_code = color_mapping.get(color, '0;39')

    if bold:
        stream.write('\033[1m')
    if italic:
        stream.write('\033[3m')
    stream.write('\033[%sm' % color_code)
    stream.write(s)
    stream.write('\033[0m')
    if newline:
        stream.write('\n')


if IS_PY3K:
    def is_callable(o):
        """
        Abstracts away the different ways to test for a callable object in
        Python 2.x and 3.x.
        """
        return isinstance(o, collections.Callable)
else:
    def is_callable(o):
        """
        Abstracts away the different ways to test for a callable object in
        Python 2.x and 3.x.
        """
        return callable(o)


class ProgressBar:
    """
    A class to display a progress bar in the terminal.

    It is designed for use with the `with` statement::

       with ProgressBar(len(items)) as bar:
           for i, item in enumerate(items):
               bar.update(i)
    """
    def __init__(self, total, stream=sys.stdout):
        """
        *total* is the number of steps in the process.
        """
        self._total = total
        self._stream = stream
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
            rows, cols = stdout.split()
            terminal_width = int(cols)
        self._bar_length = terminal_width - 29 - (num_length * 2)
        self._num_format = '| %%%dd/%%%dd ' % (num_length, num_length)
        self.update(0)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        if exc_type is None:
            self.update(self._total)
        self._stream.write('\n')
        self._stream.flush()

    def update(self, value):
        """
        Update the progress bar to the given value (out of the total
        given to the constructor.
        """
        if self._total == 0:
            frac = 1.0
        else:
            frac = float(value) / float(self._total)
        bar_fill = int(float(self._bar_length) * frac)
        self._stream.write('|')
        color_print('blue', '=' * bar_fill,
                    stream=self._stream, newline=False)
        if bar_fill < self._bar_length:
            color_print('green', '>',
                        stream=self._stream, newline=False)
            color_print('default', '-' * (self._bar_length - bar_fill - 1),
                        stream=self._stream, newline=False)
        if value >= self._total:
            t = time.time() - self._start_time
            eta = '    %02d:%02d' % (t / 60, t % 60)
        elif value <= 0:
            eta = ''
        else:
            t = ((time.time() - self._start_time) * (1.0 - frac)) / frac
            eta = 'ETA %02d:%02d' % (t / 60, t % 60)
        self._stream.write(self._num_format % (value, self._total))
        self._stream.write('(%6s%%) ' % ('%.2f' % (frac * 100.0)))
        self._stream.write(eta)
        self._stream.write('\r')
        self._stream.flush()


def map_with_progress(function, items, multiprocess=False, stream=sys.stdout):
    """
    Does a *map* while displaying a progress bar with percentage
    complete.
    """
    with ProgressBar(len(items), stream=stream) as bar:
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
