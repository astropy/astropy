import difflib
import functools
import operator
from functools import reduce
from itertools import islice

import numpy as np

from .misc import indent

__all__ = ['fixed_width_indent', 'diff_values', 'report_diff_values',
           'where_not_allclose']

# Smaller default shift-width for indent:
fixed_width_indent = functools.partial(indent, width=2)


def diff_values(a, b, rtol=0.0, atol=0.0):
    """
    Diff two scalar values.  If both values are floats they are compared to
    within the given absolute and relative tolerance.
    """

    if isinstance(a, float) and isinstance(b, float):
        if np.isnan(a) and np.isnan(b):
            return False
        return not np.allclose(a, b, rtol=rtol, atol=atol)
    else:
        return a != b


def report_diff_values(fileobj, a, b, ind=0):
    """Write a diff between two values to the specified file-like object."""

    typea = type(a)
    typeb = type(b)

    if (isinstance(a, str) and not isinstance(b, str)):
        a = repr(a).lstrip('u')
    elif (isinstance(b, str) and not isinstance(a, str)):
        b = repr(b).lstrip('u')

    if isinstance(a, (int, float, complex, np.number)):
        a = repr(a)

    if isinstance(b, (int, float, complex, np.number)):
        b = repr(b)

    if isinstance(a, np.ndarray) and isinstance(b, np.ndarray):
        diff_indices = np.where(a != b)
        num_diffs = reduce(operator.mul, map(len, diff_indices), 1)
        for idx in islice(zip(*diff_indices), 3):
            fileobj.write(
                fixed_width_indent('  at {!r}:\n'.format(list(idx)), ind))
            report_diff_values(fileobj, a[idx], b[idx], ind=ind + 1)

        if num_diffs > 3:
            fileobj.write(fixed_width_indent('  ...and at {} more indices.\n'
                                             .format(num_diffs - 3), ind))
        return

    padding = max(len(typea.__name__), len(typeb.__name__)) + 3

    for line in difflib.ndiff(str(a).splitlines(), str(b).splitlines()):
        if line[0] == '-':
            line = 'a>' + line[1:]
            if typea != typeb:
                typename = '(' + typea.__name__ + ') '
                line = typename.rjust(padding) + line

        elif line[0] == '+':
            line = 'b>' + line[1:]
            if typea != typeb:
                typename = '(' + typeb.__name__ + ') '
                line = typename.rjust(padding) + line
        else:
            line = ' ' + line
            if typea != typeb:
                line = ' ' * padding + line
        fileobj.write(
            fixed_width_indent('  {}\n'.format(line.rstrip('\n')), ind))


def where_not_allclose(a, b, rtol=1e-5, atol=1e-8):
    """
    A version of numpy.allclose that returns the indices where the two arrays
    differ, instead of just a boolean value.
    """

    # Create fixed mask arrays to handle INF and NaN; currently INF and NaN
    # are handled as equivalent
    if not np.all(np.isfinite(a)):
        a = np.ma.fix_invalid(a).data
    if not np.all(np.isfinite(b)):
        b = np.ma.fix_invalid(b).data

    if atol == 0.0 and rtol == 0.0:
        # Use a faster comparison for the most simple (and common) case
        return np.where(a != b)
    return np.where(np.abs(a - b) > (atol + rtol * np.abs(b)))
