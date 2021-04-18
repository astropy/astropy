import difflib
import functools
import sys
import numbers

import numpy as np

from .misc import indent

__all__ = ['fixed_width_indent', 'diff_values', 'report_diff_values',
           'where_not_allclose']


# Smaller default shift-width for indent
fixed_width_indent = functools.partial(indent, width=2)


def diff_values(a, b, rtol=0.0, atol=0.0):
    """
    Diff two scalar values. If both values are floats, they are compared to
    within the given absolute and relative tolerance.

    Parameters
    ----------
    a, b : int, float, str
        Scalar values to compare.

    rtol, atol : float
        Relative and absolute tolerances as accepted by
        :func:`numpy.allclose`.

    Returns
    -------
    is_different : bool
        `True` if they are different, else `False`.

    """
    if isinstance(a, float) and isinstance(b, float):
        if np.isnan(a) and np.isnan(b):
            return False
        return not np.allclose(a, b, rtol=rtol, atol=atol)
    else:
        return a != b


def report_diff_values(a, b, fileobj=sys.stdout, indent_width=0):
    """
    Write a diff report between two values to the specified file-like object.

    Parameters
    ----------
    a, b
        Values to compare. Anything that can be turned into strings
        and compared using :py:mod:`difflib` should work.

    fileobj : object
        File-like object to write to.
        The default is ``sys.stdout``, which writes to terminal.

    indent_width : int
        Character column(s) to indent.

    Returns
    -------
    identical : bool
        `True` if no diff, else `False`.

    """
    if isinstance(a, np.ndarray) and isinstance(b, np.ndarray):
        if a.shape != b.shape:
            fileobj.write(
                fixed_width_indent('  Different array shapes:\n',
                                   indent_width))
            report_diff_values(str(a.shape), str(b.shape), fileobj=fileobj,
                               indent_width=indent_width + 1)
            return False

        diff_indices = np.transpose(np.where(a != b))
        num_diffs = diff_indices.shape[0]

        for idx in diff_indices[:3]:
            lidx = idx.tolist()
            fileobj.write(
                fixed_width_indent(f'  at {lidx!r}:\n', indent_width))
            report_diff_values(a[tuple(idx)], b[tuple(idx)], fileobj=fileobj,
                               indent_width=indent_width + 1)

        if num_diffs > 3:
            fileobj.write(fixed_width_indent(
                f'  ...and at {num_diffs - 3:d} more indices.\n',
                indent_width))
            return False

        return num_diffs == 0

    typea = type(a)
    typeb = type(b)

    if typea == typeb:
        lnpad = ' '
        sign_a = 'a>'
        sign_b = 'b>'
        if isinstance(a, numbers.Number):
            a = repr(a)
            b = repr(b)
        else:
            a = str(a)
            b = str(b)
    else:
        padding = max(len(typea.__name__), len(typeb.__name__)) + 3
        lnpad = (padding + 1) * ' '
        sign_a = ('(' + typea.__name__ + ') ').rjust(padding) + 'a>'
        sign_b = ('(' + typeb.__name__ + ') ').rjust(padding) + 'b>'

        is_a_str = isinstance(a, str)
        is_b_str = isinstance(b, str)
        a = (repr(a) if ((is_a_str and not is_b_str) or
                         (not is_a_str and isinstance(a, numbers.Number)))
             else str(a))
        b = (repr(b) if ((is_b_str and not is_a_str) or
                         (not is_b_str and isinstance(b, numbers.Number)))
             else str(b))

    identical = True

    for line in difflib.ndiff(a.splitlines(), b.splitlines()):
        if line[0] == '-':
            identical = False
            line = sign_a + line[1:]
        elif line[0] == '+':
            identical = False
            line = sign_b + line[1:]
        else:
            line = lnpad + line
        fileobj.write(fixed_width_indent(
            '  {}\n'.format(line.rstrip('\n')), indent_width))

    return identical


def where_not_allclose(a, b, rtol=1e-5, atol=1e-8):
    """
    A version of :func:`numpy.allclose` that returns the indices
    where the two arrays differ, instead of just a boolean value.

    Parameters
    ----------
    a, b : array-like
        Input arrays to compare.

    rtol, atol : float
        Relative and absolute tolerances as accepted by
        :func:`numpy.allclose`.

    Returns
    -------
    idx : tuple of array
        Indices where the two arrays differ.

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
