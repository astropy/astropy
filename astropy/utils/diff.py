import difflib
import functools
import sys
from textwrap import indent

import numpy as np

__all__ = [
    "diff_values",
    "report_diff_values",
    "where_not_allclose",
]


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


def _ignore_astropy_terminal_size(func):
    @functools.wraps(func)
    def inner(*args, **kwargs):
        from astropy import conf

        with conf.set_temp("max_width", -1), conf.set_temp("max_lines", -1):
            return func(*args, **kwargs)

    return inner


@_ignore_astropy_terminal_size
def report_diff_values(a, b, fileobj=sys.stdout, indent_width=0, rtol=0.0, atol=0.0):
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

    rtol, atol : float
        Relative and absolute tolerances as accepted by
        :func:`numpy.allclose`.

    Returns
    -------
    identical : bool
        `True` if no diff, else `False`.

    """
    indent_prefix = indent_width * "  "
    if isinstance(a, np.ndarray) and isinstance(b, np.ndarray):
        if a.shape != b.shape:
            fileobj.write(indent("  Different array shapes:\n", indent_prefix))
            report_diff_values(
                str(a.shape),
                str(b.shape),
                fileobj=fileobj,
                indent_width=indent_width + 1,
            )
            return False

        if np.issubdtype(a.dtype, np.floating) and np.issubdtype(b.dtype, np.floating):
            diff_indices = np.transpose(where_not_allclose(a, b, rtol=rtol, atol=atol))
        else:
            diff_indices = np.transpose(np.where(a != b))

        num_diffs = diff_indices.shape[0]

        for idx in diff_indices[:3]:
            lidx = idx.tolist()
            fileobj.write(indent(f"  at {lidx!r}:\n", indent_prefix))
            report_diff_values(
                a[tuple(idx)],
                b[tuple(idx)],
                fileobj=fileobj,
                indent_width=indent_width + 1,
                rtol=rtol,
                atol=atol,
            )

        if num_diffs > 3:
            fileobj.write(
                indent(f"  ...and at {num_diffs - 3:d} more indices.\n", indent_prefix)
            )
            return False

        return num_diffs == 0

    typea = type(a)
    typeb = type(b)

    if typea == typeb:
        lnpad = " "
        sign_a = "a>"
        sign_b = "b>"
        a = str(a)
        b = str(b)
    else:
        padding = max(len(typea.__name__), len(typeb.__name__)) + 3
        lnpad = (padding + 1) * " "
        sign_a = ("(" + typea.__name__ + ") ").rjust(padding) + "a>"
        sign_b = ("(" + typeb.__name__ + ") ").rjust(padding) + "b>"

        is_a_str = isinstance(a, str)
        is_b_str = isinstance(b, str)
        a = repr(a) if is_a_str and not is_b_str else str(a)
        b = repr(b) if is_b_str and not is_a_str else str(b)

    identical = True

    for line in difflib.ndiff(a.splitlines(), b.splitlines()):
        if line[0] == "-":
            identical = False
            line = sign_a + line[1:]
        elif line[0] == "+":
            identical = False
            line = sign_b + line[1:]
        else:
            line = lnpad + line
        fileobj.write(indent("  {}\n".format(line.rstrip("\n")), indent_prefix))

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
