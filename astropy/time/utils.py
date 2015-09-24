# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Time utilites.

In particular, routines to do basic arithmatic on numbers represented by two
doubles, using the procedure of Shewchuk, 1997, Discrete & Computational
Geometry 18(3):305-363 -- http://www.cs.berkeley.edu/~jrs/papers/robustr.pdf
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import numpy as np


def day_frac(val1, val2, factor=1., divisor=1.):
    """
    Return the sum of ``val1`` and ``val2`` as two float64s, an integer part
    and the fractional remainder.  If ``factor`` is not 1.0 then multiply the
    sum by ``factor``.  If ``divisor`` is not 1.0 then divide the sum by
    ``divisor``.

    The arithmetic is all done with exact floating point operations so no
    precision is lost to rounding error.  This routine assumes the sum is less
    than about 1e16, otherwise the ``frac`` part will be greater than 1.0.

    Returns
    -------
    day, frac : float64
        Integer and fractional part of val1 + val2.
    """
    # Add val1 and val2 exactly, returning the result as two float64s.
    # The first is the approximate sum (with some floating point error)
    # and the second is the error of the float64 sum.
    sum12, err12 = two_sum(val1, val2)

    if np.any(factor != 1.):
        sum12, carry = two_product(sum12, factor)
        carry += err12 * factor
        sum12, err12 = two_sum(sum12, carry)

    if np.any(divisor != 1.):
        q1 = sum12 / divisor
        p1, p2 = two_product(q1, divisor)
        d1, d2 = two_sum(sum12, -p1)
        d2 += err12
        d2 -= p2
        q2 = (d1 + d2) / divisor  # 3-part float fine here; nothing can be lost
        sum12, err12 = two_sum(q1, q2)

    # get integer fraction
    day = np.round(sum12)
    extra, frac = two_sum(sum12, -day)
    frac += extra + err12
    return day, frac


def two_sum(a, b):
    """
    Add ``a`` and ``b`` exactly, returning the result as two float64s.
    The first is the approximate sum (with some floating point error)
    and the second is the error of the float64 sum.

    Using the procedure of Shewchuk, 1997,
    Discrete & Computational Geometry 18(3):305-363
    http://www.cs.berkeley.edu/~jrs/papers/robustr.pdf

    Returns
    -------
    sum, err : float64
        Approximate sum of a + b and the exact floating point error
    """
    x = a + b
    eb = x - a
    eb = b - eb
    ea = x - b
    ea = a - ea
    return x, ea + eb


def two_product(a, b):
    """
    Multiple ``a`` and ``b`` exactly, returning the result as two float64s.
    The first is the approximate product (with some floating point error)
    and the second is the error of the float64 product.

    Uses the procedure of Shewchuk, 1997,
    Discrete & Computational Geometry 18(3):305-363
    http://www.cs.berkeley.edu/~jrs/papers/robustr.pdf

    Returns
    -------
    prod, err : float64
        Approximate product a * b and the exact floating point error
    """
    x = a * b
    ah, al = split(a)
    bh, bl = split(b)
    y1 = ah * bh
    y = x - y1
    y2 = al * bh
    y -= y2
    y3 = ah * bl
    y -= y3
    y4 = al * bl
    y = y4 - y
    return x, y


def split(a):
    """
    Split float64 in two aligned parts.

    Uses the procedure of Shewchuk, 1997,
    Discrete & Computational Geometry 18(3):305-363
    http://www.cs.berkeley.edu/~jrs/papers/robustr.pdf

    """
    c = 134217729. * a  # 2**27+1.
    abig = c - a
    ah = c - abig
    al = a - ah
    return ah, al
