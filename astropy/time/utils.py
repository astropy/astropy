# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Time utilities.

In particular, routines to do basic arithmetic on numbers represented by two
doubles, using the procedure of Shewchuk, 1997, Discrete & Computational
Geometry 18(3):305-363 -- http://www.cs.berkeley.edu/~jrs/papers/robustr.pdf

Furthermore, some helper routines to turn strings and other types of
objects into two values, and vice versa.
"""
import decimal

import numpy as np
import astropy.units as u


def day_frac(val1, val2, factor=None, divisor=None):
    """Return the sum of ``val1`` and ``val2`` as two float64s.

    The returned floats are an integer part and the fractional remainder,
    with the latter guaranteed to be within -0.5 and 0.5 (inclusive on
    either side, as the integer is rounded to even).

    The arithmetic is all done with exact floating point operations so no
    precision is lost to rounding error.  It is assumed the sum is less
    than about 1e16, otherwise the remainder will be greater than 1.0.

    Parameters
    ----------
    val1, val2 : array of float
        Values to be summed.
    factor : float, optional
        If given, multiply the sum by it.
    divisor : float, optional
        If given, divide the sum by it.

    Returns
    -------
    day, frac : float64
        Integer and fractional part of val1 + val2.
    """
    # Add val1 and val2 exactly, returning the result as two float64s.
    # The first is the approximate sum (with some floating point error)
    # and the second is the error of the float64 sum.
    sum12, err12 = two_sum(val1, val2)

    if factor is not None:
        sum12, carry = two_product(sum12, factor)
        carry += err12 * factor
        sum12, err12 = two_sum(sum12, carry)

    if divisor is not None:
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
    # Our fraction can now have gotten >0.5 or <-0.5, which means we would
    # loose one bit of precision. So, correct for that.
    excess = np.round(frac)
    day += excess
    extra, frac = two_sum(sum12, -day)
    frac += extra + err12
    return day, frac


def quantity_day_frac(val1, val2=None):
    """Like ``day_frac``, but for quantities with units of time.

    The quantities are separately converted to days. Here, we need to take
    care with the conversion since while the routines here can do accurate
    multiplication, the conversion factor itself may not be accurate.  For
    instance, if the quantity is in seconds, the conversion factor is
    1./86400., which is not exactly representable as a float.

    To work around this, for conversion factors less than unity, rather than
    multiply by that possibly inaccurate factor, the value is divided by the
    conversion factor of a day to that unit (i.e., by 86400. for seconds).  For
    conversion factors larger than 1, such as 365.25 for years, we do just
    multiply.  With this scheme, one has precise conversion factors for all
    regular time units that astropy defines.  Note, however, that it does not
    necessarily work for all custom time units, and cannot work when conversion
    to time is via an equivalency.  For those cases, one remains limited by the
    fact that Quantity calculations are done in double precision, not in
    quadruple precision as for time.
    """
    if val2 is not None:
        res11, res12 = quantity_day_frac(val1)
        res21, res22 = quantity_day_frac(val2)
        # This summation is can at most lose 1 ULP in the second number.
        return res11 + res21, res12 + res22

    try:
        factor = val1.unit.to(u.day)
    except Exception:
        # Not a simple scaling, so cannot do the full-precision one.
        # But at least try normal conversion, since equivalencies may be set.
        return val1.to_value(u.day), 0.

    if factor == 1.:
        return day_frac(val1.value, 0.)

    if factor > 1:
        return day_frac(val1.value, 0., factor=factor)
    else:
        divisor = u.day.to(val1.unit)
        return day_frac(val1.value, 0., divisor=divisor)


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
    eb = x - a  # bvirtual in Shewchuk
    ea = x - eb  # avirtual in Shewchuk
    eb = b - eb  # broundoff in Shewchuk
    ea = a - ea  # aroundoff in Shewchuk
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


_enough_decimal_places = 34  # to represent two doubles


def longdouble_to_twoval(val1, val2=None):
    if val2 is None:
        val2 = val1.dtype.type(0.)
    else:
        best_type = np.result_type(val1.dtype, val2.dtype)
        val1 = val1.astype(best_type, copy=False)
        val2 = val2.astype(best_type, copy=False)

    # day_frac is independent of dtype, as long as the dtype
    # are the same and no factor or divisor is given.
    i, f = day_frac(val1, val2)
    return i.astype(float, copy=False), f.astype(float, copy=False)


def decimal_to_twoval1(val1, val2=None):
    with decimal.localcontext() as ctx:
        ctx.prec = _enough_decimal_places
        d = decimal.Decimal(val1)
        i = round(d)
        f = d - i
    return float(i), float(f)


def bytes_to_twoval1(val1, val2=None):
    return decimal_to_twoval1(val1.decode('ascii'))


def twoval_to_longdouble(val1, val2):
    return val1.astype(np.longdouble) + val2.astype(np.longdouble)


def twoval_to_decimal1(val1, val2):
    with decimal.localcontext() as ctx:
        ctx.prec = _enough_decimal_places
        return decimal.Decimal(val1) + decimal.Decimal(val2)


def twoval_to_string1(val1, val2, fmt):
    if val2 == 0.:
        # For some formats, only a single float is really used.
        # For those, let numpy take care of correct number of digits.
        return str(val1)

    result = format(twoval_to_decimal1(val1, val2), fmt).strip('0')
    if result[-1] == '.':
        result += '0'
    return result


def twoval_to_bytes1(val1, val2, fmt):
    return twoval_to_string1(val1, val2, fmt).encode('ascii')


decimal_to_twoval = np.vectorize(decimal_to_twoval1)
bytes_to_twoval = np.vectorize(bytes_to_twoval1)
twoval_to_decimal = np.vectorize(twoval_to_decimal1)
twoval_to_string = np.vectorize(twoval_to_string1, excluded='fmt')
twoval_to_bytes = np.vectorize(twoval_to_bytes1, excluded='fmt')
