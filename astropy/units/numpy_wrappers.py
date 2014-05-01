import numpy as np

from .quantity import Quantity


def assert_allclose(actual, desired, rtol=1.e-7, atol=0, err_msg='', verbose=True):
    """
    Raise an assertion if two objects are not equal up to desired tolerance.

    This is a :class:`~astropy.units.Quantity`-aware version of
    :func:`numpy.testing.assert_allclose`.
    """
    if isinstance(actual, Quantity) and isinstance(desired, Quantity):
        np.testing.assert_allclose(actual.value, desired.to(actual.unit).value,
                                   rtol=rtol, atol=atol, err_msg=err_msg, verbose=verbose)
    elif isinstance(actual, Quantity):
        raise TypeError("If `actual` is a Quantity, `desired` should also be a Quantity")
    elif isinstance(desired, Quantity):
        raise TypeError("If `desired` is a Quantity, `actual` should also be a Quantity")
    else:
        np.testing.assert_allclose(actual, desired,
                                   rtol=rtol, atol=atol, err_msg=err_msg, verbose=verbose)


def linspace(start, stop, num=50, endpoint=True, retstep=False, dtype=None):
    """
    Return evenly spaced numbers over a specified interval.

    This is a :class:`~astropy.units.Quantity`-aware version of
    :func:`numpy.linspace`.
    """
    if isinstance(start, Quantity) and isinstance(stop, Quantity):
        return np.linspace(start.value, stop.to(start.unit).value, num=num,
                           endpoint=endpoint, retstep=retstep, dtype=dtype) * start.unit
    elif isinstance(start, Quantity):
        raise TypeError("If `start` is a Quantity, `stop` should also be a Quantity")
    elif isinstance(stop, Quantity):
        raise TypeError("If `stop` is a Quantity, `start` should also be a Quantity")
    else:
        return np.linspace(start, stop, num=num,
                           endpoint=endpoint, retstep=retstep, dtype=dtype)
