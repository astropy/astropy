# Licensed under a 3-clause BSD style license - see LICENSE.rst

import functools
from math import inf
from numbers import Number

import numpy as np

from astropy.units import Quantity
from astropy.utils import isiterable
from astropy.utils.decorators import deprecated

from . import units as cu

__all__ = []  # nothing is publicly scoped

__doctest_skip__ = ["inf_like", "vectorize_if_needed"]


def vectorize_redshift_method(func=None, nin=1):
    """Vectorize a method of redshift(s).

    Parameters
    ----------
    func : callable or None
        method to wrap. If `None` returns a :func:`functools.partial`
        with ``nin`` loaded.
    nin : int
        Number of positional redshift arguments.

    Returns
    -------
    wrapper : callable
        :func:`functools.wraps` of ``func`` where the first ``nin``
        arguments are converted from |Quantity| to :class:`numpy.ndarray`.
    """
    # allow for pie-syntax & setting nin
    if func is None:
        return functools.partial(vectorize_redshift_method, nin=nin)

    @functools.wraps(func)
    def wrapper(self, *args, **kwargs):
        """
        :func:`functools.wraps` of ``func`` where the first ``nin``
        arguments are converted from |Quantity| to `numpy.ndarray` or scalar.
        """
        # process inputs
        # TODO! quantity-aware vectorization can simplify this.
        zs = [
            z if not isinstance(z, Quantity) else z.to_value(cu.redshift)
            for z in args[:nin]
        ]
        # scalar inputs
        if all(isinstance(z, (Number, np.generic)) for z in zs):
            return func(self, *zs, *args[nin:], **kwargs)
        # non-scalar. use vectorized func
        return wrapper.__vectorized__(self, *zs, *args[nin:], **kwargs)

    wrapper.__vectorized__ = np.vectorize(func)  # attach vectorized function
    # TODO! use frompyfunc when can solve return type errors

    return wrapper


@deprecated(
    since="5.0",
    message=(
        "vectorize_if_needed has been removed because it constructs a new ufunc on each"
        " call"
    ),
    alternative="use a pre-vectorized function instead for a target array 'z'",
)
def vectorize_if_needed(f, *x, **vkw):
    """Helper function to vectorize scalar functions on array inputs.

    Parameters
    ----------
    f : callable
        'f' must accept positional arguments and no mandatory keyword
        arguments.
    *x
        Arguments into ``f``.
    **vkw
        Keyword arguments into :class:`numpy.vectorize`.

    Examples
    --------
    >>> func = lambda x: x ** 2
    >>> vectorize_if_needed(func, 2)
    4
    >>> vectorize_if_needed(func, [2, 3])
    array([4, 9])
    """
    return np.vectorize(f, **vkw)(*x) if any(map(isiterable, x)) else f(*x)


@deprecated(
    since="5.0",
    message=(
        "inf_like has been removed because it duplicates "
        "functionality provided by numpy.full_like()"
    ),
    alternative="Use numpy.full_like(z, numpy.inf) instead for a target array 'z'",
)
def inf_like(x):
    """Return the shape of x with value infinity and dtype='float'.

    Preserves 'shape' for both array and scalar inputs.
    But always returns a float array, even if x is of integer type.

    Parameters
    ----------
    x : scalar or array-like
        Must work with functions `numpy.isscalar` and `numpy.full_like` (if `x`
        is not a scalar`

    Returns
    -------
    `math.inf` or ndarray[float] thereof
        Returns a scalar `~math.inf` if `x` is a scalar, an array of floats
        otherwise.

    Examples
    --------
    >>> inf_like(0.)  # float scalar
    inf
    >>> inf_like(1)  # integer scalar should give float output
    inf
    >>> inf_like([0., 1., 2., 3.])  # float list
    array([inf, inf, inf, inf])
    >>> inf_like([0, 1, 2, 3])  # integer list should give float output
    array([inf, inf, inf, inf])
    """
    return inf if np.isscalar(x) else np.full_like(x, inf, dtype=float)


def aszarr(z):
    """
    Redshift as a `~numbers.Number` or `~numpy.ndarray` / |Quantity| / |Column|.
    Allows for any ndarray ducktype by checking for attribute "shape".
    """
    if isinstance(z, (Number, np.generic)):  # scalars
        return z
    elif hasattr(z, "shape"):  # ducktypes NumPy array
        if hasattr(z, "unit"):  # Quantity Column
            return (z << cu.redshift).value  # for speed only use enabled equivs
        return z
    # not one of the preferred types: Number / array ducktype
    return Quantity(z, cu.redshift).value
