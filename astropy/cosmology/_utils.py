# Licensed under a 3-clause BSD style license - see LICENSE.rst

import functools
from numbers import Number

import numpy as np

from astropy.units import Quantity

from . import units as cu

__all__ = []  # nothing is publicly scoped


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
        """Wrapper converting arguments to numpy-compatible inputs.

        :func:`functools.wraps` of ``func`` where the first ``nin`` arguments are
        converted from |Quantity| to `numpy.ndarray` or scalar.
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


def aszarr(z):
    """Redshift as a `~numbers.Number` or |ndarray| / |Quantity| / |Column|.

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
