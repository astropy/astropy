# Licensed under a 3-clause BSD style license - see LICENSE.rst

__all__: list[str] = []  # nothing is publicly scoped

import functools
from collections.abc import Callable
from numbers import Number
from typing import Any, Final, ParamSpec, Protocol, TypeAlias, TypeVar

import numpy as np
from numpy.typing import ArrayLike, NDArray

from astropy.units import Quantity
from astropy.utils.compat import COPY_IF_NEEDED

# isort: split
import astropy.cosmology._src.units as cu

from .signature_deprecations import _depr_kws_wrap

P = ParamSpec("P")
R = TypeVar("R")


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


# ===================================================================

ScalarTypes: TypeAlias = Number | np.generic
SCALAR_TYPES: Final = (float, int, np.generic, Number)  # arranged for speed


class HasShape(Protocol):
    shape: tuple[int, ...]


def aszarr(
    z: Quantity | NDArray[Any] | ArrayLike | ScalarTypes | HasShape, /
) -> NDArray[Any]:
    """Redshift as an Array duck type.

    Allows for any ndarray ducktype by checking for attribute "shape".
    """
    # Scalars
    if isinstance(z, SCALAR_TYPES):
        return np.asarray(z)

    # Quantities. We do this before checking for normal ndarray because Quantity is a
    # subclass of ndarray.
    elif isinstance(z, Quantity):
        return z.to_value(cu.redshift)[...]

    # Arrays
    elif isinstance(z, np.ndarray):
        return z

    return Quantity(z, cu.redshift, copy=COPY_IF_NEEDED, subok=True).view(np.ndarray)


# ===================================================================


def deprecated_keywords(
    *kws: str, since: str | float | tuple[str | float, ...]
) -> Callable[[Callable[P, R]], Callable[P, R]]:
    """Deprecate calling one or more arguments as keywords.

    Parameters
    ----------
    *kws: str
        Names of the arguments that will become positional-only.

    since : str or number or sequence of str or number
        The release at which the old argument became deprecated.
    """
    return functools.partial(_depr_kws, kws=kws, since=since)


def _depr_kws(
    func: Callable[P, R],
    /,
    kws: tuple[str, ...],
    since: str | float | tuple[str | float, ...],
) -> Callable[P, R]:
    wrapper = _depr_kws_wrap(func, kws, since)
    functools.update_wrapper(wrapper, func)
    return wrapper
