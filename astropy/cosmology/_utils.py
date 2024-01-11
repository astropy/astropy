# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import annotations

__all__: list[str] = []  # nothing is publicly scoped

import functools
import operator
from dataclasses import Field
from numbers import Number
from typing import TYPE_CHECKING

import numpy as np

from astropy.units import Quantity

from . import units as cu

if TYPE_CHECKING:
    from typing import Any

    from astropy.cosmology import Parameter


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
        if getattr(z, "__module__", "").startswith("pandas"):
            # See https://github.com/astropy/astropy/issues/15576. Pandas does not play
            # well with others and will ignore unit-ful calculations so we need to
            # convert to it's underlying value.
            z = z.values
        if hasattr(z, "unit"):  # Quantity Column
            return (z << cu.redshift).value  # for speed only use enabled equivs
        return z
    # not one of the preferred types: Number / array ducktype
    return Quantity(z, cu.redshift).value


def all_cls_vars(obj: object | type, /) -> dict[str, Any]:
    """Return all variables in the whole class hierarchy."""
    cls = obj if isinstance(obj, type) else obj.__class__
    return functools.reduce(operator.__or__, map(vars, cls.mro()[::-1]))


def all_parameters(obj: object, /) -> dict[str, Field | Parameter]:
    """Get all fields of a dataclass, including those not-yet finalized.

    Parameters
    ----------
    obj : object | type
        A dataclass.

    Returns
    -------
    dict[str, Field | Parameter]
        All fields of the dataclass, including those not yet finalized in the class, if
        it's still under construction, e.g. in ``__init_subclass__``.
    """
    from astropy.cosmology.parameter import Parameter

    return {
        k: (v if isinstance(v, Parameter) else v.default)
        for k, v in all_cls_vars(obj).items()
        if (
            isinstance(v, Parameter)
            or (isinstance(v, Field) and isinstance(v.default, Parameter))
        )
    }
