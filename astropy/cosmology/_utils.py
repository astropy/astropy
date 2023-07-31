# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import annotations

import functools
import inspect
import operator
from collections import OrderedDict
from copy import deepcopy
from dataclasses import dataclass
from numbers import Number
from typing import Any, Mapping

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


# TODO: "upstream" this to `astropy.utils.metadata`?
@dataclass(frozen=True)
class MetaData:
    """Metadata for a cosmology."""

    default: Mapping[str, Any] | None = None
    copy: bool = True

    def __get__(self, instance, owner):
        if instance is None:
            return self.default
        return instance._meta

    def __set__(self, instance, value: Mapping[str, Any] | None):
        if not hasattr(instance, "_meta"):
            object.__setattr__(instance, "_meta", OrderedDict())
        instance._meta.update(
            (deepcopy(value) if self.copy else value) if value is not None else {}
        )


def all_fields(cls):
    """Get all fields of a dataclass, including in ``__init_subclass__``."""
    return functools.reduce(
        operator.__or__,
        (getattr(c, "__dataclass_fields__", {}) for c in cls.mro()[::-1] + [cls]),
    )


def all_cls_vars(cls):
    """Return all variables in the whole class hierarchy."""
    cls = cls if isinstance(cls, type) else cls.__class__
    return functools.reduce(operator.__or__, map(vars, cls.mro()[::-1] + [cls]))


def _init_signature(cls):
    """Initialization signature (without 'self')."""
    # get signature, dropping "self" by taking arguments [1:]
    sig = inspect.signature(cls.__init__)
    sig = sig.replace(parameters=list(sig.parameters.values())[1:])
    return sig
