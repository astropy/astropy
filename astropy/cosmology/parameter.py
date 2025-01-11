# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Astropy cosmology parameter module.

.. deprecated:: 7.1

    This module is deprecated and will be removed in Astropy v8.
    All the public classes and functions have been and will be available
    in the :mod:`~astropy.cosmology` module.

"""

__all__ = ["Parameter"]  # noqa: F822


import sys
import warnings
from typing import Any


def __getattr__(name: str) -> Any:
    if name != "Parameter":
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}.")

    from ._src.parameter import Parameter

    setattr(sys.modules[__name__], name, Parameter)

    warnings.warn(
        "The module `astropy.cosmology.parameter` is deprecated since v7.1 and will be "
        "removed in a future version. Import from `astropy.cosmology` instead.",
        category=DeprecationWarning,
    )

    return Parameter
