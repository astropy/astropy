# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Astropy cosmology core module.

.. deprecated:: 7.1

    This module is deprecated and will be removed in a future version. All the public
    classes and functions have been and will continue to be available in the
    :mod:`~astropy.cosmology` module.

"""

__all__ = ["Cosmology", "CosmologyError", "FlatCosmologyMixin"]  # noqa: F822


import sys
import warnings
from typing import Any


def __getattr__(name: str) -> Any:
    if name not in __all__ + ["_COSMOLOGY_CLASSES", "dataclass_decorator"]:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}.")

    from ._src import core

    obj = getattr(core, name)

    setattr(sys.modules[__name__], name, obj)

    warnings.warn(
        "The module `astropy.cosmology.core` is deprecated since v7.1 and will be "
        "removed in a future version. Import from `astropy.cosmology` instead.",
        category=DeprecationWarning,
    )

    return obj
