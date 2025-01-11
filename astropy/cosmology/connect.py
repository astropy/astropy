# Licensed under a 3-clause BSD style license - see LICENSE.rst

__all__ = [  # noqa: F822
    "CosmologyFromFormat",
    "CosmologyRead",
    "CosmologyToFormat",
    "CosmologyWrite",
]

import sys
import warnings
from typing import Any


def __getattr__(name: str) -> Any:
    if name not in __all__ + ["convert_registry", "readwrite_registry"]:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}.")

    from ._src.io import connect

    obj = getattr(connect, name)

    setattr(sys.modules[__name__], name, obj)

    warnings.warn(
        "The module `astropy.cosmology.connect` is deprecated since v7.1 and will be "
        "removed in a future version. Import from `astropy.cosmology.io` instead.",
        category=DeprecationWarning,
    )

    return obj
