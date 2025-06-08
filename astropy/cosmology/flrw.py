# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Astropy cosmology flrw module.

.. deprecated:: 7.1

    This module is deprecated and will be removed in a future version. All the public
    classes and functions have been and will continue to be available in the
    :mod:`~astropy.cosmology` module.

"""

__all__ = [  # noqa: RUF100, RUF022, F822
    # base
    "FLRW",
    "FlatFLRWMixin",
    # lambdacdm
    "LambdaCDM",
    "FlatLambdaCDM",
    # w0cdm
    "wCDM",
    "FlatwCDM",
    # w0wacdm
    "w0waCDM",
    "Flatw0waCDM",
    # w0wzcdm
    "w0wzCDM",
    "Flatw0wzCDM",
    # wpwazpcdm
    "wpwaCDM",
    "FlatwpwaCDM",
]

import sys
import warnings
from typing import Any


def __getattr__(name: str) -> Any:
    if name not in __all__:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}.")

    from ._src import flrw

    obj = getattr(flrw, name)

    setattr(sys.modules[__name__], name, obj)

    warnings.warn(
        "The module `astropy.cosmology.flrw` is deprecated since v7.1 and will be "
        "removed in a future version. Import from `astropy.cosmology` instead.",
        category=DeprecationWarning,
    )

    return obj
