# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import annotations

__all__: list[str] = []

import warnings

from astropy.utils.exceptions import AstropyDeprecationWarning


def __getattr__(name):
    """For backward compatibility (see #14982) allow imports from deprecated module.

    Raises
    ------
    AttributeError
        If "name" is not in :mod:`astropy.cosmology.io`
    """
    if name not in [
        "CosmologyRead",
        "CosmologyWrite",
        "readwrite_registry",
        "CosmologyFromFormat",
        "CosmologyToFormat",
        "convert_registry",
    ]:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}.")

    from astropy.cosmology import io

    obj = getattr(io, name)

    # Raise deprecation warning
    warnings.warn(
        f"astropy.cosmology.connect.{name} is deprecated (since v6.0)"
        " and will be removed in a future version. "
        f"Use astropy.cosmology.io.{name} instead.",
        AstropyDeprecationWarning,
        stacklevel=1,
    )

    return obj
