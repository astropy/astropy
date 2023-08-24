# Licensed under a 3-clause BSD style license - see LICENSE.rst

__all__ = []  # nothing is publicly scoped

import sys

from astropy.utils.decorators import deprecated


# TODO: complete the deprecation for v6.2 / v7
def __getattr__(name):
    """Get realizations using lazy import from ``PEP 562``.

    Raises
    ------
    AttributeError
        If "name" is not in :mod:`astropy.cosmology.realizations`
    """
    if name not in ("vectorize_redshift_method", "aszarr"):
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}.")

    from . import _utils

    func = deprecated(
        since="v6.0",
        message=(
            "this private function has been moved to the private module"
            " `astropy.cosmology._utils`"
        ),
    )(getattr(_utils, name))

    sys.modules[__name__].__dict__[name] = func  # cache for next time

    return func
