# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Support for the Python Array API standard.

Centralized all lookups of an array's namespace via `array_namespace`.
"""

from array_api_compat import array_namespace as _array_namespace
from array_api_compat import is_array_api_obj
from array_api_compat import numpy as _numpy_namespace

__all__ = ["array_namespace"]


def array_namespace(*values):
    """Return the Array API namespace of the arrays among ``values``.

    Values that are not arrays (Python scalars, lists, `None`, ...) are
    ignored. If none of the values is an array, the NumPy namespace is
    returned.

    Parameters
    ----------
    *values
        Arrays and other values.

    Returns
    -------
    namespace : module
        A namespace that follows the Array API standard (as provided by
        ``array-api-compat`` for libraries that do not fully follow it).

    Raises
    ------
    TypeError
        If the arrays belong to different, incompatible namespaces.
    """
    arrays = [value for value in values if is_array_api_obj(value)]
    if not arrays:
        return _numpy_namespace
    return _array_namespace(*arrays)
