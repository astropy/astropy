# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Metadata exceptions and warnings."""

from astropy.utils.exceptions import AstropyWarning

__all__ = ["MergeConflictError", "MergeConflictWarning"]


class MergeConflictError(TypeError):
    pass


class MergeConflictWarning(AstropyWarning):
    pass
