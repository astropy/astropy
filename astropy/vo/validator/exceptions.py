# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Exceptions related to Virtual Observatory (VO) validation."""
from __future__ import absolute_import, division, print_function, unicode_literals


__all__ = ['BaseVOValidationError', 'ValidationMultiprocessingError']


class BaseVOValidationError(Exception):  # pragma: no cover
    """Base class for VO validation exceptions."""
    pass


class ValidationMultiprocessingError(BaseVOValidationError):  # pragma: no cover
    """Validation using multiprocessing failed."""
    pass


class InvalidValidationAttribute(BaseVOValidationError):  # pragma: no cover
    """Invalid validation attribute."""
    pass
