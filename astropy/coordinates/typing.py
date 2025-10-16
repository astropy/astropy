"""Experimental module for supporting type annotations in :mod:`~astropy.coordinates`."""

__all__ = ["Matrix3x3"]


from typing import Literal, TypeAlias

import numpy as np

Matrix3x3: TypeAlias = np.ndarray[tuple[Literal[3], Literal[3]], np.dtype[np.floating]]
