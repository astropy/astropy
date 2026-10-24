# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Box Least Squares Algorithms.

AstroPy-compatible reference implementation of the transit periodogram used
to discover transiting exoplanets.
"""

__all__ = ["BoxLeastSquares", "BoxLeastSquaresResults"]

from .core import BoxLeastSquares, BoxLeastSquaresResults
