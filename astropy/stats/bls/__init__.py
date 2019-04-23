# Licensed under a 3-clause BSD style license - see LICENSE.rst

# The BoxLeastSquares periodogram functionality has been moved to
# astropy.timeseries.periodograms.bls. The purpose of this file is to provide backward-
# compatibility during a transition phase. We can't emit a deprecation warning
# simply on import of this module, since the classes are imported into the
# top-level astropy.stats, so instead we wrap the main class and emit a
# warning during initialization.

import warnings

from astropy.timeseries.periodograms.bls import (BoxLeastSquares as TimeseriesBoxLeastSquares,
                                                 BoxLeastSquaresResults as TimeseriesBoxLeastSquaresResults)
from astropy.utils.exceptions import AstropyDeprecationWarning

__all__ = ['BoxLeastSquares', 'BoxLeastSquaresResults']


class BoxLeastSquares(TimeseriesBoxLeastSquares):
    """
    Compute the box least squares periodogram.

    This class has been deprecated and will be removed in a future version.
    Use `astropy.timeseries.BoxLeastSquares` instead.
    """

    def __init__(self, *args, **kwargs):
        warnings.warn('Importing BoxLeastSquares from astropy.stats has been '
                      'deprecated and will no longer be supported in future. '
                      'Please import this class from the astropy.timeseries '
                      'module instead', AstropyDeprecationWarning)
        super().__init__(*args, **kwargs)


class BoxLeastSquaresResults(TimeseriesBoxLeastSquaresResults):
    """
    The results of a BoxLeastSquares search.

    This class has been deprecated and will be removed in a future version.
    Use `astropy.timeseries.BoxLeastSquaresResults` instead.
    """

    def __init__(self, *args, **kwargs):
        warnings.warn('Importing BoxLeastSquaresResults from astropy.stats has been '
                      'deprecated and will no longer be supported in future. '
                      'Please import this class from the astropy.timeseries '
                      'module instead', AstropyDeprecationWarning)
        super().__init__(*args, **kwargs)
