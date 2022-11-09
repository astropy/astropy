# Licensed under a 3-clause BSD style license - see LICENSE.rst

# The LombScargle periodogram functionality has been moved to
# astropy.timeseries.periodograms.bls. The purpose of this file is to provide backward-
# compatibility during a transition phase. We can't emit a deprecation warning
# simply on import of this module, since the classes are imported into the
# top-level astropy.stats, so instead we wrap the main class and emit a
# warning during initialization.

import warnings

from astropy.timeseries.periodograms.lombscargle import (
    LombScargle as TimeseriesLombScargle,
)
from astropy.utils.exceptions import AstropyDeprecationWarning

__all__ = ["LombScargle"]


class LombScargle(TimeseriesLombScargle):
    """
    Compute the Lomb-Scargle Periodogram.

    This class has been deprecated and will be removed in a future version.
    Use `astropy.timeseries.LombScargle` instead.
    """

    def __init__(self, *args, **kwargs):
        warnings.warn(
            "Importing LombScargle from astropy.stats has been "
            "deprecated and will no longer be supported in future. "
            "Please import this class from the astropy.timeseries "
            "module instead",
            AstropyDeprecationWarning,
        )
        super().__init__(*args, **kwargs)
