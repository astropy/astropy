# Licensed under a 3-clause BSD style license - see LICENSE.rst
from .binned import BinnedTimeSeries as BinnedTimeSeries
from .core import (
    BaseTimeSeries as BaseTimeSeries,
    autocheck_required_columns as autocheck_required_columns,
)
from .downsample import aggregate_downsample as aggregate_downsample
from .periodograms import (
    BasePeriodogram as BasePeriodogram,
    BoxLeastSquares as BoxLeastSquares,
    BoxLeastSquaresResults as BoxLeastSquaresResults,
    LombScargle as LombScargle,
    LombScargleMultiband as LombScargleMultiband,
    base as base,
    bls as bls,
    lombscargle as lombscargle,
    lombscargle_multiband as lombscargle_multiband,
)
from .sampled import TimeSeries as TimeSeries
from . import (
    binned as binned,
    core as core,
    downsample as downsample,
    sampled as sampled,
    io as io,
    periodograms as periodograms,
    tests as tests,
)
