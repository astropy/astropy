# Licensed under a 3-clause BSD style license - see LICENSE.rst
from .base import BaseWCSWrapper as BaseWCSWrapper
from .sliced_wcs import (
    SlicedLowLevelWCS as SlicedLowLevelWCS,
    sanitize_slices as sanitize_slices,
)
from . import (
    base as base,
    sliced_wcs as sliced_wcs,
    tests as tests,
)
