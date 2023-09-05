# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""This module contains helper functions and classes for handling metadata."""

from .core import MetaAttribute, MetaData
from .exceptions import MergeConflictError, MergeConflictWarning
from .merge import (
    MERGE_STRATEGIES,
    MergeNpConcatenate,
    MergePlus,
    MergeStrategy,
    MergeStrategyMeta,
    enable_merge_strategies,
    merge,
)
from .utils import common_dtype

__all__ = [
    # core
    "MetaData",
    "MetaAttribute",
    # exceptions
    "MergeConflictError",
    "MergeConflictWarning",
    # merge
    "MERGE_STRATEGIES",
    "MergeStrategyMeta",
    "MergeStrategy",
    "MergePlus",
    "MergeNpConcatenate",
    "enable_merge_strategies",
    "merge",
    # utils
    "common_dtype",
]
