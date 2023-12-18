# Licensed under a 3-clause BSD style license - see LICENSE.rst
from .core import (
    MetaAttribute as MetaAttribute,
    MetaData as MetaData,
)
from .exceptions import (
    MergeConflictError as MergeConflictError,
    MergeConflictWarning as MergeConflictWarning,
)
from ._merge import (
    MERGE_STRATEGIES as MERGE_STRATEGIES,
    MergeNpConcatenate as MergeNpConcatenate,
    MergePlus as MergePlus,
    MergeStrategy as MergeStrategy,
    MergeStrategyMeta as MergeStrategyMeta,
    enable_merge_strategies as enable_merge_strategies,
    merge as merge,
)
from .utils import common_dtype as common_dtype
