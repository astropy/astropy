# Licensed under a 3-clause BSD style license - see LICENSE.rst
# Includes a backport of functools.lru_cache from
# http://code.activestate.com/recipes/578078/
# (C) Raymond Hettinger 2012, under the MIT license.
from __future__ import absolute_import

try:
    from functools import lru_cache
except ImportError:
    from ._functools.lru_cache import lru_cache
