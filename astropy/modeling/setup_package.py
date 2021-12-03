# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os
from os.path import join
from collections import defaultdict

from setuptools import Extension

from extension_helpers import import_file


# This defines the set of projection functions that we want to wrap.
# The key is the projection name, and the value is the number of
# parameters.

# (These are in the order that the appear in the WCS coordinate
# systems paper).
projections = {
    'azp': 2,
    'szp': 3,
    'tan': 0,
    'stg': 0,
    'sin': 2,
    'arc': 0,
    'zea': 0,
    'air': 1,
    'cyp': 2,
    'cea': 1,
    'mer': 0,
    'sfl': 0,
    'par': 0,
    'mol': 0,
    'ait': 0,
    'cop': 2,
    'coe': 2,
    'cod': 2,
    'coo': 2,
    'bon': 1,
    'pco': 0,
    'tsc': 0,
    'csc': 0,
    'qsc': 0,
    'hpx': 2,
    'xph': 0,
}
