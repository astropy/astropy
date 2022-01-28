# Licensed under a 3-clause BSD style license - see LICENSE.rst

from astropy.io.utils import load_all_entry_points

from .core import *
from .core import _json_base_encode
from . import builtins, numpy

# After importing, load all entry points
load_all_entry_points('astropy_io_json_extensions')
