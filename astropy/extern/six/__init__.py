# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Handle loading six package from system or from the bundled copy

"""
import sys

# Trying to load alternate six packages
sys.modules['astropy.extern.six'] = None

# We have removed everything we already imported
# Importing again

import sys 

def _load_six_moves(base, dest):
    _cur_sys_modules = list(sys.modules.items())
    for i,mod in _cur_sys_modules:
        if i.startswith(base):
            pre, full, trail = i.partition(base)
            if not pre:
                modname = dest + trail
                sys.modules[modname] = mod

six_system_package = False

_dest_moves = 'astropy.extern.six.moves'
_dest_root = 'astropy.extern.six'

try:
    import six
    _system_package = True
except ImportError:
    _system_package = False

if _system_package:
    # Check six version
    _valid_version = True
    if _valid_version:
        # handle 'moves'
        _base_moves = 'six.moves'
        _load_six_moves(_base_moves, _dest_moves)
        six.system_package = True
        six.bundled_package = False
        sys.modules[_dest_root] = six

if not _system_package:
    import astropy.extern.bundled.six as six
    # handle 'moves'
    _base_moves = 'astropy.extern.bundled.six.moves'
    _load_six_moves(_base_moves, _dest_moves)
    six.system_package = False
    six.bundled_package = True
    sys.modules[_dest_root] = six

