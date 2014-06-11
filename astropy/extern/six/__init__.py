
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
    # handle 'moves'
    _base_moves = 'six.moves'
    _load_six_moves(_base_moves, _dest_moves)
    sys.modules[_dest_root] = six
    six_system_package = True
except ImportError:
    pass

if not six_system_package:
    import astropy.extern.bundled.six as six
    # handle 'moves'
    _base_moves = 'astropy.extern.bundled.six.moves'
    _load_six_moves(_base_moves, _dest_moves)
    sys.modules[_dest_root] = six

