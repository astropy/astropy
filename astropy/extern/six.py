# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Handle loading six package from system or from the bundled copy
"""

import imp
from distutils.version import StrictVersion


_SIX_MIN_VERSION = StrictVersion('1.7.3')

# Update this to prevent Astropy from using its bundled copy of six
# (but only if some other version of at least _SIX_MIN_VERSION can
# be provided)
_SIX_SEARCH_PATH = ['astropy.extern.bundled.six', 'six']


def _find_module(name, path=None):
    """
    Alternative to `imp.find_module` that can also search in subpackages.
    """

    parts = name.split('.')

    for part in parts:
        if path is not None:
            path = [path]

        fh, path, descr = imp.find_module(part, path)

    return fh, path, descr


for mod_name in _SIX_SEARCH_PATH:
    try:
        mod_info = _find_module(mod_name)
    except ImportError:
        continue

    mod = imp.load_module(__name__, *mod_info)

    try:
        if StrictVersion(mod.__version__) >= _SIX_MIN_VERSION:
            break
    except (AttributeError, ValueError):
        # Attribute error if the six module isn't what it should be and doesn't
        # have a .__version__; ValueError if the version string exists but is
        # somehow bogus/unparseable
        continue
else:
    raise ImportError(
        "Astropy requires the 'six' module of minimum version {0}; "
        "normally this is bundled with the astropy package so if you get "
        "this warning consult the packager of your Astropy "
        "distribution.".format(_SIX_MIN_VERSION))
