# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Handle loading six package from system or from the bundled copy
"""

import imp


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


try:
    six_info = _find_module('six')
except ImportError:
    six_info = _find_module('astropy.extern.bundled.six')


imp.load_module(__name__, *six_info)
