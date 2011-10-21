"""
This module contains a number of utilities for use during
setup/build/packaging that are useful to astropy as a whole.
"""

import os.path
import sys


def write_if_different(filename, data):
    """
    Write *data* to *filename*, if the content of the file is
    different.
    """
    assert isinstance(data, bytes)

    if os.path.exists(filename):
        with open(filename, 'rb') as fd:
            original_data = fd.read()
    else:
        original_data = None

    if original_data != data:
        with open(filename, 'wb') as fd:
            fd.write(data)


def check_numpy():
    """
    Check that Numpy is installed and it is of the minimum version we
    require.
    """
    import numpy

    major, minor, rest = numpy.__version__.split(".", 2)
    if (int(major), int(minor)) < (1, 3):
        raise ImportError("numpy version 1.3 or later must be installed to build astropy")


def get_numpy_include_path():
    """
    Gets the path to the numpy headers
    """
    import numpy

    try:
        numpy_include = numpy.get_include()
    except AttributeError:
        numpy_include = numpy.get_numpy_include()
    return numpy_include



################################################################################
# Backport of importlib.import_module from 3.x.  This backport was provided by
# Brett Cannon and downloaded from here:
#    http://pypi.python.org/pypi/importlib/1.0.1
# While not critical (and in no way guaranteed!), it would be nice to keep this
# code compatible with Python 2.3.

try:
    from importlib import import_module
except ImportError:

    def _resolve_name(name, package, level):
        """Return the absolute name of the module to be imported."""
        if not hasattr(package, 'rindex'):
            raise ValueError("'package' not set to a string")
        dot = len(package)
        for x in xrange(level, 1, -1):
            try:
                dot = package.rindex('.', 0, dot)
            except ValueError:
                raise ValueError("attempted relative import beyond top-level "
                                  "package")
        return "%s.%s" % (package[:dot], name)


    def import_module(name, package=None):
        """Import a module.

        The 'package' argument is required when performing a relative
        import. It specifies the package to use as the anchor point
        from which to resolve the relative import to an absolute
        import.
        """
        if name.startswith('.'):
            if not package:
                raise TypeError(
                    "relative imports require the 'package' argument")
            level = 0
            for character in name:
                if character != '.':
                    break
                level += 1
            name = _resolve_name(name[level:], package, level)
        __import__(name)
        return sys.modules[name]
