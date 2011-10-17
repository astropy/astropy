"""
This module contains a number of utilities for use during
setup/build/packaging that are useful to astropy as a whole.
"""

import os.path


def write_if_different(filename, data):
    """
    Write *data* to *filename*, if the content of the file is
    different.
    """
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
