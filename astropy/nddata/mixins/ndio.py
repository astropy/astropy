# Licensed under a 3-clause BSD style license - see LICENSE.rst
# This module implements the I/O mixin to the NDData class.

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from ...io import registry as io_registry

__all__ = ['NDIOMixin']


class NDIOMixin(object):
    """
    Mixin class to connect NDData to the astropy input/output registry.

    This mixin adds two methods to its subclasses, ``read`` and ``write``.
    """

    @classmethod
    def read(cls, *args, **kwargs):
        """
        Read and parse gridded N-dimensional data and return as an
        NDData-derived object.

        This function provides the NDDataBase interface to the astropy unified
        I/O layer.  This allows easily reading a file in the supported data
        formats.
        """
        return io_registry.read(cls, *args, **kwargs)

    def write(self, *args, **kwargs):
        """
        Write a gridded N-dimensional data object out in specified format.

        This function provides the NDDataBase interface to the astropy unified
        I/O layer.  This allows easily writing a file in the supported data
        formats.
        """
        io_registry.write(self, *args, **kwargs)
