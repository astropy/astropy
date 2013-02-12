# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)


class Base(object):
    """
    The abstract base class of all unit formats.
    """

    def parse(self, s):
        """
        Convert a string to a unit object.  Must be overridden by the
        subclass.
        """

        raise NotImplementedError(
            "Can not parse {0}".format(self.__class__.__name__))

    def to_string(self, u):
        """
        Convert a unit object to a string.  Must be overridden by the
        subclass.
        """

        raise NotImplementedError(
            "Can not output in {0} format".format(self.__class__.__name__))
