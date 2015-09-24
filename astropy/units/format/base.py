# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)


from ...utils.misc import InheritDocstrings
from ...extern import six


class _FormatterMeta(InheritDocstrings):
    registry = {}

    def __new__(mcls, name, bases, members):
        if 'name' in members:
            formatter_name = members['name'].lower()
        else:
            formatter_name = members['name'] = name.lower()

        cls = super(mcls, _FormatterMeta).__new__(mcls, name, bases, members)

        mcls.registry[formatter_name] = cls

        return cls


@six.add_metaclass(_FormatterMeta)
class Base(object):
    """
    The abstract base class of all unit formats.
    """

    def __new__(cls, *args, **kwargs):
        # This __new__ is to make it clear that there is no reason to
        # instantiate a Formatter--if you try to you'll just get back the
        # class
        return cls

    @classmethod
    def parse(cls, s):
        """
        Convert a string to a unit object.
        """

        raise NotImplementedError(
            "Can not parse {0}".format(cls.__name__))

    @classmethod
    def to_string(cls, u):
        """
        Convert a unit object to a string.
        """

        raise NotImplementedError(
            "Can not output in {0} format".format(cls.__name__))
