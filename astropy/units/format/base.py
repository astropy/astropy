# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os

from ...utils.misc import InheritDocstrings


class _FormatterMeta(InheritDocstrings):
    registry = {}

    def __new__(mcls, name, bases, members):
        if 'name' in members:
            formatter_name = members['name'].lower()
        else:
            formatter_name = members['name'] = name.lower()

        cls = super().__new__(mcls, name, bases, members)

        mcls.registry[formatter_name] = cls

        return cls


TAB_HEADER = """# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

# This file was automatically generated from ply. To re-generate this file,
# remove it from this folder, then build astropy and run the tests in-place:
#
#   python setup.py build_ext --inplace
#   pytest astropy/units
#
# You can then commit the changes to this file.

"""


class Base(metaclass=_FormatterMeta):
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

    @classmethod
    def _add_tab_header(cls, name):

        lextab_file = os.path.join(os.path.dirname(__file__), name + '.py')

        with open(lextab_file, 'r') as f:
            contents = f.read()

        with open(lextab_file, 'w') as f:
            f.write(TAB_HEADER)
            f.write(contents)
