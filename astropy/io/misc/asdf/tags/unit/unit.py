# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-

import six

from astropy.units import Unit, UnitBase
from astropy.io.misc.asdf.types import AstropyAsdfType


class UnitType(AstropyAsdfType):
    name = 'unit/unit'
    types = ['astropy.units.UnitBase']
    requires = ['astropy']

    @classmethod
    def to_tree(cls, node, ctx):
        if isinstance(node, six.string_types):
            node = Unit(node, format='vounit', parse_strict='warn')
        if isinstance(node, UnitBase):
            return node.to_string(format='vounit')
        raise TypeError(f"'{node}' is not a valid unit")

    @classmethod
    def from_tree(cls, node, ctx):
        return Unit(node, format='vounit', parse_strict='silent')
