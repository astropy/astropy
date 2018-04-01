# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-

from numpy import isscalar
from astropy.units import Quantity

from asdf.yamlutil import custom_tree_to_tagged_tree
from asdf.tags.core import NDArrayType

from ...types import AstropyAsdfType
from .unit import UnitType


class QuantityType(AstropyAsdfType):
    name = 'unit/quantity'
    types = ['astropy.units.Quantity']
    requires = ['astropy']
    version = '1.1.0'

    @classmethod
    def to_tree(cls, quantity, ctx):
        node = {}
        if isinstance(quantity, Quantity):
            node['value'] = custom_tree_to_tagged_tree(quantity.value, ctx)
            node['unit'] = custom_tree_to_tagged_tree(quantity.unit, ctx)
            return node
        raise TypeError("'{0}' is not a valid Quantity".format(quantity))

    @classmethod
    def from_tree(cls, node, ctx):
        if isinstance(node, Quantity):
            return node

        unit = UnitType.from_tree(node['unit'], ctx)
        value = node['value']
        if isinstance(value, NDArrayType):
            value = value._make_array()
        return Quantity(value, unit=unit)
