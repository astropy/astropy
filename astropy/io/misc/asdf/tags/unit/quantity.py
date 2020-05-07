# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-

from astropy.units import Quantity

from asdf.tags.core import NDArrayType

from astropy.io.misc.asdf.types import AstropyAsdfType


class QuantityType(AstropyAsdfType):
    name = 'unit/quantity'
    types = ['astropy.units.Quantity']
    requires = ['astropy']
    version = '1.1.0'

    @classmethod
    def to_tree(cls, quantity, ctx):
        node = {}
        if isinstance(quantity, Quantity):
            node['value'] = quantity.value
            node['unit'] = quantity.unit
            return node
        raise TypeError(f"'{quantity}' is not a valid Quantity")

    @classmethod
    def from_tree(cls, node, ctx):
        if isinstance(node, Quantity):
            return node

        unit = node['unit']
        value = node['value']
        if isinstance(value, NDArrayType):
            value = value._make_array()
        return Quantity(value, unit=unit)
