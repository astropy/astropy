# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-

from asdf.yamlutil import custom_tree_to_tagged_tree
from astropy.coordinates import Angle, Latitude, Longitude

from ..unit.quantity import QuantityType


__all__ = ['AngleType', 'LatitudeType', 'LongitudeType']


class AngleType(QuantityType):
    name = "coords/angle"
    types = ['astropy.coordinates.angles.Angle']
    requires = ['astropy']
    version = "1.0.0"
    organization = 'astropy.org'
    standard = 'astropy'

    @classmethod
    def from_tree(cls, node, ctx):
        return Angle(super().from_tree(node, ctx))


class LatitudeType(AngleType):
    name = "coords/latitude"
    types = ['astropy.coordinates.angles.Latitude']

    @classmethod
    def from_tree(cls, node, ctx):
        return Latitude(super().from_tree(node, ctx))


class LongitudeType(AngleType):
    name = "coords/longitude"
    types = ['astropy.coordinates.angles.Longitude']

    @classmethod
    def from_tree(cls, node, ctx):
        wrap_angle = node['wrap_angle']
        return Longitude(super().from_tree(node, ctx), wrap_angle=wrap_angle)

    @classmethod
    def to_tree(cls, longitude, ctx):
        tree = super().to_tree(longitude, ctx)
        tree['wrap_angle'] = custom_tree_to_tagged_tree(longitude.wrap_angle, ctx)

        return tree
