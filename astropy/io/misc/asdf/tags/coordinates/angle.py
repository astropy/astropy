# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
from astropy.coordinates import Angle, Latitude, Longitude

from astropy.io.misc.asdf.tags.unit.quantity import QuantityType


__all__ = ['AngleType', 'LatitudeType', 'LongitudeType']


class AngleType(QuantityType):
    name = "coordinates/angle"
    types = [Angle]
    requires = ['astropy']
    version = "1.0.0"
    organization = 'astropy.org'
    standard = 'astropy'

    @classmethod
    def from_tree(cls, node, ctx):
        return Angle(super().from_tree(node, ctx))


class LatitudeType(AngleType):
    name = "coordinates/latitude"
    types = [Latitude]

    @classmethod
    def from_tree(cls, node, ctx):
        return Latitude(super().from_tree(node, ctx))


class LongitudeType(AngleType):
    name = "coordinates/longitude"
    types = [Longitude]

    @classmethod
    def from_tree(cls, node, ctx):
        wrap_angle = node['wrap_angle']
        return Longitude(super().from_tree(node, ctx), wrap_angle=wrap_angle)

    @classmethod
    def to_tree(cls, longitude, ctx):
        tree = super().to_tree(longitude, ctx)
        tree['wrap_angle'] = longitude.wrap_angle

        return tree
