# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
from astropy.coordinates import EarthLocation
from astropy.io.misc.asdf.types import AstropyType


class EarthLocationType(AstropyType):
    name = 'coordinates/earthlocation'
    types = [EarthLocation]
    version = '1.0.0'

    @classmethod
    def to_tree(cls, obj, ctx):
        return obj.info._represent_as_dict()

    @classmethod
    def from_tree(cls, node, ctx):
        return EarthLocation.info._construct_from_dict(node)

    @classmethod
    def assert_equal(cls, old, new):
        return (old == new).all()
