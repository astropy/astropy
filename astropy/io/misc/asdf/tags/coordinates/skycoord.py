# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
from astropy.coordinates import SkyCoord
from astropy.coordinates.tests.helper import skycoord_equal
from astropy.io.misc.asdf.types import AstropyType


class SkyCoordType(AstropyType):
    name = 'coordinates/skycoord'
    types = [SkyCoord]
    version = "1.0.0"

    @classmethod
    def to_tree(cls, obj, ctx):
        return obj.info._represent_as_dict()

    @classmethod
    def from_tree(cls, tree, ctx):
        return SkyCoord.info._construct_from_dict(tree)

    @classmethod
    def assert_equal(cls, old, new):
        assert skycoord_equal(old, new)
